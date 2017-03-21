Base.max{T}(A::Array{T}) = mapslices(x->max(x...),A,1)
Base.isless{T,U}(A::Array{T},r::U) = reduce((a,b)-> a && (b < r), true, A)
"""
   matlab return [VV,Uforward,Ubackward]
   Julia conversion: Michael Alexander Ewert, Digital Domain 3.0, 2017
   UNTANGLE Given a possible self-intersecting mesh (V,F) _untangle_
   self-intersections by running a (conformalized) mean curvature flow
   (guaranteed to resolve self-intersections if converges for sphere topology
   surfaces) and then re-inflating with collision response. This is most
   similar to the method described in "Consistent Volumetric Discretizations
   Inside Self-Intersecting Surfaces" [Sacht et al. 2013], which itself is an
   expansion of the idea in "Interference-Aware Geometric Modeling" [Harmon et
   al. 2011].
  
   VV = untangle(V,F)
   
   Inputs:
     V  #V by dim list of vertex positions
     F  #F by 3 list of triangle indices into V
     Optional:
       "ExpansionEnergy" Followed by one of the following energy types:
         {"none"}  just let eltopo find a valid state and move on
         "arap"  gradient descent on ARAP energy after each reverse step.
       "FinalEnergy" Followed by one of the following energy types:
         {"none"}  just let eltopo find a valid state and move on
         "arap"  gradient descent on ARAP energy after each reverse step.
       "Delta" followed by delta value, should roughly be in range
         [1e-13,1e13] {1}
       "FinalMaxIter" followed by maximum number of iterations {100}
       "ExpansionMaxIter" followed by maximum number of iterations {10}
       "MaxEnergyIter" followed by maximum number of iterations for energy
         minimization {100}
       "LaplacianType" followed by {"cotangent"} of "uniform".
   Outputs:
     VV  #V by dim list of new vertex positions
     U  #V by dim by steps list of steps
"""
function untangle(V::Array{Float64}, F::Array{Int64}, selfintersect ; 
  delta = 3e-5,
  max_iter = 100,
  laplacian_type = "cotangent",
  f_max_eit = 100,
  e_max_eit = 10,
  expansion_energy = "none",
  final_energy = "none" )
#==
  delta = 3e-5
  max_iter = 100
  laplacian_type = "cotangent"
  f_max_eit = 100
  e_max_eit = 10
  expansion_energy = "none"
  final_energy = "none" 
  ==#

print("BEGIN untangle>")
  # matlab fakage
  Ubackward = []
  Uforward = []
  nargout = 2
  VV = []

  function append_Ubackward()
    if main_nargout>=3 
      if isempty(Ubackward)
        Ubackward = zeros([size(VV) 1+(size(Uforward,3)-1)*e_max_eit+f_max_eit]);
      end
      Ubackward[:,:,bit] = VV;
      bit = bit+1;
    end
  end

  V0 = V;
  rescale_output = false;

  if expansion_energy == "none"
    e_max_eit = 1;
  end
  if final_energy == "none"
    f_max_eit = 1;
  end

print("BEGIN cMCF>")

  (flowed,Uforward) = conformalized_mean_curvature_flow( 
    V,F,Delta = delta, MaxIter = max_iter, 
    LaplacianType = laplacian_type, 
    UntilSelfIntersectionFree = true, 
    SelfIntersect = selfintersect,
    RescaleOutput = true)
println("END cMCF.")

  Ubackward = [];
  bit = 1;
  main_nargout = nargout;
  VV = Uforward[:,:,end];
  append_Ubackward();

  (_,_,IF) = selfintersect(VV,F,DetectOnly = true, FirstOnly = true);

  if ~isempty(IF)
    println("cMCF didn't untangle mesh")
  else
    for it = size(Uforward,3)-1:-1:0
      if it == 0
        energy  = final_energy;
        max_eit = f_max_eit;
        it = 1;
      else
        energy  = expansion_energy;
        max_eit = e_max_eit;
      end
      if energy == "arap-local-global"
        # This energy option also seems useless
        VV_prev = VV;
        tol = 1e-7;
        arap_data = [];
        for eit = 1:max_eit
          # find a best fit rigid transformation for each component
          (G,E_prev,R) = arap_gradient(Uforward[:,:,it],F,VV_prev);
          (VV,arap_data) = arap(Uforward[:,:,it],F,[],[],"V0",VV,"Dynamic",zeros(size(V)),"MaxIter",1,"Data",arap_data);
          (VV,_) = eltopo(VV_prev,F,VV);
          append_Ubackward();
          C = normrow(arap_gradient(Uforward[:,:,it],F,VV));
          dVV = max(abs(VV-VV_prev));
          dVV
          if dVV < tol
            break;
          end
        end
      elseif energy == "arap"
print("BEGIN arap>")
        tol = 1e-7;
        VV_prev = VV;
        delta_t = 0.3;
        for eit = 1:max_eit
          (G,E_prev,R) = arap_gradient(Uforward[:,:,it],F,VV_prev);
          dVV = G;
          while true
            VV = VV_prev - delta_t * dVV;
            (_,E) = arap_gradient(Uforward[:,:,it],F,VV);
            if E<E_prev
              break;
            end
            delta_t = delta_t*0.9;
            assert(delta_t>0);
          end
          (VV,_) = eltopo(VV_prev,F,VV);
          append_Ubackward();
          C = normrow(G);
          dVV = max(abs(VV-VV_prev));
          if dVV < tol
            break;
          end
          (_,E) = arap_gradient(Uforward[:,:,it],F,VV);
          if E_prev < E
            delta_t = 0.5*delta_t;
          end
          VV_prev = VV;
        end
println("END arap.")
      elseif energy == "none"
        (VV,_) = eltopo(VV,F,Uforward[:,:,it]);
        append_Ubackward();
      else 
        throw(string("Unknown energy type ", energy));
      end

      (_,_,IF) = selfintersect(VV,F,DetectOnly = true, FirstOnly = true);
      if ~isempty(IF)
        println("eltopo failed");
      end
    end
  end

println("END untangle.")
  (VV,Uforward,Ubackward, flowed)
end
