"""
   CONFORMALIZED_MEAN_CURVATURE_FLOW Flow a surface according to "Can mean
   curvature flow be made non-singular?" [Kazhdan et al. 2012]
  
   Inputs:
     V  #V by dim list of vertex positions
     F  #F by 3 list of triangle indices into V
     Optional:
       "MaxIter" followed by maximum number of iterations {100}
       "Conformalize" followed by whether to rebuild _just_ the mass
         matrix each step {true}
       "delta" followed by delta value, should roughly be in range
         [1e-13,1e13] {1}
       "LaplacianType" followed by "cotangent" of "uniform".
       "V0" followed by #V by 3 mesh positions to treat as initial mesh (to
         build laplacian from)
       "RescaleOutput" followed by whether to scale output to match input
         (otherwise scaled to unit surface area and moved to origin for
         numerical robustness) {false}
   Outputs:
     U  #V by dim list of new vertex positions
     Usteps  #V by dim by iterations list of vertex positions during flow
"""
function conformalized_mean_curvature_flow{T,T2}(V::Array{T},F::Array{T2} ;
  Delta=1, MaxIter=100, LaplacianType= "cotangent", UntilSelfIntersectionFree=false,
  V0 = [], RescaleOutput=false, MinDiff=1e-15, Conformalize=true, SelfIntersect = nothing)

  U = []
  Usteps = []
  nargout = 3

  # default values
  delta = Delta
  max_iter = MaxIter
  laplacian_type = LaplacianType
  until_self_intersection_free = UntilSelfIntersectionFree
  V0 = isempty(V0) ? V : V0;
  rescale_output = RescaleOutput
  min_diff = MinDiff
  conformalize = Conformalize
  selfintersect = SelfIntersect
  
  # Map of parameter names to variable names

  function laplacian(Vl,Fl)
    if laplacian_type == "cotangent"
      Lap = cotmatrix(Vl,Fl);
    elseif laplacian_type == "uniform"
      A = adjacency_matrix(Fl);
      Lap = A - diag(sparse(sum(A,2)));
    end
    Lap
  end

  L = laplacian(V,F);

  if nargout > 1
    Usteps = zeros(size(V)..., max_iter);
  end
  U = V;
  iter = 1;
  while true
    if nargout > 1
      Usteps[:,:,iter] = U;
    end
    if until_self_intersection_free
      (_,_,IF) = selfintersect(U,F,DetectOnly=true,FirstOnly=true);
      if isempty(IF)
        break
      end
    end
    U_prev = U;
    # "full" seems slight more stable than "barycentric" which is more stable
    # than "voronoi"
    M = massmatrix(U,F,"barycentric");
    if ~conformalize
      L = laplacian(V,F)
    end
    U = (M-delta*L)\(M*U);
    area = sum(doublearea(U,F)*0.5);
    c = sum(broadcast(.*,0.5*doublearea(U,F)/area,barycenter(U,F)),1);
    U = broadcast(.-,U,c);
    U = U/sqrt(area);
    # Use difference from previous as stopping criterion
    # Better would be to look for convergence while factoring out M??bius
    # transformation.
    # Q: Stop when no change in angles?
    if ~until_self_intersection_free 
      d = trace(((U-U_prev)'*M*(U-U_prev)).^2);
      if d < min_diff
        println("converged...");
        break;
      end
    else
      #!me
      if (iter % 10) == 1
        println(trace(((U-U_prev)'*M*(U-U_prev)).^2))
      end
    end

    if iter >= max_iter
      println("Max iterations (", max_iter, ") exceeded without convergence: ",trace(((U-U_prev)'*M*(U-U_prev)).^2));
      break;
    end
    iter = iter + 1;
  end

  if nargout > 1
    Usteps = Usteps[:,:,1:iter];
  end

  if rescale_output
    area = sum(doublearea(V,F)*0.5);
    c = sum(broadcast(.*,0.5*doublearea(V,F)/area,barycenter(V,F)),1);
    U = U*sqrt(area);
    U = broadcast(.+,U,c);
    if nargout > 1
      for iter in 2:size(Usteps,3)
        Usteps[:,:,iter] = Usteps[:,:,iter]*sqrt(area);
        Usteps[:,:,iter] = broadcast(.+,Usteps[:,:,iter],c);
      end
    end
  end

  (U,Usteps)
end
