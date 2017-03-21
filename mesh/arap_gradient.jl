"""
  # ARAP_GRADIENT Compute the Gradient of an as-rigid-as-possible for a mesh in
  # rest configuration (V,F) with vertices now at U, according to "A Simple
  # Geometric Model for Elastic Deformations" [Chao et al. 2010]
  #
  # G = arap_gradient(V,F,U)
  # [G,E,R] = arap_gradient(V,F,U,"ParameterName",ParameterValue, ...)
  #
  # Inputs:
  #   V  #V by dim rest vertex positions
  #   F  #F by simplex-size simplex indices into V
  #   U  #V by dim deformed vertex positions
  #     Optional:
  #       "Rotations" followed by R  dim by dim by #R list of best fit
  #         rotations
  #       "Energy" followed by either "spokes","spokes-and-rims","elements"
  # Outputs:
  #   G  #V by dim list of gradient vectors per vertex
  #   E  arap energy
  #   R  dim by dim by #R list of best fit rotations
  #   data 
  #
  # See also: arap, arap_hessian
  #
  # Example:
  # % Given a mesh (V,F) and deformed positions U0, flow to energy minimum
  # % using Newton"s method.
  # clf;
  # hold on;
  #   tsurf(F,V,fphong,"FaceColor","r","SpecularStrength",0,"AmbientStrength",0.5);
  #   t = tsurf(F,U0,fphong,"FaceColor","b","SpecularStrength",0,"AmbientStrength",0.5);
  # hold off;
  # axis equal;
  # view(2);
  # camlight;
  # U = U0;
  # delta_t = 1e-1;
  # while true
  #   [G,E] = arap_gradient(V,F,U);
  #   U = U - delta_t * G;
  #   U = bsxfun(@plus,U,mean(V)-mean(U)+[max(V(:,1))-min(V(:,1)) 0 0]);
  #   t.Vertices = U;
  #   title(sprintf("E = %g\n",E),"FontSize",20);
  #   drawnow;
  # end
  #
"""
type ArapData
  L::SparseMatrixCSC{Float64,Int64}
  CSM::SparseMatrixCSC{Float64,Int64}
  K::SparseMatrixCSC{Float64,Int64}
end

ArapData() = ArapData( sparse(Array{Float64,2}()),sparse(Array{Float64,2}()),sparse(Array{Float64,2}()) )

function arap_gradient(V::Array{Float64}, F::Array{Int64}, U::Array{Float64};
                       energy="elements",R=[],data::ArapData=ArapData(),single_precision=true)
  nargout = 4 
  G = []
  E = []

  if size(F,2) == 4
    energy = "elements"
  elseif size(F,2) == 3 
    energy = "spokes-and-rims"
  end

  # TODO: implement "flat" arap like `arap.m`, for now use placeholders with
  # non-flat defaults.
  ref_V = V;
  ref_F = F;
  dim = size(ref_V,2);
  flat = false;

  if isempty(data.L)
    ss = size(F,2);
    data.L = cotmatrix(V,F);
    if isempty(R)
      data.CSM = covariance_scatter_matrix(ref_V,ref_F,energy=energy)
    end
    (_,data.K) = arap_rhs(ref_V,ref_F,Array{Float64,3}(0,0,0),energy=energy);
  end

  # compute covariance matrix elements
  S = zeros(size(data.CSM,1),dim);
  S[:,1:dim] = data.CSM*repmat(U,dim,1);
  # dim by dim by n list of covariance matrices
  SS = permutedims(reshape(S,(Int64(size(data.CSM,1)/dim), dim, dim)),[2 3 1]);
  # fit rotations to each deformed vertex
  (R,SS) = fit_rotations(SS,single_precision=single_precision);

  nr = size(R,3);
  Rcol = reshape(permutedims(R,[3 1 2]),nr*dim*dim,1);
  dV = data.K * Rcol;
  dV = reshape(dV,(size(V,1), dim));

  G = -(data.L*U + dV);
  if nargout > 1
    E = trace(-U'*0.5*data.L*U - U'*dV - V'*0.5*data.L*V);
  end

  (G,E,R,data)
end

