"""
  # FIT_ROTATIONS Given an input mesh and new positions find rotations for
  # every vertex that best maps its one ring to the new one ring
  # 
  # R = fit_rotations(S,'ParamName',ParamValue)
  #
  # Inputs:
  #   S  dim by dim by #rotations list of covariance matrices to fit rotations
  #     to
  #   Optional parameters
  #     'AllowFlips'  optionally followed by true or false, find best fitting
  #       rotation OR reflection
  #     'SinglePrecision'  use single precision if available.
  # Outputs:
  #   R  dim by dim by #F list of rotations
  #   SS  dim by dim by #rotations list of svd diagonals
  #
"""
function fit_rotations(S;allow_flips=false,single_precision=true)

  """
    # DET3 compute the determinant of a 3x3 matrix
    # Input:
    #   M  3 by 3 matrix
    # Output:
    #   D  determinant
    #
    # http://en.wikipedia.org/wiki/Determinant#3-by-3_matrices
  """
  # julia has a det function.  check speed.  ddmath has fast one I think
  function det3(M)
      M[1,1] * M[2,2] * M[3,3] + 
      M[1,2] * M[2,3] * M[3,1] + 
      M[1,3] * M[2,1] * M[3,2] - 
      M[1,3] * M[2,2] * M[3,1] - 
      M[1,2] * M[2,1] * M[3,3] - 
      M[1,1] * M[2,3] * M[3,2]
  end

  # Even faster way to check if mex exists
  fit_rotations_mex = false
  if fit_rotations_mex
    nr = size(S,3);
    dim = size(S,1);
    SS = reshape(permute(S,[3 1 2]),[nr*dim dim]);
    R = fit_rotations_mex(SS,varargin{:});
    R = reshape(R,[dim dim nr]);
    return;
  end

  dim = size(S,1);
  assert(dim == size(S,2));
  nr = size(S,3);

  #R = cellfun(@fit_rotation,S,'UniformOutput',false);
  #R = reshape(cell2mat(R'),[dim dim nr]);
  # For loop is faster
  R = zeros(dim, dim, nr);
  SS = zeros(dim, dim, nr);
  for ii = 1:nr
    # svd 
    (su,ss,sv)=svd(S[:,:,ii]);
    Ri = sv*su';
    SS[:,:,ii] = diagm(ss)
    # if reflection then flip last column
    if(~allow_flips && det(Ri) < 0 )
      su[:,end] = -su[:,end];
      Ri = sv*su';
    end
    # should definitely be rotation now
    #assert( det(Ri) >= 0 );
    R[:,:,ii] = Ri;
  end

  (R,SS)
end
