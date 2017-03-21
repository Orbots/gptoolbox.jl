"""
  % NORMROW  Compute l2 row vector norms
  %
  % B = normrow( A )
  %
  % Input:
  %  A  #A by D list of row vectors of dimension D
  % Output:
  %  B  #A list of norms of row vectors in A
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch), Daniele Panozzo
  %
"""
function normrow( A )

  if size(A,2) == 2 
    hypot(A[:,1],A[:,2]);
  else
    sqrt(sum(A.^2,2));
  end
end

