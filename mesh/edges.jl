using Combinatorics
"""
   EDGES Compute the unique undireced edges of a simplicial complex
   
   E = edges(F)
  
   Input:
    F #F x simplex-size  matrix of indices of simplex corners
   Output:
    E edges in sorted order, direction of each is also sorted
  
   Example:
     # get unique undirected edges
     E = edges(F);
     # get unique directed edges
     E = [E ; E(:,2) E(:, 1)];
   
   Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
"""
function edges(F)

  # all combinations of edges
  n = size(F,2)

  eges = reduce(hcat, permutations(1:n,2))'
  m = reduce(max,F[:]) 
  A = sparse(vec(F[:,eges[:,1]]),vec(F[:,eges[:,2]]),1.0,m,m)
  (EI,EJ) = findn(tril(A+A'))
  E = [EJ EI]
end

