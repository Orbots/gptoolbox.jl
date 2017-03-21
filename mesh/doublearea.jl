function doublearea(V,F)
  # DOUBLEAREA Compute the double area of the triangles of a mesh
  #
  # dblA = doublearea(V,F)
  #
  # Inputs:
  #  V #V x dim matrix of vertex coordinates
  #  F #F x 3  matrix of indices of triangle corners
  # Outputs:
  #  dblA   #F list of twice the area of each corresponding face. For dim = 2
  #    this is a signed quantity
  #
  # Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch), and Daniele Panozzo
  # Julia conversion: Michael Alexander Ewert, Digital Domain 3.0, 2017
  #
  dblA = 0.0

  # "Lecture Notes on Geometric Robustness" Shewchuck 09, Section 3.1
  # http://www.cs.berkeley.edu/~jrs/meshpapers/robnotes.pdf
  s_which = size(V,2)
  if s_which == 2
    r = V[F[:,1],:]-V[F[:,3],:];
    s = V[F[:,2],:]-V[F[:,3],:];
    dblA = r[:,1].*s[:,2] - r[:,2].*s[:,1];
  elseif s_which == 3
    dblA = sqrt( 
      doublearea(V[:,[2, 3]],F).^2 + 
      doublearea(V[:,[3, 1]],F).^2 + 
      doublearea(V[:,[1, 2]],F).^2);
  else
    # For arbitrary dimension use Kahan's heron's formula
    # triangles
    # edge lengths numbered same as opposite vertices
    #!me we don't care about this. no hyperspace meshes please
    throw("hyperspace is not real")
    #==
    l = [ 
      sqrt(sum((V[F[:,2],:]-V[F[:,3],:]).^2,2)) 
      sqrt(sum((V[F[:,3],:]-V[F[:,1],:]).^2,2)) 
      sqrt(sum((V[F[:,1],:]-V[F[:,2],:]).^2,2)) 
      ];
    dblA = doublearea_intrinsic(l);==#
  end
  dblA

end
