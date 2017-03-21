"""
  % VOLUME Compute volumes of tets T defined over vertices V
  %
  % v = volume(V,T)
  % 
  % Inputs:
  %   V  #V by dim>=3 list of vertex positions
  %   T  #T by 4 list of tetrahedra indices
  % Ouputs:
  %   v  #T list of tet volumes. Signed if dim = 3
  %
"""
function  volume(V,T)

  a = V[T[:,1],:];
  b = V[T[:,2],:];
  c = V[T[:,3],:];
  d = V[T[:,4],:];
  # http://en.wikipedia.org/wiki/Tetrahedron#Volume
  # volume for each tetrahedron

  function cross2(a,b)
    # Optimizes r = cross(a,b,2), that is it computes cross products per row
    # Faster than cross if I know that I'm calling it correctly
    r =hcat(a[:,2].*b[:,3]-a[:,3].*b[:,2], 
        a[:,3].*b[:,1]-a[:,1].*b[:,3], 
        a[:,1].*b[:,2]-a[:,2].*b[:,1]);
  end

  # Minus sign so that typical tetgen mesh has positive volume
  #v = -dot((a-d),cross2(b-d,c-d),2)./6./4;
  # Not sure where that ./4 came from...
  dot2(a,b) = [dot(a[i,:],b[i,:]) for i in 1:size(a,1)]
  v = -dot2((a-d),cross2(b-d,c-d))./6.;
end
