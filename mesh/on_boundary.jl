function on_boundary{IntT<:Integer}(F::Array{IntT,2})
  # ON_BOUNDARY  determines if each face of a manifold mesh is on the boundary
  # (contains at least on boundary edge)
  #
  # [I,C] = on_boundary(F)
  #
  # Inputs:
  #   F  #F by simplex-size list of element indices
  # Outputs:
  #   I  #F list of bools, whether on boundary
  #   C  #F by simplex size matrix of bools, whether opposite edge on boundary
  simplex_size = size(F,2);
  if simplex_size == 3
    E = [F[:,2] F[:,3]; F[:,3] F[:,1]; F[:,1] F[:,2]]
    sortedE = mapslices(r->tuple(r...), sort(E,2), [2])

    counts = Dict{Tuple{IntT,IntT},IntT}()
    for k in sortedE 
      counts[k] = 0
    end
    for k in sortedE 
      counts[k] = counts[k] + 1 
    end
   
    C = reshape( map(k->counts[k]==1,sortedE), size(F) )
    I = mapslices( any, C, [2] )

    [I, C]
  elseif simplex_size == 4
    T = F
    allF = [ T[:,2] T[:,4] T[:,3]; T[:,1] T[:,3] T[:,4]; T[:,1] T[:,4] T[:,2]; T[:,1] T[:,2] T[:,3]; ]
    # sort rows so that faces are reorder in ascending order of indices
    sortedE = mapslices(r->tuple(r...), sort(allF,2), [2])

    counts = Dict{Tuple{IntT,IntT,IntT},IntT}()
    for k in sortedE 
      counts[k] = 0
    end
    for k in sortedE 
      counts[k] = counts[k] + 1 
    end
   
    C = reshape( map(k->counts[k]==1,sortedE), size(F) )
    I = mapslices( any, C, [2] )

    [I,C]
  else
    throw(str("Unsupported simplex size", simplex_size))
    [Array{Bool,2}(),Array{Bool,2}()]
  end
end

