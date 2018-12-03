using CGPrism3D
using Test
using LinearAlgebra

const K = 3
const x = K+1
const y = K/2


@time dataA = dataPrsmA();
@time dataB = dataPrsmB();
#@time dataB2 = dataPrsmB2();
# Split the prism into tetrahedron and pyramid along the triangular face C = spins(jd1,jd2,je1) = (j4,j5,j6)
#@time for ja in 0.:0.5:y, jb in 0.:0.5:y, jc in 0.:0.5:y
#    sing = tensorSplitSVD(dataB,[1,2,3],[7,8,9,10,11,12],[4,5,6],[ja,jb,jc])[2]
#    if sing[1] != 0
#        println(sing[1:3])
#    end
#end

t6j =[]
#for ja in 0.:0.5:y, jb in 0.:0.5:y, jc in 0.:0.5:y, jd in 0.:0.5:y, je in 0.:0.5:y, jf in 0.:0.5:y
#    if delta(ja,jb,jc) != 0 
#        #dims = prod(visqrt.([ja,jb,jc,jd,je,jf]))
#        sol = numchop(Tetra6j(ja,jb,jc,jd,je,jf))
#        if sol != 0
#            push!(t6j,sol)
#        end
#    end
#end
#@show t6j
#@show sortperm(dataB[1]) == 1:488

@time tta = fullSplitPrism3D(dataA,[1,2,3],[7,8,9,10,11,12],[4,5,6]) # split into U, V
@time ttb = fullSplitPrism3D(dataB,[1,2,3],[7,8,9,10,11,12],[4,5,6]) # split into U, V
#@time ttb2 = fullSplitPrism3D(dataB2,[1,2,3],[7,8,9,10,11,12],[4,5,6]) # split into U, V

#@time tensorGlue(tta[2],ttb[1],[1,7,8],[5,3,2]) # glue along different faces to get prsmUAVB
#@time tensorGlue(ttb[2],tta[1],[2,5,6],[4,2,1]) # glue along different faces to get prsmUBVA

@time uavb = tensorGluePrism3D(ttb[2],tta[1],[1,7,8],[5,3,2])
@time ubva = tensorGluePrism3D(tta[2],ttb[1],[2,5,6],[4,2,1])



#@show sortperm(ttb[1][1])
#@show sortperm(tta[1][1])
#@show   tta[1][1] ==  ttb[1][1]
#@show sort(tta[1][1])[15:20]  #sort(ttb[2][1])
#@show (tta[1][2][sortperm(tta[1][1])][15:20], ttb[1][2][sortperm(ttb[1][1])][15:20] )
#@show  numchop.( tta[1][2][sortperm(tta[1][1])] - ttb[1][2][sortperm(ttb[1][1])] )

#@time uv = tensorGluePrism3D(ttb[2],ttb[1],[1,2,3],[4,5,6])
#@show length(uv[2])


#@show numchop.(sort(uv[2]) - sort(dataA[2]))


#@time uv = tensorGlue(ttb[2],tta[1],[1,2,3],[4,5,6])
@show (length(uavb[2]),length(ubva[2]) )

#@time @show length(tensorGlue(uavb,ubva,[5,7,9,10,12],[7,5,9,10,11])[1]) # full effective prism including bulk variable
#@time @show length(tensorGlue(dataA,dataA,[2,3,5,8,9],[2,1,4,10,11])[1])

#@time @show length(tensorGluePrism3D(uavb,ubva,[1,7,8,5,6],[2,6,5,7,8])[1])
@time preff = tensorGlue(uavb,ubva,[5,7,9,10,12],[7,5,9,10,11])
@time preff = tensorGluePrism3D(uavb,ubva,[5,7,9,10,12],[7,5,9,10,11])
#@time preffA = tensorGluePrism3D(dataA,dataA,[2,3,5,8,9],[2,1,4,10,11])
#@show (length(preff[1]),length(preffA[1]) )
@time @show  length(preff[1])
@time @show  length(tensorSum(preff,[7])[1])
#@time @show length(tensorSum(preffA,[5])[1]) 

#@time prUAVB = generalGlue(dataA,dataA,[2,3,5,8,9],[1,2,4,10,11])
#@time @show length(fullSplitPrism3D(dataA,[1,2,3],[7,8,9,10,11,12],[4,5,6])[1][2])
