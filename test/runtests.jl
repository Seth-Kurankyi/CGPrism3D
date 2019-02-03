using CGPrism3D
using Test
using LinearAlgebra
#using StatsBase

#const K = 3
#const x = K+1
#const y = K/2



@time dataT = dataTet();
@time dataF = dataFsymb();
@time dataA = dataPrsmA();

# Testing Penstagon Identity -- related to 3-2 Pachner move

@time f2 = tensorGlue(dataF,dataF,[2,4,6],[2,3,1])
@time f3 = tensorGlue(f2,dataF,[5,8,7,3,4],[1,2,4,5,6])
@time fs1 = tensorSum(f3,[4])
@show sum(f2[2]-fs1[2])

@time t2 = tensorGlueTet3D(dataT,dataT,[1,2,3],[1,2,3])
@time t3 = tensorGlueTet3D(t2,dataT,[1,5,6,9,8],[1,2,3,5,6])
@time ts1 = tensorSumTet3D(t3,[1])
@show sum(t2[2]-ts1[2])

# Tests on Prisms -- Full Algorithm 

# Create a prism from three tetrahedra -- equivalent to dataA
@time tet2 = tensorGlueTet3D(dataT,dataT,[1,2,3],[1,2,3])
@time prA = tensorGlueTet3D(tet2,dataT,[2,7,9],[1,2,3])
@show sum(prA[2]-dataA[2])

## 2-2 move on faces [5,6,9,8,1] and [4,6,11,12,2] to switch diagonals -- gets prism B

@time fa1 = tensor22moveA(prA,[5,6,8,9,1]) # face [i,l,j,k,n]
@time prB = tensor22moveA(fa1,[4,6,12,11,2])
@show sum(prB[2]-prA[2])


# Split prisms into tetrahedra and pyramid 

@time spltA = fullSplitTet3D(prA,[4,5,6],[7,8,9,10,11,12],[1,2,3]) # split into U, V
ua, va = spltA[1], spltA[2]

@time spltB = fullSplitTet3D(prB,[6,9,11],[3,4,5,7,8,12],[10,2,1]) # split into U, V
ub, vb = spltB[1], spltB[2]

@show (sum(ub[2]-dataT[2]),sum(vb[2]-tet2[2]))
@show (sum(ua[2]-dataT[2]),sum(va[2]-tet2[2]))

# Glue tetrahedra and pyramids back into effective prism

@time uavb = tensorGlueTet3D(ua,vb,[1,3,5],[2,6,8])
@time ubva = tensorGlueTet3D(ub,va,[1,3,5],[6,5,8])
@show (sum(uavb[2]-prA[2]),sum(ubva[2]-prA[2]))


#Using only Prism A to coarse grain 

@time pruvA = tensorGlueTet3D(prA,prA,[1,5,6,8,9],[2,4,6,12,11])
@time pruvAA = tensorSumTet3D(pruvA,[1])

# 2-2 move /F-move on faces 0,1,2,(3,4)- not necessary ?
#@time fA0FN = tensor22moveA(pruvAA,[18,15,6,2,7])
@time fA0F = tensor22moveA(pruvAA,[18,15,6,2,7])
#@time fA1F = tensor22move(pruvAA,[17,18,10,9,8])
#@time fA2F = tensor22move(pruvAA,[14,13,3,2,4])
#@time fA3F = tensor22move(pruvAA,[17,16,14,5,12])
#@time fA4F = tensor22move(pruvAA,[10,11,5,3,1])
#@show (length(fA0F[2]),length(fA1F[2]),length(fA2F[2]),length(fA3F[2]),length(fA4F[2]))


@time fA01F = tensor22moveA(fA0F,[17,18,10,9,8])
@time fA012F = tensor22moveA(fA01F,[14,13,3,2,4])
#@time fA0123F = tensor22move(fA012F,[17,16,14,5,12])
#@time fA01234F = tensor22move(fA0123F,[10,11,5,3,1])
#@show (length(fA012F[2]),length(fA0123F[2]),length(fA01234F[2]))

@time f133A = fullSplitTet3D(fA012F,[9,6,18],[1,2,3,4,5,10,12,13,14,15,16,17],[7,8,11])
#@time f233A = fullSplitTet3D(fA012F,[13,15,2],[1,3,5,6,8,9,10,11,12,14,17,18],[7,4,16])
#@show (length(f133A[1][2]),length(f133A[2][2]),length(f233A[1][2]),length(f233A[2][2]))

@time CGPrsmA1 = fullSplitTet3D(f133A[2],[8,10,2],[1,3,5,6,7,9,12,14,15],[13,4,11])
#@time CGPrsmA2 = fullSplitTet3D(f233A[2],[6,4,12],[1,2,3,7,9,10,11,14,15],[13,5,8])
#@show (length(CGPrsmA1[1][2]),length(CGPrsmA1[2][2]),length(CGPrsmA2[1][2]),length(CGPrsmA2[2][2]))
@show ( sum( CGPrsmA1[2][2] - prA[2]) ) #, sum( CGPrsmA2[2][2] - prA[2]) )

# Now Using prsm A and B to glue and coarsegrain

#@time pruv = tensorGlueTet3D(uavb,ubva,[1,2,6,7,8],[3,2,4,10,9]) # get effective prism
#@time pruvAB = tensorSumTet3D(pruv,[1])# sum over bulk index 

# Coarse graining procedures 
# 2-2 move /F-move on faces 0,1,2

#@time Face0 = tensor22move(pruvAB,[12,15,8,9,6])
#@time Face1 = tensor22move(pruvAB,[14,12,3,2,1])
#@time Face2 = tensor22move(pruvAB,[17,16,11,9,7])
#@show (length(Face0[2]),length(Face1[2]),length(Face2[2]))

# It doesn't matter the order of moves 
#@time Face01 = tensor22move(Face0,[14,12,3,2,1])
#@time Face10 = tensor22move(Face1,[12,15,8,9,6])
#@time Face02 = tensor22move(Face0,[17,16,11,9,7])
#@time Face20 = tensor22move(Face2,[12,15,8,9,6])
#@time Face12 = tensor22move(Face1,[17,16,11,9,7])
#@time Face21 = tensor22move(Face2,[14,12,3,2,1])
#@show (length(Face01[2]),length(Face10[2]),length(Face02[2]),length(Face20[2]),length(Face12[2]),length(Face21[2]))

#@time Face012 = tensor22move(Face01,[17,16,11,9,7])
#@time Face102 = tensor22move(Face10,[17,16,11,9,7])
#@show (length(Face012[2]),length(Face102[2]))

# Final Pachner moves-- 3-1 pachner moves using SVD 

#@time CGtet1 = fullSplitTet3D(Face012,[8,12,2],[3,4,5,7,9,11,13,14,15,16,17,18],[1,6,10])
#@time CGtet2 = fullSplitTet3D(Face012,[15,16,9],[1,2,3,4,5,8,10,11,12,13,14,17],[7,6,18])
#@show (length(CGtet1[1][1]),length(CGtet1[2][1]),length(CGtet2[1][1]),length(CGtet2[2][1]))


#@time CGtetfull1 = fullSplitTet3D(CGtet2[2],[6,2,9],[3,4,5,8,10,11,12,13,15],[1,14,7])
#@time CGtetfull2 = fullSplitTet3D(CGtet1[2],[9,10,5],[1,2,3,6,7,8,11,13,15],[4,14,12])
#@show (length(CGtetfull1[2][1]),length(CGtetfull1[2][2]),length(CGtetfull2[2][1]),length(CGtetfull2[2][2]) )
#@show fans = sum( CGtetfull2[2][2] - dataA[2])
#if numchop(fans) == 0
#    println("SUCCESS -- :)")
#end