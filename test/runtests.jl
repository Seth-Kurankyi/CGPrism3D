using CGPrism3D
using Test
using LinearAlgebra

J,K,alpha = 2,2,.1
#function testA(J::Int64,K::Int64,alpha::Float64)
@time dataT = dataTet(J,K,alpha)
#dependence of tensorGlueTet3D on J,K,alpha is for dimension factors
@time tet2 = tensorGlueTet3D(dataT,dataT,[1,2,3],[1,2,3],dataT)
@time prA = tensorGlueTet3D(tet2,dataT,[2,7,9],[1,2,3],dataT)
#Using only Prism A to coarse grain 
@time pruvA = tensorGlueTet3D(prA,prA,[1,5,6,8,9],[2,4,6,12,11],dataT)
@time pruvAA = tensorSumTet3D(pruvA,[1],dataT)
#@time pruvAAn = tensorSumTet3D(pruvA,[1],dataT)
#choose Fsymb for higher level to do 2-2 move
@time fA0F = tensor22moveF(pruvAA,[18,15,6,2,7],K,K,alpha)
@time fA01F = tensor22moveF(fA0F,[17,18,10,9,8],K,K,alpha)
@time fA012F = tensor22moveF(fA01F,[14,13,3,2,4],K,K,alpha)
#dependence of fullSplitTet3D on J,K,alpha is for dimension factors
@time f133A = fullSplitTet3D(fA012F,[9,6,18],[1,2,3,4,5,10,12,13,14,15,16,17],[7,8,11],dataT,J,K,alpha)
@time CGPrsmA1 = fullSplitTet3D(f133A[2],[8,10,2],[1,3,5,6,7,9,12,14,15],[13,4,11],dataT,J,K,alpha)
    #length(CGPrsmA1[1][2]),length(CGPrsmA1[2][2])
#end
tp = tensorPermute(CGPrsmA1[2],[5,1,11,2,6,3,10,12,7,8,4,9])
@show length(prA[2]),length(CGPrsmA1[2][2]), sort(tp[1]) == sort(prA[1])
#function testB(J::Int64,K::Int64,alpha::Float64)
#    dataT = dataTet(J,K,alpha)
    #dependence of tensorGlueTet3D on J,K,alpha is for dimension factors
#    tet2 = tensorGlueTet3D(dataT,dataT,[1,2,3],[1,2,3],J,K,alpha)
#    prA = tensorGlueTet3D(tet2,dataT,[2,7,9],[1,2,3],J,K,alpha)
    #Using only Prism A to coarse grain 
#    pruvA = tensorGlueTet3D(prA,prA,[1,5,6,8,9],[2,4,6,12,11],J,K,alpha)
#    pruvAA = tensorSumTet3D(pruvA,[1],J,K,alpha)
    #choose Fsymb for higher level to do 2-2 move
#    fA0F = tensor22move(pruvAA,[18,15,6,2,7],K,K,alpha)
#    fA01F = tensor22move(fA0F,[17,18,10,9,8],K,K,alpha)
#    fA012F = tensor22move(fA01F,[14,13,3,2,4],K,K,alpha)
    #dependence of fullSplitTet3D on J,K,alpha is for dimension factors
#    f133A = fullSplitTet3D(fA012F,[9,6,18],[1,2,3,4,5,10,12,13,14,15,16,17],[7,8,11],J,K,alpha)
#    CGPrsmA1 = fullSplitTet3D(f133A[2],[8,10,2],[1,3,5,6,7,9,12,14,15],[13,4,11],J,K,alpha)
    #length(CGPrsmA1[1][2]),length(CGPrsmA1[2][2])
#end
#for i in 1:5
#    J,K,alpha = i,i,.1
#    @time dataT = dataTet(J,K,alpha)
#    @time tet2 = tensorGlueTet3D(dataT,dataT,[1,2,3],[1,2,3],dataT)
#    @time prA = tensorGlueTet3D(tet2,dataT,[2,7,9],[1,2,3],dataT)
#    @time cgprismA = testA(J,K,alpha)[2]
#    @show length(prA[2]),length(cgprismA[2])
#end
#@time cgprismB = testB(2,2,1.)[2]
#length(cgprismB[2])