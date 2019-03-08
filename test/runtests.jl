using CGPrism3D
using Test
using LinearAlgebra


function testA(J::Int64,K::Int64,alpha::Float64)
    dataT = dataTet(J,K,alpha)
    #dependence of tensorGlueTet3D on J,K,alpha is for dimension factors
    tet2 = tensorGlueTet3D(dataT,dataT,[1,2,3],[1,2,3],J,K,alpha)
    prA = tensorGlueTet3D(tet2,dataT,[2,7,9],[1,2,3],J,K,alpha)
    #Using only Prism A to coarse grain 
    pruvA = tensorGlueTet3D(prA,prA,[1,5,6,8,9],[2,4,6,12,11],J,K,alpha)
    pruvAA = tensorSumTet3D(pruvA,[1],J,K,alpha)
    #choose Fsymb for higher level to do 2-2 move
    fA0F = tensor22moveA(pruvAA,[18,15,6,2,7],K,K,alpha)
    fA01F = tensor22moveA(fA0F,[17,18,10,9,8],K,K,alpha)
    fA012F = tensor22moveA(fA01F,[14,13,3,2,4],K,K,alpha)
    #dependence of fullSplitTet3D on J,K,alpha is for dimension factors
    f133A = fullSplitTet3D(fA012F,[9,6,18],[1,2,3,4,5,10,12,13,14,15,16,17],[7,8,11],J,K,alpha)
    CGPrsmA1 = fullSplitTet3D(f133A[2],[8,10,2],[1,3,5,6,7,9,12,14,15],[13,4,11],J,K,alpha)
    #length(CGPrsmA1[1][2]),length(CGPrsmA1[2][2])
end

function testB(J::Int64,K::Int64,alpha::Float64)
    dataT = dataTet(J,K,alpha)
    #dependence of tensorGlueTet3D on J,K,alpha is for dimension factors
    tet2 = tensorGlueTet3D(dataT,dataT,[1,2,3],[1,2,3],J,K,alpha)
    prA = tensorGlueTet3D(tet2,dataT,[2,7,9],[1,2,3],J,K,alpha)
    #Using only Prism A to coarse grain 
    pruvA = tensorGlueTet3D(prA,prA,[1,5,6,8,9],[2,4,6,12,11],J,K,alpha)
    pruvAA = tensorSumTet3D(pruvA,[1],J,K,alpha)
    #choose Fsymb for higher level to do 2-2 move
    fA0F = tensor22move(pruvAA,[18,15,6,2,7],K,K,alpha)
    fA01F = tensor22move(fA0F,[17,18,10,9,8],K,K,alpha)
    fA012F = tensor22move(fA01F,[14,13,3,2,4],K,K,alpha)
    #dependence of fullSplitTet3D on J,K,alpha is for dimension factors
    f133A = fullSplitTet3D(fA012F,[9,6,18],[1,2,3,4,5,10,12,13,14,15,16,17],[7,8,11],J,K,alpha)
    CGPrsmA1 = fullSplitTet3D(f133A[2],[8,10,2],[1,3,5,6,7,9,12,14,15],[13,4,11],J,K,alpha)
    #length(CGPrsmA1[1][2]),length(CGPrsmA1[2][2])
end

@time cgprismA = testA(2,2,1.)[2]
#@time cgprismB = testB(2,2,1.)[2]
@show length(cgprismA[2])#length(cgprismB[2])