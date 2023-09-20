using CGPrism3D
using Test
using LinearAlgebra

J,K,alpha = 2,3,0.5
dataT = dataTet(J,K,alpha)
@show unique(dataT[2])
#dependence of tensorGlueTet3D on J,K,alpha is for dimension factors
#tet2 = tensorGlueTet3DN(dataT,dataT,[1,2,3],[1,2,3],dataT)
tb2 = tensorBlock(dataT,[],[2,3,4,5,6],[1],[0.,0.5,0.5])[1]
#U,s,V =svd(tb2)
#@show s
#@time prA = tensorGlueTet3DN(tet2,dataT,[2,7,9],[1,2,3],dataT)
#@show length(prA[1])

#@time prA2 = tensorGlueTet3DN2(tet2,dataT,[2,7,9],[1,2,3],dataT)
#@show length(prA[1]), length(prA2[1])
#@show sum(prA[2]-prA2[2])
function testA(J::Int64,K::Int64,alpha::Float64,loop)
    dataT = dataTet(J,K,alpha)
    #dependence of tensorGlueTet3D on J,K,alpha is for dimension factors
    tet2 = tensorGlueTet3DN(dataT,dataT,[1,2,3],[1,2,3],dataT)
    prA = tensorGlueTet3DN(tet2,dataT,[2,7,9],[1,2,3],dataT)

    #Using only Prism A to coarse grain
    for l in 1:loop
        println("loop=$l")
        #dataT,dataV = fullSplitTet3DN(prA,[4,5,6],[7,8,9,10,11,12],[1,2,3],dataT,J,K,alpha)
        #dataT = tensorPermuteN(dataT,[4,5,6,1,2,3])
        #prA =
        pruvA = tensorGlueTet3DN(prA,prA,[1,5,6,8,9],[2,4,6,12,11],dataT)
        pruvAA = tensorSumTetO(pruvA,[1],dataT)
        #@time pruvAAn = tensorSumTet3D(pruvA,[1],dataT)
        #choose Fsymb for higher level to do 2-2 move
        #println("Starting 2-2 move on face 1...")
        #fA0F = tensor22moveT(pruvAA,[18,15,6,2,7],dataT)
        fA0F = tensor22moveF(pruvAA,[18,15,6,2,7],K,K,alpha)
        #fA0F = tensor22moveT2(pruvAA,[18,15,6,2,7],dataT,[7,1,17],K)
        #println("2-2 move on face 1 done")
        #println("Starting 2-2 move on face 2...")
        #fA01F = tensor22moveT(fA0F,[17,18,10,9,8],dataT)
        fA01F = tensor22moveF(fA0F,[17,18,10,9,8],K,K,alpha)
        #fA01F = tensor22moveT2(fA0F,[17,18,10,9,8],dataT,[8,7,11],K)
        #println("2-2 move on face 2 done")
        #println("Starting 2-2 move on face 3...")
        #fA012F = tensor22moveT(fA01F,[14,13,3,2,4],dataT)
        fA012F = tensor22moveF(fA01F,[14,13,3,2,4],K,K,alpha)
        #fA012F = tensor22moveT2(fA01F,[14,13,3,2,4],dataT,[4,7,16],K)
        #println("2-2 move on face 3 done")
        #println("Starting 3-1 move on tet1...")
        f133A = fullSplitTet3DN(fA012F,[9,6,18],[1,2,3,4,5,10,12,13,14,15,16,17],[7,8,11],dataT,J,K,alpha)
        #println("3-1 move on tet1 done")
        #println("Starting 3-1 move on tet2...")
        CGPrsmA1 = fullSplitTet3DN(f133A[2],[8,10,2],[1,3,5,6,7,9,12,14,15],[13,4,11],dataT,J,K,alpha)[2]
        tperm = tensorPermuteN(CGPrsmA1,[5,1,11,2,6,3,10,12,7,8,4,9])
        prA = tensorProjectPrism(tperm,K)
        dataT = tensorPermuteN(fullSplitTet3DN(prA,[7,8,9,10,11,12],[4,5,6],[1,2,3],dataT,J,K,alpha)[2],[4,5,6,1,2,3])
        @show length(prA[1])
        #@show println(dataT1[1])
        #@show println(dataT[1])
        #@show sort(dataT1[1]) == sort(dataT[1])
        tb = tensorBlock(prA,[4,5,6],[7,8,9,10,11,12],[1,2,3],[0.,0.5,0.5])[1]
        @show sum(prA[2])
        U,s,V = svd(tb)
        @show s
    end
    return prA,dataT
end



#for i in 1:3
@time ff = testA(J,K,alpha,30)
#fullSplitTet3DN(ff,)
#@time ff1 = testA(J,K,alpha)
#tb = tensorBlock(ff[1],[4,5,6],[7,8,9,10,11,12],[1,2,3],[0.,0.5,0.5])[1]
#@show sum(ff[2][2])
#U,s,V = svd(tb)
#@show s
#end

#@show length(ff[1])
