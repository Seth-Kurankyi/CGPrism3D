module CGPrism3D
using LinearAlgebra
#using StatsBase

# Define chop function 
export numchop
numchop(num::Real)=(abs(num) >= 1000*eps() ? num : zero(Real))

function numchop(num::Complex)
    numchop(imag(num)) == 0. ? numchop(real(num)) : Complex(numchop(real(num)), numchop(imag(num)))
end

# Define k-level of SU(2)_k
#Define all functions needed for G-symbol
# Define delta{ijk} -> coupling rules
function delta(i::Float64,j::Float64,k::Float64,K::Int64)::Int64
    sol = 0
    if i <=(j+k) && j<=(i+k) && k<=(i+j) && i+j+k <= K && 2*(i+j+k)%2 ==0
        sol = 1
    end
    return sol
end

#Define quantum numbers qn (this is real)
function qn(n::Float64,K::Int64)
    sol = (exp(pi*n*im/(K+2)) - exp(-pi*n*im/(K+2))) / (exp(pi*im/(K+2)) - exp(-pi*im/(K+2))) 
    return real(sol)
end

#Define qn factorial
function qnfact(n::Float64,K::Int64)
    sol = 1.
    for i in 1:n
        sol *= qn(i,K)
    end
    return sol
end

#Define square root of quantum dimension
function visqrt(i::Float64,K::Int64)
    sol = ((-1+0im)^i )*sqrt(qn(2*i+1,K))
    return sol
end

#Define triangle equality 
function trian(i::Float64,j::Float64,k::Float64,K::Int64)
    sol = 0.
    if delta(i,j,k,K) == 1
        sol = delta(i,j,k,K)*sqrt(qnfact(i+j-k,K)*qnfact(i-j+k,K)*qnfact(-i+j+k,K)/qnfact(i+j+k+1,K))
    end
    return sol
end

# Define Racah-Wigner six-j symbol
function RacWig6j(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64,K::Int64)
    a = i+j+m; b = i+l+n; c = k+j+n; d = k+l+m;  e = i+j+k+l; f = i+k+m+n; g = j+l+m+n
    sol = 0.
    if delta(i,j,m,K) != 0 && delta(i,l,n,K) != 0 && delta(k,j,n,K) != 0 && delta(k,l,m,K) != 0  
        sumz = 0.
        for z in max(a,b,c,d):min(e,f,g)
            sumz += (-1)^z *qnfact(z+1,K)/
                ((qnfact(e-z,K)*qnfact(f-z,K)*qnfact(g-z,K))* (qnfact(z-a,K)*qnfact(z-b,K)*qnfact(z-c,K)*qnfact(z-d,K)))
        end
        sol = trian(i,j,m,K)*trian(i,l,n,K)*trian(k,j,n,K)*trian(k,l,m,K)*sumz
    end
    return sol
end

#Define F-symbol
function Fsymb(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64,K::Int64)
    sol = 0.
    if delta(i,j,m,K) != 0 && delta(i,l,n,K) != 0 && delta(k,j,n,K) != 0 && delta(k,l,m,K) != 0
        sol = (-1+0im)^(i+j+k+l)*sqrt(qn(2*m+1,K)*qn(2*n+1,K)) * RacWig6j(i,j,m,k,l,n,K)
    end
    return sol
end

#Define G-symbol
function Gsymb(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64,K::Int64)
    sol = 0.
    if delta(i,j,m,K) != 0 && delta(i,l,n,K) != 0 && delta(k,j,n,K) != 0 && delta(k,l,m,K) != 0
        sol = Fsymb(i,j,m,k,l,n,K) /(visqrt(m,K)*visqrt(n,K))
    end
    return sol
end

export Tetra6j, delta, visqrt
#Define G-symbol times dimensions
function Tetra6j(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64,K::Int64)
    sol = 0.
    if delta(i,j,m,K) != 0 && delta(i,l,n,K) != 0 && delta(k,j,n,K) != 0 && delta(k,l,m,K) != 0
        dims = prod(visqrt.([i,j,m,k,l,n],K))
        sol =  numchop(dims*Gsymb(i,j,m,k,l,n,K))
    end
    return sol
end

#Define G-symbol 2
function Gsymb2(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64,K::Int64)
    sol = 0.
    if delta(i,j,m,K) != 0 && delta(i,l,n,K) != 0 && delta(k,j,n,K) != 0 && delta(k,l,m,K) != 0
        sol = (-1+0im)^(i+j+k+l+m+n)*RacWig6j(i,j,m,k,l,n,K)
    end
    return sol
end

#Define square root of quantum dimension
function visqrtT(a::Float64,J::Int64,K::Int64,alpha::Float64)
    ((-1+0im)^a )*sqrt(abs(TetraJK(0.,0.,0.,a,a,a,J,K,alpha)))
end

function TetraJK(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64,J::Int64,K::Int64,alpha::Float64)
    #@assert J <= K
    alpha*Tetra6j(i,j,m,k,l,n,J) +(1-alpha)*Tetra6j(i,j,m,k,l,n,K)
end

function visqrtU(Utensor,a::Float64)
    indx = Utensor[1]
    amps = Utensor[2]
    fc = findfirst(x-> x == [0.,0.,0.,a,a,a], indx)
    ((-1+0im)^a )*sqrt(abs(amps[fc]))
end

export dataTet, dataFsymb

# get all non-zero amplitudes and spins 

# get all non-zero tetrahedron amplitude and spins (also for F-symbol) 

function dataTet(J::Int64,K::Int64,alpha::Float64)
    @assert J <= K
    @assert 0. <= alpha <= 1.
    x = K+1
    y = K/2
    t6j = Array{Float64,1}[]
    t6s = Float64[]
    for ja in 0.:0.5:y, jb in 0.:0.5:y, jc in 0.:0.5:y
        if delta(ja,jb,jc,K) != 0
            for jd in 0.:0.5:y, je in 0.:0.5:y
                if delta(jc,jd,je,K) != 0
                    for jf in 0.:0.5:y
                        if delta(ja,je,jf,K) != 0 && delta(jb,jd,jf,K) != 0
            #dims = prod(visqrt.([ja,jb,jc,jd,je,jf]))
                            sol = numchop(TetraJK(ja,jb,jc,jd,je,jf,J,K,alpha))
                            if sol != 0
                                push!(t6s,sol)
                                push!(t6j,[ja,jb,jc,jd,je,jf])
                            end
                        end
                    end
                end
            end
        end
    end
    return t6j,t6s
end

function dataFsymb(J::Int64,K::Int64,alpha::Float64)
    @assert J <= K
    @assert 0. <= alpha <= 1.
    x = K+1
    y = K/2
    t6j = Array{Float64,1}[]
    t6s = Float64[]
    for ja in 0.:0.5:y, jb in 0.:0.5:y, jc in 0.:0.5:y
        if delta(ja,jb,jc,K) != 0
            for jd in 0.:0.5:y, je in 0.:0.5:y
                if delta(jc,jd,je,K) != 0
                    for jf in 0.:0.5:y
                        if delta(ja,je,jf,K) != 0 && delta(jb,jd,jf,K) != 0
            #dims = prod(visqrt.([ja,jb,jc,jd,je,jf]))
                            sol = numchop(alpha*Fsymb(ja,jb,jc,jd,je,jf,J) +(1-alpha)*Fsymb(ja,jb,jc,jd,je,jf,K))
                            if sol != 0
                                push!(t6s,sol)
                                push!(t6j,[ja,jb,jc,jd,je,jf])
                            end
                        end
                    end
                end
            end
        end
    end
    return t6j,t6s
end

# Operation on tensors
# 1. Gluing tensors, 2. Splitting tensors, 3. Summing over tensor indices 

export tensorBlock, fullSplitTet3D, fullSplit3D, tensorGlue, tensorSum, tensorGlueTet3D, tensor22move, tensorSumTet3D

# a general function to split a tensor M^C_{AB} to U^C_{Ai} and  V^C_{iB} with C the shared face.. for particular values of C
# first reduce to matrix form and then perform svd to get U, s, and V

#This reduces the tensor M into matrix form M^C_{AB} with AB labelling matric indices. 
function tensorBlock(tensorM,posnA,posnB,posnC,spinC)
    indx = tensorM[1]
    ampvals = tensorM[2]
    fc = findall(x-> x[posnC] == spinC, indx) #find spinsC in all spins length(posnC) = length(spinC)
    ampInfo = indx[fc] # this contains all spin configuration containing ampsC
    amps = ampvals[fc] # this contains the corresponding spin amplitude values for ampsC
    Acol = sort(unique(getindex.(ampInfo, [posnA])))
    Brow = sort(unique(getindex.(ampInfo, [posnB])))
    indxU = vcat.(Acol,repeat([spinC],length(Acol)))
    indxV = vcat.(Brow,repeat([spinC],length(Brow)))
    #indxU = []
    #indxV = []
    matM = zeros(length(Acol),length(Brow))
    if length(matM) != 0
        for i in 1:length(ampInfo)
            qu = findfirst(x-> x == ampInfo[i][posnA], Acol)[1]# [1] since it returns an array [1]
            qv = findfirst(x-> x == ampInfo[i][posnB], Brow)[1]
            matM[qu,qv] = amps[i]
            #push!(indxU,vcat(Acol[qu],spinC))
            #push!(indxV,vcat(Brow[qv],spinC))
        end
        #U, s, V = svd(matM)
        return matM,indxU,indxV
    else
        return 0,0,0
    end
end

## Glue any two tensors A, B along shared faces .. 
function tensorGlue(tensorB,tensorA,posnB,posnA)# glue two tensors TA,TB along posnA from TA and posnB from TB 
    indxA = tensorA[1] # spins of TA
    indxB = tensorB[1] # spins of TB
    ampsA = tensorA[2] # amplitudes of TA
    ampsB = tensorB[2] # amplitudes of TB
    shrdA = getindex.(indxA,[posnA])
    shrdB = getindex.(indxB,[posnB])
    lena = collect(1:length(indxA[1]))
    deleteat!(lena,sort(posnA)) # assert needs posnA to be a vector
    #if length(posnA) > 1 deleteat!(lena,sort(posnA,dims=2) else  deleteat!(lena,sort(posnA)) end
    shrdAB = []#Tuple{Array{Float64,1},Int64,Int64}[]
    shrddB = []
    shrdABn = []
    qq = Array{Array{Float64,1},1}[]
    amps = Array{Float64,1}[]
    sa = sortperm(shrdA)
    sb = sortperm(shrdB)
    shrdA = shrdA[sa]
    shrdB = shrdB[sb]
    sindxA = indxA[sa]
    sindxB = indxB[sb]
    sampsA = ampsA[sa]
    sampsB = ampsB[sb]
    shrdAu = unique(shrdA)
    shrdBu = unique(shrdB)
    shrd = shrdAu
    #println(shrdAu == shrdBu)
    if shrdAu == shrdBu
        nothing
    else
        for i in shrdBu
            #shrd = shrdAu
            if !(i in shrdAu)
                push!(shrd,i)
            end
        end
    end
    for i in shrd
        cna = count(x->x==i,shrdA)
        cnb = count(x->x==i,shrdB)
        push!(shrdAB, [i, cna,cnb ])
        #println("Ca,Cna = ",ca,"\t",cna)
        #println("Cb,Cnb = ",cb,"\t",cnb)
    end
    for i in shrdAB
        if i[2] == 0 || i[3] == 0 # Put condition for when one of a or b is zero
            nothing
        else
            for j in 1:i[2]
                ans = @. sampsA[j]*sampsB[1:i[3]]
                push!(amps,ans)
                indxa = repeat([sindxA[j][lena]],i[3]) 
                #for j in indxA[i][lena]
                indxa1 = vcat.(sindxB[1:i[3]],indxa)
                push!(qq,indxa1)
            end
        end
        deleteat!(sindxA,1:i[2])
        deleteat!(sindxB,1:i[3])
        deleteat!(sampsA,1:i[2])
        deleteat!(sampsB,1:i[3])
        #cntb = cnt
    end
    indxsol = collect(Iterators.flatten(qq))
    sol = collect(Iterators.flatten(amps))
    return indxsol,sol
end

    for i in 1:length(face)
function tensorGlueTet3D(tensorB,tensorA,posnB,posnA,Utensor)# glue two tensors TA,TB along posnA from TA and posnB from TB 
    indxA = tensorA[1] # spins of TA
    indxB = tensorB[1] # spins of TB
    ampsA = tensorA[2] # amplitudes of TA
    ampsB = tensorB[2] # amplitudes of TB
    shrdA = getindex.(indxA,[posnA])
    shrdB = getindex.(indxB,[posnB])
    lena = collect(1:length(indxA[1]))
    deleteat!(lena,sort(posnA)) # assert needs posnA to be a vector
    #if length(posnA) > 1 deleteat!(lena,sort(posnA,dims=2) else  deleteat!(lena,sort(posnA)) end
    shrdAB = []#Tuple{Array{Float64,1},Int64,Int64}[]
    qq = Array{Array{Float64,1},1}[]
    amps = Array{Float64,1}[]
    sa = sortperm(shrdA)
    sb = sortperm(shrdB)
    shrdA = shrdA[sa]
    shrdB = shrdB[sb]
    sindxA = indxA[sa]
    sindxB = indxB[sb]
    sampsA = ampsA[sa]
    sampsB = ampsB[sb]
    shrdAu = unique(shrdA)
    shrdBu = unique(shrdB)
    shrd = shrdAu
    #println(shrdAu == shrdBu)
    if shrdAu == shrdBu
        nothing
    else
        for i in shrdBu
            #shrd = shrdAu
            if !(i in shrdAu)
                push!(shrd,i)
            end
        end
    end
    for i in shrd
        cna = count(x->x==i,shrdA)
        cnb = count(x->x==i,shrdB)
        push!(shrdAB, [i, cna,cnb ])
    end
    for i in shrdAB
        if i[2] == 0 || i[3] == 0 # Put condition for when one of a or b is zero
            nothing
        else
            for j in 1:i[2]
                ans = @. sampsA[j]*sampsB[1:i[3]]
                pp = 1
                for k in i[1]
                    pp *= visqrtU(Utensor,k)
                end
                ans /= pp
                ampsN = real(ans)+imag(ans)
                push!(amps,ampsN)
                indxa = repeat([sindxA[j][lena]],i[3]) 
                #for j in indxA[i][lena]
                indxa1 = vcat.(sindxB[1:i[3]],indxa)
                push!(qq,indxa1)
            end
        end
        deleteat!(sindxA,1:i[2])
        deleteat!(sindxB,1:i[3])
        deleteat!(sampsA,1:i[2])
        deleteat!(sampsB,1:i[3])
        #cntb = cnt
    end
    indxsol = collect(Iterators.flatten(qq))
    sol = collect(Iterators.flatten(amps))
    return indxsol,sol
end



#A general function to splits fully the tensor M^C_{AB} to tensors U^C_{Ai} and V^C_{iB}. 
# It uses tensorBlock, 'fixes' the sign problem from svd, and normalize 
# We keep only the first singular value (assuming a gemetric split)
# splits prism into tetrahedron and a  pyramid
# Multiply back by missing dimension factors 
function fullSplitTet3D(dataM,posnA,posnB,posnC,Utensor,J,K,alpha) # posnA,B,C must be vectors @assert length(dataM) = vcat(A,B,C)
    indx = dataM[1] 
    ampvals = dataM[2]
    #posAC = vcat(posnA,posnC) # assert posnA,posnC must be both vectors
    #posBC = vcat(posnB,posnC) # assert posnB,posnC must be both vectors
    #indxU = unique(getindex.(indx, [posAC]))
    #indxV = unique(getindex.(indx, [posBC]))
    ampsC = unique(getindex.(indx, [posnC])) # get all unique spins C (jd1,jd2,je1) 
    ampsU = [] # this will store all amplitudes for U^C_{A}
    ampsV = [] # this will store all amplitudes for V^C_{B}
    indxsU = []
    indxsV = []
    #indxUV = [] # want this to store all spins in the right order
    for i in 1:length(ampsC) # loop over the unique spins
        jd1 = ampsC[i][1]; jd2 = ampsC[i][2]; je1 = ampsC[i][3]  #jd1,jd2,je1 in spin C form a triangle
        if delta(jd1,jd2,je1,K) != 0
            spinC = [jd1,jd2,je1]
            block = tensorBlock(dataM,posnA,posnB,posnC,spinC)
            mat = block[1]
            blkU = block[2]
            blkV = block[3]
            U, s, V = svd(mat)
            #truncate and use only first singular value
            U1 = U[:,1] 
            V1 = V[:,1]
            s1 = s[1]
            pp = 1
            for j in ampsC[i]
                pp *= visqrtU(Utensor,j)
            end
            valU = U1*sqrt(s1)*sqrt(pp) # multiply by leftover dimension factors
            valV = V1*sqrt(s1)*sqrt(pp)
            # valU,valV are all either real or purely imaginary
            valU = real(valU) - imag(valU) # take -ve of imaginary part if it gives only imag
            valV = real(valV) + imag(valV)
            ja1 = blkU[1][1]; jb1 = blkU[1][2]; jg = blkU[1][3] # pick the first spins in U
            #Fixing sign problem -- keep the same sign for U and Tetra6j symbol 
            if sign(valU[1]) == sign(TetraJK(jd1,jd2,je1,ja1,jb1,jg,J,K,alpha) )
                push!(ampsU, valU)#/valU[1])
                push!(ampsV, valV)#/valV[1])
                push!(indxsU,blkU)
                push!(indxsV,blkV)
            else
                push!(ampsU, -valU)#/valU[1])
                push!(ampsV, -valV)#/valV[1])
                push!(indxsU,blkU)
                push!(indxsV,blkV)
            end    
        end
    end
    solU = collect(Iterators.flatten(ampsU))
    solV = collect(Iterators.flatten(ampsV))
    ansU = solU/solU[1] # normalization condition
    ansV = solV/solV[1]
    indxUs = collect(Iterators.flatten(indxsU))
    indxVs = collect(Iterators.flatten(indxsV))
    return (indxUs,ansU),(indxVs,ansV) #,indxUV
end



function fullSplit3D(dataM,posnA,posnB,posnC) # posnA,B,C must be vectors @assert length(dataM) = vcat(A,B,C)
    indx = dataM[1] 
    ampvals = dataM[2]
    #posAC = vcat(posnA,posnC) # assert posnA,posnC must be both vectors
    #posBC = vcat(posnB,posnC) # assert posnB,posnC must be both vectors
    #indxU = unique(getindex.(indx, [posAC]))
    #indxV = unique(getindex.(indx, [posBC]))
    ampsC = unique(getindex.(indx, [posnC])) # get all unique spins C (jd1,jd2,je1) 
    ampsU = [] # this will store all amplitudes for U^C_{A}
    ampsV = [] # this will store all amplitudes for V^C_{B}
    indxsU = []
    indxsV = []
    #indxUV = [] # if we want this will store all spins in the right order
    for i in 1:length(ampsC) # loop over the unique spins
        #jd1 = ampsC[i][1]; jd2 = ampsC[i][2]; je1 = ampsC[i][3]  #jd1,jd2,je1 in spin C form a triangle
        #if delta(jd1,jd2,je1) != 0
        #spinC = [jd1,jd2,je1]
        block = tensorBlock(dataM,posnA,posnB,posnC,ampsC[i])
        mat = block[1]
        blkU = block[2]
        blkV = block[3]
        if mat != 0
            U, s, V = svd(mat)
            #truncate and use only first singular value
            U1 = U[:,1] 
            V1 = V[:,1]
            s1 = s[1]
            valU = U1*sqrt(s1)#*sqrt(prod(visqrt.(ampsC[i]))) # multiply by leftover dimension factors
            valV = V1*sqrt(s1)#*sqrt(prod(visqrt.(ampsC[i])))
            # valU,valV are all either real or purely imaginary
            valU = real(valU) - imag(valU) # take -ve of imaginary part if it gives only imag
            valV = real(valV) + imag(valV)
            #ja1 = blkU[1][1]; jb1 = blkU[1][2]; jg = blkU[1][3] # pick the first spins in U
            #Fixing sign problem -- keep the same sign for U and Tetra6j symbol 
            #if sign(valU[1]) == sign(Tetra6j(jd1,jd2,je1,ja1,jb1,jg) )
            push!(ampsU, valU)#/valU[1])
            push!(ampsV, valV)#/valV[1])
            push!(indxsU,blkU)
            push!(indxsV,blkV)
            #else
                #push!(ampsU, -valU)#/valU[1])
                #push!(ampsV, -valV)#/valV[1])
                #push!(indxsU,blkU)
                #push!(indxsV,blkV)
            #end    
        end
    end
    solU = collect(Iterators.flatten(ampsU))
    solV = collect(Iterators.flatten(ampsV))
    ansU = solU/solU[1] # normalization condition
    ansV = solV/solV[1]
    indxUs = collect(Iterators.flatten(indxsU))
    indxVs = collect(Iterators.flatten(indxsV))
    return (indxUs,ansU),(indxVs,ansV) #,indxUV
end

# Summing a tensor  over indices labelled posnN
# Summing a tensor  over indices labelled posnN
function tensorSumN(ampJ,posnN,Utensor)
    indx = ampJ[1]
    amps = ampJ[2]
    lenS = collect(1:length(indx[1]))
    deleteat!(lenS,sort(posnN))
    #if length(posnN) > 1   deleteat!(lenS,sort(posnN,dims=2)) else  deleteat!(lenS,sort(posnN)) end
    indxef = getindex.(indx,[lenS]) # get all spins without N 
    #indxefN = getindex.(indx,[posnN])
    #indxefu = unique(indxef)
    tt = sortperm(indxef) # sort them spins 
    indxef = indxef[tt] # apply sort on spins 
    amps = amps[tt]  # apply sort on amps
    #fc = findall(x-> x[lenS]==ans ,indx)
    qq = []
    indq = []
    fc = sortperm(indxef)
    indxN = indxef[fc]
    amps = amps[fc]
    #indxN = indxefN[fc]
    if length(posnN) == 1
        blkdims = getindex.(indx,posnN)[fc]
        sqdm = []
        for j in blkdims
            push!(sqdm,visqrtU(Utensor,j))
        end
        ampsn = amps .*sqdm
        ampsN = real(ampsn) - imag(ampsn)
        #data = indx, ampsN2
        #ans = tensorSum(data,posnN)
    else
        blkdims = getindex.(indx,[posnN])[fc]
        sqdm = []
        for i in blkdims
            pp = 1
            for j in i
                pp *= numchop(visqrtU(Utensor,j))
            end
            push!(sqdm,pp)
        end
        #println(sqdm)
        ampsn = amps .*sqdm
        ampsN = real(ampsn) - imag(ampsn)
        #data = indx, ampsN2
        #ans = tensorSum(data,posnN)
    end
    
    for i in 2:length(indxN)
        if indxN[i] == indxN[i-1]
            ampsN[i] += ampsN[i-1]
        elseif numchop(ampsN[i-1]) != 0
            push!(indq,indxN[i-1])
            push!(qq,ampsN[i-1])
        end
    end
    push!(indq,indxN[end])
    push!(qq,ampsN[end])
    return indq, qq
end

export tensorSumN 

function tensorSum(ampJ,posnN)
    indx = ampJ[1]
    amps = ampJ[2]
    lenS = collect(1:length(indx[1]))
    deleteat!(lenS,sort(posnN))
    #if length(posnN) > 1   deleteat!(lenS,sort(posnN,dims=2)) else  deleteat!(lenS,sort(posnN)) end
    indxef = getindex.(indx,[lenS]) # get all spins without N 
    #indxefu = unique(indxef)
    tt = sortperm(indxef) # sort them spins 
    indxef = indxef[tt] # apply sort on spins 
    amps = amps[tt]  # apply sort on amps
    #fc = findall(x-> x[lenS]==ans ,indx)
    qq = []
    indq = []
    fc = sortperm(indxef)
    indxN = indxef[fc]
    ampsN = amps[fc]
    
    for i in 2:length(indxN)
        if indxN[i] == indxN[i-1]
            ampsN[i] += ampsN[i-1]
        elseif numchop(ampsN[i-1]) != 0
            push!(indq,indxN[i-1])
            push!(qq,ampsN[i-1])
        end
    end
    push!(indq,indxN[end])
    push!(qq,ampsN[end])
    return indq, qq
end

# Take care of dimension factors when summing. This makes sure Pachner 3-2 move works 
# For now this only works for summing over one index-- we can change this easily
function tensorSumTet3D(tensor,posnN,Utensor)
    indx = tensor[1]
    amps = tensor[2]
    if length(posnN) == 1
        blkdims = getindex.(indx,posnN)
        sqdm = []
        for j in blkdims
            push!(sqdm,visqrtU(Utensor,j))
        end
        ampsN = amps .*sqdm
        ampsN2 = real(ampsN) - imag(ampsN)
        data = indx, ampsN2
        ans = tensorSum(data,posnN)
    else
        blkdims = getindex.(indx,[posnN])
        sqdm = []
        for i in blkdims
            pp = 1
            for j in i
                pp *= numchop(visqrtU(Utensor,j))
            end
            push!(sqdm,pp)
        end
        #println(sqdm)
        ampsN = amps .*sqdm
        ampsN2 = real(ampsN) - imag(ampsN)
        data = indx, ampsN2
        ans = tensorSum(data,posnN)
    end
    return ans
end

export swap2, permuteInd, tensorPermute, tensor22move, tensor22moveF, tensor22moveT

# Performs an F-move on a face.. equivalent to 2-2 Pachner move -- this uses Fsymbol
function swap2(vec,a,b)
    vec[a],vec[b] = vec[b], vec[a]
    return vec
end

function tensor22move(tensor,face,J::Int64,K::Int64,alpha::Float64)# place middle edge at last position 
    # assert face contains 5 elements
    indx = tensor[1] # spins of TB
    amps = tensor[2]
    flast = face[end]
    y = K/2
    indxlast = getindex.(indx, [flast])
    famps = []
    findx = []
    for i in 1:length(indx)
        #iMod = copy(indx[i])
        j1 = indx[i][face[1]]; j2 = indx[i][face[2]]; j3 = indx[i][face[3]]; j4 = indx[i][face[4]]; j5 = indx[i][face[5]]
    #    #ans = 0#amps[i]
    #    ind = zeros(length(indx[1]))
        for j in 0.:0.5:y
            #iMod = deepcopy(indxf1u)
            if delta(j1,j3,j,K) != 0 && delta(j2,j4,j,K) != 0 && Fsymb(j1,j3,j,j4,j2,j5,K) != 0
                iMod = copy(indx[i])
                ampsA = numchop(amps[i] * Fsymb(j1,j3,j,j4,j2,j5,K))
                push!(famps,ampsA)
                iMod[flast] = j
                push!(findx,iMod)
                #push!(findx,insert!(deleteat!(iMod,flast),flast,j))
                #fc = findall(x-> x == iMod,indxf1 )
            end
        end
    end
    #indxf1 = getindex.(indx, [ff1]) # get all spins without spin N to be summed over 
    tt = sortperm(findx) # sort them spins 
    indxf1 = findx[tt] # apply sort on spins 
    famps = famps[tt]  # apply sort on amps
    #fc = findall(x-> x[lenS]==ans ,indx)
    qq = []
    indxeff = []
    for i in 2:length(indxf1)
        if indxf1[i] == indxf1[i-1]
            famps[i] += famps[i-1]
        elseif numchop(famps[i-1]) != 0
            push!(indxeff,indxf1[i-1])
            push!(qq,famps[i-1])
        end
    end
    push!(indxeff,indxf1[end])
    push!(qq,famps[end])
    #indxeff = unique(indxf1)
    return indxeff, qq
end

function tensor22moveT(tensor,posF,Utensor) # J,K,alpha determines what Fsymbol to use
    glu = tensorGlueTet3D(tensor,Utensor,posF,[1,5,2,4,6],Utensor)
    n = posF[end]
    m = length(glu[1][1])
    swp = @. swap2(glu[1],n,m)
    tensorN = swp, glu[2]
    ans = tensorSumTet(tensorN,[length(glu[1][1])],Utensor)
    return ans
end

function tensor22moveF(tensor,posF,J::Int64,K::Int64,alpha::Float64)
    glu = tensorGlue(tensor,dataFsymb(J,K,alpha),posF,[1,5,2,4,6])
    n = posF[end]
    m = length(glu[1][1])
    swp = @. swap2(glu[1],n,m)
    tensorN = swp, glu[2]
    ans = tensorSum(tensorN,[m])
    return ans
end

# For 3-1 move, we will use SVD -- i.e. use the function fullSplitTet3D or fullSplit3D

function permuteInd(vec,perm)
    #perm = [5,1,11,2,6,3,10,12,7,8,4,9]
    nvec = vec[perm]
    return nvec
end


function tensorPermute(tensor,perm)
    indx = tensor[1]
    indxN = []
    @assert length(indx[1]) == length(perm)
    for i in indx
        push!(indxN,i[perm])
    end
    return indxN , tensor[2]
end

end # module
