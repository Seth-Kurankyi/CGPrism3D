module CGPrism3D
using LinearAlgebra
#using StatsBase

# Define chop function
export numchop
numchop(num::Real)=(abs(num) >= 1000*eps() ? num : zero(Real))
function numchop(num::Complex)
    numchop(imag(num)) == 0 ? numchop(real(num)) : Complex(numchop(real(num)), numchop(imag(num)))
end
function numchop(x::Array)
    for i in 1:length(x)
        if x[i] != 0
            x[i] = numchop(x[i])
        end
    end
    return x
end

#Define all functions needed for G-symbol
# Define delta{ijk} -> coupling rules
function delta(i::Float64,j::Float64,k::Float64,K)
    sol = 0
    if i <=(j+k) && j<=(i+k) && k<=(i+j) && i+j+k <= K && 2*(i+j+k)%2 ==0
        sol = 1
    end
    return sol
end

#Define quantum numbers qn (this is real)
function qn(n::Float64,K)
    sol = (exp(pi*n*im/(K+2)) - exp(-pi*n*im/(K+2))) / (exp(pi*im/(K+2)) - exp(-pi*im/(K+2)))
    return real(sol)
end

#Define qn factorial
function qnfact(n::Float64,K)
    sol = 1
    for i in 1:n
        sol *= qn(i,K)
    end
    return sol
end

#Define square root of quantum dimension
function visqrt(i::Float64,K)
    sol = ((-1+0im)^i )*sqrt(qn(2*i+1,K))
    return sol
end

#Define triangle equality
function trian(i::Float64,j::Float64,k::Float64,K)
    sol = 0
    if delta(i,j,k,K) == 1
        sol = delta(i,j,k,K)*sqrt(qnfact(i+j-k,K)*qnfact(i-j+k,K)*qnfact(-i+j+k,K)/qnfact(i+j+k+1,K))
    end
    return sol
end

# Define Racah-Wigner six-j symbol
function RacWig6j(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64,K)
    a = i+j+m; b = i+l+n; c = k+j+n; d = k+l+m;  e = i+j+k+l; f = i+k+m+n; g = j+l+m+n
    sol = 0
    if delta(i,j,m,K) != 0 && delta(i,l,n,K) != 0 && delta(k,j,n,K) != 0 && delta(k,l,m,K) != 0
        sumz = 0
        for z in max(a,b,c,d):min(e,f,g)
            sumz += (-1)^z *qnfact(z+1,K)/
                ((qnfact(e-z,K)*qnfact(f-z,K)*qnfact(g-z,K))* (qnfact(z-a,K)*qnfact(z-b,K)*qnfact(z-c,K)*qnfact(z-d,K)))
        end
        sol = trian(i,j,m,K)*trian(i,l,n,K)*trian(k,j,n,K)*trian(k,l,m,K)*sumz
    end
    return sol
end

#Define F-symbol
function Fsymb(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64,K)
    sol = 0
    if delta(i,j,m,K) != 0 && delta(i,l,n,K) != 0 && delta(k,j,n,K) != 0 && delta(k,l,m,K) != 0
        sol = (-1+0im)^(i+j+k+l)*sqrt(qn(2*m+1,K)*qn(2*n+1,K)) * RacWig6j(i,j,m,k,l,n,K)
    end
    return sol
end
export Gsymb
#Define G-symbol
function Gsymb(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64,K)
    sol = 0
    if delta(i,j,m,K) != 0 && delta(i,l,n,K) != 0 && delta(k,j,n,K) != 0 && delta(k,l,m,K) != 0
        sol = Fsymb(i,j,m,k,l,n,K) /(visqrt(m,K)*visqrt(n,K))
    end
    return sol
end

export Tetra6j, delta, visqrt
#Define G-symbol times dimensions
function Tetra6j(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64,K)
    sol = 0
    if delta(i,j,m,K) != 0 && delta(i,l,n,K) != 0 && delta(k,j,n,K) != 0 && delta(k,l,m,K) != 0
        dims = prod(visqrt.([i,j,m,k,l,n],K))
        sol =  numchop(dims*Gsymb(i,j,m,k,l,n,K))
    end
    return sol
end

function dataFsymb(K)
    t6j = []
    t6s = []
    x = K+1
    y = K/2
    for ja in 0.:0.5:y, jb in 0.:0.5:y, jc in 0.:0.5:y
        if delta(ja,jb,jc,K) != 0
            for jd in 0.:0.5:y, je in 0.:0.5:y
                if delta(jc,jd,je,K) != 0
                    for jf in 0.:0.5:y
                        if delta(ja,je,jf,K) != 0 && delta(jb,jd,jf,K) != 0
            #dims = prod(visqrt.([ja,jb,jc,jd,je,jf]))
                            sol = numchop(Fsymb(ja,jb,jc,jd,je,jf,K) )
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
function prismA(jd1::Float64,jd2::Float64,je1::Float64,jb1::Float64,jg::Float64,ja1::Float64,jc::Float64,
                jb2::Float64,jf1::Float64,ja2::Float64,jf2::Float64,je2::Float64,K)
    #Prism A as gluing together three tetrahedra ;  Prism B changing diagonals using F-symbol
    sol = 0
    if delta(jd1,jd2,je1,K) != 0 && delta(jd1,jb1,jg,K) != 0 && delta(ja1,jd2,jg,K) != 0 && delta(ja1,jb1,je1,K) != 0 && delta(jd1,jc,jb2,K) != 0 && delta(jf1,jd2,jb2,K) != 0 && delta(jf1,jc,je1,K) != 0 && delta(jd1,ja2,jf2,K) != 0 && delta(jc,je2,jf2,K) != 0 && delta(jb2,je2,ja2,K) != 0
        dims = visqrt(ja1,K)*visqrt(jb1,K)*visqrt(jg,K)*visqrt(jf1,K)*visqrt(jc,K)*visqrt(jb2,K)*visqrt(je2,K)*visqrt(ja2,K)*visqrt(jf2,K)*visqrt(jd1,K)*visqrt(jd2,K)*visqrt(je1,K)
        sol =  dims*Gsymb(jd1,jd2,je1,ja1,jb1,jg,K) * Gsymb(jd1,jd2,je1,jf1,jc,jb2,K) * Gsymb(jd1,jc,jb2,je2,ja2,jf2,K)
    end
    return sol
end
function dataPrsmA(K)
    ampsInfo = Array{Float64,1}[]
    amps = Float64[]
    x = K+1
    y = K/2
    for jd1 in 0.:0.5:y, jd2 in 0.:0.5:y, je1 in 0.:0.5:y
        if delta(jd1,jd2,je1,K) != 0
            for jb1 in 0.:0.5:y, jg in 0.:0.5:y
                if delta(jd1,jb1,jg,K) != 0
                    for ja1 in 0.:0.5:y
                        if delta(ja1,jd2,jg,K) != 0 && delta(ja1,jb1,je1,K) != 0
                            for jc in 0.:0.5:y, jb2 in 0.:0.5:y
                                if delta(jd1,jc,jb2,K) != 0
                                    for jf1 in 0.:0.5:y
                                        if delta(jf1,jd2,jb2,K) != 0 && delta(jf1,jc,je1,K) != 0
                                            for ja2 in 0.:0.5:y, jf2 in 0.:0.5:y
                                                if delta(jd1,ja2,jf2,K) != 0
                                                    for je2 in 0.:0.5:y
                                                        if delta(jc,je2,jf2,K) != 0 && delta(jb2,je2,ja2,K) != 0 #
                                                            sol = numchop(prismA(jd1,jd2,je1,jb1,jg,ja1,jc,jb2,jf1,ja2,jf2,je2,K))
                                                            if sol != 0
                                                                indx = [jb1,jg,ja1,jd1,jd2,je1,jc,jb2,jf1,ja2,jf2,je2]
                                                                push!(ampsInfo,indx)
                                                                push!(amps,sol)
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return ampsInfo, amps
end

### General function for block for svd, for gluing two tensor and for summing over an index
function tensorBlock(tensorM,posnA,posnB,posnC,spinC)
    # a general function to split a tensor M^C_{AB} to U^C_{Ai} and  V^C_{iB} with C the shared face.. for particular values of C
    # first reduce to matrix form and then perform svd to get U, s, and V
    #This reduces the tensor M into matrix form M^C_{AB} with AB labelling matric indices.
    indx = tensorM[1]
    ampvals = tensorM[2]
    fc = findall(x-> x[posnC] == spinC, indx) #find spinsC in all spins length(posnC) = length(spinC)
    ampInfo = indx[fc] # this contains all spin configuration containing ampsC
    amps = ampvals[fc] # this contains the corresponding spin amplitude values for ampsC
    Acol = unique(getindex.(ampInfo, [posnA]))
    Brow = unique(getindex.(ampInfo, [posnB]))
    indxU = vcat.(Acol,repeat([spinC],length(Acol)))
    indxV = vcat.(Brow,repeat([spinC],length(Brow)))
    matM = zeros(length(Acol),length(Brow))
    if length(matM) != 0
        for i in 1:length(ampInfo)
            qu = findfirst(x-> x == ampInfo[i][posnA], Acol)[1]# [1] since it returns an array [1]
            qv = findfirst(x-> x == ampInfo[i][posnB], Brow)[1]
            matM[qu,qv] = amps[i]
        end
        #U, s, V = svd(matM)
        return matM,indxU,indxV
    else
        return 0,0,0
    end
end

function tensorGlue(tensorB,tensorA,posnB,posnA)
    # glue two tensors TA,TB along posnA from TA and posnB from TB
    indxA = tensorA[1] # spins of TA
    indxB = tensorB[1] # spins of TB
    ampsA = tensorA[2] # amplitudes of TA
    ampsB = tensorB[2] # amplitudes of TB
    shrdA = getindex.(indxA,[posnA])
    shrdB = getindex.(indxB,[posnB])
    lena = collect(1:length(indxA[1]))
    deleteat!(lena,sort(posnA)) # assert needs posnA to be a vector
    #if length(posnA) > 1 deleteat!(lena,sort(posnA,dims=2) else  deleteat!(lena,sort(posnA)) end
    shrdAB = []
    qq = []
    amps = []
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
    if shrdAu == shrdBu
        shrd = shrdAu
    else
        for i in shrdBu
            shrd = shrdAu
            if !(i in shrdAu)
                push!(shrd,i)
            end
        end
    end
    #cshrdA = countmap(shrdA)
    #cshrdB = countmap(shrdB)
    for i in shrd
        ca = count(x->x==i,shrdA)
        cb = count(x->x==i,shrdB)
        push!(shrdAB, (i, ca,cb ))
    end
    cnt = 1
    cntb = 1
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

function tensorSum(ampJ,posnN)
    # Summing a tensor  over indices labelled posnN
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



# For 3-1 move, we will use SVD -- i.e. use the function fullSplitTet3D or fullSplit3D


export dataTet
function dataTet(K)
    t6j = []
    t6s = []
    x = K+1
    y = K/2
    for ja in 0.:0.5:y, jb in 0.:0.5:y, jc in 0.:0.5:y
        if delta(ja,jb,jc,K) != 0
            for jd in 0.:0.5:y, je in 0.:0.5:y
                if delta(jc,jd,je,K) != 0
                    for jf in 0.:0.5:y
                        if delta(ja,je,jf,K) != 0 && delta(jb,jd,jf,K) != 0
            #dims = prod(visqrt.([ja,jb,jc,jd,je,jf]))
                            sol = numchop(Tetra6j(ja,jb,jc,jd,je,jf,K) )
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
export sumT
function sumT(Kdata,Kdistr)
    n = length(Kdata)
    m = length(Kdistr)
    #@assert n!=m "Kdata and Kdistr must have the same length"

    dataCreation = []
    for i in 1:n
        push!(dataCreation,dataTet(Kdata[i]))
    end

    posmax = argmax(Kdata)
    Kmax = maximum(Kdata)

    setspin = dataCreation[posmax][1]
    setAmp = zeros(length(setspin))

    for i in 1:length(setspin)
        for j in 1:n
            pos = findfirst(x->x==setspin[i],dataCreation[j][1])
            if typeof(pos) == Int
                setAmp[i] += Kdistr[j]*dataCreation[j][2][pos]
            end
        end
    end

    return setspin,setAmp
end

function fullSplitTet3D(dataM,posnA,posnB,posnC,K,visqrt,Tetra6j,DimData) # posnA,B,C must be vectors @assert length(dataM) = vcat(A,B,C)
    #A general function to splits fully the tensor M^C_{AB} to tensors U^C_{Ai} and V^C_{iB}.
    # It uses tensorBlock, 'fixes' the sign problem from svd, and normalize
    # We keep only the first singular value (assuming a geometric split)
    # splits prism into tetrahedron and a  pyramid
    # Multiply back by missing dimension factors
    # Level K specification only here because of the fixing of the sign
    indx = dataM[1]
    ampvals = dataM[2]
    ampsC = unique(getindex.(indx, [posnC])) # get all unique spins C (jd1,jd2,je1)
    ampsU = [] # this will store all amplitudes for U^C_{A}
    ampsV = [] # this will store all amplitudes for V^C_{B}
    indxsU = []
    indxsV = []
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
            trunc = 1
            U1 = U[:,1:trunc]
            V1 = V[:,1:trunc]
            s1 = s[1:trunc]
            valU = U1*sqrt.(s1)*sqrt(prod(visqrt(DimData,ampsC[i][j]) for j in 1:length(ampsC[i]))) # multiply by leftover dimension factors
            valV = V1*sqrt.(s1)*sqrt(prod(visqrt(DimData,ampsC[i][j]) for j in 1:length(ampsC[i])))
            # valU,valV are all either real or purely imaginary
            valU = real(valU) - imag(valU) # take -ve of imaginary part if it gives only imag
            valV = real(valV) + imag(valV)
            ja1 = blkU[1][1]; jb1 = blkU[1][2]; jg = blkU[1][3] # pick the first spins in U
            #Fixing sign problem -- keep the same sign for U and Tetra6j symbol
            if sign(valU[1]) == sign(Tetra6j(DimData,jd1,jd2,je1,ja1,jb1,jg))
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
export tensorGlueTet3D
function tensorGlueTet3D(tensorB,tensorA,posnB,posnA,visqrt,DimData)
    #Take care of dimesion factors after gluing
    indx, amps = tensorGlue(tensorB,tensorA,posnB,posnA)
    face = getindex.(indx, [posnB])
    ans = complex(ones(length(face)))
    for i in 1:length(face)
        ans[i] = prod(visqrt(DimData,face[i][j]) for j in 1:length(face[i]))
    end
    ampsN = @. amps / ans
    ampsN = real(ampsN)+imag(ampsN)
    return indx, ampsN
end
function tensorSumTet3D(tensor,posnN,visqrt,DimData)
    # Take care of dimension factors when summing. This makes sure Pachner 3-2 move works
    # For now this only works for summing over one index-- we can change this easily
    indx = tensor[1]
    amps = tensor[2]
    blkdims = getindex.(indx,posnN)
    ampsN = amps .*visqrt(DimData,blkdims[1])
    ampsN2 = real(ampsN) - imag(ampsN)
    data = indx, ampsN2
    ans = tensorSum(data,posnN)
    return ans
end
function tensor22move(tensor,face,K,Fsymb)# place middle edge at last position, triangle are: 125 345
    # Performs an F-move on a face.. equivalent to 2-2 Pachner move -- this uses Fsymbol
    x = K+1
    y = K/2
    # assert face contains 5 elements
    indx = tensor[1] # spins of TB
    amps = tensor[2]
    flast = face[end]
    indxlast = getindex.(indx, [flast])
    famps = []
    findx = []
    for i in 1:length(indx)
        j1 = indx[i][face[1]]; j2 = indx[i][face[2]]; j3 = indx[i][face[3]]; j4 = indx[i][face[4]]; j5 = indx[i][face[5]]
        for j in 0.:0.5:y
            #iMod = deepcopy(indxf1u)
            if delta(j1,j3,j,K) != 0 && delta(j2,j4,j,K) != 0 && Fsymb(j1,j3,j,j4,j2,j5,K) != 0
                iMod = copy(indx[i])
                ampsA = numchop(amps[i] * Fsymb(j1,j3,j,j4,j2,j5,K))
                push!(famps,ampsA)
                iMod[flast] = j
                push!(findx,iMod)
            end
        end
    end
    tt = sortperm(findx) # sort them spins
    indxf1 = findx[tt] # apply sort on spins
    famps = famps[tt]  # apply sort on amps
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
    return indxeff, qq
end
##### FULL ALGO
export newTetra6j
function newTetra6j(DimData,i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64)
    #data: pos1: spin config, pos2: amplitude
    pos = findfirst(x->x == [i,j,m,k,l,n],DimData[1])
    sol = 0
    if typeof(pos) == Int
        sol = DimData[2][pos]
    else
        print([i,j,m,k,l,n])
        @show pos
    end
    return sol
end
function qdim(DimData,n::Float64)
    return abs(newTetra6j(DimData,n,n,0.,n,n,0.))
end
export newvisqrt
function newvisqrt(DimData,i::Float64)
    return ((-1+0im)^i )*sqrt(qdim(DimData,i))
end

#    perm = [8,1,6,3,5,2,7,9,12,11,10,4]
function permuteInd(vec,perm)
    nvec = vec[perm]
    return nvec
end
function tensorPermute(tensor,perm)
    indx = tensor[1]
    indxN = []
    for i in 1:length(indx)
        indxN =push!(indxN,permuteInd(indx[i],perm))
    end
    return indxN , tensor[2]
end

export Algo
function Algo(Kdata,Kdistr,loop::Int64,K2move)
    #initialization
    data = sumT(Kdata,Kdistr) #setspin,setAmp
    Kemb = maximum(Kdata)

    tet2 = tensorGlueTet3D(data,data,[1,2,3],[1,2,3],newvisqrt,data)
    prA = tensorGlueTet3D(tet2,data,[2,7,9],[1,2,3],newvisqrt,data)
    #####
    # for i in 1:loop
    u,v = fullSplitTet3D(prA,[10,11,12],[1,3,4,5,6,8],[2,7,9],Kemb,newvisqrt,newTetra6j,data)
    perm = [4,5,6,1,2,3]
    u= tensorPermute(u,perm)
    # data = u
    prAuv =tensorGlueTet3D(u,v,[3,4,5],[9,1,6],newvisqrt,data)
    #
    pruvEff = tensorGlueTet3D(prAuv,prAuv,[3,2,1,12,11],[4,2,6,9,10],newvisqrt,data)
    pruvEffSum = tensorSumTet3D(pruvEff,[3],newvisqrt,data)
    # #
    fA0F = tensor22move(pruvEffSum,[15,16,4,6,11],K2move,Fsymb)
    fA01F = tensor22move(fA0F,[16,17,7,9,10],K2move,Fsymb)
    fA012F = tensor22move(fA01F,[14,12,4,5,1],K2move,Fsymb)
    # #
    f133A = fullSplitTet3D(fA012F,[4,15,14],[2,3,5,6,7,8,9,10,12,13,16,17],[18,1,11],Kemb,newvisqrt,newTetra6j,data)
    CGPrsmA1 = fullSplitTet3D(f133A[2],[11,4,5],[1,2,3,7,9,10,12,13,14],[6,8,15],Kemb,newvisqrt,newTetra6j,data)
    perm = [9,2,6,1,5,3,7,8,12,11,10,4]
    prA = tensorPermute(CGPrsmA1[2],perm)
    # end

    return fA012F
end

# function permuteInd(vec)
#     perm = [5,1,11,2,6,3,10,12,7,8,4,9]
#     nvec = vec[perm]
#     return nvec
# end
# function tensorPermute(tensor)
#     indx = tensor[1]
#     indxN = permuteInd.(indx)
#     return indxN , tensor[2]
# end
# export AlgoSethnot
# function AlgoSethnot(Kdata,Kdistr,loop::Int64,K2move)
#     #initialization
#     data = sumT(Kdata,Kdistr) #setspin,setAmp,dimfact
#     Kemb = maximum(Kdata)
#
#     function qdim(n::Float64,Kemb)
#         return data[3][Int(2*n+1)]
#     end
#     function newvisqrt(i::Float64,Kemb)
#         return ((-1+0im)^i )*sqrt(qdim(i,Kemb))
#     end
#     function newTetra6j(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64,Kemb)
#         pos = findfirst(x->x == [i,j,m,k,l,n],data[1])
#         return data[2][pos]
#     end
#
#     #inital prism
#     tet2 = tensorGlueTet3D(data[1:2],data[1:2],[1,2,3],[1,2,3],Kemb,newvisqrt)
#     prA = tensorGlueTet3D(tet2,data[1:2],[2,7,9],[1,2,3],Kemb,newvisqrt)
#     ######
#
#     for i in 1:loop
#         pruvA = tensorGlueTet3D(prA,prA,[1,5,6,8,9],[2,4,6,12,11],Kemb,newvisqrt)
#         pruvAA = tensorSumTet3D(pruvA,[1],Kemb,newvisqrt)
#
#         fA0F = tensor22move(pruvAA,[18,15,6,2,7],Kemb,Fsymb)
#         fA01F = tensor22move(fA0F,[17,18,10,9,8],Kemb,Fsymb)
#         fA012F = tensor22move(fA01F,[14,13,3,2,4],Kemb,Fsymb)
#
#         f133A = fullSplitTet3D(fA012F,[9,6,18],[1,2,3,4,5,10,12,13,14,15,16,17],[7,8,11],Kemb,newvisqrt,newTetra6j)
#         CGPrsmA1 = fullSplitTet3D(f133A[2],[8,10,2],[1,3,5,6,7,9,12,14,15],[13,4,11],Kemb,newvisqrt,newTetra6j)
#
#         prA = tensorPermute(CGPrsmA1[2])
#     end
#
#     return prA
# end

end # module
