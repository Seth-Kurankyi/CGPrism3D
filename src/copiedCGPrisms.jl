module CGPrism3D
using LinearAlgebra
#using StatsBase

# Define chop function 
export numchop
numchop(num::Real)=(abs(num) >= 1000*eps() ? num : zero(Real))

function numchop(num::Complex)
    numchop(imag(num)) == 0 ? numchop(real(num)) : Complex(numchop(real(num)), numchop(imag(num)))
end

# Define k-level of SU(2)_k
const K = 3 # K-level
const x = K+1 #number of spins
const y = K/2 # maximum spin

#Define all functions needed for G-symbol
# Define delta{ijk} -> coupling rules
function delta(i::Float64,j::Float64,k::Float64)
    sol = 0
    if i <=(j+k) && j<=(i+k) && k<=(i+j) && i+j+k <= K && 2*(i+j+k)%2 ==0
        sol = 1
    end
    return sol
end

#Define quantum numbers qn (this is real)
function qn(n::Float64)
    sol = (exp(pi*n*im/(K+2)) - exp(-pi*n*im/(K+2))) / (exp(pi*im/(K+2)) - exp(-pi*im/(K+2))) 
    return real(sol)
end

#Define qn factorial
function qnfact(n::Float64)
    sol = 1
    for i in 1:n
        sol *= qn(i)
    end
    return sol
end

#Define square root of quantum dimension
function visqrt(i::Float64)
    sol = ((-1+0im)^i )*sqrt(qn(2*i+1))
    return sol
end

#Define triangle equality 
function trian(i::Float64,j::Float64,k::Float64)
    sol = 0
    if delta(i,j,k) == 1
        sol = delta(i,j,k)*sqrt(qnfact(i+j-k)*qnfact(i-j+k)*qnfact(-i+j+k)/qnfact(i+j+k+1))
    end
    return sol
end

# Define Racah-Wigner six-j symbol
function RacWig6j(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64)
    a = i+j+m; b = i+l+n; c = k+j+n; d = k+l+m;  e = i+j+k+l; f = i+k+m+n; g = j+l+m+n
    sol = 0
    if delta(i,j,m) != 0 && delta(i,l,n) != 0 && delta(k,j,n) != 0 && delta(k,l,m) != 0  
        sumz = 0
        for z in max(a,b,c,d):min(e,f,g)
            sumz += (-1)^z *qnfact(z+1)/
                ((qnfact(e-z)*qnfact(f-z)*qnfact(g-z))* (qnfact(z-a)*qnfact(z-b)*qnfact(z-c)*qnfact(z-d)))
        end
        sol = trian(i,j,m)*trian(i,l,n)*trian(k,j,n)*trian(k,l,m)*sumz
    end
    return sol
end

#Define F-symbol
function Fsymb(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64)
    sol = 0
    if delta(i,j,m) != 0 && delta(i,l,n) != 0 && delta(k,j,n) != 0 && delta(k,l,m) != 0
        sol = (-1+0im)^(i+j+k+l)*sqrt(qn(2*m+1)*qn(2*n+1)) * RacWig6j(i,j,m,k,l,n)
    end
    return sol
end

#Define G-symbol
function Gsymb(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64)
    sol = 0
    if delta(i,j,m) != 0 && delta(i,l,n) != 0 && delta(k,j,n) != 0 && delta(k,l,m) != 0
        sol = Fsymb(i,j,m,k,l,n) /(visqrt(m)*visqrt(n))
    end
    return sol
end

export Tetra6j, delta, visqrt
#Define G-symbol times dimensions
function Tetra6j(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64)
    sol = 0
    if delta(i,j,m) != 0 && delta(i,l,n) != 0 && delta(k,j,n) != 0 && delta(k,l,m) != 0
        dims = prod(visqrt.([i,j,m,k,l,n]))
        sol =  numchop(dims*Gsymb(i,j,m,k,l,n))
    end
    return sol
end

#Define G-symbol 2
function Gsymb2(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64)
    sol = 0 
    if delta(i,j,m) != 0 && delta(i,l,n) != 0 && delta(k,j,n) != 0 && delta(k,l,m) != 0
        sol = (-1+0im)^(i+j+k+l+m+n)*RacWig6j(i,j,m,k,l,n)
    end
    return sol
end


## Define Prisms 
#Prism A as gluing together three tetrahedra ;  Prism B changing diagonals using F-symbol
function prismA(jd1::Float64,jd2::Float64,je1::Float64,jb1::Float64,jg::Float64,ja1::Float64,jc::Float64,
                jb2::Float64,jf1::Float64,ja2::Float64,jf2::Float64,je2::Float64)
    sol = 0 
    if delta(jd1,jd2,je1) != 0 && delta(jd1,jb1,jg) != 0 && delta(ja1,jd2,jg) != 0 && delta(ja1,jb1,je1) != 0 && delta(jd1,jc,jb2) != 0 && delta(jf1,jd2,jb2) != 0 && delta(jf1,jc,je1) != 0 && delta(jd1,ja2,jf2) != 0 && delta(jc,je2,jf2) != 0 && delta(jb2,je2,ja2) != 0  
        dims = visqrt(ja1)*visqrt(jb1)*visqrt(jg)*visqrt(jf1)*visqrt(jc)*visqrt(jb2)*visqrt(je2)*visqrt(ja2)*visqrt(jf2)*visqrt(jd1)*visqrt(jd2)*visqrt(je1)
        sol =  dims*Gsymb(jd1,jd2,je1,ja1,jb1,jg) * Gsymb(jd1,jd2,je1,jf1,jc,jb2) * Gsymb(jd1,jc,jb2,je2,ja2,jf2)
    end
    return sol
end

function prismB(jd2p::Float64,jd1p::Float64,je2::Float64,jb2::Float64,jg::Float64,ja2::Float64,jc::Float64,
                jb1::Float64,jf2::Float64,ja1::Float64,jf1::Float64,je1::Float64)
    sol = 0
    if delta(jd1p,jd2p,je2) != 0 && delta(jd2p,jb2,jg) != 0 && delta(ja2,jd1p,jg) != 0 && delta(ja2,jb2,je2) != 0 && delta(jd2p,jc,jb1) != 0 && delta(jf2,jd1p,jb1) != 0 && delta(jf2,jc,je2) != 0 && delta(jd2p,ja1,jf1) != 0 && delta(jc,je1,jf1) != 0 && delta(jb1,je1,ja1) != 0
        for jd1 in 0.:0.5:y
            if delta(jb1,jg,jd1) != 0 && delta(ja2,jf2,jd1) != 0
                for jd2 in 0.:0.5:y
                    if delta(ja1,jg,jd2) != 0 && delta(jb2,jf1,jd2) != 0
                        sol += numchop(Fsymb(ja1,jg,jd2,jb2,jf1,jd2p)*Fsymb(jb1,jg,jd1,ja2,jf2,jd1p)*prismA(jd1,jd2,je1,jb1,jg,ja1,jc,jb2,jf1,ja2,jf2,je2))
                    end
                end
            end
        end
    end
    return sol
end

# This is the relabelling vertices. Does it work in general after one step in algorithm?
function prismB2(jd2::Float64,jd1::Float64,je2::Float64,jb2::Float64,jg::Float64,ja2::Float64,jc::Float64,
                jb1::Float64,jf2::Float64,ja1::Float64,jf1::Float64,je1::Float64)
    return prismA(jd1,jd2,je1,jb1,jg,ja1,jc,jb2,jf1,ja2,jf2,je2)
end



export dataPrsmA, dataPrsmB, dataPrsmB2, dataTet, dataTetF

# get all non-zero amplitudes and spins 
function dataPrsmA()
    ampsInfo = Array{Float64,1}[]
    amps = Float64[]
    for jd1 in 0.:0.5:y, jd2 in 0.:0.5:y, je1 in 0.:0.5:y
        if delta(jd1,jd2,je1) != 0 
            for jb1 in 0.:0.5:y, jg in 0.:0.5:y
                if delta(jd1,jb1,jg) != 0 
                    for ja1 in 0.:0.5:y
                        if delta(ja1,jd2,jg) != 0 && delta(ja1,jb1,je1) != 0
                            for jc in 0.:0.5:y, jb2 in 0.:0.5:y
                                if delta(jd1,jc,jb2) != 0
                                    for jf1 in 0.:0.5:y
                                        if delta(jf1,jd2,jb2) != 0 && delta(jf1,jc,je1) != 0
                                            for ja2 in 0.:0.5:y, jf2 in 0.:0.5:y
                                                if delta(jd1,ja2,jf2) != 0
                                                    for je2 in 0.:0.5:y
                                                        if delta(jc,je2,jf2) != 0 && delta(jb2,je2,ja2) != 0 #
                                                            sol = numchop(prismA(jd1,jd2,je1,jb1,jg,ja1,jc,jb2,jf1,ja2,jf2,je2))
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


function dataPrsmB()
    ampsInfo = Array{Float64,1}[]
    amps = Float64[]
    for jd2 in 0.:0.5:y, jd1 in 0.:0.5:y, je2 in 0.:0.5:y
        if delta(jd1,jd2,je2) != 0 
            for jb2 in 0.:0.5:y, jg in 0.:0.5:y
                if delta(jd2,jb2,jg) != 0 
                    for ja2 in 0.:0.5:y
                        if delta(ja2,jd1,jg) != 0 && delta(ja2,jb2,je2) != 0
                            for jc in 0.:0.5:y, jb1 in 0.:0.5:y
                                if delta(jd2,jc,jb1) != 0
                                    for jf2 in 0.:0.5:y
                                        if delta(jf2,jd1,jb1) != 0 && delta(jf2,jc,je2) != 0
                                            for ja1 in 0.:0.5:y, jf1 in 0.:0.5:y
                                                if delta(jd2,ja1,jf1) != 0
                                                    for je1 in 0.:0.5:y
                                                        if delta(jc,je1,jf1) != 0 && delta(jb1,je1,ja1) != 0 #
                                                            sol = numchop(prismB(jd2,jd1,je2,jb2,jg,ja2,jc,jb1,jf2,ja1,jf1,je1))
                                                            if sol != 0
                                                                indx = [jb2,jg,ja2,jd2,jd1,je2,jc,jb1,jf2,ja1,jf1,je1]
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


function dataPrsmB2()
    ampsInfo = Array{Float64,1}[]
    amps = Float64[]
    for jd2 in 0.:0.5:y, jd1 in 0.:0.5:y, je2 in 0.:0.5:y
        if delta(jd2,jd1,je2) != 0 
            for jb2 in 0.:0.5:y, jg in 0.:0.5:y
                if delta(jd2,jb2,jg) != 0 
                    for ja2 in 0.:0.5:y
                        if delta(ja2,jd1,jg) != 0 && delta(ja2,jb2,je2) != 0
                            for jc in 0.:0.5:y, jb1 in 0.:0.5:y
                                if delta(jd2,jc,jb1) != 0
                                    for jf2 in 0.:0.5:y
                                        if delta(jf2,jd1,jb1) != 0 && delta(jf2,jc,je2) != 0
                                            for ja1 in 0.:0.5:y, jf1 in 0.:0.5:y
                                                if delta(jd2,ja1,jf1) != 0
                                                    for je1 in 0.:0.5:y
                                                        if delta(jc,je1,jf1) != 0 && delta(jb1,je1,ja1) != 0 #
                                                            sol = numchop(prismB2(jd2,jd1,je2,jb2,jg,ja2,jc,jb1,jf2,ja1,jf1,je1))
                                                            if sol != 0
                                                                indx = [jb2,jg,ja2,jd2,jd1,je2,jc,jb1,jf2,ja1,jf1,je1]
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
    Indxspins = sort(ampsInfo)
    ampvals = amps[sortperm(ampsInfo)]
    return Indxspins, ampvals
end

function dataTet()
    t6j = []
    t6s = []
    for ja in 0.:0.5:y, jb in 0.:0.5:y, jc in 0.:0.5:y
        if delta(ja,jb,jc) != 0
            for jd in 0.:0.5:y, je in 0.:0.5:y
                if delta(jc,jd,je) != 0
                    for jf in 0.:0.5:y
                        if delta(ja,je,jf) != 0 && delta(jb,jd,jf) != 0
            #dims = prod(visqrt.([ja,jb,jc,jd,je,jf]))
                            sol = numchop(Tetra6j(ja,jb,jc,jd,je,jf) )
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

function dataTetF()
    t6j = []
    t6s = []
    for ja in 0.:0.5:y, jb in 0.:0.5:y, jc in 0.:0.5:y
        if delta(ja,jb,jc) != 0
            for jd in 0.:0.5:y, je in 0.:0.5:y
                if delta(jc,jd,je) != 0
                    for jf in 0.:0.5:y
                        if delta(ja,je,jf) != 0 && delta(jb,jd,jf) != 0
            #dims = prod(visqrt.([ja,jb,jc,jd,je,jf]))
                            sol = numchop(Fsymb(ja,jb,jc,jd,je,jf) )
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

export tensorSplitSVD, fullSplitPrism3D, tensorGlue, tensorSum, tensorGluePrism3D, tensor22move, tensor31move

# a general function to split a tensor M^C_{AB} to U^C_{Ai} and  V^C_{iB} with C the shared face.. for particular values of C
# first reduce to matrix form and then perform svd to get U, s, and V
#This returns the svd of tensorM,i.e., U^C_{Ai}, s_i, and V^C_{iB}. 
function tensorSplitSVD(tensorM,posnA,posnB,posnC,spinC)
    indx = tensorM[1]
    ampvals = tensorM[2]
    fc = findall(x-> x[posnC] == spinC, indx) #find spinsC in all spins length(posnC) = length(spinC)
    ampInfo = indx[fc] # this contains all spin configuration containing ampsC
    amps = ampvals[fc] # this contains the corresponding spin amplitude values for ampsC
    Ucol = unique(getindex.(ampInfo, [posnA]))
    Vrow = unique(getindex.(ampInfo, [posnB]))
    matM = zeros(length(Ucol),length(Vrow))
    if length(matM) != 0
        for i in 1:length(ampInfo)
            qu = findall(x-> x == ampInfo[i][posnA], Ucol)[1]# [1] since it returns an array [1]
            qv = findall(x-> x == ampInfo[i][posnB], Vrow)[1]
            matM[qu,qv] = amps[i]
        end
        U, s, V = svd(matM)
        return U,s,V
    else
        return 0,0,0
    end
end

#A general function to splits fully the tensor M^C_{AB} to tensors U^C_{Ai} and V^C_{iB}. 
# It uses tensorSplitSVD,'fixes' the sign problem from svd, and normalize 
#This is designed specifically for the prisms we are working with. It keeps only the first singular value (assuming a gemetric split)
# splits prism into tetrahedron and a  pyramid
function fullSplitPrism3D(dataM,posnA,posnB,posnC) # posnA,B,C must be vectors @assert length(dataM) = vcat(A,B,C)
    indx = dataM[1] 
    ampvals = dataM[2]
    posAC = vcat(posnA,posnC) # assert posnA,posnC must be both vectors
    posCB = vcat(posnC,posnB) # assert posnB,posnC must be both vectors
    indxU = unique(getindex.(indx, [posAC]))
    indxV = unique(getindex.(indx, [posCB]))
    ampsC = unique(getindex.(indx, [posnC])) # get all unique spins C (jd1,jd2,je1) 
    ampsU = [] # this will store all amplitudes for U^C_{A}
    ampsV = [] # this will store all amplitudes for V^C_{B}
    #indxUV = [] # if we want this will store all spins in the right order
    for i in 1:length(ampsC) # loop over the unique spins
        jd1 = ampsC[i][1]; jd2 = ampsC[i][2]; je1 = ampsC[i][3]  #jd1,jd2,je1 in spin C form a triangle
        if delta(jd1,jd2,je1) != 0
            spinC = [jd1,jd2,je1]
            U, s, V = tensorSplitSVD(dataM,posnA,posnB,posnC,spinC)
            #truncate and use only first singular value
            U1 = U[:,1] 
            V1 = V[:,1]
            s1 = s[1]
            valU = U1*sqrt(s1)*sqrt(prod(visqrt.(ampsC[i]))) # multiply by leftover dimension factors
            valV = V1*sqrt(s1)*sqrt(prod(visqrt.(ampsC[i])))
            # valU,valV are all either real or purely imaginary
            valU = real(valU) - imag(valU) # take -ve of imaginary part if it gives only imag
            valV = real(valV) + imag(valV)
            ja1 = indxU[1][1]; jb1 = indxU[1][2]; jg = indxU[1][3] # pick the first spins in U
            #Fixing sign problem -- keep the same sign for U and Tetra6j symbol 
            if sign(valU[1]) == sign(Tetra6j(jd1,jd2,je1,jb1,jg,ja1) )
                push!(ampsU, valU)#/valU[1])
                push!(ampsV, valV)#/valV[1])
            else
                push!(ampsU, -valU)#/valU[1])
                push!(ampsV, -valV)#/valV[1])
            end    
        end
    end
    solU = collect(Iterators.flatten(ampsU))
    solV = collect(Iterators.flatten(ampsV))
    ansU = solU/solU[1] # normalization condition
    ansV = solV/solV[1]
    return (indxU,ansU),(indxV,ansV) #,indxUV
end

#its better if TA the 'smallest' tensor, posnA and posnB should be of the same length
function tensorGlue(tensorB,tensorA,posnB,posnA)# glue two tensors TA,TB along posnA from TA and posnB from TB 
    indxA = tensorA[1] # spins of TA
    indxB = tensorB[1] # spins of TB
    ampsA = tensorA[2] # amplitudes of TA
    ampsB = tensorB[2] # amplitudes of TB
    lena = collect(1:length(indxA[1]))
    deleteat!(lena,sort(posnA)) # assert needs posnA to be a vector
    #if length(posnA) > 1 deleteat!(lena,sort(posnA,dims=2) else  deleteat!(lena,sort(posnA)) end
    amps = []
    qq = []
    for i in 1:length(indxA)
        fc = findall(x-> x[posnB] == indxA[i][posnA], indxB) # needs length(posnA) = length(posnB)
        ans = @. ampsA[i]*ampsB[fc]
        push!(amps,ans)
        indxa = repeat([indxA[i][lena]],length(indxB[fc])) 
        #for j in indxA[i][lena]
        indxa1 = vcat.(indxB[fc],indxa)
        push!(qq,indxa1)
        #push!(qq,[collect(Iterators.flatten([indxB[fc][j],indxA[i][lena]])) for j in 1:length(indxB[fc])])
    end
    sol = collect(Iterators.flatten(amps))
    indxsol = collect(Iterators.flatten(qq))
    return indxsol, sol
end
export tensorGlueN

function tensorGlueN(tensorB,tensorA,posnB,posnA)# glue two tensors TA,TB along posnA from TA and posnB from TB 
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
    cshrdA = countmap(shrdA)
    cshrdB = countmap(shrdB)
    for i in shrd
        ca = count(x->x==i,shrdA)
        cb = count(x->x==i,shrdB)
        push!(shrdAB, (i, ca,cb ))
    end
    cnt = 1
    cntb = 1
    for i in shrdAB
        # Put condition for when one of a or b is zero
        for j in 1:i[2]
            ans = @. sampsA[j]*sampsB[1:i[3]]
            push!(amps,ans)
            indxa = repeat([sindxA[j][lena]],i[3]) 
            #for j in indxA[i][lena]
            indxa1 = vcat.(sindxB[1:i[3]],indxa)
            push!(qq,indxa1)
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

export visqrt, tensorGluePrism3D2, tensorGluePrism3DN

function tensorGluePrism3D(tensorB,tensorA,posnB,posnA)
    indx, amps = tensorGlue(tensorB,tensorA,posnB,posnA)
    face = getindex.(indx, [posnB])
    ans = complex(zeros(length(face)))
    for i in 1:length(face)
        ans[i] = prod(visqrt.(face[i]))
    end
    ampsN = @. amps / ans
    ampsN = real(ampsN)+imag(ampsN)
    return indx, ampsN
end

function tensorGluePrism3DN(tensorB,tensorA,posnB,posnA)
    indx, amps = tensorGlueN(tensorB,tensorA,posnB,posnA)
    face = getindex.(indx, [posnB])
    ans = complex(zeros(length(face)))
    for i in 1:length(face)
        ans[i] = prod(visqrt.(face[i]))
    end
    ampsN = @. amps / ans
    ampsN = real(ampsN)+imag(ampsN)
    return indx, ampsN
end

## This is to try to take care of bulk spins 
function tensorGluePrism3D2(tensorB,tensorA,posnB,posnA)
    indx, amps = tensorGlue(tensorB,tensorA,posnB,posnA)
    face = getindex.(indx, [posnB])
    ans = complex(ones(length(face)))
    for i in 2:length(face)
        ans[i] = prod(visqrt.(face[i]))
    end
    ampsN = @. amps / ans
    ampsN = real(ampsN)+imag(ampsN)
    return indx, ampsN
end

function tensorSum(ampJ,posnN)
    indx = ampJ[1]
    amps = ampJ[2]
    lenS = collect(1:length(indx[1]))
    deleteat!(lenS,sort(posnN))
    #if length(posnN) > 1   deleteat!(lenS,sort(posnN,dims=2)) else  deleteat!(lenS,sort(posnN)) end
    indxef = getindex.(indx,[lenS]) # get all spins without N 
    indxefu = unique(indxef)
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

function tensor22move(tensor,face)# place middle edge at last position # first 2 position form a triangle
    # assert face contains 5 elements
    indx = tensor[1] # spins of TB
    amps = tensor[2]
    flast = face[end]
    ff1 = collect(1:length(indx[1]))
    deleteat!(ff1,flast)
    indxf1 = getindex.(indx, [ff1])
    indxf1u = unique(indxf1)
    #indxf1uMod = unique(indxf1)
    indxlast = getindex.(indx, [flast])
    for i in 1:4
        if face[i]>flast
            face[i] = face[i]-1
        end
    end
    sol = []
    findx = []
    for i in 1:length(indxf1u)
        iMod = copy(indxf1u[i])
        j1 = iMod[face[1]]; j2 = iMod[face[2]]; j3 = iMod[face[3]]; j4 = iMod[face[4]]; #j5 = ampsN[i][face[5]]
    #    #ans = 0#amps[i]
    #    ind = zeros(length(indx[1]))
        for j in 0.:0.5:y
            #iMod = deepcopy(indxf1u)
            if delta(j1,j3,j) != 0 && delta(j2,j4,j) != 0
                fc = findall(x-> x == iMod,indxf1 )
                ans = 0
                for t in fc
                    ju = indxlast[t]
                    ans += numchop(amps[t] * Fsymb(j1,j3,j,j4,j2,ju)) #* ampls[t] 
                end
                if numchop(ans) !=0
                    push!(sol,ans)
                    iMods = copy(indxf1u[i])
                    push!(findx,insert!(iMods,flast,j))
                end
    #            end
            end
        end
    end
    return findx,sol#indx,amps
end


function tensor31move(tensor,tet)
    # assert tet contains 6 elements
    indx = tensor[1] # spins of TB
    amps = tensor[2] # amps
    indxout = getindex.(indx, [tet[4:6]])
    ff1 = collect(1:length(indx[1]))
    deleteat!(ff1,sort(tet[4:6]))
    indxkeep = getindex.(indx, [ff1])
    indxkeepu = unique(getindex.(indx, [ff1]))
    sol = []
    indxn = []
    for i in 1:length(indxkeepu)
        ff = findfirst(x-> x == indxkeepu[i], indxkeep)
        j1 = indx[ff][tet[1]]; j2 = indx[ff][tet[2]]; j3 = indx[ff][tet[3]]
        j4 = indx[ff][tet[4]]; j5 = indx[ff][tet[5]]; j6 = indx[ff][tet[6]]
        dim = visqrt(j4)*visqrt(j5)/visqrt(j3) 
        fac = numchop(dim*Fsymb(j1,j2,j3,j4,j5,j6)) 
        ans = numchop(amps[ff]/fac )
        push!(sol,ans)
        push!(indxn,indxkeepu[i])
        #global cou = 0
        #if sign(fac) == sign(amps[i] )
        #    ans = numchop(amps[i] / fac)
        #    push!(sol,ans)
        #else
        #    ans = numchop(amps[i] / fac)
        #    push!(sol,-ans)
        #    #cou += 1
        #end
    end
    return indxn,sol#indx,sol#,couff1
end


end # module

