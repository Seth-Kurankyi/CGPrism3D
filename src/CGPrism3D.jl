module CGPrism3D
using LinearAlgebra
#using StatsBase

numchop(num::Real)=(abs(num) >= 1000*eps() ? num : zero(Real))

function numchop(num::Complex)
    numchop(imag(num)) == 0 ? numchop(real(num)) : Complex(numchop(real(num)), numchop(imag(num)))
end

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

function visqrtU(Utensor,a::Float64)
    indx = Utensor[1]
    amps = Utensor[2]
    fc = findfirst(x-> x == [0.,0.,0.,a,a,a], indx)
    ans = 0
    if typeof(fc) == Int64
        ans = ((-1+0im)^a )*sqrt(abs(amps[fc]))
    end
    ans
end

function TetraJK(i::Float64,j::Float64,m::Float64,k::Float64,l::Float64,n::Float64,J::Int64,K::Int64,alpha::Float64)
    #@assert J <= K
    alpha*Tetra6j(i,j,m,k,l,n,J) +(1-alpha)*Tetra6j(i,j,m,k,l,n,K)
end

function prismA(j1::Float64,j2::Float64,j3::Float64,j4::Float64,j5::Float64,j6::Float64,j7::Float64,
                j8::Float64,j9::Float64,j10::Float64,j11::Float64,j12::Float64,K::Int64)
    sol = 0 
    if delta(j1,j2,j3,K) != 0 && delta(j1,j5,j6,K) != 0 && delta(j2,j4,j6,K) != 0 && delta(j3,j4,j5,K) != 0 && delta(j1,j8,j9,K) != 0 && delta(j3,j7,j8,K) != 0 && delta(j2,j7,j9,K) != 0 && delta(j9,j10,j11,K) != 0 && delta(j7,j10,j12,K) != 0 && delta(j2,j11,j12,K) != 0  
        dims = prod(visqrt.([j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12],K))
        sol =  dims*Gsymb(j1,j2,j3,j4,j5,j6,K) * Gsymb(j1,j2,j3,j7,j8,j9,K) * Gsymb(j2,j7,j9,j10,j11,j12,K)
    end
    return sol
end

export dataTet, dataFsymb, dataPrsmA, dataGsymb

function dataPrsmA(J::Int64,K::Int64,alpha::Float64)
    @assert J <= K
    @assert 0. <= alpha <= 1.
    x = K+1
    y = K/2
    ampsInfo = Array{Float64,1}[]
    amps = Float64[]
    for j1 in 0.:0.5:y, j2 in 0.:0.5:y, j3 in 0.:0.5:y
        if delta(j1,j2,j3,K) != 0 
            for j4 in 0.:0.5:y, j5 in 0.:0.5:y
                if delta(j3,j4,j5,K) != 0 
                    for j6 in 0.:0.5:y
                        if delta(j1,j5,j6,K) != 0 && delta(j2,j4,j6,K) != 0
                            for j7 in 0.:0.5:y, j8 in 0.:0.5:y
                                if delta(j3,j7,j8,K) != 0
                                    for j9 in 0.:0.5:y
                                        if delta(j1,j8,j9,K) != 0 && delta(j2,j7,j9,K) != 0
                                            for j10 in 0.:0.5:y, j11 in 0.:0.5:y
                                                if delta(j9,j10,j11,K) != 0
                                                    for j12 in 0.:0.5:y
                                                        if delta(j7,j10,j12,K) != 0 && delta(j2,j11,j12,K) != 0 #
                                                            sol = numchop(alpha*prismA(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,J)+(1-alpha)*prismA(j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,K))
                                                            if sol != 0
                                                                indx = [j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12]
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
                            #if sol != 0
                            push!(t6s,sol)
                            push!(t6j,[ja,jb,jc,jd,je,jf])
                            #end
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

function dataGsymb(J::Int64,K::Int64,alpha::Float64)
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
                            sol = numchop(alpha*Gsymb(ja,jb,jc,jd,je,jf,J) +(1-alpha)*Gsymb(ja,jb,jc,jd,je,jf,K))
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

function tensorBlockOld(tensorM,posnA,posnB,posnC,spinC)
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
        println(matM)
        return matM,indxU,indxV
    else
        return 0,0,0
    end
end

export tensorBlock

function tensorBlock(tensorM,posnA,posnB,posnC,spinC)
    indx = tensorM[1]
    ampvals = tensorM[2]
    fc = findall(x-> x[posnC] == spinC, indx) #find spinsC in all spins length(posnC) = length(spinC)
    indxI = indx[fc] # this contains all spin configuration containing ampsC
    ampsI = ampvals[fc] # this contains the corresponding spin amplitude values for ampsC
    Acol = sort(unique(getindex.(indxI, [posnA])))
    Brow = sort(unique(getindex.(indxI, [posnB])))
    indxU = vcat.(Acol,repeat([spinC],length(Acol)))
    indxV = vcat.(Brow,repeat([spinC],length(Brow)))
    #indxU = []
    #indxV = []
    #println(ampInfo)
    la = length(Acol)
    lb = length(Brow)
    lpa = length(posnA)
    lpb = length(posnB)
    matM = zeros(la,lb)
    if la != 0 || lb !=0
        #if length(ampsI) == la*lb
        indxM = []
        #pos = collect(Iterators.flatten(([posnA,posnB,posnC])))
        pos = vcat(posnA,posnB,posnC)
        for i in indxI
            push!(indxM,i[pos])
        end
        sp = sortperm(indxM)
        ampInfo = indxM[sp]
        amps = ampsI[sp]
        for i in 1:la, j in 1:lb
            if ampInfo[1][1:lpa] == Acol[i] && ampInfo[1][lpa+1:lpb+lpa] == Brow[j]
                #println("sucess")
                matM[i,j] = amps[1]
                deleteat!(amps,1)
                deleteat!(ampInfo,1)
            end
        end
        return matM,indxU,indxV#,ampInfo
    else
        return 0,0,0
    end
    #return ampInfo, matM
end

export tensorGlueTet3DO, tensorGlueTet3DN, tensorGlue

function tensorGlueTet3DO(tensorB,tensorA,posnB,posnA,Utensor)
    indx, amps = tensorGlue(tensorB,tensorA,posnB,posnA)
    face = getindex.(indx, [posnB])
    ans = complex(ones(length(face)))
    #visqt = visqJK(a,J,K,alpha)
    for i in 1:length(face)
        #ans[i] = prod(visqrtT.(face[i],J,K,alpha)) 
        pp = 1
        for j in face[i]
            pp *= visqrtU(Utensor,j)
        end
        ans[i] = pp
    end
    ampsN = @. amps / ans
    ampsN = real(ampsN)+imag(ampsN)
    return indx, ampsN
end

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
    shrdAB = Tuple{Array{Float64,1},Int64,Int64}[]
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
    #println(shrdAB)
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

export tensorGlueTet3DN2

function tensorGlueTet3DN2(tensorB,tensorA,posnB,posnA,Utensor)# glue two tensors TA,TB along posnA from TA and posnB from TB 
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
    qq = []#Array{Array{Float64,1},1}[]
    amps = [] #Array{Float64,1}[]
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
        push!(shrdAB, [i, cna,cnb])
        #println("Ca,Cna = ",ca,"\t",cna)
        #println("Cb,Cnb = ",cb,"\t",cnb)
    end
    for i in shrdAB
        if i[2] == 0 || i[3] == 0 # Put condition for when one of a or b is zero
            nothing
        else
            for j in 1:i[2]
                for k in 1:i[3]
                    indx1 = vcat(sindxB[k],sindxA[j][lena])
                    pp=1
                    for k in i[1]
                        pp *= visqrtU(Utensor,k)
                    end
                    amps1 = numchop(sampsA[j]*sampsB[k]/pp)
                    push!(qq,indx1)
                    push!(amps,amps1)
                    #ans = @. sampsA[j]*sampsB[1:i[3]]
                end
                #pp = 1
                #
                #ans /= pp
                #ampsN = numchop.(real(ans)+imag(ans))
                #push!(amps,ampsN)
                #indxa = repeat([sindxA[j][lena]],i[3]) 
                #for j in indxA[i][lena]
                #indxa1 = vcat.(sindxB[1:i[3]],indxa)
                #push!(qq,indxa1)
            end
        end
        deleteat!(sindxA,1:i[2])
        deleteat!(sindxB,1:i[3])
        deleteat!(sampsA,1:i[2])
        deleteat!(sampsB,1:i[3])
        #cntb = cnt
    end
    #indxsol = collect(Iterators.flatten(qq))
    #sol = collect(Iterators.flatten(amps))
    return qq, amps #indxsol,sol
end

function tensorGlueTet3DN(tensorB,tensorA,posnB,posnA,Utensor)# glue two tensors TA,TB along posnA from TA and posnB from TB 
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
        push!(shrdAB, [i, cna,cnb])
        #println("Ca,Cna = ",ca,"\t",cna)
        #println("Cb,Cnb = ",cb,"\t",cnb)
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
                ampsN = numchop.(real(ans)+imag(ans))
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

# can also improve tensorGlueTet3D maybe?
export fullSplitTet3DN

function fullSplitTet3DN(dataM,posnA,posnB,posnC,Utensor,J::Int64,K::Int64,alpha::Float64) # posnA,B,C must be vectors @assert length(dataM) = vcat(A,B,C)
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
        #if mat != 0
            U, s, V = svd(mat)
            #truncate and use only first singular value
            trunc = 1#length(s)
            U1 = U[:,1:trunc] 
            V1 = V[:,1:trunc]
            s1 = numchop.(s[1:trunc])
            pp = 1
            #println(ampsC[i])
            for j in ampsC[i]
                #println(j)
                pp *= visqrtU(Utensor,j)
            end
            #@show trunc
            valU = U1*sqrt.(s1)*sqrt(pp) # multiply by leftover dimension factors
            valV = V1*sqrt.(s1)*sqrt(pp)
            # valU,valV are all either real or purely imaginary
            valU = real(valU) - imag(valU) # take -ve of imaginary part if it gives only imag
            valV = real(valV) + imag(valV)
            ja1 = blkU[1][1]; jb1 = blkU[1][2]; jg = blkU[1][3] # pick the first spins in U
            #Fixing sign problem -- keep the same sign for U and Tetra6j symbol 
            #fc = findfirst(x->x==[jd1,jd2,je1,ja1,jb1,jg],Utensor[1])
            #println(sign(Utensor[2][fc]),sign(valU[1]))
            if sign(valU[1]) == sign(TetraJK(jd1,jd2,je1,ja1,jb1,jg,J,K,alpha))
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
    #println("done")
    solU = collect(Iterators.flatten(ampsU))
    solV = collect(Iterators.flatten(ampsV))
    ansU = solU/solU[1] # normalization condition
    ansV = solV/solV[1]
    indxUs = collect(Iterators.flatten(indxsU))
    indxVs = collect(Iterators.flatten(indxsV))
    return (indxUs,ansU),(indxVs,ansV) #,indxUV
end


export tensorSumTetO, tensorSumTetN, tensorSum

function tensorSumTetN(ampJ,posnN,Utensor)
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

# Summing a tensor  over indices labelled posnN
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

#maybe can also improve tensorSumTet

function tensorSumTetO(tensor,posnN,Utensor)
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

export tensor22moveT, tensor22moveF

function swap2(vec,a,b)
    vec[a],vec[b] = vec[b], vec[a]
    return vec
end

#Same tensor for gluing and adding dimesion factors
function tensor22moveT(tensor,posF,Utensor) # J,K,alpha determines what Fsymbol to use
    glu = tensorGlueTet3DN(tensor,Utensor,posF,[1,5,2,4,6],Utensor)
    n = posF[end]
    m = length(glu[1][1])
    swp = @. swap2(glu[1],n,m)
    tensorN = swp, glu[2]
    ans = tensorSumTetO(tensorN,[m],Utensor)
    return ans
end

export tensor22moveT2, tensor22moveT3


function tensor22moveF(tensor,posF,J::Int64,K::Int64,alpha::Float64)
    glu = tensorGlue(tensor,dataFsymb(J,K,alpha),posF,[1,5,2,4,6])
    n = posF[end]
    m = length(glu[1][1])
    swp = @. swap2(glu[1],n,m)
    tensorN = swp, glu[2]
    ans = tensorSum(tensorN,[m])
    return ans
end


function tensorSum2(ampJ,posnN,tri,K)
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
        if delta(indxN[i-1][tri[1]],indxN[i-1][tri[2]],indxN[i-1][tri[3]],K) != 0.
            if indxN[i] == indxN[i-1]
                #println([indxN[i][tri[1]],indxN[i][tri[2]],indxN[i][tri[3]],K])
                ampsN[i] += ampsN[i-1]
            elseif numchop(ampsN[i-1]) != 0.
                push!(indq,indxN[i-1])
                push!(qq,ampsN[i-1])
            end
        end
    end
    push!(indq,indxN[end])
    push!(qq,ampsN[end])
    return indq, qq
end

function tensorSumTetO2(tensor,posnN,Utensor,tri,K)
    indx = tensor[1]
    amps = tensor[2]
    if length(posnN) == 1
        blkdims = getindex.(indx,posnN)
        sqdm = []
        for j in blkdims
            push!(sqdm,visqrtU(Utensor,j))
        end
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
    end
    ampsN = amps .*sqdm
    ampsN2 = real(ampsN) - imag(ampsN)
    data = indx, ampsN2
    ans = tensorSum2(data,posnN,tri,K)
    return ans
end

function tensor22moveT2(tensor,posF,Utensor,tri,K::Int64) # J,K,alpha determines what Fsymbol to use
    glu = tensorGlueTet3DN(tensor,Utensor,posF,[1,5,2,4,6],Utensor)
    n = posF[end]
    m = length(glu[1][1])
    swp = @. swap2(glu[1],n,m)
    tensorN = swp, glu[2]
    #if length(tri) == 3
    ans = tensorSumTetO2(tensorN,[m],Utensor,tri,K)
    #elseif length(tri) == 0
    #    ans = tensorSumTetO(tensorN,[m],Utensor)
    #end
    return ans
end

function tensor22moveT3(tensor,posF,Utensor,tri,K::Int64) # J,K,alpha determines what Fsymbol to use
    glu = tensorGlueTet3DN(tensor,Utensor,posF,[1,5,2,4,6],Utensor)
    n = posF[end]
    m = length(glu[1][1])
    swp = @. swap2(glu[1],n,m)
    tensorN = swp, glu[2]
    indx, amps = tensorSumTetO(tensorN,[m],Utensor)
    ls = []
    la = []
    #if length(tri) == 3
    for i in 1:length(indx)
        if delta(indx[i][tri[1]],indx[i][tri[2]],indx[i][tri[3]],K) != 0
            push!(ls,indx[i])
            push!(la,amps[i])
        end
    end
    return ls,la
    #elseif length(tri) == 0
    #    return indx, amps
    #end
end

export permuteInd, permuteIndN, tensorPermute, tensorPermuteN

function permuteInd(vec)
    perm = [5,1,11,2,6,3,10,12,7,8,4,9]
    nvec = vec[perm]
    return nvec
end


function tensorPermute(tensor)
    indx = tensor[1]
    indxN = permuteInd.(indx)
    return indxN , tensor[2]
end

function permuteIndN(vec,perm)
    #perm = [5,1,11,2,6,3,10,12,7,8,4,9]
    nvec = vec[perm]
    return nvec
end


function tensorPermuteN(tensor,perm)
    indx = tensor[1]
    indxN = []
    @assert length(indx[1]) == length(perm)
    for i in indx
        push!(indxN,i[perm])
    end
    return indxN , tensor[2]
end

export tensorProjectPrism

function tensorProjectPrism(tensor,K)
    indx = tensor[1]
    amps = tensor[2]
    indxP = []
    ampsP = []
    for i in 1:length(indx)
        if delta(indx[i][1],indx[i][2],indx[i][3],K) != 0 
            if delta(indx[i][3],indx[i][4],indx[i][5],K) != 0 
                if delta(indx[i][1],indx[i][5],indx[i][6],K) != 0 && delta(indx[i][2],indx[i][4],indx[i][6],K) != 0
                    if delta(indx[i][3],indx[i][7],indx[i][8],K) != 0
                        if delta(indx[i][1],indx[i][8],indx[i][9],K) != 0 && delta(indx[i][2],indx[i][7],indx[i][9],K) != 0
                            if delta(indx[i][9],indx[i][10],indx[i][11],K) != 0
                                if delta(indx[i][7],indx[i][10],indx[i][12],K) != 0 && delta(indx[i][2],indx[i][11],indx[i][12],K) != 0
                                    push!(indxP,indx[i])
                                    push!(ampsP,amps[i])
                                end
                            end
                        end
                    end
                end
            end
        end                         
    end
    return indxP , ampsP
end

end # module
