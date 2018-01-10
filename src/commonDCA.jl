#P is the phylogenetic matrix which is M x N matrix, M=num species(configuration), N=num of domains (spins)

##################################################
# Evaluate hamming distance
##################################################

@inline eval_hamm(ai,bi)= ai != bi ? 1 : 0

function hamming_distance(a::AbstractArray, b::AbstractArray)
    if(length(a) != length(b))
        throw( DimensionMismatch("First array has length $(length(a)), while second array has length $(length(b))"))
    end
    if length(a) == length(b)
        dist=0.0
        @simd for I in eachindex(a,b)
            @inbounds ai=a[I]
            @inbounds bi=b[I]
            dist+=eval_hamm(ai,bi)
        end
        dist/=length(a)
    end

    return dist

end


            

##################################################
# compute weights
##################################################

function compute_weigths(P::Matrix{Int8}, theta::Real)

    M,N=size(P)
    W=ones(M)

    theta= Float64(theta)
    @printf("%s%0.2f ",  "theta=",theta)

    if theta==0
        @printf("%s%d%s%d%s%d\n","M=",M," N=",N," Meff=",M)
        return W, Float64(M)
    end


    ZZ=Vector{Int8}[vec(P[i,:]) for i=1:M]
    @simd for i=1:M
        Zi=ZZ[i]
        @simd for j=(i+1):M
            Zj=ZZ[j]
            #if sequence has distance < threshold -> increase weight
            if hamming_distance(Zi,Zj)<theta 
                @inbounds W[i] += 1
                @inbounds W[j] += 1
            end
        end
        #weight
    end
    for i = 1:M
        W[i] = 1 / W[i]
    end 

    Meff=sum(W)

    @printf("%s%d%s%d%s%0.2f\n","M=",M," N=",N," Meff=",Meff)

    return W, Meff
end

##################################################
# compute frequencies
##################################################

# W is the vector of weights
function compute_freqs(P::Matrix{Int8}, W::Vector{Float64}, Meff::Float64)
    M,N=size(P)
    ZZ=Vector{Int8}[vec(P[:,i]) for i=1:N]

    Pij=zeros(N,N)
    Pi=zeros(N)
    for i=1:N
        Zi=ZZ[i]
        for j=i:N
            Zj=ZZ[j]
            for k=1:M
                Pij[i,j]+=W[k]*Zi[k]*Zj[k]
            end
        end
    end

    
    for i=1:N 
        Pij[i,i]/=Meff
        for j=(i+1):N
            Pij[i,j]/=Meff
            Pij[j,i]=Pij[i,j]
        end
    end
        
     Pi=diag(Pij)
    return Pi,Pij
end

function compute_new_frequencies(P::Matrix{Int8},theta::Real)

    W,Meff=compute_weigths(P,theta)
    Pi_t,Pij_t=compute_freqs(P,W,Meff)

    return Pi_t,Pij_t,Meff,W
end

##################################################
# add pseudocount
##################################################
function add_pseudocount(Pi_true::Vector{Float64},Pij_true::Matrix{Float64},Meff::Float64,pc::Float64)
   
    N=length(Pi_true)
    Pij=(1-pc) * Pij_true .+ pc/4
    for i = 1:N
        Pij[i,i] += pc/2
    end

    Pi = ( 1 - pc ) * Pi_true .+ pc/2
    return Pi, Pij
end

##################################################
# compute empirical correlation matrix C
##################################################
compute_C(Pi::Vector{Float64}, Pij::Matrix{Float64}) = Pij - Pi * Pi'

##################################################
# APC correction (minimize backgroung influence like phylogeny and site entropy) 
##################################################
function correct_APC(S::Matrix)
    N = size(S,1)
    Si = sum(S,1)
    Sj = sum(S,2)
    Sa = sum(S) * (1 - 1/N)

    S -= (Sj * Si) / Sa
    return S
end



