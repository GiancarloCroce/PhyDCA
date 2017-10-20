

############################
## TYPE DEFINITION ####
############################
abstract PhylogeneticDistance
type Hamming<:PhylogeneticDistance end
type Correlation<: PhylogeneticDistance end
type MutualInfo <: PhylogeneticDistance end
type pValue <: PhylogeneticDistance end
type mfDCA <: PhylogeneticDistance end
type plmDCA<: PhylogeneticDistance end
UnionDistances=Union{Hamming,Correlation,MutualInfo,pValue,mfDCA,plmDCA}

immutable PhyloOut 
	PhyloProfile::Array{Int8,2}
	list_domains::Array{String,1}
	list_species::Array{String,1}
	PhyloDistance::Array{Float64,2}
	result::Array{Any,2}
end



################################################
# METHODS DEFINITION 
###############################################

function evaluate_distance(::Hamming,
                          P::Matrix{Int8})
	M,N=size(P)
	hamm_matrix=Matrix{Float64}(N,N)
	for i1=1:N
		for i2=(i1+1):N
			hamm_matrix[i1,i2]=sum(P[:,i1].!=P[:,i2])/N
		end
	end
	hamm_matrix+=hamm_matrix'
	return hamm_matrix 
end

function evaluate_distance(::Correlation,
			   P::Matrix{Int8})

	return cor(P)
end

########################################
# check if MI is good...
function evaluate_distance(::MutualInfo,
			   P::Matrix{Int8})
	P=ones(P)*2-P
	Z=P'
	N, M = size(Z)
	q = Int(maximum(Z))
	theta=:auto
	pseudocount::Real=0.25 #for the 0,1 Ising model 
	#frequencies of one(s)
	Pi_true, Pij_true, Meff, _ = GaussDCA.compute_new_frequencies(Z, q, theta)
	Pi, Pij= GaussDCA.add_pseudocount(Pi_true, Pij_true, Float64(pseudocount), N, q)

	
	#compute mutual information
	mi=zeros(N,N)
	for i=1:N
		for j=(i+1):N
			Pij_00=1+Pij[i,j]-Pi[i]-Pi[j]
			Pij_01=Pi[j]-Pij[i,j]
			Pij_10=Pi[i]-Pij[i,j]
			Pij_11=Pij[i,j]
			Pi_0=1-Pi[i]
			Pj_0=1-Pi[j]
			Pi_1=Pi[i]
			Pj_1=Pi[j]
			mi[i,j]=Pij_00*log(Pij_00/Pi_0*Pj_0)+Pij_01*log(Pij_01/Pi_0*Pj_1)+Pij_10*log(Pij_10/Pi_1*Pj_0)+Pij_11*log(Pij_11/Pi_1*Pj_1)
		end
	end

	mi+=mi'
	return mi 
end
###############################

function evaluate_distance(::pValue,
                          P::Matrix{Int8})
	M,N=size(P)

	d_fisher=zeros(N,N)
	m11=0
	m12=0
	m21=0
	m22=0
	for i=1:N
		for j=(i+1):N
			for lin=1:M
				m11=m11+(Int(P[lin,i]==0 && P[lin,j]==0));
				m12=m12+(Int(P[lin,i]==1 && P[lin,j]==0));
				m21=m21+(Int(P[lin,i]==0 && P[lin,j]==1));
				m22=m22+(Int(P[lin,i]==1 && P[lin,j]==1));
			end
			res=FisherExactTest(m11,m12,m21,m22);
			pvalue=HypothesisTests.pvalue(res);
			d_fisher[i,j]=pvalue;
		end
	end
	d_fisher+=d_fisher'
	return d_fisher

end

function evaluate_distance(::mfDCA,
                           P::Matrix{Int8};
                           pseudocount::Real=0.25, #for the 0,1 Ising model 
                           theta::Float64=0.2)
    

	M, N = size(P)

	Pi_true, Pij_true, Meff, _ = compute_new_frequencies(P, theta)
	Pi, Pij = add_pseudocount(Pi_true, Pij_true, Meff, Float64(pseudocount))
	C = compute_C(Pi, Pij)

	mJ = -inv(cholfact(C))

	S = correct_APC(mJ)

    #remove diagonal elements (fields) 
	S-=Diagonal(S)
	return S
end

##MY PLM DCA
#function evaluate_distance(::plmDCA,
#                           P::Matrix{Int8};
#                           lambdaJ::Real=0.02,
#                           lambdaH::Real=0.05,
#                           theta::Float64=0.2)
#
#    M,N=size(P)
#    W,Meff=compute_weights(P,theta)
#    #PlmAlg is a type
#	plmalg =PlmDCA.PlmAlg(method,verbose, epsconv ,maxit, boolmask)
#	plmvar = PlmDCA.PlmVar(N,M,q,q*q,gaugecol,lambdaJ,lambdaH,Z,W)
#	Jmat, pslike = PlmDCA.MinimizePLAsym(plmalg,plmvar)                
#end
#

   

#TO IMPLEMENT MY PLMDCA
#
#function evaluate_distance(::plmDCA,
#                          P::Matrix{Int8})
#
#    lambdaJ::Real=0.02
#    lambdaH::Real=0.05
#    theta=:auto
#    epsconv::Real=1.0e-5
#    maxit::Int=1000
#    min_separation::Int = 1
#    verbose::Bool=true
#    method::Symbol=:LD_LBFGS #(questo metodo mi da un errore! ERROR: nlopt failure: -1... Capisci cos'e)
#    boolmask::Union{Array{Bool,2},Void}=nothing
#    gaugecol::Int=-1
#
#	P=ones(P)*2-P
#	Z=P'
#
#	N, M = size(Z)
#	q = Int(maximum(Z))
#	J=extract_J(Jmat, plmvar, min_separation)
#	
#	S = GaussDCA.correct_APC(J)
#	#remove diagonal elements (fields...)
#	S-=Diagonal(S)
#
#	return S 
#
#end
#
function extract_J(Jmat::Array{Float64,2}, var, min_separation::Int)

    q = var.q
    N = var.N

    JJ=reshape(Jmat[1:end-q,:], q,q,N-1,N)
    Jtemp1=zeros( q,q,Int(N*(N-1)/2))
    Jtemp2=zeros( q,q,Int(N*(N-1)/2))
    
    l = 1

    for i=1:(N-1)
        for j=(i+1):N
            Jtemp1[:,:,l]=JJ[:,:,j-1,i]; #J_ij as estimated from from g_i.
            Jtemp2[:,:,l]=JJ[:,:,i,j]'; #J_ij as estimated from from g_j.
            l=l+1;
        end
    end

    J1=zeros(q,q,Int(N*(N-1)/2))
    J2=zeros(q,q,Int(N*(N-1)/2))


    #switch to Ising gauge
    for l=1:Int(N*(N-1)/2)
        J1[:,:,l] = Jtemp1[:,:,l]-repmat(mean(Jtemp1[:,:,l],1),q,1)-repmat(mean(Jtemp1[:,:,l],2),1,q) .+ mean(Jtemp1[:,:,l])
        J2[:,:,l] = Jtemp2[:,:,l]-repmat(mean(Jtemp2[:,:,l],1),q,1)-repmat(mean(Jtemp2[:,:,l],2),1,q) .+ mean(Jtemp2[:,:,l])
    end
    J = 0.5 * ( J1 + J2 )


    #switch to lattice-gas model
    lattice_gas_J=zeros(N,N)
    l = 1
    for i = 1:N-1
        for j=i+1:N
		lattice_gas_J[i,j] = J[1,1,l]-J[1,2,l]-J[2,1,l]+J[2,2,l]
		l += 1
        end
    end
    

    lattice_gas_J+=lattice_gas_J'

    return lattice_gas_J
end


