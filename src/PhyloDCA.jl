module PhyloDCA


using GaussDCA
using PlmDCA
using PyPlot

using HypothesisTests 

export 
	#phylogenetic distance class
	Hamming,
	Correlation,
	pValue,
	mfDCA,
	plmDCA,
	PhyloOut,

	#function
    phylodca,
	make_phylo_profile,
	evaluate_distance,
	compute_distance_matrix

include("io.jl")
include("types.jl")
include("commonDCA.jl") #tools for DCA inference, both for mfDCA and plmDCA


########################################
# MAIN FUNCTIONS
########################################
#method 1, no file containing knwon functional relations is given
function phylodca(filename_phylo::AbstractString,
		  dist::UnionDistances)

	@printf("%s", "Computing phylogenetic matrix...")
	P,domains,species=make_phylo_profile(filename_phylo)
	@printf("%s\n", "done")

	@printf("%s", "Computing distance matrix...")
	dist_matr=compute_distance_matrix(P,dist)
	@printf("%s\n", "done")

	@printf("%s\n", "Sorting domain-pairs by their phylogenetic distance")
	final=compute_final_matrix(dist_matr,domains, dist)

	print_result(final)

	x=PhyloOut(P,domains,species,dist_matr,final)
	return x 
end


#method 2: one (or more) file(s) containing kwnon functional relation is given
function phylodca(filename_phylo::AbstractString,
		  dist::UnionDistances,
		  known_relations::Vararg{AbstractString}) #VarArg allows variable number of inputs

	@printf("%s", "Computing phylogenetic matrix...")
	P,domains,species=make_phylo_profile(filename_phylo)
	@printf("%s\n", "done")

	@printf("%s\n", "Computing distance matrix...")
	dist_matr=compute_distance_matrix(P,dist)

	@printf("%s\n", "Sorting domain-pairs by their phylogenetic couplings")
	final=compute_final_matrix(dist_matr,domains, dist,known_relations)
	print_result(final)

	x=PhyloOut(P,domains,species,dist_matr,final)
	return x
end







##############################
#1) construct phylogenetic matrix
##############################
function make_phylo_profile(filename::AbstractString)

	list_species,list_domains=readData(filename)
	P=makePhyloMatrix(filename, list_species,list_domains)

	#PLOT phylogenetic matrix.... not a good idea...

	#plot_matrix(P,list_domains,list_species)
	#plot_hist(P)
	#implementa plot_distr_domains (per vedere se mettere cut_off)

	return P,list_species,list_domains

end

##############################
#2) construct distance matrix
##############################
function compute_distance_matrix(P::Matrix{Int8},
				 dist::UnionDistances)

	!isa(P,Matrix) && error("Phlogenetic Matrix not given")

	dist_matrix=evaluate_distance(dist,P)
end

############################################################
#3) construct result matrix sorted by phylogenetic distances 
############################################################
#method 1, no file containing knwown functional relations is given
function compute_final_matrix(distance_matrix::Matrix{Float64},
			      list_domains::Array{String},
			      dist::UnionDistances)

	res_matr=make_result_matrix(distance_matrix,list_domains)
	sort_res_matr=sort_matrix(res_matr,dist)
	return sort_res_matr
end

#method 2: one (or more) file(s) containing kwnown functional relations is given
function compute_final_matrix(distance_matrix::Matrix{Float64},
			      list_domains::Array{String},
			      dist::UnionDistances,
			      bench_filenames::Vararg{AbstractString}) #VarArg allows variable number of inputs

	res_matr=make_result_matrix(distance_matrix,list_domains)

	for i=1:length(bench_filenames)
		println(bench_filenames[i])
		res_matr=add_benchmark(res_matr,list_domains,bench_filenames[i])
	end

	sort_res_matr=sort_matrix(res_matr,dist)
	return sort_res_matr
end


######################################################################
# plot the PPV (positive predictive value)
######################################################################
function  plot_PPV(final_matrix::Matrix{Any},
		   ppv_index::Vararg{Int})
	ss=sum(final_matrix[:,k] for k in ppv_index)
	num_pred,ppv_rate=compute_tp_rate(Array{Int}(ss))
	xscale("log")
	xlabel("Number of predictions")
	ylabel("PPV")
	plot(num_pred,ppv_rate)
end
















####################################################################################################
#plotting the matrix 
####################################################################################################
function plot_matrix(P::Matrix{Int8},
		     x_label,
		     y_label)
    fig = figure("Phylogenetic Profile", figsize=(10,10))
    ax = gca() 
    x_label=union(0,x_label)
    y_label=union(0,y_label)
    ax[:set_xticklabels](x_label)
    ax[:set_yticklabels](y_label)
    ax[:matshow](P) #cmap=plt.cm.gray 
end

function plot_hist(P::Matrix{Int8})
    fig = figure("Perc of batt", figsize=(10,10))
    ax = gca() 
    s=sum(P,1)
    ax[:hist](s) #cmap=plt.cm.gray 
    
end

####################################################################################################

###############
#compare with set of known relations 
##############

#for j>i, the pair {i,j} in the list {{1,2}, {1,3}, {1,4},...,{1,N},{2,3},... } is located in offset(i,j,N)
function offset(i::Int,j::Int,N::Int)
	off=(i-1)*N - i*(i-1)/2 +(j-i)
	return Int(off)
end

function make_result_matrix(Jij::Matrix{Float64},
			    list_domains::Array{String})
	N=length(list_domains)
	num=Int(N*(N-1)/2)
	result_matrix=Matrix{Any}(num,3)

	cont=1
	for i=1:N
		for j=(i+1):N
			result_matrix[cont,1]=list_domains[i]
			result_matrix[cont,2]=list_domains[j]
			result_matrix[cont,3]=Jij[i,j]
			cont+=1
		end
	end
	return result_matrix
end

function add_benchmark(result_matrix::Matrix,
		       list_domains::Array{String},
		       filename::AbstractString)
	
		bench=read_benchmark(filename)
		tot=size(result_matrix)[1]
		vec_bench=zeros(Int,tot)
		Nb=size(bench)[1]
		num_domains=size(list_domains)[1]
		for i=1:Nb
			dom1=bench[i,1]
			dom2=bench[i,2]
			pos1=find(dom1.==list_domains)[1]
			pos2=find(dom2.==list_domains)[1]
            if(length(pos1)>0 && length(pos2)>0)
                vec_bench[offset(pos1,pos2,num_domains)]=1
            end
		end

		result_matrix=hcat(result_matrix,vec_bench)

	return result_matrix
end

#N.b sort matrix must depend on the phylogenetic distance:
#Hamming, pValue distance-> minimum 
#correlation, mfDCA,plmDCA-> maximum
#
function sort_matrix(final_matrix::Matrix,
		     ::Union{Hamming,pValue})
	sorted_matr=sortrows(final_matrix, by=x->x[3])
	return sorted_matr
end
function sort_matrix(final_matrix::Matrix,
		     ::Union{Correlation,mfDCA,plmDCA,MutualInfo})
	sorted_matr=sortrows(final_matrix, by=x->x[3],rev=true)
	return sorted_matr
end




#################################
##plot tp rate
###################

function compute_tp_rate(tp::Array{Int})
    m=length(tp)
    tp_check=zeros(m)
    for a=1:m
        tp_check[a]=(tp[a]!=0)
    end
    tp_cumsum=cumsum(tp_check)

    num_pred=zeros(m)
    ppv_rate=zeros(m)
    for a=1:m
	    num_pred[a]=a
	    ppv_rate[a]=tp_cumsum[a]/a
    end
    return num_pred,ppv_rate
end





end


