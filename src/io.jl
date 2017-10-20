function readData(filename::AbstractString)

    !isfile(filename) && error("Input file not found")

	phylo_profile=open(filename)
	list_species=String[]
	list_domains=String[]
	for line in eachline(phylo_profile)
		data=(split(chomp(line)))
		N=length(data)
		for k=1:N
			if(k==1)
				if(in(data[k],list_species)==false)
					push!(list_species,data[k])
				end
				
			else
				if(in(data[k],list_domains)==false)
					push!(list_domains,data[k])
				end
			end
		end
	end
	close(phylo_profile)
	sort!(list_species)
	sort!(list_domains)
	return list_species,list_domains
end

function makePhyloMatrix(filename::AbstractString,
			 list_species_ord,
			 list_domains_ord)

	#N=number of domains, M=number of species
	M=length(list_species_ord)
	N=length(list_domains_ord)

	#make phylogenetic profile
	#P=phylogenetic matrix, M times N matrix
	P=zeros(Int8,M,N)

    
    #map species/domains to their index in the list
    dict_species=Dict(list_species_ord[i] => i for i=1:length(list_species_ord))
    dict_domains=Dict(list_domains_ord[i] => i for i=1:length(list_domains_ord))

	phylo_profile=open(filename)
	for line in eachline(phylo_profile)
		data=(split(chomp(line))) 
		len=length(data)
		index_spec=0
		index_dom=0
        index_spec=dict_species[data[1]]
		for k=2:len
			index_dom=dict_domains[data[k]]
			P[index_spec,index_dom]=1 #binary matrix
		end
	end
	return P
end

#############################################################
#TO be implemented well, that's just for now
function read_benchmark(filename::AbstractString)
	
	!isfile(filename) && error("Benchamark file not found")

	bench_matr=readdlm(filename)
	return bench_matr
end
##################################################

function print_result(final_matrix::Matrix)
	result_file=open("results.ccs", "w")

	lungh,num=size(final_matrix)
	if(num>3) 
		for i=1:lungh
			@printf(result_file,"%s %s %f", final_matrix[i,1], final_matrix[i,2], final_matrix[i,3])
			for k=3:num
				@printf(result_file, " %d", final_matrix[i,k])
			end
			@printf(result_file, "\n")
		end
	else
		for i=1:lungh
			@printf(result_file,"%s %s %f \n", final_matrix[i,1], final_matrix[i,2], final_matrix[i,3])
		end
	end
	close(result_file)
end


