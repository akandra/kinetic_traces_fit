function get_data(data_path, filenames)::Vector{Matrix{Float64}}

    data = Vector{Matrix{Float64}}
    data= []
    for f in filenames
        push!(data, readdlm(joinpath(data_path,f), Float64, comments=true))
    end
    return data
end

function get_data_cov(data_path, filenames)::Tuple{String,Vector{Matrix{Float64}}}

    cov_species = fill("",size(filenames,1))

    for (i,f) in enumerate(filenames)
        for s in eachline(joinpath(data_path,f))
            if !isnothing(match(r"species:",s)) 
                cov_species[i] = split(s)[end]
            end
        end
        if cov_species[i] == ""
            println(f)
        end
    end

    if any(x-> x == "", cov_species)
        error("get_data_cov: Put species name into above coverage files, please")
    end

    if any(x-> x!=cov_species[1],cov_species)
        error("get_data_cov: multiple coverage species names")
    end
    
    return cov_species[1], get_data(data_path, filenames)
end