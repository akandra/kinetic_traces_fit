function get_data(data_path, filenames)::Vector{Matrix{Float64}}

    data = Vector{Matrix{Float64}}
    data= []
    for f in filenames
        push!(data, readdlm(data_path*f, Float64, comments=true))
    end
    return data
end