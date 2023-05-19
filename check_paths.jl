function check_paths(path_list)
    
    # Check if data and model directories exist
path_error_str = ""
if !isdir(data_path)
    path_error_str = path_error_str * "data path"
end
if !isdir(model_path)
    path_error_str = path_error_str == "" ? path_error_str * "model path" : path_error_str * " and model path"
end
if !isdir(data_path) || !isdir(model_path)
    error("Please define "*path_error_str)
end

end