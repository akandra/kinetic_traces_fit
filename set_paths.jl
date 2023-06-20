# path variables
path       = "./"
data_path  = "./"
model_path = "./"
fit_path   = "./"
results_path = "./"

fn_delim = "-"
fn_pump_sfx = "beam"
fn_cov_sfx ="Oini"


working_dir(str::String)   = (global path = str)

path_to_data(str::String)  = (global data_path = joinpath(path, str))
path_to_data(str::String, folder::String)  = (global data_path = joinpath(folder, str))

path_to_model(str::String) = (global model_path = joinpath(path, str))
path_to_model(str::String, folder::String) = (global model_path = joinpath(folder, str))

path_to_fit_function(str::String) = (global fit_path = joinpath(path, str))
path_to_fit_function(str::String, folder::String) = (global fit_path = joinpath(folder, str))

path_to_output(str::String, m) = (global results_path = splitext(m)[1]*"_"*str )

model(str::String) = include(joinpath(model_path,str)*".jl")
fit_function(str::String) = include(joinpath(fit_path,str)*".jl")

delim(str) = (global fn_delim = str)
pump_sfx(str) = (global fn_pump_sfx = str)
cov_sfx(str) = (global fn_cov_sfx = str)