# path variables
path       = "./"
working_dir(str::String)   = (global path = str)

data_path  = "./"
path_to_data(str::String)  = (global data_path = joinpath(path, str))
path_to_data(str::String, folder::String)  = (global data_path = joinpath(folder, str))

model_path = "./"
path_to_model(str::String) = (global model_path = joinpath(path, str))
path_to_model(str::String, folder::String) = (global model_path = joinpath(folder, str))
rate_constants = String[]

function model(str::String)
    # load user's model file 
    include(joinpath(model_path,str)*".jl")
    # get rate constant names
    global rate_constants = rate_constants_base .* rate_constants_sfx
end

fit_path   = "./"
path_to_fit_function(str::String) = (global fit_path = joinpath(path, str))
path_to_fit_function(str::String, folder::String) = (global fit_path = joinpath(folder, str))
fit_function(str::String) = include(joinpath(fit_path,str)*".jl")

output_path = "./"
path_to_output(str::String, m) = (global output_path = splitext(m)[1]*"_"*str )

local_fit_path    = "./"
path_to_local_fit_results(str::String)   = (global local_fit_path = str)

fn_delim = "-"
delim(str) = (global fn_delim = str)

fn_pump_sfx = "beam"
pump_sfx(str) = (global fn_pump_sfx = str)

fn_cov_sfx ="Oini"
cov_sfx(str) = (global fn_cov_sfx = str)

# conditions for selecting data
conds = []
conditions(args...) = push!(conds, args)

# what to do variable 
wtd = []
what_to_do(args...) = (global wtd = [args...])

# maximum temperature for Arrhenius fits
T_cutoff::Number = Inf
Arrh_Tmax(x::Number) = (global T_cutoff = x)

# Set kinetic trace data cutoff fraction
data_cutoff_fraction  = "off"
cutoff_fraction(x) = (global data_cutoff_fraction = x)

# initial guesses

guess_Arrh_local = []
rate_constant_Arrh_local(; kargs...) = push!(guess_Arrh_local, kargs)

guess_Arrh_global = []
rate_constant_Arrh_global(; kargs...) = push!(guess_Arrh_global, kargs)

guess_hTST_local = []
rate_constant_hTST_local(; kargs...) = push!(guess_hTST_local, kargs)

guess_hTST_global = []
rate_constant_hTST_global(; kargs...) = push!(guess_hTST_global, kargs)

guesses = []
guess(; kargs...) = push!(guesses, kargs)

par_guesses = []
guess_par(; kargs...) = push!(par_guesses, kargs)

# get step densities

θs = []
step_density(; kargs...) = push!(θs, kargs)
