# path variables
path       = "./"
working_dir(str::String)   = (global path = str)

data_path  = "./"
path_to_data(str::String)  = (global data_path = joinpath(path, str))
path_to_data(str::String, folder::String)  = (global data_path = joinpath(folder, str))

model_path = "./"
path_to_model(str::String) = (global model_path = joinpath(path, str))
path_to_model(str::String, folder::String) = (global model_path = joinpath(folder, str))
model(str::String) = include(joinpath(model_path,str)*".jl")

fit_path   = "./"
path_to_fit_function(str::String) = (global fit_path = joinpath(path, str))
path_to_fit_function(str::String, folder::String) = (global fit_path = joinpath(folder, str))
fit_function(str::String) = include(joinpath(fit_path,str)*".jl")

output_path = "./"
path_to_output(str::String, m) = (global output_path = splitext(m)[1]*"_"*str )

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
wtd = ()
what_to_do(args...) = (global wtd = args)

# maximum temperature for Arrhenius fits
T_cutoff::Number = Inf
Arrh_Tmax(x::Number) = (global T_cutoff = x)

# Set kinetic trace data cutoff fraction
data_cutoff_fraction  = "off"
cutoff_fraction(x) = (global data_cutoff_fraction = x)

# initial guesses
rate_guesses = []
guess_rate(; kargs...) = push!(rate_guesses, kargs)

par_guesses = []
guess_par(; kargs...) = push!(par_guesses, kargs)

occ_factors = []
occupancies(; kargs...) = push!(occ_factors, kargs)

θs = []
step_density(; kargs...) = push!(θs, kargs)
