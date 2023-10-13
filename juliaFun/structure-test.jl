using Parameters

# structure to define species
@with_kw struct adsorbate

    name::String = "" # "O"
    site::String = "" # "s", "t", "c", "g"
    σ::Int       = 0  # occupancies
    type::String = "" # "beam", "coverage", "intermediate", "product"

end

species = [
    adsorbate(name="H₂", site="g", type = "beam"),
    adsorbate(name="O", site="s", σ = 2, type = "coverage"),
    adsorbate(name="O", site="t", σ = 4, type = "coverage")
    ]

findall(x-> x.name == "O" ,species)

# Creating and assigning variables
[@eval $(Symbol("i",species[i].name,species[i].site)) = ($i) for i in 1:length(species)]
iH₂g
iOs
iOt