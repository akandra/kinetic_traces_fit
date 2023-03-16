module m1
    using DelimitedFiles

    function get_data(filename)
   
        dropbox_path = "../../"
        data_path = dropbox_path * "Kinetics of Surface Reactions/data/"
        data_fn = data_path * filename
        println("hey: data_fn = ", data_fn)
        readdlm(data_fn, Float64)
    end
end

println("invoking m1.get_data")
data = m1.get_data("20170825-beam.dat");
println("first three rows")
data[1:3,:]
# println("typeof(data)", typeof(data))
# data
x = 1:10; y = rand(10); # These are the plotting data
plot(x, y)


