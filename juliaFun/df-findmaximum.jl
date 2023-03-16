df = DataFrame(Data1 = rand(10), Data2 = rand(10));

max_value = maximum(df.Data1)
findfirst(==(max_value), df.Data1)