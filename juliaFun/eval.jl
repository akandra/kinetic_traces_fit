
# Creates a symbol corresponding to the variable text
s=Symbol("text")

# Sets the variable text to some value

# with function eval()
eval(:( $s = 1 ))
text

# with macro eval
@eval($s = 2)
text
