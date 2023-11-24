from numpy import loadtxt

# algorithm to find root mean square fluctiation of peptide using xvg file of rmsf output by simulation
lines = loadtxt("rmsf_new.xvg", comments="#", unpack=False)
a,b = zip(*lines)
c = []
for item in lines:
	if item[1] > 1.75:
		c.append(item) # flag larger rmsf values

Hyd = peptide.select_atoms("type H") # select hydrogen atoms
Hyd = list(Hyd)

# convert lists to strings
listToStr = '\n '.join([str(elem) for elem in Hyd])
listToStr2 = ' '.join([str(elem) for elem in listToStr])

# s was a pre defined list within the terminal output from HPC Rosalind
# s was list of clean (non noisy) rmsf values
k = []
for elem in s:
	k.append(elem[6:9]) # interested in middle few vectors of each vector in s
k = [s.strip(': ') for s in k] # remove white spaces
print(k)

# finding list of ints using large flagged rmsf values from c
list_of_integers = [int(i) for i in k] 
x, y = zip(*c)

# find overlapping rmsf's from our list of integers and list of y values i.e numbers from c
w = set(x) & set(list_of_integers)




