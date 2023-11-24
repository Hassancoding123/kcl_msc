# algorithm to find the minimum distance between carbon alpha atoms in two different peptides
# aim = to see if peptides would bond by remaining under a threshold minimum distance

# in trimer system, minimum distance between peptide 1 and peptide 2
min_1_2 = []

# only selecting carbon alphas as otherwise computationally expensive
CA1 = peptide1.select_atoms("name CA")
CA2 = peptide2.select_atoms("name CA")
CA3 = peptide3.select_atoms("name CA")

# for each time step in the trajectory 
for ts in trajectory:
	# find carbon alpha positions in vector form for each monomer
	pos_1 = CA1.positions
	pos_2 = CA2.positions
	pos_3 = CA3.positions
	
	# find all possible position combinations between two carbon alpha atoms in peptide 1 and 2 
	all_possible_pos = product(pos_1, pos_2)
	
	# find the changes between positions of all possible position combinations
	norms = []
	changes = [abs(a-b) for a,b in all_possible_pos]
	changes_x = [abs(item[0]) for item in changes]
	changes_y = [abs(item[1]) for item in changes]
	changes_z = [abs(item[2]) for item in changes]

	#account for the periodic box i.e non erroneous changes are ignored => cycle through each x,y and z of each vector 
	non_pbc_x = []
	non_pbc_y = []
	non_pbc_z = []

	# if vector element larger than period box then this is diregarded, else added to non pbc element list
	for i in changes_x:
		if i > pbc_scale:
			a = dimensions - i
	        else:
			a = i
	        non_pbc_x.append(a)
	    
	for i in changes_y:
		if i > pbc_scale:
			a = dimensions - i
	        else:
			a = i
	        non_pbc_y.append(a)
		
	for i in changes_z:
		if i > pbc_scale:
			a = dimensions - i
	        else:
			a = i
	        non_pbc_z.append(a)

	# concatenate all non pbc elements into vectors again
	non_pbc_vecs = zip(non_pbc_x, non_pbc_y, non_pbc_z)
	non_pbc_vecs = list(non_pbc_vecs)

	# find norm of each vector of non pbc changes
	for vector in non_pbc_vecs:
		normer = norm(vector)
	        norms.append(normer)

	# find minimum of these distances for each time step then add onto final list to plot graph 
	minimum_dist = min(norms)
	min_1_2.append((trajectory.time, minimum_dist))




