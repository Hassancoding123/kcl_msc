# algorithm to find the minimum distance between carbon alpha atoms in two different peptides
# aim = to see if peptides would bond by remaining under a threshold minimum distance

min_1_2 = []
CA1 = peptide1.select_atoms("name CA")
CA2 = peptide2.select_atoms("name CA")
CA3 = peptide3.select_atoms("name CA")
for ts in trajectory:
    pos_1 = CA1.positions
    pos_2 = CA2.positions
    pos_3 = CA3.positions
    all_possible_pos = product(pos_1, pos_2)
    norms = []
    changes = [abs(a-b) for a,b in all_possible_pos]
    changes_x = [abs(item[0]) for item in changes]
    changes_y = [abs(item[1]) for item in changes]
    changes_z = [abs(item[2]) for item in changes]
    non_pbc_x = []
    non_pbc_y = []
    non_pbc_z = []
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
    non_pbc_vecs = zip(non_pbc_x, non_pbc_y, non_pbc_z)
    non_pbc_vecs = list(non_pbc_vecs)
    for vector in non_pbc_vecs:
        normer = norm(vector)
        norms.append(normer)
    minimum_dist = min(norms)
    min_1_2.append((trajectory.time, minimum_dist))

res1 = "ALA"
res2 = ["ALA", "GLY", "LEU", "ILE", "PHE", "MET", "TRP", "VAL"]
a = 7
for item in res2:
	saltbridge1 = "(resname %s)"%res1
	saltbridge2 = "(resname %s)"%item
	bridge2 = peptide2.select_atoms(saltbridge1)
	bridge3 = peptide3.select_atoms(saltbridge2)
	ca = contacts_within_cutoff(u, bridge2, bridge3, radius=a)
	x, y = zip(*ca)
	p = sum(y)
	m = p/len(y)
	print("N2_{a}-N3_{b} + thresh = {c}".format(a=res1, b=item, c=a))
	s = []
	for i in y:
		if i > 0:
			s.append(i)
	print("frac =", len(s)/len(y))
	print("av =", m)

