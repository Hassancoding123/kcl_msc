# similar to peptide min dist script except now using carbon alphas AND residue names to find residue-residue contacts during simulation between 2 peptides

res_min_dists = []
res_no = 3
for ts in trajectory:
    res_peptide2 = []
    res_peptide3 = []
    for i in range (2,50):
        a = peptide2.select_atoms("name CA and resid %d"%i).positions
        b = peptide3.select_atoms("name CA and resid %d"%i).positions
        res_peptide3.append(b)
        res_peptide2.append(a)
    real_res = res_peptide2[res_no - 2]
    CA_3 = peptide3.select_atoms("name CA")
    pos_3 = CA_3.positions
    m = product(pos_3, real_res)
    m = list(m)
    changes = [abs(a-b) for a,b in m]
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
    norms = []
    for vector in non_pbc_vecs:
        normer = norm(vector)
        norms.append(normer)
    minimum_dist = min(norms)
    res_min_dists.append((trajectory.time, minimum_dist))

np.savetxt("p2_res%d_boundto_p3.txt"%res_no, res_min_dists)
