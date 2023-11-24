# similar to contact_frac script - focused on lysine and glutamine residue contacts

saltbridge = "(resname LYS GLU) and (name CA CA)"
bridge1 = peptide1.select_atoms(saltbridge)
bridge2 = peptide2.select_atoms(saltbridge)
bridge3 = peptide3.select_atoms(saltbridge)

a = 10
ca = contacts_within_cutoff(u, bridge2, bridge3, radius=a)
np.savetxt("contacts_LYS+GLU_%dA.txt"%a, ca)
for item in ca:
	if item[1] > 0:
		print(item[1])
y = []
for item in ca:
	if item[1] > 0:
		y.append(item[1])
print(len(ca))
print(len(y))

