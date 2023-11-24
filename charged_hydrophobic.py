# algorithm made to find the charged and hydrophobic residue percentages within a peptide

w, e = [], []
for item in q:
    if 0 <= float(item[1]) <= thresh:
        w.append(item)
    else:
        e.append(item)
print(w)
u = []
for item in w:
    u.append(item[0])
print(u)
hydro = ["A", "V", "L", "W", "G", "I", "M", "F"]
charged = ["K", "E"]
b = set(u) & set(hydro)
r = set(hydro)
n = [x for x in u if x in r]
print(len(n))
d = set(u) & set(charged)
v = set(charged)
c = [x for x in u if x in v]
a = len(c)
if a == 0:
    c = 1
else:
    c = c
print("hydrophobic percentage =", len(u) / len(n) * 100)
print("set of residues =", b)
print("residues", n)
print("charged percentage =", len(u) / len(c) * 100)
print("set of residues =", d)
print("residues", c)
