def transcription(seq):

	#converts DNA Sequence to mRNA sequence
	#The conversion is as follows:
	#A -> U
	#T -> A
	#G -> C
	#C -> G

	mrna_seq = ""
	dna = 'ATGC'
	rna = 'UACG'
	for i in seq:
		if i not in "ATGC":
			return "INVALID SEQUENCE"
		else:
			mrna_seq+=rna[dna.index(i)]
	return mrna_seq

'''

Standard RNA codon table

the leader sequence ends before the codon AUG which corresponds to the amino acid Methionine
The stop codons are UAA, UGA, UAG
the other codons are 
A / GCU, GCC, GCA, GCG -> Alanine
I / AUU, AUC, AUA -> Isoleucine
R / CGU, CGC, CGA, CGG, AGA, AGG -> Arginine
L / CUU, CUC, CUA, CUG, UUA, UUG -> Leucine
N / AAU, AAC -> Asparagine
K / AAA, AAG -> Lycine
D / GAU, GAC -> Aspartic Acid
M / AUG -> Methionine (again, starter codon)
F / UUU, UUC -> Phenylalanine
C / UGU, UGC -> Cysteine
P / CCU, CCC, CCA, CCG -> Proline
Q / CAA, CAG -> Glutamine
S / UCU, UCC, UCA, UCG, AGU, AGC -> Serine
E / GAA, GAG -> Glutamic Acid
T / ACU, ACC, ACA, ACG -> Threonine
W / UGG -> Tryptophan
G / GGU, GGC, GGA, GGG -> Glycine
Y / UAU, UAC -> Tyrosine
H / CAU, CAC -> Histidine
V / GUU, GUC, GUA, GUG -> Valine
'''

codons_text = """A / GCU, GCC, GCA, GCG
I / AUU, AUC, AUA
R / CGU, CGC, CGA, CGG, AGA, AGG
L / CUU, CUC, CUA, CUG, UUA, UUG
N / AAU, AAC
K / AAA, AAG
D / GAU, GAC
M / AUG
F / UUU, UUC
C / UGU, UGC
P / CCU, CCC, CCA, CCG
Q / CAA, CAG
S / UCU, UCC, UCA, UCG, AGU, AGC
E / GAA, GAG
T / ACU, ACC, ACA, ACG
W / UGG
G / GGU, GGC, GGA, GGG
Y / UAU, UAC
H / CAU, CAC
V / GUU, GUC, GUA, GUG
X / UAA, UGA, UAG
""".strip()

# creating the RNA Codon dictionary to look up amino acids during translation
codons_text = codons_text.split('\n')
codons_text = [i.split('/') for i in codons_text]
codons_text = [[i[0][:-1],i[1][1:].split(', ')] for i in codons_text]
codon = {}
for i in codons_text:
	for j in i[1]:
		codon[j] = i[0]
#print(codon)

def translation(seq):
	#roughly works, still have to code in how to identify starting codon and stop codon
	polypeptide_sequence = ''
	seq = [seq[i:i+3] for i in range(0,len(seq),3)]
	for i in seq:
		polypeptide_sequence+=codon[i]
	return polypeptide_sequence

print(translation(transcription("TACTAG")))

