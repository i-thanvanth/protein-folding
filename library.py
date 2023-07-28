def transcription(seq):

	#converts DNA Sequence to mRNA sequence
	#The conversion is basically that thymine is replaced with uracil during the transcription process

	mrna_seq = ""
	seq = seq.upper()
	for i in seq:
		if i not in "ATGC":
			return "INVALID SEQUENCE"
		else:
			if i == 'T':
				mrna_seq+='U'
			else:
				mrna_seq+=i
			#mrna_seq+=rna[dna.index(i)]
	return mrna_seq

def invert(seq): #basically gives the antisense strand for a said sequence
	dna,a_s = 'ATGC','TACG'
	return "".join([a_s[dna.index(i)] for i in seq[::-1]]) if len([i for i in seq if i in "ATGC"]) == len(seq) else "INVALID SEQUENCE"
#print(invert('ATGC'))

'''
Standard RNA codon table

the leader sequence ends before the codon AUG which corresponds to the amino acid Methionine
The stop codons are UAA, UGA, UAG (marked by X in this code, might change in the future)
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

	#the translation starts at codon AUG.
	if seq == "INVALID SEQUENCE":
		return "INVALID SEQUENCE"
	seq = seq.upper()
	start = seq.find("AUG")
	if start != -1:
		polypeptide_sequence = ''
		seq = [seq[i:i+3] for i in range(start,len(seq),3)]
		for i in seq:
			# stop sequences, stop the reading if it encounters a stop protein, also could just add a * but for now it just stops
			if i in ['UAA','UAG','UGA']:
				break
			if len(i) == 3:
				polypeptide_sequence+=codon[i]
		return polypeptide_sequence

#print(translation((transcription("AGGACGGGCTAACTCCGCTCGTCACAAAGCGCAATGCAGCTATGGCAGATGTTCATGCCG"))))
#print(translation((transcription("TACATGCCATACGAGACGAGCGCGCCTAAGCGGCGCAGACTCATGGTCATT"))))
#print(transcription("ATGGCCGGTTATTAAGCA"))
#print(translation("AUGUUUUGG"))
#print(translation("ggaugcccaaauaa"))
#print(transcription("TTATGCATC"))
#print(translation('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'))

adenovirus5 = open('porcine_adenovirus_5.txt')
av5genome  = "".join([i for i in adenovirus5.read()  if i in 'ATGC'])
#print(translation(transcription(av5genome[418-1:939])))