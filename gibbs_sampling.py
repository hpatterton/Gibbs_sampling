import random
import copy
import matplotlib.pyplot as plt
import numpy.random


def generate_random_sequences(number_of_sequences, shortest, longest, symbol_alphabet):
	symbol_max = len(symbol_alphabet)
	sequences = []
	for i in range(number_of_sequences):
		random.seed()
		length_of_sequence = random.randint(shortest,longest)
		temp_seq = ''
		for j in range(length_of_sequence):
			temp_seq += symbol_alphabet[random.randint(0, symbol_max - 1)]
		sequences.append(temp_seq)
	return sequences

def add_defined_site_randomly_to_sequences(sequences, sequence_of_site):
	number_of_sequences = len(sequences)
	length_of_site = len(sequence_of_site)
	sequence_with_sites = []
	site_positions = []
	random.seed()
	for i in range(number_of_sequences):
		index = random.randint(0,len(sequences[i]))
		sequence_with_sites.append(sequences[i][:index]+sequence_of_site+sequences[i][index+length_of_site:])
		site_positions.append(index)
	return sequence_with_sites,site_positions

def select_sites_randomly(sequences,number_of_sites, site_width):
	random.seed()
	site_sequences =  []
	site_indexes = []
	background_sequences = []
	for i in range(number_of_sites):
		for j in range(len(sequences)):
			index = random.randint(0,len(sequences[j])-site_width)
			site_indexes.append(index)
			site_sequences.append(sequences[j][index:index+site_width])
			background_sequences.append(sequences[j][:index] + sequences[j][index + site_width:])
	return site_sequences,site_indexes,background_sequences

def remove_one_sequence(list_of_sequences, list_of_site_indexes, site_width):
	number_of_sequences = len(list_of_sequences)
	random.seed()
	index_to_remove = random.randint(0,number_of_sequences-1)
	removed_sequence = list_of_sequences[index_to_remove]
	site_start = list_of_site_indexes[index_to_remove]
	site_end = list_of_site_indexes[index_to_remove]+site_width
	removed_site_sequence = removed_sequence[site_start:site_end]
	list_of_sequences.pop(index_to_remove)
	list_of_site_indexes.pop(index_to_remove)
	return list_of_sequences,list_of_site_indexes,removed_sequence,removed_site_sequence

def calculate_q_matrix(site_sequences,pseudocount,symbol_alphabet):
	number_of_sequences = len(site_sequences)
	number_of_symbols = len(symbol_alphabet)
	site_width = len(site_sequences[0])
	count_matrix = []
	count_row = [0]*site_width
	for i in range(number_of_symbols):
		count_matrix.append(count_row.copy())
	q_matrix = copy.deepcopy(count_matrix)
	# calculate count matrix
	for position in range(site_width):
		for symbol_index in range(number_of_symbols):
			count = 0
			for index in range(number_of_sequences):
				if site_sequences[index][position] == symbol_alphabet[symbol_index]:
					count += 1
			count_matrix[symbol_index][position] = count
	#use count matrix to calculate q_matrix
	for position in range(site_width):
		for symbol_index in range(number_of_symbols):
			q_matrix[symbol_index][position] = count_matrix[symbol_index][position]+pseudocount/(number_of_sequences+number_of_sequences*pseudocount)
	return count_matrix, q_matrix

def calculate_background_matrix(background_sequences,symbol_alphabet):
	total_length = 0
	number_of_symbols = [0]*len(symbol_alphabet)
	number_of_sequences = len(background_sequences)
	for i in range(number_of_sequences):
		total_length += len(background_sequences[i])
		for j in range(len(symbol_alphabet)):
			number_of_symbols[j] += background_sequences[i].count(symbol_alphabet[j])
	b_matrix = []
	for index in range(len(symbol_alphabet)):
		b_matrix.append(number_of_symbols[index]/total_length)
	return b_matrix

def calculate_Ax(site_sequence, q_matrix, background_matrix, symbol_alphabet):
	site_length = len(site_sequence)
	symbol_indexes = []
	for i in range(site_length):
		symbol_indexes.append(symbol_alphabet.index(site_sequence[i]))
	numerator = 1
	demoninator = 1
	for i in range(site_length):
		numerator *= q_matrix[symbol_indexes[i]][i]
		demoninator *= background_matrix[symbol_indexes[i]]
	Ax = numerator/demoninator
	return Ax

def calculate_Ax_matrix(sequence, site_width, q_matrix, b_matrix, symbol_alphabet):
	seq_length = len(sequence)
	Ax = []
	for i in range(0, seq_length - site_width + 1):
		Ax.append(calculate_Ax(sequence[i:i + site_width], q_matrix, b_matrix, symbol_alphabet))
	return Ax

def Ax_to_probability_distribution(Ax):
	number_of_entries = len(Ax)
	sum = 0
	for i in range(number_of_entries):
		sum += Ax[i]
	Ax_normalised = []
	for i in range(number_of_entries):
		Ax_normalised.append(Ax[i]/sum)
	return Ax_normalised

def sample_from_probability_distribution(Ax_normalised):
	numpy.random.seed()
	a = []
	for i in range(len(Ax_normalised)):
		a.append(i)
	new_index = numpy.random.choice(a, 1,p=Ax_normalised)
	return new_index

def get_site_and_background_sequences(list_of_sequences, list_of_site_indexes, width_of_site):
	number_of_sequences = len(list_of_sequences)
	background_sequences = []
	site_sequences = []
	for i in range(number_of_sequences):
		start = list_of_site_indexes[i]
		end = list_of_site_indexes[i]+width_of_site
		site_sequences.append(list_of_sequences[i][start:end])
		background_sequences.append(list_of_sequences[i][:start]+list_of_sequences[i][end:])
	return site_sequences,background_sequences

def set_random_site_positions(sequences, width_of_site):
	number_of_sequences = len(sequences)
	site_positions = []
	random.seed()
	for i in range(number_of_sequences):
		sequence_length = len(sequences[i])
		site_positions.append(random.randint(0,sequence_length-width_of_site))
	return site_positions

def sample_from_distribution(distribution):
	length_of_distribution = len(distribution)
	cummulative = [0]*length_of_distribution
	cummulative[0] = distribution[0]
	for i in range(1,length_of_distribution):
		cummulative[i] = cummulative[i-1] + distribution[i]
	value = random.randint(0,cummulative[-1])
	index = 0
	while value > cummulative[index]:
		index += 1
	return index

def get_consensus_sequence(motifs, site_width, alphabet):
	number_of_symbols = len(alphabet)
	number_of_motifs = len(motifs)
	row = [0]*site_width
	count = []
	for i in range(number_of_symbols):
		count.append(row.copy())
	for position in range(site_width):
		for symbol_index in range(number_of_symbols):
			for motif_index in range(number_of_motifs):
				if motifs[motif_index][position] == alphabet[symbol_index]:
					count[symbol_index][position] += 1
	consensus = ''
	for position in range(site_width):
		max_count = max([count[index][position] for index in range(number_of_symbols)])
		symbol_index = [count[index][position] for index in range(number_of_symbols)].index(max_count)
		consensus += alphabet[symbol_index]
	return consensus

def score_motifs(list_of_sites, symbol_alphabet):
	length_of_site = len(list_of_sites[0])
	number_of_motifs = len(list_of_sites)
	number_of_symbols = len(symbol_alphabet)
	row = [0]*length_of_site
	counts = []
	for i in range(number_of_symbols):
		counts.append(row.copy())
	for position in range(length_of_site):
		for symbol_index in range(number_of_symbols):
			for motif in list_of_sites:
				if motif[position] == symbol_alphabet[symbol_index]:
					counts[symbol_index][position] += 1
	# get maximum counts in each column
	score = 0
	for motif_position in range(length_of_site):
		max_symbol = max([counts[symbol_index][motif_position]] for symbol_index in range(number_of_symbols))
		# the score is the number of symbols in a column that is not the maximum symbol
		score += (number_of_motifs-max_symbol[0])
	return score

####################################################################################################################
#
#
#
####################################################################################################################

# variables to use

site_width = 8
number_of_sites = 1
number_of_sequences = 6
shortest_sequence = 200
longest_sequence = 200
number_of_sequences_to_remove = 1
pseudocount = 0.25
number_of_iterations = 10000
sequence_of_site = 'GGCTAGCC'
score_array = []	# array of best motif list scores

molecule_type = 'D'
if molecule_type == 'D':
	symbol_alphabet = ['G', 'A', 'T', 'C']
else:
	symbol_alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
					   'W', 'Y']

# Generate random sequences, and insert the defined site randomly in each
sequences = generate_random_sequences(number_of_sequences, shortest_sequence, longest_sequence, symbol_alphabet)
sequence_with_sites,real_site_positions = add_defined_site_randomly_to_sequences(sequences, sequence_of_site)

#select random positions, and recover the selected site sequences
site_positions = set_random_site_positions(sequence_with_sites, site_width)
site_sequences = []
for i in range(len(sequence_with_sites)):
	site_sequences.append(sequence_with_sites[i][site_positions[i]:site_positions[i]+site_width])

# calculate the score for the list of site sequences
best_score = score_motifs(site_sequences, symbol_alphabet)
best_site_sequences = site_sequences.copy()


# iterate to find the list of sites with lowest score
iteration = 0
while iteration < number_of_iterations:

	#remove one sequence randomly
	shortened_list_of_sequences,shortened_list_of_sites,removed_sequence,removed_site_sequence = remove_one_sequence(sequence_with_sites, site_positions,site_width)
	site_sequences,background_sequences = get_site_and_background_sequences(shortened_list_of_sequences, shortened_list_of_sites, site_width)

	# calculate count and qij matrix
	c_matrix, q_matrix = calculate_q_matrix(site_sequences,pseudocount,symbol_alphabet)
	b_matrix = calculate_background_matrix(background_sequences,symbol_alphabet)

	# Calculate the A = Q/P value for each setting in the removed sequence
	Ax = calculate_Ax_matrix(removed_sequence, site_width, q_matrix, b_matrix, symbol_alphabet)
	Ax_normalised = Ax_to_probability_distribution(Ax)

	# Sample positions from the normalised Ax distribution
	position_of_good_site = numpy.random.choice(len(Ax), size=None, replace=True, p=Ax_normalised)

	# Get the sequence of the high probability site in the removed sequence
	good_site = removed_sequence[position_of_good_site:position_of_good_site+site_width]

	# Add the removed sequence back to the list of sequences
	shortened_list_of_sequences.append(removed_sequence)
	sequence_with_sites = shortened_list_of_sequences.copy()

	# randomly select new sites for the list of sequences, but retaining the good site
	# for the sequence added back
	site_positions = set_random_site_positions(sequence_with_sites, site_width)
	site_positions[-1] = position_of_good_site

	# Add the good site sequence from the removed sequence to the other site sequences
	# and score the list of site sequences
	site_sequences.append(good_site)
	score = score_motifs(site_sequences, symbol_alphabet)

	# if the score is better than the current best score
	# set the list of best site sequences to the current list of site sequences
	if score < best_score:
		best_score = score
		best_site_sequences = site_sequences.copy()

		# copy the current best score to a list
		score_array.append(best_score)

	print(iteration,best_site_sequences)

	iteration += 1
print('consensus =',get_consensus_sequence(best_site_sequences, site_width, symbol_alphabet))


# plot best_score
fig, ax = plt.subplots()
ax.plot(score_array, linewidth=2.0)
plt.show()


