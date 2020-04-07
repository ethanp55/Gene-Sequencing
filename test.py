from Node import *

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

cost_matrix = []
sequence1 = 'exponential'
sequence2 = 'polynomial'

length1 = len(sequence1)
length2 = len(sequence2)

if abs(length1 - length2) > MAXINDELS:
	exit()

offset = 0
jmax = 4
joffset = 0

for i in range(0, length1 + 1):
	cost_matrix.append([])

	if i > MAXINDELS:
		offset = offset + 1

	if offset > 1:
		joffset = joffset + 1

	if jmax > length2 + 1:
		jmax = length2 + 1

	for j in range(0 + offset, jmax):

		if i == 0:
			if j == 0:
				# corner
				value = 0
				previous = None
				edit_type = None
				char_1 = None
				char_2 = None

			else:
				# 1st row
				value = cost_matrix[i][j - 1].value + INDEL
				previous = cost_matrix[i][j - 1]
				edit_type = "insert"
				char_1 = None
				char_2 = sequence2[j - 1]

		elif j == 0:
			# 1st column
			value = cost_matrix[i - 1][j].value + INDEL
			previous = cost_matrix[i - 1][j]
			edit_type = "delete"
			char_1 = sequence1[i - 1]
			char_2 = None

		else:
			# anything else
			char_1 = sequence1[i - 1]
			char_2 = sequence2[j - 1]

			if sequence1[i - 1] == sequence2[j - 1] :
				dval = MATCH

			else:
				dval = SUB

			if i - j == MAXINDELS:
				# bottom-left corners
				value = min(cost_matrix[i - 1][j - joffset].value + INDEL, cost_matrix[i - 1][j - 1 - joffset].value + dval)

				if value == cost_matrix[i - 1][j - joffset].value + INDEL:
					previous = cost_matrix[i - 1][j - joffset]
					edit_type = "delete"

				else:
					previous = cost_matrix[i - 1][j - 1 - joffset]
					edit_type = "sub"

			elif j - i == MAXINDELS:
				# top-right corner
				value = min(cost_matrix[i - 1][j - 1 - joffset].value + dval, cost_matrix[i][j - 1 - offset].value + INDEL)

				if value == cost_matrix[i - 1][j - 1 - joffset].value + dval:
					previous = cost_matrix[i - 1][j - 1 - joffset]
					edit_type = "sub"

				else:
					previous = cost_matrix[i][j - 1 - offset]
					edit_type = "insert"

			else:
				# middle values
				value = min(cost_matrix[i - 1][j - joffset].value + INDEL, cost_matrix[i - 1][j - 1 - joffset].value + dval, cost_matrix[i][j - 1 - offset].value
				+ INDEL)

				if value == cost_matrix[i - 1][j - joffset].value + INDEL:
					previous = cost_matrix[i - 1][j - offset]
					edit_type = "delete"

				elif value == cost_matrix[i - 1][j - 1 - joffset].value + dval:
					previous = cost_matrix[i - 1][j - 1 - joffset]
					edit_type = "sub"

				else:
					previous = cost_matrix[i][j - 1 - offset]
					edit_type = "insert"

		cost_matrix[i].append(Node(value, previous, edit_type, char_1, char_2))
	jmax = jmax + 1

for i in range(0, len(cost_matrix)):
	for j in range(0, len(cost_matrix[i])):
		print(str(cost_matrix[i][j].value) + ' ', end = '')

	print()
