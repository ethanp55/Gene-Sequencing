#!/usr/bin/python3

from which_pyqt import PYQT_VER
from Node import *
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass

# This is the method called by the GUI.  _sequences_ is a list of the ten sequences, _table_ is a
# handle to the GUI so it can be updated as you find results, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def align( self, sequences, table, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		results = []

		# Time and space complexity of either O(nm) (if running the unrestricted alignment) or O(kn), where k = 7 (if
		# running the banded alignment)
		# This outer loop iterates 10 times
		for i in range(len(sequences)):
			jresults = []

			# This inner loop iterates 10 times
			for j in range(len(sequences)):

				if j < i:
					s = {}

				# Use the banded alignment
				elif banded:
					result = self.banded_align(sequences[i], sequences[j])

				# Use the unrestricted alignment
				else:
					result = self.unrestricted_align(sequences[i], sequences[j])

				# Send the score/cost and the alignments to the GUI
				score = result[0]
				alignment_1 = result[1]
				alignment_2 = result[2]

				s = {'align_cost':score, 'seqi_first100':alignment_1, 'seqj_first100':alignment_2}
				table.item(i,j).setText('{}'.format(int(score) if score != math.inf else score))
				table.repaint()
				jresults.append(s)
			results.append(jresults)
		return results

	# Helper function for performing an unrestricted alignment of the two gene sequences
	# O(nm) space complexity
	def unrestricted_align(self, sequence1, sequence2):
		cost_matrix = []

		# Check if the align length is too large for the sequences
		# Get the appropriate align lengths
		if len(sequence1) < self.MaxCharactersToAlign or len(sequence2) < self.MaxCharactersToAlign:
			length1 = len(sequence1)
			length2 = len(sequence2)

		else:
			length1 = self.MaxCharactersToAlign
			length2 = self.MaxCharactersToAlign

		# Time complexity of O(nm) because this outer for loop iterates n times and the inner for loop iterates m times.
		# Space complexity of O(nm) because there are nm iterations, and at each iteration we add a new value/node to
		# the cost matrix
		for i in range(0, length1 + 1):
			# Append a new list/row to the cost matrix
			cost_matrix.append([])

			# This loop will iterate m times
			for j in range(0, length2 + 1):
				# There are various cases for determining the next value/node to put into the cost matrix
				# Check the various cases and determine which node to add to the matrix
				if i == 0:
					if j == 0:
						# Corner
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
					# Anything else
					char_1 = sequence1[i - 1]
					char_2 = sequence2[j - 1]

					# Check if there is a match
					if sequence1[i - 1] == sequence2[j - 1]:
						dval = MATCH

					else:
						dval = SUB

					value = min(cost_matrix[i - 1][j].value + INDEL, cost_matrix[i - 1][j - 1].value + dval, cost_matrix[i][j - 1].value + INDEL)

					if value == cost_matrix[i - 1][j].value + INDEL:
						previous = cost_matrix[i - 1][j]
						edit_type = "delete"

					elif value == cost_matrix[i - 1][j - 1].value + dval:
						previous = cost_matrix[i - 1][j - 1]
						edit_type = "sub"

					else:
						previous = cost_matrix[i][j - 1]
						edit_type = "insert"

				# Once we have determined which node to add to the cost matrix, add it
				cost_matrix[i].append(Node(value, previous, edit_type, char_1, char_2))

		# Get the final node (the bottom-right node) from the matrix
		final_node = cost_matrix[length1][length2]

		# Get the final cost and extract the alignments of the two sequences
		cost = final_node.value
		alignment_1, alignment_2 = self.extract_alignments(final_node)

		return cost, alignment_1, alignment_2

	# Helper function for performing a banded alignment of the two gene sequences
	# O(kn) time and space complexity
	def banded_align(self, sequence1, sequence2):
		cost_matrix = []

		# Check if the align length is too large for the sequences
		# Get the appropriate align lengths
		if len(sequence1) < self.MaxCharactersToAlign or len(sequence2) < self.MaxCharactersToAlign:
			length1 = len(sequence1)
			length2 = len(sequence2)

		else:
			length1 = self.MaxCharactersToAlign
			length2 = self.MaxCharactersToAlign

		# Ensure that the two lengths do not differ by a significant amount (if so, we cannot perform the banded
		# alignment)
		if abs(length1 - length2) > MAXINDELS:
			return math.inf, "No Alignment Possible", "No Alignment Possible"

		# Initialize offsets (used for running in linear time and space complexities)
		offset = 0
		jmax = MAXINDELS + 1
		joffset = 0

		# Time complexity of O(kn) because the outer for loop iterates n times and the inner for loop will iterate a
		# maximum of 7, or k, times
		# Space complexity of O(kn) because we iterate kn times, and at each iteration we add a new value/node to the
		# cost matrix
		for i in range(0, length1 + 1):
			# Append a new list/row to the cost matrix
			cost_matrix.append([])

			# Logic for handling the offsets
			if i > MAXINDELS:
				offset = offset + 1

			if offset > 1:
				joffset = joffset + 1

			if jmax > length2 + 1:
				jmax = length2 + 1

			# This loop will iterate a maximum of 7, or k, times because the maximum difference between offset and jmax
			# is 7.
			for j in range(0 + offset, jmax):
				# There are various cases for determining the next value/node to put into the cost matrix
				# Check the various cases and determine which node to add to the matrix
				if i == 0:
					if j == 0:
						# Corner
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
					# Anything else
					char_1 = sequence1[i - 1]
					char_2 = sequence2[j - 1]

					# Check if there is a match
					if sequence1[i - 1] == sequence2[j - 1]:
						dval = MATCH

					else:
						dval = SUB

					if i - j == MAXINDELS:
						# Bottom-left corners
						value = min(cost_matrix[i - 1][j - joffset].value + INDEL, cost_matrix[i - 1][j - 1 - joffset].value + dval)

						if value == cost_matrix[i - 1][j - joffset].value + INDEL:
							previous = cost_matrix[i - 1][j - joffset]
							edit_type = "delete"

						else:
							previous = cost_matrix[i - 1][j - 1 - joffset]
							edit_type = "sub"

					elif j - i == MAXINDELS:
						# Top-right corners
						value = min(cost_matrix[i - 1][j - 1 - joffset].value + dval, cost_matrix[i][j - 1 - offset].value + INDEL)

						if value == cost_matrix[i - 1][j - 1 - joffset].value + dval:
							previous = cost_matrix[i - 1][j - 1 - joffset]
							edit_type = "sub"

						else:
							previous = cost_matrix[i][j - 1 - offset]
							edit_type = "insert"

					else:
						# Middle values
						value = min(cost_matrix[i - 1][j - joffset].value + INDEL, cost_matrix[i - 1][j - 1 - joffset].value + dval, cost_matrix[i][j - 1 - offset].value + INDEL)

						if value == cost_matrix[i - 1][j - joffset].value + INDEL:
							previous = cost_matrix[i - 1][j - offset]
							edit_type = "delete"

						elif value == cost_matrix[i - 1][j - 1 - joffset].value + dval:
							previous = cost_matrix[i - 1][j - 1 - joffset]
							edit_type = "sub"

						else:
							previous = cost_matrix[i][j - 1 - offset]
							edit_type = "insert"

				# Once we have determined which node to add to the cost matrix, add it
				cost_matrix[i].append(Node(value, previous, edit_type, char_1, char_2))

				# If we are at the end of the matrix, extract the final node
				if i == length1 and j == jmax - 1:
					final_node = Node(value, previous, edit_type, char_1, char_2)

			# Increment jmax
			jmax = jmax + 1

		# Get the final cost and extract the alignments of the two sequences
		cost = final_node.value
		alignment_1, alignment_2 = self.extract_alignments(final_node)

		return cost, alignment_1, alignment_2

	# Helper function for extracting the optimal alignments of the two gene sequences
	# O(n) time and space complexity
	def extract_alignments(self, final_node):
		# Initialize empty lists for the two alignments
		r_alignment_1 = []
		r_alignment_2 = []

		# Use an iterator and start at the final node in the matrix
		node_iter = final_node

		# Time complexity of O(n) since we start at the bottom-right node in the matrix and iterate backwards up until
		# the top-left node
		# Space complexity of O(n) since we are building the alignments, and they will have length n
		while node_iter is not None:
			# Determine the type of the node and append the appropriate characters to the alignments
			if node_iter.edit_type == "match" or node_iter.edit_type == "sub":
				r_alignment_1.append(node_iter.char_1)
				r_alignment_2.append(node_iter.char_2)

			elif node_iter.edit_type == "insert":
				r_alignment_1.append("-")
				r_alignment_2.append(node_iter.char_2)

			elif node_iter.edit_type == "delete":
				r_alignment_1.append(node_iter.char_1)
				r_alignment_2.append("-")

			else:
				break

			# Iterate to the next node
			node_iter = node_iter.previous

		# Convert the alignments to strings and reverse them
		string_1 = "".join(r_alignment_1)
		string_2 = "".join(r_alignment_2)

		alignment_1 = string_1[::-1]
		alignment_2 = string_2[::-1]

		# Return the first 100 characters of the alignments
		return alignment_1[0:100], alignment_2[0:100]
