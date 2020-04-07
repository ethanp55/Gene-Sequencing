# Node class
# Nodes are added to the cost matrix and give information needed for extracting the alignments
class Node:
    def __init__(self, value, previous, edit_type, char_1, char_2):
        # The cost
        self.value = value
        # The previous node
        self.previous = previous
        # The edit type (insert, delete, sub, or match)
        self.edit_type = edit_type
        # The corresponding character of sequence 1
        self.char_1 = char_1
        # The corresponding character of sequence 2
        self.char_2 = char_2
