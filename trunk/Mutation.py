class Mutation:
   def __init__(self, wildtype_AN, sequence_AN, old_aa, position_1D, new_aa,country):
      self.wildtype_AN=wildtype_AN	# wildtype Accesion Number
      self.sequence_AN=sequence_AN	# sequence to analize Accesion Number
      self.old_aa=old_aa				# old amino acid
      self.position_1D=position_1D	# position in the sequence in 1D notation
      self.new_aa=new_aa				# new amino acid (mutation)
      self.country=country				# found country
   def __repr__(self):
      return repr((self.wildtype_AN, self.sequence_AN, self.old_aa, self.position_1D, self.new_aa, self.country))
