class Mutation:
   def __init__(self, wildtype_AN, secuence_AN, old_aa, position_1D, new_aa):
      self.wildtype_AN=wildtype_AN
      self.secuence_AN=secuence_AN
      self.old_aa=old_aa
      self.position_1D=position_1D
      self.new_aa=new_aa
   def __repr__(self):
      return repr((self, wildtype_AN, secuence_AN, old_aa, position_1D, new_aa))
