class NotYetImplementedError(Exception):
  def __init__(self, message):
    super().__init__(message)

class ValueMissingError(Exception):
  def __init__(self, message):
    super().__init__(message)
