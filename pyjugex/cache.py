from .error import ValueMissingError

class MemoryCache():
  def __init__(self):
    self.store = dict()

  def get_from_key(self, key=None):
    if key is None:
      raise ValueMissingError('key is required')
    return self.store.get(key, None)

  def store_key_value(self, key=None, value=None):
    if key is None:
      raise ValueMissingError('key is required')
    self.store[key] = value