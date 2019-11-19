from main import PyjugexAnalysis
from util import ValueMissingError

import pytest

def test_PyjugexAnalysis():
  analysis=PyjugexAnalysis()
  with pytest.raises(ValueMissingError):
    analysis.differential_analysis()
