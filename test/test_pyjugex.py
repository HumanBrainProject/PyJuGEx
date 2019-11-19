from pyjugex import PyjugexAnalysis, ValueMissingError

import pytest

def test_PyjugexAnalysis():
  analysis=PyjugexAnalysis()
  with pytest.raises(ValueMissingError):
    analysis.differential_analysis()
