import sys
sys.path.append("..")

import pyjugex
import pytest

def test_analysis():
  analysis=pyjugex.analysis()
  with pytest.raises(pyjugex.ValueMissingError):
    analysis.run()
