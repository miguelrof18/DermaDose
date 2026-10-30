import numpy as np
import pytest
from skindose_verifier.factor_tables import RectilinearFactorTable

def test_combined_interpolation_preserves_weighted_product():
    rows=[dict(row_id='a',x=0,k_bs=1,k_f=2),dict(row_id='b',x=1,k_bs=3,k_f=4)]
    r=RectilinearFactorTable(rows,['x'],interpolate_combined=True).evaluate(x=np.array([0,.5,1]))
    assert r.k_bs*r.k_f==pytest.approx([2,7,12])
    assert r.audit['interpolation_quantity']=='combined_C'
    legacy=RectilinearFactorTable(rows,['x']).evaluate(x=.5)
    assert float(legacy.k_bs*legacy.k_f)==6
