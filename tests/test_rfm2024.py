"""Changes motivated by the RFM 2024 verification (Gañán Mora et al., Rev Fis Med 2024)."""
import math

import numpy as np
import pytest

from skindose_verifier.domain import InputError, IrradiationEvent
from skindose_verifier.surface_factors import make_event_provider

AL, CU = ("SRT", "C-120F9"), ("SRT", "C-127F9")


def event(*, kvp=70., cu=0.4, al=1.0, width=300., height=389., sid=1181.4, primary=0., secondary=0.,
          area_derived=False, uid="e1"):
    stack = [{"material_scheme": s, "material_code": c, "minimum_mm": t, "maximum_mm": t, "selected_mm": t}
             for (s, c), t in ((AL, al), (CU, cu)) if t]
    quality = {"filtration_mm_cu": "copper" if cu else "no_copper_filter"}
    if area_derived:
        quality["field_size"] = "square_equivalent_from_collimated_field_area"
    return IrradiationEvent(uid=uid, plane="113622", ka_r_gy=0.01, kap_gy_m2=3e-4, kvp=kvp, filtration_mm_cu=cu,
                            primary_deg=primary, secondary_deg=secondary, table_mm=(0., 0., 0.),
                            table_system="isocenter", sid_mm=sid, sod_mm=None, srp_mm=None,
                            field_width_mm=width, field_height_mm=height, quality=quality,
                            original={"filter_stack": stack}, field_reference="detector")


def provider(**changes):
    spec = {"provider": "benmakhlouf_2011", "inherent_al_mm": 3.1, "target": "W", "field_shape_policy": "equal_area",
            "out_of_domain": "error", "incident_kerma_basis": "free_in_air_without_backscatter",
            "quality_basis": "event_kvp_filtration_proxy_at_surface", "provenance": "test",
            "allow_water_phantom_approximation": True, "aluminium_policy": "fixed_plus_event", **changes}
    return make_event_provider(spec)


def test_aluminium_interpolation_is_linear_between_adjacent_families_only():
    ev = event(kvp=60., cu=0., al=1.8, width=300., height=300.)      # Canon case 6: 3.1 + 1.8 = 4.9 mm Al
    with pytest.raises(InputError, match="linear_between_families"):
        provider().evaluate_event(ev, 270.)
    pair = provider(aluminium_interpolation="linear_between_families").evaluate_event(ev, 270.)
    data = provider().dataset
    side = math.sqrt(270.)
    lo = data.evaluate("clinical", inherent_al_mm=4.1, kvp=60., added_cu_mm=0., field_side_cm=side)
    hi = data.evaluate("clinical", inherent_al_mm=5.1, kvp=60., added_cu_mm=0., field_side_cm=side)
    np.testing.assert_allclose(pair.k_bs, 0.2*lo.k_bs + 0.8*hi.k_bs)
    np.testing.assert_allclose(pair.k_f, 0.2*lo.k_f + 0.8*hi.k_f)
    audit = pair.audit["aluminium_interpolation"]
    assert audit["families_mm"] == [4.1, 5.1] and np.allclose(audit["weights"], [0.2, 0.8])
    with pytest.raises(InputError, match="no tabulado"):    # 3.1 + 2.5 = 5.6 mm: beyond the last family
        provider(aluminium_interpolation="linear_between_families").evaluate_event(event(kvp=60., cu=0., al=2.5,
                                                                                       width=300., height=300.), 270.)
    with pytest.raises(InputError, match="aluminium_interpolation"):
        provider(aluminium_interpolation="nearest")
