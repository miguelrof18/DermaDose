"""Verification-set support (Gañán Mora et al., Rev Fis Med 2024) and the changes it motivated."""
import math
from pathlib import Path

import numpy as np
import pytest

from skindose_verifier.config import Equipment
from skindose_verifier.domain import InputError, IrradiationEvent
from skindose_verifier.engine import calculate, normalize
from skindose_verifier.surface_factors import make_event_provider
from skindose_verifier import table_transmission as tt
from skindose_verifier import verification_rfm2024 as rfm

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


def test_collimated_area_reference_rescales_only_area_derived_fields():
    cfg = rfm.study_configuration("Canon", -150.)
    eq = Equipment.from_dict(cfg)
    area = normalize(event(width=200., height=200., sid=910., area_derived=True), eq)
    np.testing.assert_allclose([area.field_width_mm, area.field_height_mm], [182., 182.])
    assert area.quality["field_size"].endswith("at_fixed_distance")
    explicit = normalize(event(width=200., height=200., sid=910.), eq)
    assert explicit.field_width_mm == 200.
    # Canon case 9: 200.73 mm at 1 m gives the 11 x 11 cm field of Table 2 at the IRP (546 mm).
    assert math.isclose(200.73*546/1000, 109.6, abs_tol=0.1)


def test_delorenzo_terms_match_the_skindose2d_parametrisation():
    assert math.isclose(tt.kvp_filtration_term(80., 0.2), 1 + 0.1321*math.log(1.5892*80 + 0.1321)
                        + 0.05043*math.log(1.5892*0.2 + 0.05043))
    assert math.isclose(tt.transmission(80., 0.2, 273.), 0.8365, abs_tol=5e-4)
    assert math.isclose(tt.angular_term(30.), 1 - 4.9e-5*900 + 3.07e-4*30)
    assert tt.field_area_term(2000.) == tt.field_area_term(580.)
    assert tt.angular_term(85.) == tt.angular_term(80.)
    blocks = tt.delorenzo_factor_blocks()
    eq = Equipment.from_dict(rfm.study_configuration("Philips", -150., table_factors=blocks))
    value = eq.factors["table_normal"].evaluate(kvp=np.array([70.]), filtration_mm_cu=np.array([0.4]),
                                                field_area_cm2=np.array([364.]))
    np.testing.assert_allclose(value, tt.kvp_filtration_term(70., .4)*tt.field_area_term(364.), rtol=2e-4)


def test_reference_set_is_complete_and_consistent():
    assert sorted(rfm.REFERENCE) == list(range(1, 28))
    for case, (device, _, dw, (lo, hi), kf, sk2d, above) in rfm.REFERENCE.items():
        assert device == ("Canon" if case <= 11 else "Philips")
        assert lo < dw < hi and 1.0 < kf < 1.06 and sk2d > 0
        assert above == (110. if case >= 24 else 0.)
    assert rfm.verdict(6.4) == "sin discrepancia"
    assert rfm.verdict(-6.6).startswith("discrepancia al 68")
    assert rfm.verdict(13.1).startswith("discrepancia significativa")
    assert rfm.case_from_name(Path("SKIN DOSE, TEST_CASE_9_MeasDose=12.44mGy_FS(11.0,11.0)")) == 9
    with pytest.raises(InputError, match="no coincide"):
        rfm.case_from_name(Path("SKIN DOSE, TEST_CASE_9_MeasDose=15.17mGy"))


def test_study_model_equals_closed_form_at_the_field_centre():
    """Engine dose at the centre = Ka (d_IRP/d)^2 T k_BS k_f with CF = 1 (Philips base, PA)."""
    cfg = rfm.study_configuration("Philips", -150.)
    eq = Equipment.from_dict(cfg)
    ev = event()
    result = calculate([ev], eq, "skin_dose", detailed=False)
    pts = result.surface.points_mm
    centre = np.hypot(pts[:, 0], pts[:, 1]) <= 3.
    dose = float(np.mean(result.total_gy[centre]))
    n = normalize(ev, eq)
    area = n.field_width_mm*n.field_height_mm*(660/1181.4)**2/100
    pair = make_event_provider(cfg["surface_factors"]).evaluate_event(n, area)
    expected = 0.01*tt.transmission(70., .4, area)*float(pair.k_bs)*float(pair.k_f)
    assert math.isclose(dose, expected, rel_tol=2e-3)
