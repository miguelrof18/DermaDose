"""End-to-end: renamed XCAT and Slicer sidecar beside an equipment JSON."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from skindose_verifier.example_runner import main
from skindose_verifier.synthetic import make_synthetic_rdsr


class ExampleIntegrationTest(unittest.TestCase):
    def test_relocated_mrc_profile_keeps_spectral_dataset(self):
        root=Path(__file__).resolve().parents[1]
        source=root/'config'/'equipo_mrc200_conjunto_investigacion.json'
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            shutil.copyfile(source,folder/'equipo.json')
            shutil.copyfile(root/'skindose_verifier'/'data'/'xcat_wholebody_body.stl',folder/'otra malla.stl')
            make_synthetic_rdsr(folder/'informe.dcm')
            config=json.loads(source.read_text(encoding='utf-8'))
            self.assertEqual(main([str(folder),'--prepare-only']),0)
            effective=json.loads((folder/'resultados'/'configuracion_efectiva.json').read_text(encoding='utf-8'))
            self.assertEqual(effective['surface_factors']['inherent_al_mm'],2.5)
            self.assertEqual(effective['surface_factors']['expected_anode_angle_deg'],11)
            self.assertEqual(Path(effective['surface_factors']['dataset_directory']).name,'derived_mrc200_nominal')
            self.assertEqual(Path(effective['patient']['mesh_path']).name,'otra malla.stl')

    def test_full_calculation_with_local_renamed_mesh(self):
        data_dir=Path(__file__).resolve().parents[1]/'skindose_verifier'/'data'
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            config=json.loads((data_dir/'equipo_sintetico.json').read_text(encoding='utf-8'))
            config['patient'].update(type='mesh',mesh_path='missing.stl',mesh_unit='mm',
                mesh_axes_matrix=[[1,0,0],[0,0,1],[0,-1,0]],mesh_origin_mm=[0,-500,40])
            config_file=folder/'equipo.json'
            config_file.write_text(json.dumps(config),encoding='utf-8')
            shutil.copyfile(data_dir/'xcat_wholebody_body.stl',folder/'mi superficie.stl')
            (folder/'mi superficie.json').write_text(json.dumps({'coordinate_system':'RAS','length_unit':'mm'}))
            make_synthetic_rdsr(folder/'informe.dcm')
            before=config_file.read_bytes()
            self.assertEqual(main([str(folder)]),0)
            result=json.loads((folder/'resultados'/'resultados.json').read_text(encoding='utf-8'))
            self.assertGreater(result['psd_point_gy'],0)
            self.assertEqual(result['event_count'],2)
            self.assertEqual(config_file.read_bytes(),before)
            effective=json.loads((folder/'resultados'/'configuracion_efectiva.json').read_text(encoding='utf-8'))
            self.assertEqual(Path(effective['patient']['mesh_path']).name,'mi superficie.stl')
            self.assertTrue((folder/'resultados'/'mapa_3d_interactivo.html').is_file())
            self.assertTrue((folder/'resultados'/'colocacion_malla.html').is_file())
