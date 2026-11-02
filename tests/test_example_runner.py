import json
from pathlib import Path
import tempfile
import unittest
from hashlib import sha256

from skindose_verifier.example_runner import prepare_example, ExampleError


class ExampleRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.folder = self.root/'ejemplo07'
        self.folder.mkdir()
        self.resources = self.root/'data'
        (self.resources/'derived_mrc200_nominal').mkdir(parents=True)
        (self.resources/'derived_mrc200_nominal'/'metadata.json').write_text('{}')
        (self.resources/'xcat_wholebody_body.stl').write_bytes(b'XCAT fixture')
        (self.folder/'informe.dcm').write_bytes(b'fixture')
        (self.folder/'nombre libre.STL').write_bytes(b'mesh fixture')
        self.config = dict(version='test', manufacturer='Philips', model='Azurion',
            geometry={}, factors={}, patient=dict(mesh_path='../invalid/xcat.stl', mesh_unit='mm',
            mesh_axes_matrix=[[1,0,0],[0,1,0],[0,0,1]],mesh_origin_mm=[1,2,3],origin_on_table_mm=[0,0,0]),
            surface_factors=dict(dataset_directory='../skindose_verifier/data/derived_mrc200_nominal'))
        self.write('equipo.json', self.config)
        self.write('nombre libre.json', dict(schema_version=4, coordinate_system='RAS',length_unit='mm'))

    def write(self, name, value):
        (self.folder/name).write_text(json.dumps(value),encoding='utf-8')

    def prepare(self, **kwargs):
        return prepare_example(self.folder,data_directory=self.resources,**kwargs)

    def test_local_stl_and_equipment_detection_and_dataset_repair(self):
        before = {p.name:p.read_bytes() for p in self.folder.iterdir()}
        data,audit = self.prepare(use_configured_placement=True)
        self.assertEqual(Path(data['patient']['mesh_path']).name,'nombre libre.STL')
        self.assertEqual(Path(data['surface_factors']['dataset_directory']),self.resources/'derived_mrc200_nominal')
        self.assertEqual(data['patient']['mesh_origin_mm'],[1,2,3])
        self.assertEqual(before,{p.name:p.read_bytes() for p in self.folder.iterdir()})
        self.assertEqual(audit['surface_metadata']['coordinate_system'],'RAS')

    def test_no_implicit_xcat_placement_for_custom_mesh(self):
        with self.assertRaisesRegex(ExampleError,'registro'):
            self.prepare()

    def test_separate_patient_overrides_xcat_and_legacy_phantom(self):
        self.config['phantom']={'name':'ICRP_AM','mesh_path':'bad.stl'}
        self.write('equipo.json',self.config)
        self.write('paciente.json',{'patient':dict(mesh_unit='mm',mesh_axes_matrix=[[1,0,0],[0,1,0],[0,0,1]],mesh_origin_mm=[20,30,40],origin_on_table_mm=[0,0,0])})
        data,audit=self.prepare()
        self.assertNotIn('phantom',data)
        self.assertEqual(data['patient']['mesh_origin_mm'],[20,30,40])
        self.assertTrue(audit['placement_source'].endswith('paciente.json'))

    def test_multiple_meshes_require_selection(self):
        (self.folder/'otra.stl').write_bytes(b'other')
        with self.assertRaisesRegex(ExampleError,'varios'):
            self.prepare(use_configured_placement=True)
        data,_=self.prepare(mesh='otra.stl',use_configured_placement=True)
        self.assertEqual(Path(data['patient']['mesh_path']).name,'otra.stl')

    def test_multiple_equipment_files_not_first_alphabetically(self):
        self.write('otro_equipo.json',self.config)
        with self.assertRaisesRegex(ExampleError,'varios'):
            self.prepare(use_configured_placement=True)
        self.prepare(config='otro_equipo.json',use_configured_placement=True)

    def test_xcat_explicit_choice(self):
        data,_=self.prepare(mesh='xcat',use_configured_placement=True)
        self.assertEqual(Path(data['patient']['mesh_path']),self.resources/'xcat_wholebody_body.stl')

    def test_hash_mismatch_is_rejected(self):
        self.write('nombre libre.json',{'mesh_sha256':'incorrect'})
        with self.assertRaisesRegex(ExampleError,'hash'):
            self.prepare(use_configured_placement=True)

    def test_unit_mismatch_is_rejected(self):
        self.write('nombre libre.json',{'length_unit':'cm'})
        with self.assertRaisesRegex(ExampleError,'unidades'):
            self.prepare(use_configured_placement=True)

    def test_custom_missing_dataset_not_replaced(self):
        self.config['surface_factors']['dataset_directory']='mis_tablas'
        self.write('equipo.json',self.config)
        with self.assertRaisesRegex(ExampleError,'metadata'):
            self.prepare(use_configured_placement=True)

    def test_interactive_choice_and_placement(self):
        (self.folder/'otra.stl').write_bytes(b'other')
        answers=iter(['2','s'])
        data,_=self.prepare(interactive=True,input_fn=lambda _:next(answers))
        self.assertEqual(Path(data['patient']['mesh_path']).name,'otra.stl')


if __name__=='__main__':
    unittest.main()
