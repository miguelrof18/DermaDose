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
        self.root = Path(self.temp.name).resolve()
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
        data,audit = self.prepare()
        self.assertEqual(Path(data['patient']['mesh_path']).name,'nombre libre.STL')
        self.assertEqual(Path(data['surface_factors']['dataset_directory']),self.resources/'derived_mrc200_nominal')
        self.assertEqual(data['patient']['mesh_origin_mm'],[1,2,3])
        self.assertEqual(before,{p.name:p.read_bytes() for p in self.folder.iterdir()})
        self.assertEqual(audit['registration_status'], 'not_established')

    def test_no_orientation_confirmation_or_patient_file_required(self):
        self.write('paciente.json', {'patient': {'mesh_origin_mm': [99,99,99]}})
        data, audit = self.prepare(interactive=True, input_fn=lambda _: '2')
        self.assertEqual(data['patient']['mesh_origin_mm'], [1,2,3])
        self.assertEqual(audit['placement_source'], 'equipment_patient_assumed')

    def test_defaults_without_patient_block(self):
        self.config.pop('patient')
        self.write('equipo.json', self.config)
        data, audit = self.prepare()
        self.assertEqual(data['patient']['mesh_unit'], 'mm')
        self.assertEqual(audit['placement_source'], 'provisional_defaults_assumed')

    def test_interactive_xcat_needs_no_local_stl(self):
        (self.folder/'nombre libre.STL').unlink()
        data, audit = self.prepare(interactive=True, input_fn=lambda _: '1')
        self.assertEqual(Path(data['patient']['mesh_path']), self.resources/'xcat_wholebody_body.stl')
        self.assertEqual(audit['mesh_choice'], 'xcat')

    def test_legacy_phantom_does_not_override_selection(self):
        self.config['phantom'] = {'mesh_path': 'bad.stl'}
        self.write('equipo.json', self.config)
        data, _ = self.prepare()
        self.assertNotIn('phantom', data)

    def test_multiple_meshes_require_selection(self):
        (self.folder/'otra.stl').write_bytes(b'other')
        with self.assertRaisesRegex(ExampleError,'varios'):
            self.prepare()
        data,_=self.prepare(mesh='otra.stl',)
        self.assertEqual(Path(data['patient']['mesh_path']).name,'otra.stl')

    def test_multiple_equipment_files_not_first_alphabetically(self):
        self.write('otro_equipo.json',self.config)
        with self.assertRaisesRegex(ExampleError,'varios'):
            self.prepare()
        self.prepare(config='otro_equipo.json',)

    def test_xcat_explicit_choice(self):
        data,_=self.prepare(mesh='xcat',)
        self.assertEqual(Path(data['patient']['mesh_path']),self.resources/'xcat_wholebody_body.stl')

    def test_custom_missing_dataset_not_replaced(self):
        self.config['surface_factors']['dataset_directory']='mis_tablas'
        self.write('equipo.json',self.config)
        with self.assertRaisesRegex(ExampleError,'metadata'):
            self.prepare()

    def test_interactive_choice_and_placement(self):
        (self.folder/'otra.stl').write_bytes(b'other')
        answers=iter(['2','2'])
        data,_=self.prepare(interactive=True,input_fn=lambda _:next(answers))
        self.assertEqual(Path(data['patient']['mesh_path']).name,'otra.stl')


if __name__=='__main__':
    unittest.main()
