import math
import unittest
from skindose_verifier.placement_summary import placement_rows, html_summary


class PlacementSummaryTests(unittest.TestCase):
    def patient(self):
        return dict(mesh_unit='mm', mesh_axes_matrix=[[1,0,0],[0,0,1],[0,-1,0]],
                    mesh_origin_mm=[0,-500,40], origin_on_table_mm=[10,0,-20])

    def test_column_directions_and_displacements(self):
        rows = dict(placement_rows(self.patient()))
        self.assertEqual(rows['Eje +X del STL'], 'izquierda')
        self.assertEqual(rows['Eje +Y del STL'], 'posterior')
        self.assertEqual(rows['Eje +Z del STL'], 'cabeza')
        self.assertIn('50 cm hacia pies', str(rows))
        self.assertIn('4 cm hacia anterior', str(rows))
        self.assertIn('2 cm hacia posterior', str(rows))

    def test_oblique_rotation_is_not_reported_as_axis_aligned(self):
        patient = self.patient()
        c = math.sqrt(.5)
        patient['mesh_axes_matrix'] = [[c,-c,0],[c,c,0],[0,0,1]]
        rows = dict(placement_rows(patient))
        self.assertEqual(rows['Eje +X del STL'], '0.707 hacia izquierda + 0.707 hacia cabeza')

    def test_html_escapes_mesh_names_and_labels_assumption(self):
        html = html_summary(self.patient(), {'mesh': '<script>malla.stl</script>'})
        self.assertNotIn('<script>', html)
        self.assertIn('&lt;script&gt;', html)
        self.assertIn('No confirma la posición real', html)


if __name__ == '__main__':
    unittest.main()
