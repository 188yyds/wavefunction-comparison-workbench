import copy
import os
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(os.path.abspath(__file__)).parent.parent
sys.path.insert(0,str(ROOT))
import model

class ComparisonModelTests(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory()
        self.paths=[]
        for name in ['z 波函数.h5','a 波函数.hdf5','第三个.h5']:
            path=Path(self.directory.name)/name;path.write_bytes(b'fixture')
            self.paths.append(str(path))
        self.project=model.new_project();self.group=self.project['groups'][0]
    def tearDown(self):self.directory.cleanup()
    def test_decimal_grid_is_exact(self):
        self.assertEqual(model.grid('.1','.3','.1'),['0.1','0.2','0.3'])
        with self.assertRaises(ValueError):model.grid('0','1','.3')
    def test_import_preserves_passed_order(self):
        model.assign_files(self.project,self.group,self.paths)
        self.assertEqual([self.group['files'][u] for u in ['1','1.5','2']],self.paths)
    def test_overflow_preserves_existing_bindings(self):
        self.group['files']['1.5']=self.paths[1];before=copy.deepcopy(self.group)
        with self.assertRaises(ValueError):model.assign_files(self.project,self.group,self.paths)
        self.assertEqual(self.group,before)
    def test_invalid_input_rejects_entire_batch(self):
        before=copy.deepcopy(self.group)
        with self.assertRaises(ValueError):model.assign_files(self.project,self.group,[self.paths[0],os.path.join(self.directory.name,'missing.h5')])
        self.assertEqual(self.group,before)
    def test_changed_grid_keeps_matching_u(self):
        model.assign_files(self.project,self.group,self.paths)
        model.update_grid(self.project,['1.5','2','.5'])
        self.assertEqual(self.group['files'],{'1.5':self.paths[1],'2':self.paths[2]})
        self.assertEqual(self.group['pending'],[self.paths[0]])
    def test_example_has_no_private_file_bindings(self):
        example=model.validate_project(model.example_project())
        self.assertEqual([g['name'] for g in example['groups']],['MTBL','DMRG'])
        for group in example['groups']:
            self.assertEqual(group['files'],{});self.assertEqual(group['pending'],[])
    def test_physics_modules_are_included(self):
        for name in ['backend.jl','core/体区间累积量核心.jl','core/scaling物理量定义.jl']:
            self.assertTrue((ROOT/name).is_file(),name)

if __name__=='__main__':unittest.main()
