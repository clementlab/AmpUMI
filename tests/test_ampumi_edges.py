import unittest
import subprocess
import os
import tempfile
import sys

# Get the path to AmpUMI.py
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
AMPUMI_SCRIPT = os.path.join(SCRIPT_DIR, '..', 'AmpUMI', 'AmpUMI.py')

def write_fastq(path, reads):
    with open(path, 'w') as f:
        for r in reads:
            f.write(f"{r['id']}\n{r['seq']}\n+\n{r['qual']}\n")

class TestAmpUMISimulated(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.r1 = os.path.join(self.temp_dir.name, "r1.fq")
        self.r2 = os.path.join(self.temp_dir.name, "r2.fq")
        self.out1 = os.path.join(self.temp_dir.name, "out1.fq")
        self.out2 = os.path.join(self.temp_dir.name, "out2.fq")
        
    def tearDown(self):
        self.temp_dir.cleanup()

    def run_ampumi(self, *args):
        cmd = [sys.executable, AMPUMI_SCRIPT, 'Process'] + list(args)
        return subprocess.run(cmd, capture_output=True, text=True)

    def test_empty_fastq(self):
        write_fastq(self.r1, [])
        res = self.run_ampumi('--fastq', self.r1, '--fastq_out', self.out1, '--umi_regex', '^IIIIIIII')
        self.assertNotEqual(res.returncode, 0)
        self.assertIn("UMI command failed. Got no reads", res.stderr)

    def test_single_end_dedup(self):
        # Two identical UMIs + Seq -> Should dedup to 1
        reads = [
            {'id':'@Read1', 'seq':'AAAAAAAAACGTACGT', 'qual':'IIIIIIIIIIIIIIII'},
            {'id':'@Read2', 'seq':'AAAAAAAAACGTACGT', 'qual':'JJJJJJJJJJJJJJJJ'}, # Better qual
            {'id':'@Read3', 'seq':'TTTTTTTTACGTACGT', 'qual':'IIIIIIIIIIIIIIII'}  # Different UMI
        ]
        write_fastq(self.r1, reads)
        res = self.run_ampumi('--fastq', self.r1, '--fastq_out', self.out1, '--umi_regex', '^IIIIIIII')
        self.assertEqual(res.returncode, 0)
        
        with open(self.out1) as f:
            out_lines = f.readlines()
        self.assertEqual(len(out_lines), 8) # 2 reads output
        # Read2 should win over Read 1
        self.assertIn("@Read2\n", out_lines)
        self.assertNotIn("@Read1\n", out_lines)

    def test_umi_collision(self):
        # Same UMI (AAAAAAAA), two different sequences
        # Seq1 has 2 reads, Seq2 has 1 read -> Seq1 is majority, Seq2 is collision
        reads = [
            {'id':'@Read1', 'seq':'AAAAAAAAACGTACGT', 'qual':'IIIIIIIIIIIIIIII'},
            {'id':'@Read2', 'seq':'AAAAAAAAACGTACGT', 'qual':'IIIIIIIIIIIIIIII'},
            {'id':'@Read3', 'seq':'AAAAAAAATGCATGCA', 'qual':'IIIIIIIIIIIIIIII'} 
        ]
        write_fastq(self.r1, reads)
        res = self.run_ampumi('--fastq', self.r1, '--fastq_out', self.out1, '--umi_regex', '^IIIIIIII', '--write_alleles_with_multiple_UMIs')
        self.assertEqual(res.returncode, 0)
        self.assertIn("Observed 1 UMI collisions (1 reads)", res.stdout)
        
        # Check collision output file
        multi_out = self.out1 + ".AmpUMI.multi.out"
        with open(multi_out) as f:
            lines = f.readlines()
        self.assertTrue(any('k\t2\tAAAAAAAA\tACGTACGT' in l for l in lines))
        self.assertTrue(any('d\t1\tAAAAAAAA\tTGCATGCA' in l for l in lines))

    def test_min_umi_to_keep(self):
        # Read1 appears once, Read2 appears 3 times. If min_umi_to_keep=2, only Read2's UMI survives.
        reads = [
            {'id':'@Read1', 'seq':'AAAAAAAAACGTACGT', 'qual':'IIIIIIIIIIIIIIII'},
            {'id':'@Read2', 'seq':'TTTTTTTTACGTACGT', 'qual':'IIIIIIIIIIIIIIII'},
            {'id':'@Read3', 'seq':'TTTTTTTTACGTACGT', 'qual':'IIIIIIIIIIIIIIII'},
            {'id':'@Read4', 'seq':'TTTTTTTTACGTACGT', 'qual':'IIIIIIIIIIIIIIII'}
        ]
        write_fastq(self.r1, reads)
        res = self.run_ampumi('--fastq', self.r1, '--fastq_out', self.out1, '--umi_regex', '^IIIIIIII', '--min_umi_to_keep', '2')
        self.assertEqual(res.returncode, 0)
        self.assertIn("Observed 1 UMIs (1 reads) with too few reads", res.stdout)
        
        with open(self.out1) as f:
            out_lines = f.read()
        self.assertNotIn("@Read1", out_lines)
        self.assertIn("@Read2", out_lines)

    def test_paired_end(self):
        # Paired end, dedup relies on exact match of both R1 and R2
        r1_reads = [
            {'id':'@R1_1', 'seq':'AAAAAAAAACGTACGT', 'qual':'IIIIIIIIIIIIIIII'},
            {'id':'@R1_2', 'seq':'AAAAAAAAACGTACGT', 'qual':'JJJJJJJJJJJJJJJJ'}, # Better pair
            {'id':'@R1_3', 'seq':'AAAAAAAAACGTACGT', 'qual':'IIIIIIIIIIIIIIII'}  # R2 differs -> Collision
        ]
        r2_reads = [
            {'id':'@R1_1', 'seq':'GATCGATC', 'qual':'IIIIIIII'},
            {'id':'@R1_2', 'seq':'GATCGATC', 'qual':'JJJJJJJJ'},
            {'id':'@R1_3', 'seq':'CGATCGAT', 'qual':'IIIIIIII'}
        ]
        write_fastq(self.r1, r1_reads)
        write_fastq(self.r2, r2_reads)
        res = self.run_ampumi('--fastq', self.r1, '--fastq2', self.r2, '--fastq_out', self.out1, '--fastq_out2', self.out2, '--umi_regex', '^IIIIIIII')
        self.assertEqual(res.returncode, 0)
        
        with open(self.out1) as f1, open(self.out2) as f2:
            out1_lines = f1.read()
            out2_lines = f2.read()
            
        self.assertIn("@R1_2", out1_lines)  # Pair 2 wins over Pair 1
        self.assertNotIn("@R1_1", out1_lines)
        # Pair 3 is a collision to pair 1/2 UMI (1 majority pair, 1 minor)
        self.assertNotIn("@R1_3", out1_lines)
        self.assertIn("Observed 1 UMI collisions (1 reads)", res.stdout)
        
    def test_mean_quality(self):
        # Even though these would map to different clusters if sequence varies,
        # we can verify that the average vs sum executes without crashing
        reads = [
            {'id':'@Read1', 'seq':'AAAAAAAAACGTACGT', 'qual':'IIIIIIIIIIIIIIII'}
        ]
        write_fastq(self.r1, reads)
        res = self.run_ampumi('--fastq', self.r1, '--fastq_out', self.out1, '--umi_regex', '^IIIIIIII', '--use_sum_quality')
        self.assertEqual(res.returncode, 0)
        self.assertIn("Printed 1 deduplicated sequences", res.stdout)
        
    def test_truncate_length_varying_reads(self):
        # Two reads with same first 10 bp, but one is shorter due to trimming
        reads = [
            {'id':'@Read1', 'seq':'AAAAAAAAACGTACGT', 'qual':'IIIIIIIIIIIIIIII'}, # 16bp seq
            {'id':'@Read2', 'seq':'AAAAAAAAACGTAC', 'qual':'IIIIIIIIIIIIII'}     # 14bp seq
        ]
        write_fastq(self.r1, reads)
        # Without truncate_length, exact match fails, leading to a UMI collision (1 kept, 1 dropped = 1 printed)
        res1 = self.run_ampumi('--fastq', self.r1, '--fastq_out', self.out1, '--umi_regex', '^IIIIIIII')
        self.assertIn("WARNING: Reads of varying lengths were detected", res1.stdout)
        self.assertIn("Observed 1 UMI collisions", res1.stdout)
        self.assertIn("Printed 1 deduplicated sequences", res1.stdout)
        
        # With truncate_length=6, the 8bp and 6bp trimmed sequences match exactly and deduplicate into 1 common read
        res2 = self.run_ampumi('--fastq', self.r1, '--fastq_out', self.out1, '--umi_regex', '^IIIIIIII', '--truncate_length', '6')
        self.assertNotIn("WARNING: Reads of varying lengths were detected", res2.stdout)
        self.assertIn("Observed 0 UMI collisions", res2.stdout)
        self.assertIn("Printed 1 deduplicated sequences", res2.stdout)

    def test_truncate_length_paired_end(self):
        # Paired reads: R1 is exact, R2 is differently trimmed
        r1_reads = [
            {'id':'@R1_1', 'seq':'AAAAAAAAACGTACGT', 'qual':'IIIIIIIIIIIIIIII'},
            {'id':'@R1_2', 'seq':'AAAAAAAAACGTACGT', 'qual':'IIIIIIIIIIIIIIII'}
        ]
        r2_reads = [
            {'id':'@R1_1', 'seq':'GATCGATCGATC', 'qual':'IIIIIIIIIIII'}, # 12bp R2
            {'id':'@R1_2', 'seq':'GATCGATCGA',   'qual':'IIIIIIIIII'}   # 10bp R2
        ]
        write_fastq(self.r1, r1_reads)
        write_fastq(self.r2, r2_reads)
        
        # Without truncate, sequences differ, causing a collision
        res1 = self.run_ampumi('--fastq', self.r1, '--fastq2', self.r2, '--fastq_out', self.out1, '--fastq_out2', self.out2, '--umi_regex', '^IIIIIIII')
        self.assertIn("Observed 1 UMI collisions", res1.stdout)
        self.assertIn("Printed 1 deduplicated sequences", res1.stdout)
        
        # With truncate=10, they match exactly, preventing the collision
        res2 = self.run_ampumi('--fastq', self.r1, '--fastq2', self.r2, '--fastq_out', self.out1, '--fastq_out2', self.out2, '--umi_regex', '^IIIIIIII', '--truncate_length', '10')
        self.assertIn("Observed 0 UMI collisions", res2.stdout)
        self.assertIn("Printed 1 deduplicated sequences", res2.stdout)

if __name__ == '__main__':
    unittest.main()
