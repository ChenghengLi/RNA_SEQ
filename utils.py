import pysam
from collections import defaultdict

class Isoform:
    def __init__(self, name, sequence, id):
        self.name = name
        self.sequence = sequence
        self.id = id
    
    def __str__(self):
        return self.name
    
    def __len__(self):
        return len(self.sequence)
    
    def get_sequence(self):
        return self.sequence

class Read:
    def __init__(self, query_name, query_sequence):
        self.query_name = query_name
        self.query_sequence = query_sequence
        self.isoforms = []

    def __str__(self):
        return self.query_name
    
    def __len__(self):
        return len(self.query_sequence)

class DataLoader:
    def __init__(self, bam_file_path, fasta_file_path):
        self.bam_file_path = bam_file_path
        self.fasta_file_path = fasta_file_path
        self.isoforms = self.load_isoforms(fasta_file_path)
        self.reads = {}

    def load_isoforms(self, path):
        isoforms = {}
        i = 0
        with pysam.FastaFile(path) as fasta:
            for seq_name in fasta.references:
                sequence = fasta.fetch(seq_name)
                isoforms[seq_name] = Isoform(seq_name, sequence, i)
                i += 1
        return isoforms

    def _map_reads(self):
        with pysam.AlignmentFile(self.bam_file_path, "rb") as bamfile:
            for read in bamfile:
                if not read.is_unmapped:
                    transcript_id = bamfile.get_reference_name(read.reference_id)
                    if transcript_id in self.isoforms:
                        if read.query_name not in self.reads:
                            read_instance = Read(
                                read.query_name,
                                read.query_sequence
                            )
                            self.reads[read.query_name] = read_instance
                        self.reads[read.query_name].isoforms.append(self.isoforms[transcript_id])

    def get_mapped_reads(self):
        self._map_reads()
        return list(self.reads.values())
