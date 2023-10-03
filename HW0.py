from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Function to find ORFs in a sequence
def find_orfs(sequence, min_length=100):
    orfs = []
    for strand, nuc in [(+1, sequence), (-1, sequence.reverse_complement())]:
        for frame in range(3):
            length = len(nuc)
            for start in range(frame, length - 2, 3):
                if str(nuc[start:start+3]) in ["ATG", "GTG", "TTG"]:
                    for end in range(start + 3, length - 2, 3):
                        if str(nuc[end:end+3]) in ["TAA", "TAG", "TGA"]:
                            if end - start >= min_length:
                                orf_start = start if strand == 1 else length - end
                                orf_end = end if strand == 1 else length - start
                                orfs.append((orf_start, orf_end))
                            break
    return orfs

# Read the FASTA file
fasta_file = "your_input.fasta"
sequences = list(SeqIO.parse(fasta_file, "fasta"))

# Initialize BED output
bed_records = []

# Iterate through sequences and find ORFs
for seq_record in sequences:
    orfs = find_orfs(seq_record.seq)
    for orf_start, orf_end in orfs:
        # Convert ORF coordinates to BED format (0-based)
        bed_start = orf_start - 1
        bed_end = orf_end
        # Create a BED record
        bed_record = SeqFeature(FeatureLocation(bed_start, bed_end), type="CDS")
        # Append to the BED records list
        bed_records.append(bed_record)

# Write the BED file
bed_file = "output.bed"
with open(bed_file, "w") as output:
    for bed_record in bed_records:
        output.write(f"{bed_record.location.start} {bed_record.location.end} {bed_record.type}\n")

print(f"BED file generated: {bed_file}")
