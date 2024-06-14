import csv
from Bio import SeqIO
import argparse
import re

def convert_amino_acid(sequence, to_three_letter=False):
    """
    Convert between single-letter and three-letter amino acid abbreviations.
    Args:
        sequence (str): The amino acid sequence to be converted.
        to_three_letter (bool, optional): If True, convert from single-letter to three-letter abbreviations.
                                          If False, convert from three-letter to single-letter abbreviations.
    Returns:
        str: The converted amino acid sequence.
    """
    # Dictionary mapping single-letter to three-letter amino acid abbreviations
    single_to_three = {
        'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
        'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
        'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
        'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
    }
    # Dictionary mapping three-letter to single-letter amino acid abbreviations
    three_to_single = {v: k for k, v in single_to_three.items()}
    if to_three_letter:
        return ''.join(single_to_three.get(aa, aa) for aa in sequence)
    else:
        converted = []
        for i in range(0, len(sequence), 3):
            three_letter = sequence[i:i+3]
            single_letter = three_to_single.get(three_letter, three_letter)
            converted.append(single_letter)
        return ''.join(converted)

def convert_hgvs_to_non_standard(hgvs):
    # 提取变异类型和位置
    match = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvs)
    
    if match:
        # 解析提取的信息
        ref_aa = match.group(1)
        position = match.group(2)
        alt_aa = match.group(3)
        
        # 将三字母代码转换为单字母代码
        ref_aa_single = convert_amino_acid(ref_aa, to_three_letter=False)
        alt_aa_single = convert_amino_acid(alt_aa, to_three_letter=False)
        
        # 生成非标准命名
        non_standard = f"p.{ref_aa_single}{position}{alt_aa_single}"
        return non_standard
    else:
        # 处理其它类型的变异描述
        ref_aa = re.findall(r'([A-Z][a-z]{2})', hgvs)
        ref_aa_single = [convert_amino_acid(aa, to_three_letter=False) for aa in ref_aa]
        non_standard = hgvs
        for aa, aa_single in zip(ref_aa, ref_aa_single):
            non_standard = non_standard.replace(aa, aa_single)
        return non_standard

def get_frameshift_seq(sequence, position):
    """
    Get the frameshift sequence from a given position until the first '*' character.
    Args:
        sequence (str): The input amino acid sequence.
        position (int): The starting position to search from.
    Returns:
        str: The subsequence from the given position until the first '*' character.
             If the position is out of range, an empty string is returned.
    """
    if position < 0 or position >= len(sequence):
        return ""
    for i in range(position, len(sequence)):
        if sequence[i] == '*':
            return sequence[position:i]
    return sequence[position:]

class MutationPeptideGenerator:
    def __init__(self, csv_file, fasta_file, out_path, prefix, flank_length=15, encoding='utf-8'):
        self.csv_file = csv_file
        self.fasta_file = fasta_file
        self.out_path = out_path
        self.prefix = prefix
        self.flank_length = flank_length
        self.encoding = encoding
        self.mutations = []
        self.fasta_seqs = {}

    def read_csv(self):
        with open(self.csv_file, 'r', encoding=self.encoding) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if  row['CHECK'] == 'Y':
                    self.mutations.append({
                        'transcript_id': row['Feature_ID'],
                        'hgvs_p': row['HGVS.p'],
                        'variant_type': row['Annotation'],
                        'gene_name': row['Gene_Name']
                    })

    def read_fasta(self):
        ref_seq = None
        for record in SeqIO.parse(self.fasta_file, 'fasta'):
            if 'Ref' == record.description.split()[1]:
                ref_seq = str(record.seq)
                
            else:
                stranscript_id = str(record.id).split()[0]
                hgvs_p = record.description.split("HGVS.p:")[1]
                variant_id = "_".join([stranscript_id,hgvs_p])
                self.fasta_seqs[variant_id] = {
                    'ref': ref_seq,
                    'variant': str(record.seq)
                }
                ref_seq = None
                
    def dedup_var_peptides(self, ref_peptides, var_peptides):
        ref_peptides_dedup = []
        var_peptides_dedup = []
        var_pep_uniq = []
        for ref_peptide, var_peptide in zip(ref_peptides, var_peptides):
            var_pep = var_peptide.split("\n")[1]
            if var_pep not in  var_pep_uniq:
                ref_peptides_dedup.append(ref_peptide)
                var_peptides_dedup.append(var_peptide)
                var_pep_uniq.append(var_pep)
        return ref_peptides_dedup, var_peptides_dedup

        
    def get_var_peptides(self):
        ref_peptides = []
        var_peptides = []
        for mutation in self.mutations:
            transcript_id = mutation['transcript_id']
            hgvs_p = mutation['hgvs_p']
            gene_name = mutation['gene_name']
            variant_id =  "_".join([transcript_id, hgvs_p])
            variant_type = mutation['variant_type'].split("&")

            try:
                ref_seq = self.fasta_seqs[variant_id]['ref']
                var_seq = self.fasta_seqs[variant_id]['variant']
            except:
                print(f'{variant_id} not in fasta_seqs')
                
            ref_peptide,var_peptide = "",""
            
            if len(variant_type) == 1:
                var_type = variant_type[0]
                if var_type == "missense_variant": ## eg: p.Gly13Asp
                    
                    ref_aa = hgvs_p[2:5]
                    alt_aa = hgvs_p[-3:]
                    mut_pos = int(hgvs_p[5:].rstrip(alt_aa))
                    left_pos = mut_pos - self.flank_length -1
                    left_pos = left_pos if left_pos > 0 else 0
                    right_pos = mut_pos + self.flank_length
                    
                    ref_peptide = ref_seq[left_pos : right_pos]
                    var_peptide = var_seq[left_pos: mut_pos-1] + convert_amino_acid(alt_aa) + var_seq[mut_pos : right_pos]
                    
                elif var_type == "frameshift_variant":  ## eg: p.Ser878fs
                    ref_aa = hgvs_p[2:5]
                    mut_pos = int(hgvs_p[5:].rstrip("fs"))
                    left_pos = mut_pos - self.flank_length -1
                    left_pos = left_pos if left_pos > 0 else 0
                    right_pos = mut_pos + self.flank_length
                    ref_peptide = ref_seq[left_pos  : right_pos]
                    var_peptide = get_frameshift_seq(var_seq, left_pos)
                    
                elif var_type == "disruptive_inframe_insertion":  ## eg: p.Thr882_Glu906dup
                    left_aa = hgvs_p[2:5]
                    right_aa = hgvs_p.split("_")[1][:3]
                    left_pos = int(hgvs_p[5:].split("_")[0]) -1 
                    left_pos = left_pos if left_pos_a > 0 else 0
                    right_pos = int(hgvs_p.split("_")[1].rstrip("dup")[2:])
                    dup_peptide = var_seq[left_pos: right_pos+1]
                    if len(dup_peptide) >= self.flank_length:
                        var_peptide = dup_peptide[-15:]*2
                        ref_peptide = ref_seq[left_pos:(right_pos + self.flank_length)]
                    else:
                        extend_length = self.flank_length - len(dup_peptide)
                        left_pos_a = left_pos - extend_length
                        left_pos_a = left_pos_a if left_pos_a > 0 else 0
                        
                        var_peptide = var_seq[( left_pos-extend_length if (left_pos-extend_length) >0 else None):(right_pos+ len(dup_peptide)+extend_length)]
                        ref_peptide = ref_seq[left_pos : (right_pos+ self.flank_length)]
                        
                elif var_type == "conservative_inframe_deletion": ## eg1: p.Gln94_Gln95del eg2: p.Asp31del
                    if "_" in hgvs_p:
                        print(hgvs_p,"yesdf")
                        left_pos = int(hgvs_p[5:].split("_")[0]) - 1
                        right_pos = int(hgvs_p.split("_")[1].rstrip("del")[2:])
                        del_base_num = right_pos - left_pos + 1
                    else:
                        left_pos = int(hgvs_p[5:].rstrip("del")) -1 
                        right_pos = int(hgvs_p[5:].rstrip("del"))
                        del_base_num = 1
                        
                    left_pos_shift = left_pos - self.flank_length
                    left_pos_shift = left_pos_shift if left_pos_shift > 0 else 0
                    right_pos_shift = right_pos + self.flank_length
                    ref_peptide = ref_seq[left_pos_shift : right_pos_shift]
                    var_peptide = var_seq[left_pos_shift : (right_pos_shift - del_base_num)]
                elif var_type == "synonymous_variant":
                    pass
                elif var_type == "stop_gained":
                    pass
                    
            elif  len(variant_type) > 1:
                if "frameshift_variant" in variant_type:  ## eg: p.Gly51fs
                    ref_aa = hgvs_p[2:5]
                    mut_pos = int(hgvs_p[5:].rstrip("fs"))
                    left_pos = mut_pos - self.flank_length -1
                    left_pos = left_pos if left_pos > 0 else 0
                    right_pos = mut_pos + self.flank_length
                    ref_peptide = ref_seq[left_pos  : right_pos]
                    var_peptide = get_frameshift_seq(var_seq, left_pos)
                    
                elif "conservative_inframe_deletion" in variant_type: ## eg: p.Glu40_Asn51del
                    """
                    annotation: splice_acceptor_variant&splice_donor_variant&conservative_inframe_deletion&splice_region_variant&intron_variant
                    eg: p.Glu40_Asn51del
                    """
                    left_pos = int(hgvs_p[5:].split("_")[0])  
                    right_pos = int(hgvs_p.split("_")[1].rstrip("del")[2:])
                    ref_peptide = ref_seq[(left_pos - self.flank_length):(right_pos + self.flank_length)]
                    var_peptide =  var_seq[(left_pos - self.flank_length):(left_pos + self.flank_length)]

            if (ref_peptide != "") and (var_peptide != ""):
                non_standard_hgvs_p = convert_hgvs_to_non_standard(hgvs_p)
                ref_peptides.append(f'>Ref {transcript_id}\n{ref_peptide}')
                var_peptides.append(f'>Var {transcript_id}({gene_name}):{hgvs_p} {non_standard_hgvs_p} \n{var_peptide}')
            else:
                print(f'{hgvs_p} was filtered, please check it' )
        return ref_peptides, var_peptides

    def write_fasta(self, ref_peptides, var_peptides):
        with open(f'{self.out_path}/{self.prefix}.ref_variants.faa', 'w') as outfile:
            for ref_peptide, var_peptide  in zip(ref_peptides, var_peptides):
                outfile.write(ref_peptide + "\n")
                outfile.write(var_peptide + "\n")

    def run(self):
        self.read_csv()
        
        self.read_fasta()
        ref_peptides, var_peptides = self.get_var_peptides()
        ref_peptides, var_peptides = self.dedup_var_peptides(ref_peptides, var_peptides)
        self.write_fasta(ref_peptides, var_peptides)
        print(f'已将{len(ref_peptides)}条Ref多肽段写入{self.out_path}/{self.prefix}.ref.fasta文件。')
        print(f'已将{len(var_peptides)}条Variant多肽段写入{self.out_path}{self.prefix}.variants.fasta文件。')

if __name__ == '__main__':
    # 创建ArgumentParser对象
    parser = argparse.ArgumentParser(description='Generate mutation peptide sequences from CSV and FASTA files.')

    # 添加命令行参数
    parser.add_argument('-m','--mutation_csv', help='Path to the CSV file containing mutation information.')
    parser.add_argument('-f','--aa_fasta', help='Path to the FASTA file containing amino acid sequences.')
    parser.add_argument('-o','--outpath', help='Path to the output directory.')
    parser.add_argument('-p','--prefix', help='Prefix for output file names.')
    parser.add_argument('-l','--flank_length', type=int, default=15, help='Length of flanking sequences (default: 15).')
    parser.add_argument('-e','--encoding', default='utf-8', help='Encoding for reading files (default: utf-8).')

    # 解析命令行参数
    args = parser.parse_args()

    # 创建MutationPeptideGenerator实例并运行
    generator = MutationPeptideGenerator(args.mutation_csv, args.aa_fasta, args.outpath, args.prefix, args.flank_length, args.encoding)
    generator.run()