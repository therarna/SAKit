import pandas as pd
import gzip
import argparse

class VCFParser:
    def __init__(self, vcf_file, vcf_type="vardict", is_annotated=False):
        self.vcf_file = vcf_file
        self.vcf_type = vcf_type
        self.is_annotated = is_annotated
        self.df = None

    def parse_ann(self, ann_str):
        ann_list = []
        for ann in ann_str.split(","):
            ann_fields = ann.split("|")
            if len(ann_fields) >= 12:
                ann_dict = {
                    "Gene_Name": ann_fields[3],
                    "Feature_ID": ann_fields[6],
                    "Transcript_BioType": ann_fields[7],
                    "Annotation": ann_fields[1],
                    "HGVS.c": ann_fields[9],
                    "HGVS.p": ann_fields[10]
                }
                ann_list.append(ann_dict)
        return ann_list

    def parse_info(self, info_str):
        info_dict = {}
        for item in info_str.split(";"):
            key, value = item.split("=", 1)
            info_dict[key] = value
        return info_dict

    def parse_sample(self, sample_str, format_fields):
        sample_values = sample_str.split(":")
        sample_dict = dict(zip(format_fields, sample_values))
        return sample_dict

    def _read_vcf_lines(self):
        with gzip.open(self.vcf_file, "rt") if self.vcf_file.endswith(".gz") else open(self.vcf_file, "r") as file:
            lines = file.readlines()
        data_lines = []
        for line in lines:
            if not line.startswith("#"):
                data_lines.append(line.strip().split("\t"))
        return data_lines, lines

    def _get_header(self, lines):
        header = [h.strip("#") for h in lines if h.startswith("#CHROM")][0].strip().split("\t")
        return header

    def _get_num_samples(self, header):
        sample_names = header[9:]
        num_samples = len(sample_names)
        return num_samples

    def _get_new_cols(self, num_samples):
        cols_base = ["CHROM", "POS", "REF", "ALT", "FILTER"]
        cols_snvcall = ["STATUS", "TYPE", "MSI", "MSILEN"]
        cols_annotate = ["Gene_Name", "Feature_ID", "Transcript_BioType", "Annotation", "HGVS.c", "HGVS.p"]
        cols_paired = ["case_DP", "case_VD", "case_AF", "case_QUAL", "ctrl_DP", "ctrl_VD", "ctrl_AF", "ctrl_QUAL"]
        cols_single = ["DP", "VD", "AF", "QUAL"]

        if (self.is_annotated == True) and (num_samples == 2):
            new_cols = cols_base+cols_snvcall+cols_annotate+cols_paired
        elif (self.is_annotated == True) and (num_samples == 1):
            new_cols = cols_base+cols_snvcall+cols_annotate+cols_single
        elif  (self.is_annotated == False) and (num_samples == 2):
            new_cols = cols_base+cols_snvcall+cols_paired
        elif  (self.is_annotated == False) and (num_samples == 1):
            new_cols = cols_base+cols_snvcall+cols_single
        return new_cols
            
    def _process_vcf(self, df, num_samples):
        df["INFO"] = df["INFO"].apply(self.parse_info)

        df["STATUS"] = df["INFO"].apply(lambda x: x.get("STATUS", ""))
        df["TYPE"] = df["INFO"].apply(lambda x: x.get("TYPE", ""))
        df["MSI"] = df["INFO"].apply(lambda x: float(x.get("MSI", 0)))
        df["MSILEN"] = df["INFO"].apply(lambda x: float(x.get("MSILEN", 0)))
        format_fields = df["FORMAT"].iloc[0].split(":")
        
        if self.is_annotated == True:
            df["ANN"] = df["INFO"].apply(lambda x: x.get("ANN", ""))
            df["ANN"] = df["ANN"].apply(self.parse_ann)
            print(len(df))
            df = df.explode("ANN", ignore_index=True)
            print(len(df))
            df["Gene_Name"] = df["ANN"].apply(lambda x: x.get("Gene_Name"))
            df["Feature_ID"] = df["ANN"].apply(lambda x: x.get("Feature_ID"))
            df["Transcript_BioType"] = df["ANN"].apply(lambda x: x.get("Transcript_BioType"))
            df["Annotation"] = df["ANN"].apply(lambda x: x.get("Annotation"))
            df["HGVS.c"] = df["ANN"].apply(lambda x: x.get("HGVS.c"))
            df["HGVS.p"] = df["ANN"].apply(lambda x: x.get("HGVS.p"))

        if num_samples == 2:
            df = self._process_paired_samples(df, format_fields)
        elif num_samples == 1:
            df = self._process_single_sample(df, format_fields)
        
        new_cols = self._get_new_cols(num_samples)
        df = df[new_cols]
        
        return df

    def _process_paired_samples(self, df, format_fields):
        sample_names = df.columns[9:]
        df["case"] = df.apply(lambda row: self.parse_sample(row[sample_names[0]], format_fields), axis=1)
        df["ctrl"] = df.apply(lambda row: self.parse_sample(row[sample_names[1]], format_fields), axis=1)
        df["case_DP"] = df["case"].apply(lambda x: int(x.get("DP", 0)))
        df["case_VD"] = df["case"].apply(lambda x: int(x.get("VD", 0)))
        df["case_AF"] = df["case"].apply(lambda x: float(x.get("AF", 0)))
        df["case_QUAL"] = df["case"].apply(lambda x: float(x.get("QUAL", 0)))
        df["ctrl_DP"] = df["ctrl"].apply(lambda x: int(x.get("DP", 0)))
        df["ctrl_VD"] = df["ctrl"].apply(lambda x: int(x.get("VD", 0)))
        df["ctrl_AF"] = df["ctrl"].apply(lambda x: float(x.get("AF", 0)))
        df["ctrl_QUAL"] = df["ctrl"].apply(lambda x: float(x.get("QUAL", 0)))

        return df

    def _process_single_sample(self, df, format_fields):
        sample_name = df.columns[9]
        df["sample"] = df.apply(lambda row: self.parse_sample(row[sample_name], format_fields), axis=1)
        df["DP"] = df["sample"].apply(lambda x: int(x.get("DP", 0)))
        df["VD"] = df["sample"].apply(lambda x: int(x.get("VD", 0)))
        df["AF"] = df["sample"].apply(lambda x: float(x.get("AF", 0)))
        df["QUAL"] = df["INFO"].apply(lambda x: x.get("QUAL"))

        return df

    def vcf2dataframe(self):
        data_lines, lines = self._read_vcf_lines()
        header = self._get_header(lines)
        num_samples = self._get_num_samples(header)

        # Creating the DataFrame
        df = pd.DataFrame(data_lines, columns=header)

        if  self.vcf_type == "vardict":
            df = self._process_vcf(df, num_samples)
            # For non-annotated VCF, keep only the ID columns
        elif self.vcf_type == "dbsnp":
            cols_base = ["CHROM", "POS", "REF", "ALT", "ID"]
            df = df[cols_base]
            df["CHROM"] = "chr"+df["CHROM"].astype(str)
            
        else:
            print("vcf must be  format of vardict or dbsnp ")

        df["CHROM_POS"] = df["CHROM"] + ":" + df["POS"].astype(str)
        df["REF_ALT"] = df["REF"] + "_" + df["ALT"]
        df = df.set_index(["CHROM_POS", "REF_ALT"])

        return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse VCF file to DataFrame')
    parser.add_argument('-v', '--vcf_file', required=True, help='Path to the input VCF file')
    parser.add_argument('-t', '--vcf_type', default='vardict', choices=['vardict', 'dbsnp'],
                        help='Type of the input VCF file (default: vardict)')
    parser.add_argument('-a', '--is_annotated', action='store_true', help='Whether the input VCF file is annotated')
    parser.add_argument('-o', '--output_path', default='.', help='Path to save the output file (default: current directory)')
    parser.add_argument('-p', '--prefix', default='output', help='Prefix for the output file name (default: output)')

    args = parser.parse_args()

    vcf_parser = VCFParser(args.vcf_file, args.vcf_type, args.is_annotated)
    df = vcf_parser.vcf2dataframe()

    output_file = f"{args.output_path}/{args.prefix}.csv"
    df.to_csv(output_file, index=True)

    print(f"Parsed VCF file saved to {output_file}")
