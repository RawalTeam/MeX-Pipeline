import os


class NgsTeMapper:

    def __init__(self, chromosome, start, end, strand, family,
                 supporting_reads, typ):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.strand = strand
        self.family = family
        self.supporting_reads = supporting_reads
        self.typ = typ

    @staticmethod
    def get_inserts(raw_file):
        inserts = []
        with open(raw_file) as f:
            for line in f:
                line = line.replace(";", "\t").split("\t")
                chromosome = line[0].strip()
                start = line[1].strip()
                end = line[2].strip()
                strand = line[4].strip()
                family = line[5].strip()
                supporting_read = line[-2].strip()
                typ = line[-1].strip()
                clazz = NgsTeMapper(chromosome, start, end, strand, family, supporting_read, typ)
                inserts.append(clazz)
        return inserts


def make_vcf_files(inserts, output_dir):
    malformed_inserts = []
    perfect_inserts = []
    non_reference_inserts = []
    for insert in inserts:
        if insert.start > insert.end:
            malformed_inserts.append(insert)
        else:
            perfect_inserts.append(insert)
            if insert.typ == "non-reference":
                non_reference_inserts.append(insert)

    with open(os.path.join(output_dir, "complete.vcf"), "w") as fw:
        for insert in perfect_inserts:
            chromosome = insert.chromosome
            start = insert.start
            end = insert.end
            strand = insert.strand
            family = insert.family
            supporting_reads = insert.supporting_reads
            typ = insert.typ
            record = f"{chromosome}\t{start}\t{end}\tINS\t{strand}\t" \
                     f"FAMILY={family};SUPPORTING_READS={supporting_reads};TYPE={typ};\n"
            fw.write(record)

    with open(os.path.join(output_dir, "non_reference.vcf"), "w") as fw:
        for insert in non_reference_inserts:
            chromosome = insert.chromosome
            start = insert.start
            end = insert.end
            strand = insert.strand
            family = insert.family
            supporting_reads = insert.supporting_reads
            typ = insert.typ
            record = f"{chromosome}\t{start}\t{end}\tINS\t{strand}\t" \
                     f"FAMILY={family};SUPPORTING_READS={supporting_reads};TYPE={typ};\n"
            fw.write(record)

    if len(malformed_inserts) > 0:
        with open(os.path.join(output_dir, "malformed.vcf"), "w") as fw:
            for insert in malformed_inserts:
                chromosome = insert.chromosome
                start = insert.start
                end = insert.end
                strand = insert.strand
                family = insert.family
                supporting_reads = insert.supporting_reads
                typ = insert.typ
                record = f"{chromosome}\t{start}\t{end}\tINS\t{strand}\t" \
                         f"FAMILY={family};SUPPORTING_READS={supporting_reads};TYPE={typ};\n"
                fw.write(record)
