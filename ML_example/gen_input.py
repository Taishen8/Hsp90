import os, subprocess
import pandas
import numpy
import pickle


class Seq:
    def __init__(self, seq_id, sequence, describe=None) -> None:
        self.id = seq_id
        self.seq = sequence
        self.describe = describe
        self.pssm = None
        self.ss = None  # secondary structure
        self.ids_long = None  # intrinsically disordered score
        self.ids_short = None  # intrinsically disordered score
        self.ids_anchor = None  # intrinsically disordered score

    @staticmethod
    def read_fasta(filepath) -> str:
        """
        Read the fasta containing only one seqence
        """
        with open(filepath, "r") as file:
            return file.readlines()[1].strip()

    def to_fasta(self, append=False, todir="../data/dataprepare_log") -> None:
        to = os.path.join(todir, f"{self.id}.fasta")
        # default not append to exists fasta, only one sequence in one fasta
        if append == False and os.path.exists(to):
            return
        with open(to, "a") as file:
            file.write(f">{self.id}\n{self.seq}")

    def run_psiblast(self, rewrite=False, workdir="../data/dataprepare_log") -> None:
        outfilename = f"{os.path.join(workdir, self.id)}.uniref50.pssm"
        if rewrite == False and os.path.exists(outfilename):
            return
        cline = f"/opt/ncbi-blast-2.13.0+/bin/psiblast \
            -db /tmp/local_opt/BLAST_data/uniref50 \
            -num_iterations 3 \
            -evalue 0.001 \
            -save_pssm_after_last_round \
            -query {os.path.join(workdir, self.id)}.fasta \
            -out_ascii_pssm {outfilename} \
            -out {os.path.join(workdir, self.id)}.psiblast.uniref50.log"
        # run PSI-BLAST
        result = subprocess.run(
            cline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        # skip ignorable errors
        for line in result.stderr.decode().rstrip("\n").splitlines():
            if "composition-based score adjustment is not supported with PSSMs" in line:
                continue
            print(line)

    def run_scratch(
        self,
        rewrite=False,
        bindir="/tmp/local_opt/scratch-1d-1.2",
        workdir="../data/dataprepare_log",
        scratch_version="1.2",
    ) -> None:
        """
        SCRATCH will only use one thread for one sequence
        """
        if scratch_version == "2.0":
            cline = f"/tmp/local_opt/scratch-1d-2.0/bin/run_scratch1d_predictors.sh \
                --input_fasta {os.path.join(workdir, self.id)}.fasta \
                --output_prefix {os.path.join(workdir, self.id)}.scratch20 \
                --num_threads 1 \
                > {os.path.join(workdir, self.id)}.scratch12.log"
        elif scratch_version == "1.2":
            outfilename = f"{os.path.join(workdir, self.id)}.scratch12.ss"
            if rewrite == False and os.path.exists(outfilename):
                return
            binpath = os.path.join(bindir, "bin/run_SCRATCH-1D_predictors.sh")
            cline = f"bash {binpath} \
                {os.path.join(workdir, self.id)}.fasta \
                {os.path.join(workdir, self.id)}.scratch12 \
                1 \
                > {os.path.join(workdir, self.id)}.scratch12.log"
        else:
            print("scratch_version should be 1.2 or 2.0")
            return
        result = subprocess.run(
            cline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        print(result.stderr.decode().rstrip("\n"))

    def run_iupred2a(self, rewrite=False, workdir="../data/dataprepare_log"):
        if rewrite == False and os.path.exists(
            os.path.join(workdir, f"{self.id}.IUPred2A.long.out")
        ):
            return
        # IUPred2A can run in two mode: long or short
        # run in long mode
        cline = f"python3 /opt/iupred2a/iupred2a.py \
            -a {os.path.join(workdir, self.id)}.fasta \
            long \
            > {os.path.join(workdir, self.id)}.IUPred2A.long.out"
        result = subprocess.run(
            cline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        print(result.stderr.decode().rstrip("\n"))
        # run in short mode
        cline = f"python3 /opt/iupred2a/iupred2a.py \
            -a {os.path.join(workdir, self.id)}.fasta \
            short \
            > {os.path.join(workdir, self.id)}.IUPred2A.short.out"
        result = subprocess.run(
            cline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        print(result.stderr.decode().rstrip("\n"))

    # read the output of IUPred2A
    # Intrinsically disordered region (IDR)
    def read_idr(self, indir="../data/dataprepare_log"):
        path = os.path.join(indir, f"{self.id}.IUPred2A.long.out")
        df = pandas.read_csv(path, sep="\s+", comment="#", header=None)
        self.ids_long = df.values[:, 2].astype(float)

        path = os.path.join(indir, f"{self.id}.IUPred2A.short.out")
        df = pandas.read_csv(path, sep="\s+", comment="#", header=None)
        self.ids_short = df.values[:, 2].astype(float)

        self.ids_anchor = df.values[:, 3].astype(float)

    # read the output of SCRATCH
    # Secondary structure (SS)
    def read_ss(self, indir="../data/dataprepare_log"):
        """
        Read the SCRATCH output in SINGLE-sequence fasta form
        Actually the output can be a multi-sequence fasta
        So this method should only be used for the current case
        """
        path = os.path.join(indir, f"{self.id}.scratch12.ss")
        with open(path, "r") as file:
            ss_str = file.readlines()[1].strip()
            self.ss = numpy.array(list(ss_str))

    # read the output of psi-blast
    # PSSM
    def read_pssm(self, indir="../data/dataprepare_log"):
        path = os.path.join(indir, f"{self.id}.uniref50.pssm")
        df = pandas.read_csv(
            path, sep="\s+", skiprows=3, skipfooter=6, header=None, engine="python"
        )
        self.pssm = df.iloc[:, 2:22].T.values.astype(int)


class CAMPSeq(Seq):
    # consider non-standard residues
    amino_acid_dict = {
        "A": 1,
        "C": 2,
        "B": 3,
        "E": 4,
        "D": 5,
        "G": 6,
        "F": 7,
        "I": 8,
        "H": 9,
        "K": 10,
        "M": 11,
        "L": 12,
        "O": 13,
        "N": 14,
        "Q": 15,
        "P": 16,
        "S": 17,
        "R": 18,
        "U": 19,
        "T": 20,
        "W": 21,
        "V": 22,
        "Y": 23,
        "X": 24,
        "Z": 25,
    }
    # revise order, not necessary but used by CAMP
    ss_dict = {"H": 1, "C": 2, "E": 3}
    seq_ss_array = []
    for seq in amino_acid_dict.keys():
        for ss in ss_dict.keys():
            seq_ss_array.append(seq + ss)
    # seq_ss_array = [seq + ss for seq in amino_acid_dict.keys() for ss in ss_dict.keys()]
    seq_ss_dict = {seq_ss: index + 1 for index, seq_ss in enumerate(seq_ss_array)}
    physicochemical_dict = {
        "A": 1,
        "C": 3,
        "B": 7,
        "E": 5,
        "D": 5,
        "G": 2,
        "F": 1,
        "I": 1,
        "H": 6,
        "K": 6,
        "M": 1,
        "L": 1,
        "O": 7,
        "N": 4,
        "Q": 4,
        "P": 1,
        "S": 4,
        "R": 6,
        "U": 7,
        "T": 4,
        "W": 2,
        "V": 1,
        "Y": 4,
        "X": 7,
        "Z": 7,
    }

    def __init__(self, seq_id, sequence, describe=None) -> None:
        super().__init__(seq_id, sequence, describe)
        self.sigmoided_pssm = None
        self.mapped_seq = None
        self.seq_2 = None
        self.seq_ss = None  # concat seq and ss, eg, res A & ss H -> AH

    @staticmethod
    def sigmoid(x: numpy.ndarray) -> numpy.ndarray:
        return 1 / (1 + numpy.exp(-x))

    @staticmethod
    def padding_serial(serial: numpy.ndarray, max_len) -> numpy.ndarray:
        """
        padding the serials, eg. seq, ss, ids, pssm
        serial.shape should be 1d or 2d. if 1d, reshape to 2d first
        if 2d, [channels, seq_len]
        """
        # reshape 1d serial to 2d
        reshaped = False
        if serial.ndim == 1:
            serial = serial.reshape(1, -1)
            reshaped = True
        print(serial.shape)
        padding_array = numpy.zeros([serial.shape[0], max_len])
        padding_len = min(serial.shape[1], max_len)
        padding_array[:, :padding_len] = serial[:, :padding_len]
        if reshaped is True:
            padding_array = padding_array[0]
        return padding_array

    @staticmethod
    def array_map(serial: numpy.ndarray, map_dict):
        vectorized_fun = numpy.vectorize(map_dict.get)
        return vectorized_fun(serial)

    def seq_map(self):
        seq_array = numpy.array(list(self.seq))
        self.mapped_seq = self.array_map(seq_array, self.amino_acid_dict)

    def concat_map_seq_ss(self):
        seq_array = numpy.array(list(self.seq))
        seq_ss = numpy.array([seq + ss for seq, ss in zip(seq_array, self.ss)])
        self.seq_ss = self.array_map(seq_ss, self.seq_ss_dict)

    def seq_2_map(self):
        seq_array = numpy.array(list(self.seq))
        self.seq_2 = self.array_map(seq_array, self.physicochemical_dict)


def preparedata(seq_id, sequence=None, indir="data/dataprepare_log") -> CAMPSeq:
    """
    init a CAMPSeq
    """
    if sequence is None:
        fasta_path = os.path.join(indir, f"{seqid}.fasta")
        if not os.path.exists(fasta_path):
            print(f"ERROR: {seq_id}: no fasta or given sequence")
            return None
        sequence = CAMPSeq.read_fasta(fasta_path)

    myseq = CAMPSeq(seq_id, sequence)
    myseq.to_fasta(todir=indir)

    is_prot = len(sequence) > 50

    # PSSM
    if is_prot:
        print("running pssm")
        myseq.run_psiblast(workdir=indir)
        myseq.read_pssm(indir=indir)
        myseq.sigmoided_pssm = myseq.sigmoid(myseq.pssm)

    print("running scratch")
    # Secondary structure
    myseq.run_scratch(workdir=indir)
    myseq.read_ss(indir=indir)
    print("running iupred2a")
    # Intrinsically disordered region (IDR)
    myseq.run_iupred2a(workdir=indir)
    myseq.read_idr(indir=indir)
    print("mapping")
    # sequence encode
    myseq.seq_map()
    # secondary structure encode with sequence
    myseq.concat_map_seq_ss()
    # aminos physicochemical codes
    myseq.seq_2_map()

    if is_prot and len(myseq.seq) != myseq.pssm.shape[1]:
        print(f"len of seq {len(myseq.seq)} != PSSM {myseq.pssm.shape[1]}")

    return myseq


def dump_seq(seq: Seq, filepath):
    """dump class Seq to pickle
    using pickle or pandas
    """
    pass


def read_seq(filepath) -> Seq:
    pass


def padx(x, N):
    if len(x.shape) == 1:
        x = numpy.expand_dims(x, axis=1)
    padding_array = numpy.zeros([N, x.shape[1]])
    if x.shape[0] >= N:  # sequence is longer than N
        padding_array[:N, : x.shape[1]] = x[:N, :]
    else:
        padding_array[: x.shape[0], : x.shape[1]] = x
    return padding_array.flatten()


def gen_row(pep, prot):
    pep_seq = padx(pep.mapped_seq, 50)
    prot_seq = padx(prot.mapped_seq, 589)
    pep_ss = padx(pep.seq_ss, 50)
    prot_ss = padx(prot.seq_ss, 589)
    pep_2 = padx(pep.seq_2, 50)
    prot_2 = padx(prot.seq_2, 589)
    pep_dense = padx(
        numpy.stack([pep.ids_long, pep.ids_short, pep.ids_anchor]).transpose(1, 0), 50
    )

    prot_pssm = prot.sigmoided_pssm.transpose(1, 0)
    prot_iud = numpy.stack([prot.ids_long, prot.ids_short, prot.ids_anchor]).transpose(
        1, 0
    )
    try:
        prot_dense = padx(numpy.concatenate([prot_pssm, prot_iud], axis=1), 589)
    except:
        return None

    return numpy.concatenate(
        [pep_seq, prot_seq, pep_ss, prot_ss, pep_2, prot_2, pep_dense, prot_dense],
        axis=0,
    )


def update_gen(seqs, name):
    try:
        infodict = pickle.load(open(f"{name}s.p", "rb"))
    except:
        infodict = {}
    need = set()
    for seq in seqs:
        if seq not in infodict.keys():
            need.add(seq)

    for i, seq in enumerate(need):
        infodict[seq] = preparedata(f"{name}{i}", seq)

    pickle.dump(infodict, open(f"{name}s.p", "+wb"))
    return infodict


import shutil

if __name__ == "__main__":
    prots = open("prots.txt").read().split("\n")
    peps = open("peps.txt").read().split("\n")

    os.makedirs("./data/dataprepare_log", exist_ok=True)

    prot_seqs, pep_seqs = [], []
    prot_dict = update_gen(prots, "prot")
    pep_dict = update_gen(peps, "pep")

    for prot in prots:
        prot_seqs.append(prot_dict[prot])
    for pep in peps:
        pep_seqs.append(pep_dict[pep])

    shutil.rmtree("./data")

    rows = []

    for pep_seq in pep_seqs:
        for prot_seq in prot_seqs:
            rows.append(gen_row(pep_seq, prot_seq))

    with open("input.p", "+wb") as f:
        pickle.dump(numpy.stack(rows), f)
