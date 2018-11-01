
try:
  from string import maketrans
  COMPL = maketrans("ATGC","TACG")
except: # python 3
  COMPL = str.maketrans("ATGC","TACG")

def rc(seq):
  return seq.translate(COMPL)[::-1]

def read_fasta(fasta, split_at_space=False):
  data = open(fasta).read().strip().split('\n')
  reads = {}
  names = []
  name = None
  seq = ""
  for d in data:
    if len(d) == 0:
      continue
    if d[0] == '>':
      if name is not None:
        reads[name] = seq
      name = d[1:]
      if split_at_space and ' ' in name:
        name = name[:name.index(' ')]
      names.append(name)
      seq = ""
    else:
      seq += d
  reads[name] = seq
  return reads, names

def iter_fasta(fasta, split_at_space=False):
  fin = open(fasta, 'r')
  name = None
  seq = ""
  for d in fin:
    if len(d) == 0:
      continue
    if d[0] == '>':
      if name is not None:
        yield name, seq
      name = d[1:]
      if split_at_space and ' ' in name:
        name = name[:name.index(' ')]
      names.append(name)
      seq = ""
    else:
      seq += d
  fin.close()
  if name is not None:
    yield name, seq
  return

def write_fasta(fname, reads, names=None):
  fout = open(fname, 'w')
  if names is None:
    names = list(reads.keys())
  for r in range(len(names)):
    if not names[r] in reads:
      print("WARNING: Writing fasta -- %s missing, skipped" % names[r])
      continue
    fout.write(">%s\n" % names[r])
    seq = reads[names[r]]
    fout.write("%s\n" % '\n'.join([seq[i:min(len(seq),i+80)] for i in range(0, len(seq), 80)]))
  fout.close()

def read_fastq(fastq, split_at_space=False):
  data = open(fastq).read().strip().split('\n')
  reads = {}
  quals = {}
  names = []
  i = 0
  row = 0
  while i < len(data):
    if len(data[i].strip()) == 0:
      i += 1
      continue
    if row%4 == 0:
      name = data[i][1:]
      if split_at_space and ' ' in name:
        name = name[:name.index(' ')]
      names.append(name)
    if row%4 == 1:
      reads[name] = data[i]
    if row%4 == 3:
      quals[name] = data[i]
    i += 1
    row += 1
  return reads, quals, names

def iter_fastq(fastq, split_at_space=False):
  fin = open(fastq, 'r')
  while True:
    name = fin.readline().strip()
    if split_at_space and ' ' in name:
      name = name[:name.index(' ')]
    if name is None or len(name) == 0:
      return
    seq = fin.readline().strip()
    fin.readline() # qual header
    qual = fin.readline().strip()
    yield name, seq, qual

def write_fastq(fname, reads, quals, names=None):
  fout = open(fname, 'w')
  if names is None:
    names = list(reads.keys())
  for r in range(len(names)):
    if not (names[r] in reads) or not (names[r] in quals):
      print("WARNING: Writing fasta -- %s missing sequence or quality, skipped" % names[r])
      continue
    fout.write("@{}\n".format(names[r]))
    fout.write("{}\n".format(reads[names[r]]))
    fout.write("+\n")
    fout.write("{}\n".format(quals[names[r]]))
  fout.close()
