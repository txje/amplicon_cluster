import argparse
try:
  from string import maketrans
  COMPL = maketrans("ATGC","TACG")
except: # python 3
  COMPL = str.maketrans("ATGC","TACG")

class Alignment:
  def __init__(self, query=None, target=None, alignment=None, m=None, x=None, i=None, d=None, acc=None):
     self.query = query
     self.target = target
     self.alignment = alignment
     self.m = m # num matches
     self.x = x # num mismatches
     self.i = i # num insertions
     self.d = d # num deletions
     self.acc = acc

  def accuracy(self):
    if self.acc is not None:
      pass
    elif self.m is not None and self.x is not None and self.i is not None and self.d is not None:
      self.acc = float(self.m) / (self.m + self.x + self.i + self.d)
    elif self.alignment is not None:
      self.acc = float(self.alignment.count('|')) / len(self.alignment)
    else:
      raise ValueError("Accuracy cannot be computed")
    return self.acc

  def filter(self, length_threshold, alignment_threshold, edge_dist_threshold):
    # filter out self-alignments
    if self.query.name == self.target.name:
      return False
    # filter by length
    if self.query.end - self.query.start < length_threshold:
      return False
    # filter by accuracy
    if self.accuracy() < alignment_threshold:
      return False
    # filter by query edge overlap
    if self.query.start > edge_dist_threshold and self.query.end < self.query.length - edge_dist_threshold:
      return False
    # filter by target edge overlap
    if self.target.start > edge_dist_threshold and self.target.end < self.target.length - edge_dist_threshold:
      return False
    return True

  def invert(self):
    self.query.invert()
    self.target.invert()
    if self.alignment is not None:
      self.alignment = self.alignment[::-1]

  def make_alignment_string(self):
    if self.alignment is None:
      self.alignment = ''.join('|' if self.target.alignment[i] == self.query.alignment[i] else '*' for i in range(len(self.query.alignment)))

  def get_errors(self):
    if self.m is None:
      self.m = 0
      self.x = 0
      self.i = 0
      self.d = 0
      for i in range(len(self.query.alignment)):
        if self.query.alignment[i] == '-':
          self.d += 1
        elif self.target.alignment[i] == '-':
          self.i += 1
        elif self.target.alignment[i] == self.query.alignment[i]:
          self.m += 1
        else:
          self.x += 1

  def norm(self):
    # normalize to target + strand
    # make locations of insertions and deletions consistent relative to the target strand
    # There may be (in fact, we are optimized for) adjacent insertions and deletions

    # convert to lists for modification
    if self.target.reverse:
      self.invert()

    query = [a for a in self.query.alignment]
    target = [a for a in self.target.alignment]

    # covert mismatches into adjacent ins+del
    for i in range(len(query)):
      if query[i] != '-' and target[i] != '-' and query != target:
        query.insert(i+1, '-')
        target.insert(i, '-')

    # move deletions AFTER match, if equivalent
    trailing = None
    for i in range(len(query)):
      if query[i] == '-':
        if trailing is None or target[i] != target[trailing]:
          trailing = i
      elif trailing is not None and target[i] == target[trailing]:
        if query[i] != '-':
          query[trailing] = query[i]
          query[i] = '-'
          trailing += 1
      else:
        trailing = None

    # move insertion BEFORE match, if equivalent
    trailing = None
    for i in range(len(target)-1, -1, -1):
      if target[i] == '-':
        if trailing is None or query[i] != query[trailing]:
          trailing = i
      elif trailing is not None and query[i] == query[trailing]:
        if target[i] != '-':
          target[trailing] = target[i]
          target[i] = '-'
          trailing -= 1
      else:
        trailing = None

    # move inserts AFTER deletions
    trailing = None
    for i in range(len(query)):
      if target[i] == '-':
        if trailing is None:
          trailing = i
      elif trailing is not None and query[i] == '-':
        query.insert(trailing, query.pop(i))
        target.insert(trailing, target.pop(i))
        trailing += 1
      else:
        trailing = None

    # convert back to strings
    self.query.alignment = ''.join(query)
    self.target.alignment = ''.join(target)

  def __str__(self):
    return "%s %i:%i (%s) <-> %s %i:%i (%s) (~%ibp x %.2f%%)" % (self.query.name, self.query.start, self.query.end, '-' if self.query.reverse else '+', self.target.name, self.target.start, self.target.end, '-' if self.target.reverse else '+', self.query.end - self.query.start, self.accuracy() * 100)

# 0   qName qSeqLength qStart qEnd qStrand
# 5   tName tSeqLength tStart tEnd tStrand
# 10  score numMatch numMismatch numIns numDel
# 15  mapQV qAlignedSeq matchPattern tAlignedSeq
  def m5line(self):
    self.make_alignment_string()
    self.get_errors()
    return "%s %i %i %i %s %s %i %i %i %s %i %i %i %i %i %i %s %s %s" % (
        self.query.name, self.query.length, self.query.start, self.query.end, '-' if self.query.reverse else '+',
        self.target.name, self.target.length, self.target.start, self.target.end, '-' if self.target.reverse else '+',
        0, self.m, self.x, self.i, self.d,
        0, self.query.alignment, self.alignment, self.target.alignment)


class AlignmentRead:
  # doubles as a read in a contig, where end and alignment are not used
  def __init__(self, name=None, length=None, start=None, end=None, reverse=None, alignment=None):
    self.name = name
    self.length = length
    self.start = start
    self.end = end
    self.reverse = reverse
    self.alignment = alignment

  def invert_positions(self):
    tmp = self.start
    self.start = self.length - self.end # this works, ends are exclusive
    self.end = self.length - tmp # and beginnings are inclusive

  def invert(self):
    self.reverse = not self.reverse
    self.invert_positions()
    if self.alignment is not None:
      self.alignment = self.alignment.translate(COMPL)[::-1]

  def __str__(self):
    return "{} {}:{} ({})".format(self.name, self.start, self.end, '-' if self.reverse else '+')

  def __repr__(self):
    return str(self)


class Contig:
  def __init__(self):
    self.reads = []
    self.read_map = {}

  def flip(self):
    for r in range(len(self.reads)):
      self.reads[r].reverse = not self.reads[r].reverse # switch strand
      self.reads[r].start = -1 * self.reads[r].start - self.reads[r].length # flip over start

  def add_read(self, read, start, length):
    contig_read = AlignmentRead(read.name, length, start, None, read.reverse)
    self.reads.append(contig_read)
    self.read_map[read.name] = contig_read

  def read(self, name):
    return self.read_map[name]

  def __str__(self):
    reads = sorted(self.reads, key = lambda a: a.start)
    xoffset = min(a.start for a in reads)
    scalar = 1000
    namepad = 150
    s = ""
    for r in self.reads:
      line = "%s (%i)" % (r.name, r.start)
      for i in range(0, (namepad - len(line))*scalar + r.start - xoffset, scalar):
        line += " "
      for i in range(r.start/scalar, (r.start + r.length)/scalar + 1):
        line += '-' if r.reverse else '+'
      s += line + "\n"
    return s


# m5 format:
# 0   qName qSeqLength qStart qEnd qStrand
# 5   tName tSeqLength tStart tEnd tStrand
# 10  score numMatch numMismatch numIns numDel
# 15  mapQV qAlignedSeq matchPattern tAlignedSeq

def iter_m5(m5_file, length_threshold=0, alignment_threshold=0, edge_dist_threshold=10000000000, by_query=False):
  if by_query:
    primary_index = 0
  else:
    primary_index = 5
  alignments = []
  template_name = None
  for line in open(m5_file):
    if line[0] == '#': # comment
      continue
    data = line.strip().split()
    if template_name == None:
      template_name = data[primary_index]
    if template_name != data[primary_index]:
      yield (template_name, alignments)
      template_name = data[primary_index]
      alignments = []

    # get and fix query name (blasr appends '/start_end' to aligned query names)
    query_name = data[0]
    if '/' in query_name:
      query_name = query_name[:query_name.rindex('/')]

    # all positions are relative to the forward strand
    query = AlignmentRead(query_name, int(data[1]), int(data[2]), int(data[3]), data[4] == '-', data[16])
    template = AlignmentRead(data[5], int(data[6]), int(data[7]), int(data[8]), data[9] == '-', data[18])
    aln = Alignment(query, template, data[17], int(data[11]), int(data[12]), int(data[13]), int(data[14]))
    if aln.filter(length_threshold, alignment_threshold, edge_dist_threshold):
      alignments.append(aln)

  yield (template_name, alignments)


# .m4 Interval Alignment format:

# 0   qName tName score pctSimilarity
# 4   qStrand qAlignStart qAlignEnd qLength
# 8   tStrand tAlignStart tAlignEnd tLength
# 12  mapQV
#
# qStrand should always be 0, tStrand may be 0 or 1
# q and t positions are strand-adjusted

def iter_m4(m4_file, length_threshold, alignment_threshold, edge_dist_threshold, by_query=False):
  if by_query:
    primary_index = 0 # query name
  else:
    primary_index = 1 # target name

  alignments = []
  template_name = None
  for line in open(m4_file):
    if line[0] == '#': # comment
      continue
    data = line.strip().split()
    if template_name == None:
      template_name = data[primary_index]
    if template_name != data[primary_index]:
      yield (template_name, alignments)
      template_name = data[primary_index]
      alignments = []

    # get and fix query name (blasr appends '/start_end' to aligned query names)
    query_name = data[0]
    if '/' in query_name:
      query_name = query_name[:query_name.rindex('/')]

    # by default, positions are relative to the aligned strand, so we'll change them
    # to match the forward sequence if necessary
    query = AlignmentRead(query_name, int(data[7]), int(data[5]), int(data[6]), data[4] == '1')
    if query.reverse:
      query.invert_positions()

    target = AlignmentRead(data[1], int(data[11]), int(data[9]), int(data[10]), data[8] == '1')
    if target.reverse:
      target.invert_positions()

    aln = Alignment(query, target)
    aln.acc = float(data[3]) / 100.0

    if aln.filter(length_threshold, alignment_threshold, edge_dist_threshold):
      alignments.append(aln)

  yield (template_name, alignments)


# MAF format:
#   0     1     2        3         4     5        6
# {a,s} qName qStart qAlnLength qStrand qLen qAlignedSeq
#
# s: sequence
# a: annotation (incl. scores)

# "target" is the first listed sequence

def iter_maf(maf_file, length_threshold=100, alignment_threshold=0.7, edge_dist_threshold=10000000000, by_query=False):
  alignments = []
  template_name = None
  lastread = None
  for line in open(maf_file):
    data = line.strip().split()
    '''
    if data[0] == 'a':
      for d in range(1, len(data)):
        k,v = data[d].split('=')
        if k == "score":
          score = float(v)
    '''
    if len(data) > 0 and data[0] == 's':
      # by default, positions are relative to the aligned strand, so we'll change them
      # to match the forward sequence if necessary
      read = AlignmentRead(data[1], int(data[5]), int(data[2]), int(data[2])+int(data[3]), data[4] == '-', data[6])
      if read.reverse:
        read.invert_positions()

      if ((by_query and lastread is not None) or (not by_query and lastread is None)) and data[1] != template_name:
        if template_name is not None:
          yield (template_name, alignments)
        template_name = data[1]
        alignments = []

      if lastread is None:
        lastread = read
      else:
        aln = Alignment(read, lastread)
        aln.make_alignment_string()
        if aln.filter(length_threshold, alignment_threshold, edge_dist_threshold):
          alignments.append(aln)
        lastread = None

  yield (template_name, alignments)


def maf2m5(maf_file, m5_file):
  fout = open(m5_file, 'w')
  lines = 0
  for template, alignments in iter_maf(maf_file, 0, 0, 5000000000): # no restrictions
    for aln in alignments:
      if lines > 0:
        fout.write('\n')
      fout.write(aln.m5line())
      lines += 1
  fout.close()


def iter_bam(bam_file, length_threshold=100, alignment_threshold=0.7, edge_dist_threshold=10000000000, by_query=False, ref_name=None, st=None, en=None, ref_seqs=None, query_seqs=None):
  import pysam
  bam = pysam.AlignmentFile(bam_file, 'rb')

  if by_query:
    alignments = {}

  for r in range(len(bam.references) if ref_name is None else 1):
    template_name = bam.references[r] if ref_name is None else ref_name
    if not by_query:
      alignments = []
    for read in (bam.fetch(template_name) if st is None else bam.fetch(template_name, st, en)):
      query = AlignmentRead(read.query_name, read.query_length, read.query_alignment_start, read.query_alignment_end, False, alignment=None)
      ref = AlignmentRead(template_name, bam.lengths[r], read.reference_start, read.reference_end, read.is_reverse, alignment=None)

      if ref_seqs is None: # or query_seqs is None:
        stats = read.get_cigar_stats()[0] # a list of nucleotides per operation in the following order: match, ins, del, refskip, softclip, hardclip, pad, equal, diff, back, NM tag
        acc = stats[0] / float(sum(stats))
        aln = Alignment(query, ref, m=stats[0], x=stats[8], i=stats[1], d=stats[2], acc=acc)

      else:
        # apply cigar string to query and target to draw actual aligned bp
        #query_seq = query_seqs[query.name][query.start:query.end]
        query_seq = read.query_alignment_sequence
        if query_seq is None: # this does happen (from minimap2) - not sure why
          continue
        # ref seq doesn't need to be normalized by strand, the query_alignment_sequence has already been
        ref_seq = ref_seqs[ref.name][ref.start:ref.end] #if not ref.reverse else ref_seqs[ref.name][ref.start:ref.end].translate(COMPL)[::-1] # you might think it would work to do the substring and reverse at the same time, but it doesn't if the ref starts at 0, because it tries to backtrack up to (but not including) -1, which marks the *end* of the string... TLDR; it ends up empty

        aln_string = []
        query_string = []
        target_string = []
        q = 0
        t = 0
        for a in read.cigartuples:
          if a[0] == 0: # MATCH
            for i in range(a[1]):
              if query_seq[q+i] == ref_seq[t+i]:
                aln_string.append('|')
              else:
                aln_string.append('*')
            query_string.extend(b for b in query_seq[q:q+a[1]])
            target_string.extend(b for b in ref_seq[t:t+a[1]])
            q += a[1]
            t += a[1]
          elif a[0] in [1,3]: # INS or REFSKIP
            aln_string.extend(['*']*a[1])
            query_string.extend(b for b in query_seq[q:q+a[1]])
            target_string.extend(['-']*a[1])
            q += a[1]
          elif a[0] == 2: # DEL
            aln_string.extend(['*']*a[1])
            target_string.extend(b for b in ref_seq[t:t+a[1]])
            query_string.extend(['-']*a[1])
            t += a[1]
          # soft and hard-clipped bases are not included in the query_alignment_sequence or reference positions
          '''
          elif a[0] in [4,5]: # SOFT,HARDCLIP
            aln_string.extend(['*']*a[1])
            query_string.extend(b for b in query_seq[q:q+a[1]])
            target_string.extend(['-']*a[1])
            q += a[1]
          '''
        query_string = ''.join(query_string)
        aln_string = ''.join(aln_string)
        target_string = ''.join(target_string)

        query.alignment = query_string
        ref.alignment = target_string
        aln = Alignment(query, ref)
        aln.alignment = aln_string
        aln.get_errors() # compute i,d,m,x from alignment strings

      if aln.filter(length_threshold, alignment_threshold, edge_dist_threshold):
        if by_query:
          #yield(query.name, aln)
          if query.name not in alignments:
            alignments[query.name] = []
          alignments[query.name].append(aln)
        else:
          alignments.append(aln)

    if not by_query:
      yield (template_name, alignments)

  if by_query:
    for q in alignments:
      yield (q, alignments[q])

def iter_any(infile, *args, **kwargs):
  if '.' not in infile:
    return None
  ext = infile[infile.rindex('.')+1:].lower()
  if ext == "maf":
    return iter_maf(infile, *args, **kwargs)
  if ext == "m4":
    return iter_m4(infile, *args, **kwargs)
  if ext == "m5":
    return iter_m5(infile, *args, **kwargs)
  if ext == "bam":
    return iter_bam(infile, *args, **kwargs)
  return None


if __name__ == "__main__":
  parser = argparse.ArgumentParser("Play with your formats")
  parser.add_argument("command", help="convert (maf -> m5), possibly others")
  parser.add_argument("--maf", help="MAF file")
  parser.add_argument("--m5", help="m5 file")
  args = parser.parse_args()
  if args.command == "convert":
    maf2m5(args.maf, args.m5)
  else:
    raise Exception("Command not recognized")
