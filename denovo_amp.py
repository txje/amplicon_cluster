import matplotlib as mpl
mpl.use("Agg")
from matplotlib import pyplot as plt
import seaborn as sn
sn.set(style="whitegrid")

import argparse
import math
import fasta
import sys, os
import getopt
import numpy as np
import aln_formats as aln

import scipy
import scipy.stats
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as spadist

from skbio.stats.ordination import pcoa
from skbio.stats.distance import DistanceMatrix

from functools import reduce

try:
  from string import maketrans
  COMPL = maketrans("ATGC","TACG")
except: # python 3
  COMPL = str.maketrans("ATGC","TACG")


# compute the consensus among alignments given against a reference template
def aln_consensus(ref_seq, alns, st=2250, en=4450, verbosity=0):
  alleles = [[] for i in range(st,en)]

  for a in alns:

    # use on those that completely cover the target region
    if a.target.start > st or a.target.end < en:
      continue

    if a.target.reverse:
      a.alignment = a.alignment.translate(COMPL)[::-1]
      a.query.alignment = a.query.alignment.translate(COMPL)[::-1]
      a.target.alignment = a.target.alignment.translate(COMPL)[::-1]

    # advance through alignment until we reach the beginning of the target region
    t = a.target.start # target position
    al_pos = 0
    while al_pos < len(a.target.alignment) and t < st:
      if a.target.alignment[al_pos] != '-':
        t += 1
      al_pos += 1

    # walk through alignment and keep list of the "allele" corresponding to each target locus (may be empty [deletion] or multi-character insertion)
    current_al = ""
    while al_pos < len(a.target.alignment) and t < en:
      if a.query.alignment[al_pos] != '-':
        current_al += a.query.alignment[al_pos]
      if a.target.alignment[al_pos] != '-':
        alleles[t - st].append(current_al)
        current_al = ""
        t += 1
      al_pos += 1

  conseq = ""
  i = 0
  for al in alleles:
    al_frequencies = sorted([(a, al.count(a)) for a in set(al)], key=lambda a: -a[1])
    if verbosity > 1:
      print("pos {}, alleles:".format(i), al_frequencies)
    cons_al = al_frequencies[0][0]
    conseq += cons_al
    i += 1

  return conseq


# read_names is used to make sure the feature rows are in the same order as the fasta file
def compute_features(aln_file, read_names, ref_seq, ref_name, feature_file, st=2250, en=4450, binary=False):
  '''
  AAV:
  8260 aligned reads
  5175 reads cover 2250-4450
  '''

  mult = 5 if binary else 1
  # use the entire target genome
  l = en - st + 1
  variants = np.zeros((len(read_names), mult*l), dtype='u1') # a mask of the variants in all reads
  name_idx_map = {read_names[i]:i for i in range(len(read_names))}

  ct = 0
  for query, alignments in aln.iter_any(aln_file, 0, 0, 1000000000, by_query=True): # go by query just to save memory, they are all on the same target
    ct += 1
    if ct % 10000 == 0:
      print("{} alignments processed.".format(ct))
    if '/' in query:
      query = query[:query.rindex('/')]
    q = name_idx_map[query]

    if len(alignments) == 0:
      continue
    a = alignments[0] # we only reported one each (-bestn 1)

    # sometimes I want to use alignments to multiple targets, so we just ignore ones we aren't interested in right now
    if a.target.name != ref_name:
      continue

    # use on those that completely cover the target region
    if a.target.start > st or a.target.end < en + 1:
      continue

    t = a.target.start # target position

    if a.target.reverse:
      a.alignment = a.alignment.translate(COMPL)[::-1]
      a.query.alignment = a.query.alignment.translate(COMPL)[::-1]
      a.target.alignment = a.target.alignment.translate(COMPL)[::-1]

    for i in range(len(a.alignment)):
      if a.target.alignment[i] != '-':
        if t - st >= 0:
          if a.query.alignment[i] == '-': # deletion
            variants[q, mult*(t-st) + (4 if binary else 0)] = 1
          else:
            if a.query.alignment[i] == 'A':
              variants[q, mult*(t-st)] = (1 if binary else 2)
            elif a.query.alignment[i] == 'C':
              variants[q, mult*(t-st) + (1 if binary else 0)] = (1 if binary else 3)
            elif a.query.alignment[i] == 'G':
              variants[q, mult*(t-st) + (2 if binary else 0)] = (1 if binary else 4)
            elif a.query.alignment[i] == 'T':
              variants[q, mult*(t-st) + (3 if binary else 0)] = (1 if binary else 5)
        assert a.target.alignment[i] == ref_seq[t]
        t += 1
        if t > en:
          break

  np.save(feature_file, np.array(variants))
  return variants, [r for r in read_names if np.count_nonzero(variants[name_idx_map[r],:]) > 0] # aligned read names (ordered)


def distance(variants):
  #dist = spadist.squareform(spadist.pdist(variants, 'hamming')).astype('u2')
  dist = np.zeros((variants.shape[0], variants.shape[0]), dtype='i4')
  for i in range(variants.shape[0]):
    for j in range(i+1, variants.shape[0]):
      iz = np.count_nonzero(variants[i,:]) == 0
      jz = np.count_nonzero(variants[j,:]) == 0
      if (iz and not jz) or (not iz and jz):
        dist[i,j] = -1
      elif iz and jz:
        dist[i,j] = -2
      else:
        dist[i,j] = (variants[i,:] ^ variants[j,:]).sum()
      dist[j,i] = dist[i,j]

  return dist


# the ordering the dendrogram is used only to sort the rows/columns
def draw_heatmap(matrix, linkage, out_prefix, cutoff=0.7, color_gradient="white_blue"):

  if color_gradient == 'white_blue':
    cmap=BlueWhite()

  # color scale
  vmin = matrix.min()
  vmax = matrix.max()
  norm = mpl.colors.Normalize(vmin, vmax)

  default_window_width = 8
  default_window_height = 11 # 10.52

  plt.clf()
  fig = plt.figure(figsize=(default_window_width,default_window_height), dpi=100)

  # calculate positions for all elements
  margin = 0.02

  # axm, placement of heatmap for the data matrix
  [axm_x, axm_y, axm_w, axm_h] = [margin, 0.06+margin*2, 1-margin*2, 1-0.20-margin*4]

  # ax2, placement of dendrogram 2, on the top of the heatmap
  [ax2_x, ax2_y, ax2_w, ax2_h] = [margin, 0.82 + margin, 1-margin*2, 0.14]

  # axcb - placement of the color legend
  [axcb_x, axcb_y, axcb_w, axcb_h] = [margin, margin + 0.03, 1-margin*2, 0.04]

  # Compute and plot top dendrogram
  ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=False)
  # just plots the hierarchical clustering as a dendrogram, on the 'current' axis
  dendrogram = sch.dendrogram(linkage, color_threshold=cutoff)
  plt.axhline(y=cutoff, label="y = {}".format(cutoff), c='black', linestyle='--')
  ax2.set_xticks([]) # no ticks
  ax2.set_yticks([])

  # Plot distance matrix.
  axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h], frame_on=False)  # axes for the data matrix

  idx2 = dendrogram['leaves'] # apply the clustering for the array-dendrograms to the actual matrix data
  # order the matrix and flat clusters the same way as the plotted dendrogram
  matrix = matrix[idx2,:][:,idx2] # sort the same way on both dimensions

  ### taken from http://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python/3011894#3011894
  im = axm.matshow(matrix, aspect='auto', origin='upper', cmap=cmap, norm=norm)
  axm.set_xticks([]) # no ticks
  axm.set_yticks([])

  # color legend
  axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)  # axes for colorbar
  cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, orientation='horizontal')
  axcb.set_xlabel("Edit distance")

  plt.savefig("{}.heatmap_dendrogram.png".format(out_prefix), dpi=100)


def plot_distance_distr(dist_matrix, out_prefix):
  plt.clf()
  hist = np.zeros(dist_matrix.max()+1, dtype='u4')
  for i in range(dist_matrix.shape[0]):
    hist[dist_matrix[i,i+1:]] += 1
  plt.plot(hist)
  plt.savefig("{}.distance_hist.png".format(out_prefix))


def ndarray2list(nda):
    l=[]
    for line in nda.tolist():
        #for a in line:
        l.append(line[2])
    return l

def get_cutoff(nda, dist, out_prefix, threshold):
  dists = []
  groups = [[i] for i in range(dist.shape[0])]
  for a in nda:
    groupdist = dist[groups[int(a[0])],:][:,groups[int(a[1])]]
    avgdist = np.median(groupdist)
    dists.append(avgdist)
    if avgdist >= threshold:
      groups.append([])
      continue
    new_group = groups[int(a[0])] + groups[int(a[1])]
    groups.append(new_group)
    groups[int(a[0])] = []
    groups[int(a[1])] = []

  lim = min(50, len(dists))
  #lim = len(dists)
  l = ndarray2list(nda)
  l = [x for x in reversed(sorted(l))]

  lines = []

  # ------ use avg precomputed distance ------
  dists.sort(key=lambda a:-a)
  plt.clf()
  fig, ax1 = plt.subplots()
  ax2 = ax1.twinx()
  #lines.extend(ax2.plot(dists[:lim], label="Edit distance", color='b'))
  lines.extend(ax2.plot(range(lim-1), -np.diff(dists[:lim]), label="Edit dy/dx", color='c'))
  lines.extend(ax2.plot(range(lim-2), -np.diff(np.diff(dists[:lim])), label="Edit 2nd dy/dx", color='g'))
  der = -1 * np.diff(np.diff(dists))  # derivative
  der = der[~np.isnan(der)]
  #print(dists)
  print("-- edit dist --")
  mn = min(der[:min(lim, len(der)/2)])
  print("2nd derivative minimum: {}".format(mn))
  if list(der).index(mn)+1 < min(lim, len(der)/2): # there are no additional good mins
    mn2 = min(der[list(der).index(mn)+1:min(lim, len(der)/2)])
    print("2nd derivative 2nd minimum: {}".format(mn2))
    if mn2 < mn/2:
      mn = mn2
  print("minimum used: {}".format(mn))
  breakpoint = list(der).index(mn)
  print("breakpoint (position of minimum): {}".format(breakpoint))
  est_ncluster = breakpoint + 2
  print("edit dist clusters: {}".format(est_ncluster))
  print("stdv: {}".format(np.std(der)))
  print("min / stdv: {}".format(mn/np.std(der)))
  if abs(mn) / np.std(der) < 3:
    print("breakpoint reset to 0 since |{}| < 3".format(mn/np.std(der)))
    breakpoint = 0

  # ------ get cutoff height for this # of clusters from ward dendrogram ------
  lines.extend(ax1.plot(l[:lim], label="Ward distance", color='r'))
  lines.extend(ax1.plot(range(lim-1), -np.diff(l[:lim]), label="Ward dy/dx", color='orange'))
  lines.extend(ax1.plot(range(lim-2), -np.diff(np.diff(l[:lim])), label="Ward 2nd dy/dx", color='y'))
  der = -1 * np.diff(np.diff(l))  # derivative
  print("-- ward dist --")
  mn = min(der[:min(lim, len(der)/2)])
  print("2nd derivative minimum: {}".format(mn))
  mn2 = min(der[list(der).index(mn)+1:min(lim, len(der)/2)])
  print("2nd derivative 2nd minimum: {}".format(mn2))
  if mn2 < mn/2:
    mn = mn2
  print("minimum used: {}".format(mn))
  breakpoint2 = list(der).index(mn)
  print("breakpoint (position of minimum): {}".format(breakpoint2))
  est_ncluster2 = breakpoint2 + 2
  print("ward clusters: {}".format(est_ncluster2))

  if breakpoint2 > breakpoint:
    est_ncluster = est_ncluster2
    breakpoint = breakpoint2
  print("using {} breakpoint ({})".format("ward" if breakpoint2>breakpoint else "edit dist", breakpoint))

  y = l[breakpoint]
  print("dendrogram cutoff height: {}".format(y))

  #ax1.scatter(breakpoint, dists[breakpoint], c="g", s=20, label="Cutoff")
  #ax2.scatter(breakpoint, l[breakpoint], c="g", s=20, label="Cutoff")
  ax1.set_xlabel("Linkage graph depth (cutoff)")
  ax1.set_ylabel("Distance")
  plt.legend(lines, [l.get_label() for l in lines])
  plt.savefig("{}.cutoff.png".format(out_prefix))
  plt.clf()

  #groups = [g for g in groups if len(g) > 0]
  #print("{} clusters".format(len(groups)))

  return est_ncluster, y

# build a simple custom colormap
def BlueWhite():
  cdict = {'red':   ((0.0, 0.0, 0.0),
                     (1.0, 1.0, 1.0)),
           'green': ((0.0, 0.0, 0.0),
                     (1.0, 1.0, 1.0)),
           'blue':  ((0.0, 0.0, 1.0),
                     (1.0, 1.0, 1.0))
          }
  my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
  return my_cmap


def cluster(read_fa, ref_fa, aln_file, out_prefix, verbosity=0, st=None, en=None):
  # load read and reference sequences
  reads, read_names = fasta.read_fasta(read_fa, split_at_space=True)
  ref, ref_names = fasta.read_fasta(ref_fa)
  ref_name = ref_names[0]
  ref_seq = ref[ref_name] # we assume the ref has only one sequence, or the first is the primary
  ref_name = ref_name.split()[0]

  if st is None:
    st = 0
  if en is None:
    en = len(ref_seq)-1
  if verbosity > 0:
    print("Assessing {}-{} of {} bp in {}".format(st, en, len(ref_seq), ref_names[0]))

  # load or build distance matrix
  dist_matrix_file = "{}.pairwise_distance.npy".format(out_prefix)
  read_name_file = "{}.aligned_reads.txt".format(out_prefix)
  feature_file = "{}.features.npy".format(out_prefix)
  features = None
  try:
    print("Trying to load features and distance matrix...")
    dist = np.load(dist_matrix_file).astype('i4')
    np.save(dist_matrix_file, dist)
    aligned_read_names = open(read_name_file).read().strip().split('\n')
    features = np.load(feature_file)
  except Exception as e:
    print("Missing.")
    print("Computing features...")
    features, aligned_read_names = compute_features(aln_file, read_names, ref_seq, ref_name, feature_file, binary=True, st=st, en=en)
    open(read_name_file, 'w').write('\n'.join(aligned_read_names) + '\n')
    print("Computing distances...")
    dist = distance(features)
    np.save(dist_matrix_file, dist)

  # only distances between aligned reads (filter out any rows/cols with any -2)
  aln_indices = [d for d in range(dist.shape[0]) if -2 not in dist[d,:]]
  aln_dist = dist[aln_indices,:][:,aln_indices]
  if aln_dist.shape[0] == 0:
    print("No reads aligned to {} ({} bp) from {} - {}".format(ref_name, len(ref_seq), st, en))
    return

  print("Plotting distance distribution...")
  plot_distance_distr(aln_dist, out_prefix)

  compressed_dist_matrix = spadist.squareform(aln_dist)

  print("Agglomerative clustering (linkage)...")
  # does pretty simple agglomerative hierarchical clustering (think neighbor-joining)
  linkage = sch.linkage(compressed_dist_matrix, method="ward", metric="euclidean") # same as ward(compressed_dist_matrix)
  np.save("{}.linkage.npy".format(out_prefix), linkage)

  # convert hierarchical clustering (from linkage) to flat clustering:
  n_clusters, cutoff = get_cutoff(linkage, aln_dist, out_prefix, threshold=1000)

  print("Cutoff: {}".format(cutoff))
  cluster_indices = sch.fcluster(linkage, cutoff-1, 'distance') # this is the default behavior of dendrogram
  #print(list(cluster_indices))
  ai = np.array(aln_indices)
  np.save("{}.aligned_indices.npy".format(out_prefix), ai)
  np.save("{}.cluster_indices.npy".format(out_prefix), cluster_indices)

  print("Drawing heatmap...")
  draw_heatmap(aln_dist, linkage, out_prefix, cutoff)

  n_indices = len(set(cluster_indices))
  print("{} clusters (indices) found".format(n_indices))

  # ------ PCoA and plot colored by cluster_indices ------
  d = DistanceMatrix(aln_dist)
  pcoa_result = pcoa(d)
  if verbosity > 1:
    print("Proportion explained:", pcoa_result.proportion_explained)
    print("Eigenvalues:", pcoa_result.eigvals)
    print("Samples:", pcoa_result.samples)
    print("Features:", pcoa_result.features)

  x = pcoa_result.samples["PC1"]
  for pc in [2,3,4]:
    y = pcoa_result.samples["PC{}".format(pc)]
    plt.clf()
    f, ax = plt.subplots(figsize=(8,8))
    sn.despine(f)
    sn.scatterplot(x, y, hue=cluster_indices, palette=sn.color_palette("husl", n_indices))
    lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.xlabel("PC1: {:.2f}%".format(pcoa_result.proportion_explained["PC1"]*100))
    plt.ylabel("PC{}: {:.2f}%".format(pc, pcoa_result.proportion_explained["PC{}".format(pc)]*100))
    plt.savefig("{}.pcoa_1_{}.png".format(out_prefix, pc), bbox_extra_artists=(lgd,), bbox_inches='tight')

  # ------ Generate cluster consensus seqs ------
  generate_cluster_seqs(ref_seq, ref_name, aligned_read_names, aln_file, cluster_indices, out_prefix, st, en, verbosity)


def generate_cluster_seqs(ref_seq, ref_name, read_names, aln_file, cluster_indices, out_prefix, st, en, verbosity=0):
  cluster_ids = list(set(cluster_indices)) # start at 1
  if verbosity > 0:
    print("{} clusters".format(len(cluster_ids)))

  # build read name lists for each cluster
  clusters = [[] for i in range(len(cluster_ids))]
  read_cluster_map = {}
  #print(clusters, cluster_indices, read_names)
  for i in range(len(cluster_indices)):
    clusters[cluster_indices[i]-1].append(read_names[i])
    read_cluster_map[read_names[i]] = cluster_indices[i]-1

  # gather alignments belonging to each cluster
  alns = [[] for i in range(len(cluster_ids))]
  for query, alignments in aln.iter_any(aln_file, 0, 0, 1000000000, by_query=True):
    if alignments[0].target.name != ref_name:
      continue
    if len(alignments) > 0 and query in read_cluster_map:
      alns[read_cluster_map[query]].append(alignments[0])

  # compute consensus of each cluster
  consensi = {}
  for c in range(len(clusters)):
    if verbosity > 1:
      print("Cluster {} has {} reads".format(c, len(clusters[c])))
    conseq = aln_consensus(ref_seq, alns[c], st, en, verbosity)
    consensi["cluster{}_{}reads".format(c, len(clusters[c]))] = conseq

  # write 1-round consensus sequences
  fasta.write_fasta("{}.cluster_consensus.fasta".format(out_prefix), consensi)
  for name in consensi:
    fasta.write_fasta("{}.{}.fasta".format(out_prefix, name[:name.index('_')]), {name:consensi[name]})


if __name__ == '__main__':
  sys.setrecursionlimit(10000) # the default is 1000
  parser = argparse.ArgumentParser("De novo clustering and consensus of long (full-length) amplicon sequences")
  parser.add_argument("fa", help="Read FASTA file")
  parser.add_argument("ref", help="Initial consensus or ref FASTA")
  parser.add_argument("aln", help="Initial consensus alignment (m5)")
  parser.add_argument("out", help="Output prefix (incl. path)")
  parser.add_argument("--verbosity", help="Level of output (0-3)", type=int, default=0)
  parser.add_argument("--start", help="Start position of ref to consider", type=int)
  parser.add_argument("--end", help="End position of ref to consider", type=int)
  args = parser.parse_args()
  cluster(args.fa, args.ref, args.aln, args.out, args.verbosity, args.start, args.end)
