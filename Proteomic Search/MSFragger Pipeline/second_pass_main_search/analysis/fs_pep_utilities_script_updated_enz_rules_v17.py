#!/usr/bin/env python3
"""
pep_util_v101f.py — Peptide evidence utilities aligned to v101f peptide2parent

Works with peptide2parent TSVs from v101f, where each row is:
    <PEPTIDE>\t<header_with_pipe_separated_key=val_pairs>

Key updates vs older utility:
- Robust header parsing (order-agnostic; survives extra fields in v101f).
- Understands v101f novelty classes:
    frameshift, transitional, zero_frame,
    zero_frame_with_minus_one_stop, zero_frame_with_minus_two_stop
- Supports frameshift-only exports (zero_frame may be absent).
- Accepts both `enz` or `enzyme`, and `digest` or `digest_span`.
- Keeps your set logic, including the roll-up set
  `zero_frame_with_frameshift_stop_only`.
"""

import argparse, gzip, os, sys, re
import pandas as pd
import itertools
from collections import defaultdict, Counter

# ----------------------------
# Constants & helpers
# ----------------------------

PEP_SETS = (
    'frameshift_only',
    'transitional_only',
    'zero_frame_only',
    'transitional_and_frameshift',
    'zero_frame_with_frameshift_stop_only',  # aggregate of evidence-only stop classes
    'mixture_with_zero_frame',
    'not_found',
    'no_zero_frame_header',                   # convenience filter (no canonical support)
)

# v101f header novelty tokens we expect to see
NOVELTY_KEYS = {
    'frameshift',
    'transitional',
    'zero_frame',
    'zero_frame_with_minus_one_stop',
    'zero_frame_with_minus_two_stop',
}

# A regex to pull novelty quickly when we don't want to fully parse
NOV_RX = re.compile(
    r'novelty=(frameshift|transitional|zero_frame|zero_frame_with_minus_one_stop|zero_frame_with_minus_two_stop)',
    re.I
)

EVIDENCE_STOP_CLASSES = {
    'zero_frame_with_minus_one_stop',
    'zero_frame_with_minus_two_stop',
}


# --- new helpers for v2 headers ---
INT_FIELDS = {
    'tmd_start','tmd_end',
    'frame_transition_aa','frame_transition_rel_aa',
    'anchor_strict_aa','anchor_final_aa',
}

KEY_ALIASES = {
    # tolerate both names everywhere
    'enz': 'enzyme',
    'digest': 'digest_span',
}

def _normalize_keys(d):
    """Apply KEY_ALIASES to a parsed header dict (non-destructive for originals)."""
    out = dict(d)
    for k_src, k_dst in KEY_ALIASES.items():
        if k_src in out and k_dst not in out:
            out[k_dst] = out[k_src]
    return out

def _parse_int_fields(d):
    """Best-effort int cast for known numeric fields."""
    for k in list(d.keys()):
        if k in INT_FIELDS:
            try:
                d[k] = int(d[k])
            except Exception:
                pass
    return d

def _parse_per_enzyme_upstream(d):
    """
    Convert 'per_enzyme_upstream' like 'aspn:12,chymo:20,...' into
    a dict under 'per_enzyme_upstream_dict', keep original string.
    """
    s = d.get('per_enzyme_upstream')
    if not s:
        return d
    out = {}
    for tok in s.split(','):
        if ':' in tok:
            k, v = tok.split(':', 1)
            k = k.strip(); v = v.strip()
            try:
                out[k] = int(v)
            except Exception:
                continue
    d['per_enzyme_upstream_dict'] = out
    return d

def _novelty_from(hdr, hdict):
    nov = (hdict.get('novelty') or '').lower()
    if not nov:
        m = NOV_RX.search(hdr)
        nov = m.group(1).lower() if m else ''
    return nov


def _zopen(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)

def load_query_peptides(path):
    """Read text or FASTA; return set of peptide strings (uppercased, deduped)."""
    peptides = set()
    with _zopen(path, 'rt') as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith('>'):
                seq_chunks = []
                for l in itertools.chain([line], fh):
                    l = l.strip()
                    if not l:
                        continue
                    if l.startswith('>'):
                        if seq_chunks:
                            peptides.add(''.join(seq_chunks).upper())
                        seq_chunks = []
                    else:
                        seq_chunks.append(l)
                if seq_chunks:
                    peptides.add(''.join(seq_chunks).upper())
            else:
                peptides.add(s.upper())
                for l in fh:
                    l = l.strip()
                    if l:
                        peptides.add(l.upper())
            break
    return peptides


def load_query_with_sources(path):
    """
    Returns (query_set, pep2sources)
      - query_set: set[str] of peptides (uppercased)
      - pep2sources: dict[str, set[str]] mapping peptide -> {source,...}
        (empty dict if the query isn't a TSV with sources)
    Recognizes TSV header with columns: 'source' and one of
    'peptide_sequence'|'sequence'|'peptide'.
    """
    qset = set()
    pep2src = defaultdict(set)

    with _zopen(path, 'rt') as fh:
        first = fh.readline()
        if not first:
            return qset, pep2src

        # TSV with header?
        hdr = [c.strip().lower() for c in first.rstrip('\n').split('\t')]
        if '\t' in first and 'source' in hdr and any(c in hdr for c in ('peptide_sequence','sequence','peptide')):
            # column indices
            i_src = hdr.index('source')
            for cand in ('peptide_sequence','sequence','peptide'):
                if cand in hdr:
                    i_pep = hdr.index(cand)
                    break
            for ln in fh:
                if not ln.strip(): continue
                cols = ln.rstrip('\n').split('\t')
                if i_pep >= len(cols) or i_src >= len(cols): continue
                pep = cols[i_pep].strip().upper()
                src = cols[i_src].strip()
                if not pep: continue
                qset.add(pep)
                if src:
                    pep2src[pep].add(src)
            return qset, pep2src

        # otherwise: FASTA or plain list
        # reuse existing loader’s logic
        # include first line back into the stream
    # Fall back to your existing loader
    qset = load_query_peptides(path)
    return qset, pep2src


def parse_header_to_dict(hdr):
    out = {}
    for tok in hdr.split('|'):
        if '=' not in tok:
            continue
        k, v = tok.split('=', 1)
        k = k.strip(); v = v.strip()
        if k and k not in out:
            out[k] = v
    # v2 enrichments
    out = _normalize_keys(out)
    out = _parse_int_fields(out)
    out = _parse_per_enzyme_upstream(out)
    return out

def parse_header_selective(hdr, select_keys=None, do_int=False, parse_upstream_dict=False):
    """
    Cheap header parser for streaming:
      - If select_keys is provided, only keep those keys (if present).
      - Optionally int-cast known INT_FIELDS.
      - Optionally parse per_enzyme_upstream into a dict (off by default).
    """
    d = {}
    if select_keys:
        want = set(k.strip() for k in select_keys.split(',') if k.strip())
    else:
        want = None

    for tok in hdr.split('|'):
        if '=' not in tok:
            continue
        k, v = tok.split('=', 1)
        k = k.strip(); v = v.strip()
        if not k or k in d:
            continue
        if want is None or k in want:
            # normalize common aliases
            if k in KEY_ALIASES:
                d[KEY_ALIASES[k]] = v
            d[k] = v

    if do_int:
        for k in list(d.keys()):
            if k in INT_FIELDS:
                try:
                    d[k] = int(d[k])
                except Exception:
                    pass

    if parse_upstream_dict and 'per_enzyme_upstream' in d:
        out = {}
        for tok in d['per_enzyme_upstream'].split(','):
            if ':' in tok:
                kk, vv = tok.split(':', 1)
                kk = kk.strip(); vv = vv.strip()
                try:
                    out[kk] = int(vv)
                except Exception:
                    continue
        d['per_enzyme_upstream_dict'] = out

    # ensure aliases present if only one form was seen
    d = _normalize_keys(d)
    return d


def peptide2parent_iter(map_tsv):
    with _zopen(map_tsv) as fh:
        first = fh.readline()
        if first and first.count('\t') >= 1 and not first.startswith('#'):
            first_pep = first.rstrip('\n').split('\t', 1)[0]
            if first_pep and not any(c.isspace() for c in first_pep):
                pep, hdr = first.rstrip('\n').split('\t', 1)
                hdict = parse_header_to_dict(hdr)
                yield pep.upper(), hdr, hdict
        for ln in fh:
            ln = ln.rstrip('\n')
            if not ln:
                continue
            try:
                pep, hdr = ln.split('\t', 1)
            except ValueError:
                continue
            hdict = parse_header_to_dict(hdr)
            yield pep.upper(), hdr, hdict

def peptide2parent_iter_raw(map_tsv):
    """
    Yield (pep_upper, hdr_str) without parsing the header into a dict.
    Much cheaper for streaming use-cases (fasta/dedup).
    """
    with _zopen(map_tsv) as fh:
        first = fh.readline()
        if first and first.count('\t') >= 1 and not first.startswith('#'):
            first_pep = first.rstrip('\n').split('\t', 1)[0]
            if first_pep and not any(c.isspace() for c in first_pep):
                pep, hdr = first.rstrip('\n').split('\t', 1)
                yield pep.upper(), hdr
        for ln in fh:
            ln = ln.rstrip('\n')
            if not ln:
                continue
            try:
                pep, hdr = ln.split('\t', 1)
            except ValueError:
                continue
            yield pep.upper(), hdr

# ----------------------------
# Build peptide -> set mapping
# ----------------------------

def build_peptide_set_map(map_tsv):
    # Count by novelty per peptide
    counts = defaultdict(lambda: {
        'frameshift': 0,
        'transitional': 0,
        'zero_frame': 0,
        'zero_frame_with_minus_one_stop': 0,
        'zero_frame_with_minus_two_stop': 0
    })

    for pep, hdr, hdict in peptide2parent_iter(map_tsv):
        nov = _novelty_from(hdr, hdict)
        if nov in counts[pep]:
            counts[pep][nov] += 1

    pep2set = {}
    for pep, c in counts.items():
        fs  = c['frameshift']
        ts  = c['transitional']
        ca  = c['zero_frame']
        z1s = c['zero_frame_with_minus_one_stop']
        z2s = c['zero_frame_with_minus_two_stop']

        if fs and not ts and not ca and not (z1s or z2s):
            label = 'frameshift_only'
        elif ts and not fs and not ca and not (z1s or z2s):
            label = 'transitional_only'
        elif ca and not fs and not ts and not (z1s or z2s):
            label = 'zero_frame_only'
        elif fs and ts and not ca and not (z1s or z2s):
            label = 'transitional_and_frameshift'
        elif (z1s + z2s) > 0 and fs == 0 and ts == 0 and ca == 0:
            label = 'zero_frame_with_frameshift_stop_only'
        else:
            total = fs + ts + ca + z1s + z2s
            label = 'not_found' if total == 0 else 'mixture_with_zero_frame'
        pep2set[pep] = label

    return pep2set

def _allowed_set(args_set):
    """Return the set of labels to include for a given args.set."""
    if args_set == 'no_zero_frame_header':
        # also include the new evidence-only stop set
        return {
            'frameshift_only',
            'transitional_only',
            'transitional_and_frameshift',
            'zero_frame_with_frameshift_stop_only'
        }
    else:
        return {args_set}

# ----------------------------
# headers subcommand
# ----------------------------

def cmd_headers(args, pep2set):
    """
    Streaming, low-memory headers writer.
    - Filters early by novelty (zero_frame hidden unless --include-zero-frame).
    - If --query is provided, only emit peptides in that set (early skip).
    - If --set is provided, pep2set must have been prebuilt in main() (lazy).
    - If --lean is provided, writes directly to TSV (no DataFrame in RAM).
    - --cols lets you limit the header keys to include.
    """
    query_set = load_query_peptides(args.query) if args.query else None
    keep_zero = bool(args.include_zero_frame)

    # Column planning: determine header columns from the first kept record.
    header_cols = None  # ['peptide_seq'] + ordered header keys

    # When streaming, we open output and write header once.
    outfh = None
    writer = None  # we just use .write on outfh

    def maybe_open_and_write_header(hdict):
        nonlocal outfh, writer, header_cols
        if outfh is not None:
            return
        # pick columns
        base_keys = list(hdict.keys())

        if args.cols:
            wanted = [k.strip() for k in args.cols.split(',') if k.strip()]
            ordered = [k for k in wanted if k in hdict]
        else:
            prefer_first = [
                'novelty','enzyme','digest_span','enz','digest',
                'ensembl_transcript_id','ensembl_gene_id','gene_name',
                'tmd_start','tmd_end','tmd_seq','ss_seq','tmd-SS_gap_size',
                'fs_type','frame_transition','terminal_stop_codon',
                'frame_transition_aa','frame_transition_rel_aa',
                'anchor_rules','upstream_pad','anchor_strict_aa','anchor_final_aa',
                'per_enzyme_upstream'
            ]
            ordered = [k for k in prefer_first if k in hdict] + [k for k in base_keys if k not in prefer_first]

        header_cols = ['peptide_seq'] + ordered

        # open writer
        if args.out:
            outfh = gzip.open(args.out, 'wt') if args.out.endswith('.gz') else open(args.out, 'wt')
        else:
            outfh = sys.stdout
        writer = outfh
        writer.write('\t'.join(header_cols) + '\n')

    # iterate & stream
    n_written = 0
    try:
        for pep, hdr, _ in peptide2parent_iter(args.map):
            # early filter by query
            if query_set and pep not in query_set:
                continue

            # novelty filter (cheap): we only need novelty and a subset of keys
            nov = _novelty_from(hdr, {})  # NOV_RX on the raw header
            if (not keep_zero) and nov == 'zero_frame':
                continue

            # set filter (requires pep2set only if you asked for it)
            if args.set:
                allowed = _allowed_set(args.set)
                if pep2set.get(pep, 'not_found') not in allowed:
                    continue

            # parse header minimally
            if args.lean:
                hdict = parse_header_selective(hdr, select_keys=args.cols or None,
                                               do_int=False, parse_upstream_dict=False)
                # ensure we include 'novelty' if present in header (not required)
                if 'novelty' not in hdict and nov:
                    hdict['novelty'] = nov
            else:
                hdict = parse_header_to_dict(hdr)

            # plan columns and open output on first kept record
            if header_cols is None:
                maybe_open_and_write_header(hdict)

            # emit row
            row_vals = [pep] + [str(hdict.get(k, '')) for k in header_cols[1:]]
            writer.write('\t'.join(row_vals) + '\n')
            n_written += 1

    finally:
        if outfh is not None and outfh is not sys.stdout:
            outfh.close()

    if args.out and outfh is not sys.stdout:
        print(f"Wrote {n_written:,} rows → {args.out}")


# ----------------------------
# classify subcommand
# ----------------------------

def cmd_classify(args, pep2set):
    counts = defaultdict(lambda: {
        'frameshift': 0,
        'transitional': 0,
        'zero_frame': 0,
        'zero_frame_with_minus_one_stop': 0,
        'zero_frame_with_minus_two_stop': 0
    })

    # load query + optional sources
    if args.query:
        query_set, pep2sources = load_query_with_sources(args.query)
    else:
        query_set, pep2sources = None, {}

    for pep, hdr, hdict in peptide2parent_iter(args.map):
        if query_set and pep not in query_set:
            continue
        nov = _novelty_from(hdr, hdict)
        if nov in counts[pep]:
            counts[pep][nov] += 1

    all_peps = sorted(query_set) if query_set else sorted(counts)
    if args.set:
        allowed = _allowed_set(args.set)
        all_peps = [p for p in all_peps if pep2set.get(p) in allowed]

    rows = []
    summary = Counter()
    not_found = []

    for pep in all_peps:
        ca  = counts[pep]['zero_frame']
        fs  = counts[pep]['frameshift']
        ts  = counts[pep]['transitional']
        z1s = counts[pep]['zero_frame_with_minus_one_stop']
        z2s = counts[pep]['zero_frame_with_minus_two_stop']
        tot = ca + fs + ts + z1s + z2s

        if tot == 0:
            pset = 'not_found'; not_found.append(pep)
        elif fs and not ts and not ca and not (z1s or z2s):
            pset = 'frameshift_only'
        elif ts and not fs and not ca and not (z1s or z2s):
            pset = 'transitional_only'
        elif ca and not fs and not ts and not (z1s or z2s):
            pset = 'zero_frame_only'
        elif fs and ts and not ca and not (z1s or z2s):
            pset = 'transitional_and_frameshift'
        elif (z1s + z2s) > 0 and fs == 0 and ts == 0 and ca == 0:
            pset = 'zero_frame_with_frameshift_stop_only'
        else:
            pset = 'mixture_with_zero_frame'

        summary[pset] += 1
        rows.append({
            'peptide':        pep,
            'n_frameshift':   fs,
            'n_transitional': ts,
            'n_zero_frame':   ca,
            'n_zero_frame_with_minus_one_stop': z1s,
            'n_zero_frame_with_minus_two_stop': z2s,
            'n_total':        tot,
            'category':       pset
        })

    # write main classify table
    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)
    print(f"Wrote {len(rows):,} peptides to {args.out}")

    # print summary
    print("\nPeptide-set summary")
    for k in ('frameshift_only','transitional_only','zero_frame_only',
              'transitional_and_frameshift','zero_frame_with_frameshift_stop_only',
              'mixture_with_zero_frame','not_found'):
        print(f"{k:<36}: {summary[k]:>6}")

    # optional: write not_found table
    if args.not_found_table:
        # If no query was provided, not_found is usually empty (since we only see peps from map).
        # Still write a headered file for consistency.
        has_sources = bool(pep2sources)
        with open(args.not_found_table, 'w') as fh:
            if has_sources:
                fh.write('peptide\tsource\n')
            else:
                fh.write('peptide\n')
            for pep in not_found:
                srcs = sorted(pep2sources.get(pep, []))
                if has_sources and srcs:
                    # one row per (peptide, source)
                    for s in srcs:
                        fh.write(f"{pep}\t{s}\n")
                else:
                    fh.write(f"{pep}\n")
        print(f"Wrote not_found table to {args.not_found_table}")

    if args.summary_out:
        with open(args.summary_out,'w') as fh:
            fh.write(f"Wrote {len(rows):,} peptides to {args.out}\n\n")
            fh.write("Peptide-set summary:\n------------------------------------\n")
            for k in ('frameshift_only','transitional_only','zero_frame_only',
                      'transitional_and_frameshift','zero_frame_with_frameshift_stop_only',
                      'mixture_with_zero_frame','not_found'):
                fh.write(f"{k:<36}: {summary[k]:>6}\n")
            if not_found:
                fh.write("\npeptides not found:\n")
                for p in not_found:
                    fh.write(f"  {p}\n")
        print(f"Wrote summary to {args.summary_out}")

# ----------------------------
# motif subcommand
# ----------------------------

def cmd_motif(args, pep2set):
    motif_df = pd.read_csv(args.motif, sep=',', comment='#')
    want_keys = {
        (row['transcript_id'], row['gene_id'], row['gene_name'],
         int(row['tm_start']), int(row['tm_end']), row['tm_domain'],
         row['slipsite_seq'], int(row['gap_size']))
        for _, row in motif_df.iterrows()
    }

    if args.filter:
        pep_to_all_keys = defaultdict(set)
        for pep, hdr, hdict in peptide2parent_iter(args.map):
            try:
                key = (
                    hdict.get('ensembl_transcript_id',''),
                    hdict.get('ensembl_gene_id',''),
                    hdict.get('gene_name',''),
                    int(hdict.get('tmd_start','0')),
                    int(hdict.get('tmd_end','0')),
                    hdict.get('tmd_seq',''),
                    hdict.get('ss_seq',''),
                    int(hdict.get('tmd-SS_gap_size','0'))
                )
            except Exception:
                continue
            pep_to_all_keys[pep].add(key)
        filtered_peps = {pep for pep, keys in pep_to_all_keys.items() if keys <= want_keys}
    else:
        filtered_peps = None

    out_rows = []
    pep2hdr_first = {}
    pep2hdr_all = defaultdict(list)
    cols = None

    for pep, hdr, hdict in peptide2parent_iter(args.map):
        if filtered_peps is not None and pep not in filtered_peps:
            continue
        try:
            key = (
                hdict.get('ensembl_transcript_id',''),
                hdict.get('ensembl_gene_id',''),
                hdict.get('gene_name',''),
                int(hdict.get('tmd_start','0')),
                int(hdict.get('tmd_end','0')),
                hdict.get('tmd_seq',''),
                hdict.get('ss_seq',''),
                int(hdict.get('tmd-SS_gap_size','0'))
            )
        except Exception:
            continue
        if key not in want_keys:
            continue
        if args.set:
            allowed = _allowed_set(args.set)
            if pep2set.get(pep, 'not_found') not in allowed:
                continue

        if cols is None:
            cols = ['peptide_seq'] + list(hdict.keys())
        vals = [hdict.get(k, '') for k in cols[1:]]
        out_rows.append([pep] + vals)

        if pep not in pep2hdr_first:
            pep2hdr_first[pep] = hdr
        pep2hdr_all[pep].append(hdr)

    if not out_rows:
        print("No motif-driven peptides found.")
        return

    df = pd.DataFrame(out_rows, columns=cols)

    if args.out:
        ext = os.path.splitext(args.out)[1].lower()
        if ext in ('.tsv','.txt'):
            df.to_csv(args.out, sep='\t', index=False)
        elif ext == '.csv':
            df.to_csv(args.out, index=False)
        else:
            sys.exit(f"[ERROR] motif --out must be .tsv/.txt or .csv, got {args.out!r}")
        print(f"Wrote {len(df):,} rows → {args.out}")

    if args.fa:
        opener = gzip.open if args.fa.endswith('.gz') else open
        with opener(args.fa, 'wt') as fh:
            if args.dedup:
                for pep, hdr in pep2hdr_first.items():
                    fh.write(f">{hdr}\n{pep}\n")
            else:
                for pep, hdrs in pep2hdr_all.items():
                    for hdr in hdrs:
                        fh.write(f">{hdr}\n{pep}\n")
        print(f"Wrote FASTA to {args.fa}")

# ----------------------------
# fasta subcommand
# ----------------------------

def cmd_fasta(args, pep2set):
    """
    Stream out FASTA records without building an in-RAM peptide->headers map.

    - Reads the map TSV line-by-line (raw).
    - Early-filters by query, set, novelty.
    - Supports --dedup (first header per peptide wins).
    - If not --dedup, emits all headers that match filters; optionally cap per-peptide.
    """
    default_keep = {
        'frameshift','transitional','zero_frame',
        'zero_frame_with_minus_one_stop','zero_frame_with_minus_two_stop'
    }
    keep = set(args.novelty) if args.novelty else default_keep

    query_set = load_query_peptides(args.query) if args.query else None
    allowed_sets = _allowed_set(args.set) if args.set else None

    seen = set()
    counts_per_pep = defaultdict(int)   # <-- move OUTSIDE the loop

    opener = gzip.open if args.out.endswith('.gz') else open
    with opener(args.out, 'wt') as out:
        for pep, hdr in peptide2parent_iter_raw(args.map):
            if query_set and pep not in query_set:
                continue

            if allowed_sets is not None and pep2set.get(pep, 'not_found') not in allowed_sets:
                continue

            nov = _novelty_from(hdr, {})  # regex on raw header
            if nov and nov not in keep:
                continue

            if args.dedup:
                # --dedup takes precedence over --max-per-peptide
                if pep in seen:
                    continue
                seen.add(pep)
                out.write(f">{hdr}\n{pep}\n")
            else:
                if args.max_per_peptide is not None:
                    if counts_per_pep[pep] >= args.max_per_peptide:
                        continue
                    counts_per_pep[pep] += 1
                out.write(f">{hdr}\n{pep}\n")

    print(f"Wrote FASTA to {args.out}")

# ----------------------------
# summarize subcommand
# ----------------------------

def cmd_summarize(args):
    """
    One-row summary per peptide:
      n_headers, n_frameshift, n_transitional, n_zero_frame,
      n_zero_frame_with_minus_one_stop, n_zero_frame_with_minus_two_stop,
      transcript_ids, gene_ids, tmd_ranges, slipsite_seqs,
      enzymes, digest_coords, novelty_types[, set]
    """
    hdr_counts = Counter()
    nov_counts = defaultdict(lambda: Counter())
    tx_sets    = defaultdict(set)
    gene_sets  = defaultdict(set)
    tmd_ranges = defaultdict(list)
    slip_sets  = defaultdict(set)
    enz_sets   = defaultdict(set)
    dig_lists  = defaultdict(list)
    nov_types  = defaultdict(set)

    qset = load_query_peptides(args.query) if args.query else None

    for pep, hdr, hdict in peptide2parent_iter(args.map):
        if qset and pep not in qset:
            continue
        hdr_counts[pep] += 1

        nov = _novelty_from(hdr, hdict)
        if nov:
            nov_counts[pep][nov] += 1
            nov_types[pep].add(nov)

        tx = hdict.get('ensembl_transcript_id')
        if tx: tx_sets[pep].add(tx)

        gid = hdict.get('ensembl_gene_id')
        if gid: gene_sets[pep].add(gid)

        ts, te = hdict.get('tmd_start'), hdict.get('tmd_end')
        if ts is not None and te is not None and str(ts).isdigit() and str(te).isdigit():
            tmd_ranges[pep].append(f"{ts}-{te}")

        ss = hdict.get('ss_seq')
        if ss: slip_sets[pep].add(ss)

        # Accept both 'enz' and 'enzyme'
        enz = hdict.get('enzyme')
        if enz: enz_sets[pep].add(enz)

        # Accept both 'digest' and 'digest_span'
        dig = hdict.get('digest_span')
        if dig: dig_lists[pep].append(dig)

    rows = []
    peptides = sorted(qset) if qset else sorted(hdr_counts)
    for pep in peptides:
        total = hdr_counts.get(pep, 0)
        if total == 0:
            continue
        ca  = nov_counts[pep]['zero_frame']
        fs  = nov_counts[pep]['frameshift']
        ts  = nov_counts[pep]['transitional']
        z1s = nov_counts[pep]['zero_frame_with_minus_one_stop']
        z2s = nov_counts[pep]['zero_frame_with_minus_two_stop']

        # classify into one of our sets
        if ca==0 and fs>0 and ts==0 and z1s==0 and z2s==0:
            pset = 'frameshift_only'
        elif ca==0 and ts>0 and fs==0 and z1s==0 and z2s==0:
            pset = 'transitional_only'
        elif ca>0 and fs==0 and ts==0 and z1s==0 and z2s==0:
            pset = 'zero_frame_only'
        elif ca==0 and fs>0 and ts>0 and z1s==0 and z2s==0:
            pset = 'transitional_and_frameshift'
        elif (z1s+z2s)>0 and fs==0 and ts==0 and ca==0:
            pset = 'zero_frame_with_frameshift_stop_only'
        else:
            pset = 'mixture_with_zero_frame'

        no_ca = (ca == 0)
        if args.set:
            if args.set == 'no_zero_frame_header' and not no_ca:
                continue
            if args.set != 'no_zero_frame_header' and pset != args.set:
                continue

        rows.append({
            'peptide':        pep,
            'n_headers':      total,
            'n_frameshift':   fs,
            'n_transitional': ts,
            'n_zero_frame':   ca,
            'n_zero_frame_with_minus_one_stop': z1s,
            'n_zero_frame_with_minus_two_stop': z2s,
            'transcript_ids': ';'.join(sorted(tx_sets[pep])),
            'gene_ids':       ';'.join(sorted(gene_sets[pep])),
            'tmd_ranges':     ';'.join(tmd_ranges[pep]),
            'slipsite_seqs':  ';'.join(sorted(slip_sets[pep])),
            'enzymes':        ';'.join(sorted(enz_sets[pep])),
            'digest_coords':  ';'.join(dig_lists[pep]),
            'novelty_types':  ';'.join(sorted(nov_types[pep])),
            'set':            'no_zero_frame_header' if no_ca else pset
        })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.out, sep='\t', index=False)
    print(f"Wrote {len(rows):,} peptides to {args.out}")

# ----------------------------
# CLI
# ----------------------------

def main():
    p = argparse.ArgumentParser(description="Multi-purpose peptide utility (v101f-aware)")
    sub = p.add_subparsers(dest='cmd')

    ph = sub.add_parser('headers', help='explode headers (non-zero-frame by default)')
    ph.add_argument('-m','--map', required=True)
    ph.add_argument('-q','--query')
    ph.add_argument('-o','--out')
    ph.add_argument('--set', choices=PEP_SETS, help="filter by peptide set")
    ph.add_argument('--include-zero-frame', action='store_true',
                help='Also include zero_frame headers (default is to hide them).')

    ph.add_argument('--lean', action='store_true',
                help='Use a streaming, low-memory mode (no DataFrame; writes TSV directly).')
    ph.add_argument('--cols', help='Comma-separated list of header keys to include (default = all keys seen).')

    pc = sub.add_parser('classify', help='classify peptides into evidence sets')
    pc.add_argument('-m','--map', required=True)
    pc.add_argument('-q','--query')
    pc.add_argument('-o','--out', required=True)
    pc.add_argument('--summary-out', help='write summary + missing list here')
    pc.add_argument('--set', choices=PEP_SETS, help="filter by peptide set")
    pc.add_argument('--not-found-table', help="Write a TSV of not_found peptides; includes 'source' if --query is a TSV with source/peptide columns")

    pm = sub.add_parser('motif', help='join with motif CSV and optionally emit FASTA')
    pm.add_argument('-m','--map', required=True)
    pm.add_argument('--motif', required=True)
    pm.add_argument('-o','--out', help="TSV/CSV output")
    pm.add_argument('--fa', help="write a FASTA of peptides (.fa/.fa.gz/.fasta)")
    pm.add_argument('--dedup', action='store_true', help="emit ≤1 FASTA record per peptide")
    pm.add_argument('--set', choices=PEP_SETS, help="filter by peptide set")
    pm.add_argument('--filter', action='store_true', help="only keep peptides exclusive to these motifs")

    pf = sub.add_parser('fasta', help='emit FASTA for given peptide list and novelty filters')
    pf.add_argument('-m','--map', required=True)
    pf.add_argument('-q','--query', required=True)
    pf.add_argument('-o','--out', required=True)
    pf.add_argument('--novelty', nargs='+', choices=[
        'frameshift','transitional','zero_frame',
        'zero_frame_with_minus_one_stop','zero_frame_with_minus_two_stop'
    ])
    pf.add_argument('--set', choices=PEP_SETS, help="filter by peptide set")
    pf.add_argument('--dedup', action='store_true', help="emit ≤1 FASTA record per peptide")

    pf.add_argument('--max-per-peptide', type=int, default=None,
                help='Emit at most N headers per peptide (implies non-dedup).')


    ps = sub.add_parser('summarize', help='one-row per peptide summary')
    ps.add_argument('-m','--map', required=True)
    ps.add_argument('-q','--query')
    ps.add_argument('--set', choices=PEP_SETS, help="filter by peptide set")
    ps.add_argument('-o','--out', required=True)

    args = p.parse_args()
    if not args.cmd:
        p.print_help(); sys.exit(1)

    # Build once; re-use
    need_pep2set = (args.cmd in {'classify','motif','fasta','summarize'}) or (args.cmd == 'headers' and args.set)
    pep2set = build_peptide_set_map(args.map) if need_pep2set else {}

    if args.cmd == 'headers':
        cmd_headers(args, pep2set)
    elif args.cmd == 'classify':
        cmd_classify(args, pep2set)
    elif args.cmd == 'motif':
        cmd_motif(args, pep2set)
    elif args.cmd == 'fasta':
        cmd_fasta(args, pep2set)
    elif args.cmd == 'summarize':
        cmd_summarize(args)
    else:
        p.print_help(); sys.exit(1)

if __name__ == '__main__':
    main()
