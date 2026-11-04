# Observations: scripts and run records

The scripts and data behind the section "Observations" of the paper (A. Canli, Improved spherical codes in dimensions
3 to 32). Each script prints, byte for byte, the output kept in `expected/`, and every number in that section is a line
of one of these outputs. Counts are cut at 2 November 2026 18:46 (+04), when the table's page used by the paper was
saved (`../records/live_table.tsv`): runs and filings after that are left out.

## Running

Python 3.9 or later; the last six scripts below also need numpy (we used Python 3.9.6 and 3.12.13, numpy 2.5.3, on
macOS). From this folder:

    for s in lanes coverage history neighbours random_starts structure same_basin deform bench_recheck gpu_check gaps; do
      python3 -I $s.py | diff -q - expected/$s.out && echo "$s: same output"
    done

Each takes a few seconds. On stderr the numpy scripts say which values they recomputed and which they read from
`data/measured/` (see "What is not here" below).

## Which script gives which number

- 1,431 entries improved (opening): `history.py` and `coverage.py`.
- Methods: random starts 4.9e-5 to 1.3e-2 behind after ten minutes, `random_starts.py`; the 120 annealing runs that
  returned their start, `bench_recheck.py`; the 360 entries run both ways (303, 57, 88, four, 0.4%, the 10 entries
  with gains of at least 2.5e-4, d = 26, N = 121), `same_basin.py`; the table of runs, `lanes.py`.
- Which entries: `coverage.py`.
- Structure: moves after the best rotation (0.37%, 34%, rank correlation 0.90, the 36 entries, 16 in dimension 8),
  `deform.py`; the 433 codes that passed the strict test (377 of them improvements; contacts, rattlers, isostatic
  counts, antipodal pairs), `structure.py`.
- Convergence: 200 of the 1,242 sweep runs, `lanes.py`; everything else, `history.py`.
- No number in the paper: `neighbours.py` (entries whose value nearly equals that for N+1), `gpu_check.py` (the GPU
  annealer test of 26 October) and `gaps.py` (used to choose the contact tolerance).

## Data

- `../records/`: `IMPROVED.tsv` (every filing of a code that beat the table of 25 October, with its time), `UPLOADED.tsv`
  and `UPLOADED.log` (every confirmed upload, with its time and the value the table then served), `LIVE_SEEN.tsv`, the
  table of 25 October (`targets_20261025.csv`) and the page saved on 2 November at 18:46 (`live_table.tsv`). These are
  the same files `make_table.py` uses for the rest of the paper.
- `data/run_records.jsonl`: one line per completed run (4,543): the text of its `.done` file (the summary line
  `d n best table improved secs polish=converged|cut` and the independent checker's lines) and the file's modification
  time, which the scripts use to cut at the page time and to tell the programs apart. Git keeps no file times, hence
  one file. 564 records of runs that never finished (no summary or checker lines, mostly the checker's error message
  because the run wrote no code) are left out; the scripts skip such records anyway.
- Lanes (each one a queue of runs): `hitrate` (25-26 Oct test, annealing from the table's code; results in
  `data/runs/hitrate/results.txt`), `lowtest` (the same on low dimensions), `sweep` (marked entries with d 18-32: anneal
  the table's code with small kicks, then the LP polish), `sweep3` (the same, d 3-17), `sweep4` (the same, largest N·d
  first; no run finished by the page time), `deep`, `repolish` and `converge` (further polishing of our codes), `fast`
  (polish of the table's code with no annealing: the LP polish for 60 s on 31 Oct, the active-set polish `finish.py`
  from 1 Nov 00:30), `finish` (`finish.py` on our improved codes from 1 Nov 01:35), `unflag` and `uall` (the sweep's
  method on entries not marked "needs more optimization").
- `data/runs/<lane>/dD_nN.txt`: the output codes that `structure.py`, `same_basin.py`, `gaps.py` and `gpu_check.py` read,
  where they beat the table of 25 October by at least 1e-10 (fast 302, finish 333, sweep 269, gpu 1).
- `data/best_codes.tsv` and `data/best/`: our best code filed by the page time for each of the 1,431 improved entries,
  which `deform.py` compares with the table's code. For 562 entries it is the code we uploaded last
  (`../codes/dD_nN.txt`); the other 869 are in `data/best/` (389 entries never uploaded, 480 where a better code was
  filed after the last upload). For 7 of the 480 the file of that filing was later overwritten, and `data/best/` holds
  the closest code that survives (its largest cosine is within 2e-17 of the filed one; the manifest gives both).
- `data/runs/tests/*.log`: the from-scratch tests of 25 October (`d N best table improved secs`).
- `data/bench/picks.txt`: the eight entries of the annealing benchmark of 31 October.

## What is not here, and how the numbers still come out

Two kinds of file are listed but not included:

1. The table's own codes as we downloaded them on 25-26 October (31 October for unmarked entries): the 1,487 codes
   that `deform.py` (start codes), `structure.py` (comparison) and `bench_recheck.py` compare with. They are Henry
   Cohn's data, and neither the table's site nor its archive (https://hdl.handle.net/1721.1/153543) gives a licence to
   redistribute them. They can no longer be downloaded in that version either: the site serves only the current code
   of each entry, and its page of 4 November shows a different value for 1,428 of the 1,431 entries `deform.py`
   compares (they were improved since, by us or by the table's own polishing); the archive is a snapshot of
   31 July 2026 and holds the same file for one of them.
2. 293 of our run outputs that did not beat the table of 25 October by 1e-10: the 120 benchmark outputs, 77 from the
   fast lane, 91 from the sweep, 2 from the finish lane and 3 GPU outputs. They are the table's codes unchanged or
   polished in the last digits, so they are left out for the same reason.

`data/table_codes.tsv` and `data/unpublished_outputs.tsv` list these files with their size and SHA-256 (and each table
code with its URL and download time). The values the scripts compute from them are in `data/measured/`: per entry
for `deform.py` (largest cosine of our code and the moves), per code for `structure.py`, per pair for
`same_basin.py`, and per output for `bench_recheck.py` and `gpu_check.py`. These were written on 4 November 2026 by
the same scripts with `--measure`, run against the original files. Each script recomputes every value whose files
are here and checks it against `data/measured/`. Without the missing files, `structure.py` recomputes 377 of its 433
codes, `same_basin.py` 264 of its 360 moves, `deform.py` all 1,431 of our codes, and `gpu_check.py` one of its four
outputs. The comparisons with the table's codes come from `data/measured/`.

Whoever holds the original files, such as the table's maintainer, can check the rest. `python3 -I table_codes.py copy
DIR` puts every file under DIR whose SHA-256 matches the lists into `cache/` (not in git). `python3 -I table_codes.py
fetch` does the same for table codes the site still serves unchanged (one request a second; most will report
"changed"). The scripts then recompute those values too and report any difference from `data/measured/`. With all the
files, every value is recomputed and agrees exactly.

Left out: `anneal_check.py`, which tried to tell whether a code came from an annealing pass or from polishing steps
by the rotation between start and output. It was inconclusive, because the LP polish rotates codes too, so no number
in the paper uses it, and it would need the table's codes for all 1,349 runs. Also left out is a helper's speed
summary of the 31 October benchmark (no number in the paper uses it).

## Licence

CC BY 4.0, like the rest of this repository (`../LICENSE`). The values in `../records/live_table.tsv` and
`../records/targets_20261025.csv` are read from Henry Cohn's table at spherical-codes.org and remain its data.
