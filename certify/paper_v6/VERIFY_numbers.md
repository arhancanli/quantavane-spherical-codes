# Every macro of v6_numbers.tex: value and source

Written by make_numbers_v6.py; do not edit.

| macro | value | source |
|---|---|---|
| `\AboveCapCPUh` | 10 | extrapolation of the fit CPUSlope to the skipped systems |
| `\AboveCapGB` | 6.6 | 8 AboveCapMaxU^2 bytes |
| `\AboveCapMaxU` | 28{,}801 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper): smallest skipped system, max over entries |
| `\ArbMaxU` | 400 | recount/arbcheck.jsonl (recount/arbcheck.py, Arb ball arithmetic): max unknowns |
| `\BigCPU` | 3 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status largest certified system: cpu_seconds |
| `\BigEntry` | $(3,1020)$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status largest certified system |
| `\BigHyperAll` | 9 | census/exact_v4/big/*.json (exact_big.py): mode hyper |
| `\BigHyperGood` | 9 | census/exact_v4/big/*.json (exact_big.py): mode hyper, confirmed |
| `\BigHyperMaxU` | 3{,}061 | census/exact_v4/big/*.json (exact_big.py): mode hyper, max unknowns |
| `\BigIsoAll` | 12 | census/exact_v4/big/*.json (exact_big.py): mode iso |
| `\BigIsoGood` | 12 | census/exact_v4/big/*.json (exact_big.py): mode iso, confirmed |
| `\BigIsoMaxU` | 2{,}974 | census/exact_v4/big/*.json (exact_big.py): mode iso, max unknowns |
| `\BigSecondsMax` | 92 | census/exact_v4/big/*.json (exact_big.py): max seconds |
| `\BigStrictAll` | 9 | census/exact_v4/big/*.json (exact_big.py): mode strict |
| `\BigStrictGood` | 9 | census/exact_v4/big/*.json (exact_big.py): mode strict, confirmed |
| `\BigStrictMaxU` | -- | census/exact_v4/big/*.json (exact_big.py): mode strict, max unknowns |
| `\BigU` | 3{,}061 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status largest certified system |
| `\CPUHoursCensus` | 3.9 | census/summary_v4.json macro CPUHoursCensus (census/summarize.py from the per-entry outputs) |
| `\CPUHoursStrict` | 0.8 | census/summary_v4.json macro CPUHoursStrict (census/summarize.py from the per-entry outputs) |
| `\CPUMaxLoc` | 51 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: cpu_seconds max |
| `\CPUMaxPass` | 2 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: cpu_seconds max |
| `\CPUMedLoc` | 0.08 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: cpu_seconds median |
| `\CPUMedPass` | 0.07 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: cpu_seconds median |
| `\CPUSlope` | 2.7 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS/LOCALMIN with >=1000 unknowns: least-squares slope of log cpu on log U |
| `\CertShaCq` | c12bf1\ldots | census/results_v4/*.json (top-level record): code_sha256['cq.py'] (one value in all entries) |
| `\ChiMax` | $1.2\cdot10^{-27}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: max tau_extra |
| `\DiagGapRatioMed` | 1.01 | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): over the unseparated proposals, median of (largest cosine - next pair) / tau |
| `\DiagRowsumMed` | $1.1\cdot10^{10}$ | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): median |
| `\DiagRowsumMin` | $1.4$ | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted) (census record krawczyk_attempts[0].max_rowsum): min |
| `\DiagSpreadRatioMed` | 0.99 | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): over the unseparated proposals, median of (largest - smallest contact cosine) / tau |
| `\DriverShaA` | 4a7435\ldots | census/results_v4/*.json (top-level record): code_sha256['run_census.py'], the more frequent |
| `\DriverShaB` | de9662\ldots | census/results_v4/*.json (top-level record): code_sha256['run_census.py'], the other |
| `\ExOnlyDistCoordMax` | $1.3\cdot10^{-15}$ | v5/out/summary_v5.json (v5/summarize_v5.py) macro VExOnlyDistCoordMax |
| `\ExOnlyDistLSMax` | $1.2\cdot10^{-14}$ | v5/out/summary_v5.json (v5/summarize_v5.py) macro VExOnlyDistLSMax |
| `\ExOnlyN` | 18 | v5/out/summary_v5.json (v5/summarize_v5.py) macro VExOnlyN |
| `\ExactCPUMax` | 679 | census/exact_v4/strict/*.json (exact_strict.py): max cpu_seconds of STRICT |
| `\FamilyList` | $(22,47)$, $(22,48)$, $(23,50)$ | v5/out/summary_v5.json (v5/summarize_v5.py) macro VFamilyList |
| `\FlatColoursMax` | 52 | v5/out/summary_v5.json (v5/summarize_v5.py) macro VFlatColoursMax |
| `\FlatColoursMin` | 26 | v5/out/summary_v5.json (v5/summarize_v5.py) macro VFlatColoursMin |
| `\FlatDistMax` | $2.8\cdot10^{-14}$ | v5/out/summary_v5.json (v5/summarize_v5.py) macro VFlatDistMax |
| `\FlatList` | $(16,117)$, $(17,162)$, $(20,92)$, $(22,79)$, $(22,114)$, $(22,127)$, $(23,102)$, $(23,129)$, $(23,147)$, $(24,121)$, $(24,149)$, $(26,57)$, $(27,96)$ | v5/out/summary_v5.json (v5/summarize_v5.py) macro VFlatList |
| `\FlexDeltaMax` | $2.7\cdot10^{-5}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): delta, max |
| `\FlexDeltaMed` | $2.5\cdot10^{-9}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): delta, median |
| `\FlexDeltaMin` | $1.4\cdot10^{-15}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): delta = min(h/sqrt d, tau'/4), min |
| `\FlexEpsMax` | $6.4\cdot10^{-15}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): max of the Gram error bound eps |
| `\FlexHMed` | $1.2\cdot10^{-2}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): h/sqrt d, median |
| `\FlexTauRatioMin` | 0.88 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): min of tau'/tau |
| `\FlowAbovecapFlex` | 117 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) x class |
| `\FlowAbovecapNone` | 271 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) x class |
| `\FlowAbovecapVerdict` | 3 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) x class |
| `\FlowAttemptedFlex` | 302 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) x class |
| `\FlowAttemptedNone` | 79 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) x class |
| `\FlowAttemptedVerdict` | 1{,}795 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) x class |
| `\FlowNeverFlex` | 12{,}063 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) x class |
| `\FlowNeverNone` | 1{,}635 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) x class |
| `\ImpLowerA` | $6.4\cdot10^{-12}$ | v5/out/summary_v5.json (v5/summarize_v5.py) macro VImpLowerA |
| `\ImpLowerB` | $1.8\cdot10^{-13}$ | v5/out/summary_v5.json (v5/summarize_v5.py) macro VImpLowerB |
| `\InvShareMed` | 54 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS/LOCALMIN >=1000 unknowns: seconds_inverse/wall_seconds, median (wall-clock share) |
| `\LocDistMax` | $4.0\cdot10^{-6}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: max dist_ls (to P') |
| `\LocDistTubeMax` | $4.0\cdot10^{-6}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: max dist_ls + localmin_displacement |
| `\LocExtraMax` | 6{,}227 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: max |
| `\LocExtraMed` | 62 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: median |
| `\LocExtraMin` | 1 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: min contacts - isostatic_count |
| `\MAnyVerdict` | 93 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' (marked) |
| `\MPropNever` | 12{,}647 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) (marked) |
| `\MaxUnknownsLoc` | 3{,}061 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: max unknowns |
| `\MaxUnknownsPass` | 2{,}974 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: max unknowns |
| `\NAboveCapExact` | 3 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) |
| `\NAboveCapOpen` | 388 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) |
| `\NAnyVerdict` | 1{,}798 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NArbHyper` | 31 | recount/arbcheck.jsonl (recount/arbcheck.py, Arb ball arithmetic): census LOCALMIN |
| `\NArbIso` | 15 | recount/arbcheck.jsonl (recount/arbcheck.py, Arb ball arithmetic): census PASS |
| `\NArbOverlap` | 46 | recount/arbcheck.jsonl (recount/arbcheck.py, Arb ball arithmetic): cosine ball overlaps the census 30-digit interval |
| `\NArbPass` | 46 | recount/arbcheck.jsonl (recount/arbcheck.py, Arb ball arithmetic): Krawczyk test passes in Arb |
| `\NArbRun` | 46 | recount/arbcheck.jsonl (recount/arbcheck.py, Arb ball arithmetic) |
| `\NAttFailExist` | 348 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper): first attempt failed existence |
| `\NAttFailNoStress` | 18 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper): first attempt found no positive stress |
| `\NAttFailOther` | 15 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper): first attempt failed otherwise |
| `\NAttemptFail` | 381 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) without a verdict |
| `\NCPUSlope` | 526 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS/LOCALMIN with >=1000 unknowns |
| `\NCensus` | 16{,}265 | census/index.tsv (snapshot.py from the archived data set): listed, 3<=d<=32, N<=1024 |
| `\NClassVaries` | 9{,}385 | census/summary_v4.json macro NClassVaries (census/summarize.py from the per-entry outputs) |
| `\NCodesData` | 16{,}332 | census/index.tsv (snapshot.py from the archived data set): codes in data.txt |
| `\NCtlAgree` | 27 | phase2/out_controls/*.json (symstrict.py on entries certified strict by recognition): gamma* agrees with strict.py's c* to 35 digits |
| `\NCtlRun` | 152 | phase2/out_controls/*.json (symstrict.py on entries certified strict by recognition) |
| `\NCtlStrict` | 27 | phase2/out_controls/*.json (symstrict.py on entries certified strict by recognition): STRICT |
| `\NDecidedExact` | 17 | v5/out/summary_v5.json (v5/summarize_v5.py) macro VNDecidedExact minus (22,116) |
| `\NDiag` | 348 | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted) |
| `\NDiagClose` | 270 | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): next pair within 10 tau |
| `\NDiagCloseSing` | 183 | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): not separated and numerically singular |
| `\NDiagConv` | 72 | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): census Newton residual below 1e-25 |
| `\NDiagCorankOne` | 56 | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): exactly one singular value below 1e-12 relative |
| `\NDiagNoConv` | 276 | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): census Newton residual not below 1e-25 |
| `\NDiagNonsing` | 90 | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): smallest/largest singular value >= 1e-12 |
| `\NDiagRowsumBelowOne` | 0 | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): contraction factor below 1 |
| `\NDiagSep` | 78 | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): next pair more than 10 tau below the largest cosine |
| `\NDiagSepNonsing` | 3 | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): separated and nonsingular |
| `\NDiagSepSing` | 75 | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): separated and numerically singular |
| `\NDiagSingular` | 258 | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): smallest/largest singular value < 1e-12 |
| `\NDispBelow` | 80 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS/LOCALMIN: certified value (+ tau_extra) rounded up to 12 decimals below the displayed cosine |
| `\NDispBelowLiveAuthor` | 0 | v5/live/live_20261107.tsv (the live table read on 7 November 2026): entries of NDispBelow credited to the author |
| `\NDispBelowLiveImproved` | 77 | v5/live/live_20261107.tsv (the live table read on 7 November 2026): entries of NDispBelow whose live cosine is below the archived one |
| `\NDriverA` | 15{,}992 | census/results_v4/*.json (top-level record): entries with DriverShaA |
| `\NDriverB` | 273 | census/results_v4/*.json (top-level record): entries with DriverShaB |
| `\NEntriesStoppedMemory` | 2 | census/summary_v4.json macro NEntriesStoppedMemory (census/summarize.py from the per-entry outputs) |
| `\NExactRecognised` | 62 | census/exact_v4/strict/*.json (exact_strict.py): recognised |
| `\NExactRootInside` | 337 | census/validate_v4.tsv (validate_poly.py): 'root inside' |
| `\NExactStrict` | 39 | census/exact_v4/strict/*.json (exact_strict.py): STRICT |
| `\NExactTested` | 337 | census/validate_v4.tsv (validate_poly.py): rows with a verdict |
| `\NExactTestedHyper` | 298 | census/validate_v4.tsv (validate_poly.py): kind hyperstatic |
| `\NExactTestedHyperRecog` | 228 | census/validate_v4.tsv (validate_poly.py) and class strict_hyper (recognition with the table's polynomial) |
| `\NExactTestedHyperStrict` | 234 | census/validate_v4.tsv (validate_poly.py) and class |
| `\NExactTestedIso` | 39 | census/validate_v4.tsv (validate_poly.py): kind isostatic |
| `\NExactTried` | 942 | census/exact_v4/strict/*.json (exact_strict.py) |
| `\NFailedAll` | 82 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NFlagAll` | 13{,}244 | census/index.tsv (snapshot.py from the archived data set): census entries marked 'needs more optimization' |
| `\NFlatLive` | 6 | v5/out/summary_v5.json (v5/summarize_v5.py) macro VNFlatLive |
| `\NFlex` | 12{,}482 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NFlexDeltaAboveNine` | 8{,}066 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): delta >= 1e-9 |
| `\NFlexDeltaAtTauQuarter` | 12{,}374 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): delta attained by tau'/4 |
| `\NFlexTauBelowTau` | 904 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): entries with tau' < tau |
| `\NFlexTheorem` | 12{,}482 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): count |
| `\NFlexTolEight` | 3{,}003 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): tol 1e-08 |
| `\NFlexTolNine` | 2{,}071 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): tol 1e-09 |
| `\NFlexTolSeven` | 5{,}121 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): tol 1e-07 |
| `\NFlexTolSmaller` | 1{,}353 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): tol < 1e-10 |
| `\NFlexTolTen` | 934 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): tol 1e-10 |
| `\NFlexUnmarked` | 1{,}091 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): unmarked entries |
| `\NHyperTen` | 1{,}364 | census/summary_v4.json NClassHyperTen |
| `\NHypoEmptyTen` | 14{,}162 | census/summary_v4.json NClassEmptyTen + NClassHypoTen |
| `\NInside` | 668 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: code_inside_rho |
| `\NIsoTen` | 739 | census/summary_v4.json NClassIsoTen |
| `\NListed` | 16{,}269 | census/index.tsv (snapshot.py from the archived data set): codes listed in table.txt |
| `\NLiveAuthor` | 443 | v5/live/live_20261107.tsv (the live table read on 7 November 2026): credited to the author |
| `\NLocTolOther` | 31 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: tol != 1e-10 |
| `\NLocal` | 630 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NLocalAll` | 937 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN |
| `\NLocalAtLeast` | 1{,}783 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NLocalVFive` | 679 | v5/out/summary_v5.json VNLocal |
| `\NNoRho` | 13 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: rho absent |
| `\NNonStrict` | 2 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NNoneTried` | 12{,}564 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NNonoptThenCert` | 0 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status an attempt NONOPT and a later PASS/LOCALMIN |
| `\NNotListed` | 63 | census/index.tsv (snapshot.py from the archived data set): codes not listed |
| `\NNotLocal` | 15 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NNotLocalFlat` | 13 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NNotLocalNonopt` | 2 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NOpenDegree` | 29 | census/strict_v4/*.json (strict.py) status of the entries still open |
| `\NOpenExact` | 22 | v5/out/summary_v5.json (v5/summarize_v5.py) macro VNOpenExact |
| `\NOpenNoPoly` | 509 | census/strict_v4/*.json (strict.py) status of the entries still open |
| `\NOpenNotInField` | 49 | census/strict_v4/*.json (strict.py) status of the entries still open |
| `\NOpenTooMany` | 43 | census/strict_v4/*.json (strict.py) status of the entries still open |
| `\NOptAll` | 25 | census/index.tsv (snapshot.py from the archived data set): census entries marked 'Provably optimal' |
| `\NOptStrict` | 23 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' and index 'optimal' |
| `\NOutside` | 158 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: code_inside_rho false |
| `\NPassDistAbove` | 97 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: dist_ls > 1e-8 |
| `\NPassRattlers` | 631 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: rattler_placement.rattlers > 0 |
| `\NPassTolOther` | 106 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: tol != 1e-10 |
| `\NPolyCensus` | 517 | census/index.tsv (snapshot.py from the archived data set): census entries with a minimal polynomial |
| `\NPolyListed` | 521 | census/index.tsv (snapshot.py from the archived data set): listed codes with a minimal polynomial |
| `\NPropAboveCap` | 391 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper): none attempted, some skipped |
| `\NPropAttVerdict` | 1{,}795 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) and class |
| `\NPropAttempted` | 2{,}176 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper): some attempted |
| `\NPropNever` | 13{,}698 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper): none, none |
| `\NRattlerMoved` | 278 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: rattler_placement moved_first_order + moved_second_order > 0 |
| `\NRattlerStuck` | 0 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: rattler_placement.stuck > 0 |
| `\NRedModelRun` | 49 | phase2/out_red_model/*.json (symstrict.py --red model) |
| `\NRedModelStrict` | 0 | phase2/out_red_model/*.json (symstrict.py --red model): STRICT (must be 0) |
| `\NRedPermRun` | 49 | phase2/out_red_perm/*.json (symstrict.py --red perm) |
| `\NRedPermStrict` | 0 | phase2/out_red_perm/*.json (symstrict.py --red perm): STRICT (must be 0) |
| `\NRedSwapRun` | 49 | phase2/out_red_swap/*.json (symstrict.py --red swap) |
| `\NRedSwapStrict` | 0 | phase2/out_red_swap/*.json (symstrict.py --red swap): STRICT (must be 0) |
| `\NRhoDefined` | 826 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: rho present |
| `\NStoppedMemory` | 8 | census/summary_v4.json macro NStoppedMemory (census/summarize.py from the per-entry outputs) |
| `\NStrict` | 1{,}151 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NStrictExact` | 3 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NStrictHyper` | 258 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NStrictHyperFieldHigh` | 96 | census/strict_v4/*.json (strict.py): field_degree >= 3 |
| `\NStrictHyperFieldOne` | 73 | census/strict_v4/*.json (strict.py): field_degree 1 |
| `\NStrictHyperFieldTwo` | 89 | census/strict_v4/*.json (strict.py): field_degree 2 |
| `\NStrictIso` | 839 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NStrictNew` | 2 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NStrictOpenSec` | 2 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' (exact core, bar-rigid) |
| `\NStrictRotExact` | 258 | census/strict_v4/*.json (strict.py): rotation_exactly_orthogonal |
| `\NStrictSym` | 49 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NStrictVFive` | 1{,}103 | v5/out/summary_v5.json VNStrictAll (version 5 count, for the comparison) |
| `\NSurvey` | 679 | phase2/survey_*.jsonl (phase2/sym_survey.py, untrusted) |
| `\NSymCap` | 0 | phase2/out_open/*.json (phase2/symstrict.py): group above the enumeration cap |
| `\NSymDimThree` | 22 | phase2/out_open/*.json (phase2/symstrict.py) (d = 3) |
| `\NSymFieldQ` | 36 | phase2/out_open/*.json (phase2/symstrict.py): field_D = 1 |
| `\NSymFieldQuad` | 13 | phase2/out_open/*.json (phase2/symstrict.py): field_D > 1 |
| `\NSymModelCyc` | 3 | phase2/out_open/*.json (phase2/symstrict.py): model cyclic |
| `\NSymModelInv` | 14 | phase2/out_open/*.json (phase2/symstrict.py): model involutions |
| `\NSymModelIso` | 13 | phase2/out_open/*.json (phase2/symstrict.py): model isotypic |
| `\NSymModelLines` | 19 | phase2/out_open/*.json (phase2/symstrict.py): model lines |
| `\NSymNoModel` | 5 | phase2/out_open/*.json (phase2/symstrict.py): status 'no exact model' |
| `\NSymNoSym` | 531 | phase2/out_open/*.json (phase2/symstrict.py): status 'no symmetry' (search complete, only the identity) |
| `\NSymNontrivial` | 136 | phase2/survey_*.jsonl (phase2/sym_survey.py, untrusted): an automorphism other than the identity found |
| `\NSymNotSquare` | 82 | phase2/out_open/*.json (phase2/symstrict.py) |
| `\NSymOther` | 1 | phase2/out_open/*.json (phase2/symstrict.py): other failures |
| `\NSymRun` | 679 | phase2/out_open/*.json (phase2/symstrict.py) |
| `\NSymSquare` | 51 | phase2/survey_*.jsonl (phase2/sym_survey.py, untrusted): reduced system square |
| `\NSymSquareNoModel` | 2 | phase2/out_open/*.json (phase2/symstrict.py): square in the survey, no exact model |
| `\NSymSquareStrict` | 49 | phase2/out_open/*.json (phase2/symstrict.py) and phase2/survey_*.jsonl (phase2/sym_survey.py, untrusted): square in the survey and STRICT |
| `\NSymStrict` | 49 | phase2/out_open/*.json (phase2/symstrict.py): status STRICT |
| `\NSymTimeout` | 11 | phase2/out_open/*.json (phase2/symstrict.py): no automorphism found within the time limit |
| `\NSymUnmarked` | 46 | phase2/out_open/*.json (phase2/symstrict.py) (unmarked) |
| `\NTooLarge` | 1{,}903 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NTooLargeNever` | 1{,}632 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) |
| `\NUnflagAll` | 3{,}021 | census/index.tsv (snapshot.py from the archived data set): census entries not so marked |
| `\NUnion` | 7 | v5/out/summary_v5.json (v5/summarize_v5.py) macro VNUnion |
| `\NUniqueCoversCode` | 839 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: r_unique > dist_nonrattlers |
| `\NWatchdogEntries` | 273 | census/summary_v4.json macro NWatchdogEntries (census/summarize.py from the per-entry outputs) |
| `\NewColsA` | 1{,}760 | v5/out/rigidity_*.json (v5/rigidity.py) (10,238): cols |
| `\NewColsB` | 2{,}508 | v5/out/rigidity_*.json (v5/rigidity.py) (12,210): cols |
| `\NewDistA` | $1.1\cdot10^{-14}$ | v5/out/rigidity_*.json (v5/rigidity.py) (10,238): dist_ls_upper |
| `\NewDistB` | $1.6\cdot10^{-14}$ | v5/out/rigidity_*.json (v5/rigidity.py) (12,210): dist_ls_upper |
| `\NewDistMax` | $1.6\cdot10^{-14}$ | v5/out/rigidity_*.json (v5/rigidity.py): max dist_ls_upper of (10,238), (12,210) |
| `\NewGapA` | $1.7\cdot10^{-7}$ | v5/out/rigidity_*.json (v5/rigidity.py) (10,238): -rattler_max_minus_c_upper |
| `\NewGapB` | $1.2\cdot10^{-1}$ | v5/out/rigidity_*.json (v5/rigidity.py) (12,210): -rattler_max_minus_c_upper |
| `\NewRankA` | 1{,}715 | v5/out/rigidity_*.json (v5/rigidity.py) (10,238): rank_mod_p |
| `\NewRankB` | 2{,}442 | v5/out/rigidity_*.json (v5/rigidity.py) (12,210): rank_mod_p |
| `\NonoptDistA` | $1.4\cdot10^{-8}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status NONOPT (3, 576): dist_ls |
| `\NonoptDistB` | $7.4\cdot10^{-5}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status NONOPT (3, 791): dist_ls |
| `\NonoptList` | $(3,576)$, $(3,791)$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status NONOPT |
| `\NonoptStressA` | $-1.1\cdot10^{-6}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status NONOPT (3, 576): smallest stress over mean |
| `\NonoptStressB` | $-6.2$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status NONOPT (3, 791): smallest stress over mean |
| `\OpenContacts` | 4{,}875 | v5/out/rigidity_*.json (v5/rigidity.py) (19,166) |
| `\OpenLPMaxMin` | $4.6\cdot10^{-17}$ | v5/out/rigidity_*.json (v5/rigidity.py) (19,166) |
| `\OpenRankNeeded` | 2{,}394 | v5/out/rigidity_*.json (v5/rigidity.py) (19,166): cols - g |
| `\OpenRankStressed` | 2{,}022 | v5/out/rigidity_*.json (v5/rigidity.py) (19,166) |
| `\OpenStressed` | 3{,}810 | v5/out/rigidity_*.json (v5/rigidity.py) (19,166) |
| `\OptMissingList` | $(22,891)$, $(23,552)$ | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' and index 'optimal' |
| `\PassDistMax` | $8.9\cdot10^{-7}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: max dist_ls |
| `\PassDistMed` | $9.0\cdot10^{-14}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: dist_ls, median |
| `\PctAnyVerdict` | 11.1 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\PctPropAttVerdict` | 82 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) and class |
| `\PctPropNever` | 84 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) |
| `\PctUAnyVerdict` | 56 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' restricted to unmarked entries |
| `\PctUPropAttVerdict` | 95 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) and class (unmarked) |
| `\QCErrMax` | $2.1\cdot10^{-30}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: max c_err |
| `\QKappaHalfMed` | $6.7\cdot10^{-10}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: kappa_point/2, median |
| `\QRUniqueMed` | $2.4\cdot10^{-5}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: r_unique, median |
| `\QRhoMed` | $1.9\cdot10^{-10}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: rho, median |
| `\RattlerMovedMax` | $1.01\cdot10^{-6}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: max dist_rattlers over moved (2 digits, rounded up) |
| `\SizeCap` | 3{,}100 | census/summary_v4.json macro SizeCap (census/summarize.py from the per-entry outputs) |
| `\SplitCoreA` | 72 | v5/out/rigidity_*.json (v5/rigidity.py) (8,72): core_size |
| `\SplitCoreB` | 240 | v5/out/rigidity_*.json (v5/rigidity.py) (9,330): core_size |
| `\SplitDistMax` | $1.3\cdot10^{-14}$ | v5/out/summary_v5.json (v5/summarize_v5.py) macro VSplitDistMax |
| `\SplitRankA` | 520 | v5/out/rigidity_*.json (v5/rigidity.py) (8,72): rank_mod_p |
| `\SplitRankB` | 2{,}088 | v5/out/rigidity_*.json (v5/rigidity.py) (9,330): rank_mod_p |
| `\StrictDistMax` | $5.7\cdot10^{-10}$ | census/strict_v4/*.json (strict.py): max dist_ls_upper (sqrt(n d) times dist_coord_upper where NaN) |
| `\StrictDistMed` | $5.5\cdot10^{-15}$ | census/strict_v4/*.json (strict.py): dist_ls_upper, median (order statistic) |
| `\StrictExactList` | $(16,288)$, $(21,162)$, $(22,275)$ | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\StrictNaNList` | $(7,182)$ | census/strict_v4/*.json (strict.py): dist_ls_upper NaN |
| `\StrictNewList` | $(10,238)$, $(12,210)$ | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\SymCPUMax` | 7 | phase2/out_open/*.json (phase2/symstrict.py): max cpu_seconds |
| `\SymDistMax` | $1.3\cdot10^{-11}$ | phase2/out_open/*.json (phase2/symstrict.py): max dist_ls_upper |
| `\SymExEntry` | $(7,359)$ | phase2/out_open/*.json (phase2/symstrict.py): the entry with the largest ratio of unknowns to reduced unknowns |
| `\SymExExtra` | 714 | phase2/out_open/*.json (phase2/symstrict.py): its extra contacts |
| `\SymExGroup` | 5{,}040 | phase2/out_open/*.json (phase2/symstrict.py): its group order |
| `\SymExRed` | 12 | phase2/out_open/*.json (phase2/symstrict.py): its reduced system |
| `\SymExU` | 1{,}772 | phase2/out_open/*.json (phase2/symstrict.py): its full system |
| `\SymFieldList` | $\mathbb Q(\sqrt{3})$, $\mathbb Q(\sqrt{5})$ | phase2/out_open/*.json (phase2/symstrict.py) |
| `\SymGammaRadMax` | $1.5\cdot10^{-95}$ | phase2/out_open/*.json (phase2/symstrict.py): max radius of gamma* |
| `\SymGroupMax` | 5{,}040 | phase2/out_open/*.json (phase2/symstrict.py): group_order_exact |
| `\SymGroupMin` | 2 | phase2/out_open/*.json (phase2/symstrict.py): group_order_exact |
| `\SymOtherStatuses` | 1 not LOCALMIN on re-run | phase2/out_open/*.json (phase2/symstrict.py): statuses of the other failures (for the text, not printed if empty) |
| `\SymRedMax` | 283 | phase2/out_open/*.json (phase2/symstrict.py): reduced_unknowns max |
| `\SymSquareNoModelList` | $(4,50)$, $(12,26)$ | phase2/out_open/*.json (phase2/symstrict.py) |
| `\SymStressMinRel` | $2.9\cdot10^{-4}$ | phase2/out_open/*.json (phase2/symstrict.py): min stress_min_rel |
| `\SymStrictList` | $(3,17)$, $(3,34)$, $(3,36)$, $(3,38)$, $(3,50)$, $(3,61)$, $(3,72)$, $(3,78)$, $(3,93)$, $(3,120)$, $(3,153)$, $(3,180)$, $(3,190)$, $(3,205)$, $(3,300)$, $(3,360)$, $(3,480)$, $(3,620)$, $(3,692)$, $(3,732)$, $(3,852)$, $(3,1020)$, $(4,17)$, $(4,40)$, $(4,65)$, $(4,72)$, $(5,25)$, $(5,70)$, $(6,49)$, $(6,73)$, $(6,74)$, $(6,82)$, $(6,83)$, $(6,88)$, $(6,89)$, $(6,92)$, $(6,104)$, $(6,108)$, $(6,123)$, $(7,82)$, $(7,127)$, $(7,147)$, $(7,148)$, $(7,359)$, $(8,43)$, $(8,46)$, $(8,47)$, $(9,35)$, $(10,44)$ | phase2/out_open/*.json (phase2/symstrict.py) |
| `\SymUMax` | 3{,}061 | phase2/out_open/*.json (phase2/symstrict.py): unknowns of the full system, max |
| `\TiltMax` | $4.8\cdot10^{-42}$ | v5/out/subspace_*.json: -core_max_minus_c, max |
| `\TiltMaxEntry` | $(16,117)$ | v5/out/subspace_*.json |
| `\TiltMin` | $6.3\cdot10^{-94}$ | v5/out/subspace_*.json: -core_max_minus_c, min |
| `\TiltMinEntry` | $(23,102)$ | v5/out/subspace_*.json |
| `\TubeMax` | $1.6\cdot10^{-14}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: max localmin_displacement |
| `\UAnyVerdict` | 1{,}705 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' restricted to unmarked entries |
| `\UFlex` | 1{,}091 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' (unmarked) |
| `\UPropAttVerdict` | 1{,}702 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) and class (unmarked) |
| `\UPropAttempted` | 1{,}786 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) (unmarked) |
| `\UPropNever` | 1{,}051 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) (unmarked) |
| `\UStrict` | 1{,}126 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' (unmarked) |
| `\ValChiMax` | $1.9\cdot10^{-28}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN of the validated hyperstatic entries: max tau_extra |
