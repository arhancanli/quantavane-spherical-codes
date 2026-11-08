# Every macro of v8_numbers.tex: value and source

Written by make_numbers_v8.py; do not edit.

| macro | value | source |
|---|---|---|
| `\AboveCapCPUh` | 4 | extrapolation of the fit CPUSlope to the skipped systems |
| `\AboveCapGB` | 6.2 | 8 AboveCapMaxU^2 bytes |
| `\AboveCapMaxU` | 27{,}901 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper): smallest skipped system, max over entries |
| `\AboveFlatDistMax` | $5.4\cdot10^{-14}$ | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): flat cores, max dist_ls_upper |
| `\AboveGanContacts` | 20{,}670 | above/gan510.json: contacts |
| `\AboveGanDegs` | 80 to 82 | above/gan510.json: contacts at a point, min and max |
| `\AboveGanDist` | $3.9\cdot10^{-14}$ | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): (10,510) dist_ls_upper |
| `\AboveGanMB` | 220 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): (10,510) peak resident memory, MB (sampled) |
| `\AboveGanNext` | $1/3$ | above/gan510.json: largest inner product below 1/2 |
| `\AboveGanRank` | 5{,}055 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): (10,510) rank mod p |
| `\AboveGanSec` | 19 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): (10,510) wall-clock seconds |
| `\AboveGanzhinovVerdict` | jammed (a strict local minimum with no rattlers) | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): (10,510) |
| `\AboveHours` | 8 | above/out and above/out_more: total wall-clock hours of run_above.py and above_more.py |
| `\AboveMaxCols` | 20{,}000 | run_above.py --max-cols default |
| `\AboveMaxColsDone` | 19{,}680 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): most columns of a completed rank test |
| `\AboveMaxSec` | 1768 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): largest wall-clock seconds |
| `\AboveNonStrictClause` |  | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): local minima, not strict (clause omitted when there are none) |
| `\AboveStrictDistMax` | $7.6\cdot10^{-14}$ | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): strict or not strict, max dist_ls_upper |
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
| `\CJKTList` | $(4,24)$, $(5,40)$, $(6,72)$, $(7,126)$, $(8,240)$, $(9,306)$ | kissing configurations of [CJKT11] that are census entries (all strict) |
| `\CPUHoursCensus` | 3.9 | census/summary_v4.json macro CPUHoursCensus (census/summarize.py from the per-entry outputs) |
| `\CPUHoursStrict` | 0.8 | census/summary_v4.json macro CPUHoursStrict (census/summarize.py from the per-entry outputs) |
| `\CPUMaxLoc` | 51 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: cpu_seconds max |
| `\CPUMaxPass` | 2 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: cpu_seconds max |
| `\CPUMedLoc` | 0.08 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: cpu_seconds median |
| `\CPUMedPass` | 0.07 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: cpu_seconds median |
| `\CPUSlope` | 2.7 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS/LOCALMIN with >=1000 unknowns: least-squares slope of log cpu on log U |
| `\CertShaCq` | c12bf1\ldots | census/results_v4/*.json (top-level record): code_sha256['cq.py'] (one value in all entries) |
| `\ChiMax` | $1.2\cdot10^{-27}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: max tau_extra |
| `\ClosePos` | 22 | recount/tolscan_close.jsonl (tolscan.py close): exact PASS/LOCALMIN with the next pairs added |
| `\ClosePosList` | $(3,275)$, $(3,296)$, $(3,334)$, $(3,347)$, $(3,359)$, $(3,406)$, $(3,492)$, $(3,532)$, $(3,584)$, $(3,627)$, $(3,688)$, $(3,698)$, $(3,758)$, $(3,956)$, $(4,139)$, $(4,218)$, $(6,84)$, $(7,86)$, $(10,74)$, $(19,75)$, $(25,56)$, $(29,61)$ | recount/tolscan_close.jsonl (tolscan.py close) |
| `\CloseRun` | 43 | recount/tolscan_close.jsonl (tolscan.py close): entries |
| `\CloseSentence` | 22 of them give a positive certificate for the larger contact graph as well (16 in dimensions 3 and 4, and $(6,84)$, $(7,86)$, $(10,74)$, $(19,75)$, $(25,56)$, $(29,61)$); for these the stored digits do not decide whether the configuration the table means has the extra contacts, since the certified configurations with and without them both lie within $5.2\cdot10^{-8}$ of the stored code. For the others the larger proposal is not certified | recount/tolscan_close.jsonl (tolscan.py close) |
| `\DiagGapRatioMed` | 1.01 | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): over the unseparated proposals, median of (largest cosine - next pair) / tau |
| `\DiagRowsumMed` | $1.1\cdot10^{10}$ | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): median |
| `\DiagRowsumMin` | $1.4$ | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted) (census record krawczyk_attempts[0].max_rowsum): min |
| `\DiagSpreadRatioMed` | 0.99 | recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted): over the unseparated proposals, median of (largest - smallest contact cosine) / tau |
| `\DriverShaA` | 4a7435\ldots | census/results_v4/*.json (top-level record): code_sha256['run_census.py'], the more frequent |
| `\DriverShaB` | de9662\ldots | census/results_v4/*.json (top-level record): code_sha256['run_census.py'], the other |
| `\ExBBetaK` | $1.3\cdot10^{-14}$ | example/example_d5_n16.json (example.py 5 16): beta_K |
| `\ExBChi` | $8.0\cdot10^{-30}$ | example/example_d5_n16.json (example.py 5 16): chi |
| `\ExBDist` | $1.7\cdot10^{-29}$ | example/example_d5_n16.json (example.py 5 16): dist_ls |
| `\ExBExtra` | 25 | example/example_d5_n16.json (example.py 5 16) |
| `\ExBIso` | 55 | example/example_d5_n16.json (example.py 5 16) |
| `\ExBKappa` | $0.00122$ | example/example_d5_n16.json (example.py 5 16): kappa |
| `\ExBNContacts` | 80 | example/example_d5_n16.json (example.py 5 16) |
| `\ExBRho` | $1.14\cdot10^{-4}$ | example/example_d5_n16.json (example.py 5 16): rho |
| `\ExBStressRel` | $0.999$ | example/example_d5_n16.json (example.py 5 16): smallest relative stress |
| `\ExBTube` | $2.4\cdot10^{-25}$ | example/example_d5_n16.json (example.py 5 16): tube |
| `\ExBU` | 81 | example/example_d5_n16.json (example.py 5 16) |
| `\ExBetaK` | $1.4\cdot10^{-15}$ | example/example_d3_n7.json (example.py 3 7): beta_K |
| `\ExCDigits` | 21013\,83127\,30603\,08486\,57653\,0161 | example/example_d3_n7.json (example.py 3 7): common digits of the 30-digit interval, in groups of five |
| `\ExCPU` | 0.6 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): (3,7) cpu_seconds |
| `\ExDegrees` | three points of degree 4 and four points of degree 3 | example/example_d3_n7.json (example.py 3 7): degrees |
| `\ExDist` | $6.8\cdot10^{-17}$ | example/example_d3_n7.json (example.py 3 7): dist_ls |
| `\ExEpsGamma` | $1.4\cdot10^{-31}$ | example/example_d3_n7.json (example.py 3 7): eps_gamma |
| `\ExEpsMax` | $3.8\cdot10^{-31}$ | example/example_d3_n7.json (example.py 3 7): max eps |
| `\ExGap` | $0.643$ | example/example_d3_n7.json (example.py 3 7): gap_lower |
| `\ExKappa` | $0.0176$ | example/example_d3_n7.json (example.py 3 7): kappa |
| `\ExKappaHalf` | $0.00510$ | example/example_d3_n7.json (example.py 3 7): kappa/(2 sqrt d) |
| `\ExLambda` | $1.22$ | example/example_d3_n7.json (example.py 3 7): Lambda |
| `\ExMuP` | $3.11$ | example/example_d3_n7.json (example.py 3 7): mu_p |
| `\ExNContacts` | 12 | example/example_d3_n7.json (example.py 3 7) |
| `\ExNextGap` | $0.643$ | example/example_d3_n7.json (example.py 3 7): next pair below the largest cosine (double) |
| `\ExNormC` | $3.49$ | example/example_d3_n7.json (example.py 3 7): ||C||_inf |
| `\ExOnlyDistCoordMax` | $1.3\cdot10^{-15}$ | v5/out/summary_v5.json (v5/summarize_v5.py) macro VExOnlyDistCoordMax |
| `\ExOnlyDistLSMax` | $1.2\cdot10^{-14}$ | v5/out/summary_v5.json (v5/summarize_v5.py) macro VExOnlyDistLSMax |
| `\ExOnlyN` | 18 | v5/out/summary_v5.json (v5/summarize_v5.py) macro VExOnlyN |
| `\ExR` | $7.9\cdot10^{-31}$ | example/example_d3_n7.json (example.py 3 7): r |
| `\ExRho` | $0.00210$ | example/example_d3_n7.json (example.py 3 7): rho |
| `\ExRu` | $0.0262$ | example/example_d3_n7.json (example.py 3 7): r_unique |
| `\ExSpan` | $4.4\cdot10^{-30}$ | example/example_d3_n7.json (example.py 3 7): d r ||B^-1|| |
| `\ExStressMin` | $0.0519$ | example/example_d3_n7.json (example.py 3 7): smallest stress, lower bound |
| `\ExStressRel` | $0.623$ | example/example_d3_n7.json (example.py 3 7): smallest stress times |E| |
| `\ExTableCos` | $0.210138312731$ | example/example_d3_n7.json (example.py 3 7): table cosine |
| `\ExTheta` | $6.1\cdot10^{-15}$ | example/example_d3_n7.json (example.py 3 7): vartheta |
| `\ExU` | 22 | example/example_d3_n7.json (example.py 3 7) |
| `\ExZeta` | $0.322$ | example/example_d3_n7.json (example.py 3 7): zeta |
| `\ExactCPUMax` | 679 | census/exact_v4/strict/*.json (exact_strict.py): max cpu_seconds of STRICT |
| `\FXEDenseOps` | $2.2\cdot10^{14}$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): a dense product: sum rows^2 columns |
| `\FXEHours` | 2.3 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): cpu hours, all runs |
| `\FXELimbOps` | $4.7\cdot10^{12}$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): upper bound of the limb products, 4 per scalar product, 2d nonzeros per row of L: sum 8 d rows^2 |
| `\FXEMatchGramSeven` | 12{,}482 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): same E_sha256 as flexgram7 |
| `\FXEMaxCols` | 29{,}280 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): max columns |
| `\FXEMaxRows` | 3{,}100 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): max rows |
| `\FXEMaxSec` | 7 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): max cpu seconds |
| `\FXEQR` | 130 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): QR right inverses |
| `\FXERatioMax` | $1.641$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): max of that ratio |
| `\FXERatioMed` | $1.000$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): median of that ratio |
| `\FXERatioMin` | $0.213$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): min exact delta / flexgram7 (double kernel) delta |
| `\FXERefineList` | $(17,317)$ | recount/out_flex8r: the entries tried by flexrefine8.py |
| `\FXERefineRun` | 1 | recount/out_flex8r: entries where flexexact8.py failed, tried by flexrefine8.py |
| `\FXESVD` | 1 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): SVD right inverses (flexrefine8.py) |
| `\FamilyList` | $(22,47)$, $(22,48)$, $(23,50)$ | v5/out/summary_v5.json (v5/summarize_v5.py) macro VFamilyList |
| `\FlatColoursMax` | 52 | v5/out/summary_v5.json (v5/summarize_v5.py) macro VFlatColoursMax |
| `\FlatColoursMin` | 26 | v5/out/summary_v5.json (v5/summarize_v5.py) macro VFlatColoursMin |
| `\FlatDistAllMax` | $5.4\cdot10^{-14}$ | max of FlatDistMax (v5, census flat cores) and AboveFlatDistMax |
| `\FlatDistMax` | $2.8\cdot10^{-14}$ | v5/out/summary_v5.json (v5/summarize_v5.py) macro VFlatDistMax |
| `\FlatList` | $(16,117)$, $(17,162)$, $(20,92)$, $(22,79)$, $(22,114)$, $(22,127)$, $(23,102)$, $(23,129)$, $(23,147)$, $(24,121)$, $(24,149)$, $(26,57)$, $(27,96)$ | v5/out/summary_v5.json (v5/summarize_v5.py) macro VFlatList |
| `\FlatLiveImpRange` | $7.8\cdot10^{-11}$ to $1.1\cdot10^{-5}$ | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): their improvements |
| `\FlatOtherList` | $(16,117)$, $(17,162)$, $(20,92)$, $(22,127)$, $(23,147)$, $(24,121)$ | the other flat cores |
| `\FlatUnionList` | $(22,79)$, $(22,114)$, $(23,102)$, $(23,129)$, $(24,149)$, $(26,57)$, $(27,96)$ | census/index.tsv notes mention an orthogonal union (flat cores) |
| `\FlexCatAboveCapMed` | $2.5\cdot10^{-13}$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): category above-cap, median delta |
| `\FlexCatAboveCapN` | 117 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): category above-cap |
| `\FlexCatAfterFailedMed` | $2.6\cdot10^{-11}$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): category after-failed, median delta |
| `\FlexCatAfterFailedN` | 302 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): category after-failed |
| `\FlexCatNoProposalFlaggedMed` | $2.5\cdot10^{-9}$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): category no-proposal-flagged, median delta |
| `\FlexCatNoProposalFlaggedN` | 11{,}019 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): category no-proposal-flagged |
| `\FlexCatNoProposalUnflaggedMed` | $2.8\cdot10^{-8}$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): category no-proposal-unflagged, median delta |
| `\FlexCatNoProposalUnflaggedN` | 1{,}044 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): category no-proposal-unflagged |
| `\FlexDeltaMax` | $2.7\cdot10^{-5}$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): max delta |
| `\FlexDeltaMed` | $2.5\cdot10^{-9}$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): median delta |
| `\FlexDeltaMedRes` | $2.5\cdot10^{-9}$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): median delta, resolved |
| `\FlexDeltaMin` | $1.2\cdot10^{-15}$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): min delta |
| `\FlexDeltaMinEntry` | $(11,105)$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): entry of the min delta |
| `\FlexDeltaMinH` | $4.4\cdot10^{-15}$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): its h |
| `\FlexEpsMax` | $6.4\cdot10^{-15}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram7_*.jsonl (recount/flexgram7.py): max of the Gram error bound eps |
| `\FlexHLimitedList` | $(3,314)$, $(11,105)$, $(26,81)$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): the three h-limited entries with the smallest delta |
| `\FlexHMed` | $1.2\cdot10^{-2}$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): h/sqrt(d), median |
| `\FlexHMinEntry` | $(11,105)$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): entry of the smallest h |
| `\FlexHMinH` | $4.4\cdot10^{-15}$ | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): its h |
| `\FlexTauRatioMin` | 1.00 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): min tau'/tol |
| `\FlowAbovecapFlex` | 50 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) x class |
| `\FlowAbovecapNone` | 214 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) x class |
| `\FlowAbovecapVerdict` | 127 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) x class |
| `\FlowAttemptedFlex` | 202 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) x class |
| `\FlowAttemptedNone` | 179 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) x class |
| `\FlowAttemptedVerdict` | 1{,}795 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) x class |
| `\FlowNeverFlex` | 11{,}337 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) x class |
| `\FlowNeverNone` | 2{,}361 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) x class |
| `\GapCertDistMax` | -- | recount/gapcut_cert.jsonl (gapcut.py certify): max dist_ls |
| `\GapCertExistFail` | 50 | recount/gapcut_cert.jsonl (gapcut.py certify): existence test failed |
| `\GapCertList` | none | recount/gapcut_cert.jsonl (gapcut.py certify) |
| `\GapCertMaxU` | 1{,}500 | gapsample.py: at most 1,500 unknowns |
| `\GapCertPos` | 0 | recount/gapcut_cert.jsonl (gapcut.py certify): exact PASS or LOCALMIN |
| `\GapCertRun` | 54 | recount/gapcut_cert.jsonl (gapcut.py certify): entries |
| `\GapFound` | 1{,}336 | recount/gapcut_*.jsonl (gapcut.py scan): a non-hypostatic cut found |
| `\GapFoundSmall` | 224 | recount/gapcut_*.jsonl (gapcut.py scan): with spread <= 1e-7 |
| `\GapScanN` | 12{,}583 | recount/gapcut_*.jsonl (gapcut.py scan): entries scanned |
| `\ImpLowerA` | $6.4\cdot10^{-12}$ | v5/out/summary_v5.json (v5/summarize_v5.py) macro VImpLowerA |
| `\ImpLowerB` | $1.8\cdot10^{-13}$ | v5/out/summary_v5.json (v5/summarize_v5.py) macro VImpLowerB |
| `\InvShareMed` | 54 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS/LOCALMIN >=1000 unknowns: seconds_inverse/wall_seconds, median (wall-clock share) |
| `\LiveCertBeaten` | 597 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): of these, beaten |
| `\LiveCertListed` | 1{,}836 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): certified local minima still listed |
| `\LiveCertNotBestPct` | 33 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): of these, beaten, percent |
| `\LiveCertStillBest` | 1{,}239 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): of these, not beaten (no listed cosine below the certified value) |
| `\LiveDate` | 7 November 2026 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026) |
| `\LiveDistLowMed` | $3.4\cdot10^{-3}$ | live/out/*.json (run_live.py on codes fetched by fetch_live.py): positive-class, median dist_lower_point |
| `\LiveDistLowMedBeaten` | $3.9\cdot10^{-3}$ | live/out/*.json (run_live.py on codes fetched by fetch_live.py): beaten ones, median dist_lower_point |
| `\LiveDistLowMin` | $5.5\cdot10^{-14}$ | live/out/*.json (run_live.py on codes fetched by fetch_live.py): positive-class, min dist_lower_point |
| `\LiveDistLowMinBeaten` | $5.2\cdot10^{-5}$ | live/out/*.json (run_live.py on codes fetched by fetch_live.py): beaten ones, min dist_lower_point |
| `\LiveFetchFlex` | 175 | live/out/*.json (run_live.py on codes fetched by fetch_live.py): sampled flex entries |
| `\LiveFetchPos` | 640 | live/out/*.json (run_live.py on codes fetched by fetch_live.py): positive-class entries |
| `\LiveFlexImpMed` | $2.3\cdot10^{-4}$ | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): flex entries, median improvement |
| `\LiveFlexImproved` | 9{,}159 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): flex entries listed with a lower cosine |
| `\LiveFlexListed` | 9{,}166 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): flex entries listed |
| `\LiveFlexPct` | 99.9 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): flex entries listed with a lower cosine, percent |
| `\LiveListed` | 13{,}277 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): census entries still listed |
| `\LiveMarkedFlexPct` | 99.9 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): flex, marked and listed: lower cosine, percent |
| `\LiveMarkedLocalPct` | 73.8 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): local, marked and listed: lower cosine, percent |
| `\LiveMarkedNonePct` | 99.0 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): none, marked and listed: lower cosine, percent |
| `\LiveMarkedStrictPct` | 80.8 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): strict, marked and listed: lower cosine, percent |
| `\LiveMarkedUnresPct` | 98.7 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): unres, marked and listed: lower cosine, percent |
| `\LiveRerunBeatenFlex` | 27 | live/out/*.json (run_live.py on codes fetched by fetch_live.py): beaten ones whose new code gets FLEX |
| `\LiveRerunBeatenN` | 596 | live/out/*.json (run_live.py on codes fetched by fetch_live.py): fetched positive-class entries beaten (live cosine below the certified value) |
| `\LiveRerunBeatenPos` | 534 | live/out/*.json (run_live.py on codes fetched by fetch_live.py): beaten ones whose new code gets an exact PASS/LOCALMIN |
| `\LiveRerunFlexFlex` | 64 | live/out/*.json (run_live.py on codes fetched by fetch_live.py): flex sample, FLEX |
| `\LiveRerunFlexNone` | 48 | live/out/*.json (run_live.py on codes fetched by fetch_live.py): flex sample, no verdict |
| `\LiveRerunFlexPos` | 63 | live/out/*.json (run_live.py on codes fetched by fetch_live.py): flex sample, exact PASS/LOCALMIN |
| `\LiveRerunFlexPosPct` | 36 | live/out/*.json (run_live.py on codes fetched by fetch_live.py): flex sample, positive, percent |
| `\LiveRerunPosFlex` | 30 | live/out/*.json (run_live.py on codes fetched by fetch_live.py): positive-class, FLEX |
| `\LiveRerunPosOther` | 37 | live/out/*.json (run_live.py on codes fetched by fetch_live.py): positive-class, other |
| `\LiveRerunPosPos` | 573 | live/out/*.json (run_live.py on codes fetched by fetch_live.py): positive-class, exact PASS/LOCALMIN |
| `\LiveRerunPosStrictIso` | 381 | live/out/*.json (run_live.py on codes fetched by fetch_live.py): positive-class, exact PASS |
| `\LiveRerunUnconfirmed` | 0 | live/out/*.json (run_live.py on codes fetched by fetch_live.py): first pass positive, exact not |
| `\LiveStrictAuthor` | 4 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): beaten strict entries credited to the author |
| `\LiveStrictBeaten` | 455 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): strict entries listed with a cosine below the certified value |
| `\LiveStrictImpMed` | $8.7\cdot10^{-6}$ | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): beaten strict entries, median improvement of the displayed cosine |
| `\LiveStrictListed` | 1{,}214 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): strict entries listed |
| `\LiveStrictLower` | 463 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): strict entries listed with a lower displayed cosine |
| `\LiveStrictPct` | 37 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): strict entries beaten, percent of listed |
| `\LiveStrictRefined` | 8 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): strict entries with a lower displayed cosine not below the certified value |
| `\LiveUnmarkedFlexImp` | 881 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): flex, unmarked, listed, lower cosine |
| `\LiveUnmarkedFlexN` | 882 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): flex, unmarked and listed |
| `\LiveUnmarkedFlexPct` | 99.9 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): flex, unmarked and listed: lower cosine, percent |
| `\LiveUnmarkedLocalImp` | 124 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): local, unmarked, listed, lower cosine |
| `\LiveUnmarkedLocalN` | 555 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): local, unmarked and listed |
| `\LiveUnmarkedLocalPct` | 22.3 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): local, unmarked and listed: lower cosine, percent |
| `\LiveUnmarkedNoneImp` | 35 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): none, unmarked, listed, lower cosine |
| `\LiveUnmarkedNoneN` | 91 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): none, unmarked and listed |
| `\LiveUnmarkedNonePct` | 38.5 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): none, unmarked and listed: lower cosine, percent |
| `\LiveUnmarkedStrictImp` | 442 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): strict, unmarked, listed, lower cosine |
| `\LiveUnmarkedStrictN` | 1{,}188 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): strict, unmarked and listed |
| `\LiveUnmarkedStrictPct` | 37.2 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): strict, unmarked and listed: lower cosine, percent |
| `\LiveUnmarkedUnresImp` | 71 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): unres, unmarked, listed, lower cosine |
| `\LiveUnmarkedUnresN` | 82 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): unres, unmarked and listed |
| `\LiveUnmarkedUnresPct` | 86.6 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): unres, unmarked and listed: lower cosine, percent |
| `\LocDistMax` | $4.0\cdot10^{-6}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: max dist_ls (to P') |
| `\LocDistTubeMax` | $4.0\cdot10^{-6}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: max dist_ls + localmin_displacement |
| `\LocExtraMax` | 6{,}227 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: max |
| `\LocExtraMed` | 62 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: median |
| `\LocExtraMin` | 1 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: min contacts - isostatic_count |
| `\LowCondFamily` | $N=2d+1$, from $(5,11)$ to $(32,65)$ | census/results_v4: PASS with cond_inf_upper < 600 and at least 50 unknowns, all but two with N = 2d+1 |
| `\MAnyVerdict` | 94 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' (marked) |
| `\MPropNever` | 12{,}647 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) (marked) |
| `\MaxUnknownsLoc` | 3{,}061 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: max unknowns |
| `\MaxUnknownsPass` | 2{,}974 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: max unknowns |
| `\NAboveAll` | 155 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): entries |
| `\NAboveCapExact` | 89 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) |
| `\NAboveCapOpen` | 264 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) |
| `\NAboveCore` | 5 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): core |
| `\NAboveCoreFail` | 5 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): exact core not built |
| `\NAboveDecided` | 124 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): decided |
| `\NAboveDeg` | 2 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): degree |
| `\NAboveField` | 3 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): rational cosine, core recognised over Q(sqrt D) |
| `\NAboveFlat` | 38 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): flat |
| `\NAboveFlatUnion` | 25 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py) and census/index.tsv notes: flat cores that the table describes as orthogonal unions |
| `\NAboveLPNone` | 2 | above/out_more: bar-rigid, and the linear program finds no strictly positive stress (untrusted) |
| `\NAboveMoreRun` | 19 | above/out_more: entries of the second pass |
| `\NAboveMoreStopped` | 11 | above/out_more: second-pass entries stopped by the time or memory limit |
| `\NAboveNoStress` | 16 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): nostress |
| `\NAboveNoStressRest` | 15 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): no stress found, the others |
| `\NAboveNoStressZero` | 1 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): no stress found, and the stress closest to the uniform one vanishes on all contacts of some core point |
| `\NAboveNonStrict` | 0 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): nonstrict |
| `\NAboveNotRigid` | 3 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): notrigid |
| `\NAboveNotRigidList` | $(11,604)$, $(22,166)$, $(30,656)$ | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): rank short |
| `\NAboveNotRigidStress` | 3 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): rank short, with an exact strictly positive stress |
| `\NAboveOther` | 0 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): other |
| `\NAboveSecondPass` | 3 | above/out_more (above_more.py): entries decided by the second pass |
| `\NAboveStrict` | 86 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): strict |
| `\NAboveStrictSub` | 0 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): strict with (iii), (iv) on the support of a stress |
| `\NAboveTooLarge` | 5 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): toolarge |
| `\NAboveTried` | 153 | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): cosine of degree <= 2 |
| `\NAnyVerdict` | 1{,}922 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
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
| `\NFlatLiveImproved` | 4 | v5/live/live_20261107.tsv (the live table's index read on 7 November 2026): flat-core entries with a lower cosine now |
| `\NFlex` | 11{,}589 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NFlexCensus` | 12{,}482 | census/results_v4: entries with the negative certificate (FLEX) in the census |
| `\NFlexDeltaAboveNine` | 8{,}065 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): delta >= 1e-9 |
| `\NFlexDeltaAtTauQuarter` | 12{,}374 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): delta set by tau' |
| `\NFlexDeltaBelowNine` | 4{,}417 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): delta < 1e-9 |
| `\NFlexDeltaBelowTwelve` | 797 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): delta < 1e-12 |
| `\NFlexExactFail` | 0 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NFlexFlagged` | 11{,}391 | census/index.tsv: FLEX entries marked |
| `\NFlexHLimited` | 108 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): delta set by h |
| `\NFlexResolved` | 11{,}589 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NFlexTauBelowTau` | 55 | recount/out_flex8/*.json (flexexact8.py, the exact decision) and recount/out_flex8r/*.json (flexrefine8.py): tau' < tol |
| `\NFlexTheorem` | 12{,}482 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram7_*.jsonl (recount/flexgram7.py): count |
| `\NFlexTolEight` | 3{,}003 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram7_*.jsonl (recount/flexgram7.py): tol 1e-08 |
| `\NFlexTolNine` | 2{,}071 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram7_*.jsonl (recount/flexgram7.py): tol 1e-09 |
| `\NFlexTolSeven` | 5{,}121 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram7_*.jsonl (recount/flexgram7.py): tol 1e-07 |
| `\NFlexTolSmaller` | 1{,}353 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram7_*.jsonl (recount/flexgram7.py): tol < 1e-10 |
| `\NFlexTolTen` | 934 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram7_*.jsonl (recount/flexgram7.py): tol 1e-10 |
| `\NFlexUnmarked` | 1{,}091 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status FLEX (cq.certify_flex) and recount/flexgram7_*.jsonl (recount/flexgram7.py): unmarked entries |
| `\NFlexUnresolved` | 893 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NHyperTen` | 1{,}364 | census/summary_v4.json NClassHyperTen |
| `\NHypoEmptyTen` | 14{,}162 | census/summary_v4.json NClassEmptyTen + NClassHypoTen |
| `\NInside` | 668 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: code_inside_rho |
| `\NIsoTen` | 739 | census/summary_v4.json NClassIsoTen |
| `\NListed` | 16{,}269 | census/index.tsv (snapshot.py from the archived data set): codes listed in table.txt |
| `\NLiveAuthor` | 443 | v5/live/live_20261107.tsv (the live table read on 7 November 2026): credited to the author |
| `\NLocTolOther` | 31 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: tol != 1e-10 |
| `\NLocal` | 630 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NLocalAll` | 937 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN |
| `\NLocalAtLeast` | 1{,}869 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NLocalVFive` | 679 | v5/out/summary_v5.json VNLocal |
| `\NNoRho` | 13 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: rho absent |
| `\NNoVerdict` | 1{,}861 | NFailedAll + NTooLarge |
| `\NNoVerdictAll` | 2{,}754 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NNoVerdictOther` | 1{,}861 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NNonStrict` | 2 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NNoneTried` | 12{,}564 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NNonoptThenCert` | 0 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status an attempt NONOPT and a later PASS/LOCALMIN |
| `\NNotListed` | 63 | census/index.tsv (snapshot.py from the archived data set): codes not listed |
| `\NNotLocal` | 53 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NNotLocalFlat` | 13 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NNotLocalFlatAll` | 51 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NNotLocalNonopt` | 2 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NOpenDegree` | 29 | census/strict_v4/*.json (strict.py) status of the entries still open |
| `\NOpenExact` | 22 | v5/out/summary_v5.json (v5/summarize_v5.py) macro VNOpenExact |
| `\NOpenNoPoly` | 509 | census/strict_v4/*.json (strict.py) status of the entries still open |
| `\NOpenNotInField` | 49 | census/strict_v4/*.json (strict.py) status of the entries still open |
| `\NOpenTooMany` | 43 | census/strict_v4/*.json (strict.py) status of the entries still open |
| `\NOptAll` | 25 | census/index.tsv (snapshot.py from the archived data set): census entries marked 'Provably optimal' |
| `\NOptStrict` | 25 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' and index 'optimal' |
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
| `\NStrict` | 1{,}237 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\NStrictAbove` | 86 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
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
| `\NTooLarge` | 1{,}779 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
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
| `\OptAboveText` | they are above the size limit of the floating-point pass, and the exact route of Section~\ref{sec:above} certifies both strict | above/out/*.json (run_above.py) and above/out_more/*.json (above_more.py): (22,891) and (23,552) |
| `\OptMissingList` |  | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' and index 'optimal' |
| `\OptSentence` | all of them are certified strict (Section~\ref{sec:valid}) | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' and index 'optimal' |
| `\PassDistMax` | $8.9\cdot10^{-7}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: max dist_ls |
| `\PassDistMed` | $9.0\cdot10^{-14}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: dist_ls, median |
| `\PctAnyVerdict` | 11.8 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' |
| `\PctPropAttVerdict` | 82 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) and class |
| `\PctPropNever` | 84 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) |
| `\PctUAnyVerdict` | 61 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' restricted to unmarked entries |
| `\PctUPropAttVerdict` | 95 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) and class (unmarked) |
| `\QCErrMax` | $2.1\cdot10^{-30}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: max c_err |
| `\QKappaHalfMed` | $6.7\cdot10^{-10}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: kappa_point/2, median |
| `\QRUniqueMed` | $2.4\cdot10^{-5}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: r_unique, median |
| `\QRhoMed` | $1.9\cdot10^{-10}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: rho, median |
| `\RattlerMovedMax` | $1.01\cdot10^{-6}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status PASS: max dist_rattlers over moved (2 digits, rounded up) |
| `\SecondFlexGood` | 8 | census/exact_v4/big/flex_*.json (exact_big.py --flex): FLEX confirmed |
| `\SecondFlexN` | 8 | census/exact_v4/big/flex_*.json (exact_big.py --flex) |
| `\SecondFlexRowsMax` | 771 | census/exact_v4/big/flex_*.json (exact_big.py --flex): max rows |
| `\SecondFlexRowsMin` | 300 | census/exact_v4/big/flex_*.json (exact_big.py --flex): min rows |
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
| `\StructIsoPct` | 45 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided), above/out (core_size, contacts), census/index.tsv: isostatic share of the local minima |
| `\StructNLarge` | 11{,}102 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided), above/out (core_size, contacts), census/index.tsv: entries with 16 < N/d <= 1000000000 |
| `\StructNMid` | 4{,}288 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided), above/out (core_size, contacts), census/index.tsv: entries with 4 < N/d <= 16 |
| `\StructNSmall` | 875 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided), above/out (core_size, contacts), census/index.tsv: entries with 0 < N/d <= 4 |
| `\StructNdNeg` | 38.0 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided), above/out (core_size, contacts), census/index.tsv: median N/d of the negative entries |
| `\StructNdPos` | 13.5 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided), above/out (core_size, contacts), census/index.tsv: median N/d of the local minima |
| `\StructRatPct` | 62 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided), above/out (core_size, contacts), census/index.tsv: local minima with rattlers, share |
| `\StructShareLarge` | 8 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided), above/out (core_size, contacts), census/index.tsv: share of certified local minima for 16 < N/d <= 1000000000 |
| `\StructShareMid` | 14 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided), above/out (core_size, contacts), census/index.tsv: share of certified local minima for 4 < N/d <= 16 |
| `\StructShareSmall` | 45 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided), above/out (core_size, contacts), census/index.tsv: share of certified local minima for 0 < N/d <= 4 |
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
| `\TolConflictSentence` | no entry has positive certificates of two different kinds | recount/tolscan_other_*.jsonl (tolscan.py other): entries with positive certificates of two kinds |
| `\TolOtherConflict` | 0 | recount/tolscan_other_*.jsonl (tolscan.py other): entries with positive certificates of two kinds |
| `\TolOtherEntries` | 146 | recount/tolscan_other_*.jsonl (tolscan.py other): entries |
| `\TolOtherNone` | 230 | recount/tolscan_other_*.jsonl (tolscan.py other): other outcomes |
| `\TolOtherPos` | 0 | recount/tolscan_other_*.jsonl (tolscan.py other): exact PASS/LOCALMIN at another tolerance |
| `\TolOtherProposals` | 230 | recount/tolscan_other_*.jsonl (tolscan.py other): proposals |
| `\TolOtherSamePos` | 0 | recount/tolscan_other_*.jsonl (tolscan.py other): of the same kind |
| `\TolOtherSentence` | Running the certificates on all of them gives no further positive certificate: all 230 attempts fail | recount/tolscan_other_*.jsonl (tolscan.py other) |
| `\TubeMax` | $1.6\cdot10^{-14}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN: max localmin_displacement |
| `\UAnyVerdict` | 1{,}828 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' restricted to unmarked entries |
| `\UFlex` | 1{,}006 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' (unmarked, resolved) |
| `\UPropAttVerdict` | 1{,}702 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) and class (unmarked) |
| `\UPropAttempted` | 1{,}786 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) (unmarked) |
| `\UPropNever` | 1{,}051 | census/results_v4/*.json (cq.py via run_census.py), the attempt with status (attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper) (unmarked) |
| `\UStrict` | 1{,}211 | class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict' (unmarked) |
| `\ValChiMax` | $1.9\cdot10^{-28}$ | census/results_v4/*.json (cq.py via run_census.py), the attempt with status LOCALMIN of the validated hyperstatic entries: max tau_extra |
| `\XCErrMax` | $9\cdot10^{-30}$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): widest 30-digit interval |
| `\XCPUSlope` | 2.8 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): least-squares slope of log cpu_seconds on log unknowns, unknowns >= 1000 |
| `\XCSame` | 1{,}776 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided) vs census/results_v4: identical 30-digit intervals (all others overlap; asserted) |
| `\XChiMax` | $3.1\cdot10^{-27}$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): LOCALMIN max chi |
| `\XDistRatioMax` | $1.000$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): max of exact/census dist_ls |
| `\XDistRatioMed` | $0.166$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): median of exact/census dist_ls |
| `\XDistTinyN` | 21 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): certificates with dist_ls < 1e-19 (drawn on the lower edge of Figure 3) |
| `\XGapHypEleven` | 0 | census/results_v4: gaps.min_gap_float of the deciding attempt (first pass, rattlers included): LOCALMIN < 1e-11 |
| `\XGapHypNine` | 13 | census/results_v4: gaps.min_gap_float of the deciding attempt (first pass, rattlers included): LOCALMIN < 1e-9 |
| `\XGapIsoEleven` | 3 | census/results_v4: gaps.min_gap_float of the deciding attempt (first pass, rattlers included): PASS < 1e-11 |
| `\XGapIsoNine` | 30 | census/results_v4: gaps.min_gap_float of the deciding attempt (first pass, rattlers included): PASS < 1e-9 |
| `\XGapMin` | $2.27\cdot10^{-13}$ | census/results_v4: gaps.min_gap_float of the deciding attempt (first pass, rattlers included): smallest |
| `\XGapMinEntry` | $(15,63)$ | census/results_v4: gaps.min_gap_float of the deciding attempt (first pass, rattlers included): its entry |
| `\XGapMinExact` | $2.32\cdot10^{-13}$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): its exact gap_lower (non-rattlers) |
| `\XGapMinTol` | $10^{-13}$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): its tolerance |
| `\XHours` | 2.7 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): sum of cpu_seconds / 3600 |
| `\XInsideHyp` | 909 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): LOCALMIN, code_inside_rho |
| `\XInsideIso` | 668 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): PASS, code_inside_rho |
| `\XKappaHalfMed` | $6.7\cdot10^{-10}$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): PASS kappa_half median |
| `\XKappaRatioMax` | $1.000$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): max of exact/census kappa |
| `\XKappaRatioMin` | $0.720$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): min of exact/census kappa |
| `\XLocDistTubeMax` | $4.0\cdot10^{-6}$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): open LOCALMIN max dist_ls_plus_tube |
| `\XMaxU` | 3{,}061 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): max unknowns |
| `\XNChecked` | 1{,}778 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): files |
| `\XNDisagree` | 0 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): exact status differs from the census status ([]) |
| `\XNLocOK` | 937 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): census LOCALMIN, exact LOCALMIN |
| `\XNNonoptOK` | 2 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): census NONOPT, exact NONOPT |
| `\XNPassDistAbove` | 97 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): PASS dist_ls > 1e-8 |
| `\XNPassOK` | 839 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): census PASS, exact PASS |
| `\XNRhoDefined` | 826 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): PASS with rho |
| `\XNoRhoIso` | 13 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): PASS without rho |
| `\XNonoptDistMax` | $7.4\cdot10^{-5}$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): NONOPT max dist_ls |
| `\XOutsideIso` | 158 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): PASS, not code_inside_rho |
| `\XPassDistMax` | $8.9\cdot10^{-7}$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): PASS max dist_ls |
| `\XRUniqueMed` | $2.4\cdot10^{-5}$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): PASS r_unique median |
| `\XRepeatList` | $(19,102)$ | strictness re-checks repeated with version 2 after the untrusted LP hit its time limit (stress_lp_status 1) |
| `\XRepeatN` | 1 | their number |
| `\XRhoMed` | $1.9\cdot10^{-10}$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): PASS rho median |
| `\XSecMax` | 89 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): cpu_seconds, max |
| `\XSecMaxEntry` | $(18,163)$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): entry of the max cpu_seconds |
| `\XSecMaxU` | 2{,}917 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): its unknowns |
| `\XSecMed` | 1.2 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): cpu_seconds, median |
| `\XStageTwoList` | $(3,613)$, $(3,791)$, $(3,877)$ | their entries |
| `\XStageTwoN` | 3 | certificates whose exact stress test needed the second enclosure (stress_enclosure 'both') |
| `\XStressRatioMax` | $1.000$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): max of exact/census stress_min_lower |
| `\XStressRatioMin` | $0.720$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): min of exact/census stress_min_lower |
| `\XStrictHours` | 1.5 | verify/out_strict7/*.json (strict7.py: strict.py with exact (S') and span) and verify/out_sym7/*.json (symstrict7.py: symstrict.py with exact (S') and span): cpu_seconds / 3600 |
| `\XStrictOK` | 258 | verify/out_strict7/*.json (strict7.py: strict.py with exact (S') and span): STRICT |
| `\XStrictRun` | 258 | verify/out_strict7/*.json (strict7.py: strict.py with exact (S') and span): files |
| `\XSymOK` | 49 | verify/out_sym7/*.json (symstrict7.py: symstrict.py with exact (S') and span): STRICT |
| `\XSymRun` | 49 | verify/out_sym7/*.json (symstrict7.py: symstrict.py with exact (S') and span): files |
| `\XTubeMax` | $4.8\cdot10^{-14}$ | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): LOCALMIN max tube |
| `\XUniqueCovers` | 839 | verify/out_exact/*.json (exactcheck.py) and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided): PASS, r_unique > dist_coord |
| `\ZetaN` | 693 | recount/zetasample.jsonl (zetasample.py; double precision, untrusted): certificates |
| `\ZetaOutside` | 44 | recount/zetasample.jsonl (zetasample.py; double precision, untrusted): stored code outside rho |
| `\ZetaOutsideIn` | 0 | recount/zetasample.jsonl (zetasample.py; double precision, untrusted): of these, inside rho2 |
| `\ZetaRatioMax` | $0.60$ | recount/zetasample.jsonl (zetasample.py; double precision, untrusted): max zeta2/zeta |
| `\ZetaRatioMed` | $0.18$ | recount/zetasample.jsonl (zetasample.py; double precision, untrusted): median zeta2/zeta |
| `\ZetaRatioMin` | $0.07$ | recount/zetasample.jsonl (zetasample.py; double precision, untrusted): min zeta2/zeta |
