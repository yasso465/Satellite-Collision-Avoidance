%% Verify published C63 result files (does not rerun simulations)
% Run from the repository root in MATLAB R2026b.
clear; clc;
assert(isfile('results/early_burn_full_catalog_summary.csv'), 'Missing nominal summary CSV.');
assert(isfile('results/early_burn_mc_10000_summary.csv'), 'Missing Monte Carlo summary CSV.');
assert(isfile('results/C63_excluded_validation.csv'), 'Missing excluded-object CSV.');
T = readtable('results/early_burn_full_catalog_summary.csv');
MC = readtable('results/early_burn_mc_10000_summary.csv');
E = readtable('results/C63_excluded_validation.csv');
candidate = string(T.Candidate);
row = find(candidate == "C63_060",1);
assert(~isempty(row),'C63 missing from nominal summary');
assert(T.ObjectsOK(row)==10866 && T.ObjectsError(row)==0, 'Unexpected retained-catalog counts');
assert(abs(T.GlobalMin_01s_km(row)-11.7346843073154)<1e-4,'Unexpected C63 nominal minimum');
rowMC=find(string(MC.Candidate)=="C63_060",1);
assert(~isempty(rowMC),'C63 missing from MC summary');
assert(MC.Trials(rowMC)==10000 && MC.Breaches10p5(rowMC)==0, 'Unexpected MC results');
ok = string(E.Status)=="OK";
assert(sum(ok)==6014 && sum(~ok)==3,'Unexpected excluded-object counts');
assert(abs(min(E.Min_1s_km(ok))-60.872769)<1e-4,'Unexpected excluded-object minimum');
fprintf('PASS: C63 saved-results consistency checks.\n');
fprintf('Nominal minimum: %.6f km; MC threshold crossings: %d/%d\n',T.GlobalMin_01s_km(row),MC.Breaches10p5(rowMC),MC.Trials(rowMC));
fprintf('Excluded propagated: %d; errors: %d; minimum: %.6f km\n',sum(ok),sum(~ok),min(E.Min_1s_km(ok)));
fprintf('NOTE: This checks CSV consistency only; it does NOT rerun orbital simulations.\n');
