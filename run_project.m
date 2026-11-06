function run_project(mode)
% RUN_PROJECT  Main entry point for the Satellite Collision Avoidance study.
%
% Usage from the repository root:
%   run_project()          % quick, self-contained spacecraft demonstration
%   run_project("demo")    % same as above
%   run_project("full")    % full retained-catalog, Monte Carlo and excluded-set runs
%
% The full mode can take a long time and requires MATLAB R2026b with
% Aerospace Toolbox. It uses a temporary working directory so the source
% scripts can read their expected filenames without moving repository data.

if nargin < 1 || isempty(mode)
    mode = "demo";
end
mode = lower(string(mode));

projectRoot = fileparts(mfilename('fullpath'));
sourceDir = fullfile(projectRoot, 'src');
dataDir = fullfile(projectRoot, 'data');
resultsDir = fullfile(projectRoot, 'results');

if ~isfolder(sourceDir)
    error('run_project:MissingSource', 'Cannot find src folder: %s', sourceDir);
end
if ~isfolder(resultsDir)
    mkdir(resultsDir);
end
addpath(sourceDir);

switch mode
    case "demo"
        fprintf('\n=== Satellite Collision Avoidance: C63 demo ===\n');
        fprintf('This quick demo propagates the spacecraft only; it does NOT screen debris.\n\n');
        executeScript(fullfile(sourceDir, 'demo_C63_two_body.m'));
        fprintf('\nDemo finished. For the catalog-based study, run run_project("full").\n');

    case "full"
        required = ["debris_matlab.tle", "debris_leo.tle", "debris_filtered.tle"];
        missing = strings(0,1);
        for k = 1:numel(required)
            if ~isfile(fullfile(dataDir, required(k)))
                missing(end+1,1) = required(k); %#ok<AGROW>
            end
        end
        if ~isempty(missing)
            error('run_project:MissingData', ...
                'Missing required files in data/: %s', strjoin(missing, ', '));
        end

        fprintf('\n=== FULL RESEARCH PIPELINE ===\n');
        fprintf('MATLAB version: %s\n', version);
        if ~license('test','Aerospace_Toolbox')
            error('run_project:MissingToolbox', ...
                'Aerospace Toolbox license is not available in this MATLAB session.');
        end

        workDir = tempname;
        originalDirectory = pwd;
        mkdir(workDir);
        cleanupWork = onCleanup(@() cleanupFullRun(workDir, originalDirectory)); %#ok<NASGU>
        for k = 1:numel(required)
            copyfile(fullfile(dataDir, required(k)), workDir);
        end

        fprintf('\n[1/3] Nominal retained-catalog validation...\n');
        cd(workDir);
        executeScript(fullfile(sourceDir, 'validate_early_burn_full_catalog.m'));
        cd(projectRoot);
        copyResult(workDir, resultsDir, 'early_burn_full_catalog_results.csv');
        copyResult(workDir, resultsDir, 'early_burn_full_catalog_summary.csv');

        fprintf('\n[2/3] Paired 10,000-trial Monte Carlo...\n');
        cd(workDir);
        executeScript(fullfile(sourceDir, 'monte_carlo_early_burn_10000.m'));
        cd(projectRoot);
        copyResult(workDir, resultsDir, 'early_burn_mc_10000_trials.csv');
        copyResult(workDir, resultsDir, 'early_burn_mc_10000_summary.csv');

        fprintf('\n[3/3] Excluded-object C63 validation...\n');
        cd(workDir);
        executeScript(fullfile(sourceDir, 'validate_C63_excluded.m'));
        cd(projectRoot);
        copyResult(workDir, resultsDir, 'C63_excluded_validation.csv');

        fprintf('\n=== PIPELINE FINISHED ===\n');
        fprintf('Generated result files are in: %s\n', resultsDir);
        fprintf('Review the logs and CSVs before making scientific claims.\n');

    otherwise
        error('run_project:UnknownMode', ...
            'Unknown mode "%s". Use run_project("demo") or run_project("full").', mode);
end
end

function executeScript(scriptPath)
if ~isfile(scriptPath)
    error('run_project:MissingScript', 'Cannot find script: %s', scriptPath);
end
run(scriptPath);
end

function copyResult(sourceDirectory, destinationDirectory, fileName)
sourceFile = fullfile(sourceDirectory, fileName);
if ~isfile(sourceFile)
    error('run_project:MissingOutput', ...
        'Expected output was not generated: %s', sourceFile);
end
copyfile(sourceFile, fullfile(destinationDirectory, fileName), 'f');
end

function cleanupFullRun(workDirectory, originalDirectory)
if isfolder(originalDirectory)
    cd(originalDirectory);
end
if isfolder(workDirectory)
    rmdir(workDirectory, 's');
end
end
