% run_all_tests.m
% Run every test_*.m in this folder and report a grep-able summary. Each test
% prints its own <NAME>_PASS / <NAME>_FAIL sentinel; this aggregator captures
% stdout, greps for _FAIL, and prints ALL_TESTS_PASS / ALL_TESTS_FAIL.
setup_paths();
this_dir = fileparts(mfilename('fullpath'));
files = dir(fullfile(this_dir, 'test_*.m'));
all_ok = true;
for k = 1:numel(files)
    name = files(k).name;
    fprintf('\n===== %s =====\n', name);
    out = evalc('run(fullfile(this_dir, name))');
    disp(out);
    if contains(out, '_FAIL')
        all_ok = false;
        fprintf('  -> %s reported a FAILURE\n', name);
    end
    close all;
end
if all_ok, disp('ALL_TESTS_PASS'); else, disp('ALL_TESTS_FAIL'); end
