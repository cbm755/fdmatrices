files = dir('.');

totaltime = cputime();
num_tests = 0;
num_failed = 0;
for i=1:length(files)
  fname = files(i).name;
  % tests should be files that start with 'test_'
  if ( (~files(i).isdir) && strncmp(fname, 'test_', 5) )
    testtime = cputime();
    f = str2func(fname(1:end-2));
    num_tests = num_tests + 1;
    disp(['** Running test(s) in: ' fname ]);
    pass = f();
    testtime = cputime() - testtime;
    if all(pass)
      fprintf('** PASS: %s  [%g sec]\n', fname, testtime);
    else
      fprintf('** FAIL: %s  [%g sec]\n', fname, testtime);
      num_failed = num_failed + 1;
      pass
    end
  end
end

totaltime = cputime() - totaltime;
fprintf('\n***** Passed %d/%d tests passed (%g seconds) *****\n', ...
        num_tests-num_failed, num_tests, totaltime);
if (num_failed > 0)
  disp('***** WARNING: some tests failed *****');
end
