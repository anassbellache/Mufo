function run_tests
    suite = testsuite('mufo/tests');
    results = run(suite);
    disp(table(results));
    if any([results.Failed])
        error('One or more tests failed');
    end
end
