# Using the runTests driver to run the test

The runTests helper script can be used to simplify creation and running of new
tests. E.g.:

    runTests gun

will read the `test.cfg` configuration file where it will find the specification
for the "gun" test and it will expand it to a shell script and run it. You can
see what runTests will do by supplying the `--dry-run` / `-n` option.

Each test consists of a variable number of steps, e.g.:

    example:
      steps:
        - name: SIM 
        - name: REC
    
By default a step will expand to:

    # Running Test example
    cd test/<test name> && rm -rf *.root *.dat *.log fort* hlt hough raw* recraw/*.root recraw/*.log
    # Step: <step name>
    # Variant default

    aliroot -b -q <lowercase step name>.C 2>&1 | tee  <lowercase step name>.log
    mv syswatch.log <lowercase step name>watch.log

Where in this case:

- <test name> is `example`
- <step name> is either SIM or REC
- <lowercase step name> is either `sim` or `rec`

Steps of a same test are guaranteed to be performed in order, as specified in
the configuration file. E.g. SIM will always preceed REC for the above test.

# Customization of the steps.

By default the steps will do the typical aliroot macro execution, as
illustrated above, simply because that's the most common usecase.

It is however possible to customize a step by supplying the following
parameters:

- pre: What to execute before the step
- driver: what to execute as a given step
- post: what to execute after the step
- setup: this is a per-test customization which allows to specify
  what to do before all steps.

For example the `fastjetload` test, also defined in test.cfg, looks like:

    loadfastjet:
      setup: ""
      steps:
        - name: LOAD
          pre: "export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/lib:/usr/lib64"
          driver: "root -l -b -q ${ALICE_PHYSICS}/PWGJE/macros/TestLoadFastJet.C 2>&1"
          post: ""

This will expand to simply:

   # Running Test loadfastjet

   # Step: LOAD
   # Variant default
   export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/lib:/usr/lib64
   root -l -b -q ${ALICE_PHYSICS}/PWGJE/macros/TestLoadFastJet.C 2>&1 

## Variants

One of the main reason to have a script like `runTests`, apart from easier
parsing of the test details, is to allow for what we call "variants": simple
changes of the default behavior of the test which share most of the
configuration, but then extend it in order to use special tools, for example
valgrind or igprof. Variants are defined per step and can be specified via a
set of additional configurable options: 

- prefix: a prefix to the driver part.

For example if you want to define an igprof_performance and igprof_memory
variant for the `gun` test you can do:

    gun:
      steps:
        - name: SIM
          variants:
            default: {}
            igprof_memory:
              prefix: "igprof -mp -o igprof.%(test_name)s_%(step_name)s_MEMORY.gz "
            igprof_performance:
              prefix: "igprof -pp -o igprof.%(test_name)s_%(step_name)s_PERFORMANCE.gz "

You can specify which variant to run by passing them to `runTests` via the
`--variant` option as comma separated labels.  E.g. `test/runTests gun
--variant default,igprof_memory`.

FIXME: add the ability to specify a label for a given test, so that tests can be
       selected by label (e.g. run all the tests in the `pullrequest` label).
FIXME: add the ability to specify "cost" of a given test, so that we can select only
       those which add up to a maximum amount of time.
FIXME: "default" variant should always be there.
FIXME: "default" variant should be run if the requested variant is not found in
       a given step.
FIXME: add ability to have a step depend on the step of a different test, for example
       have a "PbPbsimulation" test which has a SIM step and then have multiple "XYZreco"
       tests which have a REC step which depends on "PbPbsimulation/SIM" so that simulation
       is done only once.
