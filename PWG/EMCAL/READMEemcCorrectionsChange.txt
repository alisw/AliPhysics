/*! \page READMEemcCorrectionsChange Switch to the EMCal Correction Framework

Swtiching to the EMCal correction framework is a straightforward process. It should be possible to configure with
your desired settings in less than an hour.

Before using the EMCal correction framework, it is exteremly important that you are updated enough to be using
EMCal containers with your analysis. See [here](\ref READMEchangefw) for instructions to update, if you need. 

# Transition your correction settings to the EMCal Corrections Framework

If you plan to test and verify that everything works the same (which we strongly encourage!) be certain not to delete your previous corrections!

# Configure your run macro or wagon

Follow the same procedure as described on \ref READMEemcCorrections.

# Test and verifying the changes                     {#emcalCorrectionsVerifyChanges}

To test and verify the changes, we have a general procedure, as well as tools to help verify automatically. **NOTE: This procedure only works for analyses already using EMCal Containers! Older tasks must first update to at least use EMCal containers to use this tool!**

The general procedure is as follows:

- Configure the EMCal Corrections Framework to run side-by-side with the previous corrections
- Add two copies of your %Analysis Task to your run macro to compare the result of using output from both sets of corrections
- %Compare and verify the result using an automatic comparison tool

## Testing the changes

First, you must add some additional options to your run macro. In particular, we need to configure the EMCal Corrections Framework to copy cells, clusters, and tracks to ensure that the two sets of corrections do not interfere with each other. This has two simple steps:

- A change to the AddTask configuration
- A change to the corrections configuration

### Configure the EMCal Correction Framework AddTask for side-by-side testing

To enable side-by-side testing, simply add the following line after you setup the EMCal Correction Framework AddTask:

~~~{.cxx}
// Create a copy of cells, clusters, and tracks to compare against the current correction framework
correctionTask->SetCreateNewObjectBranches(true);
~~~

NOTE: It is imperative that you setup the EMCal Corrections Framework **before** the EMCal tender task!

If you configured it properly, it should look something like:

~~~{.cxx}
AliEmcalCorrectionTask * correctionTask = AddTaskEmcalCorrectionTask();
correctionTask->SelectCollisionCandidates(kPhysSel);
// Set the run period, same as the track container
// If you derived from the file "runEMCalJetSampleTask.C", then it is likely stored under "sRunPeriod.Data()"
correctionTask->SetRunPeriod("LHC11h");
// Create a copy of cells, clusters, and tracks to compare against the current correction framework
correctionTask->SetCreateNewObjectBranches(true);
// Set the user configuration file, assuming that your file is called "userTestConfiguration.yaml" and is located in
// the current directory. This also supports alien:// paths!
correctionTask->SetUserConfigurationFilename("userTestConfiguration.yaml");
// Initialize the configuration files in the correction task
// It is EXTREMELY important to run this function in your run macro!
correctionTask->InitializeConfiguration();

// Configure the EMCal tender and other current corrections below here!
~~~

### Configure the corrections for side-by-side testing

Before beginning, be certain to read the introduction to YAML and configuring corrections available [here](\ref READMEemcCorrections).

You can map the current corrections to the new ones with the following table. Note that each new correction is preceded by the name `AliEmcalCorrection`:

| Current Correction Name | New Correction Name       |
| ----------------------- | ------------------------- |
| Tender                  | CellBadChannel            |
|                         | CellEnergy                |
|                         | CellBadChannel            |
| Clusterizer             | Clusterizer               |
| ClusterMaker            | ClusterExotics            |
|                         | ClusterNonLinearity       |
| ClusTrackMatcher        | ClusterTrackMatcher       |
| HadCorr                 | ClusterHadronicCorrection |

Using the table above, we need to make sure that the current and the new corrections are configured the same. Compare the settings in your run macro with those in ``$ALICE_PHYSICS/PWG/EMCAL/config/AliEmcalCorrectionConfiguration.yaml``, known as the default file. In particular, check that each parameter value in the YAML file matches the value for that variable in your run macro (and not the other way around -- it is okay if there is no matching field for every field in the AddTask). The variable names should often be the same. If values are different in the YAML file, then modify the values in **``userTestConfiguration.yaml``**, not in ``AliEmcalCorrectionConfiguration.yaml``! To make the change, write it in yourself or copy the relevant structure of the default file to your user file and modify there. For instance, in the default file, the ``cellE`` parameter in the Clusterizer is to ``0.05``. If I wanted to change that in the user file, you would add:

~~~
Clusterizer:
    cellE: 0.1    # Changed this value
~~~

Any setting that you have in the user file will override the default file!

## Configure your task to run twice

Configure your analysis task to run twice - once with the current corrections and once with the new corrections. To do so, configure your task once with different cluster and track (if needed) containers. Note that if you use jets, you need to run the jet finder twice - once for the old framework, and once for the new framework. For example, if you run on AODs, and set your cluster and track containers in run macro using "usedefault", then it would look something like:

~~~{.cxx}
AliAnalysisTaskMyAnalysis * task =
new AliAnalysisTaskMyAnalysis("usedefault", // clusters
                              "usedefault", // tracks
                              WhateverOptionsYouNormallyUse);
task->OptionOne();
// etc...

AliAnalysisTaskMyAnalysis * taskNew =
new AliAnalysisTaskMyAnalysis("caloClustersNew", // clusters
                              "tracksNew", // tracks
                              WhateverOptionsYouNormallyUse);
taskNew->OptionOne();
// etc...
~~~

**Note**: If you don't use the "usedefault" pattern, then you'll have to set the names manually. You should just set the tracks and clusters names as usual.

### Special note on the track containers

As noted above, if you use jets, you will need to configure two jet finders. Due to using nonstandard track and cluster names to compare the two frameworks, the jet finder will not be configured properly if you just pass the names. Instead, you will need to remove the particle container and add a new track container based on the "tracksNew" branch. In code, it will look like,

~~~{.cxx}
// Remove wrong particle container
pFullJet02TaskNew->RemoveParticleContainer(0);
// Add track container. The macro won't do this automatically, because the branch name is not "tracks" or "Tracks"
AliTrackContainer * newTracks = new AliTrackContainer(newFrameworkTracks.Data());
// Set the cuts for the new track container
newTracks->SetParticlePtCut(minTrackPt);
// Add it to the jet finder
pFullJet02TaskNew->AdoptTrackContainer(newTracks);
~~~

If your analysis uses a track container, and you use the "usedefault" pattern in your analysis task AddTask, repeat the above lines (adapted for your task) for your task's track container.

## Verify the changes

You are all set! Now run the run macro as normal. Once it is finished, we have a python script to compare the output histograms automatically. It is available <a href="https://gitlab.cern.ch/ALICEYale/alice-yale-dev/raw/master/analyses/utilities/compareHistos.py" download>here</a> (you may need to right click -> Save Link As..). To run it, you need to pass the path of the ``AnalysisResults.root`` file, as well as the name of the output list from your task that was generated with the new corrections. For example,

~~~{.sh}
python compareHistos.py -f ../exampleDir/AnalysisResults.root -n MyAnalysisOutputListFromNewCorrections
~~~

If you don't know the name of your output list, just pass an invalid name to the ``-n`` argument and it will show you the available options:

~~~{.sh}
python compareHistos.py -f ../exampleDir/AnalysisResults.root -n invalidName
~~~

In exceptional cases, you will need to pass the name corresponding to the name of the output list from your task that was run with the current corrections. If needed, the script will inform you. For example,

~~~{.sh}
python compareHistos.py -f ../exampleDir/AnalysisResults.root -n MyAnalysisOutputListFromNewCorrections -o MyAnalysisOutputListFromOldCorrections
~~~

If the outputs do not match initially, please double check your settings! If they still do not match, then please let us know!

Help for the script is available with `python compareHistos.py --help`. (If you are using aliBuild, you'll need to set your ``$PYTHONPATH`` variable. In bash, you can do this by setting ``export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH``. Other shells may vary.)

*/
