/*! \page READMEemcCorrectionsChange Switch to the EMCal Correction Framework

Switching to the EMCal correction framework is a straightforward process. It should be possible to configure with
your desired settings in less than an hour.

# Transition your correction settings to the EMCal Corrections Framework

If you plan to test and verify that everything works the same (which we **strongly encourage!**) be certain not to delete your previous corrections yet!

# Configure your run macro or wagon

Follow the same procedure as described on \ref READMEemcCorrections.

# Test and verify the changes                     {#emcalCorrectionsVerifyChanges}

To test and verify the changes, we have a general procedure, as well as tools to help verify automatically. **NOTE: This procedure is oriented towards using EMCal Containers. If you do not use EMCal containers, this will only be a rough guide!**

The general procedure is as follows:

- Configure the EMCal Corrections Framework to run side-by-side with the previous corrections
- Add two copies of your %Analysis Task to your run macro to compare the result of using output from both sets of corrections
- %Compare and verify the result using an automatic comparison tool

## Testing the changes

First, you must add some additional options to your run macro. In particular, we need to configure the EMCal Corrections Framework to copy cells, clusters, and tracks to ensure that the two sets of corrections do not interfere with each other. This has two simple steps:

- A change to the AddTask configuration
- A change to the corrections configuration

### Configure the EMCal Correction Framework AddTask for side-by-side testing

To enable side-by-side testing, we will need to setup the copy of branches before setting up the EMCal Correction Framework AddTask and the old correction framework AddTask (ie. this code must be **executed before the Correction Framework AddTask** and **before the old correction framework AddTask** in your run macro or LEGO train). This is required to ensure that the two correction frameworks do not interfere with each other. To copy the proper input objects, use something like the code below (be certain to set ``IsEsd`` as appropriate):

~~~{.cxx}
// Cells
bool IsEsd = (iDataType == kEsd);
AliEmcalContainerUtils::InputObject_t inputObject = AliEmcalContainerUtils::kCaloCells;
TString inputObjectBranchName = AliEmcalContainerUtils::DetermineUseDefaultName(inputObject, IsEsd);
TString newBranchName = inputObjectBranchName;
newBranchName += "New";
AliEmcalCopyCollection * copyTaskCells = AddTaskEmcalCopyCollection(inputObject, inputObjectBranchName.Data(), newBranchName.Data());

// Clusters
inputObject = AliEmcalContainerUtils::kCluster;
inputObjectBranchName = AliEmcalContainerUtils::DetermineUseDefaultName(inputObject, IsEsd);
newBranchName = inputObjectBranchName;
newBranchName += "New";
AliEmcalCopyCollection * copyTaskClusters = AddTaskEmcalCopyCollection(inputObject, inputObjectBranchName.Data(), newBranchName.Data());

// Tracks
inputObject = AliEmcalContainerUtils::kTrack;
inputObjectBranchName = AliEmcalContainerUtils::DetermineUseDefaultName(inputObject, IsEsd);
newBranchName = inputObjectBranchName;
newBranchName += "New";
AliEmcalCopyCollection * copyTaskTracks = AddTaskEmcalCopyCollection(inputObject, inputObjectBranchName.Data(), newBranchName.Data());
~~~

The names in ``newBranchName`` determine the name of the new branches. In principle, this could be anything, but we strongly recommend using the "usedefault" name and then adding "New" onto the end. So for AODs, it would be "emcalCellsNew" for cells, "caloClustersNew" for clusters, and "tracksNew" for tracks.

### Configure the corrections for side-by-side testing

Before beginning, be certain to read the introduction to YAML and configuring corrections available [here](\ref READMEemcCorrections).

You can map the current corrections to the new ones with the following table. Note that each new correction is preceded by the name `AliEmcalCorrection`:

| Current Correction Name | New Correction Name       |
| ----------------------- | ------------------------- |
| Tender                  | CellBadChannel            |
|                         | CellEnergy                |
|                         | CellTimeCalib             |
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

### Configure the input objects for the corrections

In addition to determining the proper correction settings, the input objects need to be changed slightly. In parituclar, we need to change the branch names to the names defined when copying objects (described above)!

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

### Special note on track containers

When running the tasks side-by-side, we are forced to use non-standard input object names (such as "tracksNew"). If a task implements the "usedefault" pattern then it will likely not handle the non-standard track names properly - it will create a track container for the standard name and a particle container for the non-standard name. This will almost certainly lead to disagreements, so it needs to be fixed! To fix it, you will need to remove the particle container and add a new track container based on the (for example) "tracksNew" branch. In code, it will look like,

~~~{.cxx}
// Remove wrong particle container
taskName->RemoveParticleContainer(0);
// Add track container. The macro won't do this automatically, because the branch name is not "tracks" or "Tracks"
AliTrackContainer * newTracks = new AliTrackContainer(newFrameworkTracks.Data());
// Set the cuts for the new track container
// min pt, for example
newTracks->SetParticlePtCut(minTrackPt);
// Add it to the jet finder
taskName->AdoptTrackContainer(newTracks);
~~~

Some example tasks for which this applies include the jet finder, the EMCal sample tasks, and perhaps even your own analysis task!

## Verify the changes

You are all set! Now run the run macro as normal. Once it is finished, we have a python script to compare the output histograms automatically. It is available <a href="https://raw.githubusercontent.com/ALICEYale/alice-yale-dev/master/analyses/utilities/compareHistos.py" download>here</a> (you may need to right click -> Save Link As..). To run it, you need to pass the path of the ``AnalysisResults.root`` file, as well as the name of the output list from your task that was generated with the new corrections. For example,

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

If the outputs do not match initially, please double check your settings! If they still do not match, then consider running a statistical bin-by-bin comparsion. To do so, use the ``-s`` option:

~~~{.sh}
python compareHistos.py -f ../exampleDir/AnalysisResults.root -s
~~~

When running this statistical test, you can set the threshold for the fractional disagreement using the ``-t`` option:

~~~{.sh}
python compareHistos.py -f ../exampleDir/AnalysisResults.root -s -t 0.01
~~~

After all that, if they still do not match, then please let us know!

Help for the script is available with `python compareHistos.py --help`. (If you are using aliBuild, you may need to
set your ``$PYTHONPATH`` variable. In bash, you can do this by setting
``export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH``. Other shells may vary.)

# Some suggestions and tips if the tests disagree

Testing is very important, but it can take **a few iterations** to get all of your settings right! If they
don't match the first time, check your settings again closely. It is often helpful to test locally with a
relatively small number of events (perhaps ~1000) to allow rapid iteration. Anecdotally, we have observed the
time cuts, non-linearity function, and clusterizer type to cause many of the issues. In addition, a few defaults
were updated in the EMCal and EMCal-Jet sample tasks that are work checking. Those values include:

| Settings               | Previous value     | New default          | Reason for change                         |
| ---------------------- | ------------------ | -------------------- | ----------------------------------------- |
| cell time cuts         | +/- 50 us          | +/- 1                | Clarify they are off by default. This can make a difference in some periods! |
| Non-linearity function | kBeamTestCorrected | kBeamTestCorrectedv3 | Update in EMCal %Analysis Recommendations |
| Mass for track prop    | Pion mass          | PID mass hypothesis  | Experts decided this was the best setting |
| Point for track prop   | Vertex             | DCA                  | Experts decided this was the best setting |

Lastly, please note that due to small differences in the Tender clusterizer and the New Corrections Framework
clusterizer, exact agreement may not be possible. However, the differences should be very small. Note that this
does **not** apply to the standalone clusterizer (where perfect agreement is expected)!

*/
