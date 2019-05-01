/*! \page READMEemcCorrections EMCal Correction Framework

# EMCal Correction Framework

The new EMCal correction framework centralizes all of the analysis-level EMCal corrections into a single task:
AliEmcalCorrectionTask. This task coordinates the running of a configurable list of individual "correction components", where a
correction component means a single correction to be performed: cell-level energy calibration, or cluster-level non-linearity
correction, or cluster-track matching, etc. The set of corrections to run, as well as their order and their configurable
parameters, is set in a %YAML configuration file. A default config file containing the centralized standard set of parameters is
shared by all users, and each user additionally writes their own user config file, overwriting desired variables.

The motivation of this new approach is the following:
- Centralize code so that there are not many similar versions of code used by different people
- Provide a simple and unified interface to make life easier for both users and developers
- Reduce likelihood of mistakes for users (make it obvious to a user if they are deviating from a recommended config, and avoid setting arguments through ``AddTask`` macros)

The following correction components are available:
- [CellBadChannel](\ref AliEmcalCorrectionCellBadChannel) -- Sets cells marked as bad to E = 0, using OADB bad channel map.
- [CellEnergy](\ref AliEmcalCorrectionCellEnergy) -- Performs energy calibration of cells, using OADB calibration.
- [CellTimeCalib](\ref AliEmcalCorrectionCellTimeCalib) -- Performs time calibration of cells, using OADB calibration.
- [EmulateCrosstalk](\ref AliEmcalCorrectionCellEmulateCrosstalk) -- Emulate cell-level crosstalk.
- [Clusterizer](\ref AliEmcalCorrectionClusterizer) -- Clusterizes a collection of cells into a collection of clusters.
- [ClusterExotics](\ref AliEmcalCorrectionClusterExotics) -- Flags exotic clusters for removal from the cluster collection.
- [ClusterNonLinearity](\ref AliEmcalCorrectionClusterNonLinearity) -- Corrects cluster energy for non-linear response.
- [ClusterNonLinearityMCAfterburner](\ref AliEmcalCorrectionClusterNonLinearityMCAfterburner.h) -- This is an additional correction for MC. The effect is that the pi0 mass position is the same for data and MC. It is only defined for specific periods.
- [ClusterTrackMatcher](\ref AliEmcalCorrectionClusterTrackMatcher) -- Matches each track to a single cluster, if they are in close enough proximity.
- [ClusterHadronicCorrection](\ref AliEmcalCorrectionClusterHadronicCorrection) -- For clusters that have one or more matched tracks, reduces the cluster energy in order to avoid overestimating the particle's energy.
- [PHOSCorrection](\ref AliEmcalCorrectionPHOSCorrections) -- Perform PHOS correction via an interface to the PHOS tender.

This new correction task unifies what was previously done by tasks such as:
 - [EMCal Tender](\ref AliEmcalTenderTask)
 - [ClusterizerFast](\ref AliAnalysisTaskEMCALClusterizeFast)
 - [ClusterMaker](\ref AliEmcalClusterMaker)
 - [ClusTrackMatcher](\ref AliEmcalClusTrackMatcherTask)
 - [HadCorr](\ref AliHadCorrTask)

...and other similar tasks.

# Using the EMCal Corrections Framework

Below, follow the instructions to implement the new Correction Task in your analysis. Further down the page, you will find
instructions on how to transition from the old correction method(s) to the new correction method, including how to configure the
corrections, and how to **test your analysis** results in the new vs. old correction method. Testing your corrections configuration
is an extremely important step, as failure to do so could lead to incorrect results!

# %Setup the Correction Task 

There are two steps:
- [Configure the Correction Add Task](\ref configureEMCalCorrectionRunMacro)
- [Configure the Corrections](\ref configureEMCalCorrections)

# Switching to the EMCal Corrections Framework

For those who are switching to this framework, a special page has been prepared to explain the process.
This page explains how to configure the correction framework alongside the previous set of corrections, allowing
you to test and show that you get the same results with the new and old corrections!
For the instructions, please see \subpage READMEemcCorrectionsChange.

# Configure the Correction Add Task (or train wagon)            {#configureEMCalCorrectionRunMacro}

To enable the correction task, add the following lines to your run macro:

~~~{.cxx}
AliEmcalCorrectionTask * correctionTask = AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask();
correctionTask->SelectCollisionCandidates(kPhysSel);
// Set the user configuration file, assuming that your file is called "userConfiguration.yaml" and is located in
// the current directory. This also supports alien:// paths!
correctionTask->SetUserConfigurationFilename("userConfiguration.yaml");
// Initialize the configuration and corrections in the correction task
// It is EXTREMELY important to run this function in your run macro!
correctionTask->Initialize();
~~~

Note that this should be added after `CDBconnect`, which is required to use the correction framework.

## LEGO Train Wagon                                             {#emcCorrectionsLEGOTrainWagon}

There are some special procedures for the LEGO train. There will be a centralized EMCal Correction Task wagon which
contains standard settings such as physics selection, etc. Then the user will create a wagon which calls
`PWG/EMCAL/macros/ConfigureEmcalCorrectionTaskOnLEGOTrain(std::string suffix)`. This will return an EMCal Correction Task
for you to configure with your particular %YAML configuration file and then initialize.

To setup this wagon, use the macro listed above. One argument is required: the suffix. This suffix should be some unique
identifier for your analysis - it could be your name, your analysis, whatever you like. The precise value only matters
if you are using [specialization](\ref emcCorrectionsSpecialization) - in such a case, your suffix for the configuration
wagon must match the specialization suffix.

Then in the macro customization of your configuration wagon, you specify your %YAML configuration file and initialize
the configuration. It should looks something like the following:

~~~{.cxx}
// Set your user configuration
__R_ADDTASK__->SetUserConfigurationFilename("userConfiguration.yaml");
// It is extremely important to pass "true" to Initialize().
__R_ADDTASK__->Initialize(true);
~~~

Note that your configuration wagon should depend on the centralized correction task wagon, but your tasks (such as jet
finders, user tasks, etc) should depend **only** on the centralized correction task wagon.
They should not depend on your configuration wagon! One side effect of this configuration is that intermediate train
test may fail with an error about the correction task not being configured. Intermediate tests failing is not necessarily
a problem in itself - instead, check on the result of the full train test. If that test was successful, the train is fine
and can be started. For more, see the [FAQ answer](\ref emcCorrectionsLEGOTrainTests).

# Configuring Corrections                                       {#configureEMCalCorrections}

The corrections configuration file is specified via a file written in the %YAML markup language. %YAML is quite
readable and is commonly used for configuration files. Although it should be straightforward to read and write,
more information on its use can be found through the examples available [here](https://en.wikipedia.org/wiki/YAML).
**All** configuration and options are set through the %YAML file, from setting the clusterizer type to the name of
the tracks!

The configuration is determined by two files: the base default file (located in ``$ALICE_PHYSICS/PWG/EMCAL/config/``),
as well as the user configuration file. Any changes that you need to make should be done in the user configuration file!
Note that any setting that you have in the user file will override the default file!

## Features of the configuration file

We will discuss a number of features including:
- Setting cells, clusters, and tracks (input objects)
- Comments
- Enumerations
- Advanced usage options

To discuss the configuration and introduce, consider the configuration excerpt below as an example:

~~~
# This is a comment! Comments can also be on inline - everything after the hash is a comment
inputObjects:                               # Define all of the input objects for the corretions
    cells:                                  # Configure cells
        # The user can select the name of each container
        # The names don't need to correspond to any particular scheme
        defaultCells:
            # Note that we support the "usedefault" pattern of setting branch names.
            branchName: "usedefault"
    clusterContainers:                      # Configure clusters
        # The user can select the name of each container
        # The names don't need to correspond to any particular scheme
        defaultClusterContainer:            # Name of a cluster input (corresponds to a cluster container)
            # Sets the branch name
            branchName: "usedefault"
            # Takes all default cuts!
        defaultClusterContainer_1:          # Name of another cluster input which inherits from defaultClusterContainer
            # The branch name is inherited from defaultClusterContainer!
            minE: 0.0                       # Cluster min E
            minPt: 0.0                      # Cluster min pt
        defaultClusterContainer_2:          # Name of another cluster input which inherits from defaultClusterContainer
            # The branch name is inherited from defaultClusterContainer!
            minE: 0.0                       # Cluster min E
            minPt: 0.0                      # Cluster min pt
            # Handled particularly for cluster containers
            clusNonLinCorrEnergyCut: 0.15   # formerly "minPt" in the non-linearity AddTask()
    trackContainers:
        # The user can select the name of each container
        # The names don't need to correspond to any particular scheme
        defaultTrackContainer:              # Name of a track input (corresponds to a track container)
            # Sets the branch name
            branchName: "usedefault"
            # Takes all default cuts!
        defaultTrackContainer_1:            # Name of another track input which inherits from defaultTrackContainer
            # The branch name is inherited from defaultClusterContainer!
            minPt: 0.15                    # Track min pt
            trackFilterType: kHybridTracks          # Set the track filter type. Check the documentation for possible values
            #aodFilterBits:                         # Can also set AOD filter bits. Check the docuemntation for more information
            #    - 16
            #    - 1
sharedParameters:
    # You can define here any parameters that you would like to share between tasks. See the documentation for more information
# This next section defines the settings for the "AliEmcalCorrectionCellEnergy" correction.
# Note how the "AliEmcalCorrection" part of the name is not required when defining the settings
CellEnergy:
    # All tasks can be enabled or disabled by setting the "enabled" property.
    enabled: true
    # Most tasks can create histograms, which are useful for QA or comparison with other corrections.
    createHistos: true
# This defines the settings for the "AliEmcalCorrectionCellBadChannel" correction.
CellBadChannel:
    enabled: true
    createHistos: true
# This defines the settings for the "AliEmcalCorrectionCellTimeCalibration" correction.
CellTimeCalib:
    enabled: true
    createHistos: true
# This defines the settings for the "AliEmcalCorrectionClusterizer" correction.
Clusterizer:
    enabled: true
    createHistos: true
    # Sets the cluster branch name for this particular task.
    clusterBranchName: "sharedParameters:clusterBranchName"
    # Sets the Clusterizer type based on the same familiar enumeration.
    # Note that you should _not_ include the prefix as you usually would for setting an enumeration. Only list the value.
    clusterizer: kClusterizerv3
    cellE: 0.05
    seedE: 0.1
    cellTimeMin: -1                    # Min cell time (s)
    cellTimeMax: +1                    # Max cell time (s)
    clusterTimeLength: 1               # Maximum time difference between the digits inside EMC cluster (s)
    w0: 4.5
~~~

### Setting cells, clusters, and tracks (input objects)

To properly configured the Correction Tasks, we need to indicate which cells, clusters, and tracks are needed (collectively referred to as input objects). These are configured in the ``inputObjects`` section. Within this section, we configure ``cells`` for input cells, ``clusterContainers`` for input cluster containers, and ``trackContainers`` for input tracks. Note that each subsection can contain multiple objects.

Each input object follows the same pattern:

~~~
objectName:                         # Used to refer to this object in the YMAL configuration
    branchName: "usedefault"        # The name of the branch that we want. The "usedefualt" pattern is supported here
    option: value                   # We support a wide variety of options, listed below.
~~~

Although it is not required for your task to use EMCal containers, they are used internally in the Correction Task, so it is important to understand a bit about how they work. In particular, any filtering of the collections (such as AOD Track Filtering) is applied for the corrections, but the filtering does not remove anything from the collections!
In addition, we always write cell/cluster/track corrections in place - it is generally not possible to write out new collections. Thus, it is then the users' responsibility to filter the collection using the information available! One particularly easy way to do so is to use EMCal containers, but it is not required! In the case that the user wants to manually filter, for example, tracks before running the Correction Framework, one can run AOD Track Filtering and then select no additional track filtering using ``kNoTrackFilter``.

The available configuration options are listed below. Direct configuration through the %YAML file is limited to a subset of all available configuration options. If you require additional options, please contact the developers - such additions are usually trivial. Alternatively, after calling ``Initialize()``, all containers are available for configuration as usual, so any additional options can be set by hand until implemented by the developers.

__Cells__:
| Property              | Value                             |
| --------------------- | --------------------------------- |
| branchName            | String selecting the branch name  |
| embedding             | True if the cells branch should be taken from an embedded (external) event |

__EMCal Containers__:
| Property              | Value                             |
| --------------------- | --------------------------------- |
| Container name        | This is set by the object name in the %YAML config |
| branchName            | String selecting the branch name  |
| embedding             | True if the branch should be taken from an embedded (external) event |
| minPt                 | Double setting the minimum pT of the container |
| minE                  | Double setting the minimum energy of the container |
| (minEta, maxEta)      | A pair of doubles to set the min and max eta of the container. Note that they must be set as a pair! |
| (minPhi, maxPhi)      | A pair of doubles to set the min and max phi of the container. Note that they must be set as a pair! |

__Clusters__ (Includes all options for EMCal Containers):
| Property              | Value                             |
| --------------------- | --------------------------------- |
| clusNonLinCorrEnergyCut | Double determining the cluster non-linearity energy cut |
| clusHadCorrEnergyCut  | Double determining the cluster hadronic cluster energy cut |
| defaultClusterEnergy  | String defining the default energy type of the container (kNonLinCorr, kHadCorr, ...) |
| includePHOS           | True if PHOS cluster should be included |

__Tracks__ (Includes all options for EMCal Containers):
| Property              | Value                             |
| --------------------- | --------------------------------- |
| trackFilterType       | Enumerated value determining the track filtering to apply to the container.  The enumerations are defined in AliEmcalTrackSelection and include: kNoTrackFilter, kCustomTrackFilter, kHybridTracks, kTPCOnlyTracks |
| aodFilterBits         | Sets values based on the evaluated bits (ie to set bit 3, pass 2^3 = 8). To set it properly, see below. To use aodFilterBits, Be sure to set the trackFilterType to kCustomTrackFilter! |
| trackCutsPeriod       | Sets the track cuts period for the track container. A default can also be set for all track containers using static methods of AliTrackContainer independent of the correction task. See AliTrackContainer for more information |

The aodFilterBits are set as an array. This is possible to express in two ways in YAML:

~~~
aodFilterBits: [1, 256]
# Alternatively
aodFilterBits:
 - 1
 - 256
~~~

### Comments

Comments include anything after a hash ("#"). Of course, they are not required, but are highly recommended! You can
add a comment on the same line as where you define a property - just add the "#" after you are done defining the property.

### Enumerations

Enumerations are supported as expected. Just set the value as you normally would. However, be sure that you don't
include class where it is defined!

### Example configuration and further information

For a survey of the available configuration options and information about the meaning of each, see the default
configuration file located in ``$ALICE_PHYSICS/PWG/EMCAL/config/``. This also serves as an example configuration,
Note that once the Correction Task is initialized using the ``Initialize()`` function, all corrections are
configured and created. Consequently, any additional configuration can be done manually if desired. However,
this approach is strongly discouraged. Instead, it is better to change the %YAML configuration file so that
there is a record of settings.

### Advanced usage options

There are a number of useful advanced options to make the Corrections Framework simpler and more pleasant to use. These options are described below.

#### Shared parameters

Often, a user will want to change some parameters in unison. Say, if a pt cut is changed, it should be changed
everywhere. In such a case, it is useful to able to define a variable so that one change will change things everything.
This can be accomplished by defining a parameters in the "shared parameters" section of the %YAML file. The name of the
parameter defined in the shared parameters section can be referenced in other areas of the file by prepending
``sharedParameters:`` to the parameter name. Consider the example below:

~~~
sharedParameters:
    aMinimumValue: 3
Correction1:
    exampleValue: "sharedParameters:aMinimumValue"
    ...
Correction2:
    anotherExample: "sharedParameters:aMinimumValue"
~~~

In the example, any change to ``aMinimumValue`` will be propagated to ``exampleValue`` in ``Correction1`` and
``anotherExample`` in ``Correction2``. Note that the parameter name (here, ``aMinimumValu``) can be anything that
the user desires. When setting the value, don't forget to prepend "sharedParameters:" (in our example,
"sharedParameters:aMinimumValu")!

Note that shared parameters are inherently somewhat limited. It can only retrieve values where the requested type
is arithmetic, string, or bool. The retrieved shared parameters value can only be of those same types. Note that
the shared parameters correspond to each configuration file.  ie. If ``sharedParameters:test`` is requested in the
first configuration file, then it will only look for the sharedParameters value in that configuration. Thus, if a
sharedParameter is requested in a later configuration file, the earlier configuration shared parameter values
**will not** be considered.

As an explicit example, to change the recalculation of the MC labels including the fractional energy on cell level
to ``false`` while leaving the cluster-track matcher setting the same in your user config, it should look like
(omitting unrelated values for brevity):

~~~
sharedParameter:
    enableFracEMCRecalc: false
Clusterizer:
    enableFracEMCRecalc: true
ClusterTrackMatcher:
    inducedTCardMinimum: "sharedParameter:enableFracEMCRecalc"
~~~

To change the ``enableFracEMCRecalc`` for both in your user config, you'll need to do:

~~~
sharedParameter:
    enableFracEMCRecalc: false
Clusterizer:
    enableFracEMCRecalc: "sharedParameters:enableFracEMCRecalc"
ClusterTrackMatcher:
    inducedTCardMinimum: "sharedParameters:enableFracEMCRecalc"
~~~

**It isn't enough to just change the value in the shared parameter of the user config!**. This is because the
``sharedParameters`` field of each configuration do not override each other.

## Running multiple corrections at once ("specializing")                      {#emcCorrectionsSpecialization}

Often, a user would like to run two nearly identical corrections. For instance, one could run two clusterizers with the same
configuration, but perhaps different input cells and output clusters. In such a case, the clusterizer can be "specialized", such that
each of the two clusterizers will inherit the same settings except for the cells. Consider the following example %YAML configuration
file:

~~~
inputObjects:
    cells:
        defaultCells:
            branchName: "usedefault"
        anotherCells:
            branchName: "otherCellsBranch"
    clusterContainers:
        defaultClusterContainer:
            branchName: "usedefault"
            minPt: 3
        anotherClusterContainer:
            branchName: "someOtherClustersBranch"
# Cell corrections ...
Clusterizer:
    enabled: true
    clusterizer: kClusterizerv2
    cellE: 0.05
    seedE: 0.1
    ...
    cellsNames:
        - defaultCells
    clusterContainersNames:
        - defaultClusterContainer
Clusterizer_mySpecialization:
    cellsNames:
        - anotherCells
    clusterContainersNames:
        - anotherClusterContainer
~~~

In this example, two clusterizers are configured: ``Clusterizer`` and ``Clusterizer_mySpecializedClusterizer``. Both have the exact same
configuration, with the exception that ``Clusterizer`` uses ``defaultCells`` for cells and ``defaultClusterContainer`` for clusters,
while ``Clusterizer_mySpecializedClusterizer`` uses ``anotherCells`` for cells and ``anotherClusterContainer`` for clusters (note: the
cells branch would have to be created elsewhere, say with ``AliEmcalCopyCollection``). 

To execute these corrections, we must select them in the AddTask of the Correction Task. To do so, set the suffix argument of the AddTask
to correspond with the specialization suffix that was set in the %YAML file. Continuing with the previous example, we would need to add:

~~~{.cxx}
// The default argument is an empty string, so we don't have to set it here.
AliEmcalCorrectionTask * correctionTask = AddTaskEmcalCorrectionTask(); 
// Create the specialized correction task
AliEmcalCorrectionTask * correctionTaskSpecialized = AddTaskEmcalCorrectionTask("mySpecialization");
~~~

Note that in doing this, it is then required that all additional corrections after the first specialized component are also specialized.
This is necessary to make the user's intentions clear. Continuing with the two clusterizers example, if we then want to run the cluster-track
matcher, then we must created both ``ClusterTrackMatcher`` and ``ClusterTrackMatcher_mySpecialization`` with the proper configuration.

It is extremely important to be careful to avoid apply corrections multiple times to the same collections! For instance, if running two
clusterizers on the same cells collection, then the cell corrections must be disabled for one of the two corrections! If the above example
had used the same cells, then it would have been required to disable them in one correction task (say, the "mySpecialization" task).
Note that the corrections being disabled by default does not automatically avoid this issue! Consider the following __WRONG__ example:

~~~
CellEnergy:
    enabled: true
    cellsNames:
        - defaultCells
Clusterizer:
    enabled: true
Clusterizer_example1:
    someSettings: true
Clusterizer_example2:
    someSettings: false
~~~

For specialized tasks "example1" and "example2", the correction chains are the following:

~~~
CorrectionTask_example1:
- CellEnergy
- Clusterizer_example1

CorrectionTask_example2:
- CellEnergy
- Clusterizer_example2
~~~

So now you've applied the Cell Energy correction twice! Instead, your configuration should look like:

~~~
CellEnergy:
    enabled: false
    cellsNames:
        - defaultCells
# This will cause the Cell Energy correction to only be performed in the task "example1".
CellEnergy_example1:
    enabled: true
Clusterizer:
    enabled: true
    # .. Can include any other settings
Clusterizer_example1:
    someSettings: true
Clusterizer_example2:
    someSettings: false
~~~

Now the correction chain will look like as intended:

~~~
CorrectionTask_example1:
- CellEnergy
- Clusterizer_example1

CorrectionTask_example2:
- Clusterizer_example2
~~~

### Note on specialization and LEGO trains

A practical way to use specialization on the LEGO trains is to implement a correction task wagon with sub wagons for
each specialization. By naming the sub wagon suffix the same as the specialization name, a specialized correction
task will be created for each suffix. Each specialized correction task will extract the appropriate components for
each correction chain from a single %YAML file. For an example, see train 2860 on the Jets_EMC_PbPb train. However,
the naming and dependencies of specializations with LEGO train subwagons requires special attention. In a normal run
macro, the user can always select the order of running correction task specializations. However, in the case of LEGO
train subwagons, the subwagons will always be run in **alphabetical order**. This can cause the correction tasks to
be executed in an unexpected order. For a concrete example, consider three specializations, named A, B, and C, where
A depends on B and C running first. If subwagons of a single correction task wagon are used to create these
specializations, then the "A" specialization will run **before** B and C. To run them in the proper order, an
additional correction task wagon must be created to specialize "A". It should look like:

~~~
FirstCorrectionTaskWagon
 -> Subwagon "B"
 -> Subwagon "C"

SecondCorrectionTaskWagon (depends on FirstCorrectionTaskWagon)
 -> Subwagon "A" (does not necessarily require a subwagon to create a specialization here, but it could be used)
~~~

# Using the output of the Correction Task                                    {#emcalCorrectionsOutput}

The correction generated by each component of the Correction Framework is written to the input objects TClonesArray **in place**. This means that all corrected values are immediately available to the user. How the user accesses those corrected values depends on whether their user task utilizes EMCal Containers. Both scenarios will be addressed. For both examples, it will involve retrieving clusters from an AOD with the branch name "caloClusters". More on branch names can be [here](\ref emcalContainerBranchNames).

## User task uses EMCal Containers

In your run macro or add task,
~~~{.cxx}
// If your user task in your AddTask is called "myTask"
// We add a container exactly as usual
AliClusterContainer * clusCont = task->AddClusterContainer("caloClusters");
// Configure your container as usual...
clusCont->SetClusECut(0.5);
// Note that by default, the cluster container will automatically remove clusters marked as exotic when
// you iterate over it (but only if you have run the Correction Framework to mark the exotic clusters!)
// We can also make the cut explicit here
clusCont->SetExoticCut(kTRUE);
~~~

For more information, see the [EMCal Containers documentation](\ref READMEcontainers).

## User task does not use EMCal Containers

In your main event loop (perhaps ``UserExec()``), you access the input object TClonesArray exactly as normal. For example,
~~~{.cxx}
// Retrieve the TClonesArray as usual
TClonesArray * clusters = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("caloClusters"));
// The corrected values are available after the Correction Framework has been run, but it is the _users_
// responsibility to filter the clusters.
// For example, to ensure that only EMCal clusters are considered, it would look something like
for (int i = 0; i < clusters->GetEntries(); i++) {
    AliVCluster * cluster = dynamic_cast<AliVCluster *>(clusters->At(i));
    if (cluster->IsEMCAL()) {
        Printf("EMCal cluster!")
    }
    else {
        Printf("Not an EMCal cluster");
    }
}
~~~

# FAQ                                                                       {#emcCorrectionsFAQ}

## I am seeing an error related to the EMCal geometry - what is wrong?      {#emcCorrectionsGeometryError}

For example, you may see an error similar to below.

~~~
F-AliEMCALGeometry::GetMatrixForSuperModule: Cannot find EMCAL misalignment matrices! Recover them either:
- importing TGeoManager from file geometry.root or
- from OADB in file OADB/EMCAL/EMCALlocal2master.root or
- from OCDB in directory OCDB/EMCAL/Align/Data/ or
- from AliESDs (not in AliAOD) via AliESDRun::GetEMCALMatrix(Int_t superModIndex).
Store them via AliEMCALGeometry::SetMisalMatrix(Int_t superModIndex)
~~~

This is due to the OCDB being inaccessible. Please ensure that you ran the `CDBconnet` task.

## Some of the LEGO train tests failed and I can't launch my train! Help!   {#emcCorrectionsLEGOTrainTests}

When a configure wagon is used, it is out of the main dependency chain, which can cause some LEGO train tests to fail.
In particular, this happens when the configure wagon is not run during a test, which is especially common during intermediate
tests (ie running some, but not all wagons). In the case that some of these intermediate tests fail, check the full train test.
If that test was successful, then those failed tests can safely be ignored.

However, train operators should be aware that the LEGO train framework will not allow a train to be launched if any test is
reported as failed. To workaround this issue, the operator should run a **fast test**, which will only run the baseline and
full train tests. As noted above, the full test should run successfully with the configure wagon, so assuming that both the
baseline and full train tests succeed, it is then possible to launch the train.

# Details on the framework and the corrections

You can see the code at ``$ALICE_PHYSICS/PWG/EMCAL/EMCALtasks``. The steering class is ``AliEmcalCorrectionTask``. The
individual corrections inherit from ``AliEmcalCorrectionComponent``, and are labeled ``AliEmcalCorrectionXXXX``. Note
that neither ``AliEmcalCorrectionTask`` nor ``AliEmcalCorrectionComponent`` inherit from ``AliAnalysisTaskEmcal``.
However, they provide similar functionality. The default configuration file is at ``$ALICE_PHYSICS/PWG/EMCAL/config/AliEmcalConfiguration.yaml``.

NOTE: If you are interested in how a particular correction works, you only need to look at the particular correction and its configuration!
There are many other details in the base and steering classes, but they are almost certainly not relevant!

# Developing a correction

If you are interested in developing a task that is shared by analyses using the EMCal, then the correction
framework is a great place to deploy it! The general approach is very similar to the steering in
``AliAnalysisTaskEmcal``. To create your task, your task name should start with ``AliEmcalCorrection``.
It should derive from ``AliEmcalCorrectionComponent`` and implement:

~~~{.cxx}
// Run once to initialize the component. These are run independent initializations.
virtual Bool_t Initialize();
// Execute in the first event for run dependent initialization. Same as in AliAnalysisTaskEmcal
virtual void ExecOnce();
// Executed before the first event. Same as in AliAnalysisTaskSE
virtual void UserCreateOutputObjects();
// Called each event. Same as in AliAnalysisTaskEmcal
virtual Bool_t Run();
// UserNotify() from AliAnalysisTaskSE
virtual Bool_t UserNotify();
~~~

To configure your task in the ``Initialize()`` function with the %YAML configuration, you must define the
object and then get it via the ``GetProperty("propertyName", property)`` function defined in
``AliEmcalCorrectionComponent``. For an example, see ``Initialize()`` in any of the correction components.

For your component to be accessible to the correction task, it must be registered within the framework.
To do so, add the following to your component header:

```{.cxx}
// // Allows the registration of the class so that it is availble to be used by the correction task.
static RegisterCorrectionComponent<AliEmcalCorrectionYourNewComponent> reg;
```

and then the following to your component implementation file (cxx):

```{.cxx}
// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionYourNewComponent> AliEmcalCorrectionYourNewComponent::reg("AliEmcalCorrectionYourNewComponent");
```

Lastly, the task needs to be added to the default %YAML configuration file. To do so, add the name of the
component (removing "AliEmcalCorrection") to the file, and then set default values for each parameter.
Be certain to document what each parameter means! In addition to the component configuration, be certain
to add it into the executionOrder node! Otherwise, your task will never be executed! Note that the order
of your task name in this list will determine when it is executed relative to the other components.

# Additional Info

You can find an analysis tutorial presentation [here](https://indico.cern.ch/event/586577/).

For more information, contact the developers <james.mulligan@yale.edu> and <raymond.ehlers@yale.edu>.
 
*/
