/*! \page READMEemcCorrections EMC Correction Framework

# EMC Correction Framework

The new EMCal correction framework centralizes all of the analysis-level EMCal corrections into a single task:
AliEmcalCorrectionTask. This task coordinates the running of a configurable list of individual "correction components", where a
correction component means a single correction to be performed: cell-level energy calibration, or cluster-level non-linearity
correction, or cluster-track matching, etc. The set of corrections to run, as well as their order and their configurable
parameters, is set in a YAML configuration file. A default config file containing the centralized standard set of parameters is
shared by all users, and each user additionally writes their own user config file, overwriting desired variables.

This new correction task unifies what was previously done by tasks such as:
 - EMCal Tender
 - ClusterMaker
 - ClusterizerFast
 - ClusTrackMatcher
 - HadCorr

...and other similar tasks.

The motivation of this new approach is the following:
 - Centralize code so that there are not many similar versions of code used by different people
 - Provide a simple and unified interface to make life easier for both users and developers
 - Reduce likelihood of mistakes for users (make it obvious to a user if they are deviating from a recommended config, and avoid setting arguments through ``AddTask`` macros)

# Prerequisites

Before using the EMCal correction framework, it is extremely important that your code is updated enough that it is using
EMCal containers! In particular, you need to be using the "new" EMCal framework developed by Salvatore.
See [here](\ref READMEchangefw) for instructions to update, if you need.
There are a few special (and rare) exceptions - contact the developers for further information if you think this
applies to your case.

This matters in particular because we always write cell/cluster/track corrections in place - it is generally not possible to write out new collections.

# Using the EMCal Corrections Framework

Below, follow the instructions to implement the new Correction Task in your analysis. Further down the page, you will find
instructions on how to transition from the old correction method(s) to the new correction method, including how to configure the
corrections, and how to test your analysis results in the new vs. old correction method.

# %Setup the Correction Task 

There are two steps:
- [Configure the Correction Add Task](\ref configureEMCalCorrectionRunMacro)
- [Configure the corrections](\ref configureEMCalCorrections)

# Configure the Correction Add Task (or train wagon)            {#configureEMCalCorrectionRunMacro}

To enable the correction task, add the following lines to your run macro:

~~~{.cxx}
AliEmcalCorrectionTask * correctionTask = AddTaskEmcalCorrectionTask();
correctionTask->SelectCollisionCandidates(kPhysSel);
// Set the run period, same as the track container
// If you derived from the file "runEMCalJetSampleTask.C", then it is likely stored under "sRunPeriod.Data()"
correctionTask->SetRunPeriod("LHC11h");
// Set the user configuration file, assuming that your file is called "userConfiguration.yaml" and is located in
// the current directory. This also supports alien:// paths!
correctionTask->SetUserConfigurationFilename("userConfiguration.yaml");
// Initialize the configuration files in the correction task
// It is EXTREMELY important to run this function in your run macro!
correctionTask->InitializeConfiguration();
~~~

Don't forget to also load the macro with:

~~~{.cxx}
gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalCorrectionTask.C");
~~~

# Configuring Corrections                                      {#configureEMCalCorrections}

The corrections configuration file is specified via a file written in the YAML markup language. YAML is quite
readable and is commonly used for configuration files. Although it should be straightforward to read and write,
more information on its use can be found through the examples available [here](https://en.wikipedia.org/wiki/YAML). **All**
configuration and options are set through the YAML file, from setting the clusterizer type to the name of the tracks!

The configuration is determined by two files: the base default file (located in ``$ALICE_PHYSICS/PWG/EMCAL/config/``),
as well as the user configuration file. Any changes that you need to make should be done in the user configuration file!
Note that any setting that you have in the user file will override the default file!

## Introducing features of the configuration file

We will discuss a number of features including:
- Comments
- Shared parameters
- Setting cells, clusters, and tracks
- Enumerations

To discuss the configuration and introduce, consider the configuration excerpt below as an example:

~~~
# This is a comment!
# The next line sets the cell branch name - note that we support the "usedefault" pattern of setting containers branch names.
cellBranchName: "usedefault"   # Anything after a "hash" symbol is a comment.
# This defines the "sharedParamters" section.
sharedParameters:
    # These are parameters shared by multiple correction components.
    clusterBranchName: "usedefault"
    trackBranchName: "usedefault"
# This defines the settings for the "AliEmcalCorrectionCellEnergy" correction.
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
    clusterizer: kClusterizerv2
    cellE: 0.05
    seedE: 0.1
    cellTimeMin: -1                    # Min cell time (s)
    cellTimeMax: +1                    # Max cell time (s)
    clusterTimeLength: 1               # Maximum time difference between the digits inside EMC cluster (s)
    w0: 4.5
~~~

### Comments

Comments include anything after a hash ("#"). Of course, they are not required, but are highly recommended! You can add a
comment on the same line as where you define a property - just add the "#" after you are done defining the property.

### Shared properties

If you have some parameters with the same value and you want to change them in unison (say, if you wanted to change the
track name, or a min pt cut), you can set a parameter in the ``sharedParameters`` section, and then set the value in each
component to that parameter. For example, if I wanted to set the min pt in the ``ClusterExotics`` and ``ClusterNonLinearity``,
I would have:

~~~
sharedParameters:
    ptMin: 0.15
ClusterExotics:
    clusterPtMin: "sharedParameters:ptMin"
ClusterNonLinearity:
    clusterPtMin: "sharedParameters:ptMin"
~~~

This is just provided for your convenience!

### Enumerations

Enumerations are supported as expected. Just set the value as you normally would. However, be sure that you don't include
class where it is defined!

### Example configuration and further information

For a survery of the available configuration options and information about the meaning of each, see the default configuration
file located in ``$ALICE_PHYSICS/PWG/EMCAL/config/``. This also serves as an example configuration, as it contains all possible
types of settings.

# Switching to the EMCal Corrections Framework

For those who are switching to this framework, a special page has been prepared to explain the process.
This page explains how to configure the correction framework alongside the previous set of corrections, allowing
you to test and show that you get the same results with the new and old corrections!
For the instructions, please see \subpage READMEemcCorrectionsChange.

# Details on the framework and the corrections

You can see the code at ``$ALICE_PHYSICS/PWG/EMCAL/EMCALtasks``. The steering class is ``AliEmcalCorrectionTask``. The
individual corrections inherit from ``AliEmcalCorrectionComponent``, and are labeled ``AliEmcalCorrectionXXXX``. Note that
neither ``AliEmcalCorrectionTask`` nor ``AliEmcalCorrectionComponent`` inherit from ``AliAnalysisTaskEmcal``. However, they
provide similar functionality. The default configuration file is at ``$ALICE_PHYSICS/PWG/EMCAL/config/AliEmcalConfiguration.yaml``.

NOTE: If you are interested in how a particular correction works, you only need to look at the particular correction and its configuration!
There are many other details in the base and steering classes, but they are almost certainly not relevant!

# Developing a correction

If you are interested in developing a task that is shared by analyses using the EMCal, then the correction framework is a great place to deploy it!
The general approach is very similar to the steering in ``AliAnalysisTaskEmcal``. Roughly, your task should derive from ``AliEmcalCorrectionComponent``
and implement:

~~~{.cxx}
// Run once to initialize the component. These are run independent initializations.
virtual Bool_t Initialize();
// Execute in the first event for run dependent initialization. Same as in AliAnalysisTaskEmcal
virtual void ExecOnce();
// Called each event. Same as in AliAnalysisTaskEmcal
virtual Bool_t Run();
// UserNotify() from AliAnalysisTaskSE
virtual Bool_t UserNotify();
~~~

Until this section is expanded further, for more information, contact the developers <james.mulligan@yale.edu> and <raymond.ehlers@yale.edu>.
 
*/
