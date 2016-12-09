/*! \page READMEjetfw Jet finding framework

# Jet Framework Users

For users more familiar with the framework, see the [Jet Framework Topics](\ref jetFrameworkTopics) section below.

For those who are less familiar, please see the introduction [below](\ref jetFrameworkIntroduction).

An example can be found at ``<your-path-to-AliPhysics>/PWGJE/EMCalJetTasks/macros/runEMCalJetSampleTask.C``

# Introduction to the Jet Framework              {#jetFrameworkIntroduction}

## Framework Philosophy

The philosophy of the Jet Framework is to provide access to common objects and functions, thereby allowing the user to focus on the particulars of their task, as well as reducing duplicate boilerplate code from user tasks, helping to reduce the inevitable bugs introduced by many copies. To achieve these goals, the Jet Framework uses the EMCal Framework defined in ``$ALICE_PHYSICS/PWG/EMCAL/``.

The EMCal framework extends the standard ALICE %Analysis Framework, ``AliAnalysisTaskSE`` to provide commonly used objects, such as access to clusters and tracks. In the standard framework, the user is responsible for tasks such as loading the event properly. With the EMCal framework, such tasks are already taken care of for the user, along with providing basic QA histograms and other variables to simplify tasks.

The Jet Framework extends this concept, introducing support for collections of jets in ``AliAnalysisTaskEmcalJet``. Note that despite the name including EMCal, the framework works equally well for both charged and full (charged + neutral) jets!

The EMCAL framework base class is `AliAnalysisTaskEmcal`. It contains the cluster and track collections that are used by the containers. The Jet framework base class is ``AliAnalysisTaskEmcalJet`` and it provides the objects for jet containers. Since these are used across multiple PWGs, they are located in the ``EMCAL`` and ``JETFW`` folders, respectively, in ``$ALICE_PHYSICS/PWG``.

To use the Jet Framework, the user implements an analysis task, usually named ``%AliAnalysisTaskEmcalJet{AnalysisName}``, which includes a main implementation file (.cxx) and a header file (.h). In addition, the user implements a macro to call their analysis task with the proper options, usually named ``AddTaskEmcalJet{AnalysisName}``. Lastly, the user must have a way to call the macro - this is usually done either via a run macro (which could be used to submit to the grid), or via a wagon on a LEGO train. Frequently, the analysis task is placed in ``$ALICE_PHYSICS/PWGJE/EMCALJetTasks/UserTasks/`` or ``$ALICE_PHYSICS/PWGJE/FlavourJetTasks/``, while the run macros is placed in ``$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/`` or ``$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/``, respectively.

### Containers

Containers are a central concept in the EMCal framework. To conduct an analysis, a user needs access to objects such as clusters, tracks, or jets. Containers are wrappers around these collections, allowing a consistent interface regardless of the underlying object. Furthermore, the containers help apply and manage cuts and selections.

To learn about containers, including a usage example, looking the EMCal framework page on containers, available [here](\ref READMEcontainers).

For more information on the containers, see the base class, ``AliEmcalContainer``, as well as the particular containers, ``AliClusterContainer``, ``AliParticleContainer``, ``AliTrackContainer``, ``AliMCParticleContainer``, and ``AliJetContainer``.

### Jet Finding

Since Jet Finding is such an integral part of the framework, it is taken care of with an straightforward and simple interface. To perform jet finding, the user specifies the tracks and clusters, as well as the algorithm, R parameter, jet type (charged or full), and a few additional options. These options are passed to ``AddTaskEmcalJet``, and it performs all necessary actions for the user, including loading the jets into the jet container and making them available to the user. This means that adding one task to a run macro can take care of most basic jet finding needs.

While basic jet finding is straightforward, this does not preclude additional tasks. For example, using another macro (``AddTaskRhoNew``), the average background can be calculated. The actual implementation of background subtraction is performed via another class, with information available [here](\ref READMEjetfwUtilities). In general, the user will utilize various tasks that depend on functionality implemented in ``fastjet`` (via ``AliFJWrapper``). However, a user will rarely, if ever, need to modify the wrapper! (Only in the case of adding new functionality).

For further information on jet finding, see ``AliEmcalJetTask``. For further information on background calculation, see ``AliAnalysisTaskRho``. In general, further jet documentation can be viewed at the [ALICE FJ Utilities Wrapper documentation](\ref READMEjetfwUtilities). For more advanced jet finding details and features, see the interface to ``fastjet``, ``AliFJWrapper``.

## %Analysis Chain

Using the philosophy, each user can create their own analysis task, and then run it either using a Run Macro (locally to test, and then on the grid), or create a wagon to run it on the LEGO train.

### Your %Analysis Task

Inherit from ``AliAnalysisTaskEmcalJet`` for jet focused analyses. If jet features are not required, one may instead inherit from ``AliAnalysisTaskEmcal``.

Both of the base classes inherit from ``AliAnalysisTaskSE``, but additional functionality is implemented (such as automatic loading of data into containers for each event), so the user functions change slightly. Instead of implementing ``UserExec()``, the user should implement ``Run()`` and ``FillHistograms()``.

For more information, look at the base classes, ``AliAnalysisTaskEmcal`` and ``AliAnalysisTaskEmcalJet``, as well as the sample task, ``AliAnalysisTaskEmcalJetSample``.

### Your Add Task

The ``AddTask`` is used to setup your analysis task, allowing it to run. This code is interpreted through CINT (ie through ROOT), so it is best to avoid complicated code, as debugging can be rather difficult.

Often, the ``AddTask`` is used to set properties such as the tracks or clusters. To handle this in a flexible manner, it is recommended to utilize the "usedefault" approach as demonstrated in many ``AddTask`` macros, such as ``AddTaskEmcalJet``. Doing so ensures that the proper collection of objects is loaded regardless of the type of file that is used. Alternatively, the containers can be setup in your run macro or wagon. To see the most up to date names for these collections, see ``AddTaskEmcalJetSample``.

For examples, see ``AddTaskEmcalJetSample`` and ``AddTaskEmcalJet``.

### Running Your task

Once your have created your task, the next step is to run it. There are a few different options:

- Attaching a wagon to the LEGO train. For more information, see the <a href="https://alimonitor.cern.ch/trains/">LEGO Trains</a> (requires a grid certificate) and the <a href="https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AnalysisTrains">LEGO Train TWiki Page</a>.
- Run via a Run Macro. For an example, see ``runEMCalJetSampleTask.C``. It runs the primary tasks necessary for running an analysis, so simple modification of it should be sufficient to get started. To try running it, ``runEMCalJetSampleTask.sh`` provides a slightly more use friendly interface.

  A Run Macro enables two modes of operation:
    - Run locally to test code
    - Submit to the grid

  NOTE: Since you will need to modify the run macro, it is highly recommended to copy the run macro to another location or to make a new branch for it!

Regardless of method, there are a number of considerations to take into account. A few are listed below.

#### Physics Selection

When running your task, it is necessary to specify the physics data that you are interested in analyzing. This should be specified using ``SelectCollisionCandidate()`` in your run macro. The function is defined in ``AliAnalysisTaskSE``.

Note that you only need a physics selection task when using ESDs. In that case, as per Salvatore's email on Feb 18, 2016, use the physics selection task located at ``OADB/macros/AddTaskPhysicsSelection.C`` for the old centraltiy framework, or ``OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C`` for the new centrality framework.

For more information, see ``AliAnalysisTaskSE`` and the possible selections in ``AliVEvent``. For an example of an implementation, see ``runEMCalJetAnalysisNew.C``.

#### Necessary Tasks

There are a number of tasks which are required to run before running your analysis in the proper order. These tasks include the EMCal corrections (cell corrections, cluster corrections, hadronic corrections - see [here](\ref READMEclustcorr)) and jet finding. For an up to date task, see ``<your-path-to-AliPhysics>/PWGJE/EMCalJetTasks/macros/runEMCalJetSampleTask.C``.

### Additional notes

A lot of event properties are also automatically available when you derive from ``AliAnalysisTaskEmcal``. For example, the centrality percentile of the current event is available in ``fCent``, and the bin as ``fCentBin`` (see the source for the precise binning options). Note that if you are using the new centrality framework (OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C), you need to set task->SetUseNewCentralityEstimation(kTRUE) 

Other additional features include:

- Calculating \f$\Delta\phi\f$ in any range using ``DeltaPhi(elementOne, elementTwo, minRange, maxRange)``. 
- Calculating \f$\Delta\eta,\Delta\phi\f$ for a track and cluster.
- Generating a fixed bin array.

In general, it is extremely helpful to be familiar with the base classes. It can save you a tremendous amount of work!

---------------------

# Jet Framework Topics              {#jetFrameworkTopics}

## Before looking for jets: Track and Cluster selection

Before running any jet finder a proper track and cluster samples have to be set/prepared for use.

To select tracks you can find instructions [here] (\ref READMEtracks).

Clusters require a bit more levels of corrections, see [here] (\ref READMEclustcorr). The different stages of cluster corrections are accessible via different data members of AliVCluster in your task.

## Basic jet finding

Basic jet finding is provided by the AliEmcalJetTask class which is found in the library libPWGJEEMCALJetTasks (source code in PWGJE/EMCALJetTasks). An add task macro is provided in PWGJE/EMCALJetTasks/AddTaskEmcalJet.C:

~~~{.cxx}
AliEmcalJetTask* AddTaskEmcalJet(
  const char *nTracks        = "usedefault",
  const char *nClusters      = "usedefault",
  const Int_t algo           = AliJetContainer::antikt_algorithm,
  const Double_t radius      = 0.4,
  const Int_t type           = AliJetContainer::kFullJet,
  const Double_t minTrPt     = 0.15,
  const Double_t minClPt     = 0.30,
  const Double_t ghostArea   = 0.005,
  const Int_t recombScheme   = AliJetContainer::pt_scheme,
  const char *tag            = "Jet",
  const Double_t minJetPt    = 0.,
  const Bool_t lockTask       = kTRUE,
  const Bool_t bFillGhosts    = kFALSE
)
~~~

If you use "usedefault" for nTrack and nCluster the AliParticleContainer and AliClusterContainer will be created using the default names, namely "tracks" and "caloClusters" for AOD and "Tracks" and "CaloClusters" for ESD.

### Charged jets
For charged jets, take care of your [track selection] (\ref READMEtracks) and run the jet finder e.g. like this:
~~~
AliEmcalJetTask* jetFinderCh = AliEmcalJetTask("usedefault", "",  AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kChargedJet)
~~~

### Full jets
EMCal/DCal cluster corrections have to be applied beforehand as explained [here](\ref READMEclustcorr).

EMCal/DCal cluster objects contain different fields to accomodate different "levels" of corrections to the energy deposition. **The user must be careful in selecting the level of corrections required for his/her analysis**. The various corrected energies are selected via
~~~{.cxx}
pFuJetTask->GetClusterContainer(0)->SetDefaultClusterEnergy(energyType);
~~~
where pFuJetTask is the pointer to the AliEmcalJetTask object and energyType is an integer number. It can be either -1 (no additional corrections to the cluster energy) or one of the enum constants defined in AliVCluster:

~~~{.cxx}
enum VCluUserDefEnergy_t {
    kNonLinCorr          = 0,
    kHadCorr             = 1,
    kUserDefEnergy1      = 2,
    kUserDefEnergy2      = 3,
    kLastUserDefEnergy   = 4
  };
~~~

### Particle level jets (MC)
For particle level jets it is usually enough to filter primary particles (see \ref READMEtracks).

## Jet containers

For an introduction to containers, see [here](\ref READMEcontainers). The jet container ``AliJetContainer`` allows you to apply a variety of basic cuts to your jet collection.

### Acceptance cut

For example, you can set the geometrical jet acceptance selection you would like to consider -- the allowed options are listed in AliEmcalJet::JetAcceptanceType. The user can select a single type (e.g. kEMCAL), or a bitwise combination (e.g. kEMCAL | kDCAL). The container can be configured via AliJetContainer::SetJetAcceptanceType() or when adding a jet container via one of the AliAnalysisTaskEmcalJet::AddJetContainer functions. The cut is implemented in AliJetContainer by comparing a jet's bits (set automatically in the jet finder) to the container's bits (set by user).

Additionally, since the jet acceptance bits are stored in the jet, you can manipulate them in your analysis using the AliEmcalJet::GetJetAcceptanceType() function. This may be useful when studying jets in the DCal region, since there are several different partially overlapping bits.

## Utilities (e.g. FJ contribs)

For information on utilities such as ``fastjet`` contrib, see \subpage READMEjetfwUtilities.

## Embedding
Embedding here means to combine two events at the level of reconstructed tracks and EMCal cells or random tracks or clusters.

For more information, see the (new) [Embedding Framework](\ref READMEemcEmbedding).

There is also an old embedding framework which uses an older version of the framework that requires manual filtering of the tracks. See \subpage READMEembedding for more information on the old embedding classes.

## Jet tagger
The task AliAnalysisTaskEmcalJetTagger allows to tag a jets as "close" and/or sharing a minimun fraction of constituent p<sub>T</sub>. This is useful for:
-# Matching PYTHIA jet before and after embedding
-# Match jets with "signal" jets. A signal jet can be defined e.g. as a jet with a minimum cut p<sub>T</sub> > 4 GeV/c on the constituents

See \subpage READMEtagging for more details.

## Unfolding
\subpage READMEunfolding

## User tasks
\subpage READMEJEtasks

## Jet QA				{#jetQA}
\subpage READMEjetQA

### Appendix

For clarifications, corrections or improvements, please contact [Raymond Ehlers](mailto:raymond.ehlers@cern.ch).

*/
