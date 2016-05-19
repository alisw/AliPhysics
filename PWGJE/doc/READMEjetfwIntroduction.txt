/*! \page READMEjetfwIntroduction Introduction to the Jet Framework

# Introduction to the Jet Framework

## Framework Philosophy

The philosophy of the Jet Framework is to provide access to common objects and functions, thereby allowing the user to focus on the particulars of their task, as well as reducing duplicate boilerplate code from user tasks, helping to reduce the inevitable bugs introduced by many copies. To achieve these goals, the Jet Framework uses the EMCal Framework defined in ``$ALICE_PHYSICS/PWG/EMCAL/``.

The EMCal framework extends the standard ALICE %Analysis Framework, ``AliAnalysisTaskSE`` to provide commonly used objects, such as access to clusters and tracks. In the standard framework, the user is responsible for tasks such as loading the event properly. With the EMCal framework, such tasks are already taken care of for the user, along with providing basic QA histograms and other variables to simplify tasks.

The Jet Framework extends this concept, introducing support for collections of jets in ``AliAnalysisTaskEmcalJet``. Note that despite the name including EMCal, the framework works equally well for both charged and full (charged + neutral) jets!

The EMCAL framework base class is `AliAnalysisTaskEmcal`. It contains the cluster and track collections that are used by the containers. The Jet framework base class is ``AliAnalysisTaskEmcalJet`` and it provides the objects for jet containers. Since these are used across multiple PWGs, they are located in the ``EMCAL`` and ``JETFW`` folders, respectively, in ``$ALICE_PHYSICS/PWG``.

To use the Jet Framework, the user implements an analysis task, usually named ``%AliAnalysisTaskEmcalJet{AnalysisName}``, which includes a main implementation file (.cxx) and a header file (.h). In addition, the user implements a macro to call their analysis task with the proper options, usually named ``AddTaskEmcalJet{AnalysisName}``. Lastly, the user must have a way to call the macro - this is usually done either via a run macro (which could be used to submit to the grid), or via a wagon on a LEGO train. Frequently, the analysis task is placed in ``$ALICE_PHYSICS/PWGJE/EMCALJetTasks/UserTasks/`` or ``$ALICE_PHYSICS/PWGJE/FlavourJetTasks/``, while the run macros is placed in ``$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/`` or ``$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/``, respectively.

### Containers

Containers are a central concept in the EMCal framework. To conduct an analysis, a user needs access to objects such as clusters, tracks, or jets. Containers are wrappers around these collections, allowing a consistent interface regardless of the underlying object. Furthermore, the containers help apply and manage cuts and selections.

Usage of a container can be best summarized using the specific example of accessing tracks:

 - Create a container in your ``AddTask`` macro using the name of the object of interest:

 ~~~{.cxx}
 AliTrackContainer * exampleTracks = new AliTrackContainer("tracks", "myTrackContainer");
 // Make the tracks available in your task which inherits from AliAnalysisTaskEmcal
 myTask->AdoptParticleContainer(exampleTracks);
 // Can alternatively be created using AddTrackContainer(). See below.
 ~~~
 
 Note that the user may create as many containers as desired, and each one will be treated independently. In general, the container interface is extremely flexible, allowing the user to use the contained objects however is desired. As an alternative, the user can create the container and add it to the task in one step by using AliAnalysisTaskEmcal::AddTrackContainer(). Either approach yields equivalent results.
 
 - A user may also apply any desired cuts:
 
 ~~~{.cxx}
 exampleTracks->SetMinPt(2);
 exampleTracks->SetEtaLimits(-0.9, 0.9);
 // ...
 ~~~
 
 A wide variety of cuts are available, reducing code duplication and allowing the analyzer to focus on their analysis. One may also apply track selection. See \ref READMEtracks for more information.

 - The containers are automatically loaded and ready to use in the ``Run()`` method in your task. All you need to do is retrieve the tracks and iterate through the available tracks subject to the specified cuts:

 ~~~{.cxx}
 // Retrieve the tracks
 // Can also retrieve tracks by name
 AliTrackContainer * tracks = dynamic_cast<AliTrackContainer *>(GetParticleContainer("myTrackContainer"));
 
 // Iterable approach (using C++11)
 AliTLorentzVector track;
 for (auto trackIterator : tracks->accepted_momentum() )
 {
     // trackIterator is a std::map of AliTLorentzVector and AliVTrack
     // Get the proper track kinematics
     track.Clear();
     track = trackIterator.first;
     // Full access to the full track is also available with:
     AliVTrack * fullTrack = trackIterator.second;
     // However, you need to be careful with this object to ensure that you get the proper values!
     // See the note below.
 }
 ~~~
 
 Users can also access all tracks using AliEmcalContainer::all_momentum(), or just the ``AliVTrack`` object with AliEmcalContainer::accepted() and AliEmcalContainer::all(). If you access the objects directly, be careful! The object will not necessarily give you the values! To get the properly corrected values, either use the ``AliTLorentzVector`` or explicitly call the corrected value (for example, AliVCluster::GetHadCorrEnergy()). For more details on this issue, as well as more generally on iteration techniques, see ``AliEmcalContainer``.

For more information on the containers, see the base class, ``AliEmcalContainer``, as well as the particular containers, ``AliClusterContainer``, ``AliParticleContainer``, ``AliTrackContainer``, ``AliMCParticleContainer``, and ``AliJetContainer``.

### Jet Finding

Since Jet Finding is such an integral part of the framework, it is taken care of with an straightforward and simple interface. To perform jet finding, the user specifies the tracks and clusters, as well as the algorithm, R parameter, jet type (charged or full), and a few additional options. These options are passed to ``AddTaskEmcalJet``, and it performs all necessary actions for the user, including loading the jets into the jet container and making them available to the user. This means that adding one task to a run macro can take care of most basic jet finding needs.

While basic jet finding is straightforward, this does not preclude additional tasks. For example, using another macro (``AddTaskRhoNew``), the average background can be calculated. The actual implementation of background subtraction is performed via another class, with information available [here](\ref fjContribUtilities). In general, the user will utilize various tasks that depend on functionality implemented in ``fastjet`` (via ``AliFJWrapper``). However, a user will rarely, if ever, need to modify the wrapper! (Only in the case of adding new functionality).

For further information on jet finding, see ``AliEmcalJetTask``. For further information on background calculation, see ``AliAnalysisTaskRho``. In general, further jet documentation can be viewed at the [ALICE FJ Utilities Wrapper documentation](\ref fjContribUtilities). For more advanced jet finding details and features, see the interface to ``fastjet``, ``AliFJWrapper``.

## %Analysis Chain

Using the philosophy, each user can create their own analysis task, and then run it either using a Run Macro (locally to test, and then on the grid), or create a wagon to run it on the LEGO train.

### Your %Analysis Task

Inherit from ``AliAnalysisTaskEmcalJet`` for jet focused analyses. If jet features are not required, one may instead inherit from ``AliAnalysisTaskEmcal``.

Both of the base classes inherit from ``AliAnalysisTaskSE``, but additional functionality is implemented (such as automatic loading of data into containers for each event), so the user functions change slightly. Instead of implementing ``UserExec()``, the user should implement ``Run()`` and ``FillHistograms()``.

For more information, look at the base classes, ``AliAnalysisTaskEmcal`` and ``AliAnalysisTaskEmcalJet``, as well as the sample task, ``AliAnalysisTaskEmcalJetSample``.

### Your Add Task

The ``AddTask`` is used to setup your analysis task, allowing it to run. This code is interpreted through CINT (ie through ROOT), so it is best to avoid complicated code, as debugging can be rather difficult.

Often, the ``AddTask`` is used to set properties such as the tracks or clusters. To handle this in a flexible manner, it is recommended to utilize the "usedefault" approach as demonstrated in many ``AddTask`` macros, such as ``AddTaskEmcalJet``. Doing so ensures that the proper collection of objects is loaded regardless of the type of file that is used. Alternatively, the containers can be setup in your run macro or wagon. To see the most up to date names for these collections, see ``AddTaskEmcalJetSample``.

For examples, see ``AddTaskEmcalJetSample`` and ``AddTaskEmcalJet``

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

Note that you only need a physics selection task when using ESDs. In that case, as per Salvatore's email on Feb 18, 2016, use the physics selection task located at ``OADB/macros/AddTaskPhysicsSelection.C``

For more information, see ``AliAnalysisTaskSE`` and the possible selections in ``AliVEvent``. For an example of an implementation, see ``runEMCalJetAnalysisNew.C``.

#### Necessary Tasks

There are a number of tasks which are required to run before running your analysis in the proper order. These tasks include the EMCal corrections (cell corrections, cluster corrections, hadronic corrections) and jet finding. For an up to date task, see ``runEMCalJetAnalysisNew.C``.

## Appendix

For clarifications or corrections, please contact [raymond.ehlers@yale.edu](mailto:raymond.ehlers@yale.edu)

*/
