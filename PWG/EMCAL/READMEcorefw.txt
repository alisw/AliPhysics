/*! \page READMEcorefw The EMCAL core framework

The EMCAL framework consists of 

- An analysis task AliAnalysisTaskEmcal extending the functionality of the standard analysis task

- A set of containers inheriting from AliEmcalContainer
  - AliParticleContainer
  - AliClusterContainer
  - AliTrackContainer
  
- A virtual track selection AliEmcal track selection with their implementations
  - AliEmcalTrackSelectionESD for ESD analysis
  - AliEmcalTrackSelectionAOD for AOD analysis
  
# Writing an analysis using the AliAnalysisTaskEmcal   {#EMCALAnalysisTask}

The core of the EMCAL core framework is the AliAnalysisTaskEmcal. It should be
used for EMCAL or jet related analyses. The AliAnalysisTaskEmcal

- Handles cluster-/particle-/track-containers (see below)
- Performs event selection
- Fills QA histograms
- Runs the user code

In order to run the user code, users have to implement either of the
two functions:

- bool FillHistograms()
- bool Run()

The Run function is indended to contain the actual user analysis, while the FillHistograms
function should fill user histograms (for example on selected particles or jets). Input for
the analysis comes from the different containers attached to the task, connecting to data
from the input event. Users must not overwrite the function UserExec as it contains further
functionality essential for the AliAnalysisTaskEmcal to work.

In case the AliAnalysisTaskEmcal should handle the output as well, which is the normal 
behaviour, the task should be configured to also fill histograms. In this case several
general histograms monitoring the event selection are filled as well, among them the 
vertex distribution before and after selection, but also Monte-Carlo related information
like the cross section or the number of trials from the event generator. Furthermore
the list fOutput can handle all the user histograms. For this the function 
UserCreateOutputObjects needs to be implemented by the user and the corresponding function
of the base class has to be called at the beginning.

In the Add macro use the named constructor in order to setup the task and attach it to
the analysis manager.

# Using EMCAL containers              {#EMCALContainers}

EMCAL containers are used to handle arrays of objects (particles, clusters, jets) shared
among different tasks within the input ESD or AOD event via an easy interface hiding the
direct access from the user. Each task can handle several containers as needed. Once they
are defined, the AliAnalysisTaskEmcal will automatically connect them to the content they
are supposed to handle, so the content can be used from that moment on. Normally EMCAL 
containers are defined in the add macro of the task.

To better understand containers, consider the specific example of a track container:

- Create a container in your ``AddTask`` macro using the name of the object of interest:
 
  ~~~{.cxx}
  // The name of the tracks provided by the user (in this case, "tracks") must match the name in the input event!
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

- The containers are automatically loaded and ready to use in the ``Run()`` method in the user task. All you need to do is retrieve the tracks and iterate through the available tracks subject to the specified cuts:

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
  
  Users can also access all tracks using ``AliEmcalContainer::all_momentum()``, or just the ``AliVTrack`` object with ``AliEmcalContainer::accepted()`` and ``AliEmcalContainer::all()``. If you access the objects directly, be careful! The object will not necessarily give you the values! To get the properly corrected values, either use the ``AliTLorentzVector`` or explicitly call the corrected value (for example, ``AliVCluster::GetHadCorrEnergy()``). For more details on this issue, as well as more generally on iteration techniques, see ``AliEmcalIterableContainer``.

# Track Selection

Track selection is normally applied in the track selection task as part of the 
framework. For more information, including on virtual track selection, see the
page on \subpage READMEtracks

# EMCal/DCal Cluster Corrections

Information about the cell and cluster corrections is availabe on the \subpage READMEclustcorr page.
 
*/
