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
  
\section EMCALAnalysisTask Writing an analysis using the AliAnalysisTaskEmcal

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

The Run function is indended to contain the actual user Analysis, while the FillHistograms
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

\section EMCALcontainers Using the EMCAL containers

EMCAL containers are used to handle arrays of objects (particles, clusters, jets) shared
among different tasks within the input ESD or AOD event via an easy interface hiding the
direct access from the user. Each task can handle several containers as needed. Once they
are defined, the AliAnalysisTaskEmcal will automatically connect them to the content they
are supposed to handle, so the content can be used from that moment on. Normally EMCAL 
containers are defined in the add macro of the task.

- In the add task users tasks need to create a container which can connect to the array inside the 
  input event. It is important that the name provided matches the name of the content in the input event.
  <pre>
  AliParticleContainer *cont = task->AddParticleContainer("tracks");
  cont->SetName("MyTrackContainer");
  </pre>
  If one likes, one can apply also kinematical cuts
  <pre>
  ...
  cont->SetMinPt(0.5);
  cont->SetEtaRange(-0.8, 0.8);
  ...
  </pre> 
- In the user task one needs to get the container from the task. Afterwards one can iterate over the content
  <pre>
  ...
  AliParticleContainer *trackcont = GetParticleContainer("MyTrackContainer");
  trackcont->ResetCurrentID(-1);
  AliVParticle *track;
  while((track = trackcont->GetNextAcceptParticle())){
    // Do something with the particle
  }
  </pre>
  The loop is done only on particles surviving the kinematical selection.

\section VirtualTrackSelection Using the virtual track selection

Track selection is normally applied in the track selection task as part of the 
framework. In case a more strict track selection needs to be applied, users 
should use the virtual track selection. The virtual track selection provides
an interface handling the track selection in a seemless way both for ESDs and
AODs. The selection criteria have to be implemented by the user.

<pre>
AliEmcalTrackSelectionESD *sel = new AliEmcalTrackSelectionESD;
sel->AddTrackCuts(AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 1));
</pre>

Inside the task the virtual track selection can be used in two ways: Either one checks for
every track separate, via the function IsTrackAccepted, or one requests a TObjArray of all 
accepted tracks via GetAcceptedTracks. Everything needed at this step is provided by
AliEmcalTrackSelection, the virtual base class. 
*/