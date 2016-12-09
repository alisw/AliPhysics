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

The Run function is intended to contain the actual user analysis, while the FillHistograms
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

# EMCal Containers

Containers are a central concept in the EMCal framework. To conduct an analysis, a user needs access to objects such as clusters, tracks, or jets. Containers are wrappers around these collections, allowing a consistent interface regardless of the underlying object. Furthermore, the containers help apply and manage cuts and selections. For some basics of working with containers, see \subpage READMEcontainers. 

# Track Selection

Track selection is normally applied in the track selection task as part of the 
framework. For more information, including on virtual track selection, see \subpage READMEtracks.

# EMCal Corrections Framework

Information about the corrections framework is available at \subpage READMEemcCorrections.

# EMCal/DCal Cluster Corrections

Information about the cell and cluster corrections is available at \subpage READMEclustcorr.

# EMCal Embedding Framework

Information about embedding in the EMCal framework is available at \subpage READMEemcEmbedding.

# How to transition from old to new framework

Here you find detailed information which compares old and new framework and tells you the things you have to keep in mind when using the new EMCal framework. Please read \subpage READMEchangefw.  
*/
