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
*/