/*! \page READMEjetfw The jet finding framework

# Basic jet finding

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

## Charged jets
For charged jets, running the AliEmcalJetTask is sufficient. The track selection is described in \subpage READMEtracks.

## Full jets
EMCal/DCal cluster corrections have to be applied beforehand as explained here \subpage READMEclustcorr.

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

## Particle level jets (MC)
For particle level jets it is usually enough to filter primary particles (see \subpage READMEtracks).

# Utilities (e.g. FJ contribs)
Additional utilities can be attached to the AliEmcalJetTask object, in a similar fashion as it is done for the AliTender class. The utility classes have to derive from the abstract class AliEmcalJetUtility. An EMCal jet utility class can implement any of the following four methods declared as virtual in AliEmcalJetUtility:

~~~{.cxx}
virtual void Init() = 0;                                                        // Executed only once at the end of AliEmcalJetTask::DoInit()
virtual void Prepare(AliFJWrapper& fjw) = 0;                                    // Executed for each event at the beginning of AliEmcalJetTask::FillJetBranch()
virtual void ProcessJet(AliEmcalJet* jet, Int_t ij, AliFJWrapper& fjw) = 0;     // Executed for each jet in the loop in AliEmcalJetTask::FillJetBranch()
virtual void Terminate(AliFJWrapper& fjw) = 0;                                  // Executed for each event at the end of AliEmcalJetTask::FillJetBranch()
~~~

At the moment two utilities are available, which make use of the FastJet contribs: the generic subtractor implemented in AliEmcalJetUtilityGenSubtractor and the constituent subtractor implemented in AliEmcalJetUtilityConstSubtractor. For example, to add the generic subtractor to a previously defined AliEmcalJetTask object named jetTask one can use the following code.

~~~{.cxx}
AliEmcalJetUtilityGenSubtractor* genSub = jetTask->AddUtility(new AliEmcalJetUtilityGenSubtractor("GenSubtractor"));
genSub->SetUseExternalBkg(kTRUE);
genSub->SetRhoName(rhoname);
genSub->SetRhomName(rhomname);
genSub->SetGenericSubtractionJetMass(kTRUE);
~~~

or for the constituent subtraction:

~~~{.cxx}
AliEmcalJetUtilityConstSubtractor* constUtil = jetTask->AddUtility(new AliEmcalJetUtilityConstSubtractor("ConstSubtractor"));
constUtil->SetUseExternalBkg(kTRUE);
constUtil->SetRhoName(rhoname);
constUtil->SetRhomName(rhomname);
constUtil->SetJetsSubName(Form("%sConstSub",jetTask->GetName()));
constUtil->SetParticlesSubName(Form("%s_%sConstSub",kTracksName.Data(),jetTask->GetName()));
~~~

_Note:_ A description of the AliFWWrapper should be added

# User task

To write your first "EMCal framework task", start from AliAnalysisTaskEmcalJetSpectraQA (in $ALICE_PHYSICS/PWGJE/EMCALJetTasks) and corresponding AddTaskEmcalJetSpectraQA (in $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros). If you do not need to access jets, you may use instead AliAnalysisTaskEmcalJetQA and corresponding AddTaskEmcalJetQA.C.

- Inherit from AliAnalysisTaskEmcalJet (if you access jets) or AliAnalysisTaskEmcal (only tracks and clusters)
- In the AddTask macro add as many containers as you need: 

~~~{.cxx}
AliTrackContainer* trackCont = task->AddTrackContainer("tracks"); // hybrid track (data and MC detector level)
AliParticleContainer* partCont = task->AddParticleContainer("PicoTracks"); // if you still use the old framework (e.g. embedding)
AliMCParticleContainer* mcPartCont = task->AddMCParticleContainer("mcparticles"); // MC particles at generator level (physical primaries are selected by default)
AliClusterContainer* clusCont = task->AddClusterContainer("caloClusters"); // EMCal/DCal clusters
AliJetContainer *jetCont = task->AddJetContainer(AliJetContainer::kFullJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, 
radius, AliJetContainer::kEMCALfid, trackCont, clusCont);
~~~

where `task` is a valid pointer to an analysis task object derived from AliAnalysisTaskEmcal or AliAnalysisTaskEmcalJet and radius is a `Double_t` variable with the resolution parameter of the jets.

- In your task get the containers: 

~~~{.cxx}
AliJetContainer *jetCont = GetJetContainer(0);
AliTrackContainer *trackCont = GetTrackContainer(0);
~~~

- Use the iterator methods in the containers to loop over objects (e.g. AliJetContainer::GetNextJet(), AliJetContainer::GetJet(i)) and to apply cuts (e.g. AliJetContainer::AcceptJet(i)).
- A lot of event properties are also automatically available when you derive from AliAnalysisTaskEmcal. For example, the centrality percentile of the current event is available in fCent. Note that if you are using the new centrality framework (OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C), you need to set task->SetUseNewCentralityEstimation(kTRUE) 

_Note_: $ALICE_PHYSICS/PWGJE/EMCALJetTasks/AliAnalysisTaskEmcalJetSample.h/cxx is currently outdated and should not be used as an example. Use instead the QA tasks as mentioned above.

# Embedding
Embedding here means to combine two events at the level of reconstructed tracks and EMCal cells or random tracks or clusters.

Currently the embedding framework uses an older version of the framework that requires manual filtering of the tracks as explained in **find doc and add here link**.

The embedding classes are described here \subpage READMEembedding

# Jet tagger
The task AliAnalysisTaskEmcalJetTagger allows to tag a jets as "close" and/or sharing a minimun fraction of constituent p<sub>T</sub>. This is useful for:
-# Matching PYTHIA jet before and after embedding
-# Match jets with "signal" jets. A signal jet can be defined e.g. as a jet with a minimum cut p<sub>T</sub> > 4 GeV/c on the constituents

See \subpage READMEtagging for more details.

# Unfolding
\subpage READMEunfolding

*/
