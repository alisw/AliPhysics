/*! \page READMEembedding Old Embedding Framework

__IMPORTANT__: This page describes the old embedding framework and has been superseded by the new EMCal Embedding Framework! Please see documentation on the [new framework](\ref READMEemcEmbedding).

# Embedding

__Page under construction__

The basic functionalities are in the class AliJetModelBaseTask, in particular the methods to add tracks, clusters, or MC particles to the corresponding track or clusters array. Those methods are AliJetModelBaseTask::AddTrack, AliJetModelBaseTask::AddCluster, AliJetModelBaseTask::AddMCParticle. The embedded objects can be added to the original array or into a copy of it (`SetCopyArray(kTRUE)` and `SetSuffix(newName)`).

There are currently few types of embedding available:
-# Embedding of AOD events (AliJetEmbeddingFromAODTask and AliJetEmbeddingFromPYTHIATask with additional functionalities for our general productions, e.g. p<sub>T</sub> hard bins)
-# Embedding of generated PYTHIA events (AliJetEmbeddingFromGenTask)
-# Embedding of single particles (AliJetEmbeddingTask)
-# If you need specific features that are not implemented in the available classes, you can write your own class inheriting from AliJetModelBaseTask

_Note_: The Embedding framework in going to be restructured and updated for the use within the new jet framework. It is currently mandatory to use the old jet framwork since there's an explicit use of the AliPicoTrack in AddTrack.

_More info_: JIRA ticket for the new embedding development [ALPHY-53](https://alice.its.cern.ch/jira/browse/ALPHY-53)

## Embedding of reconstructed (detector level) PYTHIA

Use the AddTaskJetEmbeddingFromPYTHIA.C macro that runs the task AliJetEmbeddingFromPYTHIATask. Also particle level objects will be available. To distinguish the the embedded objects get the label of the track (`track->GetLabel()`), it'll be !=0, differently from data tracks. If you need to use the EMCal you have to run the clusterizer after embedding.

To perform the analysis of the embedded event proceed normally by running the jet finder(s). You probably want to run one on the total event and one on the PYTHIA only event:

- PYTHIA+data jets:

~~~{.cxx}
AliEmcalJetTask *jetTask= AddTaskEmcalJet("PicoTracks","", 1, 0.2, 1, 0.15, 0.30, 0.005, "Jet");
~~~

- PYTHIA-only jets:
~~~{.cxx}
AliEmcalJetTask *jetTaskMC= AddTaskEmcalJet("PicoTracks","", 1, 0.2, 1, 0.15, 0.30, 0.005, "JetMConly");
jetTaskMC->SelectConstituents(TObject::kBitMask, 0);
~~~

## Embedding of generated (particle level) PYTHIA

It is possible also to embed PYTHIA events from a generator. The relevant task is AliJetEmbeddingFromGenTask.

- Generate Pythia events on-the-fly
~~~{.cxx}
gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddMCGenPythia.C");
Double_t e_cms = 5020.;
Double_t ptminHard = 30.;
Double_t ptmaxHard = 100.;
AliGenPythia *generPythia = AddMCGenPythia(e_cms,ptminHard,ptmaxHard,2);
~~~

- Emded the tracks into the data event

~~~{.cxx}
Int_t labPythia = 99999.; // a large number, to give the tracks a label that is larger than that of the tracks already present in the event
TString kTracksName = "PicoTracks";
gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetEmbeddingFromGen.C");
AliJetEmbeddingFromGenTask *genTask = AddTaskJetEmbeddingFromGen(generPythia,
                             kTracksName,
                             "EmbFromGenTask",
                             0.15, 1e6,
                             -0.9,0.9,
                             0., TMath::TwoPi(),
                             kTRUE,kTRUE);
genTask->SetSuffix("");  //this will change the name of the track array
genTask->SetChargedOnly(kTRUE);
genTask->SetMarkMC(labPythia);
genTask->SetGeometryName("EMCAL_COMPLETE12SMV1");  //need to set the Geometry name, not so important which one at generated level
~~~

<strong>Important!</strong>
The embedding of on-the-fly generated PYTHIA events needs an empty ESD event. You need to create a list (e.g. `TString localFiles("AliESDs.list")`) that contains the path to the `esdempty.root` file that can be downloaded here: [AliESDs.list](https://twiki.cern.ch/twiki/pub/ALICE/EMCalJetEmbedding/AliESDs.list), [esdempty.root](https://twiki.cern.ch/twiki/pub/ALICE/EMCalJetEmbedding/esdempty.root)

## Embedding of single track

The task providing this functionality is AliJetEmbeddingTask and the corresponding AddTaskJetEmbedding.C.

The p<sub>T</sub> (and mass) of the tracks to be embedded can be drawn randomly
   - From a flat distribution
   - From a function or histogram given as input
   - From a TTree: in this case the tracks are taken in order from the tree (it is suggested to randomize the tree entries first)


For the following examples these variable are used:
~~~{.cxx}
Bool_t copyTracks = kTRUE; // whether the track array is copied to a new one before embedding the new track(s)
TString newName = "Emb"; // remember to set a suffix to be appended to the name of the copy of the track array
Int_t NTrEmb = 1; // number of tracks to be embedded per event
Int_t NClEmb = 0; // number of clusters to be embedded per event
Int_t labEmb = 99999; // a large number that will be the label of the first embedded track/event
~~~

- Example 1: embedding of massless tracks with flat p<sub>T</sub> distribution
~~~{.cxx}
AliJetEmbeddingTask *taskEmb = AddTaskJetEmbedding(tracksName.Data(),"","SingleTrackEmbedding",40.,120.,-0.5,0.5,0., TMath::TwoPi(), NTrEmb, NClEmb, copyTracks);
  taskEmb->SelectCollisionCandidates(pSel);
  taskEmb->SetMarkMC(labEmb);
  taskEmb->SetMasslessParticles(kTRUE);
  taskEmb->SetSuffix(newName);
  TString trackEmb = taskEmb->GetOutTrackName();
~~~

- Example 2: embedding of massive tracks with flat p<sub>T</sub> distribution. As Example 1 substituting `taskEmb->SetMasslessParticles(kTRUE);` with

~~~{.cxx}
  Int_t mass = 8; //GeV
  taskEmbM->SetMasslessParticles(kFALSE);
  taskEmbM->SetMass(mass);
~~~

- Example 3: embedding of tracks from a p<sub>T</sub> and a mass distribution

~~~{.cxx}
TString pathFileMass, pathFilepT, histonameM, histonamepT; //set appropriately
AliJetEmbeddingTask *taskEmbM = AddTaskJetEmbedding(tracksName.Data(),"","SingleTrackEmbedding",0.,0.,-0.5,0.5,0., TMath::TwoPi(), NTrEmb, NClEmb, copyTracks);
  taskEmbDistr->SelectCollisionCandidates(pSel);
  taskEmbDistr->SetMarkMC(labEmb);
  taskEmbDistr->SetMasslessParticles(kFALSE);
  taskEmbDistr->SetSuffix(newName);
  taskEmbDistr->SetMassAndPtDistributionFromFile(pathFileMass, pathFilepT, histonameM, histonamepT); // the file can be stored on alien
  TString trackEmbDistr = taskEmbDistr->GetOutTrackName();
~~~

Several methods avaiable to give the input distributions, check AliJetModelBaseTask and AliJetEmbeddingTask.

- Example 4: embedding with input TTree containing TLorentzVector of the tracks to be embedded. This allows e.g. to keep track of the particle level TLorentzVector of each track in another TTree. The TTree can be prepared with a task like AliAnalysisTaskPrepareInputForEmbedding.

By default, the embedding is performed using detector level tracks, but applying a cut at particle level p<sub>T</sub>.

~~~{.cxx}
TString pathFileTree, treeName, branchName; //set appropriately
AliJetEmbeddingTask *taskEmb4Vect = AddTaskJetEmbedding(tracksName.Data(),"","SingleTrackEmbedding",0.,0.,0,0,0., TMath::TwoPi(), NTrEmb, NClEmb, copyTracks);
  taskEmb4Vect->SelectCollisionCandidates(pSel);
  taskEmb4Vect->SetMarkMC(labEmb);
  taskEmb4Vect->SetNamesForTree(pathFileTree, treeName, branchName);
  taskEmb4Vect->SetSuffix(newName);
  TString trackEmb4Vect = taskEmb4Vect->GetOutTrackName(); // gives the final name of the track array with embedded track
~~~

# Example of analysis

A task to perform the matching between jets in different collections is available: AliAnalysisTaskEmcalJetTagger (see also https://twiki.cern.ch/twiki/bin/view/ALICE/EMCalJetTagger):


~~~{.cxx}
AliAnalysisTaskEmcalJetTagger *taskTagPy =  AddTaskEmcalJetTagger( jetTask->GetName(), jetTaskMC->GetName(), 0.4, "", "", "PicoTracks","","TPC", "V0M", pSel);
taskTagPy->SetNCentBins(1);
taskTagPy->SetIsPythia(kTRUE);
taskTagPy->SetJetTaggingType(AliAnalysisTaskEmcalJetTagger::kClosest);
taskTagPy->SetMinFractionShared(0.5);
~~~

Run a task to match the two jet collections. You can use PWGJE/EMCALJetTasks/AliJetResponseMaker.cxx\h and the add task macro PWGJE / EMCALJetTasks / macros / AddTaskJetResponseMaker.C

~~~{.cxx}
AliJetResponseMaker *taskRM = AddTaskJetResponseMaker("PicoTracks","", jetTaskMC->GetName(),0.2,"PicoTracks","", jetTask->GetName(),"",0.2);
taskRM->GetJetContainer(0)->SetIsParticleLevel(kFALSE);
taskRM->GetJetContainer(1)->SetIsParticleLevel(kFALSE);
taskRM->SetIsPythia(kFALSE);
taskRM->SetIsEmbedded(kTRUE);
~~~



*/
