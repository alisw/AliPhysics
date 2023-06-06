// Analysis task to perform embedding of the simulated AOD
// into real data (AOD)
// Output container contain CaloCluster branches with
//   0. Signal
//   1. Background (Data before embedding)
//   2. Merged (Data+with embedded signal)
//   3. MC information from signal
// CaloClusters are produced as if PHOS Tender was run, no explicit Tender wagon is necessary
// Data for embedding are obtained using standard Manager
// Signal is read from TChain with aodTree
//
// Authors: Dmitri Peressounko
// Date   : 28.05.2011
// updated to support Run2, Daiki Sekihata

#include "AliPHOSEmbeddingRun2.h"

#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBPath.h"
#include "AliESDEvent.h"
#include "AliGRPObject.h"
#include "AliGeomManager.h"
#include "AliInputEventHandler.h"
#include "AliMagF.h"
#include "AliMultSelection.h"
#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliPHOSAodCluster.h"
#include "AliPHOSCalibData.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSCpvRecPoint.h"
#include "AliPHOSDigit.h"
#include "AliPHOSEmcRecPoint.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSReconstructor.h"
#include "AliVEventHandler.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TGeoGlobalMagField.h"
#include "TGeoManager.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TRandom.h"

ClassImp(AliPHOSEmbeddingRun2)

  //________________________________________________________________________
  AliPHOSEmbeddingRun2::AliPHOSEmbeddingRun2(const char* name)
  : AliAnalysisTaskSE(name),
    fAODChain(nullptr),
    fSignal(nullptr),
    fDigitsTree(nullptr),
    fClustersTree(nullptr),
    fTreeOut(nullptr),
    fDigitsArr(nullptr),
    fmcparticles(nullptr),
    fSignalClusters(nullptr),
    fEmbeddedClusters(nullptr),
    fEmbeddedCells(nullptr),
    fCellsPHOS(nullptr),
    fEmcRecPoints(nullptr),
    fCpvRecPoints(nullptr),
    fPHOSGeo(nullptr),
    fClusterizer(nullptr),
    fPHOSReconstructor(nullptr),
    fCalibMC(nullptr),
    fCalibData(nullptr),
    fPathPrivateOADBMC("$ALICE_PHYSICS/OADB/PHOS/PHOSMCCalibrations.root"),
    fSignalECorrection(1.),
    fNSignal(0),
    fNCaloClustersOld(0),
    fRunNumber(0),
    fMF(0),
    fInitialized(0),
    fVtx(0., 0., 0.)
{
  // Constructor
  for (int i = 0; i < 5; i++) {
    fOldPHOSCalibration[i] = nullptr;
    fRunByRunCorr[i] = 0.136;
  }
}
//________________________________________________________________________
void AliPHOSEmbeddingRun2::UserCreateOutputObjects()
{
  // prepare output containers

  fDigitsArr = new TClonesArray("AliPHOSDigit", 4 * 56 * 64);

  AliAODHandler* handler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
  if (handler) {
    fTreeOut = handler->GetTree();
  } else {
    AliWarning("No output AOD Event Handler connected.");
  }
  PostData(0, fTreeOut);
}
//______________________________________________________________________________
void AliPHOSEmbeddingRun2::Init()
{
  if (fInitialized)
    return;

  AliAODEvent* event = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!event) {
    AliError("Can not obtain Bg AOD InputEvent!");
    return;
  }
  fRunNumber = event->GetRunNumber();

  int nEvents =
    ((AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler()))->GetTree()->GetEntriesFast();
  nEvents *= fSelEventsPart; // Only 20% of events will be useful
  RunSignalSimulation(nEvents);
  // Connect simulated signal
  TChain* chainAOD = new TChain("aodTree");
  chainAOD->AddFile("AliAOD.root");
  SetSignalChain(chainAOD);

  AliCDBManager::Instance()->SetDefaultStorage("raw://");
  AliCDBManager::Instance()->SetRun(fRunNumber);
  AliCDBPath path("PHOS", "Calib", "RecoParam");
  AliCDBEntry* entry = AliCDBManager::Instance()->Get(path.GetPath());
  if (!entry) {
    AliError(Form("Can not get OCDB entry %s", path.GetPath().Data()));
    return;
  }
  TObjArray* recoParamArray = (TObjArray*)entry->GetObject();
  AliPHOSRecoParam* recoParam = (AliPHOSRecoParam*)recoParamArray->At(2);

  fPHOSGeo = AliPHOSGeometry::GetInstance("Run2");

  fPHOSReconstructor = new AliPHOSReconstructor("Run2"); // Name of geometry!
  fPHOSReconstructor->SetRecoParam(recoParam);

  InitMF();

  // Geometry
  AliOADBContainer geomContainer("phosGeo");
  geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSMCGeometry.root", "PHOSMCRotationMatrixes");
  TObjArray* matrixes = (TObjArray*)geomContainer.GetObject(fRunNumber, "PHOSRotationMatrixes");
  for (int mod = 0; mod < 6; mod++) {
    if (!matrixes->At(mod)) {
      continue;
    }
    fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)), mod);
    AliInfo(Form(".........Adding Matrix(%d), geo=%p\n", mod, fPHOSGeo));
    ((TGeoHMatrix*)matrixes->At(mod))->Print();
  }
  // Now geometry for Clusterizer
  //  Import ideal TGeo geometry and apply misalignment
  if (!gGeoManager) {
    AliGeomManager::LoadGeometry("geometry.root");
    if (!gGeoManager) {
      AliFatal("Can not load geometry");
    }
    if (!AliGeomManager::CheckSymNamesLUT("PHOS")) {
      AliFatal("CheckSymNamesLUT");
    }
  }
  TString detStr = "PHOS";
  TString loadAlObjsListOfDets = "PHOS";
  if (AliGeomManager::GetNalignable("GRP") != 0)
    loadAlObjsListOfDets.Prepend("GRP "); // add alignment objects for non-sensitive modules
  //  AliGeomManager::ApplyAlignObjsFromCDB(loadAlObjsListOfDets.Data());

  AliCDBManager::Instance()->UnloadFromCache("*/Align/*");
  AliCDBManager::Instance()->UnloadFromCache("GRP/Geometry/Data");

  // Get AODB MC calibration for this run
  AliOADBContainer calibContainerMC("phosRecalibration");
  calibContainerMC.InitFromFile(fPathPrivateOADBMC, "phosRecalibration");
  TObjArray* recalibMC = (TObjArray*)calibContainerMC.GetObject(fRunNumber, "PHOSRecalibration", "Run2NoNcellCut");
  if (!recalibMC) {
    AliFatal(Form("Can not read calibrations for run %d, do not apply OADB de-calibration\n", fRunNumber));
  }
  AliInfo(Form("Reading recalibration for MC from pass %s", fPathPrivateOADBMC.Data()));
  const int recoPass = 1;
  fCalibMC = (AliPHOSCalibData*)recalibMC->At(recoPass - 1);

  // Calibration of data (different from MC)
  // Check the pass1-pass2-pass3 reconstruction
  AliOADBContainer calibContainerRaw("phosRecalibration");
  calibContainerRaw.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSCalibrations.root", "phosRecalibration");
  TObjArray* recalibRaw = (TObjArray*)calibContainerRaw.GetObject(fRunNumber, "PHOSRecalibration");
  if (!recalibRaw) {
    AliFatal(Form(
      "Can not read calibrations for run %d\n. You may choose your specific calibration with ForceUsingCalibration()\n",
      fRunNumber));
  } else {
    fCalibData = (AliPHOSCalibData*)recalibRaw->At(recoPass - 1);
    if (!fCalibData) {
      AliFatal(Form("Can not find calibration for run %d, pass %d \n", fRunNumber, recoPass));
    }
  }

  // Calibration used for Clusterization
  fgCalibData = new AliPHOSCalibData(-1); // use AliCDBManager's run number

  // L1phase and run-by-run correction for Run2
  // L1phase for Run2
  AliOADBContainer L1Container("phosL1Calibration");
  L1Container.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSL1Calibrations.root", "phosL1Calibration");
  TNamed* a = (TNamed*)L1Container.GetObject(fRunNumber);
  if (!a) {
    AliError(Form("L1phase for run %d was not found, time calibration will be wrong!\n", fRunNumber));
    for (int ii = 0; ii < 15; ii++)
      fL1phase[ii] = 0;
  } else {
    const char* c = a->GetName();
    for (int ii = 0; ii < 15; ii++)
      fL1phase[ii] = c[ii] - '0';
  }

  // Run-by-run correction
  AliOADBContainer runByRunContainer("phosRunByRunCalibration");
  runByRunContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSRunByRunCalibrations.root", "phosRunByRunCalibration");
  TNamed* rbr = (TNamed*)runByRunContainer.GetObject(fRunNumber);
  if (rbr) {
    sscanf(rbr->GetName(), "%f,%f,%f,%f", &fRunByRunCorr[1], &fRunByRunCorr[2], &fRunByRunCorr[3], &fRunByRunCorr[4]);
  }
  // In any case correction should not be zero
  // If it is zero, set default and write warning
  for (int mod = 1; mod < 5; mod++) {
    if (fRunByRunCorr[mod] == 0.) {
      fRunByRunCorr[mod] = 0.136;
      AliWarning(Form("Run-by-Run correction for mod. %d is zero in run %d", mod, fRunNumber));
    }
  }

  AliOADBContainer badmapContainer(Form("phosBadMap"));
  // if(fPrivateOADBBadMap.Length()!=0){
  //    //Load standard bad maps file if no OADB file is force loaded
  //    AliInfo(Form("using custom bad channel map from %s\n",fPrivateOADBBadMap.Data()));
  //     badmapContainer.InitFromFile(fPrivateOADBBadMap.Data(),"phosBadMap");
  //  } else {
  //    //Load force loaded OADB file
  //    AliInfo("using standard bad channel map from $ALICE_PHYSICS/OADB/PHOS/PHOSBadMaps.root\n");
  badmapContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSBadMaps.root", "phosBadMap");
  // }
  TObjArray* maps = (TObjArray*)badmapContainer.GetObject(fRunNumber, "phosBadMap");
  if (!maps) {
    AliError(
      Form("Can not read Bad map for run %d. \n You may choose to use your map with ForceUsingBadMap()\n", fRunNumber));
  } else {
    AliInfo(Form("Setting PHOS bad map with name %s \n", maps->GetName()));
    for (int mod = 1; mod < 5; mod++) {
      TH2I* h = (TH2I*)maps->At(mod);
      fPHOSBadMap[mod] = new TH2I(*h);
    }
  }

  fInitialized = kTRUE;
}

//________________________________________________________________________
void AliPHOSEmbeddingRun2::UserExec(Option_t*)
{
  // Main loop, called for each event
  // Perform embedding
  // To avoid double counting of evens if new was not embedded
  Init();

  if (fmcparticles)
    fmcparticles->Clear();
  if (fSignalClusters)
    fSignalClusters->Clear();
  if (fEmbeddedClusters)
    fEmbeddedClusters->Clear();

  AliAODEvent* event = dynamic_cast<AliAODEvent*>(InputEvent()); // for general information before embedding
  if (!event) {
    AliError("ERROR: Could not retrieve Bg AOD event");
    PostData(0, fTreeOut);
    return;
  }

  TString trigClasses = event->GetFiredTriggerClasses();
  AliInfo(Form("trig.class %s, period %d, bc %d, orbit %d", trigClasses.Data(), event->GetPeriodNumber(),
               event->GetBunchCrossNumber(), event->GetOrbitNumber()));

  float acentrality = -1;
  AliMultSelection* fMultSelection = (AliMultSelection*)event->FindListObject("MultSelection");
  if (!fMultSelection) {
    AliError("AliMultSelection object not found!");
    return;
  } else {
    acentrality = fMultSelection->GetMultiplicityPercentile("V0M");
  }

  if (acentrality <= 0. || acentrality > 101.) {
    AliError(Form("Centrality = %f", acentrality));
    PostData(0, fTreeOut);
    return;
  }

  const AliAODVertex* vtx = event->GetPrimaryVertex();
  fVtx.SetXYZ(vtx->GetX(), vtx->GetY(), vtx->GetZ());
  if (TMath::Abs(vtx->GetZ()) > 10.) {
    PostData(0, fTreeOut);
    return;
  }

  // Create output branch if necessary
  //  Note that AODevent will keepo pointer to TClonesArray, but Embedding should own it
  //   MC particles
  TObject* br0 = event->FindListObject(AliAODMCParticle::StdBranchName());
  if (!br0) {
    fmcparticles = new TClonesArray("AliAODMCParticle", 0);
    fmcparticles->SetName(AliAODMCParticle::StdBranchName());
    event->AddObject(fmcparticles);
  }

  TObject* br1 = event->FindListObject("SignalCaloClusters");
  if (!br1) {
    if (!fSignalClusters) {
      fSignalClusters = new TClonesArray("AliAODCaloCluster", 0);
    }
    fSignalClusters->SetName("SignalCaloClusters");
    event->AddObject(fSignalClusters);
  }

  TObject* br2 = event->FindListObject("EmbeddedCaloClusters");
  if (!br2) {
    if (!fEmbeddedClusters) {
      fEmbeddedClusters = new TClonesArray("AliAODCaloCluster", 0);
    }
    fEmbeddedClusters->SetName("EmbeddedCaloClusters");
    event->AddObject(fEmbeddedClusters);
  }

  //Check if run number changed. If changed, do nothing, PHOS reconstruction will fail
  if(fRunNumber != event->GetRunNumber()){
    PostData(0, fTreeOut);
    return;
  }


  // Read next AOD event
  // If necesary method checks if AOD event is good to embed
  // e.g. there are PHOS clusters etc.
  if (!GetNextSignalEvent()) {
    AliError("ERROR: Could not retrieve signal event");
    PostData(0, fTreeOut);
    return;
  }

  CopyRecalibrateSignal();

  CopyRecalibrateBackground();

  // Now perform real embedding
  // Embed signal clusters into event
  MakeEmbedding();

  PostData(0, fTreeOut);
}
//______________________________________________________________________________
void AliPHOSEmbeddingRun2::CopyRecalibrateSignal()
{
  // Apply recalibration a-la Tender and copy Signal event to output
  // Copy MC particles
  ConvertMCParticles();

  // For re-calibration
  const double logWeight = 4.5;

  TClonesArray* clusters = fSignal->GetCaloClusters();
  AliAODCaloCells* cells = fSignal->GetPHOSCells();

  int multClust = clusters->GetEntriesFast();
  // Add noise
  if (fAddNoiseMC) {
    short ncell = cells->GetNumberOfCells();
    short cellNumber;
    double amplitude = 0., time = 0., efrac = 0.;
    int mclabel;
    for (short pos = 0; pos < ncell; pos++) {
      cells->GetCell(pos, cellNumber, amplitude, time, mclabel, efrac);
      amplitude = TMath::Max(0., amplitude + gRandom->Gaus(0, fNoiseMC));
      // Apply ZeroSuppression if necessary
      if (amplitude < fZScut) {
        amplitude = 1.e-6;
      }
      bool isHG = cells->GetHighGain(pos);
      cells->SetCell(pos, cellNumber, amplitude, time, mclabel, efrac, isHG);
    }
  }

  int iCluNew = 0;
  for (int i = 0; i < multClust; i++) {
    AliAODCaloCluster* cluOld = static_cast<AliAODCaloCluster*>(clusters->At(i));
    if (cluOld->GetType() != AliVCluster::kPHOSNeutral)
      continue;

    float positionOld[3];
    cluOld->GetPosition(positionOld);
    TVector3 globalOld(positionOld);

    // Apply re-Calibreation
    AliPHOSAodCluster cluPHOS(*cluOld);
    cluPHOS.Recalibrate(fCalibMC, cells); // modify the cell energies
    cluPHOS.EvalAll(logWeight, fVtx);     // recalculate the cluster parameters

    float position[3];
    cluPHOS.GetPosition(position);
    TVector3 global(position);
    int relId[4];
    fPHOSGeo->GlobalPos2RelId(global, relId);
    int mod = relId[0];
    int cellX = relId[2];
    int cellZ = relId[3];
    if (!IsGoodChannel(mod, cellX, cellZ)) {
      continue;
    }

    float energy = CorrectNonlinearityMC(cluPHOS.E());
    AliAODCaloCluster* clu = new ((*fSignalClusters)[iCluNew]) AliAODCaloCluster(
      iCluNew, cluOld->GetNLabels(), cluOld->GetLabels(), energy, position, nullptr, AliVCluster::kPHOSNeutral, 0);
    iCluNew++;

    double ecore = CoreEnergy(&cluPHOS);
    ecore = CorrectNonlinearity(ecore);
    clu->SetCoreEnergy(ecore); // core particle energy

    // Eval FullDispersion
    clu->SetDispersion(TestFullLambda(clu->E(), cluPHOS.GetM20(), cluPHOS.GetM02()));
    // Eval CoreDispersion
    double m02 = 0., m20 = 0.;
    EvalLambdas(&cluPHOS, m02, m20);
    clu->SetChi2(TestCoreLambda(clu->E(), m20, m02)); // not yet implemented
    clu->SetM02(m02);                                 // second moment M2x
    clu->SetM20(m20);                                 // second moment M2z
    clu->SetNCells(cluOld->GetNCells());

    // double minDist=clu->GetDistanceToBadChannel() ;//Already calculated
    //    DistanceToBadChannel(mod,&locPos,minDist);
    //    clu->SetDistanceToBadChannel(minDist) ;

    // double ecross = EvalEcross(&cluPHOS);
    // clu->SetMCEnergyFraction(ecross) ;
  }
}
//______________________________________________________________________________
void AliPHOSEmbeddingRun2::CopyRecalibrateBackground()
{
  // Recalibrate background a-la Tender
  const double logWeight = 4.5;

  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(InputEvent()); // for general information before embedding
  TClonesArray* clusters = aod->GetCaloClusters();
  AliAODCaloCells* cells = aod->GetPHOSCells();
  int multClust = clusters->GetEntriesFast();
  for (int i = 0; i < multClust; i++) {
    AliAODCaloCluster* clu = static_cast<AliAODCaloCluster*>(clusters->At(i));
    if (clu->GetType() != AliVCluster::kPHOSNeutral)
      continue;

    float positionOld[3];
    clu->GetPosition(positionOld);
    TVector3 globalOld(positionOld);

    // Apply re-Calibreation
    AliPHOSAodCluster cluPHOS(*clu);
    cluPHOS.Recalibrate(fCalibData, cells); // modify the cell energies
    cluPHOS.EvalAll(logWeight, fVtx);       // recalculate the cluster parameters

    float position[3];
    cluPHOS.GetPosition(position);
    clu->SetPosition(position); // rec.point position in MARS
    TVector3 global(position);
    int relId[4];
    fPHOSGeo->GlobalPos2RelId(global, relId);
    int mod = relId[0];
    int cellX = relId[2];
    int cellZ = relId[3];
    if (!IsGoodChannel(mod, cellX, cellZ)) {
      clu->SetE(0.);
      continue;
    }
    cluPHOS.SetE(0.136 / fRunByRunCorr[mod] * CorrectNonlinearity(cluPHOS.E())); // Users's nonlinearity
    TVector3 locPosOld; // Use it to re-calculate distance to track
    fPHOSGeo->Global2Local(locPosOld, globalOld, mod);

    double ecore = CoreEnergy(&cluPHOS);
    ecore = 0.136 / fRunByRunCorr[mod] * CorrectNonlinearity(ecore);

    clu->SetE(cluPHOS.E());    // total particle energy
    clu->SetCoreEnergy(ecore); // core particle energy

    // Eval FullDispersion
    clu->SetDispersion(TestFullLambda(clu->E(), cluPHOS.GetM20(), cluPHOS.GetM02()));
    // Eval CoreDispersion
    double m02 = 0., m20 = 0.;
    EvalLambdas(&cluPHOS, m02, m20);
    clu->SetChi2(TestCoreLambda(clu->E(), m20, m02)); // not yet implemented
    clu->SetM02(m02);                                 // second moment M2x
    clu->SetM20(m20);                                 // second moment M2z
    clu->SetNCells(cluPHOS.GetNCells());

    // correct distance to track
    double pttrack = 0.;
    int charge = 0;
    double dx = 999., dz = 999.;
    TVector3 locPos;
    fPHOSGeo->Global2Local(locPos, global, mod);
    int itr = FindTrackMatching(mod, &locPos, dx, dz, pttrack, charge);
    clu->SetTrackDistance(dx, dz);
    double r = 999.; // Big distance
    //      int nTracksMatched = clu->GetNTracksMatched();
    if (itr > 0) {
      r = TestCPV(dx, dz, pttrack, charge);
    }
    clu->SetEmcCpvDistance(r); // Distance in sigmas
    if (itr >= 0) {            // there is a track
      // Remove existing
      AliAODTrack* tr = (AliAODTrack*)aod->GetTrack(itr);
      int ntrM = clu->GetNTracksMatched();
      if (ntrM > 0) {
        AliAODTrack* trStored = (AliAODTrack*)clu->GetTrackMatched(0);
        if (trStored != tr) {
          clu->RemoveTrackMatched(trStored);
          clu->AddTrackMatched(tr);
        }
      } else {
        clu->AddTrackMatched(tr);
      }
    }

    double tof = EvalTOF(&cluPHOS, cells);
    clu->SetTOF(tof);
    // double minDist=clu->GetDistanceToBadChannel() ;//Already calculated
    // DistanceToBadChannel(mod,&locPos,minDist);
    // clu->SetDistanceToBadChannel(minDist) ;

    // double ecross = EvalEcross(&cluPHOS);
    // clu->SetMCEnergyFraction(ecross) ;
  }
}
//________________________________________________________________________
void AliPHOSEmbeddingRun2::MakeEmbedding()
{
  // Perform embedding of the signal to the event

  const AliAODEvent* bgevent = dynamic_cast<AliAODEvent*>(InputEvent());
  // Prepare format for Clusterizer/Reconstruturer
  gROOT->cd(); // make sure that the digits and RecPoints Trees are memory resident
  if (fDigitsTree) {
    delete fDigitsTree;
  }
  fDigitsTree = new TTree("digitstree", "digitstree");
  fDigitsTree->Branch("PHOS", "TClonesArray", &fDigitsArr, 32000);

  if (fClustersTree) {
    delete fClustersTree;
  }
  fClustersTree = new TTree("clustertree", "clustertree");

  // Remember number of Clusters before we added new ones...
  fNCaloClustersOld = bgevent->GetNumberOfCaloClusters();

  int nPHOSBefore = 0;
  int nCPVBefore = 0;
  int nEMCALBefore = 0;

  for (int iClust = 0; iClust < bgevent->GetNumberOfCaloClusters(); ++iClust) {
    AliAODCaloCluster* cluster = bgevent->GetCaloCluster(iClust);

    if (cluster->IsPHOS()) {
      if (cluster->GetType() == AliVCluster::kPHOSNeutral)
        nPHOSBefore++;
      if (cluster->GetType() == AliVCluster::kPHOSCharged)
        nCPVBefore++;
    } else {
      nEMCALBefore++;
    }
  }
  AliInfo(Form("Before embedding: Nall = %d , nPHOS = %d , nCPV = %d , nEMCAL = %d.", fNCaloClustersOld, nPHOSBefore,
               nCPVBefore, nEMCALBefore));

  // create digits
  MakeDigits(bgevent, fSignal);
  // clusterize and make tracking
  fPHOSReconstructor->Reconstruct(fDigitsTree, fClustersTree);
  ConvertEmbeddedClusters();
}
//________________________________________________________________________
bool AliPHOSEmbeddingRun2::GetNextSignalEvent()
{
  // Read signal AOD event from the chain

  if (fAODChain == nullptr) {
    AliError(Form("No chain to read signal events: "));
    return false;
  }

  if (fSignal) {
    delete fSignal;
  }
  fSignal = new AliAODEvent;
  fSignal->ReadFromTree(fAODChain);
  return fAODChain->GetEvent(fNSignal++) > 0; // event read without errors
}
//______________________________________________________________________________
void AliPHOSEmbeddingRun2::ConvertEmbeddedClusters()
{
  // Copy PHOS clusters and cells after embedding
  TBranch* emcbranch = fClustersTree->GetBranch("PHOSEmcRP");
  if (!emcbranch) {
    AliError("can't get the branch with the PHOS EMC clusters !");
    return;
  }
  emcbranch->SetAddress(&fEmcRecPoints);
  emcbranch->GetEntry(0);

  TBranch* cpvbranch = fClustersTree->GetBranch("PHOSCpvRP");
  if (cpvbranch) {
    cpvbranch->SetAddress(&fCpvRecPoints);
    cpvbranch->GetEntry(0);
  }

  // Access to the AOD container of clusters
  fEmbeddedClusters->Clear();
  fEmbeddedClusters->Expand(fEmcRecPoints->GetEntriesFast() + fCpvRecPoints->GetEntriesFast());

  int jClusters(0);
  for (int iClust = 0; iClust < fEmcRecPoints->GetEntriesFast(); ++iClust) {
    AliPHOSEmcRecPoint* emcRP = static_cast<AliPHOSEmcRecPoint*>(fEmcRecPoints->At(iClust));

    // Primaries
    int nLabel = 0;
    int* labels = emcRP->GetPrimaries(nLabel);

    float energy = emcRP->GetEnergy();
    energy = 0.136 / fRunByRunCorr[emcRP->GetPHOSMod()] * CorrectNonlinearity(energy);
    float ecore = emcRP->GetCoreEnergy();
    ecore = 0.136 / fRunByRunCorr[emcRP->GetPHOSMod()] * CorrectNonlinearity(ecore);
    if (energy < 0.001) { // skip zero energy clusters
      continue;
    }

    TVector3 local;
    emcRP->GetLocalPosition(local);
    // Correct for the non-perpendicular incidence
    //  Correction for the depth of the shower starting point (TDR p 127)
    const float para = 0.925;
    const float parb = 6.52;

    // Remove Old correction (vertex at 0,0,0)
    TVector3 vtxOld(0., 0., 0.);
    TVector3 vInc;
    float x = local.X();
    float z = local.Z();
    fPHOSGeo->GetIncidentVector(vtxOld, emcRP->GetPHOSMod(), x, z, vInc);
    float depthxOld = 0.;
    float depthzOld = 0.;
    if (energy > 0 && vInc.Y() != 0.) {
      depthxOld = (para * TMath::Log(energy) + parb) * vInc.X() / TMath::Abs(vInc.Y());
      depthzOld = (para * TMath::Log(energy) + parb) * vInc.Z() / TMath::Abs(vInc.Y());
    } else {
      AliError("Cluster with zero energy \n");
    }
    // Apply Real vertex
    fPHOSGeo->GetIncidentVector(fVtx, emcRP->GetPHOSMod(), x, z, vInc);
    float depthx = 0.;
    float depthz = 0.;
    if (energy > 0 && vInc.Y() != 0.) {
      depthx = (para * TMath::Log(energy) + parb) * vInc.X() / TMath::Abs(vInc.Y());
      depthz = (para * TMath::Log(energy) + parb) * vInc.Z() / TMath::Abs(vInc.Y());
    }

    // Correct for the vertex position and shower depth
    double xd = x + (depthxOld - depthx);
    double zd = z + (depthzOld - depthz);
    TVector3 dir(0, 0, 0);
    fPHOSGeo->Local2Global(emcRP->GetPHOSMod(), xd, zd, dir);
    //  account this in GetMomentum function
    //    dir-=fVtx ;
    //    dir.SetMag(1.) ; //Momentum direction
    float posF[3] = { float(dir.x()), float(dir.y()), float(dir.z()) };

    AliAODCaloCluster* caloCluster = new ((*fEmbeddedClusters)[jClusters])
      AliAODCaloCluster(jClusters, nLabel, labels, energy, posF, nullptr, AliVCluster::kPHOSNeutral, 0);
    jClusters++;
    caloCluster->SetCoreEnergy(ecore);
    // fills the ESDCaloCluster
    caloCluster->SetDispersion(emcRP->GetDispersion()); // cluster dispersion
    caloCluster->SetM02(emcRP->GetM2x());               // second moment M2x
    caloCluster->SetM20(emcRP->GetM2z());               // second moment M2z
    caloCluster->SetNExMax(emcRP->GetNExMax());         // number of local maxima
    caloCluster->SetChi2(-1);                           // not yet implemented
    caloCluster->SetNCells(emcRP->GetMultiplicity());

    // recalibrate time?
    caloCluster->SetTOF(emcRP->GetTime()); // Time of flight - already calibrated in EMCRecPoint

    // Cells contributing to clusters
    int cellMult = emcRP->GetDigitsMultiplicity();
    int* digitsList = emcRP->GetDigitsList();
    float* rpElist = emcRP->GetEnergiesList();
    unsigned short absIdList[cellMult];
    double fracList[cellMult];
    for (int iCell = 0; iCell < cellMult; iCell++) {
      AliPHOSDigit* digit = static_cast<AliPHOSDigit*>(fDigitsArr->At(digitsList[iCell]));
      absIdList[iCell] = (unsigned short)(digit->GetId());
      fracList[iCell] = rpElist[iCell];
    }
    caloCluster->SetNCells(cellMult);
    caloCluster->SetCellsAbsId(absIdList);
    caloCluster->SetCellsAmplitudeFraction(fracList);

    // a-la Tender distance to track
    double dx = 999., dz = 999.;
    double pttrack = 0.;
    int charge = 0;
    int itr = FindTrackMatching(emcRP->GetPHOSMod(), &local, dx, dz, pttrack, charge);
    caloCluster->SetTrackDistance(dx, dz);
    double r = TestCPV(dx, dz, pttrack, charge);
    caloCluster->SetEmcCpvDistance(r);
    if (itr >= 0) { // there is a track
      caloCluster->AddTrackMatched((AliAODTrack*)dynamic_cast<AliAODEvent*>(InputEvent())->GetTrack(itr));
    }
  }

  // Fill CPV clusters
  int nOfCPVclu = fCpvRecPoints->GetEntriesFast();
  for (int recpoint = 0; recpoint < nOfCPVclu; recpoint++) {
    AliPHOSCpvRecPoint* cpvRP = static_cast<AliPHOSCpvRecPoint*>(fCpvRecPoints->At(recpoint));

    TVector3 pos;
    fPHOSGeo->GetGlobalPHOS(cpvRP, pos);
    double xyz[3] = { pos.x(), pos.y(), pos.z() };

    // Create cell lists
    int cellMult = cpvRP->GetDigitsMultiplicity();
    // int    *digitsList = cpvRP->GetDigitsList();
    // float  *rpElist    = cpvRP->GetEnergiesList() ;
    // unsigned short *absIdList  = new unsigned short[cellMult];
    // double *fracList   = new double[cellMult];

    // Primaries
    int primMult = 0;
    int* primList = cpvRP->GetPrimaries(primMult);
    float energy = cpvRP->GetEnergy();

    AliAODCaloCluster* caloCluster = new ((*fEmbeddedClusters)[jClusters++])
      AliAODCaloCluster(jClusters, primMult, primList, energy, xyz, nullptr, AliVCluster::kPHOSCharged, 0);

    caloCluster->SetNCells(cellMult);
    caloCluster->SetDispersion(cpvRP->GetDispersion()); // cluster dispersion
    caloCluster->SetM02(cpvRP->GetM2x());               // second moment M2x
    caloCluster->SetM20(cpvRP->GetM2z());               // second moment M2z
    caloCluster->SetNExMax(cpvRP->GetNExMax());         // number of local maxima
    caloCluster->SetChi2(-1);                           // not yet implemented
  }
}
//______________________________________________________________________________
void AliPHOSEmbeddingRun2::ConvertMCParticles()
{
  // Copy MC branches to new AOD

  TClonesArray* mcArray = (TClonesArray*)fSignal->FindListObject(AliAODMCParticle::StdBranchName());

  fmcparticles->Clear();
  for (int i = 0; i < mcArray->GetEntriesFast(); i++) {
    AliAODMCParticle* aodpart = (AliAODMCParticle*)mcArray->At(i);
    new ((*fmcparticles)[i]) AliAODMCParticle(*aodpart);
  }
}
//__________________________________________________________________________________
void AliPHOSEmbeddingRun2::MakeDigits(const AliAODEvent* bg, const AliAODEvent* signal)
{
  //-------------------------------------------------------------------------------------
  // Transform CaloCells into Digits which can be used for standard reconstruction
  // Add signal digits to the event
  // Clusterizser uses calibration parameters fgCalibData
  fCellsPHOS = bg->GetPHOSCells();

  fDigitsArr->Clear();

  // First copy data digits
  int ndigit = 0;
  for (short icell = 0; icell < fCellsPHOS->GetNumberOfCells(); icell++) {
    short id = 0;
    double time = 0., amp = 0.;
    int mclabel;
    double efrac = 0.;
    if (fCellsPHOS->GetCell(icell, id, amp, time, mclabel, efrac) != kTRUE)
      break;
    bool isHG = fCellsPHOS->GetHighGain(icell);
    int idLong = id;
    if (id < 0)
      idLong = -id + 56 * 64 * 5; // CPV digits
    int relId[4];
    fPHOSGeo->AbsToRelNumbering(idLong, relId);
    int module = relId[0];
    int row = relId[2];
    int column = relId[3];
    // Change Amp back from Energy to ADC counts
    float calibration = fgCalibData->GetADCchannelEmc(module, column, row);
    // Apply fine tuning. Final energy should be multiplied by coeff.
    double c2 = fCalibData->GetADCchannelEmc(module, column, row);
    if (c2 > 0) {
      calibration /= c2;
    } else {
      calibration = 0;
    }
    amp /= calibration;
    // Un-calibrate time
    double dt = 0;
    if (isHG) {
      dt = fgCalibData->GetTimeShiftEmc(module, column, row);
      dt += fCalibData->GetTimeShiftEmc(module, column, row);
    } else {
      dt = fgCalibData->GetLGTimeShiftEmc(module, column, row);
      dt += fCalibData->GetLGTimeShiftEmc(module, column, row);
    }
    time -= dt;

    // Calibration in Clusterizser with fgCalibData :
    // else{ //EMC
    //   if(isLG)
    //     time += fgCalibData->GetLGTimeShiftEmc(module,column,row);
    //   else
    //     time += fgCalibData->GetTimeShiftEmc(module,column,row);
    // In Tender with fCalibData
    //     tof-=fPHOSCalibData->GetTimeShiftEmc(module, column, row);
    //   else{
    //     tof-=fPHOSCalibData->GetLGTimeShiftEmc(module, column, row);
    //   }

    AliPHOSDigit* d = new ((*fDigitsArr)[ndigit]) AliPHOSDigit(-1, idLong, float(amp), float(time), ndigit);
    ndigit++;
    d->SetLG(!isHG);
  }

  // Add Digits from Signal
  TClonesArray sdigits("AliPHOSDigit", 0);
  int isdigit = 0;
  if (signal) {
    AliAODCaloCells* cellsS = signal->GetPHOSCells();
    int cellLabels[1000] = { 0 }; // 1000 should be enough for simulated
    // int cellSecondLabels[1000]={0} ; //low-statistics event.
    for (int i = 0; i < cellsS->GetNumberOfCells(); i++) {
      cellLabels[i] = -1;
      // cellSecondLabels[i]=-1;
    }
    //------------------------------------------------------------------------------------
    // Ancestry information
    // Celect digits contributing to signal clusters and add primary information
    //(it is not stored in CaloCells)
    //------------------------------------------------------------------------------------
    sdigits.Expand(cellsS->GetNumberOfCells());
    for (int i = 0; i < signal->GetNumberOfCaloClusters(); i++) {
      // cluster from embedded signal
      AliVCluster* clus = signal->GetCaloCluster(i);

      if (!clus->IsPHOS())
        continue;

      int label = clus->GetLabel();
      // int label2 = -1 ;
      // if (clus->GetNLabels()>=2) label2 = clus->GetLabelAt(1) ;

      unsigned short* index = clus->GetCellsAbsId();
      for (int ic = 0; ic < clus->GetNCells(); ic++) {
        for (int icell = 0; icell < cellsS->GetNumberOfCells(); icell++) {
          short cellNumber;
          double cellAmplitude = 0., cellTime = 0.;
          int mclabel;
          double efrac = 0.;
          cellsS->GetCell(icell, cellNumber, cellAmplitude, cellTime, mclabel, efrac);
          int longCellNumber = cellNumber;
          if (cellNumber < 0)
            longCellNumber = -cellNumber + 56 * 64 * 5; // CPV digits
          if (longCellNumber == index[ic]) {
            cellLabels[icell] = label;
            // cellSecondLabels[icell]=label2;
            break;
          }
        }
      }
    }

    for (int icell = 0; icell < cellsS->GetNumberOfCells(); icell++) {
      short cellNumber;
      double cellAmplitude = 0., cellTime = 0.;
      int mclabel;
      double efrac = 0.;
      if (cellsS->GetCell(icell, cellNumber, cellAmplitude, cellTime, mclabel, efrac) != kTRUE)
        break;

      bool isHG = fCellsPHOS->GetHighGain(icell);
      int longCellNumber = cellNumber;
      if (cellNumber < 0)
        longCellNumber = -cellNumber + 56 * 64 * 5; // CPV digits
      if (cellNumber > 0) {
        int relId[4];
        fPHOSGeo->AbsToRelNumbering(longCellNumber, relId);
        int module = relId[0];
        int row = relId[2];
        int column = relId[3];
        // Change Amp back from Energy to ADC counts
        float calibration = fgCalibData->GetADCchannelEmc(module, column, row);
        // Apply fine tuning. Final energy should be multiplied by coeff.
        double c2 = fCalibMC->GetADCchannelEmc(module, column, row);
        if (c2 > 0) {
          calibration /= c2;
        } else {
          calibration = 0;
        }
        cellAmplitude /= calibration;
        // Un-calibrate time
        double dt = 0;
        if (isHG) {
          dt = fgCalibData->GetTimeShiftEmc(module, column, row);
          dt += fCalibMC->GetTimeShiftEmc(module, column, row);
        } else {
          dt = fgCalibData->GetLGTimeShiftEmc(module, column, row);
          dt += fCalibMC->GetLGTimeShiftEmc(module, column, row);
        }
        cellTime -= dt;
      }

      AliPHOSDigit* d = new (sdigits[isdigit])
        AliPHOSDigit(cellLabels[icell], longCellNumber, float(cellAmplitude), float(cellTime), isdigit);
      isdigit++;
      d->SetLG(!isHG);
    }
  }

  // Merge digits
  int icurrent = 0; // index of the last used digit in underlying event
  fDigitsArr->Expand(ndigit + isdigit);
  for (int i = 0; i < isdigit; i++) {
    AliPHOSDigit* sdigit = static_cast<AliPHOSDigit*>(sdigits.At(i));
    bool added = false;
    for (int id = icurrent; id < ndigit; id++) {
      AliPHOSDigit* digit = static_cast<AliPHOSDigit*>(fDigitsArr->At(id));
      if (sdigit->GetId() == digit->GetId()) {
        *digit += *sdigit; // add energies
        icurrent = id + 1;
        added = kTRUE;
        break; // no more digits with same ID in the list
      }
      if (sdigit->GetId() < digit->GetId()) {
        icurrent = id;
        break; // no more digits with same ID in the list
      }
    }
    if (!added) {
      new ((*fDigitsArr)[ndigit]) AliPHOSDigit(*sdigit);
      ndigit++;
    }
  }

  fDigitsArr->Sort();
  for (int i = 0; i < ndigit; i++) {
    AliPHOSDigit* digit = static_cast<AliPHOSDigit*>(fDigitsArr->At(i));
    digit->SetIndexInList(i);
  }
  fDigitsTree->Fill();
}
//____________________________________________________________________________
void AliPHOSEmbeddingRun2::InitMF()
{
  AliInfo("............Init MF \n");
  //------------------------------------
  // Initialization of the Mag.Fiels from GRP entry
  // Copied from AliReconstruction::InitGRP()
  //------------------------------------
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  AliGRPObject* aGRPData = 0;
  if (entry) {
    TMap* m = dynamic_cast<TMap*>(entry->GetObject()); // old GRP entry

    if (m) {
      AliInfo("Found a TMap in GRP/GRP/Data, converting it into an AliGRPObject");
      m->Print();
      aGRPData = new AliGRPObject();
      aGRPData->ReadValuesFromMap(m);
    }

    else {
      AliInfo("Found an AliGRPObject in GRP/GRP/Data, reading it");
      aGRPData = dynamic_cast<AliGRPObject*>(entry->GetObject()); // new GRP entry
      entry->SetOwner(0);
    }
  }

  if (!aGRPData) {
    AliError("No GRP entry found in OCDB!");
    return;
  }
  //*** Dealing with the magnetic field map

  TString lhcState = aGRPData->GetLHCState();
  if (lhcState == AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the LHC state ! Using UNKNOWN");
    lhcState = "UNKNOWN";
  }

  TString beamType = aGRPData->GetBeamType();
  if (beamType == AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
    beamType = "UNKNOWN";
  }

  float beamEnergy = aGRPData->GetBeamEnergy();
  if (beamEnergy == AliGRPObject::GetInvalidFloat()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
    beamEnergy = 0;
  }

  TString runType = aGRPData->GetRunType();
  if (runType == AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the run type ! Using UNKNOWN");
    runType = "UNKNOWN";
  }

  if (TGeoGlobalMagField::Instance()->IsLocked()) {
    if (TGeoGlobalMagField::Instance()->GetField()->TestBit(AliMagF::kOverrideGRP)) {
      AliInfo("PHOSEmbedding: MF is locked - GRP information will be ignored !");
      AliInfo("Running with the externally locked B field !");
    } else {
      AliInfo("Destroying existing B field instance!");
      delete TGeoGlobalMagField::Instance();
    }
  }
  if (!TGeoGlobalMagField::Instance()->IsLocked()) {
    // Construct the field map out of the information retrieved from GRP.
    bool ok = kTRUE;
    // L3
    float l3Current = aGRPData->GetL3Current((AliGRPObject::Stats)0);
    if (l3Current == AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 current !");
      ok = false;
    }

    Char_t l3Polarity = aGRPData->GetL3Polarity();
    if (l3Polarity == AliGRPObject::GetInvalidChar()) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 polarity !");
      ok = false;
    }
    fMF = (int)l3Polarity;

    // Dipole
    float diCurrent = aGRPData->GetDipoleCurrent((AliGRPObject::Stats)0);
    if (diCurrent == AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole current !");
      ok = false;
    }

    Char_t diPolarity = aGRPData->GetDipolePolarity();
    if (diPolarity == AliGRPObject::GetInvalidChar()) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole polarity !");
      ok = false;
    }

    // read special bits for the polarity convention and map type
    int polConvention = aGRPData->IsPolarityConventionLHC() ? AliMagF::kConvLHC : AliMagF::kConvDCS2008;
    bool uniformB = aGRPData->IsUniformBMap();

    if (ok) {
      AliMagF* fld = AliMagF::CreateFieldMap(TMath::Abs(l3Current) * (l3Polarity ? -1 : 1),
                                             TMath::Abs(diCurrent) * (diPolarity ? -1 : 1), polConvention, uniformB,
                                             beamEnergy, beamType.Data());
      if (fld) {
        TGeoGlobalMagField::Instance()->SetField(fld);
        TGeoGlobalMagField::Instance()->Lock();
        AliInfo("Running with the B field constructed out of GRP !");
      } else
        AliFatal("Failed to create a B field map !");
    } else
      AliFatal("B field is neither set nor constructed from GRP ! Exitig...");
  }
}
//___________________________________________________________________________________________________
int AliPHOSEmbeddingRun2::FindTrackMatching(int mod, TVector3* locpos, double& dx, double& dz, double& pt, int& charge)
{
  // Find track with closest extrapolation to cluster
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(InputEvent());

  double magF = aod->GetMagneticField();

  double magSign = 1.0;
  if (magF < 0)
    magSign = -1.0;

  if (!TGeoGlobalMagField::Instance()->GetField()) {
    AliError("Margnetic filed was not initialized, use default");
    AliMagF* field = new AliMagF("Maps", "Maps", magSign, magSign, AliMagF::k5kG);
    TGeoGlobalMagField::Instance()->SetField(field);
  }

  // *** Start the matching
  int nt = aod->GetNumberOfTracks();

  // Calculate actual distance to PHOS module
  TVector3 globaPos;
  fPHOSGeo->Local2Global(mod, 0., 0., globaPos);
  const double rPHOS = globaPos.Pt();               // Distance to center of  PHOS module
  const double kYmax = 72. + 10.;                   // Size of the module (with some reserve) in phi direction
  const double kZmax = 64. + 10.;                   // Size of the module (with some reserve) in z direction
  const double kAlpha0 = 330. / 180. * TMath::Pi(); // First PHOS module angular direction
  const double kAlpha = 20. / 180. * TMath::Pi();   // PHOS module angular size
  double minDistance = 1.e6;

  double gposTrack[3];

  double bz = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->SolenoidField();
  bz = TMath::Sign(0.5 * kAlmost0Field, bz) + bz;

  double b[3];
  int itr = -1;
  AliAODTrack* aodTrack = 0x0;
  double xyz[3] = { 0 }, pxpypz[3] = { 0 }, cv[21] = { 0 };
  for (int i = 0; i < nt; i++) {
    aodTrack = (AliAODTrack*)aod->GetTrack(i);

    // Continue extrapolation from TPC outer surface
    AliExternalTrackParam outerParam;
    aodTrack->GetPxPyPz(pxpypz);
    aodTrack->GetXYZ(xyz);
    aodTrack->GetCovarianceXYZPxPyPz(cv);
    outerParam.Set(xyz, pxpypz, cv, aodTrack->Charge());

    double z;
    if (!outerParam.GetZAt(rPHOS, bz, z))
      continue;

    if (TMath::Abs(z) > kZmax)
      continue; // Some tracks miss the PHOS in Z

    // Direction to the current PHOS module
    double phiMod = kAlpha0 - kAlpha * mod;
    if (!outerParam.RotateParamOnly(phiMod))
      continue; // RS use faster rotation if errors are not needed

    double y; // Some tracks do not reach the PHOS
    if (!outerParam.GetYAt(rPHOS, bz, y))
      continue; //    because of the bending

    if (TMath::Abs(y) < kYmax) {
      outerParam.GetBxByBz(b);
      outerParam.PropagateToBxByBz(rPHOS, b); // Propagate to the matching module
      // outerParam.CorrectForMaterial(...); // Correct for the TOF material, if needed
      outerParam.GetXYZ(gposTrack);
      TVector3 globalPositionTr(gposTrack);
      TVector3 localPositionTr;
      fPHOSGeo->Global2Local(localPositionTr, globalPositionTr, mod);
      double ddx = locpos->X() - localPositionTr.X();
      double ddz = locpos->Z() - localPositionTr.Z();
      double d2 = ddx * ddx + ddz * ddz;
      if (d2 < minDistance) {
        dx = ddx;
        dz = ddz;
        minDistance = d2;
        itr = i;
        pt = aodTrack->Pt();
        charge = aodTrack->Charge();
      }
    }
  } // Scanned all tracks

  return itr;
}
//____________________________________________________________________________
float AliPHOSEmbeddingRun2::TestCPV(double dx, double dz, double pt, int charge)
{
  // Parameterization of LHC10h period
  //_true if neutral_

  double meanX = 0;
  double meanZ = 0.;
  double sx = TMath::Min(5.4, 2.59719e+02 * TMath::Exp(-pt / 1.02053e-01) +
                                6.58365e-01 * 5.91917e-01 * 5.91917e-01 /
                                  ((pt - 9.61306e-01) * (pt - 9.61306e-01) + 5.91917e-01 * 5.91917e-01) +
                                1.59219);
  double sz = TMath::Min(2.75, 4.90341e+02 * 1.91456e-02 * 1.91456e-02 / (pt * pt + 1.91456e-02 * 1.91456e-02) + 1.60);

  if (fMF < 0.) { // field --
    meanZ = -0.468318;
    if (charge > 0)
      meanX = TMath::Min(7.3, 3.89994 * 1.20679 * 1.20679 / (pt * pt + 1.20679 * 1.20679) + 0.249029 +
                                2.49088e+07 * TMath::Exp(-pt * 3.33650e+01));
    else
      meanX = -TMath::Min(7.7, 3.86040 * 0.912499 * 0.912499 / (pt * pt + 0.912499 * 0.912499) + 1.23114 +
                                 4.48277e+05 * TMath::Exp(-pt * 2.57070e+01));
  } else { // Field ++
    meanZ = -0.468318;
    if (charge > 0)
      meanX = -TMath::Min(8.0, 3.86040 * 1.31357 * 1.31357 / (pt * pt + 1.31357 * 1.31357) + 0.880579 +
                                 7.56199e+06 * TMath::Exp(-pt * 3.08451e+01));
    else
      meanX = TMath::Min(6.85, 3.89994 * 1.16240 * 1.16240 / (pt * pt + 1.16240 * 1.16240) - 0.120787 +
                                 2.20275e+05 * TMath::Exp(-pt * 2.40913e+01));
  }

  double rz = (dz - meanZ) / sz;
  double rx = (dx - meanX) / sx;
  return TMath::Sqrt(rx * rx + rz * rz);
}
//____________________________________________________________________________
float AliPHOSEmbeddingRun2::TestCPVRun2(double dx, double dz, double pt, int charge)
{
  // Parameterization of LHC15o period
  //_true if neutral_

  double meanX = 0.;
  double meanZ = 0.;

  double sx = TMath::Min(5.2, 1.160 + 0.52 * TMath::Exp(-0.042 * pt * pt) + 5.1 / TMath::Power(pt + 0.62, 3));
  double sz = TMath::Min(3.3, 1.10 + 0.39 * TMath::Exp(-0.027 * pt * pt) + 0.70 / TMath::Power(pt + 0.223, 3));

  AliESDEvent* event = dynamic_cast<AliESDEvent*>(InputEvent());
  double mf = event->GetMagneticField(); // Positive for ++ and negative for --

  if (mf < 0.) { // field --
    meanZ = 0.077;
    if (charge > 0)
      meanX = TMath::Min(5.8, 0.2 + 0.7 * TMath::Exp(-0.019 * pt * pt) + 34. / TMath::Power(pt + 1.39, 3));
    else
      meanX = -TMath::Min(5.8, 0.1 + 0.7 * TMath::Exp(-0.014 * pt * pt) + 30. / TMath::Power(pt + 1.36, 3));
  } else { // Field ++
    meanZ = 0.077;
    if (charge > 0)
      meanX = -TMath::Min(5.8, 0.3 + 0.7 * TMath::Exp(-0.012 * pt * pt) + 35. / TMath::Power(pt + 1.43, 3));
    else
      meanX = TMath::Min(5.8, 0.2 + 0.6 * TMath::Exp(-0.014 * pt * pt) + 28. / TMath::Power(pt + 1.27, 3));
  }

  double rz = (dz - meanZ) / sz;
  double rx = (dx - meanX) / sx;
  return TMath::Sqrt(rx * rx + rz * rz);
}
//________________________________________________________________________
bool AliPHOSEmbeddingRun2::IsGoodChannel(int mod, int ix, int iz)
{
  // Check if this channel belogs to the good ones

  if (!fPHOSBadMap[mod]) {
    AliError(Form("No Bad map for PHOS module %d", mod));
    return kFALSE;
  }
  if (fPHOSBadMap[mod]->GetBinContent(ix, iz) > 0)
    return kFALSE;
  else
    return kTRUE;
}
//________________________________________________________________________
double AliPHOSEmbeddingRun2::CorrectNonlinearity(double en)
{
  // Non-linearity correction in real data
  if (en <= 0.)
    return 0.;
  const double xMin = 0.36; // low part of the param (optimized from pi0 peak)
  const double xMax = 5.17; // Upper part of the param (optimized from pi0 peak)

  // middle part param
  const double a = 1.02165;
  const double b = -2.548e-01;
  const double c = 6.483e-01;
  const double d = -0.4805;
  const double e = 0.1275;

  double ecorr = 0.;
  if (en < xMin) {
    const double beta = 2. * a * sqrt(xMin) + b - d / (xMin)-2. * e / (xMin * sqrt(xMin));
    const double alpha = a * xMin + b * sqrt(xMin) + c + d / sqrt(xMin) + e / xMin - beta * sqrt(xMin);
    ecorr = 1.0312526 * (alpha + beta * sqrt(en));
  } else {
    if (en < xMax) {
      ecorr = 1.0312526 * (a * en + b * sqrt(en) + c + d / sqrt(en) + e / en);
    } else {
      const double beta = b + 2. * c / sqrt(xMax) + 3. * d / xMax + 4. * e / xMax / sqrt(xMax);
      const double alpha =
        a + b / sqrt(xMax) + c / xMax + d / xMax / sqrt(xMax) + e / (xMax * xMax) - beta / sqrt(xMax);
      ecorr = 1.0312526 * (alpha * en + beta * sqrt(en));
    }
  }
  if (ecorr > 0) {
    return ecorr;
  } else {
    return 0.;
  }
}
//________________________________________________________________________
double AliPHOSEmbeddingRun2::CorrectNonlinearityMC(double en)
{
  // Nonlinearity correction of MC
  if (en <= 0.)
    return 0.;

  const double xMin = 0.850; // low part of the param (optimized from pi0 peak)
  const double xMax = 5.17;  // Upper part of the param (optimized from pi0 peak)

  // middle part param
  const double a = 1.02165;
  const double b = -2.548e-01;
  const double c = 0.6483;
  const double d = -0.4980;
  const double e = 0.1245;

  double ecorr = 0.;
  if (en < xMin) {
    const double gamma = 0.150;
    const double beta = 0.5 * (0.5 * b * sqrt(xMin) + c + 1.5 * d / sqrt(xMin) + 2. * e / xMin) *
                        TMath::Power((xMin * xMin + gamma * gamma), 2) / (xMin * xMin * xMin);
    const double alpha =
      (a * xMin + b * sqrt(xMin) + c + d / sqrt(xMin) + e / xMin - beta * xMin / (xMin * xMin + gamma * gamma)) / xMin;
    ecorr = 1.0328783 * (alpha * en + beta * en / (en * en + gamma * gamma));
  } else {
    if (en < xMax) {
      ecorr = 1.0328783 * (a * en + b * sqrt(en) + c + d / sqrt(en) + e / en);
    } else {
      const double beta = b + 2. * c / sqrt(xMax) + 3. * d / xMax + 4. * e / xMax / sqrt(xMax);
      const double alpha =
        a + b / sqrt(xMax) + c / xMax + d / xMax / sqrt(xMax) + e / (xMax * xMax) - beta / sqrt(xMax);
      ecorr = 1.0328783 * (alpha * en + beta * sqrt(en));
    }
  }
  if (ecorr < 0) {
    return 0.;
  }

  return ecorr;
}
//____________________________________________________________________________
double AliPHOSEmbeddingRun2::CoreEnergy(AliVCluster* clu)
{
  // calculate energy of the cluster in the circle with radius distanceCut around the maximum

  // Can not use already calculated coordinates?
  // They have incidence correction...
  const double distanceCut = 3.5;
  const double logWeight = 4.5;

  Double32_t* elist = clu->GetCellsAmplitudeFraction();
  // Calculates the center of gravity in the local PHOS-module coordinates
  float wtot = 0;
  double xc[100] = { 0 };
  double zc[100] = { 0 };
  double x = 0;
  double z = 0;
  int mulDigit = TMath::Min(100, clu->GetNCells());
  for (int iDigit = 0; iDigit < mulDigit; iDigit++) {
    int relid[4];
    float xi;
    float zi;
    fPHOSGeo->AbsToRelNumbering(clu->GetCellAbsId(iDigit), relid);
    fPHOSGeo->RelPosInModule(relid, xi, zi);
    xc[iDigit] = xi;
    zc[iDigit] = zi;
    if (clu->E() > 0 && elist[iDigit] > 0) {
      float w = TMath::Max(0., logWeight + TMath::Log(elist[iDigit] / clu->E()));
      x += xc[iDigit] * w;
      z += zc[iDigit] * w;
      wtot += w;
    }
  }
  if (wtot > 0) {
    x /= wtot;
    z /= wtot;
  }
  double coreE = 0.;
  for (int iDigit = 0; iDigit < mulDigit; iDigit++) {
    double distance = TMath::Sqrt((xc[iDigit] - x) * (xc[iDigit] - x) + (zc[iDigit] - z) * (zc[iDigit] - z));
    if (distance < distanceCut)
      coreE += elist[iDigit];
  }
  // Apply non-linearity correction
  return coreE;
}
//________________________________________________________________________
void AliPHOSEmbeddingRun2::EvalLambdas(AliVCluster* clu, double& m02, double& m20)
{
  // calculate dispecrsion of the cluster in the circle with radius distanceCut around the maximum

  const double rCut = 4.5;

  Double32_t* elist = clu->GetCellsAmplitudeFraction();
  // Calculates the center of gravity in the local PHOS-module coordinates
  float wtot = 0;
  double xc[100] = { 0 };
  double zc[100] = { 0 };
  double x = 0;
  double z = 0;
  int mulDigit = TMath::Min(100, clu->GetNCells());
  const double logWeight = 4.5;
  for (int iDigit = 0; iDigit < mulDigit; iDigit++) {
    int relid[4];
    float xi;
    float zi;
    fPHOSGeo->AbsToRelNumbering(clu->GetCellAbsId(iDigit), relid);
    fPHOSGeo->RelPosInModule(relid, xi, zi);
    xc[iDigit] = xi;
    zc[iDigit] = zi;
    if (clu->E() > 0 && elist[iDigit] > 0) {
      float w = TMath::Max(0., logWeight + TMath::Log(elist[iDigit] / clu->E()));
      x += xc[iDigit] * w;
      z += zc[iDigit] * w;
      wtot += w;
    }
  }
  if (wtot > 0) {
    x /= wtot;
    z /= wtot;
  }

  wtot = 0.;
  double dxx = 0.;
  double dzz = 0.;
  double dxz = 0.;
  double xCut = 0.;
  double zCut = 0.;
  for (int iDigit = 0; iDigit < mulDigit; iDigit++) {
    if (clu->E() > 0 && elist[iDigit] > 0.) {
      double w = TMath::Max(0., logWeight + TMath::Log(elist[iDigit] / clu->E()));
      double xi = xc[iDigit];
      double zi = zc[iDigit];
      if ((xi - x) * (xi - x) + (zi - z) * (zi - z) < rCut * rCut) {
        xCut += w * xi;
        zCut += w * zi;
        dxx += w * xi * xi;
        dzz += w * zi * zi;
        dxz += w * xi * zi;
        wtot += w;
      }
    }
  }
  if (wtot > 0) {
    xCut /= wtot;
    zCut /= wtot;
    dxx /= wtot;
    dzz /= wtot;
    dxz /= wtot;
    dxx -= xCut * xCut;
    dzz -= zCut * zCut;
    dxz -= xCut * zCut;

    m02 = 0.5 * (dxx + dzz) + TMath::Sqrt(0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz);
    m20 = 0.5 * (dxx + dzz) - TMath::Sqrt(0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz);
  } else {
    m20 = m02 = 0.;
  }
}
//_____________________________________________________________________________
double AliPHOSEmbeddingRun2::TestCoreLambda(double pt, double l1, double l2)
{
  // Parameterization for core dispersion
  // For R=4.5
  const double l1Mean = 1.150200 + 0.097886 / (1. + 1.486645 * pt + 0.000038 * pt * pt);
  const double l2Mean = 1.574706 + 0.997966 * exp(-0.895075 * pt) - 0.010666 * pt;
  const double l1Sigma = 0.100255 + 0.337177 * exp(-0.517684 * pt) + 0.001170 * pt;
  const double l2Sigma = 0.232580 + 0.573401 * exp(-0.735903 * pt) - 0.002325 * pt;
  const double c = -0.110983 - 0.017353 / (1. - 1.836995 * pt + 0.934517 * pt * pt);

  double R2 = 0.5 * (l1 - l1Mean) * (l1 - l1Mean) / l1Sigma / l1Sigma +
              0.5 * (l2 - l2Mean) * (l2 - l2Mean) / l2Sigma / l2Sigma +
              0.5 * c * (l1 - l1Mean) * (l2 - l2Mean) / l1Sigma / l2Sigma;
  return R2;
}
//_____________________________________________________________________________
double AliPHOSEmbeddingRun2::TestFullLambda(double pt, double l1, double l2)
{
  // Parameterization for full dispersion
  // Parameterizatino for full dispersion
  const double l2Mean = 1.53126 + 9.50835e+06 / (1. + 1.08728e+07 * pt + 1.73420e+06 * pt * pt);
  const double l1Mean = 1.12365 + 0.123770 * TMath::Exp(-pt * 0.246551) + 5.30000e-03 * pt;
  const double l2Sigma = 6.48260e-02 + 7.60261e+10 / (1. + 1.53012e+11 * pt + 5.01265e+05 * pt * pt) + 9.00000e-03 * pt;
  const double l1Sigma = 4.44719e-04 + 6.99839e-01 / (1. + 1.22497e+00 * pt + 6.78604e-07 * pt * pt) + 9.00000e-03 * pt;
  const double c = -0.35 - 0.550 * TMath::Exp(-0.390730 * pt);

  double R2 = 0.5 * (l1 - l1Mean) * (l1 - l1Mean) / l1Sigma / l1Sigma +
              0.5 * (l2 - l2Mean) * (l2 - l2Mean) / l2Sigma / l2Sigma +
              0.5 * c * (l1 - l1Mean) * (l2 - l2Mean) / l1Sigma / l2Sigma;
  return R2;
}
//________________________________________________________________________
double AliPHOSEmbeddingRun2::EvalTOF(AliVCluster* clu, AliVCaloCells* cells)
{
  // Evaluate TOF of the cluster after re-calibration
  // TOF here is weighted average of digits
  //  -within 50ns from the most energetic cell
  //  -not too soft.

  Double32_t* elist = clu->GetCellsAmplitudeFraction();
  int mulDigit = clu->GetNCells();

  // Slewing correction from LHC16eghklmn
  const float sA = -2.57668e-09;
  const float sB = 8.19737e-09;
  const float sC = -3.16538e-09;
  const float sD = 1.02124e-09;
  const float sE = -1.11128e-10;

  // Slewing correction LHC15xx  2 values does not depend on trigger.
  const double saturate = -4.42; //-4.42  ± 0.02  for pp at 5 TeV LHC15n
  const double slope1 = 4.82;    // 4.82  ± 0.02  for pp at 5 TeV LHC15n
  const double slope2 = -0.100;  //-0.100 ± 0.006 for pp at 5 TeV LHC15n

  float tMax = 0.; // Time at the maximum
  float eMax = 0.;
  float eMaxHG = 0.; // High Gain only
  // bool isHGMax = kTRUE;
  int absIdMax = -1, absIdMaxHG = -1;
  for (int iDigit = 0; iDigit < mulDigit; iDigit++) {
    int absId = clu->GetCellAbsId(iDigit);
    bool isHG = cells->GetCellHighGain(absId);
    float ei = elist[iDigit];
    if (isHG && (ei > eMaxHG)) { // only HG channels
      eMaxHG = ei;
      absIdMaxHG = absId;
    }
    if (ei > eMax) {
      eMax = ei;
      absIdMax = absId;
      // isHGMax = isHG;
    }
  }
  // Use LG for time
  tMax = cells->GetCellTime(absIdMax);
  // tMax=CalibrateTOF(cells->GetCellTime(absIdMax),absIdMax,isHGMax) ;
  // Slewing correction
  if (eMax > 0 && fRunNumber > 209122) { // Run2
    if (fRunNumber < 252603)             // before LHC16e
      tMax -= (saturate + slope1 / eMax + slope2 / (eMax * eMax)) * 1e-9;
    //          tMax-= sA15+sB15/clu->E();
    else
      tMax -= sA + sB / eMax + sC / eMax / eMax + sD / eMax / eMax / eMax + sE / eMax / eMax / eMax / eMax;
  }

  tMax = cells->GetCellTime(absIdMaxHG);
  // tMax=CalibrateTOF(cells->GetCellTime(absIdMaxHG),absIdMaxHG,kTRUE) ;
  // Slewing correction
  if (eMaxHG > 0 && fRunNumber > 209122) {
    if (fRunNumber < 252603) // before LHC16e
      tMax -= (saturate + slope1 / eMax + slope2 / (eMax * eMax)) * 1e-9;
    //          tMax-= sA15+sB15/clu->E();
    else
      tMax -= sA + sB / eMaxHG + sC / eMaxHG / eMaxHG + sD / eMaxHG / eMaxHG / eMaxHG +
              sE / eMaxHG / eMaxHG / eMaxHG / eMaxHG;
  }

  return tMax;
}
//________________________________________________________________________
void AliPHOSEmbeddingRun2::RunSignalSimulation(int nEvents)
{
  // run MC simulation
  AliInfo(Form("Simulating %d events\n", nEvents));
  // set the seed environment variable
  gSystem->Setenv("CONFIG_SEED", "0");
  gSystem->Setenv("CONFIG_RUN_TYPE", "kPythia6"); // kPythia6 or kPhojet
  gSystem->Setenv("CONFIG_FIELD", "k5kG");        // kNoField or k5kG
  gSystem->Setenv("CONFIG_ENERGY", "5020");
  gSystem->Setenv("DC_RUN", Form("%d", fRunNumber));
  gSystem->Setenv("SIM_EVENTS", Form("%d", nEvents));
  gSystem->Exec("cp $ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_Run2embedding/macros/* .");
  gSystem->Exec("chmod u+x aliroot_gghbtsim.sh");
  gSystem->Exec("alien_cp alien:///alice/cern.ch/user/p/prsnko/Tagging/macros.LHC15o/alidpg.tgz file:");
  gSystem->Exec(
    Form("./aliroot_gghbtsim.sh --run %d --mode simrec  "
         "--generator Custom --detector PhosOnly --simulation PhosOnly --reconstruction PhosOnly --nevents %d --system "
         "PbPb --ocdb alien",
         fRunNumber, nEvents));
  gSystem->Exec("root -q -b CreateAOD.C");
}
