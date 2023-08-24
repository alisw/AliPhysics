/*************** Tree Maker Class **********************
*                                                      *
* Created: 05.10.2016                                  *
* Authors: Aaron Capon      (aaron.capon@cern.ch)      *
*          Sebastian Lehner (sebastian.lehner@cern.ch) *
*          Elisa Meninno (elisa.meninno@cern.ch)       *
*******************************************************/
/**********************************************************************************
*This analysis task is designed to create simple flat structured TTrees which     *
*can then be worked on locally. The motivation for this was to create TTrees that *
*could be easily used within TMVA.                                                *
*                                                                                 *
*The standard TTree relies on the dielectron framework and sources cuts from a    *
"CutLibrary" config file. An example can be found in macrosLMEE/LMEECuLib_acapon.C*
*The desired cut is then specified in the AddTaskSimpleTreeMaker.C.               *
*                                                                                 *
*The GRID PID number for each job is stored with each event along with a simple   *
*event counter ID, so that after merging each event still has a unique label.     *
*                                                                                 *
*By default, the class will return a TTree focused on selecting electrons from    *
*primary decays. There are setter functions which can be used to instead create a *
*TTree filled with tracks originating from V0 decays. In this case the DCA cuts   *
*are removed and standard methods for applying cuts are used.                     *
**********************************************************************************/
#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"

#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"

#include "AliVEvent.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "AliAODMCParticle.h"

#include "TParticle.h"
#include "AliTrackReference.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#include "AliVVertex.h"

#include "AliESDv0.h"
#include "AliDielectronTrackCuts.h"
#include "AliDielectronVarManager.h"
#include "AliDielectronV0Cuts.h"

#include "TClonesArray.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"

#include "TSystem.h"
#include <string>
#include "AliAnalysisTaskSimpleTreeMaker.h"



ClassImp(AliAnalysisTaskSimpleTreeMaker)

AliAnalysisTaskSimpleTreeMaker::AliAnalysisTaskSimpleTreeMaker():
AliAnalysisTaskSE(),
  eventNum(0),
  hasMC(kFALSE),
  fESDtrackCuts(0),
  fPIDResponse(0),
  fMCevent(0),
  fTree(0x0),
  eventCuts(0),
  eventFilter(0),
  varCuts(0),
  trackCuts(0),
  pidCuts(0),
  cuts(0),
  trackFilter(0),
  varManager(0),
  primaryVertex{0,0,0},
  multiplicityV0A(0),
  multiplicityV0C(0),
  multiplicityCL1(0),
  runNumber(0),
  event(0),
  pt(0),
  eta(0),
  phi(0),
  nTPCclusters(0),
  nTPCcrossed(0),
  fTPCcrossOverFind(0),
  nTPCfindable(0),
  tpcSharedMap(0),
  nTPCshared(0),
  chi2TPC(0),
  DCA{0,0},
  DCAsigma{0,0},
  DCA3{0,0,0},
  nITS(0),
  chi2ITS(0),
  fITSshared(0),
  SPDfirst(0),
  charge(0),
  EnSigmaITS(0),
  EnSigmaITScorr(0),
  EnSigmaTPC(0),
  EnSigmaTPCcorr(0),
  EnSigmaTOF(0),
  EnSigmaTOFcorr(0),
  maxPtPIDcorrection(5),
  PnSigmaITS(0),
  PnSigmaTPC(0),
  PnSigmaTOF(0),
  storeKaonPID(0),
  KnSigmaITS(0),
  KnSigmaTPC(0),
  KnSigmaTOF(0),
  storeProtonPID(0),
  PrnSigmaITS(0),
  PrnSigmaTPC(0),
  PrnSigmaTOF(0),
  ITSsignal(0),
  TPCsignal(0),
  TOFsignal(0),
  goldenChi2(0),
  mcEta(0),
  mcPhi(0),
  mcPt(0),
  mcVert{0,0,0},
  iPdg(0),
  iPdgMother(0),
  HasMother(0),
  motherLabel(0),
  isInj(0),
  isPhysPrimary(0),
  iPdgFirstMother(0),
  gLabelFirstMother(0),
  gLabelMinFirstMother(0),
  gLabelMaxFirstMother(0),
  pointingAngle(0),
  daughtersDCA(0),
  decayLength(0),
  v0mass(0),
  ptArm(0),
  alpha(0),
  fQAhist(0),
  fCentralityPercentileMin(0),
  fCentralityPercentileMax(100),
  fPtMin(0.2),
  fPtMax(10),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fESigTPCMin(-4.),
  fESigTPCMax(4.),
  fHasSDD(kTRUE),
  fExtraDCA(kFALSE),
  fIsV0tree(kFALSE),
  fArmPlot(0),
  fIsAOD(kTRUE),
  fFilterBit(16),
  fIsGRIDanalysis(kTRUE),
  fGridPID(-1),
  fUseTPCcorr(kFALSE),
  fWidthTPC(0),
  fMeanTPC(0),
  fUseITScorr(kFALSE),
  fWidthITS(0),
  fMeanITS(0),
  fUseTOFcorr(kFALSE),
  fWidthTOF(0),
  fMeanTOF(0),
  TOFstartMask(0),
  fGeneratorHashes(0)
{

}

AliAnalysisTaskSimpleTreeMaker::AliAnalysisTaskSimpleTreeMaker(const char *name, Bool_t addDCA) :
  AliAnalysisTaskSE(name),
  eventNum(0),
  hasMC(kFALSE),
  fESDtrackCuts(0),
  fPIDResponse(0),
  fMCevent(0),
  fTree(0),
  eventCuts(0),
  eventFilter(0),
  varCuts(0),
  trackCuts(0),
  pidCuts(0),
  cuts(0),
  trackFilter(0),
  varManager(0),
  primaryVertex{0,0,0},
  multiplicityV0A(0),
  multiplicityV0C(0),
  multiplicityCL1(0),
  runNumber(0),
  event(0),
  pt(0),
  eta(0),
  phi(0),
  nTPCclusters(0),
  nTPCcrossed(0),
  fTPCcrossOverFind(0),
  nTPCfindable(0),
  tpcSharedMap(0),
  nTPCshared(0),
  chi2TPC(0),
  DCA{0,0},
  DCAsigma{0,0},
  DCA3{0,0,0},
  nITS(0),
  chi2ITS(0),
  fITSshared(0),
  SPDfirst(0),
  charge(0),
  EnSigmaITS(0),
  EnSigmaITScorr(0),
  EnSigmaTPC(0),
  EnSigmaTPCcorr(0),
  EnSigmaTOF(0),
  EnSigmaTOFcorr(0),
  maxPtPIDcorrection(5),
  PnSigmaITS(0),
  PnSigmaTPC(0),
  PnSigmaTOF(0),
  storeKaonPID(0),
  KnSigmaITS(0),
  KnSigmaTPC(0),
  KnSigmaTOF(0),
  storeProtonPID(0),
  PrnSigmaITS(0),
  PrnSigmaTPC(0),
  PrnSigmaTOF(0),
  ITSsignal(0),
  TPCsignal(0),
  TOFsignal(0),
  goldenChi2(0),
  mcEta(0),
  mcPhi(0),
  mcPt(0),
  mcVert{0,0,0},
  iPdg(0),
  iPdgMother(0),
  HasMother(0),
  motherLabel(0),
  isInj(0),
  isPhysPrimary(0),
  iPdgFirstMother(0),
  gLabelFirstMother(0),
  gLabelMinFirstMother(0),
  gLabelMaxFirstMother(0),
  pointingAngle(0),
  daughtersDCA(0),
  decayLength(0),
  v0mass(0),
  ptArm(0),
  alpha(0),
  fQAhist(0),
  fCentralityPercentileMin(0),
  fCentralityPercentileMax(100),
  fPtMin(0.2),
  fPtMax(10),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fESigTPCMin(-4.),
  fESigTPCMax(4.),
  fHasSDD(kTRUE),
  fExtraDCA(addDCA),
  fIsV0tree(kFALSE),
  fArmPlot(0),
  fIsAOD(kTRUE),
  fFilterBit(16),
  fIsGRIDanalysis(kTRUE),
  fGridPID(0),
  fUseTPCcorr(kFALSE),
  fWidthTPC(0),
  fMeanTPC(0),
  fUseITScorr(kFALSE),
  fWidthITS(0),
  fMeanITS(0),
  fUseTOFcorr(kFALSE),
  fWidthTOF(0),
  fMeanTOF(0),
  TOFstartMask(0),
  fGeneratorHashes(0)

{
  if(fIsV0tree){
    fESDtrackCuts = AliESDtrackCuts::GetStandardV0DaughterCuts();
  }

  // Create hashes of generators used for injected signals in MC
  TString generatorNames = "Pythia CC_1;Pythia BB_1;Pythia B_1;Jpsi2ee_1;B2Jpsi2ee_1";
  TObjArray arr = *(generatorNames.Tokenize(";"));
  for(Int_t i = 0; i < arr.GetEntries(); i++){
    TString temp = arr.At(i)->GetName();
    std::cout << "---" << temp << std::endl;
    fGeneratorHashes.push_back(temp.Hash());
  }

  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineOutput(1, TTree::Class()); //will be connected to fTree
  DefineOutput(2, TH1F::Class());
  DefineOutput(3, TH2F::Class());

}

//________________________________________________________________________

AliAnalysisTaskSimpleTreeMaker::~AliAnalysisTaskSimpleTreeMaker(){

  delete eventCuts;
  delete eventFilter;

  delete varCuts;
  delete trackCuts;
  delete pidCuts;
  delete cuts;
  delete trackFilter;

}


//________________________________________________________________________

void AliAnalysisTaskSimpleTreeMaker::UserCreateOutputObjects(){

  AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  inputHandler->SetNeedField();

  fPIDResponse = inputHandler->GetPIDResponse();
  if(!fPIDResponse){
    AliFatal("This task needs the PID response attached to the inputHandler");
    return;
  }

  fTree = new TTree("tracks", "tracks");

  // Track variables
  fTree->Branch("pt",                &pt,                "pt/F");
  fTree->Branch("eta",               &eta,               "eta/F");
  fTree->Branch("phi",               &phi,               "phi/F");
  fTree->Branch("nTPCclusters",      &nTPCclusters,      "nTPCclusters/F");
  fTree->Branch("nTPCcrossed",       &nTPCcrossed,       "nTPCcrossed/F");
  fTree->Branch("nTPCfindable",      &nTPCfindable,      "nTPCfindable/F");
  fTree->Branch("nTPCshared",        &nTPCshared,        "nTPCshared/F");
  fTree->Branch("fTPCcrossOverFind", &fTPCcrossOverFind, "fTPCcrossOverFind/F");
  fTree->Branch("chi2TPC",           &chi2TPC,           "chi2TPC/F");
  fTree->Branch("nITS",              &nITS,              "nITS/I");
  fTree->Branch("fITSshared",        &fITSshared,        "fITSshared/F");
  fTree->Branch("chi2ITS",           &chi2ITS,           "chi2ITS/F");
  fTree->Branch("SPDfirst",          &SPDfirst,          "SPDfirst/O");
  fTree->Branch("DCAxy",             &DCA[0],            "DCAxy/F");
  fTree->Branch("DCAz",              &DCA[1],            "DCAz/F");
  if(fExtraDCA){
    fTree->Branch("DCAsigmaxy",        &DCAsigma[0],    "DCAsigmaxy/F");
    fTree->Branch("DCAsigmaz",         &DCAsigma[1],    "DCAsigmaz/F");
    fTree->Branch("DCA3x",             &DCA3[0],        "DCA3x/F");
    fTree->Branch("DCA3y",             &DCA3[1],        "DCA3y/F");
    fTree->Branch("DCA3z",             &DCA3[2],        "DCA3z/F");
  }
  fTree->Branch("goldenChi2",        &goldenChi2,     "goldenChi2/F");
  fTree->Branch("charge",            &charge,         "charge/I");
  fTree->Branch("EsigITS",           &EnSigmaITS,     "EsigITS/F");
  if(fUseITScorr){
    fTree->Branch("EsigITScorr",     &EnSigmaITScorr, "EsigITScorr/F");
  }
  fTree->Branch("EsigTPC",           &EnSigmaTPC,     "EsigTPC/F");
  if(fUseTPCcorr){
    fTree->Branch("EsigTPCcorr",     &EnSigmaTPCcorr, "EsigTPCcorr/F");
  }
  fTree->Branch("EsigTOF",           &EnSigmaTOF,     "EsigTOF/F");
  if(fUseTOFcorr){
    fTree->Branch("EsigTOFcorr",     &EnSigmaTOFcorr, "EsigTOFcorr/F");
  }
  fTree->Branch("PsigITS",           &PnSigmaITS, "PsigITS/F");
  fTree->Branch("PsigTPC",           &PnSigmaTPC, "PsigTPC/F");
  fTree->Branch("PsigTOF",           &PnSigmaTOF, "PsigTOF/F");
  if(storeKaonPID){
    fTree->Branch("KsigITS",           &KnSigmaITS, "KsigITS/F");
    fTree->Branch("KsigTPC",           &KnSigmaTPC, "KsigTPC/F");
    fTree->Branch("KsigTOF",           &KnSigmaTOF, "KsigTOF/F");
  }
  if(storeProtonPID){
    fTree->Branch("PrsigITS",           &PrnSigmaITS, "PrsigITS/F");
    fTree->Branch("PrsigTPC",           &PrnSigmaTPC, "PrsigTPC/F");
    fTree->Branch("PrsigTOF",           &PrnSigmaTOF, "PrsigTOF/F");
  }
  fTree->Branch("ITSsignal",         &ITSsignal,  "ITSsignal/F");
  fTree->Branch("TPCsignal",         &TPCsignal,  "TPCsignal/F");
  fTree->Branch("TOFsignal",         &TOFsignal,  "TOFsignal/F");
  // V0 tree variables
  if(fIsV0tree){
    fTree->Branch("pointingAngle", &pointingAngle, "pointingAngle/F");
    fTree->Branch("daughtersDCA",  &daughtersDCA,  "daughtersDCA/F");
    fTree->Branch("decayLength",   &decayLength,   "decayLength/F");
    fTree->Branch("v0mass",        &v0mass,        "v0mass/F");
    fTree->Branch("ptArm",         &ptArm,         "ptArm/F");
    fTree->Branch("alpha",         &alpha,         "alpha/F");
  }
  // MC variables
  if(hasMC){
    fTree->Branch("mcPt",              &mcPt,                 "mcPt/F");
    fTree->Branch("mcEta",             &mcEta,                "mcEta/F");
    fTree->Branch("mcPhi",             &mcPhi,                "mcPhi/F");
    fTree->Branch("mcVertx",           &mcVert[0],            "mcVertx/F");
    fTree->Branch("mcVerty",           &mcVert[1],            "mcVerty/F");
    fTree->Branch("mcVertz",           &mcVert[2],            "mcVertz/F");
    fTree->Branch("pdg",               &iPdg,                 "pdg/I");
    fTree->Branch("pdgMother",         &iPdgMother,           "pdgMother/I");
    fTree->Branch("HasMother",         &HasMother,            "HasMother/O");
    fTree->Branch("motherLabel",       &motherLabel,          "motherLabel/I");
    fTree->Branch("pdgFirstMother",    &iPdgFirstMother,      "pdgFirstMother/I");
    fTree->Branch("gLabelFirstMother", &gLabelFirstMother,    "gLabelFirstMother/I");
    fTree->Branch("labelMinInitial",   &gLabelMinFirstMother, "labelMinInitial/I");
    fTree->Branch("labelMaxInitial",   &gLabelMaxFirstMother, "labelMaxInitial/I");
    fTree->Branch("isInjected",        &isInj,                "isInjected/I");
    fTree->Branch("isPhysPrim",        &isPhysPrimary,        "isPhysPrim/O");
  }
  // Event variables
  fTree->Branch("vertexX",         &primaryVertex[0], "vertexX/F");
  fTree->Branch("vertexY",         &primaryVertex[1], "vertexY/F");
  fTree->Branch("vertexZ",         &primaryVertex[2], "vertexZ/F");
  fTree->Branch("runNumber",       &runNumber,        "runNumber/I");
  fTree->Branch("eventNum",        &eventNum,         "eventNum/I");
  fTree->Branch("multiplicityV0A", &multiplicityV0A,  "multiplicityV0A/F");
  fTree->Branch("multiplicityV0C", &multiplicityV0C,  "multiplicityV0C/F");
  fTree->Branch("multiplicityCL1", &multiplicityCL1,  "multiplicityCL1/F");
  fTree->Branch("gridPID",         &fGridPID,         "gridPID/I");
  fTree->Branch("TOFstartMask",    &TOFstartMask,     "TOFstartMask/I");

  // Get grid PID which can be used later to assign unique event numbers
  if(fIsGRIDanalysis){
    const char* gridIDchar = gSystem->Getenv("ALIEN_PROC_ID");
    std::string str(gridIDchar);
    SetGridPID(str);
  }
  else{
    fGridPID = -1;
  }

  // Setup correction maps
  if(fUseTPCcorr){
    AliDielectronPID::SetCentroidCorrFunction( (TH1*)fMeanTPC->Clone() );
    AliDielectronPID::SetWidthCorrFunction( (TH1*)fWidthTPC->Clone() );
    ::Info("AliAnalysisTaskSimpleTreeMaker::UserExec","Setting TPC Correction Histos");
  }

  if(fUseITScorr){
    AliDielectronPID::SetCentroidCorrFunctionITS( (TH1*)fMeanITS->Clone() );
    AliDielectronPID::SetWidthCorrFunctionITS( (TH1*)fWidthITS->Clone() );
    ::Info("AliAnalysisTaskSimpleTreeMaker::UserExec","Setting ITS Correction Histos");
  }

  if(fUseTOFcorr){
    AliDielectronPID::SetCentroidCorrFunctionTOF( (TH1*)fMeanTOF->Clone() );
    AliDielectronPID::SetWidthCorrFunctionTOF( (TH1*)fWidthTOF->Clone() );
    ::Info("AliAnalysisTaskSimpleTreeMaker::UserExec","Setting TOF Correction Histos");
  }

  // Needed by the dielectron framework
  varManager->SetPIDResponse(fPIDResponse);

  // Create TH2F for armenteros plot. Filled if creating v0 tree
  // NOTE: Class designed to study electrons with weak EsigTPC cuts
  // applied. Therefore, the number of tracks will not necessarily be twice the
  // number of V0 particles (as one might expect).
  fArmPlot = new TH2F("ArmPlot", "Armenteros Plot", 100, -1, 1, 100, 0, 0.4);
  if(fIsV0tree){
    fArmPlot->GetXaxis()->SetTitle("#alpha = (p^{+}-p^{-})/(p^{+}+p^{-})");
    fArmPlot->GetYaxis()->SetTitle("p_{T}");
  }

  fQAhist = new TH1F("h1", "Event and track QA", 8, 0, 8);
  // Fill each bin with nothing to ensure correct ordering
  fQAhist->Fill("Events_check",0);
  fQAhist->Fill("Events_accepted",0);
  fQAhist->Fill("Events_MCcheck",0);
  fQAhist->Fill("Tracks_all",0);
  fQAhist->Fill("Tracks_MCcheck",0);
  fQAhist->Fill("Tracks_Cuts",0);

  PostData(1, fTree);
  PostData(2, fQAhist);
  PostData(3, fArmPlot);
}

//________________________________________________________________________
// Main loop: called for each event
void AliAnalysisTaskSimpleTreeMaker::UserExec(Option_t *){

  AliVEvent* event = 0x0;
  AliESDEvent* esdEvent = 0x0;
  if(!fIsV0tree){
    event = dynamic_cast<AliVEvent*>(InputEvent());
  }else{
    esdEvent = dynamic_cast<AliESDEvent*>(InputEvent());
    event = esdEvent;
  }
  if(!event){
    AliError("No event");
    return;
  }
  if(!fIsV0tree){
    fQAhist->Fill("Events_check",1);
  }

  // Check event cuts
  if(IsEventAccepted(event) == 0){
    return;
  }
  if(!fIsV0tree){
    fQAhist->Fill("Events_accepted",1);
  }

  // Check if running on MC files
  fMCevent = MCEvent();
  if(fMCevent){
    hasMC = kTRUE;
    if(!fIsV0tree){
      fQAhist->Fill("Events_MCcheck",1);
    }
  }
  else{
    hasMC = kFALSE;
  }

  eventNum += 1;

  // PID Response task active?
  fPIDResponse = (dynamic_cast<AliInputEventHandler*>((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())))->GetPIDResponse();

  if(!fPIDResponse){ AliFatal("This task needs the PID response attached to the inputHandler"); }

  AliVVertex* vertex = const_cast<AliVVertex*>(event->GetPrimaryVertex());

  primaryVertex[0] = vertex->GetX();
  primaryVertex[1] = vertex->GetY();
  primaryVertex[2] = vertex->GetZ();

  // Get Multiplicity
  AliMultSelection* multSelection = dynamic_cast<AliMultSelection*>(event->FindListObject("MultSelection"));
  if(!multSelection){
    AliWarning("AliMultSelection object not found");
  }
  else{
    multiplicityV0A = multSelection->GetMultiplicityPercentile("V0A");
    multiplicityV0C = multSelection->GetMultiplicityPercentile("V0C");
    multiplicityCL1 = multSelection->GetMultiplicityPercentile("CL1");
  }

  // Check if ESD or AOD analysis
  // Get ClassName from file and set appropriate ESD or AOD flag
  TString className = static_cast<TString>(event->ClassName());
  if(className.Contains("AOD")){
    fIsAOD = kTRUE;
  }
  else if(className.Contains("ESD")){
    fIsAOD = kFALSE;
  }
  else{
    AliError("-----!! Analysis type unknown !!--------");
    return;
  }


  runNumber         = event->GetRunNumber();
  Int_t eventTracks = event->GetNumberOfTracks();
  Int_t numV0s      = event->GetNumberOfV0s();

  AliVParticle* mcTrack = 0x0;
  AliVParticle* motherMCtrack = 0x0;
  //Printf("myclass!!");
  // Loop over tracks for event
  if(!fIsV0tree){
    for(Int_t iTrack = 0; iTrack < eventTracks; iTrack++){

      AliVTrack* track = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
      if(!track){
        AliError(Form("Could not receive track %d", iTrack));
        continue;
      }

      fQAhist->Fill("Tracks_all",1);

      // Get MC information
      if(hasMC){

        //Printf("--- Particle: %i", iTrack);

        // Set MC features to dummy values
        mcEta           = -99;
        mcPhi           = -99;
        mcPt            = -99;
        // Vert values that will be written to TTree as floats
        mcVert[0]       = {-99};
        mcVert[1]       = {-99};
        mcVert[2]       = {-99};
        // Vert values that are passed to mcTrack as doubles
        Double_t mcVertD[3] = {-99, -99, -99};
        iPdg            = -9999;
        iPdgMother      = -9999;
        HasMother       = kFALSE;
        motherLabel     = -9999;
        isInj = 0; // Default is "not injected" (needed because DPMJET)
        // First mother information
        iPdgFirstMother      = -9999;
        gLabelFirstMother    = -9999;
        gLabelMinFirstMother = -9999;
        gLabelMaxFirstMother = -9999;

        if(fIsAOD){
          mcTrack = dynamic_cast<AliAODMCParticle*>(fMCevent->GetTrack(TMath::Abs(track->GetLabel())));
          if(!mcTrack){
            continue;
          }
        }else{
          mcTrack = dynamic_cast<AliMCParticle*>(fMCevent->GetTrack(TMath::Abs(track->GetLabel())));
          if(!mcTrack){
            continue;
          }
        }

        fQAhist->Fill("Tracks_MCcheck", 1);

        // Get basic MC information
        iPdg  = mcTrack->PdgCode();
        // Skip lead ion (present in some MC productions....)
        if(iPdg == 1000822080){
          continue;
        }
        mcEta = mcTrack->Eta();
        mcPhi = mcTrack->Phi();
        mcPt  = mcTrack->Pt();
        mcTrack->XvYvZv(mcVertD);

        for(Int_t i = 0; i < 3; ++i){
          mcVert[i] = (Float_t)mcVertD[i];
        }

        // Check which generator was used
        isInj = CheckGenerator(track->GetLabel());
        // Flag to indicated physical primary or not
        isPhysPrimary = mcTrack->IsPhysicalPrimary();

        // Get label of mother particle
        // Will return -1 if no mother particle
        Int_t gMotherIndex = mcTrack->GetMother();
        if(!(gMotherIndex < 0)){
          if(fIsAOD){
            motherMCtrack = dynamic_cast<AliAODMCParticle*>((fMCevent->GetTrack(gMotherIndex)));
            // Check for mother particle.
            if(!motherMCtrack){
              continue;
            }
          }else{
            motherMCtrack = dynamic_cast<AliMCParticle*>((fMCevent->GetTrack(gMotherIndex)));
            // Check for mother particle.
            if(!motherMCtrack){
              continue;
            }
          }
          // Fill mother particle information
          HasMother   = kTRUE;
          iPdgMother  = motherMCtrack->PdgCode();
          // Get mother label so tracks can be correctly paired
          motherLabel = TMath::Abs(motherMCtrack->GetLabel());

          // If index is minus then particle has no mother particle
          // Otherwise, begin search for original particle
          Int_t gFirstMotherIndex = motherMCtrack->GetMother();

          if(gFirstMotherIndex != -1){

            // Retreive grandmother particle of trackx
            AliMCParticle* firstMotherTrack = (AliMCParticle*)(fMCevent->GetTrack(gFirstMotherIndex));
            // Scan down decay chain until a negative index is returned
            // I.e. first particle in decay is found
            while(gFirstMotherIndex > 0){

              gLabelFirstMother = gFirstMotherIndex; // Use label as temp. index storage
              firstMotherTrack  = (AliMCParticle*)(fMCevent->GetTrack(gLabelFirstMother));
              gFirstMotherIndex = firstMotherTrack->GetMother();
            }

            // If greatgrand-mother (etc) was found, store pdg code
            // Otherwise, grandmother was already primary
            if(gLabelFirstMother != -1){
              iPdgFirstMother = firstMotherTrack->PdgCode();
            }
            else{
              gLabelFirstMother = gFirstMotherIndex; // set mother to first mother
              iPdgFirstMother   = iPdgMother;
            }

            Int_t nParticles = fMCevent->GetNumberOfTracks();
            // Ignore initial lead ion
            if(iPdgFirstMother == 1000822080){
              continue;
            }
            // Needed for HIJING....
            // find range of -1 - minimum
            gLabelMinFirstMother = gLabelFirstMother;
            while(gFirstMotherIndex < 0){
              gLabelMinFirstMother--;
              if(gLabelMinFirstMother < 0){
                gFirstMotherIndex = 0;
              }
              else{
                firstMotherTrack = (AliMCParticle*)(fMCevent->GetTrack(gLabelMinFirstMother));
                gFirstMotherIndex = firstMotherTrack->GetMother();
              }
            }
            gLabelMinFirstMother ++; // set back by one
            gFirstMotherIndex = -1; // set back to -1

            // find range of -1 - maximum
            gLabelMaxFirstMother = gLabelFirstMother;

            while(gFirstMotherIndex < 0){
              gLabelMaxFirstMother++;
              if(gLabelMaxFirstMother > nParticles){
                gFirstMotherIndex = 0;
              }
              else{
                firstMotherTrack = (AliMCParticle*)(fMCevent->GetTrack(gLabelMaxFirstMother));
                gFirstMotherIndex = firstMotherTrack->GetMother();
              }
            }
            gLabelMaxFirstMother --; // set back by one

                  }
        }
      }// End if(hasMC)

      // Apply global track trackFilter
      if(!fIsAOD){
        if(!(fESDtrackCuts->AcceptTrack(dynamic_cast<const AliESDtrack*>(track)))){
          continue;
        }
      }
      else{
        UInt_t selectedMask = (1<<trackFilter->GetCuts()->GetEntries())-1;
        if( selectedMask != (trackFilter->IsSelected((AliVParticle*)track)) ){
          continue;
        }
      }

      fQAhist->Fill("Tracks_Cuts", 1);

      pt  = track->Pt();
      eta = track->Eta();
      phi = track->Phi();

      // Get PID response values from available detectors
      EnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
      EnSigmaTPCcorr = EnSigmaTPC;
      if(fUseTPCcorr && pt < maxPtPIDcorrection){
        EnSigmaTPCcorr -= AliDielectronPID::GetCntrdCorr(track);
        EnSigmaTPCcorr /= AliDielectronPID::GetWdthCorr(track);
      }
      if(fHasSDD){
        EnSigmaITS = fPIDResponse->NumberOfSigmasITS(track, AliPID::kElectron);
        EnSigmaITScorr = EnSigmaITS;
        if(fUseITScorr && (EnSigmaITS != -999.) && pt < maxPtPIDcorrection){
          EnSigmaITScorr -= AliDielectronPID::GetCntrdCorrITS(track);
          EnSigmaITScorr /= AliDielectronPID::GetWdthCorrITS(track);
        }
      }
      EnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
      EnSigmaTOFcorr = EnSigmaTOF;
      if(fUseTOFcorr && (EnSigmaTOF != -999.) && pt < maxPtPIDcorrection){
        EnSigmaTOFcorr -= AliDielectronPID::GetCntrdCorrTOF(track);
        EnSigmaTOFcorr /= AliDielectronPID::GetWdthCorrTOF(track);
      }

      PnSigmaITS = fPIDResponse->NumberOfSigmasITS(track, AliPID::kPion);
      PnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
      PnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);

      if(storeKaonPID){
        KnSigmaITS = fPIDResponse->NumberOfSigmasITS(track, AliPID::kKaon);
        KnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        KnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
      }
      
      if(storeProtonPID){
        PrnSigmaITS = fPIDResponse->NumberOfSigmasITS(track, AliPID::kProton);
        PrnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        PrnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
      }
      // Get TOF start time mask
      AliTOFPIDResponse TOFresponse = fPIDResponse->GetTOFResponse();
      TOFstartMask = TOFresponse.GetStartTimeMask((Float_t)track->GetP());

      // Get TPC information
      // kNclsTPC
      nTPCclusters = track->GetTPCNcls();

      // kNFclsTPCr
      nTPCcrossed = track->GetTPCClusterInfo(2,1);

      fTPCcrossOverFind = 0;
      nTPCfindable = track->GetTPCNclsF();
      // kNFclsTPCfCross
      if(nTPCfindable > 0){
        fTPCcrossOverFind = nTPCcrossed/nTPCfindable;
      }

      tpcSharedMap = 0;
      if(fIsAOD){
        tpcSharedMap = (dynamic_cast<AliAODTrack*>(track))->GetTPCSharedMap();
      }else{
        tpcSharedMap = (dynamic_cast<AliESDtrack*>(track))->GetTPCSharedMap();
      }

      nTPCshared = -1;
      if(fIsAOD){
        nTPCshared = tpcSharedMap.CountBits(0) - tpcSharedMap.CountBits(159);
      }else{
        nTPCshared = (dynamic_cast<AliESDtrack*>(track))->GetTPCnclsS();
      }

      chi2TPC = track->GetTPCchi2();

      // DCA values
      Float_t  DCAesd[2] = {-999, -999};
      Double_t DCAaod[2] = {-999, -999};
      Double_t DCAcov[3] = {-999., -999., -999.}; // Covariance matrix

      if(!fIsAOD){
        // Arguments: xy, z
        // 'Float_t' needed for 'virtual void AliESDtrack::GetImpactParameters(Float_t p[2], Float_t cov[3]) const'
        Float_t DCAcov1[3] = {-999., -999., -999.}; // Covariance matrix
        track->GetImpactParameters(DCAesd, DCAcov1);
        DCAcov[0] = DCAcov1[0]; DCAcov[1] = DCAcov1[1]; DCAcov[2] = DCAcov1[2];
        // Store DCA values
        DCA[0] = DCAesd[0];
        DCA[1] = DCAesd[1];
      }
      else{
        GetDCA(const_cast<const AliVEvent*>(event), dynamic_cast<const AliAODTrack*>(track), DCAaod, DCAcov);
        // Store DCA values
        DCA[0] = static_cast<Float_t>(DCAaod[0]);
        DCA[1] = static_cast<Float_t>(DCAaod[1]);
      }

      if(fExtraDCA){
        // Compute the normalized DCAs
        // Neglect the resolution mixed term DCAcov[1]
        if(DCAcov[0] > 0.){
          DCAsigma[0] = DCA[0]/TMath::Sqrt(static_cast<Float_t>(DCAcov[0]));
        }
        if(DCAcov[2] > 0.){
          DCAsigma[1] = DCA[1]/TMath::Sqrt(static_cast<Float_t>(DCAcov[2]));
        }

        // Approximation of the Impact Parameter
        // Get TVectors for primary and secondary vertex as well as for particle momentum
        // ImpPar = |p x (priVertex- secVertex)|/|p|

        TVector3 priVtx(primaryVertex[0],primaryVertex[1],primaryVertex[2]); // Primary vertex

        Double_t secvtx[3];
        track->GetXYZ(secvtx);
        TVector3 secVtx(secvtx[0],secvtx[1],secvtx[2]); // Particle position

        Double_t  mompart[3];
        track->GetPxPyPz(mompart);
        TVector3 momPart(mompart[0],mompart[1],mompart[2]); // Particle momentum

        priVtx -= secVtx;

        TVector3 numerator = momPart.Cross(priVtx);
        Double_t denominator = momPart.Mag();

        DCA3[0] = numerator.X()/denominator;
        DCA3[1] = numerator.Y()/denominator;
        DCA3[2] = numerator.Z()/denominator;
      }

      // Get ITS information
      // kNclsITS
      nITS = track->GetNcls(0);
      chi2ITS = track->GetITSchi2();

      fITSshared = 0.;
      for(Int_t d = 0; d < 6; d++){
        fITSshared += static_cast<Double_t>(track->HasSharedPointOnITSLayer(d));
      }
      fITSshared /= nITS;

      if(fIsAOD){
        SPDfirst = (dynamic_cast<AliAODTrack*>(track))->HasPointOnITSLayer(0);
      }else{
        SPDfirst = (dynamic_cast<AliESDtrack*>(track))->HasPointOnITSLayer(0);
      }

      // Get ITS and TPC signals
      ITSsignal = track->GetITSsignal();
      TPCsignal = track->GetTPCsignal();
      TOFsignal = track->GetTOFsignal();

      Int_t fCutMaxChi2TPCConstrainedVsGlobalVertexType = fESDtrackCuts->kVertexTracks | fESDtrackCuts->kVertexSPD;

      const AliVVertex* vertex = 0;

      if(fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDtrackCuts->kVertexTracks){
        vertex = track->GetEvent()->GetPrimaryVertexTracks();
      }

      if((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDtrackCuts->kVertexSPD){
        vertex = track->GetEvent()->GetPrimaryVertexSPD();
      }

      if((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDtrackCuts->kVertexTPC){
        vertex = track->GetEvent()->GetPrimaryVertexTPC();
      }

      // Get golden Chi2
      goldenChi2 = -1;
      if(vertex->GetStatus()){
        if(fIsAOD){
          goldenChi2 = dynamic_cast<AliAODTrack*>(track)->GetChi2TPCConstrainedVsGlobal();
        }
        else{
          goldenChi2 = dynamic_cast<AliESDtrack*>(track)->GetChi2TPCConstrainedVsGlobal(dynamic_cast<const AliESDVertex*>(vertex));
        }
      }

      charge = track->Charge();

      fTree->Fill();

    } // End loop over tracks
  }// End of standard TTree creation
  else if(fIsV0tree){
    for(Int_t iV0 = 0; iV0 < numV0s; iV0++){

      AliESDv0* V0vertex = esdEvent->GetV0(iV0);

      if(!V0vertex){
        AliError(Form("Could not receive V0 track %d", iV0));
        continue;
      }

      //fQAhist->Fill("Tracks_all",1);

      //Get V0 daughter tracks
      AliESDtrack* negTrack = esdEvent->GetTrack(V0vertex->GetIndex(0));
      AliESDtrack* posTrack = esdEvent->GetTrack(V0vertex->GetIndex(1));
      if(!negTrack || !posTrack){
        Printf("Daughter track of v0 not found: %p - %p \n",negTrack, posTrack);
        continue;
      }
      //Check for like-sign V0 candidates
      if(negTrack->Charge() == posTrack->Charge()){ continue; }

      //Apply kinematic and PID cuts to both legs
      if(isV0daughterAccepted(negTrack) != kTRUE){ continue; }
      if(isV0daughterAccepted(posTrack) != kTRUE){ continue; }

      pointingAngle = V0vertex->GetV0CosineOfPointingAngle();
      daughtersDCA  = V0vertex->GetDcaV0Daughters();
      decayLength   = V0vertex->GetRr();
      v0mass        = V0vertex->M();
      //Super loose cuts on V0 topological qualities(stored in Tree to be cut on later)
      if( pointingAngle < 0.8 || daughtersDCA < 0.05 || decayLength < 0.01 ){ continue; }

      //fQAhist->Fill("Arm. cuts",1);

      ptArm = V0vertex->PtArmV0();
      alpha = V0vertex->AlphaV0();
      //Armentors plots is not filled until MC checks completed



      //TODO: Improve efficiency of MC section
      //Currently: checks neg particle, then pos particle.
      //If both pass, write pos particle features, then get neg features
      //again!
      Int_t label = -9999;
      Double_t mcVertD[3] = {-99, -99, -99};
      if(hasMC){
        //-------- Check negative charge particle --------------
        //Only want to check validity of track
        label = negTrack->GetLabel();

        mcTrack = dynamic_cast<AliMCParticle*>(fMCevent->GetTrack(TMath::Abs(label)));
        if(!mcTrack){
          continue;
        }

        Int_t gMotherIndex = mcTrack->GetMother();

        if(!(gMotherIndex < 0)){
          if(fIsAOD){
            motherMCtrack = dynamic_cast<AliAODMCParticle*>((fMCevent->GetTrack(gMotherIndex)));
            //Check for mother particle.
            if(!motherMCtrack){
              continue;
            }
          }else{
            motherMCtrack = dynamic_cast<AliMCParticle*>((fMCevent->GetTrack(gMotherIndex)));
            //Check for mother particle.
            if(!motherMCtrack){
              continue;
            }
          }
        }//End loop over motherIndex
        //--------- End negative particle check  -------------

        //--------- Get postive particle features
        //Set features for writing into TTree

        label = posTrack->GetLabel();

        mcTrack = dynamic_cast<AliMCParticle*>(fMCevent->GetTrack(TMath::Abs(label)));
        if(!mcTrack){
          continue;
        }

        //Fill MC check
        //fQAhist->Fill("Tracks_MCcheck", 1);

        iPdg = mcTrack->PdgCode();

        mcEta = mcTrack->Eta();
        mcPhi = mcTrack->Phi();
        mcPt  = mcTrack->Pt();
        mcTrack->XvYvZv(mcVertD);

        for(Int_t i = 0; i < 3; ++i){
          mcVert[i] = (Float_t)mcVertD[i];
        }

        gMotherIndex = mcTrack->GetMother();

        if(!(gMotherIndex < 0)){
          if(fIsAOD){
            motherMCtrack = dynamic_cast<AliAODMCParticle*>((fMCevent->GetTrack(gMotherIndex)));
            //Check for mother particle.
            if(!motherMCtrack){
              continue;
            }
          }else{
            motherMCtrack = dynamic_cast<AliMCParticle*>((fMCevent->GetTrack(gMotherIndex)));
            //Check for mother particle.
            if(!motherMCtrack){
              continue;
            }
          }
          HasMother   = kTRUE;
          iPdgMother  = motherMCtrack->PdgCode();
          motherLabel = TMath::Abs(motherMCtrack->GetLabel());
        }//End loop over motherIndex
      }//End ifMC section

      fArmPlot->Fill(alpha, ptArm);

      //Get positive particle obsevables
      pt  = posTrack->Pt();
      eta = posTrack->Eta();
      phi = posTrack->Phi();

      EnSigmaITS = -999;
      PnSigmaITS = -999;
      KnSigmaITS = -999;

      if(fHasSDD){
        EnSigmaITS = fPIDResponse->NumberOfSigmasITS(posTrack,(AliPID::EParticleType)AliPID::kElectron);
        PnSigmaITS = fPIDResponse->NumberOfSigmasITS(posTrack,(AliPID::EParticleType)AliPID::kPion);
        KnSigmaITS = fPIDResponse->NumberOfSigmasITS(posTrack,(AliPID::EParticleType)AliPID::kKaon);
      }
      EnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(posTrack,(AliPID::EParticleType)AliPID::kElectron);
      EnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(posTrack,(AliPID::EParticleType)AliPID::kElectron);
      PnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(posTrack,(AliPID::EParticleType)AliPID::kPion);
      PnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(posTrack,(AliPID::EParticleType)AliPID::kPion);
      KnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(posTrack,(AliPID::EParticleType)AliPID::kKaon);
      KnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(posTrack,(AliPID::EParticleType)AliPID::kKaon);


      //DCA values
      //Initialized here incase goto is invoked
      Float_t DCAesd[2]  = {0.0,0.0};
      Double_t DCAaod[2] = {0.0,0.0};
      Double_t DCAcov[2] = {0.0, 0.0};

      if(TMath::Abs(EnSigmaTPC) > fESigTPCMax){
        goto negativeTrack;
      }

      //Get TPC information
      nTPCclusters = posTrack->GetTPCNcls();
      if(nTPCclusters < 70){ continue;}

      nTPCcrossed = posTrack->GetTPCClusterInfo(2,1);
      if(nTPCcrossed < 60){ continue;}

      fTPCcrossOverFind = 0;

      nTPCfindable = posTrack->GetTPCNclsF();
      if(nTPCfindable > 0){
        fTPCcrossOverFind = nTPCcrossed/nTPCfindable;
      }

      if(fTPCcrossOverFind < 0.3 || fTPCcrossOverFind > 1.1){ continue;}

      tpcSharedMap = 0;
      if(fIsAOD){
        tpcSharedMap = (dynamic_cast<AliAODTrack*>(posTrack))->GetTPCSharedMap();
      }else{
        tpcSharedMap = (dynamic_cast<AliESDtrack*>(posTrack))->GetTPCSharedMap();
      }

      nTPCshared = -1;
      if(fIsAOD){
        nTPCshared = tpcSharedMap.CountBits(0) - tpcSharedMap.CountBits(159);
      }else{
        nTPCshared = (dynamic_cast<AliESDtrack*>(posTrack))->GetTPCnclsS();
      }

      chi2TPC = posTrack->GetTPCchi2(); //Function only implemented in ESDs. Returns dumym value for AODs

      //Check for refits
      if((posTrack->GetStatus() & AliVTrack::kITSrefit) <= 0){ continue;}
      if((posTrack->GetStatus() & AliVTrack::kTPCrefit) <= 0){ continue;}

      if(!fIsAOD){
        posTrack->GetImpactParameters( &DCAesd[0], &DCAesd[1]);
      }
      else{
        GetDCA(const_cast<const AliVEvent*>(event), dynamic_cast<const AliAODTrack*>(posTrack), DCAaod, DCAcov);
      }
      //Final DCA values stored here
      if(!fIsAOD){
        DCA[0] = static_cast<Double_t>(DCAesd[0]);
        DCA[1] = static_cast<Double_t>(DCAesd[1]);
      }
      else{
        DCA[0] = static_cast<Double_t>(DCAaod[0]);
        DCA[1] = static_cast<Double_t>(DCAaod[1]);
      }

      //ITS clusters and shared clusters
      nITS = posTrack->GetNumberOfITSClusters();
      fITSshared = 0.;
      for(Int_t d = 0; d < 6; d++){
        fITSshared += static_cast<Double_t>(posTrack->HasSharedPointOnITSLayer(d));
      }
      fITSshared /= nITS;

      SPDfirst = (dynamic_cast<AliESDtrack*>(posTrack))->HasPointOnITSLayer(0);

      charge = posTrack->Charge();

      //Fill with feature from positive track
      fTree->Fill();

      //Skip to this point if positive track failed PID cut
      negativeTrack:

      //Get negative particle obsevables
      pt  = negTrack->Pt();
      eta = negTrack->Eta();
      phi = negTrack->Phi();

      if(fHasSDD){
        EnSigmaITS = fPIDResponse->NumberOfSigmasITS(negTrack,(AliPID::EParticleType)AliPID::kElectron);
        PnSigmaITS = fPIDResponse->NumberOfSigmasITS(negTrack,(AliPID::EParticleType)AliPID::kPion);
        KnSigmaITS = fPIDResponse->NumberOfSigmasITS(negTrack,(AliPID::EParticleType)AliPID::kKaon);
      }
      EnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(negTrack,(AliPID::EParticleType)AliPID::kElectron);
      EnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(negTrack,(AliPID::EParticleType)AliPID::kElectron);
      PnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(negTrack,(AliPID::EParticleType)AliPID::kPion);
      PnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(negTrack,(AliPID::EParticleType)AliPID::kPion);
      KnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(negTrack,(AliPID::EParticleType)AliPID::kKaon);
      KnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(negTrack,(AliPID::EParticleType)AliPID::kKaon);

      if(TMath::Abs(EnSigmaTPC) > fESigTPCMax){
        continue;
      }

      //Get TPC information
      nTPCclusters = negTrack->GetTPCNcls();
      if(nTPCclusters < 70){ continue;}

      nTPCcrossed = negTrack->GetTPCClusterInfo(2,1);
      if(nTPCcrossed < 60){ continue;}

      fTPCcrossOverFind = 0;

      nTPCfindable = negTrack->GetTPCNclsF();
      if(nTPCfindable > 0){
        fTPCcrossOverFind = nTPCcrossed/nTPCfindable;
      }

      if(fTPCcrossOverFind < 0.3 || fTPCcrossOverFind > 1.1){ continue;}

      tpcSharedMap = 0;
      if(fIsAOD){
        tpcSharedMap = (dynamic_cast<AliAODTrack*>(negTrack))->GetTPCSharedMap();
      }else{
        tpcSharedMap = (dynamic_cast<AliESDtrack*>(negTrack))->GetTPCSharedMap();
      }

      nTPCshared = -1;
      if(fIsAOD){
        nTPCshared = tpcSharedMap.CountBits(0) - tpcSharedMap.CountBits(159);
      }else{
        nTPCshared = (dynamic_cast<AliESDtrack*>(negTrack))->GetTPCnclsS();
      }

      chi2TPC = negTrack->GetTPCchi2(); //Function only implemented in ESDs. Returns dumym value for AODs

      //Check for refits
      if((negTrack->GetStatus() & AliVTrack::kITSrefit) <= 0){ continue;}
      if((negTrack->GetStatus() & AliVTrack::kTPCrefit) <= 0){ continue;}

      //DCA values
      if(!fIsAOD){
        negTrack->GetImpactParameters( &DCAesd[0], &DCAesd[1]);
      }
      else{
        GetDCA(const_cast<const AliVEvent*>(event), dynamic_cast<const AliAODTrack*>(negTrack), DCAaod, DCAcov);
      }
      //Final DCA values stored here
      if(!fIsAOD){
        DCA[0] = static_cast<Double_t>(DCAesd[0]);
        DCA[1] = static_cast<Double_t>(DCAesd[1]);
      }
      else{
        DCA[0] = static_cast<Double_t>(DCAaod[0]);
        DCA[1] = static_cast<Double_t>(DCAaod[1]);
      }

      //ITS clusters and shared clusters
      nITS = negTrack->GetNumberOfITSClusters();
      fITSshared = 0.;
      for(Int_t d = 0; d < 6; d++){
        fITSshared += static_cast<Double_t>(negTrack->HasSharedPointOnITSLayer(d));
      }
      fITSshared /= nITS;

      SPDfirst = (dynamic_cast<AliESDtrack*>(negTrack))->HasPointOnITSLayer(0);

      charge = negTrack->Charge();
      //Write negative observales to tree (v0 information written twice. Filter by looking at only pos or neg charge)
      if(hasMC){
        label = negTrack->GetLabel();

        mcTrack = dynamic_cast<AliMCParticle*>(fMCevent->GetTrack(TMath::Abs(label)));
        //Redundant?
        if(!mcTrack){
          continue;
        }

        //Fill MC check
        //fQAhist->Fill("Tracks_MCcheck", 1);

        iPdg = mcTrack->PdgCode();

        mcEta = mcTrack->Eta();
        mcPhi = mcTrack->Phi();
        mcPt  = mcTrack->Pt();
        mcTrack->XvYvZv(mcVertD);

        for(Int_t i = 0; i < 3; ++i){
          mcVert[i] = (Float_t)mcVertD[i];
        }
        Int_t gMotherIndex = mcTrack->GetMother();

        if(!(gMotherIndex < 0)){
          if(fIsAOD){
            motherMCtrack = dynamic_cast<AliAODMCParticle*>((fMCevent->GetTrack(gMotherIndex)));
            //Check for mother particle.
            if(!motherMCtrack){
              continue;
            }
          }else{
            motherMCtrack = dynamic_cast<AliMCParticle*>((fMCevent->GetTrack(gMotherIndex)));
            //Check for mother particle.
            if(!motherMCtrack){
              continue;
            }
          }
          HasMother = kTRUE;
          iPdgMother = motherMCtrack->PdgCode();
          motherLabel = TMath::Abs(motherMCtrack->GetLabel());
        }//End loop over motherIndex
      }


      //Fill with features from negative track
      fTree->Fill();
    }//End loop over v0's for this event
  }//End V0 code

}

//~ //________________________________________________________________________

void  AliAnalysisTaskSimpleTreeMaker::FinishTaskOutput(){
    // Finish task output

    // not implemented ...

}
//~

//~ //________________________________________________________________________

void AliAnalysisTaskSimpleTreeMaker::Terminate(Option_t *){
    // Draw result to the screen

    // Called once at the end of the query

    // not implemented ...

}
//~


//________________________________________________________________________

Int_t AliAnalysisTaskSimpleTreeMaker::IsEventAccepted(AliVEvent *event){

    if(!fIsV0tree){
      UInt_t selectedMask = (1<<eventFilter->GetCuts()->GetEntries())-1;
      varManager->SetEvent(event);
      if(selectedMask == (eventFilter->IsSelected(event))){
        return 1;
      }
    }
    else{
      if(TMath::Abs(event->GetPrimaryVertexSPD()->GetZ()) < 10){
        return 1;
      }
    }
    return 0;
}

Bool_t AliAnalysisTaskSimpleTreeMaker::isV0daughterAccepted(AliVTrack* track){

    Bool_t answer = kFALSE;
    //Kinematic cuts
    Float_t pt   = (Float_t)track->Pt();
    if( pt < fPtMin || pt > fPtMax ){ return answer; }
    Float_t eta  = (Float_t)track->Eta();
    if( eta < fEtaMin || eta > fEtaMax ){ return answer; }

    if(!fIsAOD){
        if(!(fESDtrackCuts->AcceptTrack(dynamic_cast<AliESDtrack*>(track)))){
            return answer;
        }
    }
    else{
        if(!(dynamic_cast<AliAODTrack*>(track))->TestFilterBit(fFilterBit)){
            return answer;
        }
    }
    //fQAhist->Fill("Tracks_KineCuts", 1);
    //Do not apply PID cuts
    //fQAhist->Fill("Tracks_PIDcuts",1);

    answer = kTRUE;
    return answer;
}

//
void AliAnalysisTaskSimpleTreeMaker::SetupTrackCuts(AliDielectronCutGroup* cuts){

  //Initialise track cut object and add to DielectronCutGroup
  //Track cuts include kinematic cuts, track cuts and PID cuts
  trackFilter = new AliAnalysisFilter("TrackFilter", "trackCuts");
  trackFilter->AddCuts(cuts);
}

void AliAnalysisTaskSimpleTreeMaker::SetupEventCuts(AliDielectronEventCuts* cuts){

  //Initiliase the event trackFilter object and add event cuts
  eventFilter = new AliAnalysisFilter("eventFilter", "eventFilter");
  eventFilter->AddCuts(cuts);
}

// Check if the generator is on the list of generators
// If found, assign track with integer value correspding to generator
//0 = gen purp, 1=Pythia CC_1, 2= Pythia BB_1, 3=Pythia B_1, 4=Jpsi2ee_1, 5=B2Jpsi2ee_1";
Int_t AliAnalysisTaskSimpleTreeMaker::CheckGenerator(Int_t trackID){

  if(fGeneratorHashes.size() == 0){
    return -1;
  }

  TString genname;
  Bool_t hasGenerator = fMCevent->GetCocktailGenerator(TMath::Abs(trackID), genname);
  if(!hasGenerator){
    Printf("no cocktail header list was found for this track");
    return -1;
  }
  else{
    for(UInt_t i = 0; i < fGeneratorHashes.size(); ++i){
      if(genname.Hash() == fGeneratorHashes[i]){
        return i+1;
      }
    }
    return 0;
  }
  return -1; // should not happen
}
