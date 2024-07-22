#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom2.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliMCParticle.h"
#include "AliVParticle.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliStack.h"
#include "TPDGCode.h"
#include "AliEventCuts.h"
#include "AliAnalysisTaskHypTritTestAOD.h"


using namespace std;
/// \cond CLASSIMP
ClassImp(AliAnalysisTaskHypTritTestAOD)
/// \endcond

// Default Constructor
AliAnalysisTaskHypTritTestAOD::AliAnalysisTaskHypTritTestAOD()
  :AliAnalysisTaskSE("AliAnalysisTaskHypTritTestAOD"),
  fInputHandler(0),
  fPID(0),
  fESDevent(0),
  fHistdEdx(0),
  fHistdEdxV0(0),
  fHistNumEvents(0),
  fHistTrigger(0),
  fHistV0(0),
  fTree(0),
  fTreeMCGen(0),
  fHistogramList(NULL),
  fPrimaryVertex(),
  fSecondaryVertex(),
  tMagField(),
  fMCtrue(0),
  fEventCuts(),
  fPeriod(0),
  fTriggerMask(0),
  fBetheSplines(kFALSE),
  fBetheParamsHe(),
  fBetheParamsT(),
  fUseExternalSplines(kFALSE),
  tM(0),
  tCt(0),
  tPt(0), 
  tP(0), 
  tPVx(0),
  tPVy(0),
  tPVz(0),
  tSVx(0),
  tSVy(0),
  tSVz(0),
  tPx(0),
  tPy(0),
  tPz(0),
  tEta(0),
  tPhi(0),
  tHePx(0),
  tHePy(0),
  tHePz(0),
  tPiPx(0),
  tPiPy(0),
  tPiPz(0),
  tDca(0), 
  tCosPA(0), 
  tY(0), 
  tPiDca(0), 
  tHeDca(0), 
  tHeDedxSigma(0),
  tPiDedxSigma(0),
  tTDedxSigma(0), 
  tHeTOFmass(-999), 
  tPiTOFmass(-999),   
  tHeP(0), 
  tPiP(0), 
  tHeDedx(0), 
  tPiDedx(0), 
  tTrig(0), 
  tHePt(0), 
  tPiPt(0), 
  tImpPiHe(0),
  tChi2(0),
  tNDF(0),
  tZ(0), 
  tMc(0), 
  tRunnumber(0),
  tMultV0M(-999),  
  tMultOfV0M(-999),      
  tMultSPDTracklet(-999),  
  tMultSPDCluster(-999),  
  tMultRef05(-999),      
  tMultRef08(-999),
  tKinkHe(-999),
  tTPCrefitHe(-999),
  tITSrefitHe(-999),  
  tnTPCclusterHe(-999),
  tnITSclusterHe(-999),
  tTPCchi2He(-999),  
  tKinkPi(-999),
  tTPCrefitPi(-999),
  tITSrefitPi(-999),  
  tnTPCclusterPi(-999),
  tnITSclusterPi(-999),
  tTPCchi2Pi(-999),
  tTRDvalid(0),
  tTRDtrigHNU(0),
  tTRDtrigHQU(0),
  tTRDPid(0),
  tTRDnTracklets(0),
  tTRDPt(0),
  tTRDLayerMask(0),
  tTRDSagitta(-1),
  tTRDStack(0),
  tTRDSector(0),
  tTRDPID0(0),
  tTRDPID1(0),
  tTRDPID2(0),
  tTRDPID3(0),
  tTRDPID4(0),
  tTRDPID5(0),
  fYear(0) {

  }

// Constructor
AliAnalysisTaskHypTritTestAOD::AliAnalysisTaskHypTritTestAOD(const char *name)
  :AliAnalysisTaskSE(name),
  fInputHandler(0),
  fPID(0),
  fESDevent(0),
  fHistdEdx(0),
  fHistdEdxV0(0),
  fHistNumEvents(0),
  fHistTrigger(0),
  fHistV0(0),
  fTree(0),
  fTreeMCGen(0),
  fHistogramList(NULL),
  fPrimaryVertex(),
  fSecondaryVertex(),
  tMagField(),
  fMCtrue(0),
  fEventCuts(),
  fPeriod(0),
  fTriggerMask(0),
  fBetheSplines(kFALSE),
  fBetheParamsHe(),
  fBetheParamsT(),
  fUseExternalSplines(kFALSE),
  tM(0),
  tCt(0),
  tPt(0), 
  tP(0),
  tPVx(0),
  tPVy(0),
  tPVz(0),
  tSVx(0),
  tSVy(0),
  tSVz(0),
  tPx(0),
  tPy(0),
  tPz(0),
  tEta(0),
  tPhi(0),  
  tHePx(0),
  tHePy(0),
  tHePz(0),
  tPiPx(0),
  tPiPy(0),
  tPiPz(0),
  tDca(0), 
  tCosPA(0), 
  tY(0), 
  tPiDca(0),  
  tHeDca(0), 
  tHeDedxSigma(0),
   tPiDedxSigma(0),
  tTDedxSigma(0), 
  tHeTOFmass(-999), 
  tPiTOFmass(-999),    
  tHeP(0), 
  tPiP(0), 
  tHeDedx(0), 
  tPiDedx(0), 
  tTrig(0), 
  tHePt(0), 
  tPiPt(0), 
  tImpPiHe(0),
  tChi2(0),
  tNDF(0),
  tZ(0), 
  tMc(0), 
  tRunnumber(0),
  tMultV0M(-999),  
  tMultOfV0M(-999),      
  tMultSPDTracklet(-999),  
  tMultSPDCluster(-999),  
  tMultRef05(-999),      
  tMultRef08(-999),
  tKinkHe(-999),
  tTPCrefitHe(-999),
  tITSrefitHe(-999),  
  tnTPCclusterHe(-999),
  tnITSclusterHe(-999),
  tTPCchi2He(-999),  
  tKinkPi(-999),
  tTPCrefitPi(-999),
  tITSrefitPi(-999),  
  tnTPCclusterPi(-999),
  tnITSclusterPi(-999),
  tTPCchi2Pi(-999),
  tTRDvalid(0),
  tTRDtrigHNU(0),
  tTRDtrigHQU(0),
  tTRDPid(0),
  tTRDnTracklets(0),
  tTRDPt(0),
  tTRDLayerMask(0),
  tTRDSagitta(-1),
  tTRDStack(0),
  tTRDSector(0),
  tTRDPID0(0),
  tTRDPID1(0),
  tTRDPID2(0),
  tTRDPID3(0),
  tTRDPID4(0),
  tTRDPID5(0),
  fYear(0)
  {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
    DefineOutput(3, TTree::Class());
  }

// Destructor
AliAnalysisTaskHypTritTestAOD::~AliAnalysisTaskHypTritTestAOD() {
}

void AliAnalysisTaskHypTritTestAOD::UserCreateOutputObjects() {
  fInputHandler = dynamic_cast<AliInputEventHandler*>
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!fInputHandler) {
    AliError("Could not get ESD InputHandler.\n");
    return;
  }
  fPID = fInputHandler->GetPIDResponse();
  if (!fPID) {
    AliError("Could not get PID response.\n");
    return;
  }
  fHistdEdx = new TH2F("fHistdEdX","dE/dx;#frac{#it{p}}{z} (GeV/#it{c});TPC Signal (a.u.)",1000,-5.0,5.0,1000,0.0,1500);
  fHistdEdxV0 = new TH2F("fHistdEdXV0","dE/dx;#frac{#it{p}}{z} (GeV/#it{c});TPC Signal (a.u.)",1000,-5.0,5.0,1000,0.0,1500);

  fHistNumEvents = new TH1F("fHistNumEvents","Number of Events",2,0,2);
  fHistNumEvents->GetXaxis()->SetBinLabel(1,"before PhysSel");
  fHistNumEvents->GetXaxis()->SetBinLabel(2,"after PhysSel");

  fHistTrigger = new TH1F("fHistTrigger","Trigger",8,0,8);
  fHistTrigger->GetXaxis()->SetBinLabel(1,"other");
  fHistTrigger->GetXaxis()->SetBinLabel(2,"kINT7");
  fHistTrigger->GetXaxis()->SetBinLabel(3,"kHighMultV0");
  fHistTrigger->GetXaxis()->SetBinLabel(4,"kHighMultSPD");
  fHistTrigger->GetXaxis()->SetBinLabel(5,"HNU");
  fHistTrigger->GetXaxis()->SetBinLabel(6,"HQU");
  fHistTrigger->GetXaxis()->SetBinLabel(7,"HJT");
  fHistTrigger->GetXaxis()->SetBinLabel(8,"HSE");
  fHistV0 = new TH1F("fHistV0","Trigger V0s",8,0,8);
  fHistV0->GetXaxis()->SetBinLabel(1,"other");
  fHistV0->GetXaxis()->SetBinLabel(2,"kINT7");
  fHistV0->GetXaxis()->SetBinLabel(3,"kHighMultV0");
  fHistV0->GetXaxis()->SetBinLabel(4,"kHighMultSPD");
  fHistV0->GetXaxis()->SetBinLabel(5,"HNU");
  fHistV0->GetXaxis()->SetBinLabel(6,"HQU");
  fHistV0->GetXaxis()->SetBinLabel(7,"HJT");
  fHistV0->GetXaxis()->SetBinLabel(8,"HSE");

  fHistogramList = new TList();
  fHistogramList->SetOwner(kTRUE);
  fHistogramList->SetName(GetName());
  fHistogramList->Add(fHistdEdx);
  fHistogramList->Add(fHistdEdxV0);
  fHistogramList->Add(fHistNumEvents);
  fHistogramList->Add(fHistTrigger);
  fHistogramList->Add(fHistV0);


  fEventCuts.AddQAplotsToList(fHistogramList);

  fTree = new TTree("treeKF","fTree");
  fTree->Branch("ct",               &tCt,                  "ct/D");
  fTree->Branch("m",                &tM,                   "m/D");
  fTree->Branch("pt",               &tPt,                  "pt/D");
  fTree->Branch("tPx",              &tPx,                  "tPx/D");
  fTree->Branch("tPy",              &tPy,                  "tPy/D");
  fTree->Branch("tPz",              &tPz,                  "tPz/D");  
  fTree->Branch("tHePx",            &tHePx,                "tHePx/D");
  fTree->Branch("tHePy",            &tHePy,                "tHePy/D");
  fTree->Branch("tHePz",            &tHePz,                "tHePz/D");
  fTree->Branch("tPiPx",            &tPiPx,                "tPiPx/D");    
  fTree->Branch("tPiPy",            &tPiPy,                "tPiPy/D");
  fTree->Branch("tPiPz",            &tPiPz,                "tPiPz/D");
  fTree->Branch("p",                &tP,                   "p/D");
  fTree->Branch("hePt",             &tHePt,                "hePt/D");
  fTree->Branch("piPt",             &tPiPt,                "piPt/D");
  fTree->Branch("dca",              &tDca,                 "dca/D");
  fTree->Branch("cosPA",            &tCosPA,               "cosPA/D");
  fTree->Branch("y",                &tY,                   "y/D");
  fTree->Branch("piDca",            &tPiDca,               "piDca/D");
  fTree->Branch("heNcls",           &tnTPCclusterHe,       "heNcls/D");
  fTree->Branch("piNcls",           &tnTPCclusterPi,       "piNcls/D");
  fTree->Branch("heItsNcls",        &tnITSclusterHe,       "heItsNcls/D");
  fTree->Branch("piItsNcls",        &tnITSclusterPi,       "piItsNcls/D");
  fTree->Branch("heTPCrefit",       &tTPCrefitHe,          "heTPCrefit/D");
  fTree->Branch("piTPCrefit",       &tTPCrefitPi,          "piTPCrefit/D");
  fTree->Branch("heITSrefit",       &tITSrefitHe,          "heITSrefit/D");
  fTree->Branch("piITSrefit",       &tITSrefitPi,          "piITSrefit/D");
  fTree->Branch("heKink",           &tKinkHe,              "heKink/D");
  fTree->Branch("piKink",           &tKinkPi,              "piKink/D");
  fTree->Branch("heDca",            &tHeDca,               "heDca/D");
  fTree->Branch("impPiHe",          &tImpPiHe,             "impPiHe/D");
  fTree->Branch("heDedxSigma",      &tHeDedxSigma,         "heDedxSigma/D");
  fTree->Branch("tDedxSigma",       &tTDedxSigma,          "ptDedxSigma/D");
  fTree->Branch("tPiDedxSigma",     &tPiDedxSigma,          "tPiDedxSigma/D");
  fTree->Branch("heTOFmass",        &tHeTOFmass,           "heTOFmass/D");
  fTree->Branch("piTOFmass",        &tPiTOFmass,           "piTOFmass/D");
  fTree->Branch("heP",              &tHeP,                 "heP/D");
  fTree->Branch("piP",              &tPiP,                 "piP/D");
  fTree->Branch("heDedx",           &tHeDedx,              "heDedx/D");
  fTree->Branch("piDedx",           &tPiDedx,              "piDedx/D");
  fTree->Branch("trig",             &tTrig,                "trig/D");
  fTree->Branch("chi2",             &tChi2,                "chi2/D");
  fTree->Branch("NDF",              &tNDF,                 "NDF/D");
  /*fTree->Branch("tTRDvalid",        &tTRDvalid,            "tTRDvalid/I");
  fTree->Branch("tTRDtrigHNU",      &tTRDtrigHNU,          "tTRDtrigHNU/I");
  fTree->Branch("tTRDtrigHQU",      &tTRDtrigHQU,          "tTRDtrigHQU/I");
  fTree->Branch("tTRDPid",          &tTRDPid,              "tTRDPid/I");
  fTree->Branch("tTRDnTracklets",   &tTRDnTracklets,       "tTRDnTracklets/I");
  fTree->Branch("tTRDPt",           &tTRDPt,               "tTRDPt/I");
  fTree->Branch("tTRDLayerMask",    &tTRDLayerMask,        "tTRDLayerMask/I");
  fTree->Branch("tTRDSagitta",      &tTRDSagitta,          "tTRDSagitta/D");
  fTree->Branch("tTRDStack",        &tTRDStack,            "tTRDStack/I");
  fTree->Branch("tTRDSector",       &tTRDSector,           "tTRDSector/I");
  fTree->Branch("tTRDPID0",         &tTRDPID0,             "tTRDPID0/i");
  fTree->Branch("tTRDPID1",         &tTRDPID1,             "tTRDPID1/i");
  fTree->Branch("tTRDPID2",         &tTRDPID2,             "tTRDPID2/i");
  fTree->Branch("tTRDPID3",         &tTRDPID3,             "tTRDPID3/i");
  fTree->Branch("tTRDPID4",         &tTRDPID4,             "tTRDPID4/i");
  fTree->Branch("tTRDPID5",         &tTRDPID5,             "tTRDPID5/i"); */
  fTree->Branch("tPVx",             &tPVx,                 "tPVx/D");
  fTree->Branch("tPVy",             &tPVy,                 "tPVy/D");
  fTree->Branch("tPVz",             &tPVz,                 "tPVz/D");
  fTree->Branch("tSVx",             &tSVx,                 "tSVx/D");
  fTree->Branch("tSVy",             &tSVy,                 "tSVy/D");
  fTree->Branch("tSVz",             &tSVz,                 "tSVz/D");
  fTree->Branch("magField",         &tMagField,            "magField/D");
  fTree->Branch("tMultV0M",         &tMultV0M,             "tMultV0M/D");
  fTree->Branch("tMultOfV0M",       &tMultOfV0M,           "tMultOfV0M/D");
  fTree->Branch("tMultSPDTracklet", &tMultSPDTracklet,     "tMultSPDTracklet/D");
  fTree->Branch("tMultSPDCluster",  &tMultSPDCluster,      "tMultSPDCluster/D");
  fTree->Branch("tMultRef05",       &tMultRef05,           "tMultRef05/D");
  fTree->Branch("tMultRef08",       &tMultRef08,           "tMultRef08/D");
  fTree->Branch("runnumber",        &tRunnumber,           "runnumber/I");
  fTree->Branch("z",                &tZ,                   "z/I");
  fTree->Branch("mc",               &tMc,                  "mc/D");

  fTreeMCGen = new TTree("treeKF_mc", "fTreeMCGen");
  fTreeMCGen->Branch("m",           &tM,                   "m/D");
  fTreeMCGen->Branch("ct",          &tCt,                  "ct/D");
  fTreeMCGen->Branch("pt",          &tPt,                  "pt/D");
  fTreeMCGen->Branch("p",           &tP,                   "p/D");
  fTreeMCGen->Branch("eta",         &tEta,                 "eta/D");
  fTreeMCGen->Branch("phi",         &tPhi,                 "phi/D");
  fTreeMCGen->Branch("y",           &tY,                   "y/D");
  fTreeMCGen->Branch("trig",        &tTrig,                "trig/D");
  fTreeMCGen->Branch("runnumber",   &tRunnumber,           "runnumber/I");
  fTreeMCGen->Branch("z",           &tZ,                   "z/I");
  
  PostData(1, fHistogramList);
  PostData(2, fTree);
  PostData(3, fTreeMCGen);
}

void AliAnalysisTaskHypTritTestAOD::UserExec(Option_t *) {

   // Data
  fESDevent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fESDevent) {
    AliError("Could not get ESD Event.\n");
    return;
  }
  if (!fPID) {
    AliError("Could not get PID response.\n");
    return;
  }
  
  // MC 
  fMCtrue = kTRUE;
  /*AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*>
      (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!mcEventHandler) fMCtrue = kFALSE;
  
  AliMCEvent* mcEvent = 0x0;
  if (mcEventHandler) mcEvent = mcEventHandler->MCEvent();
  if (!mcEvent && fMCtrue) return;*/
  AliMCEvent* mcEvent = 0x0;
  AliAODMCHeader* header = 0x0;
  TClonesArray* MCTrackArray = 0x0;
  if (fMCtrue) {
    //OOB pileup
    header = static_cast<AliAODMCHeader *>(fESDevent->FindListObject("mcHeader"));//AliAODMCHeader::StdBranchName()
    if (!header)
    {
      //AliWarning("No header found.");
      //PostAllData();
      fMCtrue = kFALSE;
      //return;
    }
    MCTrackArray = dynamic_cast<TClonesArray *>(fESDevent->FindListObject("mcparticles"));//AliAODMCParticle::StdBranchName()
    if (!MCTrackArray)
    {
      //AliWarning("No MC track array found.");
      //PostAllData();
      fMCtrue = kFALSE;
      //return;
    }
  }  
  if(fMCtrue) {
    mcEvent = MCEvent();
    if(!fMCEvent){ cout<<"failed"<<endl; return; }
  }

  tRunnumber = fESDevent->GetRunNumber();
 
  fHistNumEvents->Fill(0);  // Count events before event cuts
  fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kAny);
  if ((tRunnumber == 297219 || tRunnumber == 297194 || tRunnumber == 297029
		   || tRunnumber == 296890 || tRunnumber == 296849 || tRunnumber == 296750
		   || tRunnumber == 296749 || tRunnumber == 297481)) fEventCuts.UseTimeRangeCut();
  if(!fEventCuts.AcceptEvent(fESDevent)) {
    PostData(1,fHistogramList);
    return;
  }
  fHistNumEvents->Fill(1);  // Count events after event cuts

  AliMultSelection *multSelection = dynamic_cast<AliMultSelection*> (fESDevent->FindListObject("MultSelection"));
  if (multSelection) {
    tMultV0M = multSelection->GetMultiplicityPercentile("V0M");
    tMultOfV0M = multSelection->GetMultiplicityPercentile("OnlineV0M");
    tMultSPDTracklet = multSelection->GetMultiplicityPercentile("SPDClusters");
    tMultSPDCluster = multSelection->GetMultiplicityPercentile("SPDTracklets");
    tMultRef05 = multSelection->GetMultiplicityPercentile("RefMult05");
    tMultRef08 = multSelection->GetMultiplicityPercentile("RefMult08");
  }

  //AliESDtrackCuts trackCuts("AlitrackCuts", "AlitrackCuts");
  //trackCuts.SetEtaRange(-0.9,0.9);
  //trackCuts.SetAcceptKinkDaughters(kTRUE);
  //trackCuts.SetRequireTPCRefit(kFALSE);
  //trackCuts.SetMaxChi2PerClusterTPC(8);
  //trackCuts.SetMinNClustersTPC(40);   
  
  if (!fUseExternalSplines) SetBetheBlochParams(tRunnumber);
  
  fHistV0->Fill(tTrig);
  tMagField  = fESDevent->GetMagneticField();
  KFParticle3LHJ::SetField(tMagField);
  const AliAODVertex *vertex = fESDevent->GetPrimaryVertexTracks(); 
  fPrimaryVertex.SetXYZ(vertex->GetX(),vertex->GetY(),vertex->GetZ());
  tPVx = vertex->GetX();
  tPVy = vertex->GetY();
  tPVz = vertex->GetZ();

  if (fMCtrue) {
    MCStackLoop(mcEvent);
  }
		
  // Track loop - identify 3He/pi tracks and store them in std::vector
  std::vector<int> heTracks, piTracks;
  const Int_t nTracks = fESDevent->GetNumberOfTracks();
  
  for (Int_t itrack = 0; itrack < nTracks; itrack++) {
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fESDevent->GetTrack(itrack));
    if (!track) continue;
    if (TMath::Abs(track->Eta()) > 0.9) continue;
    if (track->GetTPCNcls() < 40) continue;
    if ((track->GetTPCchi2() / (float)track->GetTPCNcls())>8.) continue;
  
    fHistdEdxV0->Fill(track->GetTPCmomentum() * track->GetSign(), track->GetTPCsignal());

    if (track->GetTPCsignal() > 1500 || track->GetTPCsignal() > 1500) continue;
    if (track->GetTPCmomentum() > 5 || track->GetTPCmomentum() > 5) continue;
    
    if (TMath::Abs(fPID->NumberOfSigmasTPC(track, AliPID::kPion)) < 4) {
      piTracks.push_back(itrack);
    }
    if (TMath::Abs(Bethe(*track, AliPID::ParticleMass(AliPID::kHe3),  2, fBetheParamsHe)) < 5) {
      heTracks.push_back(itrack);
    }
    /*Int_t labelT = TMath::Abs(track->GetLabel());         
    AliAODMCParticle *daughterT = dynamic_cast<AliAODMCParticle*>(mcEvent->GetTrack(labelT));
    if (TMath::Abs(daughterT->GetPdgCode()) == 211) {
      piTracks.push_back(itrack);
    }
    if (TMath::Abs(daughterT->GetPdgCode()) == 1000020030) {
      heTracks.push_back(itrack);
      }*/
  } // end track loop

  KFVertex primKFVertex = CreateKFVertex(*vertex);
  for (Int_t idxHe : heTracks) {
    AliAODTrack* trackHe = dynamic_cast<AliAODTrack*>(fESDevent->GetTrack(idxHe));

    if(fMCtrue){
      Int_t labelHe = TMath::Abs(trackHe->GetLabel());
      Int_t labelMotherHe = TMath::Abs(mcEvent->GetLabelOfParticleMother(labelHe));
      AliAODMCParticle *particleMotherHe = dynamic_cast<AliAODMCParticle*> (mcEvent->GetTrack(labelMotherHe));     
      if(TMath::Abs(particleMotherHe->GetPdgCode()) != 1010010030) continue;
    }
    
    for (Int_t idxPi : piTracks) {

      AliAODTrack* trackPi = dynamic_cast<AliAODTrack*>(fESDevent->GetTrack(idxPi));
      if (trackHe->GetSign() == trackPi->GetSign()) continue;
      
      tZ = 0;
      if (trackHe->GetSign() == -1 && trackPi->GetSign() ==  1) tZ = -1;
      if (trackHe->GetSign() == 1  && trackPi->GetSign() == -1) tZ =  1;
      if (!tZ) continue;

      tMc = 0.;
      if (fMCtrue) {
        
        Int_t labelPi = TMath::Abs(trackPi->GetLabel());
	Int_t labelHe = TMath::Abs(trackHe->GetLabel());      

	//AliAODMCParticle *daughterHe = dynamic_cast<AliAODMCParticle*>(mcEvent->GetTrack(labelHe));
	//AliAODMCParticle *daughterPi =  dynamic_cast<AliAODMCParticle*>(mcEvent->GetTrack(labelPi));
	//if(TMath::Abs(daughterPi->GetPdgCode()) == 211 && TMath::Abs(daughterHe->GetPdgCode()) == 1000020030){
	//if ((daughterPi->GetPdgCode() == 211 && daughterHe->GetPdgCode() == -1000020030) || (daughterHe->GetPdgCode() == 1000020030 && daughterPi->GetPdgCode() == -211)) {
      
	Int_t labelMotherPi = TMath::Abs(mcEvent->GetLabelOfParticleMother(labelPi));
	Int_t labelMotherHe = TMath::Abs(mcEvent->GetLabelOfParticleMother(labelHe));
	       
	AliAODMCParticle *particleMotherPi = dynamic_cast<AliAODMCParticle*> (mcEvent->GetTrack(labelMotherPi));

	if(TMath::Abs(particleMotherPi->GetPdgCode()) != 1010010030) continue;
	
	if (labelMotherHe != labelMotherPi) continue;
	else tMc = 1.;
	//}
	//}
      }
      if(fMCtrue && tMc == 0) continue;
      
      AliExternalTrackParam exTrackHe;
      exTrackHe.CopyFromVTrack(trackHe);
      AliExternalTrackParam exTrackPi;        
      exTrackPi.CopyFromVTrack(trackPi);
      
      KFParticle3LHJ KFHelium = CreateKFParticle(exTrackHe, 2.80839f, 2);
      KFParticle3LHJ KFPion   = CreateKFParticle(exTrackPi, 0.13957f, 1);      
      

      KFParticle3LHJ KFHypertriton;     
      KFHypertriton.AddDaughters(KFPion, KFHelium, 2); // (daughter1, daughter2, constructMethod (0 or 2)

      KFHypertriton.TransportToDecayVertex();
      
      fSecondaryVertex.SetXYZ(KFHypertriton.GetX(),KFHypertriton.GetY(),KFHypertriton.GetZ());
            
      fSecondaryVertex = fSecondaryVertex - fPrimaryVertex;      

			tSVx = KFHypertriton.GetX();
			tSVy = KFHypertriton.GetY();
			tSVz = KFHypertriton.GetZ();
	      
      TVector3 momvect(KFHypertriton.GetPx(), KFHypertriton.GetPy(), KFHypertriton.GetPz());
      
      if (!KFHypertriton.GetMass() || KFHypertriton.GetMass() == 0. || KFHypertriton.GetP() == 0.) continue;
      if (((KFHypertriton.GetE() + KFHypertriton.GetPz()) / (KFHypertriton.GetE() - KFHypertriton.GetPz())) < 0) continue;
      
      
      tM                 = KFHypertriton.GetMass();
      tPt                = KFHypertriton.GetPt();      
      tP                 = KFHypertriton.GetP();
      tPx                = KFHypertriton.GetPx();
      tPy                = KFHypertriton.GetPy();
      tPz                = KFHypertriton.GetPz();
      tChi2              = KFHypertriton.GetChi2();
      tNDF               = KFHypertriton.GetNDF();
      tDca               = (Double_t) KFHelium.GetDistanceFromParticle(KFPion); 
      tCt                = fSecondaryVertex.Mag() * tM / tP;
      tY                 = KFHypertriton.GetRapidity();
      tCosPA             = TMath::Cos(momvect.Angle(fSecondaryVertex));
      tHePt              = KFHelium.GetPt(); 
      tHePx              = KFHelium.GetPx(); 
      tHePy              = KFHelium.GetPy(); 
      tHePz              = KFHelium.GetPz(); 
      tPiPt              = KFPion.GetPt(); 
      tPiPx              = KFPion.GetPx(); 
      tPiPy              = KFPion.GetPy(); 
      tPiPz              = KFPion.GetPz(); 
      tHeP               = trackHe->GetTPCmomentum();
      tPiP               = trackPi->GetTPCmomentum();
      tnTPCclusterHe     = trackHe->GetTPCNcls(); 
      tnTPCclusterPi     = trackPi->GetTPCNcls();
      tnITSclusterHe     = trackHe->GetITSNcls(); 
      tnITSclusterPi     = trackPi->GetITSNcls();
      
      tHeTOFmass         = GetTOFSignal(*trackHe);        
      tPiTOFmass         = GetTOFSignal(*trackPi);              
      tHeDedxSigma       = Bethe(*trackHe, AliPID::ParticleMass(AliPID::kHe3),    2, fBetheParamsHe);
      tTDedxSigma        = Bethe(*trackHe, AliPID::ParticleMass(AliPID::kTriton), 1, fBetheParamsT);
      tPiDedxSigma       = fPID->NumberOfSigmasTPC(trackPi, AliPID::kPion);
      tHeDedx            = trackHe->GetTPCsignal();
      tPiDedx            = trackPi->GetTPCsignal();
      tHeDca             = TMath::Abs((Double_t) KFHelium.GetDistanceFromVertex(primKFVertex));
      tPiDca             = TMath::Abs((Double_t) KFPion.GetDistanceFromVertex(primKFVertex));
      tImpPiHe           = (Double_t) KFHelium.GetDistanceFromVertex(primKFVertex) * (Double_t) KFPion.GetDistanceFromVertex(primKFVertex);

      AliAODVertex *vtx1 = (AliAODVertex*)trackHe->GetProdVertex();
      tKinkHe = Int_t(vtx1->GetType()) == AliAODVertex::kKink ? 1. : 0.;
      tTPCrefitHe = (trackHe->GetStatus() & AliAODTrack::kTPCrefit) != 0 ? 1. : 0.;
      tITSrefitHe = (trackHe->GetStatus() & AliAODTrack::kITSrefit) != 0 ? 1. : 0.;

      vtx1 = (AliAODVertex*)trackPi->GetProdVertex();
      tKinkPi = Int_t(vtx1->GetType()) == AliAODVertex::kKink ? 1. : 0.;
      tTPCrefitPi = (trackPi->GetStatus() & AliAODTrack::kTPCrefit) != 0 ? 1. : 0.;
      tITSrefitPi = (trackPi->GetStatus() & AliAODTrack::kITSrefit) != 0 ? 1. : 0.;

      if (tCosPA < 0.9)  continue;
      if (tDca > 2.5)    continue;
      if (tM > 3.1)      continue; 
      
      // MC part
      fTree->Fill();
    }
  }
  
  PostData(1, fHistogramList);
  PostData(2, fTree);
  PostData(3, fTreeMCGen);
}
//_____________________________________________________________________________
KFVertex AliAnalysisTaskHypTritTestAOD::CreateKFVertex(const AliVVertex& vertex){
  
  /// GetTrack parameters
  Double_t param[6];
  Double_t cov[6];
  
  vertex.GetXYZ(param);
  vertex.GetCovarianceMatrix(cov);
  
  KFPVertex kfpVtx;
  /// Set the values
  Float_t paramF[3] = {(Float_t) param[0],(Float_t) param[1],(Float_t) param[2]};
  kfpVtx.SetXYZ(paramF);
  Float_t covF[6] = {(Float_t) cov[0],(Float_t) cov[1],(Float_t) cov[2],
    (Float_t) cov[3],(Float_t) cov[4],(Float_t) cov[5]};
  kfpVtx.SetCovarianceMatrix(covF);
  KFVertex KFVtx(kfpVtx);
  return KFVtx;
}
//_____________________________________________________________________________
KFParticle3LHJ AliAnalysisTaskHypTritTestAOD::CreateKFParticle(AliExternalTrackParam& track, float Mass, Int_t Charge) {

  Double_t fP[6];
  track.GetXYZ(fP);
  track.PxPyPz(fP+3);
  Int_t fQ = track.Charge()*TMath::Abs(Charge);
  fP[3] *= TMath::Abs(Charge);
  fP[4] *= TMath::Abs(Charge);
  fP[5] *= TMath::Abs(Charge);

  Double_t pt=1./TMath::Abs(track.GetParameter()[4]) * TMath::Abs(Charge);
  Double_t cs=TMath::Cos(track.GetAlpha()), sn=TMath::Sin(track.GetAlpha());
  Double_t r=TMath::Sqrt((1.-track.GetParameter()[2])*(1.+track.GetParameter()[2]));

  Double_t m00=-sn, m10=cs;
  Double_t m23=-pt*(sn + track.GetParameter()[2]*cs/r), m43=-pt*pt*(r*cs - track.GetParameter()[2]*sn);
  Double_t m24= pt*(cs - track.GetParameter()[2]*sn/r), m44=-pt*pt*(r*sn + track.GetParameter()[2]*cs);
  Double_t m35=pt, m45=-pt*pt*track.GetParameter()[3];

  m43*=track.GetSign();
  m44*=track.GetSign();
  m45*=track.GetSign();

  const Double_t *cTr = track.GetCovariance();
  Double_t fC[21];
  fC[0 ] = cTr[0]*m00*m00;
  fC[1 ] = cTr[0]*m00*m10; 
  fC[2 ] = cTr[0]*m10*m10;
  fC[3 ] = cTr[1]*m00; 
  fC[4 ] = cTr[1]*m10; 
  fC[5 ] = cTr[2];
  fC[6 ] = m00*(cTr[3]*m23 + cTr[10]*m43); 
  fC[7 ] = m10*(cTr[3]*m23 + cTr[10]*m43); 
  fC[8 ] = cTr[4]*m23 + cTr[11]*m43; 
  fC[9 ] = m23*(cTr[5]*m23 + cTr[12]*m43)  +  m43*(cTr[12]*m23 + cTr[14]*m43);
  fC[10] = m00*(cTr[3]*m24 + cTr[10]*m44); 
  fC[11] = m10*(cTr[3]*m24 + cTr[10]*m44); 
  fC[12] = cTr[4]*m24 + cTr[11]*m44; 
  fC[13] = m23*(cTr[5]*m24 + cTr[12]*m44)  +  m43*(cTr[12]*m24 + cTr[14]*m44);
  fC[14] = m24*(cTr[5]*m24 + cTr[12]*m44)  +  m44*(cTr[12]*m24 + cTr[14]*m44);
  fC[15] = m00*(cTr[6]*m35 + cTr[10]*m45); 
  fC[16] = m10*(cTr[6]*m35 + cTr[10]*m45); 
  fC[17] = cTr[7]*m35 + cTr[11]*m45; 
  fC[18] = m23*(cTr[8]*m35 + cTr[12]*m45)  +  m43*(cTr[13]*m35 + cTr[14]*m45);
  fC[19] = m24*(cTr[8]*m35 + cTr[12]*m45)  +  m44*(cTr[13]*m35 + cTr[14]*m45); 
  fC[20] = m35*(cTr[9]*m35 + cTr[13]*m45)  +  m45*(cTr[13]*m35 + cTr[14]*m45);

  KFParticle3LHJ *part = new KFParticle3LHJ();
  part->Create(fP,fC,fQ,Mass);
  return *part;
}
//_____________________________________________________________________________
void AliAnalysisTaskHypTritTestAOD::Terminate(const Option_t*) {
  if (!GetOutputData(0)) return;
}
//_____________________________________________________________________________
/// Loops over MC stack and matches generated particles with reconstructed particles
/// \param stack MC stack
void AliAnalysisTaskHypTritTestAOD::MCStackLoop(AliMCEvent* mcEvent) {
  
  Int_t nV0Gen = 0;
  for (Int_t istack = 0; istack < mcEvent->GetNumberOfTracks(); istack++) {

    AliAODMCParticle *tparticleMother = (AliAODMCParticle*) mcEvent->GetTrack(istack);
    if (!tparticleMother) continue;
    Long_t pdgCodeMother = tparticleMother->GetPdgCode();
    if (TMath::Abs(pdgCodeMother) != 1010010030) continue;
    
    
    AliAODMCParticle *he3 = 0; AliAODMCParticle *pi = 0;
    for (int daughteriD = tparticleMother->GetDaughterFirst(); daughteriD <= tparticleMother->GetDaughterLast(); daughteriD++) {
      AliAODMCParticle *tparticleDaughter = (AliAODMCParticle*) mcEvent->GetTrack(TMath::Abs(daughteriD));
      //if (!(tparticleDaughter && tparticleDaughter->IsSecondaryFromWeakDecay())) continue;//mcEvent->IsSecondaryFromWeakDecay(daughteriD)
      if(!tparticleDaughter) continue;
      if (TMath::Abs(tparticleDaughter->GetPdgCode()) == 1000020030)
        he3 = tparticleDaughter;
      if (TMath::Abs(tparticleDaughter->GetPdgCode()) == 211)
        pi = tparticleDaughter;   
    }
    
    if (!he3 || !pi) continue;
    Double_t posx = he3->Xv();
    Double_t posy = he3->Yv();
    Double_t posz = he3->Zv();
    Double_t disx = posx - tparticleMother->Xv();
    Double_t disy = posy - tparticleMother->Yv();
    Double_t disz = posz - tparticleMother->Zv();
    Double_t distance = TMath::Sqrt(disx*disx + disy*disy + disz*disz );
    tM = tparticleMother->M();
    tP = tparticleMother->P();
    tEta = tparticleMother->Eta();
    tPhi = tparticleMother->Phi();        
    tPt = tparticleMother->Pt();    
    tY = tparticleMother->Y();
    tCt = distance * tparticleMother->M() / tparticleMother->P();
    tZ = 1;
    if (pdgCodeMother == -1010010030) 
      tZ = -1; 
    
    fTreeMCGen->Fill();
  }

}
//_____________________________________________________________________________
Float_t AliAnalysisTaskHypTritTestAOD::GetInvPtDevFromBC(Int_t b, Int_t c) {
  //returns d(1/Pt) in c/GeV 
  //in case of no gtu simulation -> return maximum 0.5
  if(b==0 && c==0) return 0.5;
  Int_t tmp = (((b & 0xfff) << 12) ^ 0x800000) - 0x800000;
  tmp += (c & 0xfff);
  Float_t invPtDev = tmp * 0.000001;
  return invPtDev;
}
//_____________________________________________________________________________
/// Calculates number of sigma deviation from expected dE/dx in TPC
/// \param track particle track
/// \param mass mass hypothesis of particle
/// \param charge particle charge hypothesis
/// \param params Parameters of Aleph parametrization of Bethe Energy-loss
Double_t AliAnalysisTaskHypTritTestAOD::Bethe(const AliAODTrack& track, Double_t mass, Int_t charge, Double_t* params){
  Double_t expected = charge*charge*AliExternalTrackParam::BetheBlochAleph(charge*track.GetTPCmomentum()/mass,params[0],params[1],params[2],params[3],params[4]);
  Double_t sigma = expected*params[5];
  if (TMath::IsNaN(expected)) return -999;
  return (track.GetTPCsignal() - expected) / sigma;
}
//_____________________________________________________________________________
// _________________________________________________ //
Double_t AliAnalysisTaskHypTritTestAOD::GetTOFSignal(const AliAODTrack& track) {
  Float_t mass = 0;
  Float_t time = -1;
  Float_t beta = 0;
  Float_t gamma = 0;
  Float_t length = 0;
  Float_t time0 = 0;
  length = track.GetIntegratedLength();
  time0 = fPID->GetTOFResponse().GetStartTime(track.P());//fPID->GetTOFResponse().GetTimeZero();
  time = track.GetTOFsignal() - time0;
  if (time > 0) {
    beta = length / (2.99792457999999984e-02 * time);
    if (beta * beta >= 1) return -99;
    gamma = 1 / TMath::Sqrt(1 - beta * beta);
    if (gamma * gamma <= 1) return -99;
    mass = (track.GetTPCmomentum()) / TMath::Sqrt(gamma * gamma - 1); // using inner TPC mom. as approx.
  }
  return mass;
}
// _________________________________________________ //
//_____________________________________________________________________________
void AliAnalysisTaskHypTritTestAOD::SetBetheBlochParams(Int_t runNumber) {
	// set Bethe-Bloch parameter
	if (runNumber >= 252235 && runNumber <= 265589) { // 2016 pp data
		fYear = 2016;
		// He3
		fBetheParamsHe[0] = 1.81085;
		fBetheParamsHe[1] = 29.4656;
		fBetheParamsHe[2] = 0.0458225;
		fBetheParamsHe[3] = 2.08689;
		fBetheParamsHe[4] = 2.28772;
		fBetheParamsHe[5] = 0.06;
		// Triton
		fBetheParamsT[0] = 0.427978;
		fBetheParamsT[1] = 105.46;
		fBetheParamsT[2] = -7.08642e-07;
		fBetheParamsT[3] = 2.23332;
		fBetheParamsT[4] = 18.8231;
		fBetheParamsT[5] = 0.06;
	}
	if (runNumber > 265589 && runNumber <= 267166) { // 2016 p-Pb data
		fYear = 2016;
		// He3
		fBetheParamsHe[0] = 0.715489;
		fBetheParamsHe[1] = 59.5463;
		fBetheParamsHe[2] = 4.44487e-12;
		fBetheParamsHe[3] = 2.69874;
		fBetheParamsHe[4] = 24.063;
		fBetheParamsHe[5] = 0.04725;
		// Triton
		fBetheParamsT[0] = 0.223948;
		fBetheParamsT[1] = 180.564;
		fBetheParamsT[2] = -3.03884e-10;
		fBetheParamsT[3] = 2.30095;
		fBetheParamsT[4] = 34.2269;
		fBetheParamsT[5] = 0.06517;	
	} 	
	if (runNumber >= 270581 && runNumber <= 282704) { // 2017 pp data
		fYear = 2017;
		// He3
		fBetheParamsHe[0] = 3.20025;
		fBetheParamsHe[1] = 16.4971;
		fBetheParamsHe[2] = -0.0116571;
		fBetheParamsHe[3] = 2.3152;
		fBetheParamsHe[4] = 3.11135;
		fBetheParamsHe[5] = 0.06;
		// Triton
		fBetheParamsT[0] = 0.420434;
		fBetheParamsT[1] = 106.102;
		fBetheParamsT[2] = -3.15587e-07;
		fBetheParamsT[3] = 2.32499;
		fBetheParamsT[4] = 21.3439;
		fBetheParamsT[5] = 0.06;
	}
	if (runNumber >= 285009 && runNumber <= 294925) { // 2018 pp data
		fYear = 2018;
		// He3
		fBetheParamsHe[0] = 1.81085;
		fBetheParamsHe[1] = 29.4656;
		fBetheParamsHe[2] = 0.0458225;
		fBetheParamsHe[3] = 2.08689;
		fBetheParamsHe[4] = 2.28772;
		fBetheParamsHe[5] = 0.06;
		// Triton
		fBetheParamsT[0] = 0.427978;
		fBetheParamsT[1] = 105.46;
		fBetheParamsT[2] = -7.08642e-07;
		fBetheParamsT[3] = 2.23332;
		fBetheParamsT[4] = 18.8231;
		fBetheParamsT[5] = 0.06;
	}

	/* __ 2018 PbPb pass 3__
	if (runNumber >= 295581 && runNumber <= 297624) {
	  fBetheParamsT[0] = 0.648689;
	  fBetheParamsT[1] = 56.6706;
	  fBetheParamsT[2] = -1.63243e-10;
	  fBetheParamsT[3] = 2.46921;
	  fBetheParamsT[4] = 16.8531;
	  fBetheParamsT[5] = 0.06;

	  fBetheParamsHe[0] = 1.70184;
	  fBetheParamsHe[1] = 28.4426;
	  fBetheParamsHe[2] = 3.21871e-12;
	  fBetheParamsHe[3] = 2.06952;
	  fBetheParamsHe[4] = 2.77971;
	  fBetheParamsHe[5] = 0.06;
	}*/
 //  2018 PbPb pass 1 //
  if(runNumber >= 295581 && runNumber <= 297624){
     fBetheParamsT[0] = 0.669634;
     fBetheParamsT[1] = 53.1497;
     fBetheParamsT[2] =-1.32853e-08;
     fBetheParamsT[3] = 2.5775;
     fBetheParamsT[4] = 17.7607;
     fBetheParamsT[5] = 0.06;
     fBetheParamsHe[0] = 1.50582;
     fBetheParamsHe[1] = 33.7232;
     fBetheParamsHe[2] = -0.0923749;
     fBetheParamsHe[3] = 2.00901;
     fBetheParamsHe[4] = 2.28772;
     fBetheParamsHe[5] = 0.06;
     }


}
//_____________________________________________________________________________

//_____________________________________________________________________________
