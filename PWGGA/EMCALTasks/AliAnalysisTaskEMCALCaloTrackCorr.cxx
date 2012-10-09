#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TClonesArray.h"
#include "TString.h"

#include "AliVEvent.h"
#include "AliStack.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskEMCALCaloTrackCorr.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliVCaloCells.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliVVertex.h"
#include "AliESDtrackCuts.h"
#include "AliAODTrack.h"
#include "AliCentrality.h"
#include "AliCaloTrackParticle.h"
#include "AliTriggerAnalysis.h"
#include "AliFiducialCut.h"
#include "AliEventplane.h"
#include "AliAODHeader.h"


ClassImp(AliAnalysisTaskEMCALCaloTrackCorr)

//________________________________________________________________________
AliAnalysisTaskEMCALCaloTrackCorr::AliAnalysisTaskEMCALCaloTrackCorr(const char *name) 
 :AliAnalysisTaskSE(name),
  fManager(NULL),     fInputHandler(NULL),
  fEvent(NULL),       fMCEvent(NULL),     fStack(NULL),        
  fCentrality(0x0),   fEventPlane(0x0),   
  fEMCALRecU(0x0),    fEMCALGeom(0x0),    fESDtrackCuts(0x0), /* fFidCut(0x0),*/
  outputContainer(0), 

  fEMCALGeomName("EMCAL_COMPLETEV1"),
  fCentralityClass(""),  fCentralityBinMin(0), fCentralityBinMax(0),
  fEventPlaneMethod(""), fEventTriggerMaks(0), fNCentralityBins(0),      fNEventPlaneBins(0), 
  fHistoPtBins(0),       fHistoPtMax(0.),      fHistoPtMin(0.),
  fHistoPhiBins(0),      fHistoPhiMax(0.),     fHistoPhiMin(0.),
  fHistoEtaBins(0),      fHistoEtaMax(0.),     fHistoEtaMin(0.),
  fMinNCells(0)   ,      fMinE(0),             fMinDistBad(0),
  fL0CutMin(0),          fL0CutMax(0),         fTimeCutMin(0),
  fTimeCutMax(0),        fPhotonPairTimeCut(0),fEMCALDPhiCut(0),
  fEMCALDEtaCut(0),      fZVertexCut(0),       fDebug(0), 
  fAnaMesonType(0),      fAsymmetryCut(0),     fDataType(0),     
  fTrackFilterMask(0),

  fInvMassMinCut(0),     fInvMassMaxCut(0),    fLeftBandMinCut(0),
  fLeftBandMaxCut(0),    fRightBandMinCut(0),  fRightBandMaxCut(0),

  kMC(0),                  kNeutralMesonHistos(0),   kDoMixEventsAna(0), 
  kDoPhotonCorrAna(0),     kDoAsymmetryCut(0),       kDoSelectHybridTracks(0),
  kDoMesonFill(0),         kDoMesonCorrAna(0),       kDoIsolatedAna(0),
  kDoTrackMultBins(0),     kUELeftRight(0),          kUENearAway(0),
  kDecayPhotonCorr(0),     kAnaMCTruthCorr(0),       kAnaMCPrimaryCorr(0), 
  kAnaPi0Prim(0),          kAnaEtaPrim(0),           kAnaPhotonPrim(0),
  kMakeAbsoluteLeading(0), kMakeNearSideLeading(0),  kTwoTracksCorr(0),
  kPhotonInAcceptance(0),  kAnaDecayMapping(0),      /*fCheckFidCut(0),*/
  kEventTriggerAtSE(0),    kPhotonPairTimeCut(0),    kDoPhotonIDCut(0),

  fhNEvents(0),            fnEvents(0),             fhNEventsAnalyized(0),  fEventAnalyized(0),
  fPhotonEvent(0),         fPhotonPairEvent(0),     fCTSEvent(0),           
  nPhotonsEMCAL(0),        nTracksCTS(0),           fAnaTypeInIsolated(0),
  nMixedEvents(0),       
  fSetConeR(0),            fSetPtThreshold(0),      fSetSumPtThreshold(0),
  fSetPtFraction(0),       fICMethod(0),            fParticlesInCone(0),  

  fTriggPtArray(0),        fNTriggPtBins(0),        fptTriggerBegin(0),     fptTriggerEnd(0),
  fAssocPtArray(0),        fNAssocPtBins(0),        fptAssociatedBegin(0),

  fDeltaPhiMaxCut(0.),     fDeltaPhiMinCut(0.),     fUeDeltaPhiSize(0),  
  fUeDeltaPhiFix(0),       fDeltaPhiHRSize(0),      

  fhPhotonE(0),            fhPhotonPtPhi(0),        fhPhotonPtEta(0),     fhPhotonPhiEta(0),  
  fhMesonE(0),             fhMesonPtPhi(0),         fhMesonPtEta(0),      fhMesonPhiEta(0),  

  fhAnglePairNoCut(0),     fhInvMassPairNoCut(0),   fhAsyNoCut(0),
  fhInvMassPairAsyCut(0),  fhAnglePairAsyCut(0),    fhInvMassPairPhi(0),
  fhInvMassPairEta(0),     fhInvMassPairAllCut(0),  fhAnglePairAllCut(0),
  fhAsyAllCut(0),

  fhPi0DecayPhoton1(0),    fhPi0DecayPhoton1Dphi(0),     fhDecayPhoton1Pi0Dphi(0),
  fhPi0DecayPhoton2(0),    fhPi0DecayPhoton2Dphi(0),     fhDecayPhoton2Pi0Dphi(0),
  fhDecayPhoton1Photon2(0),fhDecayPhoton1Photon2Dphi(0), fhDecayPhoton2Photon1Dphi(0),

   

  fhNtracksAll(0),         fhNtracksEMC7(0),        fhNtracksAnyINT(0),
  fhNtracksCentral(0),     fhNtracksSemiCentral(0), fhNtracksOtherTirgger(0),   
  fhNtracksCorr(0),        fhTrackPtPhi(0),          
  fhTrackPtEta(0),         fhTrackPhiEta(0),        fhPtPhiLeading(0),
  fhPtEtaLeading(0),       fhMixPtPhiLeading(0),    fhMixPtEtaLeading(0),
  fhDPhiTriggPtAssocPt(0), fhDEtaTriggPtAssocPt(0),
  fhAssocPtTriggPt(0),     fhxELogTriggPt(0),       fhpoutTriggPt(0),
  fhzTTriggPt(0),          fhxETriggPt(0),       
  fhAssocPtTriggPtHR(0),   fhxELogTriggPtHR(0),     fhpoutTriggPtHR(0),
  fhzTTriggPtHR(0),        fhxETriggPtHR(0),
  fhNUeAssocPtTriggPt(0),  fhNUepoutTriggPt(0),     fhNUexETriggPt(0),
  fhNUezTTriggPt(0),       fhNUexELogTriggPt(0),    fhNUeDPhiDEta(0),
  fhAUeAssocPtTriggPt(0),  fhAUepoutTriggPt(0),     fhAUezTTriggPt(0), 
  fhAUexETriggPt(0),       fhAUexELogTriggPt(0),    fhAUeDPhiDEta(0),

  fhMCPtPhiLeading(0),     fhMCPtEtaLeading(0),    
  fhMCAssocPtTriggPt(0),   fhMCxELogTriggPt(0),    fhMCpoutTriggPt(0),
  fhMCzTTriggPt(0),        fhMCxETriggPt(0),       
  fhMCAssocPtTriggPtHR(0), fhMCxELogTriggPtHR(0),  fhMCpoutTriggPtHR(0),
  fhMCzTTriggPtHR(0),      fhMCxETriggPtHR(0),
  fhMCNUeAssocPtTriggPt(0),fhMCNUepoutTriggPt(0),  fhMCNUexETriggPt(0),
  fhMCNUezTTriggPt(0),     fhMCNUexELogTriggPt(0),
  fhMCAUeAssocPtTriggPt(0),fhMCAUepoutTriggPt(0),  fhMCAUezTTriggPt(0), 
  fhMCAUexETriggPt(0),     fhMCAUexELogTriggPt(0),
  fhDPhiAssocPt15T(0),     fhDEtaAssocPt15T(0),    fhMixDPhiAssocPt15T(0),
  fhMixDEtaAssocPt15T(0),  fhMCDPhiAssocPt15T(0),  fhMCDEtaAssocPt15T(0)

{
  // Constructor
  for(Int_t imix=0;imix<10;imix++){
    for(Int_t jmix=0;jmix<10;jmix++){
     for(Int_t kmix=0;kmix<10;kmix++)
      fListMixEvents[imix][jmix][kmix]=0x0;
    }
  }

  // Default constructor.
  for(Int_t i = 0; i< GetNAssocPtBins(); i++){
   fhDPhiTriggPtT[i] = 0;
   fhDEtaTriggPtT[i] = 0;
   fhMixDPhiTriggPtT[i] = 0;
   fhMixDEtaTriggPtT[i] = 0;
  }
 
  for(Int_t j = 0; j < GetNTriggPtBins(); j++){
   for(Int_t k = 0; k< GetNAssocPtBins(); k++){
    fhDPhiSumPtBin[j][k]   = 0;
    fhDEtaSumPtBin[j][k]   = 0;
    fhDPhiDEtaBin[j][k]    = 0;
    fhMixDPhiDEtaBin[j][k]    = 0;
   }
  }
 
  for(Int_t j1 = 0; j1 < GetNTriggPtBins(); j1++){  
   fhDPhiAssocPtA[j1] = 0;
   fhDEtaAssocPtA[j1] = 0;
   fhMixDPhiAssocPtA[j1] = 0;
   fhMixDEtaAssocPtA[j1] = 0;
  } 

   /////MC
  for(Int_t i = 0; i< GetNAssocPtBins(); i++){
   fhMCDPhiTriggPtT[i] = 0;
   fhMCDEtaTriggPtT[i] = 0;
  }
  
  for(Int_t j = 0; j < GetNTriggPtBins(); j++){
   for(Int_t k = 0; k< GetNAssocPtBins(); k++){
    fhMCDPhiSumPtBin[j][k]   = 0;
    fhMCDEtaSumPtBin[j][k]   = 0;
    fhMCDPhiDEtaBin[j][k]    = 0;
   }
  }
 
  for(Int_t j1 = 0; j1 < GetNTriggPtBins(); j1++){  
   fhMCDPhiAssocPtA[j1] = 0;
   fhMCDEtaAssocPtA[j1] = 0;
  }

  //Initialize parameters
  InitParameters();
  
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskEMCALCaloTrackCorr::UserCreateOutputObjects()
{
 
  fEMCALGeom = AliEMCALGeometry::GetInstance(fEMCALGeomName);
  fEMCALRecU = new AliEMCALRecoUtils(); 

  Int_t nptbins  = GetHistoPtBins();
  Int_t nphibins = GetHistoPhiBins();
  Int_t netabins = GetHistoEtaBins();
  Float_t ptmax  = GetHistoPtMax();
  Float_t phimax = GetHistoPhiMax();
  Float_t etamax = GetHistoEtaMax();
  Float_t ptmin  = GetHistoPtMin();
  Float_t phimin = GetHistoPhiMin();
  Float_t etamin = GetHistoEtaMin();

  Int_t ndphibins = 350;
  Float_t dphimin = -2.;
  Float_t dphimax = 5.;
 
  Int_t ndetabins = 360;
  Float_t detamin = -1.8;
  Float_t detamax = 1.8;

  Int_t   nzTbins  = 150;
  Float_t zTmin    = 0.;
  Float_t zTmax    = 3.;

  Int_t   nxEbins  = 150;
  Float_t xEmin    = 0.;
  Float_t xEmax    = 3.;

  Int_t   nxElogbins  = 180;
  Float_t xElogmin    = 0.;
  Float_t xElogmax    = 9.;

  Int_t   nmassbins = 200;
  Float_t massmin   = 0.;
  Float_t massmax   = 1.;

  Int_t   nanglebins = 300;
  Float_t anglemin   = 0.;
  Float_t anglemax   = 0.6;

  outputContainer = new TList();
  outputContainer->SetOwner(kTRUE);

  fhNEvents = new TH1F ("hNEvents","Number of all events",1,0.5,1.5);
  outputContainer->Add(fhNEvents);

  fhNEventsAnalyized = new TH1F ("hNEventsAnalyized","Number of events analyzed Last",1,0.5,1.5);
  outputContainer->Add(fhNEventsAnalyized);	

  fhPhotonE  = new TH1F("hPhotonE","Number of #gamma over calorimeter vs energy",nptbins,ptmin,ptmax);
  fhPhotonE->SetYTitle("N");
  fhPhotonE->SetXTitle("E_{#gamma}(GeV)");
  outputContainer->Add(fhPhotonE) ;

  fhPhotonPtPhi  = new TH2F
    ("hPhotonPtPhi","#phi_{#gamma} vs p_{T}",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
  fhPhotonPtPhi->SetYTitle("#phi (rad)");
  fhPhotonPtPhi->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhPhotonPtPhi) ;

  fhPhotonPtEta  = new TH2F
    ("hPhotonPtEta","#eta_{#gamma} vs p_{T}",nptbins,ptmin,ptmax,netabins,etamin,etamax);
  fhPhotonPtEta->SetYTitle("#eta");
  fhPhotonPtEta->SetXTitle("p_{T #gamma} (GeV/c)");
  outputContainer->Add(fhPhotonPtEta) ;

  fhPhotonPhiEta  = new TH2F
  ("hPhotonPhiEta","#phi vs #eta",nphibins,phimin,phimax, netabins,etamin,etamax);
  fhPhotonPhiEta->SetXTitle("#phi (rad)");
  fhPhotonPhiEta->SetYTitle("#eta");
  outputContainer->Add(fhPhotonPhiEta);

  if(kDoMesonFill){
   fhMesonE  = new TH1F("hMesonE","Number of #pi^{0}(#eta) over calorimeter vs energy",nptbins,ptmin,ptmax);
   fhMesonE->SetYTitle("N");
   fhMesonE->SetXTitle("E_{#meson}(GeV)");
   outputContainer->Add(fhMesonE) ;

   fhMesonPtPhi  = new TH2F
    ("hMesonPtPhi","#phi vs p_{T} of #pi^{0}(#eta)",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
   fhMesonPtPhi->SetYTitle("#phi (rad)");
   fhMesonPtPhi->SetXTitle("p_{T meson} (GeV/c)");
   outputContainer->Add(fhMesonPtPhi) ;

   fhMesonPtEta  = new TH2F
    ("hMesonPtEta","#eta vs p_{T}  of #pi^{0}(#eta)",nptbins,ptmin,ptmax,netabins,etamin,etamax);
   fhMesonPtEta->SetYTitle("#eta");
   fhMesonPtEta->SetXTitle("p_{T meson} (GeV/c)");
   outputContainer->Add(fhMesonPtEta) ;

   fhMesonPhiEta  = new TH2F
    ("hMesonPhiEta","#phi vs #eta of #pi^{0}(#eta)",nphibins,phimin,phimax, netabins,etamin,etamax);
   fhMesonPhiEta->SetXTitle("#phi (rad)");
   fhMesonPhiEta->SetYTitle("#eta");
   outputContainer->Add(fhMesonPhiEta);

   if(kNeutralMesonHistos){
    fhAnglePairNoCut  = new TH2F("hAnglePairNoCut","Angle between all #gamma pair vs E_{#pi^{0}}",
                                 nptbins,ptmin,ptmax,nanglebins,anglemin,anglemax);
    fhAnglePairNoCut->SetYTitle("Angle (rad)");
    fhAnglePairNoCut->SetXTitle("E_{ #pi^{0}} (GeV)");

    fhAsyNoCut  = new TH2F("hAsymmetryNoCut","Asymmetry of all #gamma pair vs E_{#pi^{0}}",
                           nptbins,ptmin,ptmax,100,0,1);
    fhAsyNoCut->SetYTitle("Asymmetry");
    fhAsyNoCut->SetXTitle("E_{ #pi^{0}} (GeV)");

    fhInvMassPairNoCut  = new TH2F("hInvMassPairNoCut","Invariant Mass of all #gamma pair vs E_{#pi^{0}}",
                                   nptbins,ptmin,ptmax,nmassbins, massmin, massmax);
    fhInvMassPairNoCut->SetYTitle("hInvariant Mass (GeV/c^{2})");
    fhInvMassPairNoCut->SetXTitle("E_{ #pi^{0}} (GeV)");

    outputContainer->Add(fhAnglePairNoCut) ;
    outputContainer->Add(fhAsyNoCut) ;
    outputContainer->Add(fhInvMassPairNoCut) ;


    if(kDoAsymmetryCut) {
     fhAnglePairAsyCut  = new TH2F("hAnglePairAsymmetryCut","AnglePairAsymmetryCut",
                                        nptbins,ptmin,ptmax,nanglebins,anglemin,anglemax);
     fhAnglePairAsyCut->SetYTitle("Angle (rad)");
     fhAnglePairAsyCut->SetXTitle("E_{ #pi^{0}} (GeV)");

     fhInvMassPairAsyCut  = new TH2F("hInvMassPairAsymmetryCut","Invariant Mass of #gamma pair vs E_{#pi^{0}}",
                                                                  nptbins,ptmin,ptmax,nmassbins, massmin, massmax);
     fhInvMassPairAsyCut->SetYTitle("Invariant Mass (GeV/c^{2})");
     fhInvMassPairAsyCut->SetXTitle("E_{#pi^{0}}(GeV)");

     fhInvMassPairPhi = new TH2F("hInvMassPairPhi", "M_{#gamma#gamma} vs #phi",
                                nptbins,ptmin,ptmax,nphibins, phimin, phimax);
     fhInvMassPairPhi->SetYTitle("#phi (rad)");
     fhInvMassPairPhi->SetXTitle("E_{ #pi^{0}} (GeV)");

     fhInvMassPairEta = new TH2F("hInvMassPairEta","M_{#gamma#gamma} vs #eta",
                              nptbins,ptmin,ptmax, netabins, etamin, etamax);
     fhInvMassPairEta->SetYTitle("#eta");
     fhInvMassPairEta->SetXTitle("E_{ #pi^{0}} (GeV)");

     outputContainer->Add(fhAnglePairAsyCut) ;
     outputContainer->Add(fhInvMassPairAsyCut) ;
     outputContainer->Add(fhInvMassPairPhi);
     outputContainer->Add(fhInvMassPairEta);
    }

    fhAnglePairAllCut  = new TH2F("hAnglePairAllCut", "Angle between all #gamma pair (opening angle + asymmetry + inv mass cut) vs E_{#pi^{0}}",
                                 nptbins,ptmin,ptmax,nanglebins,anglemin,anglemax);
    fhAnglePairAllCut->SetYTitle("Angle (rad)");
    fhAnglePairAllCut->SetXTitle("E_{ #pi^{0}} (GeV)");

    fhInvMassPairAllCut  = new TH2F("hInvMassPairAllCut","Invariant Mass of #gamma pair (opening angle + asymmetry + invmass cut) vs E_{#pi^{0}}",
                                    nptbins,ptmin,ptmax,nmassbins, massmin, massmax);
    fhInvMassPairAllCut->SetYTitle("Invariant Mass (GeV/c^{2})");
    fhInvMassPairAllCut->SetXTitle("E_{#pi^{0}}(GeV)");

    fhAsyAllCut  = new TH2F("hAsymmetryAllCut", "Asymmetry of #gamma pair (opening angle+invmass cut) vs E_{#pi^{0}}",
                           nptbins,ptmin,ptmax,100,0,1);
    fhAsyAllCut->SetYTitle("Asymmetry");
    fhAsyAllCut->SetXTitle("E_{#pi^{0}}(GeV)");

    outputContainer->Add(fhAnglePairAllCut) ;
    outputContainer->Add(fhAsyAllCut) ;
    outputContainer->Add(fhInvMassPairAllCut) ;
   } 

   fhPi0DecayPhoton1  = new TH2F("hPi0DecayPhoton1","#pi^{0} vs its decay #gamma1 ",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
   fhPi0DecayPhoton1->SetYTitle("p_{T #gamma1}(GeV/c)");
   fhPi0DecayPhoton1->SetXTitle("p_{T #pi^{0}}(GeV/c)");
   outputContainer->Add(fhPi0DecayPhoton1);

   fhPi0DecayPhoton1Dphi  = new TH2F("hPi0DecayPhoton1Dphi","#pi^{0} vs #Delta#phi ",nptbins,ptmin,ptmax,300,0,0.6);
   fhPi0DecayPhoton1Dphi->SetYTitle("#Delta#phi(rad)");
   fhPi0DecayPhoton1Dphi->SetXTitle("p_{T #pi^{0}}(GeV/c)");
   outputContainer->Add(fhPi0DecayPhoton1Dphi);

   fhDecayPhoton1Pi0Dphi  = new TH2F("hDecayPhoton1Pi0Dphi","decay #gamma1 vs #Delta#phi",nptbins,ptmin,ptmax,300,0,0.6);
   fhDecayPhoton1Pi0Dphi->SetYTitle("#Delta#phi(rad)");
   fhDecayPhoton1Pi0Dphi->SetXTitle("p_{T #gamma1}(GeV/c)");
   outputContainer->Add(fhDecayPhoton1Pi0Dphi);

   fhPi0DecayPhoton2  = new TH2F("hPi0DecayPhoton2","#pi^{0} vs its decay #gamma2",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
   fhPi0DecayPhoton2->SetYTitle("p_{T #gamma2}(GeV/c)");
   fhPi0DecayPhoton2->SetXTitle("p_{T #pi^{0}}(GeV/c)");
   outputContainer->Add(fhPi0DecayPhoton2);

   fhPi0DecayPhoton2Dphi  = new TH2F("hPi0DecayPhoton2Dphi","#pi^{0} vs #Delta#phi",nptbins,ptmin,ptmax,300,0,0.6);
   fhPi0DecayPhoton2Dphi->SetYTitle("#Delta#phi(rad)");
   fhPi0DecayPhoton2Dphi->SetXTitle("p_{T #pi^{0}}(GeV/c)");
   outputContainer->Add(fhPi0DecayPhoton2Dphi);

   fhDecayPhoton2Pi0Dphi  = new TH2F("hDecayPhoton2Pi0Dphi","decay #gamma2 vs #Delta#phi",nptbins,ptmin,ptmax,300,0,0.6);
   fhDecayPhoton2Pi0Dphi->SetYTitle("#Delta#phi(rad)");
   fhDecayPhoton2Pi0Dphi->SetXTitle("p_{T #gamma2}(GeV/c)");
   outputContainer->Add(fhDecayPhoton2Pi0Dphi);

   fhDecayPhoton1Photon2  = new TH2F("hDecayPhoton1Photon2","#gamma one vs #gamma two",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
   fhDecayPhoton1Photon2->SetYTitle("p_{T #gamma2}(GeV/c)");
   fhDecayPhoton1Photon2->SetXTitle("p_{T #gamma1}(GeV/c)");
   outputContainer->Add(fhDecayPhoton1Photon2);

   fhDecayPhoton1Photon2Dphi  = new TH2F("hDecayPhoton1Photon2Dphi","#gamma one vs #Delta#phi",nptbins,ptmin,ptmax,300,0,0.6);
   fhDecayPhoton1Photon2Dphi->SetYTitle("#Delta#phi(rad)");
   fhDecayPhoton1Photon2Dphi->SetXTitle("p_{T #gamma1}(GeV/c)");
   outputContainer->Add(fhDecayPhoton1Photon2Dphi);

   fhDecayPhoton2Photon1Dphi  = new TH2F("hDecayPhoton2Photon1Dphi","#gamma two vs #Delta#phi",nptbins,ptmin,ptmax,300,0,0.6);
   fhDecayPhoton2Photon1Dphi->SetZTitle("#Delta#phi(rad)");
   fhDecayPhoton2Photon1Dphi->SetXTitle("p_{T #gamma2}(GeV/c)");
   outputContainer->Add(fhDecayPhoton2Photon1Dphi);
  }

  fhNtracksAll=new TH1F("hNtracksAll","Number of tracks w/o event trigger",2000,0,2000);
  outputContainer->Add(fhNtracksAll);

  fhNtracksEMC7=new TH1F("hNtracksEMC7","Number of tracks w/ event trigger kEMC7",2000,0,2000);
  outputContainer->Add(fhNtracksEMC7);

  fhNtracksAnyINT=new TH1F("hNtracksAnyINT","Number of tracks w/ event trigger kAnyINT",2000,0,2000);
  outputContainer->Add(fhNtracksAnyINT);

  fhNtracksCentral=new TH1F("hNtracksCentral","Number of tracks w/ event trigger kCentral",2000,0,2000);
  outputContainer->Add(fhNtracksCentral);

  fhNtracksSemiCentral=new TH1F("hNtracksSemiCentral","Number of tracks w/ event trigger kSemiCentral",2000,0,2000);
  outputContainer->Add(fhNtracksSemiCentral);

  fhNtracksOtherTirgger=new TH1F("hNtracksOtherTirgger","Number of tracks w/ event trigger other ",2000,0,2000);
  outputContainer->Add(fhNtracksOtherTirgger);
 
  fhNtracksCorr=new TH1F("hNtracksCorr","Number of tracks w/ event trigger same with correlation",2000,0,2000);
  outputContainer->Add(fhNtracksCorr);

  fhTrackPtPhi  = new TH2F ("hTrackPtPhi","p_{T}  and #phi distribution of tracks", nptbins,ptmin,ptmax, nphibins, phimin, phimax);
  outputContainer->Add(fhTrackPtPhi);

  fhTrackPtEta  = new TH2F ("hTrackPtEta","p_{T} and #eta distribution of tracks",nptbins,ptmin,ptmax, netabins,etamin,etamax);
  outputContainer->Add(fhTrackPtEta);

  fhTrackPhiEta  = new TH2F ("hTrackPhiEta","#phi and #eta distribution of tracks",nphibins, phimin, phimax, netabins,etamin,etamax);
  outputContainer->Add(fhTrackPhiEta);

  fhPtPhiLeading  = new TH2F ("hPtPhiLeading","p_{T}  and #phi distribution of leading particles", nptbins,ptmin,ptmax, nphibins, phimin, phimax); 
  outputContainer->Add(fhPtPhiLeading);
  
  fhPtEtaLeading  = new TH2F ("hPtEtaLeading","p_{T} and #eta distribution of leading",nptbins,ptmin,ptmax, netabins,etamin,etamax); 
  outputContainer->Add(fhPtEtaLeading);

   fhMixPtPhiLeading  = new TH2F ("hMixPtPhiLeading","p_{T}  and #phi distribution of mixed leading particles", nptbins,ptmin,ptmax, nphibins, phimin, phimax);
  outputContainer->Add(fhMixPtPhiLeading);

  fhMixPtEtaLeading  = new TH2F ("hMixPtEtaLeading","p_{T} and #eta distribution of mixed leading particles",nptbins,ptmin,ptmax, netabins,etamin,etamax);
  outputContainer->Add(fhMixPtEtaLeading);

  for(Int_t i=0; i< GetNAssocPtBins(); i++){
   fhDPhiTriggPtT[i] = new TH2F(Form("hDPhiTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fT",fptTriggerBegin, fptTriggerEnd, fAssocPtArray[i], fAssocPtArray[i+1]),
   Form("DPhiTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fT",fptTriggerBegin, fptTriggerEnd, fAssocPtArray[i], fAssocPtArray[i+1]),ndphibins,dphimin,dphimax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhDPhiTriggPtT[i]);

   fhDEtaTriggPtT[i] = new TH2F(Form("hDEtaTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fT",fptTriggerBegin, fptTriggerEnd, fAssocPtArray[i], fAssocPtArray[i+1]),
   Form("DEtaTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fT",fptTriggerBegin, fptTriggerEnd, fAssocPtArray[i], fAssocPtArray[i+1]), ndetabins,detamin,detamax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhDEtaTriggPtT[i]);

   fhMixDPhiTriggPtT[i] = new TH2F(Form("hMixDPhiTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fT",fptTriggerBegin, fptTriggerEnd, fAssocPtArray[i], fAssocPtArray[i+1]),
   Form("MixDPhiTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fT",fptTriggerBegin, fptTriggerEnd, fAssocPtArray[i], fAssocPtArray[i+1]),ndphibins,dphimin,dphimax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMixDPhiTriggPtT[i]);

   fhMixDEtaTriggPtT[i] = new TH2F(Form("hMixDEtaTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fT",fptTriggerBegin, fptTriggerEnd, fAssocPtArray[i], fAssocPtArray[i+1]),
   Form("MixDEtaTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fT",fptTriggerBegin, fptTriggerEnd, fAssocPtArray[i], fAssocPtArray[i+1]), ndetabins,detamin,detamax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMixDEtaTriggPtT[i]);
  }

  fhDPhiAssocPt15T = new TH2F(Form("hDPhiTriggPtTrigg%1.f_%1.fAssoc1.0_5.0",fptTriggerBegin, fptTriggerEnd),
  Form("DPhiTriggPtTrigg%1.f_%1.fAssoc1.0_5.0",fptTriggerBegin, fptTriggerEnd),ndphibins,dphimin,dphimax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhDPhiAssocPt15T);

  fhDEtaAssocPt15T = new TH2F(Form("hDEtaTriggPtTrigg%1.f_%1.fAssoc1.0_5.0",fptTriggerBegin, fptTriggerEnd),
  Form("DEtaTriggPtTrigg%1.f_%1.fAssoc1.0_5.0",fptTriggerBegin, fptTriggerEnd),ndetabins,detamin,detamax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhDEtaAssocPt15T);

  fhMixDPhiAssocPt15T = new TH2F(Form("hMixDPhiTriggPtTrigg%1.f_%1.fAssoc1.0_5.0",fptTriggerBegin, fptTriggerEnd),
  Form("MixDPhiTriggPtTrigg%1.f_%1.fAssoc1.0_5.0",fptTriggerBegin, fptTriggerEnd),ndphibins,dphimin,dphimax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhMixDPhiAssocPt15T);

  fhMixDEtaAssocPt15T = new TH2F(Form("hMixDEtaTriggPtTrigg%1.f_%1.fAssoc1.0_5.0",fptTriggerBegin, fptTriggerEnd),
  Form("MixDEtaTriggPtTrigg%1.f_%1.fAssoc1.0_5.0",fptTriggerBegin, fptTriggerEnd),ndetabins,detamin,detamax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhMixDEtaAssocPt15T);

  for(Int_t i = 0 ; i < GetNTriggPtBins() ; i++){
   for(Int_t j=0; j<GetNAssocPtBins(); j++){
    fhDPhiSumPtBin[i][j] = new TH2F(Form("hDPhiSumPtBinTrigg%1.f_%1.fAssoc%1.f_%1.f",fTriggPtArray[i], fTriggPtArray[i+1],fAssocPtArray[j],fAssocPtArray[j+1]),
    Form("DPhiSumPtBinTrigg%1.f_%1.fAssoc%1.f_%1.f",fTriggPtArray[i], fTriggPtArray[i+1],fAssocPtArray[j],fAssocPtArray[j+1]), ndphibins,dphimin,dphimax,nptbins, ptmin, ptmax);
    outputContainer->Add(fhDPhiSumPtBin[i][j]);
    
    fhDEtaSumPtBin[i][j] = new TH2F(Form("hDEtaSumPtBinTrigg%1.f_%1.fAssoc%1.f_%1.f",fTriggPtArray[i], fTriggPtArray[i+1],fAssocPtArray[j],fAssocPtArray[j+1]),
    Form("DEtaSumPtBinTrigg%1.f_%1.fAssoc%1.f_%1.f",fTriggPtArray[i], fTriggPtArray[i+1],fAssocPtArray[j],fAssocPtArray[j+1]), ndetabins,detamin,detamax,nptbins, ptmin, ptmax);
    outputContainer->Add(fhDEtaSumPtBin[i][j]);

    fhDPhiDEtaBin[i][j] = new TH2F(Form("hDPhiDEtaBinTrigg%1.f_%1.fAssoc%1.f_%1.f",fTriggPtArray[i], fTriggPtArray[i+1],fAssocPtArray[j],fAssocPtArray[j+1]),
    Form("hDPhiDEtaBinTrigg%1.f_%1.fAssoc%1.f_%1.f",fTriggPtArray[i], fTriggPtArray[i+1], fAssocPtArray[j],fAssocPtArray[j+1]), ndphibins, dphimin, dphimax, ndetabins,detamin,detamax);
    outputContainer->Add(fhDPhiDEtaBin[i][j]);

   fhMixDPhiDEtaBin[i][j] = new TH2F(Form("hMixDPhiDEtaBinTrigg%1.f_%1.fAssoc%1.f_%1.f",fTriggPtArray[i], fTriggPtArray[i+1],fAssocPtArray[j],fAssocPtArray[j+1]),
    Form("hMixDPhiDEtaBinTrigg%1.f_%1.fAssoc%1.f_%1.f",fTriggPtArray[i], fTriggPtArray[i+1],fAssocPtArray[j],fAssocPtArray[j+1]), ndphibins, dphimin, dphimax, ndetabins,detamin,detamax);
    outputContainer->Add(fhMixDPhiDEtaBin[i][j]);

   }/////end loop for associated bins
  }
  
  for(Int_t i=0; i<GetNTriggPtBins(); i++){
   fhDPhiAssocPtA[i] = new TH2F(Form("hDPhiAssocPtTrigg%1.f_%1.fAssoc%1.f_%1.fA",fTriggPtArray[i], fTriggPtArray[i+1],fptAssociatedBegin, 0.),
   Form("DPhiAssocPtTrigg%1.f_%1.fAssoc%1.f_%1.fA",fTriggPtArray[i], fTriggPtArray[i+1],fptAssociatedBegin, 0.),ndphibins,dphimin,dphimax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhDPhiAssocPtA[i]);

   fhDEtaAssocPtA[i] = new TH2F(Form("hDEtaAssocPtTrigg%1.f_%1.fAssoc%1.f_%1.fA",fTriggPtArray[i], fTriggPtArray[i+1],fptAssociatedBegin, 0.),
   Form("DEtaAssocPtTrigg%1.f_%1.fAssoc%1.f_%1.fA",fTriggPtArray[i], fTriggPtArray[i+1],fptAssociatedBegin, 0.), ndetabins,detamin,detamax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhDEtaAssocPtA[i]);

   fhMixDPhiAssocPtA[i] = new TH2F(Form("hMixDPhiAssocPtTrigg%1.f_%1.fAssoc%1.f_%1.fA",fTriggPtArray[i], fTriggPtArray[i+1],fptAssociatedBegin, 0.),
   Form("MixDPhiAssocPtTrigg%1.f_%1.fAssoc%1.f_%1.fA",fTriggPtArray[i], fTriggPtArray[i+1],fptAssociatedBegin, 0.),ndphibins,dphimin,dphimax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMixDPhiAssocPtA[i]);

   fhMixDEtaAssocPtA[i] = new TH2F(Form("hMixDEtaAssocPtTrigg%1.f_%1.fAssoc%1.f_%1.fA",fTriggPtArray[i], fTriggPtArray[i+1],fptAssociatedBegin, 0.),
   Form("MixDEtaAssocPtTrigg%1.f_%1.fAssoc%1.f_%1.fA",fTriggPtArray[i], fTriggPtArray[i+1],fptAssociatedBegin, 0.), ndetabins,detamin,detamax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMixDEtaAssocPtA[i]);
  }

   /////all
  fhDPhiTriggPtAssocPt = new TH2F(Form("hDPhiTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
  Form("DPhiTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),ndphibins,dphimin,dphimax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhDPhiTriggPtAssocPt);

  fhDEtaTriggPtAssocPt = new TH2F(Form("hDEtaTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
  Form("DEtaTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),ndetabins,detamin,detamax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhDEtaTriggPtAssocPt);

  fhAssocPtTriggPt = new TH2F(Form("hAssocPtTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
  Form("AssocPtTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nptbins, ptmin, ptmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhAssocPtTriggPt);

  fhpoutTriggPt = new TH2F(Form("hpoutTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
  Form("poutTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nptbins, ptmin, ptmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhpoutTriggPt);

  fhzTTriggPt = new TH2F(Form("hzTTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
  Form("zTTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nzTbins,zTmin,zTmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhzTTriggPt);

  fhxETriggPt = new TH2F(Form("hxETriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
  Form("xETriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nxEbins,xEmin,xEmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhxETriggPt);

  fhxELogTriggPt = new TH2F(Form("hxELogTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
  Form("xELogTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nxElogbins,xElogmin,xElogmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhxELogTriggPt);

  ////pout at Head Region in Away side
  fhAssocPtTriggPtHR = new TH2F(Form("hAssocPtTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
  Form("AssocPtTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nptbins, ptmin, ptmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhAssocPtTriggPtHR);    

  fhpoutTriggPtHR = new TH2F(Form("hpoutTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
  Form("poutTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nptbins, ptmin, ptmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhpoutTriggPtHR);

  fhzTTriggPtHR = new TH2F(Form("hzTTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
  Form("zTTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nzTbins,zTmin,zTmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhzTTriggPtHR);

  ////xE at Head Region in Away side
  fhxETriggPtHR = new TH2F(Form("hxETriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
  Form("xETriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nxEbins,xEmin,xEmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhxETriggPtHR);

  ////-log(xE) at Head Region in away side
  fhxELogTriggPtHR = new TH2F(Form("hxELogTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
  Form("xELogTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nxElogbins,xElogmin,xElogmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhxELogTriggPtHR);
  //////underlying observalbes
  ///////////Ue at Left side
  fhNUeAssocPtTriggPt = new TH2F(Form("hNUeAssocPtTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
  Form("NUeAssocPtTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nptbins, ptmin, ptmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhNUeAssocPtTriggPt);

  fhNUepoutTriggPt= new TH2F(Form("hNUepoutTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
  Form("NUepoutTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nptbins, ptmin, ptmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhNUepoutTriggPt);

  fhNUezTTriggPt = new TH2F(Form("hNUezTTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
  Form("NUezTTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nzTbins,zTmin,zTmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhNUezTTriggPt);

  fhNUexETriggPt = new TH2F(Form("hNUexETriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
  Form("NUexETriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nxEbins,xEmin,xEmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhNUexETriggPt);

  fhNUexELogTriggPt = new TH2F(Form("hNUexELogTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
  Form("NUexELogTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nxElogbins,xElogmin,xElogmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhNUexELogTriggPt);

  fhNUeDPhiDEta = new TH2F(Form("hNUeDPhiDEta%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
  Form("NUeDPhiDEtaTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), ndphibins, dphimin, dphimax, ndetabins,detamin,detamax);
  outputContainer->Add(fhNUeDPhiDEta);

  fhAUeAssocPtTriggPt = new TH2F(Form("hAUeAssocPtTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
  Form("AUeAssocPtTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nptbins, ptmin, ptmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhAUeAssocPtTriggPt);

  fhAUepoutTriggPt = new TH2F(Form("hAUepoutTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
  Form("AUepoutTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nptbins, ptmin, ptmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhAUepoutTriggPt);

  fhAUezTTriggPt = new TH2F(Form("hAUezTTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
  Form("AUezTTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nzTbins,zTmin,zTmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhAUezTTriggPt);

  fhAUexETriggPt = new TH2F(Form("hAUexETriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
  Form("AUexETriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nxEbins,xEmin,xEmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhAUexETriggPt);

  fhAUexELogTriggPt = new TH2F(Form("hAUexELogTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
  Form("AUexELogTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nxElogbins,xElogmin,xElogmax,nptbins, ptmin, ptmax);
  outputContainer->Add(fhAUexELogTriggPt);

  fhAUeDPhiDEta = new TH2F(Form("hAUeDPhiDEta%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
  Form("AUeDPhiDEtaTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), ndphibins, dphimin, dphimax, ndetabins,detamin,detamax);
  outputContainer->Add(fhAUeDPhiDEta);

  if(IsDataMC() && (kAnaMCTruthCorr || kAnaMCPrimaryCorr)){
   fhMCPtPhiLeading  = new TH2F ("hMCPtPhiLeading","p_{T}  and #phi distribution of leading particles", nptbins,ptmin,ptmax, nphibins, phimin, phimax);
   outputContainer->Add(fhMCPtPhiLeading);
   fhMCPtEtaLeading  = new TH2F ("hMCPtEtaLeading","p_{T} and #eta distribution of leading",nptbins,ptmin,ptmax, netabins,etamin,etamax);
   outputContainer->Add(fhMCPtEtaLeading);

   for(Int_t i=0; i< GetNAssocPtBins(); i++){
    fhMCDPhiTriggPtT[i] = new TH2F(Form("hMCDPhiTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fT",fptTriggerBegin, fptTriggerEnd, fAssocPtArray[i], fAssocPtArray[i+1]),
    Form("MCDPhiTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fT",fptTriggerBegin, fptTriggerEnd, fAssocPtArray[i], fAssocPtArray[i+1]),ndphibins,dphimin,dphimax,nptbins, ptmin, ptmax);
    outputContainer->Add(fhMCDPhiTriggPtT[i]);

    fhMCDEtaTriggPtT[i] = new TH2F(Form("hMCDEtaTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fT",fptTriggerBegin, fptTriggerEnd, fAssocPtArray[i], fAssocPtArray[i+1]),
    Form("MCDEtaTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fT",fptTriggerBegin, fptTriggerEnd, fAssocPtArray[i], fAssocPtArray[i+1]), ndetabins,detamin,detamax,nptbins, ptmin, ptmax);
    outputContainer->Add(fhMCDEtaTriggPtT[i]);
   }

   for(Int_t i = 0 ; i < GetNTriggPtBins() ; i++){
    for(Int_t j=0; j<GetNAssocPtBins(); j++){
     fhMCDPhiSumPtBin[i][j] = new TH2F(Form("hMCDPhiSumPtBinTrigg%1.f_%1.fAssoc%1.f_%1.f",fTriggPtArray[i], fTriggPtArray[i+1],fAssocPtArray[j],fAssocPtArray[j+1]),
     Form("MCDPhiSumPtBinTrigg%1.f_%1.fAssoc%1.f_%1.f",fTriggPtArray[i], fTriggPtArray[i+1],fAssocPtArray[j],fAssocPtArray[j+1]), ndphibins,dphimin,dphimax,nptbins, ptmin, ptmax);
     outputContainer->Add(fhMCDPhiSumPtBin[i][j]);
   
     fhMCDEtaSumPtBin[i][j] = new TH2F(Form("hMCDEtaSumPtBinTrigg%1.f_%1.fAssoc%1.f_%1.f",fTriggPtArray[i], fTriggPtArray[i+1],fAssocPtArray[j],fAssocPtArray[j+1]),
     Form("MCDEtaSumPtBinTrigg%1.f_%1.fAssoc%1.f_%1.f",fTriggPtArray[i], fTriggPtArray[i+1],fAssocPtArray[j],fAssocPtArray[j+1]), ndetabins,detamin,detamax,nptbins, ptmin, ptmax);
     outputContainer->Add(fhMCDEtaSumPtBin[i][j]);

     fhMCDPhiDEtaBin[i][j] = new TH2F(Form("hMCDPhiDEtaBinTrigg%1.f_%1.fAssoc%1.f_%1.f",fTriggPtArray[i], fTriggPtArray[i+1],fAssocPtArray[j],fAssocPtArray[j+1]),
     Form("hMCDPhiDEtaBinTrigg%1.f_%1.fAssoc%1.f_%1.f",fTriggPtArray[i], fTriggPtArray[i+1], fAssocPtArray[j],fAssocPtArray[j+1]), ndphibins, dphimin, dphimax, ndetabins,detamin,detamax);
     outputContainer->Add(fhMCDPhiDEtaBin[i][j]);
    }/////end loop for associated bins
   }
  
   for(Int_t i=0; i<GetNTriggPtBins(); i++){
    fhMCDPhiAssocPtA[i] = new TH2F(Form("hMCDPhiAssocPtTrigg%1.f_%1.fAssoc%1.f_%1.fA",fTriggPtArray[i], fTriggPtArray[i+1],fptAssociatedBegin, 0.),
    Form("MCDPhiAssocPtTrigg%1.f_%1.fAssoc%1.f_%1.fA",fTriggPtArray[i], fTriggPtArray[i+1],fptAssociatedBegin, 0.),ndphibins,dphimin,dphimax,nptbins, ptmin, ptmax);
    outputContainer->Add(fhMCDPhiAssocPtA[i]);

    fhMCDEtaAssocPtA[i] = new TH2F(Form("hMCDEtaAssocPtTrigg%1.f_%1.fAssoc%1.f_%1.fA",fTriggPtArray[i], fTriggPtArray[i+1],fptAssociatedBegin, 0.),
    Form("MCDEtaAssocPtTrigg%1.f_%1.fAssoc%1.f_%1.fA",fTriggPtArray[i], fTriggPtArray[i+1],fptAssociatedBegin, 0.), ndetabins,detamin,detamax,nptbins, ptmin, ptmax);
    outputContainer->Add(fhMCDEtaAssocPtA[i]);
   }

   /////all
   fhMCDPhiAssocPt15T = new TH2F(Form("hMCDPhiTriggPtTrigg%1.f_%1.fAssoc1.0_5.0",fptTriggerBegin, fptTriggerEnd),
   Form("MCDPhiTriggPtTrigg%1.f_%1.fAssoc1.0_5.0",fptTriggerBegin, fptTriggerEnd),ndphibins,dphimin,dphimax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCDPhiAssocPt15T);

   fhMCDEtaAssocPt15T = new TH2F(Form("hMCDEtaTriggPtTrigg%1.f_%1.fAssoc1.0_5.0",fptTriggerBegin, fptTriggerEnd),
   Form("MCDEtaTriggPtTrigg%1.f_%1.fAssoc1.0_5.0",fptTriggerBegin, fptTriggerEnd),ndetabins,detamin,detamax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCDEtaAssocPt15T);

   fhMCAssocPtTriggPt = new TH2F(Form("hMCAssocPtTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
   Form("MCAssocPtTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nptbins, ptmin, ptmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCAssocPtTriggPt);

   fhMCpoutTriggPt = new TH2F(Form("hMCpoutTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
   Form("MCpoutTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nptbins, ptmin, ptmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCpoutTriggPt);

   fhMCzTTriggPt = new TH2F(Form("hMCzTTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
   Form("MCzTTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nzTbins,zTmin,zTmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCzTTriggPt);

   fhMCxETriggPt = new TH2F(Form("hMCxETriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
   Form("MCxETriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nxEbins,xEmin,xEmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCxETriggPt);

   fhMCxELogTriggPt = new TH2F(Form("hMCxELogTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
   Form("MCxELogTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nxElogbins,xElogmin,xElogmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCxELogTriggPt);

   ////pout at Head Region in Away side
   fhMCAssocPtTriggPtHR = new TH2F(Form("hMCAssocPtTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
   Form("MCAssocPtTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nptbins, ptmin, ptmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCAssocPtTriggPtHR);    

   fhMCpoutTriggPtHR = new TH2F(Form("hMCpoutTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
   Form("MCpoutTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nptbins, ptmin, ptmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCpoutTriggPtHR);

   fhMCzTTriggPtHR = new TH2F(Form("hMCzTTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
   Form("MCzTTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nzTbins,zTmin,zTmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCzTTriggPtHR);

   ////xE at Head Region in Away side
   fhMCxETriggPtHR = new TH2F(Form("hMCxETriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
   Form("MCxETriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nxEbins,xEmin,xEmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCxETriggPtHR);

   ////-log(xE) at Head Region in away side
   fhMCxELogTriggPtHR = new TH2F(Form("hMCxELogTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
   Form("MCxELogTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.fHR",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nxElogbins,xElogmin,xElogmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCxELogTriggPtHR);

   //////underlying observalbes
   ///////////Ue at Left side
   fhMCNUeAssocPtTriggPt = new TH2F(Form("hMCNUeAssocPtTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
   Form("MCNUeAssocPtTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nptbins, ptmin, ptmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCNUeAssocPtTriggPt);
 
   fhMCNUepoutTriggPt= new TH2F(Form("hMCNUepoutTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
   Form("MCNUepoutTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nptbins, ptmin, ptmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCNUepoutTriggPt);
 
   fhMCNUezTTriggPt = new TH2F(Form("hMCNUezTTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
   Form("MCNUezTTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nzTbins,zTmin,zTmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCNUezTTriggPt);
 
   fhMCNUexETriggPt = new TH2F(Form("hMCNUexETriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
   Form("MCNUexETriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nxEbins,xEmin,xEmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCNUexETriggPt);

   fhMCNUexELogTriggPt = new TH2F(Form("hMCNUexELogTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
   Form("MCNUexELogTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nxElogbins,xElogmin,xElogmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCNUexELogTriggPt);
 
   fhMCAUeAssocPtTriggPt = new TH2F(Form("hMCAUeAssocPtTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
   Form("MCAUeAssocPtTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nptbins, ptmin, ptmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCAUeAssocPtTriggPt);

   fhMCAUepoutTriggPt = new TH2F(Form("hMCAUepoutTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
   Form("MCAUepoutTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nptbins, ptmin, ptmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCAUepoutTriggPt);

   fhMCAUezTTriggPt = new TH2F(Form("hMCAUezTTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
   Form("MCAUezTTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.), nzTbins,zTmin,zTmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCAUezTTriggPt);

   fhMCAUexETriggPt = new TH2F(Form("hMCAUexETriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.),
   Form("MCAUexETriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nxEbins,xEmin,xEmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCAUexETriggPt);

   fhMCAUexELogTriggPt = new TH2F(Form("hMCAUexELogTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd,fptAssociatedBegin, 0.),
   Form("MCAUexELogTriggPtTrigg%1.f_%1.fAssoc%1.f_%1.f",fptTriggerBegin, fptTriggerEnd, fptAssociatedBegin, 0.), nxElogbins,xElogmin,xElogmax,nptbins, ptmin, ptmax);
   outputContainer->Add(fhMCAUexELogTriggPt);

  }


  PostData(1,outputContainer);

}

//________________________________________________________________________
void AliAnalysisTaskEMCALCaloTrackCorr::UserExec(Option_t *) 
{
  // Execute analysis for current event
  fManager = AliAnalysisManager::GetAnalysisManager();
  fInputHandler = dynamic_cast<AliInputEventHandler*>(fManager->GetInputEventHandler());

  if(!fInputHandler) return;

  if(!FillInputEvent()) return;  //select good event

  fEventPlane = (dynamic_cast<AliVEvent*>(InputEvent()))->GetEventplane();

  UInt_t isSelectedReal = 0;
  if(!kEventTriggerAtSE) isSelectedReal = fInputHandler->IsEventSelected() & fEventTriggerMaks;

  const AliVVertex *fEventVertex = fEvent->GetPrimaryVertex();
  Double_t vtx[3];
  vtx[0] = fEventVertex->GetX();
  vtx[1] = fEventVertex->GetY();
  vtx[2] = fEventVertex->GetZ();

  if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()){
   fMCEvent = MCEvent();
   fStack = fMCEvent->Stack();
  }

  fnEvents++;
 
  if(fPhotonEvent)
   fPhotonEvent->Clear();
  else
   fPhotonEvent = new TClonesArray("AliCaloTrackParticle", InputEvent()->GetNumberOfCaloClusters());
   nPhotonsEMCAL=0;
  if(fPhotonPairEvent)
   fPhotonPairEvent->Clear();
  else
   fPhotonPairEvent = new TClonesArray("AliCaloTrackParticle", InputEvent()->GetNumberOfCaloClusters());

  if(fCTSEvent)
    fCTSEvent->Clear();
  else
    fCTSEvent = new TClonesArray("AliCaloTrackParticle", InputEvent()->GetNumberOfTracks());
    nTracksCTS = 0;

  if(kEventTriggerAtSE){
   FillInputTrack();
   FillInputPhoton();
   if(kDoMesonFill)FillInputMeson();
  }
  else {
    FillInputTrack();
    if(isSelectedReal){
     FillInputPhoton();
     if(kDoMesonFill)FillInputMeson();
   }
  }
  if(nTracksCTS==0) return;


  UInt_t isSelectedMix = 0;
  isSelectedMix =((fInputHandler->IsEventSelected() & AliVEvent::kAnyINT) ||
                  (fInputHandler->IsEventSelected() & AliVEvent::kMB) ||
                  (fInputHandler->IsEventSelected() & AliVEvent::kCentral) ||
                  (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral));

  Double_t vtx5[3] ={ fEvent->GetPrimaryVertex()->GetX(),
                      fEvent->GetPrimaryVertex()->GetY(),
                      fEvent->GetPrimaryVertex()->GetZ()
                    };
  Int_t zvtx      =  -999;
  Int_t centr     =  -999;
  Int_t replane   =  -999;
  Int_t nTracks = fCTSEvent->GetEntriesFast(); 
 
  if(vtx5[2]>=-10. && vtx5[2]<-8)
   zvtx=0;
   else if(vtx5[2]>=-8. && vtx5[2]<-6)
    zvtx=1;
    else if(vtx5[2]>=-6. && vtx5[2]<-4)
     zvtx=2;
     else if(vtx5[2]>=-4. && vtx5[2]<-2)
      zvtx=3;
      else if(vtx5[2]>=-2. && vtx5[2]<-0)
       zvtx=4;
       else if(vtx5[2]>=0. && vtx5[2]<=2)
        zvtx=5;
        else if(vtx5[2]>2. && vtx5[2]<=4)
         zvtx=6;
         else if(vtx5[2]>4. && vtx5[2]<=6)
           zvtx=7;
           else if(vtx5[2]>6. && vtx5[2]<=8)
            zvtx=8;
            else if(vtx5[2]>8. && vtx5[2]<=10)
             zvtx=9;

  if(kDoTrackMultBins){  // for pp analysis
   if(nTracks<=5)
    centr=8;
    else if(nTracks<=10)
     centr=7;
     else if(nTracks<=15)
      centr=6;
      else if(nTracks<=20)
       centr=5;
       else if(nTracks<=30)
        centr=4;
        else if(nTracks<=40)
         centr=3;
         else if(nTracks<=55)
          centr=2;
          else if(nTracks<=70)
           centr=1;
           else centr=0;

   replane = 0;
  }
  else{///for PbPb
   Float_t fCentralityPerBin = (fCentralityBinMax -fCentralityBinMin)/fNCentralityBins;
   Float_t fcentrality = 0; 
  
   if(fDataType == "ESD")fCentrality->GetCentralityPercentile(fCentralityClass);
   else if(fDataType == "AOD") fcentrality = ((AliAODHeader*)fEvent->GetHeader())->GetCentrality();
   else fcentrality =0;   

   if(fcentrality <=(fCentralityBinMin+1*fCentralityPerBin))
    centr=9;
    else if(fcentrality <=(fCentralityBinMin+2*fCentralityPerBin))
     centr=8;
     else if(fcentrality <=(fCentralityBinMin+3*fCentralityPerBin))
      centr=7;
      else if(fcentrality <=(fCentralityBinMin+4*fCentralityPerBin))
       centr=6;
       else if(fcentrality <=(fCentralityBinMin+5*fCentralityPerBin))
        centr=5;
        else if(fcentrality <=(fCentralityBinMin+6*fCentralityPerBin))
         centr=4;
         else if(fcentrality <=(fCentralityBinMin+7*fCentralityPerBin))
          centr=3;
          else if(fcentrality <=(fCentralityBinMin+8*fCentralityPerBin))
           centr=2;
           else if(fcentrality <=(fCentralityBinMin+9*fCentralityPerBin))
            centr=1;
            else if(fcentrality <=(fCentralityBinMin+10*fCentralityPerBin))
             centr=0;

   if(fEventPlane){
    Float_t feventplane = fEventPlane->GetEventplane(fEventPlaneMethod);
    if(feventplane<0)feventplane+=TMath::Pi();
    if(feventplane>TMath::Pi())feventplane-=TMath::Pi();

    replane = Int_t((fNEventPlaneBins*feventplane)/TMath::Pi());
    if(replane>(fNEventPlaneBins -1)) replane = fNEventPlaneBins -1;
   }
   else{
    replane = 0;
   }
  }

  if(!fListMixEvents[zvtx][centr][replane]) fListMixEvents[zvtx][centr][replane]=new TList() ;
  TList *pool = fListMixEvents[zvtx][centr][replane];

  TClonesArray *fAnaAODParticle = new TClonesArray();
  if(kDoPhotonCorrAna && !kDoMesonCorrAna){
   fAnaAODParticle = fPhotonEvent;
  }  
  if(kDoMesonCorrAna && !kDoPhotonCorrAna){
   fAnaAODParticle = fPhotonPairEvent;
  }

  Int_t naod = fAnaAODParticle->GetEntriesFast();
  for(Int_t iaod=0; iaod<naod; iaod++){
   AliCaloTrackParticle *aodParticle=(AliCaloTrackParticle*)fAnaAODParticle->At(iaod);
   if(!aodParticle) continue;
   if(kDoIsolatedAna) {
    if(aodParticle->IsIsolated()) 
     printf("Inclusive correlation Isolated  pt=%f, id=%d\n",aodParticle->Pt(),iaod);
     if(GetAnaTypeInIsolated() == "Iso" && !aodParticle->IsIsolated()) continue;
     if(GetAnaTypeInIsolated() == "NoIso" && aodParticle->IsIsolated())continue;
   }
   
   if(kMakeAbsoluteLeading  || kMakeNearSideLeading){
    if(!aodParticle->IsLeading())continue;
   }   

   Double_t ptTrigg0  = aodParticle->Pt();
   Double_t phiTrigg0 = aodParticle->Phi();
   if(phiTrigg0<0.) phiTrigg0+=TMath::TwoPi();
   Double_t etaTrigg0 = aodParticle->Eta();
   
   if(MakeChargedCorrelation(iaod, ptTrigg0, phiTrigg0, etaTrigg0)){
    if(kDoMixEventsAna) MakeChargedMixCorrelation(ptTrigg0, phiTrigg0, etaTrigg0, pool);
   }

   aodParticle = 0;

  } 

  if(isSelectedMix && kDoMixEventsAna){
   if(fCTSEvent->GetEntriesFast()>0){
    pool->AddFirst(fCTSEvent);
    fCTSEvent=0;
    if(pool->GetSize()>nMixedEvents){//Remove redundant events
     TClonesArray * tmp = static_cast<TClonesArray*>(pool->Last()) ;
     pool->RemoveLast() ;
     delete tmp ;
    }
   }

  }

  PostData(1, outputContainer);
  fEventAnalyized++;

}

//=============================================
void AliAnalysisTaskEMCALCaloTrackCorr::InitParameters()
{
  fHistoPtBins  = 300 ;  fHistoPtMax  = 60 ;            fHistoPtMin  = 0.  ;
  fHistoPhiBins = 360 ;  fHistoPhiMax = TMath::TwoPi(); fHistoPhiMin = 0.  ;
  fHistoEtaBins = 180 ;  fHistoEtaMax = 0.9;            fHistoEtaMin = -0.9;

  fSetConeR          = 0.4 ;
  fSetPtThreshold    = 0.5 ;
  fSetSumPtThreshold = 0.5 ;
  fSetPtFraction     = 0.1 ;
  fParticlesInCone   = kIsolatedNeutralAndCharged;
  fICMethod          = kSumPtFracationIC ;

  fCentralityClass  = "V0M";
  fCentralityBinMin = fCentralityBinMax=-1;
  fEventPlaneMethod = "Q";
  fEventTriggerMaks = AliVEvent::kAny;
  fNCentralityBins  = 5;
  fNEventPlaneBins  = 5;

  kMC = kFALSE ;
  kDoMixEventsAna    = kTRUE  ;
  kDoPhotonCorrAna   = kTRUE  ;
  kDoAsymmetryCut    = kTRUE  ;
  kDoSelectHybridTracks = kFALSE;
  kDoMesonFill        = kTRUE ;
  kNeutralMesonHistos = kTRUE  ;
  kDoMesonCorrAna     = kFALSE ;
  kDoIsolatedAna      = kFALSE ;
  kDoTrackMultBins   = kTRUE  ;
  kUELeftRight       = kFALSE ;
  kUENearAway        = kTRUE  ;
  kDecayPhotonCorr   = kFALSE ;
  kAnaMCTruthCorr    = kTRUE  ;
  kAnaMCPrimaryCorr  = kFALSE ;
  kAnaPi0Prim        = kFALSE ;
  kAnaEtaPrim        = kFALSE ;
  kAnaPhotonPrim     = kFALSE ;
  kMakeAbsoluteLeading = kFALSE ;
  kMakeNearSideLeading = kFALSE ;
  kTwoTracksCorr     = kFALSE ;
  kPhotonInAcceptance  = kTRUE;
  kAnaDecayMapping   = kTRUE  ;
  //fCheckFidCut       = kTRUE  ;
  kEventTriggerAtSE  = kTRUE  ;
  kPhotonPairTimeCut = kFALSE ;
  kDoPhotonIDCut = kTRUE;

  fMinNCells = 1 ;
  fMinE = 0.5;
  fMinDistBad = 2. ;
  
  fL0CutMin = 0.1;
  fL0CutMax = 0.27;
  fTimeCutMin = -1000;
  fTimeCutMax = 1000;
  fPhotonPairTimeCut = 100;
  fEMCALDPhiCut = 0.03;
  fEMCALDEtaCut = 0.025;
  fZVertexCut = 10.;

  fDebug = -1 ;
  fAnaMesonType = "Pi0";
  fAsymmetryCut = 0.7;
  fDataType = "ESD" ;
  fTrackFilterMask = 128;


  fAnaTypeInIsolated   = "Iso"; 
  nMixedEvents     = 200;
  fptTriggerBegin = 5.;
  fptTriggerEnd   = 25.;
  fptAssociatedBegin = 1.;

  fDeltaPhiMinCut  = TMath::Pi()/2.;
  fDeltaPhiMaxCut  = 3*TMath::Pi()/2.;
  fDeltaPhiHRSize  = TMath::Pi()/5.;
  fUeDeltaPhiSize  = 0.2;
  fUeDeltaPhiFix   = 1.3;

  fESDtrackCuts =  AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();

}
//=====================================
Bool_t AliAnalysisTaskEMCALCaloTrackCorr::FillInputEvent()
{
   Float_t  fcentrality = 0;

   if(GetDataType() == "ESD"){
    fEvent=dynamic_cast<AliESDEvent*>(InputEvent());
    if(!fEvent) {
     Printf("ERROR: Could not retrieve event");
     return kFALSE;
    }

    AliESDEvent * fESDEvent = dynamic_cast<AliESDEvent*> (InputEvent());
    if(!fESDEvent) return kFALSE; 
    if(fESDEvent->GetPrimaryVertexTracks()->GetNContributors() <= 0){
     return kFALSE;
    }

    if(!IsDataMC()){;
     if((TMath::Abs(fEvent->GetPrimaryVertex()->GetX()) < 1.e-6) &&
      (TMath::Abs(fEvent->GetPrimaryVertex()->GetY()) < 1.e-6) &&
      (TMath::Abs(fEvent->GetPrimaryVertex()->GetZ()) < 1.e-6))
     return kFALSE;
    }

    fCentrality = fEvent->GetCentrality();
    if(fCentrality)
    fcentrality=fCentrality->GetCentralityPercentile(fCentralityClass);

    if((fcentrality < fCentralityBinMin) || (fcentrality>fCentralityBinMax)) return kFALSE; 

   }
   else if(GetDataType() == "AOD"){
    fEvent=dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fEvent) {
     Printf("ERROR: Could not retrieve event");
     return kFALSE;
    }
    
    fcentrality = ((AliAODHeader*)fEvent->GetHeader())->GetCentrality();
    if((fcentrality < fCentralityBinMin) || (fcentrality>fCentralityBinMax)) return kFALSE;

   }
   else return kFALSE;

   if(TMath::Abs(fEvent->GetPrimaryVertex()->GetZ())>=fZVertexCut) {
    return kFALSE; ///cut for primary vertex
   }
  
   return kTRUE;

}

//======================================
void AliAnalysisTaskEMCALCaloTrackCorr::FillInputTrack()
{
  //select good track!

  for(Int_t itrack=0; itrack<fEvent->GetNumberOfTracks(); itrack++){
   AliVTrack* track = (AliVTrack*)fEvent->GetTrack(itrack);

   if(GetDataType() == "ESD"){
    AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*>(track);
    if(!esdTrack) {
    AliError(Form("Couldn't get ESD track %d\n", itrack));
    continue;
    }

    Bool_t trkOK = fESDtrackCuts->AcceptTrack(esdTrack);
    if (!trkOK) continue;
   
   }
   else if(GetDataType() == "AOD"){
    AliAODTrack *aodTrack = dynamic_cast <AliAODTrack*>(track);
    if(!aodTrack) {
    AliError(Form("Couldn't get AOD track %d\n", itrack));
    continue;
    }

    if(kDoSelectHybridTracks){
     if(!aodTrack->IsHybridGlobalConstrainedGlobal())  continue;
    }
    else{
     if(aodTrack->TestFilterBit(fTrackFilterMask)==kFALSE) continue;
    }

    if(aodTrack->GetType()!= AliAODTrack::kPrimary) continue; 


   }
  
   new((*fCTSEvent)[nTracksCTS])AliCaloTrackParticle(track->Px(),track->Py(),track->Pz(),0);
   AliCaloTrackParticle * tr = (AliCaloTrackParticle*)fCTSEvent->At(nTracksCTS);
   tr->SetChargedSign(track->Charge());
   nTracksCTS++;

  }///end loop for tracks

  if(nTracksCTS==0) return;

  fhNtracksAll->Fill(nTracksCTS);

  if(fInputHandler->IsEventSelected( ) & AliVEvent::kAnyINT){
   fhNtracksAnyINT->Fill(nTracksCTS);
  }
  else if(fInputHandler->IsEventSelected( ) & AliVEvent::kCentral){
   fhNtracksCentral->Fill(nTracksCTS);
  }
  else if(fInputHandler->IsEventSelected( ) & AliVEvent::kSemiCentral){
   fhNtracksSemiCentral->Fill(nTracksCTS);
  }
  else if(fInputHandler->IsEventSelected( ) & AliVEvent::kEMC7){
   fhNtracksEMC7->Fill(nTracksCTS);
  }
  else{
   fhNtracksOtherTirgger->Fill(nTracksCTS);
  }

}

//=========================================
void AliAnalysisTaskEMCALCaloTrackCorr::FillInputPhoton()
{
  Int_t nClusters = 0;
  Double_t tmpPt = 0;
  Int_t    tmpId = -1;

  nClusters=fEvent->GetNumberOfCaloClusters();
  if(nClusters==0) return;

  for (Int_t j=0; j<nClusters; j++) {
   AliVCluster *clusEMCAL = fEvent->GetCaloCluster(j);
   AliVCaloCells *cellEMCAL = fEvent->GetEMCALCells();
   if(!clusEMCAL->IsEMCAL())  continue;
   if(clusEMCAL->E() < fMinE) continue;
   if(clusEMCAL->GetNCells() < fMinNCells) continue;
   if(clusEMCAL->GetDistanceToBadChannel() < fMinDistBad) continue;
   if(!fEMCALRecU->IsGoodCluster(clusEMCAL, fEMCALGeom, cellEMCAL)) continue;

   TLorentzVector ph;
   Double_t vtx[3];
   clusEMCAL ->GetMomentum(ph, vtx);
   Int_t cellAbsId=-1;
   fEMCALGeom->GetAbsCellIdFromEtaPhi(ph.Eta(), ph.Phi(), cellAbsId);
   Int_t iSM=-1, iTower=-1, Iphi=-1, Ieta=-1;
   fEMCALGeom->GetCellIndex(cellAbsId, iSM, iTower, Iphi, Ieta);
   Int_t iPhi=-1, iEta=-1;
   fEMCALGeom->GetCellPhiEtaIndexInSModule(iSM, iTower, Iphi, Ieta, iPhi, iEta); 
   
   if(ph.Pt()>tmpPt){
    tmpPt = ph.Pt();
    tmpId = clusEMCAL->GetID();
   }

   ////put selected clusters into AOD (AliCaloTrackParticle)
   new((*fPhotonEvent)[nPhotonsEMCAL])AliCaloTrackParticle(ph.Px(), ph.Py(), ph.Pz(), ph.E());
   AliCaloTrackParticle * phcluster = (AliCaloTrackParticle*)fPhotonEvent->At(nPhotonsEMCAL);
   phcluster->SetModule(iSM);
   phcluster->SetNCells(clusEMCAL->GetNCells());
   phcluster->SetLambdas(clusEMCAL->GetM02(),clusEMCAL->GetM20());
   phcluster->SetDistBad(clusEMCAL->GetDistanceToBadChannel());
   phcluster->SetTOF(clusEMCAL->GetTOF()*1e9);
   phcluster->SetClusterID(clusEMCAL->GetID());
   phcluster->SetAODClusterID(nPhotonsEMCAL);
   phcluster->SetSSABit(kFALSE);
   phcluster->SetTOFBit(kFALSE);
   phcluster->SetTrackMatchedBit(kFALSE);

   Double_t tof = clusEMCAL->GetTOF()*1e9;
   if(tof >= fTimeCutMin && tof <= fTimeCutMax)
    phcluster->SetTOFBit(kTRUE);

   if(clusEMCAL->GetM02() > fL0CutMin && clusEMCAL->GetM02() < fL0CutMax)
    phcluster->SetSSABit(kTRUE);

   Int_t nMatches = clusEMCAL->GetNTracksMatched();
   if(nMatches>0){
    Float_t dZ  = 2000.;
    Float_t dR  = 2000.;
    fEMCALRecU->GetMatchedResiduals(clusEMCAL->GetID(), dZ, dR);
    if(TMath::Abs(dR) < fEMCALDPhiCut && TMath::Abs(dZ) < fEMCALDEtaCut){
     phcluster->SetTrackMatchedBit(kTRUE);
    }
   }
   
   nPhotonsEMCAL++;

  }///end loop for clusters

  for(Int_t j=0; j< nPhotonsEMCAL; j++){
   AliCaloTrackParticle * phCan=(AliCaloTrackParticle*)fPhotonEvent->At(j);
   Double_t phPt  = phCan->Pt();
   Double_t phPhi = phCan->Phi();
   if(phPhi<0) phPhi+=TMath::TwoPi();
   Double_t phEta = phCan->Eta();

   phCan->SetIsLeading(kFALSE);
   if(tmpId == phCan->GetClusterID())
   phCan->SetIsLeading(kTRUE);
  
   phCan->SetIsolated(kFALSE);
   phCan->SetIsolated(IsolatedPhoton(fPhotonEvent, fCTSEvent, j, phPt, phPhi, phEta));

   if(kDoPhotonIDCut && !phCan->IsPIDOK(8)) continue;
    
   fhPhotonE->Fill(phCan->E());
   fhPhotonPtPhi->Fill(phCan->Pt(), phCan->Phi());
   fhPhotonPtEta->Fill(phCan->Pt(), phCan->Eta());
   fhPhotonPhiEta->Fill(phCan->Phi(), phCan->Eta());
  
  }


}

//===================================
void AliAnalysisTaskEMCALCaloTrackCorr::FillInputMeson()
{  ////do pi0 and eta meson tag
  TLorentzVector ph12;
  Int_t nPhotonPairEMCAL = 0;

  for (Int_t i=0; i<nPhotonsEMCAL-1; i++) {
   AliCaloTrackParticle * ph1=(AliCaloTrackParticle*)fPhotonEvent->At(i);
   if(kDoPhotonIDCut && !ph1->IsPIDOK(8)) continue;
   for (Int_t j=i+1; j<nPhotonsEMCAL; j++) {
    AliCaloTrackParticle * ph2=(AliCaloTrackParticle*)fPhotonEvent->At(j);
    if(kDoPhotonIDCut && !ph2->IsPIDOK(8)) continue;
    ph12  = *ph1  + *ph2;
    Double_t phi=ph12.Phi();
    if(phi < 0) phi+=TMath::TwoPi();
    Float_t angle   = ph1->Angle(ph2->Vect());
    Float_t asy     = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy())); 
    Double_t deltaTime = TMath::Abs(ph1->GetTOF()-ph2->GetTOF());
    Int_t    deltaModule = TMath::Abs(ph1->GetModule()-ph2->GetModule());  

    new((*fPhotonPairEvent)[nPhotonPairEMCAL])AliCaloTrackParticle(ph12.Px(), ph12.Py(), ph12.Pz(), ph12.E());
    AliCaloTrackParticle * phPair = (AliCaloTrackParticle*)fPhotonPairEvent->At(nPhotonPairEMCAL);
    phPair->SetPhotonPairAsy(asy); 
    phPair->SetPhotonPairAngle(angle);
    phPair->SetPhotonPairDTime(deltaTime);
    phPair->SetPhotonPairDModule(deltaModule);
    phPair->SetPhotonPairID(ph1->GetClusterID(), ph2->GetClusterID());   
    phPair->SetAODPhotonPairID(ph1->GetAODClusterID(), ph2->GetAODClusterID());   
    nPhotonPairEMCAL++;
   }//ph2
  }//ph1

  Double_t tmpMesonPt = 0;
  Int_t    tmpMesonID = -1;
 
  for(Int_t k=0; k<nPhotonPairEMCAL; k++){
   AliCaloTrackParticle * phPairCan=(AliCaloTrackParticle*)fPhotonPairEvent->At(k);
   if(kPhotonPairTimeCut && phPairCan->GetPhotonPairDTime() > fPhotonPairTimeCut) continue;
   if(SelectPair(phPairCan)){
    AliCaloTrackParticle *photon1 = (AliCaloTrackParticle*)fPhotonEvent->At(phPairCan->GetAODPhotonPairID(0));
    AliCaloTrackParticle *photon2 = (AliCaloTrackParticle*)fPhotonEvent->At(phPairCan->GetAODPhotonPairID(1));

    fhMesonE->Fill(phPairCan->E());
    fhMesonPtPhi->Fill(phPairCan->Pt(), phPairCan->Phi());
    fhMesonPtEta->Fill(phPairCan->Pt(), phPairCan->Eta());
    fhMesonPhiEta->Fill(phPairCan->Phi(), phPairCan->Eta());

    if(phPairCan->Pt() > tmpMesonPt){
     tmpMesonPt = phPairCan->Pt();
     tmpMesonID = k;
    }

    if(kAnaDecayMapping){
     Double_t deltaphi1 = TMath::Abs(phPairCan->Phi()-photon1->Phi());
     if(deltaphi1 > TMath::TwoPi()) deltaphi1-= TMath::TwoPi();

     Double_t deltaphi2 = TMath::Abs(phPairCan->Phi()-photon2->Phi());
     if(deltaphi2 > TMath::TwoPi()) deltaphi2-= TMath::TwoPi();

     Double_t deltaphi3 = TMath::Abs(photon1->Phi()-photon2->Phi());
     if(deltaphi3 > TMath::TwoPi()) deltaphi3-= TMath::TwoPi();

     fhPi0DecayPhoton1->Fill(phPairCan->Pt(),photon1->Pt());
     fhPi0DecayPhoton1Dphi->Fill(phPairCan->Pt(),deltaphi1);
     fhDecayPhoton1Pi0Dphi->Fill(photon1->Pt(),deltaphi1);
     fhPi0DecayPhoton2->Fill(phPairCan->Pt(),photon2->Pt());
     fhPi0DecayPhoton2Dphi->Fill(phPairCan->Pt(),deltaphi2);
     fhDecayPhoton2Pi0Dphi->Fill(photon2->Pt(),deltaphi2);
     fhDecayPhoton1Photon2->Fill(photon1->Pt(),photon2->Pt());
     fhDecayPhoton1Photon2Dphi->Fill(photon1->Pt(),deltaphi3);
     fhDecayPhoton2Photon1Dphi->Fill(photon2->Pt(),deltaphi3);
    }
   }
  }//end loop for mesons

  for(Int_t kmeson=0; kmeson<nPhotonPairEMCAL; kmeson++){
   AliCaloTrackParticle * phPairCan=(AliCaloTrackParticle*)fPhotonPairEvent->At(kmeson);
   phPairCan->SetIsLeading(kFALSE);
   if(kmeson == tmpMesonID) phPairCan->SetIsLeading(kTRUE);
  }

}

//======================================
Bool_t AliAnalysisTaskEMCALCaloTrackCorr::SelectPair(AliCaloTrackParticle *mesonCandidate)
{
  Double_t fInvMassMaxCutParam[3]={0.};
  if(fAnaMesonType.Contains("Pi0")){
   fInvMassMaxCutParam[0] = 0.0   ;
   fInvMassMaxCutParam[1] =-7.e-5 ;
   fInvMassMaxCutParam[2] = 8.e-5 ;
  }
  else {
   fInvMassMaxCutParam[0] = 0.00 ;
   fInvMassMaxCutParam[1] = 0.00 ;
   fInvMassMaxCutParam[2] = 0.00 ;
  }

  if(fAnaMesonType.Contains("Pi0")){
   fInvMassMinCut = 0.12;
   fInvMassMaxCut = 0.15;
   fLeftBandMinCut = 0;
   fLeftBandMaxCut = -1;
   fRightBandMinCut = 0.17;
   fRightBandMaxCut = 0.20;
  }

  if(fAnaMesonType.Contains("Eta")){
   fInvMassMinCut = 0.52;
   fInvMassMaxCut = 0.58;
   fLeftBandMinCut = 0.38;
   fLeftBandMaxCut = 0.43;
   fRightBandMinCut = 0.65;
   fRightBandMaxCut = 0.70;
  }

  Double_t phi = mesonCandidate->Phi();
  if(phi < 0) phi+=TMath::TwoPi();
  Double_t eta = mesonCandidate->Eta();
  Double_t invmass = mesonCandidate->M();
  Double_t angle = mesonCandidate->GetPhotonPairAngle();
  Double_t e = mesonCandidate->Energy();
  Double_t asy = mesonCandidate->GetPhotonPairAsy();

  if(kNeutralMesonHistos){
   fhAnglePairNoCut  ->Fill(e,angle);
   fhInvMassPairNoCut->Fill(e,invmass);
   fhAsyNoCut  ->Fill(e,asy);
  }
  
  // Asymmetry cut
  if(kDoAsymmetryCut){
   if(asy < fAsymmetryCut){
    if(kNeutralMesonHistos){
     fhInvMassPairAsyCut->Fill(e,invmass);
     fhAnglePairAsyCut  ->Fill(e,angle);
     fhInvMassPairPhi->Fill(e, phi);
     fhInvMassPairEta->Fill(e, eta);
    }
   } else return kFALSE;
  }

  Float_t invmassmaxcut = fInvMassMaxCut;
  Float_t invmassRightBandMinCut = fRightBandMinCut;
  Float_t invmassRightBandMaxCut = fRightBandMaxCut;

  //for EMCAL, pi0s, mass depends strongly with energy for e > 6, loose max cut
  if(e > 10.){
   invmassmaxcut = (fInvMassMaxCutParam[0]+fInvMassMaxCut)+fInvMassMaxCutParam[1]*e+fInvMassMaxCutParam[2]*e*e;
   invmassRightBandMinCut = (fInvMassMaxCutParam[0]+fRightBandMinCut)+fInvMassMaxCutParam[1]*e+fInvMassMaxCutParam[2]*e*e;
   invmassRightBandMaxCut = (fInvMassMaxCutParam[0]+fRightBandMaxCut)+fInvMassMaxCutParam[1]*e+fInvMassMaxCutParam[2]*e*e;
  }

  if(!fAnaMesonType.Contains("SideBand")){
   if(invmass > fInvMassMinCut && invmass < invmassmaxcut){
    if(kNeutralMesonHistos){
     fhInvMassPairAllCut->Fill(e,invmass);
     fhAnglePairAllCut  ->Fill(e,angle);
     fhAsyAllCut  ->Fill(e,asy);
    }
    return kTRUE;
   }
   else{
    return kFALSE;
   }
  }//normal selection
  else if(fAnaMesonType.Contains("SideBand")){
   // select a band around pi0/eta
   if((invmass > fLeftBandMinCut  && invmass < fLeftBandMaxCut)||
     (invmass > invmassRightBandMinCut && invmass < invmassRightBandMaxCut)){
    if(kNeutralMesonHistos){
     fhInvMassPairAllCut->Fill(e,invmass);
     fhAnglePairAllCut  ->Fill(e,angle);
     fhAsyAllCut  ->Fill(e,asy);
    }
    return kTRUE;
   }
   else{
    return kFALSE;
   }
  }
  else return kFALSE;
}

//==============================================
Bool_t AliAnalysisTaskEMCALCaloTrackCorr::IsolatedPhoton(TClonesArray *fPhotonEventIsolated, TClonesArray *fCTSEventIsolated, Int_t fIndexPhotonCan, Double_t ptPhotonCan, Double_t phiPhotonCan, Double_t etaPhotonCan)
{
 
  if(!kDoIsolatedAna) return kFALSE;

  Double_t fCalculatedConeR   = -1.;
  Double_t ptTrack    = -1.;
  Double_t phiTrack   = -999.;
  Double_t etaTrack   = -999.;
  Double_t ptNeutral  = -1.;
  Double_t phiNeutral = -999.;
  Double_t etaNeutral = -999.;

  Float_t  fSumPtInCone = 0.;
  Int_t    fNum         = 0;
  Int_t    fNumFrac     = 0;
  Bool_t   kIsolated    = kFALSE;

   ////skip photon in the edge of r<0.1
  Double_t fPhiCut1 = TMath::Pi()*80./180.-0.1;
  Double_t fPhiCut2 = TMath::Pi()-0.1;
  Double_t fEtaCut1 = -0.6;
  Double_t fEtaCut2 = 0.6;

  if(etaPhotonCan>fEtaCut2 || etaPhotonCan<fEtaCut1 || phiPhotonCan>fPhiCut2 || phiPhotonCan<fPhiCut1) return kFALSE;

  
  if(fDebug>2) printf("candidate photon pt=%f, phi=%f, eta=%f\n",ptPhotonCan,phiPhotonCan,etaPhotonCan);
 
  if(fParticlesInCone==kIsolatedOnlyCharged ||fParticlesInCone==kIsolatedNeutralAndCharged){
   for(Int_t j1=0; j1<fCTSEventIsolated->GetEntries(); j1++){
    AliCaloTrackParticle *traForIso=(AliCaloTrackParticle*)fCTSEventIsolated->At(j1); 
    ptTrack  = traForIso->Pt();
    if(fDebug>2) printf("track pt=%f\n",ptTrack);
    phiTrack = traForIso->Phi();
    if(phiTrack<0.)  phiTrack+=TMath::TwoPi();
    etaTrack = traForIso->Eta();  

    fCalculatedConeR = TMath::Sqrt((etaPhotonCan-etaTrack)*(etaPhotonCan-etaTrack)+ (phiPhotonCan-phiTrack)*(phiPhotonCan-phiTrack));     

    if(fDebug>2) printf("track dR=%f\n",fCalculatedConeR);
 
    if(fCalculatedConeR<fSetConeR){
     if(fDebug>2) printf("track in Cone pt=%f\n",ptTrack);
     //calcuate three cases
     fSumPtInCone+=ptTrack;
     if(ptTrack>fSetPtThreshold)        fNum++;    //larger the threshold pt, fNum++
     if(ptTrack>fSetPtFraction*ptPhotonCan) fNumFrac++;  //larger the fraction of ptPhIso, fNumFrac++

    }////Inside Cone R
   }////loop for Isolated track
  }////Onlycharged || NeutralAndCharged for whether Isolated photon  
   
  if(fDebug>2) printf("track SumPtInCone=%f\n",fSumPtInCone);

  if(fParticlesInCone==kIsolatedOnlyNeutral || fParticlesInCone==kIsolatedNeutralAndCharged){
   for(Int_t i2=0; i2<fPhotonEventIsolated->GetEntries(); i2++){
    AliCaloTrackParticle *ph2=(AliCaloTrackParticle*)fPhotonEventIsolated->At(i2);
    if(!ph2->IsInTrackMatched() || i2==fIndexPhotonCan) continue;
     
    ptNeutral  = ph2->Pt();
    if(fDebug>2) printf("photon or merged pt=%f\n",ptNeutral);
    phiNeutral = ph2->Phi();
    if(phiNeutral<0.)  phiNeutral+=TMath::TwoPi();
    etaNeutral = ph2->Eta(); 
     
    fCalculatedConeR = TMath::Sqrt((etaPhotonCan-etaNeutral)*(etaPhotonCan-etaNeutral)+(phiPhotonCan-phiNeutral)*(phiPhotonCan-phiNeutral));

    if(fDebug>2) printf("photon dR=%f\n",fCalculatedConeR);

    if(fCalculatedConeR<fSetConeR){
     if(fDebug>2)printf("photon in Cone pt=%f\n",ptNeutral);
     //calcuate four cases
     fSumPtInCone+=ptNeutral;
     if(ptNeutral>fSetPtThreshold)        fNum++;      //larger the threshold pt, fNum++
     if(ptNeutral>fSetPtFraction*ptPhotonCan) fNumFrac++;  //larger the fraction of ptPhIso, fNumFrac++
    
     if(fSetPtFraction*ptPhotonCan<fSetPtThreshold) {
      if(ptNeutral>fSetPtThreshold)        fNumFrac++ ;
     }
     else {
      if(ptNeutral>fSetPtFraction*ptPhotonCan) fNumFrac++;
     }

    }////Inside Cone R
   }////loop for neutral
  }////if OnlyNeutral || NeutralAndCharged
  
  if(fDebug>2) printf("track and phton SumPtInCone=%f\n",fSumPtInCone);
  
  //check Isolation, depending on method
  if(fICMethod == kPtThresholdIC){
   if(fNum == -1) kIsolated = kTRUE ;
  }
  else if(fICMethod == kSumPtInConeIC){
   if(fSumPtInCone<fSetSumPtThreshold) kIsolated = kTRUE ;
  }
  else if(fICMethod == kPtFracationIC){
   if(fNumFrac == -1) kIsolated = kTRUE ;
  }
  else if(fICMethod == kSumPtFracationIC){
    //when the fPtFraction*ptC<fSumPtThreshold then consider the later case
   if(fSetPtFraction*ptPhotonCan<fSetSumPtThreshold && fSumPtInCone<fSetSumPtThreshold)     kIsolated = kTRUE ;
   if(fSetPtFraction*ptPhotonCan>fSetSumPtThreshold && fSumPtInCone<fSetPtFraction*ptPhotonCan) kIsolated = kTRUE ;
  }
 
  if(fDebug>2) printf("kIsolated=%d\n",kIsolated); 
  return kIsolated;

}


//=======================================
void AliAnalysisTaskEMCALCaloTrackCorr::SetTriggerBins(Float_t *ptTriggBins)
{
  fTriggPtArray=new Float_t[GetNTriggPtBins()];
  for(Int_t i=0;i<=GetNTriggPtBins(); i++){
   fTriggPtArray[i] = ptTriggBins[i];
  }
}

//========================================_
void AliAnalysisTaskEMCALCaloTrackCorr::SetAssociatedBins(Float_t *ptAssocBins)
{
  fAssocPtArray=new Float_t[GetNAssocPtBins()];
  for(Int_t i=0;i<=GetNAssocPtBins(); i++){
   fAssocPtArray[i] = ptAssocBins[i];
  }
}

//========================================
Bool_t AliAnalysisTaskEMCALCaloTrackCorr::MakeChargedCorrelation(Int_t fTrackIndex, Double_t ptTrigg,
                                                                  Double_t phiTrigg, Double_t etaTrigg)
{  
  if(GetDebug() > 1)printf("AliAnalysisTaskEMCALCaloTrackCorr::MakeChargedCorrelation() - Make trigger particle - charged hadron correlation \n");
   if(fPhotonEvent->GetEntriesFast()<=0) return kFALSE;


  Int_t nTracks = fCTSEvent->GetEntriesFast();
  fhNtracksCorr->Fill(nTracks);

  Double_t ptAssoc  = -999.;
  Double_t phiAssoc = -999. ;
  Double_t etaAssoc = -999.;
  Double_t deltaPhi = -999.;
  Double_t deltaPhiUE = -999.;
  Double_t deltaEta = -999.;
  Double_t pout  = -999.;
  Double_t zT    = -999.; 
  Double_t xE    = -999.; 
  Double_t xELog = -999.; 
  Double_t uedPhi  = -999.;
  Double_t uexE    = -999.;
  Double_t uepout  = -999.;
  Double_t uexELog = -999.; 

  for(Int_t j1 = 0;j1 < nTracks; j1++ ){
   AliCaloTrackParticle *track = (AliCaloTrackParticle *)(fCTSEvent->At(j1)) ;
   if(!track) continue;
   if(kTwoTracksCorr && fTrackIndex ==j1) continue;
   ptAssoc  = track->Pt();
   etaAssoc = track->Eta();
   phiAssoc = track->Phi() ;
   if(phiAssoc < 0) phiAssoc+=TMath::TwoPi();

   fhTrackPtPhi->Fill(ptAssoc, phiAssoc);
   fhTrackPtEta->Fill(ptAssoc, etaAssoc);
   fhTrackPhiEta->Fill(phiAssoc, etaAssoc);

   deltaPhi = phiTrigg-phiAssoc;
   if(deltaPhi < -TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
   if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();

   //jump out this event if near side associated particle pt larger than trigger
   if(kMakeNearSideLeading){
    if(ptAssoc > ptTrigg && (TMath::Abs(deltaPhi) < TMath::PiOver2()))  return kFALSE;
   }
   //jump out this event if there is any other particle with pt larger than trigger
   else if(kMakeAbsoluteLeading){
    if(ptAssoc > ptTrigg)  return kFALSE;
   }
   
   deltaEta = etaTrigg-etaAssoc;

   pout = ptAssoc*TMath::Abs(TMath::Sin(deltaPhi));
   xE   =-ptAssoc/ptTrigg*TMath::Cos(deltaPhi);
   if(xE<0) xE=-xE;
   xELog = TMath::Log(1/xE);
   zT   = ptAssoc/ptTrigg;

   uedPhi = gRandom->Uniform(TMath::Pi()/2,3*TMath::Pi()/2);
   uepout = ptAssoc*TMath::Sin(uedPhi);
   uexE   =(-1*ptAssoc/ptTrigg)*TMath::Cos(uedPhi);
   if(uexE<0) uexE=-uexE;
   uexELog = TMath::Log(1/uexE);

   if(GetDebug()>0)
    printf("AliAnalysisTaskEMCALCaloTrackCorr::MakeChargedCorrelation(): deltaPhi= %f, deltaEta=%f, pout=%f, xE=%f\n",deltaPhi, deltaEta, pout, xE);

   for(Int_t kAssoc=0;kAssoc<GetNAssocPtBins();kAssoc++){
    if(ptAssoc>fAssocPtArray[kAssoc] && ptAssoc<fAssocPtArray[kAssoc+1]){
     if(ptTrigg>fptTriggerBegin && ptTrigg<fptTriggerEnd){
      fhDPhiTriggPtT[kAssoc]->Fill(deltaPhi, ptTrigg);
      fhDEtaTriggPtT[kAssoc]->Fill(deltaEta, ptTrigg);
     }
    }
   }

   if(ptTrigg>fptTriggerBegin && ptTrigg<fptTriggerEnd) {
  
    if(ptAssoc>1. && ptAssoc<5.){
     fhDPhiAssocPt15T->Fill(deltaPhi, ptTrigg);
     fhDEtaAssocPt15T->Fill(deltaEta, ptTrigg);
    }

    if(ptAssoc>fptAssociatedBegin){
    
     fhDPhiTriggPtAssocPt->Fill(deltaPhi, ptAssoc);
     fhDEtaTriggPtAssocPt->Fill(deltaEta, ptAssoc);
    
     if((deltaPhi > fDeltaPhiMinCut) && (deltaPhi < fDeltaPhiMaxCut)){
      fhAssocPtTriggPt->Fill(ptAssoc, ptTrigg);      
      fhxELogTriggPt->Fill(xELog, ptTrigg);
      fhpoutTriggPt->Fill(pout  , ptTrigg);
      fhzTTriggPt->Fill(zT, ptTrigg);
      fhxETriggPt->Fill(xE, ptTrigg);
     }////Only analysis Away Side
    
     if(TMath::Abs(deltaPhi-TMath::Pi())<fDeltaPhiHRSize){
      fhAssocPtTriggPtHR->Fill(ptAssoc, ptTrigg);
      fhxELogTriggPtHR->Fill(xELog, ptTrigg);
      fhpoutTriggPtHR->Fill(pout  , ptTrigg);
      fhzTTriggPtHR->Fill(zT, ptTrigg);
      fhxETriggPtHR->Fill(xE, ptTrigg);
     }////Only analysis head region side in Away Side

     deltaPhiUE = deltaPhi;
     Bool_t kAnaUECorr1 = kFALSE;
     Bool_t kAnaUECorr2 = kFALSE;
     Double_t fUeAwaySide = 0.;
     if(fUeDeltaPhiFix == TMath::Pi()/2.) fUeAwaySide = 3*TMath::Pi()/2.;
     else fUeAwaySide = -1*fUeDeltaPhiFix;
     
     if(kUELeftRight && !kUENearAway) {
      if(deltaPhiUE>(-1*fUeDeltaPhiFix-fUeDeltaPhiSize) && deltaPhiUE < (-1*fUeDeltaPhiFix+fUeDeltaPhiSize))
       kAnaUECorr1 = kTRUE;
      if(deltaPhiUE>(fUeDeltaPhiFix-fUeDeltaPhiSize) && deltaPhiUE < (fUeDeltaPhiFix+fUeDeltaPhiSize))
       kAnaUECorr2 = kTRUE;
     }
     if(!kUELeftRight && kUENearAway) {
      if((deltaPhiUE>-1*fUeDeltaPhiFix && deltaPhiUE<(-1*fUeDeltaPhiFix + fUeDeltaPhiSize)) || (deltaPhiUE>(fUeDeltaPhiFix-fUeDeltaPhiSize) && deltaPhiUE<fUeDeltaPhiFix))
       kAnaUECorr1 = kTRUE;
      if((deltaPhiUE>fUeDeltaPhiFix && deltaPhiUE<(fUeDeltaPhiFix + fUeDeltaPhiSize)) || (deltaPhiUE>(fUeAwaySide - fUeDeltaPhiSize) && deltaPhiUE<fUeAwaySide))
       kAnaUECorr2 = kTRUE;
     }

     ////Ue at Near Side: pt_assoc in the side, but dphi at pi/2-3*pi/2
     if(kAnaUECorr1){
      fhNUeAssocPtTriggPt->Fill(ptAssoc, ptTrigg);
      fhNUexELogTriggPt->Fill(uexELog, ptTrigg);
      fhNUepoutTriggPt->Fill(uepout  , ptTrigg);
      fhNUexETriggPt->Fill(uexE, ptTrigg);
      fhNUezTTriggPt->Fill(zT  , ptTrigg);
      fhNUexETriggPt->Fill(uexE, ptTrigg);
      fhNUeDPhiDEta->Fill(deltaPhi, deltaEta);
     }////Define the UE and analysis some UE physics observables at left

      ////Ue at Away Side
     if(kAnaUECorr2){
      fhAUeAssocPtTriggPt->Fill(ptAssoc, ptTrigg);
      fhAUexELogTriggPt->Fill(uexELog, ptTrigg);
      fhAUepoutTriggPt->Fill(uepout, ptTrigg);
      fhAUezTTriggPt->Fill(zT, ptTrigg);
      fhAUexETriggPt->Fill(uexE, ptTrigg);
      fhAUeDPhiDEta->Fill(deltaPhi, deltaEta);
     }////Define the UE and analysis some UE physics observables at left
    }
   }
    
   for(Int_t iTrigg=0;iTrigg<GetNTriggPtBins();iTrigg++){
    if(ptTrigg>=fTriggPtArray[iTrigg] && ptTrigg<fTriggPtArray[iTrigg+1]){
     if(ptAssoc>fptAssociatedBegin){
      fhDPhiAssocPtA[iTrigg]->Fill(deltaPhi, ptAssoc);
      fhDEtaAssocPtA[iTrigg]->Fill(deltaEta, ptAssoc);
     } 

     for(Int_t jAssoc=0;jAssoc<GetNAssocPtBins();jAssoc++){
      if(ptAssoc>fAssocPtArray[jAssoc] && ptAssoc<fAssocPtArray[jAssoc+1]){
       fhDPhiSumPtBin[iTrigg][jAssoc]->Fill(deltaPhi,ptAssoc+ptTrigg);
       fhDEtaSumPtBin[iTrigg][jAssoc]->Fill(deltaEta,ptAssoc+ptTrigg);
       fhDPhiDEtaBin[iTrigg][jAssoc]->Fill(deltaPhi, deltaEta);
      } 
     }////
    }
   } 
   
   track=0;

  }//end loop associated hadrons 

  fhPtPhiLeading->Fill(ptTrigg,phiTrigg);
  fhPtEtaLeading->Fill(ptTrigg, etaTrigg);
 
  return kTRUE;

}

//=======================================
void AliAnalysisTaskEMCALCaloTrackCorr::MakeChargedMixCorrelation(Double_t ptTriggMix, Double_t phiTriggMix, Double_t etaTriggMix, TList *poolMix)
{ 
  Double_t ptAssoc  = -999.;
  Double_t phiAssoc = -999.;
  Double_t etaAssoc = -999.;
  Double_t deltaPhi = -999.;
  Double_t deltaEta = -999.;

  for(Int_t ev=0;ev<poolMix->GetSize();ev++){
   TClonesArray *bgTracks = static_cast<TClonesArray*>(poolMix->At(ev));
   Int_t nTracks=bgTracks->GetEntriesFast();
   for(Int_t j1 = 0;j1 <nTracks; j1++ ){
    AliCaloTrackParticle *track = (AliCaloTrackParticle*)bgTracks->At(j1) ;
    if(!track) continue;
    ptAssoc  = track->Pt();
    etaAssoc = track->Eta();
    phiAssoc = track->Phi() ;
    if(phiAssoc < 0) phiAssoc+=TMath::TwoPi();

    deltaPhi = phiTriggMix-phiAssoc;
    if(deltaPhi < -TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
    if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();
    deltaEta = etaTriggMix-etaAssoc;

    //jump out this event if near side associated particle pt larger than trigger
    if(kMakeNearSideLeading){
     if(ptAssoc > ptTriggMix && (TMath::Abs(deltaPhi) < TMath::PiOver2()))  return;
    }
    //jump out this event if there is any other particle with pt larger than trigger
    else if(kMakeAbsoluteLeading){
     if(ptAssoc > ptTriggMix)  return;
    }

    if(GetDebug()>0)
     printf("AliAnalysisTaskEMCALCaloTrackCorr::MakeChargedMixCorrelation(): deltaPhi= %f, deltaEta=%f\n",deltaPhi, deltaEta);
 
    for(Int_t iTrigg=0;iTrigg<GetNTriggPtBins();iTrigg++){
     if(ptTriggMix>=fTriggPtArray[iTrigg] && ptTriggMix<fTriggPtArray[iTrigg+1]){
      if(ptAssoc>fptAssociatedBegin){
       fhMixDPhiAssocPtA[iTrigg]->Fill(deltaPhi, ptAssoc);
       fhMixDEtaAssocPtA[iTrigg]->Fill(deltaEta, ptAssoc);
      }
  
      for(Int_t jAssoc=0;jAssoc<GetNAssocPtBins();jAssoc++){
       if(ptAssoc>fAssocPtArray[jAssoc] && ptAssoc<fAssocPtArray[jAssoc+1]){
        fhMixDPhiDEtaBin[iTrigg][jAssoc]->Fill(deltaPhi, deltaEta);
       }//end if Associated pt bin
      }//end loop Associated pt bins
     }//end if trigger pt bin 
    }//end loop trigger bins

    if(ptTriggMix>fptTriggerBegin && ptTriggMix<fptTriggerEnd){
     if(ptAssoc>1. && ptAssoc<5.){
      fhMixDPhiAssocPt15T->Fill(deltaPhi, ptTriggMix);
      fhMixDEtaAssocPt15T->Fill(deltaEta, ptTriggMix);
     }

     for(Int_t kAssoc=0;kAssoc<GetNAssocPtBins();kAssoc++){
      if(ptAssoc>fAssocPtArray[kAssoc] && ptAssoc<fAssocPtArray[kAssoc+1]){
       fhMixDPhiTriggPtT[kAssoc]->Fill(deltaPhi, ptTriggMix);
       fhMixDEtaTriggPtT[kAssoc]->Fill(deltaEta, ptTriggMix);
      }
     }//end loop Associated pt bin
    }//end if Associated pt bins
     
   }//end loop associated hadrons

  }//end loop Mixed event in pool 

  fhMixPtPhiLeading->Fill(ptTriggMix, phiTriggMix);
  fhMixPtEtaLeading->Fill(ptTriggMix, etaTriggMix);
  
}

//========================================


//============================================
void  AliAnalysisTaskEMCALCaloTrackCorr::SetTrackCuts(AliESDtrackCuts * cuts)
{

  if(fESDtrackCuts) delete fESDtrackCuts ;
  fESDtrackCuts = cuts ;

}

//=============================================
void AliAnalysisTaskEMCALCaloTrackCorr::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

  fhNEvents->Fill(1,fnEvents);
  Printf("fnEvents=%d", fnEvents);

  fhNEventsAnalyized->Fill(1,fEventAnalyized);
  Printf("fEventAnalyized=%d", fEventAnalyized);

}
