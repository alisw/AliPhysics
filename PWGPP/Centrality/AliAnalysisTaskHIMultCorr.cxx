//-*- Mode: C++ -*-

#if __GNUS__ >= 3
using namespace std;
#endif

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TMath.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliTracker.h" 
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliESDtrackCuts.h"
#include "AliKineTrackCuts.h"
#include "AliMCParticle.h"
#include "AliESDVZERO.h"
#include "AliAnalysisTaskHIMultCorr.h"
#include "AliCentrality.h"

#include "AliESDVZERO.h"
#include "AliTriggerAnalysis.h"
#include "AliMultiplicity.h"

#include "AliMultiplicityCorrelations.h"

// Task for HI Multiplicity correlation checks
// Author: Jochen Thaeder <jochen@thaeder.de>

ClassImp(AliAnalysisTaskHIMultCorr)

#define USE_STREAMER 0

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
AliAnalysisTaskHIMultCorr::AliAnalysisTaskHIMultCorr(const char *name) :
  AliAnalysisTaskSE(name),
  fpcstream(NULL),
  fIsMC(kFALSE),
  fHStat(NULL),
  fOutList(NULL),
  fESD(NULL), fESDTrackCuts(NULL), fESDTrackCuts2(NULL),
  fUseCentralitySel(0),
  fCentralityBin(-1),
  fCentralitySPDBin(-1),
  fCentralityVZEROBin(-1),
  fCentralitySPD(-1.),
  fCentralityVZERO(-1.),
  fMaxVertexZ(30.),
  fTriggerAnalysis(NULL),
  fCorrObj(NULL), 
  fCorrObjCent0(NULL), fCorrObjCent1(NULL), fCorrObjCent2(NULL) {
  // Constructor   

#if USE_STREAMER
  fpcstream = new TTreeSRedirector("eventInfoCorr.root");
#endif

  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskHIMultCorr::~AliAnalysisTaskHIMultCorr() {
  // Destructor

  if (fTriggerAnalysis)
    delete fTriggerAnalysis;
  fTriggerAnalysis = NULL;

#if USE_STREAMER
  if (fpcstream)
    delete fpcstream;
  fpcstream = NULL;
#endif

}

/*
 * ---------------------------------------------------------------------------------
 *                                    Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisTaskHIMultCorr::UserCreateOutputObjects() {
  // Create histograms
  // Called once

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fOutList = new TList;
  fOutList->SetName(GetName()) ;
  fOutList->SetOwner(kTRUE) ;

  fTriggerAnalysis = new AliTriggerAnalysis;

  // ------------------------------------------------------------------
  // -- Correlations
  // ------------------------------------------------------------------

  fCorrObj = new AliMultiplicityCorrelations;
  fCorrObj->Initialize("All");
  fOutList->Add(fCorrObj->GetHistList());

  fCorrObjCent0 = new AliMultiplicityCorrelations;
  fCorrObjCent0->Initialize("Cent_0-5");
  //  fCorrObjCent0->SetCleanSample(1700.,2900.);      // Cleaning of Sample
  fOutList->Add(fCorrObjCent0->GetHistList());

  fCorrObjCent1 = new AliMultiplicityCorrelations;
  fCorrObjCent1->Initialize("Cent_80-90");
  fOutList->Add(fCorrObjCent1->GetHistList());
  
  fCorrObjCent2 = new AliMultiplicityCorrelations;
  fCorrObjCent2->Initialize("Cent_70-80");
  //  fCorrObjCent2->SetCleanSample(0.,180.);         // Cleaning of Sample
  fOutList->Add(fCorrObjCent2->GetHistList());

  if (fIsMC) {
    fCorrObj->SetIsMC();
    fCorrObjCent0->SetIsMC();
    fCorrObjCent1->SetIsMC();
    fCorrObjCent2->SetIsMC();
  }

  if (fESDTrackCuts) {
    fCorrObj->SetESDTrackCuts(fESDTrackCuts);
    fCorrObjCent0->SetESDTrackCuts(fESDTrackCuts);
    fCorrObjCent1->SetESDTrackCuts(fESDTrackCuts);
    fCorrObjCent2->SetESDTrackCuts(fESDTrackCuts);
  }
  else
    AliError("No ESD trackCuts!!");

  if (fESDTrackCuts2) {
    fCorrObj->SetESDTrackCuts2(fESDTrackCuts2);
    fCorrObjCent0->SetESDTrackCuts2(fESDTrackCuts2);
    fCorrObjCent1->SetESDTrackCuts2(fESDTrackCuts2);
    fCorrObjCent2->SetESDTrackCuts2(fESDTrackCuts2);
  }

  // ------------------------------------------------------------------
  // -- Centrality Distributions
  // ------------------------------------------------------------------
  
  fOutList->Add(new TH1F("fVZEROCentrality",   "VZERO Centrality Distribution;Centrality;N_{Events}",     102,0.,101.));
  fOutList->Add(new TH1F("fVZEROCentralityBin","VZERO Centrality Bin Distribution;Centrality;N_{Events}", 102,0.,101.));

  fOutList->Add(new TH1F("fSPDCentrality",   "SPD Centrality Distribution;Centrality;N_{Events}",         102,0.,101.));
  fOutList->Add(new TH1F("fSPDCentralityBin","SPD Centrality Bin Distribution;Centrality;N_{Events}",     102,0.,101.));

  fOutList->Add(new TH2F("fCorrCentrality",
			 "Centrality - SPD vs VZERO;Centrality_{SPD};Centrality_{VZERO}",      102,0.,101., 102,0.,101.));
  fOutList->Add(new TH2F("fCorrCentralityBin",
			 "Centrality Bin -  SPD vs VZERO;Centrality_{SPD};Centrality_{VZERO}", 102,0.,101., 102,0.,101.));

  // ------------------------------------------------------------------
  // -- Vertex : TPC - Tracks
  // ------------------------------------------------------------------
  
  fOutList->Add(new TH2F("fCorrVertexTracks",
			 "Tracks Vertex - Status vs Contributers;Status;N_{contributers}",     2,-1.,3., 500,0.,3000.));
  fOutList->Add(new TH2F("fCrossCorrVertexTracks",
			 "Tracks Vertex - TPC Status vs Contributers;Status;N_{contributers}", 2,-1.,3., 500,0.,3000.));

  fOutList->Add(new TH2F("fCorrVertexTPC",
			 "TPC Vertex - Status vs Contributers;Status;N_{contributers}",        2,-1.,3., 500,0.,3000.));
  fOutList->Add(new TH2F("fCrossCorrVertexTPC",
			 "TPC Vertex - Tracks Status vs Contributers;Status;N_{contributers}", 2,-1.,3., 500,0.,3000.));

  fOutList->Add(new TH2F("fCorrVertexNCont",
			 "N Contributers - Vertex Tracks vs Vertex TPC;N_{contributers}^{Tracks};N_{contributers}^{TPC}", 
			 500,0.,3000., 500,0.,3000.));

  fOutList->Add(new TH1F("fDeltaVx","Vertex Tracks - Vertex TPC #Deltax;#Delta X [cm]", 201,-10.,10.));
  fOutList->Add(new TH1F("fDeltaVy","Vertex Tracks - Vertex TPC #Deltay;#Delta Y [cm]", 201,-10.,10.));
  fOutList->Add(new TH1F("fDeltaVz","Vertex Tracks - Vertex TPC #Deltaz;#Delta z [cm]", 201,-10.,10.));

  // ------------------------------------------------------------------
  // -- Vertex Z : TPC - Tracks - SPD
  // ------------------------------------------------------------------

  fOutList->Add(new TH2F("fVzTPCvsSPD",   
			 "Z - Vertex TPC vs Vertex SPD;Vz_{TPC} [cm];Vz_{SPD} [cm]",       401,-20.,20., 401,-20.,20.));
  fOutList->Add(new TH2F("fVzTPCvsTracks", 
			 "Z - Vertex TPC vs Vertex Tracks;Vz_{TPC} [cm];Vz_{Tracks} [cm]", 401,-20.,20., 401,-20.,20.));
  fOutList->Add(new TH2F("fVzSPDvsTracks",
			 "Z - Vertex SPD vs Vertex Tracks;Vz_{SPD} [cm];Vz_{Tracks} [cm]", 401,-20.,20., 401,-20.,20.));

  // ------------------------------------------------------------------
  // -- Cuts 
  // ------------------------------------------------------------------

  fOutList->Add(new TH1F("fHStat","Cut statistics", 20,0.,19));
  fHStat = dynamic_cast<TH1F*>(fOutList->Last());

  TH1::AddDirectory(oldStatus);
}

//________________________________________________________________________
void AliAnalysisTaskHIMultCorr::UserExec(Option_t *) {
  // Called for each event
  
#if USE_STREAMER
    Int_t   i1 = -1;
    Int_t   i2 = -1;
    Float_t f1 = -1.;
    Float_t f2 = -1.;
    Float_t f3 = -1.;
#endif

  // -- Setup Event
  // ----------------
  if ( !SetupEvent() ) {
#if USE_STREAMER
    (*fpcstream) <<"event"<<
      "nTracksTPC=" << i1 <<
      "nTracks="    << i2 <<
      "v0A="        << f1 <<
      "v0C="        << f2 <<
      "v0Corr="     << f3 << "\n";
#endif
    return;
  }
    
  // -- Centrality dependence
  // --------------------------
  if (fCentralityBin < 0) { 
#if USE_STREAMER
    (*fpcstream) <<"event"<<
      "nTracksTPC=" << i1 <<
      "nTracks="    << i2 <<
      "v0A="        << f1 <<
      "v0C="        << f2 <<
      "v0Corr="     << f3 << "\n";
#endif
    return;
  }
    
  // -- Process Event
  // ------------------
  Int_t iResult = 0;
  if (fCorrObj) iResult = fCorrObj->ProcessEvent(fESD);
  if (iResult == -3) {
#if USE_STREAMER
    (*fpcstream) <<"event"<<
      "nTracksTPC=" << i1 <<
      "nTracks="    << i2 <<
      "v0A="        << f1 <<
      "v0C="        << f2 <<
      "v0Corr="     << f3 << "\n";
#endif
    return;
  }
  
  // -- After "Jochen's cut"
  fHStat->Fill(7);
  fHStat->Fill(17);

  if      (fCentralityBin == 0  && fCorrObjCent0) fCorrObjCent0->ProcessEvent(fESD);
  else if (fCentralityBin == 80 && fCorrObjCent1) fCorrObjCent1->ProcessEvent(fESD);
  else if (fCentralityBin == 70 && fCorrObjCent2) fCorrObjCent2->ProcessEvent(fESD);

  // -- Marian's Streamer
  // ----------------------
#if USE_STREAMER
  if (fCorrObj) {
    i1 = fCorrObj->GetNTracksTPC();
    i2 = fCorrObj->GetNTracks();
    f1 = fCorrObj->GetVZEROA();
    f2 = fCorrObj->GetVZEROC();
    f3 = fCorrObj->GetVZEROCorr();
  }
  
  (*fpcstream) <<"event"<<
    "nTracksTPC=" << i1 <<
    "nTracks="    << i2 <<
    "v0A="        << f1 <<
    "v0C="        << f2 <<
    "v0Corr="     << f3 << "\n";
#endif

  // -- Post output data
  // ---------------------
  PostData(1,fOutList);

  return;
}      

/*
 * ---------------------------------------------------------------------------------
 *                            Setup Methods - private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
Bool_t AliAnalysisTaskHIMultCorr::SetupEvent() {
  // Setup Reading of event
  
  Bool_t aCuts[] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- ESD Event
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> 
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
    AliError("Could not get ESDInputHandler");
    return kFALSE;
  } 

  fESD = (AliESDEvent*)esdH->GetEvent();
  if (!fESD) {
    AliError("fESD not available");
    return kFALSE;
  }

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- Check object existence
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  const AliESDVertex*    vtxESD    = fESD->GetPrimaryVertexTracks();
  const AliESDVertex*    vtxESDTPC = fESD->GetPrimaryVertexTPC();  
  const AliESDVertex*    vtxESDSPD = fESD->GetPrimaryVertexSPD();  
  const AliMultiplicity* multESD   = fESD->GetMultiplicity();  

  if ( !vtxESD ){
    AliError("No Tracks Vertex");
    return kFALSE;
  }

  if ( !vtxESDTPC ){
    AliError("No TPC Vertex");
    return kFALSE;
  }

  if ( !vtxESDSPD ){
    AliError("No SPD Vertex");
    return kFALSE;
  }

  if ( !multESD ){
    AliError("No Multiplicity");
    return kFALSE;
  }

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- After Physics Selection;
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  aCuts[0] = kFALSE;

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- ZDC cut
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  if ( !fIsMC ) {
    if (!fTriggerAnalysis->ZDCTimeTrigger(fESD)) 
      aCuts[1] = kTRUE;
  }

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- Francesco Prino cut's
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  if (!vtxESDTPC->GetStatus()) 
    aCuts[2] = kTRUE;

  if( vtxESDTPC->GetNContributors() < (-10.+0.25*multESD->GetNumberOfITSClusters(0)) ) 
    aCuts[3] = kTRUE;

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- Vertex Cuts - Tracks 
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  //  if (!vtxESD->GetStatus())
  //    aCuts[4] = kTRUE;
  
  if(vtxESDTPC->GetZ() > fMaxVertexZ || vtxESDTPC->GetZ() < (-1.*fMaxVertexZ)) 
    aCuts[5] = kTRUE;

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- CentralityBin
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  if (GetCentralityBin() < 0)
    aCuts[6] = kTRUE;

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- Marian's debug streamer
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
#if USE_STREAMER
  //  TFile *inFile=(dynamic_cast<TChain*>(GetInputData(0)))->GetCurrentFile();
  //  TString fileName(inFile->GetName()); 
  static Int_t runNo = fESD->GetRunNumber();
  static  UInt_t tt = fESD->GetTimeStamp();
  static Float_t zvTPC = vtxESDTPC->GetZ();
  static Float_t zvSPD = vtxESDSPD->GetZ();
  static Float_t zvTra = vtxESD->GetZ();
  static Int_t ncontTPCV = vtxESDTPC->GetNContributors();
  static Int_t ncontSPDV = vtxESDSPD->GetNContributors();
  static Int_t ncontTraV = vtxESD->GetNContributors();
  static Bool_t cutPrino = aCuts[3];
  static Bool_t cutZDCTi = aCuts[1];
  static Int_t spd1 =  multESD->GetNumberOfITSClusters(0);
  static Int_t spd2 =  multESD->GetNumberOfITSClusters(1);

  (*fpcstream) << "event" <<
    "run="     << runNo                << //runNo
    "time="    << tt                   <<
    //    "fname="   << fileName             <<
    //    "eventNr=" << eventNumber<<
    "zvTPC="   << zvTPC                << //zvertex TPC
    "zvSPD="   << zvSPD                << //zvertex SPD
    "zvTra="   << zvTra                << //zvertex Tracks
    "ncontTPCV"<< ncontTPCV            << // N contributors to TPC vtx
    "ncontSPDV"<< ncontSPDV            << // N contributors to SPD vtx
    "ncontTraV"<< ncontTraV            << // N contributors to Tracks vtx
    "cutPrino="<< cutPrino             << //francecos cut 
    "cutZDCTi="<< cutZDCTi             << //zdctiming cut 
    "spd1="    << spd1                 << // N clus in SPD0
    "spd2="    << spd2                 << // N clus in SPD1
    "centSPD=" << fCentralitySPD       <<
    "centSPDB="<< fCentralitySPDBin    <<
    "centV0="  << fCentralityVZERO     <<
    "centV0B=" << fCentralityVZEROBin;
#endif
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- Fill statistics / reject event
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
 
  Bool_t isRejected = kFALSE;

  for (Int_t idx = 0; idx < 7; ++idx) {
    if (aCuts[idx])
      isRejected = kTRUE;
    else
      fHStat->Fill(idx);
  }
  
  for (Int_t idx = 0; idx < 7; ++idx) {
    if (aCuts[idx])
      break;
    fHStat->Fill(10+idx);
  }

  if (isRejected)
    return kFALSE;
  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- Fill Centrality Histograms
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  (static_cast<TH1F*>(fOutList->FindObject("fVZEROCentrality")))->Fill(fCentralityVZERO);
  (static_cast<TH1F*>(fOutList->FindObject("fVZEROCentralityBin")))->Fill(fCentralityVZEROBin);

  (static_cast<TH1F*>(fOutList->FindObject("fSPDCentrality")))->Fill(fCentralitySPD);
  (static_cast<TH1F*>(fOutList->FindObject("fSPDCentralityBin")))->Fill(fCentralitySPDBin);

  (static_cast<TH2F*>(fOutList->FindObject("fCorrCentrality")))->Fill(fCentralitySPD,fCentralityVZERO);
  (static_cast<TH2F*>(fOutList->FindObject("fCorrCentralityBin")))->Fill(fCentralitySPDBin,fCentralityVZEROBin);

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- Fill Vertex Histograms
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  
  (static_cast<TH2F*>(fOutList->FindObject("fCorrVertexTracks")))->Fill(vtxESD->GetStatus(),vtxESD->GetNContributors());
  (static_cast<TH2F*>(fOutList->FindObject("fCorrVertexTPC")))->Fill(vtxESDTPC->GetStatus(),vtxESDTPC->GetNContributors());
  (static_cast<TH2F*>(fOutList->FindObject("fCrossCorrVertexTracks")))->Fill(vtxESDTPC->GetStatus(),vtxESD->GetNContributors());
  (static_cast<TH2F*>(fOutList->FindObject("fCrossCorrVertexTPC")))->Fill(vtxESD->GetStatus(),vtxESDTPC->GetNContributors());
  (static_cast<TH2F*>(fOutList->FindObject("fCorrVertexNCont")))->Fill(vtxESD->GetNContributors(),vtxESDTPC->GetNContributors());

  (static_cast<TH1F*>(fOutList->FindObject("fDeltaVx")))->Fill(vtxESD->GetX() - vtxESDTPC->GetX());
  (static_cast<TH1F*>(fOutList->FindObject("fDeltaVy")))->Fill(vtxESD->GetY() - vtxESDTPC->GetY());
  (static_cast<TH1F*>(fOutList->FindObject("fDeltaVz")))->Fill(vtxESD->GetZ() - vtxESDTPC->GetZ());

  (static_cast<TH2F*>(fOutList->FindObject("fVzTPCvsSPD")))->Fill(vtxESDTPC->GetZ(),vtxESDSPD->GetZ());
  (static_cast<TH2F*>(fOutList->FindObject("fVzTPCvsTracks")))->Fill(vtxESDTPC->GetZ(),vtxESD->GetZ());
  (static_cast<TH2F*>(fOutList->FindObject("fVzSPDvsTracks")))->Fill(vtxESDSPD->GetZ(),vtxESD->GetZ());

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  return kTRUE;
}

//________________________________________________________________________
Int_t AliAnalysisTaskHIMultCorr::GetCentralityBin(){
  // Get Centrality bin
  
  fCentralityBin = -1;

  if (fUseCentralitySel == 0)
    return fCentralityBin;
  
  AliCentrality *esdCentrality = fESD->GetCentrality();
  fCentralityVZERO  = esdCentrality->GetCentralityPercentile("V0M");  
  fCentralitySPD    = esdCentrality->GetCentralityPercentile("CL1");

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  fCentralityVZEROBin = esdCentrality->GetCentralityClass10("V0M")*10;
  if (fCentralityVZEROBin == 0 && esdCentrality->GetCentralityClass5("V0M") == 1) 
    fCentralityVZEROBin = 5;

  fCentralitySPDBin = esdCentrality->GetCentralityClass10("CL1")*10;
  if (fCentralitySPDBin == 0 && esdCentrality->GetCentralityClass5("CL1") == 1) 
    fCentralitySPDBin = 5;

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  if ( fUseCentralitySel == 1 )
    fCentralityBin = fCentralityVZEROBin;
  else if ( fUseCentralitySel == 2 )
    fCentralityBin = fCentralitySPDBin;

  if ( fCentralityBin == 100 )
    fCentralityBin = -1;

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      
  return fCentralityBin;
}

//________________________________________________________________________
void AliAnalysisTaskHIMultCorr::Terminate(Option_t *){
  // Terminate
#if USE_STREAMER
  if (fpcstream) fpcstream->GetFile()->Write();
#endif

}
