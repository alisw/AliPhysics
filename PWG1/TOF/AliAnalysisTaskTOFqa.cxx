/*  created by fbellini@cern.ch on 14/09/2010 */
/*  last modified by fbellini   on 31/03/2010 */


#ifndef ALIANALYSISTASKTOFQA_CXX
#define ALIANALYSISTASKTOFQA_CXX

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliAnalysisTaskTOFqa.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliTOFRawStream.h"
#include "AliTOFGeometry.h"

ClassImp(AliAnalysisTaskTOFqa)

//________________________________________________________________________
AliAnalysisTaskTOFqa::AliAnalysisTaskTOFqa() :
  fRunNumber(0), 
  fESD(0x0), 
  fTrackFilter(0x0), 
  fVertex(0x0),
  fESDpid(new AliESDpid()),
  fNTOFtracks(0), 
//  fNPrimaryTracks(0), 
  fHlist(0x0),
  fHlistTimeZero(0x0),
  fHlistPID(0x0)
 {
  // Default constructor
   
   for (Int_t j=0;j<3;j++ ) {
     if (j<3) fT0[j]=0.0;
     fSigmaSpecie[j]=0.0;
     fTrkExpTimes[j]=0.0;
     fThExpTimes[j]=0.0;
   }
 }
//________________________________________________________________________
AliAnalysisTaskTOFqa::AliAnalysisTaskTOFqa(const char *name) : 
  AliAnalysisTaskSE(name), 
  fRunNumber(0), 
  fESD(0x0), 
  fTrackFilter(0x0),
  fVertex(0x0),
fESDpid(new AliESDpid()),
  fNTOFtracks(0), 
  // fNPrimaryTracks(0),  
  fHlist(0x0),
  fHlistTimeZero(0),
  fHlistPID(0x0)
 {
  // Constructor
  // Define input and output slots here
   Info("AliAnalysisTaskTOFqa","Calling Constructor");
   
   for (Int_t j=0;j<5;j++ ) {
     if (j<3) fT0[j]=0.0;
     fSigmaSpecie[j]=0.0;
     fTrkExpTimes[j]=0.0;
     fThExpTimes[j]=0.0;
   }
   // Input slot #0 works with a TChain
   DefineInput(0, TChain::Class());
   
   // Output slot #0 writes into a TH1 container
   // Output slot #1 writes into a user defined  container
   DefineOutput(1, TList::Class());
   DefineOutput(2, TList::Class());
   DefineOutput(3, TList::Class());
 }

//________________________________________________________________________
AliAnalysisTaskTOFqa::AliAnalysisTaskTOFqa(const AliAnalysisTaskTOFqa& copy) 
: AliAnalysisTaskSE(), 
  fRunNumber(copy.fRunNumber), 
  fESD(copy.fESD), 
  fTrackFilter(copy.fTrackFilter), 
  fVertex(copy.fVertex),
  fESDpid(copy.fESDpid),
  fNTOFtracks(copy.fNTOFtracks), 
  //fNPrimaryTracks(copy.fNPrimaryTracks), 
  fHlist(copy.fHlist),
  fHlistTimeZero(copy.fHlistTimeZero),
  fHlistPID(copy.fHlistPID)
{
  // Copy constructor
   for (Int_t j=0;j<5;j++ ) {
     if (j<3) fT0[j]=copy.fT0[j];
     fSigmaSpecie[j]=copy.fSigmaSpecie[j];
     fTrkExpTimes[j]=copy.fTrkExpTimes[j];
     fThExpTimes[j]=copy.fThExpTimes[j];
   }
  

}

//___________________________________________________________________________
AliAnalysisTaskTOFqa& AliAnalysisTaskTOFqa::operator=(const AliAnalysisTaskTOFqa& copy) 
{
  //
  // Assignment operator
  //
  if (this!=&copy) {
    AliAnalysisTaskSE::operator=(copy) ;
    fRunNumber=copy.fRunNumber; 
    fESD=copy.fESD;
    fTrackFilter=copy.fTrackFilter;
    fVertex=copy.fVertex;
    fESDpid=copy.fESDpid;
    fNTOFtracks=copy.fNTOFtracks; 
    //fNPrimaryTracks=copy.fNPrimaryTracks; 
    for (Int_t j=0;j<5;j++ ) {
      if (j<3) fT0[j]=copy.fT0[j];
      fSigmaSpecie[j]=copy.fSigmaSpecie[j];
      fTrkExpTimes[j]=copy.fTrkExpTimes[j];
      fThExpTimes[j]=copy.fThExpTimes[j];
    }
    fHlist=copy.fHlist;
    fHlistTimeZero=copy.fHlistTimeZero;
    fHlistPID=copy.fHlistPID;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskTOFqa::~AliAnalysisTaskTOFqa() {
  //
  //destructor
  //

  Info("~AliAnalysisTaskTOFqa","Calling Destructor");
  if (fESDpid) delete fESDpid;
  if (fVertex) delete fVertex;
  if (fTrackFilter) delete fTrackFilter;

  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if (fHlist) {
    delete fHlist;
    fHlist = 0;
  }
  if (fHlistTimeZero) {
    delete fHlistTimeZero;
    fHlistTimeZero = 0;
  }
  if (fHlistPID){
    delete fHlistPID;
    fHlistPID = 0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskTOFqa::UserCreateOutputObjects()
{
  //Defines output objects and histograms
  Info("CreateOutputObjects","CreateOutputObjects (TList) of task %s", GetName());
  OpenFile(1);
  if (!fHlist) fHlist = new TList();	
  fHlist->SetOwner(kTRUE);
  if (!fHlistTimeZero) fHlistTimeZero = new TList();	
  fHlistTimeZero->SetOwner(kTRUE);
  if (!fHlistPID) fHlistPID = new TList();	
  fHlistPID->SetOwner(kTRUE);

//0
  TH1I* hTOFmatchedESDperEvt = new TH1I("hTOFmatchedPerEvt", "Matched TOF tracks per event (|#eta| #leq 0.9 and pT #geq 0.3 GeV/c);TOF-matched ESD tracks;Events", 100, 0, 100) ;  
  hTOFmatchedESDperEvt->Sumw2() ;
  hTOFmatchedESDperEvt->SetLineWidth(2);
  hTOFmatchedESDperEvt->SetLineColor(kBlue);
  hTOFmatchedESDperEvt->SetMarkerStyle(20);
  hTOFmatchedESDperEvt->SetMarkerSize(0.8);
  hTOFmatchedESDperEvt->SetMarkerColor(kBlue);
  fHlist->AddLast(hTOFmatchedESDperEvt) ;
  //1
  TH1F* hTOFmatchedESDtime = new TH1F("hTOFmatchedESDtime", "Matched  ESDs tracks: TOF Time spectrum; t [ns];Counts", 250, 0., 610. ) ; 
  hTOFmatchedESDtime->Sumw2() ;
  hTOFmatchedESDtime->SetLineWidth(2);
  hTOFmatchedESDtime->SetLineColor(kBlue);
  hTOFmatchedESDtime->SetFillColor(kBlue);
  hTOFmatchedESDtime->SetDrawOption("BAR");
  fHlist->AddLast(hTOFmatchedESDtime) ;
  //2
  TH1F* hTOFmatchedESDrawTime = new TH1F("hTOFmatchedESDrawTime", "Matched ESDs tracks: TOF raw Time spectrum;t_{raw} [ns];Counts", 250, 0., 610.) ; 
  hTOFmatchedESDrawTime->Sumw2() ;
  hTOFmatchedESDrawTime->SetLineWidth(2);
  hTOFmatchedESDrawTime->SetLineColor(kAzure+2);
  hTOFmatchedESDrawTime->SetFillColor(kAzure+2);
  hTOFmatchedESDrawTime->SetDrawOption("BAR");
  fHlist->AddLast(hTOFmatchedESDrawTime) ;
  //3
  TH1F* hTOFmatchedESDToT = new TH1F("hTOFmatchedESDToT", "Matched ESDs tracks: TOF ToT spectrum; ToT [ns];Counts",100, 0., 48.8) ; 
  hTOFmatchedESDToT->Sumw2() ;
  hTOFmatchedESDToT->SetLineColor(kOrange+1);
  hTOFmatchedESDToT->SetMarkerColor(kOrange+1);
  hTOFmatchedESDToT->SetFillColor(kOrange+1);
  hTOFmatchedESDToT->SetDrawOption("BAR");
  fHlist->AddLast(hTOFmatchedESDToT) ;
  //4
  TH1F* hTOFmatchedESDtrkLength  = new TH1F("hTOFmatchedESDtrkLength", "Matched ESDs tracks length; Track length [cm];Counts", 1200, -400., 800) ; 
  hTOFmatchedESDtrkLength->Sumw2();
  hTOFmatchedESDtrkLength->SetLineColor(kViolet-3);
  hTOFmatchedESDtrkLength->SetMarkerColor(kViolet-3);
  hTOFmatchedESDtrkLength->SetFillColor(kViolet-3);
  hTOFmatchedESDtrkLength->SetDrawOption("BAR"); 
  fHlist->AddLast(hTOFmatchedESDtrkLength);
  //5
  TH1F* hTOFmatchedESDP  = new TH1F("hTOFmatchedESDP", "TPC-TOF matched tracks momentum distribution (GeV/c); p (GeV/c);tracks", 500,0.,5.) ;  
  hTOFmatchedESDP->Sumw2() ;
  hTOFmatchedESDP->SetLineColor(kBlue);
  hTOFmatchedESDP->SetMarkerStyle(20);
  hTOFmatchedESDP->SetMarkerSize(0.7);
  hTOFmatchedESDP->SetMarkerColor(kBlue);
  fHlist->AddLast(hTOFmatchedESDP) ; 
  //6
  TH1F* hTOFmatchedESDPt  = new TH1F("hTOFmatchedESDPt", "TPC-TOF matched tracks p_{T} distribution (GeV/c); p_{T}(GeV/c);tracks", 500,0.,5.) ;  
  hTOFmatchedESDPt->Sumw2() ;
  hTOFmatchedESDPt->SetLineColor(kBlue);
  hTOFmatchedESDPt->SetMarkerStyle(21);
  hTOFmatchedESDPt->SetMarkerSize(0.7);
  hTOFmatchedESDPt->SetMarkerColor(kBlue);
  fHlist->AddLast(hTOFmatchedESDPt) ; 

  //7
  TH1F* hTOFmatchedESDeta = new TH1F("hTOFmatchedESDeta", "Matched ESDtracks #eta (p_{T} #geq 0.5 GeV/c); #eta;Counts", 200, -1., 1.) ; 
  hTOFmatchedESDeta->Sumw2();
  hTOFmatchedESDeta->SetLineColor(kBlue);
  fHlist->AddLast(hTOFmatchedESDeta) ; 
  //8
   TH1F* hTOFmatchedESDphi = new TH1F("hTOFmatchedESDphi", "Matched ESDtracks #phi; #phi (deg);Counts", 72, 0., 360.) ; 
  hTOFmatchedESDphi->Sumw2();
  hTOFmatchedESDphi->SetLineColor(kBlue);
  fHlist->AddLast(hTOFmatchedESDphi) ; 

  //9
  TH1F* hESDprimaryTrackP = new TH1F("hESDprimaryTrackP", "All ESDs tracks p_{T} distribution (GeV/c); p_{T}(GeV/c);tracks", 500, 0., 5.0) ;  
  hESDprimaryTrackP->Sumw2();
  hESDprimaryTrackP->SetLineWidth(1);
  hESDprimaryTrackP->SetMarkerStyle(24);
  hESDprimaryTrackP->SetMarkerSize(0.7);
  hESDprimaryTrackP->SetMarkerColor(kRed);
  hESDprimaryTrackP->SetLineColor(kRed);
  fHlist->AddLast(hESDprimaryTrackP);
  //10
  TH1F* hESDprimaryTrackPt = new TH1F("hESDprimaryTrackPt", "ESDs primary tracks p_{T} distribution (GeV/c); p_{T}(GeV/c);tracks", 500, 0.0, 5.0) ;  
  hESDprimaryTrackPt->Sumw2();
  hESDprimaryTrackPt->SetLineWidth(1);
  hESDprimaryTrackPt->SetMarkerStyle(25);
  hESDprimaryTrackPt->SetMarkerSize(0.7);
  hESDprimaryTrackPt->SetLineColor(kRed);
  hESDprimaryTrackPt->SetMarkerColor(kRed);
  fHlist->AddLast(hESDprimaryTrackPt);
  //11
  TH1F* hTOFprimaryESDeta = new TH1F("hTOFprimaryESDeta", "Primary ESDtracks #eta (p_{T} #geq 0.5 GeV/c); #eta;Counts",200, -1., 1.) ; 
  hTOFprimaryESDeta->Sumw2();
  hTOFprimaryESDeta->SetLineColor(kRed);
  fHlist->AddLast(hTOFprimaryESDeta) ; 
  //12
  TH1F* hTOFprimaryESDphi = new TH1F("hTOFprimaryESDphi", "Primary ESDtracks #phi;#phi (deg);Counts", 72, 0., 360.) ; 
  hTOFprimaryESDphi->Sumw2();
  hTOFprimaryESDphi->SetLineColor(kRed);
  fHlist->AddLast(hTOFprimaryESDphi) ; 
  //13
  TH2F* hTOFmatchedDxVsPtPos = new TH2F("hTOFmatchedDxVsPtPos", "Dx vs p_{T} for positive tracks;p_{T} (GeV/c); Dx [cm]; hits", 500,0.,5.,200, -10., 10.) ; 
  hTOFmatchedDxVsPtPos->Sumw2();
  fHlist->AddLast(hTOFmatchedDxVsPtPos) ; 
 //14
  TH2F* hTOFmatchedDxVsPtNeg = new TH2F("hTOFmatchedDxVsPtNeg", "Dx vs p_{T} for negative tracks;p_{T} (GeV/c); Dx [cm]; hits", 500,0.,5.,200, -10., 10.) ; 
  hTOFmatchedDxVsPtNeg->Sumw2();
  fHlist->AddLast(hTOFmatchedDxVsPtNeg) ; 

  //15
  TH2F* hTOFmatchedDzVsStrip = new TH2F("hTOFmatchedDzVsStrip", "Dz vs strip; strip (#eta); Dz [cm]; hits", 92,0.,92.,200, -10., 10.) ; 
  hTOFmatchedDzVsStrip->Sumw2();
  fHlist->AddLast(hTOFmatchedDzVsStrip) ; 
 
 //----------------------------------------------timeZero QA plots
  //TimeZero 0
  TH1D* hEventT0DetAND = new TH1D("hEventT0DetAND", "Event timeZero from T0AC detector ; t0 [ps]; events", 1000, -25000., 25000. ) ; 
  hEventT0DetAND->Sumw2() ;
  hEventT0DetAND->SetLineWidth(2);
  hEventT0DetAND->SetLineColor(kRed);
  hEventT0DetAND->SetFillColor(kRed);
  fHlistTimeZero->AddLast(hEventT0DetAND) ;

  //TImeZero 1
  TH1D* hEventT0DetA = new TH1D("hEventT0DetA", "Event timeZero from T0A detector; t0 [ps]; events", 1000, -25000., 25000. ) ; 
  hEventT0DetA->Sumw2() ;
  hEventT0DetA->SetLineWidth(2);
  hEventT0DetA->SetLineColor(kBlue);
  hEventT0DetA->SetFillColor(kBlue);
  fHlistTimeZero->AddLast(hEventT0DetA) ;

   //TImeZero 2
  TH1D* hEventT0DetC = new TH1D("hEventT0DetC", "Event timeZero from T0C detector; t0 [ps]; events", 1000, -25000., 25000.) ; 
  hEventT0DetC->Sumw2() ;
  hEventT0DetC->SetLineWidth(2);
  hEventT0DetC->SetLineColor(kGreen);
  hEventT0DetC->SetFillColor(kGreen);
  fHlistTimeZero->AddLast(hEventT0DetC);

   //TimeZero 3
  TH1F* hT0DetRes = new TH1F("hT0DetRes", "T0 detector (T0A-T0C)/2; (T0A-T0C)/2 [ps]; events", 200, -500.,500. ) ; 
  hT0DetRes->Sumw2() ;
  hT0DetRes->SetMarkerStyle(24);
  hT0DetRes->SetMarkerSize(0.7);
  hT0DetRes->SetMarkerColor(kMagenta+2);
  hT0DetRes->SetLineColor(kMagenta+2);
  hT0DetRes->SetFillColor(kMagenta+2);  
  fHlistTimeZero->AddLast(hT0DetRes) ; 

     //timeZero 4
  TH1F* hT0fill = new TH1F("hT0fill", "Event timeZero of fill; t0 [ps]; events", 1000, -25000., 25000. ) ; 
  hT0fill->Sumw2() ;
  hT0fill->SetMarkerStyle(20);
  hT0fill->SetMarkerColor(kBlack);
  hT0fill->SetLineColor(kBlack);
  fHlistTimeZero->AddLast(hT0fill) ; 

  //TimeZero 5
  TH1F* hT0TOF = new TH1F("hT0TOF", "Event timeZero estimated by TOF; t0 [ps]; events", 1000, -25000., 25000. ) ; 
  hT0TOF->Sumw2() ;
  hT0TOF->SetMarkerStyle(20);
  hT0TOF->SetMarkerColor(kBlue);
  hT0TOF->SetLineColor(kBlue);
  hT0TOF->SetFillColor(kBlue);
  fHlistTimeZero->AddLast(hT0TOF) ;


   //timeZero 6
  TH1F* hT0T0 = new TH1F("hT0T0", "Event timeZero measured by T0 detector (best between AC, A, C); t0 [ps]; events", 1000, -25000.,25000. ) ; 
  hT0T0->Sumw2() ;
  hT0T0->SetMarkerStyle(20);
  hT0T0->SetMarkerColor(kGreen+1);
  hT0T0->SetLineColor(kGreen+1);
  hT0T0->SetFillColor(kGreen+1);
  fHlistTimeZero->AddLast(hT0T0) ; 

   //timeZero 7
  TH1F* hT0best = new TH1F("hT0best", "Event timeZero estimated as T0best; t0 [ps]; events", 1000, -25000.,25000. ) ; 
  hT0best->Sumw2() ;
  hT0best->SetMarkerStyle(20);
  hT0best->SetMarkerColor(kRed);
  hT0best->SetLineColor(kRed);
  hT0best->SetFillColor(kRed); 
  fHlistTimeZero->AddLast(hT0best) ; 

   //TimeZero 8
  TH1F* hT0fillRes = new TH1F("hT0fillRes", "Resolution of fillT0; #sigma_{bestT0} [ps];events", 250, 0.,250. ) ; 
  hT0fillRes->Sumw2() ;
  hT0fillRes->SetMarkerStyle(21);
  hT0fillRes->SetMarkerColor(kBlack);
  hT0fillRes->SetLineColor(kBlack);
  hT0fillRes->SetFillColor(kBlack); 
  fHlistTimeZero->AddLast(hT0fillRes) ; 
 
  //TimeZero 9
  TH1F* hT0TOFRes = new TH1F("hT0TOFRes", "Resolution of timeZero from TOF; #sigma_{TOFT0} [ps];events", 250, 0.,250. ) ; 
  hT0TOFRes->Sumw2() ;
  hT0TOFRes->SetLineWidth(1);
  hT0TOFRes->SetMarkerStyle(21);
  hT0TOFRes->SetMarkerColor(kBlue);
  hT0TOFRes->SetLineColor(kBlue);
  hT0TOFRes->SetFillColor(kBlue); 
  fHlistTimeZero->AddLast(hT0TOFRes) ; 

   //TimeZero 10
  TH1F* hT0T0res = new TH1F("hT0T0res", "Resolution of timeZero from T0;#sigma_{T0T0}  [ps];events", 250, -0., 250. ) ; 
  hT0T0res->Sumw2() ;
  hT0T0res->SetMarkerStyle(21);
  hT0T0res->SetMarkerColor(kGreen+1);
  hT0T0res->SetLineColor(kGreen+1);
  hT0T0res->SetFillColor(kGreen+1); 
  fHlistTimeZero->AddLast(hT0T0res) ; 

   //TimeZero 11
  TH1F* hT0bestRes = new TH1F("hT0bestRes", "Resolution of bestT0; #sigma_{bestT0} [ps];events", 250, 0.,250. ) ; 
  hT0bestRes->Sumw2() ;
  hT0fillRes->SetMarkerStyle(21);
  hT0fillRes->SetMarkerColor(kRed);
  hT0fillRes->SetLineColor(kRed);
  hT0fillRes->SetFillColor(kRed); 
  fHlistTimeZero->AddLast(hT0bestRes) ; 

  //timeZero 12
  TH2F* hT0TOFvsNtrk = new TH2F("hT0TOFvsNtrk", "Event timeZero estimated by TOF vs. number of tracks in event;TOF-matching tracks; t0 [ps]", 100, 0., 100.,1000,-25000.,25000. ) ; 
  hT0TOFvsNtrk->Sumw2() ;
  fHlistTimeZero->AddLast(hT0TOFvsNtrk) ;

//--------------------------------------------- TOF PID QA plots
  //PID 0
  TH2F* hTOFmatchedESDpVsBeta  = new TH2F("hTOFmatchedESDpVsBeta", "Matched ESDs tracks beta vs. p; p(GeV/c); beta", 500, 0.0, 5.0, 150,0., 1.5) ; 
  fHlistPID->AddLast(hTOFmatchedESDpVsBeta);
  
  //PID 1 
  TH1F* hTOFmatchedMass= new TH1F("hTOFmatchedMass","Matched ESD tracks mass distribution - (L>0); M (GeV/c^{2}); entries", 500,0., 5. );
  hTOFmatchedMass->Sumw2();
  hTOFmatchedMass->SetLineWidth(2);
  hTOFmatchedMass->SetLineColor(kBlue);
  hTOFmatchedMass->SetLineColor(kBlue);
  fHlistPID->AddLast(hTOFmatchedMass);
  
  //PID 2
  TH2F* hTOFmatchedExpTimePiVsEta = new TH2F("hTOFmatchedExpTimePiVsEta", "ESDs t_{TOF}-t_{#pi,exp} (from tracking); strip (#eta); t_{TOF}-t_{#pi,exp} [ps]",92,0,92,2000, -5000., 5000. ) ; 
  hTOFmatchedExpTimePiVsEta->Sumw2() ;
  fHlistPID->AddLast(hTOFmatchedExpTimePiVsEta) ;
  
  //PID 3
  TH1F* hTOFmatchedExpTimePi = new TH1F("hTOFmatchedExpTimePi", "ESDs t_{TOF}-t_{#pi,exp} (from tracking); t_{TOF}-t_{#pi,exp} [ps];Counts",5000, -25000., 25000. ) ; 
  hTOFmatchedExpTimePi->Sumw2() ;
  hTOFmatchedExpTimePi->SetLineWidth(1);
  hTOFmatchedExpTimePi->SetLineColor(kRed);
  hTOFmatchedExpTimePi->SetMarkerStyle(20);
  hTOFmatchedExpTimePi->SetMarkerSize(0.8); 
  hTOFmatchedExpTimePi->SetMarkerColor(kRed);
  fHlistPID->AddLast(hTOFmatchedExpTimePi) ;
  
  //PID 4
  TH2F* hTOFmatchedExpTimePiVsP = new TH2F("hTOFmatchedExpTimePiVsP", "ESDs t_{TOF}-t_{#pi,exp} (from tracking) Vs P ; p (GeV/c);t_{TOF}-t_{#pi,exp} [ps];Counts",500, 0.,5.,1000, -25000., 25000. ) ; 
  hTOFmatchedExpTimePiVsP->Sumw2() ;
  fHlistPID->AddLast(hTOFmatchedExpTimePiVsP) ;

  //PID 5
  TH1F* hTOFtheoreticalExpTimePi = new TH1F("hTOFtheoreticalExpTimePi", "ESDs t_{TOF}-t_{#pi,exp} (theoretical); t_{TOF}-t_{#pi,exp} [ps];Counts", 5000, -25000., 25000. ) ; 
  hTOFtheoreticalExpTimePi->Sumw2() ;
  hTOFtheoreticalExpTimePi->SetLineWidth(1);
  hTOFtheoreticalExpTimePi->SetLineColor(kRed);
  hTOFtheoreticalExpTimePi->SetMarkerStyle(24);
  hTOFtheoreticalExpTimePi->SetMarkerSize(0.8); 
  hTOFtheoreticalExpTimePi->SetMarkerColor(kRed);
  fHlistPID->AddLast(hTOFtheoreticalExpTimePi) ;

  //PID 6
  TH2F* hTOFtheoreticalExpTimePiVsP = new TH2F("hTOFtheoreticalExpTimePiVsP", "ESDs t_{TOF}-t_{#pi,exp} (theoretical) Vs P ; p (GeV/c);t_{TOF}-t_{#pi,exp} [ps];Counts",500, 0.,5.,1000, -25000., 25000. ) ; 
  hTOFtheoreticalExpTimePiVsP->Sumw2() ;
  fHlistPID->AddLast(hTOFtheoreticalExpTimePiVsP) ;

  //PID 7
  TH2F* hTOFExpSigmaPi = new TH2F("hTOFExpSigmaPi", "ESDs TOF n#sigma_{PID,#pi} vs p_{T}; p_{T} (GeV/c); n#sigma_{PID,#pi};Tracks", 500,0.,5.,200, -10., 10. ) ; 
  hTOFExpSigmaPi->Sumw2() ;
  fHlistPID->AddLast(hTOFExpSigmaPi) ;

  //PID 8
  TH1F* hTOFmatchedExpTimeKa = new TH1F("hTOFmatchedExpTimeKa", "ESDs t_{TOF}-t_{K,exp} (from tracking); t_{TOF}-t_{K,exp} [ps];Counts", 500, -5000., 5000. ) ; 
  hTOFmatchedExpTimeKa->Sumw2() ;
  hTOFmatchedExpTimeKa->SetLineWidth(1);
  hTOFmatchedExpTimeKa->SetLineColor(kBlue);
  hTOFmatchedExpTimeKa->SetMarkerStyle(21);
  hTOFmatchedExpTimeKa->SetMarkerSize(0.8); 
  hTOFmatchedExpTimeKa->SetMarkerColor(kBlue);
  fHlistPID->AddLast(hTOFmatchedExpTimeKa);

  //PID 9
  TH2F* hTOFmatchedExpTimeKaVsP = new TH2F("hTOFmatchedExpTimeKaVsP", "ESDs t_{TOF}-t_{K,exp} (from tracking) Vs P ; p (GeV/c);t_{TOF}-t_{K,exp} [ps];Counts",500, 0.,5.,1000, -25000., 25000. ) ; 
  hTOFmatchedExpTimeKaVsP->Sumw2() ;
  fHlistPID->AddLast(hTOFmatchedExpTimeKaVsP) ; 
  
  //PID 10
  TH1F* hTOFtheoreticalExpTimeKa = new TH1F("hTOFtheoreticalExpTimeKa", "ESDs t_{TOF}-t_{K,exp} (theoretical); t_{TOF}-t_{K,exp} [ps];Counts", 5000, -25000., 25000. ) ; 
  hTOFtheoreticalExpTimeKa->Sumw2() ;
  hTOFtheoreticalExpTimeKa->SetLineWidth(1);
  hTOFtheoreticalExpTimeKa->SetLineColor(kBlue);
  hTOFtheoreticalExpTimeKa->SetMarkerStyle(24);
  hTOFtheoreticalExpTimeKa->SetMarkerSize(0.8); 
  hTOFtheoreticalExpTimeKa->SetMarkerColor(kBlue);
  fHlistPID->AddLast(hTOFtheoreticalExpTimeKa) ;  
  
  //PID 11
  TH2F* hTOFtheoreticalExpTimeKaVsP = new TH2F("hTOFtheoreticalExpTimeKaVsP", "ESDs t_{TOF}-t_{K,exp} (theoretical) Vs P ; p (GeV/c);t_{TOF}-t_{K,exp} [ps];Counts",500, 0.,5.,1000, -25000., 25000. ) ; 
  hTOFtheoreticalExpTimeKaVsP->Sumw2() ;
  fHlistPID->AddLast(hTOFtheoreticalExpTimeKaVsP) ; 
  
  //PID 12
  TH2F* hTOFExpSigmaKa = new TH2F("hTOFExpSigmaKa", "ESDs TOF n#sigma_{PID,K} vs p_{T}; p_{T} (GeV/c);n#sigma_{PID,K};Tracks", 500, 0.,5.,200, -10., 10. ) ; 
  hTOFExpSigmaKa->Sumw2() ;
  fHlistPID->AddLast(hTOFExpSigmaKa) ;
  
  //PID 13
  TH1F* hTOFmatchedExpTimePro = new TH1F("hTOFmatchedExpTimePro", "ESDs t_{TOF}-t_{p,exp} (from tracking); t_{TOF}-t_{p,exp} [ps];Counts", 500, -5000., 5000. ) ; 
  hTOFmatchedExpTimePro->Sumw2() ;
  hTOFmatchedExpTimePro->SetLineWidth(1);
  hTOFmatchedExpTimePro->SetLineColor(kGreen+1);
  hTOFmatchedExpTimePro->SetMarkerStyle(22);
  hTOFmatchedExpTimePro->SetMarkerSize(0.8); 
  hTOFmatchedExpTimePro->SetMarkerColor(kGreen+1);
  fHlistPID->AddLast(hTOFmatchedExpTimePro) ;

   //PID 14
  TH2F* hTOFmatchedExpTimeProVsP = new TH2F("hTOFmatchedExpTimeProVsP", "ESDs t_{TOF}-t_{p,exp} (from tracking) Vs P ; p (GeV/c);t_{TOF}-t_{p,exp} [ps];Counts",500, 0.,5.,1000, -25000., 25000. ) ; 
  hTOFmatchedExpTimeProVsP->Sumw2() ;
  fHlistPID->AddLast(hTOFmatchedExpTimeProVsP) ;
  
  //PID 15
  TH1F* hTOFtheoreticalExpTimePro = new TH1F("hTOFtheoreticalExpTimePro", "ESDs t_{TOF}-t_{p,exp} (theoretical); t_{TOF}-t_{p,exp} [ps];Counts", 500, -5000., 5000. ) ; 
  hTOFtheoreticalExpTimePro->Sumw2() ;
  hTOFtheoreticalExpTimePro->SetLineWidth(1);
  hTOFtheoreticalExpTimePro->SetLineColor(kGreen+1);
  hTOFtheoreticalExpTimePro->SetMarkerStyle(26);
  hTOFtheoreticalExpTimePro->SetMarkerSize(0.8); 
  hTOFtheoreticalExpTimePro->SetMarkerColor(kGreen+1);
  fHlistPID->AddLast(hTOFtheoreticalExpTimePro) ;

  //PID 16
  TH2F* hTOFtheoreticalExpTimeProVsP = new TH2F("hTOFtheoreticalExpTimeProVsP", "ESDs t_{TOF}-t_{p,exp} (theoretical) Vs P ; p (GeV/c);t_{TOF}-t_{p,exp} [ps];Counts",500, 0.,5.,1000, -25000., 25000. ) ; 
  hTOFtheoreticalExpTimeProVsP->Sumw2() ;
  fHlistPID->AddLast(hTOFtheoreticalExpTimeProVsP) ;

  //PID 17
  TH2F* hTOFExpSigmaPro = new TH2F("hTOFExpSigmaPro", "ESDs TOF n#sigma_{PID,p} vs. p_{T}; p_{T} (GeV/c); n#sigma_{PID,p};Tracks", 500, 0.,5.,200, -10., 10. ) ; 
  hTOFExpSigmaPro->Sumw2() ;
  fHlistPID->AddLast(hTOFExpSigmaPro) ;

  PostData(1, fHlist);
  PostData(2, fHlistTimeZero);
  PostData(3, fHlistPID);

}
//________________________________________________________________________
void AliAnalysisTaskTOFqa::UserExec(Option_t *) 
{ 
  /* Main - executed for each event.
    It extracts event information and track information after selecting 
    primary tracks via standard cuts. */
  
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
    Printf("ERROR: Could not get ESDInputHandler");
    return;
  } else {
    fESD = (AliESDEvent*) esdH->GetEvent();
  } 
  
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }
  
  /* info from T0 detector QA */
  for (Int_t j=0;j<3;j++){
    fT0[j]= (Float_t) fESD->GetT0TOF(j);//ps
    if (fT0[j]>90000.) fT0[j]=99999.;//fix old default values to the new one
  }
  Float_t t0cut = 90000.; 
  //Float_t t0cut =3 * t0spread; //use this cut to check t0 used in tof response
  // if(t0cut < 500) t0cut = 500;
  
  if(TMath::Abs(fT0[1]) < t0cut && TMath::Abs(fT0[2]) < t0cut ) {
    //&& TMath::Abs(fT0[2]-fT0[1]) < 500)  //add this condition to check t0 used in tof response
    ((TH1F*)fHlistTimeZero->At(3))->Fill((fT0[2]-fT0[1])*0.5);
    ((TH1F*)fHlistTimeZero->At(0))->Fill(fT0[0]);  
  } 
  if(TMath::Abs(fT0[1]) < t0cut){
    ((TH1F*)fHlistTimeZero->At(1))->Fill(fT0[1]);   
  }
  if(TMath::Abs(fT0[2]) < t0cut){
    ((TH1F*)fHlistTimeZero->At(2))->Fill(fT0[2]);
  }
  
  /*  event timeZero QA via AliESDpid::SetTOFResponse() */
  Double_t timeZero[4]={99999.,99999.,99999.,99999.};
  Double_t timeZeroRes[4]={99999.,99999.,99999.,99999.}; 
  
  for (Int_t j=0;j<4;j++){
    fESDpid->SetTOFResponse(fESD, (AliESDpid::EStartTimeType_t) j);//(fill_t0, tof_t0, t0_t0, best_t0)
    timeZero[j]=fESDpid->GetTOFResponse().GetStartTime(10.); //timeZero for bin pT>10GeV/c
    timeZeroRes[j]=fESDpid->GetTOFResponse().GetStartTimeRes(10.); //timeZero for bin pT>10GeV/c
    ((TH1D*)fHlistTimeZero->At(4+j))->Fill(timeZero[j]);
    ((TH1D*)fHlistTimeZero->At(8+j))->Fill(timeZeroRes[j]);
  }


  /* loop over ESD tracks */
  fNTOFtracks=0;
  // fNPrimaryTracks=0;

  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }

    //primary tracks selection: kTPCrefit and std cuts
    if(!fTrackFilter->IsSelected(track)) continue;
    Double_t eta=track->Eta();
    if (TMath::Abs(eta)>0.9) continue; //cut for acceptance

    Double_t mom=track->P();
    Double_t mom2 = mom*mom;
    Double_t length=track->GetIntegratedLength();
    Double_t pT = track->Pt();
    Double_t phi=track->Phi()*TMath::RadToDeg();
    track->GetIntegratedTimes(fTrkExpTimes);
      
    ((TH1F*)fHlist->At(9))->Fill(mom); 
    ((TH1F*)fHlist->At(10))->Fill(pT); 
    if (pT>=0.5)
      ((TH1F*)fHlist->At(11))->Fill(eta);
    ((TH1F*)fHlist->At(12))->Fill(phi);
    
    //matched tracks selection: kTOFout and kTIME
    if ((track->GetStatus() & AliESDtrack::kTOFout) &&
	(track->GetStatus() & AliESDtrack::kTIME)) {      
      
      Double_t tofTime=track->GetTOFsignal();//in ps
      Double_t tofTimeRaw=track->GetTOFsignalRaw();//in ps
      Double_t tofToT=track->GetTOFsignalToT(); //in ps
      Int_t channel=track->GetTOFCalChannel(); 
      Int_t volId[5]; //(sector, plate,strip,padZ,padX)
      AliTOFGeometry::GetVolumeIndices(channel,volId);
      
      if (pT>=0.3) fNTOFtracks++; //matched counter
      Double_t tof= tofTime*1E-3; // ns, average T0 fill subtracted, no info from T0detector 	 
      ((TH1F*)fHlist->At(1))->Fill(tof); //ns
      ((TH1F*)fHlist->At(2))->Fill(tofTimeRaw*1E-3); //ns
      ((TH1F*)fHlist->At(3))->Fill(tofToT);
      ((TH1F*)fHlist->At(4))->Fill(length);  
      ((TH1F*)fHlist->At(5))->Fill(mom);
      ((TH1F*)fHlist->At(6))->Fill(pT);
      if (pT>=0.5)
	((TH1F*)fHlist->At(7))->Fill(eta);
      ((TH1F*)fHlist->At(8))->Fill(phi);
      if (track->GetSign()>0)
	((TH2F*)fHlist->At(13))->Fill(pT,track->GetTOFsignalDx());
      else ((TH2F*)fHlist->At(14))->Fill(pT,track->GetTOFsignalDx());
      ((TH2F*)fHlist->At(15))->Fill((Int_t)GetStripIndex(volId),track->GetTOFsignalDz());

      //basic PID performance check
      if (tof<=0) {
	printf("WARNING: track with negative TOF time found! Skipping this track for PID checks\n");
	continue;
      }
      if (mom2==0) {
	printf("WARNING: track with negative square momentum found! Skipping this track for PID checks\n");
	continue;
      }
      if (length<=0){
	printf("WARNING: track with negative length found!Skipping this track for PID checks\n");
	continue;
      }
      Double_t c=TMath::C()*1.E-9;// m/ns
      Double_t mass=0.; //GeV
      length =length*0.01; // in meters
      tof=tof*c;
      Double_t beta=length/tof;
      Double_t fact= (tof/length)*(tof/length) -1.;
      if(fact<=0) {
	mass = -mom*TMath::Sqrt(-fact);
      }else{ 
	mass = mom*TMath::Sqrt(fact); 
      }
      ((TH2F*)fHlistPID->At(0))->Fill(mom,beta);
      ((TH1F*) fHlistPID->At(1))->Fill(mass);
      
      //PID sigmas
      Bool_t isValidBeta[AliPID::kSPECIES]={0,0,0,0,0};
      for (Int_t specie = 0; specie < AliPID::kSPECIES; specie++){
	fSigmaSpecie[specie] = fESDpid->GetTOFResponse().GetExpectedSigma(mom, fTrkExpTimes[specie], AliPID::ParticleMass(specie));
	beta=1/TMath::Sqrt(1+AliPID::ParticleMass(specie)*AliPID::ParticleMass(specie)/(mom2));
	if (beta>0) {
	  fThExpTimes[specie]=length*1.E3/(beta*c);//ps
	  isValidBeta[specie]=kTRUE;
	} else {
	  fThExpTimes[specie]=1E-10;
	  isValidBeta[specie]=kFALSE;
	}
      }
      
      if (isValidBeta[AliPID::kPion]){
	((TH2F*)fHlistPID->At(2))->Fill((Int_t)GetStripIndex(volId),tofTime-fTrkExpTimes[AliPID::kPion]);//ps
	((TH1F*)fHlistPID->At(3))->Fill(tofTime-fTrkExpTimes[AliPID::kPion]);//ps
	((TH2F*)fHlistPID->At(4))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kPion]));
	((TH1F*)fHlistPID->At(5))->Fill(tofTime-fThExpTimes[AliPID::kPion]);//ps
	((TH2F*)fHlistPID->At(6))->Fill(mom,(tofTime-fThExpTimes[AliPID::kPion]));	
	((TH2F*)fHlistPID->At(7))->Fill(pT,(tofTime-fTrkExpTimes[AliPID::kPion])/fSigmaSpecie[AliPID::kPion]);
      }
      
      if (isValidBeta[AliPID::kKaon]){
	((TH1F*)fHlistPID->At(8))->Fill(tofTime-fTrkExpTimes[AliPID::kKaon]);//ps
	((TH2F*)fHlistPID->At(9))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kKaon]));
	((TH1F*)fHlistPID->At(10))->Fill(tofTime-fThExpTimes[AliPID::kKaon]);//ps
	((TH2F*)fHlistPID->At(11))->Fill(mom,(tofTime-fThExpTimes[AliPID::kKaon]));
	((TH2F*)fHlistPID->At(12))->Fill(pT,(tofTime-fTrkExpTimes[AliPID::kKaon])/fSigmaSpecie[AliPID::kKaon]);
      }
      if (isValidBeta[AliPID::kProton]){
	((TH1F*)fHlistPID->At(13))->Fill(tofTime-fTrkExpTimes[AliPID::kProton]);//ps
	((TH2F*)fHlistPID->At(14))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kProton]));
	((TH1F*)fHlistPID->At(15))->Fill(tofTime-fThExpTimes[AliPID::kProton]);//ps
	((TH2F*)fHlistPID->At(16))->Fill(mom,(tofTime-fThExpTimes[AliPID::kProton]));
	((TH2F*)fHlistPID->At(17))->Fill(pT,(tofTime-fTrkExpTimes[AliPID::kProton])/fSigmaSpecie[AliPID::kProton]);
      }
    }//matched
  }//end loop on tracks
  
  ((TH1F*)fHlist->At(0))->Fill(fNTOFtracks) ;
  ((TH2F*)fHlistTimeZero->At(12))->Fill(fNTOFtracks,timeZero[AliESDpid::kTOF_T0]);
  
  PostData(1, fHlist);
  PostData(2, fHlistTimeZero);
  PostData(3, fHlistPID);
  
  
}      

//________________________________________________________________________
void AliAnalysisTaskTOFqa::Terminate(Option_t *) 
{
  //check on output validity
  fHlist = dynamic_cast<TList*> (GetOutputData(1));
  if (!fHlist || !fHlistTimeZero) {
    Printf("ERROR: lists not available");
    return;   
  }   
 
}

//---------------------------------------------------------------
Int_t AliAnalysisTaskTOFqa::GetStripIndex(const Int_t * const in)
{
  /* return tof strip index between 0 and 91 */
  
  Int_t nStripA = AliTOFGeometry::NStripA();
  Int_t nStripB = AliTOFGeometry::NStripB();
  Int_t nStripC = AliTOFGeometry::NStripC();

  Int_t iplate = in[1];
  Int_t istrip = in[2];
  
  Int_t stripOffset = 0;
  switch (iplate) {
  case 0:
    stripOffset = 0;
      break;
  case 1:
    stripOffset = nStripC;
    break;
  case 2:
    stripOffset = nStripC+nStripB;
    break;
  case 3:
    stripOffset = nStripC+nStripB+nStripA;
    break;
  case 4:
    stripOffset = nStripC+nStripB+nStripA+nStripB;
    break;
  default:
    stripOffset=-1;
    break;
  };
  
  if (stripOffset<0 || stripOffset>92) return -1;
  else 
    return (stripOffset+istrip);
}

#endif
