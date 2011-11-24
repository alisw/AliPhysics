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
  fEnableAdvancedCheck(kFALSE),
  fExpTimeBinWidth(24.4),
  fExpTimeRangeMin(-25010.),
  fExpTimeRangeMax(25010.),
  fExpTimeSmallRangeMin(-5002.),
  fExpTimeSmallRangeMax(5002.),
  fHlist(0x0),
  fHlistTimeZero(0x0),
  fHlistPID(0x0),
  fHpos(0x0),
  fHneg(0x0)
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
  fEnableAdvancedCheck(kFALSE),
  fExpTimeBinWidth(24.4),
  fExpTimeRangeMin(-25010.),
  fExpTimeRangeMax(25010.),
  fExpTimeSmallRangeMin(-5002.),
  fExpTimeSmallRangeMax(5002.),
  fHlist(0x0),
  fHlistTimeZero(0),
  fHlistPID(0x0),
  fHpos(0x0),
  fHneg(0x0)
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
   DefineOutput(4, TList::Class());
   DefineOutput(5, TList::Class());
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
  fEnableAdvancedCheck(copy.fEnableAdvancedCheck),
  fExpTimeBinWidth(copy.fExpTimeBinWidth),
  fExpTimeRangeMin(copy.fExpTimeRangeMin),
  fExpTimeRangeMax(copy.fExpTimeRangeMax),
  fExpTimeSmallRangeMin(copy.fExpTimeSmallRangeMin),
  fExpTimeSmallRangeMax(copy.fExpTimeSmallRangeMax),
  fHlist(copy.fHlist),
  fHlistTimeZero(copy.fHlistTimeZero),
  fHlistPID(copy.fHlistPID),
  fHpos(copy.fHpos),
  fHneg(copy.fHneg)
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
    fEnableAdvancedCheck=copy.fEnableAdvancedCheck;
    fEnableAdvancedCheck=copy.fEnableAdvancedCheck;
    fExpTimeBinWidth=copy.fExpTimeBinWidth;
    fExpTimeRangeMin=copy.fExpTimeRangeMin;
    fExpTimeRangeMax=copy.fExpTimeRangeMax;
    fExpTimeSmallRangeMin=copy.fExpTimeSmallRangeMin;
    fExpTimeSmallRangeMax=copy.fExpTimeSmallRangeMax;
    
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
    fHpos=copy.fHpos;
    fHneg=copy.fHneg;
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
  if (fHpos){
    delete fHpos;
    fHpos = 0;
  }
  if (fHneg){
    delete fHneg;
    fHneg = 0;
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
  if (!fHpos) fHpos = new TList();	
  fHpos->SetOwner(kTRUE);
  if (!fHneg) fHneg = new TList();	
  fHneg->SetOwner(kTRUE);

  Int_t nExpTimeBins=1;
  Int_t nExpTimeSmallBins=1;
  
  if (fExpTimeRangeMax<fExpTimeRangeMin) {
    SetExpTimeHistoRange(-25010.,25010.);
  }
  nExpTimeBins = TMath::Nint((fExpTimeRangeMax - fExpTimeRangeMin)/fExpTimeBinWidth);//ps
  fExpTimeRangeMax=fExpTimeRangeMin+nExpTimeBins*fExpTimeBinWidth;//ps
  
  if (fExpTimeSmallRangeMax<fExpTimeSmallRangeMin) {
    SetExpTimeHistoSmallRange(-5002.,5002.);
  }
  nExpTimeSmallBins = TMath::Nint((fExpTimeSmallRangeMax - fExpTimeSmallRangeMin)/fExpTimeBinWidth);//ps
  fExpTimeSmallRangeMax=fExpTimeSmallRangeMin+nExpTimeSmallBins*fExpTimeBinWidth;//ps
  
  //  nExpTimeSmallBins, fExpTimeSmallRangeMin, fExpTimeSmallRangeMax
  //  nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax

//0
  TH1I* hTOFmatchedESDperEvt = new TH1I("hTOFmatchedPerEvt", "Matched TOF tracks per event (|#eta| #leq 0.8 and pT #geq 0.3 GeV/c);TOF-matched ESD tracks;Events", 100, 0, 100) ;  
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
  TH1F* hESDprimaryTrackP = new TH1F("hESDprimaryTrackP", "All ESDs tracks p distribution (GeV/c); p(GeV/c);tracks", 500, 0., 5.0) ;  
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
  TH1F* hT0fillRes = new TH1F("hT0fillRes", "Resolution of fillT0; #sigma_{fillT0} [ps];events", 250, 0.,250. ) ; 
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
  TH1F* hT0T0Res = new TH1F("hT0T0Res", "Resolution of timeZero from T0;#sigma_{T0T0}  [ps];events", 250, -0., 250. ) ; 
  hT0T0Res->Sumw2() ;
  hT0T0Res->SetMarkerStyle(21);
  hT0T0Res->SetMarkerColor(kGreen+1);
  hT0T0Res->SetLineColor(kGreen+1);
  hT0T0Res->SetFillColor(kGreen+1); 
  fHlistTimeZero->AddLast(hT0T0Res) ; 

   //TimeZero 11
  TH1F* hT0bestRes = new TH1F("hT0bestRes", "Resolution of bestT0; #sigma_{bestT0} [ps];events", 250, 0.,250. ) ; 
  hT0bestRes->Sumw2() ;
  hT0bestRes->SetMarkerStyle(21);
  hT0bestRes->SetMarkerColor(kRed);
  hT0bestRes->SetLineColor(kRed);
  hT0bestRes->SetFillColor(kRed); 
  fHlistTimeZero->AddLast(hT0bestRes) ; 

  //timeZero 12
  TH2F* hT0TOFvsNtrk = new TH2F("hT0TOFvsNtrk", "Event timeZero estimated by TOF vs. number of tracks in event;TOF-matching tracks; t0 [ps]", 100, 0., 100.,1000,-25000.,25000. ) ; 
  hT0TOFvsNtrk->Sumw2() ;
  fHlistTimeZero->AddLast(hT0TOFvsNtrk) ;

 //TimeZero 13
  TH2F* hEventT0MeanVsVtx = new TH2F("hEventT0MeanVsVtx", "T0 detector: mean vs vertex ; (t0_{A}-t0_{C})/2 [ns]; (t0_{A}+t0_{C})/2 [ns]; events", 500, -25., 25., 500, -25., 25. ) ; 
  hEventT0MeanVsVtx->Sumw2() ;
  fHlistTimeZero->AddLast(hEventT0MeanVsVtx) ;

 //TimeZero 14
  TH2F* hEventV0MeanVsVtx = new TH2F("hEventV0MeanVsVtx", "V0 detector: mean vs vertex ; (V0_{A}-V0_{C})/2 [ns]; (V0_{A}+V0_{C})/2 [ns]; events", 500, -50., 50., 500, -50., 50. ) ; 
  hEventV0MeanVsVtx->Sumw2() ;
  fHlistTimeZero->AddLast(hEventV0MeanVsVtx) ;
//--------------------------------------------- TOF PID QA plots
  //PID 0
  TH2F* hTOFmatchedESDpVsBeta  = new TH2F("hTOFmatchedESDpVsBeta", "Matched ESDs tracks beta vs. p; p(GeV/c); beta", 500, 0.0, 5.0, 150, 0., 1.5) ; 
  fHlistPID->AddLast(hTOFmatchedESDpVsBeta);
  
  //PID 1 
  TH1F* hTOFmatchedMass= new TH1F("hTOFmatchedMass","Matched ESD tracks mass distribution - (L>0); M (GeV/c^{2}); entries", 500, 0., 5. );
  hTOFmatchedMass->Sumw2();
  hTOFmatchedMass->SetLineWidth(2);
  hTOFmatchedMass->SetLineColor(kBlue);
  hTOFmatchedMass->SetLineColor(kBlue);
  fHlistPID->AddLast(hTOFmatchedMass);
  
  //PID 2
  TH2F* hTOFmatchedExpTimePiVsEta = new TH2F("hTOFmatchedExpTimePiVsEta", "ESDs t_{TOF}-t_{#pi,exp} (from tracking); strip (#eta); t_{TOF}-t_{#pi,exp} [ps]", 92, 0, 92,  nExpTimeSmallBins, fExpTimeSmallRangeMin, fExpTimeSmallRangeMax) ; 
  hTOFmatchedExpTimePiVsEta->Sumw2() ;
  fHlistPID->AddLast(hTOFmatchedExpTimePiVsEta) ;
  
  //PID 3
  TH1F* hTOFmatchedExpTimePi = new TH1F("hTOFmatchedExpTimePi", "ESDs t_{TOF}-t_{#pi,exp} (from tracking); t_{TOF}-t_{#pi,exp} [ps];Counts", nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFmatchedExpTimePi->Sumw2() ;
  hTOFmatchedExpTimePi->SetLineWidth(1);
  hTOFmatchedExpTimePi->SetLineColor(kRed);
  hTOFmatchedExpTimePi->SetMarkerStyle(20);
  hTOFmatchedExpTimePi->SetMarkerSize(0.8); 
  hTOFmatchedExpTimePi->SetMarkerColor(kRed);
  fHlistPID->AddLast(hTOFmatchedExpTimePi) ;
  
  //PID 4
  TH2F* hTOFmatchedExpTimePiVsP = new TH2F("hTOFmatchedExpTimePiVsP", "ESDs t_{TOF}-t_{#pi,exp} (from tracking) Vs P ; p (GeV/c);t_{TOF}-t_{#pi,exp} [ps];Counts",500, 0.,5., nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFmatchedExpTimePiVsP->Sumw2() ;
  fHlistPID->AddLast(hTOFmatchedExpTimePiVsP) ;

  //PID 5
  TH1F* hTOFtheoreticalExpTimePi = new TH1F("hTOFtheoreticalExpTimePi", "ESDs t_{TOF}-t_{#pi,exp} (theoretical); t_{TOF}-t_{#pi,exp} [ps];Counts", nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFtheoreticalExpTimePi->Sumw2() ;
  hTOFtheoreticalExpTimePi->SetLineWidth(1);
  hTOFtheoreticalExpTimePi->SetLineColor(kRed);
  hTOFtheoreticalExpTimePi->SetMarkerStyle(24);
  hTOFtheoreticalExpTimePi->SetMarkerSize(0.8); 
  hTOFtheoreticalExpTimePi->SetMarkerColor(kRed);
  fHlistPID->AddLast(hTOFtheoreticalExpTimePi) ;

  //PID 6
  TH2F* hTOFtheoreticalExpTimePiVsP = new TH2F("hTOFtheoreticalExpTimePiVsP", "ESDs t_{TOF}-t_{#pi,exp} (theoretical) Vs P ; p (GeV/c);t_{TOF}-t_{#pi,exp} [ps];Counts",500, 0.,5., nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFtheoreticalExpTimePiVsP->Sumw2() ;
  fHlistPID->AddLast(hTOFtheoreticalExpTimePiVsP) ;

  //PID 7
  TH2F* hTOFExpSigmaPi = new TH2F("hTOFExpSigmaPi", "ESDs TOF n#sigma_{PID,#pi} vs p_{T}; p_{T} (GeV/c); n#sigma_{PID,#pi};Tracks", 500,0.,5.,200, -10., 10. ) ; 
  hTOFExpSigmaPi->Sumw2() ;
  fHlistPID->AddLast(hTOFExpSigmaPi) ;

  //PID 8
  TH1F* hTOFmatchedExpTimeKa = new TH1F("hTOFmatchedExpTimeKa", "ESDs t_{TOF}-t_{K,exp} (from tracking); t_{TOF}-t_{K,exp} [ps];Counts", nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFmatchedExpTimeKa->Sumw2() ;
  hTOFmatchedExpTimeKa->SetLineWidth(1);
  hTOFmatchedExpTimeKa->SetLineColor(kBlue);
  hTOFmatchedExpTimeKa->SetMarkerStyle(21);
  hTOFmatchedExpTimeKa->SetMarkerSize(0.8); 
  hTOFmatchedExpTimeKa->SetMarkerColor(kBlue);
  fHlistPID->AddLast(hTOFmatchedExpTimeKa);

  //PID 9
  TH2F* hTOFmatchedExpTimeKaVsP = new TH2F("hTOFmatchedExpTimeKaVsP", "ESDs t_{TOF}-t_{K,exp} (from tracking) Vs P ; p (GeV/c);t_{TOF}-t_{K,exp} [ps];Counts",500, 0.,5.,nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFmatchedExpTimeKaVsP->Sumw2() ;
  fHlistPID->AddLast(hTOFmatchedExpTimeKaVsP) ; 
  
  //PID 10
  TH1F* hTOFtheoreticalExpTimeKa = new TH1F("hTOFtheoreticalExpTimeKa", "ESDs t_{TOF}-t_{K,exp} (theoretical); t_{TOF}-t_{K,exp} [ps];Counts", nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFtheoreticalExpTimeKa->Sumw2() ;
  hTOFtheoreticalExpTimeKa->SetLineWidth(1);
  hTOFtheoreticalExpTimeKa->SetLineColor(kBlue);
  hTOFtheoreticalExpTimeKa->SetMarkerStyle(24);
  hTOFtheoreticalExpTimeKa->SetMarkerSize(0.8); 
  hTOFtheoreticalExpTimeKa->SetMarkerColor(kBlue);
  fHlistPID->AddLast(hTOFtheoreticalExpTimeKa) ;  
  
  //PID 11
  TH2F* hTOFtheoreticalExpTimeKaVsP = new TH2F("hTOFtheoreticalExpTimeKaVsP", "ESDs t_{TOF}-t_{K,exp} (theoretical) Vs P ; p (GeV/c);t_{TOF}-t_{K,exp} [ps];Counts",500, 0.,5., nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFtheoreticalExpTimeKaVsP->Sumw2() ;
  fHlistPID->AddLast(hTOFtheoreticalExpTimeKaVsP) ; 
  
  //PID 12
  TH2F* hTOFExpSigmaKa = new TH2F("hTOFExpSigmaKa", "ESDs TOF n#sigma_{PID,K} vs p_{T}; p_{T} (GeV/c);n#sigma_{PID,K};Tracks", 500, 0.,5.,200, -10., 10. ) ; 
  hTOFExpSigmaKa->Sumw2() ;
  fHlistPID->AddLast(hTOFExpSigmaKa) ;
  
  //PID 13
  TH1F* hTOFmatchedExpTimePro = new TH1F("hTOFmatchedExpTimePro", "ESDs t_{TOF}-t_{p,exp} (from tracking); t_{TOF}-t_{p,exp} [ps];Counts", nExpTimeBins, fExpTimeRangeMin,fExpTimeRangeMax) ; 
  hTOFmatchedExpTimePro->Sumw2() ;
  hTOFmatchedExpTimePro->SetLineWidth(1);
  hTOFmatchedExpTimePro->SetLineColor(kGreen+1);
  hTOFmatchedExpTimePro->SetMarkerStyle(22);
  hTOFmatchedExpTimePro->SetMarkerSize(0.8); 
  hTOFmatchedExpTimePro->SetMarkerColor(kGreen+1);
  fHlistPID->AddLast(hTOFmatchedExpTimePro) ;

   //PID 14
  TH2F* hTOFmatchedExpTimeProVsP = new TH2F("hTOFmatchedExpTimeProVsP", "ESDs t_{TOF}-t_{p,exp} (from tracking) Vs P ; p (GeV/c);t_{TOF}-t_{p,exp} [ps];Counts",500, 0.,5., nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFmatchedExpTimeProVsP->Sumw2() ;
  fHlistPID->AddLast(hTOFmatchedExpTimeProVsP) ;
  
  //PID 15
  TH1F* hTOFtheoreticalExpTimePro = new TH1F("hTOFtheoreticalExpTimePro", "ESDs t_{TOF}-t_{p,exp} (theoretical); t_{TOF}-t_{p,exp} [ps];Counts", nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFtheoreticalExpTimePro->Sumw2() ;
  hTOFtheoreticalExpTimePro->SetLineWidth(1);
  hTOFtheoreticalExpTimePro->SetLineColor(kGreen+1);
  hTOFtheoreticalExpTimePro->SetMarkerStyle(26);
  hTOFtheoreticalExpTimePro->SetMarkerSize(0.8); 
  hTOFtheoreticalExpTimePro->SetMarkerColor(kGreen+1);
  fHlistPID->AddLast(hTOFtheoreticalExpTimePro) ;

  //PID 16
  TH2F* hTOFtheoreticalExpTimeProVsP = new TH2F("hTOFtheoreticalExpTimeProVsP", "ESDs t_{TOF}-t_{p,exp} (theoretical) Vs P ; p (GeV/c);t_{TOF}-t_{p,exp} [ps];Counts",500, 0.,5., nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFtheoreticalExpTimeProVsP->Sumw2() ;
  fHlistPID->AddLast(hTOFtheoreticalExpTimeProVsP) ;

  //PID 17
  TH2F* hTOFExpSigmaPro = new TH2F("hTOFExpSigmaPro", "ESDs TOF n#sigma_{PID,p} vs. p_{T}; p_{T} (GeV/c); n#sigma_{PID,p};Tracks", 500, 0.,5.,200, -10., 10. ) ; 
  hTOFExpSigmaPro->Sumw2() ;
  fHlistPID->AddLast(hTOFExpSigmaPro) ;

   //PID 18
  TH2F* hTOFmatchedExpTimePiVsPTRDPos = new TH2F("hTOFmatchedExpTimePiVsPTRDPos", "ESDs t_{TOF}-t_{p,exp} (from tracking) Vs P (#phi_{outerTPC}#leq 30 || 150 #leq#phi_{outerTPC}#leq 230 || #phi_{outerTPC}#geq 310 ) ; p (GeV/c);t_{TOF}-t_{p,exp} [ps];Counts",500, 0.,5., nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFmatchedExpTimePiVsPTRDPos->Sumw2() ;
  if (fEnableAdvancedCheck)
    fHlistPID->AddLast(hTOFmatchedExpTimePiVsPTRDPos) ;

   //PID 19
  TH2F* hTOFmatchedExpTimePiVsPNoTRDPos = new TH2F("hTOFmatchedExpTimePiVsPNoTRDPos", "ESDs t_{TOF}-t_{p,exp} (from tracking) Vs P (50 #leq#phi_{outerTPC}#leq 130 || 250 #leq#phi_{outerTPC}#leq 290); p (GeV/c);t_{TOF}-t_{p,exp} [ps];Counts",500, 0.,5., nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFmatchedExpTimePiVsPNoTRDPos->Sumw2() ;
  if (fEnableAdvancedCheck)
    fHlistPID->AddLast(hTOFmatchedExpTimePiVsPNoTRDPos) ;

   //PID 20
  TH2F* hTOFmatchedExpTimePiVsPTRDNeg = new TH2F("hTOFmatchedExpTimePiVsPTRDNeg", "ESDs t_{TOF}-t_{p,exp} (from tracking) Vs P (#phi_{outerTPC}#leq 30 || 150 #leq#phi_{outerTPC}#leq 230 || #phi_{outerTPC}#geq 310 ) ; p (GeV/c);t_{TOF}-t_{p,exp} [ps];Counts",500, 0.,5., nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFmatchedExpTimePiVsPTRDNeg->Sumw2() ;
  if (fEnableAdvancedCheck)
    fHlistPID->AddLast(hTOFmatchedExpTimePiVsPTRDNeg) ;

   //PID 21
  TH2F* hTOFmatchedExpTimePiVsPNoTRDNeg = new TH2F("hTOFmatchedExpTimePiVsPNoTRDNeg", "ESDs t_{TOF}-t_{p,exp} (from tracking) Vs P (50 #leq#phi_{outerTPC}#leq 130 || 250 #leq#phi_{outerTPC}#leq 290); p (GeV/c);t_{TOF}-t_{p,exp} [ps];Counts",500, 0.,5.,nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFmatchedExpTimePiVsPNoTRDNeg->Sumw2() ;
  if (fEnableAdvancedCheck)
    fHlistPID->AddLast(hTOFmatchedExpTimePiVsPNoTRDNeg) ;
  
 //PID 22
  TH2F* hTOFmatchedExpTimeKaVsPTRDPos = new TH2F("hTOFmatchedExpTimeKaVsPTRDPos", "ESDs t_{TOF}-t_{p,exp} (from tracking) Vs P (#phi_{outerTPC}#leq 30 || 150 #leq#phi_{outerTPC}#leq 230 || #phi_{outerTPC}#geq 310 ) ; p (GeV/c);t_{TOF}-t_{p,exp} [ps];Counts",500, 0.,5.,nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFmatchedExpTimeKaVsPTRDPos->Sumw2() ;
  if (fEnableAdvancedCheck)
    fHlistPID->AddLast(hTOFmatchedExpTimeKaVsPTRDPos) ;
  
   //PID 23
  TH2F* hTOFmatchedExpTimeKaVsPNoTRDPos = new TH2F("hTOFmatchedExpTimeKaVsPNoTRDPos", "ESDs t_{TOF}-t_{p,exp} (from tracking) Vs P (50 #leq#phi_{outerTPC}#leq 130 || 250 #leq#phi_{outerTPC}#leq 290); p (GeV/c);t_{TOF}-t_{p,exp} [ps];Counts",500, 0.,5.,nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFmatchedExpTimeKaVsPNoTRDPos->Sumw2() ;
    if (fEnableAdvancedCheck)
      fHlistPID->AddLast(hTOFmatchedExpTimeKaVsPNoTRDPos) ;

   //PID 24
  TH2F* hTOFmatchedExpTimeKaVsPTRDNeg = new TH2F("hTOFmatchedExpTimeKaVsPTRDNeg", "ESDs t_{TOF}-t_{p,exp} (from tracking) Vs P (#phi_{outerTPC}#leq 30 || 150 #leq#phi_{outerTPC}#leq 230 || #phi_{outerTPC}#geq 310 ) ; p (GeV/c);t_{TOF}-t_{p,exp} [ps];Counts",500, 0.,5.,nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFmatchedExpTimeKaVsPTRDNeg->Sumw2() ;
  if (fEnableAdvancedCheck)
    fHlistPID->AddLast(hTOFmatchedExpTimeKaVsPTRDNeg) ;

   //PID 25
  TH2F* hTOFmatchedExpTimeKaVsPNoTRDNeg = new TH2F("hTOFmatchedExpTimeKaVsPNoTRDNeg", "ESDs t_{TOF}-t_{p,exp} (from tracking) Vs P (50 #leq#phi_{outerTPC}#leq 130 || 250 #leq#phi_{outerTPC}#leq 290); p (GeV/c);t_{TOF}-t_{p,exp} [ps];Counts",500, 0.,5.,nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFmatchedExpTimeKaVsPNoTRDNeg->Sumw2() ;
  if (fEnableAdvancedCheck)
    fHlistPID->AddLast(hTOFmatchedExpTimeKaVsPNoTRDNeg) ;
    
 //PID 26
  TH2F* hTOFmatchedExpTimeProVsPTRDPos = new TH2F("hTOFmatchedExpTimeProVsPTRDPos", "ESDs t_{TOF}-t_{p,exp} (from tracking) Vs P (#phi_{outerTPC}#leq 30 || 150 #leq#phi_{outerTPC}#leq 230 || #phi_{outerTPC}#geq 310 ) ; p (GeV/c);t_{TOF}-t_{p,exp} [ps];Counts",500, 0.,5.,nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFmatchedExpTimeProVsPTRDPos->Sumw2() ;
  if (fEnableAdvancedCheck)
    fHlistPID->AddLast(hTOFmatchedExpTimeProVsPTRDPos) ;
  
   //PID 27
  TH2F* hTOFmatchedExpTimeProVsPNoTRDPos = new TH2F("hTOFmatchedExpTimeProVsPNoTRDPos", "ESDs t_{TOF}-t_{p,exp} (from tracking) Vs P (50 #leq#phi_{outerTPC}#leq 130 || 250 #leq#phi_{outerTPC}#leq 290); p (GeV/c);t_{TOF}-t_{p,exp} [ps];Counts",500, 0.,5.,nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFmatchedExpTimeProVsPNoTRDPos->Sumw2() ;
    if (fEnableAdvancedCheck)
      fHlistPID->AddLast(hTOFmatchedExpTimeProVsPNoTRDPos) ;
    
    //PID 28
   TH2F* hTOFmatchedExpTimeProVsPTRDNeg = new TH2F("hTOFmatchedExpTimeProVsPTRDNeg", "ESDs t_{TOF}-t_{p,exp} (from tracking) Vs P (#phi_{outerTPC}#leq 30 || 150 #leq#phi_{outerTPC}#leq 230 || #phi_{outerTPC}#geq 310 ) ; p (GeV/c);t_{TOF}-t_{p,exp} [ps];Counts",500, 0.,5.,nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
    hTOFmatchedExpTimeProVsPTRDNeg->Sumw2() ;
    if (fEnableAdvancedCheck)
      fHlistPID->AddLast(hTOFmatchedExpTimeProVsPTRDNeg) ;

   //PID 29
  TH2F* hTOFmatchedExpTimeProVsPNoTRDNeg = new TH2F("hTOFmatchedExpTimeProVsPNoTRDNeg", "ESDs t_{TOF}-t_{p,exp} (from tracking) Vs P (50 #leq#phi_{outerTPC}#leq 130 || 250 #leq#phi_{outerTPC}#leq 290); p (GeV/c);t_{TOF}-t_{p,exp} [ps];Counts",500, 0.,5.,nExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  hTOFmatchedExpTimeProVsPNoTRDNeg->Sumw2() ;
  if (fEnableAdvancedCheck)
    fHlistPID->AddLast(hTOFmatchedExpTimeProVsPNoTRDNeg) ;

  //----------------------------------------------------------POSITIVE TRACKS
  //0
  TH1F* hTOFmatchedESDtrkLengthPos  = new TH1F("hTOFmatchedESDtrkLengthPos", "Matched positive ESDs tracks length; Track length [cm];Counts", 1600, -800., 800) ; 
  hTOFmatchedESDtrkLengthPos->Sumw2();
  hTOFmatchedESDtrkLengthPos->SetLineColor(kRed);
  hTOFmatchedESDtrkLengthPos->SetMarkerColor(kRed);
  hTOFmatchedESDtrkLengthPos->SetFillColor(kRed);
  hTOFmatchedESDtrkLengthPos->SetDrawOption("BAR"); 
  if (fEnableAdvancedCheck)
    fHpos->AddLast(hTOFmatchedESDtrkLengthPos);
  //1
  TH1F* hTOFmatchedESDPPos  = new TH1F("hTOFmatchedESDPPos", "TPC-TOF matched positive tracks momentum distribution (GeV/c); p (GeV/c);tracks", 500,0.,5.) ;  
  hTOFmatchedESDPPos->Sumw2() ;
  hTOFmatchedESDPPos->SetLineColor(kRed);
  hTOFmatchedESDPPos->SetMarkerStyle(20);
  hTOFmatchedESDPPos->SetMarkerSize(0.7);
  hTOFmatchedESDPPos->SetMarkerColor(kRed);
  if (fEnableAdvancedCheck)
    fHpos->AddLast(hTOFmatchedESDPPos) ; 
  //2
  TH1F* hTOFmatchedESDPtPos  = new TH1F("hTOFmatchedESDPtPos", "TPC-TOF positive matched tracks p_{T} distribution (GeV/c); p_{T}(GeV/c);tracks", 500,0.,5.) ;  
  hTOFmatchedESDPtPos->Sumw2() ;
  hTOFmatchedESDPtPos->SetLineColor(kRed);
  hTOFmatchedESDPtPos->SetMarkerStyle(21);
  hTOFmatchedESDPtPos->SetMarkerSize(0.7);
  hTOFmatchedESDPtPos->SetMarkerColor(kRed);
  if (fEnableAdvancedCheck)
    fHpos->AddLast(hTOFmatchedESDPtPos) ; 

  //3
  TH1F* hTOFmatchedESDetaPos = new TH1F("hTOFmatchedESDetaPos", "Matched positive ESD tracks #eta (p_{T} #geq 0.5 GeV/c); eta;Counts", 200, -1., 1.) ; 
  hTOFmatchedESDetaPos->Sumw2();
  hTOFmatchedESDetaPos->SetLineColor(kRed);
    if (fEnableAdvancedCheck)
fHpos->AddLast(hTOFmatchedESDetaPos) ; 
  //4
   TH1F* hTOFmatchedESDphiPos = new TH1F("hTOFmatchedESDphiPos", "Matched Positive ESD tracks #phi; #phi (deg);Counts", 72, 0., 360.) ; 
  hTOFmatchedESDphiPos->Sumw2();
  hTOFmatchedESDphiPos->SetLineColor(kRed);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hTOFmatchedESDphiPos) ; 

  //5
  TH1F* hESDmatchedTrackPtPosTRDout = new TH1F("hESDmatchedTrackPtPosTRDout", "ESDs Pos matched tracks p_{T} with kTRDout distribution (GeV/c); p_{T}(GeV/c);tracks", 500, 0.0, 5.0) ;  
  hESDmatchedTrackPtPosTRDout->Sumw2();
  hESDmatchedTrackPtPosTRDout->SetLineWidth(1);
  hESDmatchedTrackPtPosTRDout->SetMarkerStyle(20);
  hESDmatchedTrackPtPosTRDout->SetMarkerSize(0.7);
  hESDmatchedTrackPtPosTRDout->SetLineColor(kRed);
  hESDmatchedTrackPtPosTRDout->SetMarkerColor(kRed);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hESDmatchedTrackPtPosTRDout);

  //6
  TH1F* hESDmatchedTrackPtPosNoTRDout = new TH1F("hESDmatchedTrackPtPosNoTRDout", "ESDs Pos matched tracks p_{T} with !kTRDout distribution (GeV/c); p_{T}(GeV/c);tracks", 500, 0.0, 5.0) ;  
  hESDmatchedTrackPtPosNoTRDout->Sumw2();
  hESDmatchedTrackPtPosNoTRDout->SetLineWidth(1);
  hESDmatchedTrackPtPosNoTRDout->SetMarkerStyle(25);
  hESDmatchedTrackPtPosNoTRDout->SetMarkerSize(0.7);
  hESDmatchedTrackPtPosNoTRDout->SetLineColor(kRed);
  hESDmatchedTrackPtPosNoTRDout->SetMarkerColor(kRed);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hESDmatchedTrackPtPosNoTRDout);

  //7
  TH1F* hESDprimaryTrackPPos = new TH1F("hESDprimaryTrackPPos", "All Pos ESDs tracks p distribution (GeV/c); p(GeV/c);tracks", 500, 0., 5.0) ;  
  hESDprimaryTrackPPos->Sumw2();
  hESDprimaryTrackPPos->SetLineWidth(1);
  hESDprimaryTrackPPos->SetMarkerStyle(24);
  hESDprimaryTrackPPos->SetMarkerSize(0.7);
  hESDprimaryTrackPPos->SetMarkerColor(kRed);
  hESDprimaryTrackPPos->SetLineColor(kRed);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hESDprimaryTrackPPos);
  //8
  TH1F* hESDprimaryTrackPtPos = new TH1F("hESDprimaryTrackPtPos", "ESDs Pos primary tracks p_{T} distribution (GeV/c); p_{T}(GeV/c);tracks", 500, 0.0, 5.0) ;  
  hESDprimaryTrackPtPos->Sumw2();
  hESDprimaryTrackPtPos->SetLineWidth(1);
  hESDprimaryTrackPtPos->SetMarkerStyle(25);
  hESDprimaryTrackPtPos->SetMarkerSize(0.7);
  hESDprimaryTrackPtPos->SetLineColor(kRed);
  hESDprimaryTrackPtPos->SetMarkerColor(kRed);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hESDprimaryTrackPtPos);
  //9
  TH1F* hTOFprimaryESDetaPos = new TH1F("hTOFprimaryESDetaPos", "Primary ESD Pos tracks #eta (p_{T} #geq 0.5 GeV/c); #eta;Counts", 200, -1., 1.) ; 
  hTOFprimaryESDetaPos->Sumw2();
  hTOFprimaryESDetaPos->SetLineColor(kRed);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hTOFprimaryESDetaPos) ; 
  //10
  TH1F* hTOFprimaryESDphiPos = new TH1F("hTOFprimaryESDphiPos", "Primary ESD Pos tracks #phi;#phi (deg);Counts", 72, 0., 360.) ; 
  hTOFprimaryESDphiPos->Sumw2();
  hTOFprimaryESDphiPos->SetLineColor(kRed);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hTOFprimaryESDphiPos) ; 

  //11
  TH1F* hESDprimaryTrackPtPosTRDout = new TH1F("hESDprimaryTrackPtPosTRDout", "ESDs Pos primary tracks p_{T} with kTRDout distribution (GeV/c); p_{T}(GeV/c);tracks", 500, 0.0, 5.0) ;  
  hESDprimaryTrackPtPosTRDout->Sumw2();
  hESDprimaryTrackPtPosTRDout->SetLineWidth(1);
  hESDprimaryTrackPtPosTRDout->SetMarkerStyle(20);
  hESDprimaryTrackPtPosTRDout->SetMarkerSize(0.7);
  hESDprimaryTrackPtPosTRDout->SetLineColor(kRed);
  hESDprimaryTrackPtPosTRDout->SetMarkerColor(kRed);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hESDprimaryTrackPtPosTRDout);

  //12
  TH1F* hESDprimaryTrackPtPosNoTRDout = new TH1F("hESDprimaryTrackPtPosNoTRDout", "ESDs Pos primary tracks p_{T} with !kTRDout distribution (GeV/c); p_{T}(GeV/c);tracks", 500, 0.0, 5.0) ;  
  hESDprimaryTrackPtPosNoTRDout->Sumw2();
  hESDprimaryTrackPtPosNoTRDout->SetLineWidth(1);
  hESDprimaryTrackPtPosNoTRDout->SetMarkerStyle(25);
  hESDprimaryTrackPtPosNoTRDout->SetMarkerSize(0.7);
  hESDprimaryTrackPtPosNoTRDout->SetLineColor(kRed);
  hESDprimaryTrackPtPosNoTRDout->SetMarkerColor(kRed);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hESDprimaryTrackPtPosNoTRDout);

    //13
  TH1F* hESDprimaryTrackEtaPosTRDout = new TH1F("hESDprimaryTrackEtaPosTRDout", "ESDs Pos primary tracks #eta with kTRDout distribution; #eta;tracks", 200, -1.0, 1.0) ;  
  hESDprimaryTrackEtaPosTRDout->Sumw2();
  hESDprimaryTrackEtaPosTRDout->SetLineWidth(1);
  hESDprimaryTrackEtaPosTRDout->SetMarkerStyle(20);
  hESDprimaryTrackEtaPosTRDout->SetMarkerSize(0.7);
  hESDprimaryTrackEtaPosTRDout->SetLineColor(kBlue);
  hESDprimaryTrackEtaPosTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hESDprimaryTrackEtaPosTRDout);

  //14
  TH1F* hESDprimaryTrackEtaPosNoTRDout = new TH1F("hESDprimaryTrackEtaPosNoTRDout", "ESDs Pos primary tracks #eta with !kTRDout distribution;#eta;tracks", 200, -1.0, 1.0) ;  
  hESDprimaryTrackEtaPosNoTRDout->Sumw2();
  hESDprimaryTrackEtaPosNoTRDout->SetLineWidth(1);
  hESDprimaryTrackEtaPosNoTRDout->SetMarkerStyle(25);
  hESDprimaryTrackEtaPosNoTRDout->SetMarkerSize(0.7);
  hESDprimaryTrackEtaPosNoTRDout->SetLineColor(kBlue);
  hESDprimaryTrackEtaPosNoTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hESDprimaryTrackEtaPosNoTRDout);

  //15
  TH1F* hESDprimaryTrackPhiPosTRDout = new TH1F("hESDprimaryTrackPhiPosTRDout", "ESDs Pos primary tracks #phi with kTRDout distribution; #phi(deg);tracks", 72, 0., 360.) ;  
  hESDprimaryTrackPhiPosTRDout->Sumw2();
  hESDprimaryTrackPhiPosTRDout->SetLineWidth(1);
  hESDprimaryTrackPhiPosTRDout->SetMarkerStyle(20);
  hESDprimaryTrackPhiPosTRDout->SetMarkerSize(0.7);
  hESDprimaryTrackPhiPosTRDout->SetLineColor(kBlue);
  hESDprimaryTrackPhiPosTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hESDprimaryTrackPhiPosTRDout);

  //16
  TH1F* hESDprimaryTrackPhiPosNoTRDout = new TH1F("hESDprimaryTrackPhiPosNoTRDout", "ESDs Pos primary tracks #phi with !kTRDout distribution; #phi(deg);tracks", 72, 0., 360.) ;  
  hESDprimaryTrackPhiPosNoTRDout->Sumw2();
  hESDprimaryTrackPhiPosNoTRDout->SetLineWidth(1);
  hESDprimaryTrackPhiPosNoTRDout->SetMarkerStyle(25);
  hESDprimaryTrackPhiPosNoTRDout->SetMarkerSize(0.7);
  hESDprimaryTrackPhiPosNoTRDout->SetLineColor(kBlue);
  hESDprimaryTrackPhiPosNoTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hESDprimaryTrackPhiPosNoTRDout);

  //POS 17
  TH2F* hTOFmatchedExpTimeProVsPhiPos = new TH2F("hTOFmatchedExpTimeProVsPhiPos", "ESDs t_{TOF}-t_{p,exp} vs #phi at TOF (r=378.cm) (p_{T}<1GeV/c);#phi (deg); t_{TOF}-t_{p,exp} [ps];Counts",72, 0.,360., 500, -5000., 5000. ) ; 
  hTOFmatchedExpTimeProVsPhiPos->Sumw2() ;
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hTOFmatchedExpTimeProVsPhiPos) ;

  //POS 18
  TH2F* hTOFmatchedExpTimeProVsPhiVertexPos = new TH2F("hTOFmatchedExpTimeProVsPhiVertexPos", "ESDs t_{TOF}-t_{p,exp} vs #phi at vtx (p_{T}<1GeV/c);#phi (deg); t_{TOF}-t_{p,exp} [ps];Counts",72, 0.,360., nExpTimeSmallBins, fExpTimeSmallRangeMin, fExpTimeSmallRangeMax) ; 
  hTOFmatchedExpTimeProVsPhiVertexPos->Sumw2() ;
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hTOFmatchedExpTimeProVsPhiVertexPos) ;

  //POS 19
  TH2F* hTOFmatchedExpTimeProVsPhiTPCPos = new TH2F("hTOFmatchedExpTimeProVsPhiTPCPos", "ESDs t_{TOF}-t_{p,exp} vs #phi at outer TPC (p_{T}<1GeV/c);#phi (deg); t_{TOF}-t_{p,exp} [ps];Counts",72, 0.,360., nExpTimeSmallBins, fExpTimeSmallRangeMin, fExpTimeSmallRangeMax) ; 
  hTOFmatchedExpTimeProVsPhiTPCPos->Sumw2() ;
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hTOFmatchedExpTimeProVsPhiTPCPos) ;

   //20
  TH1F* hESDmatchedTrackEtaPosTRDout = new TH1F("hESDmatchedTrackEtaPosTRDout", "ESDs Pos matched tracks #eta with kTRDout distribution; #eta;tracks", 200, -1.0, 1.0) ;  
  hESDmatchedTrackEtaPosTRDout->Sumw2();
  hESDmatchedTrackEtaPosTRDout->SetLineWidth(1);
  hESDmatchedTrackEtaPosTRDout->SetMarkerStyle(20);
  hESDmatchedTrackEtaPosTRDout->SetMarkerSize(0.7);
  hESDmatchedTrackEtaPosTRDout->SetLineColor(kBlue);
  hESDmatchedTrackEtaPosTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hESDmatchedTrackEtaPosTRDout);

  //21
  TH1F* hESDmatchedTrackEtaPosNoTRDout = new TH1F("hESDmatchedTrackEtaPosNoTRDout", "ESDs Pos matched tracks #eta with !kTRDout distribution;#eta;tracks", 200, -1.0, 1.0) ;  
  hESDmatchedTrackEtaPosNoTRDout->Sumw2();
  hESDmatchedTrackEtaPosNoTRDout->SetLineWidth(1);
  hESDmatchedTrackEtaPosNoTRDout->SetMarkerStyle(25);
  hESDmatchedTrackEtaPosNoTRDout->SetMarkerSize(0.7);
  hESDmatchedTrackEtaPosNoTRDout->SetLineColor(kBlue);
  hESDmatchedTrackEtaPosNoTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hESDmatchedTrackEtaPosNoTRDout);

  //22
  TH1F* hESDmatchedTrackPhiPosTRDout = new TH1F("hESDmatchedTrackPhiPosTRDout", "ESDs Pos matched tracks #phi with kTRDout distribution; #phi(deg);tracks", 72, 0., 360.) ;  
  hESDmatchedTrackPhiPosTRDout->Sumw2();
  hESDmatchedTrackPhiPosTRDout->SetLineWidth(1);
  hESDmatchedTrackPhiPosTRDout->SetMarkerStyle(20);
  hESDmatchedTrackPhiPosTRDout->SetMarkerSize(0.7);
  hESDmatchedTrackPhiPosTRDout->SetLineColor(kBlue);
  hESDmatchedTrackPhiPosTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hESDmatchedTrackPhiPosTRDout);

  //23
  TH1F* hESDmatchedTrackPhiPosNoTRDout = new TH1F("hESDmatchedTrackPhiPosNoTRDout", "ESDs Pos matched tracks #phi with !kTRDout distribution; #phi(deg);tracks", 72, 0., 360.) ;  
  hESDmatchedTrackPhiPosNoTRDout->Sumw2();
  hESDmatchedTrackPhiPosNoTRDout->SetLineWidth(1);
  hESDmatchedTrackPhiPosNoTRDout->SetMarkerStyle(25);
  hESDmatchedTrackPhiPosNoTRDout->SetMarkerSize(0.7);
  hESDmatchedTrackPhiPosNoTRDout->SetLineColor(kBlue);
  hESDmatchedTrackPhiPosNoTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hESDmatchedTrackPhiPosNoTRDout);

  //POS 24
  TH2F* hTOFmatchedExpTimePiVsPhiTPCPos = new TH2F("hTOFmatchedExpTimePiVsPhiTPCPos", "ESDs t_{TOF}-t_{#pi,exp} vs #phi at outer TPC (p_{T}<1GeV/c);#phi (deg); t_{TOF}-t_{p,exp} [ps];Counts",72, 0.,360., nExpTimeSmallBins, fExpTimeSmallRangeMin, fExpTimeSmallRangeMax) ; 
  hTOFmatchedExpTimePiVsPhiTPCPos->Sumw2() ;
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hTOFmatchedExpTimePiVsPhiTPCPos) ;

  //POS 25
  TH2F* hTOFmatchedExpTimeKaVsPhiTPCPos = new TH2F("hTOFmatchedExpTimeKaVsPhiTPCPos", "ESDs t_{TOF}-t_{K,exp} vs #phi at outer TPC (p_{T}<1GeV/c);#phi (deg); t_{TOF}-t_{p,exp} [ps];Counts",72, 0.,360., nExpTimeSmallBins, fExpTimeSmallRangeMin, fExpTimeSmallRangeMax) ; 
  hTOFmatchedExpTimeKaVsPhiTPCPos->Sumw2() ;
  if (fEnableAdvancedCheck)
  fHpos->AddLast(hTOFmatchedExpTimeKaVsPhiTPCPos) ;

//----------------------------------------------------------NEGATIVE TRACKS
  //0
  TH1F* hTOFmatchedESDtrkLengthNeg  = new TH1F("hTOFmatchedESDtrkLengthNeg", "Matched positive ESDs tracks length; Track length [cm];Counts", 1600, -800., 800) ; 
  hTOFmatchedESDtrkLengthNeg->Sumw2();
  hTOFmatchedESDtrkLengthNeg->SetLineColor(kBlue);
  hTOFmatchedESDtrkLengthNeg->SetMarkerColor(kBlue);
  hTOFmatchedESDtrkLengthNeg->SetFillColor(kBlue);
  hTOFmatchedESDtrkLengthNeg->SetDrawOption("BAR"); 
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hTOFmatchedESDtrkLengthNeg);
  //1
  TH1F* hTOFmatchedESDPNeg  = new TH1F("hTOFmatchedESDPNeg", "TPC-TOF matched positive tracks momentum distribution (GeV/c); p (GeV/c);tracks", 500,0.,5.) ;  
  hTOFmatchedESDPNeg->Sumw2() ;
  hTOFmatchedESDPNeg->SetLineColor(kBlue);
  hTOFmatchedESDPNeg->SetMarkerStyle(20);
  hTOFmatchedESDPNeg->SetMarkerSize(0.7);
  hTOFmatchedESDPNeg->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hTOFmatchedESDPNeg) ; 
  //2
  TH1F* hTOFmatchedESDPtNeg  = new TH1F("hTOFmatchedESDPtNeg", "TPC-TOF positive matched tracks p_{T} distribution (GeV/c); p_{T}(GeV/c);tracks", 500,0.,5.) ;  
  hTOFmatchedESDPtNeg->Sumw2() ;
  hTOFmatchedESDPtNeg->SetLineColor(kBlue);
  hTOFmatchedESDPtNeg->SetMarkerStyle(21);
  hTOFmatchedESDPtNeg->SetMarkerSize(0.7);
  hTOFmatchedESDPtNeg->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hTOFmatchedESDPtNeg) ; 

  //3
  TH1F* hTOFmatchedESDetaNeg = new TH1F("hTOFmatchedESDetaNeg", "Matched positive ESD tracks #eta (p_{T} #geq 0.5 GeV/c); eta;Counts", 200, -1., 1.) ; 
  hTOFmatchedESDetaNeg->Sumw2();
  hTOFmatchedESDetaNeg->SetLineColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hTOFmatchedESDetaNeg) ; 
  //4
   TH1F* hTOFmatchedESDphiNeg = new TH1F("hTOFmatchedESDphiNeg", "Matched Negitive ESD tracks #phi; #phi (deg);Counts", 72, 0., 360.) ; 
  hTOFmatchedESDphiNeg->Sumw2();
  hTOFmatchedESDphiNeg->SetLineColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hTOFmatchedESDphiNeg) ; 

  //5
  TH1F* hESDmatchedTrackPtNegTRDout = new TH1F("hESDmatchedTrackPtNegTRDout", "ESDs Neg matched tracks p_{T} with kTRDout distribution (GeV/c); p_{T}(GeV/c);tracks", 500, 0.0, 5.0) ;  
  hESDmatchedTrackPtNegTRDout->Sumw2();
  hESDmatchedTrackPtNegTRDout->SetLineWidth(1);
  hESDmatchedTrackPtNegTRDout->SetMarkerStyle(20);
  hESDmatchedTrackPtNegTRDout->SetMarkerSize(0.7);
  hESDmatchedTrackPtNegTRDout->SetLineColor(kBlue);
  hESDmatchedTrackPtNegTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hESDmatchedTrackPtNegTRDout);

  //6
  TH1F* hESDmatchedTrackPtNegNoTRDout = new TH1F("hESDmatchedTrackPtNegNoTRDout", "ESDs Neg matched tracks p_{T} with !kTRDout distribution (GeV/c); p_{T}(GeV/c);tracks", 500, 0.0, 5.0) ;  
  hESDmatchedTrackPtNegNoTRDout->Sumw2();
  hESDmatchedTrackPtNegNoTRDout->SetLineWidth(1);
  hESDmatchedTrackPtNegNoTRDout->SetMarkerStyle(25);
  hESDmatchedTrackPtNegNoTRDout->SetMarkerSize(0.7);
  hESDmatchedTrackPtNegNoTRDout->SetLineColor(kBlue);
  hESDmatchedTrackPtNegNoTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hESDmatchedTrackPtNegNoTRDout);

  //7
  TH1F* hESDprimaryTrackPNeg = new TH1F("hESDprimaryTrackPNeg", "All Neg ESDs tracks p distribution (GeV/c); p(GeV/c);tracks", 500, 0., 5.0) ;  
  hESDprimaryTrackPNeg->Sumw2();
  hESDprimaryTrackPNeg->SetLineWidth(1);
  hESDprimaryTrackPNeg->SetMarkerStyle(24);
  hESDprimaryTrackPNeg->SetMarkerSize(0.7);
  hESDprimaryTrackPNeg->SetMarkerColor(kBlue);
  hESDprimaryTrackPNeg->SetLineColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hESDprimaryTrackPNeg);
  //8
  TH1F* hESDprimaryTrackPtNeg = new TH1F("hESDprimaryTrackPtNeg", "ESDs Neg primary tracks p_{T} distribution (GeV/c); p_{T}(GeV/c);tracks", 500, 0.0, 5.0) ;  
  hESDprimaryTrackPtNeg->Sumw2();
  hESDprimaryTrackPtNeg->SetLineWidth(1);
  hESDprimaryTrackPtNeg->SetMarkerStyle(25);
  hESDprimaryTrackPtNeg->SetMarkerSize(0.7);
  hESDprimaryTrackPtNeg->SetLineColor(kBlue);
  hESDprimaryTrackPtNeg->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hESDprimaryTrackPtNeg);
  //9
  TH1F* hTOFprimaryESDetaNeg = new TH1F("hTOFprimaryESDetaNeg", "Primary ESD Neg tracks #eta (p_{T} #geq 0.5 GeV/c); #eta;Counts", 200, -1., 1.) ; 
  hTOFprimaryESDetaNeg->Sumw2();
  hTOFprimaryESDetaNeg->SetLineColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hTOFprimaryESDetaNeg) ; 
  //10
  TH1F* hTOFprimaryESDphiNeg = new TH1F("hTOFprimaryESDphiNeg", "Primary ESD Neg tracks #phi;#phi (deg);Counts", 72, 0., 360.) ; 
  hTOFprimaryESDphiNeg->Sumw2();
  hTOFprimaryESDphiNeg->SetLineColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hTOFprimaryESDphiNeg) ; 

 //11
  TH1F* hESDprimaryTrackPtNegTRDout = new TH1F("hESDprimaryTrackPtNegTRDout", "ESDs Neg primary tracks p_{T} with kTRDout distribution (GeV/c); p_{T}(GeV/c);tracks", 500, 0.0, 5.0) ;  
  hESDprimaryTrackPtNegTRDout->Sumw2();
  hESDprimaryTrackPtNegTRDout->SetLineWidth(1);
  hESDprimaryTrackPtNegTRDout->SetMarkerStyle(20);
  hESDprimaryTrackPtNegTRDout->SetMarkerSize(0.7);
  hESDprimaryTrackPtNegTRDout->SetLineColor(kBlue);
  hESDprimaryTrackPtNegTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hESDprimaryTrackPtNegTRDout);

  //12
  TH1F* hESDprimaryTrackPtNegNoTRDout = new TH1F("hESDprimaryTrackPtNegNoTRDout", "ESDs Neg primary tracks p_{T} with !kTRDout distribution (GeV/c); p_{T}(GeV/c);tracks", 500, 0.0, 5.0) ;  
  hESDprimaryTrackPtNegNoTRDout->Sumw2();
  hESDprimaryTrackPtNegNoTRDout->SetLineWidth(1);
  hESDprimaryTrackPtNegNoTRDout->SetMarkerStyle(25);
  hESDprimaryTrackPtNegNoTRDout->SetMarkerSize(0.7);
  hESDprimaryTrackPtNegNoTRDout->SetLineColor(kBlue);
  hESDprimaryTrackPtNegNoTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hESDprimaryTrackPtNegNoTRDout);

  //13
  TH1F* hESDprimaryTrackEtaNegTRDout = new TH1F("hESDprimaryTrackEtaNegTRDout", "ESDs Neg primary tracks #eta with kTRDout distribution; #eta;tracks", 200, -1.0, 1.0) ;  
  hESDprimaryTrackEtaNegTRDout->Sumw2();
  hESDprimaryTrackEtaNegTRDout->SetLineWidth(1);
  hESDprimaryTrackEtaNegTRDout->SetMarkerStyle(20);
  hESDprimaryTrackEtaNegTRDout->SetMarkerSize(0.7);
  hESDprimaryTrackEtaNegTRDout->SetLineColor(kBlue);
  hESDprimaryTrackEtaNegTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hESDprimaryTrackEtaNegTRDout);

  //14
  TH1F* hESDprimaryTrackEtaNegNoTRDout = new TH1F("hESDprimaryTrackEtaNegNoTRDout", "ESDs Neg primary tracks #eta with !kTRDout distribution;#eta;tracks", 200, -1.0, 1.0) ;  
  hESDprimaryTrackEtaNegNoTRDout->Sumw2();
  hESDprimaryTrackEtaNegNoTRDout->SetLineWidth(1);
  hESDprimaryTrackEtaNegNoTRDout->SetMarkerStyle(25);
  hESDprimaryTrackEtaNegNoTRDout->SetMarkerSize(0.7);
  hESDprimaryTrackEtaNegNoTRDout->SetLineColor(kBlue);
  hESDprimaryTrackEtaNegNoTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hESDprimaryTrackEtaNegNoTRDout);

  //15
  TH1F* hESDprimaryTrackPhiNegTRDout = new TH1F("hESDprimaryTrackPhiNegTRDout", "ESDs Neg primary tracks #phi with kTRDout distribution; #phi(deg);tracks", 72, 0., 360.) ;  
  hESDprimaryTrackPhiNegTRDout->Sumw2();
  hESDprimaryTrackPhiNegTRDout->SetLineWidth(1);
  hESDprimaryTrackPhiNegTRDout->SetMarkerStyle(20);
  hESDprimaryTrackPhiNegTRDout->SetMarkerSize(0.7);
  hESDprimaryTrackPhiNegTRDout->SetLineColor(kBlue);
  hESDprimaryTrackPhiNegTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hESDprimaryTrackPhiNegTRDout);

  //16
  TH1F* hESDprimaryTrackPhiNegNoTRDout = new TH1F("hESDprimaryTrackPhiNegNoTRDout", "ESDs Neg primary tracks #phi with !kTRDout distribution; #phi(deg);tracks", 72, 0., 360.) ;  
  hESDprimaryTrackPhiNegNoTRDout->Sumw2();
  hESDprimaryTrackPhiNegNoTRDout->SetLineWidth(1);
  hESDprimaryTrackPhiNegNoTRDout->SetMarkerStyle(25);
  hESDprimaryTrackPhiNegNoTRDout->SetMarkerSize(0.7);
  hESDprimaryTrackPhiNegNoTRDout->SetLineColor(kBlue);
  hESDprimaryTrackPhiNegNoTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hESDprimaryTrackPhiNegNoTRDout);

  //NEG 17
  TH2F* hTOFmatchedExpTimeProVsPhiNeg = new TH2F("hTOFmatchedExpTimeProVsPhiNeg", "ESDs t_{TOF}-t_{p,exp} vs #phi at TOF (r=378.cm) (p_{T}<1GeV/c);#phi (deg); t_{TOF}-t_{p,exp} [ps];Counts",72, 0.,360., nExpTimeSmallBins, fExpTimeSmallRangeMin, fExpTimeSmallRangeMax) ; 
  hTOFmatchedExpTimeProVsPhiNeg->Sumw2() ;
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hTOFmatchedExpTimeProVsPhiNeg) ;

  //NEG 18
  TH2F* hTOFmatchedExpTimeProVsPhiVertexNeg = new TH2F("hTOFmatchedExpTimeProVsPhiVertexNeg", "ESDs t_{TOF}-t_{p,exp} vs #phi at vtx (p_{T}<1GeV/c);#phi (deg); t_{TOF}-t_{p,exp} [ps];Counts",72, 0.,360., nExpTimeSmallBins, fExpTimeSmallRangeMin, fExpTimeSmallRangeMax) ; 
  hTOFmatchedExpTimeProVsPhiVertexNeg->Sumw2() ;
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hTOFmatchedExpTimeProVsPhiVertexNeg) ;

  //NEG 19
  TH2F* hTOFmatchedExpTimeProVsPhiTPCNeg = new TH2F("hTOFmatchedExpTimeProVsPhiTPCNeg", "ESDs t_{TOF}-t_{p,exp} vs #phi at outer TPC (p_{T}<1GeV/c);#phi (deg); t_{TOF}-t_{p,exp} [ps];Counts",72, 0.,360., nExpTimeSmallBins, fExpTimeSmallRangeMin, fExpTimeSmallRangeMax) ; 
  hTOFmatchedExpTimeProVsPhiTPCNeg->Sumw2() ;
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hTOFmatchedExpTimeProVsPhiTPCNeg) ;


  //20
  TH1F* hESDmatchedTrackEtaNegTRDout = new TH1F("hESDmatchedTrackEtaNegTRDout", "ESDs Neg matched tracks #eta with kTRDout distribution; #eta;tracks", 200, -1.0, 1.0) ;  
  hESDmatchedTrackEtaNegTRDout->Sumw2();
  hESDmatchedTrackEtaNegTRDout->SetLineWidth(1);
  hESDmatchedTrackEtaNegTRDout->SetMarkerStyle(20);
  hESDmatchedTrackEtaNegTRDout->SetMarkerSize(0.7);
  hESDmatchedTrackEtaNegTRDout->SetLineColor(kBlue);
  hESDmatchedTrackEtaNegTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hESDmatchedTrackEtaNegTRDout);

  //21
  TH1F* hESDmatchedTrackEtaNegNoTRDout = new TH1F("hESDmatchedTrackEtaNegNoTRDout", "ESDs Neg matched tracks #eta with !kTRDout distribution;#eta;tracks", 200, -1.0, 1.0) ;  
  hESDmatchedTrackEtaNegNoTRDout->Sumw2();
  hESDmatchedTrackEtaNegNoTRDout->SetLineWidth(1);
  hESDmatchedTrackEtaNegNoTRDout->SetMarkerStyle(25);
  hESDmatchedTrackEtaNegNoTRDout->SetMarkerSize(0.7);
  hESDmatchedTrackEtaNegNoTRDout->SetLineColor(kBlue);
  hESDmatchedTrackEtaNegNoTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hESDmatchedTrackEtaNegNoTRDout);

  //22
  TH1F* hESDmatchedTrackPhiNegTRDout = new TH1F("hESDmatchedTrackPhiNegTRDout", "ESDs Neg matched tracks #phi with kTRDout distribution; #phi(deg);tracks", 72, 0., 360.) ;  
  hESDmatchedTrackPhiNegTRDout->Sumw2();
  hESDmatchedTrackPhiNegTRDout->SetLineWidth(1);
  hESDmatchedTrackPhiNegTRDout->SetMarkerStyle(20);
  hESDmatchedTrackPhiNegTRDout->SetMarkerSize(0.7);
  hESDmatchedTrackPhiNegTRDout->SetLineColor(kBlue);
  hESDmatchedTrackPhiNegTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hESDmatchedTrackPhiNegTRDout);

  //23
  TH1F* hESDmatchedTrackPhiNegNoTRDout = new TH1F("hESDmatchedTrackPhiNegNoTRDout", "ESDs Neg matched tracks #phi with !kTRDout distribution; #phi(deg);tracks", 72, 0., 360.) ;  
  hESDmatchedTrackPhiNegNoTRDout->Sumw2();
  hESDmatchedTrackPhiNegNoTRDout->SetLineWidth(1);
  hESDmatchedTrackPhiNegNoTRDout->SetMarkerStyle(25);
  hESDmatchedTrackPhiNegNoTRDout->SetMarkerSize(0.7);
  hESDmatchedTrackPhiNegNoTRDout->SetLineColor(kBlue);
  hESDmatchedTrackPhiNegNoTRDout->SetMarkerColor(kBlue);
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hESDmatchedTrackPhiNegNoTRDout);
 
 //NEG 24
  TH2F* hTOFmatchedExpTimePiVsPhiTPCNeg = new TH2F("hTOFmatchedExpTimePiVsPhiTPCNeg", "ESDs t_{TOF}-t_{#pi,exp} vs #phi at outer TPC (p_{T}<1GeV/c);#phi (deg); t_{TOF}-t_{p,exp} [ps];Counts",72, 0.,360., nExpTimeSmallBins, fExpTimeSmallRangeMin, fExpTimeSmallRangeMax) ; 
  hTOFmatchedExpTimePiVsPhiTPCNeg->Sumw2() ;
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hTOFmatchedExpTimePiVsPhiTPCNeg) ;

  //NEG 25
  TH2F* hTOFmatchedExpTimeKaVsPhiTPCNeg = new TH2F("hTOFmatchedExpTimeKaVsPhiTPCNeg", "ESDs t_{TOF}-t_{K,exp} vs #phi at outer TPC (p_{T}<1GeV/c);#phi (deg); t_{TOF}-t_{p,exp} [ps];Counts",72, 0.,360., nExpTimeSmallBins, fExpTimeSmallRangeMin, fExpTimeSmallRangeMax) ; 
  hTOFmatchedExpTimeKaVsPhiTPCNeg->Sumw2() ;
  if (fEnableAdvancedCheck)
  fHneg->AddLast(hTOFmatchedExpTimeKaVsPhiTPCNeg) ;

  PostData(1, fHlist);
  PostData(2, fHlistTimeZero);
  PostData(3, fHlistPID);
  PostData(4, fHpos);
  PostData(5, fHneg);
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

  // info from V0 detector QA 
  AliESDVZERO * vzero = fESD->GetVZEROData();
  Float_t V0Atime = vzero->GetV0ATime();
  Float_t V0Ctime = vzero->GetV0CTime(); 
  ((TH2F*)fHlistTimeZero->FindObject("hEventV0MeanVsVtx"))->Fill((V0Atime-V0Ctime)*0.5,(V0Atime+V0Ctime)*0.5);

  
  // info from T0 detector QA 
  for (Int_t j=0;j<3;j++){
    fT0[j]= (Float_t) fESD->GetT0TOF(j);//ps
    if (fT0[j]>90000.) fT0[j]=99999.;//fix old default values to the new one
  }
  Float_t t0cut = 90000.; 
  //Float_t t0cut =3 * t0spread; //use this cut to check t0 used in tof response
  // if(t0cut < 500) t0cut = 500;
  
  if(TMath::Abs(fT0[1]) < t0cut && TMath::Abs(fT0[2]) < t0cut ) {
    //&& TMath::Abs(fT0[2]-fT0[1]) < 500)  //add this condition to check t0 used in tof response
    ((TH1F*)fHlistTimeZero->FindObject("hT0DetRes"))->Fill((fT0[2]-fT0[1])*0.5);
    ((TH1F*)fHlistTimeZero->FindObject("hEventT0DetAND"))->Fill(fT0[0]);  
    ((TH2F*)fHlistTimeZero->FindObject("hEventT0MeanVsVtx"))->Fill((fT0[2]-fT0[1])*0.5,(fT0[2]+fT0[1])*0.5);
  } 
  if(TMath::Abs(fT0[1]) < t0cut){
    ((TH1F*)fHlistTimeZero->FindObject("hEventT0DetA"))->Fill(fT0[1]);   
  }
  if(TMath::Abs(fT0[2]) < t0cut){
	  ((TH1F*)fHlistTimeZero->FindObject("hEventT0DetC"))->Fill(fT0[2]);
  }
  
  //  event timeZero QA via AliESDpid::SetTOFResponse() 
  Double_t timeZero[4]={99999.,99999.,99999.,99999.};
  Double_t timeZeroRes[4]={99999.,99999.,99999.,99999.}; 
  
  TString timeZeroHisto[4]={"hT0fill","hT0TOF","hT0T0","hT0best"};
  TString timeZeroHistoRes[4]={"hT0fillRes","hT0TOFRes","hT0T0Res","hT0bestRes"};
  for (Int_t j=0;j<4;j++){
    fESDpid->SetTOFResponse(fESD, (AliESDpid::EStartTimeType_t) j);//(fill_t0, tof_t0, t0_t0, best_t0)
    timeZero[j]=fESDpid->GetTOFResponse().GetStartTime(10.); //timeZero for bin pT>10GeV/c
    timeZeroRes[j]=fESDpid->GetTOFResponse().GetStartTimeRes(10.); //timeZero for bin pT>10GeV/c
    ((TH1D*)(fHlistTimeZero->FindObject(timeZeroHisto[j].Data())))->Fill(timeZero[j]);
    ((TH1D*)(fHlistTimeZero->FindObject(timeZeroHistoRes[j].Data())))->Fill(timeZeroRes[j]);
  }

  //re-set response tof_t0 for all other checks
  fESDpid->SetTOFResponse(fESD,AliESDpid::kTOF_T0);//(fill_t0, tof_t0, t0_t0, best_t0)
  
  // loop over ESD tracks 
  fNTOFtracks=0;
  // fNPrimaryTracks=0;

  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }

    //primary tracks selection: kTPCrefit and std cuts
    if (fTrackFilter){
	    if(!fTrackFilter->IsSelected(track)) continue;
    }
    else{
	    Printf("No track filter found, skipping the track loop");
	    break;
    }

    Double_t eta=track->Eta();
    if (TMath::Abs(eta)>0.8) continue; //cut for acceptance

    Double_t mom=track->P();
    Double_t mom2 = mom*mom;
    Double_t length=track->GetIntegratedLength();
    Double_t pT = track->Pt();
    Double_t phi=track->Phi()*TMath::RadToDeg();
    track->GetIntegratedTimes(fTrkExpTimes);
      
    ((TH1F*)fHlist->FindObject("hESDprimaryTrackP"))->Fill(mom); 
    ((TH1F*)fHlist->FindObject("hESDprimaryTrackPt"))->Fill(pT); 
    if (pT>=0.5)
      ((TH1F*)fHlist->FindObject("hTOFprimaryESDeta"))->Fill(eta);
    ((TH1F*)fHlist->FindObject("hTOFprimaryESDphi"))->Fill(phi);
    
     //evaluate sign  
     if (fEnableAdvancedCheck){
       if (track->GetSign()>0){
 	((TH1F*)fHpos->FindObject("hESDprimaryTrackPPos"))->Fill(mom); 
 	((TH1F*)fHpos->FindObject("hESDprimaryTrackPtPos"))->Fill(pT);	
 	if (pT>=0.5){
 	  ((TH1F*)fHpos->FindObject("hTOFprimaryESDetaPos"))->Fill(eta);
 	  ((TH1F*)fHpos->FindObject("hTOFprimaryESDphiPos"))->Fill(phi);
 	}
 	if ( (track->IsOn(AliESDtrack::kTRDout)) && (track->IsOn(AliESDtrack::kTPCout)) ) {
 	  ((TH1F*)fHpos->FindObject("hESDprimaryTrackPtPosTRDout"))->Fill(pT);
 	  if (pT>=0.5) {
 	    ((TH1F*)fHpos->FindObject("hESDprimaryTrackEtaPosTRDout"))->Fill(eta);
 	    ((TH1F*)fHpos->FindObject("hESDprimaryTrackPhiPosTRDout"))->Fill(phi);
 	  }
 	}
 	if ((!(track->IsOn(AliESDtrack::kTRDout)))  && (track->IsOn(AliESDtrack::kTPCout)) ) {
 	  ((TH1F*)fHpos->FindObject("hESDprimaryTrackPtPosNoTRDout"))->Fill(pT);
 	  if (pT>=0.5) {
 	    ((TH1F*)fHpos->FindObject("hESDprimaryTrackEtaPosNoTRDout"))->Fill(eta);
 	  ((TH1F*)fHpos->FindObject("hESDprimaryTrackPhiPosNoTRDout"))->Fill(phi);
 	  }
 	}
       } //end positive 
       else {
	       ((TH1F*)fHneg->FindObject("hESDprimaryTrackPNeg"))->Fill(mom); 
 	((TH1F*)fHneg->FindObject("hESDprimaryTrackPtNeg"))->Fill(pT); 
 	if (pT>=0.5){
 	  ((TH1F*)fHneg->FindObject("hTOFprimaryESDetaNeg"))->Fill(eta);
 	  ((TH1F*)fHneg->FindObject("hTOFprimaryESDphiNeg"))->Fill(phi);		
 	}
 	if ((track->IsOn(AliESDtrack::kTRDout))  && (track->IsOn(AliESDtrack::kTPCout)) ) {
 	  ((TH1F*)fHneg->FindObject("hESDprimaryTrackPtNegTRDout"))->Fill(pT);
 	  if (pT>0.5) {
 	    ((TH1F*)fHneg->FindObject("hESDprimaryTrackEtaNegTRDout"))->Fill(eta);
 	    ((TH1F*)fHneg->FindObject("hESDprimaryTrackPhiNegTRDout"))->Fill(phi); 
 	  }
 	}
 	if  ((!(track->IsOn(AliESDtrack::kTRDout)))  && (track->IsOn(AliESDtrack::kTPCout)) ) {
 	  ((TH1F*)fHneg->FindObject("hESDprimaryTrackPtNegNoTRDout"))->Fill(pT);
 	  if (pT>0.5){
 	    ((TH1F*)fHneg->FindObject("hESDprimaryTrackEtaNegNoTRDout"))->Fill(eta);
 	    ((TH1F*)fHneg->FindObject("hESDprimaryTrackPhiNegNoTRDout"))->Fill(phi);
 	  }
 	}
       }//end negative
     }//end flag advanced check

    //matched tracks selection: kTOFout and kTIME
    if ( (track->IsOn(AliESDtrack::kTOFout)) &&
	 (track->IsOn(AliESDtrack::kTIME)) &&
	 (track->IsOn(AliESDtrack::kTPCout))  ) {      
      
      Double_t tofTime=track->GetTOFsignal();//in ps
      Double_t tofTimeRaw=track->GetTOFsignalRaw();//in ps
      Double_t tofToT=track->GetTOFsignalToT(); //in ps
      Int_t channel=track->GetTOFCalChannel(); 
      Int_t volId[5]; //(sector, plate,strip,padZ,padX)
      AliTOFGeometry::GetVolumeIndices(channel,volId);
      
      if (pT>=0.3) fNTOFtracks++; //matched counter
      Double_t tof= tofTime*1E-3; // ns, average T0 fill subtracted, no info from T0detector 	 
      ((TH1F*)fHlist->FindObject("hTOFmatchedESDtime"))->Fill(tof); //ns
      ((TH1F*)fHlist->FindObject("hTOFmatchedESDrawTime"))->Fill(tofTimeRaw*1E-3); //ns
      ((TH1F*)fHlist->FindObject("hTOFmatchedESDToT"))->Fill(tofToT);
      ((TH1F*)fHlist->FindObject("hTOFmatchedESDtrkLength"))->Fill(length);  
      ((TH1F*)fHlist->FindObject("hTOFmatchedESDP"))->Fill(mom);
      ((TH1F*)fHlist->FindObject("hTOFmatchedESDPt"))->Fill(pT);
      if (pT>=0.5)
	((TH1F*)fHlist->FindObject("hTOFmatchedESDeta"))->Fill(eta);
      ((TH1F*)fHlist->FindObject("hTOFmatchedESDphi"))->Fill(phi);
      if (track->GetSign()>0)
	((TH2F*)fHlist->FindObject("hTOFmatchedDxVsPtPos"))->Fill(pT,track->GetTOFsignalDx());
      else ((TH2F*)fHlist->FindObject("hTOFmatchedDxVsPtNeg"))->Fill(pT,track->GetTOFsignalDx());
      ((TH2F*)fHlist->FindObject("hTOFmatchedDzVsStrip"))->Fill((Int_t)GetStripIndex(volId),track->GetTOFsignalDz());

      //evaluate sign
      if (fEnableAdvancedCheck){
        if (track->GetSign()>0){
 	  ((TH1F*)fHpos->FindObject("hTOFmatchedESDtrkLengthPos"))->Fill(length);  
 	  ((TH1F*)fHpos->FindObject("hTOFmatchedESDPPos"))->Fill(mom); 
 	  ((TH1F*)fHpos->FindObject("hTOFmatchedESDPtPos"))->Fill(pT); 
 	  if (pT>=0.5){
 	    ((TH1F*)fHpos->FindObject("hTOFmatchedESDetaPos"))->Fill(eta);
 	    ((TH1F*)fHpos->FindObject("hTOFmatchedESDphiPos"))->Fill(phi);
 	  }	
 	  if ( (track->IsOn(AliESDtrack::kTRDout))  && (track->IsOn(AliESDtrack::kTPCout)) ) {
 	    ((TH1F*)fHpos->FindObject("hESDmatchedTrackPtPosTRDout"))->Fill(pT);
 	    if (pT>=0.5){
 	      ((TH1F*)fHpos->FindObject("hESDmatchedTrackEtaPosTRDout"))->Fill(eta);
 	      ((TH1F*)fHpos->FindObject("hESDmatchedTrackPhiPosTRDout"))->Fill(phi);
 	    }
 	  } 
 	  if ((!(track->IsOn(AliESDtrack::kTRDout)))  && (track->IsOn(AliESDtrack::kTPCout)) ) {
 	    ((TH1F*)fHpos->FindObject("hESDmatchedTrackPtPosNoTRDout"))->Fill(pT);
 	    if (pT>=0.5) {
 	      ((TH1F*)fHpos->FindObject("hESDmatchedTrackEtaPosNoTRDout"))->Fill(eta);
 	      ((TH1F*)fHpos->FindObject("hESDmatchedTrackPhiPosNoTRDout"))->Fill(phi);
 	    }
 	  }
       }// END POSITIVE
       else {	
         ((TH1F*)fHneg->FindObject("hTOFmatchedESDtrkLengthNeg"))->Fill(length);  
         ((TH1F*)fHneg->FindObject("hTOFmatchedESDPNeg"))->Fill(mom); 
         ((TH1F*)fHneg->FindObject("hTOFmatchedESDPtNeg"))->Fill(pT); 
         if (pT>=0.5) {
           ((TH1F*)fHneg->FindObject("hTOFmatchedESDetaNeg"))->Fill(eta);
           ((TH1F*)fHneg->FindObject("hTOFmatchedESDphiNeg"))->Fill(phi);
         }	  
         if ((track->IsOn(AliESDtrack::kTRDout))  && (track->IsOn(AliESDtrack::kTPCout)) ) {
           ((TH1F*)fHneg->FindObject("hESDmatchedTrackPtNegTRDout"))->Fill(pT);
           if (pT>=0.5) {
             ((TH1F*)fHneg->FindObject("hESDmatchedTrackEtaNegTRDout"))->Fill(eta);
             ((TH1F*)fHneg->FindObject("hESDmatchedTrackPhiNegTRDout"))->Fill(phi);	    
           }
         }
         if ((!(track->IsOn(AliESDtrack::kTRDout)))  && (track->IsOn(AliESDtrack::kTPCout)) ) {
           ((TH1F*)fHneg->FindObject("hESDmatchedTrackPtNegNoTRDout"))->Fill(pT);
           if (pT>=0.5) {
             ((TH1F*)fHneg->FindObject("hESDmatchedTrackEtaNegNoTRDout"))->Fill(eta);
             ((TH1F*)fHneg->FindObject("hESDmatchedTrackPhiNegNoTRDout"))->Fill(phi);
           }
         }	  
       }//end negative
      }//end advanced check

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
      ((TH2F*)fHlistPID->FindObject("hTOFmatchedESDpVsBeta"))->Fill(mom,beta);
      ((TH1F*) fHlistPID->FindObject("hTOFmatchedMass"))->Fill(mass);
      
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
	((TH2F*)fHlistPID->FindObject("hTOFmatchedExpTimePiVsEta"))->Fill((Int_t)GetStripIndex(volId),tofTime-fTrkExpTimes[AliPID::kPion]);//ps
	((TH1F*)fHlistPID->FindObject("hTOFmatchedExpTimePi"))->Fill(tofTime-fTrkExpTimes[AliPID::kPion]);//ps
	((TH2F*)fHlistPID->FindObject("hTOFmatchedExpTimePiVsP"))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kPion]));
	((TH1F*)fHlistPID->FindObject("hTOFtheoreticalExpTimePi"))->Fill(tofTime-fThExpTimes[AliPID::kPion]);//ps
	((TH2F*)fHlistPID->FindObject("hTOFtheoreticalExpTimePiVsP"))->Fill(mom,(tofTime-fThExpTimes[AliPID::kPion]));	
	((TH2F*)fHlistPID->FindObject("hTOFExpSigmaPi"))->Fill(pT,(tofTime-fTrkExpTimes[AliPID::kPion])/fSigmaSpecie[AliPID::kPion]);
      }
      
      if (isValidBeta[AliPID::kKaon]){
	((TH1F*)fHlistPID->FindObject("hTOFmatchedExpTimeKa"))->Fill(tofTime-fTrkExpTimes[AliPID::kKaon]);//ps
	((TH2F*)fHlistPID->FindObject("hTOFmatchedExpTimeKaVsP"))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kKaon]));
	((TH1F*)fHlistPID->FindObject("hTOFtheoreticalExpTimeKa"))->Fill(tofTime-fThExpTimes[AliPID::kKaon]);//ps
	((TH2F*)fHlistPID->FindObject("hTOFtheoreticalExpTimeKaVsP"))->Fill(mom,(tofTime-fThExpTimes[AliPID::kKaon]));
	((TH2F*)fHlistPID->FindObject("hTOFExpSigmaKa"))->Fill(pT,(tofTime-fTrkExpTimes[AliPID::kKaon])/fSigmaSpecie[AliPID::kKaon]);
      }
      if (isValidBeta[AliPID::kProton]){
	((TH1F*)fHlistPID->FindObject("hTOFmatchedExpTimePro"))->Fill(tofTime-fTrkExpTimes[AliPID::kProton]);//ps
	((TH2F*)fHlistPID->FindObject("hTOFmatchedExpTimeProVsP"))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kProton]));
	((TH1F*)fHlistPID->FindObject("hTOFtheoreticalExpTimePro"))->Fill(tofTime-fThExpTimes[AliPID::kProton]);//ps
	((TH2F*)fHlistPID->FindObject("hTOFtheoreticalExpTimeProVsP"))->Fill(mom,(tofTime-fThExpTimes[AliPID::kProton]));
	((TH2F*)fHlistPID->FindObject("hTOFExpSigmaPro"))->Fill(pT,(tofTime-fTrkExpTimes[AliPID::kProton])/fSigmaSpecie[AliPID::kProton]);
      }

      if (fEnableAdvancedCheck && (pT<1.)) {
	Double_t pos[3]={0.,0.,0.};
	track->GetXYZAt(378.,5.,pos);
	if ((pos[0]==0.)&&(pos[1]==0.)&&(pos[2]==0.))continue;
	
	Double_t phiTOF=TMath::ATan2(pos[1],pos[0])*TMath::RadToDeg();
	if (phiTOF<0) phiTOF+= (2*TMath::Pi()*TMath::RadToDeg());
	
	if (isValidBeta[AliPID::kProton]){
	  if (track->GetSign()>0){
	    ((TH2F*)fHpos->FindObject("hTOFmatchedExpTimeProVsPhiPos"))->Fill(phiTOF,tofTime-fTrkExpTimes[AliPID::kProton]);//ps
	    ((TH2F*)fHpos->FindObject("hTOFmatchedExpTimeProVsPhiVertexPos"))->Fill(phi,tofTime-fTrkExpTimes[AliPID::kProton]);//ps
	  } else {
	    ((TH2F*)fHneg->FindObject("hTOFmatchedExpTimeProVsPhiNeg"))->Fill(phiTOF,tofTime-fTrkExpTimes[AliPID::kProton]);//ps
	    ((TH2F*)fHneg->FindObject("hTOFmatchedExpTimeProVsPhiVertexNeg"))->Fill(phi,tofTime-fTrkExpTimes[AliPID::kProton]);//ps
	  }
	}
	
	track->GetOuterXYZ(pos);
	Double_t phiOuterTPC=TMath::ATan2(pos[1],pos[0])*TMath::RadToDeg();
	if (phiOuterTPC<0) phiOuterTPC+= (2*TMath::Pi()*TMath::RadToDeg());
	
	if (track->GetSign()>0){
	  if (isValidBeta[AliPID::kProton])
	    ((TH2F*)fHpos->FindObject("hTOFmatchedExpTimeProVsPhiTPCPos"))->Fill(phiOuterTPC,tofTime-fTrkExpTimes[AliPID::kProton]);//ps
	  if (isValidBeta[AliPID::kPion])
	    ((TH2F*)fHpos->FindObject("hTOFmatchedExpTimePiVsPhiTPCPos"))->Fill(phiOuterTPC,tofTime-fTrkExpTimes[AliPID::kPion]);//ps
	  if (isValidBeta[AliPID::kKaon])
	    ((TH2F*)fHpos->FindObject("hTOFmatchedExpTimeKaVsPhiTPCPos"))->Fill(phiOuterTPC,tofTime-fTrkExpTimes[AliPID::kKaon]);//ps
	  
	  if ((phiOuterTPC<=30) || ((phiOuterTPC>=150)&&(phiOuterTPC<=230)) || (phiOuterTPC>=310) ) { //TRD sectors
	    if (isValidBeta[AliPID::kPion])
	      ((TH2F*)fHlistPID->FindObject("hTOFmatchedExpTimePiVsPTRDPos"))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kPion]));
	    if (isValidBeta[AliPID::kKaon])
	      ((TH2F*)fHlistPID->FindObject("hTOFmatchedExpTimeKaVsPTRDPos"))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kKaon]));
	    if (isValidBeta[AliPID::kProton])
	      ((TH2F*)fHlistPID->FindObject("hTOFmatchedExpTimeProVsPTRDPos"))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kProton]));
	  }
	  if ( ((phiOuterTPC>=50)&&(phiOuterTPC<=130)) || ((phiOuterTPC>=250)&&(phiOuterTPC<=290)) ) {//no TRD sectors
	    if (isValidBeta[AliPID::kPion])
	      ((TH2F*)fHlistPID->FindObject("hTOFmatchedExpTimePiVsPNoTRDPos"))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kPion]));
	    if (isValidBeta[AliPID::kKaon])
	      ((TH2F*)fHlistPID->FindObject("hTOFmatchedExpTimeKaVsPNoTRDPos"))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kKaon]));
	    if (isValidBeta[AliPID::kProton])
	      ((TH2F*)fHlistPID->FindObject("hTOFmatchedExpTimeProVsPNoTRDPos"))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kProton]));
	  }
	}else{
	  if (isValidBeta[AliPID::kPion])
	    ((TH2F*)fHneg->FindObject("hTOFmatchedExpTimePiVsPhiTPCNeg"))->Fill(phiOuterTPC,tofTime-fTrkExpTimes[AliPID::kPion]);//ps
	  if (isValidBeta[AliPID::kKaon])
	    ((TH2F*)fHneg->FindObject("hTOFmatchedExpTimeKaVsPhiTPCNeg"))->Fill(phiOuterTPC,tofTime-fTrkExpTimes[AliPID::kKaon]);//ps
	  if (isValidBeta[AliPID::kProton])
	    ((TH2F*)fHneg->FindObject("hTOFmatchedExpTimeProVsPhiTPCNeg"))->Fill(phiOuterTPC,tofTime-fTrkExpTimes[AliPID::kProton]);//ps
	  
	  if ((phiOuterTPC<=30) || ((phiOuterTPC>=150)&&(phiOuterTPC<=230)) || (phiOuterTPC>=310) ) { //TRD sectors
	    if (isValidBeta[AliPID::kPion])
	      ((TH2F*)fHlistPID->FindObject("hTOFmatchedExpTimePiVsPTRDNeg"))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kPion]));
	    if (isValidBeta[AliPID::kKaon])
	      ((TH2F*)fHlistPID->FindObject("hTOFmatchedExpTimeKaVsPTRDNeg"))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kKaon]));
	    if (isValidBeta[AliPID::kProton])
	      ((TH2F*)fHlistPID->FindObject("hTOFmatchedExpTimeProVsPTRDNeg"))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kProton]));
	  }
	  if ( ((phiOuterTPC>=50)&&(phiOuterTPC<=130)) || ((phiOuterTPC>=250)&&(phiOuterTPC<=290)) ) {//no TRD sectors
	    if (isValidBeta[AliPID::kPion])
	      ((TH2F*)fHlistPID->FindObject("hTOFmatchedExpTimePiVsPNoTRDNeg"))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kPion]));
	      if (isValidBeta[AliPID::kKaon])
		((TH2F*)fHlistPID->FindObject("hTOFmatchedExpTimeKaVsPNoTRDNeg"))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kKaon]));
	      if (isValidBeta[AliPID::kProton])
		((TH2F*)fHlistPID->FindObject("hTOFmatchedExpTimeProVsPNoTRDNeg"))->Fill(mom,(tofTime-fTrkExpTimes[AliPID::kProton]));
	  }
	}
	
      }//end advanced checks && cut on pT for t-texp phi distrib      
    }//matched
  }//end loop on tracks
  
  ((TH1F*)fHlist->FindObject("hTOFmatchedPerEvt"))->Fill(fNTOFtracks) ;
  ((TH2F*)fHlistTimeZero->FindObject("hT0TOFvsNtrk"))->Fill(fNTOFtracks,timeZero[AliESDpid::kTOF_T0]);
  
  PostData(1, fHlist);
  PostData(2, fHlistTimeZero);
  PostData(3, fHlistPID);
  PostData(4, fHpos);
  PostData(5, fHneg);
  
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
