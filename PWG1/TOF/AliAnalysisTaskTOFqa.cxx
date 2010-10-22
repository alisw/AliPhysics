/*  created by fbellini@cern.ch on 14/09/2010 */
/*  last modified by fbellini   on 21/10/2010 */


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
  fNTOFtracks(0), 
  fNPrimaryTracks(0), 
  fT0(0), 
  fHlist(0x0),
  fHlistExperts(0x0)
 {
  // Default constructor
}
//________________________________________________________________________
AliAnalysisTaskTOFqa::AliAnalysisTaskTOFqa(const char *name) : 
  AliAnalysisTaskSE(name), 
  fRunNumber(0), 
  fESD(0x0), 
  fTrackFilter(0x0),
  fVertex(0x0),
  fNTOFtracks(0), 
  fNPrimaryTracks(0), 
  fT0(0), 
  fHlist(0x0),
  fHlistExperts(0)
 {
  // Constructor
  // Define input and output slots here
   Info("AliAnalysisTaskTOFqa","Calling Constructor");
   
   // Input slot #0 works with a TChain
   DefineInput(0, TChain::Class());
   
   // Output slot #0 writes into a TH1 container
   // Output slot #1 writes into a user defined  container
   DefineOutput(1, TList::Class());
   DefineOutput(2, TList::Class());
 }

//________________________________________________________________________
AliAnalysisTaskTOFqa::AliAnalysisTaskTOFqa(const AliAnalysisTaskTOFqa& copy) 
: AliAnalysisTaskSE(), 
  fRunNumber(copy.fRunNumber), 
  fESD(copy.fESD), 
  fTrackFilter(copy.fTrackFilter), 
  fVertex(copy.fVertex),
  fNTOFtracks(copy.fNTOFtracks), 
  fNPrimaryTracks(copy.fNPrimaryTracks), 
  fT0(copy.fT0),
  fHlist(copy.fHlist),
  fHlistExperts(copy.fHlistExperts)
{
  // Copy constructor
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
    fNTOFtracks=copy.fNTOFtracks; 
    fNPrimaryTracks=copy.fNPrimaryTracks; 
    fT0=copy.fT0;
    fHlist=copy.fHlist;
    fHlistExperts=copy.fHlistExperts;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskTOFqa::~AliAnalysisTaskTOFqa() {
  //
  //destructor
  //

  Info("~AliAnalysisTaskTOFqa","Calling Destructor");
  if (fVertex) delete fVertex;
  if (fTrackFilter) delete fTrackFilter;
  if (fHlist) {
    delete fHlist;
    fHlist = 0;
  }
  if (fHlistExperts) {
    delete fHlistExperts;
    fHlistExperts = 0;
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
  if (!fHlistExperts) fHlistExperts = new TList();	
  fHlistExperts->SetOwner(kTRUE);
  //0
  TH1I* hTOFmatchedESDperEvt = new TH1I("hTOFmatchedPerEvt", "Number of matched TOF tracks per event;Number of TOF matched ESD tracks;Counts", 100, 0, 100) ;  
  hTOFmatchedESDperEvt->SetLineWidth(2);
  hTOFmatchedESDperEvt->SetLineColor(kBlue);
  fHlist->AddLast(hTOFmatchedESDperEvt) ;
  //1
  TH1F* hTOFmatchedESDtime = new TH1F("hTOFmatchedESDtime", "Matched  ESDs tracks: TOF Time spectrum; t [ns];Counts", 250, 0., 610. ) ; 
  hTOFmatchedESDtime->SetLineWidth(2);
  hTOFmatchedESDtime->SetLineColor(kBlue);
  hTOFmatchedESDtime->SetFillColor(kBlue);
  fHlist->AddLast(hTOFmatchedESDtime) ;
  //2
  TH1F* hTOFmatchedESDrawTime = new TH1F("hTOFmatchedESDrawTime", "Matched ESDs tracks: TOF raw Time spectrum;t_{raw} [ns];Counts", 250, 0., 610.) ; 
  hTOFmatchedESDrawTime->SetLineWidth(2);
  hTOFmatchedESDrawTime->SetLineColor(kGreen+2);
  hTOFmatchedESDrawTime->SetFillColor(kGreen+2);
  fHlist->AddLast(hTOFmatchedESDrawTime) ;
  //3
  TH1F* hTOFmatchedESDToT = new TH1F("hTOFmatchedESDToT", "Matched ESDs tracks: TOF ToT spectrum; ToT [ns];Counts",100, 0., 48.8) ; 
  hTOFmatchedESDToT->SetLineWidth(2);
  hTOFmatchedESDToT->SetLineColor(kBlue);
  hTOFmatchedESDToT->SetFillColor(kBlue);
  fHlist->AddLast(hTOFmatchedESDToT) ;
  //4
  TH1F* hTOFmatchedESDP  = new TH1F("hTOFmatchedESDP", "TPC-TOF matched tracks momentum distribution (GeV/c); P(GeV/c);tracks", 500,0.,5.) ;  
  hTOFmatchedESDP->SetLineWidth(1);
  hTOFmatchedESDP->SetLineColor(kBlue);
  hTOFmatchedESDP->SetMarkerStyle(21);
  hTOFmatchedESDP->SetMarkerSize(0.6);
  hTOFmatchedESDP->SetMarkerColor(kBlue);
  fHlist->AddLast(hTOFmatchedESDP) ; 
  //5
  TH1F* hTOFmatchedESDPt  = new TH1F("hTOFmatchedESDPt", "TPC-TOF matched tracks p_{T} distribution (GeV/c); p_{T}(GeV/c);tracks", 500,0.,5.) ;  
  hTOFmatchedESDPt->SetLineWidth(1);
  hTOFmatchedESDPt->SetLineColor(kBlue);
  hTOFmatchedESDPt->SetMarkerStyle(4);
  hTOFmatchedESDPt->SetMarkerSize(0.6);
  hTOFmatchedESDPt->SetMarkerColor(kBlue);
  fHlist->AddLast(hTOFmatchedESDPt) ; 
  //6
  TH1F* hTOFmatchedESDtrkLength  = new TH1F("hTOFmatchedESDtrkLength", "Matched ESDs tracks length; Track length [cm];Counts", 1600, -800., 800) ; 
  hTOFmatchedESDtrkLength->SetLineWidth(1);
  hTOFmatchedESDtrkLength->SetLineColor(kBlue);
  fHlist->AddLast(hTOFmatchedESDtrkLength) ; 
  //7
  TH2F* hTOFmatchedESDpVsBeta  = new TH2F("hTOFmatchedESDpVsBeta", "Matched ESDs tracks p vs. beta; p(GeV/c); beta", 500, 0., 5.,500, 0., 5.) ; 
  fHlist->AddLast(hTOFmatchedESDpVsBeta);
  
  //8
  TH1F* hESDprimaryTrackP = new TH1F("hESDprimaryTrackP", "All ESDs tracks Pt distribution (GeV/c); p_{T}(GeV/c);tracks", 500, 0., 5.0) ;  
  hESDprimaryTrackP->SetLineWidth(1);
  hESDprimaryTrackP->SetMarkerStyle(8);
  hESDprimaryTrackP->SetMarkerSize(0.6);
  hESDprimaryTrackP->SetLineColor(kGray+1);
  fHlist->AddLast(hESDprimaryTrackP);
  //9
  TH1F* hESDprimaryTrackPt = new TH1F("hESDprimaryTrackPt", "ESDs primary tracks Pt distribution (GeV/c); p_{T}(GeV/c);tracks", 500, 0.0, 5.0) ;  
  hESDprimaryTrackPt->SetLineWidth(1);
  hESDprimaryTrackPt->SetMarkerStyle(4);
  hESDprimaryTrackPt->SetMarkerSize(0.6);
  hESDprimaryTrackPt->SetLineColor(kBlack);
  fHlist->AddLast(hESDprimaryTrackPt);
  //10
  TH1F* hTOFmatchedESDeta = new TH1F("hTOFmatchedESDeta", "Matched ESDtracks eta; eta;Counts", 500, -2.5, 2.5) ; 
  fHlist->AddLast(hTOFmatchedESDeta) ; 
  //11
  TH1F* hTOFprimaryESDeta = new TH1F("hTOFprimaryESDeta", "Primary ESDtracks eta; eta;Counts", 500, -2.5, 2.5) ; 
  fHlist->AddLast(hTOFprimaryESDeta) ; 
  //12
  TH1F* hTOFmatchedESDphi = new TH1F("hTOFmatchedESDphi", "Matched ESDtracks Phi; Phi (rad);Counts", 65, 0., 6.5) ; 
  fHlist->AddLast(hTOFmatchedESDphi) ; 
  //13
  TH1F* hTOFprimaryESDphi = new TH1F("hTOFprimaryESDphi", "Primary ESDtracks Phi;Phi (rad);Counts", 65, 0., 6.5) ; 
  fHlist->AddLast(hTOFprimaryESDphi) ; 
  
  //Experts 0
  TH1F* hTOFESDsMatchingProb  = new TH1F("hTOFESDsMatchingProb", "TPC-TOF track-matching probability per event(|eta|<0.9 && pt>0.5GeV/c);TPC-TOF track-matching probability (%)  ;Counts",21, 0., 110.) ;  
  hTOFESDsMatchingProb->SetLineColor(kRed);
  fHlistExperts->AddLast(hTOFESDsMatchingProb) ; 
  
  //Experts 1
  TH1F* hTOFmatchedESDdiffTime  = new TH1F("hTOFmatchedESDdiffTime", "ESDs t_{TOF}-t_{pi,exp} spectrum in TOF (ps); t_{TOF}-t_{pi,exp} [ps];Counts", 4000, -50000., 50000.) ; 
  hTOFmatchedESDdiffTime->SetLineWidth(1);
  hTOFmatchedESDdiffTime->SetLineColor(kBlack);
  hTOFmatchedESDdiffTime->SetMarkerStyle(8);
  hTOFmatchedESDdiffTime->SetMarkerSize(0.8);
  hTOFmatchedESDdiffTime->SetMarkerColor(kAzure+7);
  hTOFmatchedESDdiffTime->SetFillColor(kAzure-2);
  fHlistExperts->AddLast(hTOFmatchedESDdiffTime) ; 
  
  //Experts 2
  TH2F* hTOFmatchedESDdiffTimeVsStrip = new TH2F("hTOFmatchedESDdiffTimeVsStrip", "ESDs t_{TOF}-t_{pi,exp} vs strip number; strip (Eta); t_{TOF}-t_{pi,exp} [ps]", 92,0.,92,400, -5000., 5000.) ; 
  fHlistExperts->AddLast(hTOFmatchedESDdiffTimeVsStrip) ; 
  
  //Experts 3
  TH2F* hTOFmatchedESDdxVsEta = new TH2F("hTOFmatchedESDdxVsEta", "Matched ESD tracks Dx vs eta; strip(eta); Dx (cm)", 92,0.,92., 200,-10.,10.) ; 
  fHlistExperts->AddLast(hTOFmatchedESDdxVsEta) ; 
  
  //Experts 4
  TH2F* hTOFmatchedESDdzVsEta  = new TH2F("hTOFmatchedESDdzVsEta", "Matched ESDtracks Dz vs eta; strip(eta); Dz (cm)", 92,0.,92., 200,-10.,10.) ; 
  fHlistExperts->AddLast(hTOFmatchedESDdzVsEta) ; 
  
  //Experts 5
  TH1F* hTOFmatchedMass= new TH1F("hTOFmatchedMass","Matched ESD tracks mass distribution; M (GeV/c^{2}); entries", 500,0., 5. );
  hTOFmatchedMass->SetLineWidth(2);
  hTOFmatchedMass->SetLineColor(kBlue);
  hTOFmatchedMass->SetLineColor(kBlue);
  fHlistExperts->AddLast(hTOFmatchedMass);
   
  //Experts 6
  TH1D* hEventT0DetAND = new TH1D("hEventT0DetAND", "Event T0 from T0 detector (A&C); t [ps];Counts", 4000, -50000., 50000. ) ; 
  hEventT0DetAND->SetLineWidth(2);
  hEventT0DetAND->SetLineColor(kRed);
  hEventT0DetAND->SetFillColor(kRed);
  fHlistExperts->AddLast(hEventT0DetAND) ;

  //Experts 7
  TH1D* hEventT0DetA = new TH1D("hEventT0DetA", "Event T0 from T0 detector (A side); t [ps];Counts", 4000, -50000., 50000. ) ; 
  hEventT0DetA->SetLineWidth(2);
  hEventT0DetA->SetLineColor(kBlue);
  hEventT0DetA->SetFillColor(kBlue);
  fHlistExperts->AddLast(hEventT0DetA) ;

   //Experts 8
  TH1D* hEventT0DetC = new TH1D("hEventT0DetC", "Event T0 from T0 detector (C side); t [ps];Counts", 4000, -50000., 50000. ) ; 
  hEventT0DetC->SetLineWidth(2);
  hEventT0DetC->SetLineColor(kGreen);
  hEventT0DetC->SetFillColor(kGreen);
  fHlistExperts->AddLast(hEventT0DetC);
 
  //Experts 9
  TH1F* hTOFmatchedExpTime = new TH1F("hTOFmatchedExpTime", "Matched  ESDs tracks - pions expected  time; t [ns];Counts", 4000, -50000., 50000. ) ; 
  hTOFmatchedExpTime->SetLineWidth(1);
  hTOFmatchedExpTime->SetLineColor(kBlack);
  hTOFmatchedExpTime->SetMarkerStyle(8);
  hTOFmatchedExpTime->SetMarkerSize(0.8); 
  hTOFmatchedExpTime->SetMarkerColor(kRed);
  fHlistExperts->AddLast(hTOFmatchedExpTime) ;

  
  PostData(1, fHlist);
  PostData(2, fHlistExperts);
}
//________________________________________________________________________
void AliAnalysisTaskTOFqa::UserExec(Option_t *) 
{ 
  /* Main - executed for each event.
    It extracts event information and track information after selecting 
    primary tracks via standard cuts. */
  const Double_t speedOfLight =  TMath::C()*1E2*1E-12; // cm/ps
  
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
  
  // loop over ESD tracks
  fNTOFtracks=0;
  fNPrimaryTracks=0;
  
  //info from To detector
  fT0= fESD->GetT0TOF(0);//ps
  ((TH1D*)fHlistExperts->At(6))->Fill(fT0);//ps
  ((TH1D*)fHlistExperts->At(7))->Fill((Double_t)fESD->GetT0TOF(1)); //ps
  ((TH1D*)fHlistExperts->At(8))->Fill((Double_t)fESD->GetT0TOF(2));//ps
      
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }

    //primary tracks selection: kTPCrefit and std cuts
    if(!fTrackFilter->IsSelected(track)) continue;
   
    Double_t tofTime=track->GetTOFsignal();//in ps
    Double_t tofTimeRaw=track->GetTOFsignalRaw();//in ps
    Double_t tofToT=track->GetTOFsignalToT(); //in ps
    Double_t expTimes[5];
    track->GetIntegratedTimes(expTimes);
    Double_t length =track->GetIntegratedLength();
    Double_t ptot[3];
    track->GetConstrainedPxPyPz(ptot);    
    Double_t pT = TMath::Sqrt(ptot[0]*ptot[0] + ptot[1]*ptot[1]);
    Double_t P2 = pT*pT + ptot[2]*ptot[2];
    Double_t eta=track->Eta();
    Int_t channel=track->GetTOFCalChannel(); 
    Int_t volId[5]; //(sector, plate,strip,padZ,padX)
    AliTOFGeometry::GetVolumeIndices(channel,volId);
    
    if (TMath::Abs(eta)<0.9) { //cut for acceptance
      if (P2>=0)
	((TH1F*)fHlist->At(8))->Fill(TMath::Sqrt(P2)); 
      ((TH1F*)fHlist->At(9))->Fill(track->Pt()); //all esd tracks within acceptance
      ((TH1F*)fHlist->At(11))->Fill(eta);
      ((TH1F*)fHlist->At(13))->Fill(track->Phi());
      if ((TMath::Abs(track->Eta())<0.9)&&(track->Pt()>0.5)) fNPrimaryTracks++;
      
      //matched tracks selection: kTOFout and kTIME
      if ((track->GetStatus() & AliESDtrack::kTOFout) &&
	  (track->GetStatus() & AliESDtrack::kTIME)) {      
	
	if ((TMath::Abs(track->Eta())<0.9)&&(track->Pt()>0.5)) fNTOFtracks++; //matched counter
	
	  ((TH1F*)fHlist->At(1))->Fill(tofTime*1E-3); //ns
	  ((TH1F*)fHlist->At(2))->Fill(tofTimeRaw*1E-3); //ns
	  ((TH1F*)fHlist->At(3))->Fill(tofToT);
	  if (P2>=0)
	    ((TH1F*)fHlist->At(4))->Fill(TMath::Sqrt(P2));
	  ((TH1F*)fHlist->At(5))->Fill(track->Pt());
	  ((TH1F*)fHlist->At(6))->Fill(length);  
	  ((TH1F*)fHlist->At(10))->Fill(eta);
	  ((TH1F*)fHlist->At(12))->Fill(track->Phi());
	  
	  ((TH1F*)fHlistExperts->At(1))->Fill(tofTime-expTimes[2]);//ps
	  ((TH1F*)fHlistExperts->At(2))->Fill((Int_t)GetStripIndex(volId),tofTime-expTimes[2]); //ps
	  ((TH1F*)fHlistExperts->At(3))->Fill((Int_t)GetStripIndex(volId),track->GetTOFsignalDx());
	  ((TH1F*)fHlistExperts->At(4))->Fill((Int_t)GetStripIndex(volId),track->GetTOFsignalDz());
	  ((TH1F*)fHlistExperts->At(9))->Fill(expTimes[2]);//ps

	  //basic PID performance check 
	  Double_t tof= tofTime; //average T0 fill subtracted, no info from T0detector 
	  if (length>350){
	    Double_t beta= length/(tof*speedOfLight);
	    //Double_t mass2=P2*((tof/length)*(tof/length)-(1/(speedOfLight*speedOfLight)));
	    ((TH1F*)fHlist->At(7))->Fill(TMath::Sqrt(P2),beta);
	    //if (mass2>=0)((TH1F*)fHlistExperts->At(5))->Fill(TMath::Sqrt(mass2));
	  }
      }//matched
    }//acceptance cut
  }//end loop on tracks
  ((TH1F*)fHlist->At(0))->Fill(fNTOFtracks) ;
  //if (fNTOFtracks>fNPrimaryTracks) printf("Something strange!!!\n");
  if(fNPrimaryTracks>0){
    ((TH1F*)fHlistExperts->At(0))->Fill((fNTOFtracks/(Float_t)fNPrimaryTracks)*100) ;
  }
  
 
  PostData(1, fHlist);
  PostData(2, fHlistExperts);
  
}      

//________________________________________________________________________
void AliAnalysisTaskTOFqa::Terminate(Option_t *) 
{
  //check on output validity
  fHlist = dynamic_cast<TList*> (GetOutputData(1));
  if (!fHlist || !fHlistExperts) {
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
