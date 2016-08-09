/*  created by fbellini@cern.ch on 29/04/2013 */
/*  last modified by fbellini   on 19/08/2013 */


#ifndef ALIANALYSISTASKTOFQAID_CXX
#define ALIANALYSISTASKTOFQAID_CXX

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliESDpid.h"
#include "AliTOFPIDParams.h"
#include "AliTOFChannelOnlineStatusArray.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBPath.h"
#include "AliTOFcalib.h"
#include "AliTOFT0maker.h"
#include "AliTOFT0v1.h"
#include "AliAnalysisTaskTOFqaID.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliTOFRawStream.h"
#include "AliTOFGeometry.h"


ClassImp(AliAnalysisTaskTOFqaID)

//________________________________________________________________________
AliAnalysisTaskTOFqaID::AliAnalysisTaskTOFqaID() :
  fRunNumber(0), 
  fESD(0x0), 
  fMCevent(0x0),
  fTrackFilter(0x0), 
  fVertex(0x0),
  fESDpid(new AliESDpid()),
  fTOFHeader(0x0),
  fEnableAdvancedCheck(kFALSE),
  fEnableChargeSplit(kFALSE),
  fExpTimeBinWidth(24.4),
  fExpTimeRangeMin(-25010.),
  fExpTimeRangeMax(25010.),
  fExpTimeSmallRangeMin(-5002.),
  fExpTimeSmallRangeMax(5002.),
  fnExpTimeBins(1),
  fnExpTimeSmallBins(1),  
  fMyTimeZeroTOF(1e20),
  fMyTimeZeroTOFsigma(1e20),
  fMyTimeZeroTOFtracks(-1),
  fIsMC(kFALSE),
  fSelectedPdg(0),
  fP(1e10),
  fPt(1e10),
  fEta(1e10),
  fPhi(1e10),
  fTPCOuterPhi(1e10),
  fL(1e10),
  fMatchingMomCut(0.0),
  fMatchingEtaCut(1e10),
  fTof(1e10),
  fOCDBLocation("local://$ALICE_PHYSICS/OCDB"),
  fChannelArray(0x0),
  fCalib(0x0),
  fHlist(0x0),
  fHlistTimeZero(0x0),
  fHlistPID(0x0),
  fHlistTRD(0x0),
  fHlistTrigger(0x0)
{
  // Default constructor
   
   for (Int_t j=0;j<3;j++ ) {
     if (j<3) {
       fT0[j]=0.0;
       fNTOFtracks[j]=0;
     }
     fSigmaSpecie[j]=0.0;
     fTrkExpTimes[j]=0.0;
     fThExpTimes[j]=0.0;
   }
}
//________________________________________________________________________
AliAnalysisTaskTOFqaID::AliAnalysisTaskTOFqaID(const char *name) : 
  AliAnalysisTaskSE(name), 
  fRunNumber(0), 
  fESD(0x0), 
  fMCevent(0x0),
  fTrackFilter(0x0),
  fVertex(0x0),
  fESDpid(new AliESDpid()),
  fTOFHeader(0x0),
  fEnableAdvancedCheck(kFALSE),
  fEnableChargeSplit(kFALSE),
  fExpTimeBinWidth(24.4),
  fExpTimeRangeMin(-25010.),
  fExpTimeRangeMax(25010.),
  fExpTimeSmallRangeMin(-5002.),
  fExpTimeSmallRangeMax(5002.),
  fnExpTimeBins(1),
  fnExpTimeSmallBins(1),
  fMyTimeZeroTOF(1e20),
  fMyTimeZeroTOFsigma(1e20),
  fMyTimeZeroTOFtracks(-1),
  fIsMC(kFALSE),
  fSelectedPdg(0),
  fP(1e10),
  fPt(1e10),
  fEta(1e10),
  fPhi(1e10),
  fTPCOuterPhi(1e10),
  fL(1e10),
  fMatchingMomCut(1.0),
  fMatchingEtaCut(0.8),
  fTof(1e10),
  fOCDBLocation("local://$ALICE_PHYSICS/OCDB"),
  fChannelArray(0x0),
  fCalib(0x0),
  fHlist(0x0),
  fHlistTimeZero(0x0),
  fHlistPID(0x0),
  fHlistTRD(0x0),
  fHlistTrigger(0x0)
 {
  // Constructor
  // Define input and output slots here
   Info("AliAnalysisTaskTOFqaID","Calling Constructor");
   
   for (Int_t j=0;j<5;j++ ) {
     if (j<3){ 
       fT0[j]=0.0;
       fNTOFtracks[j]=0;
     }
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
AliAnalysisTaskTOFqaID::AliAnalysisTaskTOFqaID(const AliAnalysisTaskTOFqaID& copy) 
: AliAnalysisTaskSE(), 
  fRunNumber(copy.fRunNumber), 
  fESD(copy.fESD), 
  fMCevent(copy.fMCevent),
  fTrackFilter(copy.fTrackFilter), 
  fVertex(copy.fVertex),
  fESDpid(copy.fESDpid),
  fTOFHeader(copy.fTOFHeader),
  fEnableAdvancedCheck(copy.fEnableAdvancedCheck),
  fEnableChargeSplit(copy.fEnableChargeSplit),
  fExpTimeBinWidth(copy.fExpTimeBinWidth),
  fExpTimeRangeMin(copy.fExpTimeRangeMin),
  fExpTimeRangeMax(copy.fExpTimeRangeMax),
  fExpTimeSmallRangeMin(copy.fExpTimeSmallRangeMin),
  fExpTimeSmallRangeMax(copy.fExpTimeSmallRangeMax),
  fnExpTimeBins(copy.fnExpTimeBins),
  fnExpTimeSmallBins(copy.fnExpTimeSmallBins),
  fMyTimeZeroTOF(copy.fMyTimeZeroTOF),
  fMyTimeZeroTOFsigma(copy.fMyTimeZeroTOFsigma),
  fMyTimeZeroTOFtracks(copy.fMyTimeZeroTOFtracks),
  fIsMC(copy.fIsMC),
  fSelectedPdg(copy.fSelectedPdg),
  fP(copy.fP),
  fPt(copy.fPt),
  fEta(copy.fEta),
  fPhi(copy.fPhi),
  fTPCOuterPhi(copy.fTPCOuterPhi),
  fL(copy.fL),
  fMatchingMomCut(copy.fMatchingMomCut),
  fMatchingEtaCut(copy.fMatchingEtaCut),
  fTof(copy.fTof),
  fOCDBLocation(copy.fOCDBLocation),
  fChannelArray(copy.fChannelArray),
  fCalib(copy.fCalib),
  fHlist(copy.fHlist),
  fHlistTimeZero(copy.fHlistTimeZero),
  fHlistPID(copy.fHlistPID),
  fHlistTRD(copy.fHlistTRD),
  fHlistTrigger(copy.fHlistTrigger)
{
  // Copy constructor
   for (Int_t j=0;j<5;j++ ) {
     if (j<3) { 
       fT0[j]=copy.fT0[j];
       fNTOFtracks[j]=copy.fNTOFtracks[j];
     }
     fSigmaSpecie[j]=copy.fSigmaSpecie[j];
     fTrkExpTimes[j]=copy.fTrkExpTimes[j];
     fThExpTimes[j]=copy.fThExpTimes[j];
   }
  

}

//___________________________________________________________________________
AliAnalysisTaskTOFqaID& AliAnalysisTaskTOFqaID::operator=(const AliAnalysisTaskTOFqaID& copy) 
{
  //
  // Assignment operator
  //
  if (this!=&copy) {
    AliAnalysisTaskSE::operator=(copy) ;
    fRunNumber=copy.fRunNumber; 
    fESD=copy.fESD;
    fMCevent=copy.fMCevent;
    fTrackFilter=copy.fTrackFilter;
    fVertex=copy.fVertex;
    fESDpid=copy.fESDpid;
    fTOFHeader=copy.fTOFHeader;
    fEnableAdvancedCheck=copy.fEnableAdvancedCheck;
    fEnableChargeSplit=copy.fEnableChargeSplit;
    fExpTimeBinWidth=copy.fExpTimeBinWidth;
    fExpTimeRangeMin=copy.fExpTimeRangeMin;
    fExpTimeRangeMax=copy.fExpTimeRangeMax;
    fExpTimeSmallRangeMin=copy.fExpTimeSmallRangeMin;
    fExpTimeSmallRangeMax=copy.fExpTimeSmallRangeMax;
    fnExpTimeBins=copy.fnExpTimeBins;
    fnExpTimeSmallBins=copy.fnExpTimeSmallBins;
    fMyTimeZeroTOF=copy.fMyTimeZeroTOF;
    fMyTimeZeroTOFsigma=copy.fMyTimeZeroTOFsigma;
    fMyTimeZeroTOFtracks=copy.fMyTimeZeroTOFtracks;
    for (Int_t j=0;j<5;j++ ) {
      if (j<3) {
	fT0[j]=copy.fT0[j];
        fNTOFtracks[j]=copy.fNTOFtracks[j];
      }
      fSigmaSpecie[j]=copy.fSigmaSpecie[j];
      fTrkExpTimes[j]=copy.fTrkExpTimes[j];
      fThExpTimes[j]=copy.fThExpTimes[j];
    }
    fIsMC=copy.fIsMC;
    fSelectedPdg=copy.fSelectedPdg;
    fP=copy.fP;
    fPt=copy.fPt;
    fEta=copy.fEta;
    fPhi=copy.fPhi;
    fTPCOuterPhi=copy.fTPCOuterPhi;
    fL=copy.fL;
    fMatchingMomCut=copy.fMatchingMomCut;
    fMatchingEtaCut=copy.fMatchingEtaCut;
    fTof=copy.fTof;
    fOCDBLocation=copy.fOCDBLocation;
    fChannelArray=copy.fChannelArray;
    fCalib=copy.fCalib;
    fHlist=copy.fHlist;
    fHlistTimeZero=copy.fHlistTimeZero;
    fHlistPID=copy.fHlistPID;
    fHlistTRD=copy.fHlistTRD;
    fHlistTrigger=copy.fHlistTrigger;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskTOFqaID::~AliAnalysisTaskTOFqaID() {
  //
  //destructor
  //

  Info("~AliAnalysisTaskTOFqaID","Calling Destructor");
  if (fESDpid) delete fESDpid;
  if (fTOFHeader) delete fTOFHeader;
  if (fVertex) delete fVertex;
  if (fTrackFilter) delete fTrackFilter;
  if (fChannelArray) delete fChannelArray;
  if (fCalib) delete fCalib;
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
  if (fHlistTRD){
    delete fHlistTRD;
    fHlistTRD = 0;
  }
  if (fHlistTrigger){
    delete fHlistTrigger;
    fHlistTrigger = 0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskTOFqaID::UserCreateOutputObjects()
{
  //
  //Define output objects and histograms
  //
 
  //retrieve PID response object 
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if (!man)  AliFatal("Analysis manager needed");
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if (!inputHandler) AliFatal("Input handler needed");
  //pid response object
  fESDpid=(AliESDpid*)inputHandler->GetPIDResponse();
  if (!fESDpid) AliError("PIDResponse object was not created");
  //fESDpid->SetOADBPath("$ALICE_PHYSICS/OADB");

  Info("CreateOutputObjects","CreateOutputObjects (TList) of task %s", GetName());
  OpenFile(1);

  if (!fHlist) fHlist = new TList();	
  fHlist->SetOwner(kTRUE);
  fHlist->SetName("base");

  if (!fHlistTimeZero) fHlistTimeZero = new TList();	
  fHlistTimeZero->SetOwner(kTRUE);
  fHlistTimeZero->SetName("startTime");

  if (!fHlistPID) fHlistPID = new TList();	
  fHlistPID->SetOwner(kTRUE);
  fHlistPID->SetName("pid");

  if (!fHlistTRD) fHlistTRD = new TList();	
  fHlistTRD->SetOwner(kTRUE);  
  fHlistTRD->SetName("TRD");

  if (!fHlistTrigger) fHlistTrigger = new TList();	
  fHlistTrigger->SetOwner(kTRUE);  
  fHlistTrigger->SetName("trigger");

  if (fExpTimeRangeMax<fExpTimeRangeMin) {
    SetExpTimeHistoRange(-25010.,25010.);
  }
  fnExpTimeBins = TMath::Nint((fExpTimeRangeMax - fExpTimeRangeMin)/fExpTimeBinWidth);//ps
  fExpTimeRangeMax=fExpTimeRangeMin+fnExpTimeBins*fExpTimeBinWidth;//ps
  
  if (fExpTimeSmallRangeMax<fExpTimeSmallRangeMin) {
    SetExpTimeHistoSmallRange(-5002.,5002.);
  }
  fnExpTimeSmallBins = TMath::Nint((fExpTimeSmallRangeMax - fExpTimeSmallRangeMin)/fExpTimeBinWidth);//ps
  fExpTimeSmallRangeMax=fExpTimeSmallRangeMin+fnExpTimeSmallBins*fExpTimeBinWidth;//ps
  
  //add plots for start time QA
  AddStartTimeHisto(fHlistTimeZero,"");
  
  //add plots for base TOF quantities
  if (fEnableChargeSplit) {
    AddTofBaseHisto(fHlist,  1, "");
    AddTofBaseHisto(fHlist, -1, "");
  } else {
    AddTofBaseHisto(fHlist,  0, "");
  }
  //add plots for matching efficiency
   if (fEnableChargeSplit) {
     AddMatchingEffHisto(fHlist,  1, "");
     AddMatchingEffHisto(fHlist, -1, "");
   } else {
     AddMatchingEffHisto(fHlist,  0, "");
   }
  //add plots for PID checks
   if (fEnableChargeSplit) {
     AddPidHisto(fHlistPID,  1, ""); 
     AddPidHisto(fHlistPID, -1, ""); 
   } else {
     AddPidHisto(fHlistPID,  0, ""); 
   }
  //add trd plots
  if (fEnableAdvancedCheck) {
    AddTrdHisto();
  }
  //Add trigger plots
  AddTofTrgHisto("");
  
  PostData(1, fHlist);
  PostData(2, fHlistTimeZero);
  PostData(3, fHlistPID);
  PostData(4, fHlistTRD);
  PostData(5, fHlistTrigger);
}
//________________________________________________________________________
void AliAnalysisTaskTOFqaID::UserExec(Option_t *) 
{ 
  /* Main - executed for each event.
    It extracts event information and track information after selecting 
    primary tracks via standard cuts. */
  /*
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
    AliError("ERROR: Could not get ESDInputHandler");
    return;
  } else {
    fESD = (AliESDEvent*) esdH->GetEvent();
  } 
  
  */
  fESD=(AliESDEvent*)InputEvent();
  if (!fESD) {
    AliError("fESD event not available");
    return;
  }
  
  if (!fESDpid) {
    AliError("PID object fESDpid not available");
    return;
  }
  
  //retrieve default start time type from PIDresponse
  AliPIDResponse::EStartTimeType_t startTimeMethodDefault = AliPIDResponse::kBest_T0;  
  if (fESDpid->GetTOFPIDParams()) {  // during reconstruction OADB not yet available
    startTimeMethodDefault = ((AliTOFPIDParams *)fESDpid->GetTOFPIDParams())->GetStartTimeMethod();
  }
  
  //access MC event handler for MC truth information
  if (fIsMC) {
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!mcH) {
      AliError("Cannot get MCeventHandler");
      return;
    } else {
      fMCevent = (AliMCEvent *) mcH->MCEvent();
      if (!fMCevent) {
	AliError("Trying to retrieve an invalid MC event.");
	return;
      }
    }
  }
  
  // get run number
  Int_t runNb = fESD->GetRunNumber();
  if (runNb != fRunNumber) {
    //retrieve maps of channels from OCDB
    fRunNumber = runNb;
    LoadChannelMapsFromOCDB();
  }
  
  //reset matched track counters
  for (Int_t j=0;j<3;j++){fNTOFtracks[j]=0;}  

  //Get vertex info and apply vertex cut
  if (!IsEventSelected(fESD)) return;
  
  //set response tof_t0 for all other checks
  fESDpid->SetTOFResponse(fESD,AliESDpid::kTOF_T0);//(fill_t0, tof_t0, t0_t0, best_t0)
  
  AliDebug(3, Form("Momentum cut for eta and phi distributions set: Pt>%3.2f", fMatchingMomCut));

  //check existence of track filter
  if (!fTrackFilter){
    AliInfo("No track filter found, skipping the track loop");
    return;
  }
  
  // loop over ESD tracks
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      AliInfo(Form("Cannot receive track %d", iTracks));
      continue;
    }
    
    //primary tracks selection: kTPCrefit and std cuts
    if (!fTrackFilter->IsSelected(track)) continue;
    
    //select specie if MC
    if ( fIsMC && 
	 (!SelectMCspecies(fMCevent, track))) {
      AliDebug(4, Form("MC tracks selection: Track=%i  label=%i  Not Accepted", iTracks, track->GetLabel()));
      continue;
    }
    
    //apply cut for eta acceptance
    fEta=track->Eta();
    if (TMath::Abs(fEta)>fMatchingEtaCut) continue; 
    
    //get other track variables
    fP = track->P();
    fPt = track->Pt();
    fPhi = track->Phi()*TMath::RadToDeg();
    fTPCOuterPhi = GetPhiAtTPCouterRadius(track);
    fL = track->GetIntegratedLength();
    track->GetIntegratedTimes(fTrkExpTimes);
    
    Int_t charge = 0;
    if (fEnableChargeSplit) charge = track->Charge();
    
    //Fill histograms for primary particles
    FillPrimaryTrkHisto(charge,"");
    
    if (IsTPCTOFMatched(track)) {     
      fTof=track->GetTOFsignal()*1E-3;//in ps
      //increment track counters
      fNTOFtracks[0]++;
      if (charge>0) fNTOFtracks[1]++;
      if (charge<0) fNTOFtracks[2]++;
      //fill histos
      FillTofBaseHisto(track, charge,"");
      FillMatchedTrkHisto(charge,"");
      FillPidHisto(track, charge, "");
    }    
    if (fEnableAdvancedCheck) FillTrdHisto(track, charge);
  }//end loop on tracks
  
  //fill time zero histos  
  FillStartTimeHisto("");  
  if (fEnableChargeSplit) {
    ((TH1F*)fHlist->FindObject("hTOFmulti_pos"))->Fill(fNTOFtracks[1]);
    ((TH1F*)fHlist->FindObject("hTOFmulti_neg"))->Fill(fNTOFtracks[2]);
  } else {
    ((TH1F*)fHlist->FindObject("hTOFmulti_all"))->Fill(fNTOFtracks[0]);
  }
  //fill TOF trg histos from infos in TOF header
  fTOFHeader=(AliTOFHeader*)fESD->GetTOFHeader();
  if (!fTOFHeader) {
    AliWarning("Cannot get TOF header: no TOF trigger info available");
  } else {
    FillTofTrgHisto("");
  }
  
  //restore value set by AliPIDResponseTask for subsequent wagons
  fESDpid->SetTOFResponse(fESD,startTimeMethodDefault);
  
  PostData(1, fHlist);
  PostData(2, fHlistTimeZero);
  PostData(3, fHlistPID);
  PostData(4, fHlistTRD);
  PostData(5, fHlistTrigger);

}      

//________________________________________________________________________
void AliAnalysisTaskTOFqaID::Terminate(Option_t *) 
{
  //check on output validity
  fHlist = dynamic_cast<TList*> (GetOutputData(1));
  if (!fHlist) {
    AliError("Base histograms list not available");
    return;   
  } 
  
  // TH1F*hDummy = ((TH1F*)fHlist->FindObject("hTOFmatchedESDPt"));
  // TH1F*hMatchingEff = (TH1F*) hDummy->Clone("hMatchingEff");
  // hMatchingEff->SetTitle("Matching efficiency");
  // hMatchingEff->Divide((TH1F*) fHlist->FindObject("hESDprimaryTrackPt"));
  // TCanvas *c1 = new TCanvas("AliAnalysisTaskTOFqaID","Matching vs Pt",10,10,510,510);
  // c1->cd(1)->SetLogy();
  // hMatchingEff->DrawCopy("E");
  // fHlist->AddLast(hMatchingEff);
  ComputeMatchingEfficiency(fHlist, "pt");
  ComputeMatchingEfficiency(fHlist, "eta");
  ComputeMatchingEfficiency(fHlist, "phi");

  PostData(1, fHlist);
}

//---------------------------------------------------------------
Int_t AliAnalysisTaskTOFqaID::GetStripIndex(const Int_t * in)
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

//-----------------------------------------------------------------
Double_t AliAnalysisTaskTOFqaID::GetPhiAtTPCouterRadius(AliESDtrack * track)
{
  //get track phi at TPC outer radius
  if (!track) return 1e10;
  Double_t tpcoutcoord[3]={0.,0.,0.};
  track->GetOuterXYZ(tpcoutcoord);
  Double_t phiOuterTPC=TMath::ATan2(tpcoutcoord[1],tpcoutcoord[0])*TMath::RadToDeg();
  if (phiOuterTPC<0) 
    phiOuterTPC+= (2*TMath::Pi()*TMath::RadToDeg());
  return phiOuterTPC;
}
//-----------------------------------------------------------------
Bool_t  AliAnalysisTaskTOFqaID::IsEventSelected(AliESDEvent * event)
{
  //select event based on primary vertex
  if (!event) {
    AliError("Invalid ESD event");
    return kFALSE;
  }
  fVertex = (AliESDVertex*) event->GetPrimaryVertexTracks(); 
  if(fVertex->GetNContributors()<1) { 
    // SPD vertex
    fVertex = (AliESDVertex*) event->GetPrimaryVertexSPD(); 
    if(fVertex->GetNContributors()<1) fVertex = 0x0;
  }
  if (!fVertex) return kFALSE; 
  if (TMath::Abs(fVertex->GetZ())<10.0) return kTRUE;
  else return kFALSE;
}

//-----------------------------------------------------------------
Bool_t  AliAnalysisTaskTOFqaID::IsTPCTOFMatched(AliESDtrack * track)
{
  //defines TOF matching
  if (!track){
    AliWarning("Invalid track object");
    return kFALSE;
  }
  
  if ( (track->IsOn(AliESDtrack::kTOFout)) &&
       (track->IsOn(AliESDtrack::kTIME)) &&
       (track->IsOn(AliESDtrack::kTPCout))  )
    return kTRUE;
  else 
    return kFALSE;
}
//-----------------------------------------------------------------
Bool_t  AliAnalysisTaskTOFqaID::IsInTRD(AliESDtrack * track)
{
  //defines cut to select particles in/out TRD
  if (!track){
    AliWarning("Invalid track object");
    return kFALSE;
  }
  
  if ( track->IsOn(AliESDtrack::kTPCout) 
       && track->IsOn(AliESDtrack::kTRDout) )
    return kTRUE;
  else 
    return kFALSE; 
}
//-----------------------------------------------------------------
void AliAnalysisTaskTOFqaID::FillStartTimeMaskHisto(TString suffix)
{
  /* set pid response to use best_T0 and for each
     accepted track fills the histogram with the 
     used start time 
  */

  //set response best_t0 
  //fESDpid->SetTOFResponse(fESD,AliESDpid::kBest_T0);

  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      AliInfo(Form("Cannot receive track %d", iTracks));
      continue;
    }    
    //primary tracks selection: kTPCrefit and std cuts
    if (fTrackFilter){
      if(!fTrackFilter->IsSelected(track)) continue;
    }
    else{
      AliInfo("No track filter found, skipping the track loop");
      break;
    }
    if (TMath::Abs(track->Eta())>fMatchingEtaCut) continue; //cut for acceptance  
    
    Int_t StartTimeBit = fESDpid->GetTOFResponse().GetStartTimeMask(track->P());
    ((TH2F*)fHlistTimeZero->FindObject(Form("hStartTimeMask%s",suffix.Data())))->Fill(track->P(),StartTimeBit);
    
    //matched tracks selection: kTOFout and kTIME
    if ( (track->IsOn(AliESDtrack::kTOFout)) &&
	 (track->IsOn(AliESDtrack::kTIME)) &&
	 (track->IsOn(AliESDtrack::kTPCout))  ) {
      ((TH2F*)fHlistTimeZero->FindObject(Form("hStartTimeMaskMatched%s",suffix.Data())))->Fill(track->P(),StartTimeBit);
    }
  }
  return;
}

//----------------------------------------------------
Bool_t AliAnalysisTaskTOFqaID::ComputeTimeZeroByTOF1GeV()
{
  /* compute T0-TOF for tracks within momentum range [0.95, 1.05] */
  /* init T0-TOF */
  AliTOFT0v1 *fTOFT0v1 = new AliTOFT0v1(fESDpid); // TOF-T0 v1
  fTOFT0v1->Init(fESD);
  fTOFT0v1->DefineT0("all", 0.95, 1.05);
  fMyTimeZeroTOF = -1000. * fTOFT0v1->GetResult(0);
  fMyTimeZeroTOFsigma = 1000. * fTOFT0v1->GetResult(1);
  fMyTimeZeroTOFtracks = fTOFT0v1->GetResult(3);
  Bool_t hasTimeZeroTOF = kFALSE;
  /* check T0-TOF sigma */
  if (fMyTimeZeroTOFsigma < 250.)
    hasTimeZeroTOF = kTRUE;  
  return hasTimeZeroTOF;
}

//------------------------------------------------------
TString AliAnalysisTaskTOFqaID::GetSpeciesName(Int_t absPdgCode)
{
  //returns name of selected specie 
  TString name;
  switch (absPdgCode){
  case 11:
    name = "electron";
    break;
  case 13:
    name = "muon";
    break;
  case 211:
    name = "pion";
    break;
  case 321:
    name = "kaon";
      break;
  case 2212:
    name = "proton";
    break;
  default:
    name = "noPID";
    break;
  }
  return name.Data();
}

//-----------------------------------------------
Bool_t AliAnalysisTaskTOFqaID::SelectMCspecies(AliMCEvent * ev, AliESDtrack * track)
{
  //
  //Retrieves particle true ID from MC and selects the desired species
  //
  if ((!ev) || (!track)) {
    AliError("SelectMCspecies - Invalid object set as argument");
    return kFALSE;
  }
  
  if (fSelectedPdg==0) return kTRUE;   //if fSelectedPdg==0, no species selection is applied

  Long_t label = track->GetLabel();
  if (label<0) return kFALSE;
  
  // get number of particles
  Long_t nMC = ev->GetNumberOfTracks();
  // if label too large --> failed
  if (label>= nMC) {
    AliWarning(Form("Stack overflow: track label = %li -- stack maximum = %li", label, nMC));
    return kFALSE;
  } 
  // retrieve particle
  AliMCParticle *mcPart = (AliMCParticle *)ev->GetTrack(label);
  if (!mcPart) {// if particle = NULL --> failed
    AliWarning(Form("Stack discontinuity: label %li refers to a NULL object", label));
    return kFALSE;
  }

  Int_t pdgCode = mcPart->PdgCode();
  if (!(TMath::Abs(pdgCode)==fSelectedPdg))
    return kFALSE;
  else  
    return  kTRUE;
}

//----------------------------------------------------------------------------------
Bool_t  AliAnalysisTaskTOFqaID::ComputeMatchingEfficiency(TList* list, TString variable)
{
  //computes matching efficiency from previously filled histos
  // to be called in terminate function  
  if (!list) return kFALSE;

  TString matchedName, primaryName, xAxisTitle;
  if (variable.Contains("pt")) {
    matchedName = "hTOFmatchedESDPt";
    primaryName = "hESDprimaryTrackPt";
    xAxisTitle = "p_{T} (GeV/c)";
  }
  if (variable.Contains("eta")) {
    matchedName = "hTOFmatchedESDeta";
    primaryName = "hTOFprimaryESDeta";
    xAxisTitle = "#eta";
  }
  if (variable.Contains("phi")) {
    matchedName = "hTOFmatchedESDphi";
    primaryName = "hTOFprimaryESDphi";
    xAxisTitle = "#phi_vtx (deg)";
  }
  
  TH1F*hDummy = ((TH1F*)list->FindObject(matchedName.Data()));
  if (!hDummy) return 0;

  TH1F*hMatchingEff = (TH1F*) hDummy->Clone("hMatchingEff");
  hMatchingEff->SetNameTitle(Form("hMatchingEff_%s", variable.Data()),Form("Matching efficiency vs %s", variable.Data()));
  hMatchingEff->Divide((TH1F*) list->FindObject(primaryName.Data()));
  hMatchingEff->GetXaxis()->SetTitle(xAxisTitle.Data());
  hMatchingEff->GetYaxis()->SetRangeUser(0.0,1.0);
  hMatchingEff->GetYaxis()->SetTitle("#epsilon_{match}");
  list->AddLast(hMatchingEff);
  return 1;
}
//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::HistogramMakeUp(TH1* hist, Color_t color, Int_t markerStyle, TString drawOpt, TString newName, TString newTitle, TString xTitle, TString yTitle)
{
  //set histogram style and axes style at once
  if (!hist) return;  
  if (!newName.IsNull()) hist->SetName(newName.Data());
  if (!newTitle.IsNull()) hist->SetTitle(newTitle.Data());
  if (!xTitle.IsNull()) hist->GetXaxis()->SetTitle(xTitle.Data());
  if (!yTitle.IsNull()) hist->GetYaxis()->SetTitle(yTitle.Data());
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerSize(0.7);
  hist->SetDrawOption(drawOpt.Data());
  //hist->Sumw2();
  return;
}

//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::AddTofBaseHisto(TList *list, Int_t charge, TString suffix)
{
  //Creates histograms for monitoring TOF signal, time alignement and matching-related quantities
  if (!list){
    AliError("Invalid list passed as argument.");
    return;
  }
 
  TString cLabel; 
  if (charge == 0) cLabel.Form("all");
  else 
    if (charge<0) cLabel.Form("neg"); 
    else 
      if (charge>0) cLabel.Form("pos"); 
  
  
  TH1I* hTOFmulti = new TH1I(Form("hTOFmulti%s_%s",suffix.Data(), cLabel.Data()), Form("%s matched trk per event (|#eta|#leq%3.2f, p_{T}#geq0.3 GeV/c)", cLabel.Data(), fMatchingEtaCut), 100, 0, 100);  
  HistogramMakeUp(hTOFmulti, ((charge>0)? kRed : kBlue+2), 1, "E1", "","", "N","events");
  list->AddLast(hTOFmulti);
  
  TH1F* hTOFtime = new TH1F(Form("hTime%s_%s",suffix.Data(),cLabel.Data()), Form("%s matched trk TOF signal", cLabel.Data()), 250, 0., 610. ) ; 
  HistogramMakeUp(hTOFtime,((charge>0)? kRed+2 : kBlue+2), 1, "E1", "","", "t (ns)","tracks");
  list->AddLast(hTOFtime);
  
  TH1F* hTOFrawTime = new TH1F(Form("hRawTime%s_%s",suffix.Data(),cLabel.Data()), Form("%s matched trk TOF raw signal", cLabel.Data()), 250, 0., 610. ) ; 
  HistogramMakeUp(hTOFrawTime,((charge>0)? kRed+2 : kBlue+2), 1, "E1", "","", "t_{raw} (ns)","tracks");
  list->AddLast(hTOFrawTime);

  TH1F* hTOFtot = new TH1F(Form("hTot%s_%s",suffix.Data(),cLabel.Data()), Form("%s matched trk ToT", cLabel.Data()), 50, 0., 50. ) ; 
  HistogramMakeUp(hTOFtot,((charge>0)? kRed+2 : kBlue+2), 1, "E1", "","", "ToT (ns)","tracks");
  list->AddLast(hTOFtot);
    
  TH1F* hMatchedL  = new TH1F(Form("hMatchedL%s_%s",suffix.Data(),cLabel.Data()), Form("%s matched trk lenght", cLabel.Data()), 900, -100., 800) ; 
  HistogramMakeUp(hMatchedL,((charge>0)? kRed+2 : kBlue+2), 1, "E1", "","", "L (cm)","tracks");
  list->AddLast(hMatchedL);
  
  const Int_t nBinsPt = 300;
  Double_t xBins[nBinsPt+1];
  for (Int_t j=0;j<nBinsPt+1; j++) {  
    if (j<200) xBins[j] = j*0.025;
    else xBins[j] = 5.0 + (j-200)*0.050;  
  }

  TH2F* hMatchedDxVsPt = new TH2F(Form("hMatchedDxVsPt%s_%s",suffix.Data(),cLabel.Data()), Form("%s matched trk dx vs.p_{T}", cLabel.Data()), nBinsPt, xBins, 200, -10., 10.); 
  HistogramMakeUp(hMatchedDxVsPt,((charge>0)? kRed+2 : kBlue+2), 1, "colz", "","", "p_{T} (GeV/c)","dx (cm)");
  list->AddLast(hMatchedDxVsPt); 

  TH2F* hMatchedDzVsStrip = new TH2F(Form("hMatchedDzVsStrip%s_%s",suffix.Data(),cLabel.Data()), Form("%s matched trk dz vs. strip (#eta)", cLabel.Data()), 92, 0., 92., 200, -10., 10.) ; 
  HistogramMakeUp(hMatchedDzVsStrip,((charge>0)? kRed+2 : kBlue+2), 1, "colz", "","", "strip index","dz (cm)");
  list->AddLast(hMatchedDzVsStrip) ; 

  TProfile *hMatchedDxVsCh = new TProfile(Form("hMatchedDxVsCh%s_%s",suffix.Data(),cLabel.Data()), Form("%s matched trk dx vs. channel", cLabel.Data()), 157248., 0.,157248.);
  HistogramMakeUp(hMatchedDxVsCh,((charge>0)? kRed+2 : kBlue+2), 1, "", "","", "channel index","dx (cm)");
  list->AddLast(hMatchedDxVsCh);
   
  TProfile *hMatchedDzVsCh = new TProfile(Form("hMatchedDzVsCh%s_%s",suffix.Data(),cLabel.Data()), Form("%s matched trk dz vs. channel", cLabel.Data()), 157248., 0.,157248.);
  HistogramMakeUp(hMatchedDzVsCh,((charge>0)? kRed+2 : kBlue+2), 1, "", "","", "channel index","dz (cm)");
  list->AddLast(hMatchedDzVsCh);    

  return;
}

//----------------------------------------------------------------------------------
void    AliAnalysisTaskTOFqaID::AddMatchingEffHisto(TList *list, Int_t charge, TString suffix)
{
  if (!list){
    AliError("Invalid list passed as argument.");
    return;
  }
  TString cLabel; 
  if (charge == 0) cLabel.Form("all");
  else 
    if (charge<0) cLabel.Form("neg"); 
    else 
      if (charge>0) cLabel.Form("pos"); 

  const Int_t nBinsX = 300;
  Double_t xBins[nBinsX+1];
  for (Int_t j=0;j<nBinsX+1; j++) {  
    if (j<200) xBins[j] = j*0.025;
    else xBins[j] = 5.0 + (j-200)*0.050;  
  }

  TH1F* hMatchedP  = new TH1F(Form("hMatchedP%s_%s",suffix.Data(),cLabel.Data()), Form("%s matched trk p", cLabel.Data()), nBinsX, xBins);// 1000,0.,10.) ;  
  HistogramMakeUp(hMatchedP,((charge>0)? kRed+2 : kBlue+2), 1, "E1", "","", "p (GeV/c)","tracks");
  list->AddLast(hMatchedP) ; 

  TH1F* hMatchedPt  = new TH1F(Form("hMatchedPt%s_%s",suffix.Data(),cLabel.Data()), Form("%s matched trk p_{T}", cLabel.Data()), nBinsX, xBins);// 1000,0.,10.) ;  
  HistogramMakeUp(hMatchedPt,((charge>0)? kRed+2 : kBlue+2), 1, "E1", "","", "p_{T} (GeV/c)","tracks");
  list->AddLast(hMatchedPt) ; 

  TH1F* hMatchedEta = new TH1F(Form("hMatchedEta%s_%s",suffix.Data(),cLabel.Data()), Form("%s matched trk #eta", cLabel.Data()), 200, -1., 1.) ; 
  HistogramMakeUp(hMatchedEta,((charge>0)? kRed+2 : kBlue+2), 1, "E1", "","", "#eta","tracks");
  list->AddLast(hMatchedEta) ; 

  TH1F* hMatchedPhi = new TH1F(Form("hMatchedPhi%s_%s",suffix.Data(),cLabel.Data()), Form("%s matched trk #phi_{vtx}", cLabel.Data()), 72, 0., 360.) ; 
  HistogramMakeUp(hMatchedPhi,((charge>0)? kRed+2 : kBlue+2), 1, "E1", "","", "#phi_{vtx} (deg)","tracks");
  list->AddLast(hMatchedPhi) ; 

  TH2F* hMatchedPtVsOutPhi  = new TH2F(Form("hMatchedPtVsOutPhi%s_%s",suffix.Data(),cLabel.Data()), Form("%s matched trk p_{T} vs. #phi_{TPC out}", cLabel.Data()), 72, 0.0, 360.0, nBinsX, xBins);// 1000,0.,10.) ;  
  HistogramMakeUp(hMatchedPtVsOutPhi,((charge>0)? kRed+2 : kBlue+2), 1, "colz", "","", "#phi_{TPC out} (deg)","p_{T} (GeV/c)");
  list->AddLast(hMatchedPtVsOutPhi) ;
   
  TH1F* hPrimaryP  = new TH1F(Form("hPrimaryP%s_%s",suffix.Data(),cLabel.Data()), Form("%s primary trk p", cLabel.Data()), nBinsX, xBins);// 1000,0.,10.) ;  
  HistogramMakeUp(hPrimaryP,((charge>0)? kRed+2 : kBlue+2), 1, "E1", "","", "p (GeV/c)","tracks");
  list->AddLast(hPrimaryP) ; 

  TH1F* hPrimaryPt  = new TH1F(Form("hPrimaryPt%s_%s",suffix.Data(),cLabel.Data()), Form("%s primary trk p_{T}", cLabel.Data()), nBinsX, xBins);// 1000,0.,10.) ;  
  HistogramMakeUp(hPrimaryPt,((charge>0)? kRed+2 : kBlue+2), 1, "E1", "","", "p_{T} (GeV/c)","tracks");
  list->AddLast(hPrimaryPt) ; 

  TH1F* hPrimaryEta = new TH1F(Form("hPrimaryEta%s_%s",suffix.Data(),cLabel.Data()), Form("%s primary trk #eta", cLabel.Data()), 200, -1., 1.) ; 
  HistogramMakeUp(hPrimaryEta,((charge>0)? kRed+2 : kBlue+2), 1, "E1", "","", "#eta","tracks");
  list->AddLast(hPrimaryEta) ; 

  TH1F* hPrimaryPhi = new TH1F(Form("hPrimaryPhi%s_%s",suffix.Data(),cLabel.Data()), Form("%s primary trk #phi_{vtx}", cLabel.Data()), 72, 0., 360.) ; 
  HistogramMakeUp(hPrimaryPhi,((charge>0)? kRed+2 : kBlue+2), 1, "E1", "","", "#phi_{vtx} (deg)","tracks");
  list->AddLast(hPrimaryPhi) ; 
   
  TH2F* hPrimaryPtVsOutPhi  = new TH2F(Form("hPrimaryPtVsOutPhi%s_%s",suffix.Data(),cLabel.Data()), Form("%s primary trk p_{T} vs. #phi_{TPC out}", cLabel.Data()), 72, 0.0, 360.0, nBinsX, xBins);// 1000,0.,10.) ;  
  HistogramMakeUp(hPrimaryPtVsOutPhi,((charge>0)? kRed+2 : kBlue+2), 1, "colz", "","", "#phi_{TPC out} (deg)","p_{T} (GeV/c)");
  list->AddLast(hPrimaryPtVsOutPhi) ; 
  return;
}
  
//----------------------------------------------------------------------------------
void  AliAnalysisTaskTOFqaID::AddPidHisto(TList *list, Int_t charge, TString suffix)
{
  //Creates histograms for monitoring TOF PID
  if (!list){
    AliError("Invalid list passed as argument.");
    return;
  }
  TString cLabel; 
  if (charge == 0) cLabel.Form("all");
  else 
    if (charge<0) cLabel.Form("neg"); 
    else 
      if (charge>0) cLabel.Form("pos"); 
  
  const Int_t nBinsX = 300;
  Double_t xBins[nBinsX+1];
  for (Int_t j=0;j<nBinsX+1; j++) {  
    if (j<200) xBins[j] = j*0.025;
    else xBins[j] = 5.0 + (j-200)*0.050;  
  }

  TH2F* hMatchedBetaVsP  = new TH2F(Form("hMatchedBetaVsP%s_%s",suffix.Data(),cLabel.Data()), Form("%s matched trk #beta vs. p", cLabel.Data()), nBinsX, xBins, 150, 0., 1.5) ; 
  HistogramMakeUp(hMatchedBetaVsP,((charge>0)? kRed+2 : kBlue+2), 1, "colz", "","", "p (GeV/c)","#beta");
  list->AddLast(hMatchedBetaVsP);
    
  TH1F* hMatchedMass= new TH1F(Form("hMatchedMass%s_%s",suffix.Data(),cLabel.Data()), Form("%s matched p.le M", cLabel.Data()), 500, 0., 5. );
  HistogramMakeUp(hMatchedMass,((charge>0)? kRed+2 : kBlue+2), 1, "", "","", "M (GeV/c^{2})","entries");
  list->AddLast(hMatchedMass);
    
  TH1F* hMatchedMass2= new TH1F(Form("hMatchedMass2%s_%s",suffix.Data(),cLabel.Data()), Form("%s matched p.le M^{2}", cLabel.Data()), 500, 0., 10. );
  HistogramMakeUp(hMatchedMass2,((charge>0)? kRed+2 : kBlue+2), 1, "", "","", "M^{2} (GeV^{2}/c^{4})","entries");
  list->AddLast(hMatchedMass2);

  TH2F* hExpTimePiVsStrip = new TH2F(Form("hExpTimePiVsStrip%s_%s",suffix.Data(),cLabel.Data()),Form("%s matched trk t_{TOF}-t_{#pi,exp} vs strip",cLabel.Data()), 92, 0, 92,  fnExpTimeSmallBins, fExpTimeSmallRangeMin, fExpTimeSmallRangeMax) ; 
  HistogramMakeUp(hExpTimePiVsStrip,((charge>0)? kRed+2 : kBlue+2), 1, "", "","", "strip (#eta)","t_{TOF}-t_{#pi,exp} [ps]");
  list->AddLast(hExpTimePiVsStrip);
  
  TH2F* hExpTimePiT0Sub1GeV = new TH2F(Form("hExpTimePiT0Sub1GeV%s_%s",suffix.Data(),cLabel.Data()), Form("%s trk (0.95#leq p_{T}#leq 1.05 GeV/c) t_{TOF}-t_{#pi,exp}-t_{0}^{TOF}",cLabel.Data()), 500, 0., 500., fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  HistogramMakeUp(hExpTimePiT0Sub1GeV,((charge>0)? kRed+2 : kBlue+2), 1, "colz", "","","n. tracks used for t_{0}^{TOF}","t_{TOF}-t_{#pi,exp}-t_{0}^{TOF}");    
  list->AddLast(hExpTimePiT0Sub1GeV) ;

  TH1F* hExpTimePiFillSub = new TH1F(Form("hExpTimePiFillSub%s_%s",suffix.Data(),cLabel.Data()), Form("%s trk t_{TOF}-t_{#pi,exp}-t_{0,fill}",cLabel.Data()), 6150, -75030., 75030.) ; 
  HistogramMakeUp(hExpTimePiFillSub,((charge>0)? kRed+2 : kBlue+2), 1, "", "","","t_{TOF}-t_{#pi,exp} -t_{0,fill} [ps]","entries");    
  list->AddLast(hExpTimePiFillSub) ;
  
  TH1F* hExpTimePi = new TH1F(Form("hExpTimePi%s_%s",suffix.Data(),cLabel.Data()),Form("%s matched trk t_{TOF}-t_{#pi,exp}",cLabel.Data()), 6150, -75030., 75030.) ; 
  HistogramMakeUp(hExpTimePi,((charge>0)? kRed+2 : kBlue+2), 1, "", "","", "t_{TOF}-t_{#pi,exp} [ps]","tracks");
  list->AddLast(hExpTimePi);
  
  TH2F* hExpTimePiVsP = new TH2F(Form("hExpTimePiVsP%s_%s",suffix.Data(),cLabel.Data()),Form("%s matched trk t_{TOF}-t_{#pi,exp}",cLabel.Data()), nBinsX, xBins, fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  HistogramMakeUp(hExpTimePiVsP,kRed+2, 1, "colz", "","", "p (GeV/c)","t_{TOF}-t_{#pi,exp} [ps]");
  list->AddLast(hExpTimePiVsP);
  
  TH2F* hExpTimeKaVsP = new TH2F(Form("hExpTimeKaVsP%s_%s",suffix.Data(),cLabel.Data()),Form("%s matched trk t_{TOF}-t_{K,exp}",cLabel.Data()), nBinsX, xBins, fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  HistogramMakeUp(hExpTimeKaVsP,kBlue+2, 1, "colz", "","", "p (GeV/c)","t_{TOF}-t_{K,exp} [ps]");
  list->AddLast(hExpTimeKaVsP);
  
  TH2F* hExpTimeProVsP = new TH2F(Form("hExpTimeProVsP%s_%s",suffix.Data(),cLabel.Data()),Form("%s matched trk t_{TOF}-t_{p,exp}",cLabel.Data()), nBinsX, xBins, fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  HistogramMakeUp(hExpTimeProVsP,kGreen+2, 1, "colz", "","", "p (GeV/c)","t_{TOF}-t_{p,exp} [ps]");
  list->AddLast(hExpTimeProVsP);
  
  TH2F* hTOFpidSigmaPi = new TH2F(Form("hTOFpidSigmaPi%s_%s",suffix.Data(),cLabel.Data()), Form("%s trk n#sigma^{TOF}_{#pi} vs p_{T}",cLabel.Data()), 500,0.,5.,200, -10., 10. ) ; 
  HistogramMakeUp(hTOFpidSigmaPi,kRed+2, 1, "colz", "","", "p (GeV/c)","n#sigma_{#pi,exp} [ps]");
  list->AddLast(hTOFpidSigmaPi) ;
  
  TH2F* hTOFpidSigmaKa = new TH2F(Form("hTOFpidSigmaKa%s_%s",suffix.Data(),cLabel.Data()), Form("%s trk n#sigma^{TOF}_{K} vs p_{T}",cLabel.Data()), 500, 0.,5.,200, -10., 10. ) ; 
  HistogramMakeUp(hTOFpidSigmaKa,kBlue+2, 1, "colz", "","", "p (GeV/c)","n#sigma_{K,exp} [ps]");
  list->AddLast(hTOFpidSigmaKa) ;
    
  TH2F* hTOFpidSigmaPro = new TH2F(Form("hTOFpidSigmaPro%s_%s",suffix.Data(),cLabel.Data()), Form("%s trk TOF n#sigma^{TOF}_{p} vs p_{T}",cLabel.Data()), 500, 0.,5.,200, -10., 10. ) ; 
  HistogramMakeUp(hTOFpidSigmaPro,kGreen+2, 1, "colz", "","","p (GeV/c)","n#sigma_{p,exp} [ps]");
  list->AddLast(hTOFpidSigmaPro);
    
  TH2F* hExpTimePiT0SubVsP = new TH2F(Form("hExpTimePiT0SubVsP%s_%s",suffix.Data(),cLabel.Data()), Form("%s trk t_{TOF}-t_{#pi,exp}-t_{0}^{TOF}",cLabel.Data()), nBinsX, xBins, fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  HistogramMakeUp(hExpTimePiT0SubVsP,kRed+2, 1, "colz", "","","p (GeV/c)","t_{TOF}-t_{#pi,exp}-t_{0}^{TOF}");
  list->AddLast(hExpTimePiT0SubVsP) ;
    
  TH2F* hExpTimeKaT0SubVsP = new TH2F(Form("hExpTimeKaT0SubVsP%s_%s",suffix.Data(),cLabel.Data()), Form("%s trk t_{TOF}-t_{K,exp}-t_{0}^{TOF}",cLabel.Data()), nBinsX, xBins, fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  HistogramMakeUp(hExpTimeKaT0SubVsP,kBlue+2, 1, "colz", "","","p (GeV/c)","t_{TOF}-t_{K,exp}-t_{0}^{TOF}");
  list->AddLast(hExpTimeKaT0SubVsP) ;
    
  TH2F* hExpTimeProT0SubVsP = new TH2F(Form("hExpTimeProT0SubVsP%s_%s",suffix.Data(),cLabel.Data()), Form("%s trk t_{TOF}-t_{p,exp}-t_{0}^{TOF}",cLabel.Data()), nBinsX, xBins, fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  HistogramMakeUp(hExpTimeProT0SubVsP,kGreen+2, 1, "colz", "","","p (GeV/c)","t_{TOF}-t_{p,exp}-t_{0}^{TOF}");
  list->AddLast(hExpTimeProT0SubVsP) ;
  
  TH2F* hExpTimePiVsOutPhi = new TH2F(Form("hExpTimePiVsOutPhi%s_%s",suffix.Data(),cLabel.Data()),Form("%s matched trk t_{TOF}-t_{#pi,exp} vs #phi_{TPC out}",cLabel.Data()), 72, 0.0, 360.0, fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  HistogramMakeUp(hExpTimePiVsOutPhi,kRed+2, 1, "colz", "","", "#phi_{TPC out} (deg)","t_{TOF}-t_{#pi,exp} [ps]");
  list->AddLast(hExpTimePiVsOutPhi);
  
  TH2F* hExpTimeKaVsOutPhi = new TH2F(Form("hExpTimeKaVsOutPhi%s_%s",suffix.Data(),cLabel.Data()),Form("%s matched trk t_{TOF}-t_{K,exp} vs #phi_{TPC out}",cLabel.Data()), 72, 0.0, 360.0, fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  HistogramMakeUp(hExpTimeKaVsOutPhi,kBlue+2, 1, "colz", "","", "#phi_{TPC out} (deg)","t_{TOF}-t_{K,exp} [ps]");
  list->AddLast(hExpTimeKaVsOutPhi);
  
  TH2F* hExpTimeProVsOutPhi = new TH2F(Form("hExpTimeProVsOutPhi%s_%s",suffix.Data(),cLabel.Data()),Form("%s matched trk t_{TOF}-t_{p,exp} vs #phi_{TPC out}",cLabel.Data()), 72, 0.0, 360.0, fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  HistogramMakeUp(hExpTimeProVsOutPhi,kGreen+2, 1, "colz", "","", "#phi_{TPC out} (deg)","t_{TOF}-t_{p,exp} [ps]");
  list->AddLast(hExpTimeProVsOutPhi);

  TH2F* hExpTimePiVsPgoodCh = new TH2F(Form("hExpTimePiVsPgoodCh%s_%s",suffix.Data(),cLabel.Data()),Form("%s matched trk t_{TOF}-t_{#pi,exp} - good channels",cLabel.Data()), nBinsX, xBins, fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  HistogramMakeUp(hExpTimePiVsPgoodCh,kRed+2, 1, "colz", "","", "p (GeV/c)","t_{TOF}-t_{#pi,exp} [ps]");
  list->AddLast(hExpTimePiVsPgoodCh);
  
  TH2F* hExpTimeKaVsPgoodCh = new TH2F(Form("hExpTimeKaVsPgoodCh%s_%s",suffix.Data(),cLabel.Data()),Form("%s matched trk t_{TOF}-t_{K,exp} - good channels",cLabel.Data()), nBinsX, xBins, fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  HistogramMakeUp(hExpTimeKaVsPgoodCh,kBlue+2, 1, "colz", "","", "p (GeV/c)","t_{TOF}-t_{K,exp} [ps]");
  list->AddLast(hExpTimeKaVsPgoodCh);
  
  TH2F* hExpTimeProVsPgoodCh = new TH2F(Form("hExpTimeProVsPgoodCh%s_%s",suffix.Data(),cLabel.Data()),Form("%s matched trk t_{TOF}-t_{p,exp} - good channels",cLabel.Data()), nBinsX, xBins, fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax) ; 
  HistogramMakeUp(hExpTimeProVsPgoodCh,kGreen+2, 1, "colz", "","", "p (GeV/c)","t_{TOF}-t_{p,exp} [ps]");
  list->AddLast(hExpTimeProVsPgoodCh);
  
  return;
}
//----------------------------------------------------------------------------------
void    AliAnalysisTaskTOFqaID::AddStartTimeHisto(TList *list, TString suffix)
{  
  //Creates histograms for monitoring T0 signal and start-time related quantities
  if (!list){
    AliError("Invalid list passed as argument.");
    return;
  }
  TH1F* hT0AC = new TH1F(Form("hT0AC%s",suffix.Data()), "Event timeZero from T0A&C; t_{0,AC} [ps]; events", 1000, -12500., 12500.) ; 
  HistogramMakeUp(hT0AC, kRed+2, 20, "", "","","","");    
  list->AddLast(hT0AC);

  TH1F* hT0A = new TH1F(Form("hT0A%s",suffix.Data()), "Event timeZero from T0A; t_{0,A} [ps]; events", 1000, -12500., 12500.) ; 
  HistogramMakeUp(hT0A, kBlue+2, 25, "", "","","","");    
  list->AddLast(hT0A);
    
  TH1F* hT0C = new TH1F(Form("hT0C%s",suffix.Data()), "Event timeZero from T0C; t_{0,C} [ps]; events", 1000, -12500., 12500.) ; 
  HistogramMakeUp(hT0C, kGreen+2, 28, "", "","","","");    
  list->AddLast(hT0C);
    
  TH1F* hT0DetRes = new TH1F(Form("hT0DetRes%s",suffix.Data()), "T0 detector (T0A-T0C)/2; (T0A-T0C)/2 [ps]; events", 200, -500.,500. ) ; 
  HistogramMakeUp(hT0DetRes, kMagenta+1, 1, "", "","","","");    
  list->AddLast(hT0DetRes) ; 

  TH2F* hEventT0MeanVsVtx = new TH2F(Form("hEventT0MeanVsVtx%s",suffix.Data()), "T0 detector: mean vs vertex ; (t0_{A}-t0_{C})/2 [ns]; (t0_{A}+t0_{C})/2 [ns]; events", 50, -25., 25., 50, -25., 25. );
  HistogramMakeUp(hEventT0MeanVsVtx, kBlue+2, 1, "colz", "","","","");
  list->AddLast(hEventT0MeanVsVtx) ;
    
  TH2F* hEventV0MeanVsVtx = new TH2F(Form("hEventV0MeanVsVtx%s",suffix.Data()), "V0 detector: mean vs vertex ; (V0_{A}-V0_{C})/2 [ns]; (V0_{A}+V0_{C})/2 [ns]; events", 50, -25., 25., 50, -25., 25.);
  HistogramMakeUp(hEventV0MeanVsVtx, kBlack, 1, "colz", "","","","");
  list->AddLast(hEventV0MeanVsVtx);

  TH2F* hStartTime = new TH2F(Form("hStartTime%s",suffix.Data()),"Start time for each method (mask); start time (ps); method;", 140, -700, 700, 6, -1., 5.);
  hStartTime->GetYaxis()->SetBinLabel(1,"best_t0");
  hStartTime->GetYaxis()->SetBinLabel(2,"fill_t0");
  hStartTime->GetYaxis()->SetBinLabel(3,"tof_t0");
  hStartTime->GetYaxis()->SetBinLabel(4,"T0AC");
  hStartTime->GetYaxis()->SetBinLabel(5,"T0A");
  hStartTime->GetYaxis()->SetBinLabel(6,"T0C");
  list->AddLast(hStartTime);

  TH2F* hStartTimeRes = new TH2F(Form("hStartTimeRes%s",suffix.Data()),"Start time resolution for each method (mask); resolution (ps); method;", 300, 0., 300., 6, -1., 5.);
  hStartTimeRes->GetYaxis()->SetBinLabel(1,"best_t0");
  hStartTimeRes->GetYaxis()->SetBinLabel(2,"fill_t0");
  hStartTimeRes->GetYaxis()->SetBinLabel(3,"tof_t0");
  hStartTimeRes->GetYaxis()->SetBinLabel(4,"T0AC");
  hStartTimeRes->GetYaxis()->SetBinLabel(5,"T0A");
  hStartTimeRes->GetYaxis()->SetBinLabel(6,"T0C");
  list->AddLast(hStartTimeRes);

  TH2F* hT0TOFvsNtrk = new TH2F(Form("hT0TOFvsNtrk%s",suffix.Data()), "Event timeZero estimated by TOF vs. TOF-matching tracks; N_{TOF}; t0 [ps]", 700, 0., 700., 140, -700., 700.) ; 
  HistogramMakeUp(hT0TOFvsNtrk, kTeal-5, 1, "colz", "","","","");    
  list->AddLast(hT0TOFvsNtrk) ;

  TH2F* hT0ACvsNtrk = new TH2F(Form("hT0ACvsNtrk%s",suffix.Data()), "Event timeZero estimated by T0AC vs. TOF-matching tracks;  N_{TOF}; t0 [ps]", 700, 0., 700., 140, -700., 700.) ; 
  HistogramMakeUp(hT0ACvsNtrk, kRed+2, 1, "colz", "","","","");    
  list->AddLast(hT0ACvsNtrk) ;

  TH2F* hT0AvsNtrk = new TH2F(Form("hT0AvsNtrk%s",suffix.Data()), "Event timeZero estimated by T0A vs. TOF-matching tracks;  N_{TOF}; t0 [ps]", 700, 0., 700., 140, -700., 700.) ; 
  HistogramMakeUp(hT0AvsNtrk, kBlue+2, 1, "colz", "","","","");    
  list->AddLast(hT0AvsNtrk) ; 
 
  TH2F* hT0CvsNtrk = new TH2F(Form("hT0CvsNtrk%s",suffix.Data()), "Event timeZero estimated by T0C vs. TOF-matching tracks;  N_{TOF}; t0 [ps]", 700, 0., 700., 140, -700., 700.) ; 
  HistogramMakeUp(hT0CvsNtrk, kGreen+2, 1, "colz", "","","","");    
  list->AddLast(hT0CvsNtrk) ;

  const Double_t startTimeMomBins[13]={ 0.0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.5, 2., 3., 10.};
  TH2F* hStartTimeMaskMatched = new TH2F(Form("hStartTimeMaskMatched%s",suffix.Data()),"Start Time Mask vs p bin for matched tracks; p(GeV/c);", 12, startTimeMomBins, 8,0.,8.);
  hStartTimeMaskMatched->GetYaxis()->SetBinLabel(1,"fill_t0");
  hStartTimeMaskMatched->GetYaxis()->SetBinLabel(2,"tof_t0");
  hStartTimeMaskMatched->GetYaxis()->SetBinLabel(3,"T0AC");
  hStartTimeMaskMatched->GetYaxis()->SetBinLabel(4,"T0AC & tof_t0");
  hStartTimeMaskMatched->GetYaxis()->SetBinLabel(5,"T0A");
  hStartTimeMaskMatched->GetYaxis()->SetBinLabel(6,"T0A & tof_t0");
  hStartTimeMaskMatched->GetYaxis()->SetBinLabel(7,"T0C");
  hStartTimeMaskMatched->GetYaxis()->SetBinLabel(8,"T0C & tof_t0");
  HistogramMakeUp(hStartTimeMaskMatched, kRed+2, 1, "colz", "","","","");    
  list->AddLast(hStartTimeMaskMatched);
    
  TH2F* hStartTimeMask = new TH2F(Form("hStartTimeMask%s",suffix.Data()),"Start Time Mask vs p bin for primary tracks; p(GeV/c);", 12, startTimeMomBins, 8,0.,8.);
  hStartTimeMask->GetYaxis()->SetBinLabel(1,"fill_t0");
  hStartTimeMask->GetYaxis()->SetBinLabel(2,"tof_t0");
  hStartTimeMask->GetYaxis()->SetBinLabel(3,"T0AC");
  hStartTimeMask->GetYaxis()->SetBinLabel(4,"T0AC & tof_t0");
  hStartTimeMask->GetYaxis()->SetBinLabel(5,"T0A");
  hStartTimeMask->GetYaxis()->SetBinLabel(6,"T0A & tof_t0");
  hStartTimeMask->GetYaxis()->SetBinLabel(7,"T0C");
  hStartTimeMask->GetYaxis()->SetBinLabel(8,"T0C & tof_t0");
  HistogramMakeUp(hStartTimeMask, kRed+2, 1, "colz", "","","","");    
  list->AddLast(hStartTimeMask);

  return;
}
//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::AddTrdHisto()
{
  //Creates histograms for monitoring TOF base quantities wrt TRD/no TRD selection
  if (!fHlistTRD){
    AliError("Invalid TRD list");
    return;
  }

  if (fEnableChargeSplit) {
    AddMatchingEffHisto(fHlistTRD, 1, "_noTrd");
    AddMatchingEffHisto(fHlistTRD, -1, "_noTrd");
    AddMatchingEffHisto(fHlistTRD, 1, "_Trd");
    AddMatchingEffHisto(fHlistTRD, -1, "_Trd");
    
    AddPidHisto(fHlistTRD, 1, "_noTrd");
    AddPidHisto(fHlistTRD, -1, "_noTrd");
    AddPidHisto(fHlistTRD, 1, "_Trd");
    AddPidHisto(fHlistTRD, -1, "_Trd");
  } else {
    AddMatchingEffHisto(fHlistTRD, 0, "_noTrd");
    AddMatchingEffHisto(fHlistTRD, 0, "_Trd");
    AddPidHisto(fHlistTRD, 0, "_noTrd");
    AddPidHisto(fHlistTRD, 0, "_Trd");
  }

  return;
}
 

//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::AddTofTrgHisto(TString suffix)
{
  //defines histo with trigger info
  if (!fHlistTrigger){
    AliError("Invalid TOF trigger list");
    return;
  }
  
  TH1I* hFiredMaxipad = new TH1I(Form("hFiredMaxipad%s",suffix.Data()), Form("Fired maxipad per event"), 1584, 0, 1584);  
  HistogramMakeUp(hFiredMaxipad, kBlue+2, 1, "E1", "","", "N_{maxipad}","events");
  fHlistTrigger->AddLast(hFiredMaxipad);
  
  TH1I* hFiredReadoutPad = new TH1I(Form("hFiredReadoutPad%s",suffix.Data()), Form("Fired readout pad per event"), 153000, 0, 153000);  
  HistogramMakeUp(hFiredReadoutPad, kRed+2, 1, "E1", "","", "N_{pad}","events");
  fHlistTrigger->AddLast(hFiredReadoutPad);
  
  TH1I* hFiredReadoutTrgPad = new TH1I(Form("hFiredReadoutTrgPad%s",suffix.Data()), Form("Fired readout pad in trg window"), 153000, 0, 153000);  
  HistogramMakeUp(hFiredReadoutTrgPad, kBlack, 1, "E1", "","", "N_{pad} in trg window","events");
  fHlistTrigger->AddLast(hFiredReadoutTrgPad);
  
  TH2I* hFiredMaxipadVsTrgPad = new TH2I(Form("hFiredMaxipadVsTrgPad%s",suffix.Data()), Form("Fired maxipad vs pads in trg window per event"), 100, 0, 100, 100, 0, 100);  
  HistogramMakeUp(hFiredMaxipadVsTrgPad, kBlue+2, 1, "colz", "","", "N_{pad} in trg window","N_{maxipad}");
  fHlistTrigger->AddLast(hFiredMaxipadVsTrgPad);
  
  TH2I* hTrgMap = new TH2I(Form("hTrgMap%s",suffix.Data()), Form("Map of fired maxipads"), 72, 0, 72, 23, 0, 23);  
  HistogramMakeUp(hTrgMap, kBlue+2, 1, "colz", "","", "LTM","maxipad");
  fHlistTrigger->AddLast(hTrgMap);
  
  return;
}

//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::FillTofBaseHisto(AliESDtrack * track, Int_t charge, TString suffix)
{
  //fill histo with TOF base quantities
  if (!track) return;
  
  //  Double_t tofTime=track->GetTOFsignal();//in ps
  Double_t tofTimeRaw=track->GetTOFsignalRaw();//in ps
  Double_t tofToT=track->GetTOFsignalToT(); //in ps
  Int_t channel=track->GetTOFCalChannel(); 
  Int_t volId[5]; //(sector, plate,strip,padZ,padX)
  AliTOFGeometry::GetVolumeIndices(channel,volId);
  TString cLabel;
  if (charge == 0) cLabel.Form("all");
  else 
    if (charge<0) cLabel.Form("neg"); 
    else 
      if (charge>0) cLabel.Form("pos"); 

  ((TH1F*)fHlist->FindObject(Form("hTime%s_%s",suffix.Data(),cLabel.Data())))->Fill(fTof); //ns
  ((TH1F*)fHlist->FindObject(Form("hRawTime%s_%s",suffix.Data(),cLabel.Data())))->Fill(tofTimeRaw*1E-3); //ns
  ((TH1F*)fHlist->FindObject(Form("hTot%s_%s",suffix.Data(),cLabel.Data())))->Fill(tofToT);
  ((TH1F*)fHlist->FindObject(Form("hMatchedL%s_%s", suffix.Data(), cLabel.Data())))->Fill(fL);      
  ((TH2F*)fHlist->FindObject(Form("hMatchedDxVsPt%s_%s",suffix.Data(),cLabel.Data())))->Fill(fPt,track->GetTOFsignalDx());
  ((TH2F*)fHlist->FindObject(Form("hMatchedDzVsStrip%s_%s",suffix.Data(),cLabel.Data())))->Fill((Int_t)GetStripIndex(volId),track->GetTOFsignalDz());
  ((TProfile*)fHlist->FindObject(Form("hMatchedDxVsCh%s_%s",suffix.Data(),cLabel.Data())))->Fill(channel,track->GetTOFsignalDx());
  ((TProfile*)fHlist->FindObject(Form("hMatchedDzVsCh%s_%s",suffix.Data(),cLabel.Data())))->Fill(channel,track->GetTOFsignalDz());
  
  return;
}
//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::FillPrimaryTrkHisto(Int_t charge, TString suffix)
{
  // fill histos with primary tracks info
  // => denominator for matching efficiency
  TString cLabel; 
  if (charge == 0) cLabel.Form("all");
  else 
    if (charge<0) cLabel.Form("neg"); 
    else 
      if (charge>0) cLabel.Form("pos"); 
  
  if (suffix.Contains("Trd")) {
    ((TH1F*)fHlistTRD->FindObject(Form("hPrimaryP%s_%s",suffix.Data(),cLabel.Data())))->Fill(fP); 
    ((TH1F*)fHlistTRD->FindObject(Form("hPrimaryPt%s_%s",suffix.Data(),cLabel.Data())))->Fill(fPt); 
    ((TH2F*)fHlistTRD->FindObject(Form("hPrimaryPtVsOutPhi%s_%s",suffix.Data(),cLabel.Data())))->Fill(fTPCOuterPhi,fPt);
    if (fPt>=fMatchingMomCut) {
      ((TH1F*)fHlistTRD->FindObject(Form("hPrimaryEta%s_%s",suffix.Data(),cLabel.Data())))->Fill(fEta);
      ((TH1F*)fHlistTRD->FindObject(Form("hPrimaryPhi%s_%s",suffix.Data(),cLabel.Data())))->Fill(fPhi);
    }
  } else {
    ((TH1F*)fHlist->FindObject(Form("hPrimaryP%s_%s",suffix.Data(),cLabel.Data())))->Fill(fP); 
    ((TH1F*)fHlist->FindObject(Form("hPrimaryPt%s_%s",suffix.Data(),cLabel.Data())))->Fill(fPt); 
    ((TH2F*)fHlist->FindObject(Form("hPrimaryPtVsOutPhi%s_%s",suffix.Data(),cLabel.Data())))->Fill(fTPCOuterPhi,fPt); 
    if (fPt>=fMatchingMomCut) {
      ((TH1F*)fHlist->FindObject(Form("hPrimaryEta%s_%s",suffix.Data(),cLabel.Data())))->Fill(fEta);
      ((TH1F*)fHlist->FindObject(Form("hPrimaryPhi%s_%s",suffix.Data(),cLabel.Data())))->Fill(fPhi);
    }
  }
  return;
}
//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::FillMatchedTrkHisto(Int_t charge, TString suffix)
{
  //get matched tracks variables (matching cut to be applied externally)
  //=> numerator for matching efficiency
  TString cLabel; 
  if (charge == 0) cLabel.Form("all");
  else 
    if (charge<0) cLabel.Form("neg"); 
    else 
      if (charge>0) cLabel.Form("pos"); 
  
  if (suffix.Contains("Trd")) { 
    ((TH1F*)fHlistTRD->FindObject(Form("hMatchedP%s_%s",suffix.Data(),cLabel.Data())))->Fill(fP); 
    ((TH1F*)fHlistTRD->FindObject(Form("hMatchedPt%s_%s",suffix.Data(),cLabel.Data())))->Fill(fPt); 
    ((TH2F*)fHlistTRD->FindObject(Form("hMatchedPtVsOutPhi%s_%s",suffix.Data(),cLabel.Data())))->Fill(fTPCOuterPhi,fPt);
    if (fPt>=fMatchingMomCut) {
      ((TH1F*)fHlistTRD->FindObject(Form("hMatchedEta%s_%s",suffix.Data(),cLabel.Data())))->Fill(fEta);
      ((TH1F*)fHlistTRD->FindObject(Form("hMatchedPhi%s_%s",suffix.Data(),cLabel.Data())))->Fill(fPhi);
    }
  } else {
    ((TH1F*)fHlist->FindObject(Form("hMatchedP%s_%s",suffix.Data(),cLabel.Data())))->Fill(fP); 
    ((TH1F*)fHlist->FindObject(Form("hMatchedPt%s_%s",suffix.Data(),cLabel.Data())))->Fill(fPt); 
    ((TH2F*)fHlist->FindObject(Form("hMatchedPtVsOutPhi%s_%s",suffix.Data(),cLabel.Data())))->Fill(fTPCOuterPhi,fPt);
    if (fPt>=fMatchingMomCut) {
      ((TH1F*)fHlist->FindObject(Form("hMatchedEta%s_%s",suffix.Data(),cLabel.Data())))->Fill(fEta);
      ((TH1F*)fHlist->FindObject(Form("hMatchedPhi%s_%s",suffix.Data(),cLabel.Data())))->Fill(fPhi);
    }
  }
  return;
}

//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::FillPidHisto(AliESDtrack * track, Int_t charge, TString suffix)
{
  //basic PID performance check
  if (fTof<=0) {
    printf("WARNING: track with negative TOF time found! Skipping this track for PID checks\n");
    return;
  }
  if (fL<=0){
    printf("WARNING: track with negative length found!Skipping this track for PID checks\n");
    return;
  }
  if (!track) return;
  
  TString cLabel; 
  if (charge == 0) cLabel.Form("all");
  else 
    if (charge<0) cLabel.Form("neg"); 
    else 
      if (charge>0) cLabel.Form("pos"); 
  
  //calculate beta
  Double_t c=TMath::C()*1.E-9;// m/ns
  Double_t mass=0.; //GeV
  Double_t length=fL*0.01; // in meters
  Double_t tof=fTof*c;
  Double_t beta=length/tof;
  Double_t fact= (tof/length)*(tof/length) -1.;
  Double_t fP2 = fP*fP;
  
  if(fact<=0) {
    mass = -fP*TMath::Sqrt(-fact);
  } else { 
    mass = fP*TMath::Sqrt(fact); 
  }
  
  if (suffix.Contains("Trd")) {
    ((TH2F*) fHlistTRD->FindObject(Form("hMatchedBetaVsP%s_%s",suffix.Data(),cLabel.Data())))->Fill(fP,beta);
    ((TH1F*) fHlistTRD->FindObject(Form("hMatchedMass%s_%s",suffix.Data(),cLabel.Data())))->Fill(mass);
    ((TH1F*) fHlistTRD->FindObject(Form("hMatchedMass2%s_%s",suffix.Data(),cLabel.Data())))->Fill(mass*mass);
  } else {
    ((TH2F*) fHlistPID->FindObject(Form("hMatchedBetaVsP%s_%s",suffix.Data(),cLabel.Data())))->Fill(fP,beta);
    ((TH1F*) fHlistPID->FindObject(Form("hMatchedMass%s_%s",suffix.Data(),cLabel.Data())))->Fill(mass);
    ((TH1F*) fHlistPID->FindObject(Form("hMatchedMass2%s_%s",suffix.Data(),cLabel.Data())))->Fill(mass*mass);
  }
  
  //PID sigmas
  Bool_t isValidBeta[AliPID::kSPECIES]={0,0,0,0,0};
  for (Int_t specie = 0; specie < AliPID::kSPECIES; specie++){
    fSigmaSpecie[specie] = fESDpid->GetTOFResponse().GetExpectedSigma(fP, fTrkExpTimes[specie], AliPID::ParticleMass(specie));
    beta=1/TMath::Sqrt(1+AliPID::ParticleMass(specie)*AliPID::ParticleMass(specie)/(fP2));
    if (beta>0) {
      fThExpTimes[specie]=length*1.E3/(beta*c);//ps
      isValidBeta[specie]=kTRUE;
    } else {
      fThExpTimes[specie]=1E-10;
      isValidBeta[specie]=kFALSE;
    }
  }
  Float_t timeZeroTOF = (Float_t) fESDpid->GetTOFResponse().GetStartTime(fPt);
  Double_t tofps=fTof*1E3;//ps for t-texp
  Int_t channel=track->GetTOFCalChannel(); 
  Int_t volId[5]; //(sector, plate,strip,padZ,padX)
  AliTOFGeometry::GetVolumeIndices(channel,volId);
  Char_t partName[3][4] = {"Pi","Ka","Pro"}; 
  
  if (suffix.Contains("Trd")) {
    //fill histos for pion only
    ((TH2F*)fHlistTRD->FindObject(Form("hExpTimePiVsStrip%s_%s",suffix.Data(),cLabel.Data())))->Fill((Int_t)GetStripIndex(volId),tofps-fTrkExpTimes[AliPID::kPion]);//ps
    ((TH1F*)fHlistTRD->FindObject(Form("hExpTimePi%s_%s",suffix.Data(),cLabel.Data())))->Fill(tofps-fTrkExpTimes[AliPID::kPion]);//ps	
    if (ComputeTimeZeroByTOF1GeV() && (fPt>0.95) && (fPt<1.05) ){
      ((TH2F*)fHlistTRD->FindObject(Form("hExpTimePiT0Sub1GeV%s_%s",suffix.Data(),cLabel.Data())))->Fill(fMyTimeZeroTOFtracks,tofps-fMyTimeZeroTOF-fTrkExpTimes[AliPID::kPion]);
    }
    //fill sigmas and deltas for each specie
    for (Int_t specie = AliPID::kPion; specie <= AliPID::kProton; specie++){
      if (isValidBeta[specie]){
	((TH2F*)fHlistTRD->FindObject(Form("hExpTime%sVsP%s_%s",partName[specie-2], suffix.Data(),cLabel.Data())))->Fill(fP, tofps-fTrkExpTimes[specie]);
	((TH2F*)fHlistTRD->FindObject(Form("hTOFpidSigma%s%s_%s",partName[specie-2], suffix.Data(),cLabel.Data())))->Fill(fPt, (tofps-fTrkExpTimes[specie])/fSigmaSpecie[specie]);
	((TH2F*)fHlistTRD->FindObject(Form("hExpTime%sT0SubVsP%s_%s",partName[specie-2], suffix.Data(),cLabel.Data())))->Fill(fP,tofps-fTrkExpTimes[specie]-timeZeroTOF);  
	((TH2F*)fHlistTRD->FindObject(Form("hExpTime%sVsOutPhi%s_%s",partName[specie-2], suffix.Data(),cLabel.Data())))->Fill(fTPCOuterPhi,tofps-fTrkExpTimes[specie]-timeZeroTOF);
	  
      }// end check on beta
    }
  } else { 
    
    //fill histos for pion only
    ((TH2F*)fHlistPID->FindObject(Form("hExpTimePiVsStrip%s_%s",suffix.Data(),cLabel.Data())))->Fill((Int_t)GetStripIndex(volId),tofps-fTrkExpTimes[AliPID::kPion]);//ps
    ((TH1F*)fHlistPID->FindObject(Form("hExpTimePi%s_%s",suffix.Data(),cLabel.Data())))->Fill(tofps-fTrkExpTimes[AliPID::kPion]);//ps    
    if (ComputeTimeZeroByTOF1GeV() && (fPt>0.95) && (fPt<1.05) ){
      ((TH2F*)fHlistPID->FindObject(Form("hExpTimePiT0Sub1GeV%s_%s",suffix.Data(),cLabel.Data())))->Fill(fMyTimeZeroTOFtracks,tofps-fMyTimeZeroTOF-fTrkExpTimes[AliPID::kPion]);
    }
    //fill sigmas and deltas for each specie
    for (Int_t specie = AliPID::kPion; specie <= AliPID::kProton; specie++){
      if (isValidBeta[specie]){
	((TH2F*)fHlistPID->FindObject(Form("hExpTime%sVsP%s_%s",partName[specie-2], suffix.Data(),cLabel.Data())))->Fill(fP, tofps-fTrkExpTimes[specie]);
	((TH2F*)fHlistPID->FindObject(Form("hTOFpidSigma%s%s_%s",partName[specie-2], suffix.Data(),cLabel.Data())))->Fill(fPt, (tofps-fTrkExpTimes[specie])/fSigmaSpecie[specie]);
	((TH2F*)fHlistPID->FindObject(Form("hExpTime%sT0SubVsP%s_%s",partName[specie-2], suffix.Data(),cLabel.Data())))->Fill(fP,tofps-fTrkExpTimes[specie]-timeZeroTOF);  
	((TH2F*)fHlistPID->FindObject(Form("hExpTime%sVsOutPhi%s_%s",partName[specie-2], suffix.Data(),cLabel.Data())))->Fill(fTPCOuterPhi,tofps-fTrkExpTimes[specie]-timeZeroTOF);

	//fill deltas only for good channels
	if (IsChannelGood(track->GetTOFCalChannel()))
	  ((TH2F*)fHlistPID->FindObject(Form("hExpTime%sVsPgoodCh%s_%s",partName[specie-2], suffix.Data(),cLabel.Data())))->Fill(fP, tofps-fTrkExpTimes[specie]);
	
      }// end check on beta
    } 
 
  }
  
  //re-set response kFILL_T0 to check post-alignment wih OADB
  fESDpid->SetTOFResponse(fESD,AliESDpid::kFILL_T0);//(fill_t0, tof_t0, t0_t0, best_t0)
  Float_t startTimeFill=fESDpid->GetTOFResponse().GetStartTime(fP); //timeZero for bin pT>10GeV/c
  if (suffix.Contains("Trd"))
    ((TH1F*)fHlistTRD->FindObject(Form("hExpTimePiFillSub%s_%s",suffix.Data(),cLabel.Data())))->Fill(tofps-fTrkExpTimes[AliPID::kPion]-startTimeFill);//ps
  else 
    ((TH1F*)fHlistPID->FindObject(Form("hExpTimePiFillSub%s_%s",suffix.Data(),cLabel.Data())))->Fill(tofps-fTrkExpTimes[AliPID::kPion]-startTimeFill);//ps
 
  // if (fEnableAdvancedCheck && (fPt<1.)) {
  //   Double_t pos[3]={0.,0.,0.};
  //   track->GetXYZAt(378.,5.,pos);
  //   if ((pos[0]==0.)&&(pos[1]==0.)&&(pos[2]==0.))continue;
  
  //   Double_t phiTOF=TMath::ATan2(pos[1],pos[0])*TMath::RadToDeg();
  //   if (phiTOF<0) phiTOF+= (2*TMath::Pi()*TMath::RadToDeg());
  
  //   //fill t-texp vs phi@TOF
  //    if ((phiOuterTPC<=75) || ((phiOuterTPC>=125)&&(phiOuterTPC<=235)) || (phiOuterTPC>=305) ) { //TRD sectors 2012
  //    if ( ((phiOuterTPC>=85)&&(phiOuterTPC<=115)) || ((phiOuterTPC>=245)&&(phiOuterTPC<=295)) ) {//no TRD sectors 2012
  // }
  return;      
}
//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::FillStartTimeHisto(TString suffix)
{
  //fill start time histo
  if (!fESD) {
    AliError("Invalid event object");
    return;
  }
  // info from V0 detector QA 
  AliESDVZERO * vzero = fESD->GetVZEROData();
  Float_t V0Atime = vzero->GetV0ATime();
  Float_t V0Ctime = vzero->GetV0CTime(); 
  ((TH2F*)fHlistTimeZero->FindObject(Form("hEventV0MeanVsVtx%s",suffix.Data())))->Fill((V0Atime-V0Ctime)*0.5,(V0Atime+V0Ctime)*0.5);
  
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
    ((TH1F*)fHlistTimeZero->FindObject(Form("hT0DetRes%s",suffix.Data())))->Fill((fT0[2]-fT0[1])*0.5);
    ((TH1F*)fHlistTimeZero->FindObject(Form("hT0AC%s",suffix.Data())))->Fill(fT0[0]);  
    ((TH2F*)fHlistTimeZero->FindObject(Form("hEventT0MeanVsVtx%s",suffix.Data())))->Fill( ((fT0[2]-fT0[1])*0.5e-3), ((fT0[2]+fT0[1])*0.5e-3) );
  } 
  if(TMath::Abs(fT0[1]) < t0cut){
    ((TH1F*)fHlistTimeZero->FindObject(Form("hT0A%s",suffix.Data())))->Fill(fT0[1]);   
  }
  if(TMath::Abs(fT0[2]) < t0cut){
    ((TH1F*)fHlistTimeZero->FindObject(Form("hT0C%s",suffix.Data())))->Fill(fT0[2]);
  }
  
  //Fill T0fill plots
  (fESDpid->GetTOFResponse()).ResetT0info();
  fESDpid->SetTOFResponse(fESD, AliESDpid::kFILL_T0);
  Double_t timeZero = fESDpid->GetTOFResponse().GetStartTime(10.); 
  Double_t timeZeroRes = fESDpid->GetTOFResponse().GetStartTimeRes(10.);
  ((TH1F*)(fHlistTimeZero->FindObject("hStartTime")))->Fill(timeZero, AliESDpid::kFILL_T0); //will be 0 by definition
  ((TH1F*)(fHlistTimeZero->FindObject("hStartTimeRes")))->Fill(timeZeroRes, AliESDpid::kFILL_T0);//taken from the header
  
  //fill TOF_T0 plots
  (fESDpid->GetTOFResponse()).ResetT0info();
  fESDpid->SetTOFResponse(fESD, AliESDpid::kTOF_T0);
  Int_t startTimeMaskBin = fESDpid->GetTOFResponse().GetT0binMask(9);
  if (startTimeMaskBin>0) {//timeZeroMask for 10th bin, corresponding to p>10GeV/c
    //fill plot only when there is a start time other than fill_t0
    timeZero = fESDpid->GetTOFResponse().GetStartTime(10.); //timeZero for bin p>10GeV/c
    timeZeroRes = fESDpid->GetTOFResponse().GetStartTimeRes(10.); //timeZero for bin p>10GeV/c
    ((TH1F*)(fHlistTimeZero->FindObject("hStartTime")))->Fill(timeZero, AliESDpid::kTOF_T0);
    ((TH1F*)(fHlistTimeZero->FindObject("hStartTimeRes")))->Fill(timeZeroRes, AliESDpid::kTOF_T0);
    ((TH2F*)fHlistTimeZero->FindObject("hT0TOFvsNtrk"))->Fill(fNTOFtracks[0], timeZero);
  }

  //fill T0_T0 plots
  (fESDpid->GetTOFResponse()).ResetT0info();
  fESDpid->SetTOFResponse(fESD, AliESDpid::kT0_T0);
  startTimeMaskBin = fESDpid->GetTOFResponse().GetT0binMask(9);
  if (startTimeMaskBin>0) {//timeZeroMask for 10th bin, corresponding to p>10GeV/c
    //fill plot only when there is a start time other than fill_t0
    timeZero = fESDpid->GetTOFResponse().GetStartTime(10.); //timeZero for bin p>10GeV/c
    timeZeroRes = fESDpid->GetTOFResponse().GetStartTimeRes(10.); //timeZero for bin p>10GeV/c
    if (startTimeMaskBin==2) {
      ((TH1F*)(fHlistTimeZero->FindObject("hStartTime")))->Fill(timeZero, 2);
      ((TH1F*)(fHlistTimeZero->FindObject("hStartTimeRes")))->Fill(timeZeroRes, 2);
      ((TH2F*)fHlistTimeZero->FindObject("hT0AvsNtrk"))->Fill(fNTOFtracks[0],timeZero);
    } else if (startTimeMaskBin==4) {
      ((TH1F*)(fHlistTimeZero->FindObject("hStartTime")))->Fill(timeZero, 3);
      ((TH1F*)(fHlistTimeZero->FindObject("hStartTimeRes")))->Fill(timeZeroRes, 3);
      ((TH2F*)fHlistTimeZero->FindObject("hT0CvsNtrk"))->Fill(fNTOFtracks[0],timeZero);
    } else if (startTimeMaskBin==6) {
      ((TH1F*)(fHlistTimeZero->FindObject("hStartTime")))->Fill(timeZero, 4);
      ((TH1F*)(fHlistTimeZero->FindObject("hStartTimeRes")))->Fill(timeZeroRes, 4);
      ((TH2F*)fHlistTimeZero->FindObject("hT0ACvsNtrk"))->Fill(fNTOFtracks[0],timeZero);
    }
  }
   
  //response set to best_t0
  (fESDpid->GetTOFResponse()).ResetT0info();
  fESDpid->SetTOFResponse(fESD, AliESDpid::kBest_T0);
  //fill plot only when there is a start time other than fill_t0
  timeZero = fESDpid->GetTOFResponse().GetStartTime(10.); //timeZero for bin p>10GeV/c
  timeZeroRes = fESDpid->GetTOFResponse().GetStartTimeRes(10.); //timeZero for bin p>10GeV/c
  ((TH1F*)(fHlistTimeZero->FindObject("hStartTime")))->Fill(timeZero, -1);
  ((TH1F*)(fHlistTimeZero->FindObject("hStartTimeRes")))->Fill(timeZeroRes, -1);

  //Fill mask when best_t0 is set
  FillStartTimeMaskHisto(suffix.Data());

  return;
}
//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::FillTrdHisto(AliESDtrack * track, Int_t charge)
{
  //fill histograms for TRD/noTRD
  if (!track){
    AliError("Invalid track object");
    return;
  }
  
  if (IsInTRD(track)){
    FillPrimaryTrkHisto(charge,"_Trd");
    if (IsTPCTOFMatched(track)) {      
      FillMatchedTrkHisto(charge,"_Trd");
      FillPidHisto(track,charge, "_Trd");
    }
  } else {
    FillPrimaryTrkHisto(charge,"_noTrd");    
    if (IsTPCTOFMatched(track)) {      
      FillMatchedTrkHisto(charge,"_noTrd");
      FillPidHisto(track, charge, "_noTrd");
    }
  }
  return;
}
//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::FillTofTrgHisto(TString suffix)
{
  //fills histo with trigger info
  if (!fHlistTrigger){
    AliError("Invalid TOF trigger list");
    return;
  }
  if (!fTOFHeader) {
    AliWarning("Invalid AliTOFHeader object - cannot fill trg histo");
    return;    
  }
  
  Int_t nPad = fTOFHeader->GetNumberOfTOFclusters(); //all fired readout pads
  Int_t nTrgPad = fTOFHeader->GetNumberOfTOFtrgPads();// fired readout pads in the trigger window
  Int_t nMaxiPad = fTOFHeader->GetNumberOfTOFmaxipad(); //fired maxipads
  
  // update histo with fired macropads
  AliTOFTriggerMask *trgMask = fTOFHeader->GetTriggerMask();
  for(Int_t j=0;j<72;j++){
    for(Int_t i=22;i>=0;i--){
      if(trgMask->IsON(j,i))
	((TH1I*)fHlistTrigger->FindObject(Form("hTrgMap%s",suffix.Data())))->Fill(j+1,i+1); 
    }	
  }
  ((TH1I*)fHlistTrigger->FindObject(Form("hFiredMaxipad%s",suffix.Data())))->Fill(nMaxiPad);
  ((TH1I*)fHlistTrigger->FindObject(Form("hFiredReadoutPad%s",suffix.Data())))->Fill(nPad);
  ((TH1I*)fHlistTrigger->FindObject(Form("hFiredReadoutTrgPad%s",suffix.Data())))->Fill(nTrgPad);
  ((TH2I*)fHlistTrigger->FindObject(Form("hFiredMaxipadVsTrgPad%s",suffix.Data())))->Fill(nTrgPad,nMaxiPad);
  
  return;
}

//-----------------------------------------------------------
void AliAnalysisTaskTOFqaID::LoadChannelMapsFromOCDB()
{
  //method to get the channel maps from the OCDB
  // it is called at the CreatedOutputObject stage
  // to comply with the CAF environment

  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(fOCDBLocation.Data());
  cdb->SetRun(fRunNumber);  

  if(!cdb){
    AliWarning("No CDB MANAGER, maps can not be loaded");
    return;
  }
  AliCDBPath cdbPath("TOF/Calib/Status");
  AliCDBEntry *cdbe = cdb->Get(cdbPath, fRunNumber);
  if (!cdbe) {
    printf("invalid CDB entry\n");
    fChannelArray = 0x0;
    return;
  }

  fChannelArray = (AliTOFChannelOnlineStatusArray *)cdbe->GetObject();
  fCalib = new AliTOFcalib();
  fCalib->Init();
  return;
}


//-----------------------------------------------------------
Bool_t AliAnalysisTaskTOFqaID::IsChannelGood(Int_t channel = -1)
{
  // methos to check the given channel
  // against the status maps
  // taken from the OCDB

  if (channel<0) 
    return kFALSE;
  if (!fChannelArray || !fCalib) {
    AliError("Array with TOF channel status from OCDB not available.");
    return kFALSE;
  }
  if ( (fChannelArray->GetNoiseStatus(channel) == AliTOFChannelOnlineStatusArray::kTOFNoiseBad) ||
       fCalib->IsChannelProblematic(channel) ) return kFALSE;
  
  return kTRUE;
}

#endif
