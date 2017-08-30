/**************************************************************************
 * Contributors are not mentioned at all.                                 *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright noticxse appears in all *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//-----------------------------------------------------------------
//                 AliAnalysisTaskNucleiv2pPb class
//-----------------------------------------------------------------

class TTree;
class TParticle;
class TVector3;

#include "AliAnalysisManager.h"
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0; 

#include <iostream>

#include <TGrid.h>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TChain.h"
#include "Riostream.h"
#include "AliLog.h"
#include "AliCascadeVertexer.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"
#include "AliExternalTrackParam.h"
#include "AliAODEvent.h"
#include "AliInputEventHandler.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include <TRandom3.h>
#include "TFile.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliOADBContainer.h"
#include "AliAnalysisTaskNucleiv2pPb.h"

ClassImp(AliAnalysisTaskNucleiv2pPb)

using std::cout;
using std::endl;

//_____________________________________________________________________________
AliAnalysisTaskNucleiv2pPb::AliAnalysisTaskNucleiv2pPb():
  AliAnalysisTaskSE(),                       //! 
  fAODevent(0),                         //! 
  fevent(0),   
  fRun(-1),
  fMultV0(0),
  fQxnmV0A(0),
  fQynmV0A(0),
  fQxnsV0A(0),
  fQynsV0A(0),
  fQxnmV0C(0),
  fQynmV0C(0),
  fQxnsV0C(0),
  fQynsV0C(0),
  fRecPass(0),
  fCenCalV0(0), //da qui --> centrality selection
  fFilterBit(4),
  fptc(1),     
  fVzmax(10),
  fNHarm(2),
  fListHist(0), 
  fHistEventMultiplicity(0), 
  fHistTrackMultiplicity(0),
  fhBB(0),
  fhBBDeu(0),
  fhTOF(0),
  fhMassTOF(0),
  EPVzAvsCentrality(0), 
  EPVzCvsCentrality(0), 
  hQVzAQVzCvsCentrality(0),
  hQxVzAvsCentrality(0),
  hQyVzAvsCentrality(0),
  hQxVzCvsCentrality(0),
  hQyVzCvsCentrality(0),
  eventtype(-999),
  ftree(0),           
  tCentrality(0),     
  tType(0),  
  tHasTOF(0),    
  tpT(0),  
  tMassTOF(0),
  tuqV0A(0),
  tuqV0C(0),
  tCharge(0),
  tCosdeltaphiV0A(0),
  tCosdeltaphiV0C(0),
  timpactXY(0),
  timpactZ(0),
  tpull(0),
  tphi(0),
  fPIDResponse(0)
{
  cout<<"Dummy constructor"<<endl;
}

//______________________________________________________________________________
AliAnalysisTaskNucleiv2pPb::AliAnalysisTaskNucleiv2pPb(const char *name):
    AliAnalysisTaskSE(name),                   //! 
    fAODevent(0),                         //! 
    fevent(0),   
    fRun(-1),
    fMultV0(0),
    fQxnmV0A(0),
    fQynmV0A(0),
    fQxnsV0A(0),
    fQynsV0A(0),
    fQxnmV0C(0),
    fQynmV0C(0),
    fQxnsV0C(0),
    fQynsV0C(0),
    fRecPass(0),
    fCenCalV0(0), //da qui
    fFilterBit(4),
    fptc(1),     
    fVzmax(10),
    fNHarm(2),
    fListHist(0), 
    fHistEventMultiplicity(0), 
    fHistTrackMultiplicity(0),
    fhBB(0),
    fhBBDeu(0),
    fhTOF(0),
    fhMassTOF(0),
    EPVzAvsCentrality(0), 
    EPVzCvsCentrality(0), 
    hQVzAQVzCvsCentrality(0),
    hQxVzAvsCentrality(0),
    hQyVzAvsCentrality(0),
    hQxVzCvsCentrality(0),
    hQyVzCvsCentrality(0),
    eventtype(-999),
    ftree(0),           
    tCentrality(0),     
    tType(0),  
    tHasTOF(0),    
    tpT(0),  
    tMassTOF(0),
    tuqV0A(0),
    tuqV0C(0),
    tCharge(0),
    tCosdeltaphiV0A(0),
    tCosdeltaphiV0C(0),
    timpactXY(0),
    timpactZ(0),
    tpull(0),
    tphi(0),
    fPIDResponse(0)
{
  
  //
  cout<<"Real constructor"<<endl;

  //  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class()); 
}

//________________________________________________________________________
Float_t AliAnalysisTaskNucleiv2pPb::GetPhi0Pi(Float_t phi){
  // Sets the phi angle in the range 0-pi
  Float_t result=phi;
  while(result<0){
    result=result+TMath::Pi();
  }
  while(result>TMath::Pi()){
    result=result-TMath::Pi();
  }
   return result;
}


//_____________________________________________________________________________
AliAnalysisTaskNucleiv2pPb::~AliAnalysisTaskNucleiv2pPb()
{

  
}

//______________________________________________________________________________
void AliAnalysisTaskNucleiv2pPb::UserCreateOutputObjects()
{ 
  //-------------------------------------------------------
  fListHist = new TList();
  fListHist->SetOwner();  // IMPORTANT!
  
  if(! fHistEventMultiplicity ){

    fHistEventMultiplicity   = new TH1F( "fHistEventMultiplicity" , "Nb of Events" , 12 , 0.5,12.5);
    
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(1,"All Events");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(2,"Events w/PV");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(3,"Events w/good PV");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(4,"Events wo pileup");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(5,"Events w/|Vz|<10cm");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(6,"kINT7");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(7,"HM V0");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(8,"HM SPD");
    
    fListHist->Add(fHistEventMultiplicity);
  }

  if(! fHistTrackMultiplicity ){
    fHistTrackMultiplicity  = new TH2F( "fHistTrackMultiplicity", "Nb of Tracks MB Events |Vz| < 10", 250,0, 5000,105,0,105);
    fHistTrackMultiplicity->GetXaxis()->SetTitle("Number of tracks");
    fHistTrackMultiplicity->GetYaxis()->SetTitle("Percentile");
    fListHist->Add(fHistTrackMultiplicity);
  } 
 
  if(! fhBB ){
    fhBB = new TH2F( "fhBB" , "BetheBlochTPC" , 240,-10,10,250,0,1000);
    fListHist->Add(fhBB);
  }
  
  if(! fhBBDeu ){
    fhBBDeu = new TH2F( "fhBBDeu" , "BetheBlochTPC - Deuteron" , 240,-10,10,250,0,1000);
    fListHist->Add(fhBBDeu);
  }
 
  if(! fhTOF ){
    fhTOF = new TH2F( "fhTOF" , "Scatter Plot TOF" , 240,-10,10,500,0,1.2);
    fListHist->Add(fhTOF);
  }
  if(! fhMassTOF){
    fhMassTOF=new TH1F ("fhMassTOF","Particle Mass - TOF", 100,0 ,10);
    fListHist->Add(fhMassTOF);
  }
  
  EPVzAvsCentrality  = new TH2D("EPVzAvsCentrality" , "EPVzAvsCentrality" , 80,-TMath::Pi(),TMath::Pi(), 105,0,105);
  EPVzCvsCentrality  = new TH2D("EPVzCvsCentrality" , "EPVzCvsCentrality" , 80,-TMath::Pi(),TMath::Pi(), 105,0,105);

  fListHist->Add(EPVzAvsCentrality);
  fListHist->Add(EPVzCvsCentrality);
  

  if(fNHarm < 3)
    hQVzAQVzCvsCentrality = new TH2F("hQVzAQVzCvsCentrality","hQVzAQVzCvsCentrality",1000,-100,100,105,0,105);
  else
    hQVzAQVzCvsCentrality = new TH2F("hQVzAQVzCvsCentrality","hQVzAQVzCvsCentrality",5000,-1000,1000,105,0,105);
  fListHist->Add(hQVzAQVzCvsCentrality);

  if(fNHarm < 3){
    hQxVzAvsCentrality = new TH2F("hQxVzAvsCentrality","hQxVzAvsCentrality",100,-20,20,105,0,105);
    hQyVzAvsCentrality = new TH2F("hQyVzAvsCentrality","hQyVzAvsCentrality",100,-20,20,105,0,105);
    hQxVzCvsCentrality = new TH2F("hQxVzCvsCentrality","hQxVzCvsCentrality",100,-20,20,105,0,105);
    hQyVzCvsCentrality = new TH2F("hQyVzCvsCentrality","hQyVzCvsCentrality",100,-20,20,105,0,105);
  }
  
  else{
    hQxVzAvsCentrality = new TH2F("hQxVzAvsCentrality","hQxVzAvsCentrality",2000,-500,500,105,0,105);
    hQyVzAvsCentrality = new TH2F("hQyVzAvsCentrality","hQyVzAvsCentrality",2000,-500,500,105,0,105);
    hQxVzCvsCentrality = new TH2F("hQxVzCvsCentrality","hQxVzCvsCentrality",2000,-500,500,105,0,105);
    hQyVzCvsCentrality = new TH2F("hQyVzCvsCentrality","hQyVzCvsCentrality",2000,-500,500,105,0,105);
  }

  fListHist->Add(hQxVzAvsCentrality);
  fListHist->Add(hQyVzAvsCentrality);
  fListHist->Add(hQxVzCvsCentrality);
  fListHist->Add(hQyVzCvsCentrality);
 
  if(!ftree){
   
    ftree = new TTree("ftree","ftree");
 
    ftree->Branch("tCentrality"      ,&tCentrality      ,"tCentrality/D"    );
    ftree->Branch("tType"            ,&tType            ,"tType/D"          );
    ftree->Branch("tHasTOF"          ,&tHasTOF          ,"tHasTOF/D"        );
    ftree->Branch("tpT"              ,&tpT              ,"tpT/D"            );
    ftree->Branch("tMassTOF"         ,&tMassTOF         ,"tMassTOF/D"       );
    ftree->Branch("tuqV0A"           ,&tuqV0A           ,"tuqV0A/D"         );
    ftree->Branch("tuqV0C"           ,&tuqV0C           ,"tuqV0C/D"         );
    ftree->Branch("tCharge"          ,&tCharge          ,"tCharge/D"        );
    ftree->Branch("tCosdeltaphiV0A"  ,&tCosdeltaphiV0A  ,"tCosdeltaphiV0A/D");
    ftree->Branch("tCosdeltaphiV0C"  ,&tCosdeltaphiV0C  ,"tCosdeltaphiV0C/D");
    ftree->Branch("timpactXY"        ,&timpactXY        ,"timpactXY/D"      );
    ftree->Branch("timpactZ"         ,&timpactZ         ,"timpactZ/D"       );
    ftree->Branch("tpull"            ,&tpull            ,"tpull/D"          );
    ftree->Branch("tphi"             ,&tphi             ,"tphi/D"           );

  }


  PostData(1,  fListHist);
  PostData(2,  ftree);  

}

//______________________________________________________________________________
void AliAnalysisTaskNucleiv2pPb::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event


  Info("AliAnalysisTaskNucleiv2pPb","Starting UserExec");  
  fHistEventMultiplicity->Fill(1);
  AliVEvent *event = InputEvent();

  fAODevent = dynamic_cast<AliAODEvent*>(event);

  if (!fAODevent) {
    AliError("Cannot get the AOD event");
      return;
  }  
  fevent = fAODevent;

  if(!fevent || !fevent->GetHeader()){
    Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    this->Dump();
    return;
  }
      
  Int_t run = fevent->GetRunNumber();
  if(run != fRun){
    // Load the calibrations run dependent
    OpenInfoCalbration(run);
    fRun = run;
  }

  //event cuts
  // 1. primary vertex selection
  const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(fevent->GetPrimaryVertex());
  if(!vtx || vtx->GetNContributors() < 1) return;
  fHistEventMultiplicity->Fill(2);
 
  // 2. SPD vertex selection
  const AliAODVertex* vtxSPD = dynamic_cast<const AliAODVertex*>(fevent->GetPrimaryVertexSPD());

  Double_t dMaxResol = 0.25; // suggested from DPG
  Double_t cov[6] = {0};
  vtxSPD->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if ( vtxSPD->IsFromVertexerZ() && (zRes > dMaxResol)) return;
  fHistEventMultiplicity->Fill(3);
  
  // 3. pileup rejection from multivertexer
  AliAnalysisUtils utils;
  utils.SetMinPlpContribMV(5);
  utils.SetMaxPlpChi2MV(5);
  utils.SetMinWDistMV(15);
  utils.SetCheckPlpFromDifferentBCMV(kFALSE);
  Bool_t isPileupFromMV = utils.IsPileUpMV(fevent);

  if(isPileupFromMV)return;
  fHistEventMultiplicity->Fill(4);
  
  // 4. cutting on PV z-distance
  const Double_t aodVtxZ = vtx->GetZ();
  if( TMath::Abs(aodVtxZ) >  fVzmax)
    return;
  fHistEventMultiplicity->Fill(5);
 
  // 5. Physics selection (trigger)
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  UInt_t fSelectMask = inputHandler->IsEventSelected();

  Bool_t isTriggerSelected = kFALSE;
 
  Bool_t isSelectedINT7  = fSelectMask& AliVEvent::kINT7;
  Bool_t isSelectedHMV0  = fSelectMask& AliVEvent::kHighMultV0;
  Bool_t isSelectedHMSPD = fSelectMask& AliVEvent::kHighMultSPD;
  
  //Int_t eventtype = -999;
  if(isSelectedINT7){
    eventtype = 1;
    fHistEventMultiplicity->Fill(6);
  }
  if(isSelectedHMV0){
    eventtype = 2;
    fHistEventMultiplicity->Fill(7);
  }
  if(isSelectedHMSPD){
    eventtype = 3;
    fHistEventMultiplicity->Fill(8);
  }

  if(eventtype!=1 && eventtype!=2 && eventtype!=3)return;
  
  // get the PID response
  fPIDResponse=inputHandler->GetPIDResponse(); 
 

  //Analysis

  Analyze(fevent);
  
    
}

//________________________________________________________________________
void AliAnalysisTaskNucleiv2pPb::Analyze(AliVEvent* aod)
{  
  
    //Centrality
    Float_t v0Centr    = -100.;
    Float_t v0aCentr   = -100.;
    Float_t v0EqCentr  = -100.;
    Float_t v0aEqCentr = -100.;
    Float_t cl1Centr   = -100.;
    Float_t cl0Centr   = -100.;
    
    AliMultSelection* MultSelection = 0x0;
    MultSelection = (AliMultSelection*) aod->FindListObject("MultSelection");
    if( !MultSelection) {
      AliWarning("AliMultSelection object not found!");
      return;
    } else {
      v0Centr = MultSelection->GetMultiplicityPercentile("V0M");
      v0aCentr = MultSelection->GetMultiplicityPercentile("V0A");
      v0EqCentr = MultSelection->GetMultiplicityPercentile("V0MEq");
      v0aEqCentr = MultSelection->GetMultiplicityPercentile("V0AEq");
      cl1Centr = MultSelection->GetMultiplicityPercentile("CL1");
      cl0Centr = MultSelection->GetMultiplicityPercentile("CL0");
    }
    
    
    //V0 info
    Double_t Qxan = 0, Qyan = 0;
    Double_t Qxcn = 0, Qycn = 0;
    Double_t sumMa = 0, sumMc = 0;
    
    AliAODVZERO* aodV0 = (AliAODVZERO*)aod->GetVZEROData();
    
    for (Int_t iV0 = 0; iV0 < 64; iV0++) {
        
      Double_t phiV0 = TMath::PiOver4()*(0.5 + iV0 % 8);
      
      Float_t multv0 = aodV0->GetMultiplicity(iV0);
      
      if (iV0 < 32){
	
	Double_t multCorC = -10;
        
	if (iV0 < 8)
	  multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(1);
	else if (iV0 >= 8 && iV0 < 16)
	  multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(9);
	else if (iV0 >= 16 && iV0 < 24)
	  multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(17);
	else if (iV0 >= 24 && iV0 < 32)
	  multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(25);
	
        
	if (multCorC < 0){
	  cout<<"Problem with multiplicity in V0C"<<endl;
	  continue;
	}
        
	Qxcn += TMath::Cos(fNHarm*phiV0) * multCorC;
	Qycn += TMath::Sin(fNHarm*phiV0) * multCorC;
        
	sumMc = sumMc + multCorC;
        
      } else {
	
	Double_t multCorA = -10;
        
	if (iV0 >= 32 && iV0 < 40)
	  multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(33);
	else if (iV0 >= 40 && iV0 < 48)
	  multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(41);
	else if (iV0 >= 48 && iV0 < 56)
	  multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(49);
	else if (iV0 >= 56 && iV0 < 64)
	  multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(57);
	
	if (multCorA < 0){
	  cout<<"Problem with multiplicity in V0A"<<endl;
	  continue;
	}
        
	Qxan += TMath::Cos(fNHarm*phiV0) * multCorA;
	Qyan += TMath::Sin(fNHarm*phiV0) * multCorA;
        
	sumMa = sumMa + multCorA;
      }
      
    }
    
    
    if (sumMa < 0 || sumMc < 0)
      return;
    
    Int_t iCen = -10;
    if (fCenCalV0 == 1)
      iCen = Int_t(v0Centr);
    else if (fCenCalV0 == 2)
      iCen = Int_t(v0aEqCentr);
    else if (fCenCalV0 == 3)
      iCen = Int_t(cl1Centr);
    else
      iCen = Int_t(v0aCentr); //default
    
    if (iCen < 0){
      cout<<"Problem with mult: Check the value!!!"<<endl;
    }
    
    Double_t QxanCor = Qxan;
    Double_t QyanCor = (Qyan - fQynmV0A->GetBinContent(iCen+1))/fQynsV0A->GetBinContent(iCen+1);
    Double_t QxcnCor = Qxcn;
    Double_t QycnCor = (Qycn - fQynmV0C->GetBinContent(iCen+1))/fQynsV0C->GetBinContent(iCen+1);
    
    if (fNHarm != 4.){
      QxanCor = (Qxan - fQxnmV0A->GetBinContent(iCen+1))/fQxnsV0A->GetBinContent(iCen+1);
      QxcnCor = (Qxcn - fQxnmV0C->GetBinContent(iCen+1))/fQxnsV0C->GetBinContent(iCen+1);
    }
   
    //not 100% sure this makes sense...
    Double_t evPlAngV0A = TMath::ATan2(QyanCor, QxanCor)/fNHarm;
    Double_t evPlAngV0C = TMath::ATan2(QycnCor, QxcnCor)/fNHarm;

    EPVzAvsCentrality  ->Fill(evPlAngV0A  , iCen); 
    EPVzCvsCentrality  ->Fill(evPlAngV0C  , iCen); 

    //Scalar Product
  
    Double_t  QV0AQV0C = QxanCor *  QxcnCor+ QyanCor*QycnCor;
    hQVzAQVzCvsCentrality->Fill(QV0AQV0C,iCen);
    
    //NUA correction
 
    hQxVzAvsCentrality->Fill(QxanCor,iCen);
    hQyVzAvsCentrality->Fill(QyanCor,iCen);
    hQxVzCvsCentrality->Fill(QxcnCor,iCen);
    hQyVzCvsCentrality->Fill(QycnCor,iCen);

    //----------------------------------------------------
    // from here my analysis starts
    
    Float_t  impactXY=-999., impactZ=-999.;
    Double_t TPCSignal=0.;
    
    ULong_t  status=0;
    
    Double_t pmax  = 10.;
    Double_t ptmax = 6.2;
      
    Double_t ptcExp  = -999;
    Double_t pullTPC = -999;
    Double_t expbeta = -999;
    Double_t pullTOF = -999;

    Float_t deltaphiV0A = -3;
    Float_t deltaphiV0C = -3;

    Double_t massd   = 1.875612859;
    Double_t masst   = 2.808939;
    Double_t mass3he = 2.80892;
    
    Float_t  uqV0A = -999;
    Float_t  uqV0C = -999; 

    Int_t TrackNumber = fevent->GetNumberOfTracks();
    fHistTrackMultiplicity->Fill(TrackNumber,iCen); //tracce per evento

    for (Int_t j=0; j<TrackNumber; j++) { //loop on tracks
  
      AliVTrack* atrack = (AliVTrack*) fevent->GetTrack(j);
      if (!atrack)
	continue;
    
      Bool_t trkFlag = ((AliAODTrack*)atrack)->TestFilterBit(fFilterBit);
      
      if(!trkFlag)continue;
      
        status  = (ULong_t)atrack->GetStatus();
    
	Bool_t hasTOFout  = status&AliVTrack::kTOFout; 
	Bool_t hasTOF     = kFALSE;
	if (hasTOFout) hasTOF = kTRUE;
	Float_t length = atrack->GetIntegratedLength(); 
	if (length < 350.) hasTOF = kFALSE;
	
	TPCSignal=atrack->GetTPCsignal(); 
	
	if(TPCSignal<10)continue;
	if(TPCSignal>1000)continue;
	
	Double_t ptot = atrack->GetTPCmomentum(); // momentum for dEdx determination
	Double_t pt  = atrack->Pt();
	
	if(ptot<0.60)continue;
	if(pt<0.60)continue;
	
	fhBB->Fill(ptot*atrack->Charge(),TPCSignal);

	Double_t d[2], covd[3];
	AliAODTrack* track_clone=(AliAODTrack*)atrack->Clone("track_clone"); // need to clone because PropagateToDCA updates the track parameters
	Bool_t isDCA = track_clone->PropagateToDCA(fevent->GetPrimaryVertex(),fevent->GetMagneticField(),9999.,d,covd);
	delete track_clone;
	if(!isDCA)d[0]=-999.;
	impactXY = d[0];
	impactZ  = d[1];

	if(fptc==1)
	  pullTPC  = TMath::Abs((fPIDResponse->NumberOfSigmasTPC(atrack,(AliPID::EParticleType)5)));;
	if(fptc==2)
	  pullTPC  = TMath::Abs((fPIDResponse->NumberOfSigmasTPC(atrack,(AliPID::EParticleType)6)));;
	if(fptc==3)
	  pullTPC  = TMath::Abs((fPIDResponse->NumberOfSigmasTPC(atrack,(AliPID::EParticleType)7)));;
	
	Double_t p    = atrack->P();
	Double_t tof  = atrack->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(p);
	Double_t tPhi = atrack->Phi();
	
	Float_t  beta = 0;
	Float_t  gamma = 0;
	Float_t  mass  = -99;
	
	if(fptc==1)
	  expbeta = TMath::Sqrt(1-((massd*massd)/(ptot*ptot+massd*massd))); 
	if(fptc==2)
	  expbeta = TMath::Sqrt(1-((masst*masst)/(ptot*ptot+masst*masst))); 
	if(fptc==3)
	  expbeta = TMath::Sqrt(1-((mass3he*mass3he)/(ptot*ptot+mass3he*mass3he))); 
        
	if(fptc==3)
	  pt = 2*pt;

	if(TMath::Abs(ptot) < pmax  && TMath::Abs(pt) < ptmax && TMath::Abs(pullTPC) <= 3){
	 

	  if (hasTOF) {
	    beta = length / (2.99792457999999984e-02 * tof);
	    gamma = 1/TMath::Sqrt(1 - beta*beta);
	    mass = ptot/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
	    
	    if(fptc==1){
	      if(TMath::Abs(mass) > 2.65)continue;
	      if(TMath::Abs(mass) < 1.05)continue;
	    }
	    if(fptc==2){
	      if(TMath::Abs(mass) > 5.0)continue;
	      if(TMath::Abs(mass) < 1.8 )continue;
	    }
	    if(fptc==3){
	      if(TMath::Abs(mass) > 5.0)continue;
	      if(TMath::Abs(mass) < 1.8)continue;
	    }
	    fhMassTOF->Fill(mass);
	    
	    fhTOF->Fill(ptot*atrack->Charge(),beta);

	  } //has tof loop

	  fhBBDeu->Fill(ptot*atrack->Charge(),TPCSignal);
      	    
	  deltaphiV0A=TMath::Cos(fNHarm*GetPhi0Pi(tPhi-evPlAngV0A));
	  deltaphiV0C=TMath::Cos(fNHarm*GetPhi0Pi(tPhi-evPlAngV0C));
      
	  // Scalar Product
      
	  uqV0A = TMath::Cos(fNHarm*tPhi)*QxanCor+TMath::Sin(fNHarm*tPhi)*QyanCor;
	  uqV0C = TMath::Cos(fNHarm*tPhi)*QxcnCor+TMath::Sin(fNHarm*tPhi)*QycnCor;
	  
	  tCentrality      = iCen;
	  tType            = eventtype;
	  tHasTOF          = hasTOF;
	  tpT              = pt;
	  tMassTOF         = mass;
	  tuqV0A           = uqV0A;
	  tuqV0C           = uqV0C;
	  tCharge          = atrack->Charge();
	  tCosdeltaphiV0A  = deltaphiV0A;
	  tCosdeltaphiV0C  = deltaphiV0C;
	  timpactXY        = impactXY;
	  timpactZ         = impactZ;
	  tpull            = pullTPC;
	  tphi             = tPhi;
	  
	  if(pt<1.5)   
	    ftree->Fill();
	  else
	    if(hasTOF==1)
	      ftree->Fill();
	  
	}//POI selection
	
    }//track loop

    PostData(1, fListHist);
    PostData(2, ftree);

}

//_____________________________________________________________________________
void AliAnalysisTaskNucleiv2pPb::OpenInfoCalbration(Int_t run)
{
  
  TFile* foadb = 0;
  if (!gGrid) TGrid::Connect("alien");

  if (fRecPass == 0)
    foadb = TFile::Open("alien:///alice/cern.ch/user/l/lramona/CalibpPb2016/calibV0Fast.root");
  else if (fRecPass == 1)
    foadb = TFile::Open("alien:///alice/cern.ch/user/l/lramona/CalibpPb2016/calibV0NoSDD.root");
  else if (fRecPass == 2)
    foadb = TFile::Open("alien:///alice/cern.ch/user/l/lramona/CalibpPb2016/calibV0SDD.root");
    
  if(!foadb){
    printf("OADB V0 calibration file cannot be opened\n");
    return;
  }
    
 
    
  AliOADBContainer* cont = (AliOADBContainer*) foadb->Get("hMultV0BefCorPfpx");
  if(!cont){
    printf("OADB object hMultV0BefCorr is not available in the file\n");
    return;
  }
  if(!(cont->GetObject(run))){
    printf("OADB object hMultV0BefCorPfpx is not available for run %i\n", run);
    return;
  }
  fMultV0 = ((TH1D*) cont->GetObject(run));
    
    
    
  AliOADBContainer* contQxnam = 0;
  if (fCenCalV0 == 0){
    if (fNHarm == 2.)
      contQxnam = (AliOADBContainer*) foadb->Get("fqxa2mV0A");
    else if (fNHarm == 3.)
      contQxnam = (AliOADBContainer*) foadb->Get("fqxa3mV0A");
  } else if (fCenCalV0 == 1){
    if (fNHarm == 2.)
      contQxnam = (AliOADBContainer*) foadb->Get("fqxa2mV0");
    else if (fNHarm == 3.)
      contQxnam = (AliOADBContainer*) foadb->Get("fqxa3mV0");
  } else if (fCenCalV0 == 2){
    if (fNHarm == 2.)
      contQxnam = (AliOADBContainer*) foadb->Get("fqxa2mV0AEq");
    else if (fNHarm == 3.)
      contQxnam = (AliOADBContainer*) foadb->Get("fqxa3mV0AEq");
  } else if (fCenCalV0 == 3){
    if (fNHarm == 2.)
      contQxnam = (AliOADBContainer*) foadb->Get("fqxa2mCL1");
    else if (fNHarm == 3.)
      contQxnam = (AliOADBContainer*) foadb->Get("fqxa3mCL1");
  }
  if(!contQxnam){
    printf("OADB object fqxanm is not available in the file\n");
    return;
  }
  if(!(contQxnam->GetObject(run))){
    printf("OADB object fqxanm is not available for run %i\n", run);
    return;
  }
  fQxnmV0A = ((TH1D*) contQxnam->GetObject(run));
    
    
    
  AliOADBContainer* contQynam = 0;
  if (fCenCalV0 == 0){
    if (fNHarm == 2.)
      contQynam = (AliOADBContainer*) foadb->Get("fqya2mV0A");
    else if (fNHarm == 3.)
      contQynam = (AliOADBContainer*) foadb->Get("fqya3mV0A");
  } else if (fCenCalV0 == 1){
    if (fNHarm == 2.)
      contQynam = (AliOADBContainer*) foadb->Get("fqya2mV0");
    else if (fNHarm == 3.)
      contQynam = (AliOADBContainer*) foadb->Get("fqya3mV0");
  } else if (fCenCalV0 == 2){
    if (fNHarm == 2.)
      contQynam = (AliOADBContainer*) foadb->Get("fqya2mV0AEq");
    else if (fNHarm == 3.)
      contQynam = (AliOADBContainer*) foadb->Get("fqya3mV0AEq");
  } else if (fCenCalV0 == 3){
    if (fNHarm == 2.)
      contQynam = (AliOADBContainer*) foadb->Get("fqya2mCL1");
    else if (fNHarm == 3.)
      contQynam = (AliOADBContainer*) foadb->Get("fqya3mCL1");
  }
  if(!contQynam){
    printf("OADB object fqyanm is not available in the file\n");
    return;
  }
  if(!(contQynam->GetObject(run))){
    printf("OADB object fqyanm is not available for run %i\n", run);
    return;
  }
  fQynmV0A = ((TH1D*) contQynam->GetObject(run));
    
    
    
  AliOADBContainer* contQxnas = 0;
  if (fCenCalV0 == 0){
    if (fNHarm == 2.)
      contQxnas = (AliOADBContainer*) foadb->Get("fqxa2sV0A");
    else if (fNHarm == 3.)
      contQxnas = (AliOADBContainer*) foadb->Get("fqxa3sV0A");
  } else if (fCenCalV0 == 1){
    if (fNHarm == 2.)
      contQxnas = (AliOADBContainer*) foadb->Get("fqxa2sV0");
    else if (fNHarm == 3.)
      contQxnas = (AliOADBContainer*) foadb->Get("fqxa3sV0");
  } else if (fCenCalV0 == 2){
    if (fNHarm == 2.)
      contQxnas = (AliOADBContainer*) foadb->Get("fqxa2sV0AEq");
    else if (fNHarm == 3.)
      contQxnas = (AliOADBContainer*) foadb->Get("fqxa3sV0AEq");
  } else if (fCenCalV0 == 3){
    if (fNHarm == 2.)
      contQxnas = (AliOADBContainer*) foadb->Get("fqxa2sCL1");
    else if (fNHarm == 3.)
      contQxnas = (AliOADBContainer*) foadb->Get("fqxa3sCL1");
  }
  if(!contQxnas){
    printf("OADB object fqxans is not available in the file\n");
    return;
  }
  if(!(contQxnas->GetObject(run))){
    printf("OADB object fqxans is not available for run %i\n", run);
    return;
  }
  fQxnsV0A = ((TH1D*) contQxnas->GetObject(run));
    
    
    
  AliOADBContainer* contQynas = 0;
  if (fCenCalV0 == 0){
    if (fNHarm == 2.)
      contQynas = (AliOADBContainer*) foadb->Get("fqya2sV0A");
    else if (fNHarm == 3.)
      contQynas = (AliOADBContainer*) foadb->Get("fqya3sV0A");
  } else if (fCenCalV0 == 1){
    if (fNHarm == 2.)
      contQynas = (AliOADBContainer*) foadb->Get("fqya2sV0");
    else if (fNHarm == 3.)
      contQynas = (AliOADBContainer*) foadb->Get("fqya3sV0");
  } else if (fCenCalV0 == 2){
    if (fNHarm == 2.)
      contQynas = (AliOADBContainer*) foadb->Get("fqya2sV0AEq");
    else if (fNHarm == 3.)
      contQynas = (AliOADBContainer*) foadb->Get("fqya3sV0AEq");
  } else if (fCenCalV0 == 3){
    if (fNHarm == 2.)
      contQynas = (AliOADBContainer*) foadb->Get("fqya2sCL1");
    else if (fNHarm == 3.)
      contQynas = (AliOADBContainer*) foadb->Get("fqya3sCL1");
  }
  if(!contQynas){
    printf("OADB object fqyans is not available in the file\n");
    return;
  }
  if(!(contQynas->GetObject(run))){
    printf("OADB object fqyans is not available for run %i\n", run);
    return;
  }
  fQynsV0A = ((TH1D*) contQynas->GetObject(run));
    
    
    
  AliOADBContainer* contQxncm = 0;
  if (fCenCalV0 == 0){
    if (fNHarm == 2.)
      contQxncm = (AliOADBContainer*) foadb->Get("fqxc2mV0A");
    else if (fNHarm == 3.)
      contQxncm = (AliOADBContainer*) foadb->Get("fqxc3mV0A");
  } else if (fCenCalV0 == 1){
    if (fNHarm == 2.)
      contQxncm = (AliOADBContainer*) foadb->Get("fqxc2mV0");
    else if (fNHarm == 3.)
      contQxncm = (AliOADBContainer*) foadb->Get("fqxc3mV0");
  } else if (fCenCalV0 == 2){
    if (fNHarm == 2.)
      contQxncm = (AliOADBContainer*) foadb->Get("fqxc2mV0AEq");
    else if (fNHarm == 3.)
      contQxncm = (AliOADBContainer*) foadb->Get("fqxc3mV0AEq");
  } else if (fCenCalV0 == 3){
    if (fNHarm == 2.)
      contQxncm = (AliOADBContainer*) foadb->Get("fqxc2mCL1");
    else if (fNHarm == 3.)
      contQxncm = (AliOADBContainer*) foadb->Get("fqxc3mCL1");
  }
  if(!contQxncm){
    printf("OADB object fqxcnm is not available in the file\n");
    return;
  }
  if(!(contQxncm->GetObject(run))){
    printf("OADB object fqxcnm is not available for run %i\n", run);
    return;
  }
  fQxnmV0C = ((TH1D*) contQxncm->GetObject(run));
    
    
    
  AliOADBContainer* contQyncm = 0;
  if (fCenCalV0 == 0){
    if (fNHarm == 2.)
      contQyncm = (AliOADBContainer*) foadb->Get("fqyc2mV0A");
    else if (fNHarm == 3.)
      contQyncm = (AliOADBContainer*) foadb->Get("fqyc3mV0A");
  } else if (fCenCalV0 == 1){
    if (fNHarm == 2.)
      contQyncm = (AliOADBContainer*) foadb->Get("fqyc2mV0");
    else if (fNHarm == 3.)
      contQyncm = (AliOADBContainer*) foadb->Get("fqyc3mV0");
  } else if (fCenCalV0 == 2){
    if (fNHarm == 2.)
      contQyncm = (AliOADBContainer*) foadb->Get("fqyc2mV0AEq");
    else if (fNHarm == 3.)
      contQyncm = (AliOADBContainer*) foadb->Get("fqyc3mV0AEq");
  } else if (fCenCalV0 == 3){
    if (fNHarm == 2.)
      contQyncm = (AliOADBContainer*) foadb->Get("fqyc2mCL1");
    else if (fNHarm == 3.)
      contQyncm = (AliOADBContainer*) foadb->Get("fqyc3mCL1");
  }
  if(!contQyncm){
    printf("OADB object fqyc2m is not available in the file\n");
    return;
  }
  if(!(contQyncm->GetObject(run))){
    printf("OADB object fqyc2m is not available for run %i\n", run);
    return;
  }
  fQynmV0C = ((TH1D*) contQyncm->GetObject(run));
    
    
    
  AliOADBContainer* contQxncs = 0;
  if (fCenCalV0 == 0){
    if (fNHarm == 2.)
      contQxncs = (AliOADBContainer*) foadb->Get("fqxc2sV0A");
    else if (fNHarm == 3.)
      contQxncs = (AliOADBContainer*) foadb->Get("fqxc3sV0A");
  } else if (fCenCalV0 == 1){
    if (fNHarm == 2.)
      contQxncs = (AliOADBContainer*) foadb->Get("fqxc2sV0");
    else if (fNHarm == 3.)
      contQxncs = (AliOADBContainer*) foadb->Get("fqxc3sV0");
  } else if (fCenCalV0 == 2){
    if (fNHarm == 2.)
      contQxncs = (AliOADBContainer*) foadb->Get("fqxc2sV0AEq");
    else if (fNHarm == 3.)
      contQxncs = (AliOADBContainer*) foadb->Get("fqxc3sV0AEq");
  } else if (fCenCalV0 == 3){
    if (fNHarm == 2.)
      contQxncs = (AliOADBContainer*) foadb->Get("fqxc2sCL1");
    else if (fNHarm == 3.)
      contQxncs = (AliOADBContainer*) foadb->Get("fqxc3sCL1");
  }
  if(!contQxncs){
    printf("OADB object fqxc2s is not available in the file\n");
    return;
  }
  if(!(contQxncs->GetObject(run))){
    printf("OADB object fqxc2s is not available for run %i\n", run);
    return;
  }
  fQxnsV0C = ((TH1D*) contQxncs->GetObject(run));
    
    
    
  AliOADBContainer* contQyncs = 0;
  if (fCenCalV0 == 0){
    if (fNHarm == 2.)
      contQyncs = (AliOADBContainer*) foadb->Get("fqyc2sV0A");
    else if (fNHarm == 3.)
      contQyncs = (AliOADBContainer*) foadb->Get("fqyc3sV0A");
  } else if (fCenCalV0 == 1){
    if (fNHarm == 2.)
      contQyncs = (AliOADBContainer*) foadb->Get("fqyc2sV0");
    else if (fNHarm == 3.)
      contQyncs = (AliOADBContainer*) foadb->Get("fqyc3sV0");
  } else if (fCenCalV0 == 2){
    if (fNHarm == 2.)
      contQyncs = (AliOADBContainer*) foadb->Get("fqyc2sV0AEq");
    else if (fNHarm == 3.)
      contQyncs = (AliOADBContainer*) foadb->Get("fqyc3sV0AEq");
  } else if (fCenCalV0 == 3){
    if (fNHarm == 2.)
      contQyncs = (AliOADBContainer*) foadb->Get("fqyc2sCL1");
    else if (fNHarm == 3.)
      contQyncs = (AliOADBContainer*) foadb->Get("fqyc3sCL1");
  }
  if(!contQyncs){
    printf("OADB object fqycnm is not available in the file\n");
    return;
  }
  if(!(contQyncs->GetObject(run))){
    printf("OADB object fqycns is not available for run %i\n", run);
    return;
  }
  fQynsV0C = ((TH1D*) contQyncs->GetObject(run));
    
}

//_____________________________________________________________________________
void AliAnalysisTaskNucleiv2pPb::Terminate(Option_t *)
{ 
  // Terminate loop
  Printf("Terminate()");
}
