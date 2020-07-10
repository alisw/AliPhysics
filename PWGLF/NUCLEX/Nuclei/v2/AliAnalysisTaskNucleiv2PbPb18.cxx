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
//                 AliAnalysisTaskNucleiv2PbPb18 class
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
#include "AliEventCuts.h"
#include "AliAnalysisTaskNucleiv2PbPb18.h"

ClassImp(AliAnalysisTaskNucleiv2PbPb18)

using std::cout;
using std::endl;

//_____________________________________________________________________________
AliAnalysisTaskNucleiv2PbPb18::AliAnalysisTaskNucleiv2PbPb18():
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
  fPeriod(1),
  fNsigma(3),
  fSplines(0),
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
  hQVzAQTPCvsCentrality(0),
  hQVzCQTPCvsCentrality(0),
  hQxVzAvsCentrality(0),
  hQyVzAvsCentrality(0),
  hQxVzCvsCentrality(0),
  hQyVzCvsCentrality(0),
  hCos2DeltaTPCVzAvsCentrality(0),
  hCos2DeltaTPCVzCvsCentrality(0),
  hCos2DeltaVzAVzCvsCentrality(0),
  hCos2DeltaVzATPCvsCentrality(0),
  hCos2DeltaVzCTPCvsCentrality(0),
  hCos2DeltaVzCVzAvsCentrality(0),
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
  tNpidcluster(0),
  fPIDResponse(0),
  fEventCuts(0)
{
  cout<<"Dummy constructor"<<endl;
}

//______________________________________________________________________________
AliAnalysisTaskNucleiv2PbPb18::AliAnalysisTaskNucleiv2PbPb18(const char *name):
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
    fPeriod(1),
    fNsigma(3),
    fSplines(0),
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
    hQVzAQTPCvsCentrality(0),
    hQVzCQTPCvsCentrality(0),
    hQxVzAvsCentrality(0),
    hQyVzAvsCentrality(0),
    hQxVzCvsCentrality(0),
    hQyVzCvsCentrality(0),
    hCos2DeltaTPCVzAvsCentrality(0),
    hCos2DeltaTPCVzCvsCentrality(0),
    hCos2DeltaVzAVzCvsCentrality(0),
    hCos2DeltaVzATPCvsCentrality(0),
    hCos2DeltaVzCTPCvsCentrality(0),
    hCos2DeltaVzCVzAvsCentrality(0),
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
    tNpidcluster(0),
    fPIDResponse(0),
    fEventCuts(0)
{
  
  //
  cout<<"Real constructor"<<endl;

  //  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class()); 
}

//________________________________________________________________________
Float_t AliAnalysisTaskNucleiv2PbPb18::GetPhi0Pi(Float_t phi){
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
//_________________________________________________________________________
// //Parametrization of Bethe-Block
// Double_t paramDandTdata[5] = { 6.70549, 6.11866, 8.86205e-15, 2.34059, 1.07029};
// Double_t paramHe3data[5]   = { 1.74962, 27.4992, 4.00313e-15, 2.48485, 8.31768};
// Double_t paramDandTmc[5]   = { 20.1533, 2.58127, 0.00114169,  2.03730, 0.502123};
// Double_t paramHe3mc[5]     = { 20.1533, 2.58127, 0.00114169,  2.03730, 0.502123};
//________________________________________________________________________
Float_t AliAnalysisTaskNucleiv2PbPb18::nSigmaTPC3He (AliAODTrack *track)  {
  Double_t fParamHe3[5]   = { 1.74962, 27.4992, 4.00313e-15, 2.48485, 8.31768};
  //Variables
  Double_t p = track->GetTPCmomentum();
  Double_t mass = AliPID::ParticleMass (AliPID::kHe3);
  Double_t dEdx_au = track->GetTPCsignal();
  //Expected dE/dx for 3He
  Float_t hel3Exp = 4.0*AliExternalTrackParam::BetheBlochAleph(2.0*p/mass,fParamHe3[0],fParamHe3[1],fParamHe3[2],fParamHe3[3],fParamHe3[4]);
  Double_t sigma = 0.07;//dE/dx Resolution for 3He (7%)
  Double_t nSigmaHe3  = (dEdx_au - hel3Exp)/(sigma*hel3Exp);
  return nSigmaHe3;
}
//---------------------------------------------------
Float_t AliAnalysisTaskNucleiv2PbPb18::nSigmaTPCdandt (AliAODTrack *track)  {
  Double_t paramDandTdata[5] = { 6.70549, 6.11866, 8.86205e-15, 2.34059, 1.07029};
  //Variables
  Double_t p = track->GetTPCmomentum();
  Double_t mass = 0;
  if(fptc == 1)mass = AliPID::ParticleMass (AliPID::kDeuteron);
  else if(fptc == 2) mass = AliPID::ParticleMass (AliPID::kTriton);
  Double_t dEdx_au = track->GetTPCsignal();
  //Expected dE/dx for d and t 
  Float_t dandtExp = 1.0*AliExternalTrackParam::BetheBlochAleph(2.0*p/mass,paramDandTdata[0],paramDandTdata[1],paramDandTdata[2],paramDandTdata[3],paramDandTdata[4]);
  Double_t sigma = 0.07;//dE/dx Resolution for 3He (7%)
  Double_t nSigma  = (dEdx_au - dandtExp)/(sigma*dandtExp);
  return nSigma;
}
//_____________________________________________________________________________
AliAnalysisTaskNucleiv2PbPb18::~AliAnalysisTaskNucleiv2PbPb18()
{

  
}

//______________________________________________________________________________
void AliAnalysisTaskNucleiv2PbPb18::UserCreateOutputObjects()
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
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(9,"kCentral");
    fHistEventMultiplicity->GetXaxis()->SetBinLabel(10,"kSemiCentral");

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
  

  // if(fNHarm < 3){
  //   hQVzAQVzCvsCentrality = new TH2F("hQVzAQVzCvsCentrality","hQVzAQVzCvsCentrality",1000,-100,100,105,0,105);
  //   hQVzAQTPCvsCentrality = new TH2F("hQVzAQTPCvsCentrality","hQVzAQTPCvsCentrality",1000,-100,100,105,0,105);
  //   hQVzCQTPCvsCentrality = new TH2F("hQVzCQTPCvsCentrality","hQVzCQTPCvsCentrality",1000,-100,100,105,0,105);
  // }
  // else{
  hQVzAQVzCvsCentrality = new TH2F("hQVzAQVzCvsCentrality","hQVzAQVzCvsCentrality",5000,-5000,5000,105,0,105);
  hQVzAQTPCvsCentrality = new TH2F("hQVzAQTPCvsCentrality","hQVzAQTPCvsCentrality",5000,-5000,5000,105,0,105);
  hQVzCQTPCvsCentrality = new TH2F("hQVzCQTPCvsCentrality","hQVzCQTPCvsCentrality",5000,-5000,5000,105,0,105);
  //}
  fListHist->Add(hQVzAQVzCvsCentrality);
  fListHist->Add(hQVzAQTPCvsCentrality);
  fListHist->Add(hQVzCQTPCvsCentrality);

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

  hCos2DeltaTPCVzAvsCentrality   = new TH2F("hCos2DeltaTPCVzAvsCentrality"  ,"hCos2DeltaTPCVzAvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaTPCVzCvsCentrality   = new TH2F("hCos2DeltaTPCVzCvsCentrality"  ,"hCos2DeltaTPCVzCvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaVzAVzCvsCentrality   = new TH2F("hCos2DeltaVzAVzCvsCentrality"  ,"hCos2DeltaVzAVzCvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaVzATPCvsCentrality   = new TH2F("hCos2DeltaVzATPCvsCentrality"  ,"hCos2DeltaVzATPCvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaVzCTPCvsCentrality   = new TH2F("hCos2DeltaVzCTPCvsCentrality"  ,"hCos2DeltaVzCTPCvsCentrality"  ,100,-1.1,1.1,105,0,105);
  hCos2DeltaVzCVzAvsCentrality   = new TH2F("hCos2DeltaVzCVzAvsCentrality"  ,"hCos2DeltaVzCVzAvsCentrality"  ,100,-1.1,1.1,105,0,105);

  fListHist->Add(hCos2DeltaTPCVzAvsCentrality);
  fListHist->Add(hCos2DeltaTPCVzCvsCentrality);
  fListHist->Add(hCos2DeltaVzAVzCvsCentrality);
  fListHist->Add(hCos2DeltaVzATPCvsCentrality);
  fListHist->Add(hCos2DeltaVzCTPCvsCentrality);
  fListHist->Add(hCos2DeltaVzCVzAvsCentrality);
  
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
    ftree->Branch("tNpidcluster"     ,&tNpidcluster     ,"tNpidcluster/I"   );

  }

  fEventCuts.AddQAplotsToList(fListHist);

  PostData(1,  fListHist);
  PostData(2,  ftree);  

}

//______________________________________________________________________________
void AliAnalysisTaskNucleiv2PbPb18::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
 
  Info("AliAnalysisTaskNucleiv2PbPb18","Starting UserExec");  
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

  //cout<<"-------------------->RUN number "<<run<<endl;
  
  if(run != fRun){
    // Load the calibrations run dependent
    OpenInfoCalbration(run);
    fRun = run;
  }
  
  //  cout<<"Done info calib"<<endl;

  /// Use the event cut class to apply the required selections
  if (!fEventCuts.AcceptEvent(fevent)) {
    PostData(1, fListHist);
    PostData(2, ftree);
    return;
  }
  
  fHistEventMultiplicity->Fill(2);

  // Primary vertex cut
  Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};

  const AliVVertex* vertexmain = fevent->GetPrimaryVertex();
  if (!vertexmain){
    AliWarning("No prim. vertex in ESD... return!");
    
    PostData(1, fListHist);
    PostData(2, ftree);
    return;
  }

  //  cout<<"Open vertex"<<endl;

  
  vertexmain->GetXYZ( lBestPrimaryVtxPos );
 
  if((TMath::Abs(lBestPrimaryVtxPos[2])) > fVzmax) return;
  fHistEventMultiplicity->Fill(3);

  Bool_t isPileUpSpd=kFALSE;
  isPileUpSpd=fAODevent->IsPileupFromSPD();
  
  if(isPileUpSpd){  
    PostData(1, fListHist);
    PostData(2, ftree);
    return;
  }
  
  //  cout<<"is pileup"<<endl;

  fHistEventMultiplicity->Fill(4); // analyzed events with PV w/o pile up
  
  /*
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
  */
  
  // 5. Physics selection (trigger)
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  UInt_t fSelectMask = inputHandler->IsEventSelected();
  

  Bool_t isTriggerSelected = kFALSE;
 
  Bool_t isSelectedINT7        = fSelectMask& AliVEvent::kINT7;
  Bool_t isSelectedCentral     = fSelectMask& AliVEvent::kCentral;
  Bool_t isSelectedSemiCentral = fSelectMask& AliVEvent::kSemiCentral;
  Bool_t isSelectedHMV0        = fSelectMask& AliVEvent::kHighMultV0;
  Bool_t isSelectedHMSPD       = fSelectMask& AliVEvent::kHighMultSPD;
  
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
  if(isSelectedCentral){
    eventtype = 4;
    fHistEventMultiplicity->Fill(9);
  }
  if(isSelectedSemiCentral){
    eventtype = 5;
    fHistEventMultiplicity->Fill(10);
  }

  //  if(eventtype!=1 && eventtype!=2 && eventtype!=3)return;
  if(eventtype!=1 && eventtype!=4 && eventtype!=5)return;
  
  //  cout<<"event selected"<<endl;

  // get the PID response
  fPIDResponse=inputHandler->GetPIDResponse(); 
 

  //Analysis
  //cout<<"Start event"<<endl;
  Analyze(fevent);
  
    
}

//________________________________________________________________________
void AliAnalysisTaskNucleiv2PbPb18::Analyze(AliVEvent* aod)
{  
  //  cout<<"OPen analysize"<<endl;

  //new vertex selection
  const AliAODVertex* vtTrc = (AliAODVertex*)aod->GetPrimaryVertex();
  const AliAODVertex* vtSPD = (AliAODVertex*)aod->GetPrimaryVertexSPD();
  
  double covTrc[6], covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);
  
  double dz = vtTrc->GetZ() - vtSPD->GetZ();
  
  double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
  double errTrc = TMath::Sqrt(covTrc[5]);
  double nsigTot = dz/errTot;
  double nsigTrc = dz/errTrc;
    
  if (TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20)
    return; // bad vertexing
  
  //Centrality
  Float_t v0Centr    = -100.;
  Float_t cl1Centr   = -100.;
  Float_t cl0Centr   = -100.;
    
  AliMultSelection* MultSelection = 0x0;
  MultSelection = (AliMultSelection*) aod->FindListObject("MultSelection");
  if( !MultSelection) {
    AliWarning("AliMultSelection object not found!");
    return;
  } else {
    v0Centr  = MultSelection->GetMultiplicityPercentile("V0M");
    cl1Centr = MultSelection->GetMultiplicityPercentile("CL1");
    cl0Centr = MultSelection->GetMultiplicityPercentile("CL0");
  }
  
  //  cout<<"v0Centr: "<<v0Centr<<endl;
  
  if (v0Centr >= 90. || v0Centr < 0)
    return; 
  
  AliAODVZERO* aodV0 = (AliAODVZERO*)aod->GetVZEROData();
  Float_t multV0a = aodV0->GetMTotV0A();
  Float_t multV0c = aodV0->GetMTotV0C();
  Float_t multV0Tot = multV0a + multV0c;
  
  if (fCenCalV0 == 1)
    v0Centr = cl0Centr;
  else if (fCenCalV0 == 2)
    v0Centr = cl1Centr;

  /*
    Short_t centrCode = -10;
    if ((v0Centr >= 0) && (v0Centr < 5.))
    centrCode = 0;
    else if ((v0Centr >= 5.) && (v0Centr < 10.))
    centrCode = 1;
    else if ((v0Centr >= 10.) && (v0Centr < 20.))
    centrCode = 2;
    else if ((v0Centr >= 20.) && (v0Centr < 30.))
    centrCode = 3;
    else if ((v0Centr >= 30.) && (v0Centr < 40.))
    centrCode = 4;
    else if ((v0Centr >= 40.) && (v0Centr < 50.))
    centrCode = 5;
    else if ((v0Centr >= 50.) && (v0Centr < 60.))
    centrCode = 6;
    else if ((v0Centr >= 60.) && (v0Centr < 70.))
    centrCode = 7;
    else if ((v0Centr >= 70.) && (v0Centr < 80.))
    centrCode = 8;
    else if ((v0Centr >= 80.) && (v0Centr < 90.))
    centrCode = 9;
    
    if (centrCode < 0)
    return;
  */

  Int_t iCen = Int_t(v0Centr);
  
  // Int_t iCentSPD = Int_t(cl1Centr);
  // if (iCentSPD >= 90)
  //   return;
    
  // Int_t iCen = Int_t(v0Centr);
  // if (iCen >= 90)
  //   return;
  
  fHistEventMultiplicity->Fill(11);
  
  //V0 info
  Double_t Qxan = 0, Qyan = 0;
  Double_t Qxcn = 0, Qycn = 0;
  Double_t sumMa = 0, sumMc = 0;
  
  //  cout<<"qui"<<endl;

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
      
      
      // if (fIsQ2Ana){
      //     Qxc2ese += TMath::Cos(2.*phiV0) * multCorC;
      //     Qyc2ese += TMath::Sin(2.*phiV0) * multCorC;
      // } else {
      //     Qxc3ese += TMath::Cos(3.*phiV0) * multCorC;
      //     Qyc3ese += TMath::Sin(3.*phiV0) * multCorC;
      // }
      
      
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

  //cout<<"qua"<<endl;

  if (sumMa < 0 || sumMc < 0)
    return;

  Double_t QxanCor = Qxan;
  Double_t QyanCor = (Qyan - fQynmV0A->GetBinContent(iCen+1))/fQynsV0A->GetBinContent(iCen+1);
  Double_t QxcnCor = Qxcn;
  Double_t QycnCor = (Qycn - fQynmV0C->GetBinContent(iCen+1))/fQynsV0C->GetBinContent(iCen+1);
  
  if (fNHarm != 4.){
    QxanCor = (Qxan - fQxnmV0A->GetBinContent(iCen+1))/fQxnsV0A->GetBinContent(iCen+1);
    QxcnCor = (Qxcn - fQxnmV0C->GetBinContent(iCen+1))/fQxnsV0C->GetBinContent(iCen+1);
  }

   
  // cout<<"quo"<<endl;
  
  // cout<<"iCen+1 "<<iCen+1 <<endl;
  
  
  // Double_t QxanCor = Qxan;
  // Double_t QyanCor = (Qyan - fQynmV0A->GetBinContent(iCen+1))/fQynsV0A->GetBinContent(iCen+1);
  // Double_t QxcnCor = Qxcn;
  // Double_t QycnCor = (Qycn - fQynmV0C->GetBinContent(iCen+1))/fQynsV0C->GetBinContent(iCen+1);
  
  // if (fNHarm != 4.){
  //   QxanCor = (Qxan - fQxnmV0A->GetBinContent(iCen+1))/fQxnsV0A->GetBinContent(iCen+1);
  //   QxcnCor = (Qxcn - fQxnmV0C->GetBinContent(iCen+1))/fQxnsV0C->GetBinContent(iCen+1);
  // }
  
  Double_t evPlAngV0A = TMath::ATan2(QyanCor, QxanCor)/fNHarm;
  Double_t evPlAngV0C = TMath::ATan2(QycnCor, QxcnCor)/fNHarm;
 
  EPVzAvsCentrality  ->Fill(evPlAngV0A  , iCen); 
  EPVzCvsCentrality  ->Fill(evPlAngV0C  , iCen); 
  
  const Int_t nTracks = aod->GetNumberOfTracks();
  // cout<<"TPC ev plane "<<nTracks<<endl;
  Double_t Qxtn = 0, Qytn = 0;
  
  for (Int_t it1 = 0; it1 < nTracks; it1++) {
    AliAODTrack* aodTrk1 = (AliAODTrack*)aod->GetTrack(it1);
    
    if (!aodTrk1){
      delete aodTrk1;
      continue;
    }
    
    if (aodTrk1->TestFilterBit(768) && TMath::Abs(aodTrk1->Eta()) < 0.8 && aodTrk1->GetTPCNcls() >= 70 && aodTrk1->Pt() >= 0.2 && aodTrk1->Pt() < 3.){
      
      Qxtn += TMath::Cos(fNHarm*aodTrk1->Phi());
      Qytn += TMath::Sin(fNHarm*aodTrk1->Phi());
    }
  }
  // TBC
  Double_t evPlAngTPC = TMath::ATan2(Qytn, Qxtn)/fNHarm;

  hCos2DeltaTPCVzAvsCentrality  ->Fill(TMath::Cos(fNHarm*(evPlAngTPC - evPlAngV0A)) , iCen);
  hCos2DeltaTPCVzCvsCentrality  ->Fill(TMath::Cos(fNHarm*(evPlAngTPC - evPlAngV0C)) , iCen);
  hCos2DeltaVzAVzCvsCentrality  ->Fill(TMath::Cos(fNHarm*(evPlAngV0A - evPlAngV0C)) , iCen);
  hCos2DeltaVzATPCvsCentrality  ->Fill(TMath::Cos(fNHarm*(evPlAngV0A - evPlAngTPC)) , iCen);
  hCos2DeltaVzCTPCvsCentrality  ->Fill(TMath::Cos(fNHarm*(evPlAngV0C - evPlAngTPC)) , iCen);
  hCos2DeltaVzCVzAvsCentrality  ->Fill(TMath::Cos(fNHarm*(evPlAngV0C - evPlAngV0A)) , iCen);
  //Scalar Product -- Resolutions
  
  Double_t corV0AV0Cvn = QxanCor*QxcnCor + QyanCor*QycnCor;
  Double_t corV0ATPCvn = QxanCor*Qxtn + QyanCor*Qytn;
  Double_t corV0CTPCvn = QxcnCor*Qxtn + QycnCor*Qytn;
  
  // cout<<"corV0AV0Cvn "<<corV0AV0Cvn <<endl;
  // cout<<"corV0ATPCvn "<<corV0ATPCvn <<endl;
  // cout<<"corV0CTPCvn "<<corV0CTPCvn <<endl;
  
  hQVzAQVzCvsCentrality->Fill(corV0AV0Cvn,iCen);
  hQVzAQTPCvsCentrality->Fill(corV0ATPCvn,iCen);
  hQVzCQTPCvsCentrality->Fill(corV0CTPCvn,iCen);
  
  // cout<<"corV0AV0Cvn "<<corV0AV0Cvn<<endl;

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
  //Double_t ptmax = 6.2;
  Double_t ptmax = 8.2;
  if(fptc == 2)
    ptmax = 3.5;
  Double_t ptcExp  = -999;
  Double_t pullTPC = -999;
  //  Double_t expbeta = -999;
  //  Double_t pullTOF = -999;

  Float_t deltaphiV0A = -3;
  Float_t deltaphiV0C = -3;

  Double_t massd   = 1.875612859;
  Double_t masst   = 2.808938914;
  Double_t mass3he = 2.808409385;
    
  Float_t  uqV0A = -999;
  Float_t  uqV0C = -999; 

  Int_t TrackNumber = aod->GetNumberOfTracks();
  
  // cout<<"--------------------------------------: "<<fevent->GetNumberOfTracks()<<endl;
  
  fHistTrackMultiplicity->Fill(TrackNumber,iCen); //tracce per evento

  for (Int_t j=0; j<TrackNumber; j++) { //loop on tracks
  
    AliVTrack* atrack = (AliVTrack*) aod->GetTrack(j);
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
    Bool_t isDCA = track_clone->PropagateToDCA(aod->GetPrimaryVertex(),aod->GetMagneticField(),9999.,d,covd);
    delete track_clone;
    if(!isDCA)d[0]=-999.;
    impactXY = d[0];
    impactZ  = d[1];

    if(fptc==1){
      if(fSplines == kFALSE)
	pullTPC  = (fPIDResponse->NumberOfSigmasTPC(atrack,(AliPID::EParticleType)5));
      else
	pullTPC = nSigmaTPCdandt((AliAODTrack*)atrack);
    }
    if(fptc==2){
      if(fSplines == kFALSE)
	pullTPC  = (fPIDResponse->NumberOfSigmasTPC(atrack,(AliPID::EParticleType)6));
      else
	pullTPC = nSigmaTPCdandt((AliAODTrack*)atrack);
    }
    if(fptc==3){
      //pullTPC  = (fPIDResponse->NumberOfSigmasTPC(atrack,(AliPID::EParticleType)7));
      pullTPC  = nSigmaTPC3He((AliAODTrack*)atrack);
    }
    
    
    Double_t p    = atrack->P();
    Double_t tof  = atrack->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(p);
    Double_t tPhi = atrack->Phi();

    
    Float_t  beta = 0;
    Float_t  gamma = 0;
    Float_t  mass  = -99;
	
    // if(fptc==1)
    //   expbeta = TMath::Sqrt(1-((massd*massd)/(ptot*ptot+massd*massd))); 
    // if(fptc==2)
    //   expbeta = TMath::Sqrt(1-((masst*masst)/(ptot*ptot+masst*masst))); 
    // if(fptc==3)
    //   expbeta = TMath::Sqrt(1-((mass3he*mass3he)/(ptot*ptot+mass3he*mass3he))); 
        
    if(fptc==3)
      pt = 2*pt;

    if(TMath::Abs(ptot) < pmax  && TMath::Abs(pt) < ptmax && TMath::Abs(pullTPC) <= fNsigma ){
      
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
      tNpidcluster     = atrack->GetTPCsignalN();

      if(fptc < 3){
      
	if(pt<1.5)   
	  ftree->Fill();
	else
	  if(hasTOF==1)
	    ftree->Fill();
      }
      else
	ftree->Fill();
    }//POI selection
	
  }//track loop
  
  PostData(1, fListHist);
  PostData(2, ftree);
  
}

//_____________________________________________________________________________
void AliAnalysisTaskNucleiv2PbPb18::OpenInfoCalbration(Int_t run )
{

  //foadb = TFile::Open("alien:///alice/cern.ch/user/l/lramona/CalibpPb2016/calibV0NoSDD.root");
  
  TFile* foadb = 0;
  if (!gGrid) TGrid::Connect("alien");

  if (fPeriod == 0)
    foadb = TFile::Open("alien:///alice/cern.ch/user/l/lramona/CalibPbPb2018/calibV0Run2Vtx10P118q.root");
  //TFile::Open("calibV0Run2Vtx10P118q.root");
  else
    foadb = TFile::Open("alien:///alice/cern.ch/user/l/lramona/CalibPbPb2018/calibV0Run2Vtx10P118r.root");
  //foadb = TFile::Open("calibV0Run2Vtx10P118r.root");
  
  
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
  if (fNHarm == 2.)
    contQxnam = (AliOADBContainer*) foadb->Get("fqxa2m");
  else
    contQxnam = (AliOADBContainer*) foadb->Get("fqxa3m");
    
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
  if (fNHarm == 2.)
    contQynam = (AliOADBContainer*) foadb->Get("fqya2m");
  else if (fNHarm == 3.)
    contQynam = (AliOADBContainer*) foadb->Get("fqya3m");
  else if (fNHarm == 4.)
    contQynam = (AliOADBContainer*) foadb->Get("fqya4m");
    
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
  if (fNHarm == 2.)
    contQxnas = (AliOADBContainer*) foadb->Get("fqxa2s");
  else
    contQxnas = (AliOADBContainer*) foadb->Get("fqxa3s");
    
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
  if (fNHarm == 2.)
    contQynas = (AliOADBContainer*) foadb->Get("fqya2s");
  else if (fNHarm == 3.)
    contQynas = (AliOADBContainer*) foadb->Get("fqya3s");
  else if (fNHarm == 4.)
    contQynas = (AliOADBContainer*) foadb->Get("fqya4s");
    
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
  if (fNHarm == 2.)
    contQxncm = (AliOADBContainer*) foadb->Get("fqxc2m");
  else
    contQxncm = (AliOADBContainer*) foadb->Get("fqxc3m");
    
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
  if (fNHarm == 2.)
    contQyncm = (AliOADBContainer*) foadb->Get("fqyc2m");
  else if (fNHarm == 3.)
    contQyncm = (AliOADBContainer*) foadb->Get("fqyc3m");
  else if (fNHarm == 4.)
    contQyncm = (AliOADBContainer*) foadb->Get("fqyc4m");
    
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
  if (fNHarm == 2.)
    contQxncs = (AliOADBContainer*) foadb->Get("fqxc2s");
  else
    contQxncs = (AliOADBContainer*) foadb->Get("fqxc3s");
    
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
  if (fNHarm == 2.)
    contQyncs = (AliOADBContainer*) foadb->Get("fqyc2s");
  else if (fNHarm == 3.)
    contQyncs = (AliOADBContainer*) foadb->Get("fqyc3s");
  else if (fNHarm == 4.)
    contQyncs = (AliOADBContainer*) foadb->Get("fqyc4s");
    
  if(!contQyncs){
    printf("OADB object fqycnm is not available in the file\n");
    return;
  }
  if(!(contQyncs->GetObject(run))){
    printf("OADB object fqycns is not available for run %i\n", run);
    return;
  }
  fQynsV0C = ((TH1D*) contQyncs->GetObject(run));
        
        
        
  /*    
	AliOADBContainer* contQx2cmESE = (AliOADBContainer*) foadb->Get("fqxc2m");
	if(!contQx2cmESE){
        printf("OADB object fqxc2m is not available in the file\n");
        return;
	}
	if(!(contQx2cmESE->GetObject(run))){
        printf("OADB object fqxc2m is not available for run %i\n", run);
        return;
	}
	fQx2mV0CESE = ((TH1D*) contQx2cmESE->GetObject(run));
        
        
	AliOADBContainer* contQy2cmESE = (AliOADBContainer*) foadb->Get("fqyc2m");
	if(!contQy2cmESE){
        printf("OADB object fqyc2m is not available in the file\n");
        return;
	}
	if(!(contQy2cmESE->GetObject(run))){
        printf("OADB object fqyc2m is not available for run %i\n", run);
        return;
	}
	fQy2mV0CESE = ((TH1D*) contQy2cmESE->GetObject(run));
        
        
        
        
	AliOADBContainer* contQx3cmESE = (AliOADBContainer*) foadb->Get("fqxc3m");
	if(!contQx3cmESE){
        printf("OADB object fqxc3m is not available in the file\n");
        return;
	}
	if(!(contQx3cmESE->GetObject(run))){
        printf("OADB object fqxc3m is not available for run %i\n", run);
        return;
	}
	fQx3mV0CESE = ((TH1D*) contQx3cmESE->GetObject(run));
        
        
	AliOADBContainer* contQy3cmESE = (AliOADBContainer*) foadb->Get("fqyc3m");
	if(!contQy3cmESE){
        printf("OADB object fqyc3m is not available in the file\n");
        return;
	}
	if(!(contQy3cmESE->GetObject(run))){
        printf("OADB object fqyc3m is not available for run %i\n", run);
        return;
	}
	fQy3mV0CESE = ((TH1D*) contQy3cmESE->GetObject(run));
  */

}
//_____________________________________________________________________________
void AliAnalysisTaskNucleiv2PbPb18::Terminate(Option_t *)
{ 
  // Terminate loop
  Printf("Terminate()");
}
