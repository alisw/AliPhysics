#include <Riostream.h>
#include "TTree.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TDatabasePDG.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliVEventHandler.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"

#include "AliAnalysisTaskEP.h"
#include "AliPIDResponse.h"
#include "AliPhysicsSelection.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliEventCuts.h"
#include "AliPPVsMultUtils.h"
#include "AliAnalysisUtils.h"

//_____ Event pool includes
#include "AliEventPoolManager.h"
//#include "AliPool.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEP)

AliAnalysisTaskEP::AliAnalysisTaskEP(): 
AliAnalysisTaskSE(), 
fOutput(0), 
fPoolMgr(0X0), 
fPIDResponse(0), 
fVevent(0), 
lESDevent(0),
lAODevent(0), 
fEventCuts(0), 
fESDtrackCuts(0), 
fHistVz(0),
fHistCentrality(0), 
fHisteventsummary(0),
fHisteventmult(0),
frame(2),
gPsiNSet(2),
fListTRKCorr(NULL),
fListNUACorr(NULL),
fListV0MCorr(NULL),
hAvgV0ChannelsvsVz(NULL),
hAvgQNXvsCentV0C(NULL),
hAvgQNYvsCentV0C(NULL),
hAvgQNXvsCentV0A(NULL),
hAvgQNYvsCentV0A(NULL),
fHistV0CPsiNEventPlane(NULL),
fHistV0APsiNEventPlane(NULL),
hV0CV0APsiNCorrelation(NULL),
hV0CTPCPsiNCorrelation(NULL),
hV0ATPCPsiNCorrelation(NULL),
fHCorrectV0ChWgt(NULL),
fHCorrectQNxV0C(NULL),
fHCorrectQNyV0C(NULL),
fHCorrectQNxV0A(NULL),
fHCorrectQNyV0A(NULL),
fHCorrectMCposChrg(NULL),
fHCorrectMCnegChrg(NULL),
fHCorrectNUAposChrg(NULL),  
fHCorrectNUAnegChrg(NULL), 
fHCorrectEVNTWGTChrg(NULL)
{

}

AliAnalysisTaskEP::AliAnalysisTaskEP(const char *name): 
AliAnalysisTaskSE(name), 
fOutput(0), 
fPoolMgr(0X0), 
fPIDResponse(0),
fVevent(0), 
lESDevent(0),
lAODevent(0),  
fEventCuts(0), 
fESDtrackCuts(0), 
fHistVz(0),
fHistCentrality(0), 
fHisteventsummary(0),
fHisteventmult(0),
frame(2),
gPsiNSet(2),
fListTRKCorr(NULL),
fListNUACorr(NULL),
fListV0MCorr(NULL),
hAvgV0ChannelsvsVz(NULL),
hAvgQNXvsCentV0C(NULL),
hAvgQNYvsCentV0C(NULL),
hAvgQNXvsCentV0A(NULL),
hAvgQNYvsCentV0A(NULL),
fHistV0CPsiNEventPlane(NULL),
fHistV0APsiNEventPlane(NULL),
hV0CV0APsiNCorrelation(NULL),
hV0CTPCPsiNCorrelation(NULL),
hV0ATPCPsiNCorrelation(NULL),
fHCorrectV0ChWgt(NULL),
fHCorrectQNxV0C(NULL),
fHCorrectQNyV0C(NULL),
fHCorrectQNxV0A(NULL),
fHCorrectQNyV0A(NULL),
fHCorrectMCposChrg(NULL),
fHCorrectMCnegChrg(NULL),
fHCorrectNUAposChrg(NULL),  
fHCorrectNUAnegChrg(NULL), 
fHCorrectEVNTWGTChrg(NULL)
{
  DefineInput(0, TChain::Class()); 
  DefineOutput(1, TList::Class()); 
}


AliAnalysisTaskEP::~AliAnalysisTaskEP()
{
  //------------------------------------------------
  // DESTRUCTOR
  //------------------------------------------------
    
  if (fOutput){
    delete fOutput;
    fOutput = 0x0;
  }
  
  if(fVevent) {
    delete fVevent;
    fVevent=0x0;
  }

  if(lESDevent) {
    delete lESDevent;
    lESDevent=0x0;
  }

  if(lAODevent) {
    delete lAODevent;
    lAODevent=0x0;
  }


  if(fESDtrackCuts) {
    delete fESDtrackCuts;
    fESDtrackCuts=0x0;
  }


  if(fListTRKCorr)  delete fListTRKCorr;
  if(fListNUACorr)  delete fListNUACorr;
  if(fListV0MCorr)  delete fListV0MCorr;



  if (fPoolMgr) { delete fPoolMgr; fPoolMgr = 0x0; }
  if (fPIDResponse) { delete fPIDResponse; fPIDResponse = 0x0; }
}

//________________________________________________________________________
void AliAnalysisTaskEP::UserCreateOutputObjects()
{
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  if(!fESDtrackCuts ){
    fESDtrackCuts = new AliESDtrackCuts();
    fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1); 
    fESDtrackCuts->SetPtRange(0.15,20.0);
    fESDtrackCuts->SetEtaRange(-0.8,0.8);
    //fESDtrackCuts->SetRapRange(-0.5,0.5);
  }


  EventMixing();

  fOutput = new TList();
  fOutput->SetOwner(); 
  OpenFile(1);
  fHisteventsummary            = new TH1F("fHisteventsummary", "Event cut summary", 5,0.0,5.0);       
  fHisteventmult            = new TH1F("fHisteventmult", "Event multiplicity", 10,0.0,10.0);       
  fHistVz            = new TH1F("fHistZVertex", "Z vertex distribution", 100,-10,10);       
  fHistCentrality = new TH1F("fHistCentrality", "Centrality distribution", 100,0,100);
  
  fOutput->Add(fHisteventsummary);
  fOutput->Add(fHisteventmult);
  fOutput->Add(fHistVz);
  fOutput->Add(fHistCentrality);



  ///Rihan: V0A/C histograms:
  //Int_t gPsiNSet = 2;   /// Note you can add this as configurable in Addtask, to set different harmonics.  
  hAvgV0ChannelsvsVz = new TProfile2D("hAvgV0ChannelsvsVz"," <V0 Ch(n)> ; Channel;  V_{z}(cm);",64,0,64,20,-10,10);
  fOutput->Add(hAvgV0ChannelsvsVz);
  ///  <Q> for Recentering:
  hAvgQNXvsCentV0C  = new TProfile("hAvgQNXvsCentV0C",Form("V0C <Qnx> vs Cent; Cent(%%); <q_{%d,x}>",gPsiNSet),90,0,90);
  fOutput->Add(hAvgQNXvsCentV0C);
  hAvgQNYvsCentV0C  = new TProfile("hAvgQNYvsCentV0C",Form("V0C <Qny> vs Cent; Cent(%%); <q_{%d,y}>",gPsiNSet),90,0,90);
  fOutput->Add(hAvgQNYvsCentV0C);
  hAvgQNXvsCentV0A  = new TProfile("hAvgQNXvsCentV0A",Form("V0A <Qnx> vs Cent; Cent(%%); <q_{%d,x}>",gPsiNSet),90,0,90);
  fOutput->Add(hAvgQNXvsCentV0A);
  hAvgQNYvsCentV0A  = new TProfile("hAvgQNYvsCentV0A",Form("V0A <Qny> vs Cent; Cent(%%); <q_{%d,y}>",gPsiNSet),90,0,90);
  fOutput->Add(hAvgQNYvsCentV0A);
 
  /// V0 Event Planes:
  fHistV0CPsiNEventPlane = new TH2F("fHistV0CPsiNEventPlane",Form("#Psi_{n}(V0C); Cent(%%); #Psi_{%d,V0C}(rad); events",gPsiNSet),18,0,90,50,0,3.14159);
  fOutput->Add(fHistV0CPsiNEventPlane);
  fHistV0APsiNEventPlane = new TH2F("fHistV0APsiNEventPlane",Form("#Psi_{n}(V0A); Cent(%%); #Psi_{%d,V0A}(rad); events",gPsiNSet),18,0,90,50,0,3.14159);
  fOutput->Add(fHistV0APsiNEventPlane);
  
  ///Resolution Hist for V0EP: Note: You also need TPC event plane, to get Resolution from 3 Event Plane correlation.
  hV0CV0APsiNCorrelation = new TProfile("hV0CV0APsiNCorrelation",Form("V0C-V0A Psi%d; Cent(%%); Resolution",gPsiNSet),90,0,90);
  fOutput->Add(hV0CV0APsiNCorrelation);
  hV0CTPCPsiNCorrelation = new TProfile("hV0CTPCPsiNCorrelation",Form("V0C-TPC Psi%d; Cent(%%); Resolution",gPsiNSet),90,0,90);
  fOutput->Add(hV0CTPCPsiNCorrelation);
  hV0ATPCPsiNCorrelation = new TProfile("hV0ATPCPsiNCorrelation",Form("V0A-TPC Psi%d; Cent(%%); Resolution",gPsiNSet),90,0,90);
  fOutput->Add(hV0ATPCPsiNCorrelation);








  if(fListTRKCorr){
    std::cout<<"\n UserCreateOutputObject::Info() Tlist for MC tracking Efficiency Found.!!\n"<<std::endl;
  }
  else{
    std::cout<<"\n\n ******* WARNING No TList for Trk Efficiency Correction!!\n using TrkWgt = 1.0 \n "<<std::endl;
  }

  if(fListNUACorr){
    std::cout<<"\n UserCreateOutputObject::Info() Tlist for NUA Correction Found.!!\n"<<std::endl;
    //fListNUACorr->ls(); 
  }
  else{
    std::cout<<"\n\n ******* WARNING No TList NUA Correction!!\n using NUAWgt = 1.0 \n "<<std::endl;
  }

  if(fListV0MCorr){
    std::cout<<"\n UserCreateOutputObject::Info() Tlist for V0M Correction Found.!!\n"<<std::endl;
    //fListV0MCorr->ls(); 
  }
  else{
    std::cout<<"\n\n ******* WARNING No TList EVNTWGT Correction!!\n using EVNTWGTWgt = 1.0 \n "<<std::endl;
  }
 






  PostData(1, fOutput);
}


//________________________________________________________________________
void AliAnalysisTaskEP::UserExec(Option_t *)
{
     
  // Main loop
  
  
  lAODevent = dynamic_cast <AliAODEvent*> (InputEvent());


  
  if (!(lAODevent)) {
    AliWarning("ERROR: AOD event not available \n");
    PostData(1, fOutput);
    return;
  }

  ////trigger/////////////
  UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isSelected = 0;
  isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;
 

  if ( ! isSelected)
    {
      
      PostData(1, fOutput);
      return;
    }

  fHisteventsummary->Fill(0.5);

  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!fVevent){
    printf("ERROR: fVevent not available\n");
    return;
  }
  

  //event selection  
  const AliVVertex *vertex = fVevent->GetPrimaryVertex();
  if (!GoodEvent(vertex)) return; 
  double zv=vertex->GetZ();
  ////////////////////////////////


  //////centrality selection/////////
  Float_t lV0M;
  Float_t centrCL1=0; 
  Int_t lEvSelCode = 300;
  AliMultSelection *MultSelection = (AliMultSelection*) lAODevent -> FindListObject("MultSelection");
  if( !MultSelection)
    {
      AliWarning("AliMultSelection object not found!");
      
      PostData(1, fOutput);
      return;
    }
  else
    {
      lV0M = MultSelection->GetMultiplicityPercentile("V0M"); 
      centrCL1 = MultSelection->GetMultiplicityPercentile("CL1");    
    }
    
  

  fHistVz->Fill(zv);
  fHistCentrality->Fill(lV0M);
  fHisteventsummary->Fill(1.5);
 
  //getting correction histograms for eff and NUA
  Int_t runNumber = lAODevent->GetRunNumber();
  if(fListTRKCorr) GetMCCorrectionHist(runNumber,lV0M);  //use centrality dependent MC efficiency (Temporary!!)

  

  if(fListNUACorr){
    GetNUACorrectionHist(runNumber,0);         //Charge
  }
  
  if(fListV0MCorr){
    GetV0MCorrectionHist(runNumber,0);         //Charge
  }


  
  /// ***** Prottay: You need To estimate TPC event Plane, to get resolution of V0 EP
  Double_t fPsiNTPC=0;

  /*
  // Run a TPC track loop (with global track (FB96), |eta| < 0.8,  0.2 < pT < 2.0 GeV/c)
  // and estimate fPsiNTPC  
  // If your V0 EP is super flat then you do not need a very perfect TPC EP (i.e., no need for NUA)
  // But if you would later use TPC EP as EP for analysis, then you should use NUA correction to 
  // Make TPC EP flat.
  */
  
  



  /// V0 EP implementation, Start.

  Double_t fQnxV0C=0, fQnyV0C=0, fQnxV0A=0, fQnyV0A=0; 
  Bool_t kPassV0 = GetGainCorrectedV0Qvector(lAODevent, zv, gPsiNSet, fQnxV0C, fQnyV0C, fQnxV0A, fQnyV0A); //Second Order EP
  
  if(!kPassV0) return; /// V0 does not have signal for this event.              

  //Printf("\n ::UserExec() ==> Passed V0 Selection: GetGainCorrectedV0Qvector(), fQnxV0C:%f  ...\n",fQnxV0C);
    
  ApplyV0XqVectRecenter(centrCL1, gPsiNSet, fQnxV0C, fQnyV0C, fQnxV0A, fQnyV0A); // DO the <q> recentering. 

  
  ////---- fill the <q> vector vs Cent (for recentering)
  hAvgQNXvsCentV0C->Fill(centrCL1,fQnxV0C);   /// to Avoid self correlation, V0 <Q> is filled with CL1 centrality! 
  hAvgQNYvsCentV0C->Fill(centrCL1,fQnyV0C);  
  hAvgQNXvsCentV0A->Fill(centrCL1,fQnxV0A);
  hAvgQNYvsCentV0A->Fill(centrCL1,fQnyV0A);
  
  ///--------> Now, Get V0A and V0C Event Planes 
  Double_t fPsiNV0C = 0., fPsiNV0A = 0.;
    
  fPsiNV0C = (1./gPsiNSet)*TMath::ATan2(fQnyV0C,fQnxV0C);
  if(fPsiNV0C < 0) fPsiNV0C += TMath::TwoPi()/gPsiNSet;
  fPsiNV0A = (1./gPsiNSet)*TMath::ATan2(fQnyV0A,fQnxV0A);
  if(fPsiNV0A < 0) fPsiNV0A += TMath::TwoPi()/gPsiNSet;

  fHistV0CPsiNEventPlane->Fill(lV0M, fPsiNV0C);    ///*** Define This Histogram then Uncomment!
  fHistV0APsiNEventPlane->Fill(lV0M, fPsiNV0A);    ///*** Define This Histogram then Uncomment!

  /// V0A, V0C Resolutions:
  hV0CV0APsiNCorrelation->Fill(lV0M,TMath::Cos(2*fPsiNV0A - 2*fPsiNV0C));
  hV0CTPCPsiNCorrelation->Fill(lV0M,TMath::Cos(2*fPsiNTPC - 2*fPsiNV0C));
  hV0ATPCPsiNCorrelation->Fill(lV0M,TMath::Cos(2*fPsiNTPC - 2*fPsiNV0A));

  /// V0 EP implementation, End.

  //Printf("\n ::UserExec() ==> Estimated V0 EP: fPsiNV0A: %f\n",fPsiNV0A);





  
    
  PostData(1, fOutput);
}
//________________________________________________________________________
void AliAnalysisTaskEP::Terminate(Option_t *)
{
  
}

//________________________________________________________________________

Bool_t AliAnalysisTaskEP::GoodEvent(const AliVVertex *vertex) 

{

  Bool_t EventAccepted;
  EventAccepted = fEventCuts.AcceptEvent(fVevent);
  if (!EventAccepted)
    {
      PostData(1, fOutput);
      return kFALSE;
    }   


  if (!vertex)
    {
      PostData(1, fOutput);
      return kFALSE;
    }  

  if (vertex->GetNContributors() < 1)
    {
      PostData(1, fOutput);
      return kFALSE;
    }
  
  if (TMath::Abs(vertex->GetZ()) > 10.0)
    {
      PostData(1, fOutput);
      return kFALSE;
    }
  
  
  if( lAODevent->IsIncompleteDAQ() ) 
    {
      PostData(1, fOutput);
      return kFALSE;
    }

    AliAnalysisUtils *fUtils = new AliAnalysisUtils();
    /* if(fUtils->IsSPDClusterVsTrackletBG(fVevent))
    {
      PostData(1, fOutput);
      return kFALSE;
    }
    */

    //fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE); 


  return kTRUE;

}
//-----------------------------------------------
Bool_t AliAnalysisTaskEP::IsPion(AliVTrack *vtrack)
{

  
  AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(vtrack);
  Double_t nsigmatpcpion=TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodtrack, AliPID::kPion));
  Double_t nsigmatofpion=TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodtrack, AliPID::kPion));
  Bool_t TOFHIT=kFALSE;

  TOFHIT = HasTOF(aodtrack);

  
  if(!TOFHIT)
    {
      if(nsigmatpcpion>2.0)return kFALSE;
    }

  else
    {
      if(nsigmatofpion>3.0)return kFALSE;
    }

  return kTRUE;
}

Bool_t AliAnalysisTaskEP::IsKaon(AliVTrack *vtrack)
{

  AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(vtrack);
  Double_t nsigmatpckaon=TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodtrack, AliPID::kKaon));
  Double_t nsigmatofkaon=TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodtrack, AliPID::kKaon));
  Bool_t TOFHIT=kFALSE;

  TOFHIT = HasTOF(aodtrack);

  if(!TOFHIT)
    {
      if(nsigmatpckaon>2.0)return kFALSE;
    }
  else
    {
      if(nsigmatofkaon>3.0)return kFALSE;
    }

    
  return kTRUE;
}


Bool_t AliAnalysisTaskEP::HasTOF(AliAODTrack *track)
{
  bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  const float len = track->GetIntegratedLength();
  bool hasTOF = hasTOFout && hasTOFtime && (len > 350.);
  return hasTOF;
}



//_________________________________________________________________________________________________
void AliAnalysisTaskEP::EventMixing()
{

  ////////////////////////////
  // Set-up for Mixed Event //
  ////////////////////////////
  const Int_t trackDepth = 1000000;
  const Int_t poolsize = 100;
  const Int_t nmix = 5;
  Double_t centralityBins[] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0};
  Double_t vertexBins[] = {-10.0,-8.0,-6.0,-4.0,-2.0,0.0,2.0,4.0,6.0,8.0,10.0};
  const Int_t nCentralityBins = 10;
  const Int_t nVertexBins = 10;
  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centralityBins, nVertexBins, vertexBins);
  fPoolMgr->SetTargetValues(trackDepth,0.1,nmix);
}
//___________________________________________________________





void AliAnalysisTaskEP::GetNUACorrectionHist(Int_t run, Int_t kParticleID)
{

  if(fListNUACorr){

    if(kParticleID==0){ //charge
      fHCorrectNUAposChrg = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Pos_Cent%d_Run%d",0,run)); 
      fHCorrectNUAnegChrg = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Neg_Cent%d_Run%d",0,run));
    }
    else{
      fHCorrectNUAposChrg=NULL;
      fHCorrectNUAnegChrg=NULL;
    }
    
  }//------> if list Exist

  // else {
  //   printf("\n ******** Warning: No NUA Correction File/List...!! \n Run= %d, Use NUA Wgt = 1.0 ********* \n",run);
  // }
}





void AliAnalysisTaskEP::GetV0MCorrectionHist(Int_t run, Int_t kParticleID)
{

  if(fListV0MCorr){
    //cout<<"*****************AT LAST I GOT INTO EVNTWGT CORRECTION FUNCTION******************************"<<endl;
    /*if(kParticleID==0){ //charge
      fHCorrectEVNTWGTChrg = (TH1F *) fListV0MCorr->FindObject(Form("hwgtCharge_Run%d",run));     
    }
    */
    //Rihan: V0 Channel Gains:
    fHCorrectV0ChWgt = (TH2F *) fListV0MCorr->FindObject(Form("hWgtV0ChannelsvsVzRun%d",run)); //Channel Gains
    //if(fHCorrectV0ChWgt) printf("\n ==> Info:: V0 Channel Weights Found for Run %d \n ",run);    
    fHCorrectQNxV0C = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQNxvsCentV0CRun%d",run));  ///V0A,V0C <q> Vectors..
    fHCorrectQNyV0C = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQNyvsCentV0CRun%d",run));    
    fHCorrectQNxV0A = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQNxvsCentV0ARun%d",run));
    fHCorrectQNyV0A = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQNyvsCentV0ARun%d",run));
    // if(fHCorrectQNxV0C && fHCorrectQNxV0A) printf(" ==> Info:: V0A,V0C <Q> Found for Run %d \n ",run);    
  }///if List Exist!
  else{
    // cout<<"*****************AT LAST I GOT INTO EVNTWGT CORRECTION FUNCTION BUT LIST NOT PRESENT******************************"<<endl;
    fHCorrectEVNTWGTChrg=NULL;
  
    fHCorrectV0ChWgt=NULL;
    fHCorrectQNxV0C=NULL;
    fHCorrectQNyV0C=NULL;
    fHCorrectQNxV0A=NULL;
    fHCorrectQNyV0A=NULL;
  }
}




void AliAnalysisTaskEP::GetMCCorrectionHist(Int_t run,Float_t centr){

  if(fListTRKCorr) {
    //cout<<"\n =========> Info: Found TList with MC Tracking Corr Histograms <=========== "<<endl;
    /// Default: Centrality Independent MC efficiency:
    
    fHCorrectMCposChrg =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyChrgPos");
    fHCorrectMCnegChrg =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyChrgNeg");


  }
  // else if(!fListTRKCorr){
  //   std::cout<<"\n\n !!!!**** Warning : No FB Efficiency Correction, run = "<<run<<"..!!!**** \n using MC TrkWgt = 1.0 \n"<<std::endl;
  // }
}



/// Rihan--> V0 EP function #1
Bool_t AliAnalysisTaskEP::GetGainCorrectedV0Qvector(AliAODEvent *faod,Double_t fVtxZ,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A){

  const AliAODVZERO *fAODV0 = (AliAODVZERO *) faod->GetVZEROData();
  Float_t fMultV0 = 0.;
  Float_t fPhiV0  = 0.;
  Float_t fV0chGain = 1.0;

  Double_t fQxV0CHarmN=0,fQyV0CHarmN=0,fQxV0AHarmN=0,fQyV0AHarmN=0;

  Double_t fSumMV0A = 0;
  Double_t fSumMV0C = 0;
  Int_t ibinV0=0;

  for(int iV0 = 0; iV0 < 64; iV0++) { //0-31 is V0C, 32-63 VOA

    fMultV0 = fAODV0->GetMultiplicity(iV0);

    /// V0 Channel Gain Correction:
    if(fHCorrectV0ChWgt){
      ibinV0    = fHCorrectV0ChWgt->FindBin(fVtxZ,iV0);
      fV0chGain = fHCorrectV0ChWgt->GetBinContent(ibinV0); 
    }
    
    fMultV0 = fMultV0*fV0chGain;   //Corrected Multiplicity
    
    hAvgV0ChannelsvsVz->Fill(iV0+0.5, fVtxZ, fMultV0); //1st pass

    fPhiV0  = TMath::PiOver4()*(0.5 + iV0 % 8);

    if(iV0 < 32){
      qnxV0C   += TMath::Cos(gPsiN*fPhiV0) * fMultV0;
      qnyV0C   += TMath::Sin(gPsiN*fPhiV0) * fMultV0;
      fSumMV0C += fMultV0;
    }
    else if(iV0 >= 32){
      qnxV0A   += TMath::Cos(gPsiN*fPhiV0) * fMultV0;
      qnyV0A   += TMath::Sin(gPsiN*fPhiV0) * fMultV0;
      fSumMV0A += fMultV0;
    } 
  }///V0 Channel loop

  /// Now the q vectors:
  if(fSumMV0A<=1e-4 || fSumMV0C<=1e-4){
    qnxV0C = 0;
    qnyV0C = 0;
    qnxV0A = 0;
    qnyV0A = 0;    
    return kFALSE;       
  }
  else{
    if(fSumMV0C>0){
      qnxV0C = qnxV0C/fSumMV0C;
      qnyV0C = qnyV0C/fSumMV0C;
    }
    if(fSumMV0A>0){
      qnxV0A = qnxV0A/fSumMV0A;
      qnyV0A = qnyV0A/fSumMV0A;
    }
    return kTRUE;  
  }
  
}//<--- Rihan: V0 EP function #1 end.



/// Rihan: --> V0 EP function #2
void  AliAnalysisTaskEP::ApplyV0XqVectRecenter(Float_t fCent,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A){

  //** NOTE: Proper File path which contain <q> for Desired harmonic 'N' should be set in AddTask!!
  
  Int_t icentbin = 0;
  Double_t avgqx=0,avgqy=0; 
  //cout<<"Info::V0<q>::Raw!! =>" <<" N:"<<gPsiN<<"\t"<<" qnxV0C: "<<qnxV0C<<"\tqnyV0C: "<<qnyV0C<<" qnxV0A: "<<qnxV0A<<"\tqnyV0A: "<<qnyV0A<<endl;   
  if(fHCorrectQNxV0C && fHCorrectQNyV0C){
    icentbin = fHCorrectQNxV0C->FindBin(fCent);
    avgqx = fHCorrectQNxV0C->GetBinContent(icentbin);
    avgqy = fHCorrectQNyV0C->GetBinContent(icentbin);
    qnxV0C -= avgqx;
    qnyV0C -= avgqy;      
    //cout<<" V0C PsiN: "<<gPsiN<<" Cent: "<<fCent<<"\t <qx> "<<avgqx<<"\t <qy> "<<avgqy<<endl;
  }
  if(fHCorrectQNxV0A && fHCorrectQNyV0A){
    icentbin = fHCorrectQNxV0A->FindBin(fCent);
    avgqx = fHCorrectQNxV0A->GetBinContent(icentbin);
    avgqy = fHCorrectQNyV0A->GetBinContent(icentbin);
    qnxV0A -= avgqx;
    qnyV0A -= avgqy;           
    //cout<<" V0A PsiN: "<<gPsiN<<" Cent: "<<fCent<<"\t <qx> "<<avgqx<<"\t <qy> "<<avgqy<<endl;
  }
  //cout<<"Info::V0<q>::Corrected =>" <<" N:"<<gPsiN<<"\t"<<" qnxV0C: "<<qnxV0C<<"\tqnyV0C: "<<qnyV0C<<" qnxV0A: "<<qnxV0A<<"\tqnyV0A: "<<qnyV0A<<endl;
  
  return;
}
/// <--- Rihan: V0 EP function #2 end













