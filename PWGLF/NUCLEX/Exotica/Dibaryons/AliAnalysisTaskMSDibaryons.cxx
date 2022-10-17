// std
#include <iostream>
#include <stdio.h>

// ROOT includes
#include <THashList.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1.h>
#include <TFile.h>
#include <TDatabasePDG.h>

// AliRoot includes
#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisUtils.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>
#include <AliVEvent.h>
#include <AliAODVZERO.h>
#include <AliAODVertex.h>
#include <AliAODInputHandler.h> 
#include <AliAODHandler.h> 
#include <AliAODHeader.h>
#include <AliAODTrack.h> 
#include <AliAODPid.h> 
#include <AliPIDResponse.h>
#include <AliAODv0.h>
#include <AliAODcascade.h>
#include <AliAODTracklets.h>
#include <AliEventPoolManager.h>
//#include <AliCentrality>

//20220524 Flow study
#include <AliMultSelection.h>
#include <AliEventplane.h>
#include <AliQnCorrectionsManager.h>
#include <AliAnalysisTaskFlowVectorCorrections.h>
#include <AliQnCorrectionsQnVector.h>

#include "AliAnalysisTaskMSDibaryons.h"

//class AliAnalysisTaskMSDibaryons;

using namespace std;

ClassImp(AliAnalysisTaskMSDibaryons)

//_____________________________________________________________________________
AliAnalysisTaskMSDibaryons::AliAnalysisTaskMSDibaryons() : 
AliAnalysisTaskSE(),
  fEvent(0x0),
  fESD(0x0),
  fAOD(0),
  fHeader(0),
  fPIDResponse(0),
  ftrigBit(0),
  fEstimator("V0M"),
  fMultSelection(0x0),
  fCentralityMain(0.),
  fCentralityMin(0.),
  fCentralityMax(0.),
  fPoolManagerlambda(0),
  fPoolManagerantilambda(0),
  fPoolManagerxi(0),
  fPoolManagerxip(0),
  fPoolManagerproton(0),
  fPoolManagerantiproton(0),
  fPoolManagerlam(0),
  fPoolManagerantilam(0),
fPoolManagerx(0),
fPoolManagerxp(0),
fPoolManageromega(0),
fPoolManageromegap(0),
  fPoolManagerp(0),
  fPoolManagerantip(0),
  fIsFlowTask(kFALSE),
  fFlowQnVectorMgr(0x0),
  fOutputList(0)
{
 
}
 
//_____________________________________________________________________________
AliAnalysisTaskMSDibaryons::AliAnalysisTaskMSDibaryons(const char *name) : 
  AliAnalysisTaskSE(name),
  fEvent(0x0),
  fESD(0x0),
  fAOD(0),
  fHeader(0),
  fPIDResponse(0),
  ftrigBit(0),
  fEstimator("V0M"),
  fMultSelection(0x0),
  fCentralityMain(0.),
  fCentralityMin(0.),
  fCentralityMax(0.),
  fPoolManagerlambda(0),
  fPoolManagerantilambda(0),
  fPoolManagerxi(0),
  fPoolManagerxip(0),
  fPoolManagerproton(0),
  fPoolManagerantiproton(0),
  fPoolManagerlam(0),
  fPoolManagerantilam(0),
  fPoolManagerx(0),
  fPoolManagerxp(0),
  fPoolManageromega(0),
  fPoolManageromegap(0),
  fPoolManagerp(0),
  fPoolManagerantip(0),
  fIsFlowTask(kFALSE),
  fFlowQnVectorMgr(0x0),
  fOutputList(0)
{
  // Constructor
  DefineInput(0,TChain::Class());
  DefineOutput(1,THashList::Class());


}

//________________________________________________________________________
AliAnalysisTaskMSDibaryons::~AliAnalysisTaskMSDibaryons()
{
  // destructor
  if(fOutputList){
    delete fOutputList;
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskMSDibaryons::UserCreateOutputObjects()
{

  //==========define AliEventPoolMgr                                                                                                          
  Int_t MaxNEvents = 5;
  Int_t MaxNLambda = 10;
  //Int_t nMultiBins = 12;                                                                                                                 
  //Double_t multBins[]={0,10,20,30,40,50,60,70,80,90,100,110,120};                                                                      
  //Int_t nMultiBins = 17;
  //Double_t multBins[]={1,100,200,300,400,500,600,700,800,900,1000,1100,1200,1400,1600,1800,2000,3000};
  //Int_t nZvtxBins = 4;
  //Double_t vertexBins[]={-10,-5,0,5,10};
  //Int_t nMultiBinslamxi = 15;
  //Double_t multBinslamxi[]={1,5,10,15,20,25,30,35,40,50,60,70,80,90,100,200};

  // for Pb-Pb analysis
  Int_t nZvtxBins = 10;
  Double_t vertexBins[]={-10,-8,-6,-4,-2,0,2,4,6,8,10};
  //Int_t nCentBins = 8;
  //Double_t centBins[]={0,10,20,30,40,50,60,80,100};
  Int_t nCentBins = 13;
  Double_t centBins[]={0,5,10,15,20,25,30,35,40,45,50,60,80,100};
  Int_t nEPBins = 14;
  Double_t EPBins[]={TMath::Pi()*(-7/6),TMath::Pi()*(-6/6),TMath::Pi()*(-5/6),TMath::Pi()*(-4/6),TMath::Pi()*(-3/6),TMath::Pi()*(-2/6),TMath::Pi()*(-1/6),0,
    		     TMath::Pi()*(1/6),TMath::Pi()*(2/6),TMath::Pi()*(3/6),TMath::Pi()*(4/6),TMath::Pi()*(5/6),TMath::Pi()*(6/6),TMath::Pi()*(7/6)};

  fPoolManagerlambda     =new AliEventPoolManager(MaxNEvents,MaxNLambda,nCentBins,(Double_t*)centBins,nZvtxBins,(Double_t *)vertexBins);
  fPoolManagerantilambda =new AliEventPoolManager(MaxNEvents,MaxNLambda,nCentBins,(Double_t*)centBins,nZvtxBins,(Double_t *)vertexBins);
  fPoolManagerxi         =new AliEventPoolManager(MaxNEvents,MaxNLambda,nCentBins,(Double_t*)centBins,nZvtxBins,(Double_t *)vertexBins);
  fPoolManagerxip        =new AliEventPoolManager(MaxNEvents,MaxNLambda,nCentBins,(Double_t*)centBins,nZvtxBins,(Double_t *)vertexBins);
  fPoolManagerproton     =new AliEventPoolManager(MaxNEvents,MaxNLambda,nCentBins,(Double_t*)centBins,nZvtxBins,(Double_t *)vertexBins);
  fPoolManagerantiproton =new AliEventPoolManager(MaxNEvents,MaxNLambda,nCentBins,(Double_t*)centBins,nZvtxBins,(Double_t *)vertexBins);
  fPoolManagerlam        =new AliEventPoolManager(MaxNEvents,MaxNLambda,nCentBins,(Double_t*)centBins,nZvtxBins,(Double_t *)vertexBins);
  fPoolManagerantilam    =new AliEventPoolManager(MaxNEvents,MaxNLambda,nCentBins,(Double_t*)centBins,nZvtxBins,(Double_t *)vertexBins);
  fPoolManagerx          =new AliEventPoolManager(MaxNEvents,MaxNLambda,nCentBins,(Double_t*)centBins,nZvtxBins,(Double_t *)vertexBins);
  fPoolManagerxp         =new AliEventPoolManager(MaxNEvents,MaxNLambda,nCentBins,(Double_t*)centBins,nZvtxBins,(Double_t *)vertexBins);
  //========== 

  fOutputList=new THashList();
  fOutputList->SetOwner(kTRUE);

  //========== Event cut ==========
  fOutputList->Add(new TH1F("fEvtCounter","",1,0,1));
  fOutputList->Add(new TH1F("fEvtPassCut","",1,0,1));
  fOutputList->Add(new TH1F("fEvtVtxX","",400,-20.,20.));
  fOutputList->Add(new TH1F("fEvtVtxY","",400,-20.,20.));
  fOutputList->Add(new TH1F("fEvtVtxZ","",400,-20.,20.));
  fOutputList->Add(new TH1F("fEvtVtxTrk","",1000,0.,1000.));
  fOutputList->Add(new TH1F("fSPDVtxZ","",400,-20.,20.));
  fOutputList->Add(new TH1F("fSPDVtxTrk","",1000,0.,1000.));
  fOutputList->Add(new TH2F("fSPDVtxCor","",400,-20.,20.,400,-20.,20.));
  fOutputList->Add(new TH1F("fSPDVtxDisp","",300,0.,1.5));

  fOutputList->Add(new TH1F("hCentralityV0M","centrality V0M",101,0,101));
  fOutputList->Add(new TH1F("hCentralityCL0","centrality CL0",101,0,101));
  fOutputList->Add(new TH1F("hCentralityCL1","centrality CL1",101,0,101));
  fOutputList->Add(new TH1F("hCentralityZNA","centrality ZNA",101,0,101));
  fOutputList->Add(new TH1F("hCentralityZNC","centrality ZNC",101,0,101));
  fOutputList->Add(new TH1F("hCentralityV0A","centrality V0A",101,0,101));
  fOutputList->Add(new TH1F("hCentralityV0C","centrality V0C",101,0,101));
  fOutputList->Add(new TH1F("hCentralityMain","centrality Main",101,0,101));

  fOutputList->Add(new TH1F("hCentralityV0M_cent","centrality V0M MostCent",101,0,101));
  fOutputList->Add(new TH1F("hCentralityV0M_semicent","centrality V0M SemiCent",101,0,101));
  fOutputList->Add(new TH1F("hCentralityV0M_mb","centrality V0M MB",101,0,101));

  //========== multiplicity ==========
  TH2F *fmultiplicity                      =new TH2F("fmultiplicity","",2000,0,200,40,-20,20);
  TH1F *fmultiplicity_lambdalambda         =new TH1F("fmultiplicity_lambdalambda","",2000,0,200);
  TH1F *fmultiplicity_protonxi             =new TH1F("fmultiplicity_protonxi","",2000,0,200);
  TH1F *fmultiplicity_lambdaxi             =new TH1F("fmultiplicity_lambdaxi","",2000,0,200);

  fOutputList->Add(new TH2F("fmultiplicity","",2000,0,200,40,-20,20));
  fOutputList->Add(new TH1F("fmultiplicity_lambdalambda","",2000,0,200));
  fOutputList->Add(new TH1F("fmultiplicity_protonxi","",2000,0,200));            
  fOutputList->Add(new TH1F("fmultiplicity_lambdaxi","",2000,0,200));            

  //========== Proton ==========
  fOutputList->Add(new TH1F("fTrkEta","",100,-2,2));
  fOutputList->Add(new TH1F("fTrkPt","",1000,0,5));
  fOutputList->Add(new TH1F("fTrkTPCclus","",1000,0,200));
  fOutputList->Add(new TH1F("fTrkTPCcrsR","",1000,0,200));
  fOutputList->Add(new TH1F("fTrkTPCclusF","",1000,0,200));
  fOutputList->Add(new TH1F("fTrkDCAxy","",100,-3,3));
  fOutputList->Add(new TH1F("fTrkDCAz","",100,-3,3));
  fOutputList->Add(new TH2F("fPvsTPCsignal","",2000,0,10,20000,0,1000));
  fOutputList->Add(new TH2F("fPvsTOFsignal","",500,0,5,200,0,2));
  fOutputList->Add(new TH1F("fPt_proton","",1000,0,10));
  fOutputList->Add(new TH2F("PID_cut","",2000,0,10,20000,-1000,1000));

  //========== Lambda ===========
  //daughter track select contents
  fOutputList->Add(new TH1F("fEta","",100,-2,2));
  fOutputList->Add(new TH1F("fTPC","",200,0,200));
  fOutputList->Add(new TH1F("fDCA","",100,0,1));
  fOutputList->Add(new TH1F("fPIDproton","",60,0,6));
  fOutputList->Add(new TH1F("fPIDpion","",60,0,6));  
  //lambda select contents
  fOutputList->Add(new TH1F("fV0pt","",1000,0,10));  
  fOutputList->Add(new TH1F("fV0vtxx","",500,0,500));
  fOutputList->Add(new TH1F("fV0vtxy","",500,0,500));
  fOutputList->Add(new TH1F("fV0vtxz","",500,0,500));
  fOutputList->Add(new TH1F("fTransRadius","",30000,0,300));
  fOutputList->Add(new TH1F("fDCAdaugTov0Vtx","",1000,0,10));
  fOutputList->Add(new TH1F("fCPA","",200,-1,1));
  fOutputList->Add(new TH1F("hInvMassK0s","",10000,0,10));
  //v0 selection                                                                                    
  fOutputList->Add(new TH2F("hInvMassLambda","",20000,0,10,10000,0,10));
  fOutputList->Add(new TH2F("fPtv0_lambda","",1000,0,10,10000,0,10));
  fOutputList->Add(new TH1F("fPhi-Psi_lambda","",2000,-10,10));
  fOutputList->Add(new TH1F("fPhi_lambda","",2000,-10,10));

  //=========== Xi ==========
  //cascade daughter select contents
  fOutputList->Add(new TH1F("fcascdauEta","",100,-2,2));
  fOutputList->Add(new TH1F("fcascdauEtab","",100,-2,2));
  fOutputList->Add(new TH1F("fcascdauTPC","",200,0,200));
  fOutputList->Add(new TH1F("fcascdauTPCb","",200,0,200));
  fOutputList->Add(new TH1F("fcascdauPt","",1000,0,10));
  fOutputList->Add(new TH1F("fcascdauPtb","",1000,0,10));
  fOutputList->Add(new TH1F("fcascdaudca","",100,0,1)); 
  fOutputList->Add(new TH1F("fcascdaudcab","",100,0,1));
  fOutputList->Add(new TH1F("fcascdauPIDproton","",60,0,6));
  fOutputList->Add(new TH1F("fcascdauPIDpion","",60,0,6));
  fOutputList->Add(new TH1F("fcascbachPIDpion","",60,0,6));
  //cascade V0 select contents
  fOutputList->Add(new TH1F("fcascV0Transradius","",30000,0,300));
  fOutputList->Add(new TH1F("fcascV0daugdca","",1000,0,10));    
  fOutputList->Add(new TH1F("fcascV0dca","",1000,0,10));        
  fOutputList->Add(new TH1F("fcascV0cpa","",200,-1,1));        
  fOutputList->Add(new TH1F("hInvMassXidecayLambda","",10000,0,10));
  //Xi select contents
  fOutputList->Add(new TH1F("XiTranverseradius","",30000,0,300)); 
  fOutputList->Add(new TH1F("dcaXiDaughters","",1000,0,10));    
  fOutputList->Add(new TH1F("XiCPA","",200,-1,1));      
  //fOutputList->Add(new TH1F("hInvMassOmega","",10000,0,10));       
  //invariant mass xi          
  fOutputList->Add(new TH2F("hInvMassXi","",20000,0,10,10000,0,10));
  fOutputList->Add(new TH2F("fPtcascade_xi","",1000,0,10,10000,0,10));
  fOutputList->Add(new TH1F("fPtv0casc","",1000,0,10));
  fOutputList->Add(new TH1F("fPtBachelor","",1000,0,10));
  fOutputList->Add(new TH1F("fPhiXidecayLambda","",1000,0,10));
  fOutputList->Add(new TH1F("fPhi-Psi_xi","",2000,-10,10));

  //========== baryon-baryon FG ==========
  //LambdaLambda
  fOutputList->Add(new TH2F("hInvMassLambdaLambda_onlyprompt","",1000,0,10,1000,2,3));
  fOutputList->Add(new TH1F("hInvMassLambdaLambda","",10000,0,10));
  fOutputList->Add(new TH2F("hInvMassLambdaLambdavspT","",1000,0,10,1000,0,10));
  //fOutputList->Add(new TH2F("hInvMassLambdaLambda_centvsmass","",101,0,101,10000,0,10));
  //ProtonXi
  fOutputList->Add(new TH1F("hInvMassProtonXi","",10000,0,10));
  fOutputList->Add(new TH2F("hInvMassProtonXivspT","",1000,0,10,1000,0,10));
  //fOutputList->Add(new TH2F("hInvMassProtonXi_centvsmass","",101,0,101,10000,0,10));
  //LambdaXi
  fOutputList->Add(new TH1F("hInvMassLambdaXi","",10000,0,10));
  fOutputList->Add(new TH2F("hInvMassLambdaXivspT","",1000,0,10,1000,0,10));
  //fOutputList->Add(new TH2F("hInvMassLambdaXi_centvsmass","",101,0,101,10000,0,10));
  //========== baryon-baryon BG ==========
  //LambdaLambda                                                                                                                              
  fOutputList->Add(new TH1F("hInvMassLambdaLambda_evtpool","",10000,0,10));
  fOutputList->Add(new TH2F("hInvMassLambdaLambdavspT_evtpool","",1000,0,10,1000,0,10));
  //fOutputList->Add(new TH2F("hInvMassLambdaLambda_evtpool_deltaphi","",180,0,180,10000,0,10));
  //fOutputList->Add(new TH2F("hInvMassLambdaLambda_evtpool_centvsmass","",101,0,101,10000,0,10));
  //ProtonXi                                                                                                                                     
  fOutputList->Add(new TH1F("hInvMassProtonXi_evtpool","",10000,0,10));   
  fOutputList->Add(new TH2F("hInvMassProtonXivspT_evtpool","",1000,0,10,1000,0,10));   
  //fOutputList->Add(new TH2F("hInvMassProtonXi_evtpool_deltaphi","",180,0,180,10000,0,10));
  //fOutputList->Add(new TH2F("hInvMassProtonXi_evtpool_centvsmass","",101,0,101,10000,0,10));  
  //LambdaXi                                                                                                                                
  fOutputList->Add(new TH1F("hInvMassLambdaXi_evtpool","",10000,0,10));      
  fOutputList->Add(new TH2F("hInvMassLambdaXivspT_evtpool","",1000,0,10,1000,0,10));      
  //fOutputList->Add(new TH2F("hInvMassLambdaXi_evtpool_deltaphi","",180,0,180,10000,0,10));
  //fOutputList->Add(new TH2F("hInvMassLambdaXi_evtpool_centvsmass","",101,0,101,10000,0,10));
  //===============
  fOutputList->Add(new TH1F("fPt_allLambda","",1000,0,10));
  fOutputList->Add(new TH1F("fPt_xidecaylambda","",1000,0,10));
  fOutputList->Add(new TH1F("fPt_allAntiLambda","",1000,0,10));
  fOutputList->Add(new TH1F("fPt_xidecayAntiLambda","",1000,0,10));
  fOutputList->Add(new TH2F("fKaon_pt_eta","",1000,0,10,200,-1,1));

  //========== Flow analysis ==========
  fOutputList->Add(new TH1F("hEventPlane","",2000,-10,10));

  //========== Raw v2 analysis ==========
  fOutputList->Add(new TH1F("fPsi_pion","",2000,-10,10));
  fOutputList->Add(new TH1F("fPhi_pion","",2000,-10,10));
  fOutputList->Add(new TH1F("fPhi-Psi_pion","",2000,-10,10));

  //========== Eta correction for v2 ==========
  fOutputList->Add(new TH2F("fPhivsEta_lambda","",2000,-10,10,40,-1,1));
  fOutputList->Add(new TH2F("fPhivsEta_xi","",2000,-10,10,40,-1,1));

  //===============
  //============== 20220524 Flow study
  /*
  if(fIsFlowTask){
    AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask 
      = dynamic_cast<AliAnalysisTaskFlowVectorCorrections*>(AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
    if(flowQnVectorTask != NULL){
      fFlowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
    }
    else {
      AliFatal("Flow Qn vector corrections framework needed but it is not present. ABORTING!!!");
    }
  }
  */

  PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskMSDibaryons::UserExec(Option_t *)
{
  //cout<<"++++++++++++++++++++"<<endl;
  //cout<<"++ Analysis Start ++"<<endl;
  //cout<<"++++++++++++++++++++"<<endl;

  // Event loop
  
  // Input event
  //AliAnalysisManager   *man         =AliAnalysisManager::GetAnalysisManager();
  //AliInputEventHandler *inputHandler=(AliInputEventHandler*)(man->GetInputEventHandler());

  fEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if(!fEvent){
    AliError("event is not available.");
    return;
  }
  fESD = dynamic_cast<AliESDEvent*>(fEvent);
  fAOD = dynamic_cast<AliAODEvent*>(fEvent);
 
  if(!fAOD){
    printf("ERROR: fAOD not available\n");
    return;
  }
  fHeader=(AliAODHeader*)fAOD->GetHeader();
  if(!fHeader){
    printf("ERROR: fHeader not available\n");
    return;
  }
  if(fInputHandler){
    fPIDResponse=fInputHandler->GetPIDResponse();
    if(!fPIDResponse){
      printf("ERROR: fPIDResponse not available\n");
      return;
    }
  }
  
  ////=============== Event selection
  Bool_t isEvtMB=((((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()) & ftrigBit);
  if(!isEvtMB) return;

  printf("before Event selection\n");
  dynamic_cast<TH1F*>(fOutputList->FindObject("fEvtCounter"))->Fill(0.5);
  // Event pile-up
  Bool_t isPileupEvt = fAOD->IsPileupFromSPD();
  Bool_t isGoodEvt=EventSelection(fAOD);
  if(!isGoodEvt || isPileupEvt){
    PostData(1,fOutputList);
    return; 
  }
  dynamic_cast<TH1F*>(fOutputList->FindObject("fEvtPassCut"))->Fill(0.5);
 
  Float_t centralityV0M = -1.;
  Float_t centralityCL0 = -1.;
  Float_t centralityCL1 = -1.;
  Float_t centralityV0A = -1.;
  Float_t centralityV0C = -1.;
  Float_t centralityZNA = -1.;
  Float_t centralityZNC = -1.;
  Float_t centralityMain = -1.;

  //Get Centrality
  fMultSelection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
  if(!fMultSelection){
    //If you get this warning (and fCentralityV0M 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
    return;
  }
  else{
    centralityV0M  = fMultSelection->GetMultiplicityPercentile("V0M");
    centralityCL0  = fMultSelection->GetMultiplicityPercentile("CL0");
    centralityCL1  = fMultSelection->GetMultiplicityPercentile("CL1");
    centralityV0A  = fMultSelection->GetMultiplicityPercentile("V0A");
    centralityV0C  = fMultSelection->GetMultiplicityPercentile("V0C");
    centralityZNA  = fMultSelection->GetMultiplicityPercentile("ZNA");
    centralityZNC  = fMultSelection->GetMultiplicityPercentile("ZNC");
    centralityMain = fMultSelection->GetMultiplicityPercentile(fEstimator);
  }

  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityV0M")) ->Fill(centralityV0M);
  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityV0A")) ->Fill(centralityV0A);
  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityV0C")) ->Fill(centralityV0C);
  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityZNA")) ->Fill(centralityZNA);
  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityZNC")) ->Fill(centralityZNC);
  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityCL0")) ->Fill(centralityCL0);
  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityCL1")) ->Fill(centralityCL1);
  dynamic_cast<TH1F*>(fOutputList->FindObject("hCentralityMain"))->Fill(centralityMain);

  printf("before Qn\n");
  printf("Local file confirmed!\n");

  //Bool_t isok = ExtractQnVector();
  //printf("isok = %d\n",isok);

  //============= variable definition
  Int_t    nTracks   =fAOD->GetNumberOfTracks();
  Int_t    nV0s      =fAOD->GetNumberOfV0s();
  Int_t    nCascades =fAOD->GetNumberOfCascades();
  Double_t lamntrack[6][100]    ={{0},{0},{0},{0},{0},{0}};
  Double_t lambtrack[6][100]    ={{0},{0},{0},{0},{0},{0}};
  Double_t Lambda[10000]={0};
  Double_t AntiLambda[10000]={0};

  Float_t  vecTarget[3]={0.};
  Int_t multiplicity =0;
  
  vecTarget[0]=fAOD->GetPrimaryVertex()->GetX(); // primary vertex
  vecTarget[1]=fAOD->GetPrimaryVertex()->GetY();
  vecTarget[2]=fAOD->GetPrimaryVertex()->GetZ();
  
  AliAODTracklets *tracklets=(AliAODTracklets*)fAOD->GetTracklets();
  Int_t nTracklet = tracklets->GetNumberOfTracklets();
  
  for(Int_t itracklet = 0; itracklet<nTracklet; itracklet++){
    Double_t tracklet_eta   = tracklets->GetEta(itracklet);    
    if(fabs(tracklet_eta) > 0.8) continue;
    multiplicity++;    
  }

  //PDG mass                                                                                                        
  const Double_t massProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  const Double_t massLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  const Double_t massXi     = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  const Double_t massOmega  = TDatabasePDG::Instance()->GetParticle(3334)->Mass();

  //============== cascade loop
  TObjArray* fCascadeArraym=new TObjArray();
  TObjArray* fCascadeArrayp=new TObjArray();

  TObjArray* fCascadeArrayOmegam=new TObjArray();
  TObjArray* fCascadeArrayOmegap=new TObjArray();

  Int_t    nXim=0; // Xi-/event
  Int_t    nXip=0; // Xi+/event

  Int_t    nOmegam=0; // Omega-/event
  Int_t    nOmegap=0; // Omega+/event

  for (int iCasc=0; iCasc< nCascades; iCasc++){
    AliAODcascade *casc= (AliAODcascade*)fAOD->GetCascade(iCasc);
    if(!casc) continue;

    AliAODTrack* pTrackXi   = dynamic_cast<AliAODTrack*>(casc->GetDaughter(0));                     // daughter track +
    AliAODTrack* nTrackXi   = dynamic_cast<AliAODTrack*>(casc->GetDaughter(1));                     // daughter track -
    AliAODTrack* bachTrackXi= dynamic_cast<AliAODTrack*>(casc->GetDecayVertexXi()->GetDaughter(0)); // bachlor track
    double_t dcaxyp=0.,dcazp=0.,dcap=0.;  // daughter positive track
    double_t dcaxyn=0.,dcazn=0.,dcan=0.;  // daughter negative track
    double_t dcaxyb=0.,dcazb=0.,dcab=0.;  // bachelor track
    Float_t  dDCAp[2]    ={0.};        
    Float_t  cDCAp[3]    ={0.};        
    Float_t  dDCAn[2]    ={0.}; 
    Float_t  cDCAn[3]    ={0.}; 
    Float_t  dDCAb[2]    ={0.}; 
    Float_t  cDCAb[3]    ={0.}; 
    double_t v0Vtx[3]    ={0.};
    double_t xivertex[3] ={0.};
    Double_t lambdacpa=0.,lambdatransradius=0.;
    Double_t xicpa=0.,xitransradius=0.;
    Double_t v0vtxx=0.,v0vtxy=0.,v0vtxz=0.,xivtxx=0.,xivtxy=0.,xivtxz=0.;

    if((pTrackXi->Charge()) < (nTrackXi->Charge())) continue;
    
    //========== Daughter track selection 
    pTrackXi->GetImpactParameters(dDCAp,cDCAp);
    dcaxyp          =dDCAp[0];
    dcazp           =dDCAp[1];
    dcap            =sqrt(dcaxyp*dcaxyp+dcazp*dcazp);
    nTrackXi->GetImpactParameters(dDCAn,cDCAn);
    dcaxyn          =dDCAn[0];
    dcazn           =dDCAn[1];
    dcan            =sqrt(dcaxyn*dcaxyn+dcazn*dcazn);
    bachTrackXi->GetImpactParameters(dDCAb,cDCAb);
    dcaxyb          =dDCAb[0];
    dcazb           =dDCAb[1];
    dcab            =sqrt(dcaxyb*dcaxyb+dcazb*dcazb);

    //cascade daughter select contents
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdaudca"))->Fill(pTrackXi->Pt());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdaudca"))->Fill(nTrackXi->Pt());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdaudcab"))->Fill(bachTrackXi->Pt());

    if(fabs(pTrackXi   ->Eta()) > 0.8)       continue;
    if(fabs(nTrackXi   ->Eta()) > 0.8)       continue;
    if(fabs(bachTrackXi->Eta()) > 0.8)       continue;
    if(pTrackXi   ->GetTPCNcls() < 70) continue;
    if(nTrackXi   ->GetTPCNcls() < 70) continue;
    if(bachTrackXi->GetTPCNcls() < 70) continue;
    if(pTrackXi   ->Pt() < 0.3)        continue;
    if(nTrackXi   ->Pt() < 0.3)        continue;
    if(bachTrackXi->Pt() < 0.3)        continue;
    if(dcap < 0.05) continue;
    if(dcan < 0.05) continue;
    if(dcab < 0.05) continue;

    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauEta"))->Fill(pTrackXi->Eta());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauEta"))->Fill(nTrackXi->Eta());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauEtab"))->Fill(bachTrackXi->Eta());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauTPC"))->Fill(pTrackXi->GetTPCNcls());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauTPC"))->Fill(nTrackXi->GetTPCNcls());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauTPCb"))->Fill(bachTrackXi->GetTPCNcls());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauPt"))->Fill(pTrackXi->Pt());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauPt"))->Fill(nTrackXi->Pt());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauPtb"))->Fill(bachTrackXi->Pt());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauPIDproton"))->Fill(fabs(fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kProton)));
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauPIDproton"))->Fill(fabs(fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kProton)));
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauPIDpion"))->Fill(fabs(fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kPion)));
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascdauPIDpion"))->Fill(fabs(fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kPion)));
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascbachPIDpion"))->Fill(fabs(fPIDResponse->NumberOfSigmasTPC(bachTrackXi, AliPID::kPion)));

    //========== PID
    Bool_t isProton     =(fabs(fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kProton)) < 4);  // proton
    Bool_t isAntiproton =(fabs(fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kProton)) < 4);  // proton- 
    Bool_t isPospion    =(fabs(fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kPion)) < 4);    // pion+
    Bool_t isNegpion    =(fabs(fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kPion)) < 4);    // pion-
    Bool_t isBachpion   =(fabs(fPIDResponse->NumberOfSigmasTPC(bachTrackXi, AliPID::kPion)) < 4); // bacelor (pion)
    Bool_t isBachkaon   =(fabs(fPIDResponse->NumberOfSigmasTPC(bachTrackXi, AliPID::kKaon)) < 4); // bacelor (kaon)   

    Bool_t isXi         =(isProton     && isNegpion && bachTrackXi->Charge() < 0 && isBachpion);
    Bool_t isXibar      =(isAntiproton && isPospion && bachTrackXi->Charge() > 0 && isBachpion);
    Bool_t isOmega      =(isProton     && isNegpion && bachTrackXi->Charge() < 0 && isBachkaon);
    Bool_t isOmegabar   =(isAntiproton && isPospion && bachTrackXi->Charge() > 0 && isBachkaon);

    //========== lambda selection
    v0vtxx               =casc->DecayVertexV0X();                      
    v0vtxy               =casc->DecayVertexV0Y(); 
    v0vtxz               =casc->DecayVertexV0Z();
    v0Vtx[0]             =v0vtxx;
    v0Vtx[1]             =v0vtxy;
    v0Vtx[2]             =v0vtxz;
    lambdacpa            =LambdaCosPointingAngle(casc,v0Vtx,vecTarget);
    lambdatransradius    =DecayLengthXY(v0Vtx,vecTarget);

    //cascade V0 select contents
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascV0Transradius"))->Fill(lambdatransradius);
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascV0daugdca"))->Fill(casc->DcaV0Daughters());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascV0dca"))->Fill(casc->DcaV0ToPrimVertex());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fcascV0cpa"))->Fill(lambdacpa);

    if(casc->GetOnFlyStatus())           continue;// select offline     
    if(lambdacpa < 0.97)                 continue;

    Bool_t CascV0isXi    = true;
    Bool_t CascV0isOmega = true;

    if(lambdatransradius < 1.4)          CascV0isXi = false;
    if(lambdatransradius > 200)          CascV0isXi = false;
    if(casc->DcaV0Daughters() > 1.5)     CascV0isXi = false;
    if(casc->DcaV0ToPrimVertex() < 0.07) CascV0isXi = false;

    if(lambdatransradius < 1.0)          CascV0isOmega = false;
    if(lambdatransradius > 200)          CascV0isOmega = false;
    if(casc->DcaV0Daughters() > 1.2)     CascV0isOmega = false;
    if(casc->DcaV0ToPrimVertex() < 0.06) CascV0isOmega = false;

    if(CascV0isXi && isXi){
      dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassXidecayLambda"))->Fill(casc->MassLambda());
    }
    if(CascV0isXi && isXibar){
      dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassXidecayLambda"))->Fill(casc->MassAntiLambda());
    }

    //===== xi selection ===== 
    xivtxx        =casc->DecayVertexXiX();                                     
    xivtxy        =casc->DecayVertexXiY();
    xivtxz        =casc->DecayVertexXiZ();
    xivertex[0]   =xivtxx;
    xivertex[1]   =xivtxy;
    xivertex[2]   =xivtxz;
    xicpa         =casc->CosPointingAngleXi(vecTarget[0],vecTarget[1],vecTarget[2]);
    xitransradius =xiDecayLengthXY(xivertex,vecTarget);

    //Xi select contents
    dynamic_cast<TH1F*>(fOutputList->FindObject("XiTranverseradius"))->Fill(xitransradius);
    dynamic_cast<TH1F*>(fOutputList->FindObject("dcaXiDaughters"))->Fill(casc->DcaXiDaughters());
    dynamic_cast<TH1F*>(fOutputList->FindObject("XiCPA"))->Fill(xicpa);
    //dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassOmega"))->Fill(InvMassOmega(casc));

    Bool_t CascisXi    = true;
    Bool_t CascisOmega = true;

    if(xicpa < 0.98)                 CascisXi = false;
    if(xitransradius < 0.8)          CascisXi = false;
    if(xitransradius > 200)          CascisXi = false;
    if(casc->DcaXiDaughters() > 1.6) CascisXi = false;
    if((casc->MassOmega() > 1.667) && (casc->MassOmega() < 1.677)) CascisXi = false;

    if(xicpa < 0.995)                CascisOmega = false;
    if(xitransradius < 0.2)          CascisOmega = false;
    if(xitransradius > 200)          CascisOmega = false;
    if(casc->DcaXiDaughters() > 0.8) CascisOmega = false;
    if((casc->MassXi() > 1.317) && (casc->MassXi() < 1.327)) CascisOmega = false;


    //=========== lambda + pion- -> xi-
    if(isXi && CascV0isXi && CascisXi){

      if(TMath::Abs(casc->MassLambda() - massLambda) > 0.006) continue;
      dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassXi"))->Fill(sqrt(casc->Pt2Xi()),casc->MassXi());

      if(TMath::Abs(casc->MassXi() - massXi) > 0.005) continue;
      //Side band
      //if(TMath::Abs(casc->MassXi() - massXi) < 0.005) continue;
      //if(TMath::Abs(casc->MassXi() - massXi) > 0.01)  continue;
      dynamic_cast<TH2F*>(fOutputList->FindObject("fPtcascade_xi"))->Fill(sqrt(casc->Pt2Xi()),casc->MassXi());
      dynamic_cast<TH1F*>(fOutputList->FindObject("fPtv0casc"))->Fill(sqrt(casc->Pt2V0()));
      dynamic_cast<TH1F*>(fOutputList->FindObject("fPtBachelor"))->Fill(sqrt(pow(casc->MomBachX(),2)+pow(casc->MomBachY(),2)));

      nXim++;
      fCascadeArraym->Add(casc);
    }
    //========== antilambda + pion+ -> xi+
    if(isXibar && CascV0isXi && CascisXi){
     
      if(TMath::Abs(casc->MassAntiLambda() - massLambda) > 0.006) continue;
      dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassXi"))->Fill(sqrt(casc->Pt2Xi()),casc->MassXi());

      if(TMath::Abs(casc->MassXi() - massXi) > 0.005) continue;
      //Side band
      //if(TMath::Abs(casc->MassXi() - massXi) < 0.005) continue;
      //if(TMath::Abs(casc->MassXi() - massXi) > 0.01) continue;
      dynamic_cast<TH2F*>(fOutputList->FindObject("fPtcascade_xi"))->Fill(sqrt(casc->Pt2Xi()),casc->MassXi());
      dynamic_cast<TH1F*>(fOutputList->FindObject("fPtv0casc"))->Fill(sqrt(casc->Pt2V0()));
      dynamic_cast<TH1F*>(fOutputList->FindObject("fPtBachelor"))->Fill(sqrt(pow(casc->MomBachX(),2)+pow(casc->MomBachY(),2)));

      nXip++;
      fCascadeArrayp->Add(casc);
    } 
  }
  //============== cascade loop end

  //========== v0 loop
  TObjArray* fV0Arrayn=new TObjArray();
  TObjArray* fV0Arrayb=new TObjArray();

  int nLamn=0;
  int nLamb=0;

  for (Int_t iV0=0; iV0<nV0s; iV0++){
    AliAODv0 *v0= (AliAODv0*)fAOD->GetV0(iV0);
    if(!v0) continue; 

    AliAODTrack* ptrack=dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));    
    AliAODTrack* ntrack=dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));    
    Double_t dcaxyp=0.,dcazp=0.,dcap=0.; // proton
    Double_t dcaxyn=0.,dcazn=0.,dcan=0.; // pion
    Float_t  dDCAp[2]={0.};              // DCA to the vertex xy(d) and z
    Float_t  dDCAn[2]={0.};  
    Float_t  cDCAp[3]={0.};              // convariance of impact parameters
    Float_t  cDCAn[3]={0.}; 
    Double_t v0Vtx[3]={0.};
    Double_t energy=0.; 
    Double_t fcpa=0.,ftransradius=0.;

    if((ptrack->Charge()) < (ntrack->Charge())) continue;  

    ptrack->GetImpactParameters(dDCAp,cDCAp);
    dcaxyp          =dDCAp[0];
    dcazp           =dDCAp[1];
    dcap            =sqrt(dcaxyp*dcaxyp+dcazp*dcazp);
    ntrack->GetImpactParameters(dDCAn,cDCAn);
    dcaxyn          =dDCAn[0];
    dcazn           =dDCAn[1];
    dcan            =sqrt(dcaxyn*dcaxyn+dcazn*dcazn);

    //daughter track select contents
    dynamic_cast<TH1F*>(fOutputList->FindObject("fDCA"))->Fill(dcap);
    dynamic_cast<TH1F*>(fOutputList->FindObject("fDCA"))->Fill(dcan);

    //========== v0 daughter tracks selection
    if(fabs(ptrack->Eta()) > 0.8) continue;
    if(fabs(ntrack->Eta()) > 0.8) continue;
    if(ptrack->GetTPCNcls() < 70) continue;
    if(ntrack->GetTPCNcls() < 70) continue;
    if(dcap < 0.05) continue;
    if(dcan < 0.05) continue;

    dynamic_cast<TH1F*>(fOutputList->FindObject("fEta"))->Fill(ptrack->Eta());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fEta"))->Fill(ntrack->Eta());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fTPC"))->Fill(ptrack->GetTPCNcls());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fTPC"))->Fill(ntrack->GetTPCNcls());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fPIDproton"))->Fill(fabs(fPIDResponse->NumberOfSigmasTPC(ptrack, AliPID::kProton)));
    dynamic_cast<TH1F*>(fOutputList->FindObject("fPIDproton"))->Fill(fabs(fPIDResponse->NumberOfSigmasTPC(ntrack, AliPID::kProton)));
    dynamic_cast<TH1F*>(fOutputList->FindObject("fPIDpion"))->Fill(fabs(fPIDResponse->NumberOfSigmasTPC(ptrack, AliPID::kPion)));
    dynamic_cast<TH1F*>(fOutputList->FindObject("fPIDpion"))->Fill(fabs(fPIDResponse->NumberOfSigmasTPC(ntrack, AliPID::kPion)));

    //========== Out-of-bunch pile-up removal
    if(!(ptrack->HasPointOnITSLayer(0)) && !(ptrack->HasPointOnITSLayer(1))&&
       !(ptrack->HasPointOnITSLayer(4)) && !(ptrack->HasPointOnITSLayer(5))&&
       !(ptrack->GetTOFBunchCrossing()==0)) continue;
    if(!(ntrack->HasPointOnITSLayer(0)) && !(ntrack->HasPointOnITSLayer(1))&&
       !(ntrack->HasPointOnITSLayer(4)) && !(ntrack->HasPointOnITSLayer(5))&&
       !(ntrack->GetTOFBunchCrossing()==0)) continue;

    //========== daughter track PID
    Bool_t Daugproton     =(fabs(fPIDResponse->NumberOfSigmasTPC(ptrack, AliPID::kProton)) < 5); // proton +
    Bool_t Daugantiproton =(fabs(fPIDResponse->NumberOfSigmasTPC(ntrack, AliPID::kProton)) < 5); // proton -
    Bool_t Daugpion       =(fabs(fPIDResponse->NumberOfSigmasTPC(ntrack, AliPID::kPion)) < 5);   // pion -
    Bool_t Daugnpion      =(fabs(fPIDResponse->NumberOfSigmasTPC(ptrack, AliPID::kPion)) < 5);   // pion +
    Bool_t lambda         =(Daugproton && Daugpion);
    Bool_t antilambda     =(Daugnpion && Daugantiproton);

    v0Vtx[0]     =v0->DecayVertexV0X();
    v0Vtx[1]     =v0->DecayVertexV0Y();
    v0Vtx[2]     =v0->DecayVertexV0Z();
    fcpa         =CosPointingAngle(v0,v0Vtx,vecTarget);
    ftransradius =DecayLengthXY(v0Vtx,vecTarget);

    //lambda select contents
    dynamic_cast<TH1F*>(fOutputList->FindObject("fV0pt"))->Fill(sqrt(v0->Pt2V0()));
    dynamic_cast<TH1F*>(fOutputList->FindObject("fV0vtxx"))->Fill(v0Vtx[0]);
    dynamic_cast<TH1F*>(fOutputList->FindObject("fV0vtxy"))->Fill(v0Vtx[1]);
    dynamic_cast<TH1F*>(fOutputList->FindObject("fV0vtxz"))->Fill(v0Vtx[2]);
    dynamic_cast<TH1F*>(fOutputList->FindObject("fTransRadius"))->Fill(ftransradius);
    dynamic_cast<TH1F*>(fOutputList->FindObject("fDCAdaugTov0Vtx"))->Fill(v0->DcaV0Daughters());
    dynamic_cast<TH1F*>(fOutputList->FindObject("fCPA"))->Fill(fcpa);
    dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassK0s"))->Fill(InvMassK0(v0));

    //========== v0 selection
    if(sqrt(v0->Pt2V0()) < 0.3)    continue;
    if(fabs(v0Vtx[0]) > 100)       continue;
    if(fabs(v0Vtx[1]) > 100)       continue;
    if(fabs(v0Vtx[2]) > 100)       continue;
    if(ftransradius < 0.2)         continue;
    if(ftransradius > 100)         continue;
    if(v0->DcaV0Daughters() > 1.5) continue;
    if(fcpa < 0.99)                continue;
    if((v0->MassK0Short() > 0.48) && (v0->MassK0Short() < 0.515)) continue;

    if(v0->GetOnFlyStatus()) continue; // select offline v0

    //========== lambda invariant mass
    if(lambda){

      dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambda"))->Fill(sqrt(v0->Pt2V0()),v0->MassLambda());
      if(TMath::Abs(v0->MassLambda() - massLambda) > 0.004) continue;
      //Side band
      //if(TMath::Abs(v0->MassLambda() - massLambda) < 0.004) continue;
      //if(TMath::Abs(v0->MassLambda() - massLambda) > 0.008) continue;
      dynamic_cast<TH2F*>(fOutputList->FindObject("fPtv0_lambda"))->Fill(sqrt(v0->Pt2V0()),v0->MassLambda());

      nLamn++;

      fV0Arrayn->Add(v0);
    }

    //========== antilambda invariant mass
    if(antilambda){    

      dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambda"))->Fill(sqrt(v0->Pt2V0()),v0->MassAntiLambda());

      if(TMath::Abs(v0->MassAntiLambda() - massLambda) > 0.004) continue;
      //Side band
      //if(TMath::Abs(v0->MassAntiLambda() - massLambda) < 0.004) continue;
      //if(TMath::Abs(v0->MassAntiLambda() - massLambda) > 0.008) continue;
      dynamic_cast<TH2F*>(fOutputList->FindObject("fPtv0_lambda"))->Fill(sqrt(v0->Pt2V0()),v0->MassAntiLambda());

      nLamb++; 
     
      fV0Arrayb->Add(v0);
    }    
  }
  //=============== V0 loop end
  
  //=============== Track loop
  TObjArray* fProtonArray=new TObjArray();
  TObjArray* fProtonArrayb=new TObjArray();

  int np =0;
  int npb=0;

  Int_t pid=0;
  Int_t antipid=0;

  for(Int_t i=0; i<nTracks; i++){
    AliAODTrack *track=static_cast<AliAODTrack*>(fAOD->GetTrack(i));
    if(!track) continue;

    Double_t tpcratio=0.,dcaxy=0.,dcaz=0.,tofsig=0.,trklength=0.,beta=0.;
    UShort_t crsR=0,clsF=0; 
    Float_t  sigmatpc_pr=0.,sigmatpc_ka=0.,sigmatpc_pion=0.,sigmatof=0.,sigmacomb=0.;
    Float_t  dDCA[2]    ={0.};        
    Float_t  cDCA[3]    ={0.};        

    //========== proton selection 
    crsR         =track->GetTPCNCrossedRows();
    clsF         =track->GetTPCNclsF();
    tpcratio     =(clsF>0)?(double)crsR/clsF:0;
    track->GetImpactParameters(dDCA,cDCA); 
    dcaxy        =dDCA[0];
    dcaz         =dDCA[1];
    sigmatpc_pr  =fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    sigmatpc_ka  =fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    sigmatpc_pion=fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    sigmatof     =fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
    sigmacomb    =sqrt(pow(sigmatpc_pr,2)+pow(sigmatof,2));    
    tofsig       =track->GetTOFsignal();
    trklength    =track->GetIntegratedLength();
    beta         =(tofsig>0)?(double)trklength/ (2.99792457999999984e-02 * tofsig):0;

    if(fabs(track->Eta()) > 0.8) continue;
    if(track->Pt() < 0.5) continue;
    if(track->Pt() > 4.05) continue;
    if(track->GetTPCNcls() < 80) continue;
    if(crsR < 70) continue;
    if(tpcratio < 0.83) continue;
    if(!(track->GetTPCSharedMap()==0)) continue;
    if(fabs(dcaxy) > 0.1) continue;
    if(fabs(dcaz) > 0.2) continue;
    //Bool_t isProton  =((track->P() < 0.75 && fabs(sigmatpc_pr) < 3) || (track->P() > 0.75 && sigmacomb < 3));
    Bool_t isProton  =((track->P() < 0.75 && fabs(sigmatpc_pr) < 3) || (track->P() > 0.75 && fabs(sigmatpc_pr) < 3 && fabs(sigmatof) < 3));
    Bool_t isKaon    =fabs(sigmatpc_ka)<2;
    Bool_t isPion    =fabs(sigmatpc_pion)<3;
    //Bool_t isProton  =((track->P() < 0.75 && fabs(sigmatpc_pr) < 3) || (track->P() > 0.75 && ((fabs(sigmatpc_pr) < 3)||(fabs(sigmatof) < 3)) ));
    if(isKaon){
      dynamic_cast<TH2F*>(fOutputList->FindObject("fKaon_pt_eta"))->Fill(track->Pt(),track->Eta());
    }

    if(!isProton) continue;

    //proton select contents
    dynamic_cast<TH1F*>(fOutputList->FindObject("fTrkEta"))->Fill(track->Eta()); 
    dynamic_cast<TH1F*>(fOutputList->FindObject("fTrkPt"))->Fill(track->Pt()); 
    dynamic_cast<TH1F*>(fOutputList->FindObject("fTrkTPCclus"))->Fill(track->GetTPCNcls()); 
    dynamic_cast<TH1F*>(fOutputList->FindObject("fTrkTPCcrsR"))->Fill(crsR); 
    dynamic_cast<TH1F*>(fOutputList->FindObject("fTrkTPCclusF"))->Fill(tpcratio); 
    dynamic_cast<TH1F*>(fOutputList->FindObject("fTrkDCAxy"))->Fill(dcaxy); 
    dynamic_cast<TH1F*>(fOutputList->FindObject("fTrkDCAz"))->Fill(dcaz); 
    dynamic_cast<TH2F*>(fOutputList->FindObject("fPvsTPCsignal"))->Fill(track->P(),track->GetTPCsignal()); 
    dynamic_cast<TH2F*>(fOutputList->FindObject("fPvsTOFsignal"))->Fill(track->P(),beta); 

    //========== proton
    if(track->Charge() > 0){

      dynamic_cast<TH1F*>(fOutputList->FindObject("fPt_proton"))->Fill(track->Pt()); 
      dynamic_cast<TH2F*>(fOutputList->FindObject("PID_cut"))->Fill(track->P(),track->GetTPCsignal()); 
      np++;
      fProtonArray->Add(track);
    }

    //========== anti proton
    if(track->Charge() < 0){    

      dynamic_cast<TH1F*>(fOutputList->FindObject("fPt_proton"))->Fill(track->Pt()); 
      dynamic_cast<TH2F*>(fOutputList->FindObject("PID_cut"))->Fill(track->P(),track->GetTPCsignal()); 
      npb++;
      fProtonArrayb->Add(track);
    }

  }
  //=================== track loop end

  // remove same track id 
  for(Int_t i=0; i < fProtonArray->GetEntriesFast(); i++){
    TObject *objProton1=fProtonArray->At(i);    
    if(!objProton1) continue;
    AliAODTrack *proton1=dynamic_cast<AliAODTrack*>(objProton1); 
    Int_t id1 =proton1->GetID();
    if(id1 < 0) id1 =-id1-1;

    for(Int_t j=i+1; j < fProtonArray->GetEntriesFast();){
      TObject *objProton2=fProtonArray->At(j);    
      if(!objProton2) { j++; continue;}
      AliAODTrack *proton2=dynamic_cast<AliAODTrack*>(objProton2); 
      Int_t id2 =proton2->GetID();
      if(id2 < 0) id2 =-id2-1;

      if(id1==id2) fProtonArray->RemoveAt(j);
      else j++;
    }
  }
 
  for(Int_t i=0; i < fProtonArrayb->GetEntriesFast(); i++){
    TObject *objProton1=fProtonArrayb->At(i);    
    if(!objProton1) continue;
    AliAODTrack *proton1=dynamic_cast<AliAODTrack*>(objProton1); 
    Int_t id1 =proton1->GetID();
    if(id1 < 0) id1 =-id1-1;

    for(Int_t j=i+1; j < fProtonArrayb->GetEntriesFast();){
      TObject *objProton2=fProtonArrayb->At(j);    
      if(!objProton2) { j++; continue;}
      AliAODTrack *proton2=dynamic_cast<AliAODTrack*>(objProton2); 
      Int_t id2 =proton2->GetID();
      if(id2 < 0) id2 =-id2-1;

      if(id1==id2) fProtonArrayb->RemoveAt(j);
      else j++;
    }
  }

  // decay rejection
  // Lambda
  if(nLamn>0){
    for(Int_t i=0; i < fV0Arrayn->GetEntriesFast(); i++){
      TObject *objV0=fV0Arrayn->At(i);    
      if(!objV0) continue;
      AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objV0); 

      if(nXim > 0){
	for(Int_t j=0; j < fCascadeArraym->GetEntriesFast(); j++){
	  TObject *objCascade=fCascadeArraym->At(j);    
	  if(!objCascade) continue;
	  AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 
	  Bool_t isXidecayLambda = ((lambda->GetPosID()) == xi->GetPosID()) && ((lambda->GetNegID()) == (xi->GetNegID()));
	  // xi decay lambda → 1 + 0 + 0 = 1, prompt lambda → 0 + 0 + 0 = 0
	  Lambda[i] += isXidecayLambda;
	}
      }
      if(nXim == 0){
        Lambda[i] = 0;
      }
    }
  }

  // Lambda pT
  if(nLamn > 0){
    for(Int_t i=0; i < fV0Arrayn->GetEntriesFast(); i++){
      TObject *objV0=fV0Arrayn->At(i);    
      if(!objV0) continue;
      AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objV0); 

      dynamic_cast<TH1F*>(fOutputList->FindObject("fPt_allLambda"))->Fill(sqrt(pow(lambda->MomV0X(),2)+pow(lambda->MomV0Y(),2))); 
      
      if(Lambda[i]==1){
	dynamic_cast<TH1F*>(fOutputList->FindObject("fPt_xidecaylambda"))->Fill(sqrt(pow(lambda->MomV0X(),2)+pow(lambda->MomV0Y(),2)));
      }
    }
  }

  //antiLambda
  if(nLamb>0){
    for(Int_t i=0; i < fV0Arrayb->GetEntriesFast(); i++){
      TObject *objV0=fV0Arrayb->At(i);    
      if(!objV0) continue;
      AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objV0); 
      if(nXip > 0){
	for(Int_t j=0; j < fCascadeArrayp->GetEntriesFast(); j++){
	  TObject *objCascade=fCascadeArrayp->At(j);    
	  if(!objCascade) continue;
	  AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 
	    
	  Bool_t isXidecayAntiLambda = ((lambda->GetPosID()) == xi->GetPosID()) && ((lambda->GetNegID()) == (xi->GetNegID()));
	  AntiLambda[i] += isXidecayAntiLambda;
	}
      }
      if(nXip == 0){
        AntiLambda[i] = 0;
      }

    }
  }

  // antiLambda pT
  if(nLamb > 0){
    for(Int_t i=0; i < fV0Arrayb->GetEntriesFast(); i++){
      TObject *objV0=fV0Arrayb->At(i);    
      if(!objV0) continue;
      AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objV0); 

      dynamic_cast<TH1F*>(fOutputList->FindObject("fPt_allAntiLambda"))->Fill(sqrt(pow(lambda->MomV0X(),2)+pow(lambda->MomV0Y(),2)));

      if(AntiLambda[i]==1){
	dynamic_cast<TH1F*>(fOutputList->FindObject("fPt_xidecayAntiLambda"))->Fill(sqrt(pow(lambda->MomV0X(),2)+pow(lambda->MomV0Y(),2))); 
      }
    }
  }

  //==================
  Double_t energysum=0.,energy1=0.,energy2=0.,pt=0.,l1px=0.,l1py=0.,l1pz=0.,l2px=0.,l2py=0.,l2pz=0.,
    invMass=0.,relmom=0.,phi1=0.,phi2=0.,deltaphi=0.;

  const double radian_to_degree=180./TMath::Pi();

  Double_t isXi1=0.,isXi2=0.;
  Short_t posID1=0,posID2=0,negID1=0,negID2=0;
  Short_t posID=0,negID=0;

  Double_t energyxi=0.,energypr=0.,energylam=0.,prpx=0.,prpy=0.,prpz=0.,xipx=0.,xipy=0.,xipz=0.,lampx=0.,lampy=0.,lampz=0.;
  Double_t omegapx=0.,omegapy=0.,omegapz=0.,energyomega=0.;
  Short_t bachid=0,ptrackid=0,ntrackid=0,protonid=0;

  Double_t protonpt1=0.,protonpt2=0.,xidecaypt=0.,protonpt=0.;

  Bool_t islambdalambda=false;
  Bool_t isantilambdaantilambda=false;
  Bool_t isprotonxi=false;
  Bool_t isantiprotonxi=false;
  Bool_t islambdaxi=false; 
  Bool_t isantilambdaxi=false; 
  Bool_t isprotonomega=false;
  Bool_t isantiprotonomega=false;

  // p + Xi-
  if((np > 0) && (nXim > 0)){
    for(Int_t i=0; i < fCascadeArraym->GetEntriesFast(); i++){
      TObject *objCascade=fCascadeArraym->At(i);    
      if(!objCascade) continue;
      AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 
      AliAODTrack* ptrkxi   = dynamic_cast<AliAODTrack*>(xi->GetDaughter(0)); 

      for(Int_t j=0; j < fProtonArray->GetEntriesFast(); j++){
	TObject *objProton=fProtonArray->At(j);    
	if(!objProton) continue;
	AliAODTrack *proton=dynamic_cast<AliAODTrack*>(objProton); 

	xidecaypt   =ptrkxi->Pt();
	protonpt    =proton->Pt();

	xipx        =xi->MomXiX();
	xipy        =xi->MomXiY();
	xipz        =xi->MomXiZ();
	energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));
	bachid      =xi->GetBachID();
	ptrackid    =xi->GetPosID();
	ntrackid    =xi->GetNegID();
	prpx        =proton->Px();
	prpy        =proton->Py();
	prpz        =proton->Pz();
	phi1        =atan2(prpy,prpx)*radian_to_degree;
	phi2        =atan2(xipy,xipx)*radian_to_degree;
	energypr    =sqrt(pow(0.938,2)+pow(prpx,2)+pow(prpy,2)+pow(prpz,2));
	protonid    =proton->GetID();
	if(protonid < 0) protonid =-protonid-1;
	energysum   =energyxi+energypr;

	deltaphi  =fabs(phi1-phi2);
	if(deltaphi > 180){
	  deltaphi = 360-deltaphi;
	}
	relmom    =sqrt(pow((xipx-prpx),2)+pow((xipy-prpy),2)+pow((xipz-prpz),2))/2;
	pt        =sqrt((prpx+xipx)*(prpx+xipx)+(prpy+xipy)*(prpy+xipy));
	invMass   =InvariantMass(xipx,xipy,xipz,prpx,prpy,prpz,energysum);

	if(protonid == ptrackid) continue;
	if(protonid == ntrackid) continue;
	if(protonid == bachid) continue;
	if(bachid == ptrackid) continue;
	if(bachid == ntrackid) continue;
	if(ptrackid == ntrackid) continue;

	//cout<<"---------- FG pxi-"<<j<<" ----------"<<endl;
	//cout<<"(prpx,prpy,prpz,xipx,xipy,xipz) = "<<"("<<prpx<<","<<prpy<<","<<prpz<<","<<xipx<<","<<xipy<<","<<xipz<<")"<<endl;
	dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassProtonXi"))->Fill(invMass); 
	dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassProtonXivspT"))->Fill(pt,invMass); 
	//dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassProtonXi_centvsmass"))->Fill(centralityV0M,invMass); 
	dynamic_cast<TH1F*>(fOutputList->FindObject("fmultiplicity_protonxi"))->Fill(multiplicity); 

	isprotonxi=true; 

      }
    }
  }

  // pbar + Xi+
  if((npb > 0) && (nXip > 0)){
    for(Int_t i=0; i < fCascadeArrayp->GetEntriesFast(); i++){
      TObject *objCascade=fCascadeArrayp->At(i);    
      if(!objCascade) continue;
      AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 

      for(Int_t j=0; j < fProtonArrayb->GetEntriesFast(); j++){
	TObject *objProton=fProtonArrayb->At(j);    
	if(!objProton) continue;
	AliAODTrack *proton=dynamic_cast<AliAODTrack*>(objProton); 

	xipx        =xi->MomXiX();
	xipy        =xi->MomXiY();
	xipz        =xi->MomXiZ();
	energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));
	bachid      =xi->GetBachID();
	ptrackid    =xi->GetPosID();
	ntrackid    =xi->GetNegID();
	prpx        =proton->Px();
	prpy        =proton->Py();
	prpz        =proton->Pz();
	phi1        =atan2(prpy,prpx)*radian_to_degree;
	phi2        =atan2(xipy,xipx)*radian_to_degree;
	energypr    =sqrt(pow(0.938,2)+pow(prpx,2)+pow(prpy,2)+pow(prpz,2));
	protonid    =proton->GetID();
	if(protonid < 0) protonid =-protonid-1;
	energysum   =energyxi+energypr;

	deltaphi  =fabs(phi1-phi2);
	if(deltaphi > 180){
	  deltaphi = 360-deltaphi;
	}
	relmom    =sqrt(pow((xipx-prpx),2)+pow((xipy-prpy),2)+pow((xipz-prpz),2))/2;
	pt        =sqrt((prpx+xipx)*(prpx+xipx)+(prpy+xipy)*(prpy+xipy));
	invMass   =InvariantMass(xipx,xipy,xipz,prpx,prpy,prpz,energysum);

	if(protonid == ptrackid) continue;
	if(protonid == ntrackid) continue;
	if(protonid == bachid) continue;
	if(bachid == ptrackid) continue;
	if(bachid == ntrackid) continue;
	if(ptrackid == ntrackid) continue;

	dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassProtonXi"))->Fill(invMass); 
	dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassProtonXivspT"))->Fill(pt,invMass); 
	//dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassProtonXi_centvsmass"))->Fill(centralityV0M,invMass); 

	isantiprotonxi=true; 

      }
    }
  }

  //lambda+lambda
  if(nLamn>0){
    for(Int_t i=0; i < fV0Arrayn->GetEntriesFast()-1; i++){
      TObject *objV01=fV0Arrayn->At(i);    
      if(!objV01) continue;
      AliAODv0 *lambda1=dynamic_cast<AliAODv0*>(objV01); 
      AliAODTrack* ptrack1=dynamic_cast<AliAODTrack*>(lambda1->GetDaughter(0));    
      for(Int_t j=i+1; j < fV0Arrayn->GetEntriesFast(); j++){
	TObject *objV02=fV0Arrayn->At(j);    
	if(!objV02) continue;
	AliAODv0 *lambda2=dynamic_cast<AliAODv0*>(objV02); 
	AliAODTrack* ptrack2=dynamic_cast<AliAODTrack*>(lambda2->GetDaughter(0));    

	protonpt1 = ptrack1->Pt();
	protonpt2 = ptrack2->Pt();

	l1px     =lambda1->MomV0X();
	l1py     =lambda1->MomV0Y();
	l1pz     =lambda1->MomV0Z();
	l2px     =lambda2->MomV0X();
	l2py     =lambda2->MomV0Y();
	l2pz     =lambda2->MomV0Z();
	phi1     =atan2(l1py,l1px)*radian_to_degree;
	phi2     =atan2(l2py,l2px)*radian_to_degree;
	energy1  =sqrt(pow(1.115,2)+pow(l1px,2)+pow(l1py,2)+pow(l1pz,2));  
	energy2  =sqrt(pow(1.115,2)+pow(l2px,2)+pow(l2py,2)+pow(l2pz,2));  
	posID1   =lambda1->GetPosID();
	negID1   =lambda1->GetNegID();
	posID2   =lambda2->GetPosID();
	negID2   =lambda2->GetNegID();
	energysum=energy1+energy2;

	deltaphi  =fabs(phi1-phi2);
	if(deltaphi > 180){
	  deltaphi = 360-deltaphi;
	}
	pt        =sqrt(pow(l1px+l2px,2)+pow(l1py+l2py,2));
	relmom    =sqrt(pow((l1px-l2px),2)+pow((l1py-l2py),2)+pow((l1pz-l2pz),2))/2;
	invMass   =InvariantMass(l1px,l1py,l1pz,l2px,l2py,l2pz,energysum);
	isXi1     =Lambda[i];
	isXi2     =Lambda[j];

	if(posID1==posID2) continue;  
        if(negID1==negID2) continue;    
	
	// all lambda + all lambda
	dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassLambdaLambda"))->Fill(invMass); 
	dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaLambdavspT"))->Fill(pt,invMass); 
	//dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaLambda_centvsmass"))->Fill(centralityV0M,invMass); 
	dynamic_cast<TH1F*>(fOutputList->FindObject("fmultiplicity_lambdalambda"))->Fill(multiplicity); 

	// cout<<"nLam = "<<nLamn<<endl;
	// cout<<"--------------- FG lambdalambda ---------------"<<endl;
	// cout<<"(l1px,l1py,l1pz,l2px,l2py,l2pz) = "<<"("<<l1px<<","<<l1py<<","<<l1pz<<","<<l2px<<","<<l2py<<","<<l2pz<<")"<<endl;
	islambdalambda=true; 

	// prompt lambda + prompt lambda
	if((isXi1==0) && (isXi2==0)){
	  dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaLambda_onlyprompt"))->Fill(pt,invMass); 
	}
      }
    }
  }

  // antilambda + antilambda
  if(nLamb>0){
    for(Int_t i=0; i < fV0Arrayb->GetEntriesFast()-1; i++){
      TObject *objV01=fV0Arrayb->At(i);    
      if(!objV01) continue;
      AliAODv0 *lambda1=dynamic_cast<AliAODv0*>(objV01); 
      for(Int_t j=i+1; j < fV0Arrayb->GetEntriesFast(); j++){
	TObject *objV02=fV0Arrayb->At(j);    
	if(!objV02) continue;
	AliAODv0 *lambda2=dynamic_cast<AliAODv0*>(objV02); 
	
	l1px     =lambda1->MomV0X();
	l1py     =lambda1->MomV0Y();
	l1pz     =lambda1->MomV0Z();
	l2px     =lambda2->MomV0X();
	l2py     =lambda2->MomV0Y();
	l2pz     =lambda2->MomV0Z();
	phi1     =atan2(l1py,l1px)*radian_to_degree;
        phi2     =atan2(l2py,l2px)*radian_to_degree;
	energy1  =sqrt(pow(1.115,2)+pow(l1px,2)+pow(l1py,2)+pow(l1pz,2));  
	energy2  =sqrt(pow(1.115,2)+pow(l2px,2)+pow(l2py,2)+pow(l2pz,2));  
	posID1   =lambda1->GetPosID();
	negID1   =lambda1->GetNegID();
	posID2   =lambda2->GetPosID();
	negID2   =lambda2->GetNegID();
	energysum=energy1+energy2;

	deltaphi  =fabs(phi1-phi2);
	if(deltaphi > 180){
	  deltaphi = 360-deltaphi;
	}
	pt        =sqrt(pow(l1px+l2px,2)+pow(l1py+l2py,2));
	relmom    =sqrt(pow((l1px-l2px),2)+pow((l1py-l2py),2)+pow((l1pz-l2pz),2))/2;
	invMass   =InvariantMass(l1px,l1py,l1pz,l2px,l2py,l2pz,energysum);
	isXi1     =AntiLambda[i];
	isXi2     =AntiLambda[j];

	if(posID1==posID2) continue;  
        if(negID1==negID2) continue;    
	
	// all lambda + all lambda
	dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassLambdaLambda"))->Fill(invMass); 
	dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaLambdavspT"))->Fill(pt,invMass); 
	//dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaLambda_centvsmass"))->Fill(centralityV0M,invMass); 
	dynamic_cast<TH1F*>(fOutputList->FindObject("fmultiplicity_lambdalambda"))->Fill(multiplicity); 

	isantilambdaantilambda=true; 

	// prompt lambda + prompt lambda
	if((isXi1==0) && (isXi2==0)){
	  dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaLambda_onlyprompt"))->Fill(pt,invMass); 
	}
      }
    }
  }

  // lambda + xi-
  if((nLamn > 0) && (nXim > 0)){
    for(Int_t i=0; i < fCascadeArraym->GetEntriesFast(); i++){
      TObject *objCascade=fCascadeArraym->At(i);    
      if(!objCascade) continue;
      AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 
      AliAODTrack *ptrkxi=dynamic_cast<AliAODTrack*>(xi->GetDaughter(0)); 

      for(Int_t j=0; j < fV0Arrayn->GetEntriesFast(); j++){
	TObject *objV0=fV0Arrayn->At(j);    
	if(!objV0) continue;
	AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objV0); 
	AliAODTrack *ptrack=dynamic_cast<AliAODTrack*>(lambda->GetDaughter(0));    

	xidecaypt   =ptrkxi->Pt();
	protonpt1   =ptrack->Pt();

	xipx        =xi->MomXiX();
	xipy        =xi->MomXiY();
	xipz        =xi->MomXiZ();
	energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));
	bachid      =xi->GetBachID();
	ptrackid    =xi->GetPosID();
	ntrackid    =xi->GetNegID();
	lampx       =lambda->MomV0X();
	lampy       =lambda->MomV0Y();
	lampz       =lambda->MomV0Z();
	phi1        =atan2(lampy,lampx)*radian_to_degree;
        phi2        =atan2(xipy,xipx)*radian_to_degree;
	posID       =lambda->GetPosID();
	negID       =lambda->GetNegID();
	energylam   =sqrt(pow(1.115,2)+pow(lampx,2)+pow(lampy,2)+pow(lampz,2));  
	energysum   =energyxi+energylam;

	deltaphi  =fabs(phi1-phi2);
	if(deltaphi > 180){
	  deltaphi = 360-deltaphi;
	}
	relmom    =sqrt(pow((xipx-lampx),2)+pow((xipy-lampy),2)+pow((xipz-lampz),2))/2;
	pt          =sqrt(pow((lampx+xipx),2)+pow((lampy+xipy),2));
	invMass     =InvariantMass(xipx,xipy,xipz,lampx,lampy,lampz,energysum);

	if(posID == negID) continue;
	if(posID == ptrackid) continue;
	if(posID == ntrackid) continue;
	if(posID == bachid) continue;
	if(negID == ptrackid) continue;
	if(negID == ntrackid) continue;
	if(negID == bachid) continue;
	if(ptrackid == ntrackid) continue;
	if(ptrackid == bachid) continue;
	if(ntrackid == bachid) continue;

	dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassLambdaXi"))->Fill(invMass);
	dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaXivspT"))->Fill(pt,invMass);
	//dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaXi_centvsmass"))->Fill(centralityV0M,invMass); 
	dynamic_cast<TH1F*>(fOutputList->FindObject("fmultiplicity_lambdaxi"))->Fill(multiplicity); 
	islambdaxi=true; 

      }
    }
  }

  // antilambda + xi+
  if((nLamb > 0) && (nXip > 0)){
    for(Int_t i=0; i < fCascadeArrayp->GetEntriesFast(); i++){
      TObject *objCascade=fCascadeArrayp->At(i);    
      if(!objCascade) continue;
      AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 

      for(Int_t j=0; j < fV0Arrayb->GetEntriesFast(); j++){
	TObject *objV0=fV0Arrayb->At(j);    
	if(!objV0) continue;
	AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objV0); 

	xipx        =xi->MomXiX();
	xipy        =xi->MomXiY();
	xipz        =xi->MomXiZ();
	energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));
	bachid      =xi->GetBachID();
	ptrackid    =xi->GetPosID();
	ntrackid    =xi->GetNegID();
	lampx       =lambda->MomV0X();
	lampy       =lambda->MomV0Y();
	lampz       =lambda->MomV0Z();
	phi1        =atan2(lampy,lampx)*radian_to_degree;
        phi2        =atan2(xipy,xipx)*radian_to_degree;
	posID       =lambda->GetPosID();
	negID       =lambda->GetNegID();
	energylam   =sqrt(pow(1.115,2)+pow(lampx,2)+pow(lampy,2)+pow(lampz,2));  
	energysum   =energyxi+energylam;

	deltaphi  =fabs(phi1-phi2);
	if(deltaphi > 180){
	  deltaphi = 360-deltaphi;
	}
	relmom    =sqrt(pow((xipx-lampx),2)+pow((xipy-lampy),2)+pow((xipz-lampz),2))/2;
	pt          =sqrt(pow((lampx+xipx),2)+pow((lampy+xipy),2));
	invMass     =InvariantMass(xipx,xipy,xipz,lampx,lampy,lampz,energysum);

	if(posID == negID) continue;
	if(posID == ptrackid) continue;
	if(posID == ntrackid) continue;
	if(posID == bachid) continue;
	if(negID == ptrackid) continue;
	if(negID == ntrackid) continue;
	if(negID == bachid) continue;
	if(ptrackid == ntrackid) continue;
	if(ptrackid == bachid) continue;
	if(ntrackid == bachid) continue;

	dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassLambdaXi"))->Fill(invMass); 
	dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaXivspT"))->Fill(pt,invMass); 
	//dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaXi_centvsmass"))->Fill(centralityV0M,invMass); 
	dynamic_cast<TH1F*>(fOutputList->FindObject("fmultiplicity_lambdaxi"))->Fill(multiplicity); 

	isantilambdaxi=true; 

      }
    }
  }

  //========== eventmixing with AliEventpool                                                                                             
  AliEventPool *lambdapool =fPoolManagerlambda->GetEventPool(centralityV0M,vecTarget[2]);                                          
  if(!lambdapool) return;                                                                                                          
  AliEventPool *antilambdapool =fPoolManagerantilambda->GetEventPool(centralityV0M,vecTarget[2]);                                        
  if(!antilambdapool) return;                                                                                                                
  AliEventPool *xipool =fPoolManagerxi->GetEventPool(centralityV0M,vecTarget[2]);                                                         
  if(!xipool) return;                                                                                                                 
  AliEventPool *xippool =fPoolManagerxip->GetEventPool(centralityV0M,vecTarget[2]);                                                  
  if(!xippool) return;                                                                                                                    
  AliEventPool *protonpool =fPoolManagerproton->GetEventPool(centralityV0M,vecTarget[2]);                                           
  if(!protonpool) return;                                                                                                                
  AliEventPool *antiprotonpool =fPoolManagerantiproton->GetEventPool(centralityV0M,vecTarget[2]);                                          
  if(!antiprotonpool) return;                                                                                                        
  AliEventPool *lampool =fPoolManagerlam->GetEventPool(centralityV0M,vecTarget[2]);                                                          
  if(!lampool) return;                                                                                                                       
  AliEventPool *antilampool =fPoolManagerantilam->GetEventPool(centralityV0M,vecTarget[2]);                                               
  if(!antilampool) return;                                                                                                                    
  AliEventPool *xpool =fPoolManagerx->GetEventPool(centralityV0M,vecTarget[2]);                                                           
  if(!xpool) return;                                                                                                                        
  AliEventPool *xppool =fPoolManagerxp->GetEventPool(centralityV0M,vecTarget[2]);                                                          
  if(!xppool) return;

  // lambda + lambda                                                                                                                          
  if(islambdalambda){                                                                                                                           
    for(Int_t iv0=0; iv0 < fV0Arrayn->GetEntriesFast(); iv0++){                                                                              
      TObject *objV0=fV0Arrayn->At(iv0);                                                                                                      
      if(!objV0) continue;                                                                                                                   
      AliAODv0 *lambda1=dynamic_cast<AliAODv0*>(objV0);                                                                                       
      l1px     =lambda1->MomV0X();                                                                                                          
      l1py     =lambda1->MomV0Y();                                                                                                           
      l1pz     =lambda1->MomV0Z();                                                                                                            
      phi1     =atan2(l1py,l1px)*radian_to_degree;                                                                                            
      energy1  =sqrt(pow(1.115,2)+pow(l1px,2)+pow(l1py,2)+pow(l1pz,2));                                                                  
                                                                                                                                     
      //mixed event loop starts                                                                                                            
      for(Int_t imevt=0; imevt < lambdapool->GetCurrentNEvents(); imevt++){                                                           
        TObjArray *poolv0=(TObjArray*)lambdapool->GetEvent(imevt);                                                                      
        for(Int_t imv0=0; imv0<poolv0->GetEntriesFast(); imv0++){                                                                   
          TObject *objpV0=poolv0->At(imv0);                                                                                       
          if(!objpV0) continue;                                                                                                       
          AliAODv0 *lambda2=dynamic_cast<AliAODv0*>(objpV0);                                                                              
	  l2px     =lambda2->MomV0X();                                                                                             
          l2py     =lambda2->MomV0Y();                                                                                              
          l2pz     =lambda2->MomV0Z();                                                                                                
          phi2     =atan2(l2py,l2px)*radian_to_degree;                                                                            
          energy2  =sqrt(pow(1.115,2)+pow(l2px,2)+pow(l2py,2)+pow(l2pz,2));                                                          
          energysum=energy1+energy2;                                                                                                
          deltaphi  =fabs(phi1-phi2);                                                                                                    
	  if(deltaphi > 180){                                                                                                             
	    deltaphi = 360-deltaphi;                                                                                                           
	  }                                                                                                                            
	  pt        =sqrt(pow(l1px+l2px,2)+pow(l1py+l2py,2));                                                                            
	  invMass   =InvariantMass(l1px,l1py,l1pz,l2px,l2py,l2pz,energysum);                                                                
	  relmom    =sqrt(pow((l1px-l2px),2)+pow((l1py-l2py),2)+pow((l1pz-l2pz),2))/2;
	  // cout<<"========== BG lambdalambda =========="<<endl;                                                                  
	  // cout<<"(l1px,l1py,l1pz,l2px,l2py,l2pz) = "<<"("<<l1px<<","<<l1py<<","<<l1pz<<","<<l2px<<","<<l2py<<","<<l2pz<<")"<<endl;        
	  dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassLambdaLambda_evtpool"))->Fill(invMass);     
	  dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaLambdavspT_evtpool"))->Fill(pt,invMass); 
	  //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaLambda_evtpool_deltaphi"))->Fill(deltaphi,invMass);                       
	  //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaLambda_evtpool_centvsmass"))->Fill(centralityV0M,invMass); 
	}                                                                                                                                       
      }                                                                                                                                      
    }                                                                                                                                       
  }
  // antilambda + antilambda                                                                                                                
  if(isantilambdaantilambda){                                                                                                                
    for(Int_t iv0=0; iv0 < fV0Arrayb->GetEntriesFast(); iv0++){                                                                             
      TObject *objV0=fV0Arrayb->At(iv0);                                                                                                      
      if(!objV0) continue;                                                                                                                    
      AliAODv0 *lambda1=dynamic_cast<AliAODv0*>(objV0);                                                                                      
      l1px     =lambda1->MomV0X();                                                                                                             
      l1py     =lambda1->MomV0Y();                                                                                                             
      l1pz     =lambda1->MomV0Z();                                                                                                              
      phi1     =atan2(l1py,l1px)*radian_to_degree;                                                                                            
      energy1  =sqrt(pow(1.115,2)+pow(l1px,2)+pow(l1py,2)+pow(l1pz,2));
      //mixed event loop starts                                                                                                                
      for(Int_t imevt=0; imevt < antilambdapool->GetCurrentNEvents(); imevt++){                                                              
	TObjArray *poolv0=(TObjArray*)antilambdapool->GetEvent(imevt);                                                                        
	for(Int_t imv0=0; imv0<poolv0->GetEntriesFast(); imv0++){                                                                           
	  TObject *objpV0=poolv0->At(imv0);                                                                                                  
	  if(!objpV0) continue;                                                                                                              
	  AliAODv0 *lambda2=dynamic_cast<AliAODv0*>(objpV0);                                                                                   
	  l2px     =lambda2->MomV0X();                                                                                                       
	  l2py     =lambda2->MomV0Y();                                                                                                        
	  l2pz     =lambda2->MomV0Z();                                                                                                       
	  phi2     =atan2(l2py,l2px)*radian_to_degree;                                                                                        
	  energy2  =sqrt(pow(1.115,2)+pow(l2px,2)+pow(l2py,2)+pow(l2pz,2));                                                                   
	  energysum=energy1+energy2;                                                                                                   
	  deltaphi  =fabs(phi1-phi2);                                                                                                         
	  if(deltaphi > 180){                                                                                                                
	    deltaphi = 360-deltaphi;                                                                                             
          }                                                                                                                                  
	  pt        =sqrt(pow(l1px+l2px,2)+pow(l1py+l2py,2));                                                                               
	  invMass   =InvariantMass(l1px,l1py,l1pz,l2px,l2py,l2pz,energysum);                                                                 
	  relmom    =sqrt(pow((l1px-l2px),2)+pow((l1py-l2py),2)+pow((l1pz-l2pz),2))/2;                                                        
	  dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassLambdaLambda_evtpool"))->Fill(invMass);     
	  dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaLambdavspT_evtpool"))->Fill(pt,invMass);     
	  //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaLambda_evtpool_deltaphi"))->Fill(deltaphi,invMass);
	  //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaLambda_evtpool_centvsmass"))->Fill(centralityV0M,invMass); 
        }                                                                                                                                    
      }                                                                                                                                 
    }                                                                                                                                      
  }
  // p + xi-                                                                                                                                          
  // proton                                                                                                                                           
  if(isprotonxi){                                                                                                                                     
    //cout<<"MIX : (im,iz,nevts,mult,zvtx)= ("<<im<<","<<iz<<","<<xipool->GetCurrentNEvents()<<","<<multiplicity<<","<<vecTarget[2]<<")"<<endl;       
    for(Int_t itrk=0; itrk < fProtonArray->GetEntriesFast(); itrk++){                                                                                 
      TObject *objTrack=fProtonArray->At(itrk);                                                                                                       
      if(!objTrack) continue;                                                                                                                         
      AliAODTrack *proton=dynamic_cast<AliAODTrack*>(objTrack);                                                                                       
      prpx        =proton->Px();                                                                                                                      
      prpy        =proton->Py();                                                                                                                      
      prpz        =proton->Pz();                                                                                                                      
      phi1        =atan2(prpy,prpx)*radian_to_degree;                                                                                                 
      energypr    =sqrt(pow(0.938,2)+pow(prpx,2)+pow(prpy,2)+pow(prpz,2));                                                                            
          
      //cout<<"---------- BG pxi proton1 ----------"<<endl;                                                                                          
      //cout<<"(prpx,prpy,prpz) = "<<"("<<prpx<<","<<prpy<<","<<prpz<<")"<<endl;                                                                     
                                                                                                                                                      
      //mixed event loop starts                                                                                                                       
      for(Int_t imevt=0; imevt < xipool->GetCurrentNEvents(); imevt++){                                                                               
        TObjArray *poolcascade=(TObjArray*)xipool->GetEvent(imevt);                                                                                   
        for(Int_t imcascade=0; imcascade < poolcascade->GetEntriesFast(); imcascade++){                                                               
          TObject *objCascade=poolcascade->At(imcascade);                                                                                             
          if(!objCascade) continue;                                                                                                                   
          AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade);                                                                                 
          xipx        =xi->MomXiX();                                                                                                                  
          xipy        =xi->MomXiY();                                                                                                                  
          xipz        =xi->MomXiZ();                                                                                                                  
          phi2        =atan2(xipy,xipx)*radian_to_degree;                                                                                             
          energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));                                                                        
                                                                                                                                                      
	  //cout<<"---------- BG pxi xi1 ----------"<<endl;                                                                                          
	  //cout<<"(xipx,xipy,xipz) = "<<xipx<<","<<xipy<<","<<xipz<<")"<<endl;
	  deltaphi  =fabs(phi1-phi2);                                                                                                                 
          if(deltaphi > 180){                                                                                                                         
            deltaphi = 360-deltaphi;                                                                                                                  
          }                                                                                                                                           
          pt        =sqrt((prpx+xipx)*(prpx+xipx)+(prpy+xipy)*(prpy+xipy));                                                                           
          energysum = energypr+energyxi;                                                                                                              
          invMass   =InvariantMass(xipx,xipy,xipz,prpx,prpy,prpz,energysum);                                                                          
          relmom    =sqrt(pow((xipx-prpx),2)+pow((xipy-prpy),2)+pow((xipz-prpz),2))/2;                                                                
                                                                                                                                                      
          //if(imevt==0) cout<<"========== eventmixing pxi (xi loop) =========="<<endl;                                                               
          //cout<<"(prpx,prpy,prpz,xipx,xipy,xipz) = "<<"("<<prpx<<","<<prpy<<","<<prpz<<","<<xipx<<","<<xipy<<","<<xipz<<")"<<endl;                  
	  //cout<<"Centraity = "<<centralityV0M<<endl;
	  dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassProtonXi_evtpool"))->Fill(invMass);                                                  
	  dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassProtonXivspT_evtpool"))->Fill(pt,invMass);
          //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassProtonXi_evtpool_deltaphi"))->Fill(deltaphi,invMass);
	  //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassProtonXi_evtpool_centvsmass"))->Fill(centralityV0M,invMass); 
                                                                                                                                                      
        }                                                                                                                                             
      }                                                                                                                                               
    }
    // Xi                                                                                                                                             
    for(Int_t icascade=0; icascade < fCascadeArraym->GetEntriesFast(); icascade++){                                                                   
      TObject *objCascade=fCascadeArraym->At(icascade);                                                                                               
      if(!objCascade) continue;                                                                                                                       
      AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade);                                                                                     
      xipx        =xi->MomXiX();                                                                                                                      
      xipy        =xi->MomXiY();                                                                                                                      
      xipz        =xi->MomXiZ();                                                                                                                      
      phi2        =atan2(xipy,xipx)*radian_to_degree;                                                                                                 
      energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));                                                                            
                                                                                                                                                      
      //cout<<"---------- BG pxi xi2 ----------"<<endl;                                                                                              
      //cout<<"(xipx,xipy,xipz) = "<<xipx<<","<<xipy<<","<<xipz<<")"<<endl; 
      //cout<<"Centrality = "<<centralityV0M<<endl;
      //mixed event loop starts                                                                                                                       
      for(Int_t imevt=0; imevt < protonpool->GetCurrentNEvents(); imevt++){                                                                           
        TObjArray *pooltrk=(TObjArray*)protonpool->GetEvent(imevt);                                                                                   
        for(Int_t imtrk=0; imtrk < pooltrk->GetEntriesFast(); imtrk++){                                                                               
          TObject *objProton=pooltrk->At(imtrk);                                                                                                      
          if(!objProton) continue;                                                                                                                    
          AliAODTrack *proton=dynamic_cast<AliAODTrack*>(objProton);                                                                                  
          prpx        =proton->Px();                                                                                                                  
          prpy        =proton->Py();                                                                                                                  
          prpz        =proton->Pz();                                                                                                                  
          phi1        =atan2(prpy,prpx)*radian_to_degree;                                                                                             
          energypr    =sqrt(pow(0.938,2)+pow(prpx,2)+pow(prpy,2)+pow(prpz,2));                                                                        
          energysum   =energyxi+energypr;                                                                                                             
                                                                                                                                                      
          deltaphi  =fabs(phi1-phi2);                                                                                                                 
          if(deltaphi > 180){                                                                                                                         
            deltaphi = 360-deltaphi;                                                                                                                  
          }                                                                                                                                           
          pt        =sqrt((prpx+xipx)*(prpx+xipx)+(prpy+xipy)*(prpy+xipy));                                                                           
          invMass   =InvariantMass(xipx,xipy,xipz,prpx,prpy,prpz,energysum);                                                                          
          relmom    =sqrt(pow((xipx-prpx),2)+pow((xipy-prpy),2)+pow((xipz-prpz),2))/2;                                                                
                                                                                                                                                      
          //cout<<"---------- BG pxi proton2 ----------"<<endl;                                                                                      
          //cout<<"(prpx,prpy,prpz) = "<<"("<<prpx<<","<<prpy<<","<<prpz<<")"<<endl;                                                                 
	  //cout<<"Centrality = "<<centralityV0M<<endl;
                                                                                                                                            
          //if(imevt==0) cout<<"========== eventmixing pxi (proton loop) =========="<<endl;                                                           
          //cout<<"(prpx,prpy,prpz,xipx,xipy,xipz) = "<<"("<<prpx<<","<<prpy<<","<<prpz<<","<<xipx<<","<<xipy<<","<<xipz<<")"<<endl;                  

	  dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassProtonXi_evtpool"))->Fill(invMass);  
	  dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassProtonXivspT_evtpool"))->Fill(pt,invMass);  
          //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassProtonXi_evtpool_deltaphi"))->Fill(deltaphi,invMass);                                  
	  //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassProtonXi_evtpool_centvsmass"))->Fill(centralityV0M,invMass); 
        }                                                                                                                                             
      }                                                                                                                                               
    }                                                                                                                                                 
  }
  // pbar + xi+                                                                                                                                       
  // protonbar                                                                                                                                        
  if(isantiprotonxi){                                                                                                                                 
    for(Int_t itrk=0; itrk < fProtonArrayb->GetEntriesFast(); itrk++){                                                                                
      TObject *objTrack=fProtonArrayb->At(itrk);                                                                                                      
      if(!objTrack) continue;                                                                                                                         
      AliAODTrack *proton=dynamic_cast<AliAODTrack*>(objTrack);                                                                                       
      prpx        =proton->Px();                                                                                                                      
      prpy        =proton->Py();                                                                                                                      
      prpz        =proton->Pz();                                                                                                                      
      phi1        =atan2(prpy,prpx)*radian_to_degree;                                                                                                 
      energypr    =sqrt(pow(0.938,2)+pow(prpx,2)+pow(prpy,2)+pow(prpz,2));                                                                            
                                                                                                                                                      
      //mixed event loop starts                                                                                                                       
      for(Int_t imevt=0; imevt < xippool->GetCurrentNEvents(); imevt++){                                                                              
        TObjArray *poolcascade=(TObjArray*)xippool->GetEvent(imevt);                                                                                  
        for(Int_t imcascade=0; imcascade < poolcascade->GetEntriesFast(); imcascade++){                                                               
          TObject *objCascade=poolcascade->At(imcascade);                                                                                             
          if(!objCascade) continue;                                                                                                                   
          AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade);                                                                                 
          xipx        =xi->MomXiX();                                                                                                                  
          xipy        =xi->MomXiY();                                                                                                                  
          xipz        =xi->MomXiZ();                                                                                                                  
          phi2        =atan2(xipy,xipx)*radian_to_degree;                                                                                             
          energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));                                                                        
                                                                                                                                                      
          deltaphi  =fabs(phi1-phi2);                                                                                                                 
          if(deltaphi > 180){                                                                                                                         
            deltaphi = 360-deltaphi;                                                                                                                  
          }                                                                                                                                           
          pt        =sqrt((prpx+xipx)*(prpx+xipx)+(prpy+xipy)*(prpy+xipy));                                                                           
          energysum = energypr+energyxi;                                                                                                              
          invMass   =InvariantMass(xipx,xipy,xipz,prpx,prpy,prpz,energysum);                                                                          
          relmom    =sqrt(pow((xipx-prpx),2)+pow((xipy-prpy),2)+pow((xipz-prpz),2))/2;                                                                
                                                                                                                                                      
          dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassProtonXi_evtpool"))->Fill(invMass); 
          dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassProtonXivspT_evtpool"))->Fill(pt,invMass);  
          //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassProtonXi_evtpool_deltaphi"))->Fill(deltaphi,invMass);                                  
	  //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassProtonXi_evtpool_centvsmass"))->Fill(centralityV0M,invMass); 
        }                                                                                                                                             
      }                                                                                                                                               
    }
    // Xi+                                                                                                                                            
    for(Int_t icascade=0; icascade < fCascadeArrayp->GetEntriesFast(); icascade++){                                                                   
      TObject *objCascade=fCascadeArrayp->At(icascade);                                                                                               
      if(!objCascade) continue;                                                                                                                       
      AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade);                                                                                     
      xipx        =xi->MomXiX();                                                                                                                      
      xipy        =xi->MomXiY();                                                                                                                      
      xipz        =xi->MomXiZ();                                                                                                                      
      phi2        =atan2(xipy,xipx)*radian_to_degree;                                                                                                 
      energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));
      //mixed event loop starts                                                                                                                       
      for(Int_t imevt=0; imevt < antiprotonpool->GetCurrentNEvents(); imevt++){                                                                       
        TObjArray *pooltrk=(TObjArray*)antiprotonpool->GetEvent(imevt);                                                                               
        for(Int_t imtrk=0; imtrk < pooltrk->GetEntriesFast(); imtrk++){                                                                               
          TObject *objProton=pooltrk->At(imtrk);                                                                                                      
          if(!objProton) continue;                                                                                                                    
          AliAODTrack *proton=dynamic_cast<AliAODTrack*>(objProton);                                                                                  
          prpx        =proton->Px();                                                                                                                  
          prpy        =proton->Py();                                                                                                                  
          prpz        =proton->Pz();                                                                                                                  
          phi1        =atan2(prpy,prpx)*radian_to_degree;                                                                                             
          energypr    =sqrt(pow(0.938,2)+pow(prpx,2)+pow(prpy,2)+pow(prpz,2));                                                                        
          energysum   =energyxi+energypr;                                                                                                             
                                                                                                                                                      
          deltaphi  =fabs(phi1-phi2);                                                                                                                 
          if(deltaphi > 180){                                                                                                                         
            deltaphi = 360-deltaphi;                                                                                                                  
          }                                                                                                                                           
          pt        =sqrt((prpx+xipx)*(prpx+xipx)+(prpy+xipy)*(prpy+xipy));                                                                           
          invMass   =InvariantMass(xipx,xipy,xipz,prpx,prpy,prpz,energysum);                                                                          
          relmom    =sqrt(pow((xipx-prpx),2)+pow((xipy-prpy),2)+pow((xipz-prpz),2))/2;                                                                
                                                                                                                                                      
          dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassProtonXi_evtpool"))->Fill(invMass);
          dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassProtonXivspT_evtpool"))->Fill(pt,invMass);
          //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassProtonXi_evtpool_deltaphi"))->Fill(deltaphi,invMass);                                  
	  //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassProtonXi_evtpool_centvsmass"))->Fill(centralityV0M,invMass); 
        }                                                                                                                                             
      }                                                                                                                                               
    }                                                                                                                                                 
  }
  // lambda + Xi-                                                                                                                                     
  // lambda                                                                                                                                           
  if(islambdaxi){                                                                                                                                     
    for(Int_t iv0=0; iv0 < fV0Arrayn->GetEntriesFast(); iv0++){                                                                                       
      TObject *objV0=fV0Arrayn->At(iv0);                                                                                                              
      if(!objV0) continue;                                                                                                                            
      AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objV0);                                                                                                
      lampx       =lambda->MomV0X();                                                                                                                  
      lampy       =lambda->MomV0Y();                                                                                                                  
      lampz       =lambda->MomV0Z();                                                                                                                  
      phi1        =atan2(lampy,lampx)*radian_to_degree;                                                                                               
      energylam   =sqrt(pow(1.115,2)+pow(lampx,2)+pow(lampy,2)+pow(lampz,2));
      //mixed event loop starts                                                                                                                       
      for(Int_t imevt=0; imevt < xpool->GetCurrentNEvents(); imevt++){                                                                                
        TObjArray *poolcascade=(TObjArray*)xpool->GetEvent(imevt);                                                                                    
        for(Int_t imcascade=0; imcascade < poolcascade->GetEntriesFast(); imcascade++){                                                               
          TObject *objCascade=poolcascade->At(imcascade);                                                                                             
          if(!objCascade) continue;                                                                                                                   
          AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade);                                                                                 
          xipx        =xi->MomXiX();                                                                                                                  
          xipy        =xi->MomXiY();                                                                                                                  
          xipz        =xi->MomXiZ();                                                                                                                  
          phi2        =atan2(xipy,xipx)*radian_to_degree;                                                                                             
          energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));                                                                        
          energysum   =energyxi+energylam;                                                                                                            
                                                                                                                                                      
          deltaphi  =fabs(phi1-phi2);                                                                                                                 
          if(deltaphi > 180){                                                                                                                         
            deltaphi = 360-deltaphi;                                                                                                                  
          }                                                                                                                                           
          pt          =sqrt(pow((lampx+xipx),2)+pow((lampy+xipy),2));                                                                                 
          invMass     =InvariantMass(xipx,xipy,xipz,lampx,lampy,lampz,energysum);                                                                     
          relmom    =sqrt(pow((xipx-lampx),2)+pow((xipy-lampy),2)+pow((xipz-lampz),2))/2;                                                             
                                                                                                                                                      
          dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassLambdaXi_evtpool"))->Fill(invMass);
          dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaXivspT_evtpool"))->Fill(pt,invMass); 
          //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaXi_evtpool_deltaphi"))->Fill(deltaphi,invMass);                                  
	  //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaXi_evtpool_centvsmass"))->Fill(centralityV0M,invMass); 
        }                                                                                                                                             
      }                                                                                                                                               
    } 
    // Xi-                                                                                                                                            
    for(Int_t icascade=0; icascade < fCascadeArraym->GetEntriesFast(); icascade++){                                                                   
      TObject *objCascade=fCascadeArraym->At(icascade);                                                                                               
      if(!objCascade) continue;                                                                                                                       
      AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade);                                                                                     
      xipx        =xi->MomXiX();                                                                                                                      
      xipy        =xi->MomXiY();                                                                                                                      
      xipz        =xi->MomXiZ();                                                                                                                      
      phi2        =atan2(xipy,xipx)*radian_to_degree;                                                                                                 
      energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));
      //mixed event loop starts                                                                                                                       
      for(Int_t imevt=0; imevt < lampool->GetCurrentNEvents(); imevt++){                                                                              
        TObjArray *poolv0=(TObjArray*)lampool->GetEvent(imevt);                                                                                       
        for(Int_t imv0=0; imv0<poolv0->GetEntriesFast(); imv0++){                                                                                     
          TObject *objpV0=poolv0->At(imv0);                                                                                                           
          if(!objpV0) continue;                                                                                                                       
          AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objpV0);                                                                                           
          lampx       =lambda->MomV0X();                                                                                                              
          lampy       =lambda->MomV0Y();                                                                                                              
          lampz       =lambda->MomV0Z();                                                                                                              
          phi1        =atan2(lampy,lampx)*radian_to_degree;                                                                                           
          energylam   =sqrt(pow(1.115,2)+pow(lampx,2)+pow(lampy,2)+pow(lampz,2));                                                                     
          energysum   =energyxi+energylam;                                                                                                            
                                                                                                                                                      
          deltaphi  =fabs(phi1-phi2);                                                                                                                 
          if(deltaphi > 180){                                                                                                                         
            deltaphi = 360-deltaphi;                                                                                                                  
          }                                                                                                                                           
          pt          =sqrt(pow((lampx+xipx),2)+pow((lampy+xipy),2));                                                                                 
          invMass     =InvariantMass(xipx,xipy,xipz,lampx,lampy,lampz,energysum);                                                                     
          relmom    =sqrt(pow((xipx-lampx),2)+pow((xipy-lampy),2)+pow((xipz-lampz),2))/2;                                                             
                                                                                                                                                      
          dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassLambdaXi_evtpool"))->Fill(invMass);    
          dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaXivspT_evtpool"))->Fill(pt,invMass); 
          //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaXi_evtpool_deltaphi"))->Fill(deltaphi,invMass);                                  
	  //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaXi_evtpool_centvsmass"))->Fill(centralityV0M,invMass); 
        }                                                                                                                                             
      }                                                                                                                                               
    }                                                                                                                                                 
  }
  // antilambda + Xi+                                                                                                                                 
  // antilambda                                                                                                                                       
  if(isantilambdaxi){                                                                                                                                 
    for(Int_t iv0=0; iv0 < fV0Arrayb->GetEntriesFast(); iv0++){                                                                                       
      TObject *objV0=fV0Arrayb->At(iv0);                                                                                                              
      if(!objV0) continue;                                                                                                                            
      AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objV0);                                                                                                
      lampx       =lambda->MomV0X();                                                                                                                  
      lampy       =lambda->MomV0Y();                                                                                                                  
      lampz       =lambda->MomV0Z();                                                                                                                  
      phi1        =atan2(lampy,lampx)*radian_to_degree;                                                                                               
      energylam   =sqrt(pow(1.115,2)+pow(lampx,2)+pow(lampy,2)+pow(lampz,2));                                                                         
                                                                                                                                                      
      //mixed event loop starts                                                                                                                       
      for(Int_t imevt=0; imevt < xppool->GetCurrentNEvents(); imevt++){                                                                               
        TObjArray *poolcascade=(TObjArray*)xppool->GetEvent(imevt);                                                                                   
        for(Int_t imcascade=0; imcascade < poolcascade->GetEntriesFast(); imcascade++){                                                               
          TObject *objCascade=poolcascade->At(imcascade);                                                                                             
          if(!objCascade) continue;                                                                                                                   
          AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade);                                                                                 
          xipx        =xi->MomXiX();                                                                                                                  
          xipy        =xi->MomXiY();                                                                                                                  
          xipz        =xi->MomXiZ();                                                                                                                  
          phi2        =atan2(xipy,xipx)*radian_to_degree;                                                                                             
          energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));                                                                        
          energysum   =energyxi+energylam;
	  deltaphi  =fabs(phi1-phi2);                                                                                                                 
          if(deltaphi > 180){                                                                                                                         
            deltaphi = 360-deltaphi;                                                                                                                  
          }                                                                                                                                           
          pt          =sqrt(pow((lampx+xipx),2)+pow((lampy+xipy),2));                                                                                 
          invMass     =InvariantMass(xipx,xipy,xipz,lampx,lampy,lampz,energysum);                                                                     
          relmom    =sqrt(pow((xipx-lampx),2)+pow((xipy-lampy),2)+pow((xipz-lampz),2))/2;                                                             
                                                                                                                                                      
          dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassLambdaXi_evtpool"))->Fill(invMass);                                   
          dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaXivspT_evtpool"))->Fill(pt,invMass); 
          //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaXi_evtpool_deltaphi"))->Fill(deltaphi,invMass);                                  
	  //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaXi_evtpool_centvsmass"))->Fill(centralityV0M,invMass); 
        }                                                                                                                                             
      }                                                                                                                                               
    }
    // Xi+                                                                                                                                            
    for(Int_t icascade=0; icascade < fCascadeArrayp->GetEntriesFast(); icascade++){                                                                   
      TObject *objCascade=fCascadeArrayp->At(icascade);                                                                                               
      if(!objCascade) continue;                                                                                                                       
      AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade);                                                                                     
      xipx        =xi->MomXiX();                                                                                                                      
      xipy        =xi->MomXiY();                                                                                                                      
      xipz        =xi->MomXiZ();                                                                                                                      
      phi2        =atan2(xipy,xipx)*radian_to_degree;                                                                                                 
      energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));                                                                            
                                                                                                                                                      
      //mixed event loop starts                                                                                                                       
      for(Int_t imevt=0; imevt < antilampool->GetCurrentNEvents(); imevt++){                                                                          
        TObjArray *poolv0=(TObjArray*)antilampool->GetEvent(imevt);                                                                                   
        for(Int_t imv0=0; imv0<poolv0->GetEntriesFast(); imv0++){                                                                                     
          TObject *objpV0=poolv0->At(imv0);                                                                                                           
          if(!objpV0) continue;                                                                                                                       
          AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objpV0);                                                                                           
          lampx       =lambda->MomV0X();                                                                                                              
          lampy       =lambda->MomV0Y();                                                                                                              
          lampz       =lambda->MomV0Z();                                                                                                              
          phi1        =atan2(lampy,lampx)*radian_to_degree;                                                                                           
          energylam   =sqrt(pow(1.115,2)+pow(lampx,2)+pow(lampy,2)+pow(lampz,2));                                                                     
          energysum   =energyxi+energylam; 
	  deltaphi  =fabs(phi1-phi2);                                                                                                                 
          if(deltaphi > 180){                                                                                                                         
            deltaphi = 360-deltaphi;                                                                                                                  
          }                                                                                                                                           
          pt          =sqrt(pow((lampx+xipx),2)+pow((lampy+xipy),2));                                                                                 
          invMass     =InvariantMass(xipx,xipy,xipz,lampx,lampy,lampz,energysum);                                                                     
          relmom    =sqrt(pow((xipx-lampx),2)+pow((xipy-lampy),2)+pow((xipz-lampz),2))/2;                                                             
                                                                                                                                                      
          dynamic_cast<TH1F*>(fOutputList->FindObject("hInvMassLambdaXi_evtpool"))->Fill(invMass);
          dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaXivspT_evtpool"))->Fill(pt,invMass);
          //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaXi_evtpool_deltaphi"))->Fill(deltaphi,invMass);                                  
	  //dynamic_cast<TH2F*>(fOutputList->FindObject("hInvMassLambdaXi_evtpool_centvsmass"))->Fill(centralityV0M,invMass); 
        }                                                                                                                                             
      }                                                                                                                                               
    }                                                                                                                                                 
  }

  //lambda                                                                                                                                
  if(islambdalambda){                                                                                                                   
    TObjArray* tracksClonelambda=(TObjArray*)fV0Arrayn->Clone();                                                                         
    tracksClonelambda->SetOwner();                                                                                                      
    lambdapool->UpdatePool(tracksClonelambda);                                                                                              
  }                                                                                                                                     
  //lambda bar                                                                                                                             
  if(isantilambdaantilambda){                                                                                                              
    TObjArray* tracksClonelambdab=(TObjArray*)fV0Arrayb->Clone();                                                                        
    tracksClonelambdab->SetOwner();                                                                                                        
    antilambdapool->UpdatePool(tracksClonelambdab);                                                                                         
  }
  //pXi-                                                                                                                                   
  if(isprotonxi){                                                                                                                           
    //xi                                                                                                                                 
    TObjArray* tracksClonexim=(TObjArray*)fCascadeArraym->Clone();                                                                          
    tracksClonexim->SetOwner();                                                                                                              
    xipool->UpdatePool(tracksClonexim);                                                                                                      
    //proton                                                                                                                              
    TObjArray* tracksClonep=(TObjArray*)fProtonArray->Clone();                                                                            
    tracksClonep->SetOwner();                                                                                                           
    protonpool->UpdatePool(tracksClonep);                                                                                              
    //cout<<"POOL: (im,iz,mult,zvtx,nevts)= ("<<im<<","<<iz<<","<<xipool->GetCurrentNEvents()<<","<<multiplicity<<","<<vecTarget[2]<<")"<<endl;  
  }                                                                                                                                        
  //pbarXi+                                                                                                                                  
  if(isantiprotonxi){                                                                                                                         
    //xip                                                                                                                                      
    TObjArray* tracksClonexip=(TObjArray*)fCascadeArrayp->Clone();                                                                            
    tracksClonexip->SetOwner();                                                                                                               
    xippool->UpdatePool(tracksClonexip);                                                                                                     
    //proton bar                                                                                                                              
    TObjArray* tracksClonepb=(TObjArray*)fProtonArrayb->Clone();                                                                               
    tracksClonepb->SetOwner();                                                                                                                 
    antiprotonpool->UpdatePool(tracksClonepb);                                                                                                  
  }
  //lambdaXi                                                                                                                                  
  if(islambdaxi){                                                                                                                             
  //lambda                                                                                                                                   
    TObjArray* tracksClonelam=(TObjArray*)fV0Arrayn->Clone();                                                                                 
    tracksClonelam->SetOwner();                                                                                                               
    lampool->UpdatePool(tracksClonelam);                                                                                                       
    //Xi-                                                                                                                                   
    TObjArray* tracksClonex=(TObjArray*)fCascadeArraym->Clone();                                                                             
    tracksClonex->SetOwner();                                                                                                                   
    xpool->UpdatePool(tracksClonex);                                                                                                          
  }                                                                                                                                            

  //lambdabarXip                                                                                                                             
  if(isantilambdaxi){                                                                                                                         
    //lambdabar                                                                                                                               
    TObjArray* tracksClonelamb=(TObjArray*)fV0Arrayb->Clone();                                                                                
    tracksClonelamb->SetOwner();                                                                                                               
    antilampool->UpdatePool(tracksClonelamb);                                                                                                 
    TObjArray* tracksClonexp=(TObjArray*)fCascadeArrayp->Clone();                                                                           
    tracksClonexp->SetOwner();                                                                                                         
    xppool->UpdatePool(tracksClonexp);                                                                                                        
  }

  // Post output data
  PostData(1,fOutputList);
  cout<<"+++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"+++++++++++++++ Analysis finish +++++++++++++"<<endl;
  cout<<"+++++++++++++++++++++++++++++++++++++++++++++"<<endl;

}


//________________________________________________________________________
Bool_t AliAnalysisTaskMSDibaryons::EventSelection(AliAODEvent* data)
{
  const AliVVertex *vtx   =data->GetPrimaryVertex();
  const AliVVertex *vtxSPD=data->GetPrimaryVertexSPD();  
  
  Double_t xvtx=0.,yvtx=0.,zvtx=0.;
  Double_t xvtxSPD=0.,yvtxSPD=0.,zvtxSPD=0.;
  Int_t    ncont=0,ncontSPD=0;
  Double_t vdisp=0.;
  // Event vertex
  xvtx =vtx->GetX();
  yvtx =vtx->GetY();
  zvtx =vtx->GetZ();
  ncont=vtx->GetNContributors();

  dynamic_cast<TH1F*>(fOutputList->FindObject("fEvtVtxX"))->Fill(xvtx);
  dynamic_cast<TH1F*>(fOutputList->FindObject("fEvtVtxY"))->Fill(yvtx);
  dynamic_cast<TH1F*>(fOutputList->FindObject("fEvtVtxZ"))->Fill(zvtx);
  dynamic_cast<TH1F*>(fOutputList->FindObject("fEvtVtxTrk"))->Fill(ncont);
  // SPD vertex
  xvtxSPD =vtxSPD->GetX();
  yvtxSPD =vtxSPD->GetY();
  zvtxSPD =vtxSPD->GetZ();
  ncontSPD=vtxSPD->GetNContributors();
  vdisp=fabs(zvtx-zvtxSPD);

  dynamic_cast<TH1F*>(fOutputList->FindObject("fSPDVtxZ"))->Fill(zvtxSPD);
  dynamic_cast<TH1F*>(fOutputList->FindObject("fSPDVtxTrk"))->Fill(ncontSPD);
  dynamic_cast<TH2F*>(fOutputList->FindObject("fSPDVtxCor"))->Fill(zvtx,zvtxSPD);
  dynamic_cast<TH1F*>(fOutputList->FindObject("fSPDVtxDisp"))->Fill(vdisp);

  Bool_t zvtx_cut=kFALSE;
  //zvtx_cut=(fabs(zvtx)<10. && ncont>0 && ncontSPD>0 && vdisp<0.5);
  zvtx_cut=(abs(zvtx) < 10. && ncont > 0);
  return zvtx_cut;
}
//________________________________________________________________________

//invariant mass Lambda cascade class
Double_t AliAnalysisTaskMSDibaryons::InvMassLambda(AliAODcascade *casc)
{  
  Double_t pPx=0.,pPy=0.,pPz=0.,pE=0.; //proton
  Double_t nPx=0.,nPy=0.,nPz=0.,nE=0.; //pion-
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;  
  pPx=casc->MomPosX();
  pPy=casc->MomPosY();
  pPz=casc->MomPosZ();
  nPx=casc->MomNegX();
  nPy=casc->MomNegY();
  nPz=casc->MomNegZ();
  pE =sqrt(0.938*0.938+pPx*pPx+pPy*pPy+pPz*pPz); //proton
  nE =sqrt(0.14*0.14+nPx*nPx+nPy*nPy+nPz*nPz);   //pion-
  energysum =pE+nE;
  psum2     =(pPx+nPx)*(pPx+nPx)+(pPy+nPy)*(pPy+nPy)+(pPz+nPz)*(pPz+nPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0; 
  return invMass;
}
//invariant mass anti lambda
Double_t AliAnalysisTaskMSDibaryons::InvMassAntiLambda(AliAODcascade *casc)
{  
  Double_t nPx=0.,nPy=0.,nPz=0.,nE=0.; //anti proton
  Double_t pPx=0.,pPy=0.,pPz=0.,pE=0.; //pion+
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;  
  nPx=casc->MomNegX();
  nPy=casc->MomNegY();
  nPz=casc->MomNegZ();
  pPx=casc->MomPosX();
  pPy=casc->MomPosY();
  pPz=casc->MomPosZ();
  nE =sqrt(0.938*0.938+nPx*nPx+nPy*nPy+nPz*nPz); //anti proton
  pE =sqrt(0.14*0.14+pPx*pPx+pPy*pPy+pPz*pPz);   //pion+
  energysum =nE+pE;
  psum2     =(nPx+pPx)*(nPx+pPx)+(nPy+pPy)*(nPy+pPy)+(nPz+pPz)*(nPz+pPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0; 
  return invMass;
}
//invariant mass Xi-
Double_t AliAnalysisTaskMSDibaryons::InvMassXi(AliAODcascade *casc)
{  
  Double_t lPx=0.,lPy=0.,lPz=0.,lE=0.; //lambda
  Double_t bPx=0.,bPy=0.,bPz=0.,bE=0.; //bachelor pion-
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;  
  lPx=casc->MomV0X();
  lPy=casc->MomV0Y();
  lPz=casc->MomV0Z();
  bPx=casc->MomBachX();
  bPy=casc->MomBachY();
  bPz=casc->MomBachZ();
  lE =sqrt(1.1157*1.1157+lPx*lPx+lPy*lPy+lPz*lPz); //lambda 1.115
  bE =sqrt(0.14*0.14+bPx*bPx+bPy*bPy+bPz*bPz);   //pion-
  energysum =lE+bE;
  psum2     =(lPx+bPx)*(lPx+bPx)+(lPy+bPy)*(lPy+bPy)+(lPz+bPz)*(lPz+bPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0; 
  return invMass;
}
//omega rejection
Double_t AliAnalysisTaskMSDibaryons::InvMassOmega(AliAODcascade *casc)
{  
  Double_t lPx=0.,lPy=0.,lPz=0.,lE=0.; //lambda
  Double_t bPx=0.,bPy=0.,bPz=0.,bE=0.; //bachelor kaon
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;  
  lPx=casc->MomV0X();
  lPy=casc->MomV0Y();
  lPz=casc->MomV0Z();
  bPx=casc->MomBachX();
  bPy=casc->MomBachY();
  bPz=casc->MomBachZ();
  lE =sqrt(1.1157*1.1157+lPx*lPx+lPy*lPy+lPz*lPz); //lambda 1.116
  bE =sqrt(0.494*0.494+bPx*bPx+bPy*bPy+bPz*bPz); //kaon-
  energysum =lE+bE;
  psum2     =(lPx+bPx)*(lPx+bPx)+(lPy+bPy)*(lPy+bPy)+(lPz+bPz)*(lPz+bPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0; 
  return invMass;
}
//lambda cosine pointing angle
Double_t AliAnalysisTaskMSDibaryons::LambdaCosPointingAngle(AliAODcascade *casc,const Double_t *DecayVtx,
							    const Float_t *point) const 
{  
  TVector3 v0Mom(casc->MomV0X(),casc->MomV0Y(),casc->MomV0Z());
  TVector3 fline(DecayVtx[0] - point[0], DecayVtx[1] - point[1],
                 DecayVtx[2] - point[2]);
  Double_t ptot2 = v0Mom.Mag2() * fline.Mag2();
  if (ptot2 <= 0) {
    return 0.0;
  }
  else {
    Double_t cos = v0Mom.Dot(fline) / TMath::Sqrt(ptot2);
    if (cos > 1.0)
      cos = 1.0;
    if (cos < -1.0)
      cos = -1.0;
    return cos;
  }
}
//transverse radius of the lambda decay vertex
Double_t AliAnalysisTaskMSDibaryons::DecayLengthXY(const Double_t *DecayVtx,const Float_t *point) const 
{ return TMath::Sqrt( (point[0] - DecayVtx[0]) * (point[0] - DecayVtx[0])
		      + (point[1] - DecayVtx[1]) * (point[1] - DecayVtx[1]));
}
//transverse radius of the xi decay vertex
Double_t AliAnalysisTaskMSDibaryons::xiDecayLengthXY(const Double_t *xiDecayVtx,const Float_t *point) const 
{ return TMath::Sqrt( (point[0] - xiDecayVtx[0]) * (point[0] - xiDecayVtx[0])
		      + (point[1] - xiDecayVtx[1]) * (point[1] - xiDecayVtx[1]));
}
//invariant mass Lambda v0 class
Double_t AliAnalysisTaskMSDibaryons::InvMasslambda(AliAODv0 *v0)
{  
  Double_t pPx=0.,pPy=0.,pPz=0.,pE=0.; //proton
  Double_t nPx=0.,nPy=0.,nPz=0.,nE=0.; //pion
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;
  pPx=v0->MomPosX();
  pPy=v0->MomPosY();
  pPz=v0->MomPosZ();
  nPx=v0->MomNegX();
  nPy=v0->MomNegY();
  nPz=v0->MomNegZ();
  pE =sqrt(0.938*0.938+pPx*pPx+pPy*pPy+pPz*pPz); //proton
  nE =sqrt(0.14*0.14+nPx*nPx+nPy*nPy+nPz*nPz);   //pion-
  energysum =pE+nE;
  psum2     =(pPx+nPx)*(pPx+nPx)+(pPy+nPy)*(pPy+nPy)+(pPz+nPz)*(pPz+nPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0; 
  return invMass;
}
//K0 rejection                                                                                         
Double_t AliAnalysisTaskMSDibaryons::InvMassK0(AliAODv0 *v0)
{
  Double_t pPx=0.,pPy=0.,pPz=0.,pE=0.; //pion+                                                         
  Double_t nPx=0.,nPy=0.,nPz=0.,nE=0.; //pion                                                          
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;
  pPx=v0->MomPosX();
  pPy=v0->MomPosY();
  pPz=v0->MomPosZ();
  nPx=v0->MomNegX();
  nPy=v0->MomNegY();
  nPz=v0->MomNegZ();
  pE =sqrt(0.14*0.14+pPx*pPx+pPy*pPy+pPz*pPz); //pion+                                                 
  nE =sqrt(0.14*0.14+nPx*nPx+nPy*nPy+nPz*nPz); //pion-                                                 
  energysum =pE+nE;
  psum2     =(pPx+nPx)*(pPx+nPx)+(pPy+nPy)*(pPy+nPy)+(pPz+nPz)*(pPz+nPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0;
  return invMass;
}
//invariant mass antiLambda v0 class
Double_t AliAnalysisTaskMSDibaryons::InvMassAntilambda(AliAODv0 *v0)
{
  Double_t pPx=0.,pPy=0.,pPz=0.,pE=0.; //pion+                                                                                               
  Double_t nPx=0.,nPy=0.,nPz=0.,nE=0.; //antiproton                                                                                            
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;
  pPx=v0->MomPosX();
  pPy=v0->MomPosY();
  pPz=v0->MomPosZ();
  nPx=v0->MomNegX();
  nPy=v0->MomNegY();
  nPz=v0->MomNegZ();
  pE =sqrt(0.14*0.14+pPx*pPx+pPy*pPy+pPz*pPz);     //pion+                                                                                    
  nE =sqrt(0.938*0.938+nPx*nPx+nPy*nPy+nPz*nPz);   //antiproton                                                                              
  energysum =pE+nE;
  psum2     =(pPx+nPx)*(pPx+nPx)+(pPy+nPy)*(pPy+nPy)+(pPz+nPz)*(pPz+nPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0;
  return invMass;
}
Double_t AliAnalysisTaskMSDibaryons::CosPointingAngle(AliAODv0 *v0,const Double_t *DecayVtx,
						      const Float_t *point) const {  
  /// Cosine of pointing angle in space assuming it is produced at "point"
  TVector3 v0Mom(v0->MomV0X(),v0->MomV0Y(),v0->MomV0Z());
  TVector3 fline(DecayVtx[0] - point[0], DecayVtx[1] - point[1],
                 DecayVtx[2] - point[2]);
  Double_t ptot2 = v0Mom.Mag2() * fline.Mag2();
  if (ptot2 <= 0) {
    return 0.0;
  }
  else {
    Double_t cos = v0Mom.Dot(fline) / TMath::Sqrt(ptot2);
    if (cos > 1.0)
      cos = 1.0;
    if (cos < -1.0)
      cos = -1.0;
    return cos;
  }
}

Double_t AliAnalysisTaskMSDibaryons::OpenAngle(Double_t px1,Double_t py1,Double_t pz1,
					       Double_t px2,Double_t py2,Double_t pz2){
  Double_t lScalPtot1Ptot2=0.,lPtot1xPtot2=0.;

  lScalPtot1Ptot2 = (px1*px2)+(py1*py2)+(pz1*pz2);
  lPtot1xPtot2 = (sqrt(pow(px1,2)+pow(py1,2)+pow(pz1,2)))*(sqrt(pow(px2,2)+pow(py2,2)+pow(pz2,2)));
  
  return acos(lScalPtot1Ptot2/lPtot1xPtot2);
}

Double_t AliAnalysisTaskMSDibaryons::InvariantMass(Double_t px1,Double_t py1,Double_t pz1,
						   Double_t px2,Double_t py2,Double_t pz2,Double_t energysum){

  Double_t psum2=0.,pt=0.,invMass=0.;

  psum2     =pow((px1+px2),2)+pow((py1+py2),2)+pow((pz1+pz2),2);
  pt        =sqrt(pow(px1+px2,2)+pow(py1+py2,2));
  invMass   =sqrt((energysum*energysum)-psum2);

  return invMass;

}

//Bool_t AliAnalysisTaskMSDibaryons::ExtractQnVector()
Double_t AliAnalysisTaskMSDibaryons::ExtractQnVector()
{
  Int_t fHarmonics =2;
  Int_t fQnDetectorMain =0;

  //TString fTPCEPName[0] = "TPC";
  //TString fTPCEPName[1] = "TPCNegEta";
  //TString fTPCEPName[2] = "TPCPosEta";
  TString fTPCEPName[3] = {"TPC","TPCNegEta","TPCPosEta"};

  //fV0EPName[0] = "VZERO";
  //fV0EPName[1] = "VZEROA";
  //fV0EPName[2] = "VZEROC";

  TString fV0EPName[3] = {"VZERO","VZEROA","VZEROC"};
  TString fQNormalization = "QoverM";

  if(fHarmonics < 0){
    AliError(Form("Qn Flow vector correction flag is ON, but fHarmonics is not set. (it is %d now).",fHarmonics));
    return kFALSE;
  }

  TList* qnlist = fFlowQnVectorMgr->GetQnVectorList();

  const AliQnCorrectionsQnVector *QnVectorTPCDet[3];
  Double_t TPCEP[3] = {};
  for(Int_t i=0;i<3;i++){
    QnVectorTPCDet[i] = GetQnVectorFromList(qnlist,Form("%s%s",fTPCEPName[i].Data(),fQNormalization.Data()),"latest","plain");
    if(!QnVectorTPCDet[i]){
      AliInfo("Event is rejected because event plane is not found or bad event plane quality in TPC.");
      Printf("Event is rejected because event plane is not found or bad event plane quality in TPC.");
      return kFALSE;//Qn vector correction does not exist or bad quality.
    }
    TPCEP[i] = QnVectorTPCDet[i]->EventPlane(fHarmonics);
    if(TPCEP[i] < 0) TPCEP[i] += 2./(Double_t) fHarmonics * TMath::Pi();
    AliInfo(Form("harmonics %d | TPC sub detector name %s%s : event plane = %f (rad).",fHarmonics,fTPCEPName[i].Data(),fQNormalization.Data(),TPCEP[i]));
  }

  const AliQnCorrectionsQnVector *QnVectorV0Det[3];
  Double_t V0EP[3]  = {};
  for(Int_t i=0;i<3;i++){
    QnVectorV0Det[i]  = GetQnVectorFromList(qnlist,Form("%s%s",fV0EPName[i].Data(),fQNormalization.Data()),"latest","raw");
    if(!QnVectorV0Det[i]){
      AliInfo("Event is rejected because event plane is not found or bad event plane quality in VZERO.");
      Printf("Event is rejected because event plane is not found or bad event plane quality in VZERO.");
      return kFALSE;//Qn vector correction does not exist or bad quality.
    }
    V0EP[i] = QnVectorV0Det[i]->EventPlane(fHarmonics);
    if(V0EP[i] < 0)  V0EP[i]  += 2./(Double_t) fHarmonics * TMath::Pi();
    AliInfo(Form("harmonics %d | V0  sub detector name %s%s : event plane = %f (rad).",fHarmonics,fV0EPName[i].Data(),fQNormalization.Data(),V0EP[i]));
    //Printf("harmonics %d | V0  sub detector name %s%s : event plane = %f (rad).",fHarmonics,fV0EPName[i].Data(),fQNormalization.Data(),V0EP[i]);
  }

  //0 < event plane < 2*pi/fHarmonics.
  Double_t EP2 = -999; 
  Double_t EP3 = -999; 

  Double_t Q1[2] = {};//for Main
  Double_t Q2[2] = {};//for Sub1
  Double_t Q3[2] = {};//for Sub2

  Float_t fEventPlane = 999;

  if(fQnDetectorMain == AliAnalysisTaskMSDibaryons::kFullTPC){
    Q1[0] = QnVectorTPCDet[0]->Qx(fHarmonics);//FullTPC
    Q1[1] = QnVectorTPCDet[0]->Qy(fHarmonics);//FullTPC
    Q2[0] = QnVectorV0Det[1]->Qx(fHarmonics);//V0A
    Q2[1] = QnVectorV0Det[1]->Qy(fHarmonics);//V0A
    Q3[0] = QnVectorV0Det[2]->Qx(fHarmonics);//V0C
    Q3[1] = QnVectorV0Det[2]->Qy(fHarmonics);//V0C
    fEventPlane = TPCEP[0];
    EP2 = V0EP[1];
    EP3 = V0EP[2];
  }
  else if(fQnDetectorMain == AliAnalysisTaskMSDibaryons::kTPCNegEta){
    Q1[0] = QnVectorTPCDet[1]->Qx(fHarmonics);//TPCNegEta
    Q1[1] = QnVectorTPCDet[1]->Qy(fHarmonics);//TPCNegEta
    Q2[0] = QnVectorV0Det[1]->Qx(fHarmonics);//V0A
    Q2[1] = QnVectorV0Det[1]->Qy(fHarmonics);//V0A
    Q3[0] = QnVectorV0Det[2]->Qx(fHarmonics);//V0C
    Q3[1] = QnVectorV0Det[2]->Qy(fHarmonics);//V0C
    fEventPlane = TPCEP[2];
    EP2 = V0EP[1];
    EP3 = V0EP[2];
  }
  else if(fQnDetectorMain == AliAnalysisTaskMSDibaryons::kTPCPosEta){
    Q1[0] = QnVectorTPCDet[2]->Qx(fHarmonics);//TPCPosEta
    Q1[1] = QnVectorTPCDet[2]->Qy(fHarmonics);//TPCPosEta
    Q2[0] = QnVectorV0Det[1]->Qx(fHarmonics);//V0A
    Q2[1] = QnVectorV0Det[1]->Qy(fHarmonics);//V0A
    Q3[0] = QnVectorV0Det[2]->Qx(fHarmonics);//V0C
    Q3[1] = QnVectorV0Det[2]->Qy(fHarmonics);//V0C
    fEventPlane = TPCEP[1];
    EP2 = V0EP[1];
    EP3 = V0EP[2];
  }
  else if(fQnDetectorMain == AliAnalysisTaskMSDibaryons::kFullV0){
    Q1[0] = QnVectorV0Det[0]->Qx(fHarmonics);//FullV0
    Q1[1] = QnVectorV0Det[0]->Qy(fHarmonics);//FullV0
    Q2[0] = QnVectorTPCDet[1]->Qx(fHarmonics);//TPCNegEta
    Q2[1] = QnVectorTPCDet[1]->Qy(fHarmonics);//TPCNegEta
    Q3[0] = QnVectorTPCDet[2]->Qx(fHarmonics);//TPCPosEta
    Q3[1] = QnVectorTPCDet[2]->Qy(fHarmonics);//TPCPosEta
    fEventPlane = V0EP[0];
    EP2 = TPCEP[1];
    EP3 = TPCEP[2];
  }
  else if(fQnDetectorMain == AliAnalysisTaskMSDibaryons::kV0A){
    Q1[0] = QnVectorV0Det[1]->Qx(fHarmonics);//V0A
    Q1[1] = QnVectorV0Det[1]->Qy(fHarmonics);//V0A
    Q2[0] = QnVectorV0Det[2]->Qx(fHarmonics);//V0C
    Q2[1] = QnVectorV0Det[2]->Qy(fHarmonics);//V0C
    Q3[0] = QnVectorTPCDet[0]->Qx(fHarmonics);//full acceptance of TPC
    Q3[1] = QnVectorTPCDet[0]->Qy(fHarmonics);//full acceptance of TPC
    fEventPlane = V0EP[1];
    EP2 = V0EP[2];
    EP3 = TPCEP[0];
  }
  else if(fQnDetectorMain == AliAnalysisTaskMSDibaryons::kV0C){
    Q1[0] = QnVectorV0Det[2]->Qx(fHarmonics);//V0C
    Q1[1] = QnVectorV0Det[2]->Qy(fHarmonics);//V0C
    Q2[0] = QnVectorV0Det[1]->Qx(fHarmonics);//V0A
    Q2[1] = QnVectorV0Det[1]->Qy(fHarmonics);//V0A
    Q3[0] = QnVectorTPCDet[0]->Qx(fHarmonics);//full acceptance of TPC
    Q3[1] = QnVectorTPCDet[0]->Qy(fHarmonics);//full acceptance of TPC
    fEventPlane = V0EP[2];
    EP2 = V0EP[1];
    EP3 = TPCEP[0];
  }

  //fQVector1.Set(Q1[0],Q1[1]);
  TVector2 fQVector1(Q1[0],Q1[1]);
  TVector2 QVector2(Q2[0],Q2[1]);
  TVector2 QVector3(Q3[0],Q3[1]);
  Double_t sp12 = fQVector1 *  QVector2;//scalar product between Q1 vector and Q2 vector
  Double_t sp23 =  QVector2 *  QVector3;//scalar product between Q2 vector and Q3 vector
  Double_t sp31 =  QVector3 * fQVector1;//scalar product between Q3 vector and Q1 vector

  //AliInfo(Form("Q1x = %e , Q1y = %e , Q2x = %e , Q2y = %e , Q3x = %e , Q3y = %e ,  SP12 = %e ,  SP23 = %e ,  SP31 = %e",Q1[0],Q1[1],Q2[0],Q2[1],Q3[0],Q3[1],sp12,sp23,sp31));
  //Printf("Q1x = %e , Q1y = %e , Q2x = %e , Q2y = %e , Q3x = %e , Q3y = %e ,  SP12 = %e ,  SP23 = %e ,  SP31 = %e",Q1[0],Q1[1],Q2[0],Q2[1],Q3[0],Q3[1],sp12,sp23,sp31);
  
  //dynamic_cast<TH1F*>(fOutputList->FindObject("hEventPlane"))->Fill(EP2);
  dynamic_cast<TH1F*>(fOutputList->FindObject("hEventPlane"))->Fill(fEventPlane);
  //cout<<"========== EventPlane!! =============="<<endl;
  //cout<<fEventPlane<<endl;

  //  const Double_t delta = 2. * TMath::Pi() / Double_t(fHarmonics) / 12.;
  //  fEPBin = (Int_t)((fEventPlane) / delta);//it should be 0-11.
  //  if(fEPBin < 0)  fEPBin =  0;//protection to avoid fEPBin = -1.
  //  if(fEPBin > 11) fEPBin = 11;//protection to avoid fEPBin = 12.

  //return kTRUE;
  return fEventPlane;
}
const AliQnCorrectionsQnVector *AliAnalysisTaskMSDibaryons::GetQnVectorFromList(const TList *list, const char* subdetector, const char *expcorr, const char *altcorr)
{
  AliQnCorrectionsQnVector *theQnVector = NULL;

  TList *pQvecList = dynamic_cast<TList*> (list->FindObject(subdetector));
  if(pQvecList != NULL){//sub detector is found
    if(TString(expcorr).EqualTo("latest")) theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
    else                                   theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(expcorr);
    //theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(expcorr);

    if(theQnVector == NULL || !(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0)){ //the Qn vector for the expected correction is not found
      AliInfo(Form("expected correction (%s) is not found. use %s as an alternative step in %s.",expcorr,altcorr,subdetector));
      if(TString(altcorr).EqualTo("latest")) theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
      else                                   theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(altcorr);
      //theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(altcorr);
    }

    //check the Qn vector quality
    if(!(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0)) theQnVector = NULL; //bad quality, discarded

  }
  return theQnVector;

}

