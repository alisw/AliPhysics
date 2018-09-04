#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TObjString.h"
#include "TRandom.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskRhopPb.h"
#include "AliPHOSGeometry.h"
#include "AliESDEvent.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloCluster.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliPHOSGeoUtils.h"
#include "AliAnalysisManager.h"
#include "AliPID.h"
#include <AliPID.h>
#include <AliPIDCombined.h>
#include <AliPIDResponse.h>
#include "AliCentrality.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliTender.h"
#include "AliTenderSupply.h"
#include "AliT0TenderSupply.h"   
#include "AliTOFTenderSupply.h"  
#include "AliPIDTenderSupply.h"    
#include "AliTPCTenderSupply.h"  

#include "AliEventCuts.h"


// *****************************************************
//        REAL DATA
// Analysis task to create trees for rho in pPb analysis


ClassImp(AliAnalysisTaskRhopPb)

//________________________________________________________________________
AliAnalysisTaskRhopPb::AliAnalysisTaskRhopPb(const char *name)
: AliAnalysisTaskSE(name),
  fESDtrackCuts(0),
  fOutputContainer(0x0),
  fMyTree(0x0),
  hStats(0x0),
  PIDQA(0x0),
  phitr_all(0x0),
  phitr_tpc(0x0),
  phitr_tof(0x0),
  phitr_pid(0x0),
  cosdip(0x0),
  fEventCuts(0x0),
  fPIDResponse(0x0),
  fPIDCombined(0x0)
				      
{
  // Output slots 
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());
}

//-----------------------------------------------------------------------------
AliAnalysisTaskRhopPb::~AliAnalysisTaskRhopPb() {
  Info("~AliAnalysisTaskRhopPb","Calling Destructor");
  if (fOutputContainer){
    delete fOutputContainer;
    fOutputContainer = 0;
  }
  if (fMyTree){
    delete fMyTree;
    fMyTree = 0;
  }
}
//________________________________________________________________________
Bool_t AliAnalysisTaskRhopPb::UserNotify(){
  static Int_t count = 0;
  count++;
  Printf("Processing %d. file: %s", count, ((TTree*) GetInputData(0))->GetCurrentFile()->GetName());
  return kTRUE;
}
//-----------------------------------------------------------------------------
void AliAnalysisTaskRhopPb::UserCreateOutputObjects()
{
  // Create histograms
  printf("!!!! Entered OutputObjects !!!!!! \n");
  // ESD histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
  }
  if(fMyTree != NULL){
    delete fMyTree;
  }

  hStats= new TH1I("hStats","hStats",15,0.,15.);
  hStats->GetXaxis()->SetBinLabel(1, "Total");
  hStats->GetXaxis()->SetBinLabel(2, "kINT7");
  hStats->GetXaxis()->SetBinLabel(3, "FistEvChunk");
  hStats->GetXaxis()->SetBinLabel(4, "IncomDAQ");
  hStats->GetXaxis()->SetBinLabel(5, "SPDBGR");
  hStats->GetXaxis()->SetBinLabel(6, "Centr");
  hStats->GetXaxis()->SetBinLabel(7, "Vrtx");
  hStats->GetXaxis()->SetBinLabel(8, "Vrtx10");
  hStats->GetXaxis()->SetBinLabel(9, "NtrGl>0");
  hStats->GetXaxis()->SetBinLabel(10, "Ntr>0");
  hStats->GetXaxis()->SetBinLabel(11, "Ntr>0");
  hStats->GetXaxis()->SetBinLabel(12, "Ntr>0");
  hStats->GetXaxis()->SetBinLabel(13, "Ntr>0");

  PIDQA= new TH2F("PIDQA","PIDQA",200,0.,20., 2000, -100., 100.);
  phitr_all = new TH2F("phitr_all","phitr_all",370, 0., 360., 200, 0., 20.);
  phitr_tpc = new TH2F("phitr_tpc","phitr_tpc",370, 0., 360., 200, 0., 20.);
  phitr_tof = new TH2F("phitr_tof","phitr_tof",370, 0., 360., 200, 0., 20.);
  phitr_pid = new TH2F("phitr_pid","phitr_pid",370, 0., 360., 200, 0., 20.);
  cosdip= new TH2F("cosdip","cosdip",40,0.,20., 500, -0.001, 0.05);

  fOutputContainer = new TList();
  fOutputContainer->SetOwner(kTRUE);
  fOutputContainer->Add(hStats);

  fOutputContainer->Add(PIDQA);
  fOutputContainer->Add(phitr_all);
  fOutputContainer->Add(phitr_tpc);
  fOutputContainer->Add(phitr_tof);
  fOutputContainer->Add(phitr_pid);
  fOutputContainer->Add(cosdip);

  fOutputContainer->ls();

  // *************************************************
  fMyTree = new TTree("MT","Tree");

  fMyTree->Branch("rnumb",          &rnumb,                "rnumb/I");   
  fMyTree->Branch("ntr",            &ntr,                  "ntr/I");
  fMyTree->Branch("vert_pr",        &vtx51,                "vert_pr[3]/F");
  fMyTree->Branch("centr",          &centr,                "centr/F");
  //Treki
  fMyTree->Branch("chTr",           &chTr,                   "chTr[ntr]/I");
  fMyTree->Branch("xyzImpTr",       &xyzImpTr,               "xyzImpTr[ntr][2]/F");
  fMyTree->Branch("momTr",          &momTr,                  "momTr[ntr][3]/F");
  fMyTree->Branch("chi2TPC",        &chi2TPC,                "chi2TPC[ntr]/F");
  fMyTree->Branch("chi2ITS",        &chi2ITS,                "chi2ITS[ntr]/F");
  fMyTree->Branch("nTPCC",          &nTPCC,                  "nTPCC[ntr]/F");
  fMyTree->Branch("ratTPCCF",       &ratTPCCF,               "ratTPCCF[ntr]/F");
  fMyTree->Branch("fracTPCS",       &fracTPCS,               "fracTPCS[ntr]/F");
  // PID
  fMyTree->Branch("pid_sigma1",       &pid_sigma1,              "pid_sigma1[ntr]/F");
  fMyTree->Branch("pid_sigma2",       &pid_sigma2,              "pid_sigma2[ntr]/F");

  // Create ESD track cut

  fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "Loose");
  //TPC
  fESDtrackCuts->SetMinNCrossedRowsTPC(60);//70 -def
  fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8); //0.8 -def
  fESDtrackCuts->SetMaxChi2PerClusterTPC(6); //4 - def
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  //ITS
  fESDtrackCuts->SetRequireITSRefit(kTRUE);
  fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  //DCA
  //7*(0.0026+0.0050/pt^1.01) - old
  fESDtrackCuts->SetMaxDCAToVertexXYPtDep("9*(0.0015+0.0050/pt^1.1)"); // 7 -def
  fESDtrackCuts->SetMaxDCAToVertexZ(2.0); // 2 -def
  fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(100); // 36 - def
  fESDtrackCuts->SetDCAToVertex2D(kFALSE);
  fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);
  fESDtrackCuts->SetPtRange(0.3, 1000.);
  fESDtrackCuts->SetEtaRange(-0.8,0.8); // normally, |eta|<0.8
/////////////////////////////////////////////////////////////////////////////

  PostData(1, fMyTree);
  PostData(2, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskRhopPb::UserExec(Option_t *)
{
  AliESDEvent *event = (AliESDEvent*) fInputEvent;

  if (!event) {
    Printf("ERROR: Could not retrieve event");
    PostData(1, fMyTree);
    PostData(2, fOutputContainer);
    return;
  }

  hStats->AddBinContent(1);

  AliAnalysisUtils utils;

  UInt_t fSelectMask= fInputHandler->IsEventSelected();
  Bool_t isINT7selected = fSelectMask& AliVEvent::kINT7;
  if (isINT7selected) hStats->AddBinContent(2);

  AliVEvent *ev = InputEvent();
  if (!fEventCuts.AcceptEvent(ev)) {
    PostData(1, fMyTree);    
    PostData(2, fOutputContainer);
    return;
  }

  hStats->AddBinContent(3);
  if (!isINT7selected) return;
  hStats->AddBinContent(4);

  // *************************************
  // Centrality
  // *************************************
  Float_t lPercentile = 300;
  AliMultSelection *MultSelection = 0x0;
  MultSelection = (AliMultSelection * ) event->FindListObject("MultSelection");

  if( !MultSelection) {
     //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
     AliWarning("AliMultSelection object not found!");
  }else{
     lPercentile = MultSelection->GetMultiplicityPercentile("V0A");
     centr = lPercentile;
  }

  if( centr < 0 || centr > 100 ) return;
  hStats->AddBinContent(6);

  ntr=0; rnumb=-9999; 
  for (int i=0; i<3; i++) { vtx51[i]=-9999; }

  rnumb=event->GetRunNumber();

  // ********************************************
  //            Vertex
  // ******************************************** 

  const AliVVertex* vtx = fEventCuts.GetPrimaryVertex(); /// Best primary vertex available
  hStats->AddBinContent(7);
  Double_t vtx5[3];

  vtx5[0] = vtx->GetX();
  vtx5[1] = vtx->GetY();
  vtx5[2] = vtx->GetZ();
  
  for (int iii=0; iii<3; iii++) {vtx51[iii]=float(vtx5[iii]);}
  if (TMath::Abs(vtx5[2]) > 10. ) return;

  hStats->AddBinContent(8);

///////////// END OF VERTEX ///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/////////////// PID ///////////////////////////////////////////////////////////
  skip_PID = 0;
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) 
    {
      Printf("Bad PID - skip it \n");
      skip_PID = 1;
    }

  Double_t probTPCTOF1[AliPID::kSPECIES]={0.};
  if ( skip_PID == 0 )
    {
      // ------- setup PIDCombined
      fPIDCombined=new AliPIDCombined;
      fPIDCombined->SetDefaultTPCPriors();
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
    }
///////////////////////////////////////////////////////////////////////////////

  Double_t  mom_pip[3], mom_pim[3], pt_pip, pt_pim, ptf_pip, ptf_pim, ptf, mf, E, rap;
  Float_t   mPiPDG  = 0.13957018;
  Int_t     ch1, ch2;   
  ntr=0; 
  
  // Initializing arrays for pid
  for (int iii=0; iii<5000; iii++) {
    for (int jjj=0; jjj<5; jjj++) {
      pid_comb[iii][jjj]=-999;
      if (jjj<4)  pid_sigma[iii][jjj]=-999; 
      pid_sigma1[iii]=-999; 
      pid_sigma2[iii]=-999;
      pid_sigma3[iii]=-999;
    }
  }

  if ((event->GetNumberOfTracks()) > 0.)   hStats->AddBinContent(9);
  // ************************************
  for (Int_t i=0; i<event->GetNumberOfTracks(); i++) 
    {
      AliESDtrack *tr1 = event->GetTrack(i);//new AliESDtrack(*event->GetTrack(i)) ;
      if (!tr1)
	{
	  Printf("ERROR: Could not receive track 1 - %d", i);
	  continue;
	}

      if ( !(fESDtrackCuts->AcceptTrack(tr1) &&  TMath::Abs(tr1->Eta())< 0.8) ) continue;

      tr1 -> GetImpactParameters(xyzImpTr[ntr][0],xyzImpTr[ntr][1]);
      tr1 -> GetConstrainedPxPyPz(mom_pip);
      for (int iii=0; iii<3; iii++) {momTr[ntr][iii]=float(mom_pip[iii]);}

      pt_pip = sqrt( pow(momTr[ntr][0],2) + pow(momTr[ntr][1],2) );
      ptf_pip = sqrt( pow(momTr[ntr][0],2) + pow(momTr[ntr][1],2) + pow(momTr[ntr][2],2) );

      phi_pip = -1;
      if ( momTr[ntr][0] != 0 )
	{
	  if ( momTr[ntr][0] > 0 && momTr[ntr][1] >= 0 ) phi = atan(momTr[ntr][1]/momTr[ntr][0])/3.14156*180;
	  if ( momTr[ntr][0] < 0 && momTr[ntr][1] >= 0 ) phi = atan(momTr[ntr][1]/momTr[ntr][0])/3.14156*180 + 180;
	  if ( momTr[ntr][0] < 0 && momTr[ntr][1] <  0 ) phi = atan(momTr[ntr][1]/momTr[ntr][0])/3.14156*180 + 180;
	  if ( momTr[ntr][0] > 0 && momTr[ntr][1] <  0 ) phi = atan(momTr[ntr][1]/momTr[ntr][0])/3.14156*180 + 360;
	  phi_pip = phi;
	}

      chTr[ntr]        = tr1->Charge();

      Int_t  nITS_tr1 = 0;
      UChar_t itsCluMap_tr1;
      itsCluMap_tr1 = tr1->GetITSClusterMap();
      for (int k = 2; k < 6; k++) if (itsCluMap_tr1 & (1 << k)) nITS_tr1 = nITS_tr1 + int(pow(10.0,k));
      Int_t nSPD_tr1 = 0;
      nSPD_tr1  = TESTBIT(tr1->GetITSClusterMap(), 0);
      nSPD_tr1 += TESTBIT(tr1->GetITSClusterMap(), 1)*10;

      nITSTr[ntr]      = nITS_tr1;
      nSPDTr[ntr]      = nSPD_tr1;

      // ====================== Quality cuts
      chi2PerClusterITS = -1;
      chi2PerClusterTPC = -1;
      nClustersITS = -1;
      nClustersTPC = -1;
      nClustersTPCShared = -1;
      fracClustersTPCShared = -1.;
      nCrossedRowsTPC = -1;
      ratioCrossedRowsOverFindableClustersTPC = 1.0;
      
      //
      nClustersTPC = tr1->GetTPCclusters(0);
      //
      nClustersITS = tr1->GetITSclusters(0);  
      if (nClustersITS!=0) 
	{
	  chi2PerClusterITS = tr1->GetITSchi2()/Float_t(nClustersITS);
	}
      //
      nClustersTPCShared = tr1->GetTPCnclsS();
      nClustersTPC = tr1->GetTPCclusters(0);
      if (nClustersTPC!=0) 
	{    
	  chi2PerClusterTPC = tr1->GetTPCchi2()/Float_t(nClustersTPC);    
	  fracClustersTPCShared = Float_t(nClustersTPCShared)/Float_t(nClustersTPC);
	}
      //
      nCrossedRowsTPC = tr1->GetTPCCrossedRows();
      //
      if (tr1->GetTPCNclsF()>0) 
	{
	  ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / tr1->GetTPCNclsF();
	}
      
      chi2TPC[ntr]     = chi2PerClusterTPC;
      chi2ITS[ntr]     = chi2PerClusterITS;  fracTPCS[ntr]    = fracClustersTPCShared;
      nTPCC[ntr]       = nCrossedRowsTPC;
      ratTPCCF[ntr]    = ratioCrossedRowsOverFindableClustersTPC;      
      nclTPCS[ntr]     = nClustersTPCShared;
      // ====================== End quality cuts
      
      if ( skip_PID == 0 )
	{
	  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  //                        PID
	  // *********************************************************
	  pid_sigma[ntr][0] = float(fPIDResponse->NumberOfSigmasTPC((AliVTrack*)tr1, AliPID::kPion));
	  pid_sigma[ntr][1] = float(fPIDResponse->NumberOfSigmasTPC((AliVTrack*)tr1, AliPID::kKaon));
	  pid_sigma[ntr][2] = float(fPIDResponse->NumberOfSigmasTPC((AliVTrack*)tr1, AliPID::kProton));    
	  pid_sigma[ntr][3] = float(fPIDResponse->NumberOfSigmasTOF((AliVTrack*)tr1, AliPID::kPion));
	  pid_sigma[ntr][4] = float(fPIDResponse->NumberOfSigmasTOF((AliVTrack*)tr1, AliPID::kKaon));  
	  pid_sigma[ntr][5] = float(fPIDResponse->NumberOfSigmasTOF((AliVTrack*)tr1, AliPID::kProton));    

	  pid_sigma1[ntr] = pid_sigma[ntr][0];
	  pid_sigma2[ntr] = pid_sigma[ntr][3];

	  TOF_match = 1;
	  if ((tr1->GetStatus() & AliESDtrack::kTOFout) == 0) TOF_match = 0;
	  if ((tr1->GetStatus() & AliESDtrack::kTIME  ) == 0) TOF_match = 0;

	  if ( TOF_match != 1 ) pid_sigma2[ntr] = -9995;

	  FillHistogram("PIDQA",sqrt( pow(momTr[ntr][0],2) + pow(momTr[ntr][1],2) + pow(momTr[ntr][2],2) ),pid_sigma1[ntr]);

	  // *********************************************************
	  UInt_t detUsed1;
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
	  detUsed1 = fPIDCombined->ComputeProbabilities((AliVTrack*)tr1, fPIDResponse, probTPCTOF1);
	  if (detUsed1 !=0) { for (int iii=0; iii<5; iii++) pid_comb[ntr][iii]=probTPCTOF1[iii];}
	  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}

      if (fabs( pid_sigma1[ntr] ) > 5.0 && fabs( pid_sigma2[ntr] ) > 5.0) continue;
      if ( ntr<4998 ) { ntr++;}
    }

  fMyTree->Fill();

  fEventCounter++;
  PostData(1, fMyTree);
  PostData(2, fOutputContainer);


}

//_____________________________________________________________________________
void AliAnalysisTaskRhopPb::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
   
}
//____________________________________________________________________________
void AliAnalysisTaskRhopPb::FillHistogram(const char * key,Double_t x)const{
  //FillHistogram
  TH1I * tmpI = dynamic_cast<TH1I*>(fOutputContainer->FindObject(key)) ;
  if(tmpI){
    tmpI->Fill(x) ;
    return ;
  }
  TH1F * tmpF = dynamic_cast<TH1F*>(fOutputContainer->FindObject(key)) ;
  if(tmpF){
    tmpF->Fill(x) ;
    return ;
  }
  TH1D * tmpD = dynamic_cast<TH1D*>(fOutputContainer->FindObject(key)) ;
  if(tmpD){
    tmpD->Fill(x) ;
    return ;
  }
  AliInfo(Form("can not find histogram <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskRhopPb::FillHistogram(const char * key,Double_t x,Double_t y)const{
  //FillHistogram
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH1F")){
    ((TH1F*)tmp)->Fill(x,y) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y) ;
    return ;
  }
  AliError(Form("Calling FillHistogram with 2 parameters for histo <%s> of type %s",key,tmp->IsA()->GetName())) ;
}

//_____________________________________________________________________________
void AliAnalysisTaskRhopPb::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
  //Fills 1D histograms with key
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y,z) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH3F")){
    ((TH3F*)tmp)->Fill(x,y,z) ;
    return ;
  }
}
//_____________________________________________________________________________
