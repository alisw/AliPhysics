/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h> 
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h> 
#include <TROOT.h>
#include <TVector3.h> 

#include <AliESDEvent.h> 
#include <AliLog.h>
#include <AliPID.h>

#include <AliAnalysisTask.h>    //qa()
#include <AliAnalysisManager.h> //qa()
#include <TBenchmark.h>         //qa()
#include <TProof.h>             //qa()

class AliHMPIDQaEsd : public AliAnalysisTask {

public:
           AliHMPIDQaEsd() ;
  virtual ~AliHMPIDQaEsd() ;
   
  virtual void Exec(Option_t * opt = "") ;
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "");

private:
  TTree   * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESDEvent  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer; //output data container

  ClassDef(AliHMPIDQaEsd,0); // 
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDQaEsd::AliHMPIDQaEsd():AliAnalysisTask("HmpidQaTask",""), fChain(0), fESD(0)
{
// Constructor.
  DefineInput (0,TChain::Class());           // Input slot #0 works with an Ntuple  
  DefineOutput(0,TObjArray::Class()) ;   // Output slot #0 writes into a TH1 container
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDQaEsd::~AliHMPIDQaEsd()
{
  // dtor
  fOutputContainer->Clear() ;   delete fOutputContainer ; 
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQaEsd::ConnectInputData(const Option_t*)
{
//Virtual from AliAnalysisTask invoked by AliAnalysisTask::CheckNotify() which in turn invoked by AliAnalysisDataContainer::SetData()
  fChain = dynamic_cast<TChain *>(GetInputData(0)) ;
  if (!fChain) {AliError(Form("Input 0 for %s not found\n", GetName())); return;}
  
  // One should first check if the branch address was taken by some other task
  char ** address = (char **)GetBranchAddress(0, "ESD");
  if (address) {
    fESD = (AliESDEvent*)(*address);
  } else {
    fESD = new AliESDEvent();
    fESD->ReadFromTree(fChain);                                                                   //clm: new ESD access works for local, need to test it for PROOF!
    //SetBranchAddress(0, "esdTree", &fESD);
    //fChain->SetBranchStatus("*", 1);
    //fChain->SetBranchStatus("fTracks.*", 1);
  }
}//ConnectInputData()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQaEsd::CreateOutputObjects()
{  
 

  // create output container
  fOutputContainer = new TObjArray(9) ;  fOutputContainer->SetName(GetName()) ; 

  fOutputContainer->AddAt(new TH2F("CkovP" , "#theta_{c}, [rad];P, [GeV]"   , 150,   0,  7  ,100, 0, 1)  ,      0) ; 
  fOutputContainer->AddAt(new TH2F("SigP"  ,"#sigma_{#theta_c} [mrad];[GeV]", 150,   0,  7  ,100, 0, 1)  ,      1) ; 
  fOutputContainer->AddAt(new TH2F("MipXY" ,"mip position"                  , 260,   0,130  ,252, 0,126) ,      2) ; 
  fOutputContainer->AddAt(new TH2F("DifXY" ,"diff"                          , 200, -10, 10  ,200,-10,10) ,      3) ; 
  fOutputContainer->AddAt(new TH1F("PidE" ,"PID: e yellow #mu magenta"  ,100,0,1)                        ,      4) ; 
  fOutputContainer->AddAt(new TH1F("PidMu","pid of #mu"                 ,100,0,1)                        ,      5) ; 
  fOutputContainer->AddAt(new TH1F("PidPi","PID: #pi red K green p blue",100,0,1)                        ,      6) ; 
  fOutputContainer->AddAt(new TH1F("PidK" ,"pid of K"                   ,100,0,1)                        ,      7) ; 
  fOutputContainer->AddAt(new TH1F("PidP" ,"pid of p"                   ,100,0,1)                        ,      8) ; 
  //options for drawing

}//CreateOutputObjects()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQaEsd::Exec(Option_t *) 
{
// Virtual from TTask. 
// Invoked by AliAnalysisManager::StartAnalysis()->AliAnalysisManager::ExecAnalysis()->TTask::ExecuteTask() in case of mgr->StartAnalysis("local")
// Invoked by AliAnalysisSelector::Process()->AliAnalysisManager::ExecAnalysis()->TTask::ExecuteTask() in case of mgr->StartAnalysis("local")
    
  fChain->GetReadEntry() ;
  
  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  for(Int_t iTrk = 0 ; iTrk < fESD->GetNumberOfTracks() ; iTrk++){
    AliESDtrack *pTrk = fESD->GetTrack(iTrk) ;

    ((TH2F*)fOutputContainer->At(0))->Fill( pTrk->GetP(), pTrk->GetHMPIDsignal() ) ; 
    ((TH2F*)fOutputContainer->At(1))->Fill( pTrk->GetP(), TMath::Sqrt(pTrk->GetHMPIDchi2()) ) ;
     
     Float_t xm,ym; Int_t q,np;  pTrk->GetHMPIDmip(xm,ym,q,np);                       //mip info
     ((TH2F*)fOutputContainer->At(2))->Fill(xm,ym);
     Float_t xRad,yRad,th,ph;        pTrk->GetHMPIDtrk(xRad,yRad,th,ph);              //track info at the middle of the radiator
     Float_t xPc = xRad+9.25*TMath::Tan(th)*TMath::Cos(ph); // temporar: linear extrapol (B=0!)
     Float_t yPc = yRad+9.25*TMath::Tan(th)*TMath::Sin(ph); // temporar:          "
     ((TH2F*)fOutputContainer->At(3))->Fill(xm-xPc,ym-yPc); //track info 
     
    Double_t pid[5] ;      pTrk->GetHMPIDpid(pid) ; 
    for(Int_t i = 0 ; i < 5 ; i++)       ((TH1F*)fOutputContainer->At(4+i))->Fill(pid[i]) ;
  }//tracks loop 
       
  PostData(0,fOutputContainer);
}//Exec()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQaEsd::Terminate(Option_t *)
{
//Virual from Processing when the event loop is ended
  TObjArray *out=(TObjArray*)GetOutputData(0);
  
  TH2F *hAngP    = (TH2F*)out->At(0);
  TH2F *hErrP    = (TH2F*)out->At(1);
  TH2F *hMipXY   = (TH2F*)out->At(2);
  TH2F *hDifXY   = (TH2F*)out->At(3);
  TH1F *hProE    = (TH1F*)out->At(4);
  TH1F *hProMu   = (TH1F*)out->At(5);
  TH1F *hProPi   = (TH1F*)out->At(6);
  TH1F *hProK    = (TH1F*)out->At(7);
  TH1F *hProP    = (TH1F*)out->At(8);
  
  hProE ->SetLineColor(kYellow);
  hProMu->SetLineColor(kMagenta);
  hProPi->SetLineColor(kRed);
  hProK ->SetLineColor(kGreen);
  hProP ->SetLineColor(kBlue);
  
  Float_t n = 1.292 ; //mean freon ref idx 
  AliPID dummy ;      //just to initialize AliPID to get the correct particle masses
  TF1* funPi = new TF1("RiPiTheo", "acos(sqrt(x*x+[0]*[0])/(x*[1]))", 1.2, 7); funPi->SetLineWidth(1);   funPi->SetParameter(1,n) ; 

                                                funPi->SetLineColor(kRed);     funPi->SetParameter(0,AliPID::ParticleMass(AliPID::kPion));    
  TF1* funK=static_cast<TF1*>(funPi->Clone()) ; funK ->SetLineColor(kGreen) ;   funK->SetParameter(0,AliPID::ParticleMass(AliPID::kKaon)) ; 
  TF1* funP=static_cast<TF1*>(funPi->Clone()) ; funP ->SetLineColor(kBlue) ;    funP->SetParameter(0,AliPID::ParticleMass(AliPID::kProton)) ; 

  TCanvas * can = new TCanvas("HmpidCanvas","HMPID ESD Test"); can->SetFillColor(10) ;   can->SetHighLightColor(10) ;   can->Divide(3,2) ;

  can->cd(1);hAngP->Draw();funPi->Draw("same");funK->Draw("same");funP->Draw("same"); can->cd(2);hMipXY->Draw();  can->cd(3);hProE->Draw();hProMu->Draw("same");  
  can->cd(4);hErrP->Draw();                                                           can->cd(5);hDifXY->Draw();  can->cd(6);hProPi->Draw();hProK->Draw("same");hProP->Draw("same");
}//Terminate()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void qa(Int_t mode=0)
{
  
  /*
  AliAODHandler* aodHandler   = new AliAODHandler();
    mgr->SetEventHandler(aodHandler);
  */
  
  gBenchmark->Start("HMPID QA");
  
  TChain* chain =new TChain("esdTree");  
  AliAnalysisManager *mgr=new AliAnalysisManager("FunnyName");                                                   //clm: 
  //AliAODHandler* aodHandler   = new AliAODHandler();
  //mgr->SetEventHandler(aodHandler);
    
  AliAnalysisTask *qa=new AliHMPIDQaEsd();
  qa->ConnectInput (0,mgr->CreateContainer("EsdChain",TChain::Class()   ,AliAnalysisManager::kInputContainer));
  qa->ConnectOutput(0,mgr->CreateContainer("HistLst",TObjArray::Class(),AliAnalysisManager::kOutputContainer));
    
  mgr->AddTask(qa);
  if(!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  
  switch(mode){
    case 0:   chain->Add("AliESDs.root");
              mgr->StartAnalysis("local",chain);  
              break;
              
    case 1:   if(TProof::Open("proof://hmpid@lxb6046.cern.ch")==0x0) return; 
              gProof->UploadPackage("ESD.par"); gProof->EnablePackage("ESD");
              gProof->UploadPackage("ANALYSIS.par"); gProof->EnablePackage("ANALYSIS");                
              mgr->StartAnalysis("proof",chain);  
              break;
              
    case 2:   mgr->StartAnalysis("grid" ,chain);  
              break;
  }
  gBenchmark->Show("HMPID QA");
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

