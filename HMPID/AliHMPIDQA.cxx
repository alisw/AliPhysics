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

#include <AliESD.h> 
#include <AliLog.h>
#include <AliPID.h>

#include <AliAnalysisTask.h>      //qa()
#include <AliAnalysisManager.h> //qa()
#include <TBenchmark.h>         //qa()
#include <TProof.h>             //qa()

class AliHMPIDQA : public AliAnalysisTask {

public:
           AliHMPIDQA(Int_t mode) ;
  virtual ~AliHMPIDQA() ;
   
  virtual void Exec(Option_t * opt = "") ;
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "");

private:
  TTree   * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESD  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer; //output data container

  TH2F * fhHMPIDCkovP;  //
  TH2F * fhHMPIDMipXY;  //
  TH2F * fhHMPIDDifXY;  //
  TH2F * fhHMPIDSigP;   //
  TH1F * fhHMPIDProb[5];//
  
  ClassDef(AliHMPIDQA,0); // 
};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDQA::AliHMPIDQA(Int_t mode):AliAnalysisTask("HmpidQaTask",""),  
  fChain(0),
  fESD(0), 
  fhHMPIDCkovP(0),
  fhHMPIDMipXY(0),
  fhHMPIDDifXY(0),
  fhHMPIDSigP(0)
{
// Constructor.
  DefineInput (0,TChain::Class());           // Input slot #0 works with an Ntuple  
  DefineOutput(0,TObjArray::Class()) ;   // Output slot #0 writes into a TH1 container

  for(Int_t i=0;i<5;i++) fhHMPIDProb[i]=0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDQA::~AliHMPIDQA()
{
  // dtor
  fOutputContainer->Clear() ;   delete fOutputContainer ; 
  
  delete fhHMPIDCkovP ;  
  delete fhHMPIDMipXY ;  
  delete fhHMPIDDifXY ;  
  delete fhHMPIDSigP ;   
  delete [] fhHMPIDProb ;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQA::ConnectInputData(const Option_t*)
{
//Virtual from AliAnalysisTask invoked by AliAnalysisTask::CheckNotify() which in turn invoked by AliAnalysisDataContainer::SetData()
  fChain = dynamic_cast<TChain *>(GetInputData(0)) ;
  if (!fChain) {AliError(Form("Input 0 for %s not found\n", GetName())); return;}
  
  // One should first check if the branch address was taken by some other task
  char ** address = (char **)GetBranchAddress(0, "ESD");
  if (address) {
    fESD = (AliESD*)(*address);
  } else {
    fESD = new AliESD();
    SetBranchAddress(0, "ESD", &fESD);
    fChain->SetBranchStatus("*", 1);
    fChain->SetBranchStatus("fTracks.*", 1);
  }
}//ConnectInputData()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQA::CreateOutputObjects()
{  
  fhHMPIDCkovP    = new TH2F("CkovP" , "#theta_{c}, [rad];P, [GeV]", 150,   0,  7  ,100, -3, 1); 
  fhHMPIDSigP     = new TH2F("SigP"  ,"#sigma_{#theta_c}"          , 150,   0,  7  ,100, 0, 1e20);
  fhHMPIDMipXY    = new TH2F("MipXY" ,"mip position"               , 260,   0,130  ,252,0,126); 
  fhHMPIDDifXY    = new TH2F("DifXY" ,"diff"                       , 260, -10, 10  ,252,-10,10); 
  
  fhHMPIDProb[0] = new TH1F("PidE" ,"PID: e yellow #mu magenta"  ,100,0,1);   fhHMPIDProb[0]->SetLineColor(kYellow);
  fhHMPIDProb[1] = new TH1F("PidMu","pid of #mu"                 ,100,0,1);   fhHMPIDProb[1]->SetLineColor(kMagenta);
  fhHMPIDProb[2] = new TH1F("PidPi","PID: #pi red K green p blue",100,0,1);   fhHMPIDProb[2]->SetLineColor(kRed);
  fhHMPIDProb[3] = new TH1F("PidK" ,"pid of K"                   ,100,0,1);   fhHMPIDProb[3]->SetLineColor(kGreen);
  fhHMPIDProb[4] = new TH1F("PidP" ,"pid of p"                   ,100,0,1);   fhHMPIDProb[4]->SetLineColor(kBlue);
 

  // create output container
  fOutputContainer = new TObjArray(9) ;  fOutputContainer->SetName(GetName()) ; 

  fOutputContainer->AddAt(fhHMPIDCkovP,      0) ; 
  fOutputContainer->AddAt(fhHMPIDSigP,       1) ; 
  fOutputContainer->AddAt(fhHMPIDMipXY,      2) ; 
  fOutputContainer->AddAt(fhHMPIDDifXY,      3) ; 
  fOutputContainer->AddAt(fhHMPIDProb[0],    4) ; 
  fOutputContainer->AddAt(fhHMPIDProb[1],    5) ; 
  fOutputContainer->AddAt(fhHMPIDProb[2],    6) ; 
  fOutputContainer->AddAt(fhHMPIDProb[3],    7) ; 
  fOutputContainer->AddAt(fhHMPIDProb[4],    8) ; 
}//CreateOutputObjects()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQA::Exec(Option_t *) 
{
// Virtual from TTask. 
// Invoked by AliAnalysisManager::StartAnalysis()->AliAnalysisManager::ExecAnalysis()->TTask::ExecuteTask() in case of mgr->StartAnalysis("local")
// Invoked by AliAnalysisSelector::Process()->AliAnalysisManager::ExecAnalysis()->TTask::ExecuteTask() in case of mgr->StartAnalysis("local")
    
//  Long64_t entry = fChain->GetReadEntry() ;
  
  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  for(Int_t iTrk = 0 ; iTrk < fESD->GetNumberOfTracks() ; iTrk++){
    AliESDtrack *pTrk = fESD->GetTrack(iTrk) ;

    fhHMPIDCkovP->Fill( pTrk->GetP(), pTrk->GetHMPIDsignal() ) ; 
    fhHMPIDSigP ->Fill( pTrk->GetP(), TMath::Sqrt(pTrk->GetHMPIDchi2()) ) ;
     
//     Float_t xm,ym; Int_t q,np;  pTrk->GetHMPIDmip(xm,ym,q,np);  fMipXY->Fill(xm,ym); //mip info
//     Float_t xd,yd,th,ph;        pTrk->GetHMPIDtrk(xd,yd,th,ph); fDifXY->Fill(xd,yd); //track info 
     
    Double_t pid[5] ;      pTrk->GetHMPIDpid(pid) ; 
    for(Int_t i = 0 ; i < 5 ; i++)       fhHMPIDProb[i]->Fill(pid[i]) ;
  }//tracks loop 
       
  PostData(0,fOutputContainer);
}//Exec()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDQA::Terminate(Option_t *)
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
  
  Float_t n = 1.292 ; //mean freon ref idx 
  AliPID ppp ;                 
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
  gBenchmark->Start("HMPID QA");
  
  TChain* chain =new TChain("esdTree");  
  AliAnalysisManager *mgr=new AliAnalysisManager("FunnyName");
  
  AliAnalysisTask *qa=new AliHMPIDQA(mode);
  qa->ConnectInput (0,mgr->CreateContainer("EsdChain",TChain::Class()   ,AliAnalysisManager::kInputContainer));
  qa->ConnectOutput(0,mgr->CreateContainer("HistLst",TObjArray::Class(),AliAnalysisManager::kOutputContainer));
    
  mgr->AddTask(qa);
  if(!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  
  switch(mode){
    case 0:   chain->Add("AliESDs.root");
              mgr->StartAnalysis("local");  
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

