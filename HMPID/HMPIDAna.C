//This script is a skeleton for HMPID analisys. It consequently reads events in a givent set of directories and provides
//resulting RANA.root file in home directory which contains all the requested hists. 
#if !defined( __CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TTree.h>
#include <TH2F.h>
#include <AliESD.h>
#include <AliESDtrack.h>
#include <TFile.h>
#endif


TH2F *pMassLen2;

void HistBook()
{
  pMassLen2=new TH2F("masslen","Particle mass versus track length;cm;GeV",200,0,600,200,0,1.2);
}
//__________________________________________________________________________________________________
void HistFill(AliESDtrack *pTrack)
{
  pMassLen2->Fill(pTrack->GetIntegratedLength(),pTrack->GetMass());  
  
}
//__________________________________________________________________________________________________
void HistOut()
{
  TCanvas *pC=new TCanvas("HMPID analisys");
  pMassLen2->Draw();

  TFile outFile("~/RANA.root","RECREATE");   pC->Write();   outFile.Close();
}
//__________________________________________________________________________________________________
void Analyse(char *sDirName)
{
//Analyse info from single directory   
  ::Info("","Tring to open from %s",sDirName);
  TFile *pFile=TFile::Open(Form("%s/AliESDs.root",sDirName));if(!pFile || !pFile->IsOpen()) return;//open AliESDs.root                                                                    
  TTree* pTree = (TTree*) pFile->Get("esdTree");             if(!pTree)                     return;//get ESD tree
                                                                 
  AliESD *pESD=new AliESD;  pTree->SetBranchAddress("ESD", &pESD);
  
  Int_t iNevents=pTree->GetEntries();   //how many events in this given directory
  ::Info("","have %i events",iNevents);
  for(Int_t iEventN=0;iEventN<iNevents;iEventN++){//ESD events loop
    pTree->GetEvent(iEventN);
    Int_t iNtracks=pESD->GetNumberOfTracks();    
    for(Int_t iTrackN=0;iTrackN<iNtracks;iTrackN++){//ESD tracks loop
      HistFill(pESD->GetTrack(iTrackN));
    }//ESD tracks loop
  }//ESD events loop
  
  delete pESD;  pFile->Close();//close AliESDs.root
}   
//__________________________________________________________________________________________________
void RichAna(Int_t iDirFirst=1,Int_t iDirLast=4)
{
  gBenchmark->Start("HMPIDanalisys"); 
  
  HistBook();
  
  for(Int_t iDirN=iDirFirst;iDirN<=iDirLast;iDirN++) Analyse(Form("Ev%04i",iDirN));  //analise for single directory
  
  HistOut();

  gBenchmark->Stop("HMPIDanalisys");
  gBenchmark->Show("HMPIDanalisys"); 
}   
