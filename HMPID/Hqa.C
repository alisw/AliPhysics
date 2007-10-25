#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TFile.h>
#include <TStyle.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>

#include <AliHMPIDParam.h>
#include <AliHMPIDHit.h>
#include <AliHMPIDCluster.h>
#include <AliHMPIDDigit.h>
#endif

Int_t nEntries = 0;

TObjArray *CreateContainer(const char *classname,TTree *pTree);
void Hits(Int_t mode,TTree *pTree=0x0);
void Digs(Int_t mode, TTree *pTree=0x0);
void Clus(Int_t mode, TTree *pTree=0x0);
void Summary();
  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TObjArray *CreateContainer(const char *classname,TTree *pTree)
{
  TObjArray *pOA=new TObjArray(AliHMPIDParam::kMaxCh+1);
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){
    TClonesArray *pCA=new TClonesArray(classname);
    pOA->AddAt(pCA,iCh);    
  }
  
  pOA->SetOwner(kTRUE);  
  
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){
    pTree->SetBranchAddress(Form("HMPID%i",iCh),&(*pOA)[iCh]);
  }
  return pOA;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F *hHitQdc,*hHitQdcCh[AliHMPIDParam::kMaxCh+1]; TH2F *hHitMap[AliHMPIDParam::kMaxCh+1];  
Double_t fHitMean[AliHMPIDParam::kMaxCh+1],fHitErr[AliHMPIDParam::kMaxCh+1];
void Hits(Int_t mode,TTree *pTree)
{
  switch(mode){
    case 1:
      hHitQdc=new TH1F("HitQdc","Hit Qdc all chamber;QDC",4000,0,4000);
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){
        hHitMap[iCh]=new TH2F(Form("HitMap%i",iCh),Form("Ch%i;x_{Hit};y_{Hit}",iCh),160,0,160,160,0,160);      
        hHitQdcCh[iCh]=new TH1F(Form("HMPID%i",iCh),Form("Charge for HMPID%i",iCh),4000,0,4000);
      }
      break;
    case 2:
      if(pTree==0) return;
      TClonesArray *pHits=new TClonesArray("AliHMPIDHit");  pTree->SetBranchAddress("HMPID",&pHits);  
      for(Int_t iEnt=0;iEnt<pTree->GetEntriesFast();iEnt++){//entries loop
        pTree->GetEntry(iEnt);
        for(Int_t iHit=0;iHit<pHits->GetEntriesFast();iHit++){//hits loop
          AliHMPIDHit *pHit = (AliHMPIDHit*)pHits->UncheckedAt(iHit);
          hHitMap[pHit->Ch()]->Fill(pHit->LorsX(),pHit->LorsY());
          hHitQdc->Fill(pHit->Q());
          hHitQdcCh[pHit->Ch()]->Fill(pHit->Q());
        }//hits loop      
      }//entries loop
      delete pHits;
      break;
    case 3:
      TCanvas *c1=new TCanvas("HitCan","Hits",1280,800); c1->Divide(3,3);
  
      for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){
        if(iCh==6) c1->cd(1); if(iCh==5) c1->cd(2);
        if(iCh==4) c1->cd(4); if(iCh==3) c1->cd(5); if(iCh==2) c1->cd(6);
                              if(iCh==1) c1->cd(8); if(iCh==0) c1->cd(9);
        gStyle->SetPalette(1);      
        hHitMap[iCh]->Draw("colz");
      }  
      for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){
        fHitMean[iCh] = 0;
        fHitErr[iCh]  = 0;
        if((Int_t)hHitQdcCh[iCh]->GetEntries()<nEntries) continue;
        c1->cd(3);hHitQdcCh[iCh]->Fit("expo","Q");
        TF1 *funcfit = (TF1*)hHitQdcCh[iCh]->FindObject("expo");
        fHitMean[iCh] = funcfit->GetParameter(1);
        fHitErr[iCh]  = funcfit->GetParError(1);
      }
      TPad *pad = (TPad*)c1->cd(3); hHitQdc->SetFillColor(5); pad->SetLogy(); hHitQdc->Fit("expo","Q");
      break;
  }
}//Hits()
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F     *hDigQ;
TProfile *hDigHighQ,*hDigChEvt;
void Digs(Int_t mode, TTree *pTree)
{
  switch(mode){
    case 1:
      hDigHighQ=new TProfile("hDigHighQ","Highest charge in chamber  ",AliHMPIDParam::kMaxCh+1,AliHMPIDParam::kMinCh,AliHMPIDParam::kMaxCh+1);
      hDigChEvt=new TProfile("hDigChEvt","Chamber occupancy per event",AliHMPIDParam::kMaxCh+1,AliHMPIDParam::kMinCh,AliHMPIDParam::kMaxCh+1);
      hDigQ    =new TH1F("hDigQ        ","Charge of digits (ADC)     ",3000,0,3000);
      break;
    case 2:
      if(pTree==0) return;
      TObjArray *pLst=CreateContainer("AliHMPIDDigit",pTree); pTree->GetEntry(0);
      for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){//chambers loop
        TClonesArray *pDigs=(TClonesArray *)pLst->UncheckedAt(iCh);
        hDigChEvt->Fill(iCh,pDigs->GetEntriesFast()/(48.*80.*6.)*100.);
        Double_t highQ = 0;
        for(Int_t iDig=0;iDig<pDigs->GetEntriesFast();iDig++){//digits loop
          AliHMPIDDigit *pDig=(AliHMPIDDigit*)pDigs->UncheckedAt(iDig);
          hDigQ->Fill(pDig->Q());
          if(pDig->Q()>highQ) highQ = pDig->Q();
        }
        hDigHighQ->Fill(iCh,highQ);
      }
      delete pLst;
      break;
    case 3:
      TCanvas *c1=new TCanvas("DigQa","Digit Check",1280,800); c1->Divide(2,2);
      gStyle->SetOptFit(1);
      TPad *pad = (TPad*)c1->cd(1); pad->SetLogy(); hDigQ->Fit("expo","QN"); 
      c1->cd(2); hDigHighQ->Draw();
      c1->cd(3); hDigChEvt->Draw();
      break;
  }//switch
}//Dig()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F *hCluEvt,*hCluChi2,*hCluFlg,*hCluSize,*hCluQ;
void Clus(Int_t mode, TTree *pTree)
{
  switch(mode){
    case 1:
      hCluEvt=new TH1F("CluPerEvt","Cluster multiplicity"   ,100,0,100);
      hCluChi2  =new TH1F("CluChi2"  ,"Chi2 "               ,1000,0,100);
      hCluFlg   =new TH1F("CluFlg"   ,"Cluster flag"        ,14,-1.5,12.5);                       hCluFlg->SetFillColor(5);
      hCluSize  =new TH1F("CluSize"  ,"Cluster size        ",100,0,100);
      hCluQ     =new TH1F("CluQ"     ,"Cluster charge (ADC)",1000,0,5000);
      break;
    case 2:      
      if(pTree==0) return;

      TObjArray *pLst=CreateContainer("AliHMPIDCluster",pTree); pTree->GetEntry(0);
      for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){//chambers loop
        TClonesArray *pClus=(TClonesArray *)pLst->UncheckedAt(iCh);
        hCluEvt->Fill(pClus->GetEntriesFast());
        for(Int_t iClu=0;iClu<pClus->GetEntriesFast();iClu++){//clusters loop
          AliHMPIDCluster *pClu=(AliHMPIDCluster*)pClus->UncheckedAt(iClu);
          hCluFlg->Fill(pClu->Status());
          hCluChi2->Fill(pClu->Chi2());
          hCluSize->Fill(pClu->Size());
          hCluQ->Fill(pClu->Q());
        }
      }
      delete pLst;           
      break;
    case 3:
      TCanvas *c1=new TCanvas("CluComCan","Clusters in common",1280,800); c1->Divide(3,3);
      c1->cd(1); hCluEvt->SetFillColor(5);      hCluEvt->Draw();
      c1->cd(2); hCluChi2->SetFillColor(5);     hCluChi2->Draw(); 
      c1->cd(3); hCluFlg->SetFillColor(5);      hCluFlg->Draw(); 
      c1->cd(4); hCluSize->SetFillColor(5);     hCluSize->Draw(); 
      TPad *pad = (TPad*)c1->cd(5); hCluQ->SetFillColor(5); pad->SetLogy(); hCluQ->Draw(); 
      break;
  }//switch
}//Clus()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Hqa()
{
  TFile *fh=0; if(gSystem->IsFileInIncludePath("HMPID.Hits.root"))      fh=TFile::Open("HMPID.Hits.root"     ,"read");if(fh->IsZombie()) fh=0;
  TFile *fd=0; if(gSystem->IsFileInIncludePath("HMPID.Digits.root"))    fd=TFile::Open("HMPID.Digits.root"   ,"read");if(fd->IsZombie()) fd=0;
  TFile *fc=0; if(gSystem->IsFileInIncludePath("HMPID.RecPoints.root")) fc=TFile::Open("HMPID.RecPoints.root","read");if(fc->IsZombie()) fc=0;
  if(fh==0 && fd==0 && fc==0){Printf("Nothing to do!"); return;}
  if(fh) Hits(1); if(fc) Clus(1);  if(fd) Digs(1);//book
  Int_t iEvt=0;
  while(1){
    TTree *th=0; if(fh) th=(TTree*)fh->Get(Form("Event%i/TreeH",iEvt));
    TTree *td=0; if(fd) td=(TTree*)fd->Get(Form("Event%i/TreeD",iEvt));
    TTree *tc=0; if(fc) tc=(TTree*)fc->Get(Form("Event%i/TreeR",iEvt));
    
    Hits(2,th); Digs(2,td); Clus(2,tc); //fill
    if(th==0 && td==0 && tc==0) break;
    iEvt++;
    if(!(iEvt%50)) Printf("Event %i processed",iEvt);
  }
  
  nEntries = iEvt;
  
  if(fd) Clus(3);//plot everything
  if(fc) Digs(3); 
  if(fh) Hits(3);
  Summary();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Summary()
{
  //info for hits...
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){
    Printf(" #################### Summary for HMPID %i#################### ",iCh);
    //info for hits...
    Printf("-----Summary of Hits----- ");
    if(fHitMean[iCh]==0) {
      Printf("gain %5.2f +/- %5.2f",fHitMean[iCh],fHitErr[iCh]);
    } else {   
      Double_t gain = 1./TMath::Abs(fHitMean[iCh]);
      Double_t errgain = gain*gain*fHitErr[iCh];
      Printf("gain %5.2f +/- %5.2f",gain,errgain);
    }
    //info for digits...
    Printf("-----Summary of Digits-----");
    Printf(" Chamber %d with occupancy (%) %5.2f +/- %5.2f",iCh,hDigChEvt->GetBinContent(iCh+1),hDigChEvt->GetBinError(iCh+1));
    Printf(" Chamber %d with higest Q (ADC) %7.0f +/- %7.0f",iCh,hDigHighQ->GetBinContent(iCh+1),hDigHighQ->GetBinError(iCh+1));
  }
}
