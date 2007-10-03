#include <TSystem.h>
#include <TFile.h>
#include <TStyle.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TH2F.h>
#include <TCanvas.h>

#include <AliHMPIDHit.h>
#include <AliHMPIDCluster.h>
#include <AliHMPIDDigit.h>
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TObjArray *CreateContainer(const char *classname,TTree *pTree)
{
  TObjArray *pOA=new TObjArray(7); pOA->SetOwner();
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){
    TClonesArray *pCA=new TClonesArray(classname);
    pOA->AddAt(pCA,iCh);    
    pTree->SetBranchAddress(Form("HMPID%i",iCh),&pCA);
  }
  return pOA;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F *hHitQdc; TH2F *hHitMap[7];  
void Hits(Int_t mode,TTree *pTree=0x0)
{
  switch(mode){
    case 1:
      hHitQdc=new TH1F("HitQdc","Hit Qdc all chamber;QDC",500,0,4000);
      for(Int_t iCh=0;iCh<7;iCh++) hHitMap[iCh]=new TH2F(Form("HitMap%i",iCh),Form("Ch%i;x_{Hit};y_{Hit}",iCh),160,0,160,160,0,160);      
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
        }//hits loop      
      }//entries loop
      delete pHits;
      break;
    case 3:
      TCanvas *c1=new TCanvas("HitCan","Hits",1280,800); c1->Divide(3,3);
  
      for(Int_t iCh=0;iCh<7;iCh++){
        if(iCh==6) c1->cd(1); if(iCh==5) c1->cd(2);
        if(iCh==4) c1->cd(4); if(iCh==3) c1->cd(5); if(iCh==2) c1->cd(6);
                              if(iCh==1) c1->cd(8); if(iCh==0) c1->cd(9);
        gStyle->SetPalette(1);      
        hHitMap[iCh]->Draw("colz");
      }  
      c1->cd(3); gPad->SetLogy(); hHitQdc->SetFillColor(5);hHitQdc->Draw();
      break;
  }
}//Hits()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F *hCluEvt,*hCluChi2,*hCluFlg,*hCluSize;
void Clus(Int_t mode, TTree *pTree=0x0)
{
  switch(mode){
    case 1:
      hCluEvt=new TH1F("CluPerEvt","# clusters per event",21,-0.5,20.5);
      hCluChi2  =new TH1F("CluChi2"  ,"Chi2 "               ,1000,0,100);
      hCluFlg   =new TH1F("CluFlg"   ,"Cluster flag"        ,14,-1.5,12.5);                       hCluFlg->SetFillColor(5);
      hCluSize  =new TH1F("CluSize"  ,"Raw cluster size    ",100,0,100);
      break;
    case 2:      
      if(pTree==0) return;
      TObjArray *pLst=CreateContainer("AliHMPIDCluster",pTree); pTree->GetEntry(0);
      for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){//chambers loop
        TClonesArray *pClus=(TClonesArray *)pLst->UncheckedAt(iCh);
        hCluEvt->Fill(pClus->GetEntriesFast());
        for(Int_t iClu=0;iClu<pClus->GetEntriesFast();iClu++){//clusters loop
          AliHMPIDCluster *pClu=(AliHMPIDCluster*)pClus->UncheckedAt(iClu);
          hCluFlg->Fill(pClu->Status());  hCluChi2->Fill(pClu->Chi2());  hCluSize->Fill(pClu->Size());
        }
      }
      delete pLst;           
      break;
    case 3:
      TCanvas *c1=new TCanvas("CluComCan","Clusters in common",1280,800); c1->Divide(3,3);
      c1->cd(1); hCluEvt->SetFillColor(5);    hCluEvt->Draw();
      c1->cd(2); hCluChi2->SetFillColor(5);     hCluChi2->Draw(); 
      c1->cd(3); hCluFlg->SetFillColor(5);      hCluFlg->Draw(); 
      c1->cd(4); hCluSize->SetFillColor(5);     hCluSize->Draw(); 
      break;
  }//switch
}//Clus()
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F *hDigPcEvt,*hDigQ, *hDigChEvt;
void Digs(Int_t mode, TTree *pTree=0x0)
{
  switch(mode){
    case 1:
      hDigPcEvt=new TH1F("hDigPcEvt","PC occupancy per event",156,-1,77);
      hDigChEvt=new TH1F("hDigChEvt","Chamber occupancy per event",32,-1,7);
      hDigQ  =new TH1F("Q        ","Q                     ",3000,0,3000);
      break;
    case 2:
      if(pTree==0) return;
      TObjArray *pLst=CreateContainer("AliHMPIDDigit",pTree); pTree->GetEntry(0);
      for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){//chambers loop
        TClonesArray *pDigs=(TClonesArray *)pLst->UncheckedAt(iCh);
        hDigChEvt->Fill(iCh,pDigs->GetEntriesFast()/(48.*80.*6.));
        for(Int_t iDig=0;iDig<pDigs->GetEntriesFast();iDig++){//digits loop
          AliHMPIDDigit *pDig=(AliHMPIDDigit*)pDigs->UncheckedAt(iDig);
          hDigPcEvt->Fill(10.*iCh+pDig->Pc(),1./(48.*80.));
          hDigQ->Fill(pDig->Q());
        }
      }
      delete pLst;
      break;
    case 3:
      TCanvas *c1=new TCanvas("DigQa","Digit Check",1280,800); c1->Divide(2,2);
      c1->cd(1); hDigPcEvt->Draw();       c1->cd(2); hDigQ->Draw(); c1->cd(3); hDigChEvt->Draw();
      break;
  }//switch
}//Dig()

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
    Hits(2,th);   Clus(2,tc); Digs(2,td);//fill
    if(th==0 && td==0 && tc==0) break;
    iEvt++;
    Printf("Event %i processed",iEvt);
  }
  if(fh) Hits(3);
  if(fc) Clus(3); 
  if(fd) Digs(3);//plot
}

