#include <TSystem.h>
#include <TFile.h>
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
      return;
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
      return;
    case 3:
      TCanvas *c1=new TCanvas("HitCan","Hits",1280,800); c1->Divide(3,3);
  
      for(Int_t iCh=0;iCh<7;iCh++){
        if(iCh==6) c1->cd(1); if(iCh==5) c1->cd(2);
        if(iCh==4) c1->cd(4); if(iCh==3) c1->cd(5); if(iCh==2) c1->cd(6);
                              if(iCh==1) c1->cd(8); if(iCh==0) c1->cd(9);
        hHitMap[iCh]->Draw();
      }  
      c1->cd(3); gPad->SetLogy(); hHitQdc->Draw();
      return;
  }
}//Hits()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F *hCluEvt,*hCluChi2,*hCluFlg,*hCluSize;
void Clus(Int_t mode, TTree *pTree=0x0)
{
  switch(mode){
    case 1:
      hCluEvt=new TH1F("CluPerEvt","# clusters per event",11,-0.5,10.5);
      hCluChi2  =new TH1F("CluChi2"  ,"Chi2 "               ,1000,0,100);
      hCluFlg   =new TH1F("CluFlg"   ,"Cluster flag"        ,14,-1.5,12.5);                       hCluFlg->SetFillColor(5);
      hCluSize  =new TH1F("CluSize"  ,"Raw cluster size    ",100,0,100);
      return;
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
      return;  
    case 3:
      TCanvas *c1=new TCanvas("CluComCan","Clusters in common",1280,800); c1->Divide(3,3);
      c1->cd(1); hCluEvt->Draw();       c1->cd(2); hCluChi2->Draw(); 
      c1->cd(3); hCluFlg->Draw();       c1->cd(4); hCluSize->Draw(); 
      return;
  }//switch
}//Clus()
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F *hDigPcEvt,*hDigQ, *hDigChEvt;
void Digs(Int_t mode, TTree *pTree=0x0)
{
  switch(mode){
    case 1:
      hDigPcEvt=new TH1F("OccPerEvt","PC occupancy per event",156,-1,77);
      hDigChEvt=new TH1F("OccPerEvt","Chamber occupancy per event",32,-1,7);
      hDigQ  =new TH1F("Q        ","Q                     ",3000,0,3000);
      return;
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
      return;
    case 3:
      TCanvas *c1=new TCanvas("DigQa","Digit Check",1280,800); c1->Divide(2,2);
      c1->cd(1); hDigPcEvt->Draw();       c1->cd(2); hDigQ->Draw(); c1->cd(3); hDigChEvt->Draw();
      return;
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
 // if(fh) Hits(3);  if(fc) Clus(3); 

if(fd) Digs(3);//plot
}

/*

void BlahQA()
{
  TFile *pFile = new TFile("Hqa.root","recreate");
  AliHMPIDCluster::DoCorrSin(kFALSE);
  AliHMPIDDigit::fSigmas=4;
  Int_t nEvts=10000;
  TLegend *lQ=new TLegend(0.5,0.5,0.9,0.9);

  TH1F *hQ7eV  =new TH1F("hQ7eV"  ,"" ,300,-50,2000);     hQ7eV  ->SetLineColor(kRed);      lQ->AddEntry(hQ7eV  ,"Ckov 7 eV");  hQ7eV->SetStats(0);
  TH1F *hQ200eV=new TH1F("hQ200eV","" ,300,-50,2000);     hQ200eV->SetLineColor(kBlack);    lQ->AddEntry(hQ200eV,"mip 200 eV");
  TH1F *hQ500eV=new TH1F("hQ500eV","" ,300,-50,2000);     hQ500eV->SetLineColor(kCyan);     lQ->AddEntry(hQ500eV,"mip 500 eV");
  TH1F *hQ900eV=new TH1F("hQ900eV","" ,300,-50,2000);     hQ900eV->SetLineColor(kGreen);    lQ->AddEntry(hQ900eV,"mip 900 eV");
  
  
  gStyle->SetOptStat(10);
  TH2F *pCluMapSi1  =new TH2F("cluMapSi1","Size 1 map"       ,1700,-10,160,1700,-10,160);
  TH2F *pCluMapLo0  =new TH2F("cluMNoLo0","Loc Max 0 map"    ,1700,-10,160,1700,-10,160);     
  TH2F *pCluMapLo1  =new TH2F("cluMapLo1","Loc Max 1 map"    ,1700,-10,160,1700,-10,160);        
  TH2F *pCluMapUnf  =new TH2F("cluMapUnf","Unfolded map"     ,1700,-10,160,1700,-10,160);         
  TH2F *pCluMapEdg  =new TH2F("cluMapEdg","On edge map"      ,1700,-10,160,1700,-10,160);         
  TH2F *pCluMapCoG  =new TH2F("cluMapCoG","CoG  map"         ,1700,-10,160,1700,-10,160);            
  TH2F *pCluMapEmp  =new TH2F("cluMapEmp","undefined-empty"  ,1700,-10,160,1700,-10,160);      
  TH2F *pCluMapNoLoc=new TH2F("cluMapNoLoc","no loc maxima"  ,1700,-10,160,1700,-10,160);      
  TH2F *pCluMapAbn  =new TH2F("cluMapAbn","abnormal fit"     ,1700,-10,160,1700,-10,160);      
  TH2F *pCluMapNot  =new TH2F("cluMapNot","Raw Clusters"     ,1700,-10,160,1700,-10,160);      
  TH2F *pCluMapMax  =new TH2F("cluMapMax","N. locs excceds"  ,1700,-10,160,1700,-10,160);      

  TH1F *hHitCluDifX = new TH1F("hHitCluDifX" ,";entries;x_{Hit}-x_{Clu} [cm]"   ,1000,-1,1);          hHitCluDifX->Sumw2();    hHitCluDifX->SetFillColor(kYellow);
//  TH2F *hHitCluDifXv= new TH2F("hHitCluDifXv",";x_{Hit};x_{Hit}-x_{Clu} [cm]"   ,500,-0.5,0.5,1000,-0.2,0.2);hHitCluDifXv->Sumw2();
  TProfile *hHitCluDifXv= new TProfile("hHitCluDifXv",";x_{Hit};x_{Hit}-x_{Clu} [cm]"   ,500,-0.5,0.5);
  TH1F *hHitCluDifY = new TH1F("hHitCluDifY" ,";entries;y_{Hit}-y_{Clu} [cm]"   ,1000,-1,1);          hHitCluDifY->Sumw2();    hHitCluDifY->SetFillColor(kYellow);
  TH2F *hHitCluDifXY= new TH2F("hHitCluDifXY",";x_{Hit}-x_{Clu};y_{Hit}-y_{Clu}",1000,-1,1,1000,-1,1);hHitCluDifXY->Sumw2();  
  TH1F *hHitCluDifQ = new TH1F("hHitCluDifQ" ,";entries;(Q_{Clu}-Q_{Hit})/Q_{Hit}" ,200 ,-200,200);   hHitCluDifQ->Sumw2();    hHitCluDifQ->SetFillColor(kYellow);
  
 
  Float_t e200=200e-9,e500=500e-9,e900=900e-9,e7=7e-9;//predefined  Eloss
  
  AliHMPIDHit hit(0,0,kProton,0,0,0);
  for(int i=0;i<5000;i++){
    hQ200eV->Fill(hit.QdcTot(e200));  hQ500eV->Fill(hit.QdcTot(e500));   hQ900eV->Fill(hit.QdcTot(e900));  hQ7eV->Fill(hit.QdcTot(e7));
  }  
  TClonesArray hits("AliHMPIDHit");  TClonesArray sdigs("AliHMPIDDigit");
  TObjArray digs(7); for(Int_t i=0;i<7;i++) digs.AddAt(new TClonesArray("AliHMPIDDigit"),i);
  TObjArray clus(7); for(Int_t i=0;i<7;i++) clus.AddAt(new TClonesArray("AliHMPIDCluster"),i);
  
    
  for(Int_t iEvt=0;iEvt<nEvts;iEvt++){//events loop
    if(iEvt%50==0)Printf("============> iEvt = %d ",iEvt);
    
    Int_t ch,pid; Float_t e,hitx,hity,hitq;
//    Int_t nHits=(type==999)?1:40;
    Int_t nHits=1;
    for(Int_t iHit=0;iHit<nHits;iHit++){//hits loop for the current event
      switch(iHit){
        case 0:  ch=0;pid=kProton;e=e200;
                 hitx=AliHMPIDParam::SizePadX()*(6+gRandom->Rndm());
                 hity=AliHMPIDParam::SizePadY()*(6+gRandom->Rndm());break; //mip ramdomly distributed in one pad
        case 1:  ch=0;pid=kProton;e=e200;hitx=0.4;hity=0.42;break; //mip in left-hand bottom coner of chamber 0
        case 2:  ch=0;pid=kProton;e=e200;hitx=0.4;hity=30  ;break; //mip on left edge of chamber 0
        case 3:  ch=0;pid=kProton;e=e200;hitx=40; hity=0.42;break; //mip on bottom edge of chamber 0
        default: ch=gRandom->Rndm()*6; pid=(gRandom->Rndm()>0.9)? kProton:50000050;
                  if(pid==kProton) 
                    e=gRandom->Rndm()*900e-9; 
                  else
                    e=5.5e-9+3e-9*gRandom->Rndm();
                  hitx=gRandom->Rndm()*AliHMPIDParam::SizeAllX(); hity=gRandom->Rndm()*AliHMPIDParam::SizeAllY();break; //random hit
      }
      new(hits[iHit]) AliHMPIDHit(ch,e,pid,iHit,hitx,hity);                          
      hitq=e;
    }//hits loop
    
            case AliHMPIDCluster::kNoLoc: pCluMapNoLoc->Fill(clux,cluy);break;
            case AliHMPIDCluster::kNoLoc: pCluMapNoLoc->Fill(clux,cluy);break;
            case AliHMPIDCluster::kAbn:   pCluMapAbn->Fill(clux,cluy);  break;
            case AliHMPIDCluster::kNot:   pCluMapNot->Fill(clux,cluy);  break;
            default:     pCluMapEmp->Fill(clux,cluy); break; 		
        }
        
        hHitCluDifX->Fill(hitx-clux); hHitCluDifY->Fill(hity-cluy); hHitCluDifXY->Fill(hitx-clux,hity-cluy); hHitCluDifQ->Fill(100*(cluq-hitq)/hitq);
        // distorsion due to feedback photons
        Int_t pc,px,py;
        AliHMPIDParam::Lors2Pad(hitx,hity,pc,px,py);
        Float_t padCenterX = AliHMPIDParam::LorsX(pc,px);
        if(pClu->Size()>1)hHitCluDifXv->Fill(hitx-padCenterX,(hitx-clux));        
        //
      }//clusters loop
    }//chambers loop
      
    hits.Delete();  sdigs.Delete();  for(int i = 0;i<7;i++){((TClonesArray*)digs.At(i))->Delete();((TClonesArray*)clus.At(i))->Delete();}
  }//events loop      
}



*/
