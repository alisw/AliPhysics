void Hqa()
{
  gROOT->LoadMacro("Hdisp.C");
  HitQA();
}//Hqa()
void HitQA()
{
  Int_t nEvts=100;
  TLegend *lQ=new TLegend(0.5,0.5,0.9,0.9);

  TH1F *hQ7eV  =new TH1F("hQ7eV"  ,"" ,300,-50,2000);     hQ7eV  ->SetLineColor(kRed);      lQ->AddEntry(hQ7eV  ,"Ckov 7 eV");  hQ7eV->SetStats(0);
  TH1F *hQ200eV=new TH1F("hQ200eV","" ,300,-50,2000);     hQ200eV->SetLineColor(kBlack);    lQ->AddEntry(hQ200eV,"mip 200 eV");
  TH1F *hQ500eV=new TH1F("hQ500eV","" ,300,-50,2000);     hQ500eV->SetLineColor(kCyan);     lQ->AddEntry(hQ500eV,"mip 500 eV");
  TH1F *hQ900eV=new TH1F("hQ900eV","" ,300,-50,2000);     hQ900eV->SetLineColor(kGreen);    lQ->AddEntry(hQ900eV,"mip 900 eV");
  
  TH1F *hCluPerEvt=new TH1F("hCluPerEvt","# clusters per event",11,-0.5,10.5);
  TH1F *hCluChi2  =new TH1F("hChi2","Chi2 ",1000,0,100);
  TH1F *hCluFlg   =new TH1F("hCluFlg","Cluster flag",14,-1.5,12.5);                       hCluFlg->SetFillColor(5);
  TH1F *hCluRawSize= new TH1F("hCluRawSize","Raw cluster size ",100,0,100);
  
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

  TH1F *hHitCluDifX = new TH1F("hHitCluDifX" ,";entries;x_{Hit}-x_{Clu} [cm]"   ,2000,-2,2);          hHitCluDifX->Sumw2();    hHitCluDifX->SetFillColor(kYellow);
  TH1F *hHitCluDifY = new TH1F("hHitCluDifY" ,";entries;y_{Hit}-y_{Clu} [cm]"   ,2000,-2,2);          hHitCluDifY->Sumw2();    hHitCluDifY->SetFillColor(kYellow);
  TH2F *hHitCluDifXY= new TH2F("hHitCluDifXY",";x_{Hit}-x_{Clu};y_{Hit}-y_{Clu}",2000,-2,2,2000,-2,2);hHitCluDifXY->Sumw2();  
  TH1F *hHitCluDifQ = new TH1F("hHitCluDifQ" ,";entries;(Q_{Clu}-Q_{Hit})/Q_{Hit}" ,200 ,-200,200);   hHitCluDifQ->Sumw2();    hHitCluDifQ->SetFillColor(kYellow);
  
  TH2F *hHitMap= new TH2F("hHitMap",";x_{Hit};y_{Hit}",1700,-10,160,1700,-10,160);  
 
  Float_t e200=200e-9,e500=500e-9,e900=900e-9,e7=7e-9;//predefined  Eloss
  
  AliHMPIDHit hit(0,0,kProton,0,0,0);
  for(int i=0;i<5000;i++){
    hQ200eV->Fill(hit.QdcTot(e200));  hQ500eV->Fill(hit.QdcTot(e500));   hQ900eV->Fill(hit.QdcTot(e900));  hQ7eV->Fill(hit.QdcTot(e7));
  }  
  TClonesArray hits("AliHMPIDHit");  TClonesArray sdigs("AliHMPIDDigit");
  TObjArray digs(7); for(Int_t i=0;i<7;i++) digs.AddAt(new TClonesArray("AliHMPIDDigit"),i);
  TObjArray clus(7); for(Int_t i=0;i<7;i++) clus.AddAt(new TClonesArray("AliHMPIDCluster"),i);
  
    
  for(Int_t iEvt=0;iEvt<nEvts;iEvt++){//events loop
    if(iEvt%500==0)Printf("============> iEvt = %d ",iEvt);
    
    Int_t ch,pid; Float_t e,hitx,hity,hitq;
//    Int_t nHits=(type==999)?1:40;
    Int_t nHits=1;
    for(Int_t iHit=0;iHit<nHits;iHit++){//hits loop for the current event
      switch(iHit){
        case 0:  ch=0;pid=kProton;e=e200;hitx=16.0+gRandom->Rndm()*0.8;hity= 16.8+gRandom->Rndm()*0.84;break; //mip ramdomly distributed in one pad in the middle
        case 1:  ch=0;pid=kProton;e=e200;hitx=0.4;hity=0.42;break; //mip in left-hand bottom coner of chamber 0
        case 2:  ch=0;pid=kProton;e=e200;hitx=0.4;hity=30  ;break; //mip on left edge of chamber 0
        case 3:  ch=0;pid=kProton;e=e200;hitx=40; hity=0.42;break; //mip on bottom edge of chamber 0
        default: ch=gRandom->Rndm()*6; pid=(gRandom->Rndm()>0.9)? kProton:50000050;
                  if(pid==kProton) 
                    e=gRandom->Rndm()*900e-9; 
                  else
                    e=5.5e-9+3e-9*gRandom->Rndm();
                  hitx=gRandom->Rndm()*AliHMPIDDigit::SizeAllX(); hity=gRandom->Rndm()*AliHMPIDDigit::SizeAllY();break; //random hit
      }
      new(hits[iHit]) AliHMPIDHit(ch,e,pid,iHit,hitx,hity);                          
      hitq=e;
    }//hits loop
    
    AliHMPIDv1::Hit2Sdi(&hits,&sdigs);
    AliHMPIDDigitizer::DoNoise(kFALSE);
    AliHMPIDDigitizer::Sdi2Dig(&sdigs,&digs);
    AliHMPIDReconstructor::Dig2Clu(&digs,&clus,kFALSE);

// From here normal procedure for QA

    for(Int_t iHit=0;iHit<hits.GetEntriesFast();iHit++) {
      AliHMPIDHit *pHit = (AliHMPIDHit*)hits.UncheckedAt(iHit);
      hHitMap->Fill(pHit->LorsX(),pHit->LorsY());
    }
    Int_t kMaxCh=(nHits==1)?0:AliHMPIDDigit::kMaxCh;
    for(Int_t iCh=AliHMPIDDigit::kMinCh;iCh<=kMaxCh;iCh++){//chambers loop
      TClonesArray *pDigs=(TClonesArray *)digs.UncheckedAt(iCh);
      TClonesArray *pClus=(TClonesArray *)clus.UncheckedAt(iCh);
      
//      if(pClus->GetEntriesFast()>nHits) {hits.Print();}
//      if(pClus->GetEntriesFast()>nHits) {pDigs->Print();Printf("----S D I G I T S-------------");}
//      if(pClus->GetEntriesFast()>nHits) {sdigs.Print();Printf("-------------------------------");}
      hCluPerEvt->Fill(pClus->GetEntriesFast());
      for(Int_t iClu=0;iClu<pClus->GetEntriesFast();iClu++){//clusters loop
        AliHMPIDCluster *pClu=(AliHMPIDCluster*)pClus->UncheckedAt(iClu);
        Float_t clux=pClu->X(); Float_t cluy=pClu->Y(); Float_t cluq=pClu->Q();
        hCluFlg->Fill(pClu->Status());
        hCluChi2->Fill(pClu->Chi2());
        hCluRawSize->Fill(pClu->Size());
	 
        switch(pClu->Status()){
/*
    case        kFrm  : status="formed        "   ;break;
    case        kUnf  : status="unfolded (fit)"   ;break;
    case        kCoG  : status="coged         "   ;break;
    case        kLo1  : status="locmax 1 (fit)"   ;break;
    case        kMax  : status="exceeded (cog)"   ;break;
    case        kNot  : status="not done (cog)"   ;break;
    case        kEmp  : status="empty         "   ;break;
    case        kEdg  : status="edge     (fit)"   ;break;
    case 	kSi1  : status="size 1   (cog)"   ;break;
    case 	kNoLoc: status="no LocMax(fit)"   ;break;
    case 	kAbn  : status="Abnormal fit  "   ;break;
*/
            case AliHMPIDCluster::kSi1:   pCluMapSi1->Fill(clux,cluy); break;
            case AliHMPIDCluster::kLo1:   pCluMapLo1->Fill(clux,cluy); break;
            case AliHMPIDCluster::kUnf:   pCluMapUnf->Fill(clux,cluy); break; 
            case AliHMPIDCluster::kMax:   pCluMapMax->Fill(clux,cluy); break;
            case AliHMPIDCluster::kEdg:   pCluMapEdg->Fill(clux,cluy); break;
            case AliHMPIDCluster::kCoG:   pCluMapCoG->Fill(clux,cluy); break;
            case AliHMPIDCluster::kNoLoc: pCluMapNoLoc->Fill(clux,cluy); break;
            case AliHMPIDCluster::kAbn:   pCluMapAbn->Fill(clux,cluy); break;
            case AliHMPIDCluster::kNot:   pCluMapNot->Fill(clux,cluy); break;
            default:     pCluMapEmp->Fill(clux,cluy); break; 		
        }
        
        hHitCluDifX->Fill(hitx-clux); hHitCluDifY->Fill(hity-cluy); hHitCluDifXY->Fill(hitx-clux,hity-cluy); hHitCluDifQ->Fill(100*(cluq-hitq)/hitq);
        
      }//clusters loop
    }//chambers loop
      
    hits.Delete();  sdigs.Delete();  for(int i = 0;i<7;i++){((TClonesArray*)digs.At(i))->Delete();((TClonesArray*)clus.At(i))->Delete();}
  }//events loop      
      
  gStyle->SetPalette(1);
  TCanvas *pC2=new TCanvas("Digit canvas","Digit canvas",1280,800); pC2->Divide(3,3);
  pC2->cd(1);hHitCluDifX->Draw("hist");
  pC2->cd(2);hHitCluDifY->Draw("hist");
  pC2->cd(3);hHitCluDifXY->Draw("colz");
  pC2->cd(4);hHitCluDifQ->Draw("hist");
  pC2->cd(5);gPad->SetLogy(1);hCluFlg->Draw();
  pC2->cd(6);hCluChi2->Draw();
  pC2->cd(7);                 hCluRawSize->Draw();
  pC2->cd(8);                 hCluPerEvt->Draw();
  pC2->cd(9);                 hQ7eV->Draw(); hQ200eV->Draw("same"); hQ500eV->Draw("same"); hQ900eV->Draw("same"); lQ->Draw();
  TCanvas *pC1=new TCanvas("ClusterMaps","Cluster maps",1280,800); pC1->Divide(3,3);
  pC1->cd(1);  pCluMapSi1->Draw();  DrawCh(-1);
  pC1->cd(2);  pCluMapLo1->Draw();  DrawCh(-1);
  pC1->cd(3);  pCluMapUnf->Draw();  DrawCh(-1);
  pC1->cd(4);     hHitMap->Draw();  DrawCh(-1);
  pC1->cd(5);  pCluMapMax->Draw();  DrawCh(-1);
  pC1->cd(6);  pCluMapEdg->Draw();  DrawCh(-1);
  pC1->cd(7);  pCluMapNot->Draw();  DrawCh(-1);
  pC1->cd(8);  pCluMapNoLoc->Draw();DrawCh(-1);
  pC1->cd(9);  pCluMapEmp->Draw();  DrawCh(-1);
   
  pC1->SaveAs("$HOME/HitMaps.png");  //?????
  pC2->SaveAs("$HOME/HitCluDif.gif");  
  
  Printf("Digits - raw -digits conversion...");  
  
  AliHMPIDDigit d1,d2; Int_t ddl,r,d,a;UInt_t w32;

  for(Int_t iCh=AliHMPIDDigit::kMinCh;iCh<=AliHMPIDDigit::kMaxCh;iCh++)
  for(Int_t iPc=AliHMPIDDigit::kMinPc;iPc<=AliHMPIDDigit::kMaxPc;iPc++)
  for(Int_t iPx=AliHMPIDDigit::kMinPx;iPx<=AliHMPIDDigit::kMaxPx;iPx++)
  for(Int_t iPy=AliHMPIDDigit::kMinPy;iPy<=AliHMPIDDigit::kMaxPy;iPy++){
    d1.Set(iCh,iPc,iPx,iPy,3040);   //set digit               
    d1.Raw(w32,ddl,r,d,a);          //get raw word for this digit 
    d2.Raw(w32,ddl);                //set another digit from that raw word
    if(d1.Compare(&d2)) {d1.Print(); d2.Print(); Printf("");}//compare them
  }
  Printf("OK");
}//tst()
