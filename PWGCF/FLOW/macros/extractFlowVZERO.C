TFile *fo = NULL;
TF1 *fPsi = NULL;

Bool_t kLib = kFALSE;
Bool_t kVZEROrescorr = kTRUE; // VZERO resolution corrections
Bool_t kNUAcorr = kTRUE; // NUA corrections
Bool_t kNUOcorr = kFALSE; // Not Uniform Occupancy (TOF) corrections
Bool_t kPIDcorr = kFALSE; // PID corrections
Bool_t kCleanMemory = kTRUE; // if kFALSE take care of the large numbers of TCanvases

TH1D *hNUO[8]; // NUO corrections

void LoadLib();
void extractFlowVZERO(Int_t icentr,Int_t arm=2,Float_t pTh=0.8,Bool_t isMC=kFALSE,Int_t addbin=0);
TProfile *extractFlowVZEROsingle(Int_t icentr,Int_t spec,Int_t arm=2,Bool_t isMC=kFALSE,Float_t pTh=0.9,Int_t addbin=0,const char *nameSp="",Float_t detMin=0,Float_t detMax=1,Int_t chMin=-1,Int_t chMax=1);
void compareTPCTOF(Int_t icentr,Int_t spec,Int_t arm=2,Float_t pTh=0.9,Int_t addbin=0);
void plotQApid(Int_t ic,Float_t pt,Int_t addbin=0);

void LoadLib(){
  if(! kLib){
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libAOD.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libCORRFW.so");
    gSystem->Load("libNetx.so");
    gSystem->Load("libPWGflowBase.so");
    gSystem->Load("libPWGflowTasks.so");
    
    kLib = kTRUE;
  }
}

void extractFlowVZERO(Int_t icentr,Int_t arm,Float_t pTh,Bool_t isMC,Int_t addbin){
  LoadLib();

  new TCanvas();

  Int_t cMin[] = {0,5,10,20,30,40,50,60,70};
  Int_t cMax[] = {5,10,20,30,40,50,60,70,80};

  if(kNUOcorr){ // Compute correction for NUO in TOF
    compareTPCTOF(icentr,0,arm,pTh,addbin);
//     compareTPCTOF(icentr,1,arm,pTh,addbin);
//     compareTPCTOF(icentr,2,arm,pTh,addbin);
//     compareTPCTOF(icentr,3,arm,pTh,addbin);
//     compareTPCTOF(icentr,4,arm,pTh,addbin);
  }

  TProfile *pAll;
  pAll=extractFlowVZEROsingle(icentr,0,arm,0,pTh,addbin,"all",0,1);
  pAll->SetMarkerStyle(24);
  TProfile *pPiTOF,*pPiTPC,*pPiTPC2;
  pPiTOF=extractFlowVZEROsingle(icentr,1,arm,0,pTh,addbin,"piTOF",1,1);
  pPiTPC=extractFlowVZEROsingle(icentr,1,arm,0,pTh,addbin,"piTPC",0,0);
  pPiTPC2=extractFlowVZEROsingle(icentr,1,arm,0,pTh,addbin,"piTPC2",2,2);
  pPiTPC->Add(pPiTPC2);

  TH1D *hPi = pPiTOF->ProjectionX("hPi");
  hPi->SetLineColor(4);
  hPi->SetMarkerColor(4);
  hPi->SetMarkerStyle(20);
  for(Int_t i=1;i <=hPi->GetNbinsX();i++){
    Float_t x = hPi->GetBinCenter(i);
    if(x < 0.2){
      hPi->SetBinContent(i,0);
      hPi->SetBinError(i,0);
    }
    else if(x < 0.5){
      hPi->SetBinContent(i,pPiTPC->GetBinContent(i));
      hPi->SetBinError(i,pPiTPC->GetBinError(i));      
    }
    else{
      if(kNUOcorr){
	hPi->SetBinContent(i,pPiTOF->GetBinContent(i) + hNUO[0]->GetBinContent(i));
	hPi->SetBinError(i,pPiTOF->GetBinError(i));      
      }
      else{
	hPi->SetBinContent(i,pPiTOF->GetBinContent(i));
	hPi->SetBinError(i,pPiTOF->GetBinError(i));      
      }
    }
  }

  TProfile *pElTOF,*pElTPC,*pElTPC2;
  pElTOF=extractFlowVZEROsingle(icentr,4,arm,0,pTh,addbin,"piTOF",1,1);
  pElTPC=extractFlowVZEROsingle(icentr,4,arm,0,pTh,addbin,"piTPC",0,0);
  pElTPC2=extractFlowVZEROsingle(icentr,4,arm,0,pTh,addbin,"piTPC2",2,2);
  pElTPC->Add(pElTPC2);

  TH1D *hEl = pElTOF->ProjectionX("hEl");
  hEl->SetLineColor(6);
  hEl->SetMarkerColor(6);
  hEl->SetMarkerStyle(20);
  for(Int_t i=1;i <=hEl->GetNbinsX();i++){
    Float_t x = hEl->GetBinCenter(i);
    if(x < 0.2){
      hEl->SetBinContent(i,0);
      hEl->SetBinError(i,0);
    }
    else if(x < 0.3){
      hEl->SetBinContent(i,pElTPC->GetBinContent(i));
      hEl->SetBinError(i,pElTPC->GetBinError(i));      
    }
    else{
      if(kNUOcorr){
	hEl->SetBinContent(i,pElTOF->GetBinContent(i) + hNUO[0]->GetBinContent(i));
	hEl->SetBinError(i,pElTOF->GetBinError(i));      
      }
      else{
	hEl->SetBinContent(i,pElTOF->GetBinContent(i));
	hEl->SetBinError(i,pElTOF->GetBinError(i));      
      }
    }
  }

  TProfile *pKTOF,*pKTPC,*pKTPC2;
  pKTOF=extractFlowVZEROsingle(icentr,2,arm,0,pTh,addbin,"kaTOF",1,1);
  pKTPC=extractFlowVZEROsingle(icentr,2,arm,0,pTh,addbin,"kaTPC",0,0);
  pKTPC2=extractFlowVZEROsingle(icentr,2,arm,0,pTh,addbin,"kaTPC2",2,2);
  pKTPC->Add(pKTPC2);

  TH1D *hK = pKTOF->ProjectionX("hKa");
  hK->SetLineColor(1);
  hK->SetMarkerColor(1);
  hK->SetMarkerStyle(21);
  for(Int_t i=1;i <=hK->GetNbinsX();i++){
    Float_t x = hK->GetBinCenter(i);
    if(x < 0.25){
      hK->SetBinContent(i,0);
      hK->SetBinError(i,0);
    }
    else if(x < 0.45){
      hK->SetBinContent(i,pKTPC->GetBinContent(i));
      hK->SetBinError(i,pKTPC->GetBinError(i));      
    }
    else{
      if(kNUOcorr){
	hK->SetBinContent(i,pKTOF->GetBinContent(i) + hNUO[0]->GetBinContent(i));
	hK->SetBinError(i,pKTOF->GetBinError(i));      
      }
      else{
	hK->SetBinContent(i,pKTOF->GetBinContent(i));
	hK->SetBinError(i,pKTOF->GetBinError(i));      
      }
    }
  }

  TProfile *pPrTOF,*pPrTOF2,*pPrTPC,*pPrTPC2;
  pPrTOF=extractFlowVZEROsingle(icentr,3,arm,0,pTh,addbin,"prTOF",1,1,-1,-1);
  pPrTOF2=extractFlowVZEROsingle(icentr,3,arm,0,pTh,addbin,"prTOF2",1,1,-1,1);
  pPrTPC=extractFlowVZEROsingle(icentr,3,arm,0,pTh,addbin,"prTPC",0,0,-1,-1);
  pPrTPC2=extractFlowVZEROsingle(icentr,3,arm,0,pTh,addbin,"prTPC2",2,2,-1,-1);
  pPrTPC->Add(pPrTPC2);

  TH1D *hPr = pPrTOF->ProjectionX("hPr");
  hPr->SetLineColor(2);
  hPr->SetMarkerColor(2);
  hPr->SetMarkerStyle(22);
  for(Int_t i=1;i <=hPr->GetNbinsX();i++){
    Float_t x = hPr->GetBinCenter(i);
    if(x < 0.3){
      hPr->SetBinContent(i,0);
      hPr->SetBinError(i,0);
    }
    else if(x < 1.0){
      hPr->SetBinContent(i,pPrTPC->GetBinContent(i));
      hPr->SetBinError(i,pPrTPC->GetBinError(i));      
    }
    else{
      if(x < 3){
	if(kNUOcorr){
	  hPr->SetBinContent(i,pPrTOF->GetBinContent(i) + hNUO[0]->GetBinContent(i));
	  hPr->SetBinError(i,pPrTOF->GetBinError(i));      
	}
	else{
	  hPr->SetBinContent(i,pPrTOF->GetBinContent(i));
	  hPr->SetBinError(i,pPrTOF->GetBinError(i));      
	}
      }
      else{
	if(kNUOcorr){
	  hPr->SetBinContent(i,pPrTOF2->GetBinContent(i) + hNUO[0]->GetBinContent(i));
	  hPr->SetBinError(i,pPrTOF2->GetBinError(i));      
	}
	else{
	  hPr->SetBinContent(i,pPrTOF2->GetBinContent(i));
	  hPr->SetBinError(i,pPrTOF2->GetBinError(i));      
	}
      }
    }
  }
  
  pAll->Draw();
  hPi->Draw("SAME");
  hK->Draw("SAME");
  hPr->Draw("SAME");


  char name[100];

  // PID correction
  if(kPIDcorr){
    TFile *fPidTOF = new TFile("$ALICE_ROOT/PWGCF/FLOW/other/BayesianPIDcontTPCTOF.root");
    TFile *fPidTPC = new TFile("$ALICE_ROOT/PWGCF/FLOW/other/BayesianPIDcontTPC.root");
    // pi histos
    sprintf(name,"Pi_IDas_Picentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidPiPi=(TH1D *) fPidTOF->Get(name);
    sprintf(name,"Pi_IDas_Elcentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidPiEl=(TH1D *) fPidTOF->Get(name);
    sprintf(name,"Pi_IDas_Kacentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidPiKa=(TH1D *) fPidTOF->Get(name);
    sprintf(name,"Pi_IDas_Prcentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidPiPr=(TH1D *) fPidTOF->Get(name);
    TH1D *hPidAll = new TH1D(*hPidPiPi);
    hPidAll->Add(hPidPiKa);
    hPidAll->Add(hPidPiPr);
    hPidAll->Add(hPidPiEl);
    hPidPiPi->Divide(hPidAll);
    hPidPiKa->Divide(hPidAll);
    hPidPiPr->Divide(hPidAll);
    hPidPiEl->Divide(hPidAll);

    sprintf(name,"Pi_IDas_Picentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidPiPiTPC=(TH1D *) fPidTPC->Get(name);
    sprintf(name,"Pi_IDas_Elcentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidPiElTPC=(TH1D *) fPidTPC->Get(name);
    sprintf(name,"Pi_IDas_Kacentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidPiKaTPC=(TH1D *) fPidTPC->Get(name);
    sprintf(name,"Pi_IDas_Prcentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidPiPrTPC=(TH1D *) fPidTPC->Get(name);
    hPidAll->Reset();
    hPidAll->Add(hPidPiPiTPC);
    hPidAll->Add(hPidPiKaTPC);
    hPidAll->Add(hPidPiPrTPC);
    hPidAll->Add(hPidPiElTPC);
    hPidPiPiTPC->Divide(hPidAll);
    hPidPiKaTPC->Divide(hPidAll);
    hPidPiPrTPC->Divide(hPidAll);
    hPidPiElTPC->Divide(hPidAll);

    // K histos
    sprintf(name,"Ka_IDas_Picentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidKaPi=(TH1D *) fPidTOF->Get(name);
    sprintf(name,"Ka_IDas_Elcentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidKaEl=(TH1D *) fPidTOF->Get(name);
    sprintf(name,"Ka_IDas_Kacentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidKaKa=(TH1D *) fPidTOF->Get(name);
    sprintf(name,"Ka_IDas_Prcentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidKaPr=(TH1D *) fPidTOF->Get(name);
    hPidAll->Reset();
    hPidAll->Add(hPidKaPi);
    hPidAll->Add(hPidKaKa);
    hPidAll->Add(hPidKaPr);
    hPidAll->Add(hPidKaEl);
    hPidKaPi->Divide(hPidAll);
    hPidKaKa->Divide(hPidAll);
    hPidKaPr->Divide(hPidAll);
    hPidKaEl->Divide(hPidAll);

    sprintf(name,"Ka_IDas_Picentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidKaPiTPC=(TH1D *) fPidTPC->Get(name);
    sprintf(name,"Ka_IDas_Elcentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidKaElTPC=(TH1D *) fPidTPC->Get(name);
    sprintf(name,"Ka_IDas_Kacentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidKaKaTPC=(TH1D *) fPidTPC->Get(name);
    sprintf(name,"Ka_IDas_Prcentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidKaPrTPC=(TH1D *) fPidTPC->Get(name);
    hPidAll->Reset();
    hPidAll->Add(hPidKaPiTPC);
    hPidAll->Add(hPidKaKaTPC);
    hPidAll->Add(hPidKaPrTPC);
    hPidAll->Add(hPidKaElTPC);
    hPidKaPiTPC->Divide(hPidAll);
    hPidKaKaTPC->Divide(hPidAll);
    hPidKaPrTPC->Divide(hPidAll);
    hPidKaElTPC->Divide(hPidAll);

    // pr histos
    sprintf(name,"Pr_IDas_Picentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidPrPi=(TH1D *) fPidTOF->Get(name);
    sprintf(name,"Pr_IDas_Elcentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidPrEl=(TH1D *) fPidTOF->Get(name);
    sprintf(name,"Pr_IDas_Kacentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidPrKa=(TH1D *) fPidTOF->Get(name);
    sprintf(name,"Pr_IDas_Prcentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidPrPr=(TH1D *) fPidTOF->Get(name);
    hPidAll->Reset();
    hPidAll->Add(hPidPrPi);
    hPidAll->Add(hPidPrKa);
    hPidAll->Add(hPidPrPr);
    hPidAll->Add(hPidPrEl);
    hPidPrPi->Divide(hPidAll);
    hPidPrKa->Divide(hPidAll);
    hPidPrPr->Divide(hPidAll);
    hPidPrEl->Divide(hPidAll);

    sprintf(name,"Pr_IDas_Picentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidPrPiTPC=(TH1D *) fPidTPC->Get(name);
    sprintf(name,"Pr_IDas_Elcentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidPrElTPC=(TH1D *) fPidTPC->Get(name);
    sprintf(name,"Pr_IDas_Kacentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidPrKaTPC=(TH1D *) fPidTPC->Get(name);
    sprintf(name,"Pr_IDas_Prcentr%i_pth%4.2f",icentr,pTh);
    TH1D *hPidPrPrTPC=(TH1D *) fPidTPC->Get(name);
    hPidAll->Reset();
    hPidAll->Add(hPidPrPiTPC);
    hPidAll->Add(hPidPrKaTPC);
    hPidAll->Add(hPidPrPrTPC);
    hPidAll->Add(hPidPrElTPC);
    hPidPrPiTPC->Divide(hPidAll);
    hPidPrKaTPC->Divide(hPidAll);
    hPidPrPrTPC->Divide(hPidAll);
    hPidPrElTPC->Divide(hPidAll);
 
    for(Int_t k=1;k <=hPi->GetNbinsX();k++){
      Float_t pt = hPi->GetBinCenter(k);
      Float_t xPi=hPi->GetBinContent(k)*hPidPiPi->Interpolate(pt) + hK->GetBinContent(k)*hPidPiKa->Interpolate(pt) + hPr->GetBinContent(k)*hPidPiPr->Interpolate(pt) + hEl->GetBinContent(k)*hPidPiEl->Interpolate(pt);
      if(pt < 0.5)
	xPi=hPi->GetBinContent(k)*hPidPiPiTPC->Interpolate(pt) + hK->GetBinContent(k)*hPidPiKaTPC->Interpolate(pt) + hPr->GetBinContent(k)*hPidPiPrTPC->Interpolate(pt) + hEl->GetBinContent(k)*hPidPiElTPC->Interpolate(pt);
      Float_t xKa=hPi->GetBinContent(k)*hPidKaPi->Interpolate(pt) + hK->GetBinContent(k)*hPidKaKa->Interpolate(pt) + hPr->GetBinContent(k)*hPidKaPr->Interpolate(pt) + hEl->GetBinContent(k)*hPidKaEl->Interpolate(pt);
      if(pt < 0.45)
	xKa=hPi->GetBinContent(k)*hPidKaPiTPC->Interpolate(pt) + hK->GetBinContent(k)*hPidKaKaTPC->Interpolate(pt) + hPr->GetBinContent(k)*hPidKaPrTPC->Interpolate(pt) + hEl->GetBinContent(k)*hPidKaElTPC->Interpolate(pt);
      Float_t xPr=hPi->GetBinContent(k)*hPidPrPi->Interpolate(pt) + hK->GetBinContent(k)*hPidPrKa->Interpolate(pt) + hPr->GetBinContent(k)*hPidPrPr->Interpolate(pt) + hEl->GetBinContent(k)*hPidPrEl->Interpolate(pt);
      if(pt < 1)
	xPr=hPi->GetBinContent(k)*hPidPrPiTPC->Interpolate(pt) + hK->GetBinContent(k)*hPidPrKaTPC->Interpolate(pt) + hPr->GetBinContent(k)*hPidPrPrTPC->Interpolate(pt) + hEl->GetBinContent(k)*hPidPrElTPC->Interpolate(pt);

      hPi->SetBinContent(k,hPi->GetBinContent(k)*2 - xPi);
      hK->SetBinContent(k,hK->GetBinContent(k)*2 - xKa);
      hPr->SetBinContent(k,hPr->GetBinContent(k)*2 - xPr);    
    }
  }

  // write output
  snprintf(name,100,"results%03i-%03iv%i_pTh%3.1f.root",cMin[icentr],cMax[icentr+addbin],arm,pTh);
  TFile *fout = new TFile(name,"RECREATE");
  pAll->ProjectionX()->Write();
  hPi->Write();
  hK->Write();
  hPr->Write();
  if(isMC){
    TH1D *pTmp = extractFlowVZEROsingle(icentr,0,arm,1,pTh,addbin,"allMC",1,1,-1,1)->ProjectionX();
    pTmp->SetLineColor(6);
    pTmp->SetMarkerColor(6);
    pTmp->SetMarkerStyle(24);
    pTmp->Write();
    pTmp = extractFlowVZEROsingle(icentr,1,arm,1,pTh,addbin,"piMC",1,1,-1,1)->ProjectionX();
    pTmp->SetLineColor(4);
    pTmp->SetMarkerColor(4);
    pTmp->SetMarkerStyle(24);
    pTmp->Write();
    pTmp = extractFlowVZEROsingle(icentr,2,arm,1,pTh,addbin,"kaMC",1,1,-1,1)->ProjectionX();
    pTmp->SetLineColor(1);
    pTmp->SetMarkerColor(1);
    pTmp->SetMarkerStyle(25);
    pTmp->Write();
    pTmp = extractFlowVZEROsingle(icentr,3,arm,1,pTh,addbin,"prMC",1,1,-1,-1)->ProjectionX();
    pTmp->SetLineColor(2);
    pTmp->SetMarkerColor(2);
    pTmp->SetMarkerStyle(26);
    pTmp->Write();
  }
  fout->Close();
}

TProfile *extractFlowVZEROsingle(Int_t icentr,Int_t spec,Int_t arm,Bool_t isMC,Float_t pTh,Int_t addbin,const char *nameSp,Float_t detMin,Float_t detMax,Int_t chMin,Int_t chMax){
  LoadLib();

  pTh += 0.00001;

  // NUA correction currently are missing
  char name[100];
  char stringa[200];

  snprintf(name,100,"AnalysisResults.root");
  if(!fo) fo = new TFile(name);
  snprintf(name,100,"contVZEROv%i",arm);
  TList *cont = (TList *) fo->Get(name);

  cont->ls();

  Float_t xMin[5] = {icentr/*centrality bin*/,chMin/*charge*/,pTh/*prob*/,-TMath::Pi()/arm/*Psi*/,detMin/*PID mask*/};
  Float_t xMax[5] = {icentr+addbin,chMax,1,TMath::Pi()/arm,detMax};

  cont->ls();

  TProfile *p1 = cont->At(2);
  TProfile *p2 = cont->At(3);
  TProfile *p3 = cont->At(4);

  TH2F *hPsi2DA = cont->At(5);
  TH2F *hPsi2DC = cont->At(6);

  TH1D *hPsiA = hPsi2DA->ProjectionY("PsiA",icentr+1,icentr+addbin+1);
  TH1D *hPsiC = hPsi2DC->ProjectionY("PsiC",icentr+1,icentr+addbin+1);

  if(!fPsi) fPsi = new TF1("fPsi","pol0",-TMath::Pi()/arm,TMath::Pi()/arm);
  hPsiA->Fit(fPsi,"0");
  Float_t offsetA = fPsi->GetParameter(0);
  hPsiC->Fit(fPsi,"0");
  Float_t offsetC = fPsi->GetParameter(0);

  Int_t nbinPsi = hPsiA->GetNbinsX();

  Float_t *NUAcorrA = new Float_t[nbinPsi];
  Float_t *NUAcorrC = new Float_t[nbinPsi];

  for(Int_t i=0;i < nbinPsi;i++){
    NUAcorrA[i] = offsetA/(hPsiA->GetBinContent(i+1));
    NUAcorrC[i] = offsetC/(hPsiC->GetBinContent(i+1));
  }

  Float_t res1=0,res2=0,res3=0; 
  Float_t eres1=0,eres2=0,eres3=0; 

  for(Int_t i = icentr; i <= icentr+addbin;i++){
    if(p1->GetBinError(i+1)){
      eres1 += 1./p1->GetBinError(i+1)/p1->GetBinError(i+1);
      res1 += p1->GetBinContent(i+1)/p1->GetBinError(i+1)/p1->GetBinError(i+1);      
    }
    if(p2->GetBinError(i+1)){
      eres2 += 1./p2->GetBinError(i+1)/p2->GetBinError(i+1);
      res2 += p2->GetBinContent(i+1)/p2->GetBinError(i+1)/p2->GetBinError(i+1);      
    }
    if(p3->GetBinError(i+1)){
      eres3 += 1./p3->GetBinError(i+1)/p3->GetBinError(i+1);
      res3 += p3->GetBinContent(i+1)/p3->GetBinError(i+1)/p3->GetBinError(i+1);      
    }
  }

  res1 /= eres1;
  res2 /= eres2;
  res3 /= eres3;
  
  eres1 = sqrt(1./eres1);
  eres2 = sqrt(1./eres2);
  eres3 = sqrt(1./eres3);

  AliFlowVZEROResults *a = (AliFlowVZEROResults *) cont->At(0);
  AliFlowVZEROResults *b = (AliFlowVZEROResults *) cont->At(1);
  TProfile *pp,*pp2;
  if(kNUAcorr){ // with NUA corrections
    pp = a->GetV2reweight(spec,xMin,xMax,3,NUAcorrA);
    pp2 = b->GetV2reweight(spec,xMin,xMax,3,NUAcorrC);
  }
  else{
    pp = a->GetV2(spec,xMin,xMax);
    pp2 = b->GetV2(spec,xMin,xMax);
  }
  
  Float_t scaling = sqrt(res1*res3/res2);
  if(kVZEROrescorr){
    pp->Scale(1./scaling);
  }

  Float_t err1_2 = eres1*eres1/res1/res1/4 +
    eres2*eres2/res2/res2/4 +
    eres3*eres3/res3/res3/4;
  Float_t err2_2 = err1_2;
  err1_2 /= scaling*scaling;
  printf("resolution V0A = %f +/- %f\n",scaling,err1_2);
  scaling = sqrt(res2*res3/res1);
  err2_2 /= scaling*scaling;
  if(kVZEROrescorr){
    pp2->Scale(1./scaling);
  }
  printf("resolution V0C =%f +/- %f\n",scaling,err2_2);

  pp->SetName("V0A");
  pp2->SetName("V0C");

  if(! kCleanMemory){
    new TCanvas();  
    pp->Draw();
    pp2->Draw("SAME");
  }

  TProfile *pData = new TProfile(*pp);
  pData->Add(pp2);
  snprintf(stringa,100,"%sData",nameSp);
  pData->SetName(stringa);

  TProfile *pMc = NULL;
  
  TProfile *ppMC;
  TProfile *ppMC2;
  
  if(isMC){
    snprintf(name,100,"contVZEROmc");
    cont = (TList *) fo->Get(name);
    cont->ls();
    if(arm == 2){
      AliFlowVZEROResults *c = (AliFlowVZEROResults *) cont->At(0);
      if(! kCleanMemory) c->GetV2(spec,xMin,xMax)->Draw("SAME");
    }
    AliFlowVZEROResults *cA;
    if(fo->Get("contVZEROv2")) cA = (AliFlowVZEROResults *) cont->At(1+2*(arm==3));
    else cA = (AliFlowVZEROResults *) cont->At(0);
    AliFlowVZEROResults *cC;
    if(fo->Get("contVZEROv2")) cC = (AliFlowVZEROResults *) cont->At(2+2*(arm==3));
    else cC = (AliFlowVZEROResults *) cont->At(1);

    TProfile *p1mc = cont->At(5+3*(arm==3));
    TProfile *p2mc = cont->At(6+3*(arm==3));
    TProfile *p3mc = cont->At(7+3*(arm==3));

    Float_t resMC1=0,resMC2=0,resMC3=0; 
    Float_t eresMC1=0,eresMC2=0,eresMC3=0; 

    for(Int_t i = icentr; i <= icentr+addbin;i++){
      if(p1mc->GetBinError(i+1)){
	eresMC1 += 1./p1mc->GetBinError(i+1)/p1mc->GetBinError(i+1);
	resMC1 += p1mc->GetBinContent(i+1)/p1mc->GetBinError(i+1)/p1mc->GetBinError(i+1);      
      }
      if(p2mc->GetBinError(i+1)){
	eresMC2 += 1./p2mc->GetBinError(i+1)/p2mc->GetBinError(i+1);
	resMC2 += p2mc->GetBinContent(i+1)/p2mc->GetBinError(i+1)/p2mc->GetBinError(i+1);      
      }
      if(p3mc->GetBinError(i+1)){
	eresMC3 += 1./p3mc->GetBinError(i+1)/p3mc->GetBinError(i+1);
	resMC3 += p3mc->GetBinContent(i+1)/p3mc->GetBinError(i+1)/p3mc->GetBinError(i+1);      
      }
    }
    
    resMC1 /= eresMC1;
    resMC2 /= eresMC2;
    resMC3 /= eresMC3;
    
    eresMC1 = sqrt(1./eresMC1);
    eresMC2 = sqrt(1./eresMC2);
    eresMC3 = sqrt(1./eresMC3);

    ppMC = cA->GetV2(spec,xMin,xMax);
    ppMC2 = cC->GetV2(spec,xMin,xMax);
    
    scaling = sqrt(resMC1*resMC3/resMC2);
    ppMC->Scale(1./scaling);
  
    err1_2 = eresMC1*eresMC1/resMC1/resMC1/4 +
      eresMC2*eresMC2/resMC2/resMC2/4 +
      eresMC3*eresMC3/resMC3/resMC3/4;
    err2_2 = err1_2;
    err1_2 /= scaling*scaling;
    printf("resolution V0A (MC) = %f +/- %f\n",scaling,err1_2);
    scaling = sqrt(resMC2*resMC3/resMC1);
    err2_2 /= scaling*scaling;
    ppMC2->Scale(1./scaling);
    printf("resolution V0C (MC) =%f +/- %f\n",scaling,err2_2);

    ppMC->SetName("V0Amc");
    ppMC2->SetName("V0Cmc");

    if(! kCleanMemory){
      ppMC->Draw("SAME");
      ppMC2->Draw("SAME");
    }

    pMc = new TProfile(*ppMC);
    pMc->Add(ppMC2);
    snprintf(stringa,100,"%sMC",nameSp);
    pMc->SetName(stringa);
    pMc->SetLineColor(2);
  }

  if(! kCleanMemory){
    new TCanvas();  
    pData->Draw();
  }
  if(pMc && !kCleanMemory){
    pMc->Draw("SAME");
    TH1D *hData = pData->ProjectionX();
    TH1D *hMc = pMc->ProjectionX();
    hData->Divide(hMc);
    new TCanvas();  
    hData->Draw();
  }

  delete[] NUAcorrA;
  delete[] NUAcorrC;

  if(kCleanMemory){
    if(pp) delete pp;
    if(pp2) delete pp2;
    if(ppMC) delete ppMC;
    if(ppMC2) delete ppMC2;
  }

  delete cont;

  if(isMC) return pMc;
  return pData;
}
void compareTPCTOF(Int_t icentr,Int_t spec,Int_t arm,Float_t pTh,Int_t addbin){
  LoadLib();
  Float_t ptMaxFit[8] = {20,0.7,0.55,1.2,2,2,2,2};

  TProfile *pTPC = extractFlowVZEROsingle(icentr,spec,arm,0,pTh,addbin,"TPC",0,0,-1,-1+2*(spec!=3)); // TPC && !TOF (TPC PID)
  TProfile *pTPC2 = extractFlowVZEROsingle(icentr,spec,arm,0,pTh,addbin,"TPCextra",2,2,-1,-1+2*(spec!=3)); //TPC && TOF (TPC PID)

  TProfile *pTOF = extractFlowVZEROsingle(icentr,spec,arm,0,pTh,addbin,"TOF",1,1,-1,-1+2*(spec!=3)); //TPC && TOF (TPC&TOF PID)

  if(spec > 0) pTPC->Add(pTPC2);
  else pTPC->Add(pTOF);

  char name[100];
  snprintf(name,100,"NUO_%i",spec);

  hNUO[spec] = pTPC->ProjectionX(name);
  hNUO[spec]->Add(pTOF,-1);

  if(spec && hNUO[0]) hNUO[spec]->Add(hNUO[0],-1);

  if(! kCleanMemory){
    new TCanvas();
    pTPC->Draw();
    pTOF->Draw("SAME");

    new TCanvas();  
    pTPC->SetLineColor(1);
    pTOF->SetLineColor(2);
  }

  hNUO[spec]->Draw();

  if(kCleanMemory){
    if(pTPC) delete pTPC;
    if(pTPC2) delete pTPC2;
    if(pTOF) delete pTOF;
  }
}

void plotQApid(Int_t ic,Float_t pt,Int_t addbin){
  char name[100];
  char stringa[200];
  LoadLib();

  snprintf(name,100,"AnalysisResults.root");
  if(!fo) fo = new TFile(name);
  snprintf(name,100,"contVZEROv2");
  TList *cont = (TList *) fo->Get(name);
  AliFlowVZEROResults *pidqa = cont->FindObject("qaPID");
  Float_t xval[2] = {2.,pt+0.00001};
  Float_t xval2[2] = {2.+addbin,pt+0.00001};

  TProfile *proTPCpi = pidqa->GetV2(0,xval,xval2);
  TProfile *proTOFpi = pidqa->GetV2(1,xval,xval2);
  TProfile *proTPCka = pidqa->GetV2(2,xval,xval2);
  TProfile *proTOFka = pidqa->GetV2(3,xval,xval2);
  TProfile *proTPCpr = pidqa->GetV2(4,xval,xval2);
  TProfile *proTOFpr = pidqa->GetV2(5,xval,xval2);

  proTPCpi->SetName("hPiTPC");
  proTOFpi->SetName("hPiTOF");
  proTPCka->SetName("hKaTPC");
  proTOFka->SetName("hKaTOF");
  proTPCpr->SetName("hPrTPC");
  proTOFpr->SetName("hPrTOF");
  
  proTPCpi->GetXaxis()->SetTitle("dE/dx - dE/dx_{calc}^{#pi} (a.u)");
  proTOFpi->GetXaxis()->SetTitle("t_{tof} - t_{calc}^{#pi} (ps)");
  proTPCka->GetXaxis()->SetTitle("dE/dx - dE/dx_{calc}^{K} (a.u)");
  proTOFka->GetXaxis()->SetTitle("t_{tof} - t_{calc}^{K} (ps)");
  proTPCpr->GetXaxis()->SetTitle("dE/dx - dE/dx_{calc}^{p} (a.u)");
  proTOFpr->GetXaxis()->SetTitle("t_{tof} - t_{calc}^{p} (ps)");

  TCanvas *c = new TCanvas("cVproj","cVproj");
  c->Divide(2,3);
  c->cd(1);
  proTPCpi->Draw();
  c->cd(2);
  proTOFpi->Draw();
  c->cd(3);
  proTPCka->Draw();
  c->cd(4);
  proTOFka->Draw();
  c->cd(5);
  proTPCpr->Draw();
  c->cd(6);
  proTOFpr->Draw();
}
