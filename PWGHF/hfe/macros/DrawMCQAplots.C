void DrawMCQAplots(Int_t imode){

  // if you want to draw differential yield for charm, set imode==4
  // if you want to draw differential yield for charm, set imode==5
  // if you want to draw differential yield for charm and beauty, set imode==0
  // if you want to compare differential yield of electron yield from charm and beauty w/wo acc cut(eta<0.9), set imode==11
  // if you want to compare differetail cross section together with NLO prediction, set imode==12

  // load libs
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG0base.so");
  gSystem->Load("../libPWG3hfeDevel.so");

  // put input file name
  //char filename[]="./HFEtask.root";
  char filename[]="/lustre/alice/train/V006.MC_pp/2011-01-07_1954.4718/mergedPeriods/MC_pp/7TeV/LHC10f7a_c/HFEtask.root"; // pwg3 samples for pp
  //char NLOcalcc[]="../rootfiles/cpp10CTQ6diff.root"; // p+p 10TeV for cquark (pt, rapidity for q and fragments )
  //char NLOcalcb[]="../rootfiles/bpp10CTQ6diff.root"; // p+p 10TeV for bquark (pt, rapidity for q and fragments )
  //char NLOcalccb[]="../rootfiles/pp10CTQ6_pt_yall_cb.root"; // p+p 7TeV for c and b (pT only)
  char NLOcalcc[]="../rootfiles/cpp7CTQ6diff.root"; // p+p 7TeV for cquark (pt, rapidity for q and fragments )
  char NLOcalcb[]="../rootfiles/bpp7CTQ6diff.root"; // p+p 7TeV for bquark (pt, rapidity for q and fragments )
  char NLOcalccb[]="../rootfiles/pp7CTQ6_pt_yall_cb.root"; // p+p 7TeV for c and b (pT only)


  if(imode == 4) DrawDifferentialYield(4,filename, NLOcalcc); // charm differentail yield 
  if(imode == 5) DrawDifferentialYield(5,filename, NLOcalcb); // beauty differential yield
  if(imode == 0) DrawDifferentialYieldBoth(filename, NLOcalccb); // both charm and beauty
  if(imode == 11) DrawDifferentialYieldAcc(filename); // electron yield w/wo acc cut
  //if(imode == 12) DrawDifferentialXsectionBoth(filename, NLOcalccb); //be careful, have to look again the overall normalization  

}

//--------------------------
void DrawDifferentialYield(Int_t qtype, char *filename, char *NLOcalc){

  TFile *_file[2]; 
  _file[0] = TFile::Open(filename);
  _file[1] = TFile::Open(NLOcalc);

  setGeneralStyle();
  char configname[]="mcQA";
  enum qType {kQuark=0, kantiQuark=1, kHadron=2, kElec=3, kElec2nd=4, keHadron=5, kDeHadron=6};
  enum hqType {kCharmQuark=4, kBeautyQuark=5};
  Int_t fgkQuark = qtype;
  const Int_t fgkqType=7;

  struct hists {
        TH1F *Pt;
        TH1F *Eta;
        TH1F *Y;
  };

  struct histsComm {
        TH1F *PdgCode;
        TH1F *Nq;
  };

  hists fHist[7];
  histsComm fHistComm;

  TString kqTypeLabel[fgkqType];
  if (fgkQuark == kCharmQuark){
    kqTypeLabel[kQuark]="c";
    kqTypeLabel[kantiQuark]="cbar";
    kqTypeLabel[kHadron]="cHadron";
    kqTypeLabel[kElec]="ce";
    kqTypeLabel[kElec2nd]="nulle";
    kqTypeLabel[keHadron]="ceHadron";
    kqTypeLabel[kDeHadron]="nullHadron";
  } else if (fgkQuark == kBeautyQuark){
    kqTypeLabel[kQuark]="b";
    kqTypeLabel[kantiQuark]="bbar";
    kqTypeLabel[kHadron]="bHadron";
    kqTypeLabel[kElec]="be";
    kqTypeLabel[kElec2nd]="bce";
    kqTypeLabel[keHadron]="beHadron";
    kqTypeLabel[kDeHadron]="bDeHadron";
  }
  Int_t kColorCode[fgkqType];
  kColorCode[kQuark]=2;
  kColorCode[kantiQuark]=2;
  kColorCode[kHadron]=6;
  kColorCode[kElec]=4;
  kColorCode[kElec2nd]=4;
  kColorCode[keHadron]=3;
  kColorCode[kDeHadron]=3;

  Int_t kLineStyle[fgkqType];
  kLineStyle[kQuark]=1;
  kLineStyle[kantiQuark]=1;
  kLineStyle[kHadron]=3;
  kLineStyle[kElec]=2;
  kLineStyle[kElec2nd]=1;
  kLineStyle[keHadron]=2;
  kLineStyle[kDeHadron]=1;


  cPt = new TCanvas("cPt","pT",0,0,600,500);	
  cY = new TCanvas("cY","rapidity",0,0,600,500);	
  cPdg = new TCanvas("cPdg","pdg code",0,0,600,500);	
  cNq = new TCanvas("cNq","number of quark",0,0,600,500);	


  TList *tl = (TList *)_file[0]->Get("HFE_QA")->FindObject("MCqa");
  //count # of events
  TList *tl_result = (TList *)_file[0]->Get("HFE_Results");
  AliHFEcontainer *containerhfe = (AliHFEcontainer *) tl_result->FindObject("trackContainer");
  if(!containerhfe) {
      printf("No hfe container \n");
      return;
  }
  Int_t nEvt = (Int_t) containerhfe->GetNumberOfEvents();
  cout << "# of events " << nEvt << endl;

  Double_t scalefactor = 1./float(nEvt); //normalize with total number of event 

  TString hname; 
  for (Int_t iqType = 0; iqType < fgkqType; iqType++ ){

     //pT distribution
     hname="mcqa_Pt_";
     hname=hname+kqTypeLabel[iqType];
     fHist[iqType].Pt = (TH1F*)tl->FindObject(hname);
     fHist[iqType].Pt->SetXTitle("p_{t} [GeV/c]");
     fHist[iqType].Pt->SetYTitle("1/N_{Evt}dN/dp_{t} [1/GeVc^{-1}]");
     setDataStyle(*fHist[iqType].Pt, kColorCode[iqType], 3, kLineStyle[iqType]);

     CorrectFromTheWidth(fHist[iqType].Pt); //consider pT bin size
     fHist[iqType].Pt->Scale(scalefactor); //normalize with # of events 
     if (iqType==1) {
       fHist[iqType].Pt->Add(fHist[0].Pt); // add Q and Qbar
       fHist[iqType].Pt->Scale(0.5); // get number of ccbar pair
     }  

     //eta and rapidity distribution 	
     hname="mcqa_Eta_";
     hname=hname+kqTypeLabel[iqType];
     fHist[iqType].Eta = (TH1F*)tl->FindObject(hname);
     hname="mcqa_Y_";
     hname=hname+kqTypeLabel[iqType];
     fHist[iqType].Y = (TH1F*)tl->FindObject(hname);
     fHist[iqType].Y->SetXTitle("Y");
     fHist[iqType].Y->SetYTitle("1/N_{Evt}dN/dY");
     setDataStyle(*fHist[iqType].Y, kColorCode[iqType], 3, kLineStyle[iqType]);

     CorrectFromTheWidth(fHist[iqType].Y); //consider pT bin size
     fHist[iqType].Y->Scale(scalefactor); //normalize with # of events 
     if (iqType==1) {
       fHist[iqType].Y->Add(fHist[0].Y); // add Q and Qbar
       fHist[iqType].Y->Scale(0.5); // get number of ccbar pair
     }  

     cPt->cd();
     setPadStyle(2,gPad);
     if (iqType==1) fHist[iqType].Pt->Draw();
     if(iqType==2 || iqType==3 || iqType==4 || iqType==5 || iqType==6) fHist[iqType].Pt->Draw("same");

     cY->cd();
     setPadStyle(2,gPad);
     if (iqType==1) fHist[iqType].Y->Draw();
     if(iqType==2 || iqType==3 || iqType==4 || iqType==5 || iqType==6) fHist[iqType].Y->Draw("same");
     
  }

  cPdg->cd();	
  hname="mcqa_PdgCode_";
  hname=hname+kqTypeLabel[kQuark]+"Hadron";
  fHistComm.PdgCode = (TH1F*)tl->FindObject(hname);
  fHistComm.PdgCode->SetXTitle("hadron pdg code");
  fHistComm.PdgCode->SetYTitle("Yield");
  setDataStyle(*fHistComm.PdgCode, 2, 2, 1);
  setPadStyle(2,gPad);
  fHistComm.PdgCode->Draw();

  cNq->cd();
  hname="mcqa_Nq_";
  hname=hname+kqTypeLabel[kQuark];
  fHistComm.Nq = (TH1F*)tl->FindObject(hname);
  fHistComm.Nq->SetXTitle("number of "+kqTypeLabel[kQuark]+","+kqTypeLabel[kQuark]+"-bar per event");
  fHistComm.Nq->SetYTitle("Yield per Event");
  setDataStyle(*fHistComm.Nq, 2, 4, 1);
  setPadStyle(1,gPad);
  fHistComm.Nq->Scale(scalefactor); 
  fHistComm.Nq->Draw();



  TLegend *legend1 = new TLegend(0.50,0.73,0.98,0.99,"");
  setLegendStyle(*legend1,1);
  legend1->AddEntry(fHist[0].Pt,kqTypeLabel[kQuark]+"#bar{"+kqTypeLabel[kQuark]+"}"+" pair", "l");
  legend1->AddEntry(fHist[2].Pt,kqTypeLabel[kQuark]+" hadrons from "+kqTypeLabel[kQuark]+", #bar{"+kqTypeLabel[kQuark]+"}", "l");
  legend1->AddEntry(fHist[3].Pt,"elec. from "+kqTypeLabel[kQuark]+"->e", "l");
  if(qtype==5) legend1->AddEntry(fHist[4].Pt,"elec. from "+kqTypeLabel[kQuark]+"->c->e", "l");
  legend1->AddEntry(fHist[5].Pt,kqTypeLabel[kQuark]+" hadron decaying "+kqTypeLabel[kQuark]+"->e", "l");
  if(qtype==5) legend1->AddEntry(fHist[6].Pt,kqTypeLabel[kQuark]+" hadron decaying "+kqTypeLabel[kQuark]+"->c->e", "l");

  TLegend *legend1_ = new TLegend(0.65,0.42,0.90,0.72,"");
  setLegendStyle(*legend1_,1);
  legend1_->AddEntry(fHist[0].Pt,"Yield per event", "");
  legend1_->AddEntry(fHist[0].Pt,"(0.1<p_{t}<20)", "");
  legend1_->AddEntry(fHist[0].Pt,Form("%1.4f",fHist[1].Pt->Integral("width")), "l");
  legend1_->AddEntry(fHist[2].Pt,Form("%1.4f",fHist[2].Pt->Integral("width")), "l");
  legend1_->AddEntry(fHist[3].Pt,Form("%1.4f",fHist[3].Pt->Integral("width")), "l");
  if(qtype==5) legend1_->AddEntry(fHist[4].Pt,Form("%1.4f",fHist[4].Pt->Integral("width")), "l");
  legend1_->AddEntry(fHist[5].Pt,Form("%1.4f",fHist[5].Pt->Integral("width")), "l");
  if(qtype==5) legend1_->AddEntry(fHist[6].Pt,Form("%1.4f",fHist[6].Pt->Integral("width")), "l");

  cPt->cd();
  legend1_->Draw();

  TLegend *legend2 = new TLegend(0.26,0.15,0.74,0.37,"");
  setLegendStyle(*legend2,1);
  legend2->AddEntry(fHist[0].Y,kqTypeLabel[kQuark]+"#bar{"+kqTypeLabel[kQuark]+"}"+" pair", "l");
  legend2->AddEntry(fHist[2].Y,kqTypeLabel[kQuark]+" hadrons from "+kqTypeLabel[kQuark], "l");
  legend2->AddEntry(fHist[3].Y,"elec. from "+kqTypeLabel[kQuark]+"->e", "l");
  if(qtype==5) legend2->AddEntry(fHist[4].Y,"elec. from "+kqTypeLabel[kQuark]+"->c->e", "l");
  legend2->AddEntry(fHist[5].Y,kqTypeLabel[kQuark]+" hadron decaying "+kqTypeLabel[kQuark]+"->e", "l");
  if(qtype==5) legend2->AddEntry(fHist[6].Y,kqTypeLabel[kQuark]+" hadron decaying "+kqTypeLabel[kQuark]+"->c->e", "l");

  TLegend *legend2_ = new TLegend(0.70,0.77,0.98,0.99,"");
  setLegendStyle(*legend2_,1);
  legend2_->AddEntry(fHist[0].Y,"Yield per event", "");
  legend2_->AddEntry(fHist[0].Y,"(all p_{t},-7.5<Y<7.5)", "");
  legend2_->AddEntry(fHist[0].Y,Form("%1.4f",fHist[1].Y->Integral("width")), "l");
  legend2_->AddEntry(fHist[2].Y,Form("%1.4f",fHist[2].Y->Integral("width")), "l");
  legend2_->AddEntry(fHist[3].Y,Form("%1.4f",fHist[3].Y->Integral("width")), "l");
  if(qtype==5) legend2_->AddEntry(fHist[4].Y,Form("%1.4f",fHist[4].Y->Integral("width")), "l");
  legend2_->AddEntry(fHist[5].Y,Form("%1.4f",fHist[5].Y->Integral("width")), "l");
  if(qtype==5) legend2_->AddEntry(fHist[6].Y,Form("%1.4f",fHist[6].Y->Integral("width")), "l");

  cY->cd();
  legend2_->Draw();

  TLegend *legend3 = new TLegend(0.35,0.80,0.88,0.88,"");
  setLegendStyle(*legend3,0);
  legend3->AddEntry(fHistComm.Nq,"semi-leptonic decay", "");
	      
  cNq->cd();
  legend3->Draw();
  
  Double_t totcrossNLO=0;
  Double_t scalefactorNLO=0;
  Double_t totcrossNLOY=0;
  Double_t scalefactorNLOY=0;
  if(qtype==4) {
    //TH1F *hNLO = (TH1F*)_file[1]->Get("hpt_c");
    TH1F *hNLO = (TH1F*)_file[1]->Get("hptQwkick");
    TH1F *hNLOY = (TH1F*)_file[1]->Get("hyQwkick");
  }
  if(qtype==5) {
    //if(qtype==5) TH1F *hNLO = (TH1F*)_file[1]->Get("hpt_b");
    TH1F *hNLO = (TH1F*)_file[1]->Get("hptQwkick");
    TH1F *hNLOY = (TH1F*)_file[1]->Get("hyQwkick");
  }  
  // pt
  totcrossNLO = hNLO->Integral("width");
  hNLO->SetMarkerStyle(3);
  hNLO->SetMarkerColor(1);
  hNLO->SetLineColor(1);
  hNLO->SetLineWidth(2);
  scalefactorNLO = (fHist[1].Pt->Integral("width"))/totcrossNLO; // normalize to the b-bbar total cross section 
  hNLO->Scale(scalefactorNLO);
  // rapidity
  totcrossNLOY = hNLOY->Integral("width");
  hNLOY->SetMarkerStyle(3);
  hNLOY->SetMarkerColor(1);
  hNLOY->SetLineColor(1);
  hNLOY->SetLineWidth(2);
  scalefactorNLOY = (fHist[1].Y->Integral("width"))/totcrossNLOY; // normalize to the b-bbar total cross section 
  hNLOY->Scale(scalefactorNLOY);

  cPt->cd();
  hNLO->Draw("samep");
  if(qtype==4) legend1->AddEntry(hNLO,"charm from NLO", "lp");
  if(qtype==5) legend1->AddEntry(hNLO,"bottom from NLO", "lp");
  legend1->Draw();

  cY->cd();
  hNLOY->Draw("samep");
  if(qtype==4) legend2->AddEntry(hNLOY,"charm from NLO", "lp");
  if(qtype==5) legend2->AddEntry(hNLOY,"bottom from NLO", "lp");
  legend2->Draw();
}
 
//--------------------------
void DrawDifferentialYieldBoth(char *filename, char *NLOcalc){

  TFile *_file[1]; 
  _file[0] = TFile::Open(filename);
  _file[1] = TFile::Open(NLOcalc);

  setGeneralStyle();
  char configname[]="mcQA";
  enum qType {kQuark, kantiQuark, kHadron, kElec, kElec2nd, keHadron, kDeHadron};
  const Int_t fgkqType=7;

  struct hists {
        TH1F *Pt;
        TH1F *Eta;
        TH1F *Y;
  };

  hists fHist[2][7];

  TString kqTypeLabel[2][fgkqType];
  kqTypeLabel[0][kQuark]="c";
  kqTypeLabel[0][kantiQuark]="cbar";
  kqTypeLabel[0][kHadron]="cHadron";
  kqTypeLabel[0][kElec]="ce";
  kqTypeLabel[0][kElec2nd]="nulle";
  kqTypeLabel[0][keHadron]="ceHadron";
  kqTypeLabel[0][kDeHadron]="nullHadron";
  kqTypeLabel[1][kQuark]="b";
  kqTypeLabel[1][kantiQuark]="bbar";
  kqTypeLabel[1][kHadron]="bHadron";
  kqTypeLabel[1][kElec]="be";
  kqTypeLabel[1][kElec2nd]="bce";
  kqTypeLabel[1][keHadron]="beHadron";
  kqTypeLabel[1][kDeHadron]="bDeHadron";

  Int_t kColorCode[2][fgkqType];
  for(Int_t iq=0; iq<7; iq++){
    kColorCode[0][iq]=4;
    kColorCode[1][iq]=2;
  }  

  Int_t kLineStyle[fgkqType];
  kLineStyle[kQuark]=1;
  kLineStyle[kantiQuark]=1;
  kLineStyle[kHadron]=2;
  kLineStyle[kElec]=2;
  kLineStyle[kElec2nd]=1;
  kLineStyle[keHadron]=2;
  kLineStyle[kDeHadron]=1;


  cPt = new TCanvas("cPt","pT of Quark & Hadron",0,0,600,500);	
  cY = new TCanvas("cY","rapidity of Quark & Hadron",0,0,600,500);	

  cPtElec = new TCanvas("cPtElec","pT of Electron",0,0,600,500);	
  cYQElec = new TCanvas("cYElec","rapidity of Electron",0,0,600,500);	

  TList *tl = (TList *)_file[0]->Get("HFE_QA")->FindObject("MCqa");
  //count # of events
  TList *tl_result = (TList *)_file[0]->Get("HFE_Results");
  AliHFEcontainer *containerhfe = (AliHFEcontainer *) tl_result->FindObject("trackContainer");
  if(!containerhfe) {
      printf("No hfe container \n");
      return;
  }
  Int_t nEvt = (Int_t) containerhfe->GetNumberOfEvents();
  cout << "# of events " << nEvt << endl;

  Double_t scalefactor = 1./float(nEvt); //normalize with total number of event 

  TString hname; 
  for (Int_t iHQ= 0; iHQ < 2; iHQ++ ){
    for (Int_t iqType = 0; iqType < fgkqType; iqType++ ){

      //pT distribution
      hname="mcqa_Pt_";
      hname=hname+kqTypeLabel[iHQ][iqType];
      fHist[iHQ][iqType].Pt = (TH1F*)tl->FindObject(hname);
      fHist[iHQ][iqType].Pt->SetXTitle("p_{t} [GeV/c]");
      fHist[iHQ][iqType].Pt->SetYTitle("1/N_{Evt}dN/dp_{t} [1/GeVc^{-1}]");
      setDataStyle(*fHist[iHQ][iqType].Pt, kColorCode[iHQ][iqType], 3, kLineStyle[iqType]);

      CorrectFromTheWidth(fHist[iHQ][iqType].Pt); //consider pT bin size
      fHist[iHQ][iqType].Pt->Scale(scalefactor); //normalize with # of events 
      if (iqType==1) {
        fHist[iHQ][iqType].Pt->Add(fHist[iHQ][0].Pt); // Q+Qbar
	fHist[iHQ][iqType].Pt->Scale(0.5); // pair 
      }	

      //eta and rapidity distribution 	
      hname="mcqa_Eta_";
      hname=hname+kqTypeLabel[iHQ][iqType];
      fHist[iHQ][iqType].Eta = (TH1F*)tl->FindObject(hname);
      hname="mcqa_Y_";
      hname=hname+kqTypeLabel[iHQ][iqType];
      fHist[iHQ][iqType].Y = (TH1F*)tl->FindObject(hname);
      fHist[iHQ][iqType].Y->SetXTitle("Y");
      fHist[iHQ][iqType].Y->SetYTitle("1/N_{Evt}dN/dY");
      setDataStyle(*fHist[iHQ][iqType].Y, kColorCode[iHQ][iqType], 3, kLineStyle[iqType]);

      CorrectFromTheWidth(fHist[iHQ][iqType].Y); //consider Y bin size
      fHist[iHQ][iqType].Y->Scale(scalefactor); //normalize with # of events 
      if (iqType==1) {
        fHist[iHQ][iqType].Y->Add(fHist[iHQ][0].Y); // Q+Qbar
	fHist[iHQ][iqType].Y->Scale(0.5); // pair 
      }
    }
  }

  Int_t iorder[2] = {0,1}; 
  if (fHist[0][1].Pt->GetEntries() < fHist[1][1].Pt->GetEntries()) {iorder[0]=1; iorder[1]=0;};

  cPt->cd();
  setPadStyle(2,gPad);
  fHist[iorder[0]][1].Pt->Draw();
  fHist[iorder[1]][1].Pt->Draw("same");
  fHist[iorder[0]][2].Pt->Draw("same");
  fHist[iorder[1]][2].Pt->Draw("same");

  cY->cd();
  setPadStyle(2,gPad);
  fHist[iorder[0]][1].Y->Draw();
  fHist[iorder[1]][1].Y->Draw("same");
  fHist[iorder[0]][2].Y->Draw("same");
  fHist[iorder[1]][2].Y->Draw("same");

  cPtElec->cd();
  setPadStyle(2,gPad);
  fHist[iorder[0]][3].Pt->Draw();
  fHist[iorder[1]][3].Pt->Draw("same");
  fHist[1][4].Pt->Draw("same");

  cYElec->cd();
  setPadStyle(2,gPad);
  fHist[iorder[0]][3].Y->Draw();
  fHist[iorder[1]][3].Y->Draw("same");
  fHist[1][4].Y->Draw("same");


  TLegend *legend1 = new TLegend(0.50,0.77,0.98,0.99,"");
  setLegendStyle(*legend1,1);
  for (Int_t iHQ= 0; iHQ < 2; iHQ++ ){
    legend1->AddEntry(fHist[iHQ][0].Pt,kqTypeLabel[iHQ][kQuark]+"#bar{"+kqTypeLabel[iHQ][kQuark]+"}"+" pair", "l");
    legend1->AddEntry(fHist[iHQ][2].Pt,kqTypeLabel[iHQ][kQuark]+" hadrons from "+kqTypeLabel[iHQ][kQuark], "l");
  }
  TLegend *legend1_ = new TLegend(0.65,0.42,0.90,0.72,"");
  setLegendStyle(*legend1_,1);
  legend1_->AddEntry(fHist[0][0].Pt,"Yield per event", "");
  legend1_->AddEntry(fHist[0][0].Pt,"(0.1<p_{t}<20)", "");
  for (Int_t iHQ= 0; iHQ < 2; iHQ++ ){
    legend1_->AddEntry(fHist[iHQ][0].Pt,Form("%1.4f",fHist[iHQ][1].Pt->Integral("width")), "l");
    legend1_->AddEntry(fHist[iHQ][2].Pt,Form("%1.4f",fHist[iHQ][2].Pt->Integral("width")), "l");
  }

  TLegend *legend2 = new TLegend(0.50,0.77,0.98,0.99,"");
  setLegendStyle(*legend2,1);
  for (Int_t iHQ= 0; iHQ < 2; iHQ++ ){
    legend2->AddEntry(fHist[iHQ][3].Pt,"elec. from "+kqTypeLabel[iHQ][kQuark]+"->e", "l");
    if(iHQ==1) legend2->AddEntry(fHist[iHQ][4].Pt,"elec. from "+kqTypeLabel[iHQ][kQuark]+"->c->e", "l");
  }
  TLegend *legend2_ = new TLegend(0.65,0.42,0.90,0.72,"");
  setLegendStyle(*legend2_,1);
  legend2_->AddEntry(fHist[0][0].Pt,"Yield per event", "");
  legend2_->AddEntry(fHist[0][0].Pt,"(0.1<p_{t}<20)", "");
  for (Int_t iHQ= 0; iHQ < 2; iHQ++ ){
    legend2_->AddEntry(fHist[iHQ][3].Pt,Form("%1.4f",fHist[iHQ][3].Pt->Integral("width")), "l");
    if(iHQ==1) legend2_->AddEntry(fHist[iHQ][4].Pt,Form("%1.4f",fHist[iHQ][4].Pt->Integral("width")), "l");
  }


  TLegend *legend3 = new TLegend(0.26,0.15,0.74,0.37,"");
  setLegendStyle(*legend3,1);
  for (Int_t iHQ= 0; iHQ < 2; iHQ++ ){
    legend3->AddEntry(fHist[iHQ][0].Y,kqTypeLabel[iHQ][kQuark]+"#bar{"+kqTypeLabel[iHQ][kQuark]+"}"+" pair", "l");
    legend3->AddEntry(fHist[iHQ][2].Y,kqTypeLabel[iHQ][kQuark]+" hadrons from "+kqTypeLabel[iHQ][kQuark], "l");
  }
  TLegend *legend3_ = new TLegend(0.70,0.77,0.98,0.99,"");
  setLegendStyle(*legend3_,1);
  legend3_->AddEntry(fHist[0][0].Y,"Yield per event", "");
  legend3_->AddEntry(fHist[0][0].Y,"(all p_{t},-7.5<Y<7.5)", "");
  for (Int_t iHQ= 0; iHQ < 2; iHQ++ ){
    legend3_->AddEntry(fHist[iHQ][0].Y,Form("%1.4f",fHist[iHQ][1].Y->Integral("width")), "l");
    legend3_->AddEntry(fHist[iHQ][2].Y,Form("%1.4f",fHist[iHQ][2].Y->Integral("width")), "l");
  }

  TLegend *legend4 = new TLegend(0.30,0.15,0.70,0.27,"");
  setLegendStyle(*legend4,1);
  for (Int_t iHQ= 0; iHQ < 2; iHQ++ ){
    legend4->AddEntry(fHist[iHQ][3].Y,"elec. from "+kqTypeLabel[iHQ][kQuark]+"->e", "l");
    if(iHQ==1) legend4->AddEntry(fHist[iHQ][4].Y,"elec. from "+kqTypeLabel[iHQ][kQuark]+"->c->e", "l");
  }
  TLegend *legend4_ = new TLegend(0.70,0.82,0.98,0.99,"");
  setLegendStyle(*legend4_,1);
  legend4_->AddEntry(fHist[0][0].Y,"Yield per event", "");
  legend4_->AddEntry(fHist[0][0].Y,"(all p_{t},-7.5<Y<7.5)", "");
  for (Int_t iHQ= 0; iHQ < 2; iHQ++ ){
    legend4_->AddEntry(fHist[iHQ][3].Y,Form("%1.4f",fHist[iHQ][3].Y->Integral("width")), "l");
    if(iHQ==1) legend4_->AddEntry(fHist[iHQ][4].Y,Form("%1.4f",fHist[iHQ][4].Y->Integral("width")), "l");
  }


  cPt->cd();
  legend1->Draw();
  legend1_->Draw();

  cPtElec->cd();
  legend2->Draw();
  legend2_->Draw();

  cY->cd();
  legend3->Draw();
  legend3_->Draw();

  cYElec->cd();
  legend4->Draw();
  legend4_->Draw();

  TH1F *hNLOc = (TH1F*)_file[1]->Get("hpt_c");
  TH1F *hNLOb = (TH1F*)_file[1]->Get("hpt_b");
  Double_t ctotcrossNLO = hNLOc->Integral("width");
  Double_t btotcrossNLO = hNLOb->Integral("width");
  Double_t cscalefactorNLO = (fHist[0][1].Pt->Integral("width"))/ctotcrossNLO; // normalize to the b-bbar total cross section 
  Double_t bscalefactorNLO = (fHist[1][1].Pt->Integral("width"))/btotcrossNLO; // normalize to the b-bbar total cross section 
  hNLOc->Scale(cscalefactorNLO);
  hNLOb->Scale(bscalefactorNLO);
  cPt->cd();
  hNLOc->SetMarkerStyle(28);
  hNLOb->SetMarkerStyle(3);
  hNLOc->SetMarkerColor(1);
  hNLOb->SetMarkerColor(1);
  hNLOc->Draw("samep");
  hNLOb->Draw("samep");
  legend1->AddEntry(hNLOc,"charm from NLO", "p");
  legend1->AddEntry(hNLOb,"bottom from NLO", "p");
  legend1->Draw();

}


//--------------------------
void DrawDifferentialYieldAcc(char *filename){

  TFile *_file[1]; 
  _file[0] = TFile::Open(filename);

  setGeneralStyle();
  char configname[]="mcQA";
  enum qType {kQuark, kantiQuark, kHadron, kElec, kElec2nd, keHadron, kDeHadron};
  const Int_t fgkqType=7;

  struct hists {
        TH1F *Pt;
        TH1F *Eta;
        TH1F *Y;
        TH1F *barrel_Pt;
        TH1F *barrel_Eta;
        TH1F *barrel_Y;
  };

  hists fHist[2][7];

  TString kqTypeLabel[2][fgkqType];
  kqTypeLabel[0][kQuark]="c";
  kqTypeLabel[0][kElec]="ce";
  kqTypeLabel[0][kElec2nd]="nulle";
  kqTypeLabel[1][kQuark]="b";
  kqTypeLabel[1][kElec]="be";
  kqTypeLabel[1][kElec2nd]="bce";

  Int_t kColorCode[2][fgkqType][2];
  for(Int_t iq=0; iq<7; iq++){
    kColorCode[0][iq][0]=kRed;
    kColorCode[1][iq][0]=kBlue;
  }  
  kColorCode[1][kElec2nd][0]=kBlack;

  Int_t kMarkerStyle[2][2][fgkqType];
  kMarkerStyle[0][0][kElec]=24;
  kMarkerStyle[0][1][kElec]=25;
  kMarkerStyle[1][0][kElec]=20;
  kMarkerStyle[1][1][kElec]=21;
  kMarkerStyle[0][1][kElec2nd]=25;
  kMarkerStyle[1][1][kElec2nd]=21;

  Int_t kLineStyle[fgkqType];
  kLineStyle[kElec]=1;
  kLineStyle[kElec2nd]=1;

  cPtElec = new TCanvas("cPtElec","pT of Electron",0,0,600,500);	
  cYQElec = new TCanvas("cYElec","rapidity of Electron",0,0,600,500);	

  TList *tl = (TList *)_file[0]->Get("HFE_QA")->FindObject("MCqa");
  //count # of events
  TList *tl_result = (TList *)_file[0]->Get("HFE_Results");
  AliHFEcontainer *containerhfe = (AliHFEcontainer *) tl_result->FindObject("trackContainer");
  if(!containerhfe) {
      printf("No hfe container \n");
      return;
  }
  Int_t nEvt = (Int_t) containerhfe->GetNumberOfEvents();
  cout << "# of events " << nEvt << endl;

  Double_t scalefactor = 1./float(nEvt); //normalize with total number of event 

  TString hname; 
  for (Int_t iHQ= 0; iHQ < 2; iHQ++ ){
    for (Int_t iqType = 3; iqType < 5; iqType++ ){

      //pT distribution
      hname="mcqa_Pt_";
      hname=hname+kqTypeLabel[iHQ][iqType];
      fHist[iHQ][iqType].Pt = (TH1F*)tl->FindObject(hname);
      fHist[iHQ][iqType].Pt->SetXTitle("p_{t} [GeV/c]");
      fHist[iHQ][iqType].Pt->SetYTitle("1/N_{Evt}dN/dp_{t} [1/GeVc^{-1}]");
      setDataStyle(*fHist[iHQ][iqType].Pt, kColorCode[iHQ][iqType][0], 1, kLineStyle[iqType]);
      fHist[iHQ][iqType].Pt->SetMarkerStyle(kMarkerStyle[0][iHQ][iqType]);

      fHist[iHQ][iqType].Pt->Sumw2(); 
      CorrectFromTheWidth(fHist[iHQ][iqType].Pt); //consider pT bin size
      fHist[iHQ][iqType].Pt->Scale(scalefactor); //normalize with # of events 
      fHist[iHQ][iqType].Pt->Scale(0.5); //to get (e+e-)/2 

      hname="mcqa_barrel_Pt_";
      hname=hname+kqTypeLabel[iHQ][iqType];
      fHist[iHQ][iqType].barrel_Pt = (TH1F*)tl->FindObject(hname);
      fHist[iHQ][iqType].barrel_Pt->SetXTitle("p_{t} [GeV/c]");
      fHist[iHQ][iqType].barrel_Pt->SetYTitle("1/N_{Evt}dN/dp_{t} [1/GeVc^{-1}]");
      setDataStyle(*fHist[iHQ][iqType].barrel_Pt, kColorCode[iHQ][iqType][0], 1, kLineStyle[iqType]);
      fHist[iHQ][iqType].barrel_Pt->SetMarkerStyle(kMarkerStyle[1][iHQ][iqType]);

      fHist[iHQ][iqType].barrel_Pt->Sumw2();
      CorrectFromTheWidth(fHist[iHQ][iqType].barrel_Pt); //consider pT bin size
      fHist[iHQ][iqType].barrel_Pt->Scale(scalefactor); //normalize with # of events
      fHist[iHQ][iqType].barrel_Pt->Scale(0.5); //to get (e+e-)/2 


      //eta and rapidity distribution 	
      hname="mcqa_Eta_";
      hname=hname+kqTypeLabel[iHQ][iqType];
      fHist[iHQ][iqType].Eta = (TH1F*)tl->FindObject(hname);
      hname="mcqa_Y_";
      hname=hname+kqTypeLabel[iHQ][iqType];
      fHist[iHQ][iqType].Y = (TH1F*)tl->FindObject(hname);
      fHist[iHQ][iqType].Y->SetXTitle("Y");
      fHist[iHQ][iqType].Y->SetYTitle("1/N_{Evt}dN/dY");
      setDataStyle(*fHist[iHQ][iqType].Y, kColorCode[iHQ][iqType][0], 1, kLineStyle[iqType]);
      fHist[iHQ][iqType].Y->SetMarkerStyle(kMarkerStyle[0][1][iqType]);

      fHist[iHQ][iqType].Y->Rebin(2);
      fHist[iHQ][iqType].Y->Sumw2();
      CorrectFromTheWidth(fHist[iHQ][iqType].Y); //consider Y bin size
      fHist[iHQ][iqType].Y->Scale(scalefactor); //normalize with # of events
      fHist[iHQ][iqType].Y->Scale(0.5); //to get (e+e-)/2

      hname="mcqa_barrel_Eta_";
      hname=hname+kqTypeLabel[iHQ][iqType];
      fHist[iHQ][iqType].barrel_Eta = (TH1F*)tl->FindObject(hname);
      hname="mcqa_barrel_Y_";
      hname=hname+kqTypeLabel[iHQ][iqType];
      fHist[iHQ][iqType].barrel_Y = (TH1F*)tl->FindObject(hname);
      fHist[iHQ][iqType].barrel_Y->SetXTitle("Y");
      fHist[iHQ][iqType].barrel_Y->SetYTitle("1/N_{Evt}dN/dY");
      setDataStyle(*fHist[iHQ][iqType].barrel_Y, kColorCode[iHQ][iqType][0], 1, kLineStyle[iqType]);
      fHist[iHQ][iqType].barrel_Y->SetMarkerStyle(kMarkerStyle[1][iHQ][iqType]);

      fHist[iHQ][iqType].barrel_Y->Rebin(2);
      fHist[iHQ][iqType].barrel_Y->Sumw2();
      CorrectFromTheWidth(fHist[iHQ][iqType].barrel_Y); //consider Y bin size
      fHist[iHQ][iqType].barrel_Y->Scale(scalefactor); //normalize with # of events
      fHist[iHQ][iqType].barrel_Y->Scale(0.5); //to get (e+e-)/2
    }
  }

  Int_t iorder[2] = {0,1}; 
  if (fHist[0][3].Pt->GetEntries() < fHist[1][3].Pt->GetEntries()) {iorder[0]=1; iorder[1]=0;};

  cPtElec->cd();
  setPadStyle(2,gPad);
  fHist[iorder[0]][3].Pt->SetAxisRange(0,10,"X");
  fHist[iorder[0]][3].Pt->Draw("p");
  fHist[iorder[1]][3].Pt->Draw("samep");
  fHist[1][4].Pt->Draw("samep");
  fHist[iorder[0]][3].barrel_Pt->Draw("samep");
  fHist[iorder[1]][3].barrel_Pt->Draw("samep");
  fHist[1][4].barrel_Pt->Draw("samep");

  cYElec->cd();
  setPadStyle(2,gPad);
  fHist[iorder[0]][3].Y->Draw();
  fHist[iorder[1]][3].Y->Draw("same");
  fHist[1][4].Y->Draw("same");
  fHist[iorder[0]][3].barrel_Y->Draw("samelp");
  fHist[iorder[1]][3].barrel_Y->Draw("samelp");
  fHist[1][4].barrel_Y->Draw("samelp");


  TLegend *legend2 = new TLegend(0.62,0.74,0.98,0.99,"");
  setLegendStyle(*legend2,1);
  for (Int_t iHQ= 0; iHQ < 2; iHQ++ ){
    legend2->AddEntry(fHist[iHQ][3].Pt,"(e+e^{-})/2, "+kqTypeLabel[iHQ][kQuark]+"->e", "p");
    legend2->AddEntry(fHist[iHQ][3].barrel_Pt,"(e+e^{-})/2, "+kqTypeLabel[iHQ][kQuark]+"->e |#eta<0.9|", "p");
    if(iHQ==1) legend2->AddEntry(fHist[iHQ][4].Pt,"(e+e^{-})/2, "+kqTypeLabel[iHQ][kQuark]+"->c->e", "p");
    if(iHQ==1) legend2->AddEntry(fHist[iHQ][4].barrel_Pt,"(e+e^{-})/2, "+kqTypeLabel[iHQ][kQuark]+"->c->e |#eta<0.9|", "p");
  }
  TLegend *legend2_ = new TLegend(0.62,0.42,0.98,0.72,"");
  setLegendStyle(*legend2_,1);
  legend2_->AddEntry(fHist[0][0].Pt,"(e+e^{-})/2, yield per event", "");
  legend2_->AddEntry(fHist[0][0].Pt,"(0.1<p_{t}<20)", "");
  for (Int_t iHQ= 0; iHQ < 2; iHQ++ ){
    legend2_->AddEntry(fHist[iHQ][3].Pt,Form("%1.4f",fHist[iHQ][3].Pt->Integral("width")), "p");
    legend2_->AddEntry(fHist[iHQ][3].barrel_Pt,Form("%1.4f",fHist[iHQ][3].barrel_Pt->Integral("width")), "p");
    if(iHQ==1) legend2_->AddEntry(fHist[iHQ][4].Pt,Form("%1.4f",fHist[iHQ][4].Pt->Integral("width")), "p");
    if(iHQ==1) legend2_->AddEntry(fHist[iHQ][4].barrel_Pt,Form("%1.4f",fHist[iHQ][4].barrel_Pt->Integral("width")), "p");
  }


  TLegend *legend4 = new TLegend(0.30,0.12,0.70,0.36,"");
  setLegendStyle(*legend4,1);
  for (Int_t iHQ= 0; iHQ < 2; iHQ++ ){
    legend4->AddEntry(fHist[iHQ][3].Y,"(e+e^{-})/2, "+kqTypeLabel[iHQ][kQuark]+"->e", "p");
    legend4->AddEntry(fHist[iHQ][3].barrel_Y,"(e+e^{-})/2, "+kqTypeLabel[iHQ][kQuark]+"->e |#eta<0.9|", "p");
    if(iHQ==1) legend4->AddEntry(fHist[iHQ][4].Y,"(e+e^{-})/2, "+kqTypeLabel[iHQ][kQuark]+"->c->e", "p");
    if(iHQ==1) legend4->AddEntry(fHist[iHQ][4].barrel_Y,"(e+e^{-})/2, "+kqTypeLabel[iHQ][kQuark]+"->c->e |#eta<0.9|", "p");
  }
  TLegend *legend4_ = new TLegend(0.65,0.72,0.98,0.99,"");
  setLegendStyle(*legend4_,1);
  legend4_->AddEntry(fHist[0][0].Y,"(e+e^{-})/2, yield per event", "");
  legend4_->AddEntry(fHist[0][0].Y,"(all p_{t})", "");
  for (Int_t iHQ= 0; iHQ < 2; iHQ++ ){
    legend4_->AddEntry(fHist[iHQ][3].Y,Form("%1.4f",fHist[iHQ][3].Y->Integral("width")), "p");
    legend4_->AddEntry(fHist[iHQ][3].barrel_Y,Form("%1.4f",fHist[iHQ][3].barrel_Y->Integral("width")), "p");
    if(iHQ==1) legend4_->AddEntry(fHist[iHQ][4].Y,Form("%1.4f",fHist[iHQ][4].Y->Integral("width")), "p");
    if(iHQ==1) legend4_->AddEntry(fHist[iHQ][4].barrel_Y,Form("%1.4f",fHist[iHQ][4].barrel_Y->Integral("width")), "p");
  }


  cPtElec->cd();
  legend2->Draw();
  legend2_->Draw();

  cYElec->cd();
  legend4->Draw();
  legend4_->Draw();

}

//--------------------------
void DrawDifferentialXsectionBoth(char *filename, char *NLOcalc){

  TFile *_file[1]; 
  _file[0] = TFile::Open(filename);
  _file[1] = TFile::Open(NLOcalc);

  setGeneralStyle();
  char configname[]="mcQA";
  enum qType {kQuark, kantiQuark, kHadron, kElec, kElec2nd, keHadron, kDeHadron};
  const Int_t fgkqType=7;

  struct hists {
        TH1F *Pt;
  };

  hists fHist[2][7];

  TString kqTypeLabel[2][fgkqType];
  kqTypeLabel[0][kQuark]="c";
  kqTypeLabel[0][kantiQuark]="cbar";
  kqTypeLabel[0][kHadron]="cHadron";
  kqTypeLabel[0][kElec]="ce";
  kqTypeLabel[0][kElec2nd]="nulle";
  kqTypeLabel[0][keHadron]="ceHadron";
  kqTypeLabel[0][kDeHadron]="nullHadron";
  kqTypeLabel[1][kQuark]="b";
  kqTypeLabel[1][kantiQuark]="bbar";
  kqTypeLabel[1][kHadron]="bHadron";
  kqTypeLabel[1][kElec]="be";
  kqTypeLabel[1][kElec2nd]="bce";
  kqTypeLabel[1][keHadron]="beHadron";
  kqTypeLabel[1][kDeHadron]="bDeHadron";

  Int_t kColorCode[2][fgkqType];
  for(Int_t iq=0; iq<7; iq++){
    kColorCode[0][iq]=4;
    kColorCode[1][iq]=2;
  }  

  Int_t kLineStyle[fgkqType];
  kLineStyle[kQuark]=1;
  kLineStyle[kantiQuark]=1;
  kLineStyle[kHadron]=2;
  kLineStyle[kElec]=2;
  kLineStyle[kElec2nd]=1;
  kLineStyle[keHadron]=2;
  kLineStyle[kDeHadron]=1;


  cPt = new TCanvas("cPt","pT of Quark & Hadron",0,0,600,500);	

  Double_t totcrossNLO[2];
  TH1F *hNLOc = (TH1F*)_file[1]->Get("hpt_c");
  TH1F *hNLOb = (TH1F*)_file[1]->Get("hpt_b");
  totcrossNLO[0] = hNLOc->Integral("width");
  totcrossNLO[1] = hNLOb->Integral("width");
  cPt->cd();
  setPadStyle(2,gPad);
  hNLOc->SetMarkerStyle(28);
  hNLOb->SetMarkerStyle(3);
  hNLOc->SetMarkerColor(1);
  hNLOb->SetMarkerColor(1);
  hNLOc->Draw("p");
  hNLOb->Draw("samep");


  TList *tl = (TList *)_file[0]->Get("HFE_QA")->FindObject("MCqa");
  //count # of events
  TList *tl_result = (TList *)_file[0]->Get("HFE_Results");
  AliHFEcontainer *containerhfe = (AliHFEcontainer *) tl_result->FindObject("trackContainer");
  if(!containerhfe) {
      printf("No hfe container \n");
      return;
  }
  Int_t nEvt = (Int_t) containerhfe->GetNumberOfEvents();
  cout << "# of events " << nEvt << endl;

  Double_t scalefactor = 1./float(nEvt); //normalize with total number of event 

  TString hname; 
  for (Int_t iHQ= 0; iHQ < 2; iHQ++ ){
    for (Int_t iqType = 0; iqType < 3; iqType++ ){

      //pT distribution
      hname="mcqa_Pt_";
      hname=hname+kqTypeLabel[iHQ][iqType];
      fHist[iHQ][iqType].Pt = (TH1F*)tl->FindObject(hname);
      fHist[iHQ][iqType].Pt->SetXTitle("p_{t} [GeV/c]");
      fHist[iHQ][iqType].Pt->SetYTitle("1/N_{Evt}dN/dp_{t} [1/GeVc^{-1}]");
      setDataStyle(*fHist[iHQ][iqType].Pt, kColorCode[iHQ][iqType], 3, kLineStyle[iqType]);

      CorrectFromTheWidth(fHist[iHQ][iqType].Pt); //consider pT bin size
      fHist[iHQ][iqType].Pt->Scale(scalefactor); //normalize with # of events 
      if (iqType==1) {
        fHist[iHQ][iqType].Pt->Add(fHist[iHQ][0].Pt); //Q+Qbar
        fHist[iHQ][iqType].Pt->Scale(0.5); // pari
        Double_t xfactor = fHist[iHQ][iqType].Pt->Integral("width"); 
        fHist[iHQ][iqType].Pt->Scale(totcrossNLO[iHQ]/xfactor); //consider pT bin size 
      }
    }
  }

  Int_t iorder[2] = {0,1}; 
  if (fHist[0][1].Pt->GetEntries() < fHist[1][1].Pt->GetEntries()) {iorder[0]=1; iorder[1]=0;};

  cPt->cd();
  setPadStyle(2,gPad);
  fHist[iorder[0]][1].Pt->Draw("same");
  fHist[iorder[1]][1].Pt->Draw("same");

  TLegend *legend1 = new TLegend(0.50,0.77,0.98,0.99,"");
  setLegendStyle(*legend1,1);
  for (Int_t iHQ= 0; iHQ < 2; iHQ++ ){
    legend1->AddEntry(fHist[iHQ][0].Pt,kqTypeLabel[iHQ][kQuark]+"#bar{"+kqTypeLabel[iHQ][kQuark]+"}"+" pair", "l");
  }
  legend1->AddEntry(hNLOc,"charm from NLO", "p");
  legend1->AddEntry(hNLOb,"bottom from NLO", "p");

  cPt->cd();
  legend1->Draw();

}
//--------------------------
void setDataStyle(TH1F &h, Int_t lc, Int_t lw, Int_t ls){

  h.SetLineColor(lc);
  h.SetLineWidth(lw);
  h.SetLineStyle(ls);
  h.SetMarkerColor(lc);

}

//--------------------------
void setLegendStyle(TLegend &legend, Int_t bs){

  legend.SetBorderSize(bs);
  legend.SetFillColor(0);
  legend.SetTextFont(132);
  legend.SetTextSize(0.04);
  legend.SetMargin(0.15);

}

//--------------------------
void setPadStyle(Int_t lvl, TPad *pad){

  pad->SetLogy();
  if(lvl>0) gPad->SetGridy();
  if(lvl>1) gPad->SetGridx();

}

//--------------------------
void setGeneralStyle(){

  gStyle->SetPalette(1);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderSize(10);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.04,"X");
  gStyle->SetTitleSize(0.04,"Y");
  gStyle->SetTitleFont(132,"X");
  gStyle->SetTitleFont(132,"Y");
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(0.9);
  gStyle->SetLabelFont(132,"X");
  gStyle->SetLabelFont(132,"Y");
  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");

  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetLineWidth(2);

}

void CorrectFromTheWidth(TH1F *h1) const {
  //
  // Correct from the width of the bins --> dN/dp_{T} (GeV/c)^{-1}
  //

  TAxis *axis = h1->GetXaxis();
  Int_t nbinX = h1->GetNbinsX();
  
  for(Int_t i = 1; i <= nbinX; i++) {
  
    Double_t width = axis->GetBinWidth(i);
    Double_t content = h1->GetBinContent(i);
    Double_t error = h1->GetBinError(i);
    h1->SetBinContent(i,content/width);
    h1->SetBinError(i,error/width);
  }
  
}

