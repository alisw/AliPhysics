/// \file guiAlias.C

/// Make default Aliases for guiTime:
///
/// 1. Run the guiTime
/// ~~~
/// guiTime
/// ~~~
///
/// 2. Define aliases
/// ~~~{.cpp}
/// .L $ALICE_ROOT/TPC/CalibMacros/guiAlias.C
/// guiAlias();
/// ~~~
///
/// 3. Use aliases inside the guiTime
/// You can use them as custom cuts
///
/// ~~~{.cpp}
/// // browse special streamers
/// gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
/// gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
/// AliXRDPROOFtoolkit tool;
/// TChain * chainDCS = tool.MakeChainRandom("time.txt","dcs",0,10000);
/// TChain * chainCTP = tool.MakeChainRandom("time.txt","ctp",0,10000);
/// TChain * chainAlign = tool.MakeChainRandom("time.txt","align",0,10000);
/// ~~~

TObjArray *picArray = new TObjArray;  

TTree * guiTree =  guiTime->GetChain();
void SetStyle();
void guiAlias(){
  guiTree =  guiTime->GetChain();
  MakeAliasCE(4);
  MakeAliasLT(4);
  MakeAliasCosmic(4);   
  SetGoofieAlias();
  SetStyle();
}

void SetStyle(){ 
  Float_t mx0=0.15, mx1=0.05, my0=0.15, my1=0.1;
  guiTime->GetCanvas()->SetTicks(1,1);
  guiTime->GetCanvas()->SetMargin(mx0,mx1,my0,my1);
  gStyle->SetTitleYSize(0.03);
  gStyle->SetTitleXSize(0.03);
  gStyle->SetTitleXOffset(2);
  gStyle->SetTitleYOffset(6);
}

void MakeAliasCE(Double_t deltaT){
  /// Aliases cuts for CE

  guiTree->SetAlias("ceCut0", "tdriftCE.fElements[72]>100 && tdriftCE.fElements[73]>100");
  guiTree->SetAlias("dceCutTime", Form("sqrt(dcea^2+dcec^2)<%f",deltaT*3600));
  guiTree->SetAlias("ceCut","dceCutTime&&ceCut0");
};

void MakeAliasLT(Double_t deltaT){
  guiTree->SetAlias("ltCut", Form("sqrt(dla^2+dlc^2)<%f", deltaT*3600)); 
}

void MakeAliasCosmic(Double_t deltaT){
  guiTree->SetAlias("cosmicCut", Form("abs(dcosmic)<%f", deltaT*3600));
  guiTree->SetAlias("itsCut", Form("((dits!=0)&&abs(dits)<%f)", deltaT*3600));
}


void SetGoofieAlias(){
  /// goofie aliases

  guiTree->SetAlias("ptrelG","(goofie.fElements[17]/0.3426-1)");
  guiTree->SetAlias("vdriftGN","goofie.fElements[3]/(1+ptrelG)");
  guiTree->SetAlias("goCut","goofie.fElements[3]>0");
  //
}
//
// Make default plots
//


void DrawLaserDrift(){
  /// laser calibration default picture
  /// Data are filtered

  //
  // draw laser residuals A side -C side - when it is defined
  //
  TH1 * his=0;
  guiTree->Draw("CEgrDriftA-CEgrDriftC","CEgrDriftA>0&&CEgrDriftC>0");
  his=(TH1*)htemp->Clone();
  his->SetDirectory(0);
  his->SetTitle("#Delta CE time");
  his->SetName("#Delta CE time");
  his->GetXaxis()->SetTitle("#Delta T (time bin)");
  his->Draw();
  picArray->AddLast(his);
  //
  // laser drift CE correction
  //
  guiTree->Draw("100*(vdriftCEA-vdriftCEC)","CEgrDriftA>0&&CEgrDriftC>0");
  his=(TH1*)htemp->Clone();
  his->SetDirectory(0);
  his->SetTitle("#Delta CE drift");
  his->SetName("#Delta CE drift");
  his->GetXaxis()->SetTitle("#Delta v_{drift} (%) (laser A side - C side)");
  his->Draw();
  picArray->AddLast(his);
  //
  // laser track drift correction
  //
  guiTree->Draw("100*(vdriftLTA-vdriftLTC)","abs(dlaserA)<500&&abs(dlaserC)<500&&vdriftLTA!=0&&vdriftLTC!=0");
  his=(TH1*)htemp->Clone();
  his->SetDirectory(0);
  his->SetTitle("#Delta LaserTracks drift");
  his->SetName("#Delta LaserTracks drift");
  his->GetXaxis()->SetTitle("#Delta v_{drift} (%) (laser A side - C side)");
  his->Draw();
  picArray->AddLast(his);
  //
  // laser track drift correction: time
  //
  guiTree->Draw("100*(vdriftLTA-vdriftLTC):time","abs(dlaserA)<500&&abs(dlaserC)<500&&vdriftLTA!=0&&vdriftLTC!=0","colz");
  his=(TH1*)htemp->Clone();  
  his->SetDirectory(0);
  his->Draw();
  his->SetTitle("#Delta LaserTracks drift");
  his->SetName("#Delta LaserTracks drift");
  his->GetYaxis()->SetTitle("#Delta v_{drift} (%) (laser A side - C side)");
  his->GetXaxis()->SetTitle("time");
  his->GetXaxis()->SetTimeDisplay(kTRUE);
  his->Draw();
  picArray->AddLast(his);
  //
  //
  //
  guiTree->Draw("250*(vlaserA0-vlaserC0)","abs(dlaserA)<500&&abs(dlaserC)<500&&vdriftLTA!=0&&vdriftLTC!=0");
  his=(TH1*)htemp->Clone();
  his->SetDirectory(0);
  his->SetTitle("#Delta LaserTracks offset");
  his->SetName("#Delta LaserTracks offset");
  his->GetXaxis()->SetTitle("#Delta z0 (cm) (laser A side - C side)");
  his->Draw();
  picArray->AddLast(his);
}


void DrawITSVD(){
  /// ITS/TPC drift velocity correction

  guiTree->Draw("100*(ALIGN_ITSP_TPC_DRIFTVD-ALIGN_ITSM_TPC_DRIFTVD)","abs(dits)<3000&&ALIGN_ITSM_TPC_DRIFTVD!=0");
  his=(TH1*)htemp->Clone();
  his->SetDirectory(0);
  his->SetTitle("#Delta ITS/TPC drift correction positive and negative extrapolation");
  his->SetName("#Delta ITS/TPC drift correction");
  his->GetXaxis()->SetTitle("#Delta v_{dcorr} (%)");
  his->Draw();
  picArray->AddLast(his);
  //
  // comparison with laser
  //
  guiTree->Draw("100*(vdriftLTA-vdriftITS)","abs(dlaserA)<3600&&vdriftLTA!=0&&abs(dits)<3600");
  his=(TH1*)htemp->Clone();
  his->SetDirectory(0);
  his->SetTitle("#Delta ITS/TPC drift correction - Laser track corr");
  his->SetName("#Delta ITS/TPC drift correction - laser track corr");
  his->GetXaxis()->SetTitle("#Delta v_{dcorr} (%)");
  his->Draw();
  picArray->AddLast(his);
  //
  guiTree->Draw("100*(vdriftLTA-ALIGN_ITS_TPC_DRIFTVD)","abs(dlaserA)<3600&&vdriftLTA!=0&&abs(dits)<3600");
  his=(TH1*)htemp->Clone();
  his->SetDirectory(0);
  his->SetTitle("#Delta ITS/TPC non smmoth drift correction - Laser track corr");
  his->SetName("#Delta ITS/TPC non smooth drift correction - laser track corr");
  his->GetXaxis()->SetTitle("#Delta v_{dcorr} (%)");
  his->Draw();
  picArray->AddLast(his);
  //
  guiTree->Draw("100*(vdriftLTA-vdriftP)","abs(dlaserA)<1800&&vdriftLTA!=0&&abs(dp)<1800");
  his=(TH1*)htemp->Clone();
  his->SetDirectory(0);
  his->SetTitle("#Delta TPC/TPC - Laser track corr");
  his->SetName("#Delta TPC/TPC  - Laser track corr");
  his->GetXaxis()->SetTitle("#Delta v_{dcorr} (%)");
  his->Draw();
  picArray->AddLast(his);
  
}


void CEdrift(){
  ///

  guiTree->SetAlias("vdriftCE0","250/CEgrDriftA");
  guiTree->SetAlias("vdriftCE1","vdriftCE0/(1+ptrel0)");
}
