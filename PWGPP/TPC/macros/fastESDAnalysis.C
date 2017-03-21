/*
  Fast ESD analysis
  - make standard histograms using  functionality of AliTreePlayer::MakeHistograms
  - makeReport() -  projection/fits using  AliTreePlayer::DrawHistograms
  - makeParamtrization - perforamcne parameterization to be added
  //
  .L $NOTES/JIRA/PWGPP-272/code/fastESDAnalysis.C
  .L $NOTES/aux/NimStyle.C
  NimStyleBig();
  makeReport();
  //
  aliroot -b -q  $NOTES/aux/NimStyle.C\(2\)  $NOTES/JIRA/PWGPP-272/code/fastESDAnalysis.C\(5000,0\)

  aliroot -b -q  $NOTES/aux/NimStyle.C $NOTES/JIRA/PWGPP-272/code/fastESDAnalysis.C\(5000,1\)


*/

TObjArray *hisArray=0;
TObjArray *hisV0Array=0;
TObjArray *keepArray=0;

void fastESDAnalysis(Int_t maxEntries, Int_t fillHisto=0){
  
  if (fillHisto) MakeHistograms(maxEntries);
  makeReport();
}

void makeReport(){
  TFile *f= TFile::Open("esdHisto.root");
  hisArray=(TObjArray*)f->Get("hisArray");
  hisV0Array=(TObjArray*)f->Get("hisV0Array");
  keepArray=new TObjArray;
  Chi2ReportTPC();
  Chi2ReportITS();
  Chi2ReportTRD();
  MatchingReport(); 
  PtReport();
}



void MakeHistograms(Int_t maxEntries){
  //  Int_t maxEntries=1000;
  TChain * chain = AliXRDPROOFtoolkit::MakeChainRandom("esd.list","esdTree",0,10000);
  Double_t massK0=TDatabasePDG::Instance()->GetParticle("K_S0")->Mass();
  Double_t massLambda=TDatabasePDG::Instance()->GetParticle("Lambda0")->Mass();
  chain->SetAlias("ITSOn","(Tracks[].fFlags&0x1)>0");
  chain->SetAlias("TPCOn","(Tracks[].fFlags&0x10)>0");
  chain->SetAlias("ITSRefit","((Tracks[].fFlags&0x4)>0)");
  chain->SetAlias("TPCRefit","((Tracks[].fFlags&0x40)>0)");
  chain->SetAlias("TOFOn","(Tracks[].fFlags&0x2000)>0");
  chain->SetAlias("TRDOn","(Tracks[].fFlags&0x400)>0");
  chain->SetAlias("isPrim","abs(Tracks[].fD)<3&&abs(Tracks[].fD/sqrt(Tracks[].fCdd))<6");
  chain->SetAlias("qptResol","sqrt(Tracks[].fCp.fC[14]+0)");
  chain->SetAlias("normChi2TPC","(Tracks[].fTPCchi2/Tracks[].fTPCncls)");
  chain->SetAlias("normChi2ITS","sqrt(Tracks[].fITSchi2/Tracks[].fITSncls)");
  chain->SetAlias("normChi2TRD","sqrt(Tracks[].fTRDchi2/Tracks[].fTRDntracklets)");
  chain->SetAlias("qPt","Tracks[].fP[4]");
  chain->SetAlias("tgl","Tracks[].fP[3]");
  chain->SetAlias("mdEdx","50./Tracks[].fTPCsignal");
  chain->SetAlias("alpha","Tracks[].fAlpha");
  chain->SetAlias("notKink","Tracks[].fKinkIndexes[0]==0");
  chain->SetAlias("dMassK0",TString::Format("(V0s[].GetEffMass(2,2)-%.4f)",massK0).Data());
  chain->SetAlias("dMassL",TString::Format("(V0s[].GetEffMass(4,2)-%.4f)",massLambda).Data());  
  //
  TString hisString="";
  // PT
  hisString+="Tracks.Pt():tgl:#TPCOn&&isPrim>>hisPtAll(100,0,20,10,-1,1);"; 
  hisString+="Tracks.Pt():tgl:#TPCOn&&ITSOn&&isPrim>>hisPtTPCITS(100,0,20,10,-1,1);"; 
  hisString+="Tracks.Pt():tgl:#TPCOn&&ITSOn&&TRDOn&&isPrim>>hisPtTPCITSTRD(100,0,20,10,-1,1);"; 
  hisString+="Tracks.Pt():tgl:#TPCOn&&ITSOn&&TOFOn&&isPrim>>hisPtTPCITSTOF(100,0,20,10,-1,1);"; 
  // qpt
  hisString+="qPt:tgl:alpha:#TPCOn&&isPrim>>hisQPtAll(100,-2,2,10,-1,1,18,-3.14,2.14);"; 
  hisString+="qPt:tgl:alpha:#ITSOn&&TPCOn&&isPrim>>hisQPtTPCITS(100,-2,2,10,-1,1,18,-3.14,2.14);"; 
  hisString+="qPt:tgl:alpha:#ITSOn&&TPCOn&&TRDOn&&isPrim>>hisQPtTPCITSTRD(100,-2,2,10,-1,1,18,-3.14,2.14);";   
  // DCA 
  hisString+="Tracks[].fdTPC:abs(qPt):tgl:#ITSOn&&TPCOn&&isPrim>>hisTPCDCARTPC(100,-1,1,50,0,2,10,-1,1);";
  hisString+="Tracks[].fdTPC:abs(qPt):tgl:#ITSOn&&TRDOn&&isPrim>>hisTPCDCARTPCTRD(100,-1,1,50,0,2,10,-1,1);";
  hisString+="Tracks[].fzTPC:abs(qPt):tgl:#ITSOn&&TPCOn&&isPrim>>hisTPCDCAZTPC(100,-1,1,50,0,2,10,-1,1);";
  hisString+="Tracks[].fzTPC:abs(qPt):tgl:#ITSOn&&TRDOn&&isPrim>>hisTPCDCAZTPCTRD(100,-1,1,50,0,2,10,-1,1);";
  // Pt expected resolution
  hisString+="qptResol:abs(qPt):tgl:#ITSOn&&TPCOn&&isPrim>>hisQPtResolITSTPC(400,0,0.02,50,0,2,10,-1,1);";
  hisString+="qptResol:abs(qPt):tgl:#ITSOn&&TRDOn&&isPrim>>hisQPtResolITSTPCTRD(400,0,0.02,50,0,2,10,-1,1);";
  // Chi2 histograms ITS,TPC,TRD
  hisString+=TString::Format("normChi2ITS:tgl:qPt:mdEdx:#ITSOn&&TRDOn&&isPrim>>hisChi2ITSTglQptMdEdx(50,0,10,10,-1,1,100,-5,5,10,0,1.1);");
  hisString+=TString::Format("normChi2TPC:tgl:qPt:mdEdx:#ITSOn&&TRDOn&&isPrim>>hisChi2TPCTglQptMdEdx(50,0,10,10,-1,1,100,-5,5,10,0,1.1);");
  hisString+=TString::Format("normChi2TRD:tgl:qPt:mdEdx:#ITSOn&&TRDOn&&isPrim>>hisChi2TRDTglQptMdEdx(50,0,10,10,-1,1,100,-5,5,10,0,1.1);");
  hisString+=TString::Format("normChi2ITS:tgl:qPt:mdEdx:#ITSOn&&TRDOn&&isPrim&&notKink>>hisChi2ITSNoKinkTglQptMdEdx(50,0,10,10,-1,1,100,-5,5,10,0,1.1);");
  // Refit probability in case first pass
  hisString+="ITSRefit:qPt:tgl:#ITSOn&&TPCOn&&isPrim>>hisITSRefitAll(2,-0.5,1.5,50,-2,2,10,-1,1);";
  hisString+="TPCRefit:qPt:tgl:#ITSOn&&TPCOn&&isPrim>>hisTPCRefitAll(2,-0.5,1.5,50,-2,2,10,-1,1);";
  hisString+="ITSRefit:qPt:tgl:#ITSOn&&TPCOn&&TRDOn&&isPrim>>hisITSRefitTRDOn(2,-0.5,1.5,50,-2,2,10,-1,1);";
  hisString+="TPCRefit:qPt:tgl:#ITSOn&&TPCOn&&TRDOn&&isPrim>>hisTPCRefitTRDOn(2,-0.5,1.5,50,-2,2,10,-1,1);";
  
  
  // K0s delta mass
  TString v0String="";
  v0String+="dMassK0:V0s[].fParamP.fP[4]:V0s[].fParamP.fP[3]:V0s.fRr:#V0s[].fOnFlyStatus==1&&V0s[].PtArmV0()>0.05>>hisdMassK0(100,-0.1,0.1,10,-1,1,5,-1,1,10,0,80);";
  

  TStopwatch timer;
  hisArray = AliTreePlayer::MakeHistograms(chain, hisString, "TPCOn",0,maxEntries,10000000,15);
  timer.Print();
  TStopwatch timer;
  hisV0Array = AliTreePlayer::MakeHistograms(chain, v0String, "V0s[].fOnFlyStatus==1",0,maxEntries,10000000,15);
  timer.Print();
  TFile *f= TFile::Open("esdHisto.root","recreate");
  hisArray->Write("hisArray",TObjArray::kSingleKey);
  hisV0Array->Write("hisV0Array",TObjArray::kSingleKey);
  delete f;
}




void PtReport(){
  // pt report
  TString drawExpressionPt="";
  drawExpressionPt="[1,3]:";
  drawExpressionPt+="%Ogridx,gridy,logy;hisPtAll(10,100,0,30)(0)(err);hisPtTPCITS(10,100,0,30)(0)(err);hisPtTPCITSTRD(10,100,0,30)(0)(err);hisPtTPCITSTOF(10,100,0,30)(0)(err):";
  drawExpressionPt+="%Ogridx,gridy,logy;hisPtTPCITS(10,100,0,30)(0)(err):";
  drawExpressionPt+="%Ogridx,gridy,logy;hisPtTPCITSTRD(10,100,0,30)(0)(err):";
  drawExpressionPt+="%Ogridx,gridy,logy;hisPtTPCITSTOF(10,100,0,30)(0)(err):";
  TPad * padPt = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionPt,keepArray, 1+2+4+8);
  ((TCanvas*)padPt)->SetWindowSize(1200,800);
  padPt->SaveAs("ptReport.png");
}

void QptPtReport(){
  //
  //
  TString drawExpressionQptRes="";
  drawExpressionQptRes="[1]:";
  drawExpressionQptRes+="%Ogridx,gridy;hisQPtResolITSTPC(0,400,0,10,0,10000)(0,1)(f-mean p);hisQPtResolITSTPCTRD(0,400,0,10,0,10000)(0,1)(f-mean p)";
   TPad * padQPtRes = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionQptRes,keepArray, 1+2+4+8);
  ((TCanvas*)padQPtRes)->SetWindowSize(1200,800);
}

void DCARreport(){
  TString drawExpressionTPCDCAR="";
  drawExpressionTPCDCAR="[1,1]:";
  drawExpressionTPCDCAR+="%Ogridx,gridy;hisTPCDCARTPC(0,100,0,180,0,10000)(0,1)(f-grms p);hisTPCDCARTPCTRD(0,100,0,180,0,10000)(0,1)(f-grms p);:";

  TPad * padTPCDCAR = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionTPCDCAR,keepArray, 1+2+4+8);
  //  ((TCanvas*)pad)->SetWindowSize(1200,800);

}

void Chi2ReportTPC(){
  TString drawExpressionTPCChi2="";
  drawExpressionTPCChi2="[1,1,1]:";
  drawExpressionTPCChi2+="%Ogridx,gridy;hisChi2TPCTglQptMdEdx(0,100,0,10,20,80)(0,2)(f-mean p);:";
  drawExpressionTPCChi2+="%Ogridx,gridy;hisChi2TPCTglQptMdEdx(0,100,0,10,0,100)(0,1)(f-mean p);:";
  drawExpressionTPCChi2+="%Ogridx,gridy;hisChi2TPCTglQptMdEdx(0,100,0,10,0,100)(0,3)(f-mean p);:";
  TPad * padTPCChi2 = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionTPCChi2,keepArray, 1+2+4+8);
  padTPCChi2->SaveAs("Chi2ReportTPCDefault.png")
}

void Chi2ReportITS(){
  TString drawExpressionITSChi2="";
  drawExpressionITSChi2="[1,1,1]:";
  drawExpressionITSChi2+="%Ogridx,gridy;hisChi2ITSTglQptMdEdx(0,100,0,10,20,80)(0,2)(f-mean p);hisChi2ITSNoKinkTglQptMdEdx(0,100,0,10,20,80)(0,2)(f-mean p);:";
  drawExpressionITSChi2+="%Ogridx,gridy;hisChi2ITSTglQptMdEdx(0,100,0,10,0,100)(0,1)(f-mean p);hisChi2ITSNoKinkTglQptMdEdx(0,100,0,10,0,100)(0,1)(f-mean p);:";
  drawExpressionITSChi2+="%Ogridx,gridy;hisChi2ITSTglQptMdEdx(0,100,0,10,0,100)(0,3)(f-mean p);hisChi2ITSNoKinkTglQptMdEdx(0,100,0,10,0,100)(0,3)(f-mean p);:";
  TPad * padITSChi2 = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionITSChi2,keepArray, 1+2+4+8);
  padITSChi2->SaveAs("Chi2ReportITSDefault.png");
}

void Chi2ReportTRD(){
  TString drawExpressionTRDChi2="";
  drawExpressionTRDChi2="[1,1,1]:";
  drawExpressionTRDChi2+="%Ogridx,gridy;hisChi2TRDTglQptMdEdx(0,100,0,10,20,80)(0,2)(f-mean p);:";
  drawExpressionTRDChi2+="%Ogridx,gridy;hisChi2TRDTglQptMdEdx(0,100,0,10,0,100)(0,1)(f-mean p);:";
  drawExpressionTRDChi2+="%Ogridx,gridy;hisChi2TRDTglQptMdEdx(0,100,0,10,0,100)(0,3)(f-mean p);:";
  TPad * padTRDChi2 = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionTRDChi2,keepArray, 1+2+4+8);
  padTRDChi2->SaveAs("Chi2ReportTRDDefault.png");
}


void MatchingReport(){ 
  TString drawExpressionRefit="";
  drawExpressionRefit="[2,2]:";
  drawExpressionRefit+="%Ogridx,gridy;hisITSRefitTRDOn(0,5,0,100,0,100)(0,1)(f-mean p);hisITSRefitAll(0,5,0,100,0,100)(0,1)(f-mean p);:";
  drawExpressionRefit+="%Ogridx,gridy;hisITSRefitTRDOn(0,5,0,100,0,100)(0,2)(f-mean p);hisITSRefitAll(0,5,0,100,0,100)(0,2)(f-mean p);:";
  drawExpressionRefit+="%Ogridx,gridy;hisTPCRefitTRDOn(0,5,0,100,0,100)(0,1)(f-mean p);hisTPCRefitAll(0,5,0,100,0,100)(0,1)(f-mean p);:";
  drawExpressionRefit+="%Ogridx,gridy;hisTPCRefitTRDOn(0,5,0,100,0,100)(0,2)(f-mean p);hisTPCRefitAll(0,5,0,100,0,100)(0,2)(f-mean p);:";
  TPad * padRefit = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionRefit,keepArray, 1+2+4+8);
  padRefit->SaveAs("refitDefault.png");
}

void K0Report(){ 
  TString drawExpressionK0="";
  drawExpressionK0="[4,4]:";
  drawExpressionK0+="%Ogridx,gridy;hisdMassK0(30,70,2,2,0,100)(0)(err);:";
  drawExpressionK0+="%Ogridx,gridy;hisdMassK0(30,70,3,3,0,100)(0)(err);:";
  drawExpressionK0+="%Ogridx,gridy;hisdMassK0(30,70,4,4,0,100)(0)(err);:";
  drawExpressionK0+="%Ogridx,gridy;hisdMassK0(30,70,5,5,0,100)(0)(err);:";
  drawExpressionK0+="%Ogridx,gridy;hisdMassK0(30,70,6,6,0,100)(0)(err);:";
  drawExpressionK0+="%Ogridx,gridy;hisdMassK0(30,70,7,7,0,100)(0)(err);:";
  drawExpressionK0+="%Ogridx,gridy;hisdMassK0(30,70,8,8,0,100)(0)(err);:";
  drawExpressionK0+="%Ogridx,gridy;hisdMassK0(30,70,9,9,0,100)(0)(err);:";
  TPad * padK0 = AliTreePlayer::DrawHistograms(0,hisV0Array,drawExpressionK0,keepArray, 1+2+4+8);
 
}
  
