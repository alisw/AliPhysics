/*
  wget http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/alice/data/2015/LHC15o/pass1/LHC15o.pass1.B1.Bin0/performanceHisto.root -o performanceHisto.root

   AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");


*/

TObjArray * hisArray=NULL;

void ReadInputData() {
  // read cache file from server
  TFile *finput=NULL;
  TString fname = "./performanceHisto.root";
  if (!gSystem->AccessPathName( fname )) {
    finput = TFile::Open( fname ); // check if file in local directory exists
  }
  else {
    TFile::SetCacheFileDir(".");
    finput = TFile::Open("http://aliqatrkeos.web.cern.ch/aliqatrkeos/performance/alice/data/2015/LHC15o/pass1/LHC15o.pass1.B1.Bin0/performanceHisto.root", "CACHEREAD"); // if not: download from ROOT server
  }
  // get histogram array
  hisArray = new TObjArray();
  TList *keys = finput->GetListOfKeys();
  for (Int_t iKey = 0; iKey < keys->GetEntries(); iKey++) {
    TObject *o = finput->Get(
            TString::Format("%s;%d", keys->At(iKey)->GetName(), ((TKey *) keys->At(iKey))->GetCycle()).Data());
    hisArray->AddLast(o);
  }
}

/// Example P4 report
void makeP4Report(){
  TString drawExpressionP4="";
  drawExpressionP4="[1,1,2]:";
  drawExpressionP4+="%Ogridx,gridy;hisDeltaP4_Allv_qPt_tgl(0,400,80,120,0,10)(0,1)(f-grms p);:";
  drawExpressionP4+="%Ogridx,gridy;hisPullP4_Allv_qPt_tgl(0,400,80,120,0,10)(0,1)(f-grms p);:";
  drawExpressionP4+="%Ogridx,gridy;hisDeltaP4_Allv_qPt_tgl(100,300,99,102,0,10)(0)(err);:";
  drawExpressionP4+="%Ogridx,gridy;hisPullP4_Allv_qPt_tgl(100,300,99,102,0,10)(0)(err);:";
  TPad * padP4 = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionP4,keepArray, 1+2+4+8);
  padP4->SaveAs("deltaP4TPCtoFull.png");
  padP4->SaveAs("deltaP4TPCtoFull.C");
}

/// Example P4 report
void makeP4ReportHisto(){
  TString drawExpressionP4="";
  drawExpressionP4="[1,1]:";
  drawExpressionP4+="%Ogridx,gridy;hisDeltaP4_Allv_qPt_tgl(100,300,0,40,102,0,10)(0)(err);hisDeltaP4_Allv_qPt_tgl(100,300,40,80,102,0,10)(0)(err);hisDeltaP4_Allv_qPt_tgl(100,300,80,10,102,0,10)(0)(err);:";
  drawExpressionP4+="%Ogridx,gridy;hisPullP4_Allv_qPt_tgl(100,300,99,102,0,10)(0)(err);:";
  TPad * padP4 = AliTreePlayer::DrawHistograms(0,hisArray,drawExpressionP4,keepArray, 1+2+4+8);
  padP4->SaveAs("deltaP4TPCtoFull.png");
  padP4->SaveAs("deltaP4TPCtoFull.C");
}


