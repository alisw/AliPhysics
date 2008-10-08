/*

.x ~/rootlogon.C

gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
AliXRDPROOFtoolkit tool;
gSystem->Load("libSTAT.so");


//  Load - Calibration using the CE  
.L $ALICE_ROOT/TPC/CalibMacros/AnalyzeLaser.C+
LoadViewer();
tree->SetAlias("dedgePI","tan(10*pi/180.)*lx.fElements-ly.fElements+0.5*0.4");
tree->SetAlias("dedgeMI","ly.fElements+tan(10*pi/180.)*lx.fElements-0.5*0.4");
tree->SetAlias("dedgePO","tan(10*pi/180.)*lx.fElements-ly.fElements+0.5*0.6");
tree->SetAlias("dedgeMO","ly.fElements+tan(10*pi/180.)*lx.fElements-0.5*0.6");

// Load -Calibration using laser tracks
// PROOF neccessary


TChain * chainFit=0;
TChain * chainTrack=0;
TChain * chain=0;
//
gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
AliXRDPROOFtoolkit tool;
chainTrack = tool.MakeChain("laser.txt","Track",0,10200);
chainTrack->Lookup();
chainTrack->SetProof(kTRUE);
chain = tool.MakeChain("laser.txt","Residuals",0,10200);
chain->Lookup();
chainFit = tool.MakeChain("laser.txt","FitModels",0,10200);
chainFit->Lookup();
chainFit->SetProof(kTRUE);
chain->SetProof(kTRUE);
chainDrift = tool.MakeChain("laser.txt","driftv",0,10200);
chainDrift->Lookup();
chainDrift->SetProof(kTRUE);

 TCut cutChi2YOut("sqrt(chi2y2Out*dEdx)<5");
  TCut cutChi2ZOut("sqrt(chi2z2Out*dEdx)<5");
  TCut cutChi2YIn("sqrt(chi2y2In*dEdx)<5");
  TCut cutChi2ZIn("sqrt(chi2z2In*dEdx)<5");
  //
  TCut cutdEdx("sqrt(dEdx)<30&&sqrt(dEdx)>3");
  TCut cutDY("abs(yPol2In.fElements[2]*nclO*nclO/4.)<3");
  TCut cutN("nclO>20&&nclI>20");
  TCut cutA = cutChi2YOut+cutChi2ZOut+cutChi2YIn+cutChi2ZIn+cutN+cutdEdx;
  //
  // Cluster cuts
  //
  TCut cutClY("Cl[].fY!=0&&abs(Cl[].fY-TrYpol2.fElements)<1.");
  TCut cutClZ("abs(Cl[].fZ-TrZpol2.fElements)<1.0");
  TCut cutClX("abs(Cl[].fX)>10");
  TCut cutE("abs(Cl[].fY/Cl[].fX)<0.14");
  TCut cutSY("sqrt(Cl[].fSigmaY2)>0.05");
  TCut cutSZ("sqrt(Cl[].fSigmaZ2)>0.05");
  TCut cutQ("sqrt(Cl[].fMax)>4");
  TCut cutCl=cutClY+cutClZ+cutClX+cutE+cutSY+cutSZ+cutQ;



*/




void CalibEdgeQPad(){
  //
  // Get the Mean charge on pad as function of pad number
  //
  TH2* hisQPadPShort= 0;
  TH2* hisQPadMMiddle= 0;
  TH2* hisQPadMLong= 0;
  //
  tree->Draw("qIn.fElements/qF2.fElements:pad.fElements>>hisQPadPShort(10,0,10,50,0.1,2)","abs(qIn.fElements/qF2.fElements-1)<0.9&&sector%36>17&&sector<36","");
  hisQPadPShort = (TH2*)gROOT->FindObject("hisQPadPShort");
  hisQPadPShort->FitSlicesY();
  TH1 * hisQPadPShortM =  (TH1*)(gROOT->FindObject("hisQPadPShort_1")->Clone());
  //
  tree->Draw("qIn.fElements/qF2.fElements:pad.fElements>>hisQPadPMiddle(10,0,10,50,0.1,2)","abs(qIn.fElements/qF2.fElements-1)<0.9&&sector%36>17&&sector>36&&lx.fElements<197","");
  hisQPadPMiddle = (TH2*)gROOT->FindObject("hisQPadPMiddle");
  hisQPadPMiddle->FitSlicesY();
  TH1 * hisQPadPMiddleM =  (TH1*)(gROOT->FindObject("hisQPadPMiddle_1")->Clone());
  //
  tree->Draw("qIn.fElements/qF2.fElements:pad.fElements>>hisQPadPLong(10,0,10,50,0.1,2)","abs(qIn.fElements/qF2.fElements-1)<0.9&&sector%36>17&&sector>36&&lx.fElements>197","");
  hisQPadPLong = (TH2*)gROOT->FindObject("hisQPadPLong");
  hisQPadPLong->FitSlicesY();
  TH1 * hisQPadPLongM =  (TH1*)(gROOT->FindObject("hisQPadPLong_1")->Clone());
  //
  hisQPadPShortM->SetLineColor(2);
  hisQPadPMiddleM->SetLineColor(4);
  hisQPadPLongM->SetLineColor(6);
  hisQPadPShortM->SetMarkerStyle(23);
  hisQPadPMiddleM->SetMarkerStyle(24);
  hisQPadPLongM->SetMarkerStyle(25);

  hisQPadPShortM->Draw();
  hisQPadPShortM->SetMinimum(0);
  hisQPadPShortM->SetMaximum(1.4);
  hisQPadPMiddleM->Draw("same");
  hisQPadPLongM->Draw("same");
  TLegend *legend = new TLegend(0.45,0.12,0.85,0.55, "Mean charge vs Pad number");
  legend->SetBorderSize(1);
  legend->AddEntry(hisQPadPShortM,"Short pads");
  legend->AddEntry(hisQPadPMiddleM,"Medium pads");
  legend->AddEntry(hisQPadPLongM,"Long pads");
  legend->Draw();

}

void CalibEdgeTPad(){
  //
  // Get the Mean charge on pad as function of pad number
  //
  TH2* hisTPadPShort= 0;
  TH2* hisTPadMMiddle= 0;
  TH2* hisTPadMLong= 0;
  //
  tree->Draw("(timeIn.fElements-timeF2.fElements)*2.64:pad.fElements>>hisTPadPShort(10,0,10,50,-1,1)","timeIn.fElements!=0&&sector%36>17&&sector<36","");
  hisTPadPShort = (TH2*)gROOT->FindObject("hisTPadPShort");
  hisTPadPShort->FitSlicesY();
  TH1 * hisTPadPShortM =  (TH1*)(gROOT->FindObject("hisTPadPShort_1")->Clone());
  //
  tree->Draw("(timeIn.fElements-timeF2.fElements)*2.64:pad.fElements>>hisTPadPMiddle(10,0,10,50,-1,1)","timeIn.fElements!=0&&sector%36>17&&sector>36&&lx.fElements<197","");
  hisTPadPMiddle = (TH2*)gROOT->FindObject("hisTPadPMiddle");
  hisTPadPMiddle->FitSlicesY();
  TH1 * hisTPadPMiddleM =  (TH1*)(gROOT->FindObject("hisTPadPMiddle_1")->Clone());
  //
  tree->Draw("(timeIn.fElements-timeF2.fElements)*2.64:pad.fElements>>hisTPadPLong(10,0,10,50,-1,1)","timeIn.fElements!=0&&sector%36>17&&sector>36&&lx.fElements>197","");
  hisTPadPLong = (TH2*)gROOT->FindObject("hisTPadPLong");
  hisTPadPLong->FitSlicesY();
  TH1 * hisTPadPLongM =  (TH1*)(gROOT->FindObject("hisTPadPLong_1")->Clone());
  //



  //
  hisTPadPShortM->SetLineColor(2);
  hisTPadPMiddleM->SetLineColor(4);
  hisTPadPLongM->SetLineColor(6);
  hisTPadPShortM->SetMarkerStyle(23);
  hisTPadPMiddleM->SetMarkerStyle(24);
  hisTPadPLongM->SetMarkerStyle(25);

  hisTPadPShortM->Draw();
  hisTPadPShortM->SetMinimum(-0.3);
  hisTPadPShortM->SetMaximum(0.2);
  hisTPadPMiddleM->Draw("same");
  hisTPadPLongM->Draw("same");
  TLegend *legend = new TLegend(0.45,0.12,0.85,0.35, "Mean charge vs Pad number");
  legend->SetBorderSize(1);
  legend->AddEntry(hisTPadPShortM,"Short pads");
  legend->AddEntry(hisTPadPMiddleM,"Medium pads");
  legend->AddEntry(hisTPadPLongM,"Long pads");
  legend->Draw();

}




void CalibEdgeQ(){
  //
  // 
  //
  TH1* hisQedgePInner= 0;
  TH1* hisQedgeMInner= 0;
  TH1* hisQedgePOuter= 0;
  TH1* hisQedgeMOuter= 0;

  tree->Draw("qIn.fElements/qF2.fElements:dedgePI>>hisQedgePInner(50,0,5)","abs(qIn.fElements/qF2.fElements-1)<0.9&&sector%36>17&&sector<36","prof");
  tree->Draw("qIn.fElements/qF2.fElements:dedgeMI>>hisQedgeMInner(50,0,5)","abs(qIn.fElements/qF2.fElements-1)<0.9&&sector%36>17&&sector<36","prof");
  
  hisQedgePInner = (TH1*)gROOT->FindObject("hisQedgePInner");
  hisQedgeMInner = (TH1*)gROOT->FindObject("hisQedgeMInner");
  hisQedgePInner->SetLineColor(2);
  hisQedgeMInner->SetLineColor(4);
  hisQedgePInner->Draw();
  hisQedgeMInner->Draw("same");
  //
  //
  //
  tree->Draw("qIn.fElements/qF2.fElements:dedgePO>>hisQedgePOuter(50,0,5)","abs(qIn.fElements/qF2.fElements-1)<0.9&&sector%36>17&&sector>36&&lx.fElements<200","prof");
  tree->Draw("qIn.fElements/qF2.fElements:dedgeMO>>hisQedgeMOuter(50,0,5)","abs(qIn.fElements/qF2.fElements-1)<0.9&&sector%36>17&&sector>36&&lx.fElements<200","prof");  
  hisQedgePOuter = (TH1*)gROOT->FindObject("hisQedgePOuter");
  hisQedgeMOuter = (TH1*)gROOT->FindObject("hisQedgeMOuter");
  hisQedgePOuter->SetLineColor(2);
  hisQedgeMOuter->SetLineColor(4);
  hisQedgePOuter->Draw();
  hisQedgeMOuter->Draw("same");
}






void EdgeLaserZ(){
  

  TH2* hisZedgePInnerA   = 0;
  TH2* hisZedgePMiddleA = 0;
  TH2* hisZedgePLongA   = 0;
  TObjArray farray(100);
  delete gROOT->FindObject("hisEdge");
  chain->Draw("Cl[].fZ-TrZpol1.fElements:Cl[].fPad>>hisEdge1(50,0,10,50, -0.5,0.5)",cutA+cutClZ+"Cl[].fDetector<36&&Cl[].fDetector<18&&Cl[].fPad>0","",50000);  
  hisZedgePInnerA = (TH2*)(gROOT->FindObject("hisEdge"));
  
  hisZedgePInnerA->FitSlicesY(0,0,-1,0,"QNR",&farray);
  farray->At(1)->Draw();

}
