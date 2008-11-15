/*
  .x ~/UliStyle.C
  .x ~/rootlogon.C
  gSystem->Load("libSTAT.so");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  gSystem->Load("libSTAT.so");

  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");  
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
  AliXRDPROOFtoolkit tool; 
  TChain * chainCosmic = tool.MakeChain("cosmic.txt","Track0",0,1000000);
  chainCosmic->Lookup();
  //chainCosmic->SetProof(kTRUE);


*/


void InitCuts(){
  //
  // Parameters diff
  //
  TCut cutP0("cutP0","abs(Tr0.fP[0]+Tr1.fP[0])<10");
  TCut cutP1("cutP1","abs(Tr0.fP[1]-Tr1.fP[1])<15"); 
  TCut cutP2("cutP2","abs(Tr1.fAlpha-Tr0.fAlpha+pi)<0.1")
  TCut cutP3("cutP3","abs(Tr0.fP[3]+Tr1.fP[3])<0.1");     
  TCut cutP4("cutP4","abs(Tr0.fP[4]+Tr1.fP[4])<0.5");     
  TCut cutP   = cutP0+cutP1+cutP2+cutP3+cutP4;
  //
  // Pull dif
  //
  TCut cutPull0("cutPull0","abs(Tr0.fP[0]+Tr1.fP[0])/sqrt(Tr0.fC[0]+Tr1.fC[0])<10");
  TCut cutPull2("((Tr1.fAlpha-Tr0.fAlpha+pi))/sqrt(Tr0.fC[5]+Tr1.fC[5])<10");
  TCut cutPull4("cutPull0","abs(Tr0.fP[4]+Tr1.fP[4])/sqrt(Tr0.fC[14]+Tr1.fC[14])<10");
  TCut cutPull = cutPull0+cutPull2+cutPull4;
  //
  // Geometrical cut 
  //
  TCut cutOx("Op1.fX>240&&Op0.fX>240");
  TCut cutOz("abs(Op1.fP[1])<240&&abs(Op0.fP[1])<240");
  TCut cutIx("Ip1.fX<90&&Ip0.fX<90");
  TCut cutIz("abs(Ip1.fP[1])<240&&abs(Ip0.fP[1])<240");
  TCut cutX00("abs(x00)<70");
  TCut cutX10("abs(x10)<70");
  TCut cutG = cutOx+cutOz+cutIx+cutIz+cutX00+cutX10;
  //
  //
  //
  TCut cutN("cutN","min(Orig0.fTPCncls,Orig1.fTPCncls)>80");
  TCut cutN120("cutN120","min(Orig0.fTPCncls,Orig1.fTPCncls)>120");
  TCut cutAll = cutP+cutPull+cutG+cutN+"abs(mag)>4";
  //
  TCut cuthpt("abs(Tr0.fP[4])+abs(Tr1.fP[4])<0.5");  
  TCut cutS("cutS","!crossI&&!crossO");
  TCut cutRun("run<620600");
  //
  chainCosmic->Draw(">>listELP",cutAll,"entryList");
  TEntryList *elist = (TEntryList*)gDirectory->Get("listELP");
  chainCosmic->SetEntryList(elist);
  //
  chainCosmic->Draw(">>listELFit",cutAll+cuthpt+cutS+cutRun,"entryList");
  TEntryList *elistFit = (TEntryList*)gDirectory->Get("listELFit");
  chainCosmic->SetEntryList(elistFit);
  //
  //
}

void SetAlias(){

  chainCosmic->SetAlias("dP0","(Tr0.fP[0]+Tr1.fP[0])");
  chainCosmic->SetAlias("dP1","(Tr0.fP[1]-Tr1.fP[1])");
  chainCosmic->SetAlias("dP2","(Tr1.fAlpha-Tr0.fAlpha+pi)");
  chainCosmic->SetAlias("dP3","(Tr0.fP[3]+Tr1.fP[3])");
  chainCosmic->SetAlias("dP4","(Tr0.fP[4]+Tr1.fP[4])");
  //
  chainCosmic->SetAlias("sP0","sqrt(Tr0.fC[0]+Tr1.fC[0])");
  chainCosmic->SetAlias("sP1","sqrt(Tr0.fC[2]+Tr1.fC[2])");
  chainCosmic->SetAlias("sP2","sqrt(Tr0.fC[5]+Tr0.fC[5])");
  chainCosmic->SetAlias("sP3","sqrt(Tr0.fC[9]+Tr1.fC[9])");
  chainCosmic->SetAlias("sP4","sqrt(Tr0.fC[14]+Tr1.fC[14])");
  //
  chainCosmic->SetAlias("dR","(1-abs(Tr0.fP[1]/250))");
  chainCosmic->SetAlias("side","(-1+(Tr0.fP[1]>0)*2)");
  chainCosmic->SetAlias("meanPt","((Tr0.Pt()+Tr1.Pt())/2.)");

}

void Correction(){
  TStatToolkit toolkit;
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD fitParam;
  TMatrixD covMatrix;
  //
  TString fstring="";
  fstring+="side++";
  //
  fstring+="dR++";
  fstring+="dR*dR++";
  fstring+="Tr0.fP[3]++";
  //
  fstring+="dR*side++";
  fstring+="dR*dR*side++";
  fstring+="Tr0.fP[3]*side++";
  //
  TString * strP0 = TStatToolkit::FitPlane(chainCosmic,"dP0", fstring.Data(),"1", chi2,npoints,fitParam,covMatrix);
  chainCosmic->SetAlias("corrP0",strP0->Data());

  TString * strP1 = TStatToolkit::FitPlane(chainCosmic,"dP1", fstring.Data(),"!crossI&&!crossO", chi2,npoints,fitParam,covMatrix);
  chainCosmic->SetAlias("corrP1",strP1->Data());

  TString * strP2 = TStatToolkit::FitPlane(chainCosmic,"dP2", fstring.Data(),"1", chi2,npoints,fitParam,covMatrix);
  chainCosmic->SetAlias("corrP2",strP2->Data());

  TString * strP3 = TStatToolkit::FitPlane(chainCosmic,"dP3", fstring.Data(),"!crossI&&!crossO", chi2,npoints,fitParam,covMatrix);
  chainCosmic->SetAlias("corrP3",strP3->Data());

  TString * strP4 = TStatToolkit::FitPlane(chainCosmic,"dP4", fstring.Data(),"1", chi2,npoints,fitParam,covMatrix);
  chainCosmic->SetAlias("corrP4",strP4->Data());

}



void DrawNClusterGeom(){
  //
  //
  //
  TH2F * hNd("hNd","hNd",50,0,190,100,10,160);
  chainCosmic->Draw("Orig0.fTPCncls:sqrt(x00^2+x01^2)>>hNd",cutOx+cutOz+cutIx+cutIz);
  hNd->FitSlicesY();
  hNd_1->SetXTitle("DCA_{r} (cm)");
  hNd_1->SetYTitle("Mean number of clusters");
  gPad->SaveAs("pic/NCl_Radius.eps");
  gPad->SaveAs("pic/NCl_Radius.gif");


}


void PtResolPt(){

  TH2F * hdPtPt = new TH2F("hdPtPt","hdPtPt",20,0.5,50,100,-0.6,0.6);
  TH2F * hdPtPtNoCor = new TH2F("hdPtPtNoCor","hdPtPtNoCor",20,0.5,50,100,-0.6,0.6);
  TH2F * hdPtPtCor = new TH2F("hdPtPtCor","hdPtPtCor",20,0.5,50,100,-0.6,0.6);
  //
  AliTPCcalibV0::BinLogX(&hdPtPt);
  AliTPCcalibV0::BinLogX(&hdPtPtNoCor);
  AliTPCcalibV0::BinLogX(&hdPtPtCor);

  chainCosmic->Draw("((Tr0.Pt()-Tr1.Pt())/meanPt)/sqrt(2.):meanPt>>hdPtPt",""+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("((dP4-corrP4)*meanPt)/sqrt(2.):meanPt>>hdPtPtCorr",""+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP4*meanPt)/sqrt(2.):meanPt>>hdPtPtNoCorr",""+cutRun+cutS+cutN120,"");

  hdPtPt->FitSlicesY();
  hdPtPt_2->SetXTitle("p_{t} (GeV)");
  hdPtPt_2->SetYTitle("#Deltap_{t}/p_{t} ");
  hdPtPt_2->Draw();


  TH2F * hdP4Pt = new TH2F("hdP4Pt","hdP4Pt",5,1,30,50,-0.05,0.05);
  chainCosmic->Draw("(dP4-corrP4)/sqrt(2.):Tr0.Pt()>>hdP4Pt",""+cutRun+cutS+cutN120,"");
  hdP4Pt->FitSlicesY();
  hdP4Pt_2->SetXTitle("p_{t} (GeV)");
  hdP4Pt_2->SetYTitle("#sigma 1/p_{t} (1/GeV)");
  hdP4Pt_2->Draw();
}



void P0resolZ(){
  //
  //
  //
  TH2F * hdP0Z = new TH2F("hdP0Z","hdP0Z",10,-250,250,100,-3.05,3.05);
  TH2F * hdP0ZNoCor = new TH2F("hdP0ZNoCor","hdP0ZNoCor",10,-250,250,100,-3.05,3.05);
  chainCosmic->Draw("(dP0-corrP0)/sqrt(2.):x02>>hdP0Z",""+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP0-0)/sqrt(2.):x02>>hdP0ZNoCor",""+cutRun+cutS+cutN120,"");
  hdP0Z->FitSlicesY();
  hdP0ZNoCor->FitSlicesY();
  hdP0Z_2->SetMinimum(0);
  hdP0Z_2->SetXTitle("Z position (cm)");
  hdP0Z_2->SetYTitle("#sigma_{y} (cm)");
  hdP0Z_2->SetMarkerStyle(25);  
  hdP0ZNoCor_2->SetMarkerStyle(26);
  hdP0Z_2->Draw();  
  hdP0ZNoCor_2->Draw("same");
  gPad->SaveAs("pic/SigmaY_z.gif");
  gPad->SaveAs("pic/SigmaY_z.eps");
  //
  TH2F * hdPP0Z = new TH2F("hdPP0Z","hdPP0Z",8,-200,200,50,-5.05,5.05);
  TH2F * hncdPP0Z = new TH2F("hncdPP0Z","hncdPP0Z",8,-200,200,50,-5.05,5.05);
  chainCosmic->Draw("(dP0-corrP0)/sP0:x02>>hdPP0Z",""+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP0-0)/sP0:x02>>hncdPP0Z",""+cutRun+cutS+cutN120,"");
  hdPP0Z->FitSlicesY();
  hncdPP0Z->FitSlicesY();
  hdPP0Z_2->SetMarkerStyle(25);
  hncdPP0Z_2->SetMarkerStyle(26);
  hdPP0Z_2->SetMinimum(0);
  hdPP0Z_2->SetXTitle("Z position (cm)");
  hdPP0Z_2->SetYTitle("#sigma y - Pull ");
  hdPP0Z_2->Draw();
  hncdPP0Z_2->Draw("same");
  gPad->SaveAs("pic/PullY_z.gif");
  gPad->SaveAs("pic/PullY_z.eps");

}


void P0resol1Pt(){
  //
  // P0 - Y -DCA resolution as function of the 1/pt
  //
  TH2F * hdP01Pt = new TH2F("hdP01Pt","hdP01Pt",6,0,1.0,100,-1.05,1.05);
  TH2F * hdP01PtNoCor = new TH2F("hdP01PtNoCor","hdP01PtNoCor",6,0,1.0,100,-1.05,1.05);
  chainCosmic->Draw("(dP0-corrP0)/sqrt(2.):sqrt(1/meanPt)>>hdP01Pt","side>0"+cutRun+cutS,"");
  chainCosmic->Draw("(dP0-0)/sqrt(2.):sqrt(1/meanPt)>>hdP01PtNoCor","side>0"+cutRun+cutS,"");
  hdP01Pt->FitSlicesY();
  hdP01PtNoCor->FitSlicesY();
  hdP01Pt_2->SetMinimum(0);
  hdP01Pt_2->SetXTitle("#sqrt{1/p_{t} (1/GeV)}");
  hdP01Pt_2->SetYTitle("#sigma_{y} (cm)");
  hdP01Pt_2->SetMarkerStyle(25);  
  hdP01PtNoCor_2->SetMarkerStyle(26);
  hdP01Pt_2->Draw();  
  hdP01PtNoCor_2->Draw("same");
  gPad->SaveAs("pic/SigmaY_1pt.gif");
  gPad->SaveAs("pic/SigmaY_1pt.eps");

  //
  TH2F * hPullP01Pt = new TH2F("hhPullP01Pt","hhPullP01Pt",6,0,1,50,-5.05,5.05);
  TH2F * hncPullP01Pt = new TH2F("hhncPullP01Pt","hhncPullP01Pt",6,0,1,50,-5.05,5.05);
  chainCosmic->Draw("(dP0-corrP0)/sP0:sqrt(1/meanPt)>>hhPullP01Pt","side>0"+cutRun+cutS,"");
  chainCosmic->Draw("(dP0-0)/sP0:sqrt(1/meanPt)>>hhncPullP01Pt","side>0"+cutRun+cutS,"");
  hPullP01Pt->FitSlicesY();
  hncPullP01Pt->FitSlicesY();
  hhPullP01Pt_2->SetMarkerStyle(25);
  hhncPullP01Pt_2->SetMarkerStyle(26);
  hhPullP01Pt_2->SetMinimum(0);
  hhPullP01Pt_2->SetXTitle("#sqrt{1/p_{t} (1/GeV)}");
  hhPullP01Pt_2->SetYTitle("#sigma_{y} - Pull ");
  hhPullP01Pt_2->Draw();
  hhncPullP01Pt_2->Draw("same");
  gPad->SaveAs("pic/PullY_1pt.gif");
  gPad->SaveAs("pic/PullY_1pt.eps");
}


void P1resolZ(){
  //
  //
  //
  TH2F *hdP1Z = new TH2F("hdP1Z","hdP1Z",10,-250,250,100,-1.05,1.05);
  TH2F *hdP1ZNoCor=new TH2F("hdP1ZNoCor","hdP1ZNoCor",10,-250,250,100,-1.05,1.05);
  chainCosmic->Draw("(dP1-corrP1)/sqrt(2.):x02>>hdP1Z",""+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP1-0)/sqrt(2.):x02>>hdP1ZNoCor",""+cutRun+cutS+cutN120,"");
  hdP1Z->FitSlicesY();
  hdP1ZNoCor->FitSlicesY();
  hdP1Z_2->SetMinimum(0);
  hdP1Z_2->SetXTitle("Z position (cm)");
  hdP1Z_2->SetYTitle("#sigma_{z} (cm)");
  hdP1Z_2->SetMarkerStyle(25);  
  hdP1ZNoCor_2->SetMarkerStyle(26);
  hdP1Z_2->Draw();  
  hdP1ZNoCor_2->Draw("same");
  gPad->SaveAs("pic/SigmaZ_z.gif");
  gPad->SaveAs("pic/SigmaZ_z.eps");
  //
  TH2F * hdPP1Z = new TH2F("hdPP1Z","hdPP1Z",8,-200,200,50,-5.05,5.05);
  TH2F * hncdPP1Z = new TH2F("hncdPP1Z","hncdPP1Z",8,-200,200,50,-5.05,5.05);
  chainCosmic->Draw("(dP1-corrP1)/sP1:x02>>hdPP1Z",""+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP1-0)/sP1:x02>>hncdPP1Z",""+cutRun+cutS+cutN120,"");
  hdPP1Z->FitSlicesY();
  hncdPP1Z->FitSlicesY();
  hdPP1Z_2->SetMarkerStyle(25);
  hncdPP1Z_2->SetMarkerStyle(26);
  hdPP1Z_2->SetMinimum(0);
  hdPP1Z_2->SetXTitle("Z position (cm)");
  hdPP1Z_2->SetYTitle("#sigma_{z} (Unit) ");
  hdPP1Z_2->Draw();
  hncdPP1Z_2->Draw("same");
  gPad->SaveAs("pic/PullZ_z.gif");
  gPad->SaveAs("pic/PullZ_z.eps");

}


void P1resol1Pt(){
  //
  // P1 - Y -DCA resolution as function of the 1/pt
  //
  TH2F * hdP11Pt = new TH2F("hdP11Pt","hdP11Pt",6,0,1.0,100,-1.05,1.05);
  TH2F * hdP11PtNoCor = new TH2F("hdP11PtNoCor","hdP11PtNoCor",6,0,1.0,100,-1.05,1.05);
  chainCosmic->Draw("(dP1-corrP1)/sqrt(2.):sqrt(1/meanPt)>>hdP11Pt","side>0"+cutRun+cutS,"");
  chainCosmic->Draw("(dP1-0)/sqrt(2.):sqrt(1/meanPt)>>hdP11PtNoCor","side>0"+cutRun+cutS,"");
  hdP11Pt->FitSlicesY();
  hdP11PtNoCor->FitSlicesY();
  hdP11Pt_2->SetMinimum(0);
  hdP11Pt_2->SetXTitle("#sqrt{1/p_{t} (1/GeV)}");
  hdP11Pt_2->SetYTitle("#sigma_{z} (cm)");
  hdP11Pt_2->SetMarkerStyle(25);  
  hdP11PtNoCor_2->SetMarkerStyle(26);
  hdP11Pt_2->Draw();  
  hdP11PtNoCor_2->Draw("same");
  gPad->SaveAs("pic/SigmaZ_1pt.gif");
  gPad->SaveAs("pic/SigmaZ_1pt.eps");

  //
  TH2F * hPullP11Pt = new TH2F("hhPullP11Pt","hhPullP11Pt",6,0,1,50,-5.05,5.05);
  TH2F * hncPullP11Pt = new TH2F("hhncPullP11Pt","hhncPullP11Pt",6,0,1,50,-5.05,5.05);
  chainCosmic->Draw("(dP1-corrP1)/sP1:sqrt(1/meanPt)>>hhPullP11Pt","side>0"+cutRun+cutS,"");
  chainCosmic->Draw("(dP1-0)/sP1:sqrt(1/meanPt)>>hhncPullP11Pt","side>0"+cutRun+cutS,"");
  hPullP11Pt->FitSlicesY();
  hncPullP11Pt->FitSlicesY();
  hhPullP11Pt_2->SetMarkerStyle(25);
  hhncPullP11Pt_2->SetMarkerStyle(26);
  hhPullP11Pt_2->SetMinimum(0);
  hhPullP11Pt_2->SetXTitle("#sqrt{1/p_{t} (1/GeV)}");
  hhPullP11Pt_2->SetYTitle("#sigma_{z} - Pull ");
  hhPullP11Pt_2->Draw();
  hhncPullP11Pt_2->Draw("same");
  gPad->SaveAs("pic/PullZ_1pt.gif");
  gPad->SaveAs("pic/PullZ_1pt.eps");
}




void P4resolZ(){
  //
  //
  //
  TH2F *hdP4Z = new TH2F("hdP4Z","hdP4Z",10,-250,250,100,-0.05,0.05);
  TH2F *hdP4ZNoCor=new TH2F("hdP4ZNoCor","hdP4ZNoCor",10,-250,250,100,-0.05,0.05);
  chainCosmic->Draw("(dP4-corrP4)/sqrt(2.):x02>>hdP4Z",""+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP4-0)/sqrt(2.):x02>>hdP4ZNoCor",""+cutRun+cutS+cutN120,"");
  hdP4Z->FitSlicesY();
  hdP4ZNoCor->FitSlicesY();
  hdP4Z_2->SetMinimum(0);
  hdP4Z_2->SetXTitle("Z position (cm)");
  hdP4Z_2->SetYTitle("#sigma_{1/pt} (1/GeV)");
  hdP4Z_2->SetMarkerStyle(25);  
  hdP4ZNoCor_2->SetMarkerStyle(26);
  hdP4Z_2->Draw();  
  hdP4ZNoCor_2->Draw("same");
  gPad->SaveAs("pic/Sigma1Pt_z.gif");
  gPad->SaveAs("pic/Sigma1Pt_z.eps");
  //
  TH2F * hdPP4Z = new TH2F("hdPP4Z","hdPP4Z",8,-200,200,50,-5.05,5.05);
  TH2F * hncdPP4Z = new TH2F("hncdPP4Z","hncdPP4Z",8,-200,200,50,-5.05,5.05);
  chainCosmic->Draw("(dP4-corrP4)/sP4:x02>>hdPP4Z",""+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP4-0)/sP4:x02>>hncdPP4Z",""+cutRun+cutS+cutN120,"");
  hdPP4Z->FitSlicesY();
  hncdPP4Z->FitSlicesY();
  hdPP4Z_2->SetMarkerStyle(25);
  hncdPP4Z_2->SetMarkerStyle(26);
  hdPP4Z_2->SetMinimum(0);
  hdPP4Z_2->SetXTitle("Z position (cm)");
  hdPP4Z_2->SetYTitle("#sigma_{1/pt} - Pull ");
  hdPP4Z_2->Draw();
  hncdPP4Z_2->Draw("same");
  gPad->SaveAs("pic/Pull1Pt_z.gif");
  gPad->SaveAs("pic/Pull1Pt_z.eps");

}


void P4resol1Pt(){
  //
  // P4 - 1/Pt resolution as function of the 1/pt
  //
  TH2F * hdP41Pt = new TH2F("hdP41Pt","hdP41Pt",6,0,1.0,100,-0.05,0.05);
  TH2F * hdP41PtNoCor = new TH2F("hdP41PtNoCor","hdP41PtNoCor",6,0,1.0,100,-0.05,0.05);
  chainCosmic->Draw("(dP4-corrP4)/sqrt(2.):sqrt(1/meanPt)>>hdP41Pt","side>0"+cutRun+cutS+cutN120,"");
  chainCosmic->Draw("(dP4-0)/sqrt(2.):sqrt(1/meanPt)>>hdP41PtNoCor","side>0"+cutRun+cutS+cutN120,"");
  hdP41Pt->FitSlicesY();
  hdP41PtNoCor->FitSlicesY();
  hdP41Pt_2->SetMinimum(0);
  hdP41Pt_2->SetXTitle("#sqrt{1/p_{t} (1/GeV)}");
  hdP41Pt_2->SetYTitle("#sigma_{1/pt} (1/GeV)");
  hdP41Pt_2->SetMarkerStyle(25);  
  hdP41PtNoCor_2->SetMarkerStyle(26);
  hdP41Pt_2->Draw();  
  hdP41PtNoCor_2->Draw("same");
  gPad->SaveAs("pic/Sigma1Pt_1pt.gif");
  gPad->SaveAs("pic/Sigma1Pt_1pt.eps");

  //
  TH2F * hPullP41Pt = new TH2F("hhPullP41Pt","hhPullP41Pt",6,0,1,50,-5.05,5.05);
  TH2F * hncPullP41Pt = new TH2F("hhncPullP41Pt","hhncPullP41Pt",6,0,1,50,-5.05,5.05);
  chainCosmic->Draw("(dP4-corrP4)/sP4:sqrt(1/meanPt)>>hhPullP41Pt","side>0"+cutRun+cutS,"");
  chainCosmic->Draw("(dP4-0)/sP4:sqrt(1/meanPt)>>hhncPullP41Pt","side>0"+cutRun+cutS,"");
  hPullP41Pt->FitSlicesY();
  hncPullP41Pt->FitSlicesY();
  hhPullP41Pt_2->SetMarkerStyle(25);
  hhncPullP41Pt_2->SetMarkerStyle(26);
  hhPullP41Pt_2->SetMinimum(0);
  hhPullP41Pt_2->SetXTitle("#sqrt{1/p_{t}} (1/GeV)}");
  hhPullP41Pt_2->SetYTitle("#sigma_{1/pt} (Unit) ");
  hhPullP41Pt_2->Draw();
  hhncPullP41Pt_2->Draw("same");
  gPad->SaveAs("pic/Pull1Pt_1pt.gif");
  gPad->SaveAs("pic/Pull1Pt_1pt.eps");
}












void P2resol(){
  //
  //
  //
  TH2F * hdP2Z = new TH2F("hdP2Z","hdP2Z",10,-250,250,50,-0.02,0.02);
  chainCosmic->Draw("(dP2-corrP2)/sqrt(2.):x02>>hdP2Z",""+cutRun+cutS+cutN120,"");
  hdP2Z->FitSlicesY();
  hdP2Z_2->SetXTitle("Z position (cm)");
  hdP2Z_2->SetYTitle("#sigma phi (rad)");
  hdP2Z_2->Draw();
  //
  TH2F * hdPP2Z = new TH2F("hdPP2Z","hdPP2Z",8,-200,200,50,-5.05,5.05);
  chainCosmic->Draw("(dP2-corrP2)/sP2:x02>>hdPP2Z",""+cutRun+cutS+cutN120,"");
  hdPP2Z->FitSlicesY();
  hdPP2Z_2->SetXTitle("Z position (cm)");
  hdPP2Z_2->SetYTitle("#sigma phi - Pull ");
  hdPP2Z_2->Draw();
}

void P3resol(){
  //
  //
  //
  TH2F * hdP3Z("hdP3Z","hdP3Z",10,-250,250,50,-0.005,0.005);
  chainCosmic->Draw("(dP3-corrP3)/sqrt(2.):x02>>hdP3Z",""+cutRun+cutS+cutN120,"");
  hdP3Z->FitSlicesY();
  hdP3Z_2->SetXTitle("Z position (cm)");
  hdP3Z_2->SetYTitle("#sigma  (rad)");
  hdP3Z_2->Draw();
  //
  TH2F * hdPP3Z("hdPP3Z","hdPP3Z",8,-200,200,50,-5.05,5.05);
  chainCosmic->Draw("(dP3-corrP3)/sP3:x02>>hdPP3Z",""+cutRun+cutS+cutN120,"");
  hdPP3Z->FitSlicesY();
  hdPP3Z_2->SetXTitle("Z position (cm)");
  hdPP3Z_2->SetYTitle("#sigma phi - Pull ");
  hdPP3Z_2->Draw();
}







void PtResolN(){
  //
  //
  //
  TH2F * hdP4Ncl= new TH2F("hdp4Ncl","hdp4Ncl",5,80,160,100,-0.1,0.1);
  chainCosmic->Draw("(dP4-corrP4)/sqrt(2.):min(Orig0.fTPCncls,Orig1.fTPCncls)>>hdp4Ncl","side>0"+cuthpt+cutRun+cutS,"");
  hdp4Ncl->FitSlicesY();
  hdp4Ncl_2->SetXTitle("Number of clusters");
  hdp4Ncl_2->SetYTitle("#sigma 1/p_{t} (1/GeV)");
  hdp4Ncl_2->Draw();
  gPad->SaveAs("pic/Sigma1Pt_N.gif");
  gPad->SaveAs("pic/Sigma1Pt_N.eps");

  //
  //
  TH2F * hdP4PullNcl = new TH2F("hdP4PullNcl","hdP4PullNcl",5,80,160,100,-6.1,6.1);
  chainCosmic->Draw("(Tr1.fP[4]+Tr0.fP[4])/sqrt(Tr1.fC[14]+Tr0.fC[14]):min(Orig0.fTPCncls,Orig1.fTPCncls)>>hdP4PullNcl","side>0"+cuthpt+cutRun+cutS,"");
  hdP4PullNcl->FitSlicesY();
  hdP4PullNcl_2->SetXTitle("Number of clusters");
  hdP4PullNcl_2->SetYTitle("#sigma 1/p_{t} (Unit)");
  hdP4PullNcl_2->Draw();
  gPad->SaveAs("pic/Pull1Pt_N.gif");
  gPad->SaveAs("pic/Pull1Pt_N.eps");

}





