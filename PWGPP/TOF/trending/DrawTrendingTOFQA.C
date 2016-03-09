/*
  fbellini@cern.ch - last update on 27/02/2014
  Macro to draw the TOF QA trending plots by accessing the std tree.
  To be mainly used with the automatic scripts to fill the QA repository.
  Launch with 
  aliroot -x -b -q "DrawTrendingTOFQA.C" 
  The macro produces one png file for each trending variables
  and a .root file with the histograms
*/
Int_t DrawTrendingTOFQA(TString mergedTrendFile = "trending.root", // trending tree file 
			Bool_t displayAll = kFALSE) //set to kTRUE to display trending for expert plots
{
  //
  //reads merged trending.root file and draws trending plots from tree
  //
  if (!mergedTrendFile) {
    Printf("Cannot open merged trend file with TOF QA");
    return 1;
  }
  
  char  outfilename[200]= "ProductionQA.hist.root";
  TString plotDir(".");
  // TString plotDir(Form("PlotsTrending"));
  // gSystem->Exec(Form("mkdir %s",plotDir.Data()));

  Int_t runNumber=0;
  Double_t avTime=0., peakTime=0., spreadTime=0., peakTimeErr=0., spreadTimeErr=0.,negTimeRatio=0.,
    avRawTime=0., peakRawTime=0., spreadRawTime=0., peakRawTimeErr=0., spreadRawTimeErr=0., 
    avTot=0., peakTot=0.,spreadTot=0.,  peakTotErr=0.,spreadTotErr=0.,
    orphansRatio=0., avL=0., negLratio=0.,
    effPt1=0., effPt2=0., matchEffLinFit1Gev=0.,matchEffLinFit1GevErr=0.;
  Double_t avDiffTime=0.,peakDiffTime=0., spreadDiffTime=0.,peakDiffTimeErr=0., spreadDiffTimeErr=0.,avT0fillRes=0.;
   
  Double_t avT0A=0.,peakT0A=0., spreadT0A=0.,peakT0AErr=0., spreadT0AErr=0.;
  Double_t avT0C=0.,peakT0C=0., spreadT0C=0.,peakT0CErr=0., spreadT0CErr=0.;
  Double_t avT0AC=0.,peakT0AC=0., spreadT0AC=0.,peakT0ACErr=0., spreadT0ACErr=0.;
  Double_t avT0res=0.,peakT0res=0., spreadT0res=0.,peakT0resErr=0., spreadT0resErr=0.;
  Float_t avMulti=0;
  Float_t fractionEventsWHits=0.0;
  Double_t goodChannelRatio=0.0;
  Double_t goodChannelRatioInAcc=0.0;

  
  TFile * fin = TFile::Open(mergedTrendFile.Data());
  TTree * ttree = (TTree*) fin->Get("trending");
  if (!ttree){
    Printf("Invalid trending tree.");
    return 2;
  }
  ttree->SetBranchAddress("run",&runNumber);
  ttree->SetBranchAddress("avMulti",&avMulti);
  ttree->SetBranchAddress("goodChannelsRatio",&goodChannelRatio);   
  ttree->SetBranchAddress("goodChannelsRatioInAcc",&goodChannelRatioInAcc);   
  ttree->SetBranchAddress("avTime",&avTime); //mean time
  ttree->SetBranchAddress("peakTime",&peakTime); //main peak time after fit
  ttree->SetBranchAddress("spreadTime",&spreadTime); //spread of main peak of time after fit
  ttree->SetBranchAddress("peakTimeErr",&peakTimeErr); //main peak time after fit error
  ttree->SetBranchAddress("spreadTimeErr",&spreadTimeErr); //spread of main peak of time after fit error
  ttree->SetBranchAddress("negTimeRatio",&negTimeRatio); //negative time ratio
  ttree->SetBranchAddress("avRawTime",&avRawTime); //mean raw time
  ttree->SetBranchAddress("peakRawTime",&peakRawTime); //mean peak of raw time after fit
  ttree->SetBranchAddress("spreadRawTime",&spreadRawTime); //spread of main peak of raw time after fit
  ttree->SetBranchAddress("peakRawTimeErr",&peakRawTimeErr); //main peak raw  time after fit error
  ttree->SetBranchAddress("spreadRawTimeErr",&spreadRawTimeErr); //spread of  raw main peak of time after fit error
  ttree->SetBranchAddress("avTot",&avTot); //main peak tot
  ttree->SetBranchAddress("peakTot",&peakTot); // main peak of tot after fit
  ttree->SetBranchAddress("spreadTot",&spreadTot); //spread of main peak of tot after fit
  ttree->SetBranchAddress("peakTotErr",&peakTotErr); // main peak of tot after fit
  ttree->SetBranchAddress("spreadTotErr",&spreadTotErr); //spread of main peak of tot after fit
  ttree->SetBranchAddress("orphansRatio",&orphansRatio); //orphans ratio
  ttree->SetBranchAddress("avL",&avL); //mean track length
  ttree->SetBranchAddress("negLratio",&negLratio);//ratio of tracks with track length <350 cm
  ttree->SetBranchAddress("effPt1",&effPt1);//matching eff at 1 GeV/c
  ttree->SetBranchAddress("effPt2",&effPt2); //matching eff at 2 GeV/c
  ttree->SetBranchAddress("matchEffLinFit1Gev",&matchEffLinFit1Gev);//matching eff fit param 
  ttree->SetBranchAddress("matchEffLinFit1GevErr",&matchEffLinFit1GevErr);////matching eff fit param error
  ttree->SetBranchAddress("avPiDiffTime",&avDiffTime); //mean t-texp
  ttree->SetBranchAddress("peakPiDiffTime",&peakDiffTime); //main peak t-texp after fit
  ttree->SetBranchAddress("spreadPiDiffTime",&spreadDiffTime); //spread of main peak t-texp after fit
  ttree->SetBranchAddress("peakPiDiffTimeErr",&peakDiffTimeErr); //main peak t-texp after fit error
  ttree->SetBranchAddress("spreadPiDiffTimeErr",&spreadDiffTimeErr); //spread of main peak of t-texp after fit error
  ttree->SetBranchAddress("avT0A",&avT0A); //main peak t0A
  ttree->SetBranchAddress("peakT0A",&peakT0A); // main peak of t0A after fit
  ttree->SetBranchAddress("spreadT0A",&spreadT0A); //spread of main peak of t0A after fit
  ttree->SetBranchAddress("peakT0AErr",&peakT0AErr); // main peak of t0A after fit
  ttree->SetBranchAddress("spreadT0AErr",&spreadT0AErr); //spread of main peak of t0A after fit
  ttree->SetBranchAddress("avT0C",&avT0C); //main peak t0C
  ttree->SetBranchAddress("peakT0C",&peakT0C); // main peak of t0C after fit
  ttree->SetBranchAddress("spreadT0C",&spreadT0C); //spread of main peak of t0C after fit
  ttree->SetBranchAddress("peakT0CErr",&peakT0CErr); // main peak of t0C after fit
  ttree->SetBranchAddress("spreadT0CErr",&spreadT0CErr); //spread of main peak of t0C after fit
  ttree->SetBranchAddress("avT0AC",&avT0AC); //main peak t0AC
  ttree->SetBranchAddress("peakT0AC",&peakT0AC); // main peak of t0AC after fit
  ttree->SetBranchAddress("spreadT0AC",&spreadT0AC); //spread of main peak of t0AC after fit
  ttree->SetBranchAddress("peakT0ACErr",&peakT0ACErr); // main peak of t0AC after fit
  ttree->SetBranchAddress("spreadT0ACErr",&spreadT0ACErr); //spread of main peak of t0AC after fit
  ttree->SetBranchAddress("avT0res",&avT0res); //main peak t0AC
  ttree->SetBranchAddress("peakT0res",&peakT0res); // main peak of t0AC after fit
  ttree->SetBranchAddress("spreadT0res",&spreadT0res); //spread of main peak of t0AC after fit
  ttree->SetBranchAddress("peakT0resErr",&peakT0resErr); // main peak of t0AC after fit
  ttree->SetBranchAddress("spreadT0resErr",&spreadT0resErr); //spread of main peak of t0AC after fit
  ttree->SetBranchAddress("avT0fillRes",&avT0fillRes); //t0 fill res

  Int_t nRuns=ttree->GetEntries();
  TList lista;
   
  TH1F * hAvMulti=new TH1F("hAvMulti","Average multiplicity of matched tracks <N_{TOF}>;; <N_{TOF}>", nRuns,0., nRuns);//, 600, 0. , 600.);
  hAvMulti->SetDrawOption("E1");
  hAvMulti->SetMarkerStyle(20);
  hAvMulti->SetMarkerColor(kBlue);
	
  TH1F * hAvDiffTimeVsRun=new TH1F("hAvDiffTimeVsRun","Mean t-t_{exp} (no fit);run;<t^{TOF}-t_{exp,#pi}> (ps)",nRuns,0., nRuns);//, 600, 0. , 600.);
  hAvDiffTimeVsRun->SetDrawOption("E1");
  hAvDiffTimeVsRun->SetMarkerStyle(20);
  hAvDiffTimeVsRun->SetMarkerColor(kBlue);
  //   hAvTimeVsRun->GetYaxis()->SetRangeUser(0.0, 50.0);

  TH1F * hPeakDiffTimeVsRun=new TH1F("hPeakDiffTimeVsRun","t-t_{exp} (gaussian fit) ;; <t^{TOF}-t_{exp,#pi}> (ps)",nRuns,0., nRuns);//,600, 0. , 600. );
  hPeakDiffTimeVsRun->SetDrawOption("E1");
  hPeakDiffTimeVsRun->SetMarkerStyle(20);
  hPeakDiffTimeVsRun->SetMarkerColor(kBlue);
   
  TH1F * hSpreadDiffTimeVsRun=new TH1F("hSpreadDiffTimeVsRun","#sigma(t-t_{exp}) (gaussian fit);; #sigma(t^{TOF}-t_{exp,#pi}) (ns)",nRuns,0., nRuns);//, 100, 0. , 100.);
  hSpreadDiffTimeVsRun->SetDrawOption("E1");
  hSpreadDiffTimeVsRun->SetMarkerStyle(20);
  hSpreadDiffTimeVsRun->SetMarkerColor(kBlue);

  TH1F * hAvTimeVsRun=new TH1F("hAvTimeVsRun","<t^{TOF}>;;<t^{TOF}> (ns)",nRuns,0., nRuns);//, 600, 0. , 600.);
  hAvTimeVsRun->SetDrawOption("E1");
  hAvTimeVsRun->SetMarkerStyle(20);
  hAvTimeVsRun->SetMarkerColor(kBlue);
  //   hAvTimeVsRun->GetYaxis()->SetRangeUser(0.0, 50.0);

  TH1F * hPeakTimeVsRun=new TH1F("hPeakTimeVsRun","Peak value of t^{TOF} (landau fit);;t_{peak}^{TOF} (ns)",nRuns,0., nRuns);//,600, 0. , 600. );
  hPeakTimeVsRun->SetDrawOption("E1");
  hPeakTimeVsRun->SetMarkerStyle(20);
  hPeakTimeVsRun->SetMarkerColor(kBlue);
   
  TH1F * hSpreadTimeVsRun=new TH1F("hSpreadTimeVsRun","Spread of t^{TOF} (landau fit);; #sigma(t^{TOF}) (ns)",nRuns,0., nRuns);//, 100, 0. , 100.);
  hSpreadTimeVsRun->SetDrawOption("E1");
  hSpreadTimeVsRun->SetMarkerStyle(20);
  hSpreadTimeVsRun->SetMarkerColor(kBlue);
  
  TH1F * hAvRawTimeVsRun=new TH1F("hAvRawTimeVsRun","Peak value of raw t^{TOF};;<t_{raw}^{TOF}> (ns)",nRuns,0., nRuns);//, 600, 0. , 600.);
  hAvRawTimeVsRun->SetDrawOption("E1");
  hAvRawTimeVsRun->SetMarkerStyle(21);
  hAvRawTimeVsRun->SetMarkerColor(kGreen);

  TH1F * hPeakRawTimeVsRun=new TH1F("hPeakRawTimeVsRun","Peak value of raw t^{TOF} (landau fit);;t_{peak,raw}^{TOF} (ns)",nRuns,0., nRuns);//, 600, 0. , 600.);
  hPeakRawTimeVsRun->SetDrawOption("E1");
  hPeakRawTimeVsRun->SetMarkerStyle(21);
  hPeakRawTimeVsRun->SetMarkerColor(kGreen);

  TH1F * hSpreadRawTimeVsRun=new TH1F("hSpreadRawTimeVsRun","Spread of raw t^{TOF} (landau fit);;#sigma(t_{raw}^{TOF}) (ns)",nRuns,0., nRuns);//, 100, 0. , 100.);
  hSpreadRawTimeVsRun->SetDrawOption("E1");
  hSpreadRawTimeVsRun->SetMarkerStyle(21);
  hSpreadRawTimeVsRun->SetMarkerColor(kGreen);
   
  TH1F * hAvTotVsRun=new TH1F("hAvTotVsRun","<ToT> (no fit);run;<ToT> (ns)",nRuns,0., nRuns);//, 50, 0. , 50.);
  hAvTotVsRun->SetDrawOption("E1");
  hAvTotVsRun->SetMarkerStyle(22);
   
  TH1F * hPeakTotVsRun=new TH1F("hPeakTotVsRun","<ToT> (gaussian fit);;ToT_{peak} (ns)",nRuns,0., nRuns);//, 50, 0. , 50.);
  hPeakTotVsRun->SetDrawOption("E1");
  hPeakTotVsRun->SetMarkerStyle(22);
   
  TH1F * hSpreadTotVsRun=new TH1F("hSpreadTotVsRun","#sigma(ToT) (gaussian fit);#sigma(ToT) (ns)",nRuns,0., nRuns);//, 50, 0. , 50.);
  hSpreadTotVsRun->SetDrawOption("E1");
  hSpreadTotVsRun->SetMarkerStyle(22);
   
  TH1F * hNegTimeRatioVsRun=new TH1F("hNegTimeRatioVsRun","Ratio of tracks with t^{TOF}<12.5 ns; ; ratio of tracks with t^{TOF}<12.5 ns (%)",nRuns, 0., nRuns);//, 100, 0. , 100.);
  hNegTimeRatioVsRun->SetDrawOption("E");

  TH1F * hOrphansRatioVsRun=new TH1F("hOrphansRatioVsRun","Ratio of orphans (hits with ToT=0); ; ratio of orphans (%)",nRuns, 0., nRuns);//, 1000, 0. , 100.);
  hOrphansRatioVsRun->SetDrawOption("E");

  TH1F * hMeanLVsRun=new TH1F("hMeanLVsRun","Average track length;; <L> (cm)",nRuns, 0., nRuns);//, 350, 350. , 700.);
  hMeanLVsRun->SetDrawOption("E");
  TH1F * hNegLRatioVsRun=new TH1F("hNegLRatioVsRun","Ratio of tracks with L<350 cm;; ratio of tracks with L<350 cm (%)",nRuns, 0., nRuns);//, 1000, 0. , 100.);
  hNegLRatioVsRun->SetDrawOption("E");
  TH1F * hMatchEffVsRun=new TH1F("hMatchEffVsRun","#epsilon_{match} (linear fit for p_{T}>1.0 GeV/c);;#epsilon_{match} (p_{T}>1.0 GeV/c)",nRuns, 0., nRuns);//, 100, 0. , 1.);
  hMatchEffVsRun->SetDrawOption("E");
  TH1F * hMatchEffVsRunNormToGoodCh=new TH1F("hMatchEffVsRunNormToGoodCh","#epsilon_{match} normalized to percentage of TOF good channels;;#epsilon_{match}(p_{T}>1.0 GeV/c,|#eta|<0.8)/f_{all good}",nRuns, 0., nRuns);//, 100, 0. , 1.);
  hMatchEffVsRunNormToGoodCh->SetDrawOption("E");
	
  TH1F * hMatchEffVsRunNormToGoodChInAcc=new TH1F("hMatchEffVsRunNormToGoodChInAcc","#epsilon_{match} normalized to TOF good channels in |#eta|<0.8;;#epsilon_{match}(p_{T}>1.0 GeV/c,|#eta|<0.8/f_{good}(|#eta|<0.8)",nRuns, 0., nRuns);//, 100, 0. , 1.);
  hMatchEffVsRunNormToGoodChInAcc->SetDrawOption("E");

  TH1F * hMatchEffVsRun1=new TH1F("hMatchEffVsRun1","#epsilon_{match}(p_{T}=1.0 GeV/c);;#epsilon_{match} (p_{T}=1.0 GeV/c)",nRuns, 0., nRuns);
  hMatchEffVsRun1->SetDrawOption("E");
  TH1F * hPeakT0AVsRun=new TH1F("hPeakT0AVsRun","Peak value of T0A (gaussian fit);;t0A (ps)",nRuns,0., nRuns);
  TH1F * hPeakT0CVsRun=new TH1F("hPeakT0CVsRun","Peak value of T0C (gaussian fit);;t0AC (ps)",nRuns,0., nRuns);
  TH1F * hPeakT0ACVsRun=new TH1F("hPeakT0ACVsRun","Peak value of T0AC (gaussian fit);;t0AC (ps)",nRuns,0., nRuns);
  TH1F * hT0fillResVsRun=new TH1F("hT0fillResVsRun","t0_fill spread;;t0_spread (ps)",nRuns,0., nRuns);
	
  TH1F * hGoodChannelsRatio=new TH1F("hGoodChannelsRatio","Fraction of TOF good channels;;fraction of good channels",nRuns, 0., nRuns);//, 100, 0. , 1.);
  hGoodChannelsRatio->SetDrawOption("E");

  TH1F * hGoodChannelsRatioInAcc=new TH1F("hGoodChannelsRatioInAcc","Fraction of TOF good channels in |#eta|<0.8;;fraction of good channels in |#eta|<0.8",nRuns, 0., nRuns);//, 100, 0. , 1.);
  hGoodChannelsRatioInAcc->SetDrawOption("E");
	
  lista.Add(hAvMulti);
  lista.Add(hAvDiffTimeVsRun);
  lista.Add(hPeakDiffTimeVsRun);
  lista.Add(hSpreadDiffTimeVsRun);
  lista.Add(hAvTimeVsRun);
  lista.Add(hPeakTimeVsRun);
  lista.Add(hSpreadTimeVsRun);
  lista.Add(  hAvRawTimeVsRun);
  lista.Add(  hPeakRawTimeVsRun);
  lista.Add(  hSpreadRawTimeVsRun); 
  lista.Add(  hAvTotVsRun);
  lista.Add(  hPeakTotVsRun);
  lista.Add(  hSpreadTotVsRun);
  lista.Add(  hNegTimeRatioVsRun);
  lista.Add(  hOrphansRatioVsRun);
  lista.Add( hMeanLVsRun);
  lista.Add(  hNegLRatioVsRun);
  lista.Add(  hMatchEffVsRun);
  lista.Add(hMatchEffVsRunNormToGoodCh);
  lista.Add(hMatchEffVsRunNormToGoodChInAcc);
  lista.Add(hPeakT0AVsRun);
  lista.Add(hPeakT0CVsRun);
  lista.Add(hPeakT0ACVsRun);
  lista.Add(hT0fillResVsRun);
  lista.Add(hGoodChannelsRatio);
  lista.Add(hGoodChannelsRatioInAcc);

  char runlabel[6];
   
  for (Int_t irun=0;irun<nRuns;irun++){
    ttree->GetEntry(irun);
    
    sprintf(runlabel,"%i",runNumber);
    
    hAvMulti->SetBinContent(irun+1, avMulti);
    hAvMulti->GetXaxis()->SetBinLabel(irun+1,runlabel);

    hAvDiffTimeVsRun->SetBinContent(irun+1, avDiffTime);
    hAvDiffTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

    hPeakDiffTimeVsRun->SetBinContent(irun+1,peakDiffTime);
    hPeakDiffTimeVsRun->SetBinError(irun+1,peakDiffTimeErr);
    hPeakDiffTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

    hSpreadDiffTimeVsRun->SetBinContent(irun+1,spreadDiffTime);
    hSpreadDiffTimeVsRun->SetBinError(irun+1,spreadDiffTimeErr);
    hSpreadDiffTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

    hAvTimeVsRun->SetBinContent(irun+1, avTime);
    hAvTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

    hPeakTimeVsRun->SetBinContent(irun+1,peakTime);
    hPeakTimeVsRun->SetBinError(irun+1,peakTimeErr);
    hPeakTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

    hSpreadTimeVsRun->SetBinContent(irun+1,spreadTime);
    hSpreadTimeVsRun->SetBinError(irun+1,spreadTimeErr);
    hSpreadTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

    hAvRawTimeVsRun->SetBinContent(irun+1, avRawTime);
    hAvRawTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
    hPeakRawTimeVsRun->SetBinContent(irun+1,peakRawTime);
    hPeakRawTimeVsRun->SetBinError(irun+1,peakRawTimeErr);
    hPeakRawTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

    hSpreadRawTimeVsRun->SetBinContent(irun+1,spreadRawTime);
    hSpreadRawTimeVsRun->SetBinError(irun+1,spreadRawTimeErr);
    hSpreadRawTimeVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

    hAvTotVsRun->SetBinContent(irun+1,avTot);
    hAvTotVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

    hPeakTotVsRun->SetBinContent(irun+1,peakTot);
    hPeakTotVsRun->SetBinError(irun+1,peakTotErr);
    hPeakTotVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

    hSpreadTotVsRun->SetBinContent(irun+1,spreadTot);
    hSpreadTotVsRun->SetBinError(irun+1,spreadTotErr);
    hSpreadTotVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
    hNegTimeRatioVsRun->SetBinContent(irun+1,negTimeRatio);
    hNegTimeRatioVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
    hOrphansRatioVsRun->SetBinContent(irun+1,orphansRatio);
    hOrphansRatioVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
    hMeanLVsRun->SetBinContent(irun+1,avL);
    hMeanLVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
    hNegLRatioVsRun->SetBinContent(irun+1,negLratio);
    hNegLRatioVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

    hMatchEffVsRun->SetBinContent(irun+1,matchEffLinFit1Gev);
    hMatchEffVsRun->SetBinError(irun+1,matchEffLinFit1GevErr);
    hMatchEffVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);
    hMatchEffVsRun->SetLineColor(kRed);
    hMatchEffVsRun->SetLineWidth(2);

    if (goodChannelRatio>0)
      hMatchEffVsRunNormToGoodCh->SetBinContent(irun+1,matchEffLinFit1Gev/goodChannelRatio);
    else 
      hMatchEffVsRunNormToGoodCh->SetBinContent(irun+1, 0.0);
    hMatchEffVsRunNormToGoodCh->SetBinError(irun+1,matchEffLinFit1GevErr);
    hMatchEffVsRunNormToGoodCh->GetXaxis()->SetBinLabel(irun+1,runlabel);
    hMatchEffVsRunNormToGoodCh->SetLineColor(kCyan+2);
    hMatchEffVsRunNormToGoodCh->SetLineWidth(2);
     
    hGoodChannelsRatio->SetBinContent(irun+1, goodChannelRatio);
    hGoodChannelsRatio->SetLineColor(kCyan-1);
    hGoodChannelsRatio->SetLineWidth(2);
    hGoodChannelsRatio->GetXaxis()->SetBinLabel(irun+1,runlabel);

    if (goodChannelRatioInAcc>0)
      hMatchEffVsRunNormToGoodChInAcc->SetBinContent(irun+1,matchEffLinFit1Gev/goodChannelRatioInAcc);
    else 
      hMatchEffVsRunNormToGoodChInAcc->SetBinContent(irun+1, 0.0);
    hMatchEffVsRunNormToGoodChInAcc->SetBinError(irun+1,matchEffLinFit1GevErr);
    hMatchEffVsRunNormToGoodChInAcc->GetXaxis()->SetBinLabel(irun+1,runlabel);
    hMatchEffVsRunNormToGoodChInAcc->SetLineColor(kBlue);
    hMatchEffVsRunNormToGoodChInAcc->SetLineWidth(2);

    hGoodChannelsRatioInAcc->SetBinContent(irun+1, goodChannelRatioInAcc);
    hGoodChannelsRatioInAcc->SetLineColor(kBlue+2);
    hGoodChannelsRatioInAcc->SetLineWidth(2);
    hGoodChannelsRatioInAcc->GetXaxis()->SetBinLabel(irun+1,runlabel);

    hPeakT0AVsRun->SetBinContent(irun+1,peakT0A);
    hPeakT0AVsRun->SetBinError(irun+1,spreadT0A);
    hPeakT0AVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

    hPeakT0CVsRun->SetBinContent(irun+1,peakT0C);
    hPeakT0CVsRun->SetBinError(irun+1,spreadT0C);
    hPeakT0CVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
    hPeakT0ACVsRun->SetBinContent(irun+1,peakT0AC);
    hPeakT0ACVsRun->SetBinError(irun+1,spreadT0AC);
    hPeakT0ACVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);

    hT0fillResVsRun->SetBinContent(irun+1,avT0fillRes);
    hT0fillResVsRun->GetXaxis()->SetBinLabel(irun+1,runlabel);
  }
  
  TFile * fout=new TFile(outfilename,"recreate");
  fout->cd();
  lista.Write();
  fout->Close();
    
  gStyle->SetOptStat(10);

  TCanvas* cPeakDiffTimeVsRun = new TCanvas("cPeakDiffTimeVsRun","cPeakDiffTimeVsRun", 50,50,1050, 550);
  hPeakDiffTimeVsRun->GetYaxis()->SetRangeUser(-50.,50.);
  hPeakDiffTimeVsRun->Draw();
  cPeakDiffTimeVsRun->Print(Form("%s/cPeakDiffTimeVsRun.png",plotDir.Data()));
	
  TCanvas* cSpreadDiffTimeVsRun = new TCanvas("cSpreadDiffTimeVsRun","cSpreadDiffTimeVsRun", 50,50,1050, 550);
  hSpreadDiffTimeVsRun->GetYaxis()->SetRangeUser(200.,400.);
  hSpreadDiffTimeVsRun->Draw();
  cSpreadDiffTimeVsRun->Print(Form("%s/cSpreadDiffTimeVsRun.png",plotDir.Data()));
  
  TCanvas* cMatchEffVsRun = new TCanvas("cMatchEffVsRun","cMatchEffVsRun",50, 50, 1050, 550);
  hMatchEffVsRun->GetYaxis()->SetRangeUser(0.,1.);
  hMatchEffVsRun->Draw();
  cMatchEffVsRun->Print(Form("%s/cMatchEffVsRun.png",plotDir.Data()));
  		
  TCanvas* cMatchEffNormToGoodChInAcc = new TCanvas("cMatchEffNormToGoodChInAcc","cMatchEffNormToGoodChInAcc",50, 50,1050, 550);
  hMatchEffVsRunNormToGoodChInAcc->GetYaxis()->SetRangeUser(0.,1.);
  hMatchEffVsRunNormToGoodChInAcc->Draw();
  cMatchEffNormToGoodChInAcc->Print(Form("%s/cMatchEffNormToGoodChInAcc.png",plotDir.Data()));

  TCanvas* cMatchEffNormToGoodCh = new TCanvas("cMatchEffNormToGoodCh","cMatchEffNormToGoodCh",50, 50,1050, 550);
  hMatchEffVsRunNormToGoodCh->GetYaxis()->SetRangeUser(0.,1.);
  hMatchEffVsRunNormToGoodCh->Draw();
  cMatchEffNormToGoodCh->Print(Form("%s/cMatchEffNormToGoodCh.png",plotDir.Data()));

  TCanvas* cMatchEffSummary = new TCanvas("cMatchEffSummary","cMatchEffSummary",50, 50,1050, 550);
  hMatchEffVsRun->GetYaxis()->SetRangeUser(0.,1.);
  hMatchEffVsRun->Draw();
  hMatchEffVsRunNormToGoodCh->Draw("same");
  hMatchEffVsRunNormToGoodChInAcc->Draw("same");
  cMatchEffSummary->Print(Form("%s/cMatchEffSummary.png",plotDir.Data()));

  TCanvas* cGoodCh = new TCanvas("cGoodCh","cGoodCh",50, 50,1050, 550);
  hGoodChannelsRatio->GetYaxis()->SetRangeUser(0.75,1.);
  hGoodChannelsRatio->Draw();
  cGoodCh->Print(Form("%s/cGoodCh.png",plotDir.Data()));

  TCanvas* cGoodChInAcc = new TCanvas("cGoodChInAcc","cGoodChInAcc",50, 50,1050, 550);
  hGoodChannelsRatioInAcc->GetYaxis()->SetRangeUser(0.75,1.);
  hGoodChannelsRatioInAcc->Draw();
  cGoodChInAcc->Print(Form("%s/cGoodChInAcc.png",plotDir.Data()));

  if (displayAll) {	
    TCanvas* cPeakT0AVsRun = new TCanvas("cPeakT0AVsRun","cPeakT0AVsRun", 50,50,1050, 550);
    hPeakT0AVsRun->Draw();
    cPeakT0AVsRun->Print(Form("%s/cPeakT0AVsRun.png",plotDir.Data()));
  
    TCanvas* cPeakT0CVsRun = new TCanvas("cPeakT0CVsRun","cPeakT0CVsRun", 50,50,1050, 550);
    hPeakT0CVsRun->Draw();
    cPeakT0CVsRun->Print(Form("%s/cPeakT0CVsRun.png",plotDir.Data()));
  
    TCanvas* cPeakT0ACVsRun = new TCanvas("cPeakT0ACVsRun","cPeakT0ACVsRun", 50,50,1050, 550);
    hPeakT0ACVsRun->Draw();
    cPeakT0ACVsRun->Print(Form("%s/cPeakT0ACVsRun.png",plotDir.Data()));
  
    TCanvas* cT0fillResVsRun = new TCanvas("cT0fillResVsRun","cT0fillResVsRun", 50,50,1050, 550);
    hT0fillResVsRun->Draw();
    cT0fillResVsRun->Print(Form("%s/cT0fillResVsRun.png",plotDir.Data()));

    TCanvas* cAvDiffTimeVsRun = new TCanvas("cAvDiffTimeVsRun","cAvDiffTimeVsRun",50,50,1050, 550);
    gPad->SetGridx();
    gPad->SetGridy();
    hAvDiffTimeVsRun->Draw();
    cAvDiffTimeVsRun->Print(Form("%s/cAvDiffTimeVsRun.png",plotDir.Data()));
	  
    TCanvas* cAvTimeVsRun = new TCanvas("cAvTimeVsRun","cAvTimeVsRun", 50,50,1050, 550);
    hAvTimeVsRun->Draw();
    cAvTimeVsRun->Print(Form("%s/cAvTimeVsRun.png",plotDir.Data()));
	  
    TCanvas* cPeakTimeVsRun = new TCanvas("cPeakTimeVsRun","cPeakTimeVsRun", 50,50,1050, 550);
    hPeakTimeVsRun->Draw();
    cPeakTimeVsRun->Print(Form("%s/cPeakTimeVsRun.png",plotDir.Data()));
	  
    TCanvas* cSpreadTimeVsRun = new TCanvas("cSpreadTimeVsRun","cSpreadTimeVsRun", 50,50,1050, 550);
    hSpreadTimeVsRun->Draw();
    cSpreadTimeVsRun->Print(Form("%s/cSpreadTimeVsRun.png",plotDir.Data()));
	  
    TCanvas* cAvRawTimeVsRun = new TCanvas("cAvRawTimeVsRun","cAvRawTimeVsRun", 50,50,1050, 550);
    hAvRawTimeVsRun->Draw();
    cAvRawTimeVsRun->Print(Form("%s/cAvRawTimeVsRun.png",plotDir.Data()));
	  
    TCanvas* cPeakRawTimeVsRun = new TCanvas("cPeakRawTimeVsRun","cPeakRawTimeVsRun", 50,50,1050, 550);
    hPeakRawTimeVsRun->Draw();
    cPeakRawTimeVsRun->Print(Form("%s/cPeakRawTimeVsRun.png",plotDir.Data()));
	  
    TCanvas* cSpreadRawTimeVsRun = new TCanvas("cSpreadRawTimeVsRun","cSpreadRawTimeVsRun", 50,50,1050, 550);
    hSpreadRawTimeVsRun->Draw();
    cSpreadRawTimeVsRun->Print(Form("%s/cSpreadRawTimeVsRun.png",plotDir.Data()));
	  
    TCanvas* cAvTotVsRun = new TCanvas("cAvTotVsRun","cAvTotVsRun", 50,50,1050, 550);
    hAvTotVsRun->Draw();
    cAvTotVsRun->Print(Form("%s/cAvTotVsRun.png",plotDir.Data()));
	  
    TCanvas* cPeakTotVsRun = new TCanvas("cPeakTotVsRun","cPeakTotVsRun", 50,50,1050, 550);
    hPeakTotVsRun->Draw();
    cPeakTotVsRun->Print(Form("%s/cPeakTotVsRun.png",plotDir.Data()));
	  
    TCanvas* cSpreadTotVsRun = new TCanvas("cSpreadTotVsRun","cSpreadTotVsRun", 50,50,1050, 550);
    hSpreadTotVsRun->Draw();
    cSpreadTotVsRun->Print(Form("%s/cSpreadTotVsRun.png",plotDir.Data()));
	  
    TCanvas* cNegTimeRatioVsRun = new TCanvas("cNegTimeRatioVsRun","cNegTimeRatioVsRun", 50,50,1050, 550);
    hNegTimeRatioVsRun->Draw();
    cNegTimeRatioVsRun->Print(Form("%s/cNegTimeRatioVsRun.png",plotDir.Data()));
	  
    TCanvas* cOrphansRatioVsRun = new TCanvas("cOrphansRatioVsRun","cOrphansRatioVsRun", 50,50,1050, 550);
    hOrphansRatioVsRun->Draw();
    cOrphansRatioVsRun->Print(Form("%s/cOrphansRatioVsRun.png",plotDir.Data()));
	  
    TCanvas* cMeanLVsRun = new TCanvas("cMeanLVsRun","cMeanLVsRun", 50,50,1050, 550);
    hMeanLVsRun->Draw();
    cMeanLVsRun->Print(Form("%s/cMeanLVsRun.png",plotDir.Data()));
	  
    TCanvas* cNegLRatioVsRun = new TCanvas("cNegLRatioVsRun","cNegLRatioVsRun", 50,50,1050, 550);
    hNegLRatioVsRun->Draw();
    cNegLRatioVsRun->Print(Form("%s/cNegLRatioVsRun.png",plotDir.Data()));
  }
	
  return 0;
}
