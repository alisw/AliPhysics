/*

root.exe -l $ALICE_ROOT/TPC/Upgrade/macros/{loadlibs.C,ConfigOCDB.C}
.x $ALICE_ROOT/TPC/Upgrade/macros/toyMCRecPlots.C

*/

void toyMCRecPlots(TString inFileName = "toyMC.debug.root",Bool_t doPlots = kFALSE){
  //
  // do all the plots for Toy MC Reconstruction for the TDR 
  //

  //gStyle->SetOptStat(0);
  st->SetTitleX(0.17);
  st->SetTitleW(0.73);
  
  // parameters
  const Int_t nT0          = 2;
  const Int_t nT02D        = 2;
  const Int_t nZ0          = 2;
  const Int_t nTrackParams = 9;
  const Int_t nTrackParamsITS = 9;
  const Int_t nEff         = 1;
  Int_t nTrackParamsPlot=5;
  Int_t ydiv=2;  //divisions in y for track parameters plot, usually 3
  Int_t col = kBlack;

  TString sT0[nT0] = {"fTime0-t0","fTime0-t0+z0*TMath::Sign(1,tOrig.Eta())/vDrift"};
  //TString sT02D[nT02D] = {"fTime0-t0+z0*TMath::Sign(1,tOrig.Eta())/vDrift:tOrig.fP[3]","fTime0-t0+z0*TMath::Sign(1,tOrig.Eta())/vDrift:tOrig.fP[4]"};
  TString sT02D[nT02D] = {"(fTime0-t0+z0*TMath::Sign(1,tOrig.Eta())/vDrift)*vDrift:tOrig.fP[3]","(fTime0-t0+z0*TMath::Sign(1,tOrig.Eta())/vDrift)*vDrift:tOrig.fP[4]"};
  //TString sT02D[nT02D] = {"(fTime0-t0)*vDrift:tOrig.fP[3]","(fTime0-t0)*vDrift:tOrig.fP[4]"};
  TString sTrackParams[nTrackParams] = {"track.fP[0]-tOrig.fP[0]","track.fP[1]-tOrig.fP[1]","track.fP[2]-tOrig.fP[2]","track.fP[3]-tOrig.fP[3]","track.fP[4]-tOrig.fP[4]","track.fAlpha-tOrig.fAlpha","track.fX/tOrig.fX","track.fP[0]/tOrig.fP[0]","track.fP[1]/tOrig.fP[1]"};
  TString sTrackParamsITS[nTrackParamsITS] = {"trackITS.fP[0]-tOrigITS.fP[0]","trackITS.fP[1]-tOrigITS.fP[1]","trackITS.fP[2]-tOrigITS.fP[2]","trackITS.fP[3]-tOrigITS.fP[3]","trackITS.fP[4]-tOrigITS.fP[4]","trackITS.fAlpha-tOrigITS.fAlpha","trackITS.fX/tOrigITS.fX","trackITS.fP[0]/tOrigITS.fP[0]","trackITS.fP[1]/tOrigITS.fP[1]"};
  TString sTrackParamsITS1[nTrackParamsITS] = {"trackITS1.fP[0]-tOrigITS1.fP[0]","trackITS1.fP[1]-tOrigITS1.fP[1]","trackITS1.fP[2]-tOrigITS1.fP[2]","trackITS1.fP[3]-tOrigITS1.fP[3]","trackITS1.fP[4]-tOrigITS1.fP[4]","trackITS1.fAlpha-tOrigITS1.fAlpha","trackITS1.fX/tOrigITS1.fX","trackITS1.fP[0]/tOrigITS1.fP[0]","trackITS1.fP[1]/tOrigITS1.fP[1]"};
  TString sTrackParamsITS2[nTrackParamsITS] = {"trackITS2.fP[0]-tOrigITS2.fP[0]","trackITS2.fP[1]-tOrigITS2.fP[1]","trackITS2.fP[2]-tOrigITS2.fP[2]","trackITS2.fP[3]-tOrigITS2.fP[3]","trackITS2.fP[4]-tOrigITS2.fP[4]","trackITS2.fAlpha-tOrigITS2.fAlpha","trackITS2.fX/tOrigITS2.fX","trackITS2.fP[0]/tOrigITS2.fP[0]","trackITS2.fP[1]/tOrigITS2.fP[1]"};
  TString sZ0[nZ0] = {"(fTime0-t0)*vDrift","(fTime0-t0+z0*TMath::Sign(1,tOrig.Eta())/vDrift)*vDrift"};
  TString sEff[nEff] = {"nSeedClustersID/nClustersMC:nSeedClusters/nClustersMC"};


  TString tT0[nT0] = {"T_{0} resolution","intrinsic T_{0} resolution"};
  TString tT02D[nT0] = {"intrinsic T_{0} resolution vs. tan#lambda","intrinsic T_{0} resolution vs. 1/p_{T}"};
  TString tTrackParams[nTrackParams] = {"local y","z","sin inclination angle","tan#lambda","1/p_{T}","#alpha","X","Y","Z"};
  TString tTrackParamsITS[nTrackParamsITS] = {"local y","z","sin inclination angle","tan#lambda","1/p_{T}","#alpha","X","Y","Z"};
  TString tZ0[nZ0] = {"Z_{0} resolution","intrinsic Z_{0} resolution"};
  TString tEff[nEff] = {"nSeedClustersID/nClustersMC:nSeedClusters/nClustersMC"};


  TString sSel = "fTime0>-1"; // seeding successful

  // retrieve configuration string
  TPRegexp reg(".*([0-9]_[0-9]_[0-9]_[0-9]{3}_[0-9]{2}).*debug.root");
  TObjArray *arrMatch=0x0;
  arrMatch=reg.MatchS(inFileName);
  TString sConfig = arrMatch->At(1)->GetName();
   
  // get file
  TFile *fIn  = TFile::Open(inFileName.Data(),"READ");
  if(!fIn){
    Printf("No file %s found",inFileName.Data());
    return;
  }

  // get tree
  TTree *Tracks = fIn->Get("Tracks");
  if(!Tracks){
    Printf("No TTree found");
    return;
  }

 // output canvases
  TCanvas *cT0 = new TCanvas(Form("cT0_%s",sConfig.Data()),Form("cT0_%s",sConfig.Data()),1200,500);
  cT0->Divide(2,1);

  TCanvas *cT02D = new TCanvas(Form("cT02D_%s",sConfig.Data()),Form("cT02D_%s",sConfig.Data()),1200,500);
  cT02D->Divide(2,1);

  TCanvas *cZ0 = new TCanvas(Form("cZ0_%s",sConfig.Data()),Form("cZ0_%s",sConfig.Data()),1200,500);
  cZ0->Divide(2,1);

  //TCanvas *cTrackParams = new TCanvas(Form("cTrackParams_%s",sConfig.Data()),Form("cTrackParams_%s",sConfig.Data()),1200,900/3.*ydiv);
  //cTrackParams->Divide(3,ydiv);
  TCanvas *cTrackParams = new TCanvas(Form("cTrackParams_%s",sConfig.Data()),Form("cTrackParams_%s",sConfig.Data()),4800,3600);
  cTrackParams->Divide(3,3);

  //TCanvas *cTrackParamsITS = new TCanvas(Form("cTrackParamsITS_%s",sConfig.Data()),Form("cTrackParamsITS_%s",sConfig.Data()),1200,900/3.*ydiv);
  //cTrackParamsITS->Divide(3,ydiv);
  
  TCanvas *cTrackParamsITS1 = new TCanvas(Form("cTrackParamsITS1_%s",sConfig.Data()),Form("cTrackParamsITS1_%s",sConfig.Data()),1200,900/3.*ydiv);
  cTrackParamsITS1->Divide(3,ydiv);
  TCanvas *cTrackParamsITS = new TCanvas(Form("cTrackParamsITS_%s",sConfig.Data()),Form("cTrackParamsITS_%s",sConfig.Data()),4800,3600);
  cTrackParamsITS->Divide(3,3);

  TCanvas *cTrackParamsITS2 = new TCanvas(Form("cTrackParamsITS2_%s",sConfig.Data()),Form("cTrackParamsITS2_%s",sConfig.Data()),1200,900/3.*ydiv);
  cTrackParamsITS2->Divide(3,ydiv);
  
  TCanvas *cEff = new TCanvas(Form("cEff_%s",sConfig.Data()),Form("cEff_%s",sConfig.Data()),1200,900);
  //cEff->Divide(2,1);

 
  // legends
  TLegend *l[nT0];



  // draw T0 resolution
  for(Int_t iT0 = 0; iT0 < nT0; iT0 ++){

    cT0->cd(iT0+1);
    TStatToolkit::DrawHistogram(Tracks,sT0[iT0].Data(),sSel.Data(),Form("hT0_%s_%d",sConfig.Data(),iT0),Form("%s",tT0[iT0].Data()),3);

    //hT0[iT0]->Fit("gaus","","");
    //l[iT0]= new TLegend(0.55,0.7,0.8,0.8,Form(""));
    //myLegendSetUp(l[iT0],0.03);
    //l[iT0]->AddEntry( hT0[iT0], Form("#sigma = (%.2f #pm %.2f) #mus",hT0[iT0]->GetFunction("gaus")->GetParameter(2)*1e6,hT0[iT0]->GetFunction("gaus")->GetParError(2)*1e6));
    //l[iT0]->Draw();

  }


  // draw T0 resolution (2D)
  for(Int_t iT02D = 0; iT02D < nT02D; iT02D ++){


    cT02D->cd(iT02D+1);
    TStatToolkit::DrawHistogram(Tracks,sT02D[iT02D].Data(),sSel.Data(),Form("hT02D_%s_%d",sConfig.Data(),iT02D),Form("%s",tT02D[iT02D].Data()),3);

  }

  // draw Z0 resolution
  for(Int_t iZ0 = 0; iZ0 < nZ0; iZ0 ++){

    cZ0->cd(iZ0+1);
    TStatToolkit::DrawHistogram(Tracks,sZ0[iZ0].Data(),sSel.Data(),Form("hZ0_%s_%d",sConfig.Data(),iZ0),Form("%s",tZ0[iZ0].Data()),3);

  }

  // draw track parameters
  for(Int_t iTrackParams = 0; iTrackParams < nTrackParamsPlot ; iTrackParams ++){

    cTrackParams->cd(iTrackParams+1);
    if(inFileName.Contains("allClusters") && iTrackParams==nTrackParams-1) 
      TStatToolkit::DrawHistogram(Tracks,sTrackParams[iTrackParams].Data(),sSel.Data(),Form("hTrackParams_%s_%d",sConfig.Data(),iTrackParams),Form("%s",tTrackParams[iTrackParams].Data()),5);
    else
      TStatToolkit::DrawHistogram(Tracks,sTrackParams[iTrackParams].Data(),sSel.Data(),Form("hTrackParams_%s_%d",sConfig.Data(),iTrackParams),Form("%s",tTrackParams[iTrackParams].Data()),6);

  }


  // draw track parameters at ITS outer layer
  for(Int_t iTrackParamsITS = 0; iTrackParamsITS < nTrackParamsPlot; iTrackParamsITS ++){

    cTrackParamsITS->cd(iTrackParamsITS+1);
    TStatToolkit::DrawHistogram(Tracks,sTrackParamsITS[iTrackParamsITS].Data(),sSel.Data(),Form("hTrackParamsITS_%s_%d",sConfig.Data(),iTrackParamsITS),Form("%s",tTrackParamsITS[iTrackParamsITS].Data()),6);

    cTrackParamsITS1->cd(iTrackParamsITS+1);
    TStatToolkit::DrawHistogram(Tracks,sTrackParamsITS1[iTrackParamsITS].Data(),sSel.Data(),Form("hTrackParamsITS1_%s_%d",sConfig.Data(),iTrackParamsITS),Form("%s",tTrackParamsITS[iTrackParamsITS].Data()),6);

    cTrackParamsITS2->cd(iTrackParamsITS+1);
    TStatToolkit::DrawHistogram(Tracks,sTrackParamsITS2[iTrackParamsITS].Data(),sSel.Data(),Form("hTrackParamsITS2_%s_%d",sConfig.Data(),iTrackParamsITS),Form("%s",tTrackParamsITS[iTrackParamsITS].Data()),6);
    
  }


  // draw cluster efficiency
  if(inFileName.Contains("allClusters")){  
    for(Int_t iEff = 0; iEff < nEff; iEff ++){
      
      //cEff->cd(iEff+1);
      cEff->cd();
      TStatToolkit::DrawHistogram(Tracks,sEff[iEff].Data(),sSel.Data(),Form("hEff_%s_%d",sConfig.Data(),iEff),Form("%s",tEff[iEff].Data()),3);   
      
    }
  }
  
  // plots
  if(doPlots){
    TString outFileName = gSystem->BaseName(inFileName.Data());
    outFileName.ReplaceAll(".root","");
//     cT0->SaveAs(Form("%s_T0.eps",outFileName.Data()));
//     cZ0->SaveAs(Form("%s_Z0.eps",outFileName.Data()));
//     cT02D->SaveAs(Form("%s_T02D.eps",outFileName.Data()));
//     cTrackParams->SaveAs(Form("%s_TrackParams.eps",outFileName.Data()));
//     cTrackParamsITS->SaveAs(Form("%s_TrackParamsITS.eps",outFileName.Data()));
//     if(inFileName.Contains("allClusters"))
//       cEff->SaveAs(Form("%s_Eff.eps",outFileName.Data()));

    cT0->SaveAs(Form("%s_T0.png",outFileName.Data()));
    cZ0->SaveAs(Form("%s_Z0.png",outFileName.Data()));
    cT02D->SaveAs(Form("%s_T02D.png",outFileName.Data()));
    cTrackParams->SaveAs(Form("%s_TrackParams.png",outFileName.Data()));
    cTrackParamsITS->SaveAs(Form("%s_TrackParamsITS.png",outFileName.Data()));
    cTrackParamsITS1->SaveAs(Form("%s_TrackParamsITS1.png",outFileName.Data()));
    cTrackParamsITS2->SaveAs(Form("%s_TrackParamsITS2.png",outFileName.Data()));
    if(inFileName.Contains("allClusters"))
      cEff->SaveAs(Form("%s_Eff.eps",outFileName.Data()));

  }
}

void myLegendSetUp(TLegend *currentLegend = NULL,Float_t currentTextSize=0.07){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(1);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  return;
}

void myHistoSetUp( TH1 *h = NULL, Int_t col = kBlack, Int_t mar = 20){

  h->SetLineColor(col);
  h->SetMarkerColor(col);
  h->SetMarkerStyle(mar);

}
		   
