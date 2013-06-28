/*

root.exe -l $ALICE_ROOT/TPC/Upgrade/macros/{loadlibs.C,ConfigOCDB.C}
.x $ALICE_ROOT/TPC/Upgrade/macros/toyMCRecPlots.C

*/

void toyMCRecPlots(TString inFileName = "toyMC.debug.root"){
  //
  // do all the plots for Toy MC Reconstruction for the TDR 
  //

  //gStyle->SetOptStat(0);

  // parameters
  const Int_t nT0          = 2;
  const Int_t nT02D        = 2;
  const Int_t nTrackParams = 9;

  Int_t col = kBlack;

  TString sT0[nT0] = {"fTime0-t0","fTime0-t0+z0*TMath::Sign(1,tOrig.Eta())/vdrift"};
  TString sT02D[nT02D] = {"fTime0-t0+z0*TMath::Sign(1,tOrig.Eta())/vdrift:tOrig.fP[3]","fTime0-t0+z0*TMath::Sign(1,tOrig.Eta())/vdrift:tOrig.fP[4]"};
  TString sTrackParams[nTrackParams] = {"track.fP[0]-tOrig.fP[0]","track.fP[1]-tOrig.fP[1]","track.fP[2]-tOrig.fP[2]","track.fP[3]-tOrig.fP[3]","track.fP[4]-tOrig.fP[4]","track.fAlpha-tOrig.fAlpha","track.fX/tOrig.fX","track.fY/tOrig.fY","track.fZ/tOrig.fZ"};

  TString tT0[nT0] = {"T_{0} resolution","intrinsic T_{0} resolution"};
  TString tT02D[nT0] = {"intrinsic T_{0} resolution vs. tan#lambda","intrinsic T_{0} resolution vs. 1/p_{T}"};
  TString tTrackParams[nTrackParams] = {"local y","z","sin inclination angle","tan#lambda","1/p_{T}","#alpha","X","Y","Z"};

  TString sSel = "fTime0>-1"; // seeding successful

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
  TCanvas *cT0 = new TCanvas("cT0","cT0",1200,500);
  cT0->Divide(2,1);

  TCanvas *cT02D = new TCanvas("cT02D","cT02D",1200,500);
  cT02D->Divide(2,1);

  TCanvas *cTrackParams = new TCanvas("cTrackParams","cTrackParams",1200,900);
  cTrackParams->Divide(3,3);

 
  // legends
  TLegend *l[nT0];



  // draw T0 resolution
  for(Int_t iT0 = 0; iT0 < nT0; iT0 ++){

    cT0->cd(iT0+1);
    TStatToolkit::DrawHistogram(Tracks,sT0[iT0].Data(),sSel.Data(),Form("hT0%d",iT0),Form("%s",tT0[iT0].Data()),3);

    //hT0[iT0]->Fit("gaus","","");
    //l[iT0]= new TLegend(0.55,0.7,0.8,0.8,Form(""));
    //myLegendSetUp(l[iT0],0.03);
    //l[iT0]->AddEntry( hT0[iT0], Form("#sigma = (%.2f #pm %.2f) #mus",hT0[iT0]->GetFunction("gaus")->GetParameter(2)*1e6,hT0[iT0]->GetFunction("gaus")->GetParError(2)*1e6));
    //l[iT0]->Draw();

  }


  // draw T0 resolution (2D)
  for(Int_t iT02D = 0; iT02D < nT02D; iT02D ++){


    cT02D->cd(iT02D+1);
    TStatToolkit::DrawHistogram(Tracks,sT02D[iT02D].Data(),sSel.Data(),Form("hT02D%d",iT0),Form("%s",tT02D[iT02D].Data()),3);

  }


  // draw track parameters
  for(Int_t iTrackParams = 0; iTrackParams < nTrackParams; iTrackParams ++){

    cTrackParams->cd(iTrackParams+1);
    TStatToolkit::DrawHistogram(Tracks,sTrackParams[iTrackParams].Data(),sSel.Data(),Form("hTrackParams%d",iTrackParams),Form("%s",tTrackParams[iTrackParams].Data()),3);

  }


  // plots
  TString outFileName = inFileName;
  outFileName.ReplaceAll(".root","");
  cT0->SaveAs(Form("%s_T0.eps",outFileName.Data()));
  cT02D->SaveAs(Form("%s_T02D.eps",outFileName.Data()));
  cTrackParams->SaveAs(Form("%s_TrackParams.eps",outFileName.Data()));

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
		   
