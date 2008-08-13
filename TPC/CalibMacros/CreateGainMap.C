/*

This macro creates a gain map for the TPC based on the results of the Krypton calibration.
The main steps are the following:

1. Define outlier-pads where the krypton calibration was not succesful
2. A parabolic fit for the whole chamber is performed
3. replace outliers with fitted values
4. normalize separately IROCs and OROCs

For more details see below.



example usage:
==============

TFile f("calibKr.root")
AliTPCCalPad * kryptonRaw = new AliTPCCalPad(*fitMean)
AliTPCCalPad * kryptonMean = new AliTPCCalPad(*spectrMean)
AliTPCCalPad * kryptonChi2 = new AliTPCCalPad(*fitNormChi2)
AliTPCCalPad * kryptonRMS = new AliTPCCalPad(*fitRMS)

.L CreateGainMap.C
AliTPCCalPad * final = CreateGainMap(kryptonRaw, kryptonRMS)

TFile *h = new TFile("GainMap.root", "RECREATE")
final.Write()

*/


AliTPCCalPad * CreateGainMap(AliTPCCalPad *krypFitMean, AliTPCCalPad *krypFitRMS, AliTPCCalPad *noiseMap = 0, AliTPCCalPad *krypSpectrMean = 0, 
			     AliTPCCalPad *krypChi2 = 0, AliTPCCalPad *pulser = 0, AliTPCCalPad *electrode = 0) {

  
  // Draw input map
  TCanvas *test3 = new TCanvas("ASIDE3", "Aoriginal");
  krypFitMean->MakeHisto2D()->Draw("colz");
  TCanvas *test4 = new TCanvas("CSIDE4", "Coriginal");
  krypFitMean->MakeHisto2D(1)->Draw("colz");
  TH1F * hDEBUG = new TH1F("bla", "fitRMS/fitMean",100,0,0.3);

  const Double_t kryptonMean = 0.05012;
  const Double_t kryptonSigma = 0.00386;

  TObjArray arrayRocFinal(72);

  //
  // Loop over all sectors
  //
  for(Int_t isector=0; isector < 72; isector++){

    AliTPCCalROC *rocKrypFitMean = krypFitMean->GetCalROC(isector);
    AliTPCCalROC *rocKrypFitRMS  = krypFitRMS->GetCalROC(isector);

    AliTPCCalROC *rocOutlier = new AliTPCCalROC(*rocKrypFitMean);

    // 1. define map of outliers: the 41keV peak is fitted and krypFitMean represents the peak position and krypFitRMS its sigma.
    //    In order to control the fit qualtiy krypFitRMS/krypFitMean is computed and plotted; it is roughly gaussian distributed
    //    with the mean resolution kryptonMean = 0.05012 and sigma kryptonSigma = 0.00386. If the deviation is larger than 4sigma
    //    the fit to the 41keV peak is considered as failed.

    Int_t nPointsFit = 0;
    for(Int_t iChannel=0; iChannel < rocOutlier->GetNchannels(); iChannel++){
      Double_t fitMean = rocKrypFitMean->GetValue(iChannel);
      Double_t fitRMS = rocKrypFitRMS->GetValue(iChannel);
      if (fitRMS < 0.001 || fitMean < 0.001) continue;
       hDEBUG->Fill(fitRMS/fitMean);
      if (TMath::Abs(fitRMS/fitMean - kryptonMean) < 4*kryptonSigma || fitRMS < 0.001 || fitMean < 0.001) {
	rocOutlier->SetValue(iChannel,0);
	nPointsFit++;
      }
    }
    //TCanvas *DEBUG = new TCanvas();
    //rocOutlier->MakeHisto2D()->Draw("colz");
    
    Double_t rocMedian = rocKrypFitMean->GetMedian(rocOutlier);
    

    // 2. make a parabolic fit for the whole chamber with excluded outliers
    
    TVectorD params;
    TMatrixD cov;
    Float_t chi2;
    rocKrypFitMean->GlobalFit(rocOutlier,0,params,cov,chi2);
    AliTPCCalROC *rocParabolicFit=AliTPCCalROC::CreateGlobalFitCalROC(params,isector);

    if (nPointsFit != 0) cout << "sector: "<< isector << " chi2: " << chi2/nPointsFit <<" median: "<< rocMedian <<" n: " << nPointsFit << endl;
    if (nPointsFit == 0) {
      AliTPCCalROC *rocFinal = new AliTPCCalROC(*rocParabolicFit); // What to do with sectors being switched off??
      arrayRocFinal.AddAt(rocFinal,isector);
      continue;
    }

    //
    // VERY IMPORTANT **** VERY IMPORTANT  **** VERY IMPORTANT  **** VERY IMPORTANT  **** VERY IMPORTANT  **** VERY IMPORTANT 
    //
    // if the fit is considered as failed, the mean value is taken for the whole chamber (TO AVOID THIS, RELEASE THIS CUT)
    // if cosmic data can be used for gain calibration, this cut should be strong (e.g. 5) !!!!!!!!!!

    if (chi2/nPointsFit > 5) {
      AliTPCCalROC *rocFinal = new AliTPCCalROC(*rocParabolicFit);
      if (rocMedian == 0 && isector == 51) rocMedian = 1075; // manual patch to be removed for sector 51!
      for(Int_t iChannel=0; iChannel < rocFinal->GetNchannels(); iChannel++) rocFinal->SetValue(iChannel,rocMedian);
      arrayRocFinal.AddAt(rocFinal,isector);
      continue;
    }

    // 3. replace outliers with fitted values

    AliTPCCalROC *rocFinal = new AliTPCCalROC(*rocParabolicFit);
    for(Int_t iChannel=0; iChannel < rocFinal->GetNchannels(); iChannel++){
      Double_t fitMean = rocKrypFitMean->GetValue(iChannel);
      Double_t fitRMS = rocKrypFitRMS->GetValue(iChannel);   
      if (fitRMS < 0.001 || fitMean < 0.001 || TMath::Abs(rocParabolicFit->GetValue(iChannel)/fitMean - 1) > 0.35) {
	rocFinal->SetValue(iChannel,rocParabolicFit->GetValue(iChannel));
	continue;
      }
      if (TMath::Abs(fitRMS/fitMean - kryptonMean) < 4*kryptonSigma) {
	rocFinal->SetValue(iChannel,rocKrypFitMean->GetValue(iChannel));
      }
    }
    

    // 4. Postprocessing: Set dead channels and very noisy channels (time dependent) to 0
    
    const Double_t noiseMin = 0.01;
    const Double_t noiseMax = 2;

    if (noiseMap) {
      AliTPCCalROC *rocNoise = noiseMap->GetCalROC(isector);
      for(Int_t iChannel=0; iChannel < rocFinal->GetNchannels(); iChannel++){
	Double_t noise = rocNoise->GetValue(iChannel);
	if (noise < noiseMin || noise > noiseMax) rocFinal->SetValue(iChannel, 0);
      }
    }
    // Fill an array of ROCs

    arrayRocFinal.AddAt(rocFinal,isector);

  }

  AliTPCCalPad *final = new AliTPCCalPad(&arrayRocFinal);
  
  // 4. normalize separately IROCs and OROCs to the mean of their medians

  Double_t meanMedianIROC;
  Double_t meanMedianOROC;
  Int_t n = 0;

  for(Int_t isector=0; isector < 36; isector++){ // IROCs
    AliTPCCalROC *rocFinal = final->GetCalROC(isector);
    if (rocFinal->GetMedian() != 0) {
      meanMedianIROC += rocFinal->GetMedian();
      n++;
    }
  }
  meanMedianIROC = meanMedianIROC/n;

  n = 0;
  for(Int_t isector=36; isector < 72; isector++){ // OROCs
    AliTPCCalROC *rocFinal = final->GetCalROC(isector);
    if (rocFinal->GetMedian() != 0) {
      meanMedianOROC += rocFinal->GetMedian();
      n++;
    }
  }
  meanMedianOROC = meanMedianOROC/n;

  for(Int_t isector=0; isector < 72; isector++){ // OROCs
    AliTPCCalROC *rocFinal = final->GetCalROC(isector);
    if (isector<36) rocFinal->Multiply(1./meanMedianIROC);
    if (isector>35) rocFinal->Multiply(1./meanMedianOROC);
  }

  // Draw results
  TCanvas *test = new TCanvas("ASIDE", "A");
  final->MakeHisto2D()->Draw("colz");
  TCanvas *test2 = new TCanvas("CSIDE", "C");
  final->MakeHisto2D(1)->Draw("colz");
  TCanvas *cDEBUG = new TCanvas();
  hDEBUG->Draw();

  //return results
  final->SetName("GainMap");
  return final;

}

