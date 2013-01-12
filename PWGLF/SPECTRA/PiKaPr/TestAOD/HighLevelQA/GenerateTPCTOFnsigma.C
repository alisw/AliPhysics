////////////////////////////////////////////////////////////
// GenerateTPCTOFnsigma.C                                 //
// creates a histogram of the TPC+TOF nsigma distribution //
// written by John Groh                                   //
////////////////////////////////////////////////////////////

void GenerateTPCTOFnsigma()
{
  //-------------------------------------------
  // Get TPC+TOF nsigma projections from data -
  //-------------------------------------------
  const Int_t nPart = 3;
  TString Particle[nPart] = {"Pion","Kaon","Proton"};
  const Int_t Color[nPart] = {(Int_t)kBlack, (Int_t)kRed, (Int_t)kBlue};
  const Int_t nPt = 6;
  const Float_t Pt[nPt] = {0.6, 1.0, 1.4, 1.8, 2.2, 2.6};

  // open file and get histo manager
  TFile * fin = TFile::Open("output/No_Outliers/AnResDATATrain12_NoOutliers.root");
  TDirectoryFile * dir = (TDirectoryFile*)fin->Get("OutputAODSpectraTask_Data_Cent30to40_QVec0.0to100.0_Eta-0.8to0.8_3.0SigmaPID_TrBit1024");
  AliSpectraAODHistoManager * hman = (AliSpectraAODHistoManager*)dir->Get("SpectraHistos");

  // get nsigma projections
  TH2F * hTPCTOFnsigmaDATA[nPart];
  TH1F * hTPCTOFnsigmaDATAProj[nPart][nPt];
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      // get 2d histo
      hTPCTOFnsigmaDATA[ipart] = (TH2F*)((TH2F*)hman->GetNSigHistogram(Form("hHistNSig%sPtTPCTOF",Particle[ipart].Data())))->Clone();
      for (Int_t ipt=0; ipt<nPt; ipt++)
	{
	  // project it for the given pt bin
	  hTPCTOFnsigmaDATAProj[ipart][ipt] = (TH1F*)hTPCTOFnsigmaDATA[ipart]->ProjectionY(Form("%ss - p_{T} = %.1f GeV/c",Particle[ipart].Data(),Pt[ipt]),
											   hTPCTOFnsigmaDATA[ipart]->GetXaxis()->FindBin(Pt[ipt]),
											   hTPCTOFnsigmaDATA[ipart]->GetXaxis()->FindBin(Pt[ipt]));
	}
    }

  //-----------------------------------------------
  // Create TPC+TOF nsigma distribution histogram -
  //-----------------------------------------------
  // file for saving histogram
  TFile * fout = new TFile("TPCTOFnsigma.root","RECREATE");

  // TPC nsigma distributions - just  gaussians
  TF1 * fTPCnsigma[nPart];
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      fTPCnsigma[ipart] = new TF1("fTPCnsigma","gaus",-25,25);
      fTPCnsigma[ipart]->SetParameter(0,1); // normalization
    }
  fTPCnsigma[0]->SetParameter(1,0.0); // pion mean
  fTPCnsigma[0]->SetParameter(2,0.98); // pion sigma
  fTPCnsigma[1]->SetParameter(1,-0.1); // kaon mean
  fTPCnsigma[1]->SetParameter(2,0.98); // kaon sigma
  fTPCnsigma[2]->SetParameter(1,0.2); // proton mean
  fTPCnsigma[2]->SetParameter(2,0.72); // proton sigma

  // TOF nsigma distributions
  gROOT->LoadMacro("/opt/alice/aliroot/trunk/src/PWGLF/SPECTRA/PiKaPr/TOF/PbPb276/macros/TOFsignal.C");
  TF1 *fTOFnsigma[nPart];
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      fTOFnsigma[ipart] = new TF1(Form("fTOFnsigma_%s",Particle[ipart].Data()), TOFsignal, -25, 25, 4);
      fTOFnsigma[ipart]->SetParameter(0,1.0); // normalization
      fTOFnsigma[ipart]->SetParameter(3, 1.0); // tail
    }
  fTOFnsigma[0]->SetParameter(1, 0.0); // pion mean
  fTOFnsigma[0]->SetParameter(2, 1.02); // pion sigma
  fTOFnsigma[1]->SetParameter(1, -0.04); // kaon mean
  fTOFnsigma[1]->SetParameter(2, 1.05); // kaon sigma
  fTOFnsigma[2]->SetParameter(1, -0.18); // proton mean
  fTOFnsigma[2]->SetParameter(2, 0.82); // proton sigma

  
  // combined TPC+TOF nsigma distributions
  //  // to ensure that the binning is the same I'm copying an actual nsigma projection and emptying it before doing anything
  TH1F * hTPCTOFnsigma[nPart];
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      //hTPCTOFnsigma[ipart] = (TH1F*)hTPCTOFnsigmaDATAProj[0][0]->Clone();
      //hTPCTOFnsigma[ipart]->Reset();
      hTPCTOFnsigma[ipart] = new TH1F("","",1000,0,25);
      hTPCTOFnsigma[ipart]->SetName(Form("hTPCTOFnsigma%s",Particle[ipart].Data()));
      hTPCTOFnsigma[ipart]->SetTitle("Predicted NSigma distributions (TPC+TOF);n#sigma;");
      
      for (Int_t i=0; i<1000000; i++)
	{
	  // get random numbers from the nsig distributions
	  Float_t TPCpoint = fTPCnsigma[ipart]->GetRandom(-25,25);
	  Float_t TOFpoint = fTOFnsigma[ipart]->GetRandom(-25,25);
	  
	  // add them in quadrature and divide by root 2
	  Float_t TPCTOFcomb = TMath::Sqrt( (TMath::Power(TPCpoint,2) + TMath::Power(TOFpoint,2))/2 );
	  
	  // fill hTPCTOFnsigma w/ this
	  hTPCTOFnsigma[ipart]->Fill(TPCTOFcomb);
	}
      
      // normalize it
      Float_t integral = hTPCTOFnsigma[ipart]->Integral(hTPCTOFnsigma[ipart]->FindBin(0),hTPCTOFnsigma[ipart]->FindBin(20));
      hTPCTOFnsigma[ipart]->Scale(1./integral);
    }
  
  // draw  them
  TCanvas * cTPCTOFnsigmaPredicted = new TCanvas("cTPCTOFnsigmaPredicted","cTPCTOFnsigmaPredicted");
  gPad->SetLogy();
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      hTPCTOFnsigma[ipart]->GetXaxis()->SetTitle("n#sigma");
      hTPCTOFnsigma[ipart]->GetXaxis()->SetRangeUser(-5,20);
      hTPCTOFnsigma[ipart]->SetLineColor(Color[ipart]);
      if (ipart == 0)
	{
	  hTPCTOFnsigma[ipart]->DrawCopy("hist");
	  TLegend * lTPCTOFnsigma = new TLegend(.69,.69,.99,.99);
	  lTPCTOFnsigma->SetFillColor(0);
	}     
      else hTPCTOFnsigma[ipart]->DrawCopy("histsame");
      lTPCTOFnsigma->AddEntry(hTPCTOFnsigma[ipart],Form("%ss",Particle[ipart].Data()),"l");
    }
  lTPCTOFnsigma->DrawClone();

  // do the same for MC
  // TPC nsigma distributions - just  gaussians
  TF1 * fTPCnsigmaMC[nPart];
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      fTPCnsigmaMC[ipart] = new TF1("fTPCnsigmaMC","gaus",-25,25);
      fTPCnsigmaMC[ipart]->SetParameter(0,1); // normalization
    }
  fTPCnsigmaMC[0]->SetParameter(1,-0.03); // pion mean
  fTPCnsigmaMC[0]->SetParameter(2,1.38); // pion sigma
  fTPCnsigmaMC[1]->SetParameter(1,0.14); // kaon mean
  fTPCnsigmaMC[1]->SetParameter(2,1.27); // kaon sigma
  fTPCnsigmaMC[2]->SetParameter(1,-0.17); // proton mean
  fTPCnsigmaMC[2]->SetParameter(2,1.03); // proton sigma

  // TOF nsigma distributions
  //gROOT->LoadMacro("/opt/alice/aliroot/trunk/src/PWGLF/SPECTRA/PiKaPr/TOF/PbPb276/macros/TOFsignal.C");
  TF1 *fTOFnsigmaMC[nPart];
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      fTOFnsigmaMC[ipart] = new TF1(Form("fTOFnsigmaMC_%s",Particle[ipart].Data()), TOFsignal, -25, 25, 4);
      fTOFnsigmaMC[ipart]->SetParameter(0,1.0); // normalization
      fTOFnsigmaMC[ipart]->SetParameter(3, 1.0); // tail
    }
  fTOFnsigmaMC[0]->SetParameter(1, -0.01); // pion mean
  fTOFnsigmaMC[0]->SetParameter(2, 0.94); // pion sigma
  fTOFnsigmaMC[1]->SetParameter(1, 0.02); // kaon mean
  fTOFnsigmaMC[1]->SetParameter(2, 0.84); // kaon sigma
  fTOFnsigmaMC[2]->SetParameter(1, 0.04); // proton mean
  fTOFnsigmaMC[2]->SetParameter(2, 0.66); // proton sigma

  
  // combined TPC+TOF nsigma distributions
  //  // to ensure that the binning is the same I'm copying an actual nsigma projection and emptying it before doing anything
  TH1F * hTPCTOFnsigmaMC[nPart];
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      //hTPCTOFnsigmaMC[ipart] = (TH1F*)hTPCTOFnsigmaMCDATAProj[0][0]->Clone();
      //hTPCTOFnsigmaMC[ipart]->Reset();
      hTPCTOFnsigmaMC[ipart] = new TH1F("","",1000,0,25);
      hTPCTOFnsigmaMC[ipart]->SetName(Form("hTPCTOFnsigmaMC%s",Particle[ipart].Data()));
      hTPCTOFnsigmaMC[ipart]->SetTitle("Predicted NSigma distributions (TPC+TOF);n#sigma;");
      
      for (Int_t i=0; i<1000000; i++)
	{
	  // get random numbers from the nsig distributions
	  Float_t TPCpoint = fTPCnsigmaMC[ipart]->GetRandom(-25,25);
	  Float_t TOFpoint = fTOFnsigmaMC[ipart]->GetRandom(-25,25);
	  
	  // add them in quadrature and divide by root 2
	  Float_t TPCTOFcomb = TMath::Sqrt( (TMath::Power(TPCpoint,2) + TMath::Power(TOFpoint,2))/2 );
	  
	  // fill hTPCTOFnsigmaMC w/ this
	  hTPCTOFnsigmaMC[ipart]->Fill(TPCTOFcomb);
	}
      
      // normalize it
      Float_t integral = hTPCTOFnsigmaMC[ipart]->Integral(hTPCTOFnsigmaMC[ipart]->FindBin(0),hTPCTOFnsigmaMC[ipart]->FindBin(20));
      hTPCTOFnsigmaMC[ipart]->Scale(1./integral);
    }
  
  // draw  them
  TCanvas * cTPCTOFnsigmaPredictedMC = new TCanvas("cTPCTOFnsigmaPredictedMC","cTPCTOFnsigmaPredictedMC");
  gPad->SetLogy();
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      hTPCTOFnsigmaMC[ipart]->GetXaxis()->SetTitle("n#sigma");
      hTPCTOFnsigmaMC[ipart]->GetXaxis()->SetRangeUser(-5,20);
      hTPCTOFnsigmaMC[ipart]->SetLineColor(Color[ipart]);
      if (ipart == 0)
	{
	  hTPCTOFnsigmaMC[ipart]->DrawCopy("hist");
	  TLegend * lTPCTOFnsigmaMC = new TLegend(.69,.69,.99,.99);
	  lTPCTOFnsigmaMC->SetFillColor(0);
	}     
      else hTPCTOFnsigmaMC[ipart]->DrawCopy("histsame");
      lTPCTOFnsigmaMC->AddEntry(hTPCTOFnsigmaMC[ipart],Form("%ss",Particle[ipart].Data()),"l");
    }
  lTPCTOFnsigmaMC->DrawClone();
  

  
  //---------------------------------------------------------------------------------------------------
  // Compare to the actual TPC+TOF nsigma distributions for different pt bins and different particles -
  //---------------------------------------------------------------------------------------------------
  TCanvas * cCompareToData[nPt];
  Int_t xPos = 50; // for canvas positioning
  Int_t yPos = 25; // for canvas positioning
  TH1F * hRatio[nPart][nPt];
  for (Int_t ipt=0; ipt<nPt; ipt++)
    {
      cCompareToData[ipt] = new TCanvas(Form("cCompareToData_%.1fPt",Pt[ipt]),Form("cCompareToData_%.1fPt",Pt[ipt]),xPos += 50, yPos += 25, 700, 500);
      cCompareToData[ipt]->Divide(3,2);
      
      // direct comparison of the histograms
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  cCompareToData[ipt]->cd(ipart+1);
	  gPad->SetLogy();

	  // draw histo from data
	  hTPCTOFnsigmaDATAProj[ipart][ipt]->SetLineColor(Color[ipart]);
	  hTPCTOFnsigmaDATAProj[ipart][ipt]->SetTitle(Form("TPC+TOF n#sigma for %ss at p_{T} = %.1f GeV/c;n#sigma;",Particle[ipart].Data(),Pt[ipt]));
	  hTPCTOFnsigmaDATAProj[ipart][ipt]->GetXaxis()->SetRangeUser(-5,10);
	  hTPCTOFnsigmaDATAProj[ipart][ipt]->DrawCopy("hist");

	  // scale generated histo and draw it as well
	  Float_t maxDATA = -1;

	  for (Int_t ibin=hTPCTOFnsigmaDATAProj[ipart][ipt]->FindBin(0); ibin<hTPCTOFnsigmaDATAProj[ipart][ipt]->FindBin(1.5); ibin++)
	    {
	      if (hTPCTOFnsigmaDATAProj[ipart][ipt]->GetBinContent(ibin) > maxDATA) maxDATA = hTPCTOFnsigmaDATAProj[ipart][ipt]->GetBinContent(ibin);
	    }
	  hTPCTOFnsigma[ipart]->Scale(maxDATA / hTPCTOFnsigma[ipart]->GetMaximum());
	  hTPCTOFnsigma[ipart]->SetLineColor(kGreen+2);
	  hTPCTOFnsigma[ipart]->DrawCopy("histsame");

	  // legend
	  TLegend * lCompareToData = new TLegend(.6,.8,.95,.92);
	  lCompareToData->SetFillColor(0);
	  lCompareToData->AddEntry(hTPCTOFnsigmaDATAProj[ipart][ipt],"Data","l");
	  lCompareToData->AddEntry(hTPCTOFnsigma[ipart],"Predicted","l");
	  lCompareToData->DrawClone();
	}

      // ratios of the histograms
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  cCompareToData[ipt]->cd(ipart+1+nPart);
	  
	  // compute the ratio
	  hRatio[ipart][ipt] = (TH1F*)hTPCTOFnsigmaDATAProj[ipart][ipt]->Clone();
	  hRatio[ipart][ipt]->Divide(hTPCTOFnsigma[ipart]);
	  hRatio[ipart][ipt]->SetName("Data/Predicted");
	  hRatio[ipart][ipt]->SetTitle("Data/Predicted;n#sigma;");
	  hRatio[ipart][ipt]->SetLineColor(Color[ipart]);
	  hRatio[ipart][ipt]->GetXaxis()->SetRangeUser(0,3);
	  //	  hRatio[ipart][ipt]->GetYaxis()->SetRangeUser(-2,10);
	  hRatio[ipart][ipt]->DrawCopy("hist");
	}
    }
  

  /*
  // crazy experiment: see if I can get a TF1 with the full convolution
  TF1 *fTOFnsig=new TF1("fTOFnsig","(x <= ([3]+[1])) * [0] *TMath::Gaus(x,[1],[2])+(x>([3]+[1]))*[0]*TMath::Gaus([3]+[1],[1],[2])*TMath::Exp(-([3])*(x-[3]-[1])/([2]*[2]))",0,25);
  fTOFnsig->SetParameter(0,1);
  fTOFnsig->SetParameter(1,1);
  fTOFnsig->SetParameter(2,1);
  fTOFnsig->SetParameter(3,3);

  TF1 * fTPCnsig = new TF1("fTPCnsig","gaus",0,25);
  fTPCnsig->SetParameter(0,1);
  fTPCnsig->SetParameter(1,0);
  fTPCnsig->SetParameter(2,1);

  TF1 * fTPCTOFnsig = new TF1("fTPCTOFnsig","TMath::Sqrt( (TMath::Power(fTPCnsigma,2) + TMath::Power(fTOFnsig,2))/2 )",0,25);

  //  TCanvas * c1 = new TCanvas();
  hTPCTOFnsigma->Fit(fTPCTOFnsig,"NM","",0,3);
  fTPCTOFnsig->DrawCopy("same");

  Float_t integralOfFit = fTPCTOFnsig->Integral(0,3);
  cout << "The integral is " << integralOfFit << endl;

  Float_t integralOfHist = hTPCTOFnsigma->Integral(hTPCTOFnsigma->FindBin(0), hTPCTOFnsigma->FindBin(3));
  cout << "The integral of the histo is " << integralOfHist << endl;
  */

  
  // print to a pdf for easy viewing
  cTPCTOFnsigmaPredicted->Print("TPCTOFnsigma.pdf(","pdf");
  for (Int_t ipt=0; ipt<nPt; ipt++)
    {
      if (ipt == nPt-1) cCompareToData[ipt]->Print("TPCTOFnsigma.pdf)","pdf");
      else cCompareToData[ipt]->Print("TPCTOFnsigma.pdf","pdf");
    }
  
  // save it all to file
  fout->cd();
  cTPCTOFnsigmaPredicted->Write();
  cTPCTOFnsigmaPredictedMC->Write();
  for (Int_t ipt=0; ipt<nPt; ipt++) cCompareToData[ipt]->Write();
  fout->Close();
}
