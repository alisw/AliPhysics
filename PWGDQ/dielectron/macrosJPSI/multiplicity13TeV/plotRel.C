void plotRel(TString filename = "jpsiMany.root"){

  TFile *f =  TFile::Open( filename.Data());

  TH1D* hMult = (TH1D*)f->Get("h1MCEvents")->Clone();
  TH1D* hJpsi = (TH1D*)f->Get("histJpsiVsMult")->Clone();


  hMult->GetXaxis()->SetRange( 2,  hMult->GetXaxis()->GetNbins()  );
  hJpsi->GetXaxis()->SetRange( 2,  hJpsi->GetXaxis()->GetNbins()  );
  Double_t meanMult = hMult->GetMean();
  Double_t meanYield = hJpsi->Integral() / hMult->Integral();
  if( !meanMult || !meanYield ) return;
  const Int_t nBins = 19;
   Double_t newBins[nBins] = {
    -0.5, 2.5, 6.5, 10.5, 14.5,
     18.5, 22.5, 26.5, 30.5, 34.5,
     38.5, 42.5, 50.5, 58.5, 66.5,
    74.5, 82.5 , 98.5, 120.5 };


Double_t mean[nBins];

for(int i=0; i< nBins; ++i) {
  if(i < nBins -1 ){
    hMult->GetXaxis()->SetRangeUser( newBins[i], newBins[i+1] );
    mean[i] = hMult->GetMean();
    cout << "bin " << i << "-> "<< newBins[i] << " - " << newBins[i+1] <<  "    mean: " << mean[i] <<endl;
  }
}


// hJpsi->Draw();
  hMult = (TH1D*) hMult->Rebin( nBins-1, "multNewBins", newBins );
  hJpsi = (TH1D*) hJpsi->Rebin( nBins-1, "jpsiNewBins", newBins );

  TGraphErrors* hRel = new TGraphErrors(  );
  hRel->SetName("Errors");

  for(Int_t i = 0; i< nBins-1; ++i){
    if( !hMult->GetBinContent(i+1) ) return;
  }

  TGraph * hRel2 = (TGraph *) hRel->Clone();
  
  for(Int_t i = 0; i< nBins-1; ++i){
      Double_t yield = hJpsi->GetBinContent(i+1) / hMult->GetBinContent(i+1) / meanYield;
      Double_t yieldError = hJpsi->GetBinError(i+1) / hMult->GetBinContent(i+1) / meanYield;
      hRel->SetPoint( i, mean[i] / meanMult ,yield  );
      hRel2->SetPoint( i, mean[i] / meanMult ,yield  );
      hRel->SetPointError( i, 0. ,  yieldError   );
      cout <<  (mean[i] / meanMult) <<"  0.1  0.1  "  <<  yield<<"  "  << yieldError << "  " << yieldError <<endl;
   }


  hRel->SetFillColor(kBlue-4);
  hRel2->SetName("Line");
  hRel2->SetLineColor(kBlack);
  hRel->Draw("aP");

  TFile out( Form("relYield_%s", filename.Data() ), "RECREATE");
    hRel->Write();
    hRel2->Write();
  out.Close();
  
}
