///////////////////////////////////////////////////////////////////////
// CalculateTPCTOFPIDCorrFact.C  (called by JCGAnalysis.C)           //
//                                                                   //
// the return value is an array containing the correction factors    //
//                                                                   //
// written by John Groh                                              //
///////////////////////////////////////////////////////////////////////


Float_t * CalculateTPCTOFPIDCorrFact(TH2F ** hTPCTOFnsig)
{
  cout << "\n--- Just entered CalculateTPCTOFPIDCorrFact() ---\n";

  // return value
  Float_t * TPCTOFPIDCorrFact = new Float_t[nPart];

  const Int_t fixedBin[nPart] = {25,25,36}; // pt bin to use for the projection
  const Int_t fitMax[nPart] = {153,130,126}; // ranges for integration
  
  TH1F * hTPCTOFnsigProj[nPart]; // projection histograms

  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      // project the nsigma distribution in the desired pt bin
      hTPCTOFnsigProj[ipart] = (TH1F*)hTPCTOFnsig[ipart]->ProjectionY("",fixedBin[ipart],fixedBin[ipart])->Clone();

      // integrate and calculate the correction factors
      TPCTOFPIDCorrFact[ipart] = hTPCTOFnsigProj[ipart]->Integral(hTPCTOFnsigProj[ipart]->FindBin(0), hTPCTOFnsigProj[ipart]->FindBin(3));
      TPCTOFPIDCorrFact[ipart] /= hTPCTOFnsigProj[ipart]->Integral(hTPCTOFnsigProj[ipart]->FindBin(0), fitMax[ipart]);
    }

  // draw everything to make sure it worked
  TCanvas * cTPCTOFnsigInterp = new TCanvas("cTPCTOFnsigInterp","cTPCTOFnsigInterp");
  cTPCTOFnsigInterp->Divide(3,1);
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      cTPCTOFnsigInterp->cd(ipart+1);
      hTPCTOFnsigProj[ipart]->SetTitle(Form("TPC+TOF NSigma Projection - %ss - Bin %i;n#sigma (TPC+TOF);",Particle[ipart].Data(),fixedBin[ipart]));
      hTPCTOFnsigProj[ipart]->SetLineColor(Color[ipart]);
      hTPCTOFnsigProj[ipart]->DrawCopy("hist");      
    }
  
  cout << "\n\n\n\n\n--- PID CORRECTION FACTORS:\n"
       << "Pions: " << TPCTOFPIDCorrFact[0] << endl
       << "Kaons: " << TPCTOFPIDCorrFact[1] << endl
       << "Protons: " << TPCTOFPIDCorrFact[2] << endl << endl << endl;


  cout << "\n--- Leaving CalculateTPCTOFPIDCorrFact() ---\n";

  // return the correction factors
  return TPCTOFPIDCorrFact;
}




