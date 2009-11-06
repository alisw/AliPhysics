/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$

// Macro to generate histograms from digits
// By E. Sicking, CERN

TTree* GetTreeD(Int_t ievent, const TString& detName, Int_t nfiles) 
{
  for (Int_t file =0; file<nfiles; file++) {
    TString filename(detName);
    if ( file == 0 ) {
      filename += ".Digits.root";
    }  
    else { 
      filename += TString(Form(".Digits%d.root",file));
    }

    TFile* file0 = TFile::Open(filename.Data());

    TTree* treeD = (TTree*)file0->Get(Form("Event%d/TreeD",ievent));
    if (treeD)  return treeD;
  }
  return 0;
}  
    
void digitsTRD(Int_t nevents, Int_t nfiles)
{


  TH1F * hadc = new TH1F("hadc", "hadc", 260, -99, 1200);
  TH1F * hadcLow = new TH1F("hadcLow", "hadcLow", 100, -5, 5);

  TH1F * hadclog = new TH1F("hadclog", "hadclog", 100, -1, 4);

    TTree *treeD=0x0;
    AliTRDdigitsManager manD;

    for (Int_t event=0; event<nevents; event++) {
      cout << "Event " << event << endl;
  
      treeD = GetTreeD(event, "TRD", nfiles);
      if ( ! treeD ) {
        cerr << "Event directory not found in " << nfiles <<  " files" << endl;
        exit(1);
      }      
      manD.ReadDigits(treeD);
    
      AliTRDarrayADC *digitsD = 0;

    
      Int_t maxDet = 540;
      Int_t rowMax  = 0;
      Int_t colMax  = 0;
      Int_t timeMax = 0;
      Double_t signal=0;

      for (Int_t idet = 0; idet<maxDet; idet++)
	{
	  digitsD = manD.GetDigits(idet);
	  digitsD->Expand();
	
	  rowMax = digitsD->GetNrow();
	  colMax = digitsD->GetNcol();
	  timeMax = digitsD->GetNtime();	
	
	  //cout << "Detector " << idet << endl;
	  cout << "\r Detector " << idet; cout.flush();
	
	  for (Int_t irow = 0; irow < rowMax; irow++)
	    {
	      for (Int_t icol = 0; icol < colMax; icol++)
		{
		  for (Int_t itime = 0; itime < timeMax; itime++)
		    {
		    
		      signal= digitsD->GetData(irow, icol, itime);

		      hadc-> Fill(signal);
		      hadcLow-> Fill(signal);
		      if(signal>0.) hadclog-> Fill(TMath::Log10(signal));
		    
		    
		    
		    } //time
		} //col
	    } //row
	}//detector chamber
      cout << "\n \r Event " << event << endl;
    }//event
  
  TFile fc("digits.TRD.root","RECREATE");
  fc.cd();
  hadcLow->Write();
  hadc->Write();
  hadclog->Write();

  fc.Close();

}
