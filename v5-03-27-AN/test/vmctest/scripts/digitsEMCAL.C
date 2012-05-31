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
    
void digitsEMCAL(Int_t nevents, Int_t nfiles) {

  TH1F *hadc     = new TH1F("hadc","EMCAL digit",200, -10., 200.);
  TH1F *hadclog     = new TH1F("hadclog","EMCAL digit",200, 0.5, 3.);
  
    TTree *treeD=0x0;
  
    TClonesArray *digits =0x0;
  
    for (Int_t event=0; event<nevents; event++) {
      cout << "Event " << event << endl;

      treeD = GetTreeD(event, "EMCAL", nfiles);
      if ( ! treeD ) {
        cerr << "Event directory not found in " << nfiles <<  " files" << endl;
        exit(1);
      }      

      digits = NULL;
      treeD->SetBranchAddress("EMCAL", &digits);

      for(Int_t iev=0; iev<treeD->GetEntries(); iev++){
	treeD->GetEntry(iev);

      
	for (Int_t j = 0; j < digits->GetEntries(); j++) {
	  AliDigitNew* dig = dynamic_cast<AliDigitNew*> (digits->At(j));
	  hadc->Fill(dig->GetAmp());
	  if(dig->GetAmp()>0)hadclog->Fill(TMath::Log10(dig->GetAmp()));
	}

      }
    }

   TFile fc("digits.EMCAL.root","RECREATE");
   fc.cd();
   
   hadc->Write();
   hadclog->Write();

   fc.Close();
}







