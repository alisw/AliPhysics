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

TDirectoryFile* GetDirectory(Int_t ievent, const TString& detName, Int_t nfiles) 
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

    TDirectoryFile* dir = (TDirectoryFile*)file0->Get(Form("Event%d",ievent));
    if (dir)  return dir;
  }
  return 0;
}  
    
void digitsTPC(Int_t nevents, Int_t nfiles)
{

  TH1F * hadclog = new TH1F ("hadclog", "hadclog",100, -1,3.5 );
  TH1F * hadc = new TH1F ("hadc", "hadc",200, -100,1100 );

  TH1F * hRow = new TH1F ("hRow", "hRow",120, -100, 1100);
  TH1F * hCol = new TH1F ("hCol", "hCol",100, -5, 194);
  TH1F * hndig = new TH1F ("hndig", "hndig",120, -100000, 100000);
  TH1F * hSize = new TH1F ("hSize", "hSize",100, 1000, 5000);


  TDirectoryFile *tdf[100];     
  TDirectoryFile *tdfKine[100] ;

  TTree *ttree[100];      
  TTree *ttreeKine[100];  

  AliSimDigits *digits= NULL;  


  Int_t numberhits=0;


  //Run loader------------
  TString name;
  name = "galice.root";
  AliRunLoader* rlSig = AliRunLoader::Open(name.Data());

  // gAlice
  rlSig->LoadgAlice();
  gAlice = rlSig->GetAliRun();

  // Now load kinematics and event header
  rlSig->LoadKinematics();
  rlSig->LoadHeader();
  cout <<  rlSig->GetNumberOfEvents()<< endl;
  //----------------------

    //loop over events in the files
    for(Int_t event=0; event<nevents; event++){
      printf("###event= %d\n", event);

    tdf[event] = GetDirectory(event, "TPC", nfiles);
    if ( ! tdf[event] ) {
      cerr << "Event directory not found in " << nfiles <<  " files" << endl;
      exit(1);
    }      

      ttree[event] = (TTree*)tdf[event]->Get("TreeD");
          
      digits = NULL;
      ttree[event]->SetBranchAddress("Segment", &digits);

      // Runloader -> gives particle Stack
      rlSig->GetEvent(event);
      AliStack * stack = rlSig->Stack(); 
      //stack->DumpPStack();

      Short_t digitValue=0;
      Int_t iRow=0;
      Int_t iColumn=0;
      Short_t ndig =0;
      Int_t digSize =0;

      // loop over tracks
      for(Int_t iev=0; iev<ttree[event]->GetEntries(); iev++){
	ttree[event]->GetEntry(iev);
	digits->First();    
	
	iRow=digits->CurrentRow();
	digitValue = digits->CurrentDigit();
	iColumn=digits->CurrentColumn();
	ndig=digits->GetDigits();
	digSize=digits->GetDigitSize();

	if(digitValue>0.)hadclog->Fill(TMath::Log10(digitValue));
	hadc->Fill(digitValue);
	hRow->Fill(iRow);
	hCol->Fill(iColumn);
	hndig->Fill(ndig);
	hSize->Fill(digSize);

	while (digits->Next()){

	  digitValue = digits->CurrentDigit();
	  iRow=digits->CurrentRow();
	  iColumn=digits->CurrentColumn();
	  ndig=digits->GetDigits();
	  digSize=digits->GetDigitSize();
	  
	  if(digitValue>0.)hadclog->Fill(TMath::Log10(digitValue));
	  hadc->Fill(digitValue);
	  hRow->Fill(iRow);
	  hCol->Fill(iColumn);
	  hndig->Fill(ndig);
	  hSize->Fill(digSize);
	  //cout << digSize << endl;
	}
      }
    }
  
   TFile fc("digits.TPC.root","RECREATE");
   fc.cd();

   hadclog->Write();
   hadc->Write();
   hRow->Write();
   hCol->Write();
   hndig->Write();
   hSize->Write();

   fc.Close();

}
