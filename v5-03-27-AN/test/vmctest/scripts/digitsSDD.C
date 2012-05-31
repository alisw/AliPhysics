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
    
void digitsSDD(Int_t nevents, Int_t nfiles)
{

  TH1F * hadc = new TH1F ("hadc", "hadc",200, 0, 1300);   
  TH1F * hadclog = new TH1F ("hadclog", "hadclog",200, 1, 3.5);   

  TDirectoryFile *tdf[100];     
  TDirectoryFile *tdfKine[100] ;

  TTree *ttree[100];     
  TTree *ttreeKine[100]; 

  TClonesArray *arr= NULL; // 
 
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

    tdf[event] = GetDirectory(event, "ITS", nfiles);
    if ( ! tdf[event] ) {
      cerr << "Event directory not found in " << nfiles <<  " files" << endl;
      exit(1);
    }      
    
    ttree[event] = (TTree*)tdf[event]->Get("TreeD");
          
    arr = NULL;
    ttree[event]->SetBranchAddress("ITSDigitsSDD", &arr);

    rlSig->GetEvent(event);
    AliStack * stack = rlSig->Stack(); 
   


    // loop over tracks
    Int_t NumberPrim=0;
    for(Int_t iev=0; iev<ttree[event]->GetEntries(); iev++){
      ttree[event]->GetEntry(iev);
  
	
      for (Int_t j = 0; j < arr->GetEntries(); j++) {
	AliITSdigit* digit = dynamic_cast<AliITSdigit*> (arr->At(j));
	if (digit){
	  Int_t label = digit->GetHits();
	 
	  hadc->Fill(digit->GetSignal());
	  hadclog->Fill(TMath::Log10(digit->GetSignal()));

	}
      }
    }

  }
  
  TFile fc("digits.ITS.SDD.root","RECREATE");
  fc.cd();
   
  hadc->Write();
  hadclog->Write();

  fc.Close();

}
