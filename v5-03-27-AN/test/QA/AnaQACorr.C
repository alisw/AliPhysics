/*
 *  AnaQACorr.C
 *  
 *
 *  Created by schutz on 15/08/08.
 *  Copyright 2008 CERN. All rights reserved.
 *
 */

void AnaQACorr(Int_t run)
{
  // Analyze the QA NTUPLE for correlation among detectors
  // Open the file and get the ntuple
  TFile fin(Form("CORR.QA.%d.0.root", run)) ; 
  TNtupleD * nt = dynamic_cast<TNtupleD*>(fin.FindObjectAny(AliQA::GetQACorrName())) ; 
  TObjArray * branches = nt->GetListOfBranches() ; 
  printf("The following parameters are available:\n") ; 
  for ( Int_t b = 0 ; b < branches->GetEntries() ; b++) {
    printf("  %2d ---> %s\n", b, branches->At(b)->GetName()) ; 
  }
}

