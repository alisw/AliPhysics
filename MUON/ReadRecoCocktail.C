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

/* $Id$ */

// A. De Falco, H. Woehri, INFN Cagliari, July 2006 
// base macro to read the trees generated with DecodeRecoCocktail.C 
// 

void ReadRecoCocktail(char* fname="./MuonLight.root"){ 
  TFile *file = new TFile(fname); 
  TClonesArray *muonArray   = new TClonesArray("AliMUONTrackLight",100); 
  TClonesArray *dimuonArray = new TClonesArray("AliMUONPairLight",100); 
  TTree *tree = (TTree*) file->Get("tree"); 
  tree->SetBranchAddress("muons",&muonArray); 
  tree->SetBranchAddress("dimuons",&dimuonArray); 
  Int_t nev = tree->GetEntriesFast(); 
  printf ("%d events in tree\n",nev); 
  Int_t ndimuUncorr=0; 
  for (Int_t iev=0; iev<nev; iev++) { 
    tree->GetEvent(iev); 
    Int_t nmu = muonArray->GetEntries(); 
    for (Int_t imu=0; imu<nmu; imu++) { 
      AliMUONTrackLight *mu = (AliMUONTrackLight*) muonArray->At(imu); 
      //      mu->Dump(); 
    }
    Int_t ndimu = dimuonArray->GetEntriesFast(); 
    for (Int_t idimu=0; idimu<ndimu; idimu++) { 
      AliMUONPairLight *dimu = (AliMUONPairLight*) dimuonArray->At(idimu); 
      //      dimu->Dump(); 
      if (!dimu->IsCorrelated()) { 
	printf ("Event %d / %d; dimuon %d not correlated\n",iev,nev,idimu); 
	dimu->PrintInfo("H"); 
	ndimuUncorr++; 
      }
    }
  }
  printf ("%d uncorrelated dimuons in file\n",ndimuUncorr);  
  TCanvas *c1=new TCanvas(); 
  c1->Divide(2);
  printf ("plotting the mass and pt spectra of opposite sign dimuons\n");
  c1->cd(1);
  tree->Draw("dimuons.GetPRec().M()","dimuons.GetCharge()==0");
  c1->cd(2);
  tree->Draw("dimuons.GetPRec().Pt()","dimuons.GetCharge()==0");
}
