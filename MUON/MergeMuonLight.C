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
// This macro merges several files built with DecodeRecoCocktail.C into 
// a single one. 
// Arguments:     foutname = name of the output file
//                flistname = name of a text file containing the list of files
//                            to be merged 
void MergeMuonLight(char *foutname="MuonLightMerged.root",char *flistname="lista.lis"){ 
  // up to 2000 input files 

  TFile *file[2000];
  TClonesArray *muonArray   = new TClonesArray("AliMUONTrackLight",2); 
  TClonesArray *dimuonArray = new TClonesArray("AliMUONPairLight",1); 

  // open output file 
  TFile *fout = new TFile (foutname,"recreate"); 
  TTree *treeOut = new TTree("tree","tree"); 
  treeOut->Branch("muons",&muonArray); 
  treeOut->Branch("dimuons",&dimuonArray); 
  
  char *filename = new char[100];
  // read the list of input files 
  FILE *pf= fopen(flistname,"r");
  Int_t nfiles=0; 
  Int_t outflag=0; 
  // open the n input files 
  while (!outflag) { 
    if (fscanf(pf,"%s",filename)==1) { 
      file[nfiles++] = new TFile (filename); 
      cout << "Opening for input " << filename << endl;
    }
    else outflag = 1; 
  }
  fclose(pf); 

 
  for (Int_t ifile=0; ifile<nfiles; ifile++) {
    printf ("Scanning file %d: %s\n",ifile, file[ifile]->GetName()); 
    TTree *tree = (TTree*) file[ifile]->Get("tree"); 
    tree->SetBranchAddress("muons",&muonArray); 
    tree->SetBranchAddress("dimuons",&dimuonArray); 
    Int_t nev = tree->GetEntriesFast(); 
    printf ("Scanning file %d: %s containing %d events\n",
	    ifile, file[ifile]->GetName(),nev); 
    for (Int_t iev=0; iev<nev; iev++) {
      tree->GetEvent(iev); 
      treeOut->Fill();
    }
  }    
  fout->cd(); 
  treeOut->Write(); 
}
