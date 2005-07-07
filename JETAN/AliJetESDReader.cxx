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
 
//
// Jet ESD Reader 
// ESD reader for jet analysis
// Author: Mercedes Lopez Noriega 
// mercedes.lopez.noriega@cern.ch
//


#include <Riostream.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include "AliJetESDReader.h"
#include "AliJetESDReaderHeader.h"
#include "AliESD.h"
#include "AliESDtrack.h"

ClassImp(AliJetESDReader)

AliJetESDReader::AliJetESDReader()
{
  // Constructor    

  fReaderHeader = 0x0;
  fMass = 0;
  fSign = 0;
}

//____________________________________________________________________________

AliJetESDReader::~AliJetESDReader()
{
  // Destructor
  //  delete fReaderHeader;
}

//____________________________________________________________________________

void AliJetESDReader::OpenInputFiles()

{
  // chain for the ESDs
  fChain   = new TChain("esdTree");
  fChainMC = new TChain("mcStackTree");

  // get directory and pattern name from the header
  const char* dirName=fReaderHeader->GetDirectory();
  const char* pattern=fReaderHeader->GetPattern();
    
  // Add files matching patters to the chain
  void *dir  = gSystem->OpenDirectory(dirName);
  const char *name = 0x0;
  while ((name = gSystem->GetDirEntry(dir))){
    if (strstr(name,pattern)){
      printf("Adding %s\n",name);
      char path[256];
      sprintf(path,"%s/%s",dirName,name);
      fChain->AddFile(path,-1,"esdTree");
      fChainMC->AddFile(path,-1,"mcStackTree");
     }
  }

  gSystem ->FreeDirectory(dir);
  fChain  ->SetBranchAddress("ESD",    &fESD);
  fChainMC->SetBranchAddress("Header", &fAliHeader);
  fChainMC->SetBranchAddress("Stack",  &fArrayMC);

  int nMax = fChain->GetEntries(); 
  printf("\nTotal number of events in chain= %d",nMax);

  // set number of events in header
  if (fReaderHeader->GetLastEvent() == -1)
    fReaderHeader->SetLastEvent(nMax);
  else {
    Int_t nUsr = fReaderHeader->GetLastEvent();
    fReaderHeader->SetLastEvent(TMath::Min(nMax,nUsr));
  }
}

//____________________________________________________________________________

void AliJetESDReader::FillMomentumArray(Int_t event)
{
// Fills the momentum array
  Int_t goodTrack = 0;
  Int_t nt = 0;
  Float_t pt;
  Double_t p[3]; // track 3 momentum vector

  // clear momentum array
  ClearArray();
  // get event from chain
  fChain->GetEntry(event);
  fChainMC->GetEntry(event);
  // get number of tracks in event (for the loop)
  nt = fESD->GetNumberOfTracks();
  // tmporary storage of signal flags
  Int_t* flag  = new Int_t[nt];

  // get cuts set by user
  Float_t ptMin = ((AliJetESDReaderHeader*) fReaderHeader)->GetPtCut();

  //loop over tracks
  for (Int_t it = 0; it < nt; it++) {
      AliESDtrack *track = fESD->GetTrack(it);
      UInt_t status = track->GetStatus();
      if ((status & AliESDtrack::kITSrefit) == 0) continue; // quality check

      if (((AliJetESDReaderHeader*) fReaderHeader)->ReadSignalOnly() 
	  && TMath::Abs(track->GetLabel()) > 10000)  continue;
    
      track->GetPxPyPz(p);
      pt = TMath::Sqrt(p[0] * p[0] + p[1] * p[1]); // pt of the track

      if (pt < ptMin) continue; //check  cuts 

      new ((*fMomentumArray)[goodTrack]) 
	  TLorentzVector(p[0], p[1], p[2],
			 TMath::Sqrt(pt * pt +p[2] * p[2]));
      flag[goodTrack]=0;
      if (TMath::Abs(track->GetLabel()) < 10000) flag[goodTrack]=1;
      goodTrack++;
  }
  // set the signal flags
  fSignalFlag.Set(goodTrack,flag);

  printf("\nIn event %d, number of good tracks %d \n", event, goodTrack);
}


