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

//------------------------------------------------------------------------- 
// Jet ESD Reader 
// ESD reader for jet analysis
// Author: Mercedes Lopez Noriega (mercedes.lopez.noriega@cern.ch)
//------------------------------------------------------------------------- 


#include <Riostream.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "AliJetESDReader.h"
#include "AliJetESDReaderHeader.h"
#include "AliESD.h"
#include "AliESDtrack.h"

ClassImp(AliJetESDReader)

AliJetESDReader::AliJetESDReader():
  fMass(0),
  fSign(0)
{
  // Constructor    
  fReaderHeader = 0x0;
}

//____________________________________________________________________________

AliJetESDReader::~AliJetESDReader()
{
  // Destructor
}

//____________________________________________________________________________

void AliJetESDReader::OpenInputFiles()
{
  // chain for the ESDs
  fChain   = new TChain("esdTree");

  // get directory and pattern name from the header
   const char* dirName=fReaderHeader->GetDirectory();
   const char* pattern=fReaderHeader->GetPattern();
    
//   // Add files matching patters to the chain
  
   void *dir  = gSystem->OpenDirectory(dirName);
   const char *name = 0x0;
   int nesd = fReaderHeader->GetNesd();
   int a = 0;
   while ((name = gSystem->GetDirEntry(dir))){
       if (a>=nesd) continue;
       
       if (strstr(name,pattern)){
	   char path[256];
	   sprintf(path,"%s/%s/AliESDs.root",dirName,name);
	   fChain->AddFile(path);
	   a++;
       }
   }
  
  gSystem->FreeDirectory(dir);
  

  fESD = 0;
  fChain->SetBranchAddress("ESD",    &fESD);
  
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

void AliJetESDReader::ConnectTree(TTree* tree) {
     fChain = (TChain*) tree;
     
     fChain->SetBranchAddress("ESD",    &fESD);
     Int_t nMax = fChain->GetEntries(); 
     printf("\nTotal number of events in chain= %5d", nMax);
     // set number of events in header
     if (fReaderHeader->GetLastEvent() == -1)
	 fReaderHeader->SetLastEvent(nMax);
     else {
	 Int_t nUsr = fReaderHeader->GetLastEvent();
	 fReaderHeader->SetLastEvent(TMath::Min(nMax,nUsr));
     }
}

//____________________________________________________________________________

Bool_t AliJetESDReader::FillMomentumArray(Int_t event)
{
  // Fill momentum array

  Int_t goodTrack = 0;
  Int_t nt = 0;
  Float_t pt, eta;
  TVector3 p3;
  
  // clear momentum array
   ClearArray();
  
  // get event from chain
  fChain->GetTree()->GetEntry(event);

  if (!fESD) {
      return kFALSE;
  }
  
  // get number of tracks in event (for the loop)
  nt = fESD->GetNumberOfTracks();
  // temporary storage of signal and pt cut flag
  Int_t* sflag  = new Int_t[nt];
  Int_t* cflag  = new Int_t[nt];
  
  // get cuts set by user
  Float_t ptMin =  fReaderHeader->GetPtCut();
  Float_t etaMin = fReaderHeader->GetFiducialEtaMin();
  Float_t etaMax = fReaderHeader->GetFiducialEtaMax();  
  
  //loop over tracks
  for (Int_t it = 0; it < nt; it++) {
      AliESDtrack *track = fESD->GetTrack(it);
      UInt_t status = track->GetStatus();
      
      Double_t mom[3];
      track->GetPxPyPz(mom);
      p3.SetXYZ(mom[0],mom[1],mom[2]);
      pt = p3.Pt();
      if ((status & AliESDtrack::kTPCrefit) == 0) continue;    // quality check
      if ((status & AliESDtrack::kITSrefit) == 0) continue;    // quality check
      if (((AliJetESDReaderHeader*) fReaderHeader)->ReadSignalOnly() 
	  && TMath::Abs(track->GetLabel()) > 10000)  continue;   // quality check
      if (((AliJetESDReaderHeader*) fReaderHeader)->ReadBkgdOnly() 
	  && TMath::Abs(track->GetLabel()) < 10000)  continue;   // quality check
      eta = p3.Eta();
      if ( (eta > etaMax) || (eta < etaMin)) continue;           // checking eta cut
      
      new ((*fMomentumArray)[goodTrack]) TLorentzVector(p3,p3.Mag());
      sflag[goodTrack]=0;
      if (TMath::Abs(track->GetLabel()) < 10000) sflag[goodTrack]=1;
      cflag[goodTrack]=0;
      if (pt > ptMin) cflag[goodTrack]=1;                       // pt cut
      goodTrack++;
  }
  // set the signal flags
  fSignalFlag.Set(goodTrack,sflag);
  fCutFlag.Set(goodTrack,cflag);
  return kTRUE;
}

  


