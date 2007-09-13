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
// Jet AOD Reader 
// AOD reader for jet analysis
// Author: Davide Perrino <davide.perrino@cern.ch>
//------------------------------------------------------------------------- 


#include <Riostream.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TChain.h>

#include "AliJetAODReader.h"
#include "AliJetAODReaderHeader.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"

ClassImp(AliJetAODReader)

AliJetAODReader::AliJetAODReader():
    AliJetReader(),
    fChain(0x0),
    fAOD(0x0),
    fRef(new TRefArray),
    fDebug(0),
    fOpt(0)
{
  // Constructor    
}

//____________________________________________________________________________

AliJetAODReader::~AliJetAODReader()
{
  // Destructor
    delete fChain;
    delete fAOD;
    delete fRef;
}

//____________________________________________________________________________

void AliJetAODReader::OpenInputFiles()
{
    // Open the necessary input files
    // chain for the AODs
  fChain   = new TChain("aodTree");

  // get directory and pattern name from the header
   const char* dirName=fReaderHeader->GetDirectory();
   const char* pattern=fReaderHeader->GetPattern();

//   // Add files matching patters to the chain

   void *dir  = gSystem->OpenDirectory(dirName);
   const char *name = 0x0;
   int naod = ((AliJetAODReaderHeader*) fReaderHeader)->GetNaod();
   int a = 0;
   while ((name = gSystem->GetDirEntry(dir))){
       if (a>=naod) continue;
       
       if (strstr(name,pattern)){
	   char path[256];
	   sprintf(path,"%s/%s/aod.root",dirName,name);
	   fChain->AddFile(path);
	   a++;
       }
   }
  
  gSystem->FreeDirectory(dir);
  

  fAOD = 0;
  fChain->SetBranchAddress("AOD",&fAOD);
  
  int nMax = fChain->GetEntries(); 

  printf("\n AliJetAODReader: Total number of events in chain= %d \n",nMax);
  
  // set number of events in header
  if (fReaderHeader->GetLastEvent() == -1)
    fReaderHeader->SetLastEvent(nMax);
  else {
    Int_t nUsr = fReaderHeader->GetLastEvent();
    fReaderHeader->SetLastEvent(TMath::Min(nMax,nUsr));
  }
}

//____________________________________________________________________________

void AliJetAODReader::ConnectTree(TTree* tree, TObject* /*data*/) {
    // Connect the tree
    // For AOD reader it's needed only to set the number of events
     fChain = (TChain*)      tree;
     
     Int_t nMax = fChain->GetEntries(); 
     printf("\n AliJetAODReader: Total number of events in chain= %5d \n", nMax);
     // set number of events in header
     if (fReaderHeader->GetLastEvent() == -1)
	 fReaderHeader->SetLastEvent(nMax);
     else {
	 Int_t nUsr = fReaderHeader->GetLastEvent();
	 fReaderHeader->SetLastEvent(TMath::Min(nMax,nUsr));
     }
}

//____________________________________________________________________________

Bool_t AliJetAODReader::FillMomentumArray(Int_t /*event*/)
{
  // Clear momentum array
  ClearArray();
  fDebug = fReaderHeader->GetDebug();
  
  if (!fAOD) {
      return kFALSE;
  }
  
  // get number of tracks in event (for the loop)
  Int_t nt = fAOD->GetNTracks();
  printf("Fill Momentum Array %5d ", nt);
  
  // temporary storage of signal and pt cut flag
  Int_t* sflag  = new Int_t[nt];
  Int_t* cflag  = new Int_t[nt];
  
  // get cuts set by user
  Float_t ptMin =  fReaderHeader->GetPtCut();
  Float_t etaMin = fReaderHeader->GetFiducialEtaMin();
  Float_t etaMax = fReaderHeader->GetFiducialEtaMax();  
  
  //loop over tracks
  Int_t aodTrack = 0;
  Float_t pt, eta;
  TVector3 p3;
  for (Int_t it = 0; it < nt; it++) {
      AliAODTrack *track = fAOD->GetTrack(it);

      Double_t mom[3] = {track->Px(),track->Py(),track->Pz()};
      p3.SetXYZ(mom[0],mom[1],mom[2]);
      pt = p3.Pt();
      eta = p3.Eta();
	  if ( pt < ptMin )                      continue;      // checking pt  cut
      if ( (eta > etaMax) || (eta < etaMin)) continue;      // checking eta cut
      sflag[aodTrack]=1;
      cflag[aodTrack]=1;
	  new ((*fMomentumArray)[aodTrack++]) TLorentzVector(p3,p3.Mag());
	  fRef->Add(track);
  }
  // set the signal flags
  fSignalFlag.Set(aodTrack,sflag);
  fCutFlag.Set(aodTrack,cflag);

  return kTRUE;
}
