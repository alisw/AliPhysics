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
// Authors: Mercedes Lopez Noriega (mercedes.lopez.noriega@cern.ch)
//          Magali Estienne <magali.estienne@IReS.in2p3.fr>
//------------------------------------------------------------------------- 


#include <Riostream.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TGeoManager.h>

#include "AliJetESDReader.h"
#include "AliJetESDReaderHeader.h"
#include "AliESD.h"
#include "AliESDtrack.h"
//#include "AliEMCALGeometry.h"
#include "AliJetDummyGeo.h"
#include "AliJetFillUnitArrayTracks.h"
#include "AliJetUnitArray.h"

ClassImp(AliJetESDReader)

AliJetESDReader::AliJetESDReader():
    AliJetReader(),  
    fGeom(0),
    fChain(0x0),
    fESD(0x0),
    fHadCorr(0x0),
    fTpcGrid(0x0),
    fEmcalGrid(0x0),
    fPtCut(0),
    fHCorrection(0),
    fNumUnits(0),
    fDebug(0),
    fNIn(0),
    fOpt(0),
    fNeta(0),
    fNphi(0),
    fArrayInitialised(0) 
{
  // Constructor    
}

//____________________________________________________________________________

AliJetESDReader::~AliJetESDReader()
{
  // Destructor
    delete fChain;
    delete fESD;
    delete fTpcGrid;
    delete fEmcalGrid;
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
   int nesd = ((AliJetESDReaderHeader*) fReaderHeader)->GetNesd();
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

  printf("\nTotal number of events in chain= %d \n",nMax);
  
  // set number of events in header
  if (fReaderHeader->GetLastEvent() == -1)
    fReaderHeader->SetLastEvent(nMax);
  else {
    Int_t nUsr = fReaderHeader->GetLastEvent();
    fReaderHeader->SetLastEvent(TMath::Min(nMax,nUsr));
  }
}

void AliJetESDReader::ConnectTree(TTree* tree) {
    // Connect the tree
     fChain = (TChain*) tree;
     
     fChain->SetBranchAddress("ESD",    &fESD);
     Int_t nMax = fChain->GetEntries(); 
     printf("\nTotal number of events in chain= %5d \n", nMax);
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
  fDebug = fReaderHeader->GetDebug();
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
      if ((status & AliESDtrack::kTPCrefit) == 0)    continue;      // quality check
      if ((status & AliESDtrack::kITSrefit) == 0)    continue;      // quality check
      if (((AliJetESDReaderHeader*) fReaderHeader)->ReadSignalOnly() 
	  && TMath::Abs(track->GetLabel()) > 10000)  continue;      // quality check
      if (((AliJetESDReaderHeader*) fReaderHeader)->ReadBkgdOnly() 
	  && TMath::Abs(track->GetLabel()) < 10000)  continue;      // quality check
      eta = p3.Eta();
      if ( (eta > etaMax) || (eta < etaMin))         continue;      // checking eta cut
      
      new ((*fMomentumArray)[goodTrack]) TLorentzVector(p3,p3.Mag());
      sflag[goodTrack]=0;
      if (TMath::Abs(track->GetLabel()) < 10000) sflag[goodTrack]=1;
      cflag[goodTrack]=0;
      if (pt > ptMin) cflag[goodTrack]=1;                           // pt cut
      goodTrack++;
  }
  // set the signal flags
  fSignalFlag.Set(goodTrack,sflag);
  fCutFlag.Set(goodTrack,cflag);

//
//
  if (fTpcGrid || fEmcalGrid) {
      SetEMCALGeometry();
      InitParameters();
      AliJetFillUnitArrayTracks *fillUAFromTracks = new AliJetFillUnitArrayTracks(); 
      fillUAFromTracks->SetReaderHeader(fReaderHeader);
      fillUAFromTracks->SetMomentumArray(fMomentumArray);
      fillUAFromTracks->SetTPCGrid(fTpcGrid);
      fillUAFromTracks->SetEMCalGrid(fEmcalGrid);
      fillUAFromTracks->SetHadCorrection(fHCorrection);
      fillUAFromTracks->SetHadCorrector(fHadCorr);
      fNeta = fillUAFromTracks->GetNeta();
      fNphi = fillUAFromTracks->GetNphi();
      fillUAFromTracks->SetActive(kFALSE);
      // TPC only or Digits+TPC or Clusters+TPC
      if(fOpt%2==!0 && fOpt!=0) { 
	  fillUAFromTracks->SetActive(kTRUE);
	  fillUAFromTracks->SetUnitArray(fUnitArray);
	  fillUAFromTracks->ExecuteTask("tpc");
      }
  
      delete fillUAFromTracks;
  }

  return kTRUE;
}


void AliJetESDReader::SetEMCALGeometry()
{
  // Define EMCAL geometry to be able to read ESDs
    fGeom = AliJetDummyGeo::GetInstance();
    if (fGeom == 0)
	fGeom = AliJetDummyGeo::GetInstance("SHISH_77_TRD1_2X2_FINAL_110DEG","EMCAL");
    
    // To be setted to run some AliEMCALGeometry functions
    TGeoManager::Import("geometry.root");
    fGeom->GetTransformationForSM();  
    printf("\n EMCal Geometry set ! \n");

}
  
void AliJetESDReader::InitParameters()
{
    // Initialise parameters
    fHCorrection    = 0;                 // For hadron correction
    fHadCorr        = 0;                 // For hadron correction
    fNumUnits       = fGeom->GetNCells();      // Number of cells in EMCAL
    if(fDebug>1) printf("\n EMCal parameters initiated ! \n");
}

void AliJetESDReader::InitUnitArray()
{
  //Initialises unit arrays
    Int_t nElements = fTpcGrid->GetNEntries();
    if(fArrayInitialised) delete [] fUnitArray;
    if(fTpcGrid->GetGridType()==0)
	fUnitArray = new AliJetUnitArray[nElements];
    if(fTpcGrid->GetGridType()==1)
	fUnitArray = new AliJetUnitArray[fNumUnits+nElements];
    fArrayInitialised = 1;
}


