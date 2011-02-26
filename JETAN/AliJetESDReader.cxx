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
//          Magali Estienne <magali.estienne@subatech.in2p3.fr>
//
// **February 2011
// implemented standard geometry (AliEMCALGeoUtils) instead of dummy one (AliJetDummyGeo)
// moved geometry definition in AliJetReader
// marco.bregant@subatech.in2p3.fr
//------------------------------------------------------------------------- 

// --- Standard library ---
#include <Riostream.h>

// --- ROOT system ---
#include <TSystem.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TTask.h"
#include "TTree.h"
#include "TFile.h"
#include <TGeoManager.h>
#include <assert.h>
#include <TRefArray.h>
#include <TMath.h>
#include <TProcessID.h>
#include <TRandom3.h>

// --- AliRoot header files ---
#include "AliJetESDReader.h"
#include "AliJetESDReaderHeader.h"
#include "AliESDEvent.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliEMCALGeoUtils.h"
#include "AliEMCALEMCGeometry.h"
#include "AliJetESDFillUnitArrayTracks.h"
#include "AliJetESDFillUnitArrayEMCalDigits.h"
#include "AliJetHadronCorrectionv1.h"
#include "AliJetUnitArray.h"
#include "AliAnalysisTask.h"
#include "AliOADBContainer.h"

ClassImp(AliJetESDReader)

AliJetESDReader::AliJetESDReader():
  AliJetReader(), 
  fHadCorr(0x0),
  fTpcGrid(0x0),
  fEmcalGrid(0x0),
  fGrid0(0),
  fGrid1(0),
  fGrid2(0),
  fGrid3(0),
  fGrid4(0),
  fApplyElectronCorrection(kFALSE),
  fApplyMIPCorrection(kTRUE),
  fApplyFractionHadronicCorrection(kFALSE),
  fFractionHadronicCorrection(0.3),
  fNumUnits(0),
  fDebug(0),
  fMass(0),
  fSign(0),
  fNIn(0),
  fOpt(0),
  fDZ(0),
  fNeta(0),
  fNphi(0),
  fRefArray(0x0),
  fProcId(kFALSE)
{
  // Constructor 
}

//____________________________________________________________________________
AliJetESDReader::~AliJetESDReader()
{
  // Destructor
    delete fChain;
    delete fTree;
    delete fESD;
    delete fTpcGrid;
    delete fEmcalGrid;
    if(fDZ)
      {
	delete fGrid0;
	delete fGrid1;
	delete fGrid2;
	delete fGrid3;
	delete fGrid4;
      }
}

//____________________________________________________________________________
void AliJetESDReader::OpenInputFiles()
{
  // chain for the ESDs
  fChain   = new TChain("esdTree");

  // get directory and pattern name from the header
   const char* dirName=fReaderHeader->GetDirectory();
   const char* pattern=fReaderHeader->GetPattern();
    
   // Add files matching patters to the chain
   void *dir  = gSystem->OpenDirectory(dirName);
   const char *name = 0x0;
   int nesd = ((AliJetESDReaderHeader*) fReaderHeader)->GetNesd();
   int a = 0;
   while ((name = gSystem->GetDirEntry(dir))){
       if (a>=nesd) continue;
       
       if (strstr(name,pattern)){
	 printf("Adding %s\n",name);
	 fChain->AddFile(Form("%s/%s/AliESDs.root",dirName,name));
	 a++;
       }
   }
  
  gSystem->FreeDirectory(dir);
  
  fESD = new AliESDEvent();
  fESD->ReadFromTree(fChain);

  int nMax = fChain->GetEntries(); 

  printf("\n AliJetESDReader: Total number of events in chain= %d \n",nMax);

  // set number of events in header
  if (fReaderHeader->GetLastEvent() == -1)
    fReaderHeader->SetLastEvent(nMax);
  else {
    Int_t nUsr = fReaderHeader->GetLastEvent();
    fReaderHeader->SetLastEvent(TMath::Min(nMax,nUsr));
  }

}

//__________________________________________________________
void AliJetESDReader::SetApplyMIPCorrection(Bool_t val)
{
  //
  // Set flag to apply MIP correction fApplyMIPCorrection
  // - exclusive with fApplyFractionHadronicCorrection
  //

  fApplyMIPCorrection = val;
  if(fApplyMIPCorrection == kTRUE)
    {
      SetApplyFractionHadronicCorrection(kFALSE);
      printf("Enabling MIP Correction \n");
    }
  else 
    {
      printf("Disabling MIP Correction \n");
    }
}

//__________________________________________________________
void AliJetESDReader::SetApplyFractionHadronicCorrection(Bool_t val)
{
  //
  // Set flag to apply EMC hadronic correction fApplyFractionHadronicCorrection
  // - exclusive with fApplyMIPCorrection
  //

  fApplyFractionHadronicCorrection = val;
  if(fApplyFractionHadronicCorrection == kTRUE)
    {
      SetApplyMIPCorrection(kFALSE);
      printf("Enabling Fraction Hadronic Correction \n");
    }
  else 
    {
      printf("Disabling Fraction Hadronic Correction \n");
    }
}

//__________________________________________________________
void AliJetESDReader::SetFractionHadronicCorrection(Double_t val)
{
  //
  // Set value to fFractionHadronicCorrection (default is 0.3)
  // apply EMC hadronic correction fApplyFractionHadronicCorrection
  // - exclusive with fApplyMIPCorrection
  //

  fFractionHadronicCorrection = val;
  if(fFractionHadronicCorrection > 0.0 && fFractionHadronicCorrection <= 1.0)
    {
      SetApplyFractionHadronicCorrection(kTRUE);
      printf("Fraction Hadronic Correction %1.3f \n",fFractionHadronicCorrection);
    }
  else 
    {
      SetApplyFractionHadronicCorrection(kFALSE);
    }
}


//____________________________________________________________________________
void AliJetESDReader::SetInputEvent(const TObject* esd, const TObject* /*aod*/, const TObject* /*mc*/) {
    // Connect the tree
     fESD   = (AliESDEvent*) esd;
}

//____________________________________________________________________________
Bool_t AliJetESDReader::FillMomentumArray()
{
  // Fill momentum array

  Int_t goodTrack = 0;
  Int_t nt = 0;
  Float_t pt, eta;
  TVector3 p3;
  
  // clear momentum array
  ClearArray();
  fDebug = fReaderHeader->GetDebug();
  fOpt = fReaderHeader->GetDetector();

  if (!fESD) {
      return kFALSE;
  }

  // get number of tracks in event (for the loop)
  nt = fESD->GetNumberOfTracks();
  printf("Fill Momentum Array %5d  \n", nt);
 
  // temporary storage of signal and pt cut flag
  Int_t* sflag  = new Int_t[nt];
  Int_t* cflag  = new Int_t[nt];
  
  // get cuts set by user
  Float_t ptMin =  fReaderHeader->GetPtCut();
  Float_t etaMin = fReaderHeader->GetFiducialEtaMin();
  Float_t etaMax = fReaderHeader->GetFiducialEtaMax();  
  
  //loop over tracks in ESD
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

  delete[] sflag;
  delete[] cflag;

  return kTRUE;

}

//____________________________________________________________________________
void AliJetESDReader::CreateTasks(TChain* tree)
{
  //
  // For reader task initialization
  //

  fDebug = fReaderHeader->GetDebug();
  fDZ = fReaderHeader->GetDZ();
  fTree = tree;

  // Init EMCAL geometry, if needed
   if(fGeom == 0)
	  SetEMCALGeometry();
	else Info(" SetEMCALGeometry","was already done.. it's called just once !!");

  // Init parameters
  InitParameters();
  // Create and init unit array
  InitUnitArray();

  // Create global reader task for analysis 
  fFillUnitArray = new TTask("fFillUnitArray","Fill unit array jet finder");
  // Create a task for to fill the charged particle information 
  fFillUAFromTracks = new AliJetESDFillUnitArrayTracks(); 
  fFillUAFromTracks->SetReaderHeader(fReaderHeader);
  fFillUAFromTracks->SetGeom(fGeom);
  fFillUAFromTracks->SetTPCGrid(fTpcGrid);
  fFillUAFromTracks->SetEMCalGrid(fEmcalGrid);
  if(fDZ)
    { // Calo dead zones inclusion
      fFillUAFromTracks->SetGrid0(fGrid0);
      fFillUAFromTracks->SetGrid1(fGrid1);
      fFillUAFromTracks->SetGrid2(fGrid2);
      fFillUAFromTracks->SetGrid3(fGrid3);
      fFillUAFromTracks->SetGrid4(fGrid4);
    }
  fFillUAFromTracks->SetApplyMIPCorrection(fApplyMIPCorrection);
  fFillUAFromTracks->SetHadCorrector(fHadCorr);
  // Create a task for to fill the neutral particle information 
  fFillUAFromEMCalDigits = new AliJetESDFillUnitArrayEMCalDigits();
  fFillUAFromEMCalDigits->SetReaderHeader(fReaderHeader);
  fFillUAFromEMCalDigits->SetGeom(fGeom);
  fFillUAFromEMCalDigits->SetTPCGrid(fTpcGrid);
  fFillUAFromEMCalDigits->SetEMCalGrid(fEmcalGrid);
  fFillUAFromEMCalDigits->SetApplyFractionHadronicCorrection(fApplyFractionHadronicCorrection);
  fFillUAFromEMCalDigits->SetFractionHadronicCorrection(fFractionHadronicCorrection);
   fFillUAFromEMCalDigits->SetApplyElectronCorrection(fApplyElectronCorrection);
  // Add the task to global task
  fFillUnitArray->Add(fFillUAFromTracks);
  fFillUnitArray->Add(fFillUAFromEMCalDigits);
  fFillUAFromTracks->SetActive(kFALSE);
  fFillUAFromEMCalDigits->SetActive(kFALSE);

  cout << "Tasks instantiated at that stage ! " << endl;
  cout << "You can loop over events now ! " << endl;
   
}

//____________________________________________________________________________
Bool_t AliJetESDReader::ExecTasks(Bool_t procid, TRefArray* refArray)
{
  //
  // Reader task execussion
  //

  fProcId = procid;
  fRefArray = refArray;

  // clear momentum array
  ClearArray();

  fDebug = fReaderHeader->GetDebug();
  fOpt = fReaderHeader->GetDetector();

  if(!fESD) {
    return kFALSE;
  }
  
  // TPC only or Digits+TPC or Clusters+TPC
  if(fOpt%2==!0 && fOpt!=0){ 
    fFillUAFromTracks->SetESD(fESD);
    fFillUAFromTracks->SetActive(kTRUE);
    fFillUAFromTracks->SetUnitArray(fUnitArray);
    fFillUAFromTracks->SetRefArray(fRefArray);
    fFillUAFromTracks->SetProcId(fProcId);
    //    fFillUAFromTracks->ExecuteTask("tpc"); // Temporarily changed
    fFillUAFromTracks->Exec("tpc");              // Temporarily added
    if(fOpt==1){
      fNumCandidate = fFillUAFromTracks->GetMult();
      fNumCandidateCut = fFillUAFromTracks->GetMultCut();
    }
  }

  // Digits only or Digits+TPC  
  if(fOpt>=2 && fOpt<=3){
    fFillUAFromEMCalDigits->SetESD(fESD);
    fFillUAFromEMCalDigits->SetActive(kTRUE);
    fFillUAFromEMCalDigits->SetUnitArray(fUnitArray);
    fFillUAFromEMCalDigits->SetRefArray(fRefArray);
    fFillUAFromEMCalDigits->SetProcId(fFillUAFromTracks->GetProcId());
    fFillUAFromEMCalDigits->SetInitMult(fFillUAFromTracks->GetMult());
    fFillUAFromEMCalDigits->SetInitMultCut(fFillUAFromTracks->GetMultCut());
    fFillUAFromEMCalDigits->Exec("digits");      // Temporarily added
    fNumCandidate = fFillUAFromEMCalDigits->GetMult();
    fNumCandidateCut = fFillUAFromEMCalDigits->GetMultCut();
  }

  //  fFillUnitArray->ExecuteTask();             // Temporarily commented

  return kTRUE;
}


//____________________________________________________________________________  
void AliJetESDReader::InitParameters()
{
  // Initialise parameters
  fOpt = fReaderHeader->GetDetector();
  //  fHadCorr        = 0;          // For hadron correction
  if(fEFlag==kFALSE){
    if(fOpt==0 || fOpt==1)  
      fECorrection    = 0;        // For electron correction
    else fECorrection = 1;        // For electron correction
  }
  fNumUnits       = fGeom->GetEMCGeometry()->GetNCells();      // Number of cells in EMCAL
  if(fDebug>1) printf("\n EMCal parameters initiated ! \n");
}

//____________________________________________________________________________
void AliJetESDReader::InitUnitArray()
{
  //
  // Create and Initialise unit arrays
  //

  Int_t nElements = fTpcGrid->GetNEntries();
  Float_t eta = 0., phi = 0., Deta = 0., Dphi = 0.;
  if(fArrayInitialised) fUnitArray->Delete();

  if(fTpcGrid->GetGridType()==0)
    { // Fill the following quantities :
      // Good track ID, (Eta,Phi) position ID, eta, phi, energy, px, py, pz, Deta, Dphi, 
      // detector flag, in/out jet, pt cut, mass, cluster ID)
      for(Int_t nBin = 1; nBin < nElements+1; nBin++)
	{
	  fTpcGrid->GetEtaPhiFromIndex2(nBin,phi,eta);
	  phi = ((phi < 0) ? phi + 2. * TMath::Pi() : phi);
	  Deta = fTpcGrid->GetDeta();
	  Dphi = fTpcGrid->GetDphi();
	  new ((*fUnitArray)[nBin-1]) AliJetUnitArray(nBin-1,0,eta,phi,0.,Deta,Dphi,kTpc,kOutJet,kPtSmaller,kPtSmaller,kBad,0.,-1);
	}
    }

  if(fTpcGrid->GetGridType()==1)
    {
      Int_t nGaps = 0;
      Int_t n0 = 0, n1 = 0, n2 = 0, n3 = 0, n4 = 0;

      if(fDZ)
	{
	  // Define a grid of cell for the gaps between SM
	  Double_t phimin0 = 0., phimin1 = 0., phimin2 = 0., phimin3 = 0., phimin4 = 0.;
	  Double_t phimax0 = 0., phimax1 = 0., phimax2 = 0., phimax3 = 0., phimax4 = 0.;
	  fGeom->GetEMCGeometry()->GetPhiBoundariesOfSMGap(0,phimin0,phimax0);
	  fGrid0 = new AliJetGrid(0,95,phimin0,phimax0,-0.7,0.7); // 0.015 x 0.015
	  fGrid0->SetGridType(0);
	  fGrid0->SetMatrixIndexes();
	  fGrid0->SetIndexIJ();
	  n0 = fGrid0->GetNEntries();
	  fGeom->GetEMCGeometry()->GetPhiBoundariesOfSMGap(1,phimin1,phimax1);
	  fGrid1 = new AliJetGrid(0,95,phimin1,phimax1,-0.7,0.7); // 0.015 x 0.015
	  fGrid1->SetGridType(0);
	  fGrid1->SetMatrixIndexes();
	  fGrid1->SetIndexIJ();
	  n1 = fGrid1->GetNEntries();
	  fGeom->GetEMCGeometry()->GetPhiBoundariesOfSMGap(2,phimin2,phimax2);
	  fGrid2 = new AliJetGrid(0,95,phimin2,phimax2,-0.7,0.7); // 0.015 x 0.015
	  fGrid2->SetGridType(0);
	  fGrid2->SetMatrixIndexes();
	  fGrid2->SetIndexIJ();
	  n2 = fGrid2->GetNEntries();
	  fGeom->GetEMCGeometry()->GetPhiBoundariesOfSMGap(3,phimin3,phimax3);
	  fGrid3 = new AliJetGrid(0,95,phimin3,phimax3,-0.7,0.7); // 0.015 x 0.015
	  fGrid3->SetGridType(0);  
	  fGrid3->SetMatrixIndexes();
	  fGrid3->SetIndexIJ();
	  n3 = fGrid3->GetNEntries();
	  fGeom->GetEMCGeometry()->GetPhiBoundariesOfSMGap(4,phimin4,phimax4);
	  fGrid4 = new AliJetGrid(0,95,phimin4,phimax4,-0.7,0.7); // 0.015 x 0.015
	  fGrid4->SetGridType(0);
	  fGrid4->SetMatrixIndexes();
	  fGrid4->SetIndexIJ();
	  n4 = fGrid4->GetNEntries();

	  nGaps = n0+n1+n2+n3+n4;

	}

      for(Int_t nBin = 0; nBin < fNumUnits+nElements+nGaps; nBin++) 
	{
	  if(nBin<fNumUnits)
	    {
	      fGeom->EtaPhiFromIndex(nBin, eta, phi); // From EMCal geometry 
	      // fEmcalGrid->GetEtaPhiFromIndex2(nBin,phi,eta); // My function from Grid
	      phi = ((phi < 0) ? phi + 2. * TMath::Pi() : phi);
	      Deta = fEmcalGrid->GetDeta(); // Modify with the exact detector values
	      Dphi = fEmcalGrid->GetDphi(); // Modify with the exact detector values
	      new ((*fUnitArray)[nBin]) AliJetUnitArray(nBin,0,eta,phi,0.,Deta,Dphi,kTpc,kOutJet,kPtSmaller,kPtSmaller,kBad,0.,-1);
	    } 
	  else {
	    if(nBin>=fNumUnits && nBin<fNumUnits+nElements){
	      fTpcGrid->GetEtaPhiFromIndex2(nBin+1-fNumUnits,phi,eta);
	      phi = ((phi < 0) ? phi + 2. * TMath::Pi() : phi);
	      Deta = fTpcGrid->GetDeta();
	      Dphi = fTpcGrid->GetDphi();
	      new ((*fUnitArray)[nBin]) AliJetUnitArray(nBin,0,eta,phi,0.,Deta,Dphi,kTpc,kOutJet,kPtSmaller,kPtSmaller,kBad,0.,-1);
	    }
	    else {
	      if(fDZ) {
		if(nBin>=fNumUnits+nElements && nBin<fNumUnits+nElements+nGaps){
		  if(nBin<fNumUnits+nElements+n0)
		    {
		      Float_t phi0 = eta = 0.;
		      fGrid0->GetEtaPhiFromIndex2(nBin+1-(fNumUnits+nElements),phi0,eta);
		      Deta = fGrid0->GetDeta(); 
		      Dphi = fGrid0->GetDphi(); 
		      new ((*fUnitArray)[nBin]) AliJetUnitArray(nBin,0,eta,phi0,0.,Deta,Dphi,kTpc,kOutJet,kPtSmaller,kPtSmaller,kBad,0.,-1);
		    }
		  else if(nBin>=fNumUnits+nElements+n0 && nBin<fNumUnits+nElements+n0+n1)
		    {
		      Float_t phi0 = eta = 0.;
		      fGrid1->GetEtaPhiFromIndex2(nBin+1-(fNumUnits+nElements+n0),phi0,eta);
		      Deta = fGrid1->GetDeta(); 
		      Dphi = fGrid1->GetDphi(); 
		      new ((*fUnitArray)[nBin]) AliJetUnitArray(nBin,0,eta,phi0,0.,Deta,Dphi,kTpc,kOutJet,kPtSmaller,kPtSmaller,kBad,0.,-1);
		    }
		  else if(nBin>=fNumUnits+nElements+n0+n1 && nBin<fNumUnits+nElements+n0+n1+n2)
		    {
		      Float_t phi0 = eta = 0.;
		      fGrid2->GetEtaPhiFromIndex2(nBin+1-(fNumUnits+nElements+n0+n1),phi0,eta);
		      Deta = fGrid2->GetDeta(); 
		      Dphi = fGrid2->GetDphi(); 
		      new ((*fUnitArray)[nBin]) AliJetUnitArray(nBin,0,eta,phi0,0.,Deta,Dphi,kTpc,kOutJet,kPtSmaller,kPtSmaller,kBad,0.,-1);
		    }
		  else if(nBin>=fNumUnits+nElements+n0+n1+n2 && nBin<fNumUnits+nElements+n0+n1+n2+n3)
		    {
		      Float_t phi0 = eta = 0.;
		      fGrid3->GetEtaPhiFromIndex2(nBin+1-(fNumUnits+nElements+n0+n1+n2),phi0,eta);
		      Deta = fGrid3->GetDeta(); 
		      Dphi = fGrid3->GetDphi(); 
		      new ((*fUnitArray)[nBin]) AliJetUnitArray(nBin,0,eta,phi0,0.,Deta,Dphi,kTpc,kOutJet,kPtSmaller,kPtSmaller,kBad,0.,-1);
		    }
		  else if(nBin>=fNumUnits+nElements+n0+n1+n2+n3 && nBin<fNumUnits+nElements+nGaps)
		    {
		      Float_t phi0 = eta = 0.;
		      fGrid4->GetEtaPhiFromIndex2(nBin+1-(fNumUnits+nElements+n0+n1+n2+n3),phi0,eta);
		      Deta = fGrid4->GetDeta(); 
		      Dphi = fGrid4->GetDphi(); 
		      new ((*fUnitArray)[nBin]) AliJetUnitArray(nBin,0,eta,phi0,0.,Deta,Dphi,kTpc,kOutJet,kPtSmaller,kPtSmaller,kBad,0.,-1);
		    }
		}
	      } // end if(fDZ)
	    } // end else 2
	  } // end else 1
	} // end loop on nBin
    } // end grid type == 1
  fArrayInitialised = 1;
}



