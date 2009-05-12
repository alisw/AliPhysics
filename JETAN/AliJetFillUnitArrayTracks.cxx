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


//======================================================================
// ***July 2006
// Fill Unit Array class 
// Class used by AliJetESDReader to fill a UnitArray from the information extracted 
// from the particle tracks
// Author: magali.estienne@ires.in2p3.fr
//======================================================================


// --- Standard library ---
#include <Riostream.h>

// --- ROOT system ---
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TRefArray.h> 
#include <TVector3.h>
#include "TTask.h"
#include <TGeoManager.h>
#include <TMatrixD.h>
#include <TArrayD.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TProcessID.h>

// --- AliRoot header files ---
#include "AliJetFinder.h"
#include "AliJetReaderHeader.h"
#include "AliJetReader.h"
#include "AliJetESDReader.h"
#include "AliJetESDReaderHeader.h"
#include "AliESDEvent.h"
#include "AliJetDummyGeo.h"
#include "AliJetUnitArray.h"
#include "AliJetFillUnitArrayTracks.h"
#include "AliJetHadronCorrectionv1.h"
#include "AliJetGrid.h"

ClassImp(AliJetFillUnitArrayTracks)

//_____________________________________________________________________________
AliJetFillUnitArrayTracks::AliJetFillUnitArrayTracks()
  : TTask("AliJetFillUnitArrayTracks","Fill Unit Array with tpc/its and emcal information"),
    fNumUnits(0),
    fEtaMinCal(0),
    fEtaMaxCal(0),
    fPhiMinCal(0),
    fPhiMaxCal(0),
    fHadCorr(0),
    fHCorrection(0),
    fNTracks(0),
    fNTracksCut(0),
    fOpt(0),
    fDZ(0),
    fDebug(0),
    fReaderHeader(0x0),
    fMomentumArray(0x0),
    fUnitArray(0x0),
    fRefArray(0x0),
    fProcId(kFALSE),
    fTPCGrid(0x0),
    fEMCalGrid(0x0),
    fGeom(0x0),
    fESD(0x0),
    fGrid0(0x0),
    fGrid1(0x0),
    fGrid2(0x0),
    fGrid3(0x0),
    fGrid4(0x0),
    fNphi(0),
    fNeta(0),
    fPhi2(0x0),
    fEta2(0x0),
    fPhi(0x0),
    fEta(0x0),
    fIndex(0x0),
    fParams(0x0),
    fGrid(0),
    fPhiMin(0),
    fPhiMax(0),
    fEtaMin(0),
    fEtaMax(0),
    fEtaBinInTPCAcc(0),
    fPhiBinInTPCAcc(0),
    fEtaBinInEMCalAcc(0),
    fPhiBinInEMCalAcc(0),
    fNbinPhi(0)
{
  // constructor
}

//_____________________________________________________________________________
AliJetFillUnitArrayTracks::AliJetFillUnitArrayTracks(AliESDEvent* esd)
  : TTask("AliJetFillUnitArrayTracks","Fill Unit Array with tpc/its and emcal information"),
    fNumUnits(0),
    fEtaMinCal(0),
    fEtaMaxCal(0),
    fPhiMinCal(0),
    fPhiMaxCal(0),
    fHadCorr(0),
    fHCorrection(0),
    fNTracks(0),
    fNTracksCut(0),
    fOpt(0),
    fDZ(0),
    fDebug(0),
    fReaderHeader(0x0),
    fMomentumArray(0x0),
    fUnitArray(0x0),
    fRefArray(0x0),
    fProcId(kFALSE),
    fTPCGrid(0x0),
    fEMCalGrid(0x0),
    fGeom(0x0),
    fESD(esd),
    fGrid0(0x0),
    fGrid1(0x0),
    fGrid2(0x0),
    fGrid3(0x0),
    fGrid4(0x0),// first: loop over tracks in ESD 
    fNphi(0),
    fNeta(0),
    fPhi2(0x0),
    fEta2(0x0),
    fPhi(0x0),
    fEta(0x0),
    fIndex(0x0),
    fParams(0x0),
    fGrid(0),
    fPhiMin(0),
    fPhiMax(0),
    fEtaMin(0),
    fEtaMax(0),
    fEtaBinInTPCAcc(0),
    fPhiBinInTPCAcc(0),
    fEtaBinInEMCalAcc(0),
    fPhiBinInEMCalAcc(0),
    fNbinPhi(0)
{
  // constructor 2
}

//_____________________________________________________________________________
AliJetFillUnitArrayTracks::AliJetFillUnitArrayTracks(const AliJetFillUnitArrayTracks &det)
  : TTask(det),//"AliJetFillUnitArrayTracks","Fill Unit Array with tpc/its and emcal information"),
    fNumUnits(det.fNumUnits),
    fEtaMinCal(det.fEtaMinCal),
    fEtaMaxCal(det.fEtaMaxCal),
    fPhiMinCal(det.fPhiMinCal),
    fPhiMaxCal(det.fPhiMaxCal),
    fHadCorr(det.fHadCorr),
    fHCorrection(det.fHCorrection),
    fNTracks(det.fNTracks),
    fNTracksCut(det.fNTracksCut),
    fOpt(det.fOpt),
    fDZ(det.fDZ),
    fDebug(det.fDebug),
    fReaderHeader(det.fReaderHeader),
    fMomentumArray(det.fMomentumArray),
    fUnitArray(det.fUnitArray),
    fRefArray(det.fRefArray),
    fProcId(det.fProcId),
    fTPCGrid(det.fTPCGrid),
    fEMCalGrid(det.fEMCalGrid),
    fGeom(det.fGeom),
    fESD(det.fESD),
    fGrid0(det.fGrid0),
    fGrid1(det.fGrid1),
    fGrid2(det.fGrid2),
    fGrid3(det.fGrid3),
    fGrid4(det.fGrid4),
    fNphi(det.fNphi),
    fNeta(det.fNeta),
    fPhi2(det.fPhi2),
    fEta2(det.fEta2),
    fPhi(det.fPhi),
    fEta(det.fEta),
    fIndex(det.fIndex),
    fParams(det.fParams),
    fGrid(det.fGrid),
    fPhiMin(det.fPhiMin),
    fPhiMax(det.fPhiMax),
    fEtaMin(det.fEtaMin),
    fEtaMax(det.fEtaMax),
    fEtaBinInTPCAcc(det.fEtaBinInTPCAcc),
    fPhiBinInTPCAcc(det.fPhiBinInTPCAcc),
    fEtaBinInEMCalAcc(det.fEtaBinInEMCalAcc),
    fPhiBinInEMCalAcc(det.fPhiBinInEMCalAcc),
    fNbinPhi(det.fNbinPhi)
{
  // Copy constructor 
}

//------------------------------------------------------------------
AliJetFillUnitArrayTracks& AliJetFillUnitArrayTracks::operator=(const AliJetFillUnitArrayTracks& other)
{
  // Assignment

    fNumUnits = other.fNumUnits;
    fEtaMinCal = other.fEtaMinCal;
    fEtaMaxCal = other.fEtaMaxCal;
    fPhiMinCal = other.fPhiMinCal;
    fPhiMaxCal = other.fPhiMaxCal;
    fHadCorr = other.fHadCorr;
    fHCorrection = other.fHCorrection;
    fNTracks = other.fNTracks;
    fNTracksCut = other.fNTracksCut;
    fOpt = other.fOpt;
    fDZ = other.fDZ;
    fDebug = other.fDebug;
    fReaderHeader = other.fReaderHeader;
    fMomentumArray = other.fMomentumArray;
    fUnitArray = other.fUnitArray;
    fRefArray = other.fRefArray;
    fProcId = other.fProcId;
    fTPCGrid = other.fTPCGrid;
    fEMCalGrid = other.fEMCalGrid;
    fGeom = other.fGeom;
    fESD = other.fESD;
    fGrid0 = other.fGrid0;
    fGrid1 = other.fGrid1;
    fGrid2 = other.fGrid2;
    fGrid3 = other.fGrid3;
    fGrid4 = other.fGrid4;
    fNphi = other.fNphi;
    fNeta = other.fNeta;
    fPhi2 = other.fPhi2;
    fEta2 = other.fEta2;
    fPhi = other.fPhi;
    fEta = other.fEta;
    fIndex = other.fIndex;
    fParams = other.fParams;
    fGrid = other.fGrid;
    fPhiMin = other.fPhiMin;
    fPhiMax = other.fPhiMax;
    fEtaMin = other.fEtaMin;
    fEtaMax = other.fEtaMax;
    fEtaBinInTPCAcc = other.fEtaBinInTPCAcc;
    fPhiBinInTPCAcc = other.fPhiBinInTPCAcc;
    fEtaBinInEMCalAcc = other.fEtaBinInEMCalAcc;
    fPhiBinInEMCalAcc = other.fPhiBinInEMCalAcc;
    fNbinPhi = other.fNbinPhi;

    return (*this);
}

//____________________________________________________________________________
void AliJetFillUnitArrayTracks::InitParameters()
{
  fHadCorr        = 0;             // For hadron correction
  fNumUnits = fGeom->GetNCells();  // Number of towers in EMCAL
  //fDebug = fReaderHeader->GetDebug();

  fEtaMinCal = fGeom->GetArm1EtaMin();
  fEtaMaxCal = fGeom->GetArm1EtaMax();
  fPhiMinCal = fGeom->GetArm1PhiMin();
  fPhiMaxCal = fGeom->GetArm1PhiMax()-10.;

  fTPCGrid->GetAccParam(fNphi,fNeta,fPhiMin, 
			fPhiMax,fEtaMin,fEtaMax);
  fTPCGrid->GetBinParam(fPhiBinInTPCAcc,fEtaBinInTPCAcc, 
			fPhiBinInEMCalAcc,fEtaBinInEMCalAcc,fNbinPhi);

  fEta   = fTPCGrid->GetArrayEta();
  fPhi   = fTPCGrid->GetArrayPhi();
  fIndex = fTPCGrid->GetIndexObject();

  if(fDebug>20){
    for(Int_t i=0; i<fNphi+1; i++) cout << "phi[" << i << "] : " << (*fPhi)[i] << endl;
    for(Int_t i=0; i<fNeta+1; i++) cout << "eta[" << i << "] : " << (*fEta)[i] << endl;
    
    for(Int_t i=0; i<fNphi+1; i++)
      for(Int_t j=0; j<fNeta+1; j++) {cout << "fIndex[" << i << "," << j << "] : " <<
	  (*fIndex)(i,j) << endl; }
  } 
  if(fDebug>1) printf("\n Parameters initiated ! \n");
}

//_____________________________________________________________________________
AliJetFillUnitArrayTracks::~AliJetFillUnitArrayTracks()
{
  // destructor
}

//_____________________________________________________________________________
void AliJetFillUnitArrayTracks::Exec(Option_t* /*option*/)
{
  //
  // Main method.
  //

  fDebug = fReaderHeader->GetDebug();
  
  // Set parameters
  InitParameters();

  // get number of tracks in event (for the loop)
  Int_t goodTrack = 0;
  Int_t nt = 0;
  Int_t nmax = 0;
  Float_t pt,eta,phi;
  //  Int_t sflag = 0;
  TVector3 p3;

  if(fDebug>1) cout << "fESD in Fill array : " << fESD << endl;

  nt = fESD->GetNumberOfTracks();
  if(fDebug>1) cout << "Number of Tracks in ESD : " << nt << endl;
 
  // temporary storage of signal and pt cut flag
  Int_t* sflag  = new Int_t[nt];
  Int_t* cflag  = new Int_t[nt];

  // get cuts set by user
  Float_t ptMin  = fReaderHeader->GetPtCut();
  Float_t etaMin = fReaderHeader->GetFiducialEtaMin();
  Float_t etaMax = fReaderHeader->GetFiducialEtaMax();  
  Float_t phiMin = fReaderHeader->GetFiducialPhiMin();
  Float_t phiMax = fReaderHeader->GetFiducialPhiMax();  
  fOpt = fReaderHeader->GetDetector();
  fDZ  = fReaderHeader->GetDZ();

  Int_t nTracksEmcal      = 0;
  Int_t nTracksEmcalDZ    = 0;
  Int_t nTracksTpc        = 0;
  Int_t nTracksTpcOnly    = 0;
  Int_t nTracksEmcalCut   = 0;
  Int_t nTracksEmcalDZCut = 0;
  Int_t nTracksTpcCut     = 0;
  Int_t nTracksTpcOnlyCut = 0;

  fGrid = fTPCGrid->GetGridType();

  // Loop over tracks
  nmax = nt;  
  for (Int_t it = 0; it < nmax; it++) {
    AliESDtrack *track;
    track = fESD->GetTrack(it);
    UInt_t status = track->GetStatus();
    
    Double_t mom[3];
    track->GetPxPyPz(mom);
    
    p3.SetXYZ(mom[0],mom[1],mom[2]);
    vector<Float_t> vtmp(3);
    vtmp[0] = mom[0]; vtmp[1] = mom[1]; vtmp[2] = mom[2];
    pt = p3.Pt();
    
    Float_t mass = 0.;
    mass = track->GetMass();

    if ((status & AliESDtrack::kTPCrefit) == 0)    continue;      // quality check
    if ((status & AliESDtrack::kITSrefit) == 0)    continue;      // quality check
    if (((AliJetESDReaderHeader*) fReaderHeader)->ReadSignalOnly() 
	&& TMath::Abs(track->GetLabel()) > 10000)  continue;   // quality check
    if (((AliJetESDReaderHeader*) fReaderHeader)->ReadBkgdOnly() 
	&& TMath::Abs(track->GetLabel()) < 10000)  continue;   // quality check
    eta = p3.Eta();
    phi = ( (p3.Phi()) < 0) ? (p3.Phi()) + 2. * TMath::Pi() : (p3.Phi());
       
    if ( (eta > etaMax) || (eta < etaMin)) continue;           // checking eta cut
    if ( (phi > phiMax) || (phi < phiMin)) continue;           // checking phi cut
      
    sflag[goodTrack]=0;
    if (TMath::Abs(track->GetLabel()) < 10000) sflag[goodTrack]=1;
    cflag[goodTrack]=0;
    if (pt > ptMin) cflag[goodTrack]=1;                       // pt cut

    if(fDebug>20) cout << "In Pythia: Track:" << it << ", eta: " << eta << ", phi: " << phi << ", pt: " << pt << endl;

    if(fGrid==0)
      {
	// Only TPC filled from its grid in its total acceptance
	
	Int_t idTPC = fTPCGrid->GetIndex(phi,eta);
	Bool_t ok = kFALSE;
	Bool_t ref = kFALSE;
	
	AliJetUnitArray *uArray = (AliJetUnitArray*)fUnitArray->At(idTPC-1);
	
	uArray->SetUnitTrackID(it);
	
	Float_t unitEnergy = 0.;
	unitEnergy = uArray->GetUnitEnergy();
	// nTracksTpcOnly is necessary to count the number of candidate cells
	// but it doesn't give the real multiplicity -> it will be extracted 
	// in the jet finder and for that, it is necessary to fill the Px,Py,Pz
	// information for each tracks stored in a given unitcell 
	if(unitEnergy==0.){
	  nTracksTpcOnly++;
	  ok = kTRUE;
	  ref = kTRUE;
	}

	// Fill energy in TPC acceptance
	uArray->SetUnitEnergy(unitEnergy + pt);
	uArray->SetUnitPxPyPz(kFALSE,vtmp);
	uArray->SetUnitMass(mass);

	// Pt cut flag
	if(uArray->GetUnitEnergy()<ptMin){
	  uArray->SetUnitCutFlag(kPtSmaller);
	}
	else {
	  uArray->SetUnitCutFlag(kPtHigher);
	  if(ok) nTracksTpcOnlyCut++;
	}

	// Signal flag
	if(sflag[goodTrack] == 1) { 
	  uArray->SetUnitSignalFlag(kGood);
	  uArray->SetUnitSignalFlagC(kFALSE,kGood);
	} else uArray->SetUnitSignalFlagC(kFALSE,kBad);

	if(uArray->GetUnitEnergy()>0 && ref){
	  if(!fProcId) {
	    fProcId = kTRUE;
	    // delete fRefArray;
	    new(fRefArray) TRefArray(TProcessID::GetProcessWithUID(uArray));
	  }
	  fRefArray->Add(uArray);
	}
	
      }
    
    if(fGrid==1)
      {
	Int_t nElements = fTPCGrid->GetNEntries();
	// Fill track information in EMCAL acceptance
	if((eta >= fEtaMin && eta <= fEtaMax) &&
	   (phi >= fPhiMin && phi <= fPhiMax))// &&
	  {
	    // Include dead-zones
	    if(fDZ)
	      {
		Double_t phimin0 = 0., phimin1 = 0., phimin2 = 0., phimin3 = 0., phimin4 = 0.;
		Double_t phimax0 = 0., phimax1 = 0., phimax2 = 0., phimax3 = 0., phimax4 = 0.;
		fGeom->GetPhiBoundariesOfSMGap(0,phimin0,phimax0);
		fGeom->GetPhiBoundariesOfSMGap(1,phimin1,phimax1);
		fGeom->GetPhiBoundariesOfSMGap(2,phimin2,phimax2);
		fGeom->GetPhiBoundariesOfSMGap(3,phimin3,phimax3);
		fGeom->GetPhiBoundariesOfSMGap(4,phimin4,phimax4);
		Int_t n0 = fGrid0->GetNEntries();
		Int_t n1 = fGrid1->GetNEntries();
		Int_t n2 = fGrid2->GetNEntries();
		Int_t n3 = fGrid3->GetNEntries();
		
		if(phi >= phimin0 && phi <= phimax0){
		  Int_t id0 = fGrid0->GetIndex(phi,eta)-1;
		  AliJetUnitArray *uArray0 = (AliJetUnitArray*)fUnitArray->At(id0+fNumUnits+nElements);
		  
		  uArray0->SetUnitTrackID(it);
		  
		  Float_t uEnergy0 = uArray0->GetUnitEnergy();
		  Bool_t ok0 = kFALSE;
		  Bool_t ref0 = kFALSE;
		  if(uEnergy0==0.){
		    nTracksEmcalDZ++;
		    ok0 = kTRUE;
		    ref0 = kTRUE;
		  }
		  uArray0->SetUnitEnergy(uEnergy0+pt);
		  uArray0->SetUnitPxPyPz(kFALSE,vtmp);
		  if(uArray0->GetUnitEnergy()<ptMin)
		    uArray0->SetUnitCutFlag(kPtSmaller);
		  else {
		    uArray0->SetUnitCutFlag(kPtHigher);
		    if(ok0) nTracksEmcalDZCut++;
		  }
		  if(sflag[goodTrack] == 1) {
		    uArray0->SetUnitSignalFlag(kGood);
		    uArray0->SetUnitSignalFlagC(kFALSE,kGood);
		  } else uArray0->SetUnitSignalFlagC(kFALSE,kBad);
		  
		  if(uArray0->GetUnitEnergy()>0 && ref0){
		    if(!fProcId){
		      new(fRefArray) TRefArray(TProcessID::GetProcessWithUID(uArray0));
		      fProcId = kTRUE;
		    }
		    fRefArray->Add(uArray0);
		  }
		}
		if(phi >= phimin1 && phi <= phimax1){
		  Int_t id1 = fGrid1->GetIndex(phi,eta)-1+n0;
		  AliJetUnitArray *uArray1 = (AliJetUnitArray*)fUnitArray->At(id1+fNumUnits+nElements);
		  uArray1->SetUnitTrackID(it);

		  Float_t uEnergy1 = uArray1->GetUnitEnergy();
		  Bool_t ok1 = kFALSE;
		  Bool_t ref1 = kFALSE;
		  if(uEnergy1==0.){
		    nTracksEmcalDZ++;
		    ok1 = kTRUE;
		    ref1 = kTRUE;
		  }
		  uArray1->SetUnitEnergy(uEnergy1+pt);
		  uArray1->SetUnitPxPyPz(kFALSE,vtmp);
		  if(uArray1->GetUnitEnergy()<ptMin)
		    uArray1->SetUnitCutFlag(kPtSmaller);
		  else {
		    uArray1->SetUnitCutFlag(kPtHigher);
		    if(ok1) nTracksEmcalDZCut++;
		  }
		  if(sflag[goodTrack] == 1) {
		    uArray1->SetUnitSignalFlag(kGood);
		    uArray1->SetUnitSignalFlagC(kFALSE,kGood);
		  } else uArray1->SetUnitSignalFlagC(kFALSE,kBad);
		  
		  if(uArray1->GetUnitEnergy()>0 && ref1){
		    if(!fProcId){
		      new(fRefArray) TRefArray(TProcessID::GetProcessWithUID(uArray1));
		      fProcId = kTRUE;
		    }
		    fRefArray->Add(uArray1);
		  }
		}
		if(phi >= phimin2 && phi <= phimax2){
		  Int_t id2 = fGrid2->GetIndex(phi,eta)-1+n0+n1;
		  AliJetUnitArray *uArray2 = (AliJetUnitArray*)fUnitArray->At(id2+fNumUnits+nElements);
		  uArray2->SetUnitTrackID(it);

		  Float_t uEnergy2 = uArray2->GetUnitEnergy();
		  Bool_t ok2 = kFALSE;
		  Bool_t ref2 = kFALSE;
		  if(uEnergy2==0.){
		    nTracksEmcalDZ++;
		    ok2 = kTRUE;
		    ref2 = kTRUE;
		  }
		  uArray2->SetUnitEnergy(uEnergy2+pt);
		  uArray2->SetUnitPxPyPz(kFALSE,vtmp);
		  if(uArray2->GetUnitEnergy()<ptMin)
		    uArray2->SetUnitCutFlag(kPtSmaller);
		  else {
		    uArray2->SetUnitCutFlag(kPtHigher);
		    if(ok2) nTracksEmcalDZCut++;
		  }
		  if(sflag[goodTrack] == 1) {
		    uArray2->SetUnitSignalFlag(kGood);
		    uArray2->SetUnitSignalFlagC(kFALSE,kGood);
		  } else uArray2->SetUnitSignalFlagC(kFALSE,kBad);
		  
		  if(uArray2->GetUnitEnergy()>0 && ref2){
		    if(!fProcId){
		      new(fRefArray) TRefArray(TProcessID::GetProcessWithUID(uArray2));
		      fProcId = kTRUE;
		    }
		    fRefArray->Add(uArray2);
		  }
		}
		if(phi >= phimin3 && phi <= phimax3){
		  Int_t id3 = fGrid3->GetIndex(phi,eta)-1+n0+n1+n2;
		  AliJetUnitArray *uArray3 = (AliJetUnitArray*)fUnitArray->At(id3+fNumUnits+nElements);
		  uArray3->SetUnitTrackID(it);

		  Float_t uEnergy3 = uArray3->GetUnitEnergy();
		  Bool_t ok3 = kFALSE;
		  Bool_t ref3 = kFALSE;
		  if(uEnergy3==0.){
		    nTracksEmcalDZ++;
		    ok3 = kTRUE;
		    ref3 = kTRUE;
		  }
		  uArray3->SetUnitEnergy(uEnergy3+pt);
		  uArray3->SetUnitPxPyPz(kFALSE,vtmp);
		  if(uArray3->GetUnitEnergy()<ptMin)
		    uArray3->SetUnitCutFlag(kPtSmaller);
		  else {
		    uArray3->SetUnitCutFlag(kPtHigher);
		    if(ok3) nTracksEmcalDZCut++;
		  }
		  if(sflag[goodTrack] == 1) {
		    uArray3->SetUnitSignalFlag(kGood);
		    uArray3->SetUnitSignalFlagC(kFALSE,kGood);
		  } else uArray3->SetUnitSignalFlagC(kFALSE,kBad);
		  
		  if(uArray3->GetUnitEnergy()>0 && ref3){
		    if(!fProcId){
		      new(fRefArray) TRefArray(TProcessID::GetProcessWithUID(uArray3));
		      fProcId = kTRUE;
		    }
		    fRefArray->Add(uArray3);
		  }
		}
		if(phi >= phimin4 && phi <= phimax4){
		  Int_t id4 = fGrid4->GetIndex(phi,eta)-1+n0+n1+n2+n3;
		  AliJetUnitArray *uArray4 = (AliJetUnitArray*)fUnitArray->At(id4+fNumUnits+nElements);
		  uArray4->SetUnitTrackID(it);
		  
		  Float_t uEnergy4 = uArray4->GetUnitEnergy();
		  Bool_t ok4 = kFALSE;
		  Bool_t ref4 = kFALSE;
		  if(uEnergy4==0.){
		    nTracksEmcalDZ++;
		    ok4 = kTRUE;
		    ref4 = kTRUE;
		  }
		  uArray4->SetUnitEnergy(uEnergy4+pt);
		  uArray4->SetUnitPxPyPz(kFALSE,vtmp);
		  if(uArray4->GetUnitEnergy()<ptMin)
		    uArray4->SetUnitCutFlag(kPtSmaller);
		  else {
		    uArray4->SetUnitCutFlag(kPtHigher);
		    if(ok4) nTracksEmcalDZCut++;
		  }
		  if(sflag[goodTrack] == 1) {
		    uArray4->SetUnitSignalFlag(kGood);
		    uArray4->SetUnitSignalFlagC(kFALSE,kGood);
		  } else uArray4->SetUnitSignalFlagC(kFALSE,kBad);
		  
		  if(uArray4->GetUnitEnergy()>0 && ref4){
		    if(!fProcId){
		      new(fRefArray) TRefArray(TProcessID::GetProcessWithUID(uArray4));
		      fProcId = kTRUE;
		    }
		    fRefArray->Add(uArray4);
		  }
		}
	      } // end fDZ
	    
	    Int_t towerID = 0;
	    fGeom->GetAbsCellIdFromEtaPhi(eta,phi,towerID);
	    if(towerID==-1) continue;
	      
	    AliJetUnitArray *uArray = (AliJetUnitArray*)fUnitArray->At(towerID);
	    uArray->SetUnitTrackID(it);

	    Float_t unitEnergy = uArray->GetUnitEnergy(); 
	    Bool_t ok = kFALSE;
	    Bool_t ref = kFALSE;
	    if(unitEnergy==0.){
	      nTracksEmcal++;
	      ok=kTRUE;
	      ref=kTRUE;
	    }

	    // Do Hadron Correction
	    // Parametrization to be added
	    if (fHCorrection != 0) 
	      { 
		Float_t   hCEnergy = fHadCorr->GetEnergy(p3.Mag(), (Double_t)eta);
		unitEnergy -= hCEnergy*TMath::Sin(2.0*TMath::ATan(TMath::Exp(-eta)));
		
	      } //end Hadron Correction loop
	    
	    uArray->SetUnitEnergy(unitEnergy + pt);
	    uArray->SetUnitPxPyPz(kFALSE,vtmp);
	    // Put a pt cut flag
	    if(uArray->GetUnitEnergy()<ptMin){
	      uArray->SetUnitCutFlag(kPtSmaller);
	    }
	    else {
	      uArray->SetUnitCutFlag(kPtHigher);
	      if(ok) nTracksEmcalCut++;
	    }
	    
	    // Signal flag
	    if(sflag[goodTrack] == 1) { 
	      uArray->SetUnitSignalFlag(kGood);
	      uArray->SetUnitSignalFlagC(kFALSE,kGood);
	    } else uArray->SetUnitSignalFlagC(kFALSE,kBad);
	    
	    if(uArray->GetUnitEnergy()>0 && ref){
	      if(!fProcId){
		new(fRefArray) TRefArray(TProcessID::GetProcessWithUID(uArray));
		fProcId = kTRUE;
	      }
	      fRefArray->Add(uArray);
	    }
	  } // end loop on EMCal acceptance cut + tracks quality
	else{ 
	  // Outside EMCal acceptance
	  
	  Int_t idTPC = fTPCGrid->GetIndex(phi,eta);
	  
	  AliJetUnitArray *uArray = (AliJetUnitArray*)fUnitArray->At(fNumUnits-1+idTPC);
	  uArray->SetUnitTrackID(it);

	  Float_t unitEnergy2 = uArray->GetUnitEnergy(); // check if fNumUnits or fNumUnits-1
	  Bool_t okout = kFALSE;
	  Bool_t refout = kFALSE;
	  if(unitEnergy2==0.){
	    nTracksTpc++;
	    okout=kTRUE;
	    refout=kTRUE;
	  }
	  // Fill energy outside emcal acceptance
	  uArray->SetUnitEnergy(unitEnergy2 + pt);
	  uArray->SetUnitPxPyPz(kFALSE,vtmp);
	  
	  // Pt cut flag
	  if(uArray->GetUnitEnergy()<ptMin){
	    uArray->SetUnitCutFlag(kPtSmaller);
	  }
	  else {
	    uArray->SetUnitCutFlag(kPtHigher);
	    if(okout) nTracksTpcCut++;
	  }
	  // Signal flag
	  if(sflag[goodTrack] == 1) {
	    uArray->SetUnitSignalFlag(kGood);
	    uArray->SetUnitSignalFlagC(kFALSE,kGood);
	  } else uArray->SetUnitSignalFlagC(kFALSE,kBad);
	  
	  if(uArray->GetUnitEnergy()>0 && refout){
	    if(!fProcId){
	      new(fRefArray) TRefArray(TProcessID::GetProcessWithUID(uArray));
	      fProcId = kTRUE;
	    }
	    fRefArray->Add(uArray);
	  }
	}
      } // end fGrid==1

    goodTrack++;

  } // end loop on entries (tpc tracks)

  if(fDebug>0) 
    {
      cout << "End of Tracks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
      cout << "goodTracks: " << goodTrack << endl;
    } 

  delete[] sflag;
  delete[] cflag;

  if(fGrid==0) {
    fNTracks = nTracksTpcOnly;
    fNTracksCut = nTracksTpcOnlyCut;
    if(fDebug>10){
      cout << "fNTracks : " << fNTracks << endl;
      cout << "fNTracksCut : " << fNTracksCut << endl;
    }
  }
  if(fGrid==1) {
    fNTracks = nTracksEmcal+nTracksEmcalDZ+nTracksTpc;
    fNTracksCut = nTracksEmcalCut+nTracksEmcalDZCut+nTracksTpcCut;
    if(fDebug>10){
      cout << "fNTracks : " << fNTracks << endl;
      cout << "fNTracksCut : " << fNTracksCut << endl;
    }
  }  
  
}

//__________________________________________________________
void AliJetFillUnitArrayTracks::GetEtaPhiFromIndex(Int_t index, Float_t &eta, Float_t &phi)
{
  for(Int_t j=0; j<fNphi+1; j++) 
    {
      for(Int_t i=0; i<fNeta+1; i++) 
	{

	  // TPC grid only 
	  //-------------------------------------
	  if(fGrid==0) {	
	    if(j*(fNeta+1)+i == index) {
	      eta = fEta2->At(i); 
	      phi = fPhi2->At(j);
	    }
	  }
	  
	  // TPC-EMCAL grid
	  //-------------------------------------
	  Int_t ii = 0;
	  if(i==0) ii = 0;
	  if(i>0 && i<(fEtaBinInTPCAcc-fEtaBinInEMCalAcc)/2) ii = i; 
	  if(i>=(fEtaBinInTPCAcc+fEtaBinInEMCalAcc)/2 && i<fNeta+1) ii = i-fEtaBinInEMCalAcc;
	  
	  if(fGrid==1) {
	    if(j<(fNbinPhi+1) && j*(fNeta+1)+i == index) {
	      eta = fEta2->At(i);
	      phi = fPhi2->At(j);
	    }  
	    
	    if((j>=(fNbinPhi+1) && j<(fNbinPhi+1+fPhiBinInEMCalAcc)) && 
	       ((fNbinPhi+1)*(fNeta+1) + (j-fNbinPhi-1)*(fEtaBinInTPCAcc-fEtaBinInEMCalAcc) + ii)== index ) {
	      if(ii==0) {Int_t ind = 0; eta = fEta2->At(ind);}
	      else eta = fEta2->At(i);
	      phi = fPhi2->At(j);
	    }
	    
	    if(j>=(fNbinPhi+1+fPhiBinInEMCalAcc) && ((fNbinPhi+1)*(fNeta+1)+fPhiBinInEMCalAcc*((fEtaBinInTPCAcc-fEtaBinInEMCalAcc))+(j-(fNbinPhi+1+fPhiBinInEMCalAcc))*(fNeta+1)+i == index)) {
	      eta = fEta2->At(i);
	      phi = fPhi2->At(j);
	    }
	  }
	}
    }
}












