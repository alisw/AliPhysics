
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
// *** November 2006
// Author: magali.estienne@ires.in2p3.fr
// 1) Define 2 grids and take (eta,phi) from the grid or use the grid for the TPC and  
// EtaPhiFromIndex and TowerIndexFromEtaPhi for the particles in EMCAL acceptance
//     2 options are possible : fGrid==0, work only with full TPC acceptance (for the moment)
//                              fGrid==1, work with a part of the TPC acceptance  
// 2) Try to implement 2 full grids for TPC and EMCal separately and to merge them
// 3) Need to include Dead-zone -> Wait for exact positions in the new detector geometry
// Author: Magali Estienne (magali.estienne@ires.in2p3.fr)
//======================================================================
// ***September 2006
// TTask : Fill Unit Array for the Tracks information 
// Called by ESD reader for jet analysis
// Author: Magali Estienne (magali.estienne@ires.in2p3.fr)
//======================================================================
// *** July 2006
// 1) When the tracks are in the EMCal acceptance, the functions EtaPhiFromIndex
// and TowerIndexFromEtaPhi in the AliEMCALGeometry class are used to extract the
// index or the eta, phi position of a grid.
// 2) Define a grid for TPC
// Author: Magali Estienne (magali.estienne@ires.in2p3.fr)
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
#include <TVector3.h>
//#include "Math/Vector3D.h"
//#include "Math/Vector3Dfwd.h"
#include "TTask.h"
#include <TGeoManager.h>
#include <TMatrixD.h>
#include <TArrayD.h>

// --- AliRoot header files ---
#include "AliJetFinder.h"
#include "AliJetReaderHeader.h"
#include "AliJetReader.h"
#include "AliJetESDReader.h"
#include "AliJetESDReaderHeader.h"
#include "AliESD.h"
//#include "AliEMCALGeometry.h"
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
      fEtaMinCal(0.),
      fEtaMaxCal(0.),
      fPhiMinCal(0.),
      fPhiMaxCal(0.),
      fHadCorr(0x0),
      fHCorrection(0),
      fNIn(0),
      fOpt(0),
      fDebug(0),
      fReaderHeader(0x0),
      fMomentumArray(0x0),
      fUnitArray(0x0),
      fTPCGrid(0x0),
      fEMCalGrid(0x0),
      fGeom(0x0),
      fNphi(0),
      fNeta(0),
      fPhi2(0x0),
      fEta2(0x0),
      fPhi(0x0),
      fEta(0x0),
      fIndex(0x0),
      fParams(0x0),
      fGrid(0),
      fPhiMin(0.),
      fPhiMax(0.),
      fEtaMin(0.),
      fEtaMax(0.),
      fEtaBinInTPCAcc(0),
      fPhiBinInTPCAcc(0),
      fEtaBinInEMCalAcc(0),
      fPhiBinInEMCalAcc(0),
      fNbinPhi(0)
{
  // constructor
}

//____________________________________________________________________________
void AliJetFillUnitArrayTracks::SetEMCALGeometry()
{
    // Set EMCAL geometry information
    fGeom = AliJetDummyGeo::GetInstance();
    if (fGeom == 0)
	fGeom = AliJetDummyGeo::GetInstance("SHISH_77_TRD1_2X2_FINAL_110DEG","EMCAL");
    if(fDebug>1) printf("\n EMCAL Geometry setted ! \n");
}

//____________________________________________________________________________
void AliJetFillUnitArrayTracks::InitParameters()
{
  fHCorrection    = 0;     // For hadron correction
  fHadCorr        = 0;     // For hadron correction
  fNumUnits = fGeom->GetNCells();      // Number of towers in EMCAL
  fDebug = fReaderHeader->GetDebug();

  fEtaMinCal = fGeom->GetArm1EtaMin();
  fEtaMaxCal = fGeom->GetArm1EtaMax();
  fPhiMinCal = fGeom->GetArm1PhiMin();
  fPhiMaxCal = fGeom->GetArm1PhiMax()-10.; // A verifier quelle doit etre la derniere valeur !

  if(fDebug>30){
    cout << "fEtaMinCal : " << fEtaMinCal << endl;
    cout << "fEtaMaxCal : " << fEtaMaxCal << endl;
    cout << "fPhiMinCal : " << fPhiMinCal << endl;
    cout << "fPhiMaxCal : " << fPhiMaxCal << endl;
  }

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
void AliJetFillUnitArrayTracks::Exec(Option_t* option)
{
  //
  // Main method.
  // Explain

  fDebug = fReaderHeader->GetDebug();
  if(fDebug>1) printf("In AliJetFillUnitArrayTracks !");  
  if(fDebug>3) printf("\nfGeom->GetEntries() = %d\n", fGeom->GetNCells());
  // Set EMCal Geometry
  SetEMCALGeometry();
  // Set parameters
  InitParameters();

  TClonesArray *lvArray = fMomentumArray;     // Correct checked !
  Int_t nInT =  lvArray->GetEntries();        // Correct checked !
  Float_t ptMin = fReaderHeader->GetPtCut();  // Correct checked !
  
  // sflag -> not yet implemented !!!!
  
  if(fDebug>3) cout << "nInT : " << nInT << endl;

  if (nInT == 0) return;
  
  // local arrays for input
  Float_t* enT  = new Float_t[nInT];
  Float_t* ptT  = new Float_t[nInT];
  Float_t* etaT = new Float_t[nInT];
  Float_t* phiT = new Float_t[nInT];
  Float_t* thetaT = new Float_t[nInT];
  
  Int_t trackInEmcalAcc = 0;
  Int_t trackInTpcAcc = 0;
  Int_t trackInTpcAccOnly = 0;

  Int_t nElements = fTPCGrid->GetNEntries();
  Int_t nElements2 = fEMCalGrid->GetNEntries();
  fGrid = fTPCGrid->GetGridType();

  if(fDebug>3){
    cout << "nElements : " << nElements << endl;
    cout << "nElements2 : " << nElements2 << endl;
    cout << "fNumUnits : " << fNumUnits << endl;
    cout << "sum : " << fNumUnits+nElements << endl;
  }

  // Set energy exactly to zero
  if(fGrid==0)
    for(Int_t k=0; k<nElements; k++)
      fUnitArray[k].SetUnitEnergy(0.); 
  
  if(fGrid==1)
    for(Int_t k=0; k<fNumUnits+nElements; k++)
      fUnitArray[k].SetUnitEnergy(0.);  

  // load input vectors
  for (Int_t i = 0; i < nInT; i++)  
    {
      TLorentzVector *lv = (TLorentzVector*) lvArray->At(i);
      enT[i]  = lv->Energy();
      ptT[i]  = lv->Pt();
      etaT[i] = lv->Eta();
      phiT[i] = ((lv->Phi() < 0) ? (lv->Phi()) + 2. * TMath::Pi() : lv->Phi());
      thetaT[i] = 2.0*TMath::ATan(TMath::Exp(-etaT[i]));
      
      if(fDebug>20){    
	cout << "enT[" << i << "] : " <<  enT[i] << endl;
	cout << "ptT[" << i << "] : " <<  ptT[i] << endl;
	cout << "etaT[" << i << "] : " <<  etaT[i] << endl;
	cout << "phiT[" << i << "] : " <<  phiT[i] << endl;
	cout << "thetaT[" << i << "] : " <<  thetaT[i] << endl;
	cout << "fEtaMinCal : " << fEtaMinCal << ", fEtaMaxCal : " << fEtaMaxCal << endl;
	cout << "fPhiMinCal : " << fPhiMinCal << ", fPhiMaxCal : " << fPhiMaxCal << endl;
	cout << "fEtaMin : " << fEtaMin << ", fEtaMax : " << fEtaMax << endl;
	cout << "fPhiMin : " << fPhiMin << ", fPhiMax : " << fPhiMax << endl;
      }
      
      if(fGrid==0)
	{
	  // For the moment, only TPC filled from its grid in its total acceptance
	  if(fDebug>2)
	    cout << "In total TPC acceptance +++++++++++++++++++++++++++++++++++++++++++" << endl;
	    
	    trackInTpcAccOnly += 1;
	    
	    Int_t idTPC = fTPCGrid->GetIndex(phiT[i],etaT[i]);
	    
	    Float_t unitEnergy = 0.;
	    unitEnergy = fUnitArray[idTPC].GetUnitEnergy();
	    
	    if(unitEnergy > 0. && fDebug >10){
	      cout << "§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§" << endl;
	      cout << "§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§" << endl;
	      cout << "TPC id : " << idTPC << endl;
	      cout << "etaT[" << i << "] : " << etaT[i] << endl;
	      cout << "phiT[" << i << "] : " << phiT[i] << endl;
	      cout << "unitEnergy in TPC acceptance : " << unitEnergy << endl;
	      cout << "fUnitArray[idTPC].GetUnitEnergy(): " << 
		fUnitArray[idTPC].GetUnitEnergy() << endl;
	      cout << "§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§" << endl;
	      cout << "§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§" << endl;
	    }

	    // Fill energy in TPC acceptance
	    fUnitArray[idTPC].SetUnitEnergy(unitEnergy + ptT[i]);

	    if(fDebug > 10){
	      cout << "ptT[" << i << "] : " << ptT[i] << endl;
	      cout << "unitEnergy in TPC acceptance after : " << 
		fUnitArray[idTPC].GetUnitEnergy() << endl;
	    }
	    
	    // Pt cut flag
	    if(fUnitArray[idTPC].GetUnitEnergy()<ptMin)
	      fUnitArray[idTPC].SetUnitCutFlag(kPtSmaller);
	    else fUnitArray[idTPC].SetUnitCutFlag(kPtHigher);
	    // Detector flag
	    fUnitArray[idTPC].SetUnitDetectorFlag(kTpc);
	}

      if(fGrid==1)
	{
	  // Fill track information in EMCAL acceptance
	  if((etaT[i] >= fEtaMin && etaT[i] <= fEtaMax) &&
	     (phiT[i] >= fPhiMin && phiT[i] <= fPhiMax))// &&
	    //	 GetCutFlag(i) == 1)
	    //	 ptT[i] > ptMin)
	    {
	      trackInEmcalAcc += 1;
	      
	      if(fDebug>20){ 
		cout << "before : " << endl;
		cout << "etaT[i] : " << etaT[i] << endl;
		cout << "phiT[i] : " << phiT[i]*180/TMath::Pi() << endl;
	      }
	      
	      // This function should be modified soon
//	      Int_t towerID = fGeom->TowerIndexFromEtaPhi(etaT[i],180.0/TMath::Pi()*phiT[i]); // Obsolete
	      Int_t towerID = fGeom->TowerIndexFromEtaPhi2(etaT[i],180.0/TMath::Pi()*phiT[i]); // Mine modified
//	      Int_t towerID = fEMCalGrid->GetIndexFromEtaPhi(phiT[i],etaT[i]);  // Using an EMCal grid -> calculated
//  	      Int_t towerID = fEMCalGrid->GetIndex(phiT[i],etaT[i]);  // Using an EMCal grid -> tabulated (faster)

	      Float_t unitEnergy = fUnitArray[towerID].GetUnitEnergy(); 
	      
	      if(fDebug>20) { 
		cout << "§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§" << endl;
		cout << "§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§" << endl;
		cout << "after : " << endl;
		cout << "towerID : " << towerID << endl;
		cout << "etaT[i] : " << etaT[i] << endl;
		cout << "phiT[i](rad) : " << phiT[i] << endl;
		cout << "phiT[i] : " << phiT[i]*180/TMath::Pi() << endl;
		cout << "unitEnergy in emcal acceptance : " << unitEnergy << endl;
		cout << "fHadCorr : " << fHadCorr << endl;
		cout << "fHCorrection : " << fHCorrection << endl;
		cout << "§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§" << endl;
		cout << "§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§" << endl;
	      }

	      //OLD WAY:   //Do Hadron Correction
	      if (fHCorrection != 0) 
		{ 
		  Float_t   hCEnergy = fHadCorr->GetEnergy(enT[i], (Double_t)etaT[i]);
		  unitEnergy -= hCEnergy*TMath::Sin(thetaT[i]);
		  if(fDebug>20){
		    cout << "Inside loop for hadron correction ==========================" << endl;
		    cout << "enT[" << i << "] : " << enT[i] << endl;
		    cout << "ptT[" << i << "] : " << ptT[i] << endl;
		    cout << "etaT[" << i << "] : " << etaT[i] << endl;
		    cout << "Energy correction : " << hCEnergy << endl;
		    cout << "unitEnergy : " << unitEnergy << endl;
		    cout << "fUnitArray[towerID].GetUnitEnergy() : " << 
		      fUnitArray[towerID].GetUnitEnergy() << endl;
		  }	      
		} //end Hadron Correction loop
	      
	      fUnitArray[towerID].SetUnitEnergy(unitEnergy + ptT[i]);

	      // Put a pt cut flag
	      if(fUnitArray[towerID].GetUnitEnergy()<ptMin)
		fUnitArray[towerID].SetUnitCutFlag(kPtSmaller);
	      else fUnitArray[towerID].SetUnitCutFlag(kPtHigher);
	      // Detector flag
	      fUnitArray[towerID].SetUnitDetectorFlag(kTpc);
	      
	      if(fDebug>10){
		cout << "After pT filled ===============================================" << endl;
		cout << "PtCut : " << ptMin << endl;
		cout << "ptT[" << i << "] : " << ptT[i] << endl;
		cout << "unitEnergy : " << unitEnergy << endl;
		cout << "fUnitArray[towerID].GetUnitEnergy() : " << fUnitArray[towerID].GetUnitEnergy() << endl;
		cout << "fUnitArray[towerID].GetUnitCutFlag() : " << fUnitArray[towerID].GetUnitCutFlag() << endl;
		cout << "fUnitArray[towerID].GetUnitEta() : " << fUnitArray[towerID].GetUnitEta() << endl;
		cout << "fUnitArray[towerID].GetUnitPhi() : " << fUnitArray[towerID].GetUnitPhi() << endl;
	      }
	    } // end loop on EMCal acceptance cut + tracks quality
	  else{ 
	    // Outside EMCal acceptance
	    if(fDebug>2)
	      cout << "Outside EMCal acceptance +++++++++++++++++++++++++++++++++++++++++++" << endl;
	    
	    trackInTpcAcc += 1;
	    
	    // 	Int_t idTPC = GetIndexFromEtaPhi(etaT[i],phiT[i]);
	    Int_t idTPC = fTPCGrid->GetIndex(phiT[i],etaT[i]);
	    
	    Float_t unitEnergy2 = 0.;
	    unitEnergy2 = fUnitArray[fNumUnits-1+idTPC].GetUnitEnergy(); // check if fNumUnits or fNumUnits-1
	    
	    if(unitEnergy2 > 0. && fDebug >10){
	      cout << "§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§" << endl;
	      cout << "§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§" << endl;
	      cout << "TPC id : " << idTPC << endl;
	      cout << "Total id : " << fNumUnits-1+idTPC << endl;
	      cout << "etaT[" << i << "] : " << etaT[i] << endl;
	      cout << "phiT[" << i << "] : " << phiT[i] << endl;
	      cout << "unitEnergy outside emcal acceptance : " << unitEnergy2 << endl;
	      cout << "fUnitArray[fNumUnits-1+idTPC].GetUnitEnergy(): " << 
		fUnitArray[fNumUnits-1+idTPC].GetUnitEnergy() << endl;
	      cout << "§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§" << endl;
	      cout << "§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§" << endl;
	    }

	    // Fill energy outside emcal acceptance
	    fUnitArray[fNumUnits-1+idTPC].SetUnitEnergy(unitEnergy2 + ptT[i]);

	    // Pt cut flag
	    if(fUnitArray[fNumUnits-1+idTPC].GetUnitEnergy()<ptMin)
	      fUnitArray[fNumUnits-1+idTPC].SetUnitCutFlag(kPtSmaller);
	    else fUnitArray[fNumUnits-1+idTPC].SetUnitCutFlag(kPtHigher);
	    // Detector flag
	    fUnitArray[fNumUnits-1+idTPC].SetUnitDetectorFlag(kTpc);
	  }
	} // end fGrid==1
    } // end loop on entries (tpc tracks)
      
  if(fDebug>3){
    printf("Number of tracks in EMCAL acceptance: %d\n", trackInEmcalAcc);
    printf("Number of tracks in TPC acceptance: %d\n", trackInTpcAcc);
  }
  
  if(fGrid==0) fNIn = trackInTpcAccOnly;
  if(fGrid==1) fNIn = trackInEmcalAcc+trackInTpcAcc;
  
  fOpt = fReaderHeader->GetDetector();
  
  if(fDebug>1) printf("fOpt in FillUnitArrayFromTPCTracks : %d\n", fOpt);
  
      if(fOpt==1 || option=="tpc") // if only TPC
	{ //Set all unit flags, Eta, Phi

	  if(fGrid==0)
	    {
	      if(fDebug>1) printf("In if(fOpt==1 || option==tpc)\n");
	      for(Int_t i=0; i<nElements; i++)
		{
		  Float_t eta = 0.;
		  Float_t phi = 0.;
		  if (fDebug>10) Info("FillUnitArray","Setting all units outside jets \n");
		  fUnitArray[i].SetUnitFlag(kOutJet); // returns 0, 1, 2...
		  fUnitArray[i].SetUnitEntries(nElements);
		  fUnitArray[i].SetUnitID(i);
		  fTPCGrid->GetEtaPhiFromIndex2(fUnitArray[i].GetUnitID(),phi,eta);
		  fUnitArray[i].SetUnitEta(eta);
		  fUnitArray[i].SetUnitPhi(phi);
		  fUnitArray[i].SetUnitDeta(fTPCGrid->GetDeta());
		  fUnitArray[i].SetUnitDphi(fTPCGrid->GetDphi());
		}
	    }

	  if(fGrid==1)
	    {
	      if(fDebug>1) printf("In if(fOpt==1 || option==tpc)\n");
	      for(Int_t i=0; i<fNumUnits+nElements; i++)
		{
		  Float_t eta = 0.;
		  Float_t phi = 0.;
		  if (fDebug>10) Info("FillUnitArray","Setting all units outside jets \n");
		  //Set all units to be outside a jet initially
		  fUnitArray[i].SetUnitFlag(kOutJet); // returns 0, 1, 2...
		  fUnitArray[i].SetUnitEntries(fNumUnits+nElements);
		  fUnitArray[i].SetUnitID(i);

		  if(fUnitArray[i].GetUnitID()<fNumUnits)
		    {
		      // fGeom->EtaPhiFromIndex2(fUnitArray[i].GetUnitID(), eta, phi); // My function in HEADPythia63
		      fGeom->EtaPhiFromIndex(fUnitArray[i].GetUnitID(), eta, phi); // From EMCal geometry 
		      phi = ((phi < 0) ? phi + 2. * TMath::Pi() : phi);
		      // fEMCalGrid->GetEtaPhiFromIndex2(i,phi,eta); // My function from Grid
		      // fEMCalGrid->GetEtaPhiFromIndex2(fUnitArray[i].GetUnitID(),phi,eta); // My function from Grid
		      fUnitArray[i].SetUnitEta(eta);
		      //fUnitArray[i].SetUnitPhi(phi*TMath::Pi()/180.0);
		      fUnitArray[i].SetUnitPhi(phi);
		      fUnitArray[i].SetUnitDeta(fEMCalGrid->GetDeta());
		      fUnitArray[i].SetUnitDphi(fEMCalGrid->GetDphi());
		    } 
		  else {
		    fTPCGrid->GetEtaPhiFromIndex2(fUnitArray[i].GetUnitID()+1-fNumUnits,phi,eta);
		    fUnitArray[i].SetUnitEta(eta);
		    fUnitArray[i].SetUnitPhi(phi);
		    fUnitArray[i].SetUnitDeta(fTPCGrid->GetDeta());
		    fUnitArray[i].SetUnitDphi(fTPCGrid->GetDphi());
		    
		    if(fDebug>10)
		      {
			if(fUnitArray[i].GetUnitEnergy()!=0.){
			  cout << "(fUnitArray[" << i << "].GetUnitID()+1-fNumUnits : " << 
			    fUnitArray[i].GetUnitID()+1-fNumUnits << endl;
			  cout << "(fUnitArray[" << i << "].GetUnitEnergy() : " << 
			    fUnitArray[i].GetUnitEnergy() << endl;
			  cout << "(fUnitArray[" << i << "].GetUnitEta() : " << 
			    fUnitArray[i].GetUnitEta() << endl;
			  cout << "(fUnitArray[" << i << "].GetUnitPhi() : " << 
			    fUnitArray[i].GetUnitPhi() << endl;
			}
		      } 
		  }
		  fUnitArray[i].SetUnitClusterID(0);
		}//end loop over all units in array (same as all towers in EMCAL)
	    }

	  if(fDebug>20) 
	    {
	      for(Int_t i=0; i<fNumUnits+nElements; i++)
		{
		  if(fUnitArray[i].GetUnitEnergy()!=0.){
		    cout << "######################################################### " << endl;
		    cout << "Final UnitArray filled with energy != 0" << i << endl;
		    cout << "Pointeur UnitArray : " << fUnitArray << " ID : " << 
		      fUnitArray[i].GetUnitID() << " Energy : " << fUnitArray[i].GetUnitEnergy() << 
		      " Eta : " << fUnitArray[i].GetUnitEta() << " Phi : " << fUnitArray[i].GetUnitPhi() << endl;
		  }
		}
	    }

	} // end  if(fOpt==1 || option=="tpc")

      delete[] enT;
      delete[] ptT;
      delete[] etaT;
      delete[] phiT;
      delete[] thetaT;
      
}

//__________________________________________________________
void AliJetFillUnitArrayTracks::GetEtaPhiFromIndex(Int_t index, Float_t &eta, Float_t &phi)
{
  for(Int_t j=0; j<fNphi+1; j++) {
    for(Int_t i=0; i<fNeta+1; i++) {

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












