//-*- Mode: C++ -*-
// $Id: AliHLTJETReaderHeader.cxx  $
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTJETReader.cxx
    @author Jochen Thaeder
    @date   
    @brief  Reader for jet finder
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__ >= 3
using namespace std;
#endif

#include "TLorentzVector.h"
#include "TParticle.h"
#include "TParticlePDG.h"

#include "AliHLTJETReader.h"

#include "AliHLTJETConeJetCandidate.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETReader)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */
  
// #################################################################################
AliHLTJETReader::AliHLTJETReader()
  : 
  AliJetReader(),
  fESD(NULL), 
  fMC(NULL),
  fAOD(NULL),
#ifdef HAVE_FASTJET
  fMomentumVector( new vector<fastjet::PseudoJet> ),
#endif
  fGrid(NULL),
  fNJetCandidates(0),
  fJetCandidates(NULL),
  fSeedCuts(NULL),
  fTrackCuts(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

// #################################################################################
AliHLTJETReader::~AliHLTJETReader() {
  // see header file for class documentation

#ifdef HAVE_FASTJET
  if ( fMomentumVector )
    delete fMomentumVector;
  fMomentumVector = NULL;
#endif

  if ( fGrid )
    delete fGrid;
  fGrid = NULL;

  if ( fJetCandidates ) {
    fJetCandidates->Clear();
    delete fJetCandidates;
  }
  fJetCandidates = NULL;

}

// #################################################################################
Int_t AliHLTJETReader::Initialize() {
  // see header file for class documentation

  Int_t iResult = 0;
  AliHLTJETReaderHeader* readerHeader = NULL;

  HLTInfo(" -= AliHLTJETReader =- " );

  // -- Initialize reader header
  // -----------------------------
  if ( fReaderHeader ) {
    readerHeader = GetReaderHeader();

    iResult = readerHeader->Initialize();
    if ( iResult )
      HLTError("Error initializing HLT jet reader header");
  }
  else {
    HLTError("Reader Header not present");
    iResult = -EINPROGRESS;
  }
  
  // -- Initialize grid
  // --------------------
  if ( ! iResult ) {
   
    if ( fGrid )
      delete fGrid;

    if ( ! (fGrid = new AliHLTJETConeGrid()) ) {
      HLTError("Error instanciating grid.");
      iResult = -EINPROGRESS;
    }
  }
  
  if ( ! iResult ) {
      fGrid->SetEtaRange(   readerHeader->GetFiducialEtaMin(),
			    readerHeader->GetFiducialEtaMax(),
			    readerHeader->GetGridEtaRange() );

      fGrid->SetPhiRange(   readerHeader->GetFiducialPhiMin(),
			    readerHeader->GetFiducialPhiMax(),
			    readerHeader->GetGridPhiRange() );
      
      fGrid->SetBinning(    readerHeader->GetGridEtaBinning(),
			    readerHeader->GetGridEtaBinning() );

      fGrid->SetConeRadius( readerHeader->GetConeRadius() );

      iResult = fGrid->Initialize();
  }
 
  // -- Initialize jet candidates
  // ------------------------------
  if ( ! iResult ) {
    fJetCandidates = new TClonesArray("AliHLTJETConeJetCandidate", 30);
    if ( ! fJetCandidates) {
      HLTError("Error instanciating jet candidates.");
      iResult = -EINPROGRESS;
    }
  }
 
  // -- Initialize cuts
  // --------------------
 
  // -- Seed cuts
  if ( ! iResult ) {
    fSeedCuts = readerHeader->GetSeedCuts();
    if ( ! fSeedCuts ) {
      HLTError("Error getting ptr to seed cuts.");
      iResult = -EINPROGRESS;
    }
    else {
      HLTInfo(" -= SeedCuts =- " );
    }
  }

  // -- Track cuts
  if ( ! iResult ) {
    fTrackCuts = readerHeader->GetTrackCuts();
    if ( ! fTrackCuts ) {
      HLTError("Error getting ptr to track cuts.");
      iResult = -EINPROGRESS;
    }
  }

  // -- Final check
  // ----------------
  if ( iResult )
    HLTError("Error initializing HLT jet reader");
 
  return iResult;
}

//#################################################################################
void AliHLTJETReader::ResetEvent() {
  // see header file for class documentation

  // -- clear grid
  fGrid->Reset();

  // -- clear jet candidates
  fJetCandidates->Clear();

  fNJetCandidates = 0;

  return;  
}

/*
 * ---------------------------------------------------------------------------------
 *                            Fastjet Reader functionality
 * ---------------------------------------------------------------------------------
 */
#ifdef HAVE_FASTJET
// #################################################################################
Bool_t AliHLTJETReader::FillMomentumArrayFast() {
  // see header file for class documentation

  Bool_t bResult = kFALSE;

  if ( fESD )
    bResult = FillMomentumArrayFastESD();
  else if ( fMC )
    bResult = FillMomentumArrayFastMC();
  else if ( fAOD )
    bResult = FillMomentumArrayFastAOD();
  
  return bResult;
}

// #################################################################################
Bool_t AliHLTJETReader::FillMomentumArrayFastMC() {
  // see header file for class documentation

  if ( ! fMC ) {
    HLTError( "No MC Event present!" );
    return kFALSE;
  }

  // -- Clear input vector
  if ( fMomentumVector )
    fMomentumVector->clear();
      
  Int_t nTracks = 0;

  TParticle* particle = NULL;

  // -- Loop over particles
  // ------------------------
  while ( (particle = fMC->NextParticle() ) ) {
    
    // -- Apply cuts 
    if ( ! fTrackCuts->IsSelected(particle) )
      continue;
    
    // -- Create PseudoJet object      
    fastjet::PseudoJet part( particle->Px(), particle->Py(), 
			     particle->Pz(), particle->Energy() ); 
    
    // -- label the particle into Fastjet algortihm
    part.set_user_index( fMC->GetIndex() ); 

    // -- Add to input_particles vector  
    fMomentumVector->push_back(part);  

    nTracks++;    
    
  } // while ( (particle = fMC->NextParticle() ) ) {
  
  HLTInfo(" Number of selected tracks %d \n", nTracks);

  return kTRUE;
}

// #################################################################################
Bool_t AliHLTJETReader::FillMomentumArrayFastESD() {
  // see header file for class documentation

  return kTRUE;
}

// #################################################################################
Bool_t AliHLTJETReader::FillMomentumArrayFastAOD() {
  // see header file for class documentation

  return kFALSE;
}
#endif

/*
 * ---------------------------------------------------------------------------------
 *                               Grid functionality
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Bool_t AliHLTJETReader::FillGrid() {
  // see header file for class documentation

  Bool_t bResult = kFALSE;

  if ( fESD )
    bResult = FillGridESD();
  else if ( fMC )
    bResult = FillGridMC();
  else if ( fAOD )
    bResult = FillGridAOD();
  
  return bResult;
}

// #################################################################################
Bool_t AliHLTJETReader::FillGridMC() {
  // see header file for class documentation

  if ( ! fMC ) {
    HLTError( "No MC Event present!" );
    return kFALSE;
  }

  AliHLTJETReaderHeader* readerHeader = GetReaderHeader();

  Int_t nTracks = 0;
  TParticle* particle = NULL;

  // -- Loop over particles
  // ------------------------
  while ( ( particle = fMC->NextParticle() ) ) {
    
    // -- Apply track cuts 
    if ( ! fTrackCuts->IsSelected(particle) )
      continue;

    const Float_t aEtaPhi[]  = { particle->Eta(), particle->Phi(), particle->Pt() }; 
          Int_t   aGridIdx[] = { -1, -1, -1, -1, -1 };

    fGrid->FillTrack(particle, aEtaPhi, aGridIdx);

    nTracks++;    

    // -- Apply seed cuts 
    if ( ! fSeedCuts->IsSelected(particle) )
      continue;  

    // -- Add Seed
    AddSeed(aEtaPhi, const_cast<const Int_t*> (aGridIdx),
	    readerHeader->GetConeRadius());
  
  } // while ( (particle = fMC->NextParticle() ) ) {
  
  HLTDebug(" Number of selected tracks %d", nTracks);
  HLTDebug(" Number of seeds %d", fNJetCandidates);

  return kTRUE;
}

// #################################################################################
Bool_t AliHLTJETReader::FillGridESD() {
  // see header file for class documentation

  Bool_t bResult = kTRUE;

  if ( ! fESD ) {
    HLTError( "No ESD Event present!" );
    return kFALSE;
  }

  AliHLTJETReaderHeader* readerHeader = GetReaderHeader();

  Int_t nTracks = 0;

  // -- Loop over particles
  // ------------------------
  for ( Int_t iter = 0; iter < fESD->GetNumberOfTracks() && !bResult; iter++ ) {

    AliESDtrack* esdTrack = fESD->GetTrack(iter);
    if ( ! esdTrack ) {
      HLTError("Could not read ESD track %d from %d\n", iter, fESD->GetNumberOfTracks() );
      bResult = kFALSE;
      continue;
    }

    // -- Apply track cuts 
    if ( ! fTrackCuts->IsSelected(esdTrack) )
      continue;

    const Float_t aEtaPhi[]  = { esdTrack->Eta(), esdTrack->Phi(), esdTrack->Pt() }; 
          Int_t   aGridIdx[] = { -1, -1, -1, -1, -1 };

    // -- Fill grid
	  fGrid->FillTrack(esdTrack, aEtaPhi, aGridIdx);

    nTracks++;    

    // -- Apply seed cuts 
    if ( ! fSeedCuts->IsSelected(esdTrack) )
      continue;

    // -- Add Seed
    AddSeed(aEtaPhi, const_cast<const Int_t*> (aGridIdx),
	    readerHeader->GetConeRadius());
  
  } // for ( Int_t iter = 0; iter < fESD->GetNumberOfTracks() && !iResult; iter++ ) {
  
  HLTDebug(" Number of selected tracks %d", nTracks);
  HLTDebug(" Number of seeds %d", fNJetCandidates);

  return bResult;
}

// #################################################################################
Bool_t AliHLTJETReader::FillGridAOD() {
  // see header file for class documentation

  return kFALSE;
}

/*
 * ---------------------------------------------------------------------------------
 *                                     Setter
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
void AliHLTJETReader::SetInputEvent(TObject* esd, TObject* aod, TObject* mc) {
  // see header file for class documentation

  if ( esd )
    fESD = dynamic_cast<AliESDEvent*> (esd);
  else if ( aod )
    fAOD = dynamic_cast<AliAODEvent*> (aod);
  else if ( mc )
    fMC = dynamic_cast<AliHLTMCEvent*> (mc);
  
  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                     Seeds
 * ---------------------------------------------------------------------------------
 */

//#################################################################################
void AliHLTJETReader::AddSeed( const Float_t* aEtaPhi, const Int_t* aGridIdx, 
			       const Float_t coneRadius ) {
  // see header file for class documentation

  Bool_t useWholeCell = kTRUE ; // XXXXXXXXXXXXXXXXx get reader header finder type balhh
  useWholeCell = kFALSE ;
  // -- Add track / particle

  new( (*fJetCandidates) [fNJetCandidates] ) AliHLTJETConeJetCandidate( aEtaPhi, 
									aGridIdx,
									coneRadius,
									useWholeCell );
  fNJetCandidates++;

  HLTDebug("Added Seed Pt=%f, Eta=%f, Phi=%f", aEtaPhi[kIdxPt], 
	  aEtaPhi[kIdxEta], aEtaPhi[kIdxPhi] );

  return;
}
