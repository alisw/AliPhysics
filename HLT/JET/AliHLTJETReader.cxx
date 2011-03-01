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
#include "TDatabasePDG.h"

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
  fHLTMC(NULL),
  fAOD(NULL),
#ifdef HAVE_FASTJET
  fMomentumVector(NULL),
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

  HLTInfo(" -= AliHLTJETReader =- ");

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
  
  // -- Initialize Algorithms
  // --------------------------
  if ( !iResult ) {
    if ( readerHeader->GetJetAlgorithm() >= AliHLTJETBase::kFFSCSquareCell )
      iResult = InitializeFFSC();
    else {
#ifdef HAVE_FASTJET
      iResult = InitializeFastjet();
#else  
      HLTError("Error FastJet not present.");
      iResult = -EINPROGRESS;
#endif
    }
  }

  // -- Get ptr to cuts from reader
  // --------------------------------
  if ( !iResult ) {
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

  // -- Reset FFSC algorithms
  // --------------------------
  if (  GetReaderHeader()->GetJetAlgorithm() >= AliHLTJETBase::kFFSCSquareCell ) {
    // -- clear grid
    fGrid->Reset();
    
    // -- clear jet candidates
    fJetCandidates->Clear();
    
    fNJetCandidates = 0;
  }

  // -- Reset for FastJet algorithms
  // ---------------------------------
  else {
#ifdef HAVE_FASTJET
    // -- Clear input vector
    if ( fMomentumVector )
      fMomentumVector->clear();
#endif
  }

  return;  
}

/*
 * ---------------------------------------------------------------------------------
 *                                     Setter
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
void AliHLTJETReader::SetInputEvent(const TObject* esd, const TObject* aod, const TObject* mc) {
  // see header file for class documentation

  //  Needs "useMC" flag for running in analysis task only

  AliHLTJETReaderHeader* readerHeader = GetReaderHeader();

  // -- Fill ESD
  if ( esd && !readerHeader->GetUseMC() )
    fESD = dynamic_cast<AliESDEvent*> (const_cast<TObject*>(esd));

  // -- Fill AOD
  else if ( aod && !readerHeader->GetUseMC() )
    fAOD = dynamic_cast<AliAODEvent*> (const_cast<TObject*>(aod));

  // -- Fill MC
  else if ( mc && readerHeader->GetUseMC() ) {

    // -- if off-line MC event, 
    if ( !strcmp (mc->ClassName(),"AliMCEvent") )
      fMC = dynamic_cast<AliMCEvent*> (const_cast<TObject*>(mc));
    // -- if on-line MC event
    else 
      fHLTMC = dynamic_cast<AliHLTMCEvent*> (const_cast<TObject*>(mc));
  }
  
  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                            Fastjet Reader functionality
 * ---------------------------------------------------------------------------------
 */
#ifdef HAVE_FASTJET
// #################################################################################
Bool_t AliHLTJETReader::FillVector() {
  // see header file for class documentation

  Bool_t bResult = kFALSE;

  if ( fESD )
    bResult = FillVectorESD();
  else if ( fMC )
    bResult = FillVectorMC();
  else if ( fHLTMC )
    bResult = FillVectorHLTMC();
  else if ( fAOD )
    bResult = FillVectorAOD();
  
  return bResult;
}

// #################################################################################
Bool_t AliHLTJETReader::FillVectorMC() {
  // see header file for class documentation

  Bool_t bResult = kTRUE; 

  if ( ! fMC ) {
    HLTError( "No MC Event present!" );
    return kFALSE;
  }
  
  // -- Reset Event
  ResetEvent();
   
  Int_t nTracks = 0;

  AliStack* stack = fMC->Stack();

  for (Int_t iterStack = 0; iterStack < stack->GetNtrack() && bResult; iterStack++) {

    TParticle *particle = stack->Particle(iterStack);
    if ( !particle) {
      HLTError( "Error reading particle %i out of %i", iterStack, stack->GetNtrack() );
      bResult = kFALSE;
      continue;
    }

    // ------------------------------
    // -- Basic cuts on MC particle      --> To be done better XXX
    // ------------------------------

    // -- primary
    if ( !(stack->IsPhysicalPrimary(iterStack)) )
      continue;
    
    // -- final state
    if ( particle->GetNDaughters() != 0 )
      continue;

    // -- particle in DB
    TParticlePDG * particlePDG = particle->GetPDG();
    if ( ! particlePDG ) {
      particlePDG = TDatabasePDG::Instance()->GetParticle( particle->GetPdgCode() );

      if ( ! particlePDG ) {
	HLTError("Particle %i not in PDG database", particle->GetPdgCode() );
	bResult = kFALSE;
	continue;
      }
    }

    // ------------------------
    // -- Standard track cuts
    // ------------------------

    // -- Apply track cuts     
    if ( ! fTrackCuts->IsSelected(particle) )
      continue;
    
    // -- Create PseudoJet object      
    fastjet::PseudoJet part( particle->Px(), particle->Py(), 
			     particle->Pz(), particle->Energy() ); 
    
    // -- label the particle into Fastjet algortihm
    part.set_user_index( iterStack ); 

#ifdef HAVE_FASTJET
    // -- Add to input_particles vector  
    fMomentumVector->push_back(part);  
#endif

    nTracks++;    
  } // for (Int_t iterStack = 0; iterStack < stack->GetNtrack() && !bResult; iterStack++) {

  HLTDebug(" Number of selected tracks %d", nTracks);

  return kTRUE;
}

// #################################################################################
Bool_t AliHLTJETReader::FillVectorHLTMC() {
  // see header file for class documentation

  if ( ! fHLTMC ) {
    HLTError( "No HLT MC Event present!" );
    return kFALSE;
  }
     
  // -- Reset Event
  ResetEvent();

  Int_t nTracks = 0;

  TParticle* particle = NULL;

  // -- Loop over particles
  // ------------------------
  while ( (particle = fHLTMC->NextParticle() ) ) {
    
    // -- Apply track cuts 
    if ( ! fTrackCuts->IsSelected(particle) )
      continue;
    
    // -- Create PseudoJet object      
    fastjet::PseudoJet part( particle->Px(), particle->Py(), 
			     particle->Pz(), particle->Energy() ); 
    
    // -- label the particle into Fastjet algortihm
    part.set_user_index( fHLTMC->GetIndex() ); 

#ifdef HAVE_FASTJET
    // -- Add to input_particles vector  
    fMomentumVector->push_back(part);  
#endif

    nTracks++;    
    
  } // while ( (particle = fHLTMC->NextParticle() ) ) {
  
  HLTInfo(" Number of selected tracks %d", nTracks);

  return kTRUE;
}

// #################################################################################
Bool_t AliHLTJETReader::FillVectorESD() {
  // see header file for class documentation

  Bool_t bResult = kTRUE;

  if ( ! fESD ) {
    HLTError( "No ESD Event present!" );
    return kFALSE;
  }

  // -- Reset Event
  ResetEvent();

  Int_t nTracks = 0;

  // -- Loop over particles
  // ------------------------
  for ( Int_t iter = 0; iter < fESD->GetNumberOfTracks() && bResult; iter++ ) {

    AliESDtrack* esdTrack = fESD->GetTrack(iter);
    if ( ! esdTrack ) {
      HLTError("Could not read ESD track %d from %d", iter, fESD->GetNumberOfTracks() );
      bResult = kFALSE;
      continue;
    }

    // -- Apply track cuts 
    if ( ! fTrackCuts->IsSelected(esdTrack) )
      continue;

    // -- Create PseudoJet object      
    fastjet::PseudoJet part( esdTrack->Px(), esdTrack->Py(), 
			     esdTrack->Pz(), esdTrack->E() ); 
    
    // -- label the particle into Fastjet algortihm
    part.set_user_index( iter ); 

#ifdef HAVE_FASTJET
    // -- Add to input_particles vector  
    fMomentumVector->push_back(part);  
#endif

    nTracks++; 

  } // for ( Int_t iter = 0; iter < fESD->GetNumberOfTracks() && !iResult; iter++ ) {
  
  HLTInfo(" Number of selected tracks %d", nTracks);

  return bResult;
}

// #################################################################################
Bool_t AliHLTJETReader::FillVectorAOD() {
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
  else if ( fHLTMC )
    bResult = FillGridHLTMC();
  else if ( fAOD )
    bResult = FillGridAOD();

  return bResult;
}

// #################################################################################
Bool_t AliHLTJETReader::FillGridMC() {
  // see header file for class documentation

  Bool_t bResult = kTRUE; 

  if ( ! fMC ) {
    HLTError( "No MC Event present!" );
    return kFALSE;
  }

  // -- Reset Event
  ResetEvent();

  AliHLTJETReaderHeader* readerHeader = GetReaderHeader();

  Int_t nTracks = 0;

  AliStack* stack = fMC->Stack();

  for (Int_t iterStack = 0; iterStack < stack->GetNtrack() && bResult; iterStack++) {

    TParticle *particle = stack->Particle(iterStack);
    if ( !particle) {
      HLTError( "Error reading particle %i out of %i", iterStack, stack->GetNtrack() );
      bResult = kFALSE;
      continue;
    }

    // ------------------------------
    // -- Basic cuts on MC particle      --> To be done better XXX
    // ------------------------------

    // -- primary
    if ( !(stack->IsPhysicalPrimary(iterStack)) )
      continue;

    // -- final state
    if ( particle->GetNDaughters() != 0 )
      continue;

    // -- particle in DB
    TParticlePDG * particlePDG = particle->GetPDG();
    if ( ! particlePDG ) {
      particlePDG = TDatabasePDG::Instance()->GetParticle( particle->GetPdgCode() );

      if ( ! particlePDG ) {
	HLTError("Particle %i not in PDG database", particle->GetPdgCode() );
	bResult = kFALSE;
	continue;
      }
    }

    // ------------------------
    // -- Standard track cuts
    // ------------------------

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
    AddSeed(aEtaPhi, const_cast<const Int_t*> (aGridIdx), readerHeader->GetConeRadius());
    
  } // for (Int_t iterStack = 0; iterStack < stack->GetNtrack() && !bResult; iterStack++) {

  HLTInfo(" Number of selected tracks %d", nTracks);
  HLTInfo(" Number of seeds %d", fNJetCandidates);

  return kTRUE;
}

// #################################################################################
Bool_t AliHLTJETReader::FillGridHLTMC() {
  // see header file for class documentation

  if ( ! fHLTMC ) {
    HLTError( "No HLT MC Event present!" );
    return kFALSE;
  }

  // -- Reset Event
  ResetEvent();

  AliHLTJETReaderHeader* readerHeader = GetReaderHeader();

  Int_t nTracks = 0;
  TParticle* particle = NULL;

  // -- Loop over particles
  // ------------------------
  while ( ( particle = fHLTMC->NextParticle() ) ) {
    
    //    HLTError("=== nTracks %d ===",nTracks);

    // -- Apply track cuts 
    if ( ! fTrackCuts->IsSelected(particle) )
      continue;

    const Float_t aEtaPhi[]  = { particle->Eta(), particle->Phi(), particle->Pt() }; 
          Int_t   aGridIdx[] = { -1, -1, -1, -1, -1 };

    // -- Fill grid
    fGrid->FillTrack(particle, aEtaPhi, aGridIdx);

    nTracks++;    

    // -- Apply seed cuts 
    if ( ! fSeedCuts->IsSelected(particle) )
      continue;  

    // -- Add Seed
    AddSeed(aEtaPhi, const_cast<const Int_t*> (aGridIdx),
	    readerHeader->GetConeRadius());
  
  } // while ( (particle = fHLTMC->NextParticle() ) ) {
  
  HLTInfo(" Number of selected tracks %d", nTracks);
  HLTInfo(" Number of seeds %d", fNJetCandidates);

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

  // -- Reset Event
  ResetEvent();

  AliHLTJETReaderHeader* readerHeader = GetReaderHeader();

  Int_t nTracks = 0;

  // -- Loop over particles
  // ------------------------
  //  for ( Int_t iter = 0; iter < fESD->GetNumberOfTracks() && !bResult; iter++ ) { JMT Coverity
  for ( Int_t iter = 0; iter < fESD->GetNumberOfTracks() && bResult; iter++ ) {

    AliESDtrack* esdTrack = fESD->GetTrack(iter);
    if ( ! esdTrack ) {
      HLTError("Could not read ESD track %d from %d", iter, fESD->GetNumberOfTracks() );
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
  
  HLTInfo(" Number of selected tracks %d", nTracks);
  HLTInfo(" Number of seeds %d", fNJetCandidates);

  return bResult;
}

// #################################################################################
Bool_t AliHLTJETReader::FillGridAOD() {
  // see header file for class documentation

  return kFALSE;
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

  Bool_t useWholeCell = kTRUE;

  if ( GetReaderHeader()->GetJetAlgorithm() == AliHLTJETBase::kFFSCRadiusCell )
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

/*
 * ---------------------------------------------------------------------------------
 *                         Initialize - private
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTJETReader::InitializeFFSC() {
  // see header file for class documentation

  Int_t iResult = 0;
  AliHLTJETReaderHeader* readerHeader = GetReaderHeader();

  // -- Initialize grid
  // --------------------
  if ( fGrid )
    delete fGrid;
  
  if ( ! (fGrid = new AliHLTJETConeGrid()) ) {
    HLTError("Error instanciating grid.");
    iResult = -EINPROGRESS;
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

  // -- Get ptr to cuts from reader
  // --------------------------------
 
  // -- Seed cuts
  if ( ! iResult ) {
    fSeedCuts = readerHeader->GetSeedCuts();
    if ( ! fSeedCuts ) {
      HLTError("Error getting ptr to seed cuts.");
      iResult = -EINPROGRESS;
    }
  }
 
  return iResult;
}

#ifdef HAVE_FASTJET
// #################################################################################
Int_t AliHLTJETReader::InitializeFastjet() {
  // see header file for class documentation

  Int_t iResult = 0;

  // -- Initialize Vector
  // ----------------------
  if ( fMomentumVector )
    delete fMomentumVector;
  
  if ( ! (fMomentumVector = new vector<fastjet::PseudoJet>) ) {
    HLTError("Error instanciating momentum vector.");
    iResult = -EINPROGRESS;
  }
 
  return iResult;
}
#endif
