//-*- Mode: C++ -*-
// $Id: AliHLTJETConeFinder.cxx  $
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

/** @file   AliHLTJETConeFinder.cxx
    @author Jochen Thaeder
    @date   
    @brief  Jet cone finder
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTJETReader.h"

#include "AliHLTJETConeGrid.h"
#include "AliHLTJETConeFinder.h"
#include "AliHLTJETConeHeader.h"
#include "AliHLTJETConeJetCandidate.h"
#include "AliHLTJETConeEtaPhiCell.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETConeFinder)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */
  
// #################################################################################
AliHLTJETConeFinder::AliHLTJETConeFinder()
  : 
  AliJetFinder(),
  fGrid(NULL),
  fJets(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

// #################################################################################
AliHLTJETConeFinder::~AliHLTJETConeFinder() {
  // see header file for class documentation

}

/*
 * ---------------------------------------------------------------------------------
 *                                    Initialize
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTJETConeFinder::Initialize() {
  // see header file for class documentation

  Int_t iResult = 0;
  
  if ( !fHeader || !fReader ) {
    HLTError("No header or reader set!");
    return -EINPROGRESS;
  }

  // -- Initialize Reader
  AliHLTJETReader *reader = dynamic_cast<AliHLTJETReader*> (fReader);
  if (!reader) {
    HLTError( "Casting Reader failed!");
    return -EINPROGRESS;
  }

  iResult = reader->Initialize();
  if ( iResult ) {
    HLTError( "Initializing Reader failed!");
    return -EINPROGRESS;
  }
 
  // -- Initialize Header
  AliHLTJETConeHeader *header = dynamic_cast<AliHLTJETConeHeader*> (fHeader);
  if (!header) {
    HLTError( "Casting Header failed!");
    return -EINPROGRESS;
  }

  iResult = header->Initialize();
  if ( iResult ) {
    HLTError( "Initializing Header failed!");
    return -EINPROGRESS;
  }

  // -- Set ptr to grid
  fGrid = reader->GetGrid();
  if ( ! fGrid ) {
    HLTError( "Getting ptr to grid failed!");
    return -EINPROGRESS;
  }

  // -- Check ptr to output container
  if ( !fJets ) {
    HLTError( "Ptr to output container not set!");
    return -EINPROGRESS;
  }

  return iResult;
}

// #################################################################################
void AliHLTJETConeFinder::Reset() {
  // see header file for class documentation

  // -- Reset output container
  if (fJets)
    fJets->Reset();
  
  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                      Process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Bool_t AliHLTJETConeFinder::ProcessEvent() {
  // see header file for class documentation

  // -- Pick up jet reader
  AliHLTJETReader *reader = dynamic_cast<AliHLTJETReader*> (fReader);

  // -- Reset
  Reset();

  // -- Fill Grid
  if ( !reader->FillGrid() ){
    HLTError("Error filling grid.");
    return kFALSE;  
  }

  // -- Find Leading
  if ( FindConeLeading()  ) {
    HLTError("Error finding leading.");
    return kFALSE;
  }

  // -- Find Jets 
  if ( FindConeJets()  ) {
    HLTError("Error finding jets.");
    return kFALSE;
  }

  // -- Fill Jets and apply jet cuts
  if ( FillConeJets()  ) {
    HLTError("Error filling jets.");
    return kFALSE;
  }

  return kTRUE;
}

// #################################################################################
Bool_t AliHLTJETConeFinder::ProcessHLTEvent() {
  // see header file for class documentation

  // -- Reset
  Reset();

  // -- Find Leading
  if ( FindConeLeading()  ) {
    HLTError("Error finding leading.");
    return kFALSE;
  }

  // -- Find Jets 
  if ( FindConeJets()  ) {
    HLTError("Error finding jets.");
    return kFALSE;
  }

  // -- Fill Jets and apply jet cuts
  if ( FillConeJets()  ) {
    HLTError("Error filling jets.");
    return kFALSE;
  }

  return kTRUE;
}

/*
 * ---------------------------------------------------------------------------------
 *                             Process - private
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTJETConeFinder::FindConeLeading() {
  // see header file for class documentation

  // -- Pick up jet reader
  AliHLTJETReader *reader = dynamic_cast<AliHLTJETReader*> (fReader);

  // -- Pick up jet canidates
  TClonesArray* jetCandidates = reader->GetJetCandidates();

  // -- Check for more than 1 jet candidate
  if ( reader->GetNJetCandidates() > 1 ) {
    
    // -- Sort jet candidates with seed pt 
    jetCandidates->Sort();
    
    // -- Use leading seed only
    //    Keep index 0, remove the others
    if ( (dynamic_cast<AliHLTJETConeHeader*> (fHeader))->GetUseLeading() ) {
      
      for ( Int_t iter = 1; iter < reader->GetNJetCandidates(); iter++ )
	jetCandidates->RemoveAt(iter);
      
      reader->SetNJetCandidates(1);
    }

  } // if ( reader->GetNJetCandidates() > 1 ) {

  // -- Resize the seed TClonesArray
  jetCandidates->Compress();

  return 0;
}

// #################################################################################
Int_t AliHLTJETConeFinder::FindConeJets() {
  // see header file for class documentation

  Int_t iResult = 0;

  // -- Pick up jet reader
  AliHLTJETReader *reader = dynamic_cast<AliHLTJETReader*> (fReader);
   
  // -- Pick up jet canidates
  TClonesArray* jetCandidates = reader->GetJetCandidates();

  // -- Loop over jet candidates
  for ( Int_t iter = 0; iter < reader->GetNJetCandidates(); iter++ ) {
    
    AliHLTJETConeJetCandidate* jet = reinterpret_cast<AliHLTJETConeJetCandidate*> ((*jetCandidates)[iter]);
    
    // -- Set iterator for cells around seed
    fGrid->SetCellIter( jet->GetSeedEtaIdx(), jet->GetSeedPhiIdx() ); 

    Int_t cellIdx = 0;
    
    // -- Loop over cells around ssed
    while ( (cellIdx = fGrid->NextCell() ) >= 0 && !iResult) {

      AliHLTJETConeEtaPhiCell* cell = NULL;
      if ( ! (cell = reinterpret_cast<AliHLTJETConeEtaPhiCell*>(fGrid->UncheckedAt(cellIdx))) )
	continue;
      
      if ( ( iResult = jet->AddCell(cell) ) ) {
	HLTError( "Error adding cell %d to jet candiate %d", cellIdx, iter);
	continue;
      }
    } // while ( (cellIdx = fGrid->NextCell() ) >= 0 && !iResult) {
    
  } // for ( Int_t iter = 0; iter < reader->GetNJetCandidates(); iter++ ) {
  
  return iResult;
}

// #################################################################################
Int_t AliHLTJETConeFinder::FillConeJets() {
  // see header file for class documentation

  Int_t iResult = 0;

  // -- Pick up jet reader
  AliHLTJETReader *reader = dynamic_cast<AliHLTJETReader*> (fReader);

  // -- Pick up jet header
  AliHLTJETConeHeader *header = dynamic_cast<AliHLTJETConeHeader*> (fHeader);
  
  // -- Get jet canidates
  TClonesArray* jetCandidates = reader->GetJetCandidates();
  
  // -- Get jet cuts
  AliHLTJETJetCuts* jetCuts = header->GetJetCuts();

  // -- Loop over jet candidates
  for ( Int_t iter = 0; iter < reader->GetNJetCandidates(); iter++ ) {
    
    AliHLTJETConeJetCandidate* jet = reinterpret_cast<AliHLTJETConeJetCandidate*> ((*jetCandidates)[iter]);
  
    // -- Apply jet cuts 
    if ( ! jetCuts->IsSelected(jet) )
      continue;

    // -- Add jet as AliAODJet
    fJets->AddJet(jet->GetEta(), jet->GetPhi(), jet->GetPt(), jet->GetEt());
    
  } // for ( Int_t iter = 0; iter < reader->GetNJetCandidates(); iter++ ) {
  
  // xxx  HLTDebug( "Added %d jets", fJets->GetNAODJets());
  HLTInfo( "Added %d jets", fJets->GetNAODJets());

  return iResult;
}
