//-*- Mode: C++ -*-
// $Id: AliHLTJETFastJetFinder.cxx  $
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

/** @file   AliHLTJETFastJetFinder.cxx
    @author Jochen Thaeder
    @date   
    @brief  FastJet finder interface
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__ >= 3
using namespace std;
#endif

//#include "fastjet/AreaDefinition.hh"
//#include "fastjet/JetDefinition.hh"
// get info on how fastjet was configured
//#include "fastjet/config.h"

//#include <TLorentzVector.h>
//#include <TArrayF.h>
//#include <TClonesArray.h>

//#include "AliAODJet.h"

#include "AliHLTJETReader.h"
#include "AliHLTJETJetCuts.h"

#include "AliHLTJETFastJetFinder.h"
#include "AliHLTJETFastJetHeader.h"

//#include<sstream>  // needed for internal io
//#include<vector> 
//#include <cmath> 

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETFastJetFinder)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */
  
// #################################################################################
AliHLTJETFastJetFinder::AliHLTJETFastJetFinder()
  :
  AliJetFinder(),
  fInputVector(NULL),
  fJets(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

// #################################################################################
AliHLTJETFastJetFinder::~AliHLTJETFastJetFinder() { 
  // see header file for class documentation
  
}

/*
 * ---------------------------------------------------------------------------------
 *                                    Initialize
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTJETFastJetFinder::Initialize() {
  // see header file for class documentation

  Int_t iResult = 0;
  
  if ( !fHeader || !fReader ) {
    HLTError("No header or reader set!");
    return -EINPROGRESS;
  }

  // -- Initialize Reader
  AliHLTJETReader *reader = dynamic_cast<AliHLTJETReader*> (fReader);
  
  iResult = reader->Initialize();
  if ( iResult ) {
    HLTError( "Initializing Reader failed!");
    return iResult;
  }
 
  // -- Initialize Header
  AliHLTJETFastJetHeader *header = dynamic_cast<AliHLTJETFastJetHeader*> (fHeader);

  iResult = header->Initialize();
  if ( iResult ) {
    HLTError( "Initializing Header failed!");
    return iResult;
  }

  // -- Set ptr to vector
  fInputVector = reader->GetVector();
  if ( ! fInputVector ) {
    HLTError( "Getting ptr to vector failed!");
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
void AliHLTJETFastJetFinder::Reset() {
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
Bool_t AliHLTJETFastJetFinder::ProcessEvent() {
  // see header file for class documentation

  // -- Pick up jet reader
  AliHLTJETReader *reader = dynamic_cast<AliHLTJETReader*> (fReader);

  // -- Reset
  Reset();

  // -- Fill Vector
  if ( !reader->FillVector() ){
    HLTError("Error filling vector.");
    return kFALSE;  
  }

  // -- Find Jets, fill jets and apply jet cuts 
  if ( FindFastJets() ) {
    HLTError("Error finding jets.");
    return kFALSE;
  }

  return kTRUE;
}

// #################################################################################
Bool_t AliHLTJETFastJetFinder::ProcessHLTEvent() {
  // see header file for class documentation

  // -- Reset
  Reset();

  // -- Find Jets, fill jets and apply jet cuts 
  if ( FindFastJets() ) {
    HLTError("Error finding jets.");
    return kFALSE;
  }

  return kTRUE;
}

/*
 * ---------------------------------------------------------------------------------
 *                                      Helper
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
void AliHLTJETFastJetFinder::PrintJets( vector<fastjet::PseudoJet> &jets,
				  fastjet::ClusterSequenceArea &clust_seq ) {
  // see header file for class documentation

  // -- pick up fastjet header
  AliHLTJETFastJetHeader *header = dynamic_cast<AliHLTJETFastJetHeader*> (fHeader);

  // -- print header info
  TString comment = header->GetComment();
  comment += TString(clust_seq.strategy_string());

  HLTInfo( "--------------------------------------------------------" );
  HLTInfo( "%s", comment.Data() );
  HLTInfo( "--------------------------------------------------------" );
  
  header->PrintParameters();

  // -- print found jets
  HLTInfo( "Number of unclustered particles: %i", clust_seq.unclustered_particles().size() );

  HLTInfo( "Printing inclusive sub jets with pt > %f GeV", header->GetPtMin() );
  HLTInfo( "-------------------------------------------------------" );
  HLTInfo( " ijet   rap      phi        Pt         area  +-   err" );

  for ( UInt_t jetIter = 0; jetIter < jets.size(); jetIter++ ) { 
    
    Double_t area       = clust_seq.area(jets[jetIter]);
    Double_t area_error = clust_seq.area_error(jets[jetIter]);

    HLTInfo( "%5u %9.5f %8.5f %10.3f %8.3f +- %6.3f",
	    jetIter,jets[jetIter].rap(),jets[jetIter].phi(),jets[jetIter].perp(), area, area_error );

  } // for ( Int_t jetIter = 0; jetIter < jets.size(); jetIter++ ) { 

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                             Process - private
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTJETFastJetFinder::FindFastJets() {
  // see header file for class documentation

  // -- pick up fastjet header
  AliHLTJETFastJetHeader *header = dynamic_cast<AliHLTJETFastJetHeader*> (fHeader);
  
  // -- Run the jet clustering with the jet definition from the header
  fastjet::ClusterSequenceArea clust_seq( (*fInputVector), 
					  (*header->GetJetDefinition()),
					  (*header->GetAreaDefinition()) );

  // -- Extract the inclusive jets with pt > ptmin, sorted by pt
  vector<fastjet::PseudoJet> inclusive_jets = clust_seq.inclusive_jets(header->GetPtMin());
  
  // -- Subtract background 
  vector<fastjet::PseudoJet> sub_jets = clust_seq.subtracted_jets((*header->GetRangeDefinition()),
								  header->GetPtMin());  

  // -- Sort jets into increasing pt
  vector<fastjet::PseudoJet> jets = sorted_by_pt(sub_jets);  

  // -- Fill jets in output container
  // -----------------------------------

  // -- Get jet cuts
  AliHLTJETJetCuts* jetCuts = header->GetJetCuts();

  // -- Loop over jet candidates
  for ( UInt_t jetIter = 0; jetIter < jets.size(); jetIter++ ) { 

    // -- Fill AOD jets
    AliAODJet aodjet (jets[jetIter].px(), jets[jetIter].py(), 
		      jets[jetIter].pz(), jets[jetIter].E());

    // -- Apply jet cuts 
    if ( ! jetCuts->IsSelected(&aodjet) )
      continue;

    fJets->AddJet(&aodjet);
    
  } // for ( Int_t jetIter = 0; jetIter < jets.size(); jetIter++ ) { 

#if 0
  // -- Print found jets
  PrintJets( jets, clust_seq );
#endif

  return 0;
}

