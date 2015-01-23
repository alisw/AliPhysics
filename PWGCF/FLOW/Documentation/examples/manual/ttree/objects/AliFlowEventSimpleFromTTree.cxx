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

/*
  AliFlowEventSimpleFromTTree: fill a flow event from a TTree
  example class, to be expanded by the user
  see macros in PWGCF/FLOW/Documentation/examples/manual/ttree/macros
  for info on how to use this class

  author:   Redmer Alexander Bertens (rbertens@cern.ch)

  based on code from Mikolaj, see
  $ALICE_PHYSICS/PWG/FLOW/Base/AliFlowEventStar.cxx
*/

// root includes
#include "TClonesArray.h"

// aliroot includes
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowTrackSimpleCuts.h"

// custom includes
#include "AliFlowEventSimpleFromTTree.h"
#include "AliFlowTTreeEvent.h"
#include "AliFlowTTreeTrack.h"

ClassImp(AliFlowEventSimpleFromTTree)

//------------------------------------------------------------------------
AliFlowEventSimpleFromTTree::AliFlowEventSimpleFromTTree():
  AliFlowEventSimple()
{
  // constructor for root I/O
}
//------------------------------------------------------------------------
AliFlowEventSimpleFromTTree::AliFlowEventSimpleFromTTree(const AliFlowEventSimpleFromTTree& event):
  AliFlowEventSimple(event)
{
  // copy constructor
}
//------------------------------------------------------------------------
AliFlowEventSimpleFromTTree& AliFlowEventSimpleFromTTree::operator=( const AliFlowEventSimpleFromTTree& event )
{
  // assignment operator
  AliFlowEventSimple::operator=(event);
  return *this;
}
//------------------------------------------------------------------------
AliFlowEventSimpleFromTTree::AliFlowEventSimpleFromTTree( 
        const AliFlowTTreeEvent* event,
        const TClonesArray* array,
        const AliFlowTrackSimpleCuts* rpCuts,
        const AliFlowTrackSimpleCuts* poiCuts ) :
  AliFlowEventSimple(array->GetEntries())
{
    // construct flow event from generic input
    SetReferenceMultiplicity(array->GetEntries());
    
    // place here for some event selection
    if(event && PassesCuts(event)) {
        // track loop
        for (Int_t i(0); i < array->GetEntries(); i++) {
            AliFlowTTreeTrack* track = static_cast<AliFlowTTreeTrack*>(array->At(i));
            if (!track) continue;
            AliFlowTrackSimple* flowtrack = new AliFlowTrackSimple();
            flowtrack->SetPhi(track->GetPhi());
            flowtrack->SetEta(track->GetEta());
            flowtrack->SetPt(track->GetPt());
            flowtrack->SetCharge(track->GetCharge());
            if (rpCuts) {
                Bool_t pass = rpCuts->PassesCuts(flowtrack);
                flowtrack->TagRP(pass); //tag RPs
                if (pass) IncrementNumberOfPOIs(0);
            }
            if (poiCuts) {
                flowtrack->TagPOI(poiCuts->PassesCuts(flowtrack)); //tag POIs
            }
            AddTrack(flowtrack);
        }
    }
}
//------------------------------------------------------------------------
Bool_t AliFlowEventSimpleFromTTree::PassesCuts(const AliFlowTTreeEvent* event)
{
    // example of event cuts
    // of course in practice hardcoding cuts isn't very good practice ...
    if(event->GetV0M() > 40 || event->GetV0M() < 20) return kFALSE;
    return kTRUE;
}
//------------------------------------------------------------------------

