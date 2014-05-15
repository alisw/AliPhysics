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

/*****************************************************************
  AliFlowEventStar: Event container for flow analysis

  origin:   Mikolaj Krzewicki  (mikolaj.krzewicki@cern.ch)
*****************************************************************/

#include "Riostream.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliStarTrack.h"
#include "AliStarTrackCuts.h"
#include "AliStarEvent.h"
#include "AliFlowEventStar.h"

using std::cout;
using std::endl;
ClassImp(AliFlowEventStar)

//-----------------------------------------------------------------------

AliFlowEventStar::AliFlowEventStar():
  AliFlowEventSimple()
{
  //ctor
  cout << "AliFlowEventStar: Default constructor to be used only by root for io" << endl;
}

//-----------------------------------------------------------------------
AliFlowEventStar::AliFlowEventStar(const AliFlowEventStar& event):
  AliFlowEventSimple(event)
{
  //cpy ctor
}

//-----------------------------------------------------------------------
AliFlowEventStar& AliFlowEventStar::operator=( const AliFlowEventStar& event )
{
  //assignment operator
  AliFlowEventSimple::operator=(event);
  return *this;
}

//-----------------------------------------------------------------------
AliFlowEventStar::AliFlowEventStar( const AliStarEvent* starevent,
                                    const AliStarTrackCuts* rpCuts,
                                    const AliStarTrackCuts* poiCuts ):
  AliFlowEventSimple(starevent->GetNumberOfTracks())
{
  //construct from a star event
  SetReferenceMultiplicity(starevent->GetRefMult());
  for (Int_t i=0; i<starevent->GetNumberOfTracks(); i++)
  {
    const AliStarTrack* startrack = starevent->GetTrack(i);
    if (!startrack) continue;
    AliFlowTrackSimple* flowtrack = new AliFlowTrackSimple();
    flowtrack->SetPhi(startrack->GetPhi());
    flowtrack->SetEta(startrack->GetEta());
    flowtrack->SetPt(startrack->GetPt());
    flowtrack->SetCharge(startrack->GetCharge());
    if (rpCuts)
    {
      Bool_t pass = rpCuts->PassesCuts(startrack);
      flowtrack->TagRP(pass); //tag RPs
      if (pass) IncrementNumberOfPOIs(0);
    }
    if (poiCuts)
    {
      flowtrack->TagPOI(poiCuts->PassesCuts(startrack)); //tag POIs
    }
    AddTrack(flowtrack);
  }
}
