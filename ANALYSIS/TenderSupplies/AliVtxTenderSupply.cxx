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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Redo primary vertex on the fly, using the diamond constraint              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <AliESDEvent.h>
#include <AliESDInputHandler.h>
#include <AliVertexerTracks.h>
#include <AliTender.h>
#include <AliCDBId.h>
#include <AliCDBManager.h>
#include <AliCDBEntry.h>

#include "AliVtxTenderSupply.h"

AliVtxTenderSupply::AliVtxTenderSupply() :
  AliTenderSupply(),
  fDiamond(0x0)
{
  //
  // default ctor
  //
}

//_____________________________________________________
AliVtxTenderSupply::AliVtxTenderSupply(const char *name, const AliTender *tender) :
  AliTenderSupply(name,tender),
  fDiamond(0x0)
{
  //
  // named ctor
  //
}

//_____________________________________________________
void AliVtxTenderSupply::ProcessEvent()
{
  //
  // Recalculate the Vertex with constraint
  //

  AliESDEvent *event=fTender->GetEvent();
  if (!event) return;

  //

  if (fTender->RunChanged()){
    fDiamond=0x0;
    AliCDBEntry *meanVertex=fTender->GetCDBManager()->Get("GRP/Calib/MeanVertex",fTender->GetRun());
    if (!meanVertex) {
      AliError("No new MeanVertex entry found");
      return;
    } else {
      fDiamond=(AliESDVertex*)meanVertex->GetObject();
    }
    //printf("\nRun %d, sigmaX %f, sigmaY %f\n",fTender->GetRun(),fDiamond->GetXRes(),fDiamond->GetYRes());
  }

  if (!fDiamond) return;

  // Redo the primary with the constraint ONLY if the updated mean vertex was found in the OCDB
  if ( (fDiamond->GetXRes())<2){ 
    AliVertexerTracks vertexer(event->GetMagneticField());
    vertexer.SetITSMode();
    vertexer.SetMinClusters(4);
    vertexer.SetVtxStart(fDiamond);
    AliESDVertex *pvertex = vertexer.FindPrimaryVertex(event);
    event->SetPrimaryVertexTracks(pvertex);
  }  
}
