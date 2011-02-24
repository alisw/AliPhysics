/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////
//                Dielectron EventCuts                                  //
//                                                                       //
//                                                                       //
/*
Detailed description


*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <AliTriggerAnalysis.h>
#include <AliESDVertex.h>
#include <AliESDEvent.h>
#include <AliMultiplicity.h>

#include "AliDielectronEventCuts.h"

ClassImp(AliDielectronEventCuts)

AliDielectronEventCuts::AliDielectronEventCuts() :
  AliAnalysisCuts(),
  fVtxZmin(0.),
  fVtxZmax(0.),
  fRequireVtx(kFALSE),
  fMinVtxContributors(0),
  fMultITSTPC(kFALSE),
  fVtxType(kVtxTracks),
  fRequireV0and(0),
  fTriggerAnalysis(0x0),
  fkVertex(0x0)
{
  //
  // Default Constructor
  //
  
}

//______________________________________________
AliDielectronEventCuts::AliDielectronEventCuts(const char* name, const char* title) :
  AliAnalysisCuts(name, title),
  fVtxZmin(0.),
  fVtxZmax(0.),
  fRequireVtx(kFALSE),
  fMinVtxContributors(0),
  fMultITSTPC(kFALSE),
  fVtxType(kVtxTracks),
  fRequireV0and(0),
  fTriggerAnalysis(0x0),
  fkVertex(0x0)
{
  //
  // Named Constructor
  //
}

//______________________________________________
AliDielectronEventCuts::~AliDielectronEventCuts()
{
  //
  // Default Destructor
  //
  if (fTriggerAnalysis) delete fTriggerAnalysis;
}

//______________________________________________
Bool_t AliDielectronEventCuts::IsSelected(TObject* event)
{
  //
  // check the cuts
  //
  
  AliESDEvent *ev=dynamic_cast<AliESDEvent*>(event);
  if (!ev) return kFALSE;

  fkVertex=0x0;
  switch(fVtxType){
  case kVtxTracks: fkVertex=ev->GetPrimaryVertexTracks(); break;
  case kVtxSPD:    fkVertex=ev->GetPrimaryVertexSPD(); break;
  case kVtxTPC:    fkVertex=ev->GetPrimaryVertexTPC(); break;
  case kVtxAny:    fkVertex=ev->GetPrimaryVertex(); break;
  }

  if ((fRequireVtx||fVtxZmin<fVtxZmax||fMinVtxContributors>0)&&!fkVertex) return kFALSE;
  
  if (fVtxZmin<fVtxZmax){
    Double_t zvtx=fkVertex->GetZv();
    if (zvtx<fVtxZmin||zvtx>fVtxZmax) return kFALSE;
  }

  if (fMinVtxContributors>0){
    Int_t nCtrb = fkVertex->GetNContributors();
    if (nCtrb<fMinVtxContributors) return kFALSE;
  }

  if (fRequireV0and){
    if (!fTriggerAnalysis) fTriggerAnalysis=new AliTriggerAnalysis;
    Bool_t v0AND = kFALSE;
    if (fRequireV0and==1){
      Bool_t v0A       = fTriggerAnalysis->IsOfflineTriggerFired(ev, AliTriggerAnalysis::kV0A);
      Bool_t v0C       = fTriggerAnalysis->IsOfflineTriggerFired(ev, AliTriggerAnalysis::kV0C);
      v0AND = v0A && v0C;
    }

    if (fRequireV0and==2){
      Bool_t v0AHW     = (fTriggerAnalysis->V0Trigger(ev, AliTriggerAnalysis::kASide, kTRUE) == AliTriggerAnalysis::kV0BB);
      Bool_t v0CHW     = (fTriggerAnalysis->V0Trigger(ev, AliTriggerAnalysis::kCSide, kTRUE) == AliTriggerAnalysis::kV0BB);
      v0AND = v0AHW && v0CHW;
    }

    if (!v0AND) return kFALSE;
  }

  if (fMultITSTPC){
    const AliESDVertex *vtxESDTPC=ev->GetPrimaryVertexTPC();
    const AliMultiplicity *multESD = ev->GetMultiplicity();
    if ( vtxESDTPC && multESD && vtxESDTPC->GetNContributors() < (-10.+0.25*multESD->GetNumberOfITSClusters(0)) )
      return kFALSE;
  }
  
  return kTRUE;
}

