/*
 * AliFemtoPPbpbLamBaseParticle.cxx
 *
 *  Created on: Aug 29, 2017
 *      Author: gu74req
 */

#include "AliLightNBasePart.h"
ClassImp(AliLightNBasePart)
AliLightNBasePart::AliLightNBasePart()
:fIsReset(false)
,fGTI(0)
,fTrackBufferSize(0)
,fP()
,fMCP()
,fMom()
,fTrackChi2perNDF(0)
,fPt(0)
,fMCPt(0)
,fP_TPC(0)
,fEta(0)
,fTheta(0)
,fMCTheta(0)
,fPhi(0)
,fMCPhi(0)
,fIDTracks(0)
,fCharge(0)
,fCPA(0)
,fOrigin(kUnknown)
,fPDGCode(0)
,fMCPDGCode(0)
,fPDGMotherWeak(0)
,fIsMC(false)
,fUse(true)
,fIsSet(true)
{
}

AliLightNBasePart::~AliLightNBasePart() {
}


