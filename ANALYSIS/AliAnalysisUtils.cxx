#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliLog.h"
#include "AliAODVertex.h"

#include "AliAnalysisUtils.h"

ClassImp(AliAnalysisUtils)

//______________________________________________________________________
AliAnalysisUtils::AliAnalysisUtils():TObject(),
  fisAOD(kTRUE),
  fMinVtxContr(0),
  fMaxVtxZ(10.),
  fCutOnZVertexSPD(kTRUE)
{
  // Default contructor
}

//______________________________________________________________________
Bool_t AliAnalysisUtils::IsVertexSelected2013pA(AliVEvent *event)
{
  Bool_t accept = kFALSE;

  // Check whether the event is of AOD or ESD type
  AliAODEvent *aod = 0x0;
  AliESDEvent *esd = 0x0;
  aod = dynamic_cast<AliAODEvent*>(event);
  esd = dynamic_cast<AliESDEvent*>(event);

  if(aod) { 
    fisAOD = kTRUE; 
  } else {
    fisAOD = kFALSE;
    if(!esd) {
      AliFatal("Event is neither of AOD nor ESD type");
      return accept;
    }
  }

  const AliVVertex *trkVtx = fisAOD ? 
    dynamic_cast<const AliVVertex*>(aod->GetPrimaryVertex()) : 
    dynamic_cast<const AliVVertex*>(esd->GetPrimaryVertex()) ;
  if(!trkVtx || trkVtx->GetNContributors()<=fMinVtxContr){
    accept = kFALSE;
    return accept;
  }

  TString vtxTtl = trkVtx->GetTitle();
  if (!vtxTtl.Contains("VertexerTracks")) return accept;

  Float_t zvtx = trkVtx->GetZ();
  const AliVVertex* spdVtx = fisAOD ? 
    dynamic_cast<const AliVVertex*>(aod->GetPrimaryVertexSPD()) : 
    dynamic_cast<const AliVVertex*>(esd->GetPrimaryVertexSPD()) ;
  if (spdVtx->GetNContributors()<=fMinVtxContr) return accept;

  TString vtxTyp = spdVtx->GetTitle();
  Double_t cov[6]={0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return accept;
  if (fCutOnZVertexSPD && TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return accept;

  if (TMath::Abs(zvtx) > fMaxVtxZ) return accept;

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliAnalysisUtils::IsFirstEventInChunk(AliVEvent *event)
{

  Bool_t accept = kFALSE;

  // Check whether the event is of AOD or ESD type
  AliAODEvent *aod = 0x0;
  AliESDEvent *esd = 0x0;
  aod = dynamic_cast<AliAODEvent*>(event);
  esd = dynamic_cast<AliESDEvent*>(event);

  if(aod) { 
    fisAOD = kTRUE; 
  } else {
    fisAOD = kFALSE;
    if(!esd) {
      AliFatal("Event is neither of AOD nor ESD type");
      return accept;
    }
  }

  if(fisAOD){
    AliAODHeader *aodheader = 0x0;
    aodheader = aod->GetHeader();
    if(!aodheader){
      AliFatal("AOD header not there ?!");
      return kFALSE;
    }
    if(aodheader->GetEventNumberESDFile()==0) accept = kTRUE;
  } else {
    if(esd->GetEventNumberInFile()==0) accept = kTRUE;
  }

  return accept;
}
