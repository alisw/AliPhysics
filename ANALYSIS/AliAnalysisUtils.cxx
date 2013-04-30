#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliVTrack.h"
#include "AliVEvent.h"
#include <TMatrixDSym.h>
#include <TMath.h>

#include "AliAnalysisUtils.h"

ClassImp(AliAnalysisUtils)

//______________________________________________________________________
AliAnalysisUtils::AliAnalysisUtils():TObject(),
  fisAOD(kTRUE),
  fMinVtxContr(0),
  fMaxVtxZ(10.),
  fCutOnZVertexSPD(kTRUE),
  fUseMVPlpSelection(kFALSE),
  fUseOutOfBunchPileUp(kFALSE),
  fMinPlpContribMV(5),
  fMaxPlpChi2MV(5.),
  fMinWDistMV(15.),
  fCheckPlpFromDifferentBCMV(kFALSE),
  fMinPlpContribSPD(5),
  fMinPlpZdistSPD(0.8),
  fnSigmaPlpZdistSPD(3.),
  fnSigmaPlpDiamXYSPD(2.),
  fnSigmaPlpDiamZSPD(5.),
  fUseSPDCutInMultBins(kFALSE)
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

//______________________________________________________________________
Bool_t AliAnalysisUtils::IsPileUpEvent(AliVEvent *event)
{
  Bool_t isPileUp=kFALSE;
  //check for multiple vertices
  if(fUseMVPlpSelection)isPileUp=IsPileUpMV(event);
  else isPileUp=IsPileUpSPD(event);
  //check for different BC 
  if(fUseOutOfBunchPileUp && IsOutOfBunchPileUp(event))isPileUp=kTRUE;
  
  return isPileUp;
}

//______________________________________________________________________
Bool_t AliAnalysisUtils::IsPileUpMV(AliVEvent *event)
{
  // check for multi-vertexer pile-up
  const AliAODEvent *aod = dynamic_cast<const AliAODEvent*>(event);
  const AliESDEvent *esd = dynamic_cast<const AliESDEvent*>(event);
  //
  if (!aod && !esd) {
    AliFatal("Event is neither of AOD nor ESD type");
    return kFALSE;
  }
  //
  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;
  Int_t nPlp = 0;
  //
  if (aod) {
    if ( !(nPlp=aod->GetNumberOfPileupVerticesTracks()) ) return kFALSE;
    vtPrm = aod->GetPrimaryVertex();
    if (vtPrm == aod->GetPrimaryVertexSPD()) return kTRUE; // there are pile-up vertices but no primary
  }
  else {
    if ( !(nPlp=esd->GetNumberOfPileupVerticesTracks())) return kFALSE;
    vtPrm = esd->GetPrimaryVertexTracks();
    if (((AliESDVertex*)vtPrm)->GetStatus()!=1) return kTRUE; // there are pile-up vertices but no primary
  }
  Int_t bcPrim = vtPrm->GetBC();
  //
  for (Int_t ipl=0;ipl<nPlp;ipl++) {
    vtPlp = aod ? (const AliVVertex*)aod->GetPileupVertexTracks(ipl) : (const AliVVertex*)esd->GetPileupVertexTracks(ipl);
    //
    if (vtPlp->GetNContributors() < fMinPlpContribMV) continue;
    if (vtPlp->GetChi2perNDF() > fMaxPlpChi2MV) continue;
    if(fCheckPlpFromDifferentBCMV)
      {
	Int_t bcPlp = vtPlp->GetBC();
	if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2) return kTRUE; // pile-up from other BC
      }
    //
    Double_t wDst = GetWDist(vtPrm,vtPlp);
    if (wDst<fMinWDistMV) continue;
    //
    return kTRUE; // pile-up: well separated vertices
  }
  //
  return kFALSE;
  //
}

//______________________________________________________________________
Bool_t AliAnalysisUtils::IsPileUpSPD(AliVEvent *event)
{
  // check for SPD pile-up
  const AliAODEvent *aod = dynamic_cast<const AliAODEvent*>(event);
  const AliESDEvent *esd = dynamic_cast<const AliESDEvent*>(event);
  //
  if (!aod && !esd) {
    AliFatal("Event is neither of AOD nor ESD type");
    return kFALSE;
  }
  //
  if (aod) return (fUseSPDCutInMultBins)?aod->IsPileupFromSPDInMultBins():aod->IsPileupFromSPD(fMinPlpContribSPD,fMinPlpZdistSPD,fnSigmaPlpZdistSPD,fnSigmaPlpDiamXYSPD,fnSigmaPlpDiamZSPD);
  else return (fUseSPDCutInMultBins)?esd->IsPileupFromSPDInMultBins():esd->IsPileupFromSPD(fMinPlpContribSPD,fMinPlpZdistSPD,fnSigmaPlpZdistSPD,fnSigmaPlpDiamXYSPD,fnSigmaPlpDiamZSPD);
}

//______________________________________________________________________
Bool_t AliAnalysisUtils::IsOutOfBunchPileUp(AliVEvent *event)
{
  // check for SPD pile-up
  const AliAODEvent *aod = dynamic_cast<const AliAODEvent*>(event);
  const AliESDEvent *esd = dynamic_cast<const AliESDEvent*>(event);
  //
  if (!aod && !esd) {
    AliFatal("Event is neither of AOD nor ESD type");
    return kFALSE;
  }
  Int_t bc2 = (aod)?aod->GetHeader()->GetIRInt2ClosestInteractionMap():esd->GetHeader()->GetIRInt2ClosestInteractionMap();
  if (bc2 != 0)
    return kTRUE;
  
  Int_t bc1 = (aod)?aod->GetHeader()->GetIRInt1ClosestInteractionMap():esd->GetHeader()->GetIRInt1ClosestInteractionMap();
  if (bc1 != 0)
    return kTRUE;
  
  return kFALSE;
}

//______________________________________________________________________
Double_t AliAnalysisUtils::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
{
  // calculate sqrt of weighted distance to other vertex
  if (!v0 || !v1) {
    printf("One of vertices is not valid\n");
    return 0;
  }
  static TMatrixDSym vVb(3);
  Double_t dist = -1;
  Double_t dx = v0->GetX()-v1->GetX();
  Double_t dy = v0->GetY()-v1->GetY();
  Double_t dz = v0->GetZ()-v1->GetZ();
  Double_t cov0[6],cov1[6];
  v0->GetCovarianceMatrix(cov0);
  v1->GetCovarianceMatrix(cov1);
  vVb(0,0) = cov0[0]+cov1[0];
  vVb(1,1) = cov0[2]+cov1[2];
  vVb(2,2) = cov0[5]+cov1[5];
  vVb(1,0) = vVb(0,1) = cov0[1]+cov1[1];
  vVb(0,2) = vVb(1,2) = vVb(2,0) = vVb(2,1) = 0.;
  vVb.InvertFast();
  if (!vVb.IsValid()) {printf("Singular Matrix\n"); return dist;}
  dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
    +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
  return dist>0 ? TMath::Sqrt(dist) : -1; 

}
