#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliVVertex.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliVTrack.h"
#include "AliVEvent.h"
#include <TMatrixDSym.h>
#include <TMath.h>
#include "AliVMultiplicity.h"
#include "AliPPVsMultUtils.h"

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
  fUseSPDCutInMultBins(kFALSE),
  fASPDCvsTCut(65.),
  fBSPDCvsTCut(4.),
  fPPVsMultUtils(0x0)
{
  // Default contructor
}

AliAnalysisUtils::~AliAnalysisUtils()
{
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

  Double_t cov[6]={0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if (spdVtx->IsFromVertexerZ() && (zRes>0.25)) return accept;
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
    aodheader = dynamic_cast<AliAODHeader*>(aod->GetHeader());
    if(!aodheader) AliFatal("Not a standard AOD");
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
  Int_t bc2 = (aod)?((AliVAODHeader*)aod->GetHeader())->GetIRInt2ClosestInteractionMap():esd->GetHeader()->GetIRInt2ClosestInteractionMap();
  if (bc2 != 0)
    return kTRUE;
  
  Int_t bc1 = (aod)?((AliVAODHeader*)aod->GetHeader())->GetIRInt1ClosestInteractionMap():esd->GetHeader()->GetIRInt1ClosestInteractionMap();
  if (bc1 != 0)
    return kTRUE;
  
  return kFALSE;
}


//______________________________________________________________________
Bool_t AliAnalysisUtils::IsSPDClusterVsTrackletBG(AliVEvent *event){
  Int_t nClustersLayer0 = event->GetNumberOfITSClusters(0);
  Int_t nClustersLayer1 = event->GetNumberOfITSClusters(1);
  Int_t nTracklets      = event->GetMultiplicity()->GetNumberOfTracklets();
  if (nClustersLayer0 + nClustersLayer1 > fASPDCvsTCut + nTracklets*fBSPDCvsTCut) return kTRUE;
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

//______________________________________________________________________
Float_t AliAnalysisUtils::GetMultiplicityPercentile(AliVEvent *event, TString lMethod, Bool_t lEmbedEventSelection){
  if(!fPPVsMultUtils)
    fPPVsMultUtils=new AliPPVsMultUtils();
  if( (event->InheritsFrom("AliAODEvent")) || (event->InheritsFrom("AliESDEvent")) ) return fPPVsMultUtils->GetMultiplicityPercentile(event,lMethod,lEmbedEventSelection);
  else {
    AliFatal("Event is neither of AOD nor ESD type"); 
    return -999.;
  }
}
//______________________________________________________________________
Bool_t AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(Int_t index, AliMCEvent* mcEv){
  // interface method for analyses on ESDs
  // returns kTRUE if a particle is produced in a pileup event in case of MC with pileup

  // check if the event header is that of a cocktail (AliGenPileup creates a cocktail):
  AliGenCocktailEventHeader *cocktailHeader = dynamic_cast<AliGenCocktailEventHeader *>(mcEv->GenEventHeader());
  if (cocktailHeader == nullptr) return kFALSE;
  
  Int_t totPrimaries=mcEv->GetNumberOfPrimaries();
  if(index>=totPrimaries){
    // particles from the transport, get mother
    while(index>=totPrimaries) index=mcEv->GetLabelOfParticleMother(index);
  }
  TList *lgen = cocktailHeader->GetHeaders();
  if(!lgen) return kFALSE;
  return IsParticleFromOutOfBunchPileupCollision(index,lgen);
}

//______________________________________________________________________
Bool_t AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(Int_t index, AliAODMCHeader* aodMCHeader, TClonesArray *arrayMC){
  // interface method for analyses on AODs
  // returns kTRUE if a particle is produced in a pileup event in case of MC with pileup

  TList *lgen = aodMCHeader->GetCocktailHeaders();
  if(!lgen) return kFALSE;
  Int_t nh=lgen->GetEntries();
  Int_t totPrimaries=0;
  for(Int_t i=0;i<nh;i++){
    AliGenEventHeader* gh=(AliGenEventHeader*)lgen->At(i);
    totPrimaries+=gh->NProduced();
  }
  if(index>=totPrimaries){
    while(index>=totPrimaries){
      // particles from the transport, get mother
      AliAODMCParticle* part =(AliAODMCParticle*)arrayMC->At(index);
      index=part->GetMother();
    }
  }
  return IsParticleFromOutOfBunchPileupCollision(index,lgen);
}
//______________________________________________________________________
Bool_t AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(Int_t index, TList *lgen){
  
  // returns kTRUE if a particle is produced in a pileup event in case of MC with pileup

  // retrieve the header of the generator that produced the considered particle
  Int_t nh=lgen->GetEntries();
  AliGenEventHeader* theGener=0x0;
  Int_t nsumpart=0;
  for(Int_t i=0;i<nh;i++){
    AliGenEventHeader* gh=(AliGenEventHeader*)lgen->At(i);
    if(gh->InheritsFrom(AliGenCocktailEventHeader::Class())){
      AliGenCocktailEventHeader* gc=dynamic_cast<AliGenCocktailEventHeader*>(gh);
      TList* lh2=gc->GetHeaders();
      if(lh2){
	Int_t nh2=lh2->GetEntries();
	for(Int_t i2=0;i2<nh2;i2++){
	  AliGenEventHeader* gh2=(AliGenEventHeader*)lh2->At(i2);
	  Int_t npart=gh2->NProduced();
	  if(index>=nsumpart && index<(nsumpart+npart)){
	    theGener=gh2;
	    break;
	  }
	  nsumpart+=npart;
	}
      }
    }else{
      Int_t npart=gh->NProduced();
      if(index>=nsumpart && index<(nsumpart+npart)){
	theGener=gh;
	break;
      }
      nsumpart+=npart;
    }
    if(theGener) break;
  }
  if(!theGener) return kFALSE;
  
  Double_t timeNs = theGener->InteractionTime() * 1e9;
  if (TMath::Abs(timeNs) > 3.0){
    // Out of bunch pileup according to collision time
    // use 3 ns window around the trigger time
    return kTRUE;
  }
  return kFALSE;
}
//______________________________________________________________________
Bool_t AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(AliMCEvent* mcEv){
  // Interface method for ESDs
  // returns kTRUE if there is >1 collision in the bunch crossing of the trigger
  // use 3 ns window around the trigger time

  AliGenCocktailEventHeader *cocktailHeader = dynamic_cast<AliGenCocktailEventHeader *>(mcEv->GenEventHeader());
  if (cocktailHeader == nullptr) return kFALSE;
  TList *lgen = cocktailHeader->GetHeaders();
  return IsSameBunchPileupInGeneratedEvent(lgen);
}
//______________________________________________________________________
Bool_t AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(AliAODMCHeader* aodMCHeader){
  // Interface method for AODs
  // returns kTRUE if there is >1 collision in the bunch crossing of the trigger
  // use 3 ns window around the trigger time
  
  TList *lgen = aodMCHeader->GetCocktailHeaders();
  return IsSameBunchPileupInGeneratedEvent(lgen);
}
//______________________________________________________________________
Bool_t AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(TList *lgen){
  // returns kTRUE if there is >1 collision in the bunch crossing of the trigger
  // use 3 ns window around the trigger time

  if(!lgen) return kFALSE;
  Int_t nh=lgen->GetEntries();
  Int_t nCollis=0;
  for(Int_t i=0;i<nh;i++){
    AliGenEventHeader* gh=(AliGenEventHeader*)lgen->At(i);
    if(gh->InheritsFrom(AliGenCocktailEventHeader::Class())){
      AliGenCocktailEventHeader* gc=dynamic_cast<AliGenCocktailEventHeader*>(gh);
      TList* lh2=gc->GetHeaders();
      if(lh2){
	Int_t nh2=lh2->GetEntries();
	for(Int_t i2=0;i2<nh2;i2++){
	  AliGenEventHeader* gh2=(AliGenEventHeader*)lh2->At(i2);
	  Double_t timeNs = gh2->InteractionTime() * 1e9;
	  if (TMath::Abs(timeNs) < 3.0) nCollis++;
	}
      }
    }else{
      Double_t timeNs = gh->InteractionTime() * 1e9;
      if (TMath::Abs(timeNs) < 3.0) nCollis++;
    }
  }
  if(nCollis>1) return kTRUE;
  return kFALSE;
}
