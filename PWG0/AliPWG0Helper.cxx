/* $Id$ */

#include <AliPWG0Helper.h>

#include <TParticle.h>
#include <TParticlePDG.h>

#include <AliLog.h>
#include <AliESD.h>
#include <AliESDVertex.h>

//____________________________________________________________________
ClassImp(AliPWG0Helper)

//____________________________________________________________________
Bool_t AliPWG0Helper::IsEventTriggered(AliESD* aEsd)
{
  // check if the event was triggered
  //
  // MB should be
  // ITS_SPD_GFO_L0  : 32
  // VZERO_OR_LEFT   : 1
  // VZERO_OR_RIGHT  : 2

  ULong64_t triggerMask = aEsd->GetTriggerMask();

  if (triggerMask&32 && ((triggerMask&1) || (triggerMask&2)))
    return kTRUE;

  return kFALSE;
}

//____________________________________________________________________
Bool_t AliPWG0Helper::IsVertexReconstructed(AliESD* aEsd)
{
  // checks if the vertex is reasonable

  const AliESDVertex* vtxESD = aEsd->GetVertex();

  // the vertex should be reconstructed
  if (strcmp(vtxESD->GetName(), "default")==0)
    return kFALSE;

  Double_t vtx_res[3];
  vtx_res[0] = vtxESD->GetXRes();
  vtx_res[1] = vtxESD->GetYRes();
  vtx_res[2] = vtxESD->GetZRes();

  if (vtx_res[2]==0 || vtx_res[2]>0.1)
    return kFALSE;

  return kTRUE;
}

//____________________________________________________________________
Bool_t AliPWG0Helper::IsPrimaryCharged(TParticle* aParticle, Int_t aTotalPrimaries)
{
  //
  // Returns if the given particle is a primary particle
  // This function or a equivalent should be available in some common place of AliRoot
  //

  // if the particle has a daughter primary, we do not want to count it
  if (aParticle->GetFirstDaughter() != -1 && aParticle->GetFirstDaughter() < aTotalPrimaries)
  {
    //AliDebug(AliLog::kDebug+1, "Dropping particle because it has a daughter among the primaries.");
    return kFALSE;
  }

  Int_t pdgCode = TMath::Abs(aParticle->GetPdgCode());

  // skip quarks and gluon
  if (pdgCode <= 10 || pdgCode == 21)
  {
    //AliDebug(AliLog::kDebug+1, "Dropping particle because it is a quark or gluon.");
    return kFALSE;
  }

  if (strcmp(aParticle->GetName(),"XXX") == 0)
  {
    //AliDebug(AliLog::kDebug, Form("WARNING: There is a particle named XXX."));
    return kFALSE;
  }

  TParticlePDG* pdgPart = aParticle->GetPDG();

  if (strcmp(pdgPart->ParticleClass(),"Unknown") == 0)
  {
    //AliDebug(AliLog::kDebug, Form("WARNING: There is a particle with an unknown particle class (pdg code %d).", pdgCode));
    return kFALSE;
  }

  if (pdgPart->Charge() == 0)
  {
    //AliDebug(AliLog::kDebug+1, "Dropping particle because it is not charged.");
    return kFALSE;
  }

  return kTRUE;
}
