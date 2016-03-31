 /*
 * @file   AliForwardTriggerBiasCorrection.cxx
 * @author Valentina Zaccolo
 * @date   Mon Feb  3 11:30:24 2014
 ** 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_multdist
 */

#include <TH1D.h>
#include "AliForwardTriggerBiasCorrection.h"
#include "AliForwardMultiplicityDistribution.h"
#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliAODEvent.h"
#include "AliFMDMCEventInspector.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"


ClassImp(AliForwardTriggerBiasCorrection)
#if 0
; // This is for Emacs - do not delete
#endif
//_____________________________________________________________________
AliBaseMultTask::Bin*
AliForwardTriggerBiasCorrection::MakeBin(Double_t l, Double_t h)
{
  return new Bin(l,h);
}
//_____________________________________________________________________
Bool_t AliForwardTriggerBiasCorrection::IsESDClass(AliAODForwardMult* m) const
{
  return (m->IsTriggerBits(AliAODForwardMult::kInel) ||
	  m->IsTriggerBits(AliAODForwardMult::kV0AND));
}

//=====================================================================
void AliForwardTriggerBiasCorrection::Bin::CreateOutputObjects(TList* cont,
								Int_t max)
{
  //
  // Define eta bin output histos
  //
  AliBaseMultTask::Bin::CreateOutputObjects(cont, max);
  TList* out = static_cast<TList*>(cont->FindObject(GetName()));
  if (!out) return;

  fMCClass         = new TH1D("fMCClass","fMCClass", max,-0.5,max-.5);
  fESDClass        = new TH1D("fESDClass","fESDClass", max,-0.5,max-.5);
  fMCESDClass      = new TH1D("fMCESDClass","fMCESDClass", max,-0.5,max-.5);
  
  out->Add(fMCClass);
  out->Add(fESDClass);
  out->Add(fMCESDClass);
}
 

//_____________________________________________________________________
void 
AliForwardTriggerBiasCorrection::
Bin::Process(TH1D*              dndetaForward, 
	     TH1D*              dndetaCentral,
	     TH1D*              normForward,   
	     TH1D*              normCentral, 
	     TH1D*              mc, 
	     Double_t           ipZ,
	     Bool_t             pileup,
	     Bool_t             selectedTrigger, 
	     Bool_t             isMCClass, 
	     Bool_t             isESDClass, 
	     const AliAODEvent& aodevent,
	     Double_t           minIPz,
	     Double_t           maxIPz) 
{
  //
  // Process a single eta bin
  //  
  Double_t mcMult, mcErr, statErr, sysErr;
  Double_t mult = CalcMult(dndetaForward,
			   dndetaCentral,
			   normForward,
			   normCentral,
			   mc,
			   ipZ,
			   statErr,
			   sysErr,
			   mcMult,
			   mcErr);
  Double_t trMult  = mcMult;
  Double_t mcIPz   = mc->GetBinContent(0,0); // IP from MC stored here
  // The stuff below is redundant.  We've already filled the MC
  // histogram with the exact same information, and we have the IPz
  // from MC in the under-flow bin 0,0 of the MC histogram.
  // Furthermore, we use the information stored in the MC histogram to
  // form the response matrix, so we should also use the same
  // information to build the trigger bias correction.
  //
  // In fact, this class is really not very useful, since we could
  // have the same code in the MC class for doing the response matrix.
  // That would ensure that we have the same number in.
#if 0
  // MC particles from event
  TClonesArray* mcArray =
    static_cast<TClonesArray*>(aodevent.
			       FindListObject(AliAODMCParticle::
					      StdBranchName()));
  if(!mcArray){
    AliWarning("No MC array found in AOD. Try making it again.");
    return;
  }
  AliAODMCHeader* header = 
    static_cast<AliAODMCHeader*>(aodevent.
				 FindListObject(AliAODMCHeader::
						StdBranchName()));
  if (!header) {
    AliWarning("No header found.");
    return;
  }
  // Track loop: find MC truth multiplicity in selected eta bin - this
  // is probably redundant
  trMult           = 0;
  mcIPz            = header->GetVtxZ();
  Int_t    ntracks = mcArray->GetEntriesFast();
  for (Int_t it = 0; it < ntracks; it++) {
    AliAODMCParticle* particle = (AliAODMCParticle*)mcArray->At(it);
    if (!particle) {
      AliError(Form("Could not receive track %d", it));
      continue;
    }
    if (!particle->IsPhysicalPrimary()) continue;
    if (particle->Charge() == 0) continue;
    if (particle->Eta() > fEtaLow && particle->Eta() < fEtaHigh-0.0001)
      trMult++;
  }
#endif
    
  // fill fMCClass with multiplicity of MC truth NSD events with MC
  // truth |vtxz|<4
  if (mcIPz > minIPz && mcIPz < maxIPz) {
    if(isMCClass){
      fMCClass->Fill(trMult);
    }
  }
  // fill fESDClass with multiplicity from events with a reconstructed
  // NSD trigger and reconstructed |vtxz|<4
  if (ipZ > minIPz && ipZ < maxIPz){
    if(isESDClass){
      fESDClass->Fill(trMult);
    }
  }
  // fullfilling both requirements of fMCClass and fESDClass
  if(/* mcIPz > minIPz &&
	mcIPz < maxIPz && */
     ipZ               > minIPz &&
     ipZ               < maxIPz &&
     isMCClass && isESDClass){
    fMCESDClass->Fill(trMult);
  }

  if (!selectedTrigger) return;
  fHist->Fill(mult);
  fHistMC->Fill(mcMult);
}




//_____________________________________________________________________
//
//
// EOF
