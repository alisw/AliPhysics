/**
 * @file   AliForwardMultiplicityDistribution.cxx
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Thu Feb  7 01:03:52 2013
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_multdist
 * 
 */
#include <TH1D.h>
#include "AliForwardMultiplicityDistribution.h"

ClassImp(AliForwardMultiplicityDistribution)
#if 0
; // This is for Emacs - do not delete
#endif
//______________________________________________________________________
AliBaseMultTask::Bin*
AliForwardMultiplicityDistribution::MakeBin(Double_t l, Double_t h)
{
  return new Bin(l, h);
}
//_____________________________________________________________________
Bool_t 
AliForwardMultiplicityDistribution::CheckEvent(const AliAODForwardMult& fwd)
{
  return fIsSelected = AliBaseAODTask::CheckEvent(fwd);
}

//=====================================================================
void AliForwardMultiplicityDistribution::Bin::CreateOutputObjects(TList* cont,
								  Int_t max)
{
  //
  // Define eta bin output histograms
  //
  AliBaseMultTask::Bin::CreateOutputObjects(cont, max);
  TList* out = static_cast<TList*>(cont->FindObject(GetName()));
  if (!out) return;
  fHistPileUp      = new TH1D("multPileUp", GetTitle(), max, -0.5, max-.5);
  out->Add(fHistPileUp);
}
 
//_____________________________________________________________________
void AliForwardMultiplicityDistribution::
Bin::Process(TH1D*              dndetaForward, 
	     TH1D*              dndetaCentral,
	     TH1D*              normForward,   
	     TH1D*              normCentral, 
	     TH1D*              mc, 
	     Double_t           ipZ,
	     Bool_t             pileup,
	     Bool_t             selectedTrigger, 
	     Bool_t             isMCNSD, 
	     Bool_t             isESDNSD, 
	     const AliAODEvent& aodevent,
	     Double_t           minIPz,
	     Double_t           maxIPz) 
{
  //
  // Process a single eta bin
  //
  if (!selectedTrigger) return;
  if (ipZ < minIPz || ipZ > maxIPz) return;

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
  if (pileup) fHistPileUp->Fill(mult);
  else        fHist      ->Fill(mult);
  if (mc)     fHistMC    ->Fill(mcMult);
} 



//_____________________________________________________________________
//
//
// EOF
