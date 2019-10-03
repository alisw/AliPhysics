/**
 * @file   AliForwardCreateResponseMatrices.cxx
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Thu Feb  7 01:01:24 2013
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_multdist
 */

#include <TH1D.h>
#include "AliForwardCreateResponseMatrices.h"
#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliAODEvent.h"


ClassImp(AliForwardCreateResponseMatrices)
#if 0
; // This is for Emacs - do not delete
#endif
//_____________________________________________________________________
AliBaseMultTask::Bin*
AliForwardCreateResponseMatrices::MakeBin(Double_t l, Double_t h)
{
  return new Bin(l,h);
}

//=====================================================================
void AliForwardCreateResponseMatrices::Bin::CreateOutputObjects(TList* cont,
								Int_t max)
{
  //
  // Define eta bin output histos
  //
  AliBaseMultTask::Bin::CreateOutputObjects(cont, max);
  TList* out = static_cast<TList*>(cont->FindObject(GetName()));
  if (!out) return;
  
  fResponseMatrix  = new TH2D("responseMatrix","Response Matrix;"
			      "MC_{truth};MC_{measured}",
			      max, -0.5, max-.5, max, -0.5, max-.5);
  out->Add(fResponseMatrix);
}
 

//_____________________________________________________________________
void 
AliForwardCreateResponseMatrices::
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
  // retreive MC particles from event
  Printf("Bin %s - selected: %s  IPz: %f (%f,%f)", GetName(),
	 selectedTrigger ? "yes" : " no", ipZ, minIPz, maxIPz);
  if (pileup) return;
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
  fHist->Fill(mult);
  fHistMC->Fill(mcMult);
  fResponseMatrix->Fill(mcMult, mult);
}




//_____________________________________________________________________
//
//
// EOF
