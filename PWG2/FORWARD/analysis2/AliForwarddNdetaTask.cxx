//====================================================================
#include "AliForwarddNdetaTask.h"
#include <TMath.h>
#include <TH2D.h>
#include <TH1D.h>
#include <THStack.h>
#include <TList.h>
#include <AliAnalysisManager.h>
#include <AliAODEvent.h>
#include <AliAODHandler.h>
#include <AliAODInputHandler.h>
#include "AliForwardUtil.h"
#include "AliAODForwardMult.h"

//____________________________________________________________________
AliForwarddNdetaTask::AliForwarddNdetaTask()
  : AliBasedNdetaTask(),
    fSumPrimary(0),	//  Sum of primary histograms
    fSNNString(0),
    fSysString(0)
{}

//____________________________________________________________________
AliForwarddNdetaTask::AliForwarddNdetaTask(const char* /* name */)
  : AliBasedNdetaTask("Forward"), 
    fSumPrimary(0),	//  Sum of primary histograms
    fSNNString(0),
    fSysString(0)
{
}

//____________________________________________________________________
AliForwarddNdetaTask::AliForwarddNdetaTask(const AliForwarddNdetaTask& o)
  : AliBasedNdetaTask(o),
    fSumPrimary(o.fSumPrimary),	   // Sum of primary histograms
    fSNNString(o.fSNNString),	   //  
    fSysString(o.fSysString)	   //  
{}

//____________________________________________________________________
TH2D*
AliForwarddNdetaTask::GetHistogram(AliAODEvent* aod, Bool_t mc)
{
  TObject* obj = 0;
  if (mc) obj = aod->FindListObject("ForwardMC");
  else    obj = aod->FindListObject("Forward");

  // We should have a forward object at least 
  if (!obj) {
    if (!mc) AliWarning("No Forward object found AOD");
    return 0;
  }
  AliAODForwardMult* forward = static_cast<AliAODForwardMult*>(obj);

  // Here, we update the primary stuff 
  if (mc) { 
    TObject* oPrimary   = aod->FindListObject("primary");
    if (oPrimary) {
      TH2D* primary   = static_cast<TH2D*>(oPrimary);
      // Add contribtion from MC 
      
      if(fTriggerMask == AliAODForwardMult::kInel) {
	if (primary && !fSumPrimary) fSumPrimary = CloneHist(primary, "truth");
	else if (primary) fSumPrimary->Add(primary);
      }
      
      else if(fTriggerMask == AliAODForwardMult::kNSD && forward->IsTriggerBits(AliAODForwardMult::kMCNSD)) {
	if (primary && !fSumPrimary) fSumPrimary = CloneHist(primary, "truth"); 
	else if (primary) fSumPrimary->Add(primary);
      }
      
      
      
      /*   if(fTriggerMask == AliAODForwardMult::kNSD && forward->IsTriggerBits(AliAODForwardMult::kMCNSD)) {
	if (primary) fSumPrimary->Add(primary); 
      
      }
      
      else if (primary) fSumPrimary->Add(primary);*/
      
    }    
  }
  
  // Here, we get the update 
  if (!fSNNString) { 
    UShort_t sNN = forward->GetSNN();
    fSNNString = new TNamed("sNN", "");
    fSNNString->SetTitle(AliForwardUtil::CenterOfMassEnergyString(sNN));
    fSNNString->SetUniqueID(sNN);
    fSums->Add(fSNNString);
  
    UShort_t sys = forward->GetSystem();
    fSysString = new TNamed("sys", "");
    fSysString->SetTitle(AliForwardUtil::CollisionSystemString(sys));
    fSysString->SetUniqueID(sys);
    fSums->Add(fSysString);
  }

  return &(forward->GetHistogram());
}

//________________________________________________________________________
void 
AliForwarddNdetaTask::Terminate(Option_t *) 
{
  // Draw result to screen, or perform fitting, normalizations Called
  // once at the end of the query
  AliBasedNdetaTask::Terminate("");
  if(!fSums) {
    AliError("Could not retrieve TList fSums"); 
    return; 
  }
  if (!fOutput) { 
    AliError("Could not retrieve TList fOutput"); 
    return; 
  }

  fSumPrimary     = static_cast<TH2D*>(fSums->FindObject("truth"));
  fSNNString      = static_cast<TNamed*>(fSums->FindObject("sNN"));
  fSysString      = static_cast<TNamed*>(fSums->FindObject("sys"));

  fOutput->Add(fSNNString->Clone());
  fOutput->Add(fSysString->Clone());

  if (fSumPrimary) { 
    Int_t nAll        = Int_t(fTriggers->GetBinContent(kAll));
    Int_t nNSD        = Int_t(fTriggers->GetBinContent(kMCNSD));
    
    //  for(Int_t nn =1; nn<=fSumPrimary->GetNbinsX(); nn++) {
    //  for(Int_t mm =1; mm<=fSumPrimary->GetNbinsY(); mm++) {
    //	if(fSumPrimary->GetBinContent(nn,mm) > 0) std::cout<<fSumPrimary->GetBinContent(nn,mm)<<" +/-  "<<fSumPrimary->GetBinError(nn,mm)<<std::endl;
    //  }
    // }
 
    
    //TH1D* dndetaTruth = ProjectX(fSumPrimary,"dndetaTruth",
    //			 1,fSumPrimary->GetNbinsY(),false,true);
    TH1D* dndetaTruth = fSumPrimary->ProjectionX("dndetaTruth",1,fSumPrimary->GetNbinsY(),"e");

    
    if(fTriggerMask == AliAODForwardMult::kNSD)
      dndetaTruth->Scale(1./nNSD, "width");
    else
      dndetaTruth->Scale(1./nAll, "width");
    
    //for(Int_t nn =1; nn<=dndetaTruth->GetNbinsX(); nn++) {
    //  if(dndetaTruth->GetBinContent(nn) > 0) std::cout<<dndetaTruth->GetBinContent(nn)<<" +/-  "<<dndetaTruth->GetBinError(nn)<<std::endl;
    // }

    
    SetHistogramAttributes(dndetaTruth, kBlue+3, 22, "Monte-Carlo truth");
    
    fOutput->Add(dndetaTruth);
    fOutput->Add(Rebin(dndetaTruth));
  }
  // If only there was a away to get sqrt{s_NN} and beam type 
  

  PostData(2, fOutput);
}

//________________________________________________________________________
//
// EOF
//
