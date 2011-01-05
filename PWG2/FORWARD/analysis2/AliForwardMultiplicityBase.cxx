//====================================================================
/**
 * @file 
 *
 * @ingroup pwg2_forward_tasks 
 */
#include "AliForwardMultiplicityBase.h"
#include "AliLog.h"
#include "AliAODHandler.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include <TROOT.h>
#include <iostream>
#include <iomanip>

//====================================================================
void
AliForwardMultiplicityBase::MarkEventForStore() const
{
  // Make sure the AOD tree is filled 
  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (!ah)  
    AliFatal("No AOD output handler set in analysis manager");

  ah->SetFillAOD(kTRUE);
}

//____________________________________________________________________
void
AliForwardMultiplicityBase::Print(Option_t* option) const
{
  std::cout << "AliForwardMultiplicityBase: " << GetName() << "\n" 
	    << "  Enable low flux code:   " << (fEnableLowFlux ? "yes" : "no")
	    << std::endl;
}

//
// EOF
//
