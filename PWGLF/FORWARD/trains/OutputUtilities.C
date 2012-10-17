/** 
 * @defgroup pwglf_forward_trains_util Utilities for Train setups
 *
 * @ingroup pwglf_forward_trains
 */
/**
 * @file   OutputUtilities.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 17:55:32 2012
 * 
 * @brief  Special output handling
 * 
 * @ingroup pwglf_forward_trains_util
 */
#ifndef TREEOUTPUTHELPER_C
#define TREEOUTPUTHELPER_C
#ifndef __CINT__
# include <TString.h>
# include <TError.h>
# include <AliAnalysisManager.h>
# include <AliAnalysisDataContainer.h>
# include <AliVEventHandler.h>
#else
class TString;
#endif

// ===================================================================
/**
 * Special output handling - data sets and remote storage
 *
 * @ingroup pwglf_forward_trains_util
 */
struct OutputUtilities
{
  /** 
   * Register output data set 
   * 
   * @param dsname Data set name 
   * 
   * @return true on success
   */
  static Bool_t RegisterDataset(const TString& dsname)
  {
    // Get the manager
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    
    // If we are asked to make a data-set, get the output handler and
    // common output container.
    AliVEventHandler*         handler = mgr->GetOutputEventHandler();
    if (!handler) return true;

    // Get the container 
    AliAnalysisDataContainer* cont    = mgr->GetCommonOutputContainer();
    if (!cont) { 
      Warning("OutputUtilities::RegisterDataset", 
	      "No common output container defined");
      return false;
    }

    // Make the name 
    TString nme(dsname);
    if (nme.IsNull()) nme = mgr->GetName();
    if (nme.IsNull()) { 
      Error("OutputUtilities::RegisterDataset", "No data set name specified");
      return false;
    }

    // Flag for data-set creation
    cont->SetRegisterDataset(true);

    handler->SetOutputFileName(nme);
    // cont->SetFileName(nme);

    TString base(handler->GetOutputFileName());
    base.ReplaceAll(".root","");
    Info("OutputUtilities::RegisterDataset", 
	 "Will register tree output AODs (%s%s) as dataset",
	 base.Data(), cont->GetTitle());

    return true;
  }
  /** 
   * Get the name of the registered data set
   * 
   * 
   * @return Name of the registered data set
   */
  static TString RegisteredDataset()
  {
    TString ret;

    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    AliVEventHandler*   oh  = mgr->GetOutputEventHandler();
    if (!oh) { 
      Warning("OutputUtilities::GetOutputDataSet", 
	      "No outout event handler defined");
      return ret;
    }
    AliAnalysisDataContainer* co  = mgr->GetCommonOutputContainer();
    if (!co) { 
      Warning("OutputUtilities::GetOutputDataSet", 
	      "No common output container defined");
      return ret;
    }
    if (!co->IsRegisterDataset()) { 
      Info("OutputUtilities::GetOutputDataSet", 
	   "Common output is not registered as dataset");
      return ret;
    }
    ret = oh->GetOutputFileName();
    // ret.ReplaceAll("TTree", "");
    ret.ReplaceAll(".root", "");
    // ret.Append(co->GetTitle());

    return ret;
  }
  /** 
   * Register special putput storage 
   * 
   * @param url Url (root://host/full_path)
   * 
   * @return true on success
   */
  static Bool_t RegisterStorage(const TString& url)
  {
    if (url.IsNull()) { 
      Error("OutputUtilities::RegisterStorage", "No storage URI specified");
      return false;
    }

    // Get the manager
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    
    // Get the container 
    AliAnalysisDataContainer* cont    = mgr->GetCommonOutputContainer();
    if (!cont) { 
      Warning("OutputUtilities::RegisterStorage", 
	      "No common output container defined");
      return false;
    }

    cont->SetSpecialOutput();
    mgr->SetSpecialOutputLocation(url);

    return true;
  }
};
#endif
//
// EOF
//
