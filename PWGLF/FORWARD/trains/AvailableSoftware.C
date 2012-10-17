/**
 * @file   AvailableSoftware.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 17:54:11 2012
 * 
 * @brief  Find available packages
 * 
 * @ingroup pwglf_forward_trains_util
 */

#ifndef AVAILABLESOFTWARE_C
#define AVAILABLESOFTWARE_C
#ifndef __CINT__
# include <TString.h>
# include <TSystem.h>
# include <TError.h>
# include <TObjArray.h>
#else
class TString;
#endif

/**
 * Helper code to find available packages on Grid
 * 
 *
 * @ingroup pwglf_forward_trains_util
 */
struct AvailableSoftware
{
  static Bool_t Check(TString& aliroot, TString& root)
  {
    // Figure out what to do.  
    // If mode == 0, then do nothing. 
    // If bit 0 is set in mode (0x1), then list and exit 
    // If bit 1 is set in mode (0x2), select last AliROOT/ROOT version 
    // If bit 2 is set in mode (0x4), select ROOT corresponding to AliROOT
    UShort_t mode = 0; 

    TString c("wget -q http://alimonitor.cern.ch/packages/ -O - | "
	      "sed -n -e '/<tr/,/<\\/tr>/ p' | "
	      "sed -n '/<a.*VO_ALICE@AliRoot::/,/VO_ALICE@ROOT::/ p' | "
	      "sed -n -e 's/.*VO_ALICE@AliRoot::\\([-0-9a-zA-Z]*\\).*/%\\1%/p' "
	      "  -e 's/.*VO_ALICE@ROOT::\\([-0-9a-zA-Z]*\\).*/\\1@/p' | "
	      "tr -d '\\n' | tr '@' '\\n' | tr '%' '\\t' ");

    if (aliroot.EqualTo("list", TString::kIgnoreCase) ||
	root.EqualTo("list", TString::kIgnoreCase) || 
	aliroot.IsNull()) {
      Warning("AvaliableSoftware::Check", "No AliROOT/ROOT version specified, "
	      "available packages are:\n" 
	      "\tAliROOT \tROOT:");
      gSystem->Exec(c);
      return false;
    }

    if (aliroot.EqualTo("last", TString::kIgnoreCase)) 
      mode |= 0x2;
    else if (!aliroot.IsNull()) 
      mode |= 0x4; 

    // Nothing to do 
    if (mode == 0) return true; 
    

    TString    values = gSystem->GetFromPipe(c);
    TObjArray* tokens = values.Tokenize(" \t\n");
    Int_t      n      = tokens->GetEntries();

    // If we asked to select the last possible version, do so here and get out
    if (mode & 0x2) { 
      aliroot = tokens->At(n-2)->GetName();
      root    = tokens->At(n-1)->GetName();
      Info("AvaliableSoftware::Check", 
	   "Selecting lastest possible AliROOT/ROOT: %s/%s", 
	   aliroot.Data(), root.Data());
      delete tokens;
      return true;
    }
    
    // We get here if we're asked to find a ROOT version compatible
    // with the selected AliROOT version. 
    for (Int_t i = 0; i < n; i += 2) {
      if (aliroot.EqualTo(tokens->At(i)->GetName(), 
				  TString::kIgnoreCase)) { 
	root = tokens->At(i+1)->GetName();
	Info("AvaliableSoftware::Check",
	     "Found ROOT version compatible with AliROOT %s: %s",
	     aliroot.Data(), root.Data());
	delete tokens;
	return true;
      }
    }
    // If we get here, then we didn't find a ROOT version compatible
    // with the selected AliROOT, and we should fail. 
    Warning("AvaliableSoftware::Check",
	    "Didn't find a ROOT version compatible with AliROOT %s", 
	    aliroot.Data());
    delete tokens; 
    return false;
  }
};
#endif
