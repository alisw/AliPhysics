#if !defined(__CINT__) || defined(__MAKECINT__)
#include <EveBase/AliEveEventManager.h>

#include <TString.h>
#include <TROOT.h>
#endif

void visscan_raw_raw(const TString& path = ".", Bool_t showMuon = kTRUE, Bool_t showTrd = kFALSE)
{
  
  AliEveEventManager::SearchRawForCentralReconstruction();
  gROOT->ProcessLine(TString::Format(".x visscan_init.C(\"raw://\", \"%s\", %d, %d)",
				     path.Data(), showMuon, showTrd));
}
