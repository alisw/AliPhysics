#include <TString.h>
#include <TROOT.h>

void visscan_raw(const TString& path = ".", Bool_t showMuon = kTRUE, Bool_t showTrd = kFALSE)
{
  gROOT->ProcessLine(TString::Format(".x visscan_init.C(\"raw://\", \"%s\", %d, %d)",
				     path.Data(), showMuon, showTrd));
}
