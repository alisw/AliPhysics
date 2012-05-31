#include <TString.h>
#include <TROOT.h>

void visscan_mcresidual(const TString& path = ".", Bool_t showMuon = kTRUE, Bool_t showTrd = kFALSE)
{
  gROOT->ProcessLine(TString::Format(".x visscan_init.C(\"mcresidual://\", \"%s\", %d, %d)",
				     path.Data(), showMuon, showTrd));
}
