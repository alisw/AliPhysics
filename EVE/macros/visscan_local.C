#include <TString.h>
#include <TROOT.h>

// In order to open the HLT ESD Tree, instead of the Offline ESD Tree:
// run in the current directory:
// alieve mf_fix.C visscan_local.C'(".", kTRUE)'

void visscan_local(const TString& path = ".", Bool_t showHLTESDTree=kFALSE, Bool_t showMuon = kTRUE, Bool_t showTrd = kFALSE)
{
  gROOT->ProcessLine(TString::Format(".x visscan_init.C(\"local://$ALICE_ROOT/OCDB\", \"%s\", %d, %d)",
				     path.Data(), showHLTESDTree, showMuon, showTrd));
}
