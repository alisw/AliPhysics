#include "TCollection.h"
#include "TNamed.h"
#include <stdio.h>
#include "TH1.h"

// implement this function
// arguments are owned by the framework, don't delete.
// objects added to the outputs must ne NEW, so if you want to forward
// something from the input, forward a copy.
int process(TCollection* in, TCollection* out)
{
  TH1* h1 = NULL;
  TH1* h2 = NULL;
  TH1* h3 = NULL;
  TH1* h4 = NULL;
  TH1* h5 = NULL;
  TH1* h6 = NULL;
  TH1* h7 = NULL;
  TH1* h8 = NULL;
  TH1* h9 = NULL;

  TIter next(in);
  while (TObject* object = next())
  {
    printf("Received object: %s\n", object->GetName());
    TH1* tmp = dynamic_cast<TH1*>(object);
    if (!tmp) continue;
    TString name(tmp->GetName());
    if (name.Contains("fHistITSSPDvertexZ_CINT7ZAC-B-NOPF-CENT")) {
      h1 = tmp;
    } else
    if (name.Contains("fHistITSSPDvertexZ_CTVXZAC-B-NOPF-CENT")) {
      h2 = tmp;
    } else
    if (name.Contains("fHistITSSPDvertexZ_CINT7-B-NOPF-MUFAST")) {
      h4 = tmp;
    } else
    if (name.Contains("fHistITSSPDvertexZ_C0TVX-B-NOPF-MUFAST")) {
      h5 = tmp;
    } else
    if (name.Contains("fHistITSSPDvertexZ_C0V0M-B-NOPF-MUFAST")) {
      h7 = tmp;
    } else
    if (name.Contains("fHistITSSPDvertexZ_CTVXV0M-B-NOPF-MUFAST")) {
      h8 = tmp;
    }
  }

  if (h1 && h2) {
    h3 = new TH1F("RatioITSSPDvertexZ_h3","Ratio ITSSPDvertexZ CINT7ZAC-B-NOPF-CENT/CTVXZAC-B-NOPF-CENT",200,-100,100);
    h3->Divide(h1,h2);
    out->Add(h3);
  }
  if (h4 && h5) {
    h6 = new TH1F("RatioITSSPDvertexZ_h6","Ratio ITSSPDvertexZ C0TVX-B-NOPF-MUFAST/CINT7-B-NOPF-MUFAST",200,-100,100);
    h6->Divide(h5,h4);
    out->Add(h6);
  }
  if (h7 &&h8) {
    h9 = new TH1F("RatioITSSPDvertexZ_h9","Ratio ITSSPDvertexZ CTVXV0M-B-NOPF-MUFAST/C0V0M-B-NOPF-MUFAST",200,-100,100);
    h9->Divide(h8,h7);
    out->Add(h9);
  }

  return 0;
}
