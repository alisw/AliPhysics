#include "AliAODMultEventClass.h"
#include <iostream>
#include <TString.h>
#include <TAxis.h>

//____________________________________________________________________
const Int_t*
AliAODMultEventClass::GetBins() 
{
  static const Int_t bins[] = { 1,  3,  6,  9,
				14, 19, 24, 29,
				39, 49, 59, 69,
				79, 89, 99, -1 };
  return bins;
}
//____________________________________________________________________
const TAxis*
AliAODMultEventClass::GetAxis() 
{
  static TAxis* ret = 0;
  if (!ret) {
    // Count the number of bins
    const Int_t* arr   = AliAODMultEventClass::GetBins();    
    const Int_t* tmp   = arr;
    Int_t        n     = 0;
    while ((*tmp >= 0)) { n++; tmp++; }

    // Now create temporary array 
    TArrayD bins(n+2);
    bins[0] = 0;
    Int_t i = 1;
    tmp     = arr;
    while ((*tmp >= 0)) {
      bins[i] = *tmp+(i==1 ? 0 : 1e-6);
      // Printf("bin[%2d]=%f", i, bins[i]);
      tmp++;
      i++;
    }
    bins[i] = bins[i-1]+10;
    // Printf("bin[%2d]=%f", i, bins[i]);

    ret = new TAxis(n+1, bins.GetArray());
    ret->SetTitle("Ref. Multiplicity (|#it{#eta}|<0.8)");
  }
  return ret;
}
//____________________________________________________________________
Int_t
AliAODMultEventClass::GetMultBin() const
{
  const Int_t* ptr   = GetBins();
  Int_t        bin   = 0;
  while (*ptr > 0) {
    if (fMult < *ptr) return bin;
    bin++;
    ptr++;
  }
  return bin;
}

//____________________________________________________________________
namespace {
  Float_t multCent(Int_t m)
  {
    const TAxis* a   = AliAODMultEventClass::GetAxis();
    Double_t     max = a->GetBinUpEdge(a->GetNbins());
    Float_t r = TMath::Min(max-1, Double_t(m));
    if (m >= max) Printf("mult=%4d -> %6.1f", m, r);
    return r;
  }
}
//____________________________________________________________________
Float_t
AliAODMultEventClass::GetCentrality(const TString& which) const
{
  TString::ECaseCompare k = TString::kIgnoreCase;
  if (which.EqualTo("MULT",     k)) return multCent(fMult);
  if (which.EqualTo("V0A",      k)) return fUtilV0A;
  if (which.EqualTo("V0M",      k)) return fUtilV0M;
  if (which.EqualTo("V0C",      k)) return fUtilV0C;
  if (which.EqualTo("V0AEQ",    k)) return fUtilV0AEq;
  if (which.EqualTo("V0MEQ",    k)) return fUtilV0MEq;
  if (which.EqualTo("V0CEQ",    k)) return fUtilV0CEq;
  if (which.EqualTo("MULTV0A",  k)) return fUtilV0A;
  if (which.EqualTo("MULTV0M",  k)) return fUtilV0M;
  if (which.EqualTo("MULTV0C",  k)) return fUtilV0C;
  if (which.EqualTo("MULTV0AEQ",k)) return fUtilV0AEq;
  if (which.EqualTo("MULTV0MEQ",k)) return fUtilV0MEq;
  if (which.EqualTo("MULTV0CEQ",k)) return fUtilV0CEq;
  if (which.EqualTo("CND",      k)) return fSelCND;
  if (which.EqualTo("SELCND",   k)) return fSelCND;
  if (which.EqualTo("SEL",      k)) return fSelCND;
  if (which.EqualTo("SELV0A",   k)) return fSelV0A;
  if (which.EqualTo("SELV0M",   k)) return fSelV0M;
  if (which.EqualTo("SELV0C",   k)) return fSelV0C;
  if (which.EqualTo("SELV0AEQ", k)) return fSelV0AEq;
  if (which.EqualTo("SELV0MEQ", k)) return fSelV0MEq;
  if (which.EqualTo("SELV0CEQ", k)) return fSelV0CEq;
  Warning("GetCentrality", "Unknown estimator: %s", which.Data());
  return -1;
}

//____________________________________________________________________
Float_t
AliAODMultEventClass::GetCentrality(UShort_t which, Bool_t util) const
{
  const Float_t* x[][2][2]
    = {{{ &fSelV0M,&fUtilV0M }, { &fSelV0MEq,&fUtilV0MEq }},
       {{ &fSelV0A,&fUtilV0A }, { &fSelV0AEq,&fUtilV0AEq }},
       {{ &fSelV0C,&fUtilV0C }, { &fSelV0CEq,&fUtilV0CEq }}};
  Bool_t isEq = (which & kEq);
  switch (which & 0xff) {
  case kV0M: return *(x[0][isEq][util]); break;
  case kV0A: return *(x[1][isEq][util]); break;
  case kV0C: return *(x[2][isEq][util]); break;
  case kCND: return fSelCND; break;
  default:
    Warning("GetCentrality", "Unknown estimator 0x%x for %s", which,
	    (util ? "util" : "sel"));
    break;
  }
  return -1;
}
//____________________________________________________________________
void
AliAODMultEventClass::SetCentrality(UShort_t which, Bool_t util, Float_t c)
{
  Float_t* x[][2][2] = {{{ &fSelV0M,&fUtilV0M }, { &fSelV0MEq,&fUtilV0MEq }},
			{{ &fSelV0A,&fUtilV0A }, { &fSelV0AEq,&fUtilV0AEq }},
			{{ &fSelV0C,&fUtilV0C }, { &fSelV0CEq,&fUtilV0CEq }}};
  Bool_t isEq = (which & kEq);
  switch (which & 0xff) {
  case kV0M: *(x[0][isEq][util]) = c; break;
  case kV0A: *(x[1][isEq][util]) = c; break;
  case kV0C: *(x[2][isEq][util]) = c; break;
  case kCND: fSelCND = c; break;
  default:
    Warning("SetCentrality", "Unknown estimator 0x%x for %s", which,
	    (util ? "util" : "sel"));
    break;
  }
}

//____________________________________________________________________
void
AliAODMultEventClass::Clear(Option_t*)
{
  fMult		= -1;
  fUtilV0M	= -1;
  fUtilV0A	= -1;
  fUtilV0C	= -1;
  fUtilV0MEq	= -1;
  fUtilV0AEq	= -1;
  fUtilV0CEq	= -1;
  fSelCND       = -1;
  fSelV0M	= -1;
  fSelV0A	= -1;
  fSelV0C	= -1;
  fSelV0MEq	= -1;
  fSelV0AEq	= -1;
  fSelV0CEq	= -1;
}

#define P(T,U,S) \
  Printf("%-12s:  %6.2f%% (util) %6.2f%% (sel)",T,(U),(S))
//____________________________________________________________________
void
AliAODMultEventClass::Print(Option_t*) const
{
  Printf("%-12s: %8d (bin %2d)", "Multiplicity", fMult, GetMultBin());
  P("V0M",   fUtilV0M,   fSelV0M);
  P("V0A",   fUtilV0A,   fSelV0A);
  P("V0C",   fUtilV0C,   fSelV0C);
  P("V0MEq", fUtilV0MEq, fSelV0MEq);
  P("V0AEq", fUtilV0AEq, fSelV0AEq);
  P("V0CEq", fUtilV0CEq, fSelV0CEq);
  Printf("%-12s:  %6.2f", "CND", fSelCND);
}
//
// EOF
// 
