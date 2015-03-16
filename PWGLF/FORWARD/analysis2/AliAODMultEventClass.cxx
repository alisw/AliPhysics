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
  fSelV0M	= -1;
  fSelV0A	= -1;
  fSelV0C	= -1;
  fSelV0MEq	= -1;
  fSelV0AEq	= -1;
  fSelV0CEq	= -1;
}

#define P(T,U,S) \
  Printf("%-11s:  %6.2f%% (util) %6.2f%% (sel)",T,(100*U),(100*S))
//____________________________________________________________________
void
AliAODMultEventClass::Print(Option_t*) const
{
  std::ostream& o = std::cout;
  o << "Multiplicity:  " << fMult << " (bin: " << GetMultBin() << ")\n";
  Printf("%-11s: %8d (bin %2d)", "Multiplicity", fMult, GetMultBin());
  P("V0M",   fUtilV0M,   fSelV0M);
  P("V0A",   fUtilV0A,   fSelV0A);
  P("V0C",   fUtilV0C,   fSelV0C);
  P("V0MEq", fUtilV0MEq, fSelV0MEq);
  P("V0AEq", fUtilV0AEq, fSelV0AEq);
  P("V0CEq", fUtilV0CEq, fSelV0CEq);
}
//
// EOF
// 
