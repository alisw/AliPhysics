//
// This class contains the secondary correction and the double hit
// correction used in low-flux events.
//
//
#include "AliFMDCorrDoubleHit.h"
#include <TBrowser.h>
#include <TH1D.h>
#include <AliLog.h>
#include <iostream>

//____________________________________________________________________
AliFMDCorrDoubleHit::AliFMDCorrDoubleHit()
  : fCorrections()
{
  // Constructor 
  // 
  // 
  fCorrections.SetOwner(kTRUE);
  fCorrections.SetName("doubleHit");
}
//____________________________________________________________________
AliFMDCorrDoubleHit::AliFMDCorrDoubleHit(const AliFMDCorrDoubleHit& o)
  : TObject(o), 
    fCorrections(o.fCorrections)
{
  // Copy constructor 
  // 
  // Parameters: 
  //   o  Object to copy from 
}
//____________________________________________________________________
AliFMDCorrDoubleHit::~AliFMDCorrDoubleHit()
{
  // 
  // Destructor 
  // 
  fCorrections.Clear();
}
//____________________________________________________________________
AliFMDCorrDoubleHit&
AliFMDCorrDoubleHit::operator=(const AliFMDCorrDoubleHit& o)
{
  // Assignment operator 
  // 
  // Parameters: 
  //   o   Object to assign from 
  fCorrections   = o.fCorrections;

  return *this;
}
//____________________________________________________________________
TH1D*
AliFMDCorrDoubleHit::GetCorrection(UShort_t d, Char_t r) const
{
  // 
  // Get the double hit correction @f$ h_{r}(\eta)@f$ 
  //
  // Parameters:
  //    d Detector number 
  //    r Ring identifier 
  // 
  // Return:
  //    @f$ h_{r}(\eta)@f$ 
  //
  Int_t idx = GetRingIndex(d, r);
  if (idx < 0) return 0;

  TObject* o = fCorrections.At(idx);
  if (!o) { 
    AliWarning(Form("No double hit correction found for FMD%d%c", d, r));
    return 0;
  }
  return static_cast<TH1D*>(o);
}
//____________________________________________________________________
Int_t
AliFMDCorrDoubleHit::GetRingIndex(UShort_t d, Char_t r) const
{
  // 
  // Get the index corresponding to the given ring 
  //
  // Parameters:
  //    d Detector
  //    r Ring 
  // 
  // Return:
  //    Index (0 based) or negative in case of errors
  //

  switch (d) {
  case 1:  return 0;
  case 2:  return (r == 'I' || r == 'i' ? 1 : 2); break;  
  case 3:  return (r == 'I' || r == 'i' ? 3 : 4); break;  
  }
  AliWarning(Form("Index for FMD%d%c not found", d, r));
  return -1;
}
//____________________________________________________________________
Bool_t
AliFMDCorrDoubleHit::SetCorrection(UShort_t d, Char_t r, TH1D* h) 
{
  // 
  // Set the double hit correction @f$ h_{r}(\eta)@f$. Note, that the
  // object takes ownership of the passed pointer.
  //
  // Parameters:
  //    d  Detector number (1-3)
  //    r  Ring identifier (I or O)
  //    h  @f$ h_{r}(\eta)@f$ 
  // 
  // Return:
  //    true if operation succeeded 
  //
  Int_t idx = GetRingIndex(d, r);
  if (idx < 0) return kFALSE;

  h->SetName(Form("FMD%d%c", d, r));
  h->SetTitle(Form("Double hit correction for FMD%d%c", d,r));
  h->SetXTitle("#eta");
  h->SetYTitle("#sum_{i} N_{i,strips hit}(#eta)/"
	       "#sum_{i} N_{i,total hits}(#eta)");
  // h->SetFillColor(Color(d,r));
  h->SetFillStyle(3001);
  h->SetDirectory(0);
  h->SetStats(0);
  
  fCorrections.AddAtAndExpand(h, idx);
  return kTRUE;
}
//____________________________________________________________________
void
AliFMDCorrDoubleHit::Browse(TBrowser* b)
{
  // BRowse this object 
  // 
  // Parameters:
  //  b   Browser to use
  b->Add(&fCorrections);
}
//____________________________________________________________________
void
AliFMDCorrDoubleHit::Print(Option_t* option) const
{
  // Print this object 
  // 
  // Parameters:
  //  option   Passed to TH2D::Print 
  std::cout << "Double hit correction" << std::endl;
  fCorrections.Print(option);
}
    
//____________________________________________________________________
//
// EOF
//
