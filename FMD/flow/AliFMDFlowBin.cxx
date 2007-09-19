/** @file 
    @brief implementation of a Bin in a Flow histogram */
#include "flow/AliFMDFlowBin.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <TBrowser.h>

//====================================================================
void 
AliFMDFlowBin::Begin() 
{
  // Clear event plane calculators 
  fPsi.Clear();
  fPsiA.Clear();
  fPsiB.Clear();
}

//____________________________________________________________________
void 
AliFMDFlowBin::AddToEventPlane(Double_t phi, Double_t w, Bool_t a) 
{
  fPsi.Add(phi, w);
  if (a) fPsiA.Add(phi, w);
  else   fPsiB.Add(phi, w);
}

//____________________________________________________________________
void 
AliFMDFlowBin::AddToHarmonic(Double_t phi, Double_t w)
{
  // Disregard the obervation of phi from the event plane angle. 
  Double_t psi   = fPsi.Psi(phi, w);
  fHarmonic.Add(phi, psi);
}

//____________________________________________________________________
void 
AliFMDFlowBin::End()
{
  Double_t psi_A = fPsiA.Psi();
  Double_t psi_B = fPsiB.Psi();

  // Update the resolutions 
  fRes.Add(psi_A, psi_B);
  fResStar.Add(psi_A, psi_B);
  fResTdr.Add(psi_A, psi_B);
}

//____________________________________________________________________
void 
AliFMDFlowBin::Event(Double_t* phis, Double_t* ws, UInt_t n) 
{ 
  Begin();
  
  // Calculate split. 
  UInt_t split = n / 2;
  // First sub-event. 
  for (UInt_t i = 0; i < split; i++) 
    AddToEventPlane(phis[i], (ws ? ws[i] : 1), kTRUE);
  // Second sub-event. 
  for (UInt_t i = split; i < n; i++) 
    AddToEventPlane(phis[i], (ws ? ws[i] : 1), kFALSE);
  // Add contributions to the harmonic. 
  for (UInt_t i = 0; i < n; i++)     
    AddToHarmonic(phis[i], (ws ? ws[i] : 1));

  End();
}

//____________________________________________________________________
Double_t 
AliFMDFlowBin::Value(CorType t) const
{ 
  Double_t e;
  return Value(e, t);
}

//____________________________________________________________________
Double_t 
AliFMDFlowBin::EValue(CorType t) const 
{ 
  Double_t e2;
  Value(e2, t);
  return sqrt(e2);
}

//____________________________________________________________________
Double_t 
AliFMDFlowBin::Value(Double_t& e2, CorType t) const
{ 
  Double_t r, er2;
  r = Correction(er2, t);
  return fHarmonic.Value(r, er2, e2);
}

//____________________________________________________________________
Double_t 
AliFMDFlowBin::Correction(Double_t& er2, CorType t) const
{
  Double_t r = 1;
  UShort_t k = fHarmonic.Order()/fRes.Order();
  switch (t) { 
  case naive: r = fRes.Correction(k, er2);     break;
  case star:  r = fResStar.Correction(k, er2); break;
  case tdr:   r = fResTdr.Correction(k, er2);  break;
  default:    r = 1; er2 = 0;                  break;
  }
  return r;
}

//____________________________________________________________________
void 
AliFMDFlowBin::Finish() 
{}

//____________________________________________________________________
void
AliFMDFlowBin::Browse(TBrowser* b) 
{
  b->Add(&fPsi,      "Full event plane");
  b->Add(&fPsiA,     "Sub-event A event plane");
  b->Add(&fPsiB,     "Sub-event A event plane");
  b->Add(&fRes,      "Naive resolution");
  b->Add(&fResStar,  "STAR resolution");
  b->Add(&fResTdr,   "TDR resolution");
  b->Add(&fHarmonic, "Harmonic");
}

//____________________________________________________________________
void 
AliFMDFlowBin::Print(Option_t*) const
{
  Double_t e2v[4], v[4], r[4], e2r[4];
  const char* names[] = { "Bare", "Naive", "STAR", "TDR" };
  v[0] = 100 * Value(e2v[0], AliFMDFlowBin::none);
  v[1] = 100 * Value(e2v[1], AliFMDFlowBin::naive);
  v[2] = 100 * Value(e2v[2], AliFMDFlowBin::star);
  v[3] = 100 * Value(e2v[3], AliFMDFlowBin::tdr);
  r[0] = 100 * Correction(e2r[0], AliFMDFlowBin::none);
  r[1] = 100 * Correction(e2r[1], AliFMDFlowBin::naive);
  r[2] = 100 * Correction(e2r[2], AliFMDFlowBin::star);
  r[3] = 100 * Correction(e2r[3], AliFMDFlowBin::tdr);
  
  std::streamsize         old_prec  = std::cout.precision(3);
  std::ios_base::fmtflags old_flags = std::cout.setf(std::ios_base::fixed, 
						     std::ios_base::floatfield);
  std::cout << "  v" << std::setw(1) << fHarmonic.Order() << ": ";
  for (UInt_t i = 0; i < 4; i++) 
    std::cout << std::setw(6+(i == 0 ? 0 : 6)) << names[i] << ": " 
	      << std::setw(6) << v[i] << " +/- " 
	      << std::setw(6) << 100*sqrt(e2v[i]) << " ["
	      << std::setw(7) << r[i] << " +/- " 
	      << std::setw(7) << 100*sqrt(e2r[i]) << "]\n";
  std::cout << std::flush;
  std::cout.precision(old_prec);
  std::cout.setf(old_flags, std::ios_base::floatfield);
}


//____________________________________________________________________
//
// EOF
// 
