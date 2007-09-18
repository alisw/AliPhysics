#include "flow/AliFMDFlowTrue.h"
#include "flow/AliFMDFlowUtil.h"
#include <iostream>
#include <iomanip>
#include <TString.h>

//====================================================================
void
AliFMDFlowTrueBin::End()
{
  Double_t psi  = fPsi.Psi();
  Double_t dpsi = NormalizeAngle(fPsi.Order() * (psi-fPsiR));
  fResReal.Add(cos(dpsi));
}

//____________________________________________________________________
Double_t 
AliFMDFlowTrueBin::Value(CorType t) const
{ 
  Double_t e;
  return Value(e, t);
}
//____________________________________________________________________
Double_t
AliFMDFlowTrueBin::Value(Double_t& e2, CorType) const 
{
  return fHarmonic.Value(1, 0, e2);
}
//____________________________________________________________________
Double_t
AliFMDFlowTrueBin::Correction(Double_t& e2, CorType) const 
{
  e2 = fResReal.SqVar() / fResReal.N();
  return fResReal.Average();
}

//____________________________________________________________________
void
AliFMDFlowTrueBin::Print(Option_t*) const 
{
  Double_t e2v, e2r;
  Double_t v   = 100 * Value(e2v, AliFMDFlowBin::none);
  Double_t r   = 100 * Correction(e2r, AliFMDFlowBin::none);

  std::streamsize         old_prec  = std::cout.precision(3);
  std::ios_base::fmtflags old_flags = std::cout.setf(std::ios_base::fixed, 
						     std::ios_base::floatfield);
  std::cout << "  v" << std::setw(1) << fHarmonic.Order() << ":   True: "
	    << std::setw(6) << v << " +/- " 
	    << std::setw(6) << 100*sqrt(e2v) << " ["
	    << std::setw(7) << r << " +/- " 
	    << std::setw(7) << 100*sqrt(e2r) << "]";
  std::cout << std::endl;
  std::cout.precision(old_prec);
  std::cout.setf(old_flags, std::ios_base::floatfield);			       
}

//====================================================================
AliFMDFlowTrue1D::AliFMDFlowTrue1D(UShort_t order, const AliFMDFlowAxis& xaxis)
  : AliFMDFlowBinned1D(order, xaxis)
{
  // Delete old flow objects, and make new "true" ones. 
  for (UInt_t i = 0; i < xaxis.N(); i++) { 
    delete fBins[i];
    fBins[i] = new AliFMDFlowTrueBin(order);
  }
}

//____________________________________________________________________
void
AliFMDFlowTrue1D::SetPsi(Double_t psi) 
{ 
  for (UInt_t i = 0; i < fXAxis.N(); i++) 
    static_cast<AliFMDFlowTrueBin*>(fBins[i])->SetPsi(psi);
}

//____________________________________________________________________
void 
AliFMDFlowTrue1D::Print(Option_t* option) const
{
  TString opt(option);
  opt.ToLower();
  Bool_t det = opt.Contains("d");
  Bool_t sum = opt.Contains("s");
  if (det) AliFMDFlowBinned1D::Print("d");
  if (sum) { 
    std::cout << "    x |              Real \n" 
	      << "------+-------------------" << std::endl;

    std::streamsize         old_p = std::cout.precision(2);
    std::ios_base::fmtflags old_f = std::cout.setf(std::ios_base::fixed, 
						   std::ios_base::floatfield);
    for (UShort_t i = 0; i < fXAxis.N(); i++) { 
      Double_t x   = fXAxis.BinCenter(i);
      Double_t e2v;
      Double_t v   = fBins[i]->Value(e2v, AliFMDFlowBin::none);
      std::cout << std::setprecision(2) << std::setw(5) << x << " | " 
		<< std::setprecision(3) 
		<< std::setw(6) << 100 * v << " +/- " 
		<< std::setw(6) << 100 * sqrt(e2v) 
		<< std::endl; 
    }
    std::cout.precision(old_p);
    std::cout.setf(old_f, std::ios_base::floatfield);
  }
}


//____________________________________________________________________
//
// EOF
//
