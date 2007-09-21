/* Copyright (C) 2007 Christian Holm Christensen <cholm@nbi.dk>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */
//____________________________________________________________________
//
//  AliFMDFlowTrueBin: 
//    A specialised AliFMDFlowBin of flow in case the event plane is
//    well-known.  
//  AliFMDFlowTrue1D: 
//    A specialised AliFMDFlowBinned1D histogram in case the event
//    plane is well-known.   
//
#include "flow/AliFMDFlowTrue.h"
#include "flow/AliFMDFlowUtil.h"
#include <iostream>
#include <iomanip>
#include <TString.h>

//====================================================================
void
AliFMDFlowTrueBin::End()
{
  // Called at end of event 
  // PArameters: 
  //   none
  Double_t psi  = fPsi.Psi();
  Double_t dpsi = NormalizeAngle(fPsi.Order() * (psi-fPsiR));
  fResReal.Add(cos(dpsi));
}

//____________________________________________________________________
Double_t 
AliFMDFlowTrueBin::Value(CorType t) const
{ 
  // Get value of harmonic
  // PArameters: 
  //   see AliFMDFlowBin::Value
  Double_t e;
  return Value(e, t);
}
//____________________________________________________________________
Double_t
AliFMDFlowTrueBin::Value(Double_t& e2, CorType) const 
{
  // Get value of harmonic
  // PArameters: 
  //   see AliFMDFlowBin::Value
  return fHarmonic.Value(1, 0, e2);
}
//____________________________________________________________________
Double_t
AliFMDFlowTrueBin::Correction(Double_t& e2, CorType) const 
{
  // Get value of correction
  // PArameters: 
  //   see AliFMDFlowBin::Correction
  e2 = fResReal.SqVar() / fResReal.N();
  return fResReal.Average();
}

//____________________________________________________________________
void
AliFMDFlowTrueBin::Print(Option_t*) const 
{
  // Print to standard out
  // PArameters: 
  //   see AliFMDFlowBin::Print
  Double_t e2v, e2r;
  Double_t v   = 100 * Value(e2v, AliFMDFlowBin::kNone);
  Double_t r   = 100 * Correction(e2r, AliFMDFlowBin::kNone);

  std::streamsize         oldP  = std::cout.precision(3);
  std::ios_base::fmtflags oldF = std::cout.setf(std::ios_base::fixed, 
						     std::ios_base::floatfield);
  std::cout << "  v" << std::setw(1) << fHarmonic.Order() << ":   True: "
	    << std::setw(6) << v << " +/- " 
	    << std::setw(6) << 100*sqrt(e2v) << " ["
	    << std::setw(7) << r << " +/- " 
	    << std::setw(7) << 100*sqrt(e2r) << "]";
  std::cout << std::endl;
  std::cout.precision(oldP);
  std::cout.setf(oldF, std::ios_base::floatfield);			       
}

//====================================================================
AliFMDFlowTrue1D::AliFMDFlowTrue1D(UShort_t order, const AliFMDFlowAxis& xaxis)
  : AliFMDFlowBinned1D(order, xaxis)
{
  // Constructor. 
  // Parameters: 
  //   see AliFMDFlowBinned1D::AliFMDFlowBinned1D
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
  // Set event plane 
  // Parameters. 
  //    psi   The true, well-known, event plane angle 
  for (UInt_t i = 0; i < fXAxis.N(); i++) 
    static_cast<AliFMDFlowTrueBin*>(fBins[i])->SetPsi(psi);
}

//____________________________________________________________________
void 
AliFMDFlowTrue1D::Print(Option_t* option) const
{
  // Print to standard out. 
  // Parameters 
  //   See AliFMDFlowBinned1D::Print
  TString opt(option);
  opt.ToLower();
  Bool_t det = opt.Contains("d");
  Bool_t sum = opt.Contains("s");
  if (det) AliFMDFlowBinned1D::Print("d");
  if (sum) { 
    std::cout << "    x |              Real \n" 
	      << "------+-------------------" << std::endl;

    std::streamsize         oldP = std::cout.precision(2);
    std::ios_base::fmtflags oldF = std::cout.setf(std::ios_base::fixed, 
						  std::ios_base::floatfield);
    for (UShort_t i = 0; i < fXAxis.N(); i++) { 
      Double_t x   = fXAxis.BinCenter(i);
      Double_t e2v;
      Double_t v   = fBins[i]->Value(e2v, AliFMDFlowBin::kNone);
      std::cout << std::setprecision(2) << std::setw(5) << x << " | " 
		<< std::setprecision(3) 
		<< std::setw(6) << 100 * v << " +/- " 
		<< std::setw(6) << 100 * sqrt(e2v) 
		<< std::endl; 
    }
    std::cout.precision(oldP);
    std::cout.setf(oldF, std::ios_base::floatfield);
  }
}


//____________________________________________________________________
//
// EOF
//
