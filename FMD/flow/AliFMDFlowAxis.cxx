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
/** @file 
    @brief Implementation of an Axis in a Flow "histogram" */
#include "flow/AliFMDFlowAxis.h"
#include <cmath>
#include <iostream>
#include <iomanip>
//  Axis object for the 
//  AliFMDFlowBinned1D and 
//  AliFMDFlowBinned2D 
//  "histograms" of  objects 
//  of class AliFMDFlowBin. 
//
//====================================================================
AliFMDFlowAxis::AliFMDFlowAxis(UShort_t n, Double_t* bins) 
  : fN(n), 
    fBins(0)
{
  // Constructor. 
  // Parameters 
  //   n    Number of bins 
  //   bins Bin boundaries (n+1 elements)
  fBins = new Double_t[fN+1];
  for (UInt_t i = 0; i <= fN; i++) fBins[i] = bins[i];
}
//____________________________________________________________________
AliFMDFlowAxis::AliFMDFlowAxis(UShort_t n, Double_t min, Double_t max)
  : fN(n), 
    fBins(0)
{
  // Constructor 
  //   n    Number of bins.
  //   min  Axis minimum
  //   max  Axis maximum
  fBins     = new Double_t[fN+1];
  Double_t dx = (max-min)/fN;
  for (UInt_t i = 0; i <= fN; i++) fBins[i] = min + i * dx;
}
//____________________________________________________________________
AliFMDFlowAxis::AliFMDFlowAxis(const AliFMDFlowAxis& other)
  : TObject(other), 
    fN(other.fN), 
    fBins(0)
{
  // Copy constructor 
  // Parameters 
  //   other   Object to copy from. 
  fBins = new Double_t[fN+1];
  for (UInt_t i = 0; i <= fN; i++) fBins[i] = other.fBins[i];
}
//____________________________________________________________________
AliFMDFlowAxis& 
AliFMDFlowAxis::operator=(const AliFMDFlowAxis& other)
{
  // Assignment operator
  // Parameters 
  //   other   Object to assign from. 
  if (fBins) delete [] fBins;
  fN    = other.fN;
  fBins = new Double_t[fN+1];
  for (UInt_t i = 0; i <= fN; i++) fBins[i] = other.fBins[i];
  return *this;
}
//____________________________________________________________________
AliFMDFlowAxis::~AliFMDFlowAxis()
{
  // destructor
  // Parameters 
  //   none
  if (fBins) delete [] fBins;
  fBins = 0;
}
//____________________________________________________________________
Int_t 
AliFMDFlowAxis::FindBin(Double_t x) const 
{
  // Find a bin corresponding to x 
  // Param 
  //    x   Value to find bin number for
  if (x < fBins[0])  return -1;
  if (x > fBins[fN]) return -1;
  UInt_t above, below, middle;
  above = fN+2;
  below = 0;
  while((above - below) > 1) {
    middle = (above + below) / 2;
    if (x == fBins[middle-1]) return middle-1;
    if (x  < fBins[middle-1]) above = middle;
    else                      below = middle;
  }
  return below-1;
}
//____________________________________________________________________
Double_t 
AliFMDFlowAxis::BinWidth(UShort_t i) const
{
  // Width of bin i 
  // Parameter 
  //   i   bin number to get width of 
  if (i >= fN) return 0;
  return (fBins[i+1] - fBins[i]);
}
//____________________________________________________________________
Double_t 
AliFMDFlowAxis::BinCenter(UShort_t i) const
{
  // Center of bin i 
  // Parameter 
  //   i   bin number to get center of 
  if (i >= fN) return 0;
  return fBins[i] + BinWidth(i) / 2;
}
//____________________________________________________________________
Double_t 
AliFMDFlowAxis::BinLower(UShort_t i) const
{
  // Center of bin i 
  // Parameter 
  //   i   bin number to get center of 
  if (i >= fN) return 0;
  return fBins[i];
}
//____________________________________________________________________
Double_t 
AliFMDFlowAxis::BinUpper(UShort_t i) const
{
  // Center of bin i 
  // Parameter 
  //   i   bin number to get center of 
  if (i >= fN) return 0;
  return fBins[i+1];
}

//____________________________________________________________________
void
AliFMDFlowAxis::Print(Option_t*) const
{
  // Print out
  // Parameter 
  //   none
  std::ios_base::fmtflags oldF = std::cout.setf(std::ios_base::fixed, 
						std::ios_base::floatfield);
  for (UShort_t i = 0; i < fN; i++) 
    std::cout << std::setw(5) << BinLower(i) << " - "
	      << std::setw(5) << BinUpper(i) << std::endl;
  std::cout.setf(oldF, std::ios_base::floatfield);
}


//____________________________________________________________________
//
// EOF
//
