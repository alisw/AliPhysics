/** @file 
    @brief Implementation of an Axis in a Flow "histogram" */
#include "flow/AliFMDFlowAxis.h"
#include <cmath>

//====================================================================
AliFMDFlowAxis::AliFMDFlowAxis(UShort_t n, Double_t* bins) 
  : fN(n), 
    fBins(0)
{
  fBins = new Double_t[fN+1];
  for (UInt_t i = 0; i <= fN; i++) fBins[i] = bins[i];
}
//____________________________________________________________________
AliFMDFlowAxis::AliFMDFlowAxis(UShort_t n, Double_t min, Double_t max)
  : fN(n), 
    fBins(0)
{
  fBins     = new Double_t[fN+1];
  Double_t dx = (max-min)/fN;
  for (UInt_t i = 0; i <= fN; i++) fBins[i] = min + i * dx;
}
//____________________________________________________________________
AliFMDFlowAxis::AliFMDFlowAxis(const AliFMDFlowAxis& other)
  : fN(other.fN), 
    fBins(0)
{
  fBins = new Double_t[fN+1];
  for (UInt_t i = 0; i <= fN; i++) fBins[i] = other.fBins[i];
}
//____________________________________________________________________
AliFMDFlowAxis& 
AliFMDFlowAxis::operator=(const AliFMDFlowAxis& other)
{
  if (fBins) delete [] fBins;
  fN    = other.fN;
  fBins = new Double_t[fN+1];
  for (UInt_t i = 0; i <= fN; i++) fBins[i] = other.fBins[i];
  return *this;
}
//____________________________________________________________________
AliFMDFlowAxis::~AliFMDFlowAxis()
{
  if (fBins) delete [] fBins;
  fBins = 0;
}
//____________________________________________________________________
Int_t 
AliFMDFlowAxis::FindBin(Double_t x) const 
{
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
  if (i >= fN) return 0;
  return (fBins[i+1] - fBins[i]);
}
//____________________________________________________________________
Double_t 
AliFMDFlowAxis::BinCenter(UShort_t i) const
{
  if (i >= fN) return 0;
  return fBins[i] + BinWidth(i) / 2;
}
//____________________________________________________________________
Double_t 
AliFMDFlowAxis::BinLower(UShort_t i) const
{
  if (i >= fN) return 0;
  return fBins[i];
}
//____________________________________________________________________
Double_t 
AliFMDFlowAxis::BinUpper(UShort_t i) const
{
  if (i >= fN) return 0;
  return fBins[i+1];
}

//____________________________________________________________________
//
// EOF
//
