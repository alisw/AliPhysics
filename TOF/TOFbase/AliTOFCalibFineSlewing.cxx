/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// *
// *
// *
// * this class defines the TOF object to be stored
// * in OCDB for the fine time-slewing calibration
// * 
// * 
// *
// *
// *

#include "AliTOFCalibFineSlewing.h"
#include "TObjArray.h"
#include "TGraph.h"
#include "AliLog.h"

ClassImp(AliTOFCalibFineSlewing)

//_________________________________________________________

AliTOFCalibFineSlewing::AliTOFCalibFineSlewing() :
  TObject(),
  fSize(0),
  fX(NULL),
  fY(NULL)
{
  /*
   * default constructor
   */

  for (Int_t ich = 0; ich < 157248; ich++)
    fStart[ich] = 0;
}

//_________________________________________________________

AliTOFCalibFineSlewing::AliTOFCalibFineSlewing(TObjArray *oa) :
  TObject(),
  fSize(0),
  fX(NULL),
  fY(NULL)
{
  /*
   * standard constructor
   */
  
  // count number of valid points to be stored
  Int_t nvalid = 0;
  // loop over channels
  for (Int_t ich = 0; ich < 157248; ich++) {
    TGraph *g = (TGraph *)oa->At(ich);
    if (!g) continue;
    // loop over TGraph points
    for (Int_t ipt = 0; ipt < g->GetN(); ipt++) {
      if (!IsPointValid(g->GetX()[ipt], g->GetY()[ipt])) continue;
      nvalid++;
    }
  }
  AliInfo(Form("Found %d valid points", nvalid));
  
  // alloc arrays
  fSize = nvalid;
  fX = new UShort_t[fSize];
  fY = new Short_t[fSize];
  
  // store valid points
  Int_t nstored = 0;
  // loop over channels
  for (Int_t ich = 0; ich < 157248; ich++) {
    fStart[ich] = nstored;
    TGraph *g = (TGraph *)oa->At(ich);
    if (!g) continue;
    // loop over TGraph points
    for (Int_t ipt = 0; ipt < g->GetN(); ipt++) {
      if (!IsPointValid(g->GetX()[ipt], g->GetY()[ipt])) continue;
      fX[nstored] = (UShort_t)(1000. * g->GetX()[ipt]);
      fY[nstored] = (Short_t)(g->GetY()[ipt]);
      nstored++;
    }
  }
  AliInfo(Form("Stored %d points", nstored));
  
}

//_________________________________________________________

AliTOFCalibFineSlewing::AliTOFCalibFineSlewing(const AliTOFCalibFineSlewing &src) :
  TObject(src),
  fSize(src.fSize),
  fX(new UShort_t[src.fSize]),
  fY(new Short_t[src.fSize])
{
  /*
   * copy constructor
   */
  
  // copy data
  for (Int_t ich = 0; ich < 157248; ich++)
    fStart[ich] = src.fStart[ich];
  for (Int_t ipt = 0; ipt < fSize; ipt++) {
    fX[ipt] = src.fX[ipt];
    fY[ipt] = src.fY[ipt];
  }
  
}

//_________________________________________________________

AliTOFCalibFineSlewing &
AliTOFCalibFineSlewing::operator=(const AliTOFCalibFineSlewing &src)
{
  /*
   * operator=
   */
  
  if (this == &src) return *this;
  TObject::operator=(src);
  
  // delete and realloc if needed 
  if (fSize != src.fSize) {
    if (fX) delete [] fX;
    if (fY) delete [] fY;
    fSize = src.fSize;
    fX = new UShort_t[src.fSize];
    fY = new Short_t[src.fSize];
  }
  
  // copy data
  for (Int_t ich = 0; ich < 157248; ich++)
    fStart[ich] = src.fStart[ich];
  for (Int_t ipt = 0; ipt < fSize; ipt++) {
    fX[ipt] = src.fX[ipt];
    fY[ipt] = src.fY[ipt];
  }
  
}

//_________________________________________________________

AliTOFCalibFineSlewing::~AliTOFCalibFineSlewing()
{
  /*
   * default destructor
   */

  if (fX) delete [] fX;
  if (fY) delete [] fY;
}

//_________________________________________________________

Bool_t
AliTOFCalibFineSlewing::IsPointValid(Float_t x, Float_t y)
{
  /*
   * is point valid
   */

  if (TMath::Abs(y) >= kMaxShort) return kFALSE;
  if (x < 0 || x >= kMaxUShort) return kFALSE;
  return kTRUE;
}

//_________________________________________________________

Float_t
AliTOFCalibFineSlewing::Eval(Int_t ich, Float_t tot)
{
  /*
   * eval via interpolation
   */
  
  // get start point
  Int_t st = fStart[ich];
  // get number of points
  Int_t npt = 0;
  if (ich == 157247) npt = fSize - st;
  else               npt = fStart[ich + 1] - st;

  /*
   * the following is adapted starting from TGraph implementation 
   */

  if (npt == 0) return 0.;
  if (npt == 1) return (Float_t)fY[st];

  // convert input TOT to used x value
  Int_t x = (Int_t)(tot * 1000.);
  
  //linear interpolation
  //In case x is < fX[0] or > fX[fNpoints-1] return the extrapolated point
  
  //find points in graph around x assuming points are not sorted
  // (if point are sorted could use binary search)
  
  // find neighbours simply looping  all points
  // and find also the 2 adjacent points: (low2 < low < x < up < up2 )
  // needed in case x is outside the graph ascissa interval
  Int_t low  = -1;
  Int_t up  = -1;
  Int_t low2 = -1;
  Int_t up2 = -1;
  
  for (Int_t i = st; i < st + npt; ++i) {
    if (fX[i] < x) {
      if (low == -1 || fX[i] > fX[low])  {
	low2 = low;
	low = i;
      } else if (low2 == -1) low2 = i;
    } else if (fX[i] > x) {
      if (up  == -1 || fX[i] < fX[up])  {
	up2 = up;
	up = i;
      } else if (up2 == -1) up2 = i;
    } else // case x == fX[i]
      return (Float_t)fY[i]; // no interpolation needed
  }
  
  // treat cases when x is outside graph min max abscissa
  if (up == -1)  {
    up  = low;
    low = low2;
  }
  if (low == -1) {
    low = up;
    up  = up2;
  }
 
  //  assert(low != -1 && up != -1);
  if (!(low != -1 && up != -1)) return 0.;
  
  if (fX[low] == fX[up]) return (Float_t)fY[low];
  Float_t yup = (Float_t)fY[up];
  Float_t ylow = (Float_t)fY[low];
  Float_t a = (Float_t)(x - fX[up]);
  Float_t b = (Float_t)(fX[low] - fX[up]);
    //  Float_t yn = (Float_t)fY[up] + (x - fX[up]) * (fY[low] - fY[up]) / (fX[low] - fX[up]);
  Float_t yn = yup + a * (ylow - yup) / b;
  return yn;
  
}
