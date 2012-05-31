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
// * in OCDB in order to have TOF response correction
// * and actual resolution
// * 
// *
// *
// *

#include "AliTOFResponseParams.h"
#include "TGraph.h"

ClassImp(AliTOFResponseParams)

//_________________________________________________________

AliTOFResponseParams::AliTOFResponseParams() :
  TObject()
{
  /*
   * default constructor
   */

  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
    fNPoints[ipart] = 0;
}

//_________________________________________________________

AliTOFResponseParams::AliTOFResponseParams(Int_t *nPoints) :
  TObject()
{
  /*
   * default constructor
   */

  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++)
    fNPoints[ipart] = nPoints[ipart] < fgkMaxPoints ? nPoints[ipart] : fgkMaxPoints;
}

//_________________________________________________________

AliTOFResponseParams::~AliTOFResponseParams()
{
  /*
   * default destructor
   */

}

//_________________________________________________________

AliTOFResponseParams::AliTOFResponseParams(const AliTOFResponseParams &source) :
  TObject(source)
{
  /*
   * copy constructor
   */

  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    fNPoints[ipart] = source.fNPoints[ipart];
    for (Int_t ipoint = 0; ipoint < fNPoints[ipart]; ipoint++) {
      fP[ipart][ipoint] = source.fP[ipart][ipoint];
      fTExpCorr[ipart][ipoint] = source.fTExpCorr[ipart][ipoint];
    }
  }
    
}

//_________________________________________________________

AliTOFResponseParams &
AliTOFResponseParams::operator=(const AliTOFResponseParams &source)
{
  /*
   * operator=
   */

  if (this == &source) return *this;
  TObject::operator=(source);
  
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    fNPoints[ipart] = source.fNPoints[ipart];
    for (Int_t ipoint = 0; ipoint < fNPoints[ipart]; ipoint++) {
      fP[ipart][ipoint] = source.fP[ipart][ipoint];
      fTExpCorr[ipart][ipoint] = source.fTExpCorr[ipart][ipoint];
    }
  }
    
  return *this;
}

//_________________________________________________________

TGraph *
AliTOFResponseParams::DrawGraph(Int_t ipart, Option_t* option)
{
  /*
   * draw
   */

  if (ipart >= AliPID::kSPECIES) return NULL;
  if (fNPoints[ipart] == 0) return NULL;

  TGraph *graph = new TGraph(fNPoints[ipart], fP[ipart], fTExpCorr[ipart]);
  graph->Draw(option);
  return graph;
}

//_________________________________________________________

Double_t
AliTOFResponseParams::EvalTExpCorr(Int_t ipart, Double_t p)
{
  /*
   * eval corr
   */

  if (ipart >= AliPID::kSPECIES) return 0.;
  if (fNPoints[ipart] == 0) return 0.;
  if (p < fP[ipart][0]) return fTExpCorr[ipart][0];
  if (p >= fP[ipart][fNPoints[ipart] - 1]) return fTExpCorr[ipart][fNPoints[ipart] - 1];
  
  Int_t ipoint;
  for (ipoint = 0; ipoint < fNPoints[ipart] - 1; ipoint++)
    if (p >= fP[ipart][ipoint] && p < fP[ipart][ipoint + 1]) break;
  Double_t coeff = (fTExpCorr[ipart][ipoint + 1] - fTExpCorr[ipart][ipoint]) / (fP[ipart][ipoint + 1] - fP[ipart][ipoint]);
  Double_t corr = fTExpCorr[ipart][ipoint] + coeff * (p - fP[ipart][ipoint]);
  return corr;
}
