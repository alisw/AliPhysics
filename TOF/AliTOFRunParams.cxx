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
// * in OCDB on a run-by-run basis in order to have the measurement
// * of the time evolution of T0 and of TOF resolution including
// * average T0 uncertainty
// *
// *
// *

#include "AliTOFRunParams.h"

ClassImp(AliTOFRunParams)

//_________________________________________________________

AliTOFRunParams::AliTOFRunParams() :
  TObject(),
  fNPoints(0),
  fTimestamp(NULL),
  fT0(NULL),
  fTOFResolution(NULL),
  fT0Spread(NULL)
{
  /*
   * default constructor
   */
}

//_________________________________________________________

AliTOFRunParams::AliTOFRunParams(Int_t nPoints) :
  TObject(),
  fNPoints(nPoints),
  fTimestamp(new UInt_t[nPoints]),
  fT0(new Float_t[nPoints]),
  fTOFResolution(new Float_t[nPoints]),
  fT0Spread(new Float_t[nPoints])
{
  /*
   * standard constructor
   */
}

//_________________________________________________________

AliTOFRunParams::~AliTOFRunParams()
{
  /*
   * default destructor
   */

  if (fTimestamp) delete [] fTimestamp;
  if (fT0) delete [] fT0;
  if (fTOFResolution) delete [] fTOFResolution;
  if (fT0Spread) delete [] fT0Spread;
}

//_________________________________________________________

AliTOFRunParams::AliTOFRunParams(const AliTOFRunParams &source) :
  TObject(source),
  fNPoints(source.fNPoints),
  fTimestamp(new UInt_t[source.fNPoints]),
  fT0(new Float_t[source.fNPoints]),
  fTOFResolution(new Float_t[source.fNPoints]),
  fT0Spread(new Float_t[source.fNPoints])
{
  /*
   * copy constructor
   */

  for (Int_t i = 0; i < fNPoints; i++) {
    fTimestamp[i] = source.fTimestamp[i];
    fT0[i] = source.fT0[i];
    fTOFResolution[i] = source.fTOFResolution[i];
    fT0Spread[i] = source.fT0Spread[i];
  }
  
}

//_________________________________________________________

AliTOFRunParams &
AliTOFRunParams::operator=(const AliTOFRunParams &source)
{
  /*
   * operator=
   */

  if (this == &source) return *this;
  TObject::operator=(source);
  
  if (fNPoints != source.fNPoints) {
    if (fTimestamp) delete [] fTimestamp;
    if (fT0) delete [] fT0;
    if (fTOFResolution) delete [] fTOFResolution;
    if (fT0Spread) delete [] fT0Spread;
    fNPoints = source.fNPoints;
    fTimestamp = new UInt_t[source.fNPoints];
    fT0 = new Float_t[source.fNPoints];
    fTOFResolution = new Float_t[source.fNPoints];
    fT0Spread = new Float_t[source.fNPoints];
  }

  for (Int_t i = 0; i < fNPoints; i++) {
    fTimestamp[i] = source.fTimestamp[i];
    fT0[i] = source.fT0[i];
    fTOFResolution[i] = source.fTOFResolution[i];
    fT0Spread[i] = source.fT0Spread[i];
  }

  return *this;
}

//_________________________________________________________

Float_t
AliTOFRunParams::EvalT0(UInt_t timestamp)
{
  /*
   * eval T0
   */

  /* critical cases:
     1. no measurement -> 0.
     2. single measurement -> single value
     3. timestamp before first measurement -> first value
     4. timestamp after last measurement -> last value
  */
  if (fNPoints <= 0 || !fT0 || !fTimestamp) return 0.;
  if (fNPoints == 1) return fT0[0];
  if (timestamp <= fTimestamp[0]) return fT0[0];
  if (timestamp >= fTimestamp[fNPoints - 1]) return fT0[fNPoints - 1];

  /* interpolate value */
  Int_t ipoint;
  for (ipoint = 0; ipoint < fNPoints - 1; ipoint++)
    if (timestamp >= fTimestamp[ipoint] && timestamp < fTimestamp[ipoint + 1])
      break;
  Float_t coeff = (fT0[ipoint + 1] - fT0[ipoint]) / (Float_t)(fTimestamp[ipoint + 1] - fTimestamp[ipoint]);
  Float_t t0 = fT0[ipoint] + coeff * (timestamp - fTimestamp[ipoint]);
  
  return t0;
}

//_________________________________________________________

Float_t
AliTOFRunParams::EvalTOFResolution(UInt_t timestamp)
{
  /*
   * eval TOF resolution
   */

  /* critical cases:
     1. no measurement -> 0.
     2. single measurement -> single value
     3. timestamp before first measurement -> first value
     4. timestamp after last measurement -> last value
  */
  if (fNPoints <= 0 || !fTOFResolution || !fTimestamp) return 0.;
  if (fNPoints == 1) return fTOFResolution[0];
  if (timestamp <= fTimestamp[0]) return fTOFResolution[0];
  if (timestamp >= fTimestamp[fNPoints - 1]) return fTOFResolution[fNPoints - 1];

  /* interpolate value */
  Int_t ipoint;
  for (ipoint = 0; ipoint < fNPoints - 1; ipoint++)
    if (timestamp >= fTimestamp[ipoint] && timestamp < fTimestamp[ipoint + 1])
      break;
  Float_t coeff = (fTOFResolution[ipoint + 1] - fTOFResolution[ipoint]) / (Float_t)(fTimestamp[ipoint + 1] - fTimestamp[ipoint]);
  Float_t reso = fTOFResolution[ipoint] + coeff * (timestamp - fTimestamp[ipoint]);
  
  return reso;
}

//_________________________________________________________

Float_t
AliTOFRunParams::EvalT0Spread(UInt_t timestamp)
{
  /*
   * eval T0 spread
   */

  /* critical cases:
     1. no measurement -> 0.
     2. single measurement -> single value
     3. timestamp before first measurement -> first value
     4. timestamp after last measurement -> last value
  */
  if (fNPoints <= 0 || !fT0Spread || !fTimestamp) return 0.;
  if (fNPoints == 1) return fT0Spread[0];
  if (timestamp <= fTimestamp[0]) return fT0Spread[0];
  if (timestamp >= fTimestamp[fNPoints - 1]) return fT0Spread[fNPoints - 1];

  /* interpolate value */
  Int_t ipoint;
  for (ipoint = 0; ipoint < fNPoints - 1; ipoint++)
    if (timestamp >= fTimestamp[ipoint] && timestamp < fTimestamp[ipoint + 1])
      break;
  Float_t coeff = (fT0Spread[ipoint + 1] - fT0Spread[ipoint]) / (Float_t)(fTimestamp[ipoint + 1] - fTimestamp[ipoint]);
  Float_t spread = fT0Spread[ipoint] + coeff * (timestamp - fTimestamp[ipoint]);
  
  return spread;
}

