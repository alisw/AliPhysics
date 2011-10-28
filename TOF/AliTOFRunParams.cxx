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
#include "TGraph.h"

ClassImp(AliTOFRunParams)

//_________________________________________________________

AliTOFRunParams::AliTOFRunParams() :
  TObject(),
  fNPoints(0),
  fTimestamp(NULL),
  fT0(NULL),
  fTOFResolution(NULL),
  fT0Spread(NULL),
  fNRuns(0),
  fRunNb(NULL),
  fRunFirstPoint(NULL),
  fRunLastPoint(NULL),
  fUseLHCClockPhase(kFALSE)
{
  /*
   * default constructor
   */
}

//_________________________________________________________

AliTOFRunParams::AliTOFRunParams(Int_t nPoints, Int_t nRuns) :
  TObject(),
  fNPoints(nPoints),
  fTimestamp(new UInt_t[nPoints]),
  fT0(new Float_t[nPoints]),
  fTOFResolution(new Float_t[nPoints]),
  fT0Spread(new Float_t[nPoints]),
  fNRuns(nRuns),
  fRunNb(new UInt_t[nRuns]),
  fRunFirstPoint(new UInt_t[nRuns]),
  fRunLastPoint(new UInt_t[nRuns]),
  fUseLHCClockPhase(kFALSE)
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
  if (fRunNb) delete [] fRunNb;
  if (fRunFirstPoint) delete [] fRunFirstPoint;
  if (fRunLastPoint) delete [] fRunLastPoint;
}

//_________________________________________________________

AliTOFRunParams::AliTOFRunParams(const AliTOFRunParams &source) :
  TObject(source),
  fNPoints(source.fNPoints),
  fTimestamp(new UInt_t[source.fNPoints]),
  fT0(new Float_t[source.fNPoints]),
  fTOFResolution(new Float_t[source.fNPoints]),
  fT0Spread(new Float_t[source.fNPoints]),
  fNRuns(source.fNRuns),
  fRunNb(new UInt_t[source.fNRuns]),
  fRunFirstPoint(new UInt_t[source.fNRuns]),
  fRunLastPoint(new UInt_t[source.fNRuns]),
  fUseLHCClockPhase(source.fUseLHCClockPhase)
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

  for (Int_t i = 0; i < fNRuns; i++) {
    fRunNb[i] = source.fRunNb[i];
    fRunFirstPoint[i] = source.fRunFirstPoint[i];
    fRunLastPoint[i] = source.fRunLastPoint[i];
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

  if (fNRuns != source.fNRuns) {
    if (fRunNb) delete [] fRunNb;
    if (fRunFirstPoint) delete [] fRunFirstPoint;
    if (fRunLastPoint) delete [] fRunLastPoint;
    fNRuns = source.fNRuns;
    fRunNb = new UInt_t[source.fNRuns];
    fRunFirstPoint = new UInt_t[source.fNRuns];
    fRunLastPoint = new UInt_t[source.fNRuns];
  }

  for (Int_t i = 0; i < fNRuns; i++) {
    fRunNb[i] = source.fRunNb[i];
    fRunFirstPoint[i] = source.fRunFirstPoint[i];
    fRunLastPoint[i] = source.fRunLastPoint[i];
  }

  fUseLHCClockPhase = source.fUseLHCClockPhase;

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

//_________________________________________________________

Float_t
AliTOFRunParams::Average(Float_t *data, Int_t first, Int_t last)
{
  /*
   * average
   */

  if (first < 0) first = 0;
  if (last >= fNPoints) last = fNPoints - 1;
  Float_t value = 0.;
  Int_t npt = 0;
  for (Int_t i = first; i <= last; i++) {
    value += data[i];
    npt++;
  }
  value /= npt;
  return value;

}

//_________________________________________________________

Float_t
AliTOFRunParams::Average(Float_t *data, UInt_t runNb)
{
  /*
   * average
   */

  /* critical cases:
     1. no measurement -> 0.
     2. no runNb structure -> average over all points
     3. runNb not found -> average over all points
  */
  if (fNPoints <= 0 || !fT0 || !fTimestamp) return 0.;
  if (fNRuns <= 0 || !fRunNb || !fRunFirstPoint || !fRunLastPoint) return Average(data, 0, fNPoints - 1);


  /* search for runNb */
  UInt_t runPoint = 0;
  Bool_t gotRunNb = kFALSE;
  for (Int_t irun = 0; irun < fNRuns; irun++) {
    if (fRunNb[irun] == runNb) {
      runPoint = irun;
      gotRunNb = kTRUE;
      break;
    }
  }
  if (!gotRunNb) return Average(data, 0, fNPoints - 1);

  /* average between first and last run points */
  UInt_t firstPoint = fRunFirstPoint[runPoint];
  UInt_t lastPoint = fRunLastPoint[runPoint];
  return Average(data, firstPoint, lastPoint);

}

//_________________________________________________________

TGraph *
AliTOFRunParams::DrawGraph(Float_t *data, Option_t* option)
{
  /*
   * draw
   */

  if (fNPoints == 0 || !data || !fTimestamp) return NULL;

  Float_t ts[1000000];
  for (Int_t i = 0; i < fNPoints; i++)
    ts[i] = fTimestamp[i];

  TGraph *graph = new TGraph(fNPoints, ts, data);
  graph->Draw(option);
  return graph;
}

//_________________________________________________________

TGraph *
AliTOFRunParams::DrawCorrelationGraph(Float_t *datax, Float_t *datay, Option_t* option)
{
  /*
   * draw
   */

  if (fNPoints == 0 || !datax || !datay) return NULL;

  TGraph *graph = new TGraph(fNPoints, datax, datay);
  graph->Draw(option);
  return graph;
}

