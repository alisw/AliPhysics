#ifndef ALITOFRUNPARAMS_H
#define ALITOFRUNPARAMS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

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

#include "TObject.h"

class TGraph;

class AliTOFRunParams :
public TObject
{

 public:

  AliTOFRunParams(); // default constructor
  AliTOFRunParams(Int_t nPoints, Int_t nRuns = 1); // standard constructor
  virtual ~AliTOFRunParams(); // default destructor
  AliTOFRunParams(const AliTOFRunParams &source); // copy constructor
  AliTOFRunParams &operator=(const AliTOFRunParams &source); // operator=

  Int_t GetNPoints() const {return fNPoints;}; // getter
  UInt_t GetTimestamp(Int_t i) const {return fTimestamp && i < fNPoints ? fTimestamp[i] : 0;}; // getter
  Float_t GetT0(Int_t i) const {return fT0 && i < fNPoints ? fT0[i] : 0.;}; // getter
  Float_t GetTOFResolution(Int_t i) const {return fTOFResolution && i < fNPoints ? fTOFResolution[i] : 0.;}; // getter
  Float_t GetT0Spread(Int_t i) const {return fT0Spread && i < fNPoints ? fT0Spread[i] : 0.;}; // getter
  
  Int_t GetNRuns() const {return fNRuns;}; // getter
  UInt_t GetRunNb(Int_t i) const {return fRunNb && i < fNRuns ? fRunNb[i] : 0;}; // getter
  UInt_t GetRunFirstPoint(Int_t i) const {return fRunFirstPoint && i < fNRuns ? fRunFirstPoint[i] : 0;}; // getter
  UInt_t GetRunLastPoint(Int_t i) const {return fRunLastPoint && i < fNRuns ? fRunLastPoint[i] : 0;}; // getter

  void SetTimestamp(UInt_t *value) {if (fTimestamp) for (Int_t i = 0; i < fNPoints; i++) fTimestamp[i] = value[i];}; // setter
  void SetT0(Float_t *value) {if (fT0) for (Int_t i = 0; i < fNPoints; i++) fT0[i] = value[i];}; // setter
  void SetTOFResolution(Float_t *value) {if (fTOFResolution) for (Int_t i = 0; i < fNPoints; i++) fTOFResolution[i] = value[i];}; // setter
  void SetT0Spread(Float_t *value) {if (fT0Spread) for (Int_t i = 0; i < fNPoints; i++) fT0Spread[i] = value[i];}; // setter

  void SetRunNb(UInt_t *value) {if (fRunNb) for (Int_t i = 0; i < fNRuns; i++) fRunNb[i] = value[i];}; // setter
  void SetRunFirstPoint(UInt_t *value) {if (fRunFirstPoint) for (Int_t i = 0; i < fNRuns; i++) fRunFirstPoint[i] = value[i];}; // setter
  void SetRunLastPoint(UInt_t *value) {if (fRunLastPoint) for (Int_t i = 0; i < fNRuns; i++) fRunLastPoint[i] = value[i];}; // setter

  Float_t EvalT0(UInt_t timestamp); // eval T0
  Float_t EvalTOFResolution(UInt_t timestamp); // eval TOF resolution
  Float_t EvalT0Spread(UInt_t timestamp); // eval T0 spread
  
  Float_t AverageT0(UInt_t runNb) {return Average(fT0, runNb);}; // average T0
  Float_t AverageTOFResolution(UInt_t runNb) {return Average(fTOFResolution, runNb);}; // average TOF resolution
  Float_t AverageT0Spread(UInt_t runNb) {return Average(fT0Spread, runNb);}; // average T0 spread

  TGraph *DrawGraphT0(Option_t *option = "") {return DrawGraph(fT0, option);}; // draw graph t0
  TGraph *DrawGraphTOFResolution(Option_t *option = "") {return DrawGraph(fTOFResolution, option);}; // draw graph t0
  TGraph *DrawGraphT0Spread(Option_t *option = "") {return DrawGraph(fT0Spread, option);}; // draw graph t0
  TGraph *DrawCorrelationGraphTOFResolutionT0Spread(Option_t *option = "") {return DrawCorrelationGraph(fT0Spread, fTOFResolution, option);}; // draw correlation graph


 private:
  
  /* private methods */
  
  Float_t Average(Float_t *data, UInt_t runNb); // average
  Float_t Average(Float_t *data, Int_t first, Int_t last); // average
  TGraph *DrawGraph(Float_t *data, Option_t *option = ""); // draw graph
  TGraph *DrawCorrelationGraph(Float_t *datax, Float_t *datay, Option_t *option = ""); // draw graph
  
  /* data members */

  Int_t fNPoints;
  UInt_t *fTimestamp; //[fNPoints] time stamp
  Float_t *fT0; //[fNPoints] average T0 (ps)
  Float_t *fTOFResolution; //[fNPoints] average TOF resolution (T0 uncertainty included) (ps)
  Float_t *fT0Spread; //[fNPoints] estimated T0 spread (from vertex spread z) (ps)

  Int_t fNRuns;
  UInt_t *fRunNb; //[fNRuns] run number
  UInt_t *fRunFirstPoint; //[fNRuns] run start point
  UInt_t *fRunLastPoint; //[fNRuns] run last point


  ClassDef(AliTOFRunParams, 2);
};

#endif /* ALITOFRUNPARAMS_H */
