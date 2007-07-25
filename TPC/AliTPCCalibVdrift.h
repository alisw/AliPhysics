#ifndef AliTPCCalibVdrift_H
#define AliTPCCalibVdrift_H
/* Copyright(c) 2006-07, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//              Class AliTPCCalibVdrift
////////////////////////////////////////////////////////////////////////


class TObject;
class AliTPCSensorTempArray;
class TGraph;

class AliTPCCalibVdrift : public TNamed {

public:
  AliTPCCalibVdrift(AliTPCSensorTempArray *SensTemp, TObject *SensPres, TObject *SensGasComp);
  AliTPCCalibVdrift(const AliTPCCalibVdrift& source);
  virtual ~AliTPCCalibVdrift();
  AliTPCCalibVdrift& operator=(const AliTPCCalibVdrift& source);

  Double_t VdriftLinearHyperplaneApprox(Double_t dE, Double_t dT, Double_t dP, Double_t dCco2, Double_t dCn2);
  
  Double_t GetVdriftNominal();
  Double_t GetVdriftChange(Double_t x, Double_t y, Double_t z, UInt_t timeSec);

  Double_t GetMeanZVdriftChange(Double_t x, Double_t y, UInt_t timeSec);

  TGraph *MakeGraphMeanZVdriftChange(Double_t x, Double_t y, Int_t nPoints);

protected:

  AliTPCSensorTempArray *fSensTemp;   // Temperature sensors 
  TObject *fSensPres;         // Placeholder for Pressure sensors
  TObject *fSensGasComp;      // placeholder for GasConzentration infos  
  
  ClassDef(AliTPCCalibVdrift,1);

};
#endif
