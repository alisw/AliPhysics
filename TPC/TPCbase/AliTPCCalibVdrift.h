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
class AliTPCTempMap;
class AliTPCCalibVdrift : public TNamed {

public:
  AliTPCCalibVdrift();
  AliTPCCalibVdrift(AliTPCSensorTempArray *SensTemp, AliDCSSensor *SensPres, TObject *SensGasComp);
  AliTPCCalibVdrift(const AliTPCCalibVdrift& source);
  virtual ~AliTPCCalibVdrift();
  AliTPCCalibVdrift& operator=(const AliTPCCalibVdrift& source);
  //
  // Interface for the reconstruction
  //
  Double_t GetPTRelative(UInt_t absTimeSec, Int_t side);

  //
  // Stefan interfaces - for v drift study
  //
  Double_t VdriftLinearHyperplaneApprox(Double_t dE, Double_t dT, Double_t dP, Double_t dCco2, Double_t dCn2);
  
  Double_t GetVdriftNominal();
  Double_t GetVdriftChange(Double_t x, Double_t y, Double_t z, UInt_t absTimeSec);

  Double_t GetMeanZVdriftChange(Double_t x, Double_t y, UInt_t absTimeSec);

  TGraph *MakeGraphMeanZVdriftChange(Double_t x, Double_t y, Int_t nPoints);
  Float_t GetNominalTemperature(){return fNominalTemp;}
  Float_t GetNominalPressure(){return fNominalPress;}

protected:
  //
  AliTPCSensorTempArray *fSensTemp;   // Temperature sensors 
  AliDCSSensor          *fSensPres;   // pressure sensor (cavernpress in GRP)
  AliTPCTempMap         *fTempMap;    // Temperature Map
  TObject *fSensGasComp;      // placeholder for GasConzentration infos  
  //
  // Nominal values
  //
  Float_t               fNominalTemp;    // nominal temperature in Kelvin
  Float_t               fNominalPress;    // nominal pressure    in mbar 
  ClassDef(AliTPCCalibVdrift,1);

};
#endif
