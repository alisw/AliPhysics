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
  AliTPCCalibVdrift(AliTPCSensorTempArray *SensTemp, AliDCSSensor *SensPres, TObject *SensGasComp);
  AliTPCCalibVdrift(const AliTPCCalibVdrift& source);
  virtual ~AliTPCCalibVdrift();
  AliTPCCalibVdrift& operator=(const AliTPCCalibVdrift& source);
  //
  // Interface for the reconstruction
  //
  Double_t GetPTRelative(UInt_t timeSec, Int_t side);

  //
  // Stefan interfaces - for v drift study
  //
  Double_t VdriftLinearHyperplaneApprox(Double_t dE, Double_t dT, Double_t dP, Double_t dCco2, Double_t dCn2);
  
  Double_t GetVdriftNominal();
  Double_t GetVdriftChange(Double_t x, Double_t y, Double_t z, UInt_t timeSec);

  Double_t GetMeanZVdriftChange(Double_t x, Double_t y, UInt_t timeSec);

  TGraph *MakeGraphMeanZVdriftChange(Double_t x, Double_t y, Int_t nPoints);

protected:

  AliTPCSensorTempArray *fSensTemp;   // Temperature sensors 
  AliDCSSensor          *fSensPres;   // pressure sensors
  AliTPCTempMap         *fTempMap;    // Temerature sensor map
  TObject *fSensGasComp;      // placeholder for GasConzentration infos  
  
  ClassDef(AliTPCCalibVdrift,1);

};
#endif
