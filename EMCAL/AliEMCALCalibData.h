#ifndef ALIEMCALCALIBDATA_H
#define ALIEMCALCALIBDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//
//  class for EMCAL calibration               //
//
////////////////////////////////////////////////

#include "TNamed.h"
#include "AliEMCALGeoParams.h"

class AliEMCALCalibData: public TNamed {

 public:

  AliEMCALCalibData();
  AliEMCALCalibData(const char* name);
  AliEMCALCalibData(const AliEMCALCalibData &calibda);
  AliEMCALCalibData& operator= (const AliEMCALCalibData &calibda);
  virtual ~AliEMCALCalibData() { ; }
  
  void    Reset();
  void    Print(Option_t *option = "") const;
  
  // All indexes start from 0!
  Float_t GetADCchannel      (Int_t module, Int_t column, Int_t row) const;
  Float_t GetADCchannelDecal (Int_t module, Int_t column, Int_t row) const;
  Float_t GetADCpedestal     (Int_t module, Int_t column, Int_t row) const;
  Float_t GetTimeChannelDecal(Int_t module, Int_t column, Int_t row) const;
  Float_t GetTimeChannel     (Int_t module, Int_t column, Int_t row, Int_t bc) const;
	
  Float_t GetADCchannelRef   () const { return fADCchannelRef ; }

  void    SetADCchannel      (Int_t module, Int_t column, Int_t row, Float_t value);
  void    SetADCchannelDecal (Int_t module, Int_t column, Int_t row, Float_t value);
  void    SetADCpedestal     (Int_t module, Int_t column, Int_t row, Float_t value);
  void    SetTimeChannelDecal(Int_t module, Int_t column, Int_t row, Float_t value);
  void    SetTimeChannel     (Int_t module, Int_t column, Int_t row, Int_t bc, Float_t value);

  void    SetADCchannelRef   (Float_t value) { fADCchannelRef = value ; }

 protected:

  Float_t  fADCchannelRef ;  // base value of the ADC channel set from cosmics calibration

  Float_t  fADCchannel      [AliEMCALGeoParams::fgkEMCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows] ; // width of one ADC channel in GeV ([mod][col][row])
  Float_t  fADCchannelDecal [AliEMCALGeoParams::fgkEMCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows] ; // decalibrate width of one ADC channel in GeV ([mod][col][row])
  Float_t  fADCpedestal     [AliEMCALGeoParams::fgkEMCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows] ; // value of the  ADC pedestal ([mod][col][row]), not used
  Float_t  fTimeChannelDecal[AliEMCALGeoParams::fgkEMCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows] ;    // time shift of one ADC channel ([mod][col][row])
  Float_t  fTimeChannel     [AliEMCALGeoParams::fgkEMCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows][4] ; // time shift of one ADC channel ([mod][col][row][bunch crossing number])

  ClassDef(AliEMCALCalibData,4)    // EMCAL Calibration data
};

#endif
