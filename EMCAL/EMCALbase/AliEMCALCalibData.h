#ifndef ALIEMCALCALIBDATA_H
#define ALIEMCALCALIBDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
/// \class AliEMCALCalibData
/// \brief Cell energy calibration factors container class 
///
/// Channel energy calibration factors (ADC to GeV conversion) and pedestal, 
/// An extra decalibration parameter factor foreseen.
///
/// This container also includes arrays for time calibration, but this is not
/// under use, this functionality is in AliEMCALCalibTime. 
/// It is kept for backward compatibility reasons.
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
///_________________________________________________________________________

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
  Float_t GetADCchannelOnline(Int_t module, Int_t column, Int_t row) const;
  Float_t GetADCchannelDecal (Int_t module, Int_t column, Int_t row) const;
  Float_t GetADCpedestal     (Int_t module, Int_t column, Int_t row) const;

  Float_t GetADCchannelRef   () const { return fADCchannelRef ; }

  void    SetADCchannel      (Int_t module, Int_t column, Int_t row, Float_t value);
  void    SetADCchannelOnline(Int_t module, Int_t column, Int_t row, Float_t value);
  void    SetADCchannelDecal (Int_t module, Int_t column, Int_t row, Float_t value);
  void    SetADCpedestal     (Int_t module, Int_t column, Int_t row, Float_t value);
  
  void    SetADCchannelRef   (Float_t value) { fADCchannelRef = value ; }
  
  // Do not use, please use AliEMCALTimeCalib, keep for backward compatibility reasons  
  Float_t GetTimeChannelDecal(Int_t module, Int_t column, Int_t row) const;
  Float_t GetTimeChannel     (Int_t module, Int_t column, Int_t row, Int_t bc) const;
  
  void    SetTimeChannelDecal(Int_t module, Int_t column, Int_t row, Float_t value);
  void    SetTimeChannel     (Int_t module, Int_t column, Int_t row, Int_t bc, Float_t value);


  static const int fgkECALModules  = 12;   // number of modules in EMCAL
  static const int fgkDCALModules  = 10;   // number of modules in DCAL 8+2 in possible future
  
 protected:

  Float_t  fADCchannelRef ;  /// base value of the ADC channel set from cosmics calibration, not to be used, instead use fADCchannelOnline

  Float_t  fADCchannel          [fgkECALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]    ; /// width of one ADC channel in GeV ([mod][col][row])
  Float_t  fADCchannelOnline    [fgkECALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]    ; /// width of one ADC channel in GeV obtained from the voltage settings online
  Float_t  fADCchannelDecal     [fgkECALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]    ; /// decalibrate width of one ADC channel in GeV ([mod][col][row])
  Float_t  fADCpedestal         [fgkECALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]    ; /// value of the  ADC pedestal ([mod][col][row]), not used
                                                                                                                        
  // Do not use, please use AliEMCALTimeCalib, keep for backward compatibility reasons  
  Float_t  fTimeChannelDecal    [fgkECALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]    ; /// time shift of one ADC channel ([mod][col][row])
  Float_t  fTimeChannel         [fgkECALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows][4] ; /// time shift of one ADC channel ([mod][col][row][bunch crossing number])

  // Add specific arrays for DCal to avoid backward incompatibilities,
  // dimension of DCal SM is smaller than EMCAL but assume the same to avoid complications due to partial SM
  Float_t  fADCchannelDCAL      [fgkDCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]    ; /// width of one ADC channel in GeV ([mod][col][row])
  Float_t  fADCchannelOnlineDCAL[fgkDCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]    ; /// width of one ADC channel in GeV  obtained from the voltage settings online
  Float_t  fADCchannelDecalDCAL [fgkDCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]    ; /// decalibrate width of one ADC channel in GeV ([mod][col][row])
  Float_t  fADCpedestalDCAL     [fgkDCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]    ; /// value of the  ADC pedestal ([mod][col][row]), not used
 
  // Do not use, please use AliEMCALTimeCalib, keep for backward compatibility reasons  
  Float_t  fTimeChannelDecalDCAL[fgkDCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]    ; /// time shift of one ADC channel ([mod][col][row])
  Float_t  fTimeChannelDCAL     [fgkDCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows][4] ; /// time shift of one ADC channel ([mod][col][row][bunch crossing number])

  /// \cond CLASSIMP
  ClassDef(AliEMCALCalibData,6);   
  /// \endcond

};

#endif // ALIEMCALCALIBDATA_H


