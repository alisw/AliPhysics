#ifndef ALIEMCALCALIBTIME_H
#define ALIEMCALCALIBTIME_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
/// \class AliEMCALCalibTime
/// \brief Cell time shifts container class 
///
/// Cell time shifts container class. Each channel is assigned by 4 parameters, 
/// each corresponding to a BC%4.
/// An extra decalibration parameter shift foreseen.
/// Shifts are provided in ns!!
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
///_________________________________________________________________________

#include "TNamed.h"
#include "AliEMCALGeoParams.h"

class AliEMCALCalibTime: public TNamed {

 public:

  AliEMCALCalibTime();
  
  AliEMCALCalibTime(const char* name);
  
  AliEMCALCalibTime(const AliEMCALCalibTime &calibti);
  
  AliEMCALCalibTime& operator= (const AliEMCALCalibTime &calibti);
  
  virtual ~AliEMCALCalibTime() { ; }
  
  void    Reset();
  
  void    Print(Option_t *option = "") const;
  
  // Setters and getters,
  // all indexes start from 0!

  Float_t GetTimeChannel(Int_t supermodule, Int_t column, Int_t row, Int_t bc) const
  { return fTimeChannel[supermodule][column][row][bc] ; }
  
  Float_t GetTimeChannelDecal(Int_t supermodule, Int_t column, Int_t row) const
  { return fTimeChannelDecal[supermodule][column][row] ; }
  
  void    SetTimeChannel(Int_t supermodule, Int_t column, Int_t row, Int_t bc, Float_t value)
  { fTimeChannel[supermodule][column][row][bc] = value ; }
  
  void    SetTimeChannelDecal(Int_t supermodule, Int_t column, Int_t row, Float_t value)
  { fTimeChannelDecal[supermodule][column][row] = value ; }

  /// Total number of EMCAL super-modules, 
  /// 12(EMCal) + 8 (DCal)+ 2 (DCal extension possibility in future), 
  /// define it here because for Run1 releases where we had only 12 EMCal.
  static const int fgkECALDCALModules  = 22;    
  
 protected:

  /// time shift (ns) of one ADC channel ([mod][col][row])
  Float_t  fTimeChannelDecal    [fgkECALDCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]    ; 
  
  /// time shift (ns) of one ADC channel ([mod][col][row][bunch crossing number])
  Float_t  fTimeChannel         [fgkECALDCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows][4] ; 

  /// \cond CLASSIMP
  ClassDef(AliEMCALCalibTime,1) ;   
  /// \endcond

};

#endif // ALIEMCALCALIBTIME_H
