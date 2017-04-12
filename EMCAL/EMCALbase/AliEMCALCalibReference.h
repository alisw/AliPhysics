#ifndef ALIEMCALCALIBREFERENCE_H
#define ALIEMCALCALIBREFERENCE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <TObjArray.h>
#include "AliEMCALGeoParams.h"
class TString;
class TTree;

// ******* internal class definition *************

///////////////////////////////////////////////////////////////////////////////
///
/// \class AliEMCALCalibReferenceVal
/// \ingroup EMCALbase
/// \brief Container of calibration reference values container per single APD
///
/// Calibration reference values container per single APD
///
/// \author:  David Silvermyr, <david.silvermyr@cern.ch>, ORNL
///
//////////////////////////////////////////////////////////////////////////////

class AliEMCALCalibReferenceVal : public TObject 
{

 public:

  /// Constructor
  AliEMCALCalibReferenceVal() : TObject(), // just init values
    fHighLow(0),
    fLEDAmp(0),
    fLEDAmpRMS(0)
    { }
  
  void Init() {
    fHighLow = 0;
    fLEDAmp = 0;
    fLEDAmpRMS = 0;
    return;
  }

 public:
  
  void    SetHighLow(Int_t i)     { fHighLow      = i ; }
  Int_t   GetHighLow()      const { return fHighLow   ; }
  
  void    SetLEDAmp(Float_t f)    { fLEDAmp       = f ; }
  Float_t GetLEDAmp()       const { return fLEDAmp    ; }
  
  void    SetLEDAmpRMS(Float_t f) { fLEDAmpRMS    = f ; }
  Float_t GetLEDAmpRMS()    const { return fLEDAmpRMS ; }

 private:
  Int_t   fHighLow  ; ///< 0 (low) or 1 (high) gain, used for LEDAmp info
  Float_t fLEDAmp   ; ///< LED amplitude
  Float_t fLEDAmpRMS; ///< RMS

  /// \cond CLASSIMP
  ClassDef(AliEMCALCalibReferenceVal, 1) ;
  /// \endcond

};

///////////////////////////////////////////////////////////////////////////////
///
/// \class AliEMCALSuperModuleCalibReference 
/// \ingroup EMCALbase
/// \brief Container of calibration reference values container per SM
///
/// calibration reference values container per SM
/// 1 SuperModule's worth of info
///
/// \author:  David Silvermyr, <david.silvermyr@cern.ch>, ORNL
///
//////////////////////////////////////////////////////////////////////////////

class AliEMCALSuperModuleCalibReference : public TObject 
{

 public:

  /// Constructor, just init values.
  /// \param smNum: super module number
  AliEMCALSuperModuleCalibReference(const int smNum=0) : TObject(), 
    fSuperModuleNum(smNum),
    fReferenceTime(0)
    { for (int iref=0; iref<AliEMCALGeoParams::fgkEMCALLEDRefs; iref++) 
      { fLEDRefAmp[iref] = 0; fLEDRefAmpRMS[iref] = 0; fLEDRefHighLow[iref] = 0; }
      for (int itemp=0; itemp<AliEMCALGeoParams::fgkEMCALTempSensors; itemp++) 
      { fTemperature[itemp] = 0; fTemperatureRMS[itemp] = 0; }
      for (int icol=0; icol<AliEMCALGeoParams::fgkEMCALCols; icol++) {
      for (int irow=0; irow<AliEMCALGeoParams::fgkEMCALRows; irow++) {
	  fAPDVal[icol][irow].Init(); } } }

 public:
  
  // first
  void    SetSuperModuleNum(Int_t i)               { fSuperModuleNum            = i ; }
  Int_t   GetSuperModuleNum()                const { return fSuperModuleNum         ; }
  void    SetReferenceTime(Int_t i)                { fReferenceTime             = i ; }
  Int_t   GetReferenceTime()                 const { return fReferenceTime          ; }

  // second
  void    SetLEDRefAmp    (int iLEDRef, Float_t f) { fLEDRefAmp[iLEDRef]        = f ; }
  Float_t GetLEDRefAmp    (int iLEDRef)      const { return fLEDRefAmp[iLEDRef]     ; } 
  void    SetLEDRefAmpRMS (int iLEDRef, Float_t f) { fLEDRefAmpRMS[iLEDRef]     = f ; } 
  Float_t GetLEDRefAmpRMS (int iLEDRef)      const { return fLEDRefAmpRMS[iLEDRef]  ; } 
  void    SetLEDRefHighLow(int iLEDRef, Int_t i)   { fLEDRefHighLow[iLEDRef]    = i ; } 
  Int_t   GetLEDRefHighLow(int iLEDRef)      const { return fLEDRefHighLow[iLEDRef] ; }

  void    SetTemperature   (int itemp, Float_t f)  { fTemperature[itemp]        = f ; }
  Float_t GetTemperature   (int itemp) const       { return fTemperature[itemp]     ; }
  void    SetTemperatureRMS(int itemp, Float_t f)  { fTemperatureRMS[itemp]     = f ; } 
  Float_t GetTemperatureRMS(int itemp)       const { return fTemperatureRMS[itemp]  ; } 

  // third
  AliEMCALCalibReferenceVal * GetAPDVal(int icol, int irow) 
    { return &fAPDVal[icol][irow]; }

 private:
  
  // first: overall values for the whole SuperModule
  Int_t fSuperModuleNum ;                                           ///< which SuperModule is this?
  Int_t fReferenceTime  ;                                           ///< t0, unix timestamp
  
  // second: additional info for LED Reference and SM temperature
  Float_t fLEDRefAmp     [AliEMCALGeoParams::fgkEMCALLEDRefs] ;     ///< additional info for LED Reference, LED amplitude at  t0, low gain equivalent
  Float_t fLEDRefAmpRMS  [AliEMCALGeoParams::fgkEMCALLEDRefs] ;     ///< additional info for LED Reference, RMS
  Int_t   fLEDRefHighLow [AliEMCALGeoParams::fgkEMCALLEDRefs] ;     ///< additional info for LED Reference, 0 (low) or 1 (high) gain
  
  Float_t fTemperature   [AliEMCALGeoParams::fgkEMCALTempSensors] ; ///< additional info for SM temperature, temperature at t0
  Float_t fTemperatureRMS[AliEMCALGeoParams::fgkEMCALTempSensors] ; ///< additional info for SM temperature, RMS
  
  // third: individual info for each tower
  AliEMCALCalibReferenceVal fAPDVal[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; // individual info for each tower at t0

  /// \cond CLASSIMP
  ClassDef(AliEMCALSuperModuleCalibReference, 1) ;
  /// \endcond
  
};
// ******* end of internal class definition *************

///////////////////////////////////////////////////////////////////////////////
///
/// \class AliEMCALCalibReference
/// \ingroup EMCALbase
/// \brief Container of reference calibrations
///
/// Objects of this class contain basis for reference calibrations
///
/// total calibration factor is a product of
/// a) overall calibration factor [fAbsoluteCalib]
/// b) individual gain factor per tower [fRelativeCalib]
/// c) time-dependent correction
/// In this class we store the needed static ingredients for c)
///
/// \author:  David Silvermyr, <david.silvermyr@cern.ch>, ORNL
///
//////////////////////////////////////////////////////////////////////////////

class AliEMCALCalibReference : public TObject 
{
  
public:

  enum kProblemType {kNoLED=-999}; ///< code in possible problems

  AliEMCALCalibReference(const int nSM = AliEMCALGeoParams::fgkEMCALModules);

  // Read and Write txt I/O methods are normally not used, but are useful for 
  // filling the object before it is saved in OCDB 
  void ReadTextCalibReferenceInfo (Int_t nSM, const TString &txtFileName, Bool_t swapSides=kFALSE); 
  void WriteTextCalibReferenceInfo(const TString &txtFileName , Bool_t swapSides=kFALSE); 
  void ReadRootCalibReferenceInfo (const TString &rootFileName, Bool_t swapSides=kFALSE); 
  void ReadTreeCalibReferenceInfo (TTree *tree, Bool_t swapSides=kFALSE); 
  void WriteRootCalibReferenceInfo(const TString &rootFileName, Bool_t swapSides=kFALSE); 

  virtual ~AliEMCALCalibReference();

  // pointer to stored info.
  Int_t GetNSuperModule() const { return fNSuperModule; } 

  // - via the index in the stored array:
  virtual AliEMCALSuperModuleCalibReference * GetSuperModuleCalibReferenceId(Int_t smIndex) const
   { return (AliEMCALSuperModuleCalibReference*) fSuperModuleData[smIndex]; }

  // - or via the actual SM number
  virtual AliEMCALSuperModuleCalibReference * GetSuperModuleCalibReferenceNum(Int_t smNum) const;

protected:

  Int_t 	  fNSuperModule    ; ///< Number of supermodules.
  TObjArray fSuperModuleData ; ///< SuperModule data.

private:

  AliEMCALCalibReference             (const AliEMCALCalibReference &);
  AliEMCALCalibReference &operator = (const AliEMCALCalibReference &);
 
  /// \cond CLASSIMP
  ClassDef(AliEMCALCalibReference, 1) ;
  /// \endcond

};

#endif // ALIEMCALCALIBREFERENCE_H
