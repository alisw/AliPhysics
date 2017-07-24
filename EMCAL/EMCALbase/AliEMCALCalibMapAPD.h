#ifndef ALIEMCALCALIBMAPAPD_H
#define ALIEMCALCALIBMAPAPD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <TObjArray.h>
#include "AliEMCALGeoParams.h"
#include <cmath>
class TString;
class TTree;

// ******* internal class definition *************

///////////////////////////////////////////////////////////////////////////////
///
/// \class AliEMCALCalibMapAPDVal
/// \ingroup EMCALbase
/// \brief Container of calibration values container per single APD
///
/// Calibration values container per single APD
///
/// \author:  David Silvermyr, <david.silvermyr@cern.ch>, ORNL
///
//////////////////////////////////////////////////////////////////////////////

class AliEMCALCalibMapAPDVal : public TObject
{

 public:
  
  AliEMCALCalibMapAPDVal() : TObject(), 
    fHardWareId(0),
    fAPDNum(0),
    fV30(0),
    fBreakDown(0),
    fDarkCurrent(0) 
    { Init() ; }
  
  void Init() {
    fHardWareId = 0;
    fAPDNum = 0;
    fV30 = 0;
    fBreakDown = 0;
    fDarkCurrent = 0; 
    for (int i=0; i<3; i++) 
    {
      fPar[i] = 0;
      fParErr[i] = 0;
    }
    return;
  }
  
 public:
  void    SetHardWareId(Int_t i)       { fHardWareId    = i ; }
  Int_t   GetHardWareId()        const { return fHardWareId ; }
 
  void    SetAPDNum(Int_t i)           { fAPDNum        = i ; }
  Int_t   GetAPDNum()            const { return fAPDNum     ; }
  
  void    SetV30(Float_t f)            { fV30           = f ; }
  Float_t GetV30()               const { return fV30; } 
  
  void    SetPar(int ip, Float_t f)    { fPar[ip]       = f ; }
  Float_t GetPar(int ip)         const { return fPar[ip]    ; } 
 
  void    SetParErr(int ip, Float_t f) { fParErr[ip]    = f ; }
  Float_t GetParErr(int ip)      const { return fParErr[ip] ; } 
  
  void    SetBreakDown(Int_t i)        { fBreakDown     = i ; }
  Int_t   GetBreakDown()         const { return fBreakDown  ; }
  
  void    SetDarkCurrent(Float_t f)    { fDarkCurrent   = f ; }
  Float_t GetDarkCurrent()       const { return fDarkCurrent; } 

 private:
  
  Int_t   fHardWareId;   ///< HardWare index
  
  // info from APD calibrations
  Int_t   fAPDNum;      ///< assigned APD-PA number; Catania 10000-, Houston: 20000-
  Float_t fV30;         ///< Catania/Houston Voltage V30 (V) at T = 25 deg C
  Float_t fPar[3];      ///< fit parameters, p0,p1,p2 - for ADC vs bias measurement
  Float_t fParErr[3];   ///< error on fit parameters	
  
  Int_t   fBreakDown;   ///< Hamamatsu Breakdown Voltage (V)	
  Float_t fDarkCurrent; ///< Hamamatsu Dark Current (A)	
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALCalibMapAPDVal, 2) ;
  /// \endcond

}; // AliEMCALCalibAPDVal

// 1 SuperModule's worth of info: info on where the different APDs are

///////////////////////////////////////////////////////////////////////////////
///
/// \class AliEMCALSuperModuleCalibMapAPD
/// \ingroup EMCALbase
/// \brief Container of calibration values container per SM
///
/// calibration values container per SM
/// 1 SuperModule's worth of info
///
/// \author:  David Silvermyr, <david.silvermyr@cern.ch>, ORNL
///
//////////////////////////////////////////////////////////////////////////////

class AliEMCALSuperModuleCalibMapAPD : public TObject
{
  
 public:

  /// Constructor, just init values.
  /// \param smNum: super module number
  AliEMCALSuperModuleCalibMapAPD(const int smNum=0) : TObject(), 
    fSuperModuleNum(smNum)
    { for (int icol=0; icol<AliEMCALGeoParams::fgkEMCALCols; icol++) {
	    for (int irow=0; irow<AliEMCALGeoParams::fgkEMCALRows; irow++) {
	  fAPDVal[icol][irow].Init(); } } } // compactify, more than 3 lines must go to cxx
  
 public:
  
  void  SetSuperModuleNum(Int_t i) { fSuperModuleNum = i    ; } 
  Int_t GetSuperModuleNum()  const { return fSuperModuleNum ; } 
  
  AliEMCALCalibMapAPDVal * GetAPDVal(int icol, int irow) 
  { return &fAPDVal[icol][irow] ; }

 private:
  
  Int_t fSuperModuleNum; // SuperModule index
  AliEMCALCalibMapAPDVal fAPDVal[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; // APD calibration info
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALSuperModuleCalibMapAPD, 2) ;
  /// \endcond

};
// ******* end of internal class definition *************    
    
///////////////////////////////////////////////////////////////////////////////
///
/// \class AliEMCALCalibMapAPD
/// \ingroup EMCALbase
/// \brief Container of calibration values container 
///
/// Objects of this class contain info on APD calibration and map info,
/// such as V30 and other parameters from the tests in Catania/Houston,
/// as well as info on which APD is located where.  
///
/// \author:  David Silvermyr, <david.silvermyr@cern.ch>, ORNL
///
//////////////////////////////////////////////////////////////////////////////

class AliEMCALCalibMapAPD : public TObject 
{

public:

  enum kValType {kCalibTemp=25}; ///< 25 deg C used for all APD calibrations

  AliEMCALCalibMapAPD(const int nSM = AliEMCALGeoParams::fgkEMCALModules);

  // Read and Write txt I/O methods are normally not used, but are useful for 
  // filling the object before it is saved in OCDB 
  void ReadTextCalibMapAPDInfo (Int_t nSM, const TString &txtFileName, Bool_t swapSides=kFALSE); 
  void WriteTextCalibMapAPDInfo(const TString &txtFileName , Bool_t swapSides=kFALSE); 
  void ReadRootCalibMapAPDInfo (const TString &rootFileName, Bool_t swapSides=kFALSE);
  void ReadTreeCalibMapAPDInfo (TTree *tree, Bool_t swapSides=kFALSE);
  void WriteRootCalibMapAPDInfo(const TString &rootFileName, Bool_t swapSides=kFALSE); 

  virtual ~AliEMCALCalibMapAPD();

  Int_t GetNSuperModule() const { return fNSuperModule; }; 

  virtual AliEMCALSuperModuleCalibMapAPD * GetSuperModuleCalibMapAPDId(Int_t smIndex) const
    { return (AliEMCALSuperModuleCalibMapAPD*) fSuperModuleData[smIndex]; }

  virtual AliEMCALSuperModuleCalibMapAPD * GetSuperModuleCalibMapAPDNum(Int_t smNum) const;   
  
  /// Method to calculate gain M from fit parameters, and HV value
  Float_t GetGain(Float_t fitPar[3], Float_t hv) const 
    { return (fitPar[0] + fitPar[1] * exp(fitPar[2]*hv)); }

  /// Method to calculate gain M from fit parameters, and HV value
  Float_t GetGain(Float_t fitPar0, Float_t fitPar1, Float_t fitPar2, Float_t hv) const 
    { return (fitPar0 + fitPar1 * exp(fitPar2*hv)); }

protected:

  Int_t 	  fNSuperModule;    ///< Number of supermodules.
  TObjArray fSuperModuleData; ///< SuperModule data.

private:

  AliEMCALCalibMapAPD             (const AliEMCALCalibMapAPD &);
  AliEMCALCalibMapAPD &operator = (const AliEMCALCalibMapAPD &);

  /// \cond CLASSIMP
  ClassDef(AliEMCALCalibMapAPD, 3) ;
  /// \endcond

};

#endif // ALIEMCALCALIBMAPAPD_H
