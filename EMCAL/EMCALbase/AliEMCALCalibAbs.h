#ifndef ALIEMCALCALIBABS_H
#define ALIEMCALCALIBABS_H

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
/// \class AliEMCALSuperModuleCalibAbs
/// \ingroup EMCALbase
/// \brief Container of APD absolute calibrations per SM
///
/// Objects of this class contain basis for absolute calibrations
/// 1 SuperModule's worth of info
///
/// total calibration factor is a product of
/// a) overall calibration factor [fAbsoluteCalib]
/// b) individual gain factor per tower [fRelativeCalib]
/// c) time-dependent correction
/// In this class we store a), b) 
///
/// \author:  David Silvermyr, <david.silvermyr@cern.ch>, ORNL
///
//////////////////////////////////////////////////////////////////////////////

// 1 SuperModule's worth of info: info on where the different APDs are
class AliEMCALSuperModuleCalibAbs : public TObject {

 public:
  
  /// Constructor, just init values.
  /// \param smNum: super module number
  AliEMCALSuperModuleCalibAbs(const int smNum=0) : TObject(),
    fSuperModuleNum(smNum),
    fCalibMethod(0),
    fCalibPass(0),
    fAbsoluteCalib(0)
    { for (int icol=0; icol<AliEMCALGeoParams::fgkEMCALCols; icol++) {
	    for (int irow=0; irow<AliEMCALGeoParams::fgkEMCALRows; irow++) {
        fRelativeCalib[icol][irow] = 1.0; } } }  // compactify, more than 3 lines must go to cxx

 public:
  
  void    SetSuperModuleNum(Int_t i)   { fSuperModuleNum    = i ; }
  Int_t   GetSuperModuleNum() const    { return fSuperModuleNum ; }
  
  void    SetCalibMethod   (Int_t i)   { fCalibMethod       = i ; }
  Int_t   GetCalibMethod   () const    { return fCalibMethod    ; } 
  
  void    SetCalibPass     (Int_t i)   { fCalibPass         = i ; }
  Int_t   GetCalibPass     ()    const { return fCalibPass      ; } 
  
  void    SetAbsoluteCalib (Float_t f) { fAbsoluteCalib     = f ; }
  Float_t GetAbsoluteCalib ()    const { return fAbsoluteCalib  ; }

  void    SetRelativeCalib(int icol, int irow, Float_t f) { fRelativeCalib[icol][irow] = f    ; }
  Float_t GetRelativeCalib(int icol, int irow)      const { return fRelativeCalib[icol][irow] ; }

 private:
  
  // first: overall values for the whole SuperModule
  Int_t fSuperModuleNum;   ///< which SuperModule is this?
  Int_t fCalibMethod;      ///< a la 0=cosmics, 1=pi0, 2=electrons,3=using ecore,
  Int_t fCalibPass;        ///< which analysis iteration is this.. 1,2,..N
  Float_t fAbsoluteCalib;  ///< (ADC>GeV absolute gain/conversion)
  
  ///< individual info for each tower, values around 1, if gains are well balanced
  Float_t fRelativeCalib[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; 

  /// \cond CLASSIMP
  ClassDef(AliEMCALSuperModuleCalibAbs, 3) ;
  /// \endcond

};
// ******* end of internal class definition *************

///////////////////////////////////////////////////////////////////////////////
///
/// \class AliEMCALCalibAbs
/// \ingroup EMCALbase
/// \brief Container of APD absolute calibrations
///
/// Objects of this class contain basis for absolute calibrations
///
/// total calibration factor is a product of
/// a) overall calibration factor [fAbsoluteCalib]
/// b) individual gain factor per tower [fRelativeCalib]
/// c) time-dependent correction
/// In this class we store a), b) 
///
/// \author:  David Silvermyr, <david.silvermyr@cern.ch>, ORNL
///
//////////////////////////////////////////////////////////////////////////////


class AliEMCALCalibAbs : public TObject 
{

public:

  enum kProblemType {kNoLED=-999}; ///< code in possible problems

  AliEMCALCalibAbs(const int nSM = AliEMCALGeoParams::fgkEMCALModules);

  // Read and Write txt I/O methods are normally not used, but are useful for 
  // filling the object before it is saved in OCDB 
  void ReadTextCalibAbsInfo (Int_t nSM, const TString &txtFileName, Bool_t swapSides=kFALSE);
  void WriteTextCalibAbsInfo(const TString &txtFileName , Bool_t swapSides=kFALSE); 
  void ReadRootCalibAbsInfo (const TString &rootFileName, Bool_t swapSides=kFALSE); 
  void ReadTreeCalibAbsInfo (TTree *tree, Bool_t swapSides=kFALSE); 
  void WriteRootCalibAbsInfo(const TString &rootFileName, Bool_t swapSides=kFALSE); 

  virtual ~AliEMCALCalibAbs();

  Int_t GetNSuperModule() const { return fNSuperModule; }; 

  // - via the index in the stored array:
  virtual AliEMCALSuperModuleCalibAbs * GetSuperModuleCalibAbsId(Int_t smIndex) const
   { return (AliEMCALSuperModuleCalibAbs*) fSuperModuleData[smIndex]; };

  // - or via the actual SM number
  virtual AliEMCALSuperModuleCalibAbs * GetSuperModuleCalibAbsNum(Int_t smNum) const;

protected:

  Int_t 	  fNSuperModule;    ///< Number of supermodules.
  TObjArray fSuperModuleData; ///< SuperModule data.

private:

  AliEMCALCalibAbs             (const AliEMCALCalibAbs &);
  AliEMCALCalibAbs &operator = (const AliEMCALCalibAbs &);

  /// \cond CLASSIMP
  ClassDef(AliEMCALCalibAbs, 4) ;
  /// \endcond

};

#endif // ALIEMCALCALIBABS_H
