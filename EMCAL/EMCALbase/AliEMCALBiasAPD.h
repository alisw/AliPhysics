#ifndef ALIEMCALBIASAPD_H
#define ALIEMCALBIASAPD_H

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
/// \class AliEMCALSuperModuleBiasAPD
/// \ingroup EMCALbase
/// \brief APD bias info storage per SM
///
/// Objects of this class contain info on APD bias settings/voltages
/// 1 SuperModule's worth of info
///
/// \author:  David Silvermyr, <david.silvermyr@cern.ch>, ORNL
///
//////////////////////////////////////////////////////////////////////////////

class AliEMCALSuperModuleBiasAPD : public TObject 
{

public:
  
  /// Constructor, just init values.
  /// \param smNum: super module number
  AliEMCALSuperModuleBiasAPD(const int smNum=0) : TObject(), 
  fSuperModuleNum(smNum)
  { for (int icol=0; icol<AliEMCALGeoParams::fgkEMCALCols; icol++) { 
    for (int irow=0; irow<AliEMCALGeoParams::fgkEMCALRows; irow++) {
      fElecId[icol][irow] = 0; fDAC[icol][irow] = 0; fVoltage[icol][irow] = 0; } } } // compactify, more than 3 lines must go to cxx
  
public:
  
  void    SetSuperModuleNum(Int_t i)                { fSuperModuleNum        = i ; } 
  Int_t   GetSuperModuleNum()                 const { return fSuperModuleNum     ; } 
  
  void    SetElecId (int icol, int irow, Int_t i)   { fElecId[icol][irow]    = i ; }
  Int_t   GetElecId (int icol, int irow)      const { return fElecId[icol][irow] ; }
  
  void    SetDAC    (int icol, int irow, Int_t i)   { fDAC[icol][irow]       = i ; }
  Int_t   GetDAC    (int icol, int irow)      const { return fDAC[icol][irow]    ; }
  
  void    SetVoltage(int icol, int irow, Float_t f) { fVoltage[icol][irow]   = f ; }
  Float_t GetVoltage(int icol, int irow)      const { return fVoltage[icol][irow]; }
  
private:
  
  Int_t fSuperModuleNum;                                                              ///< SM index
  Int_t fElecId   [AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; ///< ElectronicsIndex/Address - we keep this to help ensure that the column/row info matches with electronics indices
  Int_t fDAC      [AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; ///< 0-0x3ff register
  Float_t fVoltage[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows]; ///< 210 to ca 417 V. (function of DAC setting)
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALSuperModuleBiasAPD, 2) ;
  /// \endcond

};
// ******* end of internal class definition *************

///////////////////////////////////////////////////////////////////////////////
///
/// \class AliEMCALBiasAPD
/// \ingroup EMCALbase
/// \brief APD bias info storage
///
/// Objects of this class contain info on APD bias settings/voltages
///
/// \author:  David Silvermyr, <david.silvermyr@cern.ch>, ORNL
///
//////////////////////////////////////////////////////////////////////////////

class AliEMCALBiasAPD : public TObject 
{

public:

  AliEMCALBiasAPD(const int nSM = AliEMCALGeoParams::fgkEMCALModules);

  // Read and Write txt I/O methods are normally not used, but are useful for 
  // filling the object before it is saved in OCDB 
  void ReadTextBiasAPDInfo (Int_t nSM, const TString &txtFileName, Bool_t swapSides=kFALSE); 
  void WriteTextBiasAPDInfo(const TString &txtFileName , Bool_t swapSides=kFALSE); 
  void ReadRootBiasAPDInfo (const TString &rootFileName, Bool_t swapSides=kFALSE); 
  void ReadTreeBiasAPDInfo (TTree *tree, Bool_t swapSides=kFALSE);
  void WriteRootBiasAPDInfo(const TString &rootFileName, Bool_t swapSides=kFALSE); 

  virtual ~AliEMCALBiasAPD();

  Int_t GetNSuperModule() const { return fNSuperModule ; }

  // - via the index in the stored array:
  virtual AliEMCALSuperModuleBiasAPD * GetSuperModuleBiasAPDId(Int_t smIndex) const
    { return (AliEMCALSuperModuleBiasAPD*) fSuperModuleData[smIndex] ; }

  // - or via the actual SM number
  virtual AliEMCALSuperModuleBiasAPD * GetSuperModuleBiasAPDNum(Int_t smNum) const;

protected:

  Int_t 	  fNSuperModule;    ///< Number of supermodules.
  TObjArray fSuperModuleData; ///< SuperModule data.

private:

  AliEMCALBiasAPD             (const AliEMCALBiasAPD &);
  AliEMCALBiasAPD &operator = (const AliEMCALBiasAPD &);

  /// \cond CLASSIMP
  ClassDef(AliEMCALBiasAPD, 3) ;
  /// \endcond

};

#endif // ALIEMCALBIASAPD_H
