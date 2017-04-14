#ifndef ALIEMCALRAWRESPONSE_H
#define ALIEMCALRAWRESPONSE_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliEMCALRawResponse
/// \ingroup EMCALraw
/// \brief Handling of digits to raw data transformation
///
///  Utility Class for handling simulated digits to Raw data generation
///
/// \author Per Thomas Hille <p.t.hille@fys.uio.no>, Yale. 
/// \author Gustavo Conesa <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS. Adjustments.
//_________________________________________________________________________

#include "Rtypes.h"

class  AliEMCALRawResponse
{

public:
  
  AliEMCALRawResponse();
  
  virtual ~AliEMCALRawResponse();
  
  static Double_t RawResponseFunction(Double_t *x, Double_t *par); 
  
  static Bool_t   RawSampledResponse(Double_t dtime, Double_t damp, 
                                     Int_t * adcH, Int_t * adcL, 
                                     Int_t keyErr = 0);  
  
  static Double_t GetFEENoise()               { return fgFEENoise      ; } 
  static void     SetFEENoise(Double_t val)   { fgFEENoise = val       ; }

  static Int_t    GetPedestalValue()          { return fgPedestalValue ; }
  static void     SetPedestalValue(Int_t val) { fgPedestalValue = val  ; }
  
  static Double_t GetRawFormatTimeTrigger()   { return fgTimeTrigger   ; }
  static Int_t    GetRawFormatThreshold()     { return fgThreshold     ; }   
  
private:
  
  static Double_t fgTimeTrigger ;       ///< Time shift of the digit, apply only if not done at digitization level
  
  static Int_t    fgThreshold;          ///< Store ADC values avobe this limit demanded by AliAltoBuffer::WriteChannel()
  
  static Int_t    fgPedestalValue;      ///< Pedestal value, apply only if not done at digitization level 
  
  static Double_t fgFEENoise;           ///< Electronics noise in ADC units, apply only if not done at digitization level
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALRawResponse,1) ;
  /// \endcond

};

#endif
