//____________________________________________________________________
//
// $Id$
//
// Script that contains a class to compare the raw data written to the
// digits it's created from.
//
// Use the script `Compile.C' to compile this class using ACLic. 
//
#include <AliLog.h>
#include <AliFMDDigit.h>
#include <AliFMDInput.h>
#include <AliFMDUShortMap.h>
#include <AliFMDParameters.h>
#include <iostream>
#include <TStyle.h>
#include <TArrayF.h>
#include <TCanvas.h>

class CheckRaw : public AliFMDInput
{
public:
  CheckRaw()
  {
    AddLoad(kDigits);
    AddLoad(kRaw);
  }
  Bool_t Init() 
  {
    Bool_t ret = AliFMDInput::Init();
    // AliFMDGeometry* geom = AliFMDGeometry::Instance();
    // geom->Init();
    // geom->InitTransformations();
    AliFMDParameters* param = AliFMDParameters::Instance();
    param->Init();
    return ret;
  }
  Bool_t ProcessDigit(AliFMDDigit* digit)
  {
    // Cache the energy loss 
    if (!digit) return kFALSE;
    UShort_t det = digit->Detector();
    Char_t   rng = digit->Ring();
    UShort_t sec = digit->Sector();
    UShort_t str = digit->Strip();
    if (str > 511) {
      AliWarning(Form("Bad strip number %d in digit", str));
      return kTRUE;
    }
    fMap(det, rng, sec, str) = digit->Counts();
    return kTRUE;
  }
  Bool_t ProcessRawDigit(AliFMDDigit* digit)
  {
    // Cache the energy loss 
    if (!digit) return kFALSE;
    UShort_t det = digit->Detector();
    Char_t   rng = digit->Ring();
    UShort_t sec = digit->Sector();
    UShort_t str = digit->Strip();
    if (str > 511) {
      AliWarning(Form("Bad strip number %d in digit", str));
      return kTRUE;
    }
    if (digit->Counts() != fMap(det, rng, sec, str) && 
	fMap(det, rng, sec, str) != 1024) {
      AliWarning(Form("Mismatch in digit FMD%d%c[%2d,%3d] %d != %d", 
		      det, rng, sec, str, digit->Counts(), 
		      fMap(det, rng, sec, str)));
      return kTRUE;
    }
    AliDebug(1, Form("Raw digit FMD%d%c[%2d,%3D] is good", 
		     det, rng, sec, str));
    return kTRUE;
  }
protected:
  AliFMDUShortMap fMap;
};


//____________________________________________________________________
//
// EOF
//


  
  
