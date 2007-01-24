#ifndef ALI_ITS_ONLINESPDSCANINFOMEANTH_H
#define ALI_ITS_ONLINESPDSCANINFOMEANTH_H  

/////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                      //
// This class is used as a container online.                   //
// It holds additional information needed for a mean threshold //
// scan.                                                       //
// This class should only be used through the interface of the //
// AliITSOnlineSPDscanMeanTh class.                            //
/////////////////////////////////////////////////////////////////

#include "AliITSOnlineSPDscanInfoMultiple.h"

class AliITSOnlineSPDscanInfoMeanTh :  public AliITSOnlineSPDscanInfoMultiple {

 public:
  AliITSOnlineSPDscanInfoMeanTh();
  virtual ~AliITSOnlineSPDscanInfoMeanTh();

  virtual UInt_t AddScanStep();

  void     SetDacLow(UInt_t nsi, UInt_t hs, Int_t val);
  void     SetDacHigh(UInt_t nsi, UInt_t hs, Int_t val);
  void     SetTPAmp(UInt_t nsi, UInt_t hs, Int_t val);

  Int_t    GetDacLow(UInt_t nsi, UInt_t hs) const;
  Int_t    GetDacHigh(UInt_t nsi, UInt_t hs) const;
  Int_t    GetTPAmp(UInt_t nsi, UInt_t hs) const;

 protected:
  TArrayI   fDacLow[6];        // DAC low values for each step
  TArrayI   fDacHigh[6];       // DAC high values for each step
  TArrayI   fTPAmps[6];        // test pulse amplitudes, one for each step

  ClassDef(AliITSOnlineSPDscanInfoMeanTh,1)
    };
    
#endif
