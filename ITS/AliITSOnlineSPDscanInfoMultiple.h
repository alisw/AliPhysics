#ifndef ALI_ITS_ONLINESPDSCANINFOMULTIPLE_H
#define ALI_ITS_ONLINESPDSCANINFOMULTIPLE_H  

/////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                      //
// This class is used as a container online.                   //
// It holds additional information needed for a scan with      //
// multiple steps. (dac scan, min thr. mean thr. etc.          //
// This class should only be used through the interface of the //
// AliITSOnlineSPDscanMultiple class.                          //
/////////////////////////////////////////////////////////////////

#include "AliITSOnlineSPDscanInfo.h"
#include "TArrayI.h"

class AliITSOnlineSPDscanInfoMultiple :  public AliITSOnlineSPDscanInfo {

 public:
  AliITSOnlineSPDscanInfoMultiple();
  virtual ~AliITSOnlineSPDscanInfoMultiple();

  virtual UInt_t AddScanStep(); // returns the index (nsi) of the added step

  void    SetDacId(Int_t val){fDacId=val;}
  void    SetDacValue(UInt_t nsi, Int_t val);
	  
  Int_t   GetDacId() const {return fDacId;}
  Int_t   GetDacValue(UInt_t nsi) const;


 protected:
  Int_t     fDacId;         // id of DAC used for the scan
  TArrayI   fDacValues;     // DAC values for each step


  ClassDef(AliITSOnlineSPDscanInfoMultiple,1)
    };

#endif
