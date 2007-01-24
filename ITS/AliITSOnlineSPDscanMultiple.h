#ifndef ALI_ITS_ONLINESPDSCANMULTIPLE_H
#define ALI_ITS_ONLINESPDSCANMULTIPLE_H  

////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                 //
// Interface class to the containers of an online scan    //
// with multiple steps.                                   //
////////////////////////////////////////////////////////////

#include "AliITSOnlineSPDscan.h"

class AliITSOnlineSPDscanMultiple :  public AliITSOnlineSPDscan {

 public:
  AliITSOnlineSPDscanMultiple();
  AliITSOnlineSPDscanMultiple(Char_t *fileName);
  AliITSOnlineSPDscanMultiple(const AliITSOnlineSPDscanMultiple& scan);
  virtual ~AliITSOnlineSPDscanMultiple();
  AliITSOnlineSPDscanMultiple& operator=(const AliITSOnlineSPDscanMultiple& scan);

  //  virtual void   ReadFromTObjArray(TObjArray *arr);
  virtual UInt_t AddScanStep();

  void    SetDacId(Int_t val);
  void    SetDacValue(UInt_t nsi, Int_t val);
	  
  Int_t   GetDacId();
  Int_t   GetDacValue(UInt_t nsi);


 protected:

  ClassDef(AliITSOnlineSPDscanMultiple,1)
    };

#endif
