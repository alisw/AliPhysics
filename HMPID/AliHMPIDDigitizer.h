#ifndef AliHMPIDDigitizer_h
#define AliHMPIDDigitizer_h
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include <AliDigitizer.h>
class AliRunDigitizer;
class TClonesArray;
class TObjArray;

class AliHMPIDDigitizer : public AliDigitizer //TObject-TNamed-TTask-AliDigitizer-AliHMPIDDigitizer
{
public:
           AliHMPIDDigitizer()                                                {}
           AliHMPIDDigitizer(AliRunDigitizer *pRunDig):AliDigitizer(pRunDig)  {}
  virtual ~AliHMPIDDigitizer()                                                {}
  void     Exec(Option_t* option=0);                //virtual
//   
  static void    Sdi2Dig(TClonesArray *pSDigLst,TObjArray *pDigLst);
protected:
  ClassDef(AliHMPIDDigitizer,0)
};    
#endif
