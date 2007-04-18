#ifndef AliHMPIDDigitizer_h
#define AliHMPIDDigitizer_h
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//

#include <AliRun.h>
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
  void     DoNoise(Bool_t doNoise)                           {fgDoNoise=doNoise;} // Set noise or not
//   
  static void    Sdi2Dig(TClonesArray *pSDigLst,TObjArray *pDigLst);
protected:
  static  Bool_t fgDoNoise;                        // flag to switch off/on noise
  ClassDef(AliHMPIDDigitizer,0)
};    

#endif
