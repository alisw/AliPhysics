#ifndef AliHMPIDDigitizer_h
#define AliHMPIDDigitizer_h
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//.
//.
//.

#include <AliDigitizer.h>

class TClonesArray;
class TObjArray;

class AliHMPIDDigitizer : public AliDigitizer //TObject-TNamed-AliDigitizer-AliHMPIDDigitizer
{
public:
           AliHMPIDDigitizer()                                                {}
           AliHMPIDDigitizer(AliDigitizationInput *digInp):AliDigitizer(digInp)  {}
  virtual ~AliHMPIDDigitizer()                                                {}
  void     Digitize(Option_t* option=0);                //virtual
  void     DoNoise(Bool_t doNoise)                           {fgDoNoise=doNoise;} // Set noise or not
//   
  static void    Sdi2Dig(TClonesArray *pSDigLst,TObjArray *pDigLst);
protected:
  static  Bool_t fgDoNoise;                        // flag to switch off/on noise
  ClassDef(AliHMPIDDigitizer,0)
};    

#endif
