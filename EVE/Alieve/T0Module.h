#ifndef ALIEVE_T0Module_H
#define ALIEVE_T0Module_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// The main AliEVE drawing module for the T0 detector                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <Reve/QuadSet.h>
#include <AliT0digit.h>
#include <AliT0RawReader.h>

namespace Alieve {
 
class T0Module : public Reve::QuadSet
{
 
  T0Module(const T0Module&);
  T0Module& operator=(const T0Module&);

public:
 
  T0Module(const Text_t* n="T0Module", Int_t sigType=0, AliT0digit *digits=0,AliT0RawReader *start=0);
  virtual ~T0Module();

  virtual void DigitSelected(Int_t idx);

  void LoadRaw(TString fileName, Int_t ievt);

  static void MakeModules(AliT0digit *digits);

protected:
  Int_t           fSigType; // 0 ~ ADC, 1 ~ TDC
  AliT0digit     *fDigits;
  AliT0RawReader *fStart;

   ClassDef(T0Module,1); 
};

/*
 class T0ModuleTDC : public T0Module
 {
 public:
   // constructor

    virtual void QuadSelected(Int_t idx);
 };
*/

}
#endif
