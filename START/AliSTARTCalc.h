#ifndef AliSTARTCalc_H
#define AliSTARTCalc_H

#include "TNamed.h"

class AliSTARTCalc: public TNamed {

 public:
   AliSTARTCalc();
   AliSTARTCalc(const char* name);
   AliSTARTCalc(const AliSTARTCalc &calibdata);
   AliSTARTCalc& operator= (const AliSTARTCalc &calibdata);
   virtual ~AliSTARTCalc();
   void Reset();
   void Print(const Option_t* option="") const;
   Float_t GetDelay(int channel) {return fTime[channel];}
   
   void SetTime(Float_t* daqtime, Float_t* dcstime);
  
 protected:
//   TMap fGraphs;
   Float_t fTime[24];

 ClassDef(AliSTARTCalc,4)
};

#endif
