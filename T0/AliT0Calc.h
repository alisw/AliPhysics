#ifndef AliT0Calc_H
#define AliT0Calc_H

#include "TNamed.h"

class AliT0Calc: public TNamed {

 public:
   AliT0Calc();
   AliT0Calc(const char* name);
   AliT0Calc(const AliT0Calc &calibdata);
   AliT0Calc& operator= (const AliT0Calc &calibdata);
   virtual ~AliT0Calc();
   void Reset();
   void Print(const Option_t* option="") const;
   Float_t GetDelay(int channel) {return fTime[channel];}
   
   void SetTime(Float_t* daqtime, Float_t* dcstime);
  
 protected:
//   TMap fGraphs;
   Float_t fTime[24];

 ClassDef(AliT0Calc,4)
};

typedef AliT0Calc AliSTARTCalc; // for backward compatibility

#endif
