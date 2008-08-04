#ifndef ALIITSGAINSSD_H
#define ALIITSGAINSSD_H
 
#include "TObjArray.h"
#include "TArrayF.h"

//////////////////////////////////////////////
// Author: Enrico Fragiacomo
// Date: 23/08/2007
//                                          //
//////////////////////////////////////////////
class AliITSGainSSD : public TObject {

 public:

    AliITSGainSSD();
    virtual ~AliITSGainSSD();

    void SetNGainP(Int_t n) { fGainP.Set(n); }
    void AddGainP(Int_t c, Float_t n) { fGainP.AddAt(n,c);}       
    TArrayF GetGainP() const {return fGainP; }
    Float_t GetGainP(Int_t n) {return fGainP.At(n); }
    void SetNGainN(Int_t n) { fGainN.Set(n); }
    void AddGainN(Int_t c, Float_t n) { fGainN.AddAt(n,c);}
    TArrayF GetGainN() const {return fGainN; }
    Float_t GetGainN(Int_t n) {return fGainN.At(n); }

    void SetMod(UShort_t mod) {fMod = mod;}
    UShort_t GetMod() { return fMod;}

protected:

  UShort_t fMod;       // module number (from 0 to 1535). Needed by the preprocessor to 
                       //   extract the information from the Detector Algorithm
  
  TArrayF fGainP;           // Gain for P side channels
  TArrayF fGainN;           // Gain for N side channels
  
 private:
    AliITSGainSSD(const AliITSGainSSD &source); // copy constructor
    AliITSGainSSD& operator=(const AliITSGainSSD &source); // ass. op.

    ClassDef(AliITSGainSSD,1) //Response class for SSD
};
#endif
