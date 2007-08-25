#ifndef ALIITSNOISESSD_H
#define ALIITSNOISESSD_H
 
#include "TObjArray.h"
#include "TArrayF.h"

//////////////////////////////////////////////
// Author: Enrico Fragiacomo
// Date: 23/08/2007
//                                          //
//////////////////////////////////////////////
class AliITSNoiseSSD : public TObject {

 public:

    AliITSNoiseSSD();
    virtual ~AliITSNoiseSSD();
    AliITSNoiseSSD(const AliITSNoiseSSD &source); // copy constructor
    AliITSNoiseSSD& operator=(const AliITSNoiseSSD &source); // ass. op.

    void SetNNoiseP(Int_t n) { fNoisP.Set(n); }
    void AddNoiseP(Int_t c, Float_t n) { fNoisP.AddAt(n,c);}       
    TArrayF GetNoiseP() const {return fNoisP; }
    Float_t GetNoiseP(Int_t n) {return fNoisP.At(n); }
    void SetNNoiseN(Int_t n) { fNoisN.Set(n); }
    void AddNoiseN(Int_t c, Float_t n) { fNoisN.AddAt(n,c);}
    TArrayF GetNoiseN() const {return fNoisN; }
    Float_t GetNoiseN(Int_t n) {return fNoisN.At(n); }

    void SetMod(UShort_t mod) {fMod = mod;}
    UShort_t GetMod() { return fMod;}

protected:

  UShort_t fMod;       // module number (from 0 to 1535). Needed by the preprocessor to 
                       //   extract the information from the Detector Algorithm
  
  TArrayF fNoisP;           // Noise for P side channels
  TArrayF fNoisN;           // Noise for N side channels
  
 private:

    ClassDef(AliITSNoiseSSD,1) // Noise  class for SSD
};
#endif
