#ifndef ALIITSPEDESTALSSD_H
#define ALIITSPEDESTALSSD_H
 
#include "TObjArray.h"
#include "TArrayF.h"

//////////////////////////////////////////////
// Author: Enrico Fragiacomo
// Date: 12/12/2007
//                                          //
//////////////////////////////////////////////

/* $Id$ */

class AliITSPedestalSSD : public TObject {

 public:

    AliITSPedestalSSD();
    virtual ~AliITSPedestalSSD();
    AliITSPedestalSSD(const AliITSPedestalSSD &source); // copy constructor
    AliITSPedestalSSD& operator=(const AliITSPedestalSSD &source); // ass. op.

    void SetNPedestalP(Int_t n) { fPedP.Set(n); }
    void AddPedestalP(Int_t c, Float_t n) { fPedP.AddAt(n,c);}       
    TArrayF GetPedestalP() const {return fPedP; }
    Float_t GetPedestalP(Int_t n) {return fPedP.At(n); }
    void SetNPedestalN(Int_t n) { fPedN.Set(n); }
    void AddPedestalN(Int_t c, Float_t n) { fPedN.AddAt(n,c);}
    TArrayF GetPedestalN() const {return fPedN; }
    Float_t GetPedestalN(Int_t n) {return fPedN.At(n); }

    void SetMod(UShort_t mod) {fMod = mod;}
    UShort_t GetMod() { return fMod;}

protected:

  UShort_t fMod;       // module number (from 0 to 1535). Needed by the preprocessor to 
                       //   extract the information from the Detector Algorithm
  
  TArrayF fPedP;           // Pedestal for P side channels
  TArrayF fPedN;           // Pedestal for N side channels
  
 private:

    ClassDef(AliITSPedestalSSD,1) // Pedestal  class for SSD
};
#endif
