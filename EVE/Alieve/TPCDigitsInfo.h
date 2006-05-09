// $Header$

#ifndef ALIEVE_TPCDigitsInfo_H
#define ALIEVE_TPCDigitsInfo_H

#include <Reve/VSD.h>

#include <vector>

#include <TNamed.h>
#include <TArrayI.h>
#include <TTree.h>

#include <AliTPCParam.h>
#include <AliSimDigits.h>


namespace Alieve {

  class TPCSeg {
  public:   
    Float_t   fPadWidth, fPadLength,fRlow; // vertices data
    Int_t     fNRows;                      // text & vertices data 
    Int_t     fNMaxPads;                   // texture offset data
    Float_t   fStepY[64];                  // y coord wher npads has changed
    Int_t     fNsteps;                     // number of steps

    void Dump() const;
  };

  class TPCDigitsInfo : public TNamed
  {
  private:
    void Init();

  protected:
    Int_t               fRefCount;
    TString             fDataDir;  
    Int_t               fEvent;    

  public:
    AliSimDigits        fSimDigits;
    AliTPCParam*        fParameter;
    TTree*              fTree;
    std::vector<Int_t>  fSegEnt;
    TPCSeg              fInnSeg;   
    TPCSeg              fOut1Seg;  
    TPCSeg              fOut2Seg;  

    TPCDigitsInfo(const Text_t* n="TPCDigitsInfo", const Text_t* t=0) :
      TNamed(n, t) { Init(); }
    virtual ~TPCDigitsInfo();

    void SetData(AliTPCParam* par, TTree* digits);
   
    void IncRefCount() { ++fRefCount; }
    void DecRefCount() { --fRefCount; if(fRefCount <= 0) delete this; }

    virtual void Print(Option_t* opt="") const;

    ClassDef(TPCDigitsInfo, 1);
  }; // endclass TPCDigitsInfo
}
#endif
