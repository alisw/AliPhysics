//
// Class AliRsnComparisonAT
//
// TODO
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#ifndef ALIRSNMVCOMPARISONAT_H
#define ALIRSNMVCOMPARISONAT_H

#include "TObjArray.h"

#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnBaseAT.h"
#include "AliRsnCut.h"
#include "AliRsnCutSet.h"

class AliRsnComparisonAT : public AliRsnBaseAT
{
  public:

    enum
    {
      kMyInputNum=1,
      kMyPIDInputNum=5
    };

    AliRsnComparisonAT(const char*name="AliRsnComparisonAT");

    virtual ~AliRsnComparisonAT() {}

    virtual void   InitIOVars();
    virtual void   LocalInit();
    virtual void   CreateOutputObjects();
    virtual void   Exec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void        LoopOverESDtracks();
    void        LoopOverMCtracks();
    void        LoopOverRSNDaughters();
    void        PrintStat();

    void AddMyInput(AliRsnComparisonObj *obj,const Int_t &index=0);
    void AddMyPIDInput(AliRsnComparisonObj *obj,const Int_t &index=0);
    void SetMyPIDInputName(TString name="default",const Int_t &index=0);

    void SetNumberOfPIDInputs(const Int_t& theValue) { fMyPIDInputNum = theValue; }


  private:
  
    AliRsnComparisonAT(const AliRsnComparisonAT&)
        : AliRsnBaseAT(""),fOutList(0x0),fMyInputNum(1),fMyPIDInputNum(0) {}
    AliRsnComparisonAT& operator=(const AliRsnComparisonAT&) {return *this;}

    TList       *fOutList;                      // output list

    Int_t       fMyInputNum;
    TObjArray   fMyInput[kMyInputNum];
    Int_t       fMyPIDInputNum;
    TObjArray   fMyPIDInput[kMyPIDInputNum];



    ClassDef(AliRsnComparisonAT, 1)
};

#endif
