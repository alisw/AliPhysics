#ifndef ALILIGHTNEVENTCUTS_H
#define ALILIGHTNEVENTCUTS_H

/*
 * AliFemtoEventCuts.h
 *
 *  Created on: Nov 22, 2017
 *      Author: gu74req
 */

#include "Rtypes.h"
#include "AliLightNEvent.h"
#include "AliLightNEventHist.h"
class AliLightNEventCuts {
public:
    AliLightNEventCuts();
    virtual ~AliLightNEventCuts();
    bool isSelected(AliLightNEvent *evt);
    static AliLightNEventCuts* StandardCutsRun1();
    static AliLightNEventCuts* StandardCutsRun2();
    void SetCutMinContrib(int nMinContrib) {
        fCutMinContrib=true;fMinContrib=nMinContrib;
    };
    void SetZVtxPosition(double zVtxLow,double zVtxUp) {
        fzVtxLow=zVtxLow;fzVtxUp=zVtxUp;fCutZVtx=true;
    };
    void SetMVPileUpRejection(bool apply){fUseMVPileUpRej=apply;};
    bool GetMVPileUpRejection() const {return fUseMVPileUpRej;};
    void PileUpRejection(bool apply){fPileUpRejection=apply;};
    void CleanUpMult(bool SPD,bool v0A, bool v0C, bool RefMult) {
        fUseSPDMult=SPD;fUseV0AMult=v0A;fUseV0CMult=v0C;
        fUseRef08Mult=RefMult;fCleanEvtMult=true;
    }
    //Everything else disabled if you use the following option:
    void UseDontWorryEvtCuts(bool apply) {fUseAliEvtCuts=apply;};
    void InitQA();
    //Histogram things
    TList *GetHistList() const {return fHist->GetHistList();};
    void SetName(TString name){fHist->SetName(name.Data());};
    void FillV0Mlpercentile(float val) {fHist->FillV0Mpercentile(val);};
    void FillV0MlpercentileHM(float val) {fHist->FillV0MpercentileHM(val);};
private:
    void BookQA(AliLightNEvent *evt);
    void BookCuts();
    AliLightNEventHist *fHist;  //!
    bool fCutMinContrib;            //
    int fMinContrib;                //
    bool fCutZVtx;                  //
    double fzVtxLow;                //
    double fzVtxUp;                 //
    bool fPileUpRejection;          //
    bool fUseMVPileUpRej;           // Which method of Pile Up Rej should be used
    bool fCleanEvtMult;             //
    bool fUseSPDMult;               //
    bool fUseV0AMult;               //
    bool fUseV0CMult;               //
    bool fUseRef08Mult;             //
    //Use evt cuts tuned by expert(don't worry solution)
    bool fUseAliEvtCuts;            //
    ClassDef(AliLightNEventCuts,1)
};

#endif /* ALILIGHTNEVENTCUTS_H */
