#ifndef ALILIGHTNEVENT_H
#define ALILIGHTNEVENT_H

/*
 * AliLightNEvent.h
 *
 *  Created on: 22 Nov 2017
 *      Author: bernhardhohlweger
 */


#include "TList.h"
#include "AliEventCuts.h"
#include "AliAODEvent.h"
#include "AliAnalysisUtils.h"
#include "Rtypes.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
class AliLightNEvent {
public:
    AliLightNEvent();
    AliLightNEvent(bool mvPileUp,bool EvtCutQA);
    virtual ~AliLightNEvent();
    void SetEvent(AliAODEvent *evt);
    TList *GetEvtCutList() const {return fEvtCutList;};
    void SetXVertex(double xvtx){fxVtx=xvtx;};
    double GetXVertex() const {return fxVtx;};
    void SetYVertex(double yvtx){fyVtx=yvtx;};
    double GetYVertex() const {return fyVtx;};
    void SetZVertex(double zvtx){fzVtx=zvtx;};
    double GetZVertex() const {return fzVtx;};
    void SetSPDMult(int spdMult){fSPDMult=spdMult;};
    int GetSPDMult() const {return fSPDMult;};
    void SetRefMult08(int refMult){fRefMult08=refMult;};
    int GetRefMult08() const {return fRefMult08;};
    void SetV0AMult(int v0AMult){fV0AMult=v0AMult;};
    int GetV0AMult() const {return fV0AMult;};
    void SetV0CMult(int v0CMult){fV0CMult=v0CMult;};
    int GetV0CMult() const {return fV0CMult;};
    void SetNumberOfContributers(int nContrib){fnContrib=nContrib;};
    int GetNumberOfContributers() const {return fnContrib;};
    void SetPassAliEvtSelection(bool pass){fPassAliEvtSelection=pass;};
    bool PassAliEvtSelection() const {return fPassAliEvtSelection;};
    void SetIsPileUp(bool PileUp) {fisPileUp=PileUp;};
    bool GetIsPileUp() const {return fisPileUp;};
    void SetHasVertex(bool HasVtx){fHasVertex=HasVtx;};
    bool GetHasVertex() const {return fHasVertex;};
    void SetMagneticField(bool HasField){fHasMagField=HasField;};
    bool GetMagneticField() const {return fHasMagField;};
    void SetSelectionStatus(bool pass){fisSelected=pass;};
    bool GetSelectionStatus() const {return fisSelected;};
private:
    int CalculateITSMultiplicity(AliAODEvent *evt);
    AliAnalysisUtils *fUtils;   //!
    AliEventCuts *fEvtCuts;     //!
    TList *fEvtCutList;         //!
    double fxVtx;               //!
    double fyVtx;               //!
    double fzVtx;               //!
    int fSPDMult;               //!
    int fRefMult08;             //!
    int fV0AMult;               //!
    int fV0CMult;               //!
    int fnContrib;              //!
    bool fPassAliEvtSelection;  //!
    bool fisPileUp;             //!
    bool fHasVertex;            //!
    bool fHasMagField;          //!
    bool fisSelected;           //!
    ClassDef(AliLightNEvent,1)
};


#endif /* ALILIGHTNEVENT_H */
