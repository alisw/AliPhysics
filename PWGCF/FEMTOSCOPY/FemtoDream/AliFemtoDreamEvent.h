/*
 * AliFemtoDreamEvent.h
 *
 *  Created on: 22 Nov 2017
 *      Author: bernhardhohlweger
 */

#ifndef ALIFEMTODREAMEVENT_H_
#define ALIFEMTODREAMEVENT_H_
#include "AliAnalysisUtils.h"
#include "AliAODEvent.h"
#include "AliEventCuts.h"
#include "Rtypes.h"
#include "TList.h"
class AliFemtoDreamEvent {
 public:
  AliFemtoDreamEvent();
  AliFemtoDreamEvent(bool mvPileUp,bool EvtCutQA, UInt_t trigger);
  virtual ~AliFemtoDreamEvent();
  void SetEvent(AliAODEvent *evt);
  TList *GetEvtCutList() const {return fEvtCutList;};
  void SetXVertex(float xvtx){fxVtx=xvtx;};
  float GetXVertex() const {return fxVtx;};
  void SetYVertex(float yvtx){fyVtx=yvtx;};
  float GetYVertex() const {return fyVtx;};
  void SetZVertex(float zvtx){fzVtx=zvtx;};
  float GetZVertex() const {return fzVtx;};
  void SetSPDMult(int spdMult){fSPDMult=spdMult;};
  int GetSPDMult() const {return fSPDMult;};
  void SetRefMult08(int refMult){fRefMult08=refMult;};
  int GetRefMult08() const {return fRefMult08;};
  void SetV0AMult(int v0AMult){fV0AMult=v0AMult;};
  int GetV0AMult() const {return fV0AMult;};
  void SetV0CMult(int v0CMult){fV0CMult=v0CMult;};
  int GetV0CMult() const {return fV0CMult;};
  int GetV0MMult() const {return (fV0AMult+fV0CMult)/2.;};
  float GetV0MCentrality() const {return fV0MCentrality;};
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
  TString ClassName() {return "AliFemtoDreamEvent";};
 private:
  int CalculateITSMultiplicity(AliAODEvent *evt);
  AliAnalysisUtils *fUtils;   //!
  AliEventCuts *fEvtCuts;     //!
  TList *fEvtCutList;         //!
  float fxVtx;               //!
  float fyVtx;               //!
  float fzVtx;               //!
  int fSPDMult;               //!
  int fRefMult08;             //!
  int fV0AMult;               //!
  int fV0CMult;               //!
  float fV0MCentrality;       //!
  int fnContrib;              //!
  bool fPassAliEvtSelection;  //!
  bool fisPileUp;             //!
  bool fHasVertex;            //!
  bool fHasMagField;          //!
  bool fisSelected;           //!
  ClassDef(AliFemtoDreamEvent,2)
};

#endif /* ALIFEMTODREAMEVENT_H_ */
