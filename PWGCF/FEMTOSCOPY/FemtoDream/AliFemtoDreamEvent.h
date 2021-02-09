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
#include "AliESDEvent.h"
#include "AliEventCuts.h"
#include "AliESDEvent.h"
#include "Rtypes.h"
#include "TList.h"
class AliFemtoDreamEvent {
 public:
  enum MultEstimator {
    kSPD = 0,
    kRef08 = 1,
    kV0M = 2,
    kV0A = 3,
    kV0C = 4,
  };
  AliFemtoDreamEvent();
  AliFemtoDreamEvent(bool mvPileUp, bool EvtCutQA, UInt_t trigger, bool useEvtCuts = true, float LowPtSpherCalc = 0.5);
  AliFemtoDreamEvent &operator=(const AliFemtoDreamEvent &obj);
  virtual ~AliFemtoDreamEvent();
  void SetEvent(AliAODEvent *evt);
  void SetEvent(AliVEvent *evt);
  void SetEvent(AliESDEvent *evt);
  TList *GetEvtCutList() const {
    return fEvtCutList;
  }
  ;
  void SetXVertex(float xvtx) {
    fxVtx = xvtx;
  }
  ;
  float GetXVertex() const {
    return fxVtx;
  }
  ;
  void SetYVertex(float yvtx) {
    fyVtx = yvtx;
  }
  ;
  float GetYVertex() const {
    return fyVtx;
  }
  ;
  void SetZVertex(float zvtx) {
    fzVtx = zvtx;
  }
  ;
  float GetZVertex() const {
    return fzVtx;
  }
  ;
  float GetZVertexSPD() const {
    return fzVtxSPD;
  }
  ;
  float GetBField() const {
    return fBField;
  }
  ;
  float GetZVertexTracks() const {
    return fzVtxTracks;
  }
  ;
  void SetSPDMult(int spdMult) {
    fSPDMult = spdMult;
  }
  ;
  int GetSPDMult() const {
    return fSPDMult;
  }
  ;
  void SetSPDClusterLy1(int spdCluster) {
    fNSPDClusterLy0 = spdCluster;
  }
  ;
  void SetSPDClusterLy2(int spdCluster) {
    fNSPDClusterLy1 = spdCluster;
  }
  ;
  int GetSPDCluster(int ly) const {
    if (ly == 0) {
      return fNSPDClusterLy0;
    } else if (ly == 1) {
      return fNSPDClusterLy1;
    } else {
      return 0.;
    }
  }
  ;
  void SetRefMult08(int refMult) {
    fRefMult08 = refMult;
  }
  ;
  int GetRefMult08() const {
    return fRefMult08;
  }
  ;
  void SetV0AMult(int v0AMult) {
    fV0AMult = v0AMult;
  }
  ;
  int GetV0AMult() const {
    return fV0AMult;
  }
  ;
  void SetV0CMult(int v0CMult) {
    fV0CMult = v0CMult;
  }
  ;
  int GetV0CMult() const {
    return fV0CMult;
  }
  ;
  int GetV0MMult() const {
    return (fV0AMult + fV0CMult) / 2.;
  }
  ;
  float GetV0ATime() const { return fV0ATime; }
  float GetV0CTime() const { return fV0CTime; }

  float GetV0MCentrality() const {
    return fV0MCentrality;
  }
  ;
  void SetNumberOfContributers(int nContrib) {
    fnContrib = nContrib;
  }
  ;
  int GetNumberOfContributers() const {
    return fnContrib;
  }
  ;
  void SetPassAliEvtSelection(bool pass) {
    fPassAliEvtSelection = pass;
  }
  ;
  void SetMultiplicityEstimator(AliFemtoDreamEvent::MultEstimator est) {
    fEstimator = est;
  }
  bool PassAliEvtSelection() const {
    return fPassAliEvtSelection;
  }
  ;
  void SetIsPileUp(bool PileUp) {
    fisPileUp = PileUp;
  }
  ;
  bool GetIsPileUp() const {
    return fisPileUp;
  }
  ;
  void SetHasVertex(bool HasVtx) {
    fHasVertex = HasVtx;
  }
  ;
  bool GetHasVertex() const {
    return fHasVertex;
  }
  ;
  void SetMagneticField(bool HasField) {
    fHasMagField = HasField;
  }
  ;
  bool GetMagneticField() const {
    return fHasMagField;
  }
  ;
  void SetSelectionStatus(bool pass) {
    fisSelected = pass;
  }
  ;
  bool GetSelectionStatus() const {
    return fisSelected;
  }
  ;
  int GetMultiplicity();
  TString ClassName() {
    return "AliFemtoDreamEvent";
  }
  ;
  void SetfLowPtSpherCalc(float LowPtSpherCalc) {
    fLowPtSpherCalc = LowPtSpherCalc;
  }
  ;
  void SetSpher(double spher) {
    fspher = spher;
  }
  ;
  float GetSpher() const {
    return fspher;
  }
  ;
  void SetSphero(double sphero) {
    fsphero = sphero;
  }
  float GetSphero() const {
    return fsphero;
  }
  void SetCalcSpherocity(bool calcsphero) {
      fcalcsphero=calcsphero;
  }

 private:
  AliFemtoDreamEvent(const AliFemtoDreamEvent&);
  int CalculateITSMultiplicity(AliAODEvent *evt);
  double CalculateSphericityEvent(AliAODEvent *evt, float lowerPtbound=0.5);
  double CalculateSphericityEvent(AliVEvent *evt);
  double CalculateSpherocityEvent(AliAODEvent *evt);
  double CalculateSpherocityEvent(AliVEvent *evt);
  AliAnalysisUtils *fUtils;   //!
  AliEventCuts *fEvtCuts;     //!
  bool fuseAliEvtCuts;        //!
  TList *fEvtCutList;         //!
  float fxVtx;                //!
  float fyVtx;                //!
  float fzVtx;                //!
  float fzVtxTracks;          //!
  float fzVtxSPD;             //!
  float fBField;              //!
  int fSPDMult;               //!
  int fNSPDClusterLy0;        //!
  int fNSPDClusterLy1;        //!
  int fRefMult08;             //!
  int fV0AMult;               //!
  int fV0CMult;               //!
  float fV0ATime;             //!
  float fV0CTime;             //!
  float fV0MCentrality;       //!
  int fnContrib;              //!
  bool fPassAliEvtSelection;  //!
  bool fisPileUp;             //!
  bool fHasVertex;            //!
  bool fHasMagField;          //!
  bool fisSelected;           //!
  MultEstimator fEstimator;   //!
  double fspher;            //!
  float fLowPtSpherCalc;      //!
  double fsphero;            //!
  bool fcalcsphero;         //!
ClassDef(AliFemtoDreamEvent,7)
};

#endif /* ALIFEMTODREAMEVENT_H_ */
