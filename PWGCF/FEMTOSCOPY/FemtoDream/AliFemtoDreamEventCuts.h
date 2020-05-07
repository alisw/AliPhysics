/*
 * AliFemtoEventCuts.h
 *
 *  Created on: Nov 22, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMEVENTCUTS_H_
#define ALIFEMTODREAMEVENTCUTS_H_
#include "Rtypes.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamEventHist.h"
class AliFemtoDreamEventCuts {
 public:
  AliFemtoDreamEventCuts();
  AliFemtoDreamEventCuts(const AliFemtoDreamEventCuts& cuts);
  AliFemtoDreamEventCuts& operator=(const AliFemtoDreamEventCuts& cuts);
  virtual ~AliFemtoDreamEventCuts();
  void SetMinimalBooking(bool doIt) {
    fMinimalBooking = doIt;
  }
  ;
  bool GetMinimalBooking() {
    return fMinimalBooking;
  }
  ;
  bool isSelected(AliFemtoDreamEvent *evt);
  static AliFemtoDreamEventCuts* StandardCutsRun1();
  static AliFemtoDreamEventCuts* StandardCutsRun2();
  void SetCutMinContrib(int nMinContrib) {
    fCutMinContrib = true;
    fMinContrib = nMinContrib;
  }
  ;
  void SetZVtxPosition(float zVtxLow, float zVtxUp) {
    fzVtxLow = zVtxLow;
    fzVtxUp = zVtxUp;
    fCutZVtx = true;
  }
  ;
  void SetMVPileUpRejection(bool apply) {
    fUseMVPileUpRej = apply;
  }
  ;
  bool GetMVPileUpRejection() const {
    return fUseMVPileUpRej;
  }
  ;
  void PileUpRejection(bool apply) {
    fPileUpRejection = apply;
  }
  ;
  void CleanUpMult(bool SPD, bool v0A, bool v0C, bool RefMult) {
    fUseSPDMult = SPD;
    fUseV0AMult = v0A;
    fUseV0CMult = v0C;
    fUseRef08Mult = RefMult;
    fCleanEvtMult = true;
  }
  //Everything else disabled if you use the following option:
  void UseDontWorryEvtCuts(bool apply) {
    fUseAliEvtCuts = apply;
  }
  ;
  bool GetUseAliEventCuts () const {
    return fUseAliEvtCuts;
  }
  void SetMultVsCentPlots(bool doIt) {
    fCentVsMultPlots = doIt;
  }
  ;
  void InitQA();
  TList *GetHistList() const {
    return fHist->GetHistList();
  }
  ;
  void SetDoSphericityCuts(bool doIt) {
    fDoSpherCuts = doIt;
  }
  ;
  bool GetDoSphericityCuts() const {
    return fDoSpherCuts;
  }
  ;
  void SetSphericityCuts(float spherlow, float spherup, float ptLow) {
    fSpherCutsLow = spherlow;
    fSpherCutsUp = spherup;
    fLowPtSpherCalc = ptLow;
    fDoSpherCuts = true;
  }
  ;
  float GetlowerPtBoundSpherCalc() const {
    return fLowPtSpherCalc;
  }
  ;

  void SetMultiplicityPercentileMax(float multMax) {
    fUseMultPercentileCut = true;
    fMultPercentileMax = multMax;
  }
  float GetSphericityCutsLow() const {
    return fSpherCutsLow;
  }
  ;
  float GetSphericityCutsUp() const {
    return fSpherCutsUp;
  }
  ;
  void SetDoSpherocityCuts(bool doIt) {
    fDoSpheroCuts = doIt;
  }
  ;
  bool GetDoSpherocityCuts() const {
    return fDoSpheroCuts;
  }
  ;
  void SetSpherocityCuts(float spherolow, float spheroup) {
    fSpheroCutsLow = spherolow;
    fSpheroCutsUp = spheroup;
    fDoSpheroCuts = true;
  }
  ;
  float GetSpherocityCutsLow() const {
    return fSpheroCutsLow;
  }
  ;
  float GetSpherocityCutsUp() const {
    return fSpheroCutsUp;
  }
  ;

 private:
  void BookQA(AliFemtoDreamEvent *evt);
  void BookCuts();
  AliFemtoDreamEventHist *fHist;  //!
  bool fMinimalBooking;           //
  bool fCutMinContrib;            //
  int fMinContrib;                //
  bool fCutZVtx;                  //
  float fzVtxLow;                //
  float fzVtxUp;                 //
  bool fPileUpRejection;          //
  bool fUseMVPileUpRej;           // Which method of Pile Up Rej should be used
  bool fCleanEvtMult;             //
  bool fUseSPDMult;               //
  bool fUseV0AMult;               //
  bool fUseV0CMult;               //
  bool fUseRef08Mult;             //
  bool fUseMultPercentileCut;     //
  float fMultPercentileMax;	      //
  //Use evt cuts tuned by expert(don't worry solution)
  bool fUseAliEvtCuts;            //
  bool fCentVsMultPlots;    //
  bool fDoSpherCuts;       //
  float  fLowPtSpherCalc;  //  
  float fSpherCutsLow;      //
  float fSpherCutsUp;      //
  bool fDoSpheroCuts;       //
  float fSpheroCutsLow;      //
  float fSpheroCutsUp;      //
ClassDef(AliFemtoDreamEventCuts,8)
};

#endif /* ALIFEMTODREAMEVENTCUTS_H_ */
