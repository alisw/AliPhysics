// $Id$
// Main author: Davide Caffarri 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveHF_H
#define AliEveHF_H

#include "AliAODRecoDecay.h"
#include "AliEveTrack.h"

#include <TEveVSDStructs.h>
#include <TPolyLine3D.h>

#include <TPDGCode.h>


class TH1F;
class TH2F;


class AliEveHFList;

class AliEveHF : public TEvePointSet
{
  friend class AliEveHFList;
  friend class AliEveHFEditor;

public:
  AliEveHF();
  AliEveHF(TEveRecTrack* tNeg, TEveRecTrack* tPos, Double_t primVtx[3], AliAODRecoDecay* aodObj, Double_t firstPointTrack[3], TEveTrackPropagator* rs);
  virtual ~AliEveHF();

  void MakeHF();

  virtual void SetMainColor(Color_t col)
  {
    TEvePointSet::SetMainColor(col);
    fPointingLine->SetLineColor(fMarkerColor);
  }

  void SetRnrStyle(TEveTrackPropagator* rs) { fRnrStyle = rs; }

  //HF Property

  Float_t GetPhi()     { return fRecDecayHF.Phi(); }
  Float_t GetEta()     { return fRecDecayP_HF.Eta(); }
  Float_t GetRadius()  { return fRecDecayHF.Perp(); }
  Float_t GetPt()      { return fRecDecayP_HF.Perp(); }

  Double_t GetInvariantMassPart()     { CalculateInvMass(fDecay); return fInvariantMassPart; }
  Double_t GetInvariantMassAntiPart() { CalculateInvMass(fDecay); return fInvariantMassAntiPart; }

  Float_t GetChi2Vtx()  { return fChi2SecondVtx; }
  Float_t GetCosPointingAngle() {return fPointingAngleHF; }

  AliAODRecoDecay *GetAODobj() {return fAODobj; }

  Int_t GetAODIndex() const { return fAODIndex; }
  void  SetAODIndex(Int_t ind) { fAODIndex = ind;}

  virtual const Text_t* GetName() const    { return Form("AOD_HF_%i",fAODIndex); }
  virtual const Text_t* GetTitle() const   { return Form("AOD_HF_%i",fAODIndex); }

  //Prongs Property

  Double_t GetProngDCA(Int_t iProng) { return fProngDCA[iProng]; }
  void SetProngDCA() const ;

  Double_t Getd0Prong(Int_t iProng)  { return fProngd0[iProng]; }
  void Setd0Prong() const ;

  void CalculateInvMass(Int_t decay);

  Bool_t SelectInvMass(Int_t decay, Float_t decayCuts);

  void      SetMaxProbPdgPid();
  Int_t     GetPdgProngMaxProb(Int_t iProng) { return fProngMaxProbPdg[iProng]; }
  Double_t  GetPidProngMaxProb(Int_t iProng) { return fProngMaxProbPid[iProng]; }

  TEveTrackPropagator* GetPropagator() const  { return fRnrStyle; }

  TEveTrack* GetNegTrack() { return fNegTrack; }
  TEveTrack* GetPosTrack() { return fPosTrack; }

  TEveLine*  GetPointingLine() { return fPointingLine; }

protected:

  AliAODRecoDecay  *fAODobj;

  TEveVector  fRecBirthHF;    // Reconstucted birth point of neutral particle
  TEveVector  fRecDecayHF;    // Point of closest approach
  TEveVector  fRecDecayP_HF;
  Double_t    fPointingAngleHF; // Track Pointing Angle

  TEveTrack        *fNegTrack;
  TEveTrack        *fPosTrack;

  TEveTrackPropagator *fRnrStyle;

  TEveLine         *fPointingLine;

  Int_t             fnProng;      // Number of Prong.
  Int_t             fAODIndex;    // Index in HF loop array.
  Double_t          fChi2SecondVtx;      //Secondary Vertex Chi-square.

  Double_t          *fProngDCA;//[fnProng] Distance at the point of closest approach.
  Double_t          *fProngd0;//[fnProng] Impact Paramter if each prong.
  Int_t             *fProngMaxProbPdg;//[fnProng] Maximum PDG probability for the negative daughter
  Double_t          *fProngMaxProbPid;//[fnProng] Maximum PID probability for the negative daughter

  Double_t           fInvariantMassPart;
  Double_t           fInvariantMassAntiPart;

  Int_t              fDecay;

 private:
  AliEveHF(const AliEveHF&);            // Not implemented
  AliEveHF& operator=(const AliEveHF&); // Not implemented

  ClassDef(AliEveHF,0); // Visual representation of a AliEveHF.
};


/******************************************************************************/
// AliEveHFList
/******************************************************************************/
class AliEveHFList : public TEveElementList
{
  friend class AliEveHFListEditor;
  //  friend class AliEveHF;

public:
  AliEveHFList();
  AliEveHFList(TEveTrackPropagator* rs);
  AliEveHFList(const Text_t* name, TEveTrackPropagator* rs=0);
  virtual  ~AliEveHFList() {}

  virtual const Text_t* GetTitle() const { return fTitle; }
  virtual void SetTitle(const Text_t* t) { fTitle = t; }
  virtual void SetTracksColor(Color_t cNeg, Color_t cPos, Int_t ip)
  {fProngColor[ip] = cNeg; fProngColor[ip++] = cPos;}

  virtual Bool_t CanEditMainColor() const { return kTRUE; }

  void  SetRnrStyle(TEveTrackPropagator* rst) { fRnrStyle = rst; }
  TEveTrackPropagator* GetPropagator()        { return fRnrStyle; }

  Bool_t GetRnrHFvtx()     const { return fRnrHFvtx; }
  Bool_t GetRnrHFpath()    const { return fRnrHFpath; }
  Bool_t GetRnrDaughters() const { return fRnrDaughters; }

  void   MakeHFs();

  void   FilterByPt(Float_t minPt, Float_t maxPt);
  void   FilterByRadius(Float_t minR, Float_t maxR);
  void   FilterByCosPointingAngle(Float_t minCosPointingAngle, Float_t maxCosPointingAngle);
  void   FilterByDCA(Float_t minDaughterDCA, Float_t maxDaughterDCA);
  void   FilterByd0(Float_t mind0, Float_t maxd0);
  //void   FilterByCheckedPidMinProb(Int_t rFlag, Int_t rDaughter, Int_t rPid, Float_t rProb);

  // void   SetProngCheckedPid(Int_t rProngCheckedPid) const;
  Int_t  GetProngCheckedPid(Int_t iProng) {return fProngCheckedPid[iProng];}
  // void   SetProngCheckedProb(Float_t rProngCheckedProb) const;
  // Float_t  GetProngCheckedProb(Int_t iProng) const  { return fProngCheckedProb[iProng]; }

  void   FilterByInvariantMass (Int_t decay, Float_t deltaInvariantMass);

protected:
  TString              fTitle;

  TEveTrackPropagator *fRnrStyle;

  Bool_t               fRnrDaughters;
  Bool_t               fRnrHFvtx;
  Bool_t               fRnrHFpath;

  Color_t*             fProngColor;//[fnProng]

  Float_t              fMinRCut;
  Float_t              fMaxRCut;

  Float_t              fMinDaughterDCA;
  Float_t              fMaxDaughterDCA;

  Float_t              fMinPt;
  Float_t              fMaxPt;

  Float_t              fMinCosPointingAngle;
  Float_t              fMaxCosPointingAngle;

  Float_t              fMind0;
  Float_t              fMaxd0;

  Int_t*               fProngCheckedPid;//[fnProng]

  Float_t*             fProngCheckedProb;//[fnProng]

  Float_t              fDeltaInvariantMass;
  Int_t                fDecay;

private:
  void Init();

  AliEveHFList(const AliEveHFList&);            // Not implemented
  AliEveHFList& operator=(const AliEveHFList&); // Not implemented

  ClassDef(AliEveHFList,0); // A list of AliEveHF objecs.
};


#endif
