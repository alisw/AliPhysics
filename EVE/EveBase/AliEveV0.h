// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveV0_H
#define AliEveV0_H


//------------------------------------------------------------------------------
// AliEveV0
//------------------------------------------------------------------------------
//
// Representation of a reconstructed V0.
//
//
//------------------------------------------------------------------------------
// AliEveV0List
//------------------------------------------------------------------------------
//
// Container for AliEveV0s.
//
// Allows runtime selection by pT and DCA of daughters, radius of V0
// creation and PID priobabilities

//==============================================================================

#include "AliEveTrack.h"
#include <TEveVSDStructs.h>

#include <TPolyLine3D.h>

#include <TPDGCode.h>

class TH1F;
class TH2F;


class AliEveV0List;

class AliEveV0 : public TEvePointSet
{
  friend class AliEveV0List;
  friend class AliEveV0Editor;

public:
  AliEveV0();
  AliEveV0(TEveRecTrack* tNeg, TEveRecTrack* tPos, TEveRecV0* v0,
     TEveTrackPropagator* rsNeg,TEveTrackPropagator* rsPos);
  virtual ~AliEveV0();

  void MakeV0();

  virtual void  SetMainColor(Color_t col)
  {
    TEvePointSet::SetMainColor(col);
    fPointingLine->SetLineColor(fMarkerColor);
  }

  void SetRnrStyleNeg(TEveTrackPropagator* rs) { fRnrStyleNeg = rs; }
  void SetRnrStylePos(TEveTrackPropagator* rs) { fRnrStylePos = rs; }

  Float_t GetDaughterDCA() const { return fDaughterDCA; }
  void SetDaughterDCA(Float_t dca) { fDaughterDCA = dca; }

  Float_t GetPhi()    const { return fRecDecayP.Phi(); }
  Float_t GetEta()    const { return fRecDecayP.Eta(); }
  Float_t GetRadius() const { return fRecDecayV.Perp(); }
  Float_t GetPt()     const { return fRecDecayP.Perp(); }
  Float_t GetMomentum() const { return fRecDecayP.Mag(); }

  Float_t GetInvMass(Int_t nPdgCode, Int_t pPdgCode) const;
  Float_t GetK0sInvMass() const { return GetInvMass(kPiMinus,kPiPlus); }
  Float_t GetLambdaInvMass() const { return GetInvMass(kPiMinus,kProton); }
  Float_t GetAntiLambdaInvMass() const { return GetInvMass(kProton,kPiPlus); }

  Bool_t GetOnFlyStatus()    const { return fOnFlyStatus; }
  void   SetOnFlyStatus(Bool_t fs) { fOnFlyStatus = fs; }

  void    SetMaxProbPdgPid(Int_t iDaughter, Int_t rPdg, Float_t rPid);
  Int_t   GetNegMaxProbPdg() const { return fNegMaxProbPdg; }
  Int_t   GetPosMaxProbPdg() const { return fPosMaxProbPdg; }
  Float_t GetNegMaxProbPid() const { return fNegMaxProbPid; }
  Float_t GetPosMaxProbPid() const { return fPosMaxProbPid; }

  Int_t GetESDIndex() const { return fESDIndex; }
  void  SetESDIndex(Int_t ind) { fESDIndex = ind;}

  virtual const Text_t* GetName()  const { return Form("ESDv0_%i",fESDIndex); }
  virtual const Text_t* GetTitle() const { return Form("ESDv0_%i",fESDIndex); }

  TEveTrackPropagator* GetPropagatorNeg() const  { return fRnrStyleNeg; }
  TEveTrackPropagator* GetPropagatorPos() const  { return fRnrStylePos; }

  AliEveTrack* GetNegTrack() const { return fNegTrack; }
  AliEveTrack* GetPosTrack() const { return fPosTrack; }

  TEveLine*  GetPointingLine() const { return fPointingLine; }

protected:
  TEveVector fRecBirthV;    // Reconstucted birth point of neutral particle
  TEveVector fRecDecayV;    // Point of closest approach
  TEveVector fRecDecayP;    // Reconstructed momentum of decayed particle.

  AliEveTrack        *fNegTrack; // Representation of negative track.
  AliEveTrack        *fPosTrack; // Representation of positive track.

  TEveTrackPropagator *fRnrStyleNeg; // Track propagator for neg track.
  TEveTrackPropagator *fRnrStylePos; // Track propagator for pos track.

  TEveLine         *fPointingLine; // Representation of pointing line.

  Int_t             fESDIndex;    // Index in ESD V0 array.
  Bool_t            fOnFlyStatus; // Reconstructed during tracking.
  Float_t           fDaughterDCA; // Distance at the point of closest approach. 
  Float_t           fChi2V0;      // Some Chi-square.

  Int_t             fNegMaxProbPdg; // Maximum PDG probability for the negative daughter
  Int_t             fPosMaxProbPdg; // Maximum PDG probability for the positive daughter
  Float_t           fNegMaxProbPid; // Maximum PID probability for the negative daughter
  Float_t           fPosMaxProbPid; // Maximum PID probability for the positive daughter

private:
  AliEveV0(const AliEveV0&);            // Not implemented
  AliEveV0& operator=(const AliEveV0&); // Not implemented

  ClassDef(AliEveV0, 0); // Visual representation of a AliEveV0.
};


/******************************************************************************/
// AliEveV0List
/******************************************************************************/

class AliEveV0List : public TEveElementList
{
  friend class AliEveV0ListEditor;

public:
  AliEveV0List();
  AliEveV0List(TEveTrackPropagator* rsNeg,TEveTrackPropagator* rsPos);
  AliEveV0List(const Text_t* name, TEveTrackPropagator* rsNeg=0,TEveTrackPropagator* rsPos=0);
  virtual ~AliEveV0List() {}

  virtual const Text_t* GetTitle() const { return fTitle; }
  virtual void SetTitle(const Text_t* t) { fTitle = t; }
  virtual void SetTracksColor(Color_t cNeg, Color_t cPos) {
    fNegColor = cNeg; fPosColor = cPos;}

  virtual Bool_t CanEditMainColor() const { return kTRUE; }

  void  SetRnrStyleNeg(TEveTrackPropagator* rst) { fRnrStyleNeg = rst; }
  void  SetRnrStylePos(TEveTrackPropagator* rst) { fRnrStylePos = rst; }
  TEveTrackPropagator* GetPropagatorNeg()        { return fRnrStyleNeg; }
  TEveTrackPropagator* GetPropagatorPos()        { return fRnrStylePos; }

  Bool_t GetRnrV0vtx()     const { return fRnrV0vtx; }
  Bool_t GetRnrV0path()    const { return fRnrV0path; }
  Bool_t GetRnrDaughters() const { return fRnrDaughters; }

  void   MakeV0s();

  void   FilterByRadius(Float_t minR, Float_t maxR);
  void   FilterByDaughterDCA(Float_t minDaughterDCA, Float_t maxDaughterDCA);
  void   FilterByPt(Float_t minPt, Float_t maxPt);
  void   FilterByCheckedPidMinProb(Int_t rFlag, Int_t rDaughter, Int_t rPid, Float_t rProb);
  void   SetNegCheckedPid(Int_t rNegCheckedPid) {fNegCheckedPid = rNegCheckedPid;}
  void   SetPosCheckedPid(Int_t rPosCheckedPid) {fPosCheckedPid = rPosCheckedPid;}
  Int_t  GetNegCheckedPid() const { return fNegCheckedPid; }
  Int_t  GetPosCheckedPid() const { return fPosCheckedPid; }

  void   SetNegCheckedProb(Float_t rNegCheckedProb) {fNegCheckedProb = rNegCheckedProb;}
  void   SetPosCheckedProb(Float_t rPosCheckedProb) {fPosCheckedProb = rPosCheckedProb;}
  Float_t  GetNegCheckedProb() const { return fNegCheckedProb; }
  Float_t  GetPosCheckedProb() const { return fPosCheckedProb; }

  void   FilterByInvariantMass(Float_t minPt, Float_t maxPt, Int_t nPdgCode, Int_t pPdgCode);

protected:
  TString              fTitle;    // Title of the object.

  TEveTrackPropagator *fRnrStyleNeg; // Track propagator to be passed do conteined V0s.
  TEveTrackPropagator *fRnrStylePos; // Track propagator to be passed do conteined V0s.
    
  Bool_t               fRnrDaughters; // Flag - display daughter tracks.
  Bool_t               fRnrV0vtx;     // Flag - display V0 vertex.
  Bool_t               fRnrV0path;    // Flag - display V0 path.

  Color_t              fNegColor;     // Color for negative tracks.
  Color_t              fPosColor;     // Color for positive tracks.

  Float_t              fMinRCut;      // Minimum radius cut.
  Float_t              fMaxRCut;      // Maximum radius cut.

  Float_t              fMinDaughterDCA; // Minimum daughter DCA cut.
  Float_t              fMaxDaughterDCA; // Maximum daughter DCA cut.

  Float_t              fMinPt;          // Minimum Pt cut.
  Float_t              fMaxPt;          // Maximum Pt cut.

  Int_t                fNegCheckedPid;  // BORIS ?
  Int_t                fPosCheckedPid;  // BORIS ?

  Float_t              fNegCheckedProb; // BORIS ?
  Float_t              fPosCheckedProb; // BORIS ?

  Float_t              fMinInvariantMass; // Minimum invariant mass cut.
  Float_t              fMaxInvariantMass; // Maximum invariant mass cut.

private:
  void Init();

  AliEveV0List(const AliEveV0List&);            // Not implemented
  AliEveV0List& operator=(const AliEveV0List&); // Not implemented

  ClassDef(AliEveV0List, 0); // A list of AliEveV0 objecs.
};


#endif
