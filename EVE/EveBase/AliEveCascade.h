// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef ALIEVECASCADE_H
#define ALIEVECASCADE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//----------------------------------------------------------------------------
// This code defines the reconstructed cascade (ESD level) visualized with EVE
//
// Origin : Boris Hippolyte, IPHC (hippolyt@in2p3.fr)
// Modified : Antonin Maire, April 2009, IPHC (antonin.maire@cern.ch)
//----------------------------------------------------------------------------

class TH1F;
class TH2F;
class TVector3;

class TEveVSDStructs;
class TEveTrackPropagator;
class TEveLine;

#include <TEveVector.h>
#include <TEveVSDStructs.h>
#include <TEvePointSet.h>
#include <TPDGCode.h>

 
class AliEveCascadeList;

#include "AliEveTrack.h"


class AliEveCascade : public TEvePointSet
{
  friend class AliEveCascadeList;
  friend class AliEveCascadeEditor;

public:
  AliEveCascade();
  AliEveCascade(TEveRecTrack*        tBac, 
		TEveRecTrack*        tNeg, 
		TEveRecTrack*        tPos, 
		TEveRecV0*           v0, 
		TEveRecCascade*      cascade, 
		TEveTrackPropagator* rsBac,
        TEveTrackPropagator* rsNeg,
        TEveTrackPropagator* rsPos);
  virtual ~AliEveCascade();

  void MakeCascade();

  virtual void  SetMainColor(Color_t col)
  {
    TEvePointSet::SetMainColor(col);
    fPointingCurve->SetLineColor(fMarkerColor);
    fV0Path->SetLineColor(fMarkerColor);
  }

  void 		SetRnrStyleBac( TEveTrackPropagator* const rs) { fRnrStyleBac = rs; }
  void 		SetRnrStyleNeg( TEveTrackPropagator* const rs) { fRnrStyleNeg = rs; }
  void 		SetRnrStylePos( TEveTrackPropagator* const rs) { fRnrStylePos = rs; }

  Float_t 	GetDaughterDCA() const { return fDaughterDCA; }
  void 		SetDaughterDCA(Float_t dca) { fDaughterDCA = dca; }

  Float_t 	GetRadius() const { return fRecDecayV.Perp(); }
  Float_t 	GetPt()     const { return fRecDecayP.Perp(); }
  Float_t 	GetPtot()   const { return fRecDecayP.Mag(); }
  
  Float_t 	GetPhi()    const { return fRecDecayP.Phi(); }
  Float_t 	GetTheta()  const { return fRecDecayP.Theta(); }
  Float_t 	GetEta()    const { return fRecDecayP.Eta(); }
  Int_t 	GetCharge() const { return fBacTrack->GetCharge(); }
   
  Double_t 	GetInvMass(Int_t cascadePdgCodeHyp) const;
  Float_t 	GetXiMinusInvMass()    const { return GetInvMass( kXiMinus); }
  Float_t 	GetOmegaMinusInvMass() const { return GetInvMass( kOmegaMinus); }
  Float_t 	GetXiPlusInvMass()     const { return GetInvMass(-kXiMinus); }
  Float_t 	GetOmegaPlusInvMass()  const { return GetInvMass(-kOmegaMinus); }
   

  Int_t 	GetESDIndex() const { return fESDIndex; }
  void  	SetESDIndex(Int_t ind) { fESDIndex = ind;}
  
  TVector3	GetLambdaP()  const { return fLambdaP; }
  void		SetLambdaP(Double_t px, Double_t py, Double_t pz) { fLambdaP.SetXYZ(px, py, pz); }
  
  TVector3	GetBachP()  const { return fBachP; }
  void		SetBachP(Double_t px, Double_t py, Double_t pz) { fBachP.SetXYZ(px, py, pz); }
  
  virtual const Text_t* GetName()  const   { return Form("ESDcascade_%i",fESDIndex); }
  virtual const Text_t* GetTitle() const   { return Form("ESDcascade_%i",fESDIndex); }

  TEveTrackPropagator* GetPropagatorBac() const  { return fRnrStyleBac; }
  TEveTrackPropagator* GetPropagatorNeg() const  { return fRnrStyleNeg; }
  TEveTrackPropagator* GetPropagatorPos() const  { return fRnrStylePos; }

  AliEveTrack* 	GetBacTrack()      const { return fBacTrack; }
  AliEveTrack* 	GetNegTrack()      const { return fNegTrack; }
  AliEveTrack* 	GetPosTrack()      const { return fPosTrack; }

  TEveLine*  	GetPointingCurve() const { return fPointingCurve; }
  TEveLine*  	GetV0Path()        const { return fV0Path; }

   
protected:
  TEveVector 		fRecBirthV;     // Assumed birth point of cascade
  TEveVector 		fRecDecayV;     // Xi decay point : point of closest approach between the Xi daughters
  TEveVector 		fRecDecayP;     // Reconstructed momentum of the cascade, at the Xi decay
  TEveVector 		fRecDecayV0;    // Reconstructed birth point of neutral daughter
  
  
  AliEveTrack 		*fBacTrack;	 //! Eve track for the bachelor of the cascade
  AliEveTrack 		*fNegTrack;	 //! Eve track for the neg V0 dghter, within the cascade
  AliEveTrack 		*fPosTrack;	 //! Eve track for the pos V0 dghter, within the cascade

  TEveTrackPropagator 	*fRnrStyleBac;	 //! track propagator for bachelor track
  TEveTrackPropagator 	*fRnrStyleNeg;	 //! track propagator for negative track
  TEveTrackPropagator 	*fRnrStylePos;	 //! track propagator for positive track

  TEveLine         	*fPointingCurve; //! Curve meant model the Xi trajectory
  TEveLine         	*fV0Path;	 //! Line meant to model the V0 path of the cascade

  Int_t  		fESDIndex;       // Index in ESD Cascade array.
  Float_t 		fDaughterDCA;    // Distance at the point of closest approach, between both Xi daughters
  Float_t 		fChi2Cascade;    // Some Chi-square.
  TVector3		fLambdaP;	 // Momentum of Lambda (V0 in cascade), at its decay point
  TVector3		fBachP;	 	 // Momentum of Bachelor, at the Xi decay point

private:
  AliEveCascade(const AliEveCascade&);            // Not implemented
  AliEveCascade& operator=(const AliEveCascade&); // Not implemented

  ClassDef(AliEveCascade, 1); // Visual representation of a AliEveCascade.
};


/******************************************************************************/
// AliEveCascadeList
/******************************************************************************/

class AliEveCascadeList : public TEveElementList
{
  friend class AliEveCascadeListEditor;

public:
  AliEveCascadeList();
  AliEveCascadeList(TEveTrackPropagator* rsBac,TEveTrackPropagator* rsNeg,TEveTrackPropagator* rsPos);
  AliEveCascadeList(const Text_t* name, TEveTrackPropagator* rsBac=0, TEveTrackPropagator* rsNeg=0, TEveTrackPropagator* rsPos=0);
  virtual ~AliEveCascadeList() {}

  virtual const Text_t* GetTitle() const { return fTitle; }
  virtual void SetTitle(const Text_t* t) { fTitle = t; }
  virtual void SetTracksColor(Color_t cBac) { fBacColor = cBac;}

  virtual Bool_t CanEditMainColor() const { return kTRUE; }

  void  SetRnrStyleBac(TEveTrackPropagator* const rst ) { fRnrStyleBac = rst; }
  void  SetRnrStyleNeg(TEveTrackPropagator* const rst ) { fRnrStyleNeg = rst; }
  void  SetRnrStylePos(TEveTrackPropagator* const rst ) { fRnrStylePos = rst; }
  TEveTrackPropagator* GetPropagatorBac()   const       { return fRnrStyleBac; }
  TEveTrackPropagator* GetPropagatorNeg()   const       { return fRnrStyleNeg; }
  TEveTrackPropagator* GetPropagatorPos()   const       { return fRnrStylePos; }

  Bool_t GetRnrCascadevtx()     const { return fRnrCascadevtx; }
  Bool_t GetRnrCascadepath()    const { return fRnrCascadepath; }
  Bool_t GetRnrDaughters()      const { return fRnrDaughters; }

  void   MakeCascades();

  void   FilterByRadius        (Float_t minR, Float_t maxR);
  void   FilterByDaughterDCA   (Float_t minDaughterDCA, Float_t maxDaughterDCA);
  void   FilterByPt            (Float_t minPt, Float_t maxPt);
  void   FilterByInvariantMass (Float_t minInvariantMass, Float_t maxInvariantMass, Int_t cascadePdgCodeHyp);
  
  
  void   SetInvMassHyp		(Int_t rInvMassHyp) {fInvMassHyp = rInvMassHyp;}
  Int_t  GetInvMassHyp() 	const { return fInvMassHyp; }

protected:
  TString              fTitle;			// title

  TEveTrackPropagator *fRnrStyleBac;		//! Rnr Style of bachelor track
  TEveTrackPropagator *fRnrStyleNeg;		//! Rnr Style of negative track
  TEveTrackPropagator *fRnrStylePos;		//! Rnr Style of positive track

  Bool_t               fRnrDaughters;		// Render state for the cascade daughters
  Bool_t               fRnrCascadevtx;		// Render state for the cascade decay point
  Bool_t               fRnrCascadepath;		// Render state for the path between the prim. vertex and the "Xi" decay point

  Color_t              fBacColor;		// Color of the bachelor track

  Float_t              fMinRCut;		// Min transv. radius allowed for cascade selection
  Float_t              fMaxRCut;		// Max transv. radius allowed for cascade selection

  Float_t              fMinDaughterDCA;		// Min DCA between Xi daughters, allowed for cascade selection
  Float_t              fMaxDaughterDCA;		// Max DCA between Xi daughters, allowed for cascade selection

  Float_t              fMinPt;			// Min pt allowed for cascade selection
  Float_t              fMaxPt;			// Max pt allowed for cascade selection
  
  Int_t                fInvMassHyp;		// PdgCode of the inv. mass hypothesis for the cascade
  
  Float_t              fMinInvariantMass; 	// Minimum invariant mass cut.for cascade
  Float_t              fMaxInvariantMass; 	// Maximum invariant mass cut.for cascade

private:
  void Init();

  AliEveCascadeList(const AliEveCascadeList&);            // Not implemented
  AliEveCascadeList& operator=(const AliEveCascadeList&); // Not implemented

  ClassDef(AliEveCascadeList, 0); // A list of AliEveCascade objecs.
};


#endif
