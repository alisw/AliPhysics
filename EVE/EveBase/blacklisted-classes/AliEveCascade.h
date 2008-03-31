// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveCascade_H
#define AliEveCascade_H

/***********************************************************************
*  This code defines the reconstructed cascades visualized with EVE
*
* Ludovic Gaudichet (gaudichet@to.infn.it)
************************************************************************/

#include <TEveVSDStructs.h>
#include <TEveElement.h>
#include <TEveTrack.h>

#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>

class TH1F;
class TH2F;


class CascadeList;

class AliEveCascade : public TEveElement,
                      public TPolyMarker3D
{
  friend class CascadeList;

public:
  typedef std::vector<TEvePathMark*>           vpPathMark_t;


  AliEveCascade();
  AliEveCascade(TEveTrackPropagator* rs);
  virtual ~AliEveCascade();

  virtual void  SetESDIndex(Int_t ind) { fESDIndex = ind;}
  virtual void  SetMainColor(Color_t col)
  {
    fMarkerColor = col; fMainColorPtr = &fMarkerColor;
    fPolyLineV0.SetLineColor(fMarkerColor);
  }
  virtual void  SetTracksColor(Color_t cNeg, Color_t cPos, Color_t cBach)
  {
    fPolyLineNeg.SetLineColor(cNeg); fPolyLinePos.SetLineColor(cPos);
    fPolyLineBach.SetLineColor(cBach);
  }
  void        SetRnrStyle(TEveTrackPropagator* rs) { fRnrStyle = rs; }

  void  AddPathMarkPos(TEvePathMark* pm) { fPathMarksPos.push_back(pm); }
  void  AddPathMarkNeg(TEvePathMark* pm) { fPathMarksNeg.push_back(pm); }
  void  AddPathMarkBach(TEvePathMark* pm) { fPathMarksBach.push_back(pm); }

  virtual void PaintV0Daughters(Option_t* option="") {
    if(fRnrSelf) {fPolyLineNeg.Paint(option);fPolyLinePos.Paint(option); } }
  virtual void PaintBachelor(Option_t* option="") {
    if(fRnrSelf) fPolyLineBach.Paint(option); }
  virtual void Paint(Option_t* option="") {
    if(fRnrSelf) TPolyMarker3D::Paint(option);}
  virtual void PaintV0Path(Option_t* option="") {
    if(fRnrSelf) fPolyLineV0.Paint(option);}
  virtual void PaintCasPath(Option_t* option="") {
    if(fRnrSelf) fPolyLineCas.Paint(option);}

  void Reset(TPolyLine3D* polyLine);
  void MakeTrack(vpPathMark_t& pathMark, TEveVector& vtx,  TEveVector& p,
		 Int_t charge, Float_t beta, TPolyLine3D& polyLine);
  void MakeV0path();
  void MakeCasPath();
  void MakeCascade();

  void SetBeta(Float_t betaNeg, Float_t betaPos, Float_t betaBach);
  void SetDCA_v0_Bach(Float_t dca) {fDCA_v0_Bach = dca;}
  void SetCasCosPointingAngle(Float_t cos) {fCasCosPointingAngle = cos;}
  void SetNegP(Float_t px, Float_t py, Float_t pz) {fP_neg.x = px; fP_neg.y = py; fP_neg.z = pz;}
  void SetPosP(Float_t px, Float_t py, Float_t pz) {fP_pos.x = px; fP_pos.y = py; fP_pos.z = pz;}
  void SetBachP(Float_t px, Float_t py, Float_t pz) {fP_bach.x = px; fP_bach.y = py; fP_bach.z = pz;}

  void SetV0vtx(Float_t vx, Float_t vy, Float_t vz)  {
    fV_neg.x = vx; fV_neg.y = vy; fV_neg.z = vz;
    fV_pos.x = vx; fV_pos.y = vy; fV_pos.z = vz;
  }
  void SetCascadeVtx(Float_t vx, Float_t vy, Float_t vz) {
    fV_decay.x = vx; fV_decay.y = vy; fV_decay.z = vz; }

  void SetDecayLength(Float_t primx, Float_t primy, Float_t primz);

  Int_t   GetESDIndex() const { return fESDIndex; }
  virtual const Text_t* GetName()  const { return Form("ESDcascade_%i",fESDIndex); }
  virtual const Text_t* GetTitle() const { return Form("ESDcascade_%i",fESDIndex); }

  Float_t GetDCA_v0_Bach() const;
  Float_t GetCasCosPointingAngle() const;
  Float_t GetRadius() const;
  Float_t GetPseudoRapidity() const;
  Float_t GetPt2() const;
  Float_t GetPt() const;
  Float_t GetP2() const;
  Float_t GetMomentum() const;
  Float_t GetPx() const;
  Float_t GetPy() const;
  Float_t GetPz() const;
  Float_t GetCasAlphaArmenteros() const;
  Float_t GetCasPtArmenteros() const;
  Float_t GetPosP2() const;
  Float_t GetPosP() const;
  Float_t GetPosPt() const;
  Float_t GetPosPseudoRapidity() const;
  Float_t GetNegP2() const;
  Float_t GetNegP() const;
  Float_t GetNegPt() const;
  Float_t GetNegPseudoRapidity() const;
  Float_t GetBachP2() const;
  Float_t GetBachP() const;
  Float_t GetBachPt() const;
  Float_t GetBachPseudoRapidity() const;

  Float_t GetV0P2() const;
  Float_t GetLambdaE() const;
  Float_t GetXiE() const;
  Float_t GetOmegaE() const;
  Float_t GetXiMass() const;
  Float_t GetAntiXiMass() const;
  Float_t GetOmegaMass() const;
  Float_t GetAntiOmegaMass() const;


protected:
  typedef std::vector<TEvePathMark*>::iterator vpPathMark_i;

  TEveVector fV_neg;       // Vertex of negative track
  TEveVector fP_neg;       // Momentum of negative track
  TEveVector fV_pos;       // Vertex of positive track
  TEveVector fP_pos;       // Momentum of positive track
  TEveVector fV_bach;      // Vertex of positive track
  TEveVector fP_bach;      // Momentum of positive track

  TEveVector fV_decay;     //decay point of the cascade
  TEveVector fV_birth;    // Reconstructed birth point of neutral particle

  vpPathMark_t         fPathMarksNeg;
  vpPathMark_t         fPathMarksPos;
  vpPathMark_t         fPathMarksBach;
  TEveTrackPropagator *fRnrStyle;

  TPolyLine3D       fPolyLineNeg;
  TPolyLine3D       fPolyLinePos;
  TPolyLine3D       fPolyLineBach;
  TPolyLine3D       fPolyLineV0;   // line of AliEveV0 travel
  TPolyLine3D       fPolyLineCas;  // line of cascade travel

  Float_t           fBeta_neg;
  Float_t           fBeta_pos;
  Float_t           fBeta_bach;

  Int_t             fESDIndex;

  Float_t           fDCA_v0_Bach;
  Float_t           fCasCosPointingAngle;
  Float_t           fCasDecayLength;

  static const Float_t fgkMassPion2;
  static const Float_t fgkMassKaon2;
  static const Float_t fgkMassProton2;
  static const Float_t fgkMassLambda2;

private:
  AliEveCascade(const AliEveCascade&);            // Not implemented
  AliEveCascade& operator=(const AliEveCascade&); // Not implemented

  ClassDef(AliEveCascade, 0); // Visual representation of a cascade.
};



//______________________________________________________________________________

inline void AliEveCascade::SetBeta(Float_t betaNeg, Float_t betaPos, Float_t betaBach) {
   fBeta_neg = betaNeg;
   fBeta_pos = betaPos;
   fBeta_bach = betaBach;
 }


//______________________________________________________________________________

inline Float_t AliEveCascade::GetV0P2() const {
  Float_t px = fP_neg.x + fP_pos.x, py = fP_neg.y + fP_pos.y,
    pz = fP_neg.z+fP_pos.z;
  return px*px + py*py + pz*pz;
}


inline Float_t AliEveCascade::GetLambdaE() const {
  return sqrt(fgkMassLambda2+GetV0P2());
}

inline Float_t AliEveCascade::GetXiE() const {
  Float_t e = GetLambdaE() +
    sqrt(fgkMassPion2 + fP_bach.x*fP_bach.x + fP_bach.y*fP_bach.y +
	 fP_bach.z*fP_bach.z);
  return e;
}

inline Float_t AliEveCascade::GetOmegaE() const {
  Float_t e = GetLambdaE() +
    sqrt(fgkMassKaon2 + fP_bach.x*fP_bach.x + fP_bach.y*fP_bach.y +
	 fP_bach.z*fP_bach.z);
  return e;
}

inline Float_t AliEveCascade::GetXiMass() const {
  Float_t e = GetXiE();
  return sqrt(e*e - GetP2());
}

inline Float_t AliEveCascade::GetAntiXiMass() const { return GetXiMass();}

inline Float_t AliEveCascade::GetOmegaMass() const {
  Float_t e = GetOmegaE();
  return sqrt(e*e - GetP2());
}

inline Float_t AliEveCascade::GetAntiOmegaMass() const { return GetOmegaMass();}


//______________________________________________________________________________

inline Float_t AliEveCascade::GetDCA_v0_Bach() const {
  return fDCA_v0_Bach;
}

inline Float_t AliEveCascade::GetCasCosPointingAngle() const {
return fCasCosPointingAngle;
}

inline Float_t AliEveCascade::GetRadius() const {
  return sqrt(fV_birth.x*fV_birth.x + fV_birth.y*fV_birth.y);
}

//inline Float_t AliEveCascade::GetDecayLength() const {
//return fDecayLength;
//}

inline Float_t AliEveCascade::GetPseudoRapidity() const {
  Float_t theta = acos( GetPz()/GetMomentum() );
  return ( -log(tan(theta/2.)) );
}


//______________________________________________________________________________
inline Float_t AliEveCascade::GetPt2() const {
  Float_t px = GetPx(), py = GetPy();
  return (px*px+py*py);
}

inline Float_t AliEveCascade::GetP2() const {

  Float_t px = GetPx(), py = GetPy(), pz = GetPz();
  return (px*px+py*py+pz*pz);
}

inline Float_t AliEveCascade::GetPt() const {
  return sqrt(GetPt2());
}

inline Float_t AliEveCascade::GetMomentum() const {
  return sqrt(GetP2());
}

inline Float_t AliEveCascade::GetPx() const {
  return (fP_neg.x + fP_pos.x + fP_bach.x);
}

inline Float_t AliEveCascade::GetPy() const {
  return (fP_neg.y + fP_pos.y + fP_bach.y);
}

inline Float_t AliEveCascade::GetPz() const {
  return (fP_neg.z + fP_pos.z + fP_bach.z);
}

//______________________________________________________________________________

inline Float_t AliEveCascade::GetPosP2() const {
  return (fP_pos.x*fP_pos.x + fP_pos.y*fP_pos.y + fP_pos.z*fP_pos.z);
}

inline Float_t AliEveCascade::GetPosP() const {
  return sqrt(GetPosP2());
}

inline Float_t AliEveCascade::GetPosPt() const {
  return sqrt(fP_pos.x*fP_pos.x + fP_pos.y*fP_pos.y);
}

inline Float_t AliEveCascade::GetPosPseudoRapidity() const {
  Float_t theta = acos( fP_pos.z/GetPosP() );
  return ( -log(tan(theta/2.)) );
}

//______________________________________________________________________________
inline Float_t AliEveCascade::GetNegP2() const {
  return (fP_neg.x*fP_neg.x + fP_neg.y*fP_neg.y + fP_neg.z*fP_neg.z);
}

inline Float_t AliEveCascade::GetNegP() const {
  return sqrt(GetNegP2());
}

inline Float_t AliEveCascade::GetNegPt() const {
  return sqrt(fP_neg.x*fP_neg.x + fP_neg.y*fP_neg.y);
}

inline Float_t AliEveCascade::GetNegPseudoRapidity() const {
  Float_t theta = acos( fP_neg.z/GetNegP() );
  return ( -log(tan(theta/2.)) );
}


//______________________________________________________________________________
inline Float_t AliEveCascade::GetBachP2() const {
  return (fP_bach.x*fP_bach.x + fP_bach.y*fP_bach.y + fP_bach.z*fP_bach.z);
}

inline Float_t AliEveCascade::GetBachP() const {
  return sqrt(GetBachP2());
}

inline Float_t AliEveCascade::GetBachPt() const {
  return sqrt(fP_bach.x*fP_bach.x + fP_bach.y*fP_bach.y);
}

inline Float_t AliEveCascade::GetBachPseudoRapidity() const {
  Float_t theta = acos( fP_bach.z/GetBachP() );
  return ( -log(tan(theta/2.)) );
}


/***********************************************************************
*
*  CascadeList class
*
************************************************************************/

class CascadeList : public TEveElementList
{
public:
  CascadeList(TEveTrackPropagator* rs=0);
  CascadeList(const Text_t* name, TEveTrackPropagator* rs=0);

  virtual const Text_t* GetTitle() const { return fTitle; }
  virtual void SetTitle(const Text_t* t) { fTitle = t; }
  virtual void SetTracksColor(Color_t cNeg, Color_t cPos, Color_t cBach)
  { fNegColor = cNeg; fPosColor = cPos; fBachColor = cBach; }

  virtual Bool_t CanEditMainColor()  { return kTRUE; }

  virtual void Paint(Option_t* option="");

  void  SetRnrStyle(TEveTrackPropagator* rst) { fRnrStyle= rst; }
  TEveTrackPropagator* GetPropagator()          { return fRnrStyle; }

  Bool_t GetRnrCasVtx() const { return fRnrCasVtx; }
  Bool_t GetRnrCasPath() const { return fRnrCasPath; }
  Bool_t GetRnrV0vtx() const { return fRnrV0vtx; }
  Bool_t GetRnrV0path() const { return fRnrV0path; }
  Bool_t GetRnrV0Daughters() const { return fRnrV0Daughters; }
  Bool_t GetRnrBachelor() const { return fRnrBach; }

  void   SetRnrV0vtx(Bool_t);
  void   SetRnrV0path(Bool_t);
  void   SetRnrV0Daughters(Bool_t);
  void   SetRnrBachelor(Bool_t);
  void   SetRnrCasPath(Bool_t);
  void   SetRnrCasVtx(Bool_t);
  void   SetMin(Int_t i, Float_t val) {
    if ((i>=0)&&(i<fgkNcutVar)) fMin[i]=val;}
  void   SetMax(Int_t i, Float_t val) {
    if ((i>=0)&&(i<fgkNcutVar)) fMax[i]=val;}

  void   MakeCascades();

  TH1F*   GetHist(Int_t i) {
    if ((i>=0)&&(i<fgkNcutVar)) return fHist[i]; else return 0;}
  TH2F*   GetHist2D(Int_t i) {
    if ((i>=0)&&(i<fgkNcutVar2D)) return fHist2D[i]; else return 0;}
  Float_t GetMin(Int_t i) {
    if ((i>=0)&&(i<fgkNcutVar)) return fMin[i]; else return 0;}
  Float_t GetMax(Int_t i) {
    if ((i>=0)&&(i<fgkNcutVar)) return fMax[i]; else return 0;}
  void GetCasIndexRange(Int_t &imin, Int_t &imax);

  void AdjustHist(Int_t iHist);
  void UnFill(AliEveCascade* cas);
  void Filter(AliEveCascade* cas);
  void FilterAll();

  void XiMassFilter(Float_t min, Float_t max);
  void OmegaMassFilter(Float_t min, Float_t max);
  void IndexFilter(Float_t min, Float_t max);
  void CosPointingFilter(Float_t min, Float_t max);
  void BachV0DCAFilter(Float_t min, Float_t max);
  void RadiusFilter(Float_t min, Float_t max);
  void PtFilter(Float_t min, Float_t max);
  void PseudoRapFilter(Float_t min, Float_t max);
  void NegPtFilter(Float_t min, Float_t max);
  void NegEtaFilter(Float_t min, Float_t max);
  void PosPtFilter(Float_t min, Float_t max);
  void PosEtaFilter(Float_t min, Float_t max);
  void BachPtFilter(Float_t min, Float_t max);
  void BachEtaFilter(Float_t min, Float_t max);

protected:
  TString              fTitle;

  TEveTrackPropagator *fRnrStyle;

  Bool_t               fRnrBach;
  Bool_t               fRnrV0Daughters;
  Bool_t               fRnrV0vtx;
  Bool_t               fRnrV0path;
  Bool_t               fRnrCasVtx;
  Bool_t               fRnrCasPath;

  Color_t              fNegColor;
  Color_t              fPosColor;
  Color_t              fBachColor;

  static const Int_t fgkNcutVar = 14;
  TH1F *fHist[fgkNcutVar];
  Float_t fMin[fgkNcutVar];
  Float_t fMax[fgkNcutVar];

  static const Int_t fgkNcutVar2D = 1;
  TH2F *fHist2D[fgkNcutVar2D];
  Float_t fMinX[fgkNcutVar2D];
  Float_t fMinY[fgkNcutVar2D];
  Float_t fMaxX[fgkNcutVar2D];
  Float_t fMaxY[fgkNcutVar2D];

private:
  void  Init();

  CascadeList(const CascadeList&);            // Not implemented
  CascadeList& operator=(const CascadeList&); // Not implemented

  ClassDef(CascadeList, 0); // A list of AliEveCascade objects.
};

#endif
