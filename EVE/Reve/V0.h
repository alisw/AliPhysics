
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/***********************************************************************
*  This code defines the reconstructed v0 visualized with EVE
*
* Ludovic Gaudichet (gaudichet@to.infn.it)
************************************************************************/


#ifndef REVE_V0_H
#define REVE_V0_H

#include <Reve/PODs.h>
#include <Reve/RenderElement.h>
#include <Reve/Track.h>

#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>

class TH1F;
class TH2F;


namespace Reve {

class V0List;

class V0 : public RenderElement, public TPolyMarker3D
{
  friend class V0List;

  V0(const V0&);            // Not implemented
  V0& operator=(const V0&); // Not implemented

public: 
  V0();
  V0(Reve::RecTrack* tNeg, Reve::RecTrack* tPos,Reve::RecV0* v0, TrackRnrStyle* rs);
  virtual ~V0();

  typedef std::vector<Reve::PathMark*>           vpPathMark_t;
  typedef std::vector<Reve::PathMark*>::iterator vpPathMark_i;
  void  AddPathMarkPos(Reve::PathMark* pm) { fPathMarksPos.push_back(pm); }
  void  AddPathMarkNeg(Reve::PathMark* pm) { fPathMarksNeg.push_back(pm); }

  void Reset(TPolyLine3D* polyLine);

  void MakeTrack(vpPathMark_t& pathMark, Reve::Vector& vtx,  Reve::Vector& p,
		 Int_t charge, Float_t beta, TPolyLine3D& polyLine);

  void MakeV0path();
  void MakeV0();

  virtual void PaintDaughters(Option_t* option="") {
    if(fRnrElement) {fPolyLineNeg.Paint(option);fPolyLinePos.Paint(option);} }
  virtual void Paint(Option_t* option="") {
    if(fRnrElement) TPolyMarker3D::Paint(option);}
  virtual void PaintPath(Option_t* option="") {
    if(fRnrElement) fPolyLineV0.Paint(option);}

  virtual void  SetMainColor(Color_t col) {
    fMarkerColor = col; fMainColorPtr = &fMarkerColor;
    fPolyLineV0.SetLineColor(fMarkerColor);}
  virtual void  SetTracksColor(Color_t cNeg, Color_t cPos) {
    fPolyLineNeg.SetLineColor(cNeg); fPolyLinePos.SetLineColor(cPos); }
  void          SetRnrStyle(TrackRnrStyle* rs) { fRnrStyle = rs; }
  void          SetESDIndex(Int_t ind) { fESDIndex = ind;}
  void          SetDaughterDCA(Float_t dca);
  void          SetCosPointingAngle(Float_t cos);
  void          SetDecayLength(Float_t len);
  void          SetDecayLength(Float_t primx, Float_t primy, Float_t primz);

  Float_t const GetDaughterDCA();
  Float_t const GetCosPointingAngle();
  Float_t const GetRadius();
  Float_t const GetDecayLength();

  Float_t const GetP2();
  Float_t const GetMomentum();
  Float_t const GetPt();
  Float_t const GetPt2();
  Float_t const GetPx();
  Float_t const GetPy();
  Float_t const GetPz();
  Float_t const GetPseudoRapidity();
  Float_t const GetAlphaArmenteros();
  Float_t const GetPtArmenteros();

  Float_t const GetK0sE() {return 0;};
  Float_t const GetLamE() {return 0;};
  Float_t const GetAntiLamE() {return 0;};
  Float_t const GetK0mass();
  Float_t const GetLamMass();
  Float_t const GetAntiLamMass();

  Float_t const GetPionMinusE();
  Float_t const GetPionPlusE();
  Float_t const GetProtonE();
  Float_t const GetPBarE();

  Float_t const GetPosDCAtoPrim();
  Float_t const GetPosP2();
  Float_t const GetPosP();
  Float_t const GetPosPt();
  Float_t const GetPosPseudoRapidity();

  Float_t const GetNegDCAtoPrim();
  Float_t const GetNegP2();
  Float_t const GetNegP();
  Float_t const GetNegPt();
  Float_t const GetNegPseudoRapidity();
  Int_t   const GetESDIndex() {return fESDIndex;};

  virtual const Text_t* GetName() const    { return Form("ESDv0_%i",fESDIndex); }
  virtual const Text_t* GetTitle() const   { return Form("ESDv0_%i",fESDIndex); }
  Int_t          GetLabelPos() const { return fLabel_pos; }
  Int_t          GetLabelNeg() const { return fLabel_neg; }
  TrackRnrStyle* GetRnrStyle() const  { return fRnrStyle; }
  TPolyLine3D*   GetPolyLineNeg() {return &fPolyLineNeg;}
  TPolyLine3D*   GetPolyLinePos() {return &fPolyLinePos;}
  TPolyLine3D*   GetPolyLineV0() {return &fPolyLineV0;}

protected:

  Reve::Vector fV_neg;       // Vertex of negative track
  Reve::Vector fP_neg;       // Momentum of negative track
  Reve::Vector fV_pos;       // Vertex of positive track
  Reve::Vector fP_pos;       // Momentum of positive track

  Reve::Vector fV_v0;        // Point of closest approach
  Reve::Vector fV0_birth;    // Reconstucted birth point of neutral particle

  Float_t           fBeta_neg;
  Float_t           fBeta_pos;

  Int_t             fLabel_neg;
  Int_t             fLabel_pos;

  vpPathMark_t      fPathMarksNeg;
  vpPathMark_t      fPathMarksPos;
  TrackRnrStyle*    fRnrStyle;

  TPolyLine3D       fPolyLineNeg;
  TPolyLine3D       fPolyLinePos;
  TPolyLine3D       fPolyLineV0;

  Int_t             fESDIndex;
  Float_t           fDaughterDCA;
  Float_t           fCosPointingAngle;
  Float_t           fDecayLength;

  static const Float_t fgkMassPion2;
  static const Float_t fgkMassProton2;

  ClassDef(V0, 1)
}; // endclass V0



//______________________________________________________________________
inline void V0::SetDaughterDCA(Float_t dca) { 
  fDaughterDCA = dca; 
}

inline void V0::SetCosPointingAngle(Float_t cos) { 
  fCosPointingAngle = cos; 
}

inline void V0::SetDecayLength(Float_t len) {
  fDecayLength = len;
}


//______________________________________________________________________
inline Float_t const V0::GetPt2() {
  Float_t px = fP_neg.x+fP_pos.x, py = fP_neg.y+fP_pos.y;
  return (px*px+py*py);
}

inline Float_t const V0::GetP2() {

  Float_t px = fP_neg.x+fP_pos.x, py = fP_neg.y+fP_pos.y, pz = fP_neg.z+fP_pos.z;
  return (px*px+py*py+pz*pz);
}

inline Float_t const V0::GetPt() {
  return sqrt(GetPt2());
}

inline Float_t const V0::GetMomentum() {
  return sqrt(GetP2());
}

inline Float_t const V0::GetPx() {
  return (fP_neg.x+fP_pos.x);
}

inline Float_t const V0::GetPy() {
  return (fP_neg.y+fP_pos.y);
}

inline Float_t const V0::GetPz() {
  return (fP_neg.z+fP_pos.z);
}

//______________________________________________________________________

inline Float_t const V0::GetDaughterDCA() {
return fDaughterDCA;
}

inline Float_t const V0::GetCosPointingAngle() {
return fCosPointingAngle;
}

inline Float_t const V0::GetRadius() {
  return sqrt(fV_v0.x*fV_v0.x + fV_v0.y*fV_v0.y);
}

inline Float_t const V0::GetDecayLength() {
return fDecayLength;
}

inline Float_t const V0::GetPseudoRapidity() {
  Float_t theta = acos( GetPz()/GetMomentum() );
  return ( -log(tan(theta/2.)) );
}

//______________________________________________________________________

inline Float_t const V0::GetPionMinusE() {
  return sqrt(fgkMassPion2+GetNegP2());
}

inline Float_t const V0::GetPionPlusE() {
  return sqrt(fgkMassPion2+GetPosP2());
}
inline Float_t const V0::GetProtonE() {
  return sqrt(fgkMassProton2+GetPosP2());

}
inline Float_t const V0::GetPBarE() {
  return sqrt(fgkMassProton2+GetNegP2());
}

//______________________________________________________________________

inline Float_t const V0::GetPosP2() {
  return (fP_pos.x*fP_pos.x + fP_pos.y*fP_pos.y + fP_pos.z*fP_pos.z);
}

inline Float_t const V0::GetPosP() {
  return sqrt(GetPosP2());
}

inline Float_t const V0::GetPosPt() {
  return sqrt(fP_pos.x*fP_pos.x + fP_pos.y*fP_pos.y);
}

inline Float_t const V0::GetPosPseudoRapidity() {
  Float_t theta = acos( fP_pos.z/GetPosP() );
  return ( -log(tan(theta/2.)) );
}

inline Float_t const V0::GetPosDCAtoPrim() {
  return 0;
}

//______________________________________________________________________
inline Float_t const V0::GetNegP2() {
  return (fP_neg.x*fP_neg.x + fP_neg.y*fP_neg.y + fP_neg.z*fP_neg.z);
}

inline Float_t const V0::GetNegP() {
  return sqrt(GetNegP2());
}

inline Float_t const V0::GetNegPt() {
  return sqrt(fP_neg.x*fP_neg.x + fP_neg.y*fP_neg.y);
}

inline Float_t const V0::GetNegPseudoRapidity() {
  Float_t theta = acos( fP_neg.z/GetNegP() );
  return ( -log(tan(theta/2.)) );
}

inline Float_t const V0::GetNegDCAtoPrim() {
  return 0;
}

//______________________________________________________________________

inline Float_t const V0::GetK0mass() {
  Float_t energy = GetPionMinusE() + GetPionPlusE();
  return sqrt( energy*energy - GetP2() );
}

inline Float_t const V0::GetLamMass() {
  Float_t energy = GetPionMinusE() + GetProtonE();
  return sqrt( energy*energy - GetP2() );
}

inline Float_t const V0::GetAntiLamMass() {
  Float_t energy = GetPionPlusE() + GetPBarE();
  return sqrt( energy*energy - GetP2() );
}




/**************************************************************************/
// V0List
/**************************************************************************/

class V0List : public TNamed, public RenderElementListBase
{
  V0List(const V0List&);            // Not implemented
  V0List& operator=(const V0List&); // Not implemented

public:
  V0List();
  V0List(TrackRnrStyle* rs);
  V0List(const Text_t* name, TrackRnrStyle* rs=0);
  virtual ~V0List();

  virtual const Text_t* GetTitle() const { return fTitle; }
  virtual void SetTitle(const Text_t* t) { fTitle = t; }
  virtual void SetTracksColor(Color_t cNeg, Color_t cPos) {
    fNegColor = cNeg; fPosColor = cPos;}

  virtual Bool_t CanEditMainColor()  { return kTRUE; }

  virtual void Paint(Option_t* option="");

  virtual void AddElement(RenderElement* el);

  void  SetRnrStyle(TrackRnrStyle* rst) { fRnrStyle= rst; }
  TrackRnrStyle* GetRnrStyle()          { return fRnrStyle; } 

  Bool_t GetRnrV0vtx() const { return fRnrV0vtx; }
  Bool_t GetRnrV0path() const { return fRnrV0path; }
  Bool_t GetRnrDaughters() const { return fRnrDaughters; }

  void   SetRnrV0vtx(Bool_t);
  void   SetRnrV0path(Bool_t);
  void   SetRnrDaughters(Bool_t);
  void   SetMin(Int_t i, Float_t val) {
    if ((i>=0)&&(i<fgkNcutVar)) fMin[i]=val;}
  void   SetMax(Int_t i, Float_t val) {
    if ((i>=0)&&(i<fgkNcutVar)) fMax[i]=val;}

  void   MakeV0s();
  void   MakeMarkers();

  TH1F*   GetHist(Int_t i) {
    if ((i>=0)&&(i<fgkNcutVar)) return fHist[i]; else return 0;}
  TH2F*   GetHist2D(Int_t i) {
    if ((i>=0)&&(i<fgkNcutVar2D)) return fHist2D[i]; else return 0;}
  Float_t GetMin(Int_t i) {
    if ((i>=0)&&(i<fgkNcutVar)) return fMin[i]; else return 0;}
  Float_t GetMax(Int_t i) {
    if ((i>=0)&&(i<fgkNcutVar)) return fMax[i]; else return 0;}
  void GetV0IndexRange(Int_t &imin, Int_t &imax);

  Float_t GetMaxR()         const { return fRnrStyle->fMaxZ; }
  Float_t GetMaxZ()         const { return fRnrStyle->fMaxR; }
  Float_t GetMaxOrbs()      const { return fRnrStyle->fMaxOrbs; }
  Float_t GetMinAng()       const { return fRnrStyle->fMinAng; }
  Float_t GetDelta()        const { return fRnrStyle->fDelta; }
  Bool_t  GetFitDaughters() const { return fRnrStyle->fFitDaughters; }
  Bool_t  GetFitDecay()     const { return fRnrStyle->fFitDecay; }

  void AdjustHist(Int_t iHist);
  void UnFill(V0* v0);
  void Filter(V0* v0);
  void FilterAll();
  void PtFilter(Float_t min, Float_t max);

  void K0sMFilter(Float_t min, Float_t max);
  void LamMFilter(Float_t min, Float_t max);
  void ALamMFilter(Float_t min, Float_t max);
  void CosPointingFilter(Float_t min, Float_t max);
  void DaughterDCAFilter(Float_t min, Float_t max);
  void RadiusFilter(Float_t min, Float_t max);
  void EtaFilter(Float_t min, Float_t max);
  void NegPtFilter(Float_t min, Float_t max);
  void NegEtaFilter(Float_t min, Float_t max);
  void PosPtFilter(Float_t min, Float_t max);
  void PosEtaFilter(Float_t min, Float_t max);
  void IndexFilter(Float_t min, Float_t max);


private:
  void  Init();

protected:
  TString              fTitle;

  TrackRnrStyle*       fRnrStyle;

  Bool_t               fRnrDaughters;
  Bool_t               fRnrV0vtx;
  Bool_t               fRnrV0path;

  Color_t              fNegColor;
  Color_t              fPosColor;

  static const Int_t fgkNcutVar = 13;
  TH1F *fHist[fgkNcutVar];
  Float_t fMin[fgkNcutVar];
  Float_t fMax[fgkNcutVar];

  static const Int_t fgkNcutVar2D = 1;
  TH2F *fHist2D[fgkNcutVar2D];
  Float_t fMinX[fgkNcutVar2D];
  Float_t fMinY[fgkNcutVar2D];
  Float_t fMaxX[fgkNcutVar2D];
  Float_t fMaxY[fgkNcutVar2D];

  ClassDef(V0List, 1)
};


} // namespace Reve

#endif
