#ifndef ALIEVE_CASCADE_H
#define ALIEVE_CASCADE_H

/***********************************************************************
*  This code defines the reconstructed cascades visualized with EVE
*
* Ludovic Gaudichet (gaudichet@to.infn.it)
************************************************************************/

#include <Reve/PODs.h>
#include <Reve/RenderElement.h>
#include <Reve/Track.h>

#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>

class TH1F;
class TH2F;


namespace Alieve {

class CascadeList;

class Cascade : public Reve::RenderElement,
                public TPolyMarker3D
{
public:
  typedef std::vector<Reve::PathMark*>           vpPathMark_t;

private:
  friend class CascadeList;

  Cascade(const Cascade&);            // Not implemented
  Cascade& operator=(const Cascade&); // Not implemented

protected:
  typedef std::vector<Reve::PathMark*>::iterator vpPathMark_i;

  Reve::Vector fV_neg;       // Vertex of negative track
  Reve::Vector fP_neg;       // Momentum of negative track
  Reve::Vector fV_pos;       // Vertex of positive track
  Reve::Vector fP_pos;       // Momentum of positive track
  Reve::Vector fV_bach;      // Vertex of positive track
  Reve::Vector fP_bach;      // Momentum of positive track

  Reve::Vector fV_decay;     //decay point of the cascade
  Reve::Vector fV_birth;    // Reconstructed birth point of neutral particle

  vpPathMark_t         fPathMarksNeg;
  vpPathMark_t         fPathMarksPos;
  vpPathMark_t         fPathMarksBach;
  Reve::TrackRnrStyle *fRnrStyle;

  TPolyLine3D       fPolyLineNeg;
  TPolyLine3D       fPolyLinePos;
  TPolyLine3D       fPolyLineBach;
  TPolyLine3D       fPolyLineV0;   // line of V0 travel
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

public: 
  Cascade();
  Cascade(Reve::TrackRnrStyle* rs);
  virtual ~Cascade();

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
  void        SetRnrStyle(Reve::TrackRnrStyle* rs) { fRnrStyle = rs; }

  void  AddPathMarkPos(Reve::PathMark* pm) { fPathMarksPos.push_back(pm); }
  void  AddPathMarkNeg(Reve::PathMark* pm) { fPathMarksNeg.push_back(pm); }
  void  AddPathMarkBach(Reve::PathMark* pm) { fPathMarksBach.push_back(pm); }

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
  void MakeTrack(vpPathMark_t& pathMark, Reve::Vector& vtx,  Reve::Vector& p,
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

  ClassDef(Cascade, 1); // Visual representation of a cascade.
}; // endclass Cascade



//______________________________________________________________________

inline void Cascade::SetBeta(Float_t betaNeg, Float_t betaPos, Float_t betaBach) {
   fBeta_neg = betaNeg;
   fBeta_pos = betaPos;
   fBeta_bach = betaBach;
 }


//______________________________________________________________________

inline Float_t Cascade::GetV0P2() const {
  Float_t px = fP_neg.x + fP_pos.x, py = fP_neg.y + fP_pos.y,
    pz = fP_neg.z+fP_pos.z;
  return px*px + py*py + pz*pz;
}


inline Float_t Cascade::GetLambdaE() const {
  return sqrt(fgkMassLambda2+GetV0P2());
}

inline Float_t Cascade::GetXiE() const {
  Float_t e = GetLambdaE() +
    sqrt(fgkMassPion2 + fP_bach.x*fP_bach.x + fP_bach.y*fP_bach.y +
	 fP_bach.z*fP_bach.z);
  return e;
}

inline Float_t Cascade::GetOmegaE() const {
  Float_t e = GetLambdaE() +
    sqrt(fgkMassKaon2 + fP_bach.x*fP_bach.x + fP_bach.y*fP_bach.y +
	 fP_bach.z*fP_bach.z);
  return e;
}

inline Float_t Cascade::GetXiMass() const {
  Float_t e = GetXiE();
  return sqrt(e*e - GetP2());
}

inline Float_t Cascade::GetAntiXiMass() const { return GetXiMass();}

inline Float_t Cascade::GetOmegaMass() const {
  Float_t e = GetOmegaE();
  return sqrt(e*e - GetP2());
}

inline Float_t Cascade::GetAntiOmegaMass() const { return GetOmegaMass();}


//______________________________________________________________________

inline Float_t Cascade::GetDCA_v0_Bach() const {
  return fDCA_v0_Bach;
}

inline Float_t Cascade::GetCasCosPointingAngle() const {
return fCasCosPointingAngle;
}

inline Float_t Cascade::GetRadius() const {
  return sqrt(fV_birth.x*fV_birth.x + fV_birth.y*fV_birth.y);
}

//inline Float_t Cascade::GetDecayLength() const {
//return fDecayLength;
//}

inline Float_t Cascade::GetPseudoRapidity() const {
  Float_t theta = acos( GetPz()/GetMomentum() );
  return ( -log(tan(theta/2.)) );
}


//______________________________________________________________________
inline Float_t Cascade::GetPt2() const {
  Float_t px = GetPx(), py = GetPy();
  return (px*px+py*py);
}

inline Float_t Cascade::GetP2() const {

  Float_t px = GetPx(), py = GetPy(), pz = GetPz();
  return (px*px+py*py+pz*pz);
}

inline Float_t Cascade::GetPt() const {
  return sqrt(GetPt2());
}

inline Float_t Cascade::GetMomentum() const {
  return sqrt(GetP2());
}

inline Float_t Cascade::GetPx() const {
  return (fP_neg.x + fP_pos.x + fP_bach.x);
}

inline Float_t Cascade::GetPy() const {
  return (fP_neg.y + fP_pos.y + fP_bach.y);
}

inline Float_t Cascade::GetPz() const {
  return (fP_neg.z + fP_pos.z + fP_bach.z);
}

//______________________________________________________________________

inline Float_t Cascade::GetPosP2() const {
  return (fP_pos.x*fP_pos.x + fP_pos.y*fP_pos.y + fP_pos.z*fP_pos.z);
}

inline Float_t Cascade::GetPosP() const {
  return sqrt(GetPosP2());
}

inline Float_t Cascade::GetPosPt() const {
  return sqrt(fP_pos.x*fP_pos.x + fP_pos.y*fP_pos.y);
}

inline Float_t Cascade::GetPosPseudoRapidity() const {
  Float_t theta = acos( fP_pos.z/GetPosP() );
  return ( -log(tan(theta/2.)) );
}

//______________________________________________________________________
inline Float_t Cascade::GetNegP2() const {
  return (fP_neg.x*fP_neg.x + fP_neg.y*fP_neg.y + fP_neg.z*fP_neg.z);
}

inline Float_t Cascade::GetNegP() const {
  return sqrt(GetNegP2());
}

inline Float_t Cascade::GetNegPt() const {
  return sqrt(fP_neg.x*fP_neg.x + fP_neg.y*fP_neg.y);
}

inline Float_t Cascade::GetNegPseudoRapidity() const {
  Float_t theta = acos( fP_neg.z/GetNegP() );
  return ( -log(tan(theta/2.)) );
}


//______________________________________________________________________
inline Float_t Cascade::GetBachP2() const {
  return (fP_bach.x*fP_bach.x + fP_bach.y*fP_bach.y + fP_bach.z*fP_bach.z);
}

inline Float_t Cascade::GetBachP() const {
  return sqrt(GetBachP2());
}

inline Float_t Cascade::GetBachPt() const {
  return sqrt(fP_bach.x*fP_bach.x + fP_bach.y*fP_bach.y);
}

inline Float_t Cascade::GetBachPseudoRapidity() const {
  Float_t theta = acos( fP_bach.z/GetBachP() );
  return ( -log(tan(theta/2.)) );
}


/***********************************************************************
*
*  CascadeList class
*
************************************************************************/

class CascadeList : public Reve::RenderElementList
{
  CascadeList(const CascadeList&);            // Not implemented
  CascadeList& operator=(const CascadeList&); // Not implemented

private:
  void  Init();

protected:
  TString              fTitle;

  Reve::TrackRnrStyle *fRnrStyle;

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

public:
  CascadeList(Reve::TrackRnrStyle* rs=0);
  CascadeList(const Text_t* name, Reve::TrackRnrStyle* rs=0);

  virtual const Text_t* GetTitle() const { return fTitle; }
  virtual void SetTitle(const Text_t* t) { fTitle = t; }
  virtual void SetTracksColor(Color_t cNeg, Color_t cPos, Color_t cBach)
  { fNegColor = cNeg; fPosColor = cPos; fBachColor = cBach; }

  virtual Bool_t CanEditMainColor()  { return kTRUE; }

  virtual void Paint(Option_t* option="");

  void  SetRnrStyle(Reve::TrackRnrStyle* rst) { fRnrStyle= rst; }
  Reve::TrackRnrStyle* GetRnrStyle()          { return fRnrStyle; } 

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
  void UnFill(Cascade* cas);
  void Filter(Cascade* cas);
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

  //--------------------------------

  ClassDef(CascadeList, 1); // A list of Cascade objects.
};

} // namespace Alieve

#endif
