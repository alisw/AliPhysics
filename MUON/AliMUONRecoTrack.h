#ifndef ALI_MUON_RECO_TRACK_H
#define ALI_MUON_RECO_TRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup rec
/// \class AliMUONRecoTrack
/// \brief A reconstructed muon track
///
/// \author M.Gheata, A.Gheata 09/10/00

#include <TObject.h>
#include <TMath.h>

class AliMUONRecoTrack : public TObject 
{
  public:
    AliMUONRecoTrack();
    AliMUONRecoTrack(Bool_t active);
    virtual        ~AliMUONRecoTrack() { }	//desctructor
    Double_t GetChi2r() const {return fChi2r;};
    Double_t GetMomReconstr(Int_t axis) const {return fPr[axis];};
    Int_t    GetSign() const {return fSign;};
    Double_t GetPosX(Int_t chamber) const {return fPosX[chamber];};
    Double_t GetPosY(Int_t chamber) const {return fPosY[chamber];};
    Double_t GetPosZ(Int_t chamber) const {return fPosZ[chamber];};
    Double_t GetVertexPos() const { return fZvr;};
    Double_t P() const {return TMath::Sqrt(fPr[0]*fPr[0] + fPr[1]*fPr[1] + fPr[2]*fPr[2]);};
    Double_t Phi() const;
    void           SetChi2r(Double_t chi) { fChi2r = chi;};
    void           SetHitPosition(Int_t chamber, Double_t x, Double_t y, Double_t z);
    void           SetMomReconstr(Double_t px, Double_t py, Double_t pz);
    void           SetSign(Int_t sign) {fSign = sign;};
    void           SetVertexPos(Double_t zvr) {fZvr = zvr;};
    void           SetFlag(Int_t flag)  {fFlag = flag;};

    Double_t Theta() const;
    void TrackInfo() const;

  private:
    Int_t       fSign;                  ///< charge sign
    Int_t       fFlag;                  ///<  flag of reconstructed track (0-"good", >0-"bad") 
    Double_t 	fZvr;                   ///< z of track vertex point
    Double_t 	fChi2r;	                ///< chi squared for reco. track
    Double_t 	fPr[3];	                // reconstr. momentum (same as in vertex)
    Double_t 	fPosX[10];              ///< hit X position in all chambers
    Double_t 	fPosY[10];              ///< hit Y position in all chambers    
    Double_t 	fPosZ[10];              ///< hit Z position in all chambers

  ClassDef(AliMUONRecoTrack,1)	// A reconstructed muon track
};

#endif
