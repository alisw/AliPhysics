#ifndef ALIESDMUONTRACK_H
#define ALIESDMUONTRACK_H

#include "TObject.h"

class AliESDMuonTrack : public TObject {
public:
 AliESDMuonTrack(){
    // Constructor;
 } // Constructor
 virtual ~AliESDMuonTrack(){
    // Destructor;
 } // Destructor
 AliESDMuonTrack(const AliESDMuonTrack& );
 AliESDMuonTrack& operator=(const AliESDMuonTrack& );


 // Get and Set methods for data
  Double_t GetInverseBendingMomentum(void) const {return fInverseBendingMomentum;}
  void SetInverseBendingMomentum(Double_t InverseBendingMomentum) 
    {fInverseBendingMomentum = InverseBendingMomentum;}
  Double_t GetThetaX(void) const {return fThetaX;}
  void SetThetaX(Double_t ThetaX) {fThetaX = ThetaX;}
  Double_t GetThetaY(void) const {return fThetaY;}
  void SetThetaY(Double_t ThetaY) {fThetaY = ThetaY;}
  Double_t GetZ(void) const {return fZ;}
  void SetZ(Double_t Z) {fZ = Z;}
  Double_t GetBendingCoor(void) const {return fBendingCoor;}
  void SetBendingCoor(Double_t BendingCoor) {fBendingCoor = BendingCoor;}
  Double_t GetNonBendingCoor(void) const {return fNonBendingCoor;}
  void SetNonBendingCoor(Double_t NonBendingCoor) {fNonBendingCoor = NonBendingCoor;}
  Double_t GetChi2(void) const {return fChi2;}
  void SetChi2(Double_t Chi2) {fChi2 = Chi2;}
  UInt_t GetNHit(void) const {return fNHit;}
  void SetNHit(UInt_t NHit) {fNHit = NHit;}

  Float_t GetX11() const {return fX11;}
  void SetX11(Float_t X11) {fX11 = X11;}
  Float_t GetY11() const {return fY11;}
  void SetY11(Float_t Y11) {fY11 = Y11;}
  Float_t GetThetaX11() const {return fThetaX11;}
  void SetThetaX11(Float_t ThetaX) {fThetaX11 = ThetaX;}
  Float_t GetThetaY11() const {return fThetaY11;}    
  void SetThetaY11(Float_t ThetaY) {fThetaY11 = ThetaY;}

protected:
  // tracking chamber
  Double_t fInverseBendingMomentum; // Inverse bending momentum (GeV/c ** -1) times the charge 
  Double_t fThetaX;           // Angle of track at vertex in X direction (rad)
  Double_t fThetaY;           // Angle of track at vertex in Y direction (rad)
  Double_t fZ;                // Z coordinate (cm)
  Double_t fBendingCoor;      // bending coordinate (cm)
  Double_t fNonBendingCoor;   // non bending coordinate (cm)
  Double_t fChi2;             // chi2 in the MUON track fit
  UInt_t   fNHit;              // number of hit in the track

  // trigger chamber
  Float_t fX11;    // x position of fired Y strip in MC11
  Float_t fY11;    // y position of fired X strip in MC11
  Float_t fThetaX11; // track theta angle in X   
  Float_t fThetaY11; // track theta angle in Y

  ClassDef(AliESDMuonTrack,1)  //MUON ESD track class 
};

#endif 
