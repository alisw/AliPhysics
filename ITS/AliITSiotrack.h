#ifndef ALIITSIOTRACK_H
#define ALIITSIOTRACK_H 

////////////////////////////////////////////////////
//  Reconstructed space point class for set:ITS   //
////////////////////////////////////////////////////

#include <TObject.h>

class AliITSiotrack : public TObject {


 public:

  AliITSiotrack();                        // constructor
  virtual ~AliITSiotrack() {};            // distructor
  Int_t   GetLabel() const {return fLab;} // get track label
  Int_t   GetIdPoint(Int_t i) const {return fIdPoints[i];}   // get the identification number for the point
  Int_t   GetIdModule(Int_t i) const {return fIdModules[i];} // get the module number for the point
  //Double_t * GetStateVector();            // get the state vector (5 elements)
  Double_t * GetCovMatrix();              // get the covariance matrix as 15 elements array
  Double_t GetStatePhi() const {return fStateVPhi;}  // gets the Phi angle of the state vector
  Double_t GetStateZ() const {return fStateVZ;}      // gets the Z cohordinate of the state vector
  Double_t GetStateD() const {return fStateVD;}      // gets the radial impact parameter of the state vector
  Double_t GetStateC() const {return fStateVC;}      // gets the curvature C of the state vector
  Double_t GetStateTgl() const {return fStateVTgl;}  // gets the dip angle tangent of the state vector
  Double_t GetRadius() const {return fRadius;}       // gets the radius corresponding to the state vector
  Int_t     GetCharge() const {return fCharge;}      // gets the particle charge
  Float_t GetX() const {return fX;}     // gets the x cohordinate of the found vertex
  Float_t GetZ() const {return fZ;}     // gets the z cohordinate of the found vertex
  Float_t GetY() const {return fY;}     // gets the y cohordinate of the found vertex
  Float_t GetPx() const {return fPx;}   // gets the x momentum component at the found vertex 
  Float_t GetPy() const {return fPy;}   // gets the y momentum component at the found vertex 
  Float_t GetPz()const {return fPz;}    // gets the z momentum component at the found vertex 
 
  void SetCovMatrix(TMatrix *cov);      // sets the covariance matrix
  void SetLabel(Int_t lab) {fLab=lab;}  // sets the track label
  void SetIdPoint(Int_t i,Int_t pnt) {fIdPoints[i]=pnt;}   // set the identification number for the point
  void SetIdModule(Int_t i,Int_t mod) {fIdModules[i]=mod;} // set the module number for the point
   
  void SetStatePhi(Double_t phi) {fStateVPhi=phi;}   // sets the Phi angle of the state vector 
  void SetStateZ(Double_t z) {fStateVZ=z;}           // sets the Z cohordinate of the state vector
  void SetStateD(Double_t d) {fStateVD=d;}           // sets the radial impact parameter of the state vector 
  void SetStateTgl(Double_t tgl) {fStateVTgl=tgl;}   // sets the dip angle tangent of the state vector
  void SetStateC(Double_t c) {fStateVC=c;}           // sets the curvature C of the state vector
  void SetRadius(Double_t r) {fRadius= r;}           // sets the radius corresponding to the state vector
  void SetCharge(Int_t charge) {fCharge=charge;}     // sets the particle charge

  void SetX(Float_t x){fX=x;}        // sets the x cohordinate of the found vertex
  void SetZ(Float_t z){fZ=z;}        // sets the z cohordinate of the found vertex
  void SetY(Float_t y){fY=y;}        // sets the y cohordinate of the found vertex
  void SetPx(Float_t px) {fPx=px;}   // sets the x momentum component at the found vertex 
  void SetPy(Float_t py) {fPy=py;}   // sets the y momentum component at the found vertex
  void SetPz(Float_t pz) {fPz=pz;}   // sets the z momentum component at the found vertex

 private:
    
  Int_t     fLab;       // label of reconstructed track
  Float_t   fX ;        // x cohordinate of the found vertex
  Float_t   fY ;        // y cohordinate of the found vertex
  Float_t   fZ ;        // z cohordinate of the found vertex
  Float_t   fPx;        // x component of track momentum at the found vertex
  Float_t   fPy;        // y component of track momentum at the found vertex
  Float_t   fPz;        // z component of track momentum at the found vertex
   
  //
  Int_t     fIdPoints[6];   // points assigned to the track (entry # in fRecPoints is given by module #)
  Int_t     fIdModules[6];  // module # corresponding to each point belonging to the track

  Double_t  fStateVPhi;          // 
  Double_t  fStateVZ;            //  state vector components
  Double_t  fStateVD;            //
  Double_t  fStateVTgl;          //
  Double_t  fStateVC;            //
	 
  Double_t  fRadius;             //  distance of the point from the origin
  Int_t     fCharge;             //  particle charge  

//  Covariance matrix	 
  Double_t  fC00;
  Double_t  fC10, fC11;
  Double_t  fC20, fC21, fC22;
  Double_t  fC30, fC31, fC32, fC33;
  Double_t  fC40, fC41, fC42, fC43, fC44;
	 
  ClassDef(AliITSiotrack,1)  // AliITSiotrack class
};

#endif



