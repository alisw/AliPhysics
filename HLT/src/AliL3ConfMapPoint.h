#ifndef ALIL3_ConfMapPoint
#define ALIL3_ConfMapPoint

#include "AliL3RootTypes.h"

class AliL3SpacePointData;
class AliL3ConfMapTrack;
class AliL3Vertex;

class AliL3ConfMapPoint {

 private:

  Int_t fHitNumber;
  Int_t fTrackNumber;
  Int_t fNextHitNumber;
  Bool_t fUsed;
  
  Int_t fPadrow;
  Int_t fSector;

  //global coordinates and their errors
  Double_t x;
  Double_t y;
  Double_t z;
  Double_t xerr;
  Double_t yerr;
  Double_t zerr;

  Double_t fWxy;  // x-y weight on x-y
  Double_t fWz;   // z weight on z
  Float_t s;     //track trajectory
  
   // Interaction point
  Double_t   fXt;          // x-value of the interaction point
  Double_t   fYt;          // y-value of the interaction point
  Double_t   fZt;          // z-value of the interaction point
  
  Double_t   fXterr;       // error of mXt
  Double_t   fYterr;       // error of mYt
  Double_t   fZterr;       // error of mZt
  
  // conformal mapping coordinates
  Double_t   fXprime;      // transformed x
  Double_t   fYprime;      // transformed y  
  
  Double_t   fXprimeerr;   // error of mXprime
  Double_t   fYprimeerr;   // error of mYprime
  
  // coordinates with respect to the vertex
  
  // cartesian coordinates
  Double_t   fXv;          // x with respect to vertex
  Double_t   fYv;          // y with respect to vertex
  Double_t   fZv;          // z with respect to vertex
  
  Double_t   fXverr;       // error of mXv
  Double_t   fYverr;       // error of mYv
  Double_t   fZverr;       // error of mZv
  
  // spherical coordinates
  Double_t   fPhi;         // angle phi
  Double_t   fEta;         // pseudorapidity
  
  
 public:

  AliL3ConfMapPoint();
  virtual ~AliL3ConfMapPoint();

  Bool_t ReadHits(AliL3SpacePointData* hits ); //!
 
  AliL3ConfMapPoint *nextVolumeHit; //!
  AliL3ConfMapPoint *nextRowHit;  //!
  
  AliL3ConfMapPoint *nextTrackHit; //! Linked chain of points in a track
  Short_t phiIndex;
  Short_t etaIndex;
 
  Double_t xyChi2;
  Double_t szChi2;
  Int_t fMCTrackID[3]; //MClabel of tracks, may overlap

   // getter
  Double_t GetX() const {return x;}
  Double_t GetY() const {return y;}
  Double_t GetZ() const {return z;}
  Double_t GetXerr() const {return xerr;}
  Double_t GetYerr() const {return yerr;}
  Double_t GetZerr() const {return zerr;}
  Int_t GetPadRow() const {return fPadrow;}
  Int_t GetSector() const {return fSector;}
  
  Double_t GetXYWeight() const {return fWxy;}
  Double_t GetZWeight() const {return fWz;}
  Float_t GetS()        const {return s;}

  //AliL3ConfMapTrack *GetTrack(TClonesArray *tracks) const;
  
  Bool_t GetUsage() const {return fUsed;}
  Double_t   GetPhi() const          { return fPhi;        }
  Double_t   GetEta() const          { return fEta;        }
  
  Double_t   GetXprime() const       { return fXprime;     }
  Double_t   GetYprime() const       { return fYprime;     }
  Double_t   GetXprimeerr() const    { return fXprimeerr;  }
  Double_t   GetYprimeerr() const    { return fYprimeerr;  }
  
  Double_t   GetXt() const           { return fXt;         }
  Double_t   GetYt() const           { return fYt;         }
  Double_t   GetZt() const           { return fZt;         }
  Double_t   GetXterr() const        { return fXterr;      }
  Double_t   GetYterr() const        { return fYterr;      }
  Double_t   GetZterr() const        { return fZterr;      }
  
  Double_t   GetXv() const           { return fXv;         }
  Double_t   GetYv() const           { return fYv;         }
  Double_t   GetZv() const           { return fZv;         }
  Double_t   GetXverr() const        { return fXverr;      }
  Double_t   GetYverr() const        { return fYverr;      }
  Double_t   GetZverr() const        { return fZverr;      }

  Int_t GetHitNumber() const {return fHitNumber;}
  Int_t GetNextHitNumber() const {return fNextHitNumber;}
  Int_t GetTrackNumber() const {return fTrackNumber;}
  //  Int_t const *GetMCTrackID()     const {return fMCTrackID;}

  // setter
  void SetX(Double_t f) {x=f;}
  void SetY(Double_t f) {y=f;}
  void SetZ(Double_t f) {z=f;}
  void SetXerr(Double_t f) {xerr=f;}
  void SetYerr(Double_t f) {yerr=f;}
  void SetZerr(Double_t f) {zerr=f;}
  void SetPadRow(Int_t f) {fPadrow=f;}
  void SetSector(Int_t f) {fSector=f;}
  void SetMCTrackID(Int_t f,Int_t g,Int_t h) {fMCTrackID[0] = f; fMCTrackID[1]=g; fMCTrackID[2]=h;}

  void SetXYWeight(Float_t f) {fWxy = f;}
  void SetZWeight(Float_t f) {fWz = f;}
  void SetS(Float_t f) {s = f;}

  void SetUsage(Bool_t f) {fUsed=f;}
    
  void    SetPhi(Double_t f)         {           fPhi = f; }
  void    SetEta(Double_t f)         {           fEta = f; }
  
  void    SetXprime(Double_t f)      {        fXprime = f; }
  void    SetYprime(Double_t f)      {        fYprime = f; }
  void    SetXprimeerr(Double_t f)   {     fXprimeerr = f; }
  void    SetYprimeerr(Double_t f)   {     fYprimeerr = f; }
  
  void    SetXt(Double_t f)          {            fXt = f; }
  void    SetYt(Double_t f)          {            fYt = f; }
  void    SetZt(Double_t f)          {            fZt = f; }
  void    SetXterr(Double_t f)       {         fXterr = f; }
  void    SetYterr(Double_t f)       {         fYterr = f; }
  void    SetZterr(Double_t f)       {         fZterr = f; }
  
  void    SetXv(Double_t f)          {            fXv = f; }
  void    SetYv(Double_t f)          {            fYv = f; }
  void    SetZv(Double_t f)          {            fZv = f; }
  void    SetXverr(Double_t f)       {         fXverr = f; }
  void    SetYverr(Double_t f)       {         fYverr = f; }
  void    SetZverr(Double_t f)       {         fZverr = f; }
  
  void    SetHitNumber(Int_t f) {fHitNumber=f;}
  void    SetTrackNumber(Int_t f) {fTrackNumber=f;}
  void    SetNextHitNumber(Int_t f) {fNextHitNumber=f;}

  void    Setup(AliL3Vertex *vertex);// does the usual setup in the right order
  void    SetAngles();// calculate spherical angles and set values
  void    SetIntPoint(const Double_t in_x = 0.,const Double_t in_y = 0.,
		      const Double_t in_z = 0.,const Double_t in_x_err = 0., 
		      const Double_t in_y_err = 0., const Double_t in_z_err = 0.);  
  //-> set interaction point
  void    SetShiftedCoord();// set shifted coordinates  
  void    SetAllCoord(const AliL3ConfMapPoint *hit);// set conformal mapping coordinates in respect to given hit
  void    SetConfCoord();// conformal mapping
  
//  ClassDef(AliL3ConfMapPoint, 1)   //Conformal mapping hit class.
};

#endif
