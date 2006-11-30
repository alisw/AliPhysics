// @(#) $Id$

#ifndef ALIL3ConfMapPointH
#define ALIL3ConfMapPointH

#include "AliHLTRootTypes.h"

class AliHLTSpacePointData;
class AliHLTConfMapTrack;
class AliHLTVertex;

class AliHLTConfMapPoint {

 private:

  Int_t fHitNumber;     //hit number
  Int_t fTrackNumber;   //track number
  Int_t fNextHitNumber; //next hit number
  Bool_t fUsed;         //flag is used
  Int_t fPadrow;        //padrow
  Int_t fSector;        //sector

  //global coordinates and their errors
  Double_t fx;    //glob x
  Double_t fy;    //glob y
  Double_t fz;    //glob z
  Double_t fxerr; //glob xerr
  Double_t fyerr; //glob yerr
  Double_t fzerr; //glob zerr

  Double_t fWxy;  // x-y weight on x-y
  Double_t fWz;   // z weight on z
  Float_t fs;      //track trajectory
  
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
  
  AliHLTConfMapPoint *fNextVolumeHit; //!
  AliHLTConfMapPoint *fNextRowHit;    //!
  AliHLTConfMapPoint *fNextTrackHit;  //! Linked chain of points in a track
  Short_t fPhiIndex; //phi index
  Short_t fEtaIndex; //eta index
  Double_t fXYChi2; //xy chi
  Double_t fSZChi2; //z chi
  Int_t fMCTrackID[3]; //MClabel of tracks, may overlap

  static Bool_t fgDontMap; //flag to switch off mapping  

 public:

  AliHLTConfMapPoint();
  virtual ~AliHLTConfMapPoint();
  
  void Reset();
  Bool_t ReadHits(AliHLTSpacePointData* hits );
  
   // getter
  Double_t GetX() const {return fx;}
  Double_t GetY() const {return fy;}
  Double_t GetZ() const {return fz;}
  Double_t GetXerr() const {return fxerr;}
  Double_t GetYerr() const {return fyerr;}
  Double_t GetZerr() const {return fzerr;}
  Int_t GetPadRow() const {return fPadrow;}
  Int_t GetSector() const {return fSector;}
  
  Double_t GetXYWeight() const {return fWxy;}
  Double_t GetZWeight() const {return fWz;}
  Float_t GetS()        const {return fs;}

  Bool_t GetUsage() const {return fUsed;}
  Double_t GetPhi() const {return fPhi;}
  Double_t GetEta() const {return fEta;}
  
  Double_t GetXprime() const    {return fXprime;}
  Double_t GetYprime() const    {return fYprime;}
  Double_t GetXprimeerr() const {return fXprimeerr;}
  Double_t GetYprimeerr() const {return fYprimeerr;}
  
  Double_t GetXt() const {return fXt;}
  Double_t GetYt() const {return fYt;}
  Double_t GetZt() const {return fZt;}
  Double_t GetXterr() const {return fXterr;}
  Double_t GetYterr() const {return fYterr;}
  Double_t GetZterr() const {return fZterr;}
  
  Double_t GetXv() const {return fXv;}
  Double_t GetYv() const {return fYv;}
  Double_t GetZv() const {return fZv;}
  Double_t GetXverr() const {return fXverr;}
  Double_t GetYverr() const {return fYverr;}
  Double_t GetZverr() const {return fZverr;}

  Int_t GetHitNumber() const {return fHitNumber;}
  Int_t GetNextHitNumber() const {return fNextHitNumber;}
  Int_t GetTrackNumber() const {return fTrackNumber;}
  //Int_t const *GetMCTrackID()     const {return fMCTrackID;}
  
  AliHLTConfMapPoint* GetNextVolumeHit(){return fNextVolumeHit;}
  AliHLTConfMapPoint* GetNextRowHit(){return fNextRowHit;}
  AliHLTConfMapPoint* GetNextTrackHit(){return fNextTrackHit;}
  Short_t GetPhiIndex() const {return fPhiIndex;}
  Short_t GetEtaIndex() const {return fEtaIndex;}
  Double_t GetXYChi2() const {return fXYChi2;}
  Double_t GetSZChi2() const {return fSZChi2;}
  //Int_t fMCTrackID[3]; //MClabel of tracks, may overlap

  // setter
  void SetNextVolumeHit(AliHLTConfMapPoint* p){fNextVolumeHit=p;}
  void SetNextRowHit(AliHLTConfMapPoint* p){fNextRowHit=p;}
  void SetNextTrackHit(AliHLTConfMapPoint* p){fNextTrackHit=p;}

  void SetPhiIndex(Short_t p){fPhiIndex=p;}
  void SetEtaIndex(Short_t p){fEtaIndex=p;}
  void SetXYChi2(Double_t d) {fXYChi2=d;}
  void SetSZChi2(Double_t d) {fSZChi2=d;}

  static void SetDontMap(Bool_t b){fgDontMap=b;}

  void SetX(Double_t f){fx=f;}
  void SetY(Double_t f){fy=f;}
  void SetZ(Double_t f){fz=f;}
  void SetXerr(Double_t f){fxerr=f;}
  void SetYerr(Double_t f){fyerr=f;}
  void SetZerr(Double_t f){fzerr=f;}
  void SetPadRow(Int_t f){fPadrow=f;}
  void SetSector(Int_t f){fSector=f;}
  void SetMCTrackID(Int_t f,Int_t g,Int_t h){fMCTrackID[0] = f; fMCTrackID[1]=g; fMCTrackID[2]=h;}

  void SetXYWeight(Float_t f){fWxy = f;}
  void SetZWeight(Float_t f){fWz = f;}
  void SetS(Float_t f){fs = f;}
  void SetUsage(Bool_t f){fUsed=f;}
  void SetPhi(Double_t f ){fPhi = f;}
  void SetEta(Double_t f){fEta = f;}
  void SetXprime(Double_t f){fXprime = f;}
  void SetYprime(Double_t f){fYprime = f;}
  void SetXprimeerr(Double_t f){fXprimeerr = f;}
  void SetYprimeerr(Double_t f){fYprimeerr = f;}
  void SetXt(Double_t f){fXt = f;}
  void SetYt(Double_t f){fYt = f;}
  void SetZt(Double_t f){fZt = f;}
  void SetXterr(Double_t f){fXterr = f;}
  void SetYterr(Double_t f){fYterr = f;}
  void SetZterr(Double_t f){fZterr = f;}
  void SetXv(Double_t f){fXv = f;}
  void SetYv(Double_t f){fYv = f;}
  void SetZv(Double_t f){fZv = f;}
  void SetXverr(Double_t f){fXverr = f;}
  void SetYverr(Double_t f){fYverr = f;}
  void SetZverr(Double_t f){fZverr = f;}
  void SetHitNumber(Int_t f){fHitNumber=f;}
  void SetTrackNumber(Int_t f){fTrackNumber=f;}
  void SetNextHitNumber(Int_t f){fNextHitNumber=f;}

  void Setup(AliHLTVertex *vertex);// does the usual setup in the right order
  void SetAngles();               // calculate spherical angles and set values
  void SetIntPoint(Double_t inx = 0., Double_t iny = 0.,
	           Double_t inz = 0., Double_t inxerr = 0., 
		   Double_t inyerr = 0., Double_t inzerr = 0.);  
  //-> set interaction point
  void SetShiftedCoord();// set shifted coordinates  
  void SetAllCoord(const AliHLTConfMapPoint *hit);// set conformal mapping coordinates in respect to given hit
  void SetConfCoord();// conformal mapping

  ClassDef(AliHLTConfMapPoint, 1)   //Conformal mapping hit class.
};

typedef AliHLTConfMapPoint AliL3ConfMapPoint; // for backward compatibility

#endif
