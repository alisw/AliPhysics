#ifndef ALIITSUCLUSTERPIX_H
#define ALIITSUCLUSTERPIX_H

#include "AliCluster.h"

class TGeoHMatrix;
class AliITSUGeomTGeo;


class AliITSUClusterPix : public AliCluster
{
 public:
  enum { // frame in which the track is currently defined
    kFrameLoc  = BIT(16)
    ,kFrameTrk = BIT(17)
    ,kFrameGlo = BIT(18)
    ,kFrameBits = kFrameLoc|kFrameTrk|kFrameGlo
    ,kSplit    = BIT(19)
  };
  //
  enum SortMode_t { // various modes
    kSortIdLocXZ  = BIT(0)    // sort according to ID, then X,Z of local frame
    ,kSortIdTrkYZ = BIT(1)    // sort according to ID, then Y,Z of tracking frame
    ,kSortBits = kSortIdLocXZ|kSortIdTrkYZ
  };

 public:
  AliITSUClusterPix();
  AliITSUClusterPix(const AliITSUClusterPix& cluster);
  AliITSUClusterPix &operator=(const AliITSUClusterPix& cluster);
  virtual ~AliITSUClusterPix();
  //
  Bool_t  IsFrameLoc()         const {return TestBit(kFrameLoc);}
  Bool_t  IsFrameGlo()         const {return TestBit(kFrameGlo);}
  Bool_t  IsFrameTrk()         const {return TestBit(kFrameTrk);}
  //
  Bool_t  IsSplit()            const {return TestBit(kSplit);}
  //
  void    SetFrameLoc()              {ResetBit(kFrameBits); SetBit(kFrameLoc);}
  void    SetFrameGlo()              {ResetBit(kFrameBits); SetBit(kFrameGlo);}
  void    SetFrameTrk()              {ResetBit(kFrameTrk);  SetBit(kFrameTrk);}
  //
  void    SetSplit(Bool_t v=kTRUE)   {SetBit(kSplit,v);}
  //
  void    GoToFrameGlo();
  void    GoToFrameLoc();
  void    GoToFrameTrk();
  void    GetLocalXYZ(Float_t xyz[3])                       const;
  void    GetTrackingXYZ(Float_t xyz[3])                    const; 
  //
  void    SetNxNzN(UChar_t nx,UChar_t nz,UShort_t n) {fNxNzN = (n<<16) + (nx<<8) + nz;}
  Int_t   GetNx()                                           const {return (fNxNzN>>8)&0xff;}
  Int_t   GetNz()                                           const {return fNxNzN&0xff;}
  Int_t   GetNPix()                                         const {return fNxNzN>>16;}
  //
  void    SetQ(UShort_t q)                                        {fCharge = q;}
  Int_t   GetQ()                                            const {return fCharge;}
  //
  virtual void                 Print(Option_t* option = "") const;
  virtual const TGeoHMatrix*   GetTracking2LocalMatrix()           const;
  virtual TGeoHMatrix*         GetMatrix(Bool_t original = kFALSE) const;
  virtual Bool_t               GetGlobalXYZ(Float_t xyz[3]) const;
  virtual Bool_t               GetGlobalCov(Float_t cov[6]) const;
  virtual Bool_t               GetXRefPlane(Float_t &xref)  const;
  //
  virtual Bool_t               IsSortable()                 const {return kTRUE;}
  virtual Bool_t               IsEqual(const TObject* obj)  const;
  virtual Int_t	               Compare(const TObject* obj)  const;
  //
  Bool_t  HasCommonTrack(const AliCluster* cl)          const;
  //
  static  void                 SetGeom(AliITSUGeomTGeo* gm) {fgGeom = gm;}
  static  void                 SetSortMode(SortMode_t md)   {fgMode &= ~kSortBits; fgMode |= md;}
  static  UInt_t               GetSortMode()                {return fgMode|kSortBits;}
  static  UInt_t               GetMode()                    {return fgMode;}
  static  SortMode_t           SortModeIdTrkYZ()            {return kSortIdTrkYZ;}
  static  SortMode_t           SortModeIdLocXZ()            {return kSortIdLocXZ;}
  //
 protected:
  //
  UShort_t                fCharge;        //  charge (for MC studies only)
  Int_t                   fNxNzN;         //  effective cluster size in X (1st byte) and Z (2nd byte) directions 
                                          //  and total Npix(3d&4th bytes)
  static UInt_t           fgMode;         //! general mode (sorting mode etc)
  static AliITSUGeomTGeo* fgGeom;         //! pointer on the geometry data

  ClassDef(AliITSUClusterPix,1)
};

#endif
