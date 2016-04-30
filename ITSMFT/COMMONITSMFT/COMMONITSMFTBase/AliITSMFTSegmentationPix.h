#ifndef AliITSMFTSegmentationPix_H
#define AliITSMFTSegmentationPix_H

#include "TObject.h"



// segmentation and response for pixels in ITS upgrade 

class AliITSMFTSegmentationPix : public TObject {
 public:
  AliITSMFTSegmentationPix(UInt_t id=0, int nchips=0,int ncol=0,int nrow=0,
			 float pitchX=0,float pitchZ=0,
			 float thickness=0,
			 float pitchLftC=-1,float pitchRgtC=-1,
			   float edgL=0,float edgR=0,float edgT=0,float edgB=0,float thr=0);
  
  //  AliITSMFTSegmentationPix(Option_t *opt="" );
  AliITSMFTSegmentationPix(const AliITSMFTSegmentationPix &source);
  virtual ~AliITSMFTSegmentationPix();
  AliITSMFTSegmentationPix& operator=(const AliITSMFTSegmentationPix &source);
  //
  virtual void    Init();
  //  
  virtual void    SetDetSize(Float_t p1,Float_t p2,Float_t p3) 
                    {fDx=p1; fDz=p2; fDy=p3;}
  
  // Detector length
  virtual Float_t Dx() const {return fDx;}
  // Detector width
  virtual Float_t Dz() const {return fDz;}
  // Detector thickness
  virtual Float_t Dy() const {return fDy;}

  // Cell size   
  virtual void    SetPadSize(Float_t,Float_t) {MayNotUse("SetPadSize");}

  virtual void    SetNPads(Int_t, Int_t) {MayNotUse("SetPadSize");}
  virtual Int_t   GetNPads() const {return fNCol*fNRow;}
  //
  virtual void    GetPadIxz(Float_t x,Float_t z,Int_t &ix,Int_t &iz) const;
  virtual void    GetPadCxz(Int_t ix,Int_t iz,Float_t &x,Float_t &z) const;
  virtual void    GetPadTxz(Float_t &x ,Float_t &z) const;
  virtual Bool_t  LocalToDet(Float_t x,Float_t z,Int_t &ix,Int_t &iz) const;
  virtual void    DetToLocal(Int_t ix,Int_t iz,Float_t &x,Float_t &z) const;
  virtual void    CellBoundries(Int_t ix,Int_t iz,Double_t &xl,Double_t &xu,Double_t &zl,Double_t &zu) const;
  //
  virtual Int_t    GetNumberOfChips() const {return fNChips;}
  virtual Int_t    GetMaximumChipIndex() const {return fNChips-1;}
  //
  virtual Int_t    GetChipFromLocal(Float_t, Float_t zloc) const;
  virtual Int_t    GetChipsInLocalWindow(Int_t* array, Float_t zmin, Float_t zmax, Float_t, Float_t) const;
  //
  virtual Int_t    GetChipFromChannel(Int_t, Int_t iz) const;
  //
  virtual Float_t  Dpx(Int_t ix=0) const;
  virtual Float_t  Dpz(Int_t iz)   const;
  Float_t          DxActive()      const {return fDxActive;}
  Float_t          DzActive()      const {return fDzActive;}
  Float_t          GetShiftXLoc()  const {return fShiftXLoc;}
  Float_t          GetShiftZLoc()  const {return fShiftZLoc;}
  Float_t          GetGuardLft()   const {return fGuardLft;}
  Float_t          GetGuardRgt()   const {return fGuardRgt;}
  Float_t          GetGuardTop()   const {return fGuardTop;}
  Float_t          GetGuardBot()   const {return fGuardBot;}
  Float_t          GetPitchX()     const {return fPitchX;}
  Float_t          GetPitchZ()     const {return fPitchZ;}
  //
  Int_t            GetNRow()       const {return fNRow;}
  Int_t            GetNCol()       const {return fNCol;}
  //
  virtual Int_t    Npx()           const {return GetNRow();}
  virtual Int_t    Npz()           const {return GetNCol();}
  //
  virtual void Neighbours(Int_t iX,Int_t iZ,Int_t* Nlist,Int_t Xlist[10],Int_t Zlist[10]) const;
  //
  virtual void Print(Option_t* option = "") const;
  //
  virtual Int_t                    GetChipTypeID()              const {return GetUniqueID();}
  //
  void         SetDiodShiftMatrix(Int_t nrow,Int_t ncol, const Float_t *shiftX, const Float_t *shiftZ);
  void         SetDiodShiftMatrix(Int_t nrow,Int_t ncol, const Double_t *shiftX, const Double_t *shiftZ);
  void         GetDiodShift(Int_t row,Int_t col, Float_t &dx,Float_t &dz) const;
  void         GetDiodShift(Int_t row,Int_t col, Double_t &dx,Double_t &dz) const {float dxf,dzf; GetDiodShift(row,col,dxf,dzf); dx=dxf; dz=dzf; }
  //
  Bool_t                           Store(const char* outf);
  static AliITSMFTSegmentationPix*   LoadWithID(UInt_t id, const char* inpf);
  static void                      LoadSegmentations(TObjArray* dest, const char* inpf);
  //
  static UInt_t      ComposeChipTypeID(UInt_t segmId);  
    
    
 protected:
  Float_t Z2Col(Float_t z) const;
  Float_t Col2Z(Int_t col) const;
  //
 protected:
    Float_t  fDx;    //Full width of the detector (x axis)- microns
    Float_t  fDz;    //Full length of the detector (z axis)- microns
    Float_t  fDy;    //Full thickness of the detector (y axis) -um 

    Float_t  fGuardLft;        // left guard edge
    Float_t  fGuardRgt;        // right guard edge
    Float_t  fGuardTop;        // upper guard edge
    Float_t  fGuardBot;        // bottom guard edge
    Float_t  fShiftXLoc;       // shift in local X of sensitive area wrt geometry center
    Float_t  fShiftZLoc;       // shift in local Z of sensitive area wrt geometry center
    Float_t  fDxActive;        // size of active area in X
    Float_t  fDzActive;        // size of active area in Z    
    Float_t  fPitchX;          // default pitch in X
    Float_t  fPitchZ;          // default pitch in Z
    Float_t  fPitchZLftCol;    // Z pitch of left column of each chip
    Float_t  fPitchZRgtCol;    // Z pitch of right column of each chip
    Float_t  fChipDZ;          // aux: chip size along Z
    Int_t    fNChips;          // number of chips per chip
    Int_t    fNColPerChip;     // number of columns per chip
    Int_t    fNRow;            // number of rows
    Int_t    fNCol;            // number of columns (total)
    //
    Int_t    fDiodShiftMatNCol; // periodicity of diod shift in columns
    Int_t    fDiodShiftMatNRow; // periodicity of diod shift in rows
    Int_t    fDiodShiftMatDim;  // dimension of diod shift matrix
    Float_t* fDiodShidtMatX; //[fDiodShiftMatDim] diod shift in X (along column), in fraction of X pitch
    Float_t* fDiodShidtMatZ; //[fDiodShiftMatDim] diod shift in Z (along row), in fraction of Z pitch
    //
    static const char* fgkSegmListName; // pattern for segmentations list name
    //
  ClassDef(AliITSMFTSegmentationPix,4) //Segmentation class upgrade pixels 

};

#endif
