#ifndef ALIITSSEGMENTATIONSPD_H
#define ALIITSSEGMENTATIONSPD_H

#include "AliITSsegmentation.h"


class AliITSgeom;

// segmentation and response for SPD 

class AliITSsegmentationSPD :
public AliITSsegmentation {
 public:

    AliITSsegmentationSPD();
    AliITSsegmentationSPD(AliITSgeom *gm);
    AliITSsegmentationSPD(const AliITSsegmentationSPD &source);
    virtual ~AliITSsegmentationSPD(){}
    AliITSsegmentationSPD& operator=(const AliITSsegmentationSPD &source);

    // Set Detector Segmentation Parameters

    // Maximum number of pixels along the two coordinates  
    virtual void    SetNPads(Int_t p1, Int_t p2);
    // Returns the maximum number of cells (digits) posible
    virtual Int_t   GetNPads() const {return fNpx*fNpz;}
    // Set Pixel Size Array in x and z, microns.
    virtual void    SetBinSize(Float_t *x,Float_t *z);

    // Transform from real to pixel coordinates
    virtual void    GetPadIxz(Float_t x,Float_t z,Int_t &ix,Int_t &iz) const;
    // Transform from pixel to real coordinates
    virtual void    GetPadCxz(Int_t ix,Int_t iz,Float_t &x,Float_t &z) const;
    // Transform from real local to global coordinates
    //virtual void    GetGlobal(Int_t module,Float_t *l ,Float_t *g) {}
    // Local transformation of real local coordinates -
    virtual void    GetPadTxz(Float_t &x ,Float_t &z) const;
    // Transformation from Geant cm detector center local coordinates
    // to detector segmentation/cell coordiantes starting from (0,0).
    virtual Bool_t  LocalToDet(Float_t x,Float_t z,Int_t &ix,Int_t &iz) const;
    // Transformation from detector segmentation/cell coordiantes starting
    // from (0,0) to Geant cm detector center local coordinates.
    virtual void    DetToLocal(Int_t ix,Int_t iz,Float_t &x,Float_t &z) const;
    // Returns the Cell upper and lower boundries in x and y. cell indexes
    // starting from (0,0) and return Geant cm detector centered local
    // coordinates, consistant with DetToLocal and LocalToDet functions above.
    virtual void CellBoundries(Int_t ix,Int_t iz,Double_t &xl,Double_t &xu,
			       Double_t &zl,Double_t &zu) const;
    //
    // Initialisation
    virtual void Init();
    virtual void Init300();
    //
    // Get member data
    //
    // Pixel size in x
    virtual Float_t Dpx(Int_t ix) const;
    // Pixel size in z 
    virtual Float_t Dpz(Int_t iz) const;

    // Maximum number of Pixels in x
    virtual Int_t    Npx() const {return fNpx;}
    // Maximum number of Pixels in z
    virtual Int_t    Npz() const {return fNpz;}
    //
    // Get next neighbours
    virtual void Neighbours
       (Int_t iX,Int_t iZ,Int_t* Nlist,Int_t Xlist[10],Int_t Zlist[10]) const;
    // Print default parameters (static const data members, if any)
    virtual void PrintDefaultParameters() const 
            {Warning("PrintDefaultParameters","No def. parameters defined as const static data members\n");}
    
 protected:

    virtual void Copy(TObject &obj) const;
    Int_t   fNpx;           // Number of pixels in x
    Int_t   fNpz;           // Number of pixels in z
    Float_t fCellSizeX[256];// Size for each pixel in x -microns
    Float_t fCellSizeZ[280];// Size for each pixel in z -microns

 private:

    Float_t ColFromZ300(Float_t z) const;
    Float_t ZFromCol300(Int_t col) const;
    Float_t ZpitchFromCol300(Int_t col) const;
    Float_t ColFromZ(Float_t z) const;
    Float_t ZFromCol(Int_t col) const;
    Float_t ZpitchFromCol(Int_t col) const;

  ClassDef(AliITSsegmentationSPD,2) //Segmentation class for SPD 

};

#endif
