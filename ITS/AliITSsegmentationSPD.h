#ifndef ALIITSSEGMENTATIONSPD_H
#define ALIITSSEGMENTATIONSPD_H

#include "AliITSsegmentation.h"

#include <TF1.h>

class AliITSgeom;

// segmentation and response for SPD 

class AliITSsegmentationSPD :
public AliITSsegmentation {
 public:

    AliITSsegmentationSPD();
    AliITSsegmentationSPD(AliITSgeom *gm);
    AliITSsegmentationSPD(AliITSsegmentationSPD &source);
    virtual ~AliITSsegmentationSPD(){}
    AliITSsegmentationSPD& operator=(AliITSsegmentationSPD &source);

    // Set Detector Segmentation Parameters
    //
    // Detector size along x,z,y coordinates  
    virtual void    SetDetSize(Float_t Dx, Float_t Dz, Float_t Dy);

    // Maximum number of pixels along the two coordinates  
    virtual void    SetNPads(Int_t p1, Int_t p2);
    // Returns the maximum number of cells (digits) posible
    virtual Int_t   GetNPads(){return fNpx*fNpz;}
    // Set Pixel Size Array in x and z, microns.
    virtual void    SetBinSize(Float_t *x,Float_t *z);

    // Transform from real to pixel coordinates
    virtual void    GetPadIxz(Float_t x,Float_t z,Int_t &ix,Int_t &iz);
    // Transform from pixel to real coordinates
    virtual void    GetPadCxz(Int_t ix,Int_t iz,Float_t &x,Float_t &z);
    // Transform from real global to local coordinates
    //virtual void    GetLocal(Int_t module,Float_t *g ,Float_t *l) {}
    // Transform from real local to global coordinates
    //virtual void    GetGlobal(Int_t module,Float_t *l ,Float_t *g) {}
    // Local transformation of real local coordinates -
    virtual void    GetPadTxz(Float_t &x ,Float_t &z);
    // Transformation from Geant cm detector center local coordinates
    // to detector segmentation/cell coordiantes starting from (0,0).
    virtual void    LocalToDet(Float_t x,Float_t z,Int_t &ix,Int_t &iz);
    // Transformation from detector segmentation/cell coordiantes starting
    // from (0,0) to Geant cm detector center local coordinates.
    virtual void    DetToLocal(Int_t ix,Int_t iz,Float_t &x,Float_t &z);
    // Returns the Cell upper and lower boundries in x and y. cell indexes
    // starting from (0,0) and return Geant cm detector centered local
    // coordinates, consistant with DetToLocal and LocalToDet functions above.
    virtual void CellBoundries(Int_t ix,Int_t iz,Double_t &xl,Double_t &xu,
			       Double_t &zl,Double_t &zu);
    //
    // Initialisation
    virtual void Init();
    virtual void Init300();
    //
    // Get member data
    //
    // Detector Type geometry
    virtual AliITSgeom* Geometry() {return fGeom;}
    // Detector length
    virtual Float_t Dx() {return fDx;}
    // Detector width
    virtual Float_t Dz() {return fDz;}
    // Detector thickness
    virtual Float_t Dy() {return fDy;}
    // Pixel size in x
    virtual Float_t Dpx(Int_t ix);
    // Pixel size in z 
    virtual Float_t Dpz(Int_t iz);

    // Maximum number of Pixels in x
    virtual Int_t    Npx() {return fNpx;}
    // Maximum number of Pixels in z
    virtual Int_t    Npz(){return fNpz;}
    //
    // Get next neighbours
    virtual void Neighbours
       (Int_t iX,Int_t iZ,Int_t* Nlist,Int_t Xlist[10],Int_t Zlist[10]);

 private:
    Float_t ColFromZ300(Float_t z);
    Float_t ZFromCol300(Int_t col);
    Float_t ZpitchFromCol300(Int_t col);
    Float_t ColFromZ(Float_t z);
    Float_t ZFromCol(Int_t col);
    Float_t ZpitchFromCol(Int_t col);
    
  protected:

    Int_t   fNpx;           // Number of pixels in x
    Int_t   fNpz;           // Number of pixels in z
    Float_t fDx;            // Full width of the detector (x axis)- microns
    Float_t fDz;            // Full length of the detector (z axis)- microns
    Float_t fDy;            // Full thickness of the detector (y axis) -um 
    Float_t fCellSizeX[256];// Size for each pixel in x -microns
    Float_t fCellSizeZ[280];// Size for each pixel in z -microns
    TF1*    fCorr;          // correction function
    AliITSgeom *fGeom;      //! local pointer to AliITSgeom.

  ClassDef(AliITSsegmentationSPD,1) //Segmentation class for SPD 

};

#endif
