#ifndef ALIITSSEGMENTATIONSSD_H
#define ALIITSSEGMENTATIONSSD_H

#include "AliITSsegmentation.h"

// segmentation for SSD

class AliITSsegmentationSSD :
public AliITSsegmentation {
 public:

    AliITSsegmentationSSD();
    AliITSsegmentationSSD(AliITSgeom *gm);
    AliITSsegmentationSSD(AliITSsegmentationSSD &source);
    virtual ~AliITSsegmentationSSD(){}
    AliITSsegmentationSSD& operator=(AliITSsegmentationSSD &source);


    // Detector size: x,z,y 
    virtual  void   SetDetSize
          (Float_t p1=72960., Float_t p2=40000., Float_t p3= 300.) 
                        {fDx=p1; fDz=p2; fDy=p3;}

    // Strip size  
    virtual void    SetCellSize(Float_t pitch=95., Float_t dummy=1.) 
                         {fPitch=pitch;}

    // Maximum number of strips along the two coordinates  
    virtual void    SetNCells(Int_t p1=768, Int_t dummy=1) 
                         {fNstrips=p1;}


    // Set stereo angles Pside-Nside 
    virtual void    SetAngles(Float_t pa=0.0175, Float_t na=0.0175) 
                         {fStereoP=pa; fStereoN=na;}

    // Transform from real coordinates to strips
    virtual void    GetCellIxz
    (Float_t &x ,Float_t &z ,Int_t   &iP,Int_t  &iN);
    // Transform from strips to real coordinates
    virtual void    GetCellCxz
    (Int_t iP, Int_t iN, Float_t &x , Float_t &z);

    // Transform from real global to local coordinates
    virtual void    GetLocal(Int_t module,Float_t *g ,Float_t *l) {}
    // Transform from real local to global coordinates
    virtual void    GetGlobal(Int_t module,Float_t *l ,Float_t *g) {}

    virtual void Init();

    // Detector type geometry
    virtual AliITSgeom* Geometry() {return 0;}
    // Detector length
    virtual Float_t Dx() {return fDx;}
    // Detector width
    virtual Float_t Dz() {return fDz;}
    // Detector thickness
    virtual Float_t Dy() {return fDy;}
    // Strip size in x
    virtual Float_t Dpx(Int_t) {return fPitch;}
    // Strip size in z 
    virtual Float_t Dpz(Int_t) {return fDz;}
    // Maximum number of Strips in x
    virtual Int_t    Npx() {return fNstrips;}
    // Maximum number of Strips in z
    virtual Int_t    Npz(){return 1;}

    // Angles : Pside stereo angle-Nside stereo angle
    virtual void Angles(Float_t &aP,Float_t &aN) 
                     {aP=fStereoP;aN=fStereoN;}

  protected:

  Int_t      fNstrips;       // Number of strips in x 
  Float_t    fStereoP;       // Stereo angle for Pside
  Float_t    fStereoN;       // Stereo angle for Nside
  Float_t    fPitch;         // Pitch of the strips
  Float_t    fDz;            // Full width of the detector (z axis)- microns
  Float_t    fDx;            // Full length of the detector (x axis)- microns
  Float_t    fDy;            // Full thickness of the detector (y axis) -um 
  
  AliITSgeom *fGeom;         // pointer to the geometry class
  TF1*       fCorr;          // correction function
  
  ClassDef(AliITSsegmentationSSD,1) //Segmentation class for SSD 
};

#endif
