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
    virtual void    SetDetSize(Float_t p1=72960.,Float_t p2=40000.,
			       Float_t p3= 300.) 
                        {fDx=p1; fDz=p2; fDy=p3;}

    // Strip size  
    virtual void    SetPadSize(Float_t pitch=95.,Float_t d=1.0)
	{fPitch=pitch;d=1.;}

    // Maximum number of strips along the two coordinates  
    virtual void    SetNPads(Int_t p1=768,Int_t d=1){fNstrips=p1;d=1;}
    // Returns the maximum number of cells (digits) posible
    virtual Int_t   GetNPads(){return 2*fNstrips;}


    // Set stereo angles Pside-Nside 
    virtual void    SetAngles(Float_t pa=0.0175, Float_t na=0.0175) 
                         {fStereoPl5=pa; fStereoNl5=na;
			 fStereoPl6=na; fStereoNl6=pa;}

    virtual void    SetAnglesLay5(Float_t pa=0.0075, Float_t na=0.0275) 
                         {fStereoPl5=pa; fStereoNl5=na;}

    virtual void    SetAnglesLay6(Float_t pa=0.0275, Float_t na=0.0075) 
                         {fStereoPl6=pa; fStereoNl6=na;}

    // Set stereo angles Pside-Nside 
    // Transform from real coordinates to strips
    virtual void    GetPadIxz(Float_t x ,Float_t z ,Int_t   &iP,Int_t  &iN);
    // Transform from strips to real coordinates
    virtual void    GetPadCxz(Int_t iP, Int_t iN, Float_t &x , Float_t &z);
    virtual void    GetPadTxz(Float_t &x , Float_t &z);

    // Transform from real global to local coordinates
    virtual void    GetLocal(Int_t,Float_t *,Float_t *) {}
    // Transform from real local to global coordinates
    virtual void    GetGlobal(Int_t,Float_t *,Float_t *) {}
    // Transformation from Geant cm detector center local coordinates
    // to detector P and N side strip numbers..
    virtual void    LocalToDet(Float_t x,Float_t z,Int_t &iP,Int_t &iN);
    // Transformation from detector segmentation/cell coordiantes starting
    // from 0. iPN=0 for P side and 1 for N side strip. Returned is z=0.0
    // and the corresponding x value..
    virtual void    DetToLocal(Int_t ix,Int_t iPN,Float_t &x,Float_t &z);
    // Given one P side strip and one N side strip, Returns kTRUE if they
    // cross each other and the location of the two crossing strips and
    // their correxlation matrix c[2][2].
    virtual Bool_t  GetCrossing(Int_t iP,Int_t iN,Float_t &x,Float_t &z,
				Float_t c[2][2]);

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
     virtual void Angles(Float_t &aP,Float_t &aN);
     virtual void SetLayer(Int_t l);
     virtual Int_t GetLayer() const {return fLayer;}

  protected:

  Int_t      fNstrips;       // Number of strips in x 
  Float_t    fStereoP;       // Stereo angle for Pside (rad)
  Float_t    fStereoN;       // Stereo angle for Nside (rad)
  Float_t    fPitch;         // Pitch of the strips
  Float_t    fDz;            // Full width of the detector (z axis)- microns
  Float_t    fDx;            // Full length of the detector (x axis)- microns
  Float_t    fDy;            // Full thickness of the detector (y axis) -um 
  
  Float_t    fStereoPl5;       // Stereo angle for Pside
  Float_t    fStereoNl5;       // Stereo angle for Nside
  Float_t    fStereoPl6;       // Stereo angle for Pside
  Float_t    fStereoNl6;       // Stereo angle for Nside
  Int_t      fLayer;

  AliITSgeom *fGeom;         //! pointer to the geometry class
  TF1*       fCorr;          // correction function
  
  ClassDef(AliITSsegmentationSSD,1) //Segmentation class for SSD 
};

#endif
