#ifndef ALIITSSEGMENTATIONSSD_H
#define ALIITSSEGMENTATIONSSD_H

#include "AliITSsegmentation.h"

// segmentation for SSD

class AliITSsegmentationSSD :
public AliITSsegmentation {
 public:

    AliITSsegmentationSSD();
    AliITSsegmentationSSD(AliITSgeom *gm);
    AliITSsegmentationSSD(const AliITSsegmentationSSD &source);
    virtual ~AliITSsegmentationSSD(){}
    AliITSsegmentationSSD& operator=(const AliITSsegmentationSSD &source);


    // Strip size  
    virtual void    SetPadSize(Float_t pitch,Float_t /* d */)
	{fPitch=pitch;}

    // Maximum number of strips along the two coordinates  
    virtual void    SetNPads(Int_t p1,Int_t /* d */){fNstrips=p1;}
    // Returns the maximum number of cells (digits) posible
    virtual Int_t   GetNPads() const {return 2*fNstrips;}


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
    virtual void    GetPadIxz(Float_t x ,Float_t z ,Int_t   &iP,Int_t  &iN) const;
    // Transform from strips to real coordinates
    virtual void    GetPadCxz(Int_t iP, Int_t iN, Float_t &x , Float_t &z) const;
    virtual void    GetPadTxz(Float_t &x , Float_t &z) const;
    // Transformation from Geant cm detector center local coordinates
    // to detector P and N side strip numbers..
    virtual Bool_t  LocalToDet(Float_t x,Float_t z,Int_t &iP,Int_t &iN) const;
    // Transformation from detector segmentation/cell coordiantes starting
    // from 0. iPN=0 for P side and 1 for N side strip. Returned is z=0.0
    // and the corresponding x value..
    virtual void    DetToLocal(Int_t ix,Int_t iPN,Float_t &x,Float_t &z) const;
    // Given one P side strip and one N side strip, Returns kTRUE if they
    // cross each other and the location of the two crossing strips and
    // their correxlation matrix c[2][2].
    virtual Bool_t  GetCrossing(Int_t iP,Int_t iN,Float_t &x,Float_t &z,
				Float_t c[2][2]);

    virtual void Init();

    // Strip size in x
    virtual Float_t Dpx(Int_t) const {return fPitch;}
    // Strip size in z 
    virtual Float_t Dpz(Int_t) const {return fDz;}
    // Maximum number of Strips in x
    virtual Int_t    Npx() const {return fNstrips;}
    // Maximum number of Strips in z
    virtual Int_t    Npz()const {return 1;}

    // Angles : Pside stereo angle-Nside stereo angle
     virtual void Angles(Float_t &aP,Float_t &aN) const;
     virtual void SetLayer(Int_t l);
     virtual Int_t GetLayer() const {return fLayer;}
    // Print Default parameters
     virtual void PrintDefaultParameters() const;

  protected:

  virtual void Copy(TObject &obj) const; 

  Int_t      fNstrips;       // Number of strips in x 
  Float_t    fStereoP;       // Stereo angle for Pside (rad)
  Float_t    fStereoN;       // Stereo angle for Nside (rad)
  Float_t    fPitch;         // Pitch of the strips
  
  Float_t    fStereoPl5;       // Stereo angle for Pside
  Float_t    fStereoNl5;       // Stereo angle for Nside
  Float_t    fStereoPl6;       // Stereo angle for Pside
  Float_t    fStereoNl6;       // Stereo angle for Nside
  Int_t      fLayer;           //! layer number (5 or 6)
  static const Float_t fgkDxDefault;  // Default value for fDx
  static const Float_t fgkDzDefault;  // Default value for fDz
  static const Float_t fgkDyDefault;  // Default value for fDy
  static const Float_t fgkPitchDefault; //Default value for fPitch
  static const Int_t fgkNstripsDefault; //Default value for fNstrips

  ClassDef(AliITSsegmentationSSD,2) //Segmentation class for SSD 
};

#endif
