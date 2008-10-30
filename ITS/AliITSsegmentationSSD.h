#ifndef ALIITSSEGMENTATIONSSD_H
#define ALIITSSEGMENTATIONSSD_H

#include "AliITSsegmentation.h"

// segmentation for SSD

/* $Id$ */

class AliITSsegmentationSSD :
public AliITSsegmentation {
 public:

    AliITSsegmentationSSD(Option_t *opt="");
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
    virtual void    SetAngles(Float_t pa=0.0075, Float_t na=0.0275) 
                         {fStereoP=pa; fStereoN=na;}

    // Set stereo angles Pside-Nside 
    // Transform from real coordinates to strips
    virtual void    GetPadIxz(Float_t x ,Float_t z ,Int_t   &iP,Int_t  &iN) const;
    // Transform from strips to real coordinates
            void    GetPadCxz(Float_t iP, Float_t iN, Float_t &x , Float_t &z) const;
    virtual void    GetPadCxz(Int_t iP, Int_t iN, Float_t &x , Float_t &z) const { GetPadCxz((Float_t) iP, (Float_t) iN, x, z); }
    virtual void    GetPadTxz(Float_t &x , Float_t &z) const;
    // Transformation from Geant cm detector center local coordinates
    // to detector P and N side strip numbers..
    virtual Bool_t  LocalToDet(Float_t x,Float_t z,Int_t &iP,Int_t &iN) const;
    // Transformation from detector segmentation/cell coordiantes starting
    // from 0. iPN=0 for P side and 1 for N side strip. Returned is z=0.0
    // and the corresponding x value..
    virtual void    DetToLocal(Int_t ix,Int_t iPN,Float_t &x,Float_t &z) const;

    virtual Int_t    GetNumberOfChips() const {
      return fgkNchipsPerSide;
    }
    virtual Int_t    GetMaximumChipIndex() const{
      return fgkNchipsPerSide*2-1;
    }
    virtual Int_t    GetChipFromLocal(Float_t xloc, Float_t zloc) const;
    virtual Int_t    GetChipFromChannel(Int_t ix, Int_t iz) const;
    virtual Int_t    GetChipsInLocalWindow(Int_t* array, Float_t zmin, Float_t zmax, Float_t xmin, Float_t xmax) const;

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
  
  Int_t      fLayer;           //! layer number (5 or 6)
  static const Float_t fgkDxDefault;  // Default value for fDx
  static const Float_t fgkDzDefault;  // Default value for fDz
  static const Float_t fgkDyDefault;  // Default value for fDy
  static const Float_t fgkPitchDefault; //Default value for fPitch
  static const Int_t fgkNstripsDefault; //Default value for fNstrips
  static const Int_t fgkNchipsPerSide;    //number of chips per side
  static const Int_t fgkNstripsPerChip;    //number of strips per chip

  ClassDef(AliITSsegmentationSSD,4) //Segmentation class for SSD 
};

#endif
