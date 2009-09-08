#ifndef ALITPCLASERTRACK
#define ALITPCLASERTRACK
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TString.h>

#include "AliExternalTrackParam.h"

class TObjArray;



class AliTPCLaserTrack : public AliExternalTrackParam {
public:
  AliTPCLaserTrack();
  AliTPCLaserTrack(const AliTPCLaserTrack &ltr);
  AliTPCLaserTrack(const Int_t id, const Int_t side, const Int_t rod,
                   const Int_t bundle, const Int_t beam,
                   Double_t x, Double_t alpha,
                   const Double_t param[5],
                   const Double_t covar[15], const Float_t rayLength=0);
  
  AliTPCLaserTrack& operator = (const  AliTPCLaserTrack &source);
  
  static void LoadTracks();
  static TObjArray* GetTracks() {return fgArrLaserTracks;}
  
  static Int_t IdentifyTrack(AliExternalTrackParam *track, Int_t side=-1);
  
  Int_t GetId()     const {return fId;     }
  Int_t GetSide()   const {return fSide;   }
  Int_t GetRod()    const {return fRod;    }
  Int_t GetBundle() const {return fBundle; }
  Int_t GetBeam()   const {return fBeam;   }

  Float_t GetRayLength() const {return fRayLength;}
  
  
  
  static Int_t GetNLaserTracks() { return fgkNLaserTracks; }
  static Int_t GetNLaserRodsPerSide() { return fgkNRodsPerSide; }
  static Int_t GetNMirrorBundlesPerRod() { return fgkNBundlePerRod; }
  static Int_t GetNLaserRaysPerMirrorBundle() { return fgkNBeamsPerBundle; }
  
  
  void SetId    (Int_t id)    {fId     = id;    }
  void SetSide  (Int_t side)  {fSide   = side;  }
  void SetRod   (Int_t rod)   {fRod    = rod;   }
  void SetBundle(Int_t bundle){fBundle = bundle;}
  void SetBeam  (Int_t beam)  {fBeam   = beam;  }
  void SetRayLength (Float_t len) {fRayLength = len;}
  
  
private:
  Int_t fId;              //Laser beam id            (0-335)
  Int_t fSide;            //TPC side; 0:Shaft Side (A) -- 1:Muon Side (C)
  Int_t fRod;             //Laser Rod                (0-5)
  Int_t fBundle;          //Mirror bundle in the Rod (0-3)
  Int_t fBeam;            //Laser Beam in the bundle (0-6)
  
  Float_t fRayLength;     //distance from the last common point of the laser Rays
                          //(Splitter box on the A-Side at the bottom of the TPC)
                          //to each mirror (needed for an exact drift velocity estimation)
  
  
  static TObjArray* fgArrLaserTracks; //! Array of all Laser Tracks,
                                        //  keeps instances of this class;
  
  static const Int_t fgkNLaserTracks    = 336; //Number of laser tracks
  static const Int_t fgkNRodsPerSide    = 6;   //Number of laser rods on each readout side
  static const Int_t fgkNBundlePerRod   = 4;   //Number of mirror bundles per rod
  static const Int_t fgkNBeamsPerBundle = 7;   //Number of laser rays per bundle
  
//    static const char* fgkDataFileName = "$ALIC_ROOT/TPC/Calib/LaserTracks.root";  //Path to the Data File
  
  ClassDef(AliTPCLaserTrack,2)        // Laser Track positions and track identification
};

#endif

