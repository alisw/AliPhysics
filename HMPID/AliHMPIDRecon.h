#ifndef AliHMPIDRecon_h
#define AliHMPIDRecon_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliHMPIDRecon                                                         //
//                                                                      //
// HMPID class to perfom pattern recognition based on Hough transfrom    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <TTask.h>        //base class
#include <TVector3.h>     //fields 

class TClonesArray; //CkovAngle()
class AliESDtrack;  //CkovAngle()

class AliHMPIDRecon : public TTask 
{
public : 
             AliHMPIDRecon();
    virtual ~AliHMPIDRecon()                                                          {}

  
  void     CkovAngle    (AliESDtrack *pTrk,TClonesArray *pCluLst,Double_t nmean             );                           //reconstructed Theta Cerenkov
  Double_t FindPhotCkov (Double_t cluX,Double_t cluY                                        );     //find ckov angle for single photon candidate
  Double_t FindPhotPhi  (Double_t cluX,Double_t cluY                                        );     //find phi angle for single photon candidate
  Double_t FindRingCkov (Int_t iNclus                                                       );     //best ckov for ring formed by found photon candidates
  Double_t FindRingArea (Double_t ckov                                                      )const;//estimated area of ring in cm^2
  Int_t    FlagPhot     (Double_t ckov                                                      );     //is photon ckov near most probable track ckov
  Double_t HoughResponse(                                                                   );     //most probable track ckov angle
  void     Propagate    (const TVector3 &dir,      TVector3 &pos,Double_t z                 )const;//propagate photon alogn the line  
  void     Refract      (      TVector3 &dir,                    Double_t n1,    Double_t n2)const;//refract photon on the boundary
  Double_t TracePhot    (Double_t ckovTh,Double_t ckovPh,TVector2 &pos                      )const;//trace photon created by track to PC 
  void     SetTrack     (Double_t x,Double_t y,Double_t theta,Double_t phi                  ){fTrkDir.SetMagThetaPhi(1,theta,phi);  fTrkPos.Set(x,y);}//set track
  Double_t SigLoc       (Double_t ckovTh,Double_t ckovPh,Double_t beta                      )const;//error due to cathode segmetation
  Double_t SigGeom      (Double_t ckovTh,Double_t ckovPh,Double_t beta                      )const;//error due to unknown photon origin
  Double_t SigCrom      (Double_t ckovTh,Double_t ckovPh,Double_t beta                      )const;//error due to unknonw photon energy
  Double_t Sigma2       (Double_t ckovTh,Double_t ckovPh                                    )const;//photon candidate sigma
  enum ETrackingFlags {kMipDistCut=-9,kMipQdcCut=-5,kNoPhotAccept=-11};
protected:
  static const Double_t fgkRadThick;                      //radiator thickness
  static const Double_t fgkWinThick;                      //window thickness
  static const Double_t fgkGapThick;                      //proximity gap thickness
  static const Double_t fgkWinIdx;                        //mean refractive index of WIN material (SiO2) 
  static const Double_t fgkGapIdx;                        //mean refractive index of GAP material (CH4)
  Double_t fRadNmean;                          //C6F14 mean refractive index
  Int_t    fPhotCnt;                           // counter of photons candidate
  Int_t    fPhotFlag[3000];                    // flags of photon candidates
  Double_t fPhotCkov[3000];                    // Ckov angles of photon candidates, [rad]
  Double_t fPhotPhi [3000];                    // phis of photons candidates, [rad]
  Double_t fPhotWei [3000];                    // weigths of photon candidates
  Double_t fCkovSigma2;                        // sigma2 of the reconstructed ring

  Bool_t  fIsWEIGHT;                          // flag to consider weight procedure
  Float_t fDTheta;                            // Step for sliding window
  Float_t fWindowWidth;                       // Hough width of sliding window
  
  TVector3 fTrkDir;                           //track direction in LORS
  TVector2 fTrkPos;                           //track positon in LORS at the middle of radiator
  ClassDef(AliHMPIDRecon,0)
};

#endif // #ifdef AliHMPIDRecon_cxx

