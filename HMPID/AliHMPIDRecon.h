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

class TClonesArray;
class AliHMPIDRecon : public TTask 
{
public : 
             AliHMPIDRecon();
    virtual ~AliHMPIDRecon()                                                          {}

  
  Double_t CkovAngle    (TClonesArray *pCluLst,Int_t &iNaccepted);                                                         //reconstructed Theta Cerenkov
  Double_t CkovSigma2   (                                                                   )const{ return fCkovSigma2;} //track ckov angle error squared
  Double_t FindPhotCkov (Double_t cluX,Double_t cluY                                        );     //find ckov angle for single photon candidate
  Double_t FindPhotPhi  (Double_t cluX,Double_t cluY                                        );     //find phi angle for single photon candidate
  Double_t FindRingCkov (Int_t iNclus                                                       );     //best ckov for ring formed by found photon candidates
  Double_t FindRingArea (Double_t ckov                                                      )const;//estimated area of ring in cm^2
  Int_t    FlagPhot     (Double_t ckov                                                      );     //is photon ckov near most probable track ckov
  Double_t HoughResponse(                                                                   );     //most probable track ckov angle
  void     Propagate    (const TVector3 &dir,      TVector3 &pos,Double_t z                 )const;//propagate photon alogn the line  
  void     Refract      (      TVector3 &dir,                    Double_t n1,    Double_t n2)const;//refract photon on the boundary
  Double_t TracePhot    (Double_t ckovTh,Double_t ckovPh,TVector2 &pos                      )const;//trace photon created by track to PC 
  void     SetTrack     (Double_t th,Double_t ph,Double_t x,Double_t y                      ){fTrkDir.SetMagThetaPhi(1,th,ph);  fTrkPos.Set(x,y);}//set track
  Double_t SigLoc       (Double_t ckovTh,Double_t ckovPh,Double_t beta                      )const;//error due to cathode segmetation
  Double_t SigGeom      (Double_t ckovTh,Double_t ckovPh,Double_t beta                      )const;//error due to unknown photon origin
  Double_t SigCrom      (Double_t ckovTh,Double_t ckovPh,Double_t beta                      )const;//error due to unknonw photon energy
  Double_t Sigma2       (Double_t ckovTh,Double_t ckovPh                                    )const;//photon candidate sigma
  
  static void  Display  (                                                                   );    //event display  
   
  
protected:
  const static Double_t fgkRadThick;                      //radiator thickness
  const static Double_t fgkWinThick;                      //window thickness
  const static Double_t fgkGapThick;                      //proximity gap thickness
  const static Double_t fgkRadIdx;                        //mean refractive index of RAD material (C6F14)
  const static Double_t fgkWinIdx;                        //mean refractive index of WIN material (SiO2) 
  const static Double_t fgkGapIdx;                        //mean refractive index of GAP material (CH4)
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

