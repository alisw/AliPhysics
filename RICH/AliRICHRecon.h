#ifndef AliRICHRecon_h
#define AliRICHRecon_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRICHRecon                                                         //
//                                                                      //
// RICH class to perfom pattern recognition based on Hough transfrom    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <TTask.h> //base class
#include <TVector3.h>

class TClonesArray;
class AliRICHRecon : public TTask 
{
public : 
             AliRICHRecon();
    virtual ~AliRICHRecon()                                                          {}

  
  Double_t CkovAngle    (TClonesArray *pCluLst,Int_t &iNaccepted);                                                         //reconstructed Theta Cerenkov
  Double_t CkovSigma2   (                                                                     )const{ return fCkovSigma2;} //track ckov angle error squared
  Double_t FindPhotCkov (Double_t cluX,Double_t cluY                                          );     //find ckov angle for single photon candidate
  Double_t FindPhotPhi  (Double_t cluX,Double_t cluY                                          );     //find phi angle for single photon candidate
  Double_t FindRingCkov (Int_t iNclus                                                         );     //best ckov for ring formed by found photon candidates
  Double_t FindRingArea (Double_t ckov                                                        )const;//estimated area of ring in cm^2
  Int_t    FlagPhot     (Double_t ckov                                                        );     //is photon ckov near most probable track ckov
  Double_t HoughResponse(                                                                     );     //most probable track ckov angle
  void     Propagate    (const TVector3 &dir,      TVector3 &pos,Double_t z                   )const;//propagate photon alogn the line  
  void     Refract      (      TVector3 &dir,                    Double_t n1,    Double_t n2  )const;//refract photon on the boundary
  Double_t TracePhoton  (Double_t ckovTh,Double_t ckovPh,TVector2 &pos                        )const;//trace photon created by track to PC
  
  Double_t SigLoc  (Double_t ckovTh,Double_t ckovPh,Double_t trkTh,Double_t trkPh,Double_t beta)const; //localization error
  Double_t SigGeom (Double_t ckovTh,Double_t ckovPh,Double_t trkTh,Double_t trkPh,Double_t beta)const; //geometry error
  Double_t SigCrom (Double_t ckovTh,Double_t ckovPh,Double_t trkTh,Double_t trkPh,Double_t beta)const; //cromasity error
  Double_t Sigma2  (Double_t ckovTh,Double_t ckovPh,Double_t trkTh,Double_t trkPh              )const; //photon candidate sigma
  void     SetTrack(Double_t th,Double_t ph,Double_t x,Double_t y){  fTrkDir.SetMagThetaPhi(1,th,ph);  fTrkPos.Set(x,y);}//set track info
   
  const static Double_t fkRadThick;                      //radiator thickness
  const static Double_t fkWinThick;                      //window thickness
  const static Double_t fkGapThick;                      //proximity gap thickness
  const static Double_t fkRadIdx;                        //mean refractive index of RAD material (C6F14)
  const static Double_t fkWinIdx;                        //mean refractive index of WIN material (SiO2) 
  const static Double_t fkGapIdx;                        //mean refractive index of GAP material (CH4)
  
protected:
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
  ClassDef(AliRICHRecon,0)
};
    
#endif // #ifdef AliRICHRecon_cxx

