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


#include <TNamed.h>        //base class
#include <TVector3.h>     //fields 

class TClonesArray; //CkovAngle()
class AliESDtrack;  //CkovAngle()
class AliHMPIDParam;//general pourpose

class AliHMPIDRecon : public TNamed 
{
public : 
             AliHMPIDRecon();
    virtual ~AliHMPIDRecon() {;} //dtor

  void     InitVars     (Int_t n);                                                                 //init space for variables
  void     DeleteVars   ()const;                                                                   //delete variables
  void     CkovAngle    (AliESDtrack *pTrk,TClonesArray *pCluLst,Int_t index,Double_t nmean,Float_t xRa,Float_t yRa );//reconstructed Theta Cerenkov
  Bool_t   FindPhotCkov (Double_t cluX,Double_t cluY,Double_t &thetaCer,Double_t &phiCer    );     //find ckov angle for single photon candidate
  Double_t FindRingCkov (Int_t iNclus                                                       );     //best ckov for ring formed by found photon candidates
  void     FindRingGeom (Double_t ckovAng,Int_t level=1                                     );     //estimated area of ring in cm^2 and portion accepted by geometry
  TVector2 IntWithEdge  (TVector2 p1,TVector2 p2                                            )const;//find intercection between plane and lines of 2 thetaC
  Int_t    FlagPhot     (Double_t ckov,TClonesArray *pCluLst,AliESDtrack *pTrk              );     //is photon ckov near most probable track ckov
  Double_t HoughResponse(                                                                   );     //most probable track ckov angle
  void     Propagate    (const TVector3  dir,      TVector3 &pos,Double_t z                 )const;//propagate photon alogn the line  
  void     Refract      (      TVector3 &dir,                    Double_t n1,    Double_t n2)const;//refract photon on the boundary
  TVector2 TracePhot    (Double_t ckovTh,Double_t ckovPh                                    )const;//trace photon created by track to PC 
  void     AddObjectToFriends(TClonesArray *pCluLst, Int_t photonIndex, AliESDtrack *pTrk   );     // Add AliHMPIDCluster object to ESD friends
  TVector2 TraceForward (TVector3 dirCkov                                                   )const;//tracing forward a photon from (x,y) to PC
  void     Lors2Trs     (TVector3 dirCkov,Double_t &thetaCer,Double_t &phiCer               )const;//LORS to TRS 
  void     Trs2Lors     (TVector3 dirCkov,Double_t &thetaCer,Double_t &phiCer               )const;//TRS to LORS
  TVector2 GetMip       (                                                                   )const 
                        {return fMipPos;}                                                          //mip coordinates
  Double_t GetRingArea  (                                                                   )const
                        {return fRingArea;}                                                        //area of the current ring in cm^2 
  Double_t GetRingAcc   (                                                                   )const
                        {return fRingAcc;}                                                         //portion of the ring ([0,1]) accepted by geometry.To scale n. of photons 
  Double_t FindRingExt  (Double_t ckov,Int_t ch,Double_t xPc,Double_t yPc,Double_t thRa,Double_t phRa);//find ring acceptance by external parameters
  void     SetTrack     (Double_t xRad,Double_t yRad,Double_t theta,Double_t phi            )
                                {fTrkDir.SetMagThetaPhi(1,theta,phi);  fTrkPos.Set(xRad,yRad);}    //set track parameter at RAD
  void     SetImpPC     (Double_t xPc,Double_t yPc                                          )
                                {fPc.Set(xPc,yPc);}                                                //set track impact to PC 
  void     SetMip       (Double_t xmip,Double_t ymip                                        )
                                {fMipPos.Set(xmip,ymip);}                                          //set track impact to PC
  enum ETrackingFlags {kNotPerformed=-20,kMipDistCut=-9,kMipQdcCut=-5,kNoPhotAccept=-11,kNoRad = -22};
//
protected:
  Int_t     fPhotCnt;                           // counter of photons candidate
  Int_t    *fPhotFlag;                          // flags of photon candidates
  Int_t    *fPhotClusIndex;                     // cluster index of photon candidates
  Double_t *fPhotCkov;                          // Ckov angles of photon candidates, [rad]
  Double_t *fPhotPhi;                           // phis of photons candidates, [rad]
  Double_t *fPhotWei;                           // weigths of photon candidates
  Double_t  fCkovSigma2;                        // sigma2 of the reconstructed ring

  Bool_t    fIsWEIGHT;                          // flag to consider weight procedure
  Float_t   fDTheta;                            // Step for sliding window
  Float_t   fWindowWidth;                       // Hough width of sliding window
  
  Double_t  fRingArea;                          // area of a given ring
  Double_t  fRingAcc;                           // fraction of the ring accepted by geometry
  TVector3  fTrkDir;                            // track direction in LORS at RAD
  TVector2  fTrkPos;                            // track positon in LORS at RAD
  TVector2  fMipPos;                            // mip positon for a given track
  TVector2  fPc;                                // track position at PC
  
  AliHMPIDParam *fParam;                        // Pointer to AliHMPIDParam
  
private:
  AliHMPIDRecon(const AliHMPIDRecon& r);              //dummy copy constructor
  AliHMPIDRecon &operator=(const AliHMPIDRecon& r);   //dummy assignment operator
//
  ClassDef(AliHMPIDRecon,3)
};

#endif // #ifdef AliHMPIDRecon_cxx

