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


  void     CkovAngle    (AliESDtrack *pTrk,TClonesArray *pCluLst,Double_t nmean, Double_t qthre);  //reconstructed Theta Cerenkov
  Bool_t   FindPhotCkov (Double_t cluX,Double_t cluY,Double_t &thetaCer,Double_t &phiCer    );     //find ckov angle for single photon candidate
  Double_t FindRingCkov (Int_t iNclus                                                       );     //best ckov for ring formed by found photon candidates
  Double_t FindRingArea (Double_t ckovAngMin,Double_t ckovAngMax                            )const;//estimated area of delta ring in cm^2 to weight Hough Transform
  Int_t    FlagPhot     (Double_t ckov                                                      );     //is photon ckov near most probable track ckov
  Double_t HoughResponse(                                                                   );     //most probable track ckov angle
  void     Propagate    (const TVector3  dir,      TVector3 &pos,Double_t z                 )const;//propagate photon alogn the line  
  void     Refract      (      TVector3 &dir,                    Double_t n1,    Double_t n2)const;//refract photon on the boundary
  TVector2 TracePhot    (Double_t ckovTh,Double_t ckovPh                                    )const;//trace photon created by track to PC 
  TVector2 TraceForward (TVector3 dirCkov                                                   )const;//tracing forward a photon from (x,y) to PC
  void     RecPhot      (TVector3 dirCkov,Double_t &thetaCer,Double_t &phiCer               );     //theta,phi cerenkov reconstructed
  void     SetTrack     (Double_t xRad,Double_t yRad,Double_t theta,Double_t phi            )
                                {fTrkDir.SetMagThetaPhi(1,theta,phi);  fTrkPos.Set(xRad,yRad);}    //set track parameter at RAD
  void     SetImpPC     (Double_t xPc,Double_t yPc                                          )
                                {fPc.Set(xPc,yPc);}                                                //set track impact to PC 
  Double_t SigLoc       (Double_t ckovTh,Double_t ckovPh,Double_t beta                      )const;//error due to cathode segmetation
  Double_t SigGeom      (Double_t ckovTh,Double_t ckovPh,Double_t beta                      )const;//error due to unknown photon origin
  Double_t SigCrom      (Double_t ckovTh,Double_t ckovPh,Double_t beta                      )const;//error due to unknonw photon energy
  Double_t Sigma2       (Double_t ckovTh,Double_t ckovPh                                    )const;//photon candidate sigma^2
  enum ETrackingFlags {kMipDistCut=-9,kMipQdcCut=-5,kNoPhotAccept=-11};
// HTA hidden track algorithm
  Bool_t   CkovHiddenTrk    (AliESDtrack *pTrk,TClonesArray *pCluLst,Double_t nmean,Double_t qthre);//Pattern recognition without trackinf information
  Bool_t   CluPreFilter     (TClonesArray *pClu               );                                   //Pre clustering filter to cut bkg clusters
  Bool_t   DoRecHiddenTrk   (TClonesArray *pClu               );                                   //Calling to the fitted procedures
  Bool_t   FitEllipse       (Double_t &phiRec                    );                                //Fit clusters with a conical section (kTRUE only for ellipses)
  Bool_t   FitFree          (Double_t phiRec                     );                                //Fit (th,ph) of the track and ckovFit as result
  Double_t FunConSect       (Double_t *c,Double_t x,Double_t y   );                                //Function of a general conical section
  void     SetNClu          (Int_t nclu                          ) {fNClu=nclu;}                   //Setter for # of clusters
  void     SetClCk          (Int_t i,Bool_t what                 ) {fClCk[i]=what;}                //Setter for cluster flags 
  void     SetCkovFit       (Double_t ckov                       ) {fCkovFit=ckov;}                //Setter for ckov fitted
  void     SetCkovSig2      (Double_t rms                       ) {fCkovSig2=rms;}                 //Setter for sigma2 ckov fitted
  void     SetTrkFit        (Double_t th,Double_t ph             ) {fThTrkFit = th;fPhTrkFit = ph;}//Setter for (th,ph) of the track
  void     SetRadXY         (Double_t  x,Double_t y              ) {fRadX = x;fRadY = y;}          //Setter for (th,ph) of the track
  static void     FunMinEl  (Int_t&/* */,Double_t* /* */,Double_t &f,Double_t *par,Int_t /* */);   //Fit function to find ellipes parameters
  static void     FunMinPhot(Int_t&/* */,Double_t* /* */,Double_t &f,Double_t *par,Int_t iflag);   //Fit function to minimize thetaCer RMS/Sqrt(n) of n clusters
  Int_t    IdxMip       ()const {return fIdxMip;}                                                  //Getter index of MIP
  Double_t MipX         ()const {return fMipX;}                                                    //Getter of x MIP in LORS
  Double_t MipY         ()const {return fMipY;}                                                    //Getter of y MIP in LORS
  Double_t MipQ         ()const {return fMipQ;}                                                    //Getter of Q MIP
  Double_t RadX         ()const {return fRadX;}                                                    //Getter of x at RAD in LORS
  Double_t RadY         ()const {return fRadY;}                                                    //Getter of y at RAD in LORS
  Int_t    NClu         ()const {return fNClu;}                                                    //Getter of cluster multiplicity
  Double_t XClu         (Int_t i)const {return fXClu[i];}                                          //Getter of x clu
  Double_t YClu         (Int_t i)const {return fYClu[i];}                                          //Getter of y clu
  Bool_t   ClCk         (Int_t i)const {return fClCk[i];}                                          //Getter of cluster flags
  Double_t CkovFit      ()const {return fCkovFit;}                                                 //Getter of ckov angle fitted
  Double_t ThTrkFit     ()const {return fThTrkFit;}                                                //Getter of theta fitted of the track
  Double_t PhTrkFit     ()const {return fPhTrkFit;}                                                //Getter of phi fitted of the track
//
protected:
  Double_t fRadNmean;                          //C6F14 mean refractive index
  Int_t    fPhotCnt;                           // counter of photons candidate
  Int_t    fPhotFlag[3000];                    // flags of photon candidates
  Double_t fPhotCkov[3000];                    // Ckov angles of photon candidates, [rad]
  Double_t fPhotPhi [3000];                    // phis of photons candidates, [rad]
  Double_t fPhotWei [3000];                    // weigths of photon candidates
  Double_t fCkovSigma2;                        // sigma2 of the reconstructed ring

  Bool_t   fIsWEIGHT;                          // flag to consider weight procedure
  Float_t  fDTheta;                            // Step for sliding window
  Float_t  fWindowWidth;                       // Hough width of sliding window
  
  TVector3 fTrkDir;                            //track direction in LORS at RAD
  TVector2 fTrkPos;                            //track positon in LORS at RAD
  TVector2 fPc;                                //track position at PC
// HTA hidden track algorithm
  Double_t fMipX;                              //mip X position for Hidden Track Algorithm  
  Double_t fMipY;                              //mip Y position for Hidden Track Algorithm
  Double_t fMipQ;                              //mip Q          for Hidden Track Algorithm
  Double_t fRadX;                              //rad X position for Hidden Track Algorithm  
  Double_t fRadY;                              //rad Y position for Hidden Track Algorithm
  Int_t    fIdxMip;                            //mip index in the clus list
  Int_t    fNClu;                              //n clusters to fit
  Double_t fXClu[100];                         //container for x clus position
  Double_t fYClu[100];                         //container for y clus position
  Bool_t   fClCk[100];                         //flag if cluster is used in fitting
  Double_t fThTrkFit;                          //theta fitted of the track
  Double_t fPhTrkFit;                          //phi   fitted of the track
  Double_t fCkovFit;                           //estimated ring Cherenkov angle
  Double_t fCkovSig2;                          //estimated error^2 on ring Cherenkov angle
//
private:
  static const Double_t fgkRadThick;                      //radiator thickness
  static const Double_t fgkWinThick;                      //window thickness
  static const Double_t fgkGapThick;                      //proximity gap thickness
  static const Double_t fgkWinIdx;                        //mean refractive index of WIN material (SiO2) 
  static const Double_t fgkGapIdx;                        //mean refractive index of GAP material (CH4)

  ClassDef(AliHMPIDRecon,0)
};

#endif // #ifdef AliHMPIDRecon_cxx

