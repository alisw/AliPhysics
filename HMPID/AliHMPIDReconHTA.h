#ifndef AliHMPIDReconHTA_h
#define AliHMPIDReconHTA_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliHMPIDReconHTA                                                     //
//                                                                      //
// HMPID class to perfom pattern recognition based on Hough transfrom   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <TNamed.h>        //base class

class TClonesArray; //CkovAngle()
class AliESDtrack;  //CkovAngle()
class AliHMPIDParam;//General pourpose

class AliHMPIDReconHTA : public TNamed 
{
public : 
             AliHMPIDReconHTA();
    virtual ~AliHMPIDReconHTA();

// HTA hidden track algorithm
  void     InitVars         (Int_t n);                                                             //init space for variables
  void     DeleteVars       ()const;                                                               //delete variables
  void     InitDatabase     ();                                                                    //initialization of database
  void     FindBinDB        (Double_t x,Double_t y,Int_t &binX,Int_t &binY);                       //find indices of DB
  Bool_t   UniformDistrib   ();                                                                    //to check isotropy in phi Cerenkov
  void     FillZeroChan     ()const;                                                               //complete the DB
  Bool_t   CkovHiddenTrk    (AliESDtrack *pTrk,TClonesArray *pClu,Int_t index, Double_t nmean);    //Pattern recognition without trackinf information
  Bool_t   CluPreFilter     (TClonesArray *pClu               );                                   //Pre clustering filter to cut bkg clusters
  Bool_t   DoRecHiddenTrk   (                                 );                                   //Calling to the fitted procedures
  Bool_t   FindShape        (Double_t &thTrkRec,Double_t &phiTrkRec,Double_t &thetaCRec);          //Find shape of the ring
  Bool_t   ShapeModel       (Int_t np,Double_t *phiphot,Double_t *dist,Double_t &xA,Double_t &xB,Double_t &phiStart);//initial shape model for the cluster candidates
  Double_t VertParab        (Double_t x1,Double_t y1,Double_t x2, Double_t y2, Double_t x3, Double_t y3)const;//calculate the coord. of the min. for a parabole for 3 points
  Bool_t   FitFree          (Double_t thTrkRec,Double_t phiTrkRec);                                //Fit (th,ph) of the track and ckovFit as result
  Bool_t   FitRing          (Double_t thTrkRec,Double_t phiTrkRec);                                //Fit (th,ph) of the track 
  Double_t FindSimmPhi      ();                                                                    //find phi of the ring with min. dist. algorithm 
  Int_t    r2               (Double_t *coef, Double_t &x1, Double_t &x2);                          // solution of 2nd degree equation
  void     SetNClu          (Int_t nclu                          ) {fNClu=nclu;}                   //Setter for # of clusters
  void     SetClCk          (Int_t i,Bool_t what                 ) {fClCk[i]=what;}                //Setter for cluster flags 
  void     SetCkovFit       (Double_t ckov                       ) {fCkovFit=ckov;}                //Setter for ckov fitted
  void     SetCkovSig2      (Double_t rms                        ) {fCkovSig2=rms;}                //Setter for sigma2 ckov fitted
  void     SetTrkFit        (Double_t th,Double_t ph             ) {fThTrkFit = th;fPhTrkFit = ph;}//Setter for (th,ph) of the track
  void     SetNCluFit       (Int_t ncluFit                       ) {fNCluFit=ncluFit;}             //Setter for # of clusters used in the fit
  void     SetPhotAngles    (Int_t i,Double_t th, Double_t ph    ) {fPhiPhot[i]=ph;fThetaPhot[i]=th;} //Setter of the Cerenkov angles for a given photon
  void     SetFitStatus     (Bool_t status                       ) {fFitStatus = status;}          //Setter for fit status
  void     SetRadXY         (Double_t  x,Double_t y              ) {fRadX = x;fRadY = y;}          //Setter for (th,ph) of the track
  static void     FunMinPhot(Int_t&/* */,Double_t* /* */,Double_t &f,Double_t *par,Int_t iflag);   //Fit function to minimize thetaCer RMS/Sqrt(n) of n clusters
  Int_t CompactDB       (Int_t binX,Int_t binY) const {return fgDB[binX][binY];}                   //find compact word of DB
  Int_t    IdxMip       ()const {return fIdxMip;}                                                  //Getter index of MIP
  Double_t MipX         ()const {return fMipX;}                                                    //Getter of x MIP in LORS
  Double_t MipY         ()const {return fMipY;}                                                    //Getter of y MIP in LORS
  Double_t MipQ         ()const {return fMipQ;}                                                    //Getter of Q MIP
  Double_t RadX         ()const {return fRadX;}                                                    //Getter of x at RAD in LORS
  Double_t RadY         ()const {return fRadY;}                                                    //Getter of y at RAD in LORS
  Int_t    NClu         ()const {return fNClu;}                                                    //Getter of cluster multiplicity
  Int_t    NCluFit      ()const {return fNCluFit;}                                                 //Getter of n. photons used to fit the ring
  Double_t XClu         (Int_t i)const {return fXClu[i];}                                          //Getter of x clu
  Double_t YClu         (Int_t i)const {return fYClu[i];}                                          //Getter of y clu
  Bool_t   ClCk         (Int_t i)const {return fClCk[i];}                                          //Getter of cluster flags
  Double_t PhotPhi      (Int_t i)const {return fPhiPhot[i];}                                       //Getter of the Phi   Cerenkov angle for a given photon
  Double_t PhotTheta    (Int_t i)const {return fThetaPhot[i];}                                     //Getter of the Theta Cerenkov angle for a given photon
  Bool_t   FitStatus    ()const {return fFitStatus;}                                               //Getter of status of the fit
  Double_t CkovFit      ()const {return fCkovFit;}                                                 //Getter of ckov angle fitted
  Double_t ThTrkIn      ()const {return fThTrkIn;}                                                 //Getter of theta started of the track
  Double_t PhTrkIn      ()const {return fPhTrkIn;}                                                 //Getter of phi started of the track
  Double_t ThTrkFit     ()const {return fThTrkFit;}                                                //Getter of theta fitted of the track
  Double_t PhTrkFit     ()const {return fPhTrkFit;}                                                //Getter of phi fitted of the track
//
protected:
//
  Double_t fMipX;                              //mip X position for Hidden Track Algorithm  
  Double_t fMipY;                              //mip Y position for Hidden Track Algorithm
  Double_t fMipQ;                              //mip Q          for Hidden Track Algorithm
  Double_t fRadX;                              //rad X position for Hidden Track Algorithm  
  Double_t fRadY;                              //rad Y position for Hidden Track Algorithm
  Int_t    fIdxMip;                            //mip index in the clus list
  Int_t    fNClu;                              //n clusters to fit
  Double_t *fXClu;                             //container for x clus position
  Double_t *fYClu;                             //container for y clus position
  Double_t *fPhiPhot;                          //container for phi clus
  Double_t *fThetaPhot;                        //container for theta Cerenkov clus
  Bool_t   *fClCk;                             //flag if cluster is used in fitting
  Double_t fThTrkIn;                           //theta started from ShapeModel
  Double_t fPhTrkIn;                           //phi   started from ShapeModel
  Double_t fThTrkFit;                          //theta fitted of the track
  Double_t fPhTrkFit;                          //phi   fitted of the track
  Double_t fCkovFit;                           //estimated ring Cherenkov angle
  Int_t    fNCluFit;                           //n clusters used to fit the ring
  Double_t fCkovSig2;                          //estimated error^2 on ring Cherenkov angle
  Bool_t   fFitStatus;                         //status of the fit 0=ok 1=still to be optimized
  
  AliHMPIDParam *fParam;                       //Pointer to AliHMPIDParam
  static Int_t fgDB[500][150];                 //tmp DB
//
private:
  AliHMPIDReconHTA(const AliHMPIDReconHTA& r);              //dummy copy constructor
  AliHMPIDReconHTA &operator=(const AliHMPIDReconHTA& r);   //dummy assignment operator
//
  ClassDef(AliHMPIDReconHTA,5)
};

#endif // #ifdef AliHMPIDReconHTA_cxx

