#ifndef ALIPHOSCORRECTIONFW_H
#define ALIPHOSCORRECTIONFW_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis for PHOS Tagged Photons 
// marks photons making pi0 with any other photon
// and calculates necessary corrections for fake pairs and
// decay partners escaped acceptance. If MC info is present 
// fills set of controll histograms.
//
//*-- Dmitry Peresunko
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"  
class AliAODEvent ; 
class AliPHOSGeometry;
class AliPHOSDigit ;
class AliPHOSEmcRecPoint ;
class TClonesArray ;
class TObjArray;
class AliPHOSCalibData ;

class AliPHOSCorrectionFW : public AliAnalysisTaskSE {

public:
  AliPHOSCorrectionFW() ;
  AliPHOSCorrectionFW(const char *name) ;
  AliPHOSCorrectionFW(const AliPHOSCorrectionFW& ap) ;   
  AliPHOSCorrectionFW & operator = (const AliPHOSCorrectionFW & ap) { Fatal("operator =", "not implemented"); return *this ;} 

  virtual ~AliPHOSCorrectionFW() ;
   
  virtual void UserCreateOutputObjects(){} 
  virtual void Init() ; 
  virtual void LocalInit() { Init() ; }
  virtual void UserExec(Option_t * opt = "") ;

  // void SetClusterBranch(const char * branchname="PHOSCluDef") ;
  static Double_t ShowerShape(Double_t x, Double_t z) ; // Shape of EM shower used in unfolding; 
                                            //class member function (not object member function)
 static void UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)  ;
                                            // Chi^2 of the fit. Should be static to be passed to MINUIT

 protected:
  void ConvertDigits() ; //Convert PHOS AODCells to AliPHOSDigits and apply calibration/bad map

  void Custerize() ;     //Clusterize PHOS digits
 
  void FillOutput() ;        //Fill resulting branch

  Int_t AreNeighbours(AliPHOSDigit * d1, AliPHOSDigit * d2)const ;

  Double_t DistancesToBadChannels(AliPHOSEmcRecPoint * rp)const ;

  Double_t Calibrate(Double_t e, Int_t absId) ;

  Double_t CalibrateT(Double_t t, Int_t absId, Bool_t isHG) ;

  Bool_t RunChanged(); 

  void  UnfoldCluster(AliPHOSEmcRecPoint * iniEmc,Int_t Nmax, 
           AliPHOSDigit ** maxAt,Float_t * maxAtEnergy ) ; //Unfolds cluster using TMinuit package

  Double_t DistanceToBadChannels(AliPHOSEmcRecPoint *rp) ;

  Bool_t  FindFit(AliPHOSEmcRecPoint * emcRP, AliPHOSDigit ** MaxAt, Float_t * maxAtEnergy, 
      Int_t NPar, Float_t * FitParametres) const; //Used in UnfoldClusters, calls TMinuit

  Double_t CorrectNonlinearity(Double_t en) ;

  Int_t FindTrackMatching(Int_t mod,TVector3 *locpos,
              Double_t &dx, Double_t &dz,
              Double_t &pt,Int_t &charge);

  Double_t TestCoreLambda(Double_t pt,Double_t l1,Double_t l2);

  Double_t TestFullLambda(Double_t pt,Double_t l1,Double_t l2);

  Double_t TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge); 

private:

  AliAODEvent * fEvent ;
  AliPHOSGeometry  * fPHOSGeo ;
  AliPHOSCalibData * fPHOSCalibData ;

  TClonesArray * fDigits ;
  TObjArray *    fClusters ;
  TClonesArray * fAODClusterArray ;

  TVector3 fVtx ;      //Vertex in current event

  TH2I * fPHOSBadMap[6] ;  //BadMaps
  Double_t fNonlinearityParams[10] ;         // Parameters for non-linearity calculation
  Float_t fRunByRunCorr[5] ;                 // Per module run-by-run correction
  Int_t fL1phase[15] ;                       // L1phases for PHOS DDLs (run2 only)
  Bool_t  fDigitsUsed[17920];        //Mark digits as already used in cluster (EMC:5*56*64 ; CPV: 5*56*128)

  Float_t fZScut ;     // Cut for zero suppression emulation in MC (in GeV)
  Float_t fNoiseMC;    // amount of electronic noise added to digits in MC (in GeV)
  Float_t fEmcMinE ;   // Minimal digit energy
  Float_t fClusteringThreshold;  // Minimal energy to start cluster (in GeV)
  Float_t fEmcLocMaxCut;         // Miminal height of local maximum (in GeV)
  Float_t fW0;                   // Weight (def=4.5)
  Float_t fEcoreRadius;          // cire radius (in cm)

  Int_t fRunNumber ;   // Run number from data
  Int_t fDRN ;         // Fixed RunNumber if ncecessary
  Int_t fRecoPass;     // Reco pass for calibration

  Bool_t fIsMC ;
  Bool_t fUsePrivateBadMap ; 
  Bool_t fUsePrivateCalib ;
  Bool_t fApplyZS ;           // Apply zero suppression emulation in MC
  Bool_t fAddNoiseMC;         // Smear cell energy by adding electroinc noise emulation in MC
  Bool_t fToUnfold;           // Unfold clusters

  TString fPrivateOADBBadMap ;
  TString fMCProduction ;
  TString fNonlinearityVersion ;
       
  ClassDef(AliPHOSCorrectionFW, 1);   
};
#endif // ALIPHOSCORRECTIONFW_H
