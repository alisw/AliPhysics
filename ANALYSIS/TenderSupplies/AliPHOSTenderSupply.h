#ifndef ALIPHOSTENDERSUPPLY_H
#define ALIPHOSTENDERSUPPLY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  PHOS tender, apply corrections to PHOS clusters                   //
//  and do track matching                                             //
//  Author : Dmitri Peressounko (RRC KI)                              //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include <AliTenderSupply.h>

class TVector3;
class AliPHOSGeometry;
class AliPHOSCalibData ;
class TH2I ;

class AliPHOSTenderSupply: public AliTenderSupply {
  
public:
  AliPHOSTenderSupply();
  AliPHOSTenderSupply(const char *name, const AliTender *tender=NULL);
  virtual ~AliPHOSTenderSupply();

  virtual void   Init(){}
  virtual void   ProcessEvent();
  
  void  SetNonlinearityVersion(const char * ver="Gustavo2005"){fNonlinearityVersion=ver;}
  void  SetNonlinearityParams(Int_t n, const Double_t * par){
            if(n>10){printf("Only 10 parameters allowed \n"); return ;}
            for(Int_t i=0;i<n;i++)fNonlinearityParams[i]=par[i]; }
  void  SetReconstructionPass(Int_t ipass=2){fRecoPass=ipass;}
  
  //If you want to override automatic choise of bad maps and calibration
  void ForceUsingBadMap(const char * filename="alien:///alice/cern.ch/user/p/prsnko/BadMaps/BadMap_LHC10b.root") ;
  void ForceUsingCalibration(const char * filename="alien:///alice/cern.ch/user/p/prsnko/Recalibrations/LHC10b_pass1.root") ;

protected:
  AliPHOSTenderSupply(const AliPHOSTenderSupply&c);
  AliPHOSTenderSupply& operator= (const AliPHOSTenderSupply&c);
  void   InitTender();
  void   FindTrackMatching(Int_t mod,TVector3 *locpos,Double_t &dx, Double_t &dz, Double_t &pttrack, Int_t &charge); 
  Float_t CorrectNonlinearity(Float_t en) ;
  Double_t TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge) ;
  Double_t TestLambda(Double_t pt,Double_t l1,Double_t l2) ;
  Bool_t IsGoodChannel(Int_t mod, Int_t ix, Int_t iz) ;

private:

  TString fOCDBpass ;                        //! Pass to OCDB recalibration object, local or alien
  TString fNonlinearityVersion;              //! Version of non-linearity correction to aaply
  AliPHOSGeometry   *fPHOSGeo;               //! PHOS geometry
  Double_t fNonlinearityParams[10] ;         //! Parameters for non-linearity calculation
  TH2I * fPHOSBadMap[5] ;                    //! Bad channels map
  Int_t fRunNumber ;                         //! Current run number
  Int_t fRecoPass ;                          //! Reconstruction pass
  Bool_t fUsePrivateBadMap ;
  Bool_t fUsePrivateCalib ;
  
  AliPHOSCalibData *fPHOSCalibData;          //! PHOS calibration object

 
  ClassDef(AliPHOSTenderSupply, 1); // PHOS tender task
};


#endif

