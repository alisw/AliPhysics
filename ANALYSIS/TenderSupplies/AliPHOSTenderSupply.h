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
class AliPHOSTenderSupply: public AliTenderSupply {
  
public:
  AliPHOSTenderSupply();
  AliPHOSTenderSupply(const char *name, const AliTender *tender=NULL);
  virtual ~AliPHOSTenderSupply();

  virtual void   Init();
  virtual void   ProcessEvent();
  
  void  SetNonlinearityVersion(const char * ver="Gustavo2005"){fNonlinearityVersion=ver;}
  void  SetNonlinearityParams(Int_t n, const Double_t * par){
            if(n>10){printf("Only 10 parameters allowed \n"); return ;}
            for(Int_t i=0;i<n;i++)fNonlinearityParams[i]=par[i]; }
  void  SetRecalibOCDBPass(const char * pass="local://OCDB"){fOCDBpass=pass;}

protected:
  AliPHOSTenderSupply(const AliPHOSTenderSupply&c);
  AliPHOSTenderSupply& operator= (const AliPHOSTenderSupply&c);
  void FindTrackMatching(Int_t mod,TVector3 *locpos,Double_t &r,Double_t &dx, Double_t &dz); 
  Float_t CorrectNonlinearity(Float_t en) ;
private:

  TString fOCDBpass ;                        //! Pass to OCDB recalibration object, local or alien
  TString fNonlinearityVersion;              //! Version of non-linearity correction to aaply
  AliPHOSGeometry   *fPHOSGeo;               //! PHOS geometry
  Double_t fNonlinearityParams[10] ;         //! Parameters for non-linearity calculation

  AliPHOSCalibData *fPHOSCalibData;          //! PHOS calibration object

 
  ClassDef(AliPHOSTenderSupply, 1); // PHOS tender task
};


#endif

