#ifndef AliRICHv3_h
#define AliRICHv3_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliRICH.h"

class AliRICHSDigit;

class AliRICHv3 : public AliRICH 
{    
public:
    
                 AliRICHv3():AliRICH()                                {;} 
                 AliRICHv3(const char *pcName, const char *pcTitle);    
  virtual       ~AliRICHv3();                             
  virtual Int_t  IsVersion()                                     const{return 3;}
  virtual void   StepManager();     
         Int_t   Hits2SDigits(Float_t xhit,Float_t yhit,Float_t eloss, Int_t idvol, ResponseType res);
  void     SetGeometryModel(Int_t c,AliRICHGeometry *pRICHGeo)                    {C(c)->SetGeometryModel(pRICHGeo);}
  void     SetSegmentationModel(Int_t c, AliSegmentation *pAliSeg)                {C(c)->SetSegmentationModel(pAliSeg);}
  void     SetResponseModel(Int_t c, AliRICHResponse *pRICHRes)                   {C(c)->SetResponseModel(pRICHRes);}
  AliRICHGeometry* GetGeometryModel(Int_t c=1)                               const{return C(c)->GetGeometryModel();}    
  AliSegmentation* GetSegmentationModel(Int_t c=1)                           const{return C(c)->GetSegmentationModel();}
  AliRICHResponse* GetResponseModel(Int_t c=1)                               const{return C(c)->GetResponseModel();}
private:
  ClassDef(AliRICHv3,1)  //RICH full version, configurable with azimuthal rotation	
};// class AliRICHv3

#endif // AliRICHv3_h
