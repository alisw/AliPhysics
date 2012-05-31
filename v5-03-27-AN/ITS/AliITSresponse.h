#ifndef ALIITSRESPONSE_H
#define ALIITSRESPONSE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>
#include <TString.h>

class AliITSsegmentation;
class TF1;
class AliITSgeom;

////////////////////////////////////////////////////
//                                                //
// ITS base response virtual base class           //
//                                                //
////////////////////////////////////////////////////
class AliITSresponse : public TObject {
 public:
 
    AliITSresponse();
    virtual ~AliITSresponse() {;}
    
    virtual void  SetDiffCoeff(Float_t p1, Float_t p2) {
      fDiffCoeff=p1; fDiffCoeff1=p2;}
    virtual void  DiffCoeff(Float_t &diff,Float_t &diff1) const {
      diff=fDiffCoeff; diff1=fDiffCoeff1;}


 protected:

    void NotImplemented(const char *method) const {if(gDebug>0)
         Warning(method,"This method is not implemented for this sub-class");}
   
 private:
    Float_t  fDiffCoeff;      // Diffusion Coefficient (scaling the time)
    Float_t  fDiffCoeff1;     // Diffusion Coefficient (constant term)
    

    ClassDef(AliITSresponse,5) // Detector type response virtual base class 
};

#endif
