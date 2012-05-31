#ifndef AliHMPIDPIDResponse_h
#define AliHMPIDPIDResponse_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliHMPIDPIDResponse                                                          //
//                                                                      //
// HMPID class to perfom pattern recognition based on Hough transfrom   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <TNamed.h>        //base class

class AliESDtrack;

class AliHMPIDPIDResponse : public TNamed 
{
public : 
             AliHMPIDPIDResponse();    //ctor
    virtual ~AliHMPIDPIDResponse() {;} //dtor
    Double_t CosTheta(Float_t *mod, Int_t species);
    Double_t Resolution(Double_t thetaCerTh, AliESDtrack *pTrk);   //Find the sigma for a given ThetaCerTh

//
protected:
  
private:
  AliHMPIDPIDResponse(const AliHMPIDPIDResponse& r);                //dummy copy constructor
  AliHMPIDPIDResponse &operator=(const AliHMPIDPIDResponse& r);     //dummy assignment operator
//
  ClassDef(AliHMPIDPIDResponse,0)
};

#endif // #ifdef AliHMPIDPIDResponse_cxx

