#ifndef AliHMPIDPid_h
#define AliHMPIDPid_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliHMPIDPid                                                          //
//                                                                      //
// HMPID class to perfom pattern recognition based on Hough transfrom   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <TTask.h>        //base class

class AliESDtrack;

class AliHMPIDPid : public TTask 
{
public : 
             AliHMPIDPid();    //ctor
    virtual ~AliHMPIDPid() {;} //dtor
    
    void FindPid(AliESDtrack *pESD,Int_t nsp,Double_t *prob);  //Find PID for tracks

//
protected:
  
private:
  AliHMPIDPid(const AliHMPIDPid& r);                //dummy copy constructor
  AliHMPIDPid &operator=(const AliHMPIDPid& r);     //dummy assignment operator
//
  ClassDef(AliHMPIDPid,0)
};

#endif // #ifdef AliHMPIDPid_cxx

