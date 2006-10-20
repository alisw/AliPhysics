#ifndef ALIITSpIDESD2_H
#define ALIITSpIDESD2_H
/* Copyright(c) 2005-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-----------------------------------------------------------------------//
// ITS PID class --- method # 2                                          //
//                                                                       //
//                                                                       //
//The PID is based on the likelihood of all the four ITS' layers,        //
//without using the truncated mean for the dE/dx. The response           //
//functions for each layer are convoluted Landau-Gaussian functions.     // 
// Origin: Elena Bruna bruna@to.infn.it, Massimo Masera masera@to.infn.it//
//-----------------------------------------------------------------------//
#include "AliITSpidESD.h"

class AliITStrackerMI;
class AliITSLoader;
class AliITSSteerPid;

class AliITSpidESD2 : public AliITSpidESD {
public:
  AliITSpidESD2();
  AliITSpidESD2(AliITStrackerMI *tracker,AliITSLoader* loader);
  virtual ~AliITSpidESD2();
  virtual Int_t MakePID(AliESD *event);
  AliITSpidESD2(const AliITSpidESD2 &ob); // copy constructor
  AliITSpidESD2& operator=(const AliITSpidESD2 & source); // ass. op.

private:
  AliITStrackerMI *fTracker; //!tracker MI
  AliITSLoader* fLoader;     //!ITS Loader
  AliITSSteerPid* fSp;       //!pointer to AliITSSteerPid

  ClassDef(AliITSpidESD2,1)   // ITS PID class
};

#endif

