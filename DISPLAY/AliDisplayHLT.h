#ifndef AliDISPLAYHLT_H
#define AliDISPLAYHLT_H

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////
// ALICE HLT DISPLAY CLASS                                             //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#include <Rtypes.h>

class TPolyMarker3D;

class AliDisplayHLT{
  //This classes is an interface to the HLT data
  //For the moment only for TPC, for adding modules there is two choices:
  //1) add the function LoadHLT[module](Int_t) and update the function LoadHLT
  //2) or inherit your class from AliDisplayHLT and overload LoadHLT

 public:

  AliDisplayHLT();
  virtual ~AliDisplayHLT();

  virtual void  LoadHLT(const char *name,Int_t e);//Load L3 datas whose belong to detector name and from the event e
  virtual void  LoadHLTTPC(Int_t nevent);
  virtual void  Draw();

 private:
  TPolyMarker3D *fPoints; //fPoints[i]=set of cluster coordinates in detector i;
  Int_t         fNb; // Number of HLT clusters
  char          **fName; //fName[i]=name of the detector i 

 ClassDef(AliDisplayHLT,0);
};

#endif
