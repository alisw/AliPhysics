#ifndef ALILEGOGENERATOR_H
#define ALILEGOGENERATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
//    Utility class to compute and draw Radiation Length Map                 //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliGenerator.h"

class AliLegoGenerator : public AliGenerator
{

public:
  AliLegoGenerator() {}
  AliLegoGenerator(Int_t ntheta, Float_t themin, Float_t themax,
		   Int_t nphi, Float_t phimin, Float_t phimax,
		   Float_t rmin, Float_t rmax, Float_t zmax);
  void    Generate();
  Float_t CurTheta() const {return fCurTheta;}
  Int_t   ThetaBin() const {return fThetaBin;}
  Float_t CurPhi() const {return fCurPhi;}
  Float_t ZMax() const {return fZMax;}
  Float_t RadMax() const {return fRadMax;}
  Int_t   PhiBin() const {return fPhiBin;}
  Int_t   Nphi() const {return fNphi;}
  Int_t   Ntheta() const {return fNtheta;}
  Float_t       PropagateCylinder(Float_t *x, Float_t *v, Float_t r, Float_t z);
private:
   Float_t    fRadMin;          //Generation radius
   Float_t    fRadMax;          //Maximum tracking radius
   Float_t    fZMax;            //Maximum tracking Z
   Int_t      fNtheta;          //Number of bins in Theta
   Int_t      fNphi;            //Number of bins in Phi
   Int_t      fThetaBin;        //Current theta bin
   Int_t      fPhiBin;          //Current phi bin
   Float_t    fCurTheta;        //Current theta of track
   Float_t    fCurPhi;          //Current phi of track

  ClassDef(AliLegoGenerator,1) //Lego generator
};

#endif

