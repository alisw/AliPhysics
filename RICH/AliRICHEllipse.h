#ifndef ALIRICHELLIPSE_H
#define ALIRICHELLIPSE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TPolyMarker3D.h>

class AliRICHEllipse :  public TPolyMarker3D {

 public:
    AliRICHEllipse();
    AliRICHEllipse(Float_t cx, Float_t cy, Float_t omega, Float_t theta, Float_t phi, Float_t emiss);
    
    virtual          ~AliRICHEllipse();
    virtual void CerenkovRingDrawing(Int_t chamber, Int_t track);
 private:
    Float_t fOmega;                    //Cherenkov angle 
    Float_t fTheta;                    //Incidence angle (dip angle)
    Float_t fPhi;                      //Incidence angle 
    Float_t fCx;                       //Hit coordinate-x
    Float_t fCy;                       //Hit coordinate-y
    Float_t fh;                        //Distance from radiator to pads
    Float_t fEmissPoint;               //Emission point
    Float_t fCoordEllipse[2][100];     // Ellipse points to draw
    

    ClassDef(AliRICHEllipse, 1)   //Utility class to draw an ellipse
};
#endif
	
