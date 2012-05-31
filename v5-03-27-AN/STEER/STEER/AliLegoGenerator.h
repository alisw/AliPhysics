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

class AliLegoGenerator:
public AliGenerator
{

public: 
    AliLegoGenerator();
    
    AliLegoGenerator(Int_t nc1, Float_t c1min, Float_t c1max,
		     Int_t nc2, Float_t c2min, Float_t c2max,
		     Float_t rmin, Float_t rmax, Float_t zmax);
    virtual ~AliLegoGenerator() {}
    virtual void    Generate();
    virtual void    SetCoor1Range(Int_t nbin, Float_t c1min, Float_t c1max)
	{fNCoor1=nbin; fCoor1Min=c1min; fCoor1Max=c1max;}
    virtual Float_t CurCoor1() const {return fCurCoor1;}
    virtual Int_t   Coor1Bin() const {return fCoor1Bin;}
    virtual void    SetCoor2Range(Int_t nbin, Float_t c2min, Float_t c2max)
	{fNCoor2=nbin; fCoor2Min=c2min; fCoor2Max=c2max;}    
    virtual Float_t CurCoor2() const {return fCurCoor2;}
    virtual Int_t   Coor2Bin() const {return fCoor2Bin;}

    virtual void    SetRadiusRange(Float_t rmin, Float_t rmax) 
	{fRadMin=rmin; fRadMax=rmax;}
    virtual void    SetZMax(Float_t zmax) 
	{fZMax=zmax;}
    	    
    virtual Float_t ZMax()   const {return fZMax;}
    virtual Float_t RadMax() const {return fRadMax;}
    virtual Int_t   NCoor1() const {return fNCoor1;}
    virtual Int_t   NCoor2() const {return fNCoor2;}
    virtual void    Coor1Range(Float_t &c1min, Float_t &c1max) const
	{c1min =  fCoor1Min; c1max =  fCoor1Max;}
    virtual void    Coor2Range(Float_t &c2min, Float_t &c2max) const
	{c2min =  fCoor2Min; c2max =  fCoor2Max;}
	    
    Float_t         PropagateCylinder(Float_t *x, Float_t *v, Float_t r, Float_t z);
 protected:
    Float_t    fRadMin;             // Generation radius
    Float_t    fRadMax;             // Maximum tracking radius
    Float_t    fZMax;               // Maximum tracking Z
    Int_t      fNCoor1;             // Number of bins in Coor1
    Int_t      fNCoor2;             // Number of bins in Coor2

    Float_t    fCoor1Min;           // Minimum Coor1
    Float_t    fCoor1Max;           // Maximum Coor1
    Float_t    fCoor2Min;           // Minimum Coor2
    Float_t    fCoor2Max;           // Maximum Coor2 
    
    Int_t      fCoor1Bin;           //Current Coor1 bin
    Int_t      fCoor2Bin;           //Current Coor2 bin
    Float_t    fCurCoor1;           //Current Coor1 of track
    Float_t    fCurCoor2;           //Current c2 of track
    
    ClassDef(AliLegoGenerator,1) //Lego generator
};

#endif

