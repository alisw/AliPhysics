// -*- C++ -*-
// 
// 1998/10/19
// ---------------------------------------------------------------------------
//
// AliGMaterial Class
//
// This file is part of the ALICE Geometry Database .
//
// Author:  Joana E. Santo
//
// ---------------------------------------------------------------------------

#include "AliGMaterial.h"

ClassImp(AliGMaterial)


//-------------------------------------------------------------------------

AliGMaterial::AliGMaterial( Int_t imat, Text_t* name, Text_t* title, Int_t isvol, Int_t ifield, Float_t fieldm, Float_t tmaxfd, Float_t stemax, Float_t deemax, Float_t epsil, Float_t stmin, Float_t* ubuf, Int_t nbuf, Float_t a, Float_t z, Float_t dens, Float_t radl, Float_t absl, Float_t* buf, Int_t nwbuf ) : TNamed(name, title)
{
    /* VIC: Very Important Constructor */

    fImat   = imat;
    fIsvol  = isvol;
    fIfield = ifield;
    fFieldm = fieldm;
    fTmaxfd = tmaxfd;
    fStemax = stemax;
    fDeemax = deemax;
    fEpsil  = epsil;
    fStmin  = stmin;

    fUbuf   = new Float_t[nbuf];
    
    for( int i=0; i<nbuf; i++ )
        fUbuf[i] = ubuf[i];
        
    fNbuf   = nbuf;
    fA      = a;
    fZ      = z;
    fDens   = dens;
    fRadl   = radl;
    fAbsl   = absl;
    
    fBuf    = new Float_t[nwbuf];

    for( int j=0; j<nwbuf; j++ )
        fBuf[j] = buf[j];
        
    fNwbuf  = nwbuf;
}

//-------------------------------------------------------------------------

AliGMaterial::AliGMaterial(Text_t* name, Text_t* title, Float_t A, Float_t Z, Float_t rho) : TNamed(name, title)
{
    /* Constructor */
    fA   = A;
    fRho = rho;
    fZ   = Z;
}

//-------------------------------------------------------------------------

AliGMaterial::AliGMaterial( AliGMaterial* Mat )
{
    /* Copy Constructor */
    if( Mat ) {
        fImat   = Mat->fImat;
        fIsvol  = Mat->fIsvol;
        fIfield = Mat->fIfield;
        fFieldm = Mat->fFieldm;
        fTmaxfd = Mat->fTmaxfd;
        fStemax = Mat->fStemax;
        fDeemax = Mat->fDeemax;
        fEpsil  = Mat->fEpsil;
        fStmin  = Mat->fStmin;

        fUbuf   = new Float_t[Mat->fNbuf];
    
        for( int i=0; i<Mat->fNbuf; i++ )
            fUbuf[i] = Mat->fUbuf[i];
        
        fNbuf   = Mat->fNbuf;
        fA      = Mat->fA;
        fZ      = Mat->fZ;
        fDens   = Mat->fDens;
        fRadl   = Mat->fRadl;
        fAbsl   = Mat->fAbsl;
    
        fBuf    = new Float_t[Mat->fNwbuf];

        for( int j=0; j<Mat->fNwbuf; j++ )
            fBuf[j] = Mat->fBuf[j];
        
        fNwbuf  = Mat->fNwbuf;

        fName  = Mat->GetName();
        fTitle = Mat->GetTitle();
    }
    else {
        /* Default Constructor */
        fImat   = 0;
        fIsvol  = 0;
        fIfield = 0;
        fFieldm = 0.;
        fTmaxfd = 0.;
        fStemax = 0.;
        fDeemax = 0.;
        fEpsil  = 0.;
        fStmin  = 0.;

        fUbuf   = NULL;
        
        fNbuf   = 0;
        fA      = 0.;
        fZ      = 0.;
        fDens   = 0.;
        fRadl   = 0.;
        fAbsl   = 0.;
    
        fBuf    = NULL;

        fNwbuf  = 0;

        fName  = "";
        fTitle = "";
    }
}

//-------------------------------------------------------------------------

AliGMaterial::~AliGMaterial()
{   /* Destructor */

    if( fUbuf ) delete [] fUbuf;
    if( fBuf  ) delete [] fBuf;
}

//-------------------------------------------------------------------------

AliGMaterial* AliGMaterial::operator=( const AliGMaterial* Mat )
{
    /* Operator = */
    if( this == Mat ) return this; // special case.

    fA     = Mat->fA;
    fName  = Mat->GetName();
    fRho   = Mat->fRho;
    fTitle = Mat->GetTitle();
    fZ     = Mat->fZ;

    return this;
}

//-------------------------------------------------------------------------

void AliGMaterial::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliGMaterial.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> fImat;
      R__b >> fIsvol;
      R__b >> fIfield;
      R__b >> fFieldm;
      R__b >> fTmaxfd;
      R__b >> fStemax;
      R__b >> fDeemax;
      R__b >> fEpsil;
      R__b >> fStmin;
      R__b.ReadArray(fUbuf); //
      R__b >> fNbuf;
      R__b >> fA;
      R__b >> fZ;
      R__b >> fDens;
      R__b >> fRadl;
      R__b >> fAbsl;
      R__b.ReadArray(fBuf); //
      R__b >> fNwbuf;
      R__b >> fRho;
   } else {
      R__b.WriteVersion(AliGMaterial::IsA());
      TNamed::Streamer(R__b);
      R__b << fImat;
      R__b << fIsvol;
      R__b << fIfield;
      R__b << fFieldm;
      R__b << fTmaxfd;
      R__b << fStemax;
      R__b << fDeemax;
      R__b << fEpsil;
      R__b << fStmin;
      R__b.WriteArray(fUbuf, fNbuf); //
      R__b << fNbuf;
      R__b << fA;
      R__b << fZ;
      R__b << fDens;
      R__b << fRadl;
      R__b << fAbsl;
      R__b.WriteArray(fBuf, fNbuf); //
      R__b << fNwbuf;
      R__b << fRho;
   }
}

//-------------------------------------------------------------------------

