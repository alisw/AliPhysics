/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
Revision 1.2  1999/09/29 09:24:19  fca
Introduction of the Copyright and cvs Log

*/

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
