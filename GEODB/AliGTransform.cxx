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
*/

/**********************************************/
/*                                            */
/* FILE: AliGTransform.cxx                    */
/* PURPOSE: To define the relative positions  */
/*          of AliGNodes.                     */
/* LANGUAGE: C++                              */
/* COMPILER: CC for HP-UX 9.x and 10.         */
/* AUTHOR: Joana && David                     */
/* DATE: May 28, 1999                         */
/* ADDRESS: jesanto@cern.ch, dcollado@cern.ch */
/*                                            */
/**********************************************/


#include <TMath.h>
#include <TVector.h>
#include <iostream.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include "AliGTransform.h"

ClassImp(AliGTransform)

//----------------------------------------------------------------------

AliGTransform::AliGTransform()
{
    /* Default Constructor */
    fExpression = "";
    fMatrix     = NULL;
    fName       = "";
    fTitle      = "";
    fX          = 0.;
    fY          = 0.;
    fZ          = 0.;
    fTheta	= 0.;
    fPsi	= 0.;
    fPhi	= 0.;
    
}

//----------------------------------------------------------------------

//----------------------------------------------------------------------

AliGTransform::AliGTransform(AliGTransform *tra)
{
    /* Copy Constructor */
    
    fMatrix     = tra->fMatrix;
    fName       = tra->fName;
    fTitle      = tra->fTitle;
    
    
 }   
   
//----------------------------------------------------------------------}

AliGTransform::AliGTransform(Text_t* name, Text_t* title) : TNamed(name,title)
{
    /* Constructor */
    fExpression = "";
    //float matrix[16] = {1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.}
    fMatrix = new TVector(0,15,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1., "END");
    fX          = 0.;
    fY          = 0.;
    fZ          = 0.;
    fTheta	= 0.;
    fPsi	= 0.;
    fPhi	= 0.;
}

//----------------------------------------------------------------------

AliGTransform::AliGTransform(Text_t* name, Text_t* title, Text_t *expression): TNamed(name, title)
{
    /* Constructor */
    fExpression = expression;
    
    //float matrix[16] = {1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.};
    //fMatrix     = new TArrayF(16,matrix);
    fMatrix = new TVector(0,15,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,
    "END");
    CheckExpression();
    BuildMatrix(fX,fY,fZ,fTheta,fPsi,fPhi);
    
}

//----------------------------------------------------------------------

AliGTransform::AliGTransform(Text_t* name,Text_t* title, Text_t *axis, Float_t angle) : TNamed(name, title)
{
    /* Constructor */
    fX=fY=fZ=0.;
   
    
        if (!strcmp(axis,"X")) {  
	    	    
            fTheta	= 90;
            fPsi	= angle;
            fPhi	= -90;
	    }
	if (!strcmp(axis,"Y")) {   
	  	    
            fTheta	= 0.;
            fPsi	= angle;
            fPhi	= 0.;
	     }
	if (!strcmp(axis,"Z")) {  
	  	    
            fTheta	= 0.;
            fPsi	= 0.;
            fPhi	= angle;
	     }

 //cout << "fTheta" << fTheta << endl;
 //cout << "fPsi" << fPsi << endl;
 //cout << "fPhi" << fPhi << endl;

	
	BuildMatrix(fX,fY,fZ,fTheta,fPsi,fPhi);
	     
        
} 
  
//----------------------------------------------------------------------
  
AliGTransform::AliGTransform( Text_t* name, Text_t* title, Float_t theta1,Float_t phi1, 
	                                            Float_t theta2,
						    Float_t phi2, 
						    Float_t theta3,Float_t phi3
	                                            ) : TNamed(name,title)
{
   const Double_t degrad = 0.0174532925199432958;
   float a1=0.,a2=0.,a3=0.,b1=0.,b2=0.,b3=0.,c1=0.,c2=0.,c3=0.;
   fX	= 0;
   fY	= 0;
   fZ	= 0;

   a1 = TMath::Sin(theta1*degrad)*TMath::Cos(phi1*degrad);
   a2 = TMath::Sin(theta1*degrad)*TMath::Sin(phi1*degrad);
   a3 = TMath::Cos(theta1*degrad);
   b1 = TMath::Sin(theta2*degrad)*TMath::Cos(phi2*degrad);
   b2 = TMath::Sin(theta2*degrad)*TMath::Sin(phi2*degrad);
   b3 = TMath::Cos(theta2*degrad);
   c1 = TMath::Sin(theta3*degrad)*TMath::Cos(phi3*degrad);
   c2 = TMath::Sin(theta3*degrad)*TMath::Sin(phi3*degrad);
   c3 = TMath::Cos(theta3*degrad);
   
   // fMatrix = new TVector(0,15,a1,a2,a3,0.,b1,b2,b3,0.,c1,c2,c3,0.,0.,0.,0.,1., "END");
   fMatrix = new TVector(0,15,a1,b1,c1,0.,a2,b2,c2,0.,a3,b3,c3,0.,0.,0.,0.,1., "END");
}						       
//----------------------------------------------------------------------
AliGTransform::AliGTransform( Text_t* name, Text_t* title, Float_t a1,Float_t a2,Float_t a3,Float_t b1,Float_t b2,
	Float_t b3,Float_t c1,Float_t c2,Float_t c3,Float_t Dx,Float_t
	Dy,Float_t Dz) : TNamed(name,title)
{


fMatrix = new TVector(0,15,a1,a2,a3,Dx,b1,b2,b3,Dy,c1,c2,c3,Dz,0.,0.,0.,1., "END");

}

//----------------------------------------------------------------------

AliGTransform::~AliGTransform() {
    /* Destructor */
    if(fMatrix)     delete fMatrix;
}

//----------------------------------------------------------------------

void AliGTransform::CheckExpression() 
/*Extracts the transformation arguments from the expression given in the
constructor*/

{
    TString* string = new TString(fExpression);
    string->ToUpper();
    float Dx, Dy, Dz,theta, psi,phi;
  
    
    if (strstr(*string, "+")) {
         sscanf( *string, "TRA %f %f %f + ROT %f %f %f ", &Dx, &Dy, &Dz, &theta, &psi, &phi );
        
        if( sscanf(*string, "TRA %f %f %f", &Dx, &Dy, &Dz ) == EOF ) 
            printf( "Error! Must introduce 3 distances\n" );

        if( sscanf(*string, "ROT %f %f %f", &theta, &psi, &phi ) == EOF ) 
            printf( "Error! Must introduce 3 angles\n" );
	 fX	= Dx;
   	 fY	= Dy;
   	 fZ	= Dz;
   	 fTheta	= theta;
   	 fPsi	= psi;
    	 fPhi	= phi;
	 
    } else {
  
        if( strstr(*string,"TRA") ) {
            sscanf( *string, "TRA %f %f %f", &Dx, &Dy, &Dz );
            
            if( sscanf(*string, "TRA %f %f %f", &Dx, &Dy, &Dz ) == EOF ) 
                printf( "Error! Must introduce 3 distances\n" );
		
	    fX          = Dx;
   	    fY          = Dy;
   	    fZ          = Dz;
		
	} else {
		
            if( strstr(*string,"ROT") ) {
            sscanf( *string, "ROT %f %f %f", &theta, &psi, &phi );

            if( sscanf(*string, "ROT %f %f %f", &theta, &psi, &phi ) == EOF ) 
                printf( "Error! Must introduce 3 angles\n" );
		
             fTheta	= theta;
   	     fPsi	= psi;
    	     fPhi	= phi;
		
             }
	}
     }
 
 
}

//----------------------------------------------------------------------

void AliGTransform::BuildMatrix(Float_t Dx=0., Float_t Dy=0., Float_t Dz=0., Float_t
theta=0., Float_t psi=0.,Float_t phi=0.  ) 
{
/* Builds the 4X4 matrix of a transformation */

   
    float a1=0.,a2=0.,a3=0.,b1=0.,b2=0.,b3=0.,c1=0.,c2=0.,c3=0.;
    
	const Double_t degrad = 0.0174532925199432958;
	
        a1 =
	TMath::Cos(phi*degrad)*TMath::Cos(psi*degrad)*TMath::Cos(theta*degrad) -
	TMath::Sin(phi*degrad)*TMath::Sin(theta*degrad);
        a2 =
	TMath::Cos(phi*degrad)*TMath::Cos(psi*degrad)*TMath::Sin(theta*degrad) +
	TMath::Sin(phi*degrad)*TMath::Cos(theta*degrad);
        a3 = - TMath::Cos(phi*degrad)*TMath::Sin(psi*degrad);
        b1 = - TMath::Sin(phi*degrad)*TMath::Cos(psi*degrad)*TMath::Cos(theta*degrad) -
	TMath::Cos(phi*degrad)*TMath::Sin(theta*degrad);
        b2 = - TMath::Sin(phi*degrad)*TMath::Cos(psi*degrad)*TMath::Sin(theta*degrad) +
	TMath::Cos(phi*degrad)*TMath::Cos(theta*degrad);
        b3 = TMath::Sin(phi*degrad)*TMath::Sin(psi*degrad);
        c1 = TMath::Sin(psi*degrad)*TMath::Cos(theta*degrad);
        c2 = TMath::Sin(psi*degrad)*TMath::Sin(theta*degrad);
        c3 = TMath::Cos(psi*degrad);
        
	/*
	a1 =
	TMath::Cos(psi*degrad)*TMath::Cos(phi*degrad);
        a2 =
	TMath::Cos(psi*degrad)*TMath::Sin(phi*degrad); 
        a3 = - TMath::Sin(psi*degrad);
        b1 = TMath::Sin(theta*degrad)*TMath::Sin(psi*degrad)*TMath::Cos(phi*degrad) -
	TMath::Cos(theta*degrad)*TMath::Sin(phi*degrad);
        b2 = TMath::Sin(theta*degrad)*TMath::Sin(psi*degrad)*TMath::Sin(phi*degrad) +
	TMath::Cos(theta*degrad)*TMath::Cos(phi*degrad);
        b3 = TMath::Sin(theta*degrad)*TMath::Cos(psi*degrad);
        c1 =
	TMath::Cos(theta*degrad)*TMath::Sin(psi*degrad)*TMath::Cos(phi*degrad) +
	TMath::Sin(theta*degrad)*TMath::Sin(phi*degrad);
        c2 =
	TMath::Cos(theta*degrad)*TMath::Sin(psi*degrad)*TMath::Sin(phi*degrad) -
	TMath::Sin(theta*degrad)*TMath::Cos(phi*degrad);
        c3 = TMath::Cos(theta*degrad)*TMath::Sin(psi*degrad);
	*/
	fMatrix = new TVector(0,15,a1,a2,a3,Dx,b1,b2,b3,Dy,c1,c2,c3,Dz,0.,0.,0.,1., "END");
       
}

//----------------------------------------------------------------------

/*

void AliGTransform::CheckExpression( Text_t* expression ) 
{
    char axis;
    float theta, phi, psi,Dx, Dy, Dz;

    TString* string = new TString( expression );
    string->ToUpper();

    if( strstr(*string,"ROT") ) {
        sscanf( *string, "ROT%c", &axis );

        switch (axis) {
            case 'A':
	        sscanf( *string, "ROTA%f %f %f", &phi, &psi, &theta );
	        if( (!psi) || (!theta) || (!phi) )
                    printf("Error! Must introduce 3 angles\n");
	        break;

            case 'X':
                sscanf( *string, "ROTX%f", &theta );
                if( !theta ) 
                    printf("Error! Must introduce 1 angle\n");
                break;

            case 'Y':
                sscanf( *string, "ROTY%f", &theta);
                if( !theta ) 
                    printf("Error! Must introduce 1 angle\n");
                break;

            case 'Z':
                sscanf( *string, "ROTZ%f", &theta );
                if( !theta )
                    printf("Error! Must introduce 1 angle\n");
                break;

            default:
	        printf("Unrecognised rotation around axis %c\n",axis);
        } 
    }
    else {
        if( strstr(*string, "TRA") ) {
            sscanf( *string, "TRA %f %f %f", &Dx, &Dy, &Dz );
            if( sscanf(*string, "TRA %f %f %f", &Dx, &Dy, &Dz) == EOF )
                printf("Error! Must introduce 3 distances\n");
            printf( "%f%f%f\n", Dx, Dy, Dz );
        }
        else 
            printf( "ERROR!\n" );
    }

    delete string;
}

*/
