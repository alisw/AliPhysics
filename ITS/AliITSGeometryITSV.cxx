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
$Log:
$Id$
*/

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <TMath.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBE.h>
#include <TTUBS.h>
#include <TPCON.h>
#include <TFile.h>    // only required for Tracking function?
#include <TCanvas.h>
#include <TObjArray.h>
#include <TLorentzVector.h>
#include <TObjString.h>
#include <TClonesArray.h>
#include <TBRIK.h>
#include <TSystem.h>
#include <TVector3.h>
#include <AliRun.h>
#include <AliITS.h>

#include "AliITSGeometryITSV.h"

ClassImp(AliITSGeometryITSV)

//______________________________________________________________________
AliITSGeometryITSV::AliITSGeometryITSV() : AliITSBaseGeometry(){
    //Default Constructor for SSD Cone geometry
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.

    SetScalemm();
}
//______________________________________________________________________
AliITSGeometryITSV::AliITSGeometryITSV(AliITS *its,const char *moth):
    AliITSBaseGeometry(its,0){
    //Standard Constructor for SSD Cone geometry
    // Inputs:
    //   const char *moth, The volume name into which the ITS mother volume
    //                     will reside.
    // Outputs:
    //   none.
    // Return:
    //   none.

    fAir = 5; // Air material number.

    fA.Size(16,"ITS Mother Volume");
    fA.SetVid(ITSG3VnameToIndex("ITSV")); // ITSV is a special name.
    Double_t rlim  = 50.;
    Double_t zmax  = 74.;
    Double_t ztpc = 284.;
    //
    fA.P0() = 0.;
    fA.dP() = 360.;
    //
    fA.Z(0)  = -ztpc-5.0-0.1;
    fA.Rn(0) = 46.0;   
    fA.Rx(0) = 85.0;
    //
    fA.Z(1) = -ztpc;
    fA.Rn(1) = 46;   
    fA.Rx(1) = 85.;
    //
    fA.Z(2) = -ztpc;
    fA.Rn(2) = 46;  
    fA.Rx(2) = rlim+6;
    //
    fA.Z(3) = -97.5;
    fA.Rn(3) = 46;  
    fA.Rx(3) = rlim+6;
    //
    fA.Z(4) = -zmax;
    fA.Rn(4) = 46;  
    fA.Rx(4) = rlim+6;
    //
    fA.Z(5) = -48;   
    fA.Rn(5) = 6;
    fA.Rx(5) = rlim+6;
    //
    fA.Z(6) = -28.6;   
    fA.Rn(6) = 6;
    fA.Rx(6) = rlim+6;
    //
    fA.Z(7) = -27.6;  
    fA.Rn(7) = 3.295;
    fA.Rx(7) = rlim+6; 
    //
    fA.Z(8) = 27.6;
    fA.Rn(8) = 3.295;
    fA.Rx(8) = rlim+6;
    //
    fA.Z(9) = 28.6;
    fA.Rn(9) = 6;
    fA.Rx(9) = rlim+6;
    //
    fA.Z(10) = 48;   
    fA.Rn(10) = 6;
    fA.Rx(10) = rlim+6;
    //
    fA.Z(11) = zmax;
    fA.Rn(11) = 46;
    fA.Rx(11) = rlim+6;
    //
    fA.Z(12) = 97.5;
    fA.Rn(12) = 46;  
    fA.Rx(12) = rlim+6;
    //
    fA.Z(13) = ztpc;
    fA.Rn(13) = 62;     
    fA.Rx(13) = 62+4.;
    //
    fA.Z(14) = ztpc;
    fA.Rn(14) = 62;     
    fA.Rx(14) = 85.;
    //
    fA.Z(15) = ztpc+4.+0.1;
    fA.Rn(15) = 62.4;
    fA.Rx(15) = 85.;
}
//______________________________________________________________________
void AliITSGeometryITSV::CreateG3Geometry(){
    // Calls Geant 3 geometry inilization routines with the information
    // stored in this class.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.

    PolyCone(fA,fAir);
    return;
}
//______________________________________________________________________
void AliITSGeometryITSV::PositionGeometry(const char *moth,Int_t cn,
					  TVector3 &t,Int_t irot){
    // Positions the ITSV geometry in the way needed by Geant3.
    // Inputs:
    //   const char *moth, The volume name into which the ITS mother volume
    //                     will reside.
    //   Int_t copy        copy number.
    //   TVector &t        Translection vector
    //   Int_t   irot      Rotation matrix index number
    // Outputs:
    //    none.
    // Return:
    //    none.

    gMC->Gspos("ITSV",1,moth,t.X(),t.Y(),t.Z(),irot,"ONLY");
    return;
}
//______________________________________________________________________
void AliITSGeometryITSV::CreateG3Materials(){
    // Fills the Geant 3 banks with Material and Medium definisions.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Returns:
    //   none.

    // Dry Air = 78.084% N2, 20.946% O2, 0.934% Ar, 0.033% CO2 by Volume
    // "Handbook of Chemistry and Physics" 60th eddition, CRC-Press (1980)
    // page F-211. Others components <0.003% and measured in parts per 
    // million.
    Int_t i;
    Int_t    Z[4]={6,7,8,18}; //C, N, O, Ar
    Double_t W[4],f[4]={0.00033,0.78084,0.20946,0.00934};//C02,N2,O2,Ar by vol
    Double_t T = 293.15; // 20 degrees C.
    Double_t P = 760.0; // mm of Hg.
    // op. cite, page  F-9.
    Double_t dens = 1.2929*(273.13/T)*(P/760.0); // Density in grams/leter

    dens *= 1.0E-3; // Convert 1 leter = 1E3 cm^3.
    // Covert fraction of compounds by volume to graction of atoms by weight.
    W[0] = GetA(Z[0])*f[0]; // C
    W[1] = GetA(Z[1])*2.0*f[1]; // N
    W[2] = GetA(Z[2])*2.0*(f[0]+f[2]); // O
    W[3] = GetA(Z[3])*f[3];
    // Renormilize the weights.
    f[0] = 0.0;
    for(i=0;i<4;i++) f[0] += W[i];
    for(i=0;i<4;i++) W[i] /= f[0];
    MixtureByWeight(fAir,"Dry Standard ITS Air",Z,W,dens,4,0);
}
//______________________________________________________________________
void AliITSGeometryITSV::BuildDisplayGeometry(){
    // Fill Root geometry banks for fast simple ITS simulation event
    // display. See Display.C, and related code, for more details.
    // Inputs:
    //    none.
    // Outputs:
    //   none.
    // Return:
    //  none.

    // No need to display ITS cones.
}

//______________________________________________________________________
void AliITSGeometryITSV::PolyCone(AliITSPConeData &d,Int_t med){
    // Interface to TMC->Gsvolu() for ITS PCON geometry. Poly-cone It has 9 
    // parameters or more. See SetScale() for units. Default units are geant
    // 3 [cm].
    // Inputs:
    //    AliITSPConeData &d  Object with poly cone data stored in it.
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    Float_t *param;
    Int_t n,i;

    n = 3+3*d.Nz();
    param = new Float_t[n];
    param[0] = d.Phi0();
    param[1] = d.DPhi();
    param[2] = (Float_t) d.Nz();
    for(i=0;i<d.Nz();i++){
	param[3+3*i] = fScale*d.ZAt(i);
	param[4+3*i] = fScale*d.Rmin(i);
	param[5+3*i] = fScale*d.Rmax(i);
    } // end for if
    if(fVolName==0){ // must create array.
	fVolNameSize = 38624;
	fVolName = new TString[fVolNameSize];
	fVolNameLast = 0;
    } // end if
    d.SetVid(ITSG3VnameToIndex("ITSV"));
    fVolName[d.GetVid()] = (d.GetName())->Data();
    gMC->Gsvolu("ITSV","PCON",GetMed(med),param,n);

    delete[] param;
}

