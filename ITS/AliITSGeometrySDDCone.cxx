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
#include <AliITS.h>

#include "AliITSGeometrySDDCone.h"

ClassImp(AliITSGeometrySDDCone)

//______________________________________________________________________
AliITSGeometrySDDCone::AliITSGeometrySDDCone(){
    //Default Constructor for SDD Cone geometry

    SetScalemm();
}
//______________________________________________________________________
AliITSGeometrySDDCone::AliITSGeometrySDDCone(AliITS *its,TVector3 *&tran,
					     const char moth[3],Int_t mat0):
    AliITSBaseGeometry(its,0){
    //Standard Constructor for SDD Cone geometry
    // Inputs:
    //   Double_t z0  Z-axis shift of this volume
    // Outputs:
    //   none.
    // Return:
    //   none.
}
//______________________________________________________________________
void AliITSGeometrySDDCone::CreateG3Geometry(const char moth[3],
					     TVector3 &trans){
    // Calls Geant 3 geometry inilization routines with the information
    // stored in this class.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.

    return;
}
//______________________________________________________________________
void AliITSGeometrySDDCone::CreateG3Materials(){
    // Fills the Geant 3 banks with Material and Medium definisions.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Returns:
    //   none.
    Int_t Z[5];
    Double_t W[5],dens;

    Z[0] = 1; W[0] = 0.5; // Hydrogen Content
    Z[1] = 6; W[1] = 0.5; // Carbon Content
    MixtureByWeight(fSDDcf,"Carbon Fiber for SDD support cone",Z,W,dens,2,0);
    Z[0] = 1; W[0] = 0.5; // Hydrogen Content
    Z[1] = 6; W[1] = 0.5; // Carbon Content
    MixtureByWeight(fSDDfs,"Inserto stealite 4411w for SDD support cone",
		    Z,W,dens,2,0);
    Z[0] = 1; W[0] = 0.5; // Hydrogen Content
    Z[1] = 6; W[1] = 0.5; // Carbon Content
    MixtureByWeight(fSDDfo,"Foam core (Rohacell 50A) for SDD support cone",
		    Z,W,dens,2,0);
    Z[0] =  6; W[0] = 0.5; // Carbon Content
    Z[1] = 25; W[1] = 0.5; // Iron Content
    MixtureByWeight(fSDDsw,"Stainless steal screw, pin, and stud material",
		    Z,W,dens,2,0);
}
//______________________________________________________________________
void AliITSGeometrySDDCone::BuildDisplayGeometry(){
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

