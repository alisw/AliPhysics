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
Revision 1.5  2003/02/01 14:02:20  nilsen
Work continues.

Revision 1.4  2003/01/29 16:01:14  nilsen
Update today's work.

Revision 1.3  2003/01/28 17:59:54  nilsen
Work continuing.

Revision 1.2  2003/01/26 14:35:15  nilsen
Some more geometry interface functions added and a start at the SSD support
cone geometry. Committed to allow easy updates of partical work between authors.

Revision 1.1  2003/01/20 23:32:49  nilsen
New ITS geometry. Only a Skeleton for now.

$Id$
*/

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  Inner Traking System version 11                                         //
//  This class contains the base procedures for the Inner Tracking System   //
//                                                                          //
// Authors: R. Barbera                                                      //
// version 6.                                                               //
// Created  2000.                                                           //
//                                                                          //
//  NOTE: THIS IS THE  SYMMETRIC PPR geometry of the ITS.                   //
// THIS WILL NOT WORK                                                       //
// with the geometry or module classes or any analysis classes. You are     //
// strongly encouraged to uses AliITSv5.                                    //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
// See AliITSv11::StepManager().
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


#include "AliRun.h"
#include "AliMagF.h"
#include "AliConst.h"
#include "AliITSGeant3Geometry.h"
#include "AliITShit.h"
#include "AliITS.h"
#include "AliITSv11.h"
#include "AliITSgeom.h"
#include "AliITSgeomSPD.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSSD.h"
#include "AliITSDetType.h"
#include "AliITSresponseSPD.h"
#include "AliITSresponseSDD.h"
#include "AliITSresponseSSD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSsimulationSPD.h"
#include "AliITSsimulationSDD.h"
#include "AliITSsimulationSSD.h"
#include "AliITSClusterFinderSPD.h"
#include "AliITSClusterFinderSDD.h"
#include "AliITSClusterFinderSSD.h"
//
#include "AliITSGeometrySSDCone.h"


ClassImp(AliITSv11)

//______________________________________________________________________
AliITSv11::AliITSv11() : AliITS() {
    // Standard default constructor for the ITS version 11.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   A default constructed AliITSv11 class.

    fc = 0;
}
//______________________________________________________________________
AliITSv11::AliITSv11(const char *title) : AliITS("ITS", title){
    // Standard constructor for the ITS version 11.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   A Standard constructed AliITSv11 class.

    fc = 0;
}
//______________________________________________________________________
AliITSv11::~AliITSv11() {
    // Standard destructor for the ITS version 11.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   none.

    if(fc!=0) delete fc;
}
//______________________________________________________________________
void AliITSv11::BuildGeometry(){
    // This routine defines and Creates the geometry for version 11 of the ITS
    // for use in the simulation display routines. This is a very simplified
    // geometry for speed of viewing.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   none.

    if(fc==0) fc = new AliITSGeometrySSDCone(new TVector3(0.0,0.0,0.0),"TSV",0);

    fc->BuildDisplayGeometry();
}
//______________________________________________________________________
void AliITSv11::CreateGeometry(){
    // This routine defines and Creates the geometry for version 11 of the ITS.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   none.

    if(fc==0) fc = new AliITSGeometrySSDCone(new TVector3(0.0,0.0,0.0),"TSV",0);
    TVector3 t(0.0,0.0,0.0);
    fc->CreateG3Geometry(t,"ITSV",0);
}
//______________________________________________________________________
void AliITSv11::CreateMaterials(){
    // Create ITS materials
    //     This function defines the default materials used in the Geant
    // Monte Carlo simulations for the geometries AliITSv1, AliITSv3,
    // AliITSv11.
    // In general it is automatically replaced by
    // the CreatMaterials routine defined in AliITSv?. Should the function
    // CreateMaterials not exist for the geometry version you are using this
    // one is used. See the definition found in AliITSv5 or the other routine
    // for a complete definition.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   none.

    if(fc==0) fc = new AliITSGeometrySSDCone(new TVector3(0.0,0.0,0.0),"TSV",0);

    fc->CreateG3Materials();
}
//______________________________________________________________________
void AliITSv11::InitAliITSgeom(){
    // Based on the geometry tree defined in Geant 3.21, this
    // routine initilizes the Class AliITSgeom from the Geant 3.21 ITS geometry
    // sturture.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   none.
}
//______________________________________________________________________
void AliITSv11::Init(){
    // Initialise the ITS after it has been created.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   none.
}
//______________________________________________________________________
void AliITSv11::SetDefaults(){
    // sets the default segmentation, response, digit and raw cluster classes
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   none.
}
//______________________________________________________________________
void AliITSv11::DrawModule(){
    // Draw a shaded view of the ITS version 11.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return
    //   none.
}
//______________________________________________________________________
void AliITSv11::StepManager(){
    // Called for every step in the ITS, then calles the AliITShit class
    // creator with the information to be recoreded about that hit.
    //  The value of the macro ALIITSPRINTGEOM if set to 1 will allow the
    // printing of information to a file which can be used to create a .det
    // file read in by the routine CreateGeometry(). If set to 0 or any other
    // value except 1, the default behavior, then no such file is created nor
    // it the extra variables and the like used in the printing allocated.
}
 
