/* *************************************************************************
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
Revision 1.1  2000/07/13 16:19:10  fca
Mainly coding conventions + some small bug fixes

Revision 1.8  2000/07/12 08:56:32  fca
Coding convention correction and warning removal

Revision 1.7  2000/06/28 21:27:45  morsch
Most coding rule violations corrected.
Still to do: Split the file (on file per class) ? Avoid the global variables.
Copy constructors and assignment operators (dummy ?)

Revision 1.6  2000/04/14 11:07:46  morsch
Correct volume to medium assignment in case several media are asigned to the
same material.

Revision 1.5  2000/03/20 15:11:03  fca
Mods to make the code compile on HP

Revision 1.4  2000/01/18 16:12:08  morsch
Bug in calculation of number of volume divisions and number of positionings corrected
Browser for Material and Media properties added

Revision 1.3  1999/11/14 14:31:14  fca
Correct small error and remove compilation warnings on HP

Revision 1.2  1999/11/10 16:53:35  fca
The new geometry viewer from A.Morsch

*/

/* 
 *  Version: 0
 *  Written by Andreas Morsch
 *  
 * 
 *
 * For questions critics and suggestions to this part of the code
 * contact andreas.morsch@cern.ch
 * 
 **************************************************************************/


#include "AliGeant3GeometryGUI.h"
#include "AliG3Volume.h"
#include "AliG3Material.h"
#include "AliG3Medium.h"
#include "AliGuiGeomMain.h"
#include "AliRun.h"
#include "AliG3toRoot.h"

#include <TArrayF.h>
#include <TRotMatrix.h>
#include <TGeometry.h>
#include <TFile.h>
#include <TFolder.h>

AliG3Volume    *gCurrentVolume   = new AliG3Volume("NULL");
AliG3Material  *gCurrentMaterial = new AliG3Material();
AliG3Medium    *gCurrentMedium   = new AliG3Medium();

ClassImp(AliGeant3GeometryGUI)

    AliGeant3GeometryGUI::AliGeant3GeometryGUI()
{
// Constructor
    fPanel     = new AliGuiGeomMain(gClient->GetRoot(), 500, 500);
//  Store local copy of zebra bank entries
    AliG3toRoot* geometry = new AliG3toRoot();
//    geometry->SetExpandDivisions();
    geometry->G3toRoot();
    
    AliG3Volume* top = (AliG3Volume*) 
	(geometry->GetTopFolder()->FindObject("ALIC"));
    gCurrentVolume = top;
//
//  Mediate between g3 Geometry and GUI
    fPanel->SetMaterialComboEntries(geometry->GetMaterials());
    fPanel->SetMediaComboEntries(geometry->GetMedia());
    fPanel->AddFoldersRecursively(geometry->GetTopFolder());
    fPanel->Update();
}

void AliGeant3GeometryGUI::Streamer(TBuffer &)
{
// Dummy Streamer
;
}





