/***************************************************************************
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
 $Id$
*/

/*********************************************************************
 * This class defines the "Standard" reconstruction for the ITS 
 * detector.
 **********************************************************************/
#include "AliITSDetTypeRec.h"

ClassImp(AliITSDetTypeRec)

//----------------------------------------------------------------------
AliITSDetTypeRec::AliITSDetTypeRec(): TObject(),
fGeom(),        //
fReconstruction(),// [NDet]
fSegmentation(),  // [NDet]
fCalibration(),   // [NMod]
fPreProcess(),    // [] e.g. Find Calibration values
fPostProcess(),   // [] e.g. find primary vertex
fClusters(),      //! [NMod][NClusters]
fDigits(),        //! [NMod][NDigits]
fClusterClassName(), // String with Cluster class name
fDigClassName(),     // String with digit class name.
fRecPointClassName(){// String with RecPoint class name
    // Default Constructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A properly zero-ed AliITSDetTypeRec class.
}
//----------------------------------------------------------------------
AliITSDetTypeRec::~AliITSDetTypeRec(){
    // Destructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    Nothing.
}
//----------------------------------------------------------------------
