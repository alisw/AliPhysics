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

/* $Id$ */

////////////////////////////////////////////////////////////////////////
// This Class owns the classes needed to to detector simulations and
// reconstruction. This includes the detector segmentation classes,
// the detector responce classes, the detector simulatin classes, and
// the detector reconstruction (clustering) classes for all of the ITS
// detectors.
////////////////////////////////////////////////////////////////////////
#include "AliITSDetType.h"
#include "AliITSClusterFinder.h"
#include "AliITSsimulation.h"


ClassImp(AliITSDetType)	 

//______________________________________________________________________
AliITSDetType::AliITSDetType():
TObject(),
fDetType(kND),
fResponse(0),
fSegmentation(0),
fSimulation(0),
fReconst(0),
fDigClassName(""),
fClustClassName(""){
    // Default constructor.
    // Inputs:
    //   none.
    // Output:
    //   none.
    // Return:
    //   A default constructed AliITSDetType class.
}
//______________________________________________________________________
AliITSDetType::AliITSDetType(AliITSDetector det,AliITSresponse *res,
                             AliITSsegmentation *seg,AliITSsimulation *sim,
                             AliITSClusterFinder *cf,
                             const Char_t *DigClassName,
                             const Char_t *ClustClassName):
TObject(),
fDetType(det),
fResponse(res),
fSegmentation(seg),
fSimulation(sim),
fReconst(cf),
fDigClassName(DigClassName),
fClustClassName(ClustClassName){
    // Standard constructor
    // Inputs:
    //   AliITSDetector       det  Detector type
    //   AliITSresponse      *res  response class to use
    //   AliITSsegmentation  *seg  Segmentation class to use
    //   AliITSsimulation    *sim  Simulation class to use
    //   AliITSClusterFinder *cf   Cluster Finder/Reconstruction class to use
    //   const Char_t        DigClassName   Name of the digit class to be used
    //   const Char_t        ClustClassName Name of the cluster class to be 
    //                                      used
    // Output:
    //   none.
    // Return:
    //   A Standard constructed AliITSDetType class.
}
//----------------------------------------------------------------------
AliITSDetType::~AliITSDetType(){
    // destructor
    // Inputs:
    //   none.
    // Output:
    //   none.
    // Return:
    //   none.

    if(fSegmentation!=0) delete fSegmentation; fSegmentation = 0;
    if(fResponse!=0)     delete fResponse;     fResponse     = 0;
    if(fSimulation!=0)   delete fSimulation;   fSimulation   = 0;
    if(fReconst!=0)      delete fReconst;      fReconst      = 0;
}
//______________________________________________________________________
AliITSDetType::AliITSDetType(const AliITSDetType &source) : TObject(source){
    //     Copy Constructor
    // Inputs:
    //   const AliITSDetType &source  class to copy from.
    // Output:
    //   none.
    // Return:
    //   none.

    if(&source == this) return;
    this->fDetType        = source.fDetType;
    this->fReconst        = source.fReconst;
    this->fSimulation     = source.fSimulation;
    this->fResponse       = source.fResponse;
    this->fSegmentation   = source.fSegmentation;
    this->fDigClassName   = source.fDigClassName;
    this->fClustClassName = source.fClustClassName;
    return;
}
//______________________________________________________________________
AliITSDetType& AliITSDetType::operator=(const AliITSDetType &source){
    //    Assignment operator
    // Inputs:
    //   const AliITSDetType &source  class to copy from.
    // Output:
    //   none.
    // Return:
    //   a new AliITSDetType class with the same values as in source.

    if(&source == this) return *this;
    this->fDetType        = source.fDetType;
    this->fReconst        = source.fReconst;
    this->fSimulation     = source.fSimulation;
    this->fResponse       = source.fResponse;
    this->fSegmentation   = source.fSegmentation;
    this->fDigClassName   = source.fDigClassName;
    this->fClustClassName = source.fClustClassName;
    return *this;  
}
