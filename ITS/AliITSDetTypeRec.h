#ifndef ALIITSDETTYPEREC_H
#define ALIITSDETTYPEREC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
$Id$ 
*/

/**********************************************************************
 * This class contains all of the "external" information needed to do *
 * detector specific reconstruction for the ITS.                      *
 **********************************************************************/
#include <TObject.h>
#include <TObjArray.h>
#include <TClonesArray.h>

#include "AliITSgeom.h"

class AliITSClusterFinder;
class AliITSsegmentation;
class AliITSresponce;

class AliITSDetTypeRec : public TObject {
  public:
    AliITSDetTypeRec(); // Default constructor
    virtual ~AliITSDetTypeRec(); // Proper Destructor

  private:
    AliITSgeom   *fGeom;          //
    TObjArray    *fReconstruction;// [NDet]
    TObjArray    *fSegmentation;  // [NDet]
    TObjArray    *fCalibration;   // [NMod]
    TObjArray    *fPreProcess;    // [] e.g. Find Calibration values
    TObjArray    *fPostProcess;   // [] e.g. find primary vertex
    TObjArray    *fClusters;      //! [NMod][NClusters]
    TClonesArray *fDigits;        //! [NMod][NDigits]
    TString       fClusterClassName; // String with Cluster class name
    TString       fDigClassName;     // String with digit class name.
    TString       fRecPointClassName;// String with RecPoint class name

    ClassDef(AliITSDetTypeRec,1) // ITS Reconstruction structure
}
#endif
