#ifndef ALITOFTRACKER_H
#define ALITOFTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliTracker.h"
#include "AliTOFpidESD.h"

class AliTOFGeometry;


class AliTOFtracker : public AliTracker {
public:

  AliTOFtracker(AliTOFGeometry* geom, Double_t parPID[2]) 
    {fGeom = geom; fTOFpid = new AliTOFpidESD(parPID);};
  virtual ~AliTOFtracker() {delete fTOFpid;}
  virtual Int_t Clusters2Tracks(AliESD* /*event*/) {return -1;};
  virtual Int_t PropagateBack(AliESD* event) {return fTOFpid->MakePID(event);};
  virtual Int_t RefitInward(AliESD* /*event*/) {return -1;};

  virtual Int_t LoadClusters(TTree* tree) 
    {return fTOFpid->LoadClusters(tree, fGeom);};
  virtual void UnloadClusters() {fTOFpid->UnloadClusters();};
  virtual AliCluster *GetCluster(Int_t index) const {return NULL;};

private:
  AliTOFpidESD*    fTOFpid;  // TOF PID
  AliTOFGeometry*  fGeom;    // TOF geometry

  ClassDef(AliTOFtracker, 1) // TOF tracker (wrapper for AliTOFpidESD)
};

#endif


