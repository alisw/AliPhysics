#include "AliHLTTRDCluster.h"


/**
 * Default Constructor
 */
//============================================================================
AliHLTTRDCluster::AliHLTTRDCluster():
  fSize(sizeof(AliHLTTRDCluster)),
  fX(0),
  fY(0),
  fZ(0),
  fIsInChamber(kFALSE),
  fIsShared(kFALSE),
  fDetector(-1),
  fLocalTimeBin(0),
  fQ(0),
  fClusterMasking(0),
  fPadCol(0),
  fPadRow(0),
  fPadTime(0)
{
}

/**
 * Main Constructor
 */
//============================================================================
AliHLTTRDCluster::AliHLTTRDCluster(AliTRDcluster * inCluster):
  fSize (sizeof(AliHLTTRDCluster)),
  fX (inCluster->GetX()),
  fY (inCluster->GetY()),
  fZ (inCluster->GetZ()),
  fIsInChamber(inCluster->IsInChamber()),
  fIsShared (inCluster->IsShared()),
  fDetector (inCluster->GetDetector()),
  fLocalTimeBin (inCluster->GetLocalTimeBin()),
  fQ (inCluster->GetQ()),
  fClusterMasking (inCluster->IsMasked()),
  fPadCol (inCluster->GetPadCol()),
  fPadRow (inCluster->GetPadRow()),
  fPadTime ( inCluster->GetPadTime())
{
  //  fNPads = inCluster->GetNPads();
  //  fCenter = inCluster->GetCenter();
}


/**
 * Copy data to the output TRDcluster
 */
//============================================================================
void AliHLTTRDCluster::ExportTRDCluster(AliTRDcluster *outCluster)
{
  //  Print();
  outCluster->SetX(fX);
  outCluster->SetY(fY);
  outCluster->SetZ(fZ);
  outCluster->SetInChamber(fIsInChamber);
  outCluster->SetShared(fIsShared);
  outCluster->SetDetector(fDetector);
  outCluster->SetLocalTimeBin(fLocalTimeBin);
  outCluster->SetQ(fQ);
  outCluster->SetClusterMasking(fClusterMasking);

  outCluster->SetPadCol(fPadCol);
  outCluster->SetPadRow(fPadRow);
  outCluster->SetPadTime(fPadTime);
  //  outCluster->SetNPads(fNPads);
  //  outCluster->SetCenter(fCenter);
  
  
}

/**
 * Prints main info about cluster
 */
//============================================================================
void AliHLTTRDCluster::Print()
{
  //printf("   --hltCluster-- addr 0x%x(%i); fSize %i\n", this, (int)this, this->GetSize());
  //printf("     fX %f; fY %f; fZ %f\n",fX,fY,fZ);
  
}
