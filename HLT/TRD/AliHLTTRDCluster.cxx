#include "AliHLTTRDCluster.h"

/**
 * Default Constructor
 */
//============================================================================
AliHLTTRDCluster::AliHLTTRDCluster():
  fX(0),
  fY(0),
  fZ(0),
  fQ(0),
  fIsInChamber(kFALSE),
  fIsShared(kFALSE),
  fDetector(-1),
  fLocalTimeBin(0),
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
  fX (inCluster->GetX()),
  fY (inCluster->GetY()),
  fZ (inCluster->GetZ()),
  fQ (inCluster->GetQ()),
  fIsInChamber(inCluster->IsInChamber()),
  fIsShared (inCluster->IsShared()),
  fDetector (inCluster->GetDetector()),
  fLocalTimeBin (inCluster->GetLocalTimeBin()),
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
  outCluster->SetQ(fQ);
  outCluster->SetInChamber(fIsInChamber);
  outCluster->SetShared(fIsShared);
  outCluster->SetDetector(fDetector);
  outCluster->SetLocalTimeBin(fLocalTimeBin);
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
  //printf("   --hltCluster-- addr 0x%x(%i); sizeof(*this) %i\n", this, (int)this, this->GetSize());
  //printf("     fX %f; fY %f; fZ %f\n",fX,fY,fZ);
  
}
