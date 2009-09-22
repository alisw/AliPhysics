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
  fDetector(-1),
  fLocalTimeBin(0),
  fClusterMasking(0),
  fPadCol(0),
  fPadRow(0),
  fPadTime(0),
  fIsInChamber(kFALSE),
  fIsShared(kFALSE)
{
}

/**
 * Main Constructor
 */
//============================================================================
AliHLTTRDCluster::AliHLTTRDCluster(const AliTRDcluster* const inCluster):
  fX (inCluster->GetX()),
  fY (inCluster->GetY()),
  fZ (inCluster->GetZ()),
  fQ (inCluster->fQ),
  fDetector (inCluster->fDetector),
  fLocalTimeBin (inCluster->fLocalTimeBin),
  fClusterMasking (inCluster->fClusterMasking),
  fPadCol (inCluster->fPadCol),
  fPadRow (inCluster->fPadRow),
  fPadTime (inCluster->fPadTime),
  fIsInChamber(inCluster->IsInChamber()),
  fIsShared (inCluster->IsShared())
{
  for(int i=0; i<3; i++)
    fSignals[i]=inCluster->fSignals[i+2];
}


/**
 * Copy data to the output TRDcluster
 */
//============================================================================
void AliHLTTRDCluster::ExportTRDCluster(AliTRDcluster* const outCluster) const
{
  //  Print();
  outCluster->SetX(fX);
  outCluster->SetY(fY);
  outCluster->SetZ(fZ);
  outCluster->fQ=fQ;
  outCluster->SetInChamber(fIsInChamber);
  outCluster->SetShared(fIsShared);
  outCluster->fDetector=fDetector;
  outCluster->fLocalTimeBin=fLocalTimeBin;
  outCluster->fClusterMasking=fClusterMasking;
  outCluster->fPadCol=fPadCol;
  outCluster->fPadRow=fPadRow;
  outCluster->fPadTime=fPadTime;

  for(int i=0; i<3; i++)
    outCluster->fSignals[i+2]=fSignals[i];
}

/**
 * Prints main info about cluster
 */
//============================================================================
void AliHLTTRDCluster::Print() const
{
  printf("   --hltCluster-- addr %p; sizeof(*this) %i\n", (void*)this, (int)sizeof(*this));
  printf("     fX %f; fY %f; fZ %f\n",fX,fY,fZ);
}
