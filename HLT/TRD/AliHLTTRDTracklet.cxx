#include "AliHLTTRDTracklet.h"
/**
 * Default Constructor
 */
//============================================================================
AliHLTTRDTracklet::AliHLTTRDTracklet():
  fN(0),
  fdX(-1),
  fS2Y(-1),
  fPt(-1),
  fX0(-1),
  fC(-1),
  fChi2(-1),
  fDet(-1),
  fCount(0),
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  fSize(sizeof(AliHLTTRDTracklet)-sizeof(AliHLTTRDCluster)),
#else
  fSize(sizeof(AliHLTTRDTracklet))
#endif
{
  InitArrays();
}

/**
 * Main Constructor
 */
//============================================================================
AliHLTTRDTracklet::AliHLTTRDTracklet(const AliTRDseedV1* const inTracklet):
  fN(inTracklet->fN),
  fdX(inTracklet->fdX),
  fS2Y(inTracklet->fS2Y),
  fPt(inTracklet->fPt),
  fX0(inTracklet->fX0),
  fC(inTracklet->fC),
  fChi2(inTracklet->fChi2),
  fDet(inTracklet->fDet),
  fCount(0),
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  fSize(sizeof(AliHLTTRDTracklet)-sizeof(AliHLTTRDCluster)),
#else
  fSize(sizeof(AliHLTTRDTracklet))
#endif
{
  CopyDataMembers(inTracklet);
}

/**
 * Copy simple (non-pointer) data members from TRDTracklet to HLTTRDTracklet 
 */
//============================================================================  
void AliHLTTRDTracklet::CopyDataMembers(const AliTRDseedV1* const inTracklet)
{
  //fChi2Z = inTracklet->GetChi2Z();

  for (Int_t i=0; i < 2; i++){
    fYref[i]   = inTracklet->fYref[i];
    fZref[i]   = inTracklet->fZref[i];
    fYfit[i]   = inTracklet->fYfit[i];
    fZfit[i]   = inTracklet->fZfit[i];
  }

  for (Int_t i=0; i < 3; i++){
    fPad[i] = inTracklet->fPad[i];
  }
  
  for (Int_t i=0; i < AliPID::kSPECIES; i++){
    fProb[i] = inTracklet->fProb[i];
  }

  for (Int_t iTimeBin = 0; iTimeBin < AliTRDseedV1::kNclusters; iTimeBin++)
    {
      AliTRDcluster* trdCluster = inTracklet->GetClusters(iTimeBin);
      if (trdCluster){
  	new (&fClusters[fCount]) AliHLTTRDCluster(trdCluster);
  	fCount++;
  	fSize += sizeof(AliHLTTRDCluster);
      }
    }
  
  //if((void*)&fClusters[fCount].Index!=(void*)GetEndPointer()){printf("ERRR");return;}
}

/**
 * Copy data to the output TRDseedV1
 */
//============================================================================
void AliHLTTRDTracklet::ExportTRDTracklet(AliTRDseedV1* const outTracklet) const
{
  //outTracklet->Reset(); we always use a fresh trdtracklet as input, so this is useless
  outTracklet->SetBit(AliTRDseedV1::kOwner);

  outTracklet->fN      = fN;
  outTracklet->fDet    = fDet;
  outTracklet->fdX     = fdX;
  outTracklet->fX0     = fX0;
  outTracklet->fS2Y    = fS2Y;
  outTracklet->fPt     = fPt;
  outTracklet->fC      = fC;
  outTracklet->fChi2   = fChi2;

  for (Int_t i=0; i < 2; i++){
    outTracklet->fYref[i]   = fYref[i];
    outTracklet->fZref[i]   = fZref[i];
    outTracklet->fYfit[i]   = fYfit[i];
    outTracklet->fZfit[i]   = fZfit[i];
  }

  for (Int_t i=0; i < 3; i++){
    outTracklet->fPad[i] = fPad[i];
  }

  for (Int_t i=0; i < AliPID::kSPECIES; i++){
    outTracklet->fProb[i] = fProb[i];
  }

  for (UInt_t iCluster=0; iCluster < fCount; iCluster++){
    AliTRDcluster *trdCluster = new AliTRDcluster();
    fClusters[iCluster].ExportTRDCluster(trdCluster);
    outTracklet->fClusters[iCluster] = trdCluster;
    outTracklet->fIndexes[iCluster] = iCluster;
  }
}


/**
 * Init arrays
 */
//============================================================================
void AliHLTTRDTracklet::InitArrays()
{
  for (Int_t i=0; i < 2; i++){
    fYref[i] = -1;
    fZref[i] = -1;
    fYfit[i] = -1;
    fZfit[i] = -1;
  }
  for (Int_t i=0; i < AliPID::kSPECIES; i++)
    fProb[i]=0;
}

/**
 * Prints main info about tracklet
 */
//============================================================================
void AliHLTTRDTracklet::Print(Bool_t printClusters) const
{
  //printf("--hltTracklet-- addr 0x%p(%i); fSize %i\n", this, (int)this, fSize);
  printf("      fDet %i; fPt %f; fdX %f fN %i\n", fDet, fPt, fdX, fN);

  if(!printClusters) return;
  for (UInt_t iCluster = 0; iCluster < fCount; iCluster++){
    printf(" [%i] ",iCluster);
    fClusters[iCluster].Print();
  }
}

/**
 * Read clusters to TRDtracklet from the memory
 */
//============================================================================
// void AliHLTTRDTracklet::ReadClustersFromMemory(void *input)
// {
//   AliHLTUInt8_t *iterPtr = (AliHLTUInt8_t*) input;
//   AliHLTTRDCluster* hltCluster = NULL;
  
//   for (Int_t iCluster = 0; iCluster < AliTRDseedV1::kNclusters; iCluster++){
//     // if we had something in the fClusters[iCluster] before copying,
//     // then this entry in the array should not be empty. Fill it.
//     if (fClusters[iCluster]){
//       hltCluster = (AliHLTTRDCluster*) iterPtr;
//       fClusters[iCluster] = hltCluster;
//       iterPtr += hltCluster->GetSize();
//       //hltCluster->Print();
//     }
    
//   }
// }

