#include "AliHLTTRDTracklet.h"
/**
 * Default Constructor
 */
//============================================================================
AliHLTTRDTracklet::AliHLTTRDTracklet():
  fN(0),
  fDet(-1),
  fdX(-1),
  fS2Y(-1),
  fPt(-1),
  fPad3(-1),
  fPad2(-1),
  fX0(-1),
  fC(-1),
  fChi2(-1),
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  fSize(sizeof(AliHLTTRDTracklet)-sizeof(IndexAndCluster)),
#else
  fSize(sizeof(AliHLTTRDTracklet)),
#endif
  fCount(0)
{
  InitArrays();
}

/**
 * Main Constructor
 */
//============================================================================
AliHLTTRDTracklet::AliHLTTRDTracklet(const AliTRDseedV1* const inTracklet):
  fN(inTracklet->fN),
  fDet(inTracklet->fDet),
  fdX(inTracklet->fdX),
  fS2Y(inTracklet->fS2Y),
  fPt(inTracklet->fPt),
  fPad3(inTracklet->fPad[3]),
  fPad2(inTracklet->fPad[2]),
  fX0(inTracklet->fX0),
  fC(inTracklet->fC),
  fChi2(inTracklet->fChi2),
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  fSize(sizeof(AliHLTTRDTracklet)-sizeof(IndexAndCluster)),
#else
  fSize(sizeof(AliHLTTRDTracklet)),
#endif
  fCount(0)
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
  
  for (Int_t i=0; i < AliPID::kSPECIES; i++){
    fProb[i] = inTracklet->fProb[i];
  }

  for (Int_t iTimeBin = 0; iTimeBin < AliTRDseedV1::kNclusters; iTimeBin++)
    {
      AliTRDcluster* trdCluster = inTracklet->GetClusters(iTimeBin);
      if (trdCluster){
	fClusters[fCount].Index = inTracklet->fIndexes[iTimeBin];
  	new (&fClusters[fCount].Cluster) AliHLTTRDCluster(trdCluster);
  	fCount++;
  	fSize += sizeof(IndexAndCluster);
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
  outTracklet->fPad[3] = fPad3;
  outTracklet->fPad[2] = fPad2;
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

  for (Int_t i=0; i < AliPID::kSPECIES; i++){
    outTracklet->fProb[i] = fProb[i];
  }

  for (UInt_t iCluster=0; iCluster < fCount; iCluster++){
    AliTRDcluster *trdCluster = new AliTRDcluster();
    fClusters[iCluster].Cluster.ExportTRDCluster(trdCluster);
    outTracklet->fClusters[iCluster] = trdCluster;
    outTracklet->fIndexes[iCluster] = fClusters[iCluster].Index;
    iCluster++;
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
  for (UInt_t iCount=0, iCluster = 0; iCluster < fCount; iCount++){
    printf(" [%i] ",iCount);
    fClusters[iCluster].Cluster.Print();
    iCluster++;
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

