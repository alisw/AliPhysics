// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors:                                                               *
 *          for The ALICE HLT Project.                                    *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//  @file   AliHLTTRDTracklet.cxx
//  @author Theodor Rascanu
//  @date   
//  @brief  A datacontainer for tracklets for the HLT. 
// 

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
  fChi2(-1),
  // fExB(-1),
  // fVD(-1),
  // fT0(-1),
  // fS2PRF(-1),
  // fDiffL(-1),
  // fDiffT(-1),
  // fX(-1),
  // fY(-1),
  // fZ(-1),
  // fS2Z(-1),
  fDet(-1),
  fBits(0),
  fCount(0),
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  fSize(sizeof(AliHLTTRDTracklet)-sizeof(fClusters[0])),
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
  fChi2(inTracklet->fChi2),
  // fExB(inTracklet->fExB),
  // fVD(inTracklet->fVD),
  // fT0(inTracklet->fT0),
  // fS2PRF(inTracklet->fS2PRF),
  // fDiffL(inTracklet->fDiffL),
  // fDiffT(inTracklet->fDiffT),
  // fX(inTracklet->fX),
  // fY(inTracklet->fY),
  // fZ(inTracklet->fZ),
  // fS2Z(inTracklet->fS2Z),
  fDet(inTracklet->fDet),
  fBits(0),
  fCount(0),
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  fSize(sizeof(AliHLTTRDTracklet)-sizeof(fClusters[0])),
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
  for (Int_t i=0; i < 2; i++){
    fYref[i]   = inTracklet->fYref[i];
    fZref[i]   = inTracklet->fZref[i];
    fYfit[i]   = inTracklet->fYfit[i];
    fZfit[i]   = inTracklet->fZfit[i];
  }
  fC[0] = inTracklet->GetC();
#ifndef HAVE_NOT_ALITRD_SEEDV1_r39693
  fC[1] = inTracklet->GetC(1);
#endif //HAVE_NOT_ALITRD_SEEDV1_r39693
  for (Int_t i=0; i < 3; i++){
    fPad[i] = inTracklet->fPad[i];
    // fCov[i] = inTracklet->fCov[i];
  }

  // for (Int_t i=0; i < 7; i++){
  //   fRefCov[i] = inTracklet->fRefCov[i];
  // }

  // for (Int_t i=0; i < AliTRDseedV1::kNslices; i++){
  //   fdEdx[i] = inTracklet->fdEdx[i];
  // }
  
  for (Int_t i=0; i < AliPID::kSPECIES; i++){
    fProb[i] = inTracklet->fProb[i];
  }

  fBits = UInt_t(inTracklet->TestBits(-1)) >> 14;

  for (Int_t iTimeBin = 0; iTimeBin < AliTRDseedV1::kNclusters; iTimeBin++){
    AliTRDcluster* trdCluster = inTracklet->GetClusters(iTimeBin);
    if (trdCluster){
      fPos[fCount] = iTimeBin;
      new (&fClusters[fCount]) AliHLTTRDExtCluster(trdCluster);
      fCount++;
      fSize += sizeof(fClusters[0]);
    }
  }  
  //if((void*)&fClusters[fCount]!=(void*)GetEndPointer()){printf("ERRR");return;}
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
  outTracklet->fdX     = fdX;
  outTracklet->fX0     = fX0;
  outTracklet->fS2Y    = fS2Y;
  outTracklet->fPt     = fPt;
  outTracklet->SetC(fC[0]);
#ifndef HAVE_NOT_ALITRD_SEEDV1_r39693
  outTracklet->SetC(fC[1], 1);
#endif //HAVE_NOT_ALITRD_SEEDV1_r39693
  outTracklet->fChi2   = fChi2;
  // outTracklet->fExB    = fExB;
  // outTracklet->fVD     = fVD;
  // outTracklet->fT0     = fT0;
  // outTracklet->fS2PRF  = fS2PRF;
  // outTracklet->fDiffL  = fDiffL;
  // outTracklet->fDiffT  = fDiffT;
  // outTracklet->fX      = fX;
  // outTracklet->fY      = fY;
  // outTracklet->fZ      = fZ;
  // outTracklet->fS2Z    = fS2Z;
  outTracklet->fDet    = fDet;

  for (Int_t i=0; i < 2; i++){
    outTracklet->fYref[i]   = fYref[i];
    outTracklet->fZref[i]   = fZref[i];
    outTracklet->fYfit[i]   = fYfit[i];
    outTracklet->fZfit[i]   = fZfit[i];
  }

  for (Int_t i=0; i < 3; i++){
    outTracklet->fPad[i] = fPad[i];
    // outTracklet->fCov[i] = fCov[i];
  }

  // for (Int_t i=0; i < 7; i++){
  //   outTracklet->fRefCov[i] = fRefCov[i];
  // }

  // for (Int_t i=0; i < AliTRDseedV1::kNslices; i++){
  //   outTracklet->fdEdx[i] = fdEdx[i];
  // }

  for (Int_t i=0; i < AliPID::kSPECIES; i++){
    outTracklet->fProb[i] = fProb[i];
  }

  outTracklet->SetBit(UInt_t(fBits)<<14);

  for(Int_t iCluster=0; iCluster < fCount; iCluster++){
    AliTRDcluster *trdCluster = new AliTRDcluster();
    fClusters[iCluster].ExportTRDCluster(trdCluster);
    trdCluster->SetDetector(fDet);
    outTracklet->fClusters[fPos[iCluster]] = trdCluster;
    outTracklet->fIndexes[fPos[iCluster]] = iCluster;
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
  fC[0] = 0.; fC[1] = 0.; 
  for (Int_t i=0; i < AliPID::kSPECIES; i++)
    fProb[i]=0;
  for (Int_t i=0; i<AliTRDseedV1::kNclusters; i++)
    fPos[i]=0;

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
 * Save tracklet at block position
 */
//============================================================================
AliHLTUInt32_t AliHLTTRDTracklet::SaveAt(AliHLTUInt8_t *const block, const AliTRDseedV1* const inTracklet)
{
  AliHLTUInt32_t size=0;

  memcpy(block,inTracklet,sizeof(AliTRDseedV1));
  size+=sizeof(AliTRDseedV1);

  for(int i=0; i<AliTRDseedV1::kNclusters; i++){
    AliTRDcluster* inClust = inTracklet->GetClusters(i);
    if(inClust) size+=AliHLTTRDCluster::SaveAt(block+size, inClust);
  }

  return size;
}

/**
 * Read tracklet from block
 */
//============================================================================
AliHLTUInt32_t AliHLTTRDTracklet::LoadFrom(AliTRDseedV1 *const outTracklet, const AliHLTUInt8_t *const block)
{
  AliHLTUInt32_t size=0;

  memcpy(((AliHLTUInt8_t*)outTracklet)+sizeof(void*),block+sizeof(void*),sizeof(AliTRDseedV1)-sizeof(void*));
  size+=sizeof(AliTRDseedV1);

  for(int i=0; i<AliTRDseedV1::kNclusters; i++){
    if(outTracklet->GetClusters(i)){
      AliTRDcluster *const outClust = new AliTRDcluster;
      outTracklet->fClusters[i]=outClust;
      size+=AliHLTTRDCluster::LoadFrom(outClust, block+size);
    }
  }

  outTracklet->SetBit(AliTRDseedV1::kOwner);

  return size;
}
