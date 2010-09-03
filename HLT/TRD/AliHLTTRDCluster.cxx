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

//  @file   AliHLTTRDCluster.cxx
//  @author Theodor Rascanu
//  @date   
//  @brief  A datacontainer for clusters fitting component for the HLT. 
// 

#include "AliHLTTRDCluster.h"
#include <cstring>

/**
 * Default Constructor
 */
//============================================================================
AliHLTTRDCluster::AliHLTTRDCluster():
  fSignals(0),
  fPadCol(0),
  fPadRow(0),
  fPadTime(0),
  fBits(0)
{
}

/**
 * Main Constructor
 */
//============================================================================
AliHLTTRDCluster::AliHLTTRDCluster(const AliTRDcluster* const inCluster):
  fSignals(0),
  fPadCol(inCluster->fPadCol),
  fPadRow(inCluster->fPadRow),
  fPadTime(inCluster->fPadTime),
  fBits(0)
{

  fSignals = inCluster->fSignals[2];
  fSignals|= inCluster->fSignals[3]<<10;
  fSignals|= inCluster->fSignals[4]<<21;

  fBits = UInt_t(inCluster->TestBits(-1)) >> 14; 
}

/**
 * Copy data to the output TRDcluster
 */
//============================================================================
void AliHLTTRDCluster::ExportTRDCluster(AliTRDcluster* const outCluster) const
{
  outCluster->fPadCol=fPadCol;
  outCluster->fPadRow=fPadRow;
  outCluster->fPadTime=fPadTime;
  
  outCluster->fSignals[2] = 0x3ff & fSignals;
  outCluster->fSignals[3] = 0x7ff & fSignals>>10;
  outCluster->fSignals[4] = 0x3ff & fSignals>>21;

  for(int i=2; i<5; i++){
    outCluster->fQ+=outCluster->fSignals[i];
  }

  outCluster->SetBit(UInt_t(fBits)<<14);
}


/**
 * Default Constructor
 */
//============================================================================
AliHLTTRDExtCluster::AliHLTTRDExtCluster():
  AliHLTTRDCluster(),
  fX(0),
  fY(0),
  fZ(0)
{
}

/**
 * Main Constructor
 */
//============================================================================
AliHLTTRDExtCluster::AliHLTTRDExtCluster(const AliTRDcluster* const inCluster):
  AliHLTTRDCluster(inCluster),
  fX(inCluster->GetX()),
  fY(inCluster->GetY()),
  fZ(inCluster->GetZ())
{
}


/**
 * Copy data to the output TRDcluster
 */
//============================================================================
void AliHLTTRDExtCluster::ExportTRDCluster(AliTRDcluster* const outCluster) const
{
  AliHLTTRDCluster::ExportTRDCluster(outCluster);
  outCluster->SetX(fX);
  outCluster->SetY(fY);
  outCluster->SetZ(fZ);
}

/**
 * Prints main info about cluster
 */
//============================================================================
void AliHLTTRDExtCluster::Print() const
{
  printf("   --hltCluster-- addr %p; sizeof(*this) %i\n", (void*)this, (int)sizeof(*this));
  printf("     fX %f; fY %f; fZ %f\n",fX,fY,fZ);
}

/**
 * Save cluster at block position
 */
//============================================================================
AliHLTUInt32_t AliHLTTRDCluster::SaveAt(AliHLTUInt8_t *const block, const AliTRDcluster* const inClust)
{
  AliHLTUInt32_t size=0;

  memcpy(block,inClust,sizeof(AliTRDcluster));
  size+=sizeof(AliTRDcluster);

  return size;
}

/**
 * Read cluster from block
 */
//============================================================================
AliHLTUInt32_t AliHLTTRDCluster::LoadFrom(AliTRDcluster *const outClust, const AliHLTUInt8_t *const block)
{
  AliHLTUInt32_t size=0;

  memcpy(((AliHLTUInt8_t*)outClust)+sizeof(void*),block+sizeof(void*),sizeof(AliTRDcluster)-sizeof(void*));
  size+=sizeof(AliTRDcluster);

  return size;
}
