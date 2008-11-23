// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTTPCEVE.cxx
    @author Matthias Richter
    @date   2008-11-22
    @brief  AliEVE bindings for the HLT TPC.
*/

#include <cerrno>
#include <cassert>
#include "AliHLTTPCEVE.h"
#include "TEvePointSet.h"
#include "TEveElement.h"
#include "AliRawReader.h"
#include "AliHLTOUT.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCSpacePointData.h"
#include "TSystem.h"
#include "TClass.h"
#include "TString.h"
#include "TMath.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCEVE)

AliHLTTPCEVE::AliHLTTPCEVE()
  : AliHLTLogging()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCEVE::~AliHLTTPCEVE()
{
  // see header file for class documentation
}

TEvePointSet* AliHLTTPCEVE::MakePointSetFromHLTDigits(const char* /*path*/, int eventNo, TEveElement* cont, Float_t maxR) const
{
  // see header file for class documentation
  TTree* pDigitTree=NULL;
  AliHLTOUT* pHLTOUT=AliHLTOUT::New(pDigitTree, eventNo);
  TEvePointSet* pointSet=MakePointSetFromHLTOUT(pHLTOUT, cont, maxR);
  AliHLTOUT::Delete(pHLTOUT);
  return pointSet;
}

TEvePointSet* AliHLTTPCEVE::MakePointSetFromHLTOUT(AliRawReader* pRawReader, TEveElement* cont, Float_t maxR) const
{
  // see header file for class documentation
  if (!pRawReader) return NULL;
  AliHLTOUT* pHLTOUT=AliHLTOUT::New(pRawReader);
  TEvePointSet* pointSet=MakePointSetFromHLTOUT(pHLTOUT, cont, maxR);
  AliHLTOUT::Delete(pHLTOUT);
  return pointSet;
}

TEvePointSet* AliHLTTPCEVE::MakePointSetFromHLTOUT(AliHLTOUT* pHLTOUT, TEveElement* /*cont*/, Float_t maxR) const
{
  // see header file for class documentation
  if (!pHLTOUT) return NULL;

  const Int_t kMaxCl=100*160;
  TEvePointSet* clusters = new TEvePointSet(kMaxCl);
  if (clusters) {
    clusters->SetOwnIds(kTRUE);
    clusters->SetMarkerColor(2);
    clusters->SetMarkerStyle(5);
    
    if (pHLTOUT->Init()>=0) {
      for (int idx=pHLTOUT->SelectFirstDataBlock(AliHLTTPCDefinitions::fgkClustersDataType);
	   idx>=0;
	   idx=pHLTOUT->SelectNextDataBlock()) {
	const AliHLTUInt8_t* pBuffer=NULL;
	AliHLTUInt32_t size=0;
	AliHLTComponentDataType dt=kAliHLTVoidDataType;
	AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
	pHLTOUT->GetDataBlockDescription(dt, spec);

	if (pHLTOUT->GetDataBuffer(pBuffer, size)>=0) {
	  int slice=AliHLTTPCDefinitions::GetMinSliceNr(spec);
	  if (slice!=AliHLTTPCDefinitions::GetMaxSliceNr(spec)) {
	    HLTWarning("cluster data array of multiple TPC slices, can not unambiguously determine phi; skipping data block of specification 0x%08x", spec);
	    continue;
	  }
	  if (size>0 && AddClusters(clusters, reinterpret_cast<const AliHLTTPCClusterData*>(pBuffer), size, slice, maxR)<0) {
	    // action if failed
	  }
	  pHLTOUT->ReleaseDataBuffer(pBuffer);
	}
      }
    } else {
      HLTError("initialization of HLTOUT handler failed");
    }
    TString name="HLT TPC Clusters";
    clusters->SetName(name);
    name.Form("N=%d", clusters->Size());
    clusters->SetTitle(name);
  }

  return clusters;
}

int AliHLTTPCEVE::AddClusters(TEvePointSet* clusters, const AliHLTTPCClusterData* data, unsigned int sizeInByte, int slice, Float_t maxR) const
{
  // see header file for class documentation
  int iResult=0;
  if (!clusters || !data) return -EINVAL;
  Float_t phi     = ( slice + 0.5 ) * TMath::Pi() / 9.0;
  Float_t cos     = TMath::Cos( phi );
  Float_t sin     = TMath::Sin( phi );
  Float_t maxRsqr = maxR*maxR; // maxR squared for a simple geometrical cut below

  for (iResult=0; iResult<(int)data->fSpacePointCnt && iResult>=0; iResult++) {
    if (reinterpret_cast<const AliHLTUInt8_t*>(data->fSpacePoints)+(iResult+1)*sizeof(AliHLTTPCSpacePointData)>reinterpret_cast<const AliHLTUInt8_t*>(data)+sizeInByte) {
      HLTError("data missmatch: buffer of size %d does not match size of AliHLTTPCClusterData (%d) + %d tracks x AliHLTTPCSpacePointData (%d)", 
	       sizeInByte, sizeof(AliHLTTPCClusterData), data->fSpacePointCnt, sizeof(AliHLTTPCSpacePointData));
      iResult=-ENOMSG;
      break;
    }
    const AliHLTTPCSpacePointData* sp=&data->fSpacePoints[iResult];
    if (sp->fX*sp->fX+sp->fY*sp->fY<=maxRsqr) {
      clusters->SetNextPoint(cos*sp->fX - sin*sp->fY, sin*sp->fX + cos*sp->fY, sp->fZ);
    }
  }

  if (iResult<0) {
    // reset clusters
  }

  return iResult;
}
