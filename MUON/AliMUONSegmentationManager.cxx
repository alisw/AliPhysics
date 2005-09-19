/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpeateose. It is      *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/* $Id$ */

#include "AliMUONSegmentationManager.h"

#include "AliMpFiles.h"
#include "AliLog.h"
#include "AliMpSectorReader.h"
#include "AliMpSector.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpSlatSegmentation.h"
#include "AliMpSt345Reader.h"

#include "TObjString.h"
#include "TClass.h"

ClassImp(AliMUONSegmentationManager)

TExMap AliMUONSegmentationManager::fgMap;
TExMap AliMUONSegmentationManager::fgDetElemIdToSlatTypeMap;

//_____________________________________________________________________________
AliMUONSegmentationManager::AliMUONSegmentationManager() : TObject()
{
}

//_____________________________________________________________________________
AliMUONSegmentationManager::~AliMUONSegmentationManager()
{
}

//_____________________________________________________________________________
Bool_t 
AliMUONSegmentationManager::IsValidDetElemId(Int_t detElemId)
{
  return (SlatType(detElemId) != 0);
}

//_____________________________________________________________________________
AliMpVSegmentation* 
AliMUONSegmentationManager::ReadSegmentation(Int_t detElemId, AliMpPlaneType planeType)
{
  if ( detElemId >= 500 && detElemId <= 1025 )
  {
    AliMpSlat* slat = ReadSlat(detElemId,planeType);
    return new AliMpSlatSegmentation(slat);
  }	
  else if ( detElemId < 500 )
  {
    AliMpStationType station(kStation2);
    if ( detElemId < 200 )
    { 
      station = kStation1;
    }	
    AliMpSectorReader reader(station,planeType);
    AliMpSector* sector = reader.BuildSector();
    //FIXME: get this to be able to delete the sectors:		 fStore.push_back(sector);
    return new AliMpSectorSegmentation(sector);
  }
  else
  {
    return 0x0;
  }
}

//_____________________________________________________________________________
Bool_t
AliMUONSegmentationManager::ReadDetElemIdToSlatType()
{ 
  std::ifstream in(AliMpFiles::Instance()->DetElemIdToSlatTypeFilePath().Data());
  if (!in.good()) return false;
  
  fgDetElemIdToSlatTypeMap.Delete();
  
  char line[80];
  
  while ( in.getline(line,80) )
  {    
    if ( !isdigit(line[0]) ) continue;
    TString sline(line);
    
    Ssiz_t pos = sline.First(' ');
    int detelemid = TString(sline(0,pos)).Atoi();
    fgDetElemIdToSlatTypeMap.Add((Long_t)detelemid,
                                 (Long_t)(new TObjString(sline(pos+1,sline.Length()-pos).Data())));
  }
  
  in.close();
  
  return true;
}

//_____________________________________________________________________________
AliMpSlat*
AliMUONSegmentationManager::ReadSlat(Int_t detElemId, AliMpPlaneType planeType)
{
	return AliMpSt345Reader::ReadSlat(SlatType(detElemId),planeType);
}


//_____________________________________________________________________________
AliMpVSegmentation* 
AliMUONSegmentationManager::Segmentation(Int_t detElemId, AliMpPlaneType planeType)
{
	Long_t it = fgMap.GetValue((Long_t)detElemId);
	
  if ( it )
  {
    TPair* p = (TPair*)(it);
  
    if ( planeType == kBendingPlane )
    {
      return (AliMpVSegmentation*)p->Key();
    }
    else if ( planeType == kNonBendingPlane )
    {
      return (AliMpVSegmentation*)p->Value();
    }		
    else
    {
      AliFatalClass("oups");
      return 0x0;
    }
  }
  else
  {
    AliMpVSegmentation* b = ReadSegmentation(detElemId,kBendingPlane);
    AliMpVSegmentation* nb = ReadSegmentation(detElemId,kNonBendingPlane);
    if ( !b || !nb )
    {
      AliErrorClass(Form("Could not get segmentations for detElemId=%d",detElemId));
      return 0x0;
    } 
    fgMap.Add((Long_t)(detElemId),(Long_t)(new TPair(b,nb)));
    return Segmentation(detElemId,planeType);
  }
}

//_____________________________________________________________________________
const char* 
AliMUONSegmentationManager::SlatType(int detelemid)
{
  if ( ! fgDetElemIdToSlatTypeMap.GetSize() ) ReadDetElemIdToSlatType();
  
  Long_t rv = fgDetElemIdToSlatTypeMap.GetValue(detelemid);
  
  if ( rv )
  {
    return ((TObjString*)(rv))->String().Data();
  }
  else
  {
    return 0;
  }
}

