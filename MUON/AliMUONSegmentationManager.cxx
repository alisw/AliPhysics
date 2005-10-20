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

#include "AliLog.h"
#include "AliMpSectorReader.h"
#include "AliMpSector.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpSlat.h"
#include "AliMpSlatSegmentation.h"
#include "AliMpSt345Reader.h"
#include "AliMpTriggerReader.h"
#include "AliMpTriggerSegmentation.h"
#include "AliMpTrigger.h"

#include "Riostream.h"

#include "TArrayI.h"
#include "TClass.h"
#include "TList.h"
#include "TObjString.h"
#include "TSystem.h"

ClassImp(AliMUONSegmentationManager)

AliMpExMap AliMUONSegmentationManager::fgMap(kTRUE);
AliMpExMap AliMUONSegmentationManager::fgDetElemIdToNameMap(kTRUE);
AliMpExMap AliMUONSegmentationManager::fgLocalBoardMap(kTRUE);

namespace
{
  //__________________________________________________________________________
  TString DetElemIdToNamePath(AliMpStationType stationType)
  {
    /// Get the full path of the file containing the mapping detElemId <->
    /// SlatType.
    /// The bending parameter below is of no use in this case, but
    /// we use it to re-use the PlaneDataDir() method untouched.
    
    TString filename(gSystem->ExpandPathName("${ALICE_ROOT}/MUON/data/"));
    filename += "denames_";
    filename += StationTypeName(stationType);
    filename += ".dat";
    return filename;
  }
}

//_____________________________________________________________________________
AliMUONSegmentationManager::AliMUONSegmentationManager() : TObject()
{
}

//_____________________________________________________________________________
AliMUONSegmentationManager::~AliMUONSegmentationManager()
{
}

//_____________________________________________________________________________
void
AliMUONSegmentationManager::FillLocalBoardMap(AliMpTriggerSegmentation* seg)
{
  const AliMpTrigger* slat = seg->Slat();
  for ( Int_t i = 0; i < slat->GetSize(); ++i )
  {
    TArrayI lbn;
    slat->GetAllLocalBoardNumbers(lbn);
    for ( Int_t j = 0; j < lbn.GetSize(); ++j )
    {
      TList* list = (TList*)fgLocalBoardMap.GetValue(lbn[j]);
      if (!list)
      {
        list = new TList;
        fgLocalBoardMap.Add(lbn[j],list);
      }
      if ( list->FindObject(seg) == 0 )
      {
        list->Add(seg);
      }
    }
  }
}

//_____________________________________________________________________________
Bool_t 
AliMUONSegmentationManager::IsValidDetElemId(Int_t detElemId)
{
  return (DetElemName(detElemId) != 0);
}

//_____________________________________________________________________________
AliMpVSegmentation* 
AliMUONSegmentationManager::ReadSegmentation(Int_t detElemId, AliMpPlaneType planeType)
{
  AliMpStationType station = StationType(detElemId);
  
  if (station==kStation345)
  {
    AliMpSlat* slat = AliMpSt345Reader::ReadSlat(DetElemName(detElemId),planeType);
    return new AliMpSlatSegmentation(slat);
  }
  else if ( station==kStation1 || station==kStation2 )
  {
    AliMpSectorReader reader(station,planeType);
    AliMpSector* sector = reader.BuildSector();
    //FIXME: get this to be able to delete the sectors:		 fStore.push_back(sector);
    return new AliMpSectorSegmentation(sector);
  }
  else if ( station == kStationTrigger )
  {
    AliMpTrigger* slat = AliMpTriggerReader::ReadSlat(DetElemName(detElemId),planeType);
    return new AliMpTriggerSegmentation(slat);
  }
  return 0x0;
}

//_____________________________________________________________________________
Bool_t
AliMUONSegmentationManager::ReadDetElemIdToName(AliMpStationType stationType)
{ 
  std::ifstream in(DetElemIdToNamePath(stationType).Data());
  if (!in.good()) 
  {
    AliErrorClass(Form("Cannot read file %s",DetElemIdToNamePath(stationType).Data()));
    return false;
  }
  
  char line[80];
  
  while ( in.getline(line,80) )
  {    
    if ( !isdigit(line[0]) ) continue;
    TString sline(line);
    
    Ssiz_t pos = sline.First(' ');
    int detelemid = TString(sline(0,pos)).Atoi();
    TObject* o = fgDetElemIdToNameMap.GetValue(detelemid);
    if (!o)
    {
      fgDetElemIdToNameMap.Add(detelemid,
                               new TObjString(sline(pos+1,sline.Length()-pos).Data()));
    }
  }
  
  in.close();
  
  return true;
}

//_____________________________________________________________________________
AliMpVSegmentation* 
AliMUONSegmentationManager::Segmentation(Int_t detElemId, AliMpPlaneType planeType)
{
	TObject* it = fgMap.GetValue(detElemId);
	
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
    fgMap.Add(detElemId,new TPair(b,nb));
    return Segmentation(detElemId,planeType);
  }
}

//_____________________________________________________________________________
TList* 
AliMUONSegmentationManager::SegmentationList(Int_t localBoardNumber)
{
  //
  // Method specific to trigger chamber where a single local trigger board
  // spans several detelemid.
  // This method returns a list of AliMpVSegmentation that contains
  // the given local board.
  //
  // Note that the returned TList is not the owner of its AliMpVSegmentation
  // pointers.
  
  // This method can only work if ALL trigger segmentation have been read in,
  // that's for sure.
  // FIXME: now I'm not so sure the following is the best way to achieve that.
  // Maybe a global way to get the list of detelemid of a stationType would be
  // best, and would avoid to hard-code detelemid range here.
  
  if ( fgLocalBoardMap.GetSize() == 0 )
  {
    for ( Int_t detElemId = 1100; detElemId < 1500; ++detElemId )
    {
      if ( StationType(detElemId) == kStationTrigger )
      {
        AliMpTriggerSegmentation* seg = 
        (AliMpTriggerSegmentation*)Segmentation(detElemId,kNonBendingPlane);
        FillLocalBoardMap(seg);
        seg = (AliMpTriggerSegmentation*)Segmentation(detElemId,kBendingPlane);
        FillLocalBoardMap(seg);
      }
    }
  }
  
  return (TList*)fgLocalBoardMap.GetValue(localBoardNumber);
}

//_____________________________________________________________________________
const char* 
AliMUONSegmentationManager::DetElemName(int detelemid)
{
  if ( ! fgDetElemIdToNameMap.GetSize() ) 
  {
    ReadDetElemIdToName(kStation345);
    ReadDetElemIdToName(kStationTrigger);
    ReadDetElemIdToName(kStation1);
    ReadDetElemIdToName(kStation2);    
  }
  
  TObject* rv = fgDetElemIdToNameMap.GetValue(detelemid);
  
  if ( rv )
  {
    return ((TObjString*)(rv))->String().Data();
  }
  else
  {
    return 0;
  }
}

//_____________________________________________________________________________
AliMpStationType
AliMUONSegmentationManager::StationType(Int_t detelemid)
{
  if (!IsValidDetElemId(detelemid)) return kStationInvalid;
  
  Int_t i = detelemid/100;

  switch (i)
  {
    case 1:
    case 2:  
      return kStation1;
      break;
    case 3:
    case 4:  
      return kStation2;
      break;
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:  
      return kStation345;
      break;
    case 11:
    case 12:
    case 13:
    case 14:  
      return kStationTrigger;
      break;
    default:
      AliErrorClass(Form("%d is not a valid detelemeid\n",detelemid));
      return kStationInvalid;
      break;
  };
}
