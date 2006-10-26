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

// $Id$
// $MpId: AliMpSt345Reader.cxx,v 1.11 2006/05/24 13:58:50 ivana Exp $

#include "AliMpSt345Reader.h"

#include "AliLog.h"
#include "AliMpMotifReader.h"
#include "AliMpFiles.h"
#include "AliMpMotifType.h"
#include "AliMpPCB.h"
#include "AliMpSlat.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotif.h"
#include "AliMpHelper.h"
#include "AliMpConstants.h"

#include "Riostream.h"
#include "TClass.h"
#include "TObjString.h"
#include "TString.h"

#include <sstream>
#include <assert.h>

/// 
/// \class AliMpSt345Reader
//
/// Read slat and pcb ASCII files.
/// 
/// Basically this class provides 2 static methods :
/// - AliMpSlat* ReadSlat()
/// - AliMpPCB ReadPCB()
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMpSt345Reader)
/// \endcond

TMap AliMpSt345Reader::fgPCBMap;

//_____________________________________________________________________________
AliMpSt345Reader::AliMpSt345Reader() : TObject()
{
  //
  // Default ctor.
  //
} 

//_____________________________________________________________________________
AliMpSt345Reader::~AliMpSt345Reader()
{
  //
  // Dtor.
  //
  fgPCBMap.Delete();
}

//_____________________________________________________________________________
AliMpPCB*
AliMpSt345Reader::PCB(const char* pcbType)
{
  //
  // Get access to an AliMpPCB object, given its type (e.g. N1, SB2, etc...)
  //
  // Note that the returned object is either a new one (read from file) or a 
  // reused one if it is already present in the internal map.
  //
  
  TPair* pair = (TPair*)fgPCBMap.FindObject(pcbType);
  if ( pair )
  {
    AliDebugClass(2,Form("Getting pcb %s from internal map",pcbType));
    return (AliMpPCB*)pair->Value();
  }
  else
  {
    AliDebugClass(2,Form("Reading pcb %s from file",pcbType));
    return ReadPCB(pcbType);
  }
}

//_____________________________________________________________________________
AliMpPCB*
AliMpSt345Reader::ReadPCB(const char* pcbType)
{ 
  //
  // Create a new AliMpPCB object, by reading it from file.
  //
  
  std::ifstream in(AliMpFiles::SlatPCBFilePath(kStation345,pcbType).Data());
  if (!in.good()) 
  {
    AliErrorClass(Form("Cannot open file for PCB %s",pcbType));
    return 0;
  }
 
  AliMpMotifReader reader(kStation345,kNonBendingPlane); 
  // note that the nonbending
  // parameter is of no use for station345, as far as reading motif is 
  // concerned, as all motifs are supposed to be in the same directory
  // (as they are shared by bending/non-bending planes).
     
  char line[80];
  
  const TString kSizeKeyword("SIZES");
  const TString kMotifKeyword("MOTIF");
  
  AliMpPCB* pcb = 0;
  
  while ( in.getline(line,80) )
  {
    if ( line[0] == '#' ) continue;
    
    TString sline(line);
    
    if ( sline(0,kSizeKeyword.Length()) == kSizeKeyword )
    {
      std::istringstream sin(sline(kSizeKeyword.Length(),
                                   sline.Length()-kSizeKeyword.Length()).Data());
      float padSizeX = 0.0;
      float padSizeY = 0.0;
      float pcbSizeX = 0.0;
      float pcbSizeY = 0.0;
      sin >> padSizeX >> padSizeY >> pcbSizeX >> pcbSizeY;
      assert(pcb==0);
      pcb = new AliMpPCB(pcbType,padSizeX,padSizeY,pcbSizeX,pcbSizeY);
    }
    
    if ( sline(0,kMotifKeyword.Length()) == kMotifKeyword )
    {
      std::istringstream sin(sline(kMotifKeyword.Length(),
                                   sline.Length()-kMotifKeyword.Length()).Data());
      TString sMotifType;
      int ix;
      int iy;
      sin >> sMotifType >> ix >> iy;
      
      AliMpMotifType* motifType = 
        reader.BuildMotifType(sMotifType.Data());
      
      assert(pcb!=0);
      pcb->Add(motifType,ix,iy);
    }
  }
  
  in.close();
  
  fgPCBMap.Add(new TObjString(pcbType),pcb);
  return pcb;
}

//_____________________________________________________________________________
AliMpSlat*
AliMpSt345Reader::ReadSlat(const char* slatType, AliMpPlaneType planeType)
{
  //
  // Create a new AliMpSlat object, by reading it from file.
  //
  
  std::ifstream in(AliMpFiles::SlatFilePath(kStation345,slatType,
                                            planeType).Data());
  if (!in.good()) 
  {
    AliErrorClass(Form("Cannot read slat from %s",
                       AliMpFiles::SlatFilePath(kStation345,slatType,planeType).Data()));
    return 0;
  }
  
  char line[80];
  
  const TString kpcbKeyword("PCB");
  
  AliMpSlat* slat = new AliMpSlat(slatType, planeType);
  
  while ( in.getline(line,80) )
  {
    if ( line[0] == '#' ) continue;
    
    TString sline(AliMpHelper::Normalize(line));
    
    if ( sline(0,kpcbKeyword.Length()) == kpcbKeyword )
    {
      TString tmp(sline(kpcbKeyword.Length()+1,sline.Length()-kpcbKeyword.Length()));
      Ssiz_t blankPos = tmp.First(' ');
      if ( blankPos < 0 )
	    {
        AliErrorClass("Syntax error in PCB file, should get a list of "
                      "manu ids after the pcbname");
        delete slat;
        return 0;
	    }
      
      TString pcbName(tmp(0,blankPos));
      TString manus(tmp(blankPos+1,tmp.Length()-blankPos));
      
      AliMpPCB* pcbType = PCB(pcbName.Data());	  
      if (!pcbType)
	    {
        AliErrorClass(Form("Cannot read pcbType=%s",pcbName.Data()));
	      delete slat;
	      return 0;
	    }      

      TArrayI manuList;
      AliMpHelper::DecodeName(manus,';',manuList);
      if ( manuList.GetSize() != Int_t(pcbType->GetSize()) )
	    {
        AliErrorClass(Form("Wrong number of manu ids for this PCB ("
                           "%s) : %d out of %d",pcbName.Data(),
                           manuList.GetSize(),pcbType->GetSize()));
	      delete slat;
	      return 0;
      }

      for ( Int_t i = 0; i < manuList.GetSize(); ++i )
      {
        manuList[i] |= AliMpConstants::ManuMask(planeType);
      }
      slat->Add(pcbType,manuList);
    }
  }
  
  in.close();
  
  return slat;
}

