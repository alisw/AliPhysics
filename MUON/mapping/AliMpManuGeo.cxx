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
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$
// $MpId: AliMpManuGeo.cxx,v 1.5 2006/05/24 13:58:34 ivana Exp $
// Category: management

// Class AliMpManuGeo
// ---------------
// Class that manages the maps manuId<>manuSerial#<>DE 
// Needed for geometrical calibration gain for manu
// Author: Ch. Finck; Subatech Nantes

#include <TString.h>

#include "AliMpManuGeo.h"
#include "AliMpFiles.h"
#include "AliMpHelper.h"
#include "AliMpExMap.h"
#include "AliMpStationType.h"
#include "AliMpDEManager.h"
#include "AliMpIntPair.h"
#include "AliMpDEIterator.h"

#include "AliLog.h"

#include "Riostream.h"

/// \cond CLASSIMP
ClassImp(AliMpManuGeo)
/// \endcond


//_____________________________________________________________________________
AliMpManuGeo::AliMpManuGeo()
  : TObject(),
    fDeManuToSerialNb(17000),
    fSerialNbToDeManu(17000)
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpManuGeo::~AliMpManuGeo() 
{
/// Destructor

  fDeManuToSerialNb.Delete();
  fSerialNbToDeManu.Delete();

}

//____________________________________________________________________
AliMpIntPair  AliMpManuGeo::GetDetElemManu(Int_t manuSerial) 
{
 /// getting (DE, manuId) from manu serial number
  return * (AliMpIntPair*) fSerialNbToDeManu.GetValue(manuSerial);

}

//____________________________________________________________________
Int_t AliMpManuGeo::GetManuSerial(AliMpIntPair& pair) 
{
/// getting manu serial number from (DE, manuId)

  if ( ! AliMpDEManager::IsValidDetElemId(pair.GetFirst(), true) )  return -1;

  Long_t it = fDeManuToSerialNb.GetValue(AliMpExMap::GetIndex(pair));

 if ( it ) 
   return (Int_t)it;
 else 
   return -1;
}

//____________________________________________________________________
void AliMpManuGeo::ReadGeomManuFiles()
{
/// Read manu serial numbers for all detection elements

  AliMpDEIterator it;
  for ( it.First(); ! it.IsDone(); it.Next() ) {
    ReadGeomManuFile(it.CurrentDE());
  }
}


//____________________________________________________________________
void AliMpManuGeo::ReadGeomManuFile(Int_t idDE)
{
/// Read manu serial numbers for the given detection elements
  
  AliMpStationType stationType = AliMpDEManager::GetStationType(idDE);
  TString detFileName = AliMpDEManager::GetDEName(idDE);
  TString infile = AliMpFiles::ManuToSerialPath(detFileName, stationType);

  ifstream in(infile, ios::in);
  if (!in) AliError(Form("File %s not found.", infile.Data()));
       
  char line[80];

  while ( in.getline(line,80) ) {

    if ( line[0] == '#' ) continue;

    TString tmp(AliMpHelper::Normalize(line));

    Int_t blankPos  = tmp.First(' ');

    TString sManuId(tmp(0, blankPos));

    Int_t manuId = atoi(sManuId.Data());

    TString sManuSerial(tmp(blankPos + 1, tmp.Length()-blankPos));

    Int_t manuSerial = atoi(sManuSerial.Data());
      
   
    // filling (idDE, manuId) <> manuSerial
    fDeManuToSerialNb.Add(AliMpExMap::GetIndex(AliMpIntPair(idDE, manuId)), (Long_t)manuSerial); 
    fSerialNbToDeManu.Add((Long_t)manuSerial, (Long_t)new AliMpIntPair(idDE, manuId)); 

  }
   
  in.close();

}
