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
// $MpId: AliMpMotifReader.cxx,v 1.10 2006/05/24 13:58:41 ivana Exp $
// Category: sector
//
// Class AliMpMotifReader
// -------------------
// Class that takes care of reading the sector data.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpFiles.h"
#include "AliMpMotifReader.h"
#include "AliMpMotifMap.h"
#include "AliMpMotif.h"
#include "AliMpMotifSpecial.h"
#include "AliMpMotifType.h"
#include "AliMpConnection.h"
#include "AliMpIntPair.h"
#include "AliMpDirection.h"

#include "AliLog.h"

#include <TSystem.h>
#include <TMath.h>
#include <Riostream.h>
#include <Rstrstream.h>

#if !defined(__HP_aCC) && !defined(__alpha)
  #include <sstream>
#endif

/// \cond CLASSIMP
ClassImp(AliMpMotifReader)
/// \endcond

//_____________________________________________________________________________
AliMpMotifReader::AliMpMotifReader(AliMpStationType station, 
                                   AliMpPlaneType plane) 
  : TObject(),
    fStationType(station),
    fPlaneType(plane)
{
/// Standard constructor
}

//_____________________________________________________________________________
AliMpMotifReader::AliMpMotifReader() 
  : TObject(),
    fStationType(kStation1),
    fPlaneType(kBendingPlane)
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpMotifReader::~AliMpMotifReader() 
{
/// Destructor  
}

//
// public methods
//

//_____________________________________________________________________________
AliMpMotifType* AliMpMotifReader::BuildMotifType(const TString& motifTypeId)
{
/// Read the files describing a motif in the "$MINSTALL/data" directory
/// and fill the AliMpMotifType structure with.
/// The files mentioned are are named padPos<maskName>.dat
/// and connect<maskName>.dat

  AliMpMotifType*  motifType = new AliMpMotifType(motifTypeId);	

  TString padPosFileName(AliMpFiles::PadPosFilePath(fStationType, 
                                                    fPlaneType, motifTypeId));
  ifstream padPos(padPosFileName);
  AliDebugStream(2) << "Opening file " << padPosFileName << endl;

  PadMapType positions;

  char line[256];
  do {
    padPos.getline(line,255);
    if (!padPos) break;

#if defined (__HP_aCC) || (__alpha)
    strstream strline;
    strline << line;
#else
    istringstream strline(line);
#endif    
    string key;

    strline>>key;
    if ((key=="#") || (key=="") ) continue;

    int i,j;
    strline>>i>>j;
#ifdef WITH_STL
    positions[key].first=i;
    positions[key].second=j;
#endif
#ifdef WITH_ROOT
    positions.Add( AliMpExMap::GetIndex(key), 
                   AliMpExMap::GetIndex(AliMpIntPair(i,j)) ); 
#endif
  } while (!padPos.eof());

  padPos.close();

  TString bergToGCFileName
    = AliMpFiles::BergToGCFilePath(fStationType);
  AliDebugStream(2) << "Opening file " << bergToGCFileName << endl;

  ifstream bergToGCFile(bergToGCFileName);
  const Int_t knbergpins = 
    (fStationType == kStation1 || fStationType == kStation2 ) ? 80 : 100;
  // Station1 & 2 Bergstak connectors have 80 pins, while for stations
  // 3, 4 and 5 they have 100 pins.
  Int_t gassiChannel[100];
  while(1) {
    Int_t bergNum;
    TString gcStr;
    bergToGCFile>>bergNum>>gcStr;
    if (!bergToGCFile.good()) break;
    if (gcStr=="GND") continue;
    if (bergNum>knbergpins) {
        Fatal("BuildMotifType","Berg number > 80 ...");
        continue;
    }
    gassiChannel[bergNum-1]= atoi(gcStr);
  }
  bergToGCFile.close();
  
  TString motifTypeFileName(AliMpFiles::MotifFilePath(fStationType, 
                                                      fPlaneType, motifTypeId));
  ifstream motif(motifTypeFileName);
  AliDebugStream(2) << "Opening file " << motifTypeFileName << endl;

  Int_t nofPadsX=0;
  Int_t nofPadsY=0;

  do {
  
    Int_t ix,iy,numBerg,numKapton,padNum,gassiNum;

    TString lineStr,token;
    lineStr.ReadLine(motif);
    if (!motif.good()) break;
#if defined (__HP_aCC) || (__alpha)
    strstream tokenList;
    tokenList << lineStr.Data();
#else
    istringstream tokenList(lineStr.Data());
#endif    
    
    token.ReadToken(tokenList);
    if (!tokenList.good()) continue; // column is missing...
    if ( (token.Length()>0) && (token[0]=='#') ) continue; // this is a comment line
    
    numBerg = atoi(token.Data());
    if (numBerg==0) {
      AliWarning(Form("Berg number %s invalid",token.Data()));
      continue;
    }
    
    token.ReadToken(tokenList);
    if (!tokenList.good()) continue; // column is missing...
    numKapton = atoi(token.Data());
    if (numKapton==0) continue;

    
    token.ReadToken(tokenList);
    if (!tokenList.good()) continue; // column is missing...
    if (token=="GND") continue;
    string padName = token.Data();
    padNum = motifType->PadNum(token);
    
     token.ReadToken(tokenList);
     if (token.IsNull() ) continue; // column is missing...
//     if (token[0]!='E') {
//       cerr<<"Problem : gassinumber isn't begining with E:"<<token<<endl;
//       continue;
//     }  else {
//        gassiNum = atoi(token.Data() +1 )-1;
//     }
    if ( (numBerg<1) || (numBerg>knbergpins) ) {
      AliWarning(Form("Berg number %d outside range (1..%d)",numBerg,knbergpins));
      continue;
    }
    
    gassiNum  = gassiChannel[numBerg-1];

#ifdef WITH_STL
    PadMapTypeIterator iter = positions.find(padName);
    if (iter==positions.end()) {
      AliWarningStream()
        << "Problem: Pad number " << padNum
	<< " found in the file " << motifTypeFileName
	<< " but not in the file " << padPosFileName << endl;
      continue;
    }

    ix= iter->second.first;
    iy= iter->second.second;
#endif

#ifdef WITH_ROOT
    Long_t value = positions.GetValue(AliMpExMap::GetIndex(padName));
    if (!value) {
      AliWarningStream()
        << "Problem: Pad number " << padNum
	<< " found in the file " << motifTypeFileName
	<< " but not in the file " << padPosFileName << endl;
      continue;
    }

    ix = AliMpExMap::GetPair(value).GetFirst();
    iy = AliMpExMap::GetPair(value).GetSecond();
#endif

    motifType->AddConnection(AliMpIntPair(ix,iy),
                  new AliMpConnection(padNum,numBerg,numKapton,gassiNum));

    if (ix>=nofPadsX) nofPadsX=ix+1;
    if (iy>=nofPadsY) nofPadsY=iy+1;

  } while (!motif.eof());    


  motifType->SetNofPads(nofPadsX, nofPadsY);

  motif.close();

  return motifType;
}


//_____________________________________________________________________________
AliMpMotifSpecial*  
AliMpMotifReader::BuildMotifSpecial(const TString& motifID,
                                    AliMpMotifType* motifType,
                                    Double_t scale)
{
/// Build a special motif by reading the file motifSpecial<motifId>.dat
/// in the data directory

  // Open the input file
  TString motifSpecialFileName(AliMpFiles::MotifSpecialFilePath(fStationType, 
                                                                fPlaneType, motifID));
  ifstream in(motifSpecialFileName);
  if (!in) {	
     AliErrorStream() 
       << "File " << motifSpecialFileName.Data() << " not found." << endl;
     return 0;
  }

  AliMpMotifSpecial* res = new AliMpMotifSpecial(motifID,motifType);
  Int_t i,j;
  Double_t x,y;
  in >> i;
  while (!in.eof()){
    in >>j >>x >> y;
    res->SetPadDimensions(AliMpIntPair(i,j),TVector2(x*scale/2.,y*scale/2.));
    in >> i;
  }
  res->CalculateDimensions();
  
  in.close();
  return res;
}

