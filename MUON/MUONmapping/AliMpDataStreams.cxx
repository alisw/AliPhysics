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
// $MpId: AliMpDataStreams.cxx,v 1.12 2006/05/23 13:09:54 ivana Exp $
// Category: basic

//-----------------------------------------------------------------------------
// Class AliMpDataStreams
// ----------------------
// Class for providing mapping data streams
// See detailed description in the header file.
// Author: Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpDataStreams.h"
#include "AliMpDataMap.h"
#include "AliMpFiles.h"

#include "AliLog.h"

#include <TMap.h>
#include <TFile.h>
#include <TObjString.h>
#include <TString.h>
#include <Riostream.h>

#include <string>

/// \cond CLASSIMP
ClassImp(AliMpDataStreams)
/// \endcond


//______________________________________________________________________________
AliMpDataStreams::AliMpDataStreams(AliMpDataMap* map) 
  : TObject(),
    fMap(map),
    fReadFromFiles(kTRUE)
{
/// Standard and default constructor

  if ( map ) fReadFromFiles = kFALSE;
}

//______________________________________________________________________________
AliMpDataStreams::AliMpDataStreams(TRootIOCtor* /*ioCtor*/) 
  : TObject(),
    fMap(0),
    fReadFromFiles()
{
/// Root IO constructor

}

//______________________________________________________________________________
AliMpDataStreams::~AliMpDataStreams() 
{
/// Destructor

  // delete fMap;
       // Do not delete data map as it is a CDB object 
       // which is cached
}

//
// private methods
//

//______________________________________________________________________________
void AliMpDataStreams::CutDataPath(string& dataPath) const
{ 
/// Cut the path defined in AliMpFiles as Top() + one more directory

  string top = AliMpFiles::GetTop().Data();
  if ( dataPath.find(top) != string::npos ) dataPath.erase(0, top.size()+1);
  dataPath.erase(0,dataPath.find('/')+1); 
}
  


//
// public methods
//

//______________________________________________________________________________
istream& AliMpDataStreams::CreateDataStream(const TString& path) const
{
/// Return the string with data in the mapping file spcified with path.
/// Both full path in the file system and a short path (without 
/// $LICE_ROOT/mapping/data string) can be used.


  if ( fReadFromFiles ) {                                                                
    AliDebugStream(2) << "Opening file " << path.Data() << endl;
    ifstream* fileBuffer = new ifstream();
    fileBuffer->open(path.Data());
    if ( ! fileBuffer->good() ) {	
       AliErrorStream() 
         << "Cannot open file " << path.Data() << endl;
    }     
    return *fileBuffer;
  }
  else {
    AliDebugStream(2) << "Opening stream " << path.Data() << endl;

    // Cut top from the path 
    string dataPath = path.Data();
    CutDataPath(dataPath);

    istringstream* stringBuffer 
      = new istringstream(fMap->Get(dataPath).Data());
    return *stringBuffer;
   }    
}

//______________________________________________________________________________
Bool_t  AliMpDataStreams::IsDataStream(const TString& path) const
{
/// Return true, if data with given path exists

  if ( fReadFromFiles ) {                                                                
    ifstream fileBuffer(path.Data());
    return fileBuffer.good();
  }
  else {
    // Cut top from the path 
    string dataPath = path.Data();
    CutDataPath(dataPath);

    return ( fMap->Get(dataPath, kFALSE) != "" );
  }
}      

//______________________________________________________________________________
void  AliMpDataStreams::SetReadFromFiles()
{
/// Set option to read data from files

  fReadFromFiles = kTRUE;
}  

//______________________________________________________________________________
Bool_t AliMpDataStreams::GetReadFromFiles() const
{ 
/// Return the info where the data are loaded from

  return fReadFromFiles; 
}  
