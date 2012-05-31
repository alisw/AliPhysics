/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpDataStreams.h,v 1.10 2006/05/24 13:58:07 ivana Exp $

/// \ingroup basic
/// \class AliMpDataStreams
/// \brief Mapping data streams provider
///
/// The class provides input streams for mapping data;
/// the input streams can be represented either by the
/// data files or by string streams filled from string map (AliMpDataMap).
/// The string map is set from outside (AliMpCDB) and is not
/// deleted in this class.
/// The data streams returned by CreateDataStream() function
/// must be deleted by the client code.
///
/// The selection between files and string streams is 
/// done in the constructor:
/// if data map is provided, reading is performed from streams,
/// otherwise reading is performed from file.
/// User can also use the set function to select reading
/// from files also when the data map is provided: \n
/// void SetReadFromFiles();
///
/// \author Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_DATA_STREAMS_H
#define ALI_MP_DATA_STREAMS_H

#include "AliMpDataMap.h"

#include <TObject.h>
#include <TString.h>

#include <fstream>
#include <sstream>

//class TMap;
class AliMpDataMap;

class AliMpDataStreams : public TObject
{
  public:
    AliMpDataStreams(AliMpDataMap* map = 0x0);
    AliMpDataStreams(TRootIOCtor* ioCtor);
    virtual ~AliMpDataStreams();
  
    // methods
    istream& CreateDataStream(const TString& path) const; 
    Bool_t   IsDataStream(const TString& path) const; 
  
    // set methods
    void   SetReadFromFiles();
    Bool_t GetReadFromFiles() const;

  private: 
    /// Not implemented
    AliMpDataStreams(const AliMpDataStreams& right);
    /// Not implemented
    AliMpDataStreams& operator=(const AliMpDataStreams& right);    

    // methods
    void CutDataPath(string& dataPath) const;

    // data members
    AliMpDataMap*  fMap;           ///< data map
    Bool_t         fReadFromFiles; ///< option for reading data from files
    
  ClassDef(AliMpDataStreams, 1) //File names and paths 
};  

#endif //ALI_MP_DATA_STREAMS_H
