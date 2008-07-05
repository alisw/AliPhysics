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
/// data files or by string streams filled from AliMpDataMap.
/// The string map is set from outside (AliMpCDB).
/// The data stream must be deleted by client code.
///
/// The selection between files and string streams can
/// be done via function:                                  \n 
/// void   SetReadFromFiles(Bool_t readFromFiles);
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
    AliMpDataStreams();
    virtual ~AliMpDataStreams();
  
    // static methods
    static AliMpDataStreams* Instance();

    // methods
    istream& CreateDataStream(const TString& path) const; 
    Bool_t   IsDataStream(const TString& path) const; 
  
    // set methods
    void   SetDataMap(AliMpDataMap* map);
    void   SetReadFromFiles(Bool_t readFromFiles);
    Bool_t GetReadFromFiles() const;

  private: 
    /// Not implemented
    AliMpDataStreams(const AliMpDataStreams& right);
    /// Not implemented
    AliMpDataStreams& operator=(const AliMpDataStreams& right);    

    // static data members
    static AliMpDataStreams* fgInstance; ///< Singleton instance

    // data members
    // TMap* fMap;
    AliMpDataMap*  fMap;           ///< data map
    Bool_t         fReadFromFiles; ///< option for reading data from files
    
  ClassDef(AliMpDataStreams, 1) //File names and paths 
};  

#endif //ALI_MP_DATA_STREAMS_H
