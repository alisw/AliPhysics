/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpDataMap.h,v 1.4 2006/05/24 13:58:07 ivana Exp $

/// \ingroup basic
/// \class AliMpDataMap
/// \brief TObject class containing a map of strings to strings
///
/// \author Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_DATA_MAP_H
#define ALI_MP_DATA_MAP_H

#include <TObject.h>
#include <TMap.h>

using namespace std;

class AliMpDataMap : public TObject
{
  public:
    AliMpDataMap();
    virtual ~AliMpDataMap();
    
    // methods
    void    Add(const TString& path, const TString& data); 
    TString Get(const TString& path, Bool_t warn = kTRUE) const;
    
    /// Return the map constant reference 
    const TMap& GetMap() const { return fMap; }
    
  private:  
    TMap  fMap; ///< the map of strings to strings
 
  ClassDef(AliMpDataMap,1)  // TObject class containing a map of strings to strings
};

#endif //ALI_MP_EX_MAP_H

