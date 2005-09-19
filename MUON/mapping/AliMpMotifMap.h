/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpMotifMap.h,v 1.7 2005/08/26 15:43:36 ivana Exp $

/// \ingroup motif
/// \class AliMpMotifMap
/// \brief Motif map containers

/// The class defines:
/// - map of motif objects to their string IDs
/// - map of motif type objects to their string IDs
/// - map of motif position objects to their string IDs
/// - map of motif position objects to their global indices
///
/// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_MAP_H
#define ALI_MP_MOTIF_MAP_H

#ifdef WITH_STL
  #include <map>
#endif

#include <TObject.h>

#include "AliMpMotifTypes.h"

class TString;
class TVector2;

class AliMpVMotif;
class AliMpMotifType;
class AliMpMotifPosition;
class AliMpMotifMap;

class AliMpMotifMap : public TObject
{
  public:
    AliMpMotifMap();
    virtual ~AliMpMotifMap();
    
    // methods
    Bool_t  AddMotif(AliMpVMotif* motif, Bool_t warn = true);
    Bool_t  AddMotifType(AliMpMotifType* motifType, Bool_t warn = true);
    Bool_t  AddMotifPosition(AliMpMotifPosition* motifType, Bool_t warn = true);
    void   FillMotifPositionMap2();
    virtual void Print(const char* /*option*/ = "") const;
    void   PrintGlobalIndices(const char* fileName) const;
    void   UpdateGlobalIndices(const char* fileName);
   
    // find methods
    AliMpVMotif*  FindMotif(const TString& motifID) const;
    AliMpVMotif*  FindMotif(const TString& motifID, const TString& motifTypeID, 
                            const TVector2& padDimensions) const;
    AliMpMotifType*      FindMotifType(const TString& motifTypeID) const;
    AliMpMotifPosition*  FindMotifPosition(Int_t motifPositionID) const;
    // AliMpMotifPosition*  FindMotifPosition(const AliMpIntPair& indices) const;

  private:
#ifdef WITH_ROOT
    static const Int_t   fgkSeparator;  // the separator used for conversion
                                        // of TString to Int_t
    
    // methods
    Int_t  GetIndex(const TString& s) const;
    Int_t  GetIndex(const AliMpIntPair& pair) const;
    TString  GetString(Int_t index) const;
    AliMpIntPair  GetPair(Int_t index) const;
#endif
  
    // methods
    void  PrintMotif(const AliMpVMotif* motif) const;
    void  PrintMotifType(const AliMpMotifType* motifType) const;
    void  PrintMotifPosition(const AliMpMotifPosition* motifPosition) const;
    void  PrintMotifPosition2(const AliMpMotifPosition* motifPosition) const;
    void  PrintMotifs() const;
    void  PrintMotifTypes() const;
    void  PrintMotifPositions() const;
    void  PrintMotifPositions2() const;
 
#ifdef WITH_STL
#ifdef __HP_aCC
    // data members
            // EXCLUDED FOR CINT (does not compile on HP)
    MotifMap          fMotifs; //! motifs map
    MotifTypeMap      fMotifTypes; //!motif types map 
    //MotifPositionMap  fMotifPositions;  //! motif positions map 
                                          // not taken by cint
    map<Int_t, AliMpMotifPosition*>  fMotifPositions; //! motif positions map by Id
    MotifPositionMap2 fMotifPositions2; //! motif positions map
#else
    // data members
    MotifMap          fMotifs; // motifs map
    MotifTypeMap      fMotifTypes; // motif types map 
    //MotifPositionMap  fMotifPositions;  // motif positions map 
                                          // not taken by cint
    std::map<Int_t, AliMpMotifPosition*> fMotifPositions; // motif positions map by Id
    MotifPositionMap2 fMotifPositions2; // motif positions map
#endif
#endif

#ifdef WITH_ROOT
    // data members
    mutable MotifMap           fMotifs;     //  motifs map
    mutable MotifTypeMap       fMotifTypes; //  motifs types map
    mutable MotifPositionMap   fMotifPositions; //  motifs positions map
    mutable MotifPositionMap2  fMotifPositions2;//  motifs positions map
#endif

  ClassDef(AliMpMotifMap,1)  // motif map
};

#endif //ALI_MP_MOTIF_MAP_H
