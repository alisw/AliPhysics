// $Id$
// Category: motif
//
// Class AliMpMotifMap
// -------------------
// Class describing the motif map container, where motifs are
// mapped to their string IDs.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_MOTIF_MAP_H
#define ALI_MP_MOTIF_MAP_H

#include <map>

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
    AliMpMotifPosition*  FindMotifPosition(const AliMpIntPair& indices) const;

  private:
    // methods
    void  PrintMotifs() const;
    void  PrintMotifTypes() const;
    void  PrintMotifPositions() const;
    void  PrintMotifPositions2() const;
 
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
    MotifMap          fMotifs;     // motifs map
    MotifTypeMap      fMotifTypes; // motif types map 
    //MotifPositionMap  fMotifPositions;  // motif positions map 
                                          // not taken by cint
    std::map<Int_t, AliMpMotifPosition*> fMotifPositions; // motif positions map by Id
    MotifPositionMap2 fMotifPositions2; // motif positions map
#endif

  ClassDef(AliMpMotifMap,1)  // motif map
};

#endif //ALI_MP_MOTIF_MAP_H
