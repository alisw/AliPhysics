/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpExMap.h,v 1.4 2006/05/24 13:58:07 ivana Exp $

/// \ingroup core
/// \class AliMpExMap
/// \brief Helper class making Root persistent TExMap
///
/// The objects and keys from TExMap are store in additional
/// arrays which are Root persistent.
///
/// \author Ivana Hrivnacova; IPN Orsay

#ifndef ALI_MP_EX_MAP_H
#define ALI_MP_EX_MAP_H

#include <TObject.h>
#include <TObjArray.h>
#include <TArrayL.h>
#include <TExMap.h>

class AliMpExMapIterator;

class TString;

class AliMpExMap : public TObject
{
  friend class AliMpExMapIterator;
  
  public:
    AliMpExMap();
    AliMpExMap(TRootIOCtor* /*ioCtor*/);
    AliMpExMap(const AliMpExMap& rhs);
    AliMpExMap& operator=(const AliMpExMap& rhs);
    virtual ~AliMpExMap();
    
    
    // static methods
    // conversion between String and Long_t
    static Long_t  GetIndex(const TString& s);
    static TString GetString(Long_t index);

    // methods from base class
    virtual void Clear(Option_t* opt="");
    virtual void Print(Option_t* opt="") const;

    // set methods
    void Add(Int_t keyFirst, Int_t keySecond, TObject* object);
    void Add(const TString& key, TObject* object);
    void Add(Int_t key, TObject* object);

    void SetSize(Int_t size);
    void SetOwner(Bool_t owner);
    
    // get methods
    Int_t GetSize() const;
    Int_t GetCapacity() const;
    
    TObject*    GetValue(Int_t keyFirst, Int_t keySecond) const;
    TObject*    GetValue(const TString& key) const;
    TObject*    GetValue(Int_t key) const;

    AliMpExMapIterator* CreateIterator() const;
    
  private:  
    // static methods
    static const TString&  GetCharacterMap();

    // methods

    // conversion between Int_t pair key and Long_t
    static Long_t  GetIndex(Int_t first, Int_t second);
    static Int_t   GetPairFirst(Long_t index);
    static Int_t   GetPairSecond(Long_t index);

    void FillMap();
    void AddKey(Long_t key);
    void Copy(TObject& dest) const;
    
    // static data members
    static const Int_t    fgkDefaultSize;      ///< Default initial size
    static const Bool_t   fgkDefaultOwnership; ///< Default ownership

    static const Int_t    fgkSeparator1; ///< \brief the separator used for conversion
                                         ///  of Int_t pair to Int_t
    static const Int_t    fgkSeparator2; ///< \brief the separator used for conversion
                                         ///  of TString to Int_t
    
    // data members
    mutable TExMap  fMap;     //!<  Transient map class
    TObjArray       fObjects; ///<  Array of objects 
    TArrayL         fKeys;    ///<  Array of keys 

  ClassDef(AliMpExMap,1)  // Root persistent TExMap
};

#endif //ALI_MP_EX_MAP_H

