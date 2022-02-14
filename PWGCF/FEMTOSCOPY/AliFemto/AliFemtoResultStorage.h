///
/// \file AliFemtoResultStorage.h
///


#pragma once


#ifndef ALIFEMTORESULTSTORAGE_H
#define ALIFEMTORESULTSTORAGE_H

#include "AliLog.h"

#include <TList.h>
#include <TNamed.h>


class AliFemtoManager;


/// \class AliFemtoResultStorage
/// \brief Objects used as output-container
///
/// Storage is in the form of a TDirectory created upon calling `Write()`
///
class AliFemtoResultStorage : public TNamed {
public:
  /// Construct with defaults
  AliFemtoResultStorage();

  /// Construct with name
  AliFemtoResultStorage(const TString &name);

  /// Construct with name and List
  AliFemtoResultStorage(const TString &name, TList *);

  /// Construct with name and fill with objects in femto manager
  AliFemtoResultStorage(const TString &name, AliFemtoManager &);

  /// destroy the objects
  virtual ~AliFemtoResultStorage();

  /// Add an output list to storage
  void AddOutputList(TList *list)
    { fObjects->AddLast(list); }

  /// Called upon by AliAnalysisManager to save data to container
  virtual Int_t	Write(const char *name=NULL, Int_t option=0, Int_t bufsize=0);
  virtual Int_t	Write(const char *name=NULL, Int_t option=0, Int_t bufsize=0) const;

private:
  /// non-copiable
  AliFemtoResultStorage(AliFemtoResultStorage const &);
  AliFemtoResultStorage& operator=(AliFemtoResultStorage const &);

protected:

  /// objects
  TObjArray *fObjects;

  /// \cond CLASSIMP
  ClassDef(AliFemtoResultStorage, 1);
  /// \endcond
};


inline
AliFemtoResultStorage::AliFemtoResultStorage(AliFemtoResultStorage const &orig)
  : TNamed(orig)
  , fObjects(nullptr)
{
  AliError("Uncopiable");
}


inline
AliFemtoResultStorage& AliFemtoResultStorage::operator=(AliFemtoResultStorage const &rhs)
{
  AliError("Unassignable");
  return *this;
}

#endif
