///
/// \file AliFemtoResultStorage.h
///


#pragma once


#ifndef ALIFEMTORESULTSTORAGE_H
#define ALIFEMTORESULTSTORAGE_H

#include <TNamed.h>
#include "AliLog.h"


/// \class AliFemtoResultStorage
/// \brief Objects used as output-container
///
class AliFemtoResultStorage : public TNamed {
public:
  /// Construct with defaults
  AliFemtoResultStorage();

  /// Construct with name
  AliFemtoResultStorage(const TString &name);

  /// Construct
  AliFemtoResultStorage(const TString &name, TList *);

private:
  /// non-copiable
  AliFemtoResultStorage(AliFemtoResultStorage const &);
  AliFemtoResultStorage& operator=(AliFemtoResultStorage const &);

public:
  //*
  /// Called upon by AliAnalysisManager to save data to container
  virtual Int_t	Write(const char *name=0, Int_t option=0, Int_t bufsize=0);
  virtual Int_t	Write(const char *name=0, Int_t option=0, Int_t bufsize=0) const;
  //*/

// protected:

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
