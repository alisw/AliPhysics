

#pragma once

#ifndef AliForwardFlowResultStorage_H
#define AliForwardFlowResultStorage_H

#include "AliLog.h"

#include <TList.h>
#include <TNamed.h>



/// \class AliForwardFlowResultStorage
/// \brief Objects used as output-container
///
/// Storage is in the form of a TDirectory created upon calling `Write()`
///
class AliForwardFlowResultStorage : public TNamed {
public:
  /// Construct with defaults
  AliForwardFlowResultStorage();

  /// Construct with name
  AliForwardFlowResultStorage(const TString &name);

  /// Construct with name and List
  AliForwardFlowResultStorage(const TString &name, TList *);

  /// destroy the objects
  virtual ~AliForwardFlowResultStorage();

  /// Add an output list to storage
  void AddOutputList(TList *list)
    { fObjects->AddLast(list); }

  /// Called upon by AliAnalysisManager to save data to container
  virtual Int_t	Write(const char *name=NULL, Int_t option=0, Int_t bufsize=0);
  virtual Int_t	Write(const char *name=NULL, Int_t option=0, Int_t bufsize=0) const;

private:
  /// non-copiable
  AliForwardFlowResultStorage(AliForwardFlowResultStorage const &);
  AliForwardFlowResultStorage& operator=(AliForwardFlowResultStorage const &);

protected:

  /// objects
  TObjArray *fObjects;

  /// \cond CLASSIMP
  ClassDef(AliForwardFlowResultStorage, 1);
  /// \endcond
};


inline
AliForwardFlowResultStorage::AliForwardFlowResultStorage(AliForwardFlowResultStorage const &orig)
  : TNamed(orig)
  , fObjects(nullptr)
{
  AliError("Uncopiable");
}


inline
AliForwardFlowResultStorage& AliForwardFlowResultStorage::operator=(AliForwardFlowResultStorage const &rhs)
{
  AliError("Unassignable");
  return *this;
}

#endif