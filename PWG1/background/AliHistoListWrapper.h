#ifndef ALIHISTOLISTWRAPPER_H
#define ALIHISTOLISTWRAPPER_H
 
//-------------------------------------------------------------------------
//                      AliHistoListWrapper
// This class is used to contain a list of histograms to be merged
// with another list, not necessarily containing the same histos in the
// same order.
//
// The merging method checks if the lists contain the same histograms
// and, if not, adds an empty copy of the missing histograms to the
// relevant list.
//
// Can only contain objects inheriting from TH1. Can be useful if you
// want to run a task on CAF and book your histograms dinamically
// (e.g. in your UserExec method rather than in the UserCreateObject)
// 
// Author: Michele Floris, CERN
//
//-------------------------------------------------------------------------

#include <TNamed.h>
#include "TList.h"

class TList;
class TH1;
class TCollection;




class AliHistoListWrapper : public TNamed
{
public:

  AliHistoListWrapper();
  AliHistoListWrapper(const char* name, const char* title);
  AliHistoListWrapper(const AliHistoListWrapper& obj);  
  ~AliHistoListWrapper();

  void AddHistoToList(TObject* h){fList->Add(h);}
  TList * GetList(){return fList;}

  Long64_t Merge(TCollection* list);

  AliHistoListWrapper& operator=(const AliHistoListWrapper& wrap);

protected:
  TList * fList;

private:

  ClassDef(AliHistoListWrapper, 1); 
};
 
#endif
