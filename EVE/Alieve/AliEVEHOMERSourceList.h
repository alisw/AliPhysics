// $Header$

#ifndef ALIEVE_AliEVEHOMERSourceList_H
#define ALIEVE_AliEVEHOMERSourceList_H

#include <Reve/Reve.h>
#include <Reve/RenderElement.h>

#include <TObject.h>

class AliEVEHOMERSourceList : public Reve::RenderElementList
{
private:
  AliEVEHOMERSourceList(const AliEVEHOMERSourceList&);            // Not implemented
  AliEVEHOMERSourceList& operator=(const AliEVEHOMERSourceList&); // Not implemented

protected:

public:
  AliEVEHOMERSourceList(const Text_t* n="HOMER Source List", const Text_t* t="");
  virtual ~AliEVEHOMERSourceList() {}

  void SelectAll();   // *MENU*
  void DeselectAll(); // *MENU*

  ClassDef(AliEVEHOMERSourceList, 1);
}; // endclass AliEVEHOMERSourceList

#endif
