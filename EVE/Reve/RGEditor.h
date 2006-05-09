// $Header$

#ifndef REVE_RGEditor_H
#define REVE_RGEditor_H

#include <TGedEditor.h>

namespace Reve {

class RGEditor : public TGedEditor
{
protected:

public:
  RGEditor(TCanvas* canvas=0);
  virtual ~RGEditor() {}
  
  void DisplayObject(TObject* obj);
  TGCanvas * GetCan() const {return fCan;}

  ClassDef(RGEditor, 1);
}; // endclass RGEditor

}

#endif
