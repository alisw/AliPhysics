#ifndef ALINAMEDSTRING_H
#define ALINAMEDSTRING_H

// $Id$

#include <TObjString.h>

class AliNamedString : public TObjString {
 public: 
  AliNamedString();
  AliNamedString(const char *name, const char *string="");

  void Clear(Option_t* /*option*/="")           { SetString("")      ; }

  const char* GetName()               const { return fName       ; }
  void        SetName(const char* n)        { fName = n          ; }

 protected:
  TString fName; // name of the string object
  
 private:
  AliNamedString(const AliNamedString&);             // not implemented
  AliNamedString& operator=(const AliNamedString&);  // not implemented
  
  ClassDef(AliNamedString, 1); // Named string object
};
#endif
