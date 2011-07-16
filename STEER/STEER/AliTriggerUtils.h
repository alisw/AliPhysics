#ifndef ALITRIGGERUTILS_H
#define ALITRIGGERUTILS_H

class AliTriggerConfiguration;
class TString;

#include <TObject.h>

class AliTriggerUtils: public TObject {

 public:

  AliTriggerUtils() : TObject() {}
  virtual ~AliTriggerUtils() {}
  
  Bool_t CheckConfiguration( TString & configfile, AliTriggerConfiguration * cfg );


  ClassDef( AliTriggerUtils, 0 )  // Trigger utilities

    };

#endif
