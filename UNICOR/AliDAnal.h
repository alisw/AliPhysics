// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

//=============================================================================
// parent class of all analyzers
//=============================================================================

#ifndef ALIDANAL_H
#define ALIDANAL_H

#include <TNamed.h>
#include <TObjArray.h>
#include <TDatabasePDG.h>

//=============================================================================
class AliDAnal : public TNamed {
   
 public:
  AliDAnal(char *nam);                                         // constructor
  virtual ~AliDAnal()     {printf("%s object named %s deleted\n",ClassName(),GetName());}
  void Save(const char *outfil, const char *mode="update");  // save histograms 

 protected:
  static TDatabasePDG fgPDG; // particle database
  TObjArray fHistos;          // histograms

  ClassDef(AliDAnal,1)
};
//=============================================================================
#endif
