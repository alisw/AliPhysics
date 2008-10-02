// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

#ifndef ALIDANAL_H
#define ALIDANAL_H

///////////////////////////////////////////////////////////////////////////////
// Parent analysis class
///////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TObjArray.h>
#include <TDatabasePDG.h>

/*****************************************************************************/
class AliDAnal : public TNamed {
   
 public:
  AliDAnal(char *nam);                        // constructor
  virtual ~AliDAnal()                         {printf("AliDAnal object %s destroyed\n",GetName());}
  void Save(const char *outfil, const char *mode="update");  // save histograms on root file

 protected:
  static TDatabasePDG *fPDG;               // particle database
  TObjArray fHistos;                       // histograms

  ClassDef(AliDAnal,1)
};
/*****************************************************************************/
#endif
