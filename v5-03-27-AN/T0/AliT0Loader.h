#ifndef ALIT0LOADER_H
#define ALIT0LOADER_H

#include "AliLoader.h"
#include "AliObjectLoader.h"

class AliT0digit;


class AliT0Loader: public AliLoader {
 public:
   AliT0Loader() : AliLoader() {};
   AliT0Loader(const Char_t *detname,const Char_t *eventfoldername) : 
     AliLoader(detname, eventfoldername) {InitObjectLoaders();};
   AliT0Loader(const Char_t *detname,TFolder* eventfolder) :
     AliLoader(detname, eventfolder) {InitObjectLoaders();};

   // Digits
   AliT0digit*  Digits(){ return (AliT0digit*) GetDigitsDataLoader()->GetBaseDataLoader()->Get();} // returns a pointer to the tree of  RawClusters

 private:
   void InitObjectLoaders();

   ClassDef(AliT0Loader,1)
};

typedef AliT0Loader AliSTARTLoader; // for backward compatibility
 
#endif



