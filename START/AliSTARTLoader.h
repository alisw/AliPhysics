#ifndef ALISTARTLOADER_H
#define ALISTARTLOADER_H

#include "AliLoader.h"

class AliSTARTdigit;


class AliSTARTLoader: public AliLoader {
 public:
   AliSTARTLoader() : AliLoader() {};
   AliSTARTLoader(const Char_t *detname,const Char_t *eventfoldername) : 
     AliLoader(detname, eventfoldername) {InitObjectLoaders();};
   AliSTARTLoader(const Char_t *detname,TFolder* eventfolder) :
     AliLoader(detname, eventfolder) {InitObjectLoaders();};

   // Digits
   AliSTARTdigit*  Digits(){ return (AliSTARTdigit*) GetDigitsDataLoader()->GetBaseDataLoader()->Get();} // returns a pointer to the tree of  RawClusters

 private:
   void InitObjectLoaders();

   ClassDef(AliSTARTLoader,1)
};
 
#endif



