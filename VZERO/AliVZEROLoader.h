#ifndef ALIVZEROLOADER_H
#define ALIVZEROLOADER_H

//base class for loaders 
//loader is common for reading data for all detectors
//Each detector has a loader data member
//loader is accessible via folder structure as well

#include <AliLoader.h>

class AliVZEROLoader: public AliLoader
 {
   public:
    AliVZEROLoader();
    AliVZEROLoader(const Char_t *name,const Char_t *topfoldername);
    AliVZEROLoader(const Char_t *name,TFolder *topfolder);
    
    virtual ~AliVZEROLoader(){};//-----------------

   protected:


   private:
    static const TString fgkDefaultHitsFileName;
    static const TString fgkDefaultDigitsFileName;

   public:
     ClassDef(AliVZEROLoader,1)
 };
 
#endif
