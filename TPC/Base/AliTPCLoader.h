#ifndef ALITPCGETTER_H
#define ALITPCGETTER_H

//base class for loaders 
//loader is common for reading data for all detectors
//Each detector has a loader data member
//loader is accessible via folder structure as well

#include <AliLoader.h>



class AliTPCLoader: public AliLoader
 {
   public:
    AliTPCLoader();
    AliTPCLoader(const Char_t *name,const Char_t *topfoldername);
    AliTPCLoader(const Char_t *name,TFolder *topfolder);
    
    virtual ~AliTPCLoader(){};//-----------------

   protected:


   private:
    static const TString fgkDefaultHitsFileName;
    static const TString fgkDefaultSDigitsFileName;
    static const TString fgkDefaultDigitsFileName;
    static const TString fgkDefaultRecPointsFileName;
    static const TString fgkDefaultTracksFileName;
    

   public:
     ClassDef(AliTPCLoader,1)
 };
 
#endif
