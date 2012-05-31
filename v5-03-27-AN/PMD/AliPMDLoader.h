#ifndef ALIPMDGETTER_H
#define ALIPMDGETTER_H

//base class for loaders 
//loader is common for reading data for all detectors
//Each detector has a loader data member
//loader is accessible via folder structure as well

#include <AliLoader.h>



class AliPMDLoader: public AliLoader
 {
   public:
    AliPMDLoader();
    AliPMDLoader(const Char_t *name,const Char_t *topfoldername);
    AliPMDLoader(const Char_t *name,TFolder *topfolder);
    
    virtual ~AliPMDLoader(){};//-----------------

   protected:


   private:
    static const TString fgkDefaultHitsFileName;
    static const TString fgkDefaultSDigitsFileName;
    static const TString fgkDefaultDigitsFileName;
    static const TString fgkDefaultRecPointsFileName;
    static const TString fgkDefaultTracksFileName;
    

   public:
     ClassDef(AliPMDLoader,1)
 };
 
#endif
