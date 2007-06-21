#ifndef ALIBCMLOADER_H
#define ALIBCMLOADER_H

/////////////////////////////////////////////////////////////////////
//                                                                 //
// Base class for BCMloaders.                                      //
// Loader provides base I/O facilities for standard data.          //
// Each detector has a loader data member.                         //
// Loader is always accessible via folder structure as well.       // 
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliLoader.h"

class AliBCMLoader: public AliLoader
 {
   public:
    AliBCMLoader();
    AliBCMLoader(const Char_t *name,const Char_t *topfoldername);
    AliBCMLoader(const Char_t *name,TFolder *topfolder);    
    virtual ~AliBCMLoader() {};
    
    AliBCMLoader & operator = (const AliBCMLoader & ) {return *this;}
    
   private:
    static const TString fgkDefaultHitsFileName;  // Default Name for hit file
    static const TString fgkDefaultDigitsFileName;// Default Name for digit file

   ClassDef(AliBCMLoader,1)
      
 };
 
#endif
