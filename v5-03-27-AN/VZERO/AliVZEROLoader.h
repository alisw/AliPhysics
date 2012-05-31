#ifndef ALIVZEROLOADER_H
#define ALIVZEROLOADER_H

/////////////////////////////////////////////////////////////////////
//                                                                 //
// Base class for VZEROloaders.                                    //                                          
// Loader provides base I/O facilities for standard data.          //
// Each detector has a loader data member.                         //
// Loader is always accessible via folder structure as well.       // 
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliLoader.h"

class AliVZEROLoader: public AliLoader
 {
   public:
    AliVZEROLoader();
    AliVZEROLoader(const Char_t *name,const Char_t *topfoldername);
    AliVZEROLoader(const Char_t *name,TFolder *topfolder);    
    virtual ~AliVZEROLoader() {};
    
    AliVZEROLoader & operator = (const AliVZEROLoader & ) {return *this;}
    
   private:
    static const TString fgkDefaultHitsFileName;  // Default Name for hit file
    static const TString fgkDefaultDigitsFileName;// Default Name for digit file

   ClassDef(AliVZEROLoader,1)
      
 };
 
#endif
