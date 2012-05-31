#ifndef ALIACORDELOADER_H
#define ALIACORDELOADER_H

/////////////////////////////////////////////////////////////////////
//                                                                 //
// Base class for ACORDEloaders.                                    //                                          
// Loader provides base I/O facilities for standard data.          //
// Each detector has a loader data member.                         //
// Loader is always accessible via folder structure as well.       // 
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliLoader.h"

class AliACORDELoader: public AliLoader
 {
   public:
    AliACORDELoader();
    AliACORDELoader(const Char_t *name,const Char_t *topfoldername);
    AliACORDELoader(const Char_t *name,TFolder *topfolder);    
    virtual ~AliACORDELoader() {};
    
    AliACORDELoader & operator = (const AliACORDELoader & ) {return *this;}
    
   private:
    static const TString fgkDefaultHitsFileName;  // Default Name for hit file
    static const TString fgkDefaultDigitsFileName;// Default Name for digit file

   ClassDef(AliACORDELoader,1)
      
 };
 
#endif
