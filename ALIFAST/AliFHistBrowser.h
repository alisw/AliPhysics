#ifndef AliFHistBrowser_H
#define AliFHistBrowser_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFHistBrowser                                                      //
//                                                                      //
// helper class to browse AliFast Makers histograms.                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

class AliFHistBrowser : public TNamed {

public:
                     AliFHistBrowser();
   virtual          ~AliFHistBrowser() {;}
   virtual void      Browse(TBrowser *b);
   Bool_t            IsFolder() {return kTRUE;}

   ClassDef(AliFHistBrowser, 0)   //helper class to browse AliFast Makers histograms
};

#endif
