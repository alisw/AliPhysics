#ifndef AliFBigBang_H
#define AliFBigBang_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFBigBang                                                          //
//                                                                      //
// helper class to browse generated particles.                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif
#ifndef ROOT_TObjArray
#include <TObjArray.h>
#endif

class AliFBrowsable;

class AliFBigBang : public TNamed {

private:
   TObjArray     *fBrowsables;      //List of browsable particles

public:
                     AliFBigBang();
   virtual          ~AliFBigBang();
   virtual void      Browse(TBrowser *b);
   AliFBrowsable    *GetBrowsable(Int_t i);
   Bool_t            IsFolder() {return kTRUE;}

   ClassDef(AliFBigBang, 0)   //helper class to browse generated particles.
};

#endif









