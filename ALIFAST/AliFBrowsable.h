#ifndef AliFBrowsable_H
#define AliFBrowsable_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFBrowsable                                                        //
//                                                                      //
// helper class to browse generated particles.                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

class AliFBigBang;

class AliFBrowsable : public TNamed {

private:
   TObject          *fRefObject;       //Referenced object
   AliFBigBang      *fBigBang;         //Pointer to control bigbang object

public:
                     AliFBrowsable();
   virtual          ~AliFBrowsable() {;}
   virtual void      Browse(TBrowser *b);
   Bool_t            IsFolder() const {return kTRUE;}
   virtual void      SetBigBang(AliFBigBang *bigbang) {fBigBang = bigbang;}
   virtual void      SetRefObject(TObject *obj) {fRefObject = obj;}

   ClassDef(AliFBrowsable, 0)   //helper class to browse generated particles.
};

#endif











