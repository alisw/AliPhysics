#ifndef AliFFruit_H
#define AliFFruit_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFFruit                                                            //
//                                                                      //
// Utility class to draw Electrons, photons, Jets, Clusters,etc         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TPolyLine3D
#include <TPolyLine3D.h>
#endif

class AliFDisplay;

class AliFFruit : public TPolyLine3D {

private:
   TObject          *fFruit;            //Pointer to original fruit
   
public:
                     AliFFruit() {;}
                     AliFFruit(TObject *obj, Float_t eta, Float_t phi, Float_t pt, Int_t type);
   virtual          ~AliFFruit() {;}
   virtual void      Delete(Option_t *option="");
   TObject          *Fruit() {return fFruit;}
   virtual char     *GetObjectInfo(Int_t px, Int_t py);

   ClassDef(AliFFruit, 0)   //Utility class to draw Electrons, photons, Jets, Clusters,etc
};

#endif
