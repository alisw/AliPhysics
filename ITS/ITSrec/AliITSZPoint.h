#ifndef ALIITSZPOINT_H
#define ALIITSZPOINT_H

#include<TObject.h>

///////////////////////////////////////////////////////////////////
//                                                               //
// Class used by AliITSVertexerZ                                 //
// Contains Z coordinates with their error                       //
// is sortable                                                   //
//                                                               //
///////////////////////////////////////////////////////////////////

class AliITSZPoint : public TObject {

 public:

  AliITSZPoint();      
  AliITSZPoint(Float_t z, Float_t ez);
  virtual ~AliITSZPoint();
  virtual Bool_t IsEqual(const TObject *obj) const 
    {return fZ == ((AliITSZPoint*)obj)->fZ;}
  virtual Bool_t      IsSortable() const { return kTRUE; }
  virtual Int_t       Compare(const TObject *obj) const 
    {if(fZ<((AliITSZPoint*)obj)->fZ) return -1;
    else if(fZ>((AliITSZPoint*)obj)->fZ) return 1;
    else return 0; }
  virtual void Print(Option_t *option="") const;
  Float_t GetZ() const {return fZ;}
  Float_t GetErrZ() const {return fErrz;}

 protected:
  Float_t fZ;          // Z coordinate on beam axiz
  Float_t fErrz;       // error on Z

  ClassDef(AliITSZPoint,1);
};
#endif

