#ifndef ALISTARTLOADER_H
#define ALISTARTLOADER_H

#include <AliLoader.h>
#include "AliDataLoader.h"
#include "AliSTARTdigit.h"
#include "AliSTARTvertex.h"
class AliSTART;
class AliSTARTLoader: public AliLoader
 {
 public:
   AliSTARTLoader(){};
   AliSTARTLoader(const Char_t *name,const Char_t *topfoldername);
   AliSTARTLoader(const Char_t *name,TFolder *topfolder);
   
   virtual ~AliSTARTLoader();
   // Digits
   virtual void   CleanDigits() {fDigitsDataLoader.GetBaseDataLoader()->Clean();}
   Int_t          LoadDigits(Option_t* opt=""){return fDigitsDataLoader.GetBaseDataLoader()->Load(opt);}
   
   void           UnloadDigits(){fDigitsDataLoader.GetBaseDataLoader()->Unload();}
   virtual Int_t  WriteDigits(Option_t* opt=""){return fDigitsDataLoader.GetBaseDataLoader()->WriteData(opt);}
   virtual Int_t PostDigits(AliSTARTdigit *dgt){return fDigitsDataLoader.GetBaseDataLoader()->Post(dgt);}
   
   TObject*  Digits(){ return fDigitsDataLoader.GetBaseDataLoader()->Get();} // returns a pointer to the tree of  RawClusters
   
   
   //Vertices
   virtual void   CleanRecPoints() {fVertexDataLoader.GetBaseDataLoader()->Clean();}
   Int_t          LoadRecPoints(Option_t* opt=""){return fVertexDataLoader.GetBaseDataLoader()->Load(opt);}
   void           UnloadRecPoints(){fVertexDataLoader.GetBaseDataLoader()->Unload();}
   virtual Int_t  WriteRecPoints(Option_t* opt=""){return fVertexDataLoader.GetBaseDataLoader()->WriteData(opt);}
   //  AliSTARTVertex *GetVertex(){return static_cast <AliSTARTVertex*>(fVertexDataLoader.GetBaseDataLoader()->Get());}
   virtual Int_t PostRecPoints(AliSTARTvertex *vrt){return fVertexDataLoader.GetBaseDataLoader()->Post(vrt);}
   //  TObject*  Vertex(){ return fVertexDataLoader.GetBaseDataLoader()->Get();} // returns a pointer to the tree of  RawClusters
   
   protected:
   // methods
   
   Int_t     PostDigits(){return fDigitsDataLoader.GetBaseDataLoader()->Post();}
   Int_t     PostRecPoints(){return fVertexDataLoader.GetBaseDataLoader()->Post();}
   

   // DATA
   AliDataLoader fDigitsDataLoader; //digits loader
   
   AliDataLoader fVertexDataLoader;  // RecPoints (vertex position) loader
   
 public:
   ClassDef(AliSTARTLoader,1)
      };
 
#endif



