#ifndef ITS_H
#define ITS_H
////////////////////////////////////////////////
//  Manager and hits classes for set: ITS     //
////////////////////////////////////////////////

#include "AliDetector.h"
#include "AliHit.h" 
#include "AliDigit.h" 

class AliITS : public AliDetector {
  
protected:
  
  Bool_t fEuclidOut;
  
  Int_t fIdSens1;    //First  layer identifier
  Int_t fIdSens2;    //Second layer identifier
  Int_t fIdSens3;    //Third  layer identifier
  Int_t fIdSens4;    //Fourth layer identifier
  Int_t fIdSens5;    //Fifth  layer identifier
  Int_t fIdSens6;    //Sixth  layer identifier
  
public:
  AliITS();
  AliITS(const char *name, const char *title);
  virtual       ~AliITS();
  virtual void   AddHit(Int_t, Int_t*, Float_t*);
  virtual void   AddDigit(Int_t*, Int_t*);
  virtual void   BuildGeometry();
  virtual void   CreateGeometry() {}
  virtual void   CreateMaterials();
  virtual Int_t  IsVersion() const =0;
  Int_t          DistancetoPrimitive(Int_t px, Int_t py);
  virtual void   Init();
  virtual void   MakeBranch(Option_t *opt=" ");
  virtual void   SetEUCLID(Bool_t euclid=1);
  virtual void   StepManager()=0;
  
  ClassDef(AliITS,1)  //Hits manager for set:ITS
};

 
//___________________________________________
class AliITSdigit: public AliDigit  {
public:
   Int_t       fEvent;      // Event number
   Int_t       fLayer;      // Layer number
   Int_t       fDet ;       // Detector number
   Int_t       fNoverl;     // Number of overflow
 
public:
   AliITSdigit() {}
   AliITSdigit(Int_t *tracks, Int_t *digits);
   virtual ~AliITSdigit() {}
 
   ClassDef(AliITSdigit,1)  //Digit (Header) object for set:ITS
};
 
//___________________________________________
 
class AliITShit : public AliHit {
public:
   Int_t     fLayer;      // Layer number
   Int_t     fLadder;     // Ladder number
   Int_t     fDet;        // Detector number  
   Float_t   fPx  ;       //PX
   Float_t   fPy  ;       //PY
   Float_t   fPz  ;       //PZ
   Float_t   fDestep;     // Energy deposited in the current step
public:
   AliITShit() {}
   AliITShit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
   virtual ~AliITShit() {}
 
   ClassDef(AliITShit,1)  //Hits object for set:ITS
};
 
#endif
