#ifndef AliGenFLUKAsource_H
#define AliGenFLUKAsource_H
/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////
#include "AliGenerator.h"
#include "TNamed.h"
#include "TF1.h"
#include "TArrayF.h"
#include "TTree.h"

// Read background particles from a FLUKA boundary source file

class AliGenFLUKAsource : public AliGenerator
{
 
protected:

  Int_t fIkine;               //Flag to choose type of particles to be read in
                              // 6 - all particles types
                              // 7 - only gammas
                              // 8 - only neutrons
                              // 9 - only charged particles
  Float_t     fAgeMax;        //Maximum age of particle
  Float_t     fAddWeight;     //Add weight for neutrons 
  Float_t     fZshift;        //Shift the Z of impact point by this quantity
  Float_t     fFrac;
  
  const Text_t     *fFileName;          //!Choose the file
   
  TTree           *fTreeFluka;        //pointer to the TTree
//Declaration of variables read from the file -- TTree type
   Float_t         Ip;
   Float_t         Ipp;
   Float_t         Xi;
   Float_t         Yi;
   Float_t         Zi;
   Float_t         Px;
   Float_t         Py;
   Float_t         Pz;
   Float_t         Ekin;
   Float_t         Zv;
   Float_t         Rv;
   Float_t         Itra;
   Float_t         Igas;
   Float_t         Wgt;
   Float_t         Etag;
   Float_t         Ptg;
   Float_t         Age;

public:
  AliGenFLUKAsource();
  AliGenFLUKAsource(Int_t npart);
  virtual ~AliGenFLUKAsource();
  // Initialise 
  virtual void Init() {}
  // Initialise fluka data 
  virtual void FlukaInit();
  // choose particle type
  virtual void SetPartFlag(Int_t ikine) {fIkine=ikine;}
  // set time cut 
  virtual void SetAgeMax(Float_t agemax) {fAgeMax=agemax;}
  // use additional weight on neutrals
  virtual void SetAddWeight(Float_t addwgt) {fAddWeight=addwgt;}
  // z-shift of vertex
  virtual void SetZshift(Float_t zshift) {fZshift=zshift;}
  // set file name of data file
  virtual void SetFileName(const Text_t *filname) {fFileName=filname;}
  // read only fraction of data  
  virtual void SetFraction(Float_t frac=1.){fFrac=frac;}
  // generate event
  virtual void Generate();

  ClassDef(AliGenFLUKAsource,1) //Boundary source
};
#endif






