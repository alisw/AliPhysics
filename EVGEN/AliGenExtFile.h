#ifndef AliGenExtFile_H
#define AliGenExtFile_H
/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 0    //
/////////////////////////////////////////////////////////
#include "AliGenerator.h"
#include "TNamed.h"
#include "TF1.h"
#include "TArrayF.h"
#include "TTree.h"

// Read background particles from a FLUKA boundary source file

class AliGenExtFile : public AliGenerator
{
 
protected:
  const Text_t     *fFileName;         //! Choose the file
  Int_t           fNcurrent;           // points to the next entry
  TTree           *fTreeNtuple;        // pointer to the TTree
//Declaration of variables read from the file -- TTree type
  //Declaration of leaves types
   Int_t           Nihead;
   Int_t           Ihead[12];
   Int_t           Nrhead;
   Float_t         Rhead[6];
   UInt_t          Idpart;
   Float_t         Theta;
   Float_t         Phi;
   Float_t         P;
   Float_t         E;
public:
   AliGenExtFile();
  AliGenExtFile(Int_t npart);
  virtual ~AliGenExtFile();
  // Initialise 
  virtual void Init() {}
  // Initialise fluka data 
  virtual void NtupleInit();
  // set file name of data file
  virtual void SetFileName(const Text_t *filname) {fFileName=filname;}
  // generate event
  virtual void Generate();

  ClassDef(AliGenExtFile,1) //Boundary source
};
#endif






