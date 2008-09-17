#ifndef ALITRDPIDREFMAKER_H
#define ALITRDPIDREFMAKER_H

//////////////////////////////////////////////////////
//
// Task to build PID reference tree for the training
// of neural networs for the TRD PID
//
// Author : Alex Wilk <wilka@uni-muenster.de>
//
///////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class AliTRDtrackV1;
class AliTRDReconstructor;
class AliTRDpidRefMaker : public AliTRDrecoTask
{

 public:
  AliTRDpidRefMaker();
  virtual ~AliTRDpidRefMaker();
  
  void    CreateOutputObjects();
  void    Exec(Option_t *option);
  void    GetRefFigure(Int_t ifig, Int_t &first, Int_t &last, Option_t *opt);  
  Bool_t  PostProcess();
  void    Terminate(Option_t *);


 private:
  AliTRDpidRefMaker(const AliTRDpidRefMaker&);               // not implemented
  AliTRDpidRefMaker& operator=(const AliTRDpidRefMaker&);    // not implemented

  void GetV0info(AliTRDtrackV1 *TRDtrack, Float_t *v0pdg);           // get the v0 information

  AliTRDReconstructor *fReconstructor;     //! reconstructor needed for recalculation the PID

  ClassDef(AliTRDpidRefMaker, 1); // TRD reference  maker for NN
};

#endif
