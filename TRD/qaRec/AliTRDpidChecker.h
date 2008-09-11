#ifndef ALITRDPIDCHECKER_H
#define ALITRDPIDCHECKER_H

//////////////////////////////////////////////////////
//
// Task to check PID performance of the TRD
//
// Author : Alex Wilk <wilka@uni-muenster.de>
//
///////////////////////////////////////////////////////

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TObjArray;
class TList;
class TClonesArray;
class TTreeSRedirector;
class AliTRDReconstructor;
class AliTRDpidChecker : public AliTRDrecoTask 
{
public:
  AliTRDpidChecker();
  virtual ~AliTRDpidChecker();
  
  void    CreateOutputObjects();
  void    Exec(Option_t *option);
  Bool_t  PostProcess();
  void    Terminate(Option_t *);

private:
  AliTRDpidChecker(const AliTRDpidChecker&); // not implemented
  AliTRDpidChecker& operator=(const AliTRDpidChecker&); // not implemented

  Double_t GetPionEfficiency(Int_t Index1, Int_t Index2);  // calculates the pion efficiency
  Double_t GetError(Int_t Index1, Int_t Index2);           // calculates the error
  

  AliTRDReconstructor *fReconstructor;     //! reconstructor needed for recalculation the PID

  ClassDef(AliTRDpidChecker, 1); // TRD PID checker
};

#endif
