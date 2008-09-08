#ifndef AliTRDpidChecker_cxx
#define AliTRDpidChecker_cxx

// Task to check PID performance of the TRD

class AliTRDReconstructor;


#include "AliAnalysisTask.h"

class TObjArray;
class TList;
class TClonesArray;
class TTreeSRedirector;

class AliTRDpidChecker : public AliAnalysisTask {
 public:
  AliTRDpidChecker(const char *name = "AliTRDpidChecker");
  virtual ~AliTRDpidChecker();
  
  void   ConnectInputData(Option_t *);
  void   CreateOutputObjects();
  void   Exec(Option_t *option);
  void   Terminate(Option_t *);
/*   Int_t  GetDebugLevel() const {return fDebugLevel;}  */
/*   void   SetDebugLevel(Int_t debug){fDebugLevel = debug;} */
 
 private:
  AliTRDpidChecker(const AliTRDpidChecker&); // not implemented
  AliTRDpidChecker& operator=(const AliTRDpidChecker&); // not implemented

  Double_t GetPionEfficiency(Int_t Index1, Int_t Index2);  // calculates the pion efficiency
  Double_t GetError(Int_t Index1, Int_t Index2);           // calculates the error
  
  TObjArray        *fObjectContainer;       // Container
  TObjArray        *fTracks;                // Array of tracks

  AliTRDReconstructor *fReconstructor;     // reconstructor needed for recalculation the PID
/*   Int_t            fDebugLevel;         // Debug level */
/*   TTreeSRedirector *fDebugStream;       // Debug stream */

  ClassDef(AliTRDpidChecker, 1); // example of analysis
};

#endif
