#ifndef AliObservableClassifierpTPID_cxx
#define AliObservableClassifierpTPID_cxx

#include "TH3F.h"

#include "AliObservableBase.h"

class AliObservableClassifierpTPID : public AliObservableBase {
 public:
  AliObservableClassifierpTPID();
  AliObservableClassifierpTPID(AliEventClassifierBase *classifier);
  ~AliObservableClassifierpTPID() {};

  void Fill(AliMCEvent *event, AliStack *stack);
 private:
  enum {
    kPROTON,
    kANTIPROTON,
    kLAMBDA,
    kANTILAMBDA,
    kK0S,
    kKPLUS,
    kKMINUS,
    kPIPLUS,
    kPIMINUS,
    kPI0,
    kXI,
    kANTIXI,
    kOMEGAMINUS,
    kOMEGAPLUS,
    kALLCHARGED,
    kLAMBDA0B,
    kANITLAMBDA0B,
    kNPID
  };

  TH3F *fhistogram;
  AliEventClassifierBase *fclassifier;
  Int_t Pid_enum_to_pdg(Int_t pid_enum);

  ClassDef(AliObservableClassifierpTPID, 1);
};

#endif
