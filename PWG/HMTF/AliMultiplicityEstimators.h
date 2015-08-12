#ifndef AliMultiplicityEstimators_cxx
#define AliMultiplicityEstimators_cxx

//class TH1F;
//class TH1I;
//class TGraphErrors;
class AliMCEvent;
class AliMCParticle;
class AliHeader;
class AliStack;
class AliMCParticle;

#include "TNamed.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"
#include "TNtuple.h"

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
  kNPID
};

class MultiplicityEstimatorBase : public TNamed {
 public:
  MultiplicityEstimatorBase();
  MultiplicityEstimatorBase(const char* name, const char* title);
  virtual ~MultiplicityEstimatorBase() {}

  void RegisterHistograms(TList* outputList);
  //get the ending of the name common to all histograms from this estimator:
  TString GetNamePostfix() {return TString("_") + fName;};
  //get the ending of the title common to all histograms from this estimator:
  TString GetTitlePostfix() {return TString(" ") + fTitle;};
  virtual void PreEvent(AliMCEvent* event) = 0;
  virtual void ProcessTrackForMultiplicityEstimation(AliMCParticle* track) = 0;
  virtual void ProcessTrackWithKnownMultiplicity(AliMCParticle* track) = 0;
  virtual void PostEvent() = 0;
  virtual void Terminate(TList* sum) = 0;
  const Bool_t fuseWeights;   // Use event weights?

 protected:
  /*
    Header dependent logic is done here. This function is supposed to be called in the
    beginning of PreEvent().
   */
  void ReadEventHeaders(AliMCEvent* event);
  Int_t festimator_bins;
  // histograms for weighted [0] and unweighted [1] are created where appropriate
  TH2F  *feta_Nch;          // dNdEta distributions; multiplicity is on the y-axis
  TH3F  *fNch_pT_pid;  // multiplicity class; pT; pid
  TNtuple *fEventTuple;

  AliHeader *fheader;       // Event header
  AliMCEvent *fevent;       // current event
  AliStack  *fstack;
  Float_t feventWeight;        // weight of the event read from the header
  Int_t fnMPI;          // Number of parton interactions

 private:
  MultiplicityEstimatorBase(const MultiplicityEstimatorBase&);           // not implemented
  MultiplicityEstimatorBase& operator=(const MultiplicityEstimatorBase&);// not implemented


  ClassDef(MultiplicityEstimatorBase, 1); // example of analysis
};


class EtaBase : public MultiplicityEstimatorBase {
 public:
  EtaBase();
  EtaBase(const char* name, const char* title, Float_t eta_min_backwards, Float_t eta_max_backwards,
	  Float_t eta_min_forwards, Float_t eta_max_forwards);
  EtaBase(const char* name, const char* title);
 protected:
  void PreEvent(AliMCEvent* event);
  void ProcessTrackForMultiplicityEstimation(AliMCParticle* track);
  void ProcessTrackWithKnownMultiplicity(AliMCParticle *track);
  void PostEvent();
  void Terminate(TList* sum);
  Int_t fnch_in_estimator_region;   // counter for charged particles in current event
  Int_t fn_pid_in_event[kNPID];     // counter for PID'ed particles in this event
  Float_t feta_min_forwards, feta_max_forwards;  // range in eta for mult. estimation (positive values)
  Float_t feta_min_backwards, feta_max_backwards;  // range in eta for mult. estimation (positive values)
  Bool_t fbypass_eta_selection;  // bypass the eta selection (used to get full eta range
  ClassDef(EtaBase, 1)
};

#endif
