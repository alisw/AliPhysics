/// \ingroup "PWG3muon"
/// \class AliAnalysisTaskSingleMu
/// \brief Analysis task for single muons in the spectrometer
///
//  Author Diego Stocco

#include "AliAODEvent.h"
#include "AliAODTrack.h"

class AliAnalysisTaskSingleMu : public AliAnalysisTask {
 public:
  AliAnalysisTaskSingleMu(const char *name = "AliAnalysisTaskSingleMu");
  virtual ~AliAnalysisTaskSingleMu() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

 protected:
  Bool_t FillTrackVariables(AliAODTrack &muonTrack);
  
  void InitVariables();

 private:
  AliAnalysisTaskSingleMu(const AliAnalysisTaskSingleMu&);
  AliAnalysisTaskSingleMu& operator=(const AliAnalysisTaskSingleMu&);

  AliAODEvent *fAOD; //!< ESDevent object

  TTree *fResults; //!< Tree with results

  enum {
    kVarPt,     //!< Muon pt
    kVarY,      //!< Muon rapidity
    kVarPhi,    //!< Muon phi
    kVarVz,     //!< Primary vertex longitudinal position
    kVarDCA,    //!< Transverse distance at vertex
    kNfloatVars
  };

  enum {
    kVarTrig,   //!< Matched trigger
    kNintVars
  };
  

  Float_t* fVarFloat; //!< Array of float variables
  Int_t* fVarInt;     //!< Array of int variables
  
  TString* fFloatVarName; //!< Float variable names for branches
  TString* fIntVarName;   //!< Intt variable names for branches

  ClassDef(AliAnalysisTaskSingleMu, 0); // Single muon analysis
};

