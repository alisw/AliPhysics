#ifndef ALIANALYSISTASKV0MNORM_H
#define ALIANALYSISTASKV0MNORM_H


class TH1I;
class TH1F;
class TF1;
class TH2F;
class TH2D;
class TH1D;
class TLorentzVector;
class TArrayD;
class TArrayF;
class TArrayL;
class TProfile;
class TList;
class TClonesArray;
class TString;
class AliEmcalJet;
class AliVParticle;
class AliLog;
class AliAnalysisUtils;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;
class AliMultSelection;

namespace PWGJE {
  namespace EMCALJetTasks {
    class AliAnalysisEmcalJetHelperEA;
  }
}

#include <vector>
using std::vector;

#include "AliAnalysisTaskEmcalJet.h"
#include "AliEventCuts.h"
#include "AliFiducialCut.h"
#include "AliEMCALRecoUtils.h"


// ANALYSIS OF V0M 
// Author Filip Krizek   (12 MAR 2024)

namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisTaskV0Mnorm : public AliAnalysisTaskEmcalJet {
   public:


   // ######### CONTRUCTORS/DESTRUCTORS AND STD FUNCTIONS
   AliAnalysisTaskV0Mnorm();
   AliAnalysisTaskV0Mnorm(const char *name);
   virtual  ~AliAnalysisTaskV0Mnorm();
   void     UserCreateOutputObjects();
   void     Terminate(Option_t *);

   static AliAnalysisTaskV0Mnorm* AddTaskV0Mnorm(
       const char* mcpariclearraynamePartMC = "mcparticles", //name of mcparticle TClonesArray array for MC particle level jets
       Bool_t      useVertexCut       = kTRUE,  // vertex cut
       Bool_t      usePileUpCut       = kTRUE, // discard pile up event
       const char* suffix = ""                              //SUBWAGON has to be the last parameter
  );


  // ######### SETTERS/GETTERS

  void        SetMCParticleContainerName(const char* name){ fMyParticleContainerName = name;}
  void        SetUseDefaultVertexCut (Bool_t val) {fUseDefaultVertexCut = val;}
  void        SetUsePileUpCut (Bool_t val) {fUsePileUpCut = val;}


  Bool_t      Run();
  Bool_t      FillHistograms();

  void InitEventProperties();

 private:

  // ######### CHECK FUNCTIONS
  Bool_t      IsEventInAcceptance(AliVEvent* event);

  // ######### STANDARD FUNCTIONS
  void        ExecOnceLocal();



  // ########## USAGE TRIGGERS
  Bool_t      fUseDefaultVertexCut;                   // trigger if automatic vertex cut from helper class should be done
  Bool_t      fUsePileUpCut;                          // trigger if pileup cut should be done


  // ########## SOURCE INFORMATION
  TString     fMyParticleContainerName;               // name of particle level MC particle container

  AliParticleContainer *fParticleContainerPartLevel;  //! particle level container with pythia particles


  // ########## GENERAL ////VARS
  PWGJE::EMCALJetTasks::AliAnalysisEmcalJetHelperEA  *fHelperEA;                   // wrapper for  mean V0 multiplicities
  AliAnalysisUtils*   fHelperClass;                   //! Vertex selection helper                                            
  Bool_t              fInitializedLocal;              //! trigger if tracks/jets are loaded  initiates calling   ExecOnce   
 
   TH2D *hEA_correlations;                             //! EA correlations between part. level and det. level V0M/<V0M>


   AliAnalysisTaskV0Mnorm(const AliAnalysisTaskV0Mnorm&);
   AliAnalysisTaskV0Mnorm& operator=(const AliAnalysisTaskV0Mnorm&);

   ClassDef(AliAnalysisTaskV0Mnorm, 1); // Charged jet analysis for pAliAnalysisTaskHJetSpectra/home/fkrizek/z501.ALIC

};
}
}
#endif
