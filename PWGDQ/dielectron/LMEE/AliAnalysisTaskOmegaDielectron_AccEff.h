/// \classAliAnalysisTaskOmegaDielectron_AccEfff
/// \brief Task to for pT spectra vs. multiplicity analysis


#ifndef AliAnalysisTaskOmegaDielectron_AccEff_H // header guard in case of multiple includes
#define AliAnalysisTaskOmegaDielectron_AccEff_H

#define MAX_HISTO_DIM 4
#define PRECISION 1e-6


#include "TList.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"

#include "AliPID.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"


#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliEventCuts.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisCuts.h"
#include "AliDielectronVarCuts.h"
#include "AliDielectronV0Cuts.h"
#include "AliDielectronPID.h"
#include "AliDielectronVarManager.h"
#include "AliDielectronTrackCuts.h"
#include "AliDielectronCutGroup.h"
#include "AliDielectronEventCuts.h"
#include "AliAODInputHandler.h"


#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliMCParticle.h"
#include "AliESDtrackCuts.h"

#include <iostream>
#include <vector>
using namespace std;
using std::string;
using std::vector;
using std::array;

class AliAnalysisFilter;
class AliPIDResponse;


class AliAnalysisTaskOmegaDielectron_AccEff : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskOmegaDielectron_AccEff();
  AliAnalysisTaskOmegaDielectron_AccEff(const char *name);
  virtual ~AliAnalysisTaskOmegaDielectron_AccEff();

  // static AliAnalysisTaskOmegaDielectron_AccEff* AddTaskMultDepSpec(TString controlstring, Int_t cutModeLow = 100, Int_t cutModeHigh = 121);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t* option);
  virtual void   Terminate(Option_t*);

  // Track cuts setter
  void   AddTrackCut(AliAnalysisFilter* filter) {fFilter_TrackCuts.push_back(filter);}
  void   AddPIDCut(AliAnalysisFilter* filter_PID) {fFilter_PID.push_back(filter_PID);}

  //Bool flag for sys Unc. -- cut variations
  void   SetBoolsysUncOutput(Bool_t sysUncOutput) {fsysUnc = sysUncOutput;}


private:

  TList*              fOutputList;		  //!<! Output list
  TList*              fCutVariList;		  //!<! Output list
  AliVEvent*          fEvent;			      //!<! Event object
  AliMCEvent*         fMCEvent;         //!<! MC event
  AliDielectronEventCuts  *fEventCuts;       //!<! Event cuts
  AliDielectronEventCuts  *fEventCuts_VertexZ;       //!<! Event cuts
  std::vector<AliAnalysisFilter*> fFilter_TrackCuts;
  std::vector<AliAnalysisFilter*> fFilter_PID;
  AliPIDResponse*     fPIDResponse;           //!<! Analysis Filter

  // Acceptance cuts for tracks
  Double_t              fMinEta;        ///< Minimum eta cut
  Double_t              fMaxEta;        ///< Maximum eta cut
  Double_t              fMinPt;			    ///< Minimum pT cut
  Double_t              fMaxPt;			    ///< Maximum pT cut
  Bool_t                fsysUnc;        ///< sys unc cut variation

  // pdg codes:
  Int_t                felectron_pdg;   //!
  Int_t                fpositron_pdg;   //!
  Int_t                fmother_pdg;     //!

  //storage vectors:
  vector<AliAODTrack *>   v_elec_true_omega;              //! array of strings containing the electron track from a true omega dielectron decay
  vector<AliAODTrack *>   v_posi_true_omega;              //! array of strings containing the positron track from a true omega dielectron decay
  vector<AliMCParticle *> v_elec_true_omega_MCPart;       //! array of strings containing the electron MCParticle from a true omega dielectron decay
  vector<AliMCParticle *> v_posi_true_omega_MCPart;       //! array of strings containing the positron MCParticle from a true omega dielectron decay
  vector<Int_t>           v_elec_motherID_true_omega;     //! array of strings containing the mother id of the electron
  vector<Int_t>           v_posi_motherID_true_omega;     //! array of strings containing the mother id of the positron


  // Output Histograms
  std::vector<TH3D*> fHist_MC_Omegas_gen;                 //!<! Histogram of generated omegas for different Event Cut Settings
  TH3D* fHist_MC_Omegas_gen_DaughtersinAcc;               //!<! Histogram of generated omegas with dielectron daughters in Acc
  std::vector<TH3D*> fHist_Rec_Omegas_TrackPID;           //!<! Histogram of reconstructed omegas within track and PID cuts

  TH2D* fHistVertex;                                      //!<! Histogram for z vertex distribution different EventCut stages
  TH1F* fHist_MC_Omegas_Rapidity;                         //!<! Histogram for omega rapidity distribution

  TH3D* fHist_MC_elec_gen;                                //!<! Histogram of generated electrons
  TH3D* fHist_MC_posi_gen;                                //!<! Histogram of generated positrons
  TH3D* fHist_MC_elec_gen_inAcc;                          //!<! Histogram of generated electrons in Acceptance
  TH3D* fHist_MC_posi_gen_inAcc;                          //!<! Histogram of generated positrons in Acceptance

  TH3D* fHist_MC_Omegas_withoutCuts;                      //!<! Histogram of generated Omegas (MCParticles)
  TH3D* fHist_Rec_Omegas_withoutCuts;                     //!<! Histogram of Omegas (AliVTracks)

  std::vector<TH3D*> fHist_elec_rec_inAcc;                //!<! Histogram of reconstructed electrons in Acceptance
  std::vector<TH3D*> fHist_elec_rec_inAcc_Track;          //!<! Histogram of reconstructed electrons within Acceptance, TrackCuts
  std::vector<TH3D*> fHist_elec_rec_inAcc_Track_PID;      //!<! Histogram of reconstructed electrons within Acceptance, TrackCuts and PID
  std::vector<TH3D*> fHist_posi_rec_inAcc;                //!<! Histogram of reconstructed positrons within Acceptance
  std::vector<TH3D*> fHist_posi_rec_inAcc_Track;          //!<! Histogram of reconstructed positrons within Acceptance, TrackCuts
  std::vector<TH3D*> fHist_posi_rec_inAcc_Track_PID;      //!<! Histogram of reconstructed positrons within Acceptance, TrackCuts and PID

  std::vector<TH3D*> fHist_MC_Omegas_TrackCuts;           //!<! Histogram of generated Omegas within TrackCuts
  std::vector<TH3D*> fHist_MC_Omegas_TrackPID;            //!<! Histogram of generated Omegas within TrackCuts, PID
  std::vector<TH3D*> fHist_Rec_Omegas_TrackCuts;          //!<! Histogram of reconstructed Omegas within TrackCuts



  void    SetPIDResponse(AliPIDResponse *fPIDRespIn)        {fPIDResponse = fPIDRespIn;} //!<! pid response



  Bool_t AcceptKinematics(AliVParticle* particle);                            //!<! Accept kinematic
  Bool_t CheckDielectronDecay(AliMCParticle *particle, Bool_t checkacc);      //!<! Check if particle has e+e- as daughters , with bool for die kinematic acceptance check for the daughters

  AliAnalysisTaskOmegaDielectron_AccEff(const AliAnalysisTaskOmegaDielectron_AccEff&); // not implemented
  AliAnalysisTaskOmegaDielectron_AccEff& operator=(const AliAnalysisTaskOmegaDielectron_AccEff&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskOmegaDielectron_AccEff, 1); // example of analysis
  /// \endcond
};

#endif
