#ifndef ALIPARTICLETREEHANDLERAPPLY_H
#define ALIPARTICLETREEHANDLERAPPLY_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliParticleTreeHandlerApply
// \brief helper class to handle a tree for Dstar cut optimisation and MVA analyses
// \authors:
// L. Vermunt, luuk.vermunt@cern.ch
/////////////////////////////////////////////////////////////

#include <TTree.h>
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"

class AliParticleTreeHandlerApply : public TObject
{
public:

  enum selection_type {
    kSelectedProton         = BIT(0),
    kSelectedKaon           = BIT(1),
    kSelectedPion           = BIT(2),
    kSelectedNClsTPCStd     = BIT(3),
    kSelectedNClsTPCTight   = BIT(4),
    kSelectedProtonPIDLoose = BIT(5),
    kSelectedProtonPIDTight = BIT(6),
    kSelectedKaonPIDLoose   = BIT(7),
    kSelectedKaonPIDTight   = BIT(8),
    kSelectedPionPIDLoose   = BIT(9),
    kSelectedPionPIDTight   = BIT(10)
  };

  AliParticleTreeHandlerApply();
  virtual ~AliParticleTreeHandlerApply();

  AliParticleTreeHandlerApply(const AliParticleTreeHandlerApply &source) = delete;
  AliParticleTreeHandlerApply& operator=(const AliParticleTreeHandlerApply &source) = delete;

  TTree* BuildTree(TString name="tree", TString title="tree");
  bool SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, AliAODTrack* tr, AliAODTrack* trGlobal);
  bool SetMCGenVariables(AliAODMCParticle* mcpart, int mc_label);
  bool SetPIDVariables(float nsig_TPC, float nsig_TOF, float nsig_comb_TPCTOF);

  //to be called for each candidate
  void SetSelectionType(int part, bool isstd, bool ispidloose, bool ispidtight); //0=proton, 1=kaon, 2=pion, 3=NClsTPC
  void FillTree() { //to be called for each candidate!
    fTreeVar->Fill();
    fTrackSel = 0;
  }

  void SetOnlyDedicatedBranches(bool b) { fOnlyDedicatedBranches = b; }
  void SetIsMC(bool b) { fIsMC = b; }
  void SetIsDebugMode(bool b) { fDebugMode = b; }

private:

  TTree* fTreeVar;                                                     //!<! tree with variables

  int fTrackSel;                                                       /// flag for track type (bit map above)
  float fPt;                                                           /// track pt
  float fY;                                                            /// track rapidity
  float fEta;                                                          /// track pseudorapidity
  float fPhi;                                                          /// track azimuthal angle
  float fPtGen;                                                        /// track generated pt
  float fYGen;                                                         /// track generated rapidity
  float fEtaGen;                                                       /// track generated pseudorapidity
  float fPhiGen;                                                       /// track generated azimuthal angle
  int fCharge;                                                         /// track charge
  int fID;                                                             /// track ID
  int fPDG;                                                            /// PDG code track
  int fMCLabel;                                                        /// track MC label
  int fEvID;                                                           /// event ID corresponding to the one set in fTreeEvChar, first 32 bit of fEvIDLong
  int fEvIDExt;                                                        /// event ID corresponding to the one set in fTreeEvChar, second 32 bit of fEvIDLong
  Long64_t fEvIDLong;                                                  /// event ID corresponding to the one set in fTreeEvChar, full fEvIDLong
  int fRunNumber;                                                      /// run number

  bool fOnlyDedicatedBranches;                                         /// variable to disable unnecessary branches
  bool fIsMC;                                                          /// flag to include MC branches

  bool fDebugMode;                                                     /// flag to enable debug mode (to validate selection new task wrt old one)
  float fTPCCls;                                                       /// track number of TPC clusters
  float fDCAXY;                                                        /// track DCA XY (without recalculation)
  float fDCAZ;                                                         /// track DCA Z (without recalculation)
  float fDCAXYProp;                                                    /// track DCA XY (with recalculation)
  float fDCAZProp;                                                     /// track DCA Z (with recalculation)
  float fTPCCrRows;                                                    /// track TPC crossed rows
  float fTPCClsToFnd;                                                  /// track ratio TPC clusters over findable
  float fChi2;                                                         /// track chi2
  float fPTPC;                                                         /// track TPC momentum
  float fNSigTPC;                                                      /// track PID N sigma TPC
  float fNSigTOF;                                                      /// track PID N sigma TOF
  float fNSigCombTPCTOF;                                               /// track combined N sigma TPC & TOF PID (sqrt[nsig^2 + nsig^2])

  /// \cond CLASSIMP
  ClassDef(AliParticleTreeHandlerApply,1); ///
  /// \endcond
};
#endif
