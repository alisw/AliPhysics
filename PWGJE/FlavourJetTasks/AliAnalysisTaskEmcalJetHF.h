#ifndef ALIANALYSISTASKEMCALJETHF_H
#define ALIANALYSISTASKEMCALJETHF_H
/**
 * \file AliAnalysisTaskEmcalJetHF.h
 * \class AliAnalysisTaskEmcalJetHF
 * \brief Declaration of class AliAnalysisTaskEmcalJetHF.  Task to perform analysis on HFE JEts
 *
 * Task to perform HFE Jet Spectra using EMCal
 *
 * \author Andrew Castro, <andrew.john.castro@cern.ch>, University of Tennessee
 * \date June 06, 2016
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"

//ROOT
class TClonesArray;
class TH1;
class TH2;
class TH3;
class TH1F;
class TH2F;
class TH3F;
class THnSparse;
class TLorentzVector;
class TGraph;
class TClonesArray;
class TArrayI;
class TProfile;
//ALIROOT
class AliEMCALTrack;
class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliESDtrack;
class AliESDtrackCuts;
class AliAnalysisEtCuts;
class AliDetectorPID;
class AliESDCaloCluster;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;
class AliEMCALTriggerPatchInfo;
class AliVCaloTrigger;
class AliPIDResponse;
//INCLUDES
#include <TRef.h>
#include <TBits.h>
#include <TMath.h>
#include <AliEmcalJet.h>
#include <AliVEvent.h>
#include <AliVTrack.h>
#include <AliVCluster.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TRandom3.h>
#include <AliLog.h>
#include "AliAnalysisTaskEmcalJet.h"
#include "AliEMCALPIDResponse.h"
#include <AliESDCaloCluster.h>
#include <AliESDtrackCuts.h>
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliPID.h"
#include "AliTPCdEdxInfo.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliEmcalTriggerSetupInfo.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliESDpid.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisTaskSE.h"
#include "AliStack.h"


class AliAnalysisTaskEmcalJetHF : public AliAnalysisTaskEmcalJet {
public:
    
    AliAnalysisTaskEmcalJetHF();//!
    AliAnalysisTaskEmcalJetHF(const char *name);//!
    virtual ~AliAnalysisTaskEmcalJetHF() ;//!
    
    void                        SetAODAnalysis()       { SetBit(kAODanalysis, kTRUE); };//!
    void                        SetESDAnalysis()       { SetBit(kAODanalysis, kFALSE); };//!
    Bool_t                      IsAODanalysis() const  { return TestBit(kAODanalysis); };//!

    void                        UserCreateOutputObjects();//!
    void                        Terminate(Option_t *option);//!
    void                        EMCalClusterInfo();//!
    void                        SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec);//!
    
protected:
    void                        ExecOnce();//!
    Bool_t                      FillHistograms();//!
    Bool_t                      Run();//!

    void                        AllocateJetHistograms();//!
    void                        AllocateTrackHistograms();//!
    void                        AllocateClusterHistograms();//!
    void                        AllocateCellHistograms();//!

    void                        DoJetLoop();//!
    void                        DoTrackLoop();//!
    void                        DoClusterLoop();//!
    void                        DoCellLoop();//!

    THistManager                fHistManager;///< Histogram manager

private:
    enum{
        kAODanalysis = BIT(20),
    };

    // Objects
    AliVEvent                     *fVevent;//!                                 event object
    AliESDEvent                   *fESD;//!                                    ESD object
    AliAODEvent                   *fAOD;//!                                    AOD object
    AliPIDResponse                *fpidResponse;//!                            pid response

    //Consts.
    Double_t                      fEventCounter;//!                            Event Counter
    Double_t                      fInvmassCut;//!                              Invariant mass Cut for photons

    // Histograms
    TH2F                          *fdEdx;//!                                   dedx vs pt
    TH2F                          *fM20;//!                                    pt vs M20
    TH2F                          *fM02;//!                                    pt vs M20
    TH2F                          *fM20EovP;//!                                M20 vs E/p
    TH2F                          *fM02EovP;//!                                M20 vs E/p
    TH1F                          *fInvmassLS;//!                              Invmass of LS pairs
    TH1F                          *fInvmassULS;//!                             Invmass of ULS pairs
    TH1F                          *fEMCTrketa;//!                              EMCal trk eta
    TH1F                          *fEMCTrkphi;//!                              EMCal trk phi
    TH1F                          *fHistJetEovP;//!                            Jet EovP
    TH2F                          *fHistJetEovPvPt;//!                         Jet EovP vs Pt
    TH1F                          *fHistJetEovPvEP;//!                         Jet EovP vs Event Plane
    TH1F                          *fHistClusEovP;//!                           EovP using Cluster Energy
    TH1F                          *fHistClusEovPnonlin;//!                     EovP using Non linear corrected energy
    TH1F                          *fHistClusEovPHadCorr;//!                    EovP using Hadron Corected Energy
    TH1F                          *fHistClusTrackMatchdPhi;//!                 dPhi = Phi_cluster - Phi_track Matching
    TH1F                          *fHistClusTrackMatchdEta;//!                 dEta = Eta_cluster - Eta_track Matching

    AliAnalysisTaskEmcalJetHF(const AliAnalysisTaskEmcalJetHF&)           ;//! not implemented
    AliAnalysisTaskEmcalJetHF &operator=(const AliAnalysisTaskEmcalJetHF&);//! not implemented

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskEmcalJetHF, 7);
    /// \endcond
};
#endif
