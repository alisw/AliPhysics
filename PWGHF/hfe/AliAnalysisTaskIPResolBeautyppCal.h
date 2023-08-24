#ifndef AliAnalysisTaskIPResolBeautyppCal_h
#define AliAnalysisTaskIPResolBeautyppCal_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////
//                                                       //
// Task for Non-HFE reconstruction efficiency in Run 2   //
//                                                       //
//  Author: Deepa Thomas (University of Texas at Austin) //
//          Vivek Kumar Singh (VECC)                     //
///////////////////////////////////////////////////////////

class THnSparse;
class TH2F;
class TLorentzVector;

class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliESDtrack;
class AliAODTrack;
class AliCFManager;
class AliMultSelection;

#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
#include "AliSelectNonHFE.h"
#include "AliAODMCParticle.h"

class AliAnalysisTaskIPResolBeautyppCal : public AliAnalysisTaskSE {
  public:
    enum HijingOrNot {kHijing,kElse};
    enum pi0etaType {kNotIsPrimary, kNoMother, kLightMesons, kBeauty, kCharm, kNoFeedDown};

    AliAnalysisTaskIPResolBeautyppCal();
    AliAnalysisTaskIPResolBeautyppCal(const char *name);
    virtual ~AliAnalysisTaskIPResolBeautyppCal();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    Bool_t  PassEventSelect(AliVEvent *fVevent);

    Bool_t  GetNMCPartProduced();
    Int_t   GetPrimary(Int_t id);
    Int_t   GetPi0EtaType(AliAODMCParticle *part);

    Bool_t  PassTrackCuts(AliAODTrack *atrack);
    
    //void    CheckCentrality(AliAODEvent* fAOD, Bool_t &centralitypass);
    //void    SetCentralitySelection(Double_t centMin, Double_t centMax) {fCentralityMin = centMin; fCentralityMax = centMax;};
    Int_t PhiBin(Double_t phi) const;
  void    SwitchRecalImpPar(Bool_t fSwitchRIP) {fRecalIP = fSwitchRIP;};
  void    RecalImpactParam(const AliAODTrack * const track, Double_t dz[2], Double_t covar[3]);
  AliAODVertex*   RemoveDaughtersFromPrimaryVtx(const AliAODTrack * const track);
    
    Bool_t  IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromHijing, Int_t &type, Int_t &iMom, Int_t &MomPDG, Double_t &MomPt);

  private:
    AliVEvent 		    *fVevent;//!V event object
    AliAODEvent 		*fAOD;//!AOD object
    const AliVVertex    *fpVtx; //!
    AliPIDResponse      *fpidResponse; //!pid response
    AliMultSelection    *fMultSelection;//!
    AliAODMCHeader      *fMCHeader;//!
    TClonesArray        *fMCArray;//!

    Bool_t               fRecalIP;//

    Double_t            fCentrality;//!
    Double_t            fCentralityMin;//
    Double_t            fCentralityMax;//
    Double_t            fMultiplicity;//!
    Int_t               fNTotMCpart; //! N of total MC particles produced by generator
    Int_t               fNpureMC;//! N of particles from main generator (Hijing/Pythia)
    Int_t               fNembMCpi0; //! N > fNembMCpi0 = particles from pi0 generator
    Int_t               fNembMCeta; //! N > fNembMCeta = particles from eta generator
    Double_t            fTPCnSigma;//!
    Int_t               ftype;//!

    TList       	   	*fOutputList;//!output list
    TH1F                *fNevents;//!no of events
    TH1F                *fVtxZ;//!
    TH1F                *fVtxX;//!
    TH1F                *fVtxY;//!
    TH1F                *fCentralityNoPass;//!
    TH1F                *fMultiplicityNoPass;//!
    TH2F                *fCentMultiplicityNoPass;//!
    TH1F                *fCentralityPass;//!
    TH1F                *fMultiplicityPass;//!
    TH2F                *fCentMultiplicityPass;//!
    TH1F                *fNegTrkIDPt;//!
    TH1F                *fTrkPt;//!
    TH1F                *fTrketa;//!
    TH1F                *fTrkphi;//!
    TH2F                *fdEdx;//!
    TH2F                *fTPCnsig;//!
    TH1F                *fTrkPt_Ele;//!
    TH1F                *fTrkPt_HFEle;//!
    TH1F                *fTrkPt_NHFEle;//!
    TH1F                *fTrkPt_GammaE;//!
    TH1F                *fTrkPt_DalitzE;//!
    
    THnSparseF          *fImpParSprs_All; //!<! sparse
    THnSparseF          *fImpParSprs_AllE; //!<! sparse
    THnSparseF          *fImpParSprs_HFEle; //!<! sparse
    THnSparseF          *fImpParSprs_NHFEle; //!<! sparse
    THnSparseF          *fImpParSprs_GammaE; //!<! sparse
    THnSparseF          *fImpParSprs_DalitzE; //!<! sparse


    AliAnalysisTaskIPResolBeautyppCal(const AliAnalysisTaskIPResolBeautyppCal&); // not implemented
    AliAnalysisTaskIPResolBeautyppCal& operator=(const AliAnalysisTaskIPResolBeautyppCal&); // not implemented

    ClassDef(AliAnalysisTaskIPResolBeautyppCal, 1); //!example of analysis
};
#endif
