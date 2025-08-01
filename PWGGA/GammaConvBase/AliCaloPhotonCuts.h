#ifndef ALICALOPHOTONCUTS_H
#define ALICALOPHOTONCUTS_H
#include <TObjString.h>
#include "AliConversionPhotonBase.h"
#include "AliAODConversionMother.h"
#include "AliAnalysisCuts.h"
#include "AliAnalysisUtils.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliEMCALGeometry.h"
#include "AliPHOSGeometry.h"
#include "AliAODCaloCluster.h"
#include "AliEMCALRecoUtils.h"
#include "AliCalorimeterUtils.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "TProfile2D.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TF1.h"
#include "AliAnalysisManager.h"
#include "AliCaloTrackMatcher.h"
#include "AliPhotonIsolation.h"
#include "AliGammaConvEventMixing.h"
#include <vector>

class AliESDEvent;
class AliAODEvent;
class AliConversionPhotonBase;
class TH1F;
class TH2F;
class TF1;
class AliPIDResponse;
class AliAnalysisCuts;
class iostream;
class TList;
class AliAnalysisManager;
class AliAODMCParticle;
class AliCaloTriggerMimicHelper;


/**
 * @class AliCaloPhotonCuts
 * @brief Class handling all kinds of selection cuts for Gamma Calo analysis
 * @author Friederike Bock
 * @author Daniel Muehlheim
 * @ingroup GammaConv
 *
 * The cut configuration is set as a string with an 19 digit number.
 * Each digit in the string corresponds to a certain cut type, while
 * its values represent the cut values. The cut configuration is listed here:
 *
 * | Position in the cut string (from the end) | Cut type               |
 * |-------------------------------------------|------------------------|
 * |                  0                        | Cluster Type           |
 * |                  1                        | Eta Min                |
 * |                  2                        | Eta Max                |
 * |                  3                        | Phi Min                |
 * |                  4                        | Phi Max                |
 * |                  5                        | NonLinearity1          |
 * |                  6                        | NonLinearity2          |
 * |                  7                        | DistanceToBadChannel   |
 * |                  8                        | Timing                 |
 * |                  9                        | TrackMatching          |
 * |                  10                       | ExoticCluster          |
 * |                  11                       | MinEnergy              |
 * |                  12                       | MinNCells              |
 * |                  13                       | MinM02                 |
 * |                  14                       | MaxM02                 |
 * |                  15                       | MinM20                 |
 * |                  16                       | MaxM20                 |
 * |                  17                       | MaximumDispersion      |
 * |                  18                       | NML                    |
*/

class AliCaloPhotonCuts : public AliAnalysisCuts {

  public:
    enum cutIds {
      kClusterType,
      kEtaMin,
      kEtaMax,
      kPhiMin,
      kPhiMax,
      kNonLinearity1,
      kNonLinearity2,
      kDistanceToBadChannel,
      kTiming,
      kTrackMatching,
      kExoticCluster,
      kMinEnergy,
      kNMinCells,
      kMinM02,
      kMaxM02,
      kMinMaxM20,
      kRecConv,
      kDispersion,
      kNLM,
      kNCuts
    };

    enum photonCuts {
      kPhotonIn=0,
      kDetector,
      kAcceptance,
      kClusterQuality,
      kPhotonOut
    };

    enum MCSet {
      // MC data sets
      kNoMC=0,
      // pp 7 TeV 2010
      k14j4,
      // pp 7 TeV 2011
      k14b7,
      k14k1ab,
     // pp 2.76 TeV 2011
      k12f1a,
      k12f1b,
      k12i3,
      kPP2T11P4JJ,
      k15g1b,
      // PbPb 2.76 TeV 2011
      k14a1,
      // pp 8 TeV 2012
      k14e2b,
      kPP8T12P2Pyt8,
      kPP8T12P2Pho,
      kPP8T12P2JJ,
      kPP8T12P2GJLow,
      kPP8T12P2GJHigh,
      kPP8T12P2GJGeant3,
      kPP8T12P2GJGeant4,
      // pPb 5 TeV 2013
      kPPb5T13P2DPMJet,
      kPPb5T13P4JJ,
      kPPb5T13P2HIJAdd,
      kPPb5T13P4JJlow,
      kPPb5T13P4JJhigh,
      k16c3a,
      k16c3b,
      k16c3c,
      kPPb5T13P4DPMJet,
      // pp 2.76TeV 2013
      k15g2,
      kPP2T13P1JJ,
      k15a3b,
      // pp 13 TeV 2015
      kPP13T15P2Pyt8,
      kPP13T15P2EPOS,
      k15k5,
      // pp 5 TeV 2015
      k16h8a,
      k16h8b,
      k16k3a,
      k16k5a,
      k16k5b,
      k17e2,
      k18j3,
      k16h3,
      // pp 5 TeV 2017
      k17l4b,
      k17l3b,
      k18j2,
      k18b8,
      k18b10,
      k18l2,
      // PbPb 5 TeV 2015
      kPbPb5T15HIJING,
      k16k3b,
      // PbPb 5 TeV 2018
      kPbPb5T18HIJING,
      // pp 13 TeV 2016
      kPP13T16P1Pyt8,
      kPP13T16P1Pyt8LowB,
      kPP13T16P1EPOS,
      kPP13T16P1JJ,
      kPP13T16P1JJLowB,
      kPP13T16P1JJTrigger,
      k17h8a,
      k17h8b,
      k17h8c,
      k17c3b1,
      k17c3a1,
      k17c3b2,
      k17c3a2,
      k17i3a1,
      // pPb 5 TeV 2016
      kPPb5T16EPOS,
      kPPb5T16DPMJet,
      k17g8a,
      k17d2a,
      k17d2b,
      // pPb 8 TeV 2016
      k17f3a,
      k17f3b,
      k17f4a,
      k17f4b,
      k18f3bc,
      k17g6b2a,
      k17g6b2b,
      k17g6b3a,
      k17g6b3b,
      k17g8b,
      k17g8c,
      k17g6b1a,
      k18b9b,
      k18b9c,
      // pp 13 TeV 2017
      k17k1,
      kPP13T17P1Pyt8,
      kPP13T17P1Pho,
      kPP13T17P1Pyt6,
      kPP13T17P1Pyt8Str,
      kPP13T17P1Pyt8LowB,
      kPP13T17P1JJ,
      kPP13T17P1JJTrigger,
      kLHC21j8a,
      kPP13T17HERJJ,
      // pp 13 TeV 2018
      kPP13T18P1JJ,
      kPP13T18P1JJTrigger,
      // Xe-Xe MC
      kXeXe5T17HIJING,
      //
      kPP13T18P1Pyt8,
      kPP13T18P1Pyt8LowB,
      // Pb-Pb 5 TeV 2015 Gamma-Jet
      kLHC18b11c,

      // pPb 8 TeV pass2 productions
      kLHC20f11d, // MB general purpose
      kLHC21d2a, // GJ EMC+DCAL
      kLHC21d2b, // JJLow EMC+DCAL
      kLHC21d2c, // JJHigh EMC+DCAL
      // pp 13TeV eta prime biased JJ MC - test production, 20 pT-hard bins
      kLHC23a4b,
      // pp 13TeV eta prime biased JJ MC - 7 pT-hard bins
      kLHC24j3,
      // Data starts here
      k10pp7TeV,
      k10pp900GeV,
      k10PbPb2760GeV,
      k11pp2760GeV,
      k11pp7TeV,
      k11PbPb2760GeV,
      k12pp8TeV,
      k13pPb5023GeV,
      k13pp2760GeV,
      k15pp13TeV,
      k15pp5TeV,
      k15PbPb5TeV,
      k16pp13TeV,
      k16pp13TeVLow,
      k16pPb5023GeV,
      k16pPb8TeV,
      k17pp13TeV,
      k17pp13TeVLow,
      k17pp13TeVNo,
      k17XeXe5440GeV,
      k17pp5TeV,
      k18pp13TeV,
      k18pp13TeVLow,
      k18PbPb5TeV
    };


    //handeling of CutString
    static const char * fgkCutNames[kNCuts];
    Bool_t      SetCutIds(TString cutString);
    Int_t       fCuts[kNCuts];
    Bool_t      SetCut(cutIds cutID, Int_t cut);
    Bool_t      UpdateCutString();
    void        PrintCuts();
    void        PrintCutsWithValues(const TString analysisCutSelection);

    Bool_t      InitializeCutsFromCutString(const TString analysisCutSelection);
    TString     GetCutNumber();
    Int_t       GetClusterType() {return fClusterType;}
    Int_t       GetMinNLMCut() {return fMinNLM;}
    Int_t       GetMaxNLMCut() {return fMaxNLM;}
    Bool_t      IsNLMCutUsed() {return fUseNLM;}

    //Constructors
    AliCaloPhotonCuts(Int_t isJetJet=0, const char *name="ClusterCuts", const char * title="Cluster Cuts");
    AliCaloPhotonCuts(const AliCaloPhotonCuts&);
    AliCaloPhotonCuts& operator=(const AliCaloPhotonCuts&);

    //virtual destructor
    virtual     ~AliCaloPhotonCuts();

    virtual Bool_t  IsSelected(TObject* /*obj*/)               {return kTRUE;}
    virtual Bool_t  IsSelected(TList* /*list*/)                {return kTRUE;}

    Bool_t      ClusterIsSelected(AliVCluster* cluster, AliVEvent *event, AliMCEvent *mcEvent,Int_t isMC, Double_t weight=1., Long_t clusterID = -1);
    Bool_t      ClusterIsSelectedBeforeTrackMatch(){return fIsCurrentClusterAcceptedBeforeTM;}
    Bool_t      ClusterIsSelectedMC(AliMCParticle *particle,AliMCEvent *mcEvent);
    Bool_t      ClusterIsSelectedElecMC(AliMCParticle *particle,AliMCEvent *mcEvent);
    Bool_t      ClusterIsSelectedElecAODMC(AliAODMCParticle *particle,TClonesArray *aodmcArray);
    Bool_t      ClusterIsSelectedAODMC(AliAODMCParticle *particle,TClonesArray *aodmcArray);
    Bool_t      ClusterIsIsolated(Int_t clusterID, AliAODConversionPhoton *PhotonCandidate);
    Bool_t      ParticleIsSelectedAODMC(AliVEvent* event, AliAODMCParticle *particle, bool selectStable = true);

    void        SetLightOutput( Int_t flag )                  {fDoLightOutput = flag; return;}
    //correct NonLinearity
    void        SetV0ReaderName(TString name)                  {fV0ReaderName = name; return;}
    void        SetCaloTrackMatcherName(TString name)          {fCaloTrackMatcherName = name; return;}
    void        SetCaloTriggerMimicHelperName(TString name)    {fCaloTriggerMimicHelperName = name; return;}
    void        SetCaloIsolationName(TString name)             {fCaloIsolationName = name; return;}
    MCSet       FindEnumForMCSet(TString namePeriod);

    void        ApplyNonLinearity(AliVCluster* cluster, Int_t isMC, AliVEvent *event = 0x0);
    void        ApplySMWiseEnergyCorrection(AliVCluster* cluster, Int_t isMC, AliVEvent *event = 0x0);

    Float_t     FunctionNL_kPi0MC(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5, Float_t p6);
    Float_t     FunctionNL_PHOSOnlyMC(Float_t e, Float_t p0, Float_t p1, Float_t p2);

    Float_t     FunctionNL_kSDM(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3 = 1.0);
    Float_t     FunctionNL_DPOW(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5);
    Float_t     FunctionNL_SPOW(Float_t e, Float_t p0, Float_t p1, Float_t p2);
    Float_t     FunctionNL_DExp(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5, Float_t p6 = 1.0, Float_t p7 = 1.0);
    Float_t     FunctionNL_SExp(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3 = 1.0);
    Float_t     FunctionNL_ExpExp(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3);
    Float_t     FunctionNL_LinLogConst(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t const1, Float_t const2);
    //predefined functions
    Float_t     FunctionNL_kPi0MCv1(Float_t e);
    Float_t     FunctionNL_kPi0MCv2(Float_t e);
    Float_t     FunctionNL_kPi0MCv3(Float_t e);
    Float_t     FunctionNL_kPi0MCv5(Float_t e);
    Float_t     FunctionNL_kPi0MCv6(Float_t e);
    Float_t     FunctionNL_kPi0MCMod(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5, Float_t p6);
    Float_t     FunctionNL_OfficialTB_50MeV_Data(Float_t e);
    Float_t     FunctionNL_OfficialTB_100MeV_Data(Float_t e);
    Float_t     FunctionNL_OfficialTB_100MeV_Data_Sys0(Float_t e);
    Float_t     FunctionNL_OfficialTB_100MeV_Data_Sys1(Float_t e);
    Float_t     FunctionNL_OfficialTB_100MeV_Data_Sys2(Float_t e);
    Float_t     FunctionNL_OfficialTB_100MeV_Data_Sys3(Float_t e);
    Float_t     FunctionNL_OfficialTB_100MeV_Data_V2(Float_t e);
    Float_t     FunctionNL_OfficialTB_100MeV_Data_V2_NoScale(Float_t e);
    Float_t     FunctionNL_OfficialTB_150MeV_Data(Float_t e);
    Float_t     FunctionNL_OfficialTB_300MeV_Data(Float_t e);
    Float_t     FunctionNL_OfficialTB_50MeV_MC(Float_t e);
    Float_t     FunctionNL_OfficialTB_100MeV_MC(Float_t e);
    Float_t     FunctionNL_OfficialTB_100MeV_MC_V2(Float_t e);
    Float_t     FunctionNL_OfficialTB_150MeV_MC(Float_t e);
    Float_t     FunctionNL_OfficialTB_300MeV_MC(Float_t e);
    Float_t     FunctionNL_kSDMv5(Float_t e);
    Float_t     FunctionNL_kSDMv6(Float_t e);
    Float_t     FunctionNL_kTestBeamCorrected(Float_t e);
    Float_t     FunctionNL_kTestBeamv2(Float_t e);
    Float_t     FunctionNL_kTestBeamv3(Float_t e);
    Float_t     FunctionNL_kTestBeamv4(Float_t e);
    Float_t     FunctionNL_kTestBeamMod(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5, Float_t p6);

    void        InitCutHistograms(TString name="");
    void        SetFillCutHistograms(TString name="")           {if(!fHistograms){InitCutHistograms(name);} return;}
    TList*      GetCutHistograms()                              {return fHistograms;}
    TList*      GetExtQAHistograms()                            {return fHistExtQA;}
    void        FillClusterCutIndex(Int_t photoncut)            {if(fHistCutIndex)fHistCutIndex->Fill(photoncut); return;}
    void        InitializeEMCAL(AliVEvent *event);
    void        InitializePHOS(AliVEvent *event);

    void        SetExtendedMatchAndQA(Int_t extendedMatchAndQA) {fExtendedMatchAndQA = extendedMatchAndQA; return;}
    void        FillHistogramsExtendedQA(AliVEvent *event, Int_t isMC);
    Double_t    GetTotalEnergyDeposit(AliVEvent *event);
    void        SetIsPureCaloCut(Int_t merged)                  {fIsPureCalo = merged; return;}
    Int_t       GetIsPureCaloCut()                              {return fIsPureCalo;}
    Int_t       GetNactiveEmcalCells()                          {return fNactiveEmcalCells;}
    Int_t       GetIsConversionRecovery()                       {return fUseRecConv;}
    Float_t     GetInvMassConversionRecovery()                  {return fMaxMGGRecConv;}
    Bool_t      GetReduceTriggeredPhiDueBadDDLs()               {return fReduceTriggeredPhiDueBadDDLs;}
    void        FillEtaPhiMapForClusterInBg(const Double_t eta, const Double_t phi, const Double_t weight){  if(fHistClusterEtavsPhiBackground && !fDoLightOutput && fUseEtaPhiMapForBackCand) fHistClusterEtavsPhiBackground->Fill(phi, eta, weight);}
    void        SetUseEtaPhiMapForBackCand(const Bool_t tmp)    {fUseEtaPhiMapForBackCand = tmp;}

    int         GetParticleIndex(int pdgcode);
    void        FillHistEnergyFracPartInClus(AliVEvent* event, AliAODCaloCluster* clus, double weight);
    void        FillHistParticleAbundanceMC(AliVEvent* event, AliAODMCParticle* part, double weight);

    // Cut functions
    Bool_t      AcceptanceCuts(AliVCluster* cluster, AliVEvent *event, Double_t weight);
    Bool_t      ClusterQualityCuts(AliVCluster* cluster,AliVEvent *event, AliMCEvent *mcEvent, Int_t isMC, Double_t weight, Long_t clusterID);

    Bool_t      MatchConvPhotonToCluster(AliAODConversionPhoton* convPhoton, AliVCluster* cluster, AliVEvent* event, Double_t weight=1.);
    void        MatchTracksToClusters(AliVEvent* event, Double_t weight=1., Bool_t isEMCalOnly = kTRUE, AliMCEvent *mcEvent = 0x0);
    void        MatchElectronTracksToClusters(AliVEvent* event, AliMCEvent* MCevent, vector<AliVCluster*> vCluster, Int_t isMC, vector<Int_t> vElectronTracks = {}, Double_t weight = 1.);
    Bool_t      CheckClusterForTrackMatch(AliVCluster* cluster);
    Int_t       GetNumberOfLocalMaxima(AliVCluster* cluster, AliVEvent * event);
    Int_t       GetNumberOfLocalMaxima(AliVCluster* cluster, AliVEvent * event,  Int_t *absCellIdList, Float_t* maxEList);
    Bool_t      AreNeighbours(Int_t absCellId1, Int_t absCellId2);
    Int_t       GetModuleNumberAndCellPosition(Int_t absCellId, Int_t & icol, Int_t & irow);
    void        SplitEnergy(Int_t absCellId1, Int_t absCellId2, AliVCluster* cluster, AliVEvent* event,
                            Int_t isMC, AliAODCaloCluster* cluster1, AliAODCaloCluster* cluster2);
    Int_t       FindLargestCellInCluster(AliVCluster* cluster, AliVEvent* event);
    Int_t       FindSecondLargestCellInCluster(AliVCluster* cluster, AliVEvent* event);
    Bool_t      CheckDistanceToBadChannel(AliVCluster* cluster, AliVEvent* event);
    Bool_t      CheckDistanceToBadChannelSwapping(const Int_t CellID, AliVEvent* event, const Int_t DistanceToBorder = 0);
    Int_t       ClassifyClusterForTMEffi(AliVCluster* cluster, AliVEvent* event, AliMCEvent* mcEvent, Bool_t isESD);

    std::vector<Int_t> GetVectorMatchedTracksToCluster(AliVEvent* event, AliVCluster* cluster);
    std::vector<Int_t> GetVectorMatchedSecTracksToCluster(AliVEvent* event, AliVCluster* cluster);
    Bool_t      GetClosestMatchedTrackToCluster(AliVEvent* event, AliVCluster* cluster, Int_t &trackLabel);
    Bool_t      GetHighestPtMatchedTrackToCluster(AliVEvent* event, AliVCluster* cluster, Int_t &trackLabel);
    Bool_t      IsClusterPi0(AliVEvent *event, AliMCEvent *mcEvent, AliVCluster *cluster);
    Bool_t      CheckForReconstructedConversionPairs(vector<AliAODConversionPhoton*> &vecPhotons, vector<Int_t> &vecReject);
    std::pair<double, double>    GetAngleForConvReject(const TLorentzVector &photon1, const TLorentzVector &photon2);
    std::pair<double, double>    GetAngleForConvReject(const AliAODConversionPhoton *photon1, const AliAODConversionPhoton *photon2);
    Bool_t      CheckVectorForIndexAndAdd(vector<Int_t> &vec, Int_t tobechecked, Bool_t addIndex );
    void        CleanClusterLabels(AliVCluster* clus,AliMCEvent *mcEvent);
    void        CleanClusterLabels(AliVCluster* clus,TClonesArray *aodTrackArray);
    AliCaloTrackMatcher* GetCaloTrackMatcherInstance()          {return fCaloTrackMatcher;}
    AliPhotonIsolation* GetPhotonIsolationInstance()        { return fCaloIsolation; }

    Bool_t      GetIsAcceptedForBasicCounting()                 {return fIsAcceptedForBasic;}

    Bool_t      GetDoFlatEnergySubtraction()                    {return fDoFlatEnergySubtraction;}
    Bool_t      GetDoSecondaryTrackMatching()                   {return fDoSecondaryTrackMatching;}

    Bool_t      SetPoissonParamCentFunction(int isMC);
    Bool_t      SetNMatchedTracksFunc(float meanCent);
    Double_t    Sigmoid(Double_t x, Double_t a, Double_t b, Double_t c);
    Double_t    NOCFunction(Double_t x, Double_t start_y, Double_t a_drop, Double_t b_drop, Double_t c_drop, Double_t c_rise, Double_t transition_width);        // Function which describes the Neutral Overlap Correction factor
    Double_t    CorrectEnergyForOverlap(float meanCent, float E = 0);
    Int_t       GetDoEnergyCorrectionForOverlap()               {return fDoEnergyCorrectionForOverlap;}
    Double_t    GetMeanEForOverlap(Double_t cent, Double_t* par);
    Double_t    GetNOCParameter0(Double_t cent, Double_t* par);
    Double_t    GetNOCParameter1(Double_t cent, Double_t* par);
    Double_t    GetNOCParameter2(Double_t cent, Double_t* par);
    Double_t    GetNOCParameter3(Double_t cent, Double_t* par);
    Double_t    GetNOCParameter4(Double_t cent, Double_t* par);

    // modify acceptance via histogram with cellID
    void        SetHistoToModifyAcceptance(TH1S* histAcc)       {fHistoModifyAcc  = histAcc; return;}

    // Set basic merging cuts
    void        SetSeedEnergy(Double_t seed)                    {fSeedEnergy      = seed; return;}
    void        SetLocMaxCutEDiff(Double_t diffCut)             {fLocMaxCutEDiff  = diffCut; return;}

    //Set Electron Cluster calibration
    void        SetElectronClusterCalibration(Bool_t calib)     {fUseElectronClusterCalibration = calib; return;}
    Bool_t      GetElectronClusterCalibration()                 {return fUseElectronClusterCalibration;}

    // GetCellIds for EMCal PHOS
    Int_t       GetCaloCellIdFromEtaPhi(const Double_t eta, const Double_t phi);

    // Set Individual Cuts
    Bool_t      SetClusterTypeCut(Int_t);
    Bool_t      SetMinEtaCut(Int_t);
    Bool_t      SetMaxEtaCut(Int_t);
    Bool_t      SetMinPhiCut(Int_t);
    Bool_t      SetMaxPhiCut(Int_t);
    Bool_t      SetDistanceToBadChannelCut(Int_t);
    Bool_t      SetTimingCut(Int_t);
    Bool_t      SetTrackMatchingCut(Int_t);
    Bool_t      SetExoticClusterCut(Int_t);
    Bool_t      SetMinEnergyCut(Int_t);
    Bool_t      SetMinNCellsCut(Int_t);
    Bool_t      SetMaxM02(Int_t);
    Bool_t      SetMinM02(Int_t);
    Bool_t      SetMinMaxM20(Int_t);
    Bool_t      SetRecConv(Int_t);
    Bool_t      SetDispersion(Int_t);
    Bool_t      SetNLM(Int_t);
    Bool_t      SetNonLinearity1(Int_t);
    Bool_t      SetNonLinearity2(Int_t);


    Int_t       GetNonLinearity()                                {return fSwitchNonLinearity;}
    void        SetUseNonLinearitySwitch( Bool_t useNonLin)      { fUseNonLinearity = useNonLin;}

    Float_t     FunctionM02 (Float_t E, Float_t a, Float_t b, Float_t c, Float_t d, Float_t e);
    Float_t     CalculateMaxM02 (Int_t maxM02, Float_t clusEnergy);
    Float_t     CalculateMinM02 (Int_t minM02, Float_t clusEnergy);
    Double_t    GetDistanceBetweenClusters(AliVCluster* cluster1, AliVCluster* cluster2);
    void        SetLogBinningXTH1 (TH1* histoRebin);
    void        SetLogBinningXTH2 (TH2* histoRebin);
    void        SetLogBinningYTH2 (TH2* histoRebin);

    Bool_t      IsExoticCluster ( AliVCluster *cluster, AliVEvent *event, Float_t& energyStar );
    Float_t     GetECross ( Int_t absID, AliVCaloCells* cells );
    Bool_t      IsCellNextToCluster( Int_t absID, Double_t Ethresh, AliVCaloCells* cells, const std::vector<int> &vCellID = {} );
    Bool_t      AcceptCellByBadChannelMap (Int_t absID );
    Bool_t      IsAbsIDsFromTCard(Int_t absId1, Int_t absId2) const;
    void        SetExoticsMinCellEnergyCut(Double_t minE)       { fExoticMinEnergyCell = minE; return;}
    void        SetExoticsQA(Bool_t enable)                     { fDoExoticsQA         = enable; return;}

    Float_t     GetMinClusterEnergy()                           { return fMinEnergy;}

    Double_t    GetMinEtaCut() {return fMinEtaCut;}
    Double_t    GetMaxEtaCut() {return fMaxEtaCut;}
    Double_t    GetMinPhiCut() {return fMinPhiCut;}
    Double_t    GetMaxPhiCut() {return fMaxPhiCut;}
    Double_t    GetMinPhiCutDMC() {return fMinPhiCutDMC;}
    Double_t    GetMaxPhiCutDMC() {return fMaxPhiCutDMC;}
    // Function to set correction task setting
    void SetCorrectionTaskSetting(TString setting) {fCorrTaskSetting = setting;}

    AliEMCALGeometry* GetGeomEMCAL(){return fGeomEMCAL;}
    AliPHOSGeometry*  GetGeomPHOS() {return fGeomPHOS;}


  protected:
    TList      *fHistograms;
    TList      *fHistExtQA;

    AliCaloTrackMatcher* fCaloTrackMatcher;             // pointer to CaloTrackMatcher
    AliCaloTriggerMimicHelper* fCaloTriggerMimicHelper; //Array wich points to AliCaloTriggerMimicHelper for each Event Cut
    AliPhotonIsolation*  fCaloIsolation;                // pointer to PhotonIsolation
    AliEMCALGeometry*   fGeomEMCAL;                     // pointer to EMCAL geometry
    AliEMCALRecoUtils*  fEMCALRecUtils;                 // pointer to EMCAL recUtils
    Bool_t     fEMCALInitialized;                       // flag for EMCal initialization
    AliPHOSGeometry*    fGeomPHOS;                      // pointer to PHOS geometry
    Bool_t     fPHOSInitialized;                        // flag for PHOS initialization
    Int_t      fPHOSCurrentRun;                         // PHOS: current processed run for bad channel map
    TObjArray* fEMCALBadChannelsMap;                    // pointer to EMCAL bad channel map
    TH1C*      fEMCALBadChannelsMap1D;                  // pointer to EMCAL bad channel map (1D)
    TH2I**     fPHOSBadChannelsMap;                     // pointer to PHOS bad channel map
    TProfile*  fBadChannels;                            // TProfile with bad channels
    Int_t      fNMaxEMCalModules;                       // max number of EMCal Modules
    Int_t      fNMaxPHOSModules;                        // max number of PHOS Modules
    TH1S*      fHistoModifyAcc;                         // hisogram for modified acceptance, if leadCellID->1 accept cluster, if leadCellID->0 reject cluster
    TClonesArray*   fAODMCTrackArray;                   // pointer to track array

    Int_t    fDoLightOutput;                           // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
    Int_t     fIsMC;                                    // Flag for usage of JetJet MC

    Bool_t     fIsCurrentClusterAcceptedBeforeTM;       // flag if latest checked cluster would have been accepted before track matching cut

    Bool_t     fUseEtaPhiMapForBackCand;                // flag to switch on Eta Phi Map for cluster candidates when doinf background calculation (for swapping method)

    //for NonLinearity correction
    TString   fV0ReaderName;                            // Name of V0Reader
    TString   fCorrTaskSetting;                         // Name of Correction Task Setting
    TString   fCaloTrackMatcherName;                    // Name of global TrackMatching instance
    TString   fCaloTriggerMimicHelperName;              // Name of global TriggerMimicHelper instance
    TString   fCaloIsolationName;                       // Name of global Isolation instance
    TString   fPeriodName;                              // PeriodName of MC
    MCSet     fCurrentMC;                               // enum for current MC set being processed

    //cuts
    Int_t     fClusterType;                             // which cluster do we have
    Float_t   fIsolationRadius;                         // radius of isolation cone
    Float_t   fMomPercentage;                           // percentage of the isolated Photon Pt
    Bool_t    fUsePhotonIsolation;                      // is photonisolation turned on
    Double_t  fMinEtaCut;                               // min eta cut
    Double_t  fMinEtaInnerEdge;                         // min eta of inner Edge (DCal)
    Double_t  fMaxEtaCut;                               // max eta cut
    Double_t  fMaxEtaInnerEdge;                         // max eta of inner Edge (DCal)
    Bool_t    fUseEtaCut;                               // flag for switching on eta cut
    Double_t  fMinPhiCut;                               // phi cut
    Double_t  fMaxPhiCut;                               // phi cut
    Double_t  fMinPhiCutDMC;                            // phi cut
    Double_t  fMaxPhiCutDMC;                            // phi cut
    Bool_t    fUsePhiCut;                               // flag for switching on phi cut
    Bool_t    fReduceTriggeredPhiDueBadDDLs;            // Throw out DDLs in Trigger which might be bad
    Double_t  fMinDistanceToBadChannel;                 // minimum distance to bad channel
    Int_t     fUseDistanceToBadChannel;                 // flag for switching on distance to bad channel cut: 0 off, 1 on without corners, 2 on with corners included
    Double_t  fMaxTimeDiff;                             // maximum time difference to triggered collision
    Double_t  fMinTimeDiff;                             // minimum time difference to triggered collision
    Double_t  fMaxTimeDiffHighPt;                       // maximum time difference to triggered collision at high Pt
    Double_t  fMinTimeDiffHighPt;                       // minimum time difference to triggered collision at high Pt
    Double_t  fOffsetTimeMC;                            // time offset difference to triggered collision in MC, seen in PbPb 5.02 TeV production
    Bool_t    fUseTimeDiff;                             // flag for switching on time difference cut
    Bool_t    fUseTimeDiffMC;                           // flag for switching on time difference cut for MC, needed in PbPb 5.02 TeV production
    Double_t  fMaxDistTrackToClusterEta;                // minimum distance between track and cluster in eta
    Double_t  fMinDistTrackToClusterPhi;                // minimum distance between track and cluster in phi
    Double_t  fMaxDistTrackToClusterPhi;                // maximum distance between track and cluster in phi
    Bool_t    fUseDistTrackToCluster;                   // flag for switching on distance between track and cluster cut
    Int_t     fUsePtDepTrackToCluster;                  // flag for switching on pT dependent matching parameters
    Bool_t    fUseOnlyMatchedClusters;                  // flag for switching on to use only the matched clusters (only for track matching studies)
    Bool_t    fTriggerMimicHelper_found;                //Flag if TriggerMimicHelper has been found in this instance
    TF1*      fFuncPtDepEta;                            // TF1 for pT dep cutting in eta
    TF1*      fFuncPtDepPhi;                            // TF1 for pT dep cutting in phi
    TRandom3  fRandom;                                  // random for effi generation
    Int_t     fUseTimingEfficiencyMCSimCluster;         // flag for switching on TimingEfficiencyMCSimCluster
    TF1*      fFuncTimingEfficiencyMCSimCluster;        // TF1 for TimingEfficiencyMCSimCluster
    TF1*      fFuncTimingEfficiencyMCSimClusterHighPt;  // TF1 for fFuncTimingEfficiencyMCSimClusterHighPt
    Float_t   fTimingEfficiencyMCSimClusterLowPtEnd;    // fFuncTimingEfficiencyMCSimCluster is used until this Threshold; above 1 is used until fFuncTimingEfficiencyMCSimClusterHighPt is used
    Float_t   fTimingEfficiencyMCSimClusterHighPtStart; // fFuncTimingEfficiencyMCSimClusterHighPt is used above this Threshold
    TF1*      fFuncNCellCutEfficiencyEMCal;             // TF1 for NCell cut efficiency applied on EMCal MC clusters with NCell<2
    Float_t   fMinTMDistSigma;                          // number of sigma's for TM using PHOS
    Bool_t    fUseEOverPVetoTM;                         // flag for switching on E/P veto (forbidding tracks to match clusters if clusterE/trackP > someValue
    Double_t  fEOverPMax;                               // maximum value for E/P of a track to be considered for TM
    Bool_t    fUseTMMIPsubtraction;                     // flag for switching on MIP subtraction
    Bool_t    fUseElectronClusterCalibration;           // flag for switching on electron cluster calibration
    Int_t     fExtendedMatchAndQA;                      // switching on ext matching histograms (1) / ext QA_noCell (2) / ext matching + ext QA_noCell (3) / extQA + cell (4) / ext match + extQA + cell (5) or all off (0)
    Double_t  fExoticEnergyFracCluster;                 // exotic energy compared to E_cross cluster cut
    Double_t  fExoticMinEnergyTCard;                    // min exotic energy to check for all cells in the same Tcard
    Double_t  fExoticMinEnergyCell;                     // minimum energy of cell to test for exotics
    Int_t     fUseExoticCluster;                        // flag for switching on exotic cluster cut
    Bool_t    fDoExoticsQA;                             // flag for switching on exotic cluster cut
    Double_t  fMinEnergy;                               // minium energy per cluster
    Bool_t    fDoFlatEnergySubtraction;                 // enable flat energy subtraction
    Double_t  fSeedEnergy;                              // seed energy for clusters
    Double_t  fLocMaxCutEDiff;                          // cut on energy difference between two cells
    Bool_t    fUseMinEnergy;                            // flag for switching on minimum energy cut
    Int_t     fMinNCells;                               // minimum number of cells
    Int_t     fMaxNCells;                               // maximum number of cells
    Double_t  fMinENCell;                               // minimum energy for number of cells cut
    Int_t     fUseNCells;                               // flag for switching on minimum N Cells cut
    Double_t  fMaxM02;                                  // maximum M02
    Double_t  fMinM02;                                  // minimum M02
    Int_t     fUseM02;                                  // flag for switching on M02 cut
    Int_t     fMaxM02CutNr;                             // maximum M02 CutNr
    Int_t     fMinM02CutNr;                             // minimum M02 CutNr
    Double_t  fMaxM20;                                  // maximum M20
    Double_t  fMinM20;                                  // minimum M20
    Bool_t    fUseM20;                                  // flag for switching on M20 cut
    Float_t   fMaxMGGRecConv;                           // maximum invariant mass below which the 2 clusters are gonna be combined assuming they are from a photon
    Float_t   fConvRejMinAngle;                         // minimum angle with respect to z that two clusters need in order to not be rejected as conversions (conversion should only be bend by magnetic field and therefore their position in z should be the same)
    Float_t   fConvRejMaxOpenAngle;                     // maximum opening angle of a cluster pair to be considerd to be a conversion from the same photon
    Int_t     fUseRecConv;                              // flag to switch on conversion recovery
    Double_t  fMaxDispersion;                           // maximum dispersion
    int       fUseDispersion;                           // flag for switching on dispersion cut. If set to two, this acts as the cluster efficiency (cut is overloaded from cutID=10 onwards)
    Int_t     fMinNLM;                                  // minimum number of local maxima in cluster
    Int_t     fMaxNLM;                                  // maximum number of local maxima in cluster
    Bool_t    fUseNLM;                                  // flag for switching on NLM cut
    Int_t     fNonLinearity1;                           // selection of nonlinearity correction, part1
    Int_t     fNonLinearity2;                           // selection of nonlinearity correction, part2
    Int_t     fSwitchNonLinearity;                      // selection (combined) of NonLinearity
    Bool_t    fUseNonLinearity;                         // flag for switching NonLinearity correction
    Int_t     fIsPureCalo;                              // flag for MergedCluster analysis
    Int_t     fNactiveEmcalCells;                       // total number of active emcal cells
    Bool_t    fDoSecondaryTrackMatching;                // flag to switch on secondary trackmatching
    Int_t     fDoEnergyCorrectionForOverlap;            // mask to switch on a special for PbPb developed cluster energy correction, 0 = off, 1 = on with mean, 2 = on with random values
    TF1*      fFuncPoissonParamCent;                    // TF1 to describe the poisson parameter that you get from fitting a poisson dsitribution to the number of matched tracks per cluster as function of centrality
    TF1*      fFuncNMatchedTracks;                      // TF1 poisson distribution to describe the number of matched tracks per cluster for a specific centrality
    Double_t  fParamMeanTrackPt[3];                     // TF1 distribution to describe the mean pT of tracks as function of centrality. Half of this value is used as neutral energy overlap correction.
    Double_t  fNOCParam0[2];                            // parameter for function [0]*exp(x*[1]+2) to describe 1st parameter as function of cent of the 5 params for the NOC
    Double_t  fNOCParam1[3];                            // parameter for function [0]/(1+exp([1]*(x-[2]))) to describe 2nd parameter as function of cent of the 5 params for the NOC
    Double_t  fNOCParam2[4];                            // parameter for function pol3 to describe 3rd parameter as function of cent of the 5 params for the NOC
    Double_t  fNOCParam3[3];                            // parameter for function pol2 to describe 4th parameter as function of cent of the 5 params for the NOC
    Double_t  fNOCParam4[3];                            // parameter for function pol2 to describe 5th parameter as function of cent of the 5 params for the NOC
    Float_t   fMeanNMatchedTracks;                      // Mean number of matched primary tracks, stored to reduce CPU time for neutral overlap correction
    TF1*      fFuncNOCMaxBoltz;                         // TF1 Maxwell Boltzmann to describe the neutral overlap correction for a specific centrality
    // for cluster efficiency application in data
    bool      fApplyClusterEffOnData;                   // switch weather cluster efficiency should be applied on data or MC
    bool      fApplyClusterEffToEMProbesOnly;           // if true and we are in the MC case, only clusters with a leading photon or electron contribution are considered
    TF1*      fClusterEfficiencyFunc;                   // function giveing the probability a cluster should survive the cuts as function of energy
    // reject cluster next to higher energetic cluster
    double fMinEnergyClusNeighbor;                      // value of min energy that the neighboring cluster has to have


    //vector
    std::vector<Int_t> fVectorMatchedClusterIDs;        // vector with cluster IDs that have been matched to tracks in merged cluster analysis

    // CutString
    TObjString* fCutString;                             // cut number used for analysis
    TString     fCutStringRead;

    // Histograms
    TH1F*     fHistCutIndex;                            // bookkeeping for cuts
    TH1F*     fHistAcceptanceCuts;                      // bookkeeping for acceptance cuts
    TH2F*     fHistClusterIdentificationCuts;           // bookkeeping for cluster identification cuts

    TH2F*     fHistClusterEtavsPhiBeforeAcc;            // eta-phi-distribution before acceptance cuts
    TH2F*     fHistClusterEtavsPhiAfterAcc;             // eta-phi-distribution of all after acceptance cuts
    TH2F*     fHistClusterEtavsPhiAfterQA;              // eta-phi-distribution of all after cluster quality cuts
    TH2F*     fHistClusterEtavsPhiAfterQA_onlyTriggered;// eta-phi-distribution of all after cluster quality cuts
    TH2F*     fHistClusterEtavsPhiBackground;           // eta-phi-distribution of all clusters in background calculation
    TH2F*     fHistClusterEtavsPhiMatchedClusters;      // eta-phi-distribution of all track-matched clusters
    TH2F*     fHistClusterTimevsEBeforeQA;              // Cluster time vs E before cluster quality cuts
    TH2F*     fHistClusterTimevsEAfterQA;               // Cluster time vs E after cluster quality cuts
    TH2F*     fHistClusterTimevsELowGain;               // Cluster time vs E for low gain cluster
    TH2F*     fHistClusterTimevsEHighGain;              // Cluster time vs E for high gain cluster
    TH2F*     fHistClusterTimevsEHighGainAllCells;     // Cluster time vs E for high gain cluster (for all cells in the cluster)
    TH1F*     fHistEnergyOfClusterBeforeNL;             // enery per cluster before NonLinearity correction
    TH1F*     fHistEnergyOfClusterAfterNL;              // enery per cluster after NonLinearity correction
    TH1F*     fHistEnergyOfClusterBeforeQA;             // enery per cluster before acceptance cuts
    TH1F*     fHistEnergyOfClusterAfterQA;              // enery per cluster after cluster quality cuts
    TH1F*     fHistNCellsBeforeQA;                      // number of cells per cluster before acceptance cuts
    TH1F*     fHistNCellsAfterQA;                       // number of cells per cluster after cluster quality cuts
    TH1F*     fHistM02BeforeQA;                         // M02 before acceptance cuts
    TH1F*     fHistM02AfterQA;                          // M02 after cluster quality cuts
    TH1F*     fHistM20BeforeQA;                         // M20 before acceptance cuts
    TH1F*     fHistM20AfterQA;                          // M20 after cluster quality cuts
    TH1F*     fHistDispersionBeforeQA;                  // dispersion before acceptance cuts
    TH1F*     fHistDispersionAfterQA;                   // dispersion after cluster quality cuts
    TH1F*     fHistNLMBeforeQA;                         // number of local maxima in cluster before acceptance cuts
    TH1F*     fHistNLMAfterQA;                          // number of local maxima in cluster after cluster quality cuts
//     TH2F*     fHistNLMAvsNLMBBeforeQA;                  // number of local maxima in cluster after cluster quality cuts
    TH2F*     fHistNLMVsNCellsAfterQA;                  // number of local maxima vs Ncells in cluster after cluster quality cuts
    TH2F*     fHistNLMVsEAfterQA;                       // number of local maxima vs E in cluster after cluster quality cuts
    //More histograms
    TH2F*     fHistClusterEnergyvsMod;                  // Cluster Energy vs Module Number
    TH2F*     fHistNCellsBigger100MeVvsMod;             // NCells with >0.1 GeV vs Module Number
    TH2F*     fHistNCellsBigger1500MeVvsMod;            // NCells with >1.5 GeV vs Module Number
    TH2F*     fHistEnergyOfModvsMod;                    // Deposited Energy vs Module Number
    TH2F*     fHistClusterEnergyvsNCellsBeforeQA;       // Cluster Energy vs NCells before QA
    TH2F*     fHistClusterEnergyvsNCellsAfterQA;        // Cluster Energy vs NCells after QA
    TH2F*     fHistClusterEnergyvsExoticityBeforeQA;    // Cluster Energy vs Exoticity before QA
    TH2F*     fHistClusterEnergyvsExoticityAfterQA;     // Cluster Energy vs Exoticity before QA
    TH2F*     fHistCellEnergyvsCellID;                  // Cell Energy vs CellID
    TH2F*     fHistClusterEnergyvsCellID;               // Cluster Energy vs CellID
    TH1F*     fHistCellEnergyLG;                        // Cell Energy of low gain cells
    TH1F*     fHistCellEnergyHG;                        // Cell Energy of high gain cells
    TH2F*     fHistCellTimevsCellID;                    // Cell Time vs CellID
    TH2F*     fHistClusterEM02BeforeQA;                 // 2-dim plot E vs. M02
    TH2F*     fHistClusterEM02AfterQA;                  // 2-dim plot E vs. M02
    TH2F*     fHistClusterEM20BeforeQA;                 // 2-dim plot E vs. M20
    TH2F*     fHistClusterEM20AfterQA;                  // 2-dim plot E vs. M20
    TH1F*     fHistClusterIncludedCellsBeforeQA;        // CellIDs in Cluster
    TH1F*     fHistClusterIncludedCellsAfterQA;         // CellIDs in Cluster of accepted ones
    TH1F*     fHistClusterEnergyFracCellsBeforeQA;      // Energy fraction of CellIDs in Cluster
    TH1F*     fHistClusterEnergyFracCellsAfterQA;       // Energy fraction of CellIDs in Cluster of accepted ones
    TH2F*     fHistClusterIncludedCellsTimingAfterQA;   // Timing of CellIDs in Cluster of accepted ones
    TH2F*     fHistClusterIncludedCellsTimingEnergyAfterQA;   // Timing vs Energy of CellIDs in Cluster of accepted ones
    TH2F*     fHistClusterDistanceInTimeCut;            // distance of clusters: within cluster timing cut + within cluster timing cut
    TH2F*     fHistClusterDistanceOutTimeCut;           // distance of clusters: within cluster timing cut + outside cluster timing cut
    TH1F*     fHistClusterDistance1DInTimeCut;          // 1D distance of clusters: within cluster timing cut + within cluster timing cut

    //Track matching histograms
    TH1F*     fHistClusterRBeforeQA;                    // cluster position in R=SQRT(x^2+y^2) (before QA)
    TH1F*     fHistClusterRAfterQA;                     // cluster position in R=SQRT(x^2+y^2) for matched tracks (After QA)
    TH2F*     fHistClusterdEtadPhiBeforeQA;             // 2-dim plot dEta vs. dPhi
    TH2F*     fHistClusterdEtadPhiAfterQA;              // 2-dim plot dEta vs. dPhi for matched tracks (after QA)
    TH1F*     fHistDistanceTrackToClusterBeforeQA;      // distance cluster to track before acceptance cuts
    TH1F*     fHistDistanceTrackToClusterAfterQA;       // distance cluster to track after cluster quality cuts

    //Extended track matching histograms
    TH2F*     fHistClusterdEtadPhiPosTracksBeforeQA;    // 2-dim plot dEta vs. dPhi
    TH2F*     fHistClusterdEtadPhiNegTracksBeforeQA;    // 2-dim plot dEta vs. dPhi
    TH2F*     fHistClusterdEtadPhiPosTracksAfterQA;     // 2-dim plot dEta vs. dPhi
    TH2F*     fHistClusterdEtadPhiNegTracksAfterQA;     // 2-dim plot dEta vs. dPhi
    TH2F*     fHistClusterdEtadPhiPosTracksP_000_075BeforeQA; // 2-dim plot dEta vs. dPhi, positive Tracks, P < 0.75
    TH2F*     fHistClusterdEtadPhiPosTracksP_075_125BeforeQA; // 2-dim plot dEta vs. dPhi, positive Tracks, 0.75 < P < 1.25
    TH2F*     fHistClusterdEtadPhiPosTracksP_125_999BeforeQA; // 2-dim plot dEta vs. dPhi, positive Tracks, P > 1.25
    TH2F*     fHistClusterdEtadPhiNegTracksP_000_075BeforeQA; // 2-dim plot dEta vs. dPhi, negative Tracks, P < 0.75
    TH2F*     fHistClusterdEtadPhiNegTracksP_075_125BeforeQA; // 2-dim plot dEta vs. dPhi, negative Tracks, 0.75 < P < 1.25
    TH2F*     fHistClusterdEtadPhiNegTracksP_125_999BeforeQA; // 2-dim plot dEta vs. dPhi, negative Tracks, P > 1.25
    TH2F*     fHistClusterdEtadPtBeforeQA;              // 2-dim plot dEta vs. Pt
    TH2F*     fHistClusterdEtadPtAfterQA;               // 2-dim plot dEta vs. Pt
    TH2F*     fHistClusterdEtadPtTrueMatched;           // 2-dim plot dEta vs. Pt for validated match
    TH2F*     fHistClusterdPhidPtPosTracksBeforeQA;     // 2-dim plot dPhi vs. Pt, positive Tracks
    TH2F*     fHistClusterdPhidPtNegTracksBeforeQA;     // 2-dim plot dPhi vs. Pt, negative Tracks
    TH2F*     fHistClusterdPhidPtAfterQA;               // 2-dim plot dPhi vs. Pt
    TH2F*     fHistClusterdPhidPtPosTracksTrueMatched;  // 2-dim plot dPhi vs. Pt, positive Tracks, validated match
    TH2F*     fHistClusterdPhidPtNegTracksTrueMatched;  // 2-dim plot dPhi vs. Pt, negative Tracks, validated match
    TH2F*     fHistClusterM02M20BeforeQA;               // 2-dim plot M20 vs. M02
    TH2F*     fHistClusterM02M20AfterQA;                // 2-dim plot M20 vs. M20

    // QA histograms for exotic clusters
    TH2F*     fHistClusterEtavsPhiExotics;              // eta-phi-distribution of exotic clusters
    TH2F*     fHistClusterEM02Exotics;                  // 2-dim plot E vs. M02 for exotic clusters
    TH2F*     fHistClusterEnergyvsNCellsExotics;        // Cluster Energy vs NCells
    TH2F*     fHistClusterEEstarExotics;                // 2-dim plot E vs. Estar for exotic clusters

    // histograms for track matching efficiency
    TH2F*     fHistClusterTMEffiInput;                  //
    TH2F*     fHistClusterTrueElecEtaPhiBeforeTM_30_00; // electron clusters before track matching, only for E_clus > 30 GeV
    TH2F*     fHistClusterTrueElecEtaPhiAfterTM_30_00;  // electron clusters after track matching, only for E_clus > 30 GeV
    TH2F*     fHistClusterEvsTrackECharged;             //
    TH2F*     fHistClusterEvsTrackEChargedLead;         //
    TH2F*     fHistClusterEvsTrackENeutral;             //
    TH2F*     fHistClusterEvsTrackENeutralSubCharged;   //
    TH2F*     fHistClusterEvsTrackEGamma;               //
    TH2F*     fHistClusterEvsTrackEGammaSubCharged;     //
    TH2F*     fHistClusterEvsTrackEConv;                //
    TH2F*     fHistClusterENMatchesNeutral;             //
    TH2F*     fHistClusterENMatchesConv;                //
    TH2F*     fHistClusterENMatchesCharged;             //
    TH2F*     fHistClusterEvsTrackEPrimaryButNoElec;    //
    TH2F*     fHistClusterEvsTrackSumEPrimaryButNoElec; //
    TH1F*     fHistClusterNMatched;                     // Number of matched tracks per cluster

    TH1F*     fHistClusETruePi0_BeforeTM;               // for checking the false positives: how much true pi0s are matched away?
    TH1F*     fHistClusETruePi0_Matched;                //

    // histograms for studying track matching veto with track momentum vs. cluster energy
    TH2F*     fHistMatchedTrackPClusE;                  // track P vs cluster E in case of matching with a cluster
    TH2F*     fHistMatchedTrackPClusEAfterEOverPVeto;   // track P vs cluster E for matched tracks surviving the E/P veto
    TH2F*     fHistMatchedTrackPClusETruePi0Clus;       // track P vs cluster E in case of matching with a true pi0 cluster

    TH1F*     fHistElectronPositronOnEMC;                       // Electron/Positron P for all tracks in EMCal acceptance
    TH1F*     fHistElectronPositronOnEMCCell;                   // Electron/Positron P for all tracks in EMCal acceptance and that hit a good cell
    TH1F*     fHistVerifiedElectronPositronOnEMCCell;           // Electron/Positron P for all tracks in EMCal acceptance and that hit a good cell and are really an electron
    TH2F*     fHistElectronPositronClusterMatch;                // Electron/Positron P vs cluster E in case of matching with a cluster
    TH2F*     fHistElectronPositronClusterMatchSub;             // Electron/Positron P vs E - P in case of matching with a cluster
    TH2F*     fHistElectronPositronClusterMatchEoverP;          // Electron/Positron E/P vs PT of track in case of matching of Electron with cluster
    TH2F*     fHistElectronPositronClusterMatchEoverPonVtx;     // Electron/Positron E/P(at Vtx) vs PT of track in case of matching of Electron with cluster
    TH2F*     fHistElectronPositronClusterMatchEoverPVsE;       // Electron/Positron E/P vs E of cluster in case of matching of Electron with cluster
    TH2F*     fHistElectronPositronClusterMatchEoverPonVtxVsE;  // Electron/Positron E/P(at Vtx) vs E of cluster in case of matching of Electron with cluster
    TH2F*     fHistElectronClusterMatch;                        // Electron P vs cluster E in case of matching with a cluster
    TH2F*     fHistPositronClusterMatch;                        // Positron P vs cluster E in case of matching with a cluster
    TH2F*     fHistTrueElectronPositronClusterMatch;            // True Electron/Positron P vs cluster E in case of matching with a cluster
    TH2F*     fHistTrueNoElectronPositronClusterMatch;          // True No Electron/Positron P vs cluster E in case of matching with a cluster
    TH2F*     fHistElectronClusterMatchTruePID;                 // MC true histogram for purity studies of selected electrons
    TH2F*     fHistElectronClusterMatchTruePIDRecP;             // MC true histogram for purity studies of selected electrons with rec mom on EMC surface
    TH2F*     fHistTrueElectronPositronClusterMatchEoverP;      // True Electron/Positron E/P vs PT of cluster in case of matching of Electron with cluster
    TH2F*     fHistElectronClusterNCellsVsE;                    // Electron/Positron cluster number of cells vs. cluster energy
    TH2F*     fHistElectronClusterMatchEta;                     // distance between track and cluster in eta vs. track momentum
    TH2F*     fHistElectronClusterMatchPhi;                     // distance between track and cluster in phi vs. track momentum
    TH2F*     fHistElectronPositronVsSM;                        // track momentum for electron and positron tracks for each SM seperatly
    TH2F*     fHistElectronPositronVsEta;                       // track momentum for electron and positron tracks vs. eta of cluster
    TH2F*     fHistElectronPositronClusterMatchVsSM;            // track momentum for matched electron and positron tracks for each SM seperatly
    TH2F*     fHistElectronPositronClusterMatchVsEta;           // track momentum for matched electron and positron tracks vs. eta of cluster
    TH1F*     fHistElectronPositronClusterMatchDoubleCount;     // track momentum for matched electron and positron tracks if matched to multiple clusters
    TH2F*     fHistElectronClusterMatchM02;                     // track momentum for matched electron and positron tracks vs M02 of matched cluster
    TH1F*     fHistElectronPositronOnEMCCellMixing;             // Electron/Positron P for all tracks in EMCal acceptance and that hit a good cell
    TH1F*     fHistElectronPositronOnEMCCellMatchedMixing;      // Electron/Positron P for all tracks in EMCal acceptance and that hit a good cell and matched to a cluster
    TH2F*     fHistElectronPositronClusterMatchEoverPMixing;    // Electron/Positron P for all tracks in EMCal acceptance and that hit a good cell and matched to a cluster vs E/p
    std::vector<ElectronMixing> vecElectronMixing;              //! vector holding information on electron tracks for electron mixing 

    // histogram for cluster energy fraction studies
    TH3F*     fHistEnergyFracClusVsPart;                // histogram with particle index, energy fraction of particle in cluster and cluster energy
    TH3F*     fHistDepEnergyClusVsPart;                 // histogram with particle index, total deposited energy of particle and cluster energy
    TH3F*     fHistDepEnergyPartPVsPart;                // histogram with particle index, deposited energy in cluster and MC particle energy
    TH3F*     fHistDepEnergyClusEVsPart;                // histogram with particle index, deposited energy in cluster and cluster energy
    TH2F*     fHistPartPVsPartIndex;                    // histogram with MC gen particle momentum for different particles
    // histogram for conv candidate rejection
    TH2F*     fHistInvMassDiCluster;                    // histogram for monitoring di-cluster mass
    TH2F*     fHistInvMassConvFlagging;                 // histogram for monitoring rejected di-cluster mass
    TH1F*     fHistDiClusterAngle;                      // histogram for monitoring di-cluster angle with respect to z
    Int_t     fNMaxDCalModules;                         // max number of DCal Modules
    Int_t     fgkDCALCols;                              // Number of columns in DCal
    Bool_t    fIsAcceptedForBasic;                      // basic counting

  private:

    ClassDef(AliCaloPhotonCuts,148)
};

#endif
