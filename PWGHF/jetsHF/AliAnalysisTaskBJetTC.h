#ifndef ALIANALYSISTASKBJETTC_H
#define ALIANALYSISTASKBJETTC_H
#include "AliAnalysisTaskEmcalJet.h"
#include "AliV0ReaderV1.h"
#include "AliConvEventCuts.h"
#include "AliAnalysisTaskWeakDecayVertexer.h"
class AliEmcalJet;
class AliAODVertex;
class AliAODTrack;
class TList;
class TH1D;
class TH2D;
class TH3D;
class THnSparse;
class AliHFJetsTagging;
class TClonesArray;
class AliAODMCParticle;
class AliAnalysisUtils;
class TRandom3;
class AliPIDResponse;
class AliHFJetsTaggingVertex;
class AliRDHFJetsCutsVertex;
class AliVertexerTracks;

class AliAnalysisTaskBJetTC : public AliAnalysisTaskEmcalJet
{
public:
    enum TTypeImpPar
    {
        kXY,
        kXYSig,
        kXYZ,
        kXYZSig,
        kXYZSigmaOnly
    };

    AliAnalysisTaskBJetTC();
    AliAnalysisTaskBJetTC(const char *name);
    //AliAnalysisTaskBJetTC(const AliAnalysisTaskBJetTC&);	// copy constructor not implemented yet
    //AliAnalysisTaskBJetTC& operator=(const AliAnalysisTaskBJetTC&);	// assignment operator not implemented yet

    virtual ~AliAnalysisTaskBJetTC();
    virtual void UserCreateOutputObjects();
    virtual void Terminate(Option_t *);
    virtual Bool_t Run();
    virtual Double_t GetDeltaPtRandomCone();
    virtual Double_t GetDeltaPtRandomConeWithSignal();
    virtual Bool_t Notify();

    void SetPtHardBin(Int_t b) { fSelectPtHardBin = b; }
    void SetPtHardThreshold(Double_t b) { fPtHardThreshold = b; }

    Bool_t IsEventSelected();
    enum EPileup
    {
        kNoPileupSelection,
        kRejectPileupEvent,
        kRejectTracksFromPileupVertex
    };
    enum ERejBits
    {
        kNotSelTrigger,
        kNoVertex,
        kTooFewVtxContrib,
        kZVtxOutFid,
        kPileupSPD,
        kOutsideCentrality,
        kVertexZContrib,
        kPhysicsSelection,
        kNoContributors,
        kDeltaVertexZ,
        kNoVertexTracks,
        kVertexZResolution,
        kMVPileup,
        kSPDClusterCut,
        kZVtxSPDOutFid,
        kCentralityFlattening,
        kSelPtHardBin
    };
    Bool_t IsSelected(Int_t &WhyRejected, ULong_t &RejectionBits);
    void DoJetProbabilityAnalysis(Bool_t val = true) { fDoJetProbabilityAnalysis = val; }
    void DoCharmFractions(Bool_t val = true) { fDoCharmFractions = val; }
    void UsePartonDefinition(Bool_t val = true) { fUsePartonDef = val; }
    void DoJetMassAnalysis(Bool_t val = true) { fDoJetMass = val; }
    void DoSVEnergyFractionAnalysis(Bool_t val = true) { fDoSVEnergyFraction = val; }
    void DoPtRelAnalysis(Bool_t val = true) { fDoPtRelAnalysis = val; }
    void DoPtRelEventSelection(Bool_t val = true) { fDoSelectionPtRel = val; }

    void SetUseImpactParameterSignificance(bool val = true) { fUseIPs = val; }
    void SetUseNormalIPCalculation(bool val = true) { fUseNewIPsCalculation = val; }
    void SetCorrectResolution(bool val = true) { fCorrectRes = val; }

    void SetCorrectionFunctionPscat(TF1 *corrFunc, Int_t nITS) { fCorrectionFactorsPscat[nITS] = corrFunc; }
    void SetCorrectionFunctionNvtxContrib(TF1 *corrFunc, Int_t nITS) { fCorrectionFactorsNvtxContrib[nITS] = corrFunc; }
    Double_t CorrectResolution(Double_t resolution, Double_t pScat, Int_t nVtxContrib, Int_t nITS);

    void UseCorrectedRhoPt(bool val = true) { fUseCorrPt = val; }
    void UseGammaV0Rejection(bool val = true) { fEnableV0GammaRejection = val; }
    Bool_t IsV0InJet(TVector3 vV0, Double_t dJetPtMin);
    Bool_t IsElectronHF(AliAODTrack *track);
    void FillCandidates(Bool_t isK, Bool_t isL, Bool_t isAL, Int_t iCut);
    void FillCandidateJet(Int_t CutIndex, Int_t JetFlavor);

    // V0 selection
    void SetCutTPCRefit(Bool_t val = kTRUE) { fbTPCRefit = val; }
    void SetCutRejectKinks(Bool_t val = kTRUE) { fbRejectKinks = val; }
    void SetCutFindableClusters(Bool_t val = kTRUE) { fbFindableClusters = val; }
    void SetCutNCrossedRowsTPCMin(Double_t val = 70.) { fdCutNCrossedRowsTPCMin = val; }
    void SetCutCrossedRowsOverFindMin(Double_t val = 0.8) { fdCutCrossedRowsOverFindMin = val; }
    void SetCutCrossedRowsOverFindMax(Double_t val = 1e3) { fdCutCrossedRowsOverFindMax = val; }
    void SetCutPtDaughterMin(Double_t val = 0.150) { fdCutPtDaughterMin = val; }
    void SetCutDCAToPrimVtxMin(Double_t val = 0.1) { fdCutDCAToPrimVtxMin = val; }
    void SetCutDCADaughtersMax(Double_t val = 1.) { fdCutDCADaughtersMax = val; }
    void SetCutEtaDaughterMax(Double_t val = 0.8) { fdCutEtaDaughterMax = val; }
    void SetCutNSigmadEdxMax(Double_t val = 3.) { fdCutNSigmadEdxMax = val; }
    void SetPtProtonPIDMax(Double_t val = 1.) { fdPtProtonPIDMax = val; }
    void SetOnFly(Bool_t val = 0) { fbOnFly = val; }
    void SetCutCPAKMin(Double_t val = 0.998) { fdCutCPAKMin = val; }
    void SetCutCPALMin(Double_t val = 0.998) { fdCutCPALMin = val; }
    void SetCutRadiusDecayMin(Double_t val = 5.) { fdCutRadiusDecayMin = val; }
    void SetCutRadiusDecayMax(Double_t val = 100.) { fdCutRadiusDecayMax = val; }
    void SetCutEtaV0Max(Double_t val = 0.7) { fdCutEtaV0Max = val; }
    void SetCutRapV0Max(Double_t val = 0.75) { fdCutRapV0Max = val; }
    void SetCutNTauKMax(Double_t val = 5.0) { fdCutNTauKMax = val; }
    void SetCutNTauLMax(Double_t val = 5.0) { fdCutNTauLMax = val; }
    void SetCutArmPod(Bool_t val = kTRUE) { fbCutArmPod = val; }
    void SetCutCross(Bool_t val = kTRUE) { fbCutCross = val; }

    Bool_t SetResFunction(TF1 *f, Int_t j)
    {
        fResolutionFunction[j] = f;
        return kTRUE;
    }
    Bool_t SetResFunctionb(TF1 *f, Int_t j)
    {
        fResolutionFunctionb[j] = f;
        return kTRUE;
    }
    Bool_t SetResFunctionc(TF1 *f, Int_t j)
    {
        fResolutionFunctionc[j] = f;
        return kTRUE;
    }
    Bool_t SetResFunctionlf(TF1 *f, Int_t j)
    {
        fResolutionFunctionlf[j] = f;
        return kTRUE;
    }

    void ApplyV0Reconstruction(Bool_t val = kTRUE) { fApplyV0Rec = val; }
    void ApplyV0RejectionAll(Bool_t val = kTRUE) { fApplyV0RejectionAll = val; }

    void SetV0ReaderName(TString name)
    {
        fV0ReaderName = name;
        return;
    }

    void SetDoSVAnalysis(Bool_t value) { fDoSVAnalysis = value; }
    void SetDoTCAnalysis(Bool_t value) { fDoTrackCountingAnalysis = value; }
    void SetDoForthLargestIP(Bool_t value) { fDoForthIP = value; }

    void SetMinTrackProbability(Double_t value) { fMinTrackProb = value; }

    void SetDoImprovedDCACut(Bool_t value) { fDoImprovedDCACut = value; }

    void SetMaxFactorPtHardJet(Double_t value) { fMaxFactorPtHardJet = value; }

    void SetCalculateDCATruth(Bool_t value) { fCalcDCATruth = value; }

    void SetTaggerWorkingPoint(Double_t value) { fThresholdIP = value; }

    void SetDoDeltaPtWithSignal(Bool_t value) { fDoDeltaPtWithSignal = value; }

    void SetDoTaggedJetsDRM(Bool_t value) { fDoTaggedDRM = value; }

    Int_t FindVertices6Prong(const AliEmcalJet *jet,
                             TClonesArray *fTrackArrayIn,
                             AliAODEvent *aodEvent,
                             AliESDVertex *primaryESDVertex,
                             Double_t magZkG,
                             TClonesArray *arrayVtxHF,
                             Int_t &nDauRejCount);

    Int_t FindVertices5Prong(const AliEmcalJet *jet,
                             TClonesArray *fTrackArrayIn,
                             AliAODEvent *aodEvent,
                             AliESDVertex *primaryESDVertex,
                             Double_t magZkG,
                             TClonesArray *arrayVtxHF,
                             Int_t &nDauRejCount);

    Int_t FindVertices4Prong(const AliEmcalJet *jet,
                             TClonesArray *fTrackArrayIn,
                             AliAODEvent *aodEvent,
                             AliESDVertex *primaryESDVertex,
                             Double_t magZkG,
                             TClonesArray *arrayVtxHF,
                             Int_t &nDauRejCount);
    // B jet tracks selection
    void SetTrackMinPt(Double_t val) { fTCMinTrackPt = val; }
    void SetTPCClusterMin(Int_t val) { fTCMinClusTPC = val; }
    void SetITSHitsMin(Int_t val) { fTCMinHitsITS = val; }
    void SetMaxChi2perNDF(Double_t val) { fTCMaxChi2pNDF = val; }
    void SetMaxIPxy(Double_t val) { fTCMaxIPxy = val; }
    void SetMaxIPz(Double_t val) { fTCMaxIPz = val; }
    void SetMaxbDecayLength(Double_t val) { fTCMaxDecayLength = val; }
    void SetMaxDCATrackJet(Double_t val) { fTCMaxDCATrackJet = val; }

private:
    AliV0ReaderV1 *fV0Reader;    //!
    TString fV0ReaderName;       //
    TClonesArray *fReaderGammas; //!

    Bool_t SelectV0CandidateVIT();
    Bool_t IsParticleInCone(const AliVParticle *part, const AliEmcalJet *jet, Double_t dRMax) const; // decides whether a particle is inside a jet cone
    Bool_t CalculateTrackImpactParameter(AliAODTrack *track, double *impar, double *cov);            // Removes track from Vertex calculation first
    Bool_t CalculateTrackImpactParameterTruth(AliAODTrack *track, double *impar, double *cov);       // calculates DCA on MC particle/event information
    Bool_t CalculateJetSignedTrackImpactParameter(AliAODTrack *track, AliEmcalJet *jet, double *impar, double *cov, double &sign, double &dcajetrack, double &lineardecaylength);
    Double_t GetValImpactParameter(TTypeImpPar type, double *impar, double *cov);
    Bool_t IsV0PhotonFromBeamPipeDaughter(const AliAODTrack *track);
    Bool_t IsV0Daughter(const AliAODTrack *track);
    Bool_t IsTrackAccepted(AliAODTrack *track);
    Bool_t IsTrackAcceptedFidutial(AliAODTrack *track);
    Bool_t IsTrackAcceptedQuality(AliAODTrack *track, AliEmcalJet *Jet, Int_t &QualityClass, double *imp, double *cov, double &sign);
    Bool_t IsTrackAcceptedBJetCuts(AliAODTrack *track, Int_t jetFlavour);
    Bool_t MatchJetsGeometricDefault();                                                        //jet matching function 1/4
    void DoJetLoop();                                                                          //jet matching function 2/4
    void SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, int matching = 0);             //jet matching function 3/4
    void GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const; //jet matching function 4/4
    Double_t CalculateTrackProb(Double_t significance, Int_t trclass, Int_t jetFlavor);
    Double_t CalculateJetProb(AliEmcalJet *jet, Int_t jetFlavor); //!
    void FillResolutionFunctionHists(AliAODTrack *track, AliEmcalJet *jet, Int_t jetFlavor);

    AliAODMCParticle *GetMCTrack(const AliAODTrack *_track);

    Double_t GetPtRel(AliAODTrack *Lep, AliEmcalJet *jet, Bool_t AddLepToJet);
    void GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff);

    AliHFJetsTagging *fHFJetUtils; //!

    Double_t fPtHardThreshold;      //
    Bool_t fUseCorrPt;              //
    Bool_t fEnableV0GammaRejection; //
    Float_t fPythiaEventWeight;     //!
    Bool_t fDoImprovedDCACut;       //
    Bool_t fVertexConstraint;       //!
    Double_t fThresholdIP;          //
    Bool_t fDoDeltaPtWithSignal;    //

    AliESDVertex *fDiamond;       //!
    AliVertexerTracks *fVertexer; //!

    AliPIDResponse *fRespoPID; //!

    TRandom3 *fRandom; //! Random cone input

    TH1D *fh1dEventRejectionRDHFCuts; //! Store Rejection reasons and number of accepted events
    TH1D *fh1dVertexZ;                //!
    TH1D *fh1dVertexZAccepted;        //!
    TH2D *fh1dVertexR;                //!
    TH2D *fh1dVertexRAccepted;        //!

    // Bjet cuts
    Double_t fTCMinTrackPt;     //
    Int_t fTCMinClusTPC;        //
    Int_t fTCMinHitsITS;        //
    Double_t fTCMaxChi2pNDF;    //
    Double_t fTCMaxIPxy;        //
    Double_t fTCMaxIPz;         //
    Double_t fTCMaxDecayLength; //
    Double_t fTCMaxDCATrackJet; //

    Double_t fMaxFactorPtHardJet; //

    TH1D *fhistInclusiveJetCuts; //!
    TH1D *fhistbJetCuts;         //!
    TH1D *fhistcJetCuts;         //!
    TH1D *fhistlfJetCuts;        //!

    TH1D *fh1dTracksAccepeted; //!

    TH1D *fh1dJetRecPtAcceptedunCorr; //!

    TH2D *f2histRhoVsDeltaPt;       //!
    TH2D *f2histRhoVsDeltaPtFirst;  //!
    TH2D *f2histRhoVsDeltaPtSecond; //!
    TH2D *f2histRhoVsDeltaPtThird;  //!

    TH2D *f2histRhoVsDeltaPtWithSignal;       //!
    TH2D *f2histRhoVsDeltaPtWithSignalFirst;  //!
    TH2D *f2histRhoVsDeltaPtWithSignalSecond; //!
    TH2D *f2histRhoVsDeltaPtWithSignalThird;  //!

    TH1D *fh1dTracksImpParXY;              //! R Impact Parameter
    TH1D *fh1dTracksImpParXYZ;             //! R+z Impact Parameter
    TH1D *fh1dTracksImpParXYSignificance;  //! R Impact Parameter Significance
    TH1D *fh1dTracksImpParXYZSignificance; //! R+z Impact Parameter Significance

    TH1D *fh1dTracksImpParXYTruth;          //! R Impact Parameter truth
    TH1D *fh1dTracksImpParXYZTruth;         //! R+z Impact Parameter truth
    TH1D *fh1dTracksImpParXYResidualTruth;  //! Delta R Impact Parameter over uncertainty
    TH1D *fh1dTracksImpParXYZResidualTruth; //! Delta R+z Impact Parameter truth over uncertainty
    // inclusive  impact parameter distributions  for tracks with light meson (TODO D meson) correction in Monte Carlo

    TH2D *fh2dVertexChi2NDFNESDTracks; //!

    TH1D *fh1dJetGenPt;             //! Generator level jet momentum for unfolding
    TH1D *fh1dJetGenPtUnidentified; //!
    TH1D *fh1dJetGenPtudsg;         //!
    TH1D *fh1dJetGenPtc;            //!
    TH1D *fh1dJetGenPtb;            //!

    TH1D *fh1dJetRecPt;             //! Detector level jets
    TH1D *fh1dJetRecPtAccepted;     //! Detector level jets accepted
    TH2D *fh1dJetRecEtaPhiAccepted; //! Detector level jets accepted

    TH1D *fh1dJetRecPtUnidentified;         //!
    TH1D *fh1dJetRecPtudsg;                 //!
    TH1D *fh1dJetRecPtc;                    //!
    TH1D *fh1dJetRecPtb;                    //!
    TH1D *fh1dJetRecPtUnidentifiedAccepted; //!
    TH1D *fh1dJetRecPtudsgAccepted;         //!
    TH1D *fh1dJetRecPtcAccepted;            //!
    TH1D *fh1dJetRecPtbAccepted;            //!

    Bool_t fDoTaggedDRM;                //Flag whether to do tagged jets RDM
    TH2D *fh2dJetGenPtVsJetRecPt;       //! raw momentum response matrix
    TH2D *fh2dJetGenPtVsJetRecPtFirst;  //! raw momentum response matrix N=1
    TH2D *fh2dJetGenPtVsJetRecPtSecond; //! raw momentum response matrix N=2
    TH2D *fh2dJetGenPtVsJetRecPtThird;  //! raw momentum response matrix N=3
    TH2D *fh2dJetGenPtVsJetRecPtb;      //! b momentum response matrix
    TH2D *fh2dJetGenPtVsJetRecPtc;      //! c momentum response matrix
    TH2D *fh2dJetGenPtVsJetRecPtudsg;   //! udsg momentum response matrix
    // inclusive signed impact parameter distributions
    TH2D *fh2dJetSignedImpParXY;             //!
    TH2D *fh2dJetSignedImpParXYUnidentified; //!
    TH2D *fh2dJetSignedImpParXYudsg;         //!
    TH2D *fh2dJetSignedImpParXYb;            //!
    TH2D *fh2dJetSignedImpParXYc;            //!

    TH2D *fh2dJetSignedImpParXYSignificance;             //!
    TH2D *fh2dJetSignedImpParXYSignificanceUnidentified; //!
    TH2D *fh2dJetSignedImpParXYSignificanceudsg;         //!
    TH2D *fh2dJetSignedImpParXYSignificanceb;            //!
    TH2D *fh2dJetSignedImpParXYSignificancec;            //!

    TH2D *fh2dJetSignedImpParXYZ;             //!
    TH2D *fh2dJetSignedImpParXYZUnidentified; //!
    TH2D *fh2dJetSignedImpParXYZudsg;         //!
    TH2D *fh2dJetSignedImpParXYZb;            //!
    TH2D *fh2dJetSignedImpParXYZc;            //!

    TH2D *fh2dJetSignedImpParXYZSignificance;             //!
    TH2D *fh2dJetSignedImpParXYZSignificanceUnidentified; //!
    TH2D *fh2dJetSignedImpParXYZSignificanceudsg;         //!
    TH2D *fh2dJetSignedImpParXYZSignificanceb;            //!
    TH2D *fh2dJetSignedImpParXYZSignificancec;            //!

    //################################ Jet Probabilty
    TH2D *fh2dJetSignedImpParXY_Class1;              //!
    TH2D *fh2dJetSignedImpParXYSignificance_Class1;  //!
    TH2D *fh2dJetSignedImpParXYZ_Class1;             //!
    TH2D *fh2dJetSignedImpParXYZSignificance_Class1; //!

    TH2D *fh2dJetSignedImpParXYSignificanceb_Class1;  //!
    TH2D *fh2dJetSignedImpParXYSignificancec_Class1;  //!
    TH2D *fh2dJetSignedImpParXYSignificancelf_Class1; //!

    TH2D *fh2dJetSignedImpParXY_Class2;              //!
    TH2D *fh2dJetSignedImpParXYSignificance_Class2;  //!
    TH2D *fh2dJetSignedImpParXYZ_Class2;             //!
    TH2D *fh2dJetSignedImpParXYZSignificance_Class2; //!

    TH2D *fh2dJetSignedImpParXYSignificanceb_Class2;  //!
    TH2D *fh2dJetSignedImpParXYSignificancec_Class2;  //!
    TH2D *fh2dJetSignedImpParXYSignificancelf_Class2; //!

    TH2D *fh2dJetSignedImpParXY_Class3;              //!
    TH2D *fh2dJetSignedImpParXYSignificance_Class3;  //!
    TH2D *fh2dJetSignedImpParXYZ_Class3;             //!
    TH2D *fh2dJetSignedImpParXYZSignificance_Class3; //!

    TH2D *fh2dJetSignedImpParXYSignificanceb_Class3;  //!
    TH2D *fh2dJetSignedImpParXYSignificancec_Class3;  //!
    TH2D *fh2dJetSignedImpParXYSignificancelf_Class3; //!

    TH2D *fh2dJetSignedImpParXY_Class4;              //!
    TH2D *fh2dJetSignedImpParXYSignificance_Class4;  //!
    TH2D *fh2dJetSignedImpParXYZ_Class4;             //!
    TH2D *fh2dJetSignedImpParXYZSignificance_Class4; //!

    TH2D *fh2dJetSignedImpParXYSignificanceb_Class4;  //!
    TH2D *fh2dJetSignedImpParXYSignificancec_Class4;  //!
    TH2D *fh2dJetSignedImpParXYSignificancelf_Class4; //!

    //Jet Mass Histograms
    TH2D *fhistJetMass;              //!
    TH2D *fhistJetMass_Unidentified; //!
    TH2D *fhistJetMass_udsg;         //!
    TH2D *fhistJetMass_c;            //!
    TH2D *fhistJetMass_b;            //!

    TH2D *fhistJetMassFirst;  //!
    TH2D *fhistJetMassSecond; //!
    TH2D *fhistJetMassThird;  //!

    TH2D *fhistJetMass_UnidentifiedFirst; //!
    TH2D *fhistJetMass_udsgFirst;         //!
    TH2D *fhistJetMass_cFirst;            //!
    TH2D *fhistJetMass_bFirst;            //!

    TH2D *fhistJetMass_UnidentifiedSecond; //!
    TH2D *fhistJetMass_udsgSecond;         //!
    TH2D *fhistJetMass_cSecond;            //!
    TH2D *fhistJetMass_bSecond;            //!

    TH2D *fhistJetMass_UnidentifiedThird; //!
    TH2D *fhistJetMass_udsgThird;         //!
    TH2D *fhistJetMass_cThird;            //!
    TH2D *fhistJetMass_bThird;            //!

    //Secondary Vertex energy fraction
    TH2D *fhistSVEnergyFraction;              //!
    TH2D *fhistSVEnergyFraction_Unidentified; //!
    TH2D *fhistSVEnergyFraction_udsg;         //!
    TH2D *fhistSVEnergyFraction_c;            //!
    TH2D *fhistSVEnergyFraction_b;            //!

    TH2D *fhistSVEnergyFractionFirst;              //!
    TH2D *fhistSVEnergyFraction_UnidentifiedFirst; //!
    TH2D *fhistSVEnergyFraction_udsgFirst;         //!
    TH2D *fhistSVEnergyFraction_cFirst;            //!
    TH2D *fhistSVEnergyFraction_bFirst;            //!

    TH2D *fhistSVEnergyFractionSecond;              //!
    TH2D *fhistSVEnergyFraction_UnidentifiedSecond; //!
    TH2D *fhistSVEnergyFraction_udsgSecond;         //!
    TH2D *fhistSVEnergyFraction_cSecond;            //!
    TH2D *fhistSVEnergyFraction_bSecond;            //!

    TH2D *fhistSVEnergyFractionThird;              //!
    TH2D *fhistSVEnergyFraction_UnidentifiedThird; //!
    TH2D *fhistSVEnergyFraction_udsgThird;         //!
    TH2D *fhistSVEnergyFraction_cThird;            //!
    TH2D *fhistSVEnergyFraction_bThird;            //!

    TH2D *fhistSVnProngs;              //!
    TH2D *fhistSVnProngs_Unidentified; //!
    TH2D *fhistSVnProngs_udsg;         //!
    TH2D *fhistSVnProngs_c;            //!
    TH2D *fhistSVnProngs_b;            //!

    //Jet Probability Histograms
    TH2D *fhistJetProbability;              //!
    TH2D *fhistJetProbability_Unidentified; //!
    TH2D *fhistJetProbability_udsg;         //!
    TH2D *fhistJetProbability_c;            //!
    TH2D *fhistJetProbability_b;            //!

    TH2D *fhistJetProbabilityLog;              //!
    TH2D *fhistJetProbability_UnidentifiedLog; //!
    TH2D *fhistJetProbability_udsgLog;         //!
    TH2D *fhistJetProbability_cLog;            //!

    TH2D *fhistJetProbability_cLog_D0; //!
    TH2D *fhistJetProbability_cLog_Dp; //!
    TH2D *fhistJetProbability_cLog_Ds; //!
    TH2D *fhistJetProbability_cLog_Lc; //!

    TH2D *fhistJetProbability_bLog; //!

    TH2D *fhistJetProbabilityLogFirst;              //!
    TH2D *fhistJetProbability_UnidentifiedLogFirst; //!
    TH2D *fhistJetProbability_udsgLogFirst;         //!
    TH2D *fhistJetProbability_cLogFirst;            //!

    TH2D *fhistJetProbability_cLogFirst_D0; //!
    TH2D *fhistJetProbability_cLogFirst_Dp; //!
    TH2D *fhistJetProbability_cLogFirst_Ds; //!
    TH2D *fhistJetProbability_cLogFirst_Lc; //!

    TH2D *fhistJetProbability_bLogFirst; //!

    TH2D *fhistJetProbabilityLogSecond;              //!
    TH2D *fhistJetProbability_UnidentifiedLogSecond; //!
    TH2D *fhistJetProbability_udsgLogSecond;         //!
    TH2D *fhistJetProbability_cLogSecond;            //!

    TH2D *fhistJetProbability_cLogSecond_D0; //!
    TH2D *fhistJetProbability_cLogSecond_Dp; //!
    TH2D *fhistJetProbability_cLogSecond_Ds; //!
    TH2D *fhistJetProbability_cLogSecond_Lc; //!

    TH2D *fhistJetProbability_bLogSecond; //!

    TH2D *fhistJetProbabilityLogThird;              //!
    TH2D *fhistJetProbability_UnidentifiedLogThird; //!
    TH2D *fhistJetProbability_udsgLogThird;         //!
    TH2D *fhistJetProbability_cLogThird;            //!

    TH2D *fhistJetProbability_cLogThird_D0; //!
    TH2D *fhistJetProbability_cLogThird_Dp; //!
    TH2D *fhistJetProbability_cLogThird_Ds; //!
    TH2D *fhistJetProbability_cLogThird_Lc; //!

    TH2D *fhistJetProbability_bLogThird; //!

    TH2D *fhistJetProbabilityLogSVHE;              //!
    TH2D *fhistJetProbability_UnidentifiedLogSVHE; //!
    TH2D *fhistJetProbability_udsgLogSVHE;         //!
    TH2D *fhistJetProbability_cLogSVHE;            //!
    TH2D *fhistJetProbability_bLogSVHE;            //!

    TH2D *fhistJetProbabilityLogSVHP;              //!
    TH2D *fhistJetProbability_UnidentifiedLogSVHP; //!
    TH2D *fhistJetProbability_udsgLogSVHP;         //!
    TH2D *fhistJetProbability_cLogSVHP;            //!
    TH2D *fhistJetProbability_bLogSVHP;            //!

    Double_t fMinTrackProb;

    // inclusive signed impact parameter distributions
    //First
    TH2D *fh2dJetSignedImpParXYFirst;             //!
    TH2D *fh2dJetSignedImpParXYUnidentifiedFirst; //!
    TH2D *fh2dJetSignedImpParXYudsgFirst;         //!
    TH2D *fh2dJetSignedImpParXYbFirst;            //!
    TH2D *fh2dJetSignedImpParXYcFirst;            //!

    TH2D *fh2dJetSignedImpParXYSignificanceFirst;             //!
    TH2D *fh2dJetSignedImpParXYSignificanceUnidentifiedFirst; //!
    TH2D *fh2dJetSignedImpParXYSignificanceudsgFirst;         //!
    TH2D *fh2dJetSignedImpParXYSignificancebFirst;            //!
    TH2D *fh2dJetSignedImpParXYSignificancecFirst;            //!

    TH2D *fh2dJetSignedImpParXYZFirst;             //!
    TH2D *fh2dJetSignedImpParXYZUnidentifiedFirst; //!
    TH2D *fh2dJetSignedImpParXYZudsgFirst;         //!
    TH2D *fh2dJetSignedImpParXYZbFirst;            //!
    TH2D *fh2dJetSignedImpParXYZcFirst;            //!

    TH2D *fh2dJetSignedImpParXYZSignificanceFirst;             //!
    TH2D *fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst; //!
    TH2D *fh2dJetSignedImpParXYZSignificanceudsgFirst;         //!
    TH2D *fh2dJetSignedImpParXYZSignificancebFirst;            //!
    TH2D *fh2dJetSignedImpParXYZSignificancecFirst;            //!
    //Second
    TH2D *fh2dJetSignedImpParXYSecond;             //!
    TH2D *fh2dJetSignedImpParXYUnidentifiedSecond; //!
    TH2D *fh2dJetSignedImpParXYudsgSecond;         //!
    TH2D *fh2dJetSignedImpParXYbSecond;            //!
    TH2D *fh2dJetSignedImpParXYcSecond;            //!

    TH2D *fh2dJetSignedImpParXYSignificanceSecond;             //!
    TH2D *fh2dJetSignedImpParXYSignificanceUnidentifiedSecond; //!
    TH2D *fh2dJetSignedImpParXYSignificanceudsgSecond;         //!
    TH2D *fh2dJetSignedImpParXYSignificancebSecond;            //!
    TH2D *fh2dJetSignedImpParXYSignificancecSecond;            //!

    TH2D *fh2dJetSignedImpParXYZSecond;             //!
    TH2D *fh2dJetSignedImpParXYZUnidentifiedSecond; //!
    TH2D *fh2dJetSignedImpParXYZudsgSecond;         //!
    TH2D *fh2dJetSignedImpParXYZbSecond;            //!
    TH2D *fh2dJetSignedImpParXYZcSecond;            //!

    TH2D *fh2dJetSignedImpParXYZSignificanceSecond;             //!
    TH2D *fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond; //!
    TH2D *fh2dJetSignedImpParXYZSignificanceudsgSecond;         //!
    TH2D *fh2dJetSignedImpParXYZSignificancebSecond;            //!
    TH2D *fh2dJetSignedImpParXYZSignificancecSecond;            //!
    //Third
    TH2D *fh2dJetSignedImpParXYThird;             //!
    TH2D *fh2dJetSignedImpParXYUnidentifiedThird; //!
    TH2D *fh2dJetSignedImpParXYudsgThird;         //!
    TH2D *fh2dJetSignedImpParXYbThird;            //!
    TH2D *fh2dJetSignedImpParXYcThird;            //!

    TH2D *fh2dJetSignedImpParXYSignificanceThird;             //!
    TH2D *fh2dJetSignedImpParXYSignificanceUnidentifiedThird; //!
    TH2D *fh2dJetSignedImpParXYSignificanceudsgThird;         //!
    TH2D *fh2dJetSignedImpParXYSignificancebThird;            //!
    TH2D *fh2dJetSignedImpParXYSignificancecThird;            //!

    TH2D *fh2dJetSignedImpParXYZThird;             //!
    TH2D *fh2dJetSignedImpParXYZUnidentifiedThird; //!
    TH2D *fh2dJetSignedImpParXYZudsgThird;         //!
    TH2D *fh2dJetSignedImpParXYZbThird;            //!
    TH2D *fh2dJetSignedImpParXYZcThird;            //!

    TH2D *fh2dJetSignedImpParXYZSignificanceThird;             //!
    TH2D *fh2dJetSignedImpParXYZSignificanceUnidentifiedThird; //!
    TH2D *fh2dJetSignedImpParXYZSignificanceudsgThird;         //!
    TH2D *fh2dJetSignedImpParXYZSignificancebThird;            //!
    TH2D *fh2dJetSignedImpParXYZSignificancecThird;            //!

    //Forth
    Bool_t fDoForthIP;

    TH2D *fh2dJetSignedImpParXYForth;     //!
    TH2D *fh2dJetSignedImpParXYudsgForth; //!
    TH2D *fh2dJetSignedImpParXYbForth;    //!
    TH2D *fh2dJetSignedImpParXYcForth;    //!

    TH2D *fh2dJetSignedImpParXYSignificanceForth;     //!
    TH2D *fh2dJetSignedImpParXYSignificanceudsgForth; //!
    TH2D *fh2dJetSignedImpParXYSignificancebForth;    //!
    TH2D *fh2dJetSignedImpParXYSignificancecForth;    //!

    //V0
    TH2D *fh2dKshortMassVsPt;  //!
    TH2D *fh2dLamdaMassVsPt;   //!
    TH2D *fh2dAnLamdaMassVsPt; //!

    TH2D *fh2dKshortMassVsPtReal;  //!
    TH2D *fh2dLamdaMassVsPtReal;   //!
    TH2D *fh2dAnLamdaMassVsPtReal; //!

    TH2D *fh2dKshortRecPtVsGenPt;  //!
    TH2D *fh2dLamdaRecPtVsGenPt;   //!
    TH2D *fh2dAnLamdaRecPtVsGenPt; //!

    TH1D *fh1dPhotonPt; //!

    TH1D *fh1dKshortPtMC;  //!
    TH1D *fh1dLamdaPtMC;   //!
    TH1D *fh1dAnLamdaPtMC; //!

    TH2D *fh2dKshortPtVsJetPtMC;  //!
    TH2D *fh2dLamdaPtVsJetPtMC;   //!
    TH2D *fh2dAnLamdaPtVsJetPtMC; //!

    //V0 Reconstruction
    TH1D *fh1V0CounterCentK0s;     //! number of K0s candidates after various cuts
    TH1D *fh1V0CounterCentLambda;  //! number of Lambda candidates after various cuts
    TH1D *fh1V0CounterCentALambda; //! number of ALambda candidates after various cuts

    TClonesArray *fMCArray;           //!
    TClonesArray *fCaloClusters;      //! Tender cluster
    AliAnalysisUtils *fUtils;         //!
    Bool_t fDoJetProbabilityAnalysis; //
    Bool_t fDoCharmFractions;         // Flag for using different template for D0 Dp Ds Lc for reweighting
    Bool_t fUsePartonDef;             // Flag for using the parton definition to set the jet flavor
    Bool_t fUseIPs;                   // Flag for using the IPs instead of IP for tagging
    Bool_t fUseNewIPsCalculation;     // Flag whether to use to the ALICE IP calculation or keep the DELPHI way
    Bool_t fCorrectRes;               // Flag whether to correct the IP resolution interms of the primary vertex contributors and pScat
    Bool_t fDoJetMass;                //
    Bool_t fDoSVEnergyFraction;       //
    Bool_t fDoPtRelAnalysis;          //
    Bool_t fDoSelectionPtRel;         //
    Bool_t fUsePicoTracks;            //!

    TF1 *fCorrectionFactorsPscat[5];       //
    TF1 *fCorrectionFactorsNvtxContrib[5]; //

    // V0 selection
    // Daughter tracks
    Bool_t fbTPCRefit;                    // (yes) TPC refit for daughter tracks
    Bool_t fbRejectKinks;                 // (no) reject kink-like production vertices of daughter tracks
    Bool_t fbFindableClusters;            // (no) require positive number of findable clusters
    Double_t fdCutNCrossedRowsTPCMin;     // (70.) min number of crossed TPC rows
    Double_t fdCutCrossedRowsOverFindMin; // (0.8) min ratio crossed rows / findable clusters
    Double_t fdCutCrossedRowsOverFindMax; // (1e3) max ratio crossed rows / findable clusters
    Double_t fdCutPtDaughterMin;          // (0.150) [GeV/c] min transverse momentum of daughter tracks, to reject primaries which do not make it to the TPC
    Double_t fdCutDCAToPrimVtxMin;        // (0.1) [cm] min DCA of daughters to the prim vtx
    Double_t fdCutDCADaughtersMax;        // (1.) [sigma of TPC tracking] max DCA between daughters
    Double_t fdCutEtaDaughterMax;         // (0.8) max |pseudorapidity| of daughter tracks, historical reasons: tracking in MC for 2010 was restricted to 0.7
    Double_t fdCutNSigmadEdxMax;          // (3.) [sigma dE/dx] max difference between measured and expected signal of dE/dx in the TPC
    Double_t fdPtProtonPIDMax;            // (1.) [GeV/c] maxium pT of proton for applying PID cut in Pb-Pb
    // V0 candidate
    Bool_t fbOnFly;               // (0) on-the-fly (yes) or offline (no) reconstructed
    Double_t fdCutCPAKMin;        // (0.998) min cosine of the pointing angle, K0S
    Double_t fdCutCPALMin;        // (0.998) min cosine of the pointing angle, Lambda
    Double_t fdCutRadiusDecayMin; // (5.) [cm] min radial distance of the decay vertex
    Double_t fdCutRadiusDecayMax; // (100.) [cm] max radial distance of the decay vertex
    Double_t fdCutEtaV0Max;       // (0.7) max |pseudorapidity| of V0
    Double_t fdCutRapV0Max;       // (0.75) max |rapidity| of V0 (turned off)
    Double_t fdCutNTauKMax;       // (5.0) [tau] max proper lifetime in multiples of the mean lifetime, K0S
    Double_t fdCutNTauLMax;       // (5.0) [tau] max proper lifetime in multiples of the mean lifetime, Lambda
    Bool_t fbCutArmPod;           // (yes) Armenteros-Podolanski for K0S
    Bool_t fbCutCross;            // (no) cross-contamination

    THnSparse *fhnV0InJetK0s;     //! V0 in jet cones, in a centrality bin, m_V0; pt_V0; eta_V0; pt_jet
    THnSparse *fhnV0InJetLambda;  //!
    THnSparse *fhnV0InJetALambda; //!

    //########################## PtRel
    TH1D *fhistPtRelEvents;              //!
    TH2D *fhistPtRelVsJetPt;             //!
    TH2D *fhistLepIPVsJetPt;             //!
    TH2D *fhistPtRelVsJetPtUnidentified; //!
    TH2D *fhistPtRelVsJetPtudsg;         //!
    TH2D *fhistPtRelVsJetPtc;            //!
    TH2D *fhistPtRelVsJetPtb;            //!
    TH2D *fhistLepIPVsJetPtUnidentified; //!
    TH2D *fhistLepIPVsJetPtudsg;         //!
    TH2D *fhistLepIPVsJetPtc;            //!
    TH2D *fhistLepIPVsJetPtb;            //!
    TH2D *fHistMcEopEle;                 //!
    TH2D *fHistMcEopHad;                 //!
    TH2D *fTPCnsigMcEle;                 //!
    TH2D *fTPCnsigMcHad;                 //!

    TH2D *fhistPtRelVsJetPtTaggedFirst; //!
    TH2D *fhistLepIPVsJetPtTaggedFirst; //!

    TH2D *fhistPtRelVsJetPtTaggedUnidentifiedFirst; //!
    TH2D *fhistLepIPVsJetPtTaggedUnidentifiedFirst; //!

    TH2D *fhistPtRelVsJetPtTaggedudsgFirst; //!
    TH2D *fhistLepIPVsJetPtTaggedudsgFirst; //!

    TH2D *fhistPtRelVsJetPtTaggedcFirst; //!
    TH2D *fhistLepIPVsJetPtTaggedcFirst; //!

    TH2D *fhistPtRelVsJetPtTaggedbFirst; //!
    TH2D *fhistLepIPVsJetPtTaggedbFirst; //!

    TH2D *fhistPtRelVsJetPtTaggedSecond; //!
    TH2D *fhistLepIPVsJetPtTaggedSecond; //!

    TH2D *fhistPtRelVsJetPtTaggedUnidentifiedSecond; //!
    TH2D *fhistLepIPVsJetPtTaggedUnidentifiedSecond; //!

    TH2D *fhistPtRelVsJetPtTaggedudsgSecond; //!
    TH2D *fhistLepIPVsJetPtTaggedudsgSecond; //!

    TH2D *fhistPtRelVsJetPtTaggedcSecond; //!
    TH2D *fhistLepIPVsJetPtTaggedcSecond; //!

    TH2D *fhistPtRelVsJetPtTaggedbSecond; //!
    TH2D *fhistLepIPVsJetPtTaggedbSecond; //!

    TH2D *fhistPtRelVsJetPtTaggedThird; //!
    TH2D *fhistLepIPVsJetPtTaggedThird; //!

    TH2D *fhistPtRelVsJetPtTaggedUnidentifiedThird; //!
    TH2D *fhistLepIPVsJetPtTaggedUnidentifiedThird; //!

    TH2D *fhistPtRelVsJetPtTaggedudsgThird; //!
    TH2D *fhistLepIPVsJetPtTaggedudsgThird; //!

    TH2D *fhistPtRelVsJetPtTaggedcThird; //!
    TH2D *fhistLepIPVsJetPtTaggedcThird; //!

    TH2D *fhistPtRelVsJetPtTaggedbThird; //!
    TH2D *fhistLepIPVsJetPtTaggedbThird; //!

    Bool_t fApplyV0Rec;          //
    Bool_t fApplyV0RejectionAll; //

    TClonesArray *fV0CandidateArray;    //!
    AliJetContainer *fJetContainerMC;   //! Container with reconstructed jets
    AliJetContainer *fJetContainerData; //! Container with reconstructed jets
    AliAODEvent *fAODIn;                //! AOD Input Event
    AliAODVertex *fPrimaryVertex;       //! Event Primary Vertex

    TF1 *fResolutionFunction[7];   //
    TF1 *fResolutionFunctionb[7];  //
    TF1 *fResolutionFunctionc[7];  //
    TF1 *fResolutionFunctionlf[7]; //

    //Secondary Vertex
    Bool_t fDoSVAnalysis;            //
    Bool_t fDoTrackCountingAnalysis; //

    AliHFJetsTaggingVertex *fVtxTagger3Prong; //!
    AliHFJetsTaggingVertex *fVtxTagger2Prong; //!

    AliRDHFJetsCutsVertex *fjetCuts3Prong; //! SV cuts
    AliRDHFJetsCutsVertex *fjetCuts2Prong; //! SV cuts

    TObjArray *fTrackArray; //! Tracks selected for the SV reconstruction

    AliESDtrackCuts *fEsdTrackCuts; //! cuts used on the track selected for the SV reconstruction

    TH3D *fHistSV2Prong;             //!
    TH3D *fHistSV2ProngUnidentified; //!
    TH3D *fHistSV2Prongb;            //!
    TH3D *fHistSV2Prongc;            //!
    TH3D *fHistSV2Pronglf;           //!

    TH2D *fHistDispersion2Prong;             //!
    TH2D *fHistDispersion2ProngUnidentified; //!
    TH2D *fHistDispersion2Prongb;            //!
    TH2D *fHistDispersion2Prongc;            //!
    TH2D *fHistDispersion2Pronglf;           //!

    TH3D *fHistSV3Prong;             //!
    TH3D *fHistSV3ProngUnidentified; //!
    TH3D *fHistSV3Prongb;            //!
    TH3D *fHistSV3Prongc;            //!
    TH3D *fHistSV3Pronglf;           //!

    TH2D *fHistDispersion3Prong;             //!
    TH2D *fHistDispersion3ProngUnidentified; //!
    TH2D *fHistDispersion3Prongb;            //!
    TH2D *fHistDispersion3Prongc;            //!
    TH2D *fHistDispersion3Pronglf;           //!

    Double_t fInvariantMass; //!
    Double_t fJetMass;       //!
    Double_t fDispersion;    //!
    Double_t fDecayLength;   //!
    Double_t fLxySign;       //!
    Double_t fJetPt;         //!
    Int_t fJetFlavor;        //!
    Double_t fValJetProb;    //!
    Double_t fLogJetProb;    //!

    Bool_t fCalcDCATruth; //

    AliAnalysisTaskWeakDecayVertexer *fDecayVertex; //!

    static const Double_t fgkMassPion;    //
    static const Double_t fgkMassKshort;  //
    static const Double_t fgkMassLambda;  //
    static const Double_t fgkMassProton;  //
    static const Int_t fgkiNCategV0 = 18; // number of V0 selection steps

    ClassDef(AliAnalysisTaskBJetTC, 66)
};
#endif
//
