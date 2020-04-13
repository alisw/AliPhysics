#ifndef ALIANALYSISTASKSECHARMTRIGGERSTUDY_H
#define ALIANALYSISTASKSECHARMTRIGGERSTUDY_H

//**************************************************************************************
// \class AliAnalysisTaskSECharmTriggerStudy
// \brief task that produces an output tree for the charm trigger studies
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include <vector>

#include <TList.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include "AliAnalysisTaskSE.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODv0.h"
#include "AliEventCuts.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCutsLctoV0.h"

using std::vector;

struct Charm2Prong
{
    float fInvMassD0;                   /// inv mass of D0 hypothesis
    float fInvMassD0bar;                /// inv mass of D0bar hypothesis
    Double32_t fPt;                     //[0.0,65.535,16]
    Double32_t fY;                      //[-1.023,1.023,11]
    Double32_t fCosP;                   //[0.67233000,1.,15]
    Double32_t fCosPXY;                 //[0.67233000,1.,15]
    Double32_t fDecayLength;            //[0.0,6.5536,16]
    Double32_t fNormDecayLengthXY;      //[0.0,102.4,10]
    Double32_t fImpParProd;             //[-0.32768,0.32768,16]
    Double32_t fPtMinDau;               //[0.0,65.535,16]
    Double32_t fd0MinDau;               //[0.0,0.065536,16]
    int fGenLabel;                      /// label to match with MC
    int fProngIdx0;                     /// prong index 0
    int fProngIdx1;                     /// prong index 1
    unsigned char fCandType;            /// flag for cand type (signal, bkg, reflected, prompt, FD)
    unsigned char fDecay;               /// flag for decay channel
    int fSelBit;                        /// selection bit
};

struct Charm3Prong
{
    float fInvMassDplus;
    float fInvMassDstoKKpi;
    float fInvMassDstopiKK;
    float fInvMassLctopKpi;
    float fInvMassLctopiKp;
    float fInvMassPhiKKpi;
    float fInvMassPhipiKK;
    Double32_t fPt;                     //[0.0,65.535,16]
    Double32_t fYDplus;                 //[-1.023,1.023,11]
    Double32_t fYDs;                    //[-1.023,1.023,11]
    Double32_t fYLc;                    //[-1.023,1.023,11]
    Double32_t fCosP;                   //[0.67233000,1.,15]
    Double32_t fCosPXY;                 //[0.67233000,1.,15]
    Double32_t fDecayLength;            //[0.0,6.5535,16]
    Double32_t fNormDecayLengthXY;      //[0.0,102.3,10]
    Double32_t fSigmaVtx;               //[0.0,0.08190,12]
    Double32_t fPtMinDau;               //[0.0,65.535,16]
    Double32_t fd0MinDau;               //[0.0,0.065536,16]
    int fGenLabel;                      /// label to match with MC
    int fProngIdx0;                     /// prong index 0
    int fProngIdx1;                     /// prong index 1
    int fProngIdx2;                     /// prong index 2
    unsigned char fCandType;            /// flag for cand type (signal, bkg, reflected, prompt, FD)
    unsigned char fDecay;               /// flag for decay channel
    int fSelBit;                        /// selection bit
};

struct Dstar
{
    float fInvMass;
    float fInvMassD0;
    Double32_t fPt;                     //[0.0,65.535,16]
    Double32_t fY;                      //[-1.023,1.023,11]
    Double32_t fCosPD0;                 //[0.67233000,1.,15]
    Double32_t fCosPXYD0;               //[0.67233000,1.,15]
    Double32_t fDecayLengthD0;          //[0.0,6.5535,16]
    Double32_t fNormDecayLengthXYD0;    //[0.0,102.3,10]
    Double32_t fPtMinDau;               //[0.0,65.535,16]
    Double32_t fd0MinDau;               //[0.0,0.065536,16]
    int fGenLabel;                      /// label to match with MC
    int fProngIdx0;                     /// prong index 0
    int fProngIdx1;                     /// prong index 1
    int fProngIdx2;                     /// prong index 2
    unsigned char fCandType;            /// flag for cand type (signal, bkg, reflected, prompt, FD)
    unsigned char fDecay;               /// flag for decay channel
    int fSelBit;                        /// selection bit
};

struct CharmCascade
{
    float fInvMassLctopK0s;
    float fInvMassLctopiLambda;
    float fInvMassK0s;
    float fInvMassLambda;
    Double32_t fPt;                     //[0.0,65.535,16]
    Double32_t fY;                      //[-1.023,1.023,11]
    Double32_t fCosPV0;                 //[0.67233000,1.,15]
    Double32_t fCosPXYV0;               //[0.67233000,1.,15]
    Double32_t fPtMinDau;               //[0.0,65.535,16]
    Double32_t fd0MinDau;               //[0.0,0.065536,16]
    int fGenLabel;                      /// label to match with MC
    int fProngIdx0;                     /// prong index 0
    int fProngIdx1;                     /// prong index 1
    int fProngIdx2;                     /// prong index 2
    unsigned char fCandType;            /// flag for cand type (signal, bkg, reflected, prompt, FD)
    unsigned char fDecay;               /// flag for decay channel
    int fSelBit;                        /// selection bit
};

struct Beauty3Prong
{
    float fInvMassBplustoD0pi;
    float fInvMassD0;
    Double32_t fPt;                     //[0.0,65.535,16]
    Double32_t fY;                      //[-1.023,1.023,11]
    Double32_t fCosP;                   //[0.67233000,1.,15]
    Double32_t fCosPXY;                 //[0.67233000,1.,15]
    Double32_t fDecayLength;            //[0.0,6.5536,16]
    Double32_t fNormDecayLengthXY;      //[0.0,102.4,10]
    Double32_t fImpParProd;             //[-0.32768,0.32768,16]
    Double32_t fPtD0;                   //[0.0,65.535,16]
    Double32_t fCosPD0;                 //[0.67233000,1.,15]
    Double32_t fCosPXYD0;               //[0.67233000,1.,15]
    Double32_t fDecayLengthD0;          //[0.0,6.5536,16]
    Double32_t fNormDecayLengthXYD0;    //[0.0,102.4,10]
    Double32_t fImpParProdD0;           //[-0.32768,0.32768,16]
    int fGenLabel;                      /// label to match with MC
    unsigned char fCandType;            /// flag for cand type (signal, bkg, reflected, prompt, FD)
    unsigned char fDecay;               /// flag for decay channel
    int fSelBit;                        /// selection bit
};

struct Beauty4Prong
{
    float fInvMassB0toDminuspi;
    float fInvMassBstoDsminuspi;
    float fInvMassLbtoLcpluspi;
    float fInvMassDplus;
    float fInvMassDs;
    float fInvMassLc;
    Double32_t fPt;                     //[0.0,65.535,16]
    Double32_t fYB0;                    //[-1.023,1.023,11]
    Double32_t fYBs;                    //[-1.023,1.023,11]
    Double32_t fYLb;                    //[-1.023,1.023,11]
    Double32_t fCosP;                   //[0.67233000,1.,15]
    Double32_t fCosPXY;                 //[0.67233000,1.,15]
    Double32_t fDecayLength;            //[0.0,6.5536,16]
    Double32_t fNormDecayLengthXY;      //[0.0,102.4,10]
    Double32_t fImpParProd;             //[-0.32768,0.32768,16]
    Double32_t fPtD;                    //[0.0,65.535,16]
    Double32_t fCosPD;                  //[0.67233000,1.,15]
    Double32_t fCosPXYD;                //[0.67233000,1.,15]
    Double32_t fDecayLengthD;           //[0.0,6.5536,16]
    Double32_t fNormDecayLengthXYD;     //[0.0,102.4,10]
    Double32_t fSigmaVtxD;              //[0.0,0.08190,12]
    int fGenLabel;                      /// label to match with MC
    unsigned char fCandType;            /// flag for cand type (signal, bkg, reflected, prompt, FD)
    unsigned char fDecay;               /// flag for decay channel
    int fSelBit;                        /// selection bit
};

struct GenHadron
{
    float fPt;
    float fY;
    int fGenLabel;                      /// label to match with MC
    unsigned char fCandType;            /// flag for cand type (prompt, FD)
    unsigned char fDecay;               /// flag for decay channel
};

class AliAnalysisTaskSECharmTriggerStudy : public AliAnalysisTaskSE
{
public:
    enum kCandType
    {
        kSignal      = BIT(0),
        kBackground  = BIT(1),
        kPrompt      = BIT(2),
        kFeedDown    = BIT(3),
        kHasDauInAcc = BIT(4),
        kIsInFidAcc  = BIT(5)
    };

    enum kDecay
    {
        kNone,
        kDzerotoKpi,
        kDzerotopiK,
        kDplustoKpipi,
        kDstoKKpi,
        kDstopiKK,
        kDplustoKKpi,
        kDplustopiKK,
        kDstartoKpipi,
        kLctopKpi,
        kLctopiKp,
        kLctopK0s,
        kLctopiLambda,
        kBplustoD0pi,
        kB0toDminuspi,
        kBstoDsminuspi,
        kLbtoLcpluspi
    };

    enum kSelBit
    {
        kDzerotoKpiCuts      = BIT(0),
        kDzerotopiKCuts      = BIT(1),
        kDplustoKpipiCuts    = BIT(2),
        kDstartoKpipiCuts    = BIT(3),
        kDstoKKpiCuts        = BIT(4),
        kDstopiKKCuts        = BIT(5),
        kLctopKpiCuts        = BIT(6),
        kLctopiKpCuts        = BIT(7),
        kLctoV0bachCuts      = BIT(8),
        kBplustoD0piCuts     = BIT(9),
        kB0toDminuspiCuts    = BIT(10),
        kBstoDminuspiCuts    = BIT(11),
        kLbtoLcpluspiCuts    = BIT(12),
        kDzerotoKpiCutsPID   = BIT(13),
        kDzerotopiKCutsPID   = BIT(14),
        kDplustoKpipiCutsPID = BIT(15),
        kDstartoKpipiCutsPID = BIT(16),
        kDstoKKpiCutsPID     = BIT(17),
        kDstopiKKCutsPID     = BIT(18),
        kLctopKpiCutsPID     = BIT(19),
        kLctopiKpCutsPID     = BIT(20),
        kLctoV0bachCutsPID   = BIT(21),
        kDzerotoKpiFidAcc    = BIT(22),
        kDplustoKpipiFidAcc  = BIT(23),
        kDstartoKpipiFidAcc  = BIT(24),
        kDstoKKpiFidAcc      = BIT(25),
        kLctopKpiFidAcc      = BIT(26),
        kLctoV0bachFidAcc    = BIT(27),
        kBplustoD0piFidAcc   = BIT(28),
        kB0toDminuspiFidAcc  = BIT(29),
        kBstoDsminuspiFidAcc = BIT(30),
        kLbtoLcpluspiFidAcc  = BIT(31)
    };

    enum kSystem
    {
        kpp,
        kPbPb
    };

    AliAnalysisTaskSECharmTriggerStudy();
    AliAnalysisTaskSECharmTriggerStudy(const char *name, TList *cutlist);
    virtual ~AliAnalysisTaskSECharmTriggerStudy();

    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);

    void Enable2Prongs(bool enable = true)               {fEnable2Prongs = enable;}
    void Enable3Prongs(bool enableDplus = true,
                       bool enableDs = true,
                       bool enableLc = true)             {fEnable3Prongs = 0; if(enableDplus) fEnable3Prongs |= BIT(0); if(enableDs) fEnable3Prongs |= BIT(1); if(enableLc) fEnable3Prongs |= BIT(2);}
    void EnableDstars(bool enable = true)                {fEnableDstars = enable;}
    void EnableCascades(bool enable = true)              {fEnableCascades = enable;}
    void EnableBeauty3Prongs(bool enable = true)         {fEnableBeauty3Prongs = enable;}
    void EnableBeauty4Prongs(bool enableB0 = true,
                             bool enableBs = true,
                             bool enableLb = true)       {fEnableBeauty4Prongs = 0; if(enableB0) fEnableBeauty4Prongs |= BIT(0); if(enableBs) fEnableBeauty4Prongs |= BIT(1); if(enableLb) fEnableBeauty4Prongs |= BIT(2);}
    void SetFillOnlySignal(bool fillonlysignal = true)   {fFillOnlySignal = fillonlysignal;}
    void SetFillGenTree(bool fill = true)                {fFillGenTree = fill;}

    void SetSystem(int system = kpp)                     {fSystem = system;}
    void ApplyCuts(bool applycuts = true)                {fApplyCuts = applycuts;}
    void SetReadMC(bool readMC = true)                   {fReadMC = readMC;}

private:

    void FillCharm2Prong(AliAODRecoDecayHF2Prong* cand, int issel);
    void FillCharm3Prong(AliAODRecoDecayHF3Prong* cand, bool isselDplus, int isselDs, int isselLc);
    void FillDstar(AliAODRecoCascadeHF* cand, AliAODRecoDecayHF2Prong* dau, bool issel);
    void FillCharmCascade(AliAODRecoCascadeHF* cand, AliAODv0* dau, int issel);
    void FillBeauty3Prong(AliAODRecoDecayHF2Prong* cand, AliAODRecoDecayHF2Prong* dau, bool issel);
    void FillBeauty4Prong(AliAODRecoDecayHF2Prong* cand, AliAODRecoDecayHF3Prong* dau, bool isselB0, bool isselBs, bool isselLb, int isselDs, int isselLc);
    void FillGenerated(AliAODMCParticle* part, int origin, int decay, bool aredauinacc);
    bool RecalcOwnPrimaryVertex(AliAODRecoDecayHF* cand);
    void CleanOwnPrimaryVertex(AliAODRecoDecayHF* cand, AliAODVertex* origvtx);
    bool AreDauInAcc(int nProng, int *labDau);
    bool IsInFiducialAcceptance(double pt, double y);
    AliAODVertex* ReconstructDisplVertex(const AliVVertex *primary, TObjArray *tracks, double bField, double dispersion);

    TList *fOutput;                             //!<! List of output histograms
    TH1F* fHistNEvents;                         //!<! Histogram for event info
    TTree *fRecoTree;                           //!<! Output tree with reco candidates
    TTree *fGenTree;                            //!<! Output tree with generated particles

    AliEventCuts fEventCuts;                    /// object for event selection
    int fSystem;                                /// system (pp or PbPb)
    AliAODEvent *fAOD;                          //!<! AOD event
    int fAODProtection;                         /// protection for delta AOD mismatch
    TClonesArray* fMCArray;                     //!<! MC array
    float fRecoZvtx;                            /// Z of the reconstructed primary vertex
    float fGenZvtx;                             /// Z of the generated primary vertex
    int fNtracklets;                            /// number of tracklets in |eta| < 1

    vector<Charm2Prong> fCharm2Prong;           //!<! vector of charm 2 prongs
    vector<Charm3Prong> fCharm3Prong;           //!<! vector of charm 3 prongs
    vector<Dstar> fDstar;                       //!<! vector of Dstar
    vector<CharmCascade> fCharmCascade;         //!<! vector of charm cascades
    vector<Beauty3Prong> fBeauty3Prong;         //!<! vector of beauty 3 prongs
    vector<Beauty4Prong> fBeauty4Prong;         //!<! vector of beauty 4 prongs
    vector<GenHadron> fGenHadron;               //!<! vector of generated charm/beauty hadrons

    bool fEnable2Prongs;                        /// flag to enable 2-prong branch
    int fEnable3Prongs;                         /// flag to enable 3-prong branch (with D+ and/or Ds+ and/or Lc)
    bool fEnableDstars;                         /// flag to enable Dstar branch
    bool fEnableCascades;                       /// flag to enable cascade branch
    bool fEnableBeauty3Prongs;                  /// flag to enable B+
    int fEnableBeauty4Prongs;                   /// flag to enable B0 / Bs / Lb
    bool fFillOnlySignal;                       /// flag to fill only signal
    bool fFillGenTree;                          /// flag to fill tree with generated information
    bool fReadMC;                               /// flag to read MC info (if MC production)

    AliRDHFCutsD0toKpi* fCutsD0toKpi;           /// cut object for D0->Kpi
    AliRDHFCutsDplustoKpipi* fCutsDplustoKpipi; /// cut object for D+->Kpipi
    AliRDHFCutsDStartoKpipi* fCutsDstartoKpipi; /// cut object for D*+->D0pi->Kpipi
    AliRDHFCutsDstoKKpi* fCutsDstoKKpi;         /// cut object for Ds+->phipi->KKpi
    AliRDHFCutsLctopKpi* fCutsLctopKpi;         /// cut object for Lc->pKpi
    AliRDHFCutsLctoV0* fCutsLctoV0bach;         /// cut object for Lc->V0bachelor

    bool fApplyCuts;                            /// flag to enable cuts application
    TList* fListCuts;                           /// list of cut objects

    ClassDef(AliAnalysisTaskSECharmTriggerStudy, 1);
};

#endif