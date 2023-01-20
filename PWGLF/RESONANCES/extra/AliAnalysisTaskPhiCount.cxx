// !TODO LIST
// TODO: Check Overlapping histograms
// TODO: Tracklets up to 5k and rebinning (5k bins)

#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TTree.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODInputHandler.h"
#include "AliMCVertex.h"
#include "AliAODMCParticle.h"
#include "AliAODHeader.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliPPVsMultUtils.h"
#include "AliAnalysisTaskPhiCount.h"
#include "AliESDtrackCuts.h"

class AliAnalysisTaskPhiCount;
class AliPIDResponse;
 
using namespace std;

ClassImp(AliAnalysisTaskPhiCount)

            AliAnalysisTaskPhiCount::AliAnalysisTaskPhiCount() : AliAnalysisTaskSE(), kMCbool(0), kPhibool(0), kKaonbool(0), kSpherocityPTWeight(0), kComputeRT(0), kFilterBit(5), kDCAzCut(2.), kNSigmaPtDepXYDCA(7.), kMinTPCclusters(70.), kChi2TPCcluster(4.), kChi2TPCGlobal(36.), kChi2ITScluster(36.), kTPCClsOverFndbl(0), fIs_p_p(1), fIs_p_Pb(0), fIs_Pb_Pb(0), fAOD(0), fESD(0), fMCD(0), AODMCTrackArray(0), fPIDResponse(0), fPrimaryVertex(0), fMultSelection(0), fCurrent_Track(0), fCurrent_Track_Charge(0), fCurrent_Track_Momentum(0),fCurrent_Track_TransMom(0), fCurrent_Track_Eta(0), fCurrent_Track_Phi(0), fCurrent_Track_DCAXY(0), fCurrent_Track_DCAZ(0), fCurrent_SPH(0), fCurrent_RT(0), fIsTPCAvailable(0), fIsTOFAvailable(0), fBetaFromTOFSignal(0), fTPCSignal(0), fQCOutputList(0), fRunName(0), kTriggerMask(0), kSgTPC_Alone(0), kSgTPC_TOFVt(0), kSgTOF_Veto(0), fCurrent_V0M(0), fCurrent_TRK(0), fCurrent_Run(0), fnPhiRec(0), fEventMask(0), fIsINELgt0(0), fTrueEventMask(0), fIsTrueINELgt0(0), fPhiCandidate(0), fnPhi(0), fKaonCandidate(0), fnKaon(0), fPhiEfficiency(0), fnPhiTru(0), fKaonEfficiency(0), fnKaonTru(0), fAnalysisOutputList(0), kVertexCut(0), /* fQC Histograms*/ fQC_Event_Enum_FLL(0), fQC_Event_Enum_V0M(0), fQC_Event_Enum_TRK(0), fQC_Event_Vertex_Fll(0), fQC_Event_Vertex_Cut(0), fQC_Tracks_Momentum(0), fQC_Tracks_TMomentum(0), fQC_Tracks_Eta(0), fQC_Tracks_Phi(0), fQC_Tracks_V0M(0), fQC_Tracks_TRK(0), fQC_Tracks_DCAXY_P(0), fQC_Tracks_DCAXY_PT(0), fQC_Tracks_DCAZ_P(0), fQC_Tracks_DCAZ_PT(0), fQC_Tracks_TOF_P(0), fQC_Tracks_TOF_PT(0), fQC_Tracks_TPC_P(0), fQC_Tracks_TPC_PT(0), fQC_Tracks_TPC_CLS(0), fQC_Tracks_TPC_FRC(0), fQC_Tracks_TPC_CHI(0), fQC_Tracks_TPC_CNS(0), fQC_Tracks_ITS_CHI(0), fQC_Kaons_Momentum(0), fQC_Kaons_TMomentum(0), fQC_Kaons_Eta(0), fQC_Kaons_Phi(0), fQC_Kaons_V0M(0), fQC_Kaons_TRK(0), fQC_Kaons_DCAXY_P(0), fQC_Kaons_DCAXY_PT(0), fQC_Kaons_DCAZ_P(0), fQC_Kaons_DCAZ_PT(0), fQC_Kaons_TOF_P(0), fQC_Kaons_TOF_PT(0), fQC_Kaons_TPC_P(0), fQC_Kaons_TPC_PT(0), fQC_PID_TOF_Kaons_P(0), fQC_PID_TOF_Kaons_PT(0), fQC_PID_TPC_Kaons_P(0), fQC_PID_TPC_Kaons_PT(0), fQC_PID_TPC_TOF_Kaons_PT(0), fQC_PID_TOF_NSig_SEL_Kaons_P(0), fQC_PID_TOF_NSig_SEL_Kaons_PT(0), fQC_PID_TPC_NSig_SEL_Kaons_P(0), fQC_PID_TPC_NSig_SEL_Kaons_PT(0), fQC_PID_TOF_Sgnl_SEL_Kaons_P(0), fQC_PID_TOF_Sgnl_SEL_Kaons_PT(0), fQC_PID_TPC_Sgnl_SEL_Kaons_P(0), fQC_PID_TPC_Sgnl_SEL_Kaons_PT(0), fQC_Event_Enum_V0T(0), fQC_Phi_InvMass_Rec(0), fQC_Phi_InvMass_Gen(0), fQC_Phi_InvMass_Eff(0), fQC_Event_Spherocity(0), fQC_Event_FullRT(0), fQC_Event_FullRF(0), fQC_Event_FullRB(0) {
            }

//_____________________________________________________________________________

            AliAnalysisTaskPhiCount::AliAnalysisTaskPhiCount(const char* name) : AliAnalysisTaskSE(name), kMCbool(0), kPhibool(0), kKaonbool(0), kSpherocityPTWeight(0), kComputeRT(0), kFilterBit(5), kDCAzCut(2.), kNSigmaPtDepXYDCA(7.), kMinTPCclusters(70.), kChi2TPCcluster(4.), kChi2TPCGlobal(36.), kChi2ITScluster(36.), kTPCClsOverFndbl(0), fIs_p_p(1), fIs_p_Pb(0), fIs_Pb_Pb(0), fAOD(0), fESD(0), fMCD(0), AODMCTrackArray(0), fPIDResponse(0), fPrimaryVertex(0), fMultSelection(0), fCurrent_Track(0), fCurrent_Track_Charge(0), fCurrent_Track_Momentum(0),fCurrent_Track_TransMom(0), fCurrent_Track_Eta(0), fCurrent_Track_Phi(0), fCurrent_Track_DCAXY(0), fCurrent_Track_DCAZ(0), fCurrent_SPH(0), fCurrent_RT(0), fIsTPCAvailable(0), fIsTOFAvailable(0), fBetaFromTOFSignal(0), fTPCSignal(0), fQCOutputList(0), fRunName(0), kTriggerMask(0), kSgTPC_Alone(0), kSgTPC_TOFVt(0), kSgTOF_Veto(0), fCurrent_V0M(0), fCurrent_TRK(0), fCurrent_Run(0), fnPhiRec(0), fEventMask(0), fIsINELgt0(0), fTrueEventMask(0), fIsTrueINELgt0(0), fPhiCandidate(0), fnPhi(0), fKaonCandidate(0), fnKaon(0), fPhiEfficiency(0), fnPhiTru(0), fKaonEfficiency(0), fnKaonTru(0), fAnalysisOutputList(0), kVertexCut(0), /* fQC Histograms*/ fQC_Event_Enum_FLL(0), fQC_Event_Enum_V0M(0), fQC_Event_Enum_TRK(0), fQC_Event_Vertex_Fll(0), fQC_Event_Vertex_Cut(0), fQC_Tracks_Momentum(0), fQC_Tracks_TMomentum(0), fQC_Tracks_Eta(0), fQC_Tracks_Phi(0), fQC_Tracks_V0M(0), fQC_Tracks_TRK(0), fQC_Tracks_DCAXY_P(0), fQC_Tracks_DCAXY_PT(0), fQC_Tracks_DCAZ_P(0), fQC_Tracks_DCAZ_PT(0), fQC_Tracks_TOF_P(0), fQC_Tracks_TOF_PT(0), fQC_Tracks_TPC_P(0), fQC_Tracks_TPC_PT(0), fQC_Tracks_TPC_CLS(0), fQC_Tracks_TPC_FRC(0), fQC_Tracks_TPC_CHI(0), fQC_Tracks_TPC_CNS(0), fQC_Tracks_ITS_CHI(0), fQC_Kaons_Momentum(0), fQC_Kaons_TMomentum(0), fQC_Kaons_Eta(0), fQC_Kaons_Phi(0), fQC_Kaons_V0M(0), fQC_Kaons_TRK(0), fQC_Kaons_DCAXY_P(0), fQC_Kaons_DCAXY_PT(0), fQC_Kaons_DCAZ_P(0), fQC_Kaons_DCAZ_PT(0), fQC_Kaons_TOF_P(0), fQC_Kaons_TOF_PT(0), fQC_Kaons_TPC_P(0), fQC_Kaons_TPC_PT(0), fQC_PID_TOF_Kaons_P(0), fQC_PID_TOF_Kaons_PT(0), fQC_PID_TPC_Kaons_P(0), fQC_PID_TPC_Kaons_PT(0), fQC_PID_TPC_TOF_Kaons_PT(0), fQC_PID_TOF_NSig_SEL_Kaons_P(0), fQC_PID_TOF_NSig_SEL_Kaons_PT(0), fQC_PID_TPC_NSig_SEL_Kaons_P(0), fQC_PID_TPC_NSig_SEL_Kaons_PT(0), fQC_PID_TOF_Sgnl_SEL_Kaons_P(0), fQC_PID_TOF_Sgnl_SEL_Kaons_PT(0), fQC_PID_TPC_Sgnl_SEL_Kaons_P(0), fQC_PID_TPC_Sgnl_SEL_Kaons_PT(0), fQC_Event_Enum_V0T(0), fQC_Phi_InvMass_Rec(0), fQC_Phi_InvMass_Gen(0), fQC_Phi_InvMass_Eff(0), fQC_Event_Spherocity(0), fQC_Event_FullRT(0), fQC_Event_FullRF(0), fQC_Event_FullRB(0) {
// Define Input
DefineInput(0, TChain::Class());

// Define Output
DefineOutput(1, TList::Class());
DefineOutput(2, TList::Class());
DefineOutput(3, TTree::Class());
DefineOutput(4, TTree::Class());
DefineOutput(5, TTree::Class());
DefineOutput(6, TTree::Class());
}

//_____________________________________________________________________________

            AliAnalysisTaskPhiCount::~AliAnalysisTaskPhiCount()                 {
    // Deleting TLists
    if( fAnalysisOutputList )
    {
        delete fAnalysisOutputList;
    }
    if( fQCOutputList )
    {
        delete fQCOutputList;
    }
    if( fKaonCandidate )
    {
        delete fKaonCandidate;
    }
    if( fPhiCandidate )
    {
        delete fPhiCandidate;
    }
    if( fKaonEfficiency )
    {
        delete fKaonEfficiency;
    }
    if( fPhiEfficiency )
    {
        delete fPhiEfficiency;
    }
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::UserCreateOutputObjects()                  {
    //  PPvsMult
    fMultSelection          =   new AliPPVsMultUtils();
    
    // Analysis Output Histograms TList initialisation
    fAnalysisOutputList     = new TList();
    fAnalysisOutputList     ->SetOwner(kTRUE);
    PostData(1, fAnalysisOutputList);
    
    // QC utility Histograms TList initialisation
    fQCOutputList   = new TList();
    fQCOutputList   ->SetOwner(kTRUE);
    
    //_____________________________________________________________________________
    //
    //  --- Event
    //  --- --- Multiplicity
    //  --- --- --- V0M
    Int_t   kV0M_bin    =   100000;
    Float_t kV0M_min    =   -100.;
    Float_t kV0M_max    =   +900.;
    //
    //  --- --- --- Tracklets
    Int_t   kTRK_bin    =   220;
    Float_t kTRK_min    =   -20.;
    Float_t kTRK_max    =   +200.;
    //
    //  --- --- Z Vertex
    Int_t   zVTX_bin    =   300;
    Float_t zVTX_min    =   -15.;
    Float_t zVTX_max    =   +15.;
    //
    //  --- Track
    //  --- --- Transverse momentum
    Int_t   kMOM_bin    =   400;
    Float_t kMOM_min    =   +0.;
    Float_t kMOM_max    =   +10.;
    //
    //  --- --- Phi angle
    Int_t   kPhA_bin    =   370;
    Float_t kPhA_min    =   -5.;
    Float_t kPhA_max    =   +365.;
    //
    //  --- --- Eta
    Int_t   kETA_bin    =   100;
    Float_t kETA_min    =   -1.;
    Float_t kETA_max    =   +1.;
    //
    //  --- --- Nsgima PID
    Int_t   sPID_bin    =   200;
    Float_t sPID_min    =   -10.;
    Float_t sPID_max    =   +10.;
    //
    //_____________________________________________________________________________
    fQC_Event_Enum_FLL              = new TH1D("fQC_Event_Enum_FLL",       "Event Selection",                                  29, -0.5, 28.5);
    fQC_Event_Enum_FLL              ->  GetXaxis()  ->  SetTitle("");
    fQC_Event_Enum_FLL              ->  GetYaxis()  ->  SetTitle("Accepted Events");
    fSetEventCountLabels(fQC_Event_Enum_FLL);
    fQCOutputList->Add(fQC_Event_Enum_FLL);
    
    fQC_Event_Enum_V0M              = new TH1D("fQC_Event_Enum_V0M",        "Acc. Events in Mult.",                             kV0M_bin, kV0M_min, kV0M_max);
    fQC_Event_Enum_V0M              ->  GetXaxis()  ->  SetTitle("V0M Multiplicity");
    fQC_Event_Enum_V0M              ->  GetYaxis()  ->  SetTitle("Accepted Event");
    fQCOutputList->Add(fQC_Event_Enum_V0M);
    
    fQC_Event_Enum_TRK              = new TH1D("fQC_Event_Enum_TRK",        "Acc. Events in Mult.",                             kTRK_bin, kTRK_min, kTRK_max);
    fQC_Event_Enum_TRK              ->  GetXaxis()  ->  SetTitle("Tracklets Multiplicity");
    fQC_Event_Enum_TRK              ->  GetYaxis()  ->  SetTitle("Accepted Event");
    fQCOutputList->Add(fQC_Event_Enum_TRK);
    
    fQC_Event_Enum_V0T              = new TH2D("fQC_Event_Enum_V0T",        "Acc. Events in Mult.",                             kV0M_bin, kV0M_min, kV0M_max, kTRK_bin, kTRK_min, kTRK_max);
    fQC_Event_Enum_V0T              ->  GetXaxis()  ->  SetTitle("V0M Multiplicity");
    fQC_Event_Enum_V0T              ->  GetYaxis()  ->  SetTitle("Tracklets Multiplicity");
    fQC_Event_Enum_V0T              ->  GetZaxis()  ->  SetTitle("Accepted Event");
    fQCOutputList->Add(fQC_Event_Enum_V0T);
    
    fQC_Event_Vertex_Fll            = new TH1F("fQC_Event_Vertex_Fll",      "Collision Vertex (FULL)",                          zVTX_bin, zVTX_min, zVTX_max);
    fQC_Event_Vertex_Fll            ->  GetXaxis()  ->  SetTitle("Vertex z-position (cm)");
    fQC_Event_Vertex_Fll            ->  GetYaxis()  ->  SetTitle("Accepted Events");
    fQCOutputList->Add(fQC_Event_Vertex_Fll);
    
    fQC_Event_Vertex_Cut            = new TH1F("fQC_Event_Vertex_Cut",      "Collision Vertex (CUTS)",                          zVTX_bin, zVTX_min, zVTX_max);
    fQC_Event_Vertex_Cut            ->  GetXaxis()  ->  SetTitle("Vertex z-position (cm)");
    fQC_Event_Vertex_Fll            ->  GetYaxis()  ->  SetTitle("Accepted Events");
    fQCOutputList->Add(fQC_Event_Vertex_Cut);
    
    fQC_Event_Spherocity            = new TH1F("fQC_Event_Spherocity",      "Event Spherocity",                                 410, -3, 1.1);
    fQC_Event_Spherocity            ->  GetXaxis()  ->  SetTitle("Spherocity");
    fQC_Event_Spherocity            ->  GetYaxis()  ->  SetTitle("Accepted Events");
    fQCOutputList->Add(fQC_Event_Spherocity);
    
    fQC_Event_FullRF                = new TH1F("fQC_Event_FullRF",          "Event RForward",                                   503, -3, 500);
    fQC_Event_FullRF                ->  GetXaxis()  ->  SetTitle("RForward");
    fQC_Event_FullRF                ->  GetYaxis()  ->  SetTitle("Accepted Events");
    fQCOutputList->Add(fQC_Event_FullRF);
    
    fQC_Event_FullRB                = new TH1F("fQC_Event_FullRB",          "Event RBackward",                                  503, -3, 500);
    fQC_Event_FullRB                ->  GetXaxis()  ->  SetTitle("RBackward");
    fQC_Event_FullRB                ->  GetYaxis()  ->  SetTitle("Accepted Events");
    fQCOutputList->Add(fQC_Event_FullRB);
    
    fQC_Event_FullRT                = new TH1F("fQC_Event_FullRT",          "Event RTransverse",                                503, -3, 500);
    fQC_Event_FullRT                ->  GetXaxis()  ->  SetTitle("RTransverse");
    fQC_Event_FullRT                ->  GetYaxis()  ->  SetTitle("Accepted Events");
    fQCOutputList->Add(fQC_Event_FullRT);
    
    //_____________________________________________________________________________
    //
    //  TRACKS
    //_____________________________________________________________________________
    fQC_Tracks_Momentum             = new TH2F("fQC_Tracks_Momentum",       "Acc. Tracks Momentum",                             2, -1., 1., kMOM_bin, kMOM_min, kMOM_max);
    fQC_Tracks_Momentum             ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_Tracks_Momentum             ->  GetYaxis()  ->  SetTitle("Track Momentum (GeV/c)");
    fQC_Tracks_Momentum             ->  GetZaxis()  ->  SetTitle("Accepted Tracks");
    fQCOutputList->Add(fQC_Tracks_Momentum);
    
    fQC_Tracks_TMomentum            = new TH2F("fQC_Tracks_TMomentum",      "Acc. Tracks Transverse Momentum",                  2, -1., 1., kMOM_bin, kMOM_min, kMOM_max);
    fQC_Tracks_TMomentum            ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_Tracks_TMomentum            ->  GetYaxis()  ->  SetTitle("Track Transverse Momentum (GeV/c)");
    fQC_Tracks_TMomentum            ->  GetZaxis()  ->  SetTitle("Accepted Tracks");
    fQCOutputList->Add(fQC_Tracks_TMomentum);
    
    fQC_Tracks_Eta                  = new TH2F("fQC_Tracks_Eta",            "Acc. Tracks #eta",                                 2, -1., 1., kETA_bin, kETA_min, kETA_max);
    fQC_Tracks_Eta                  ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_Tracks_Eta                  ->  GetYaxis()  ->  SetTitle("Track #eta");
    fQC_Tracks_Eta                  ->  GetZaxis()  ->  SetTitle("Accepted Tracks");
    fQCOutputList->Add(fQC_Tracks_Eta);
    
    fQC_Tracks_Phi                  = new TH2F("fQC_Tracks_Phi",            "Acc. Tracks #phi",                                 2, -1., 1., kPhA_bin, kPhA_min, kPhA_max);
    fQC_Tracks_Phi                  ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_Tracks_Phi                  ->  GetYaxis()  ->  SetTitle("Track #varphi");
    fQC_Tracks_Phi                  ->  GetZaxis()  ->  SetTitle("Accepted Tracks");
    fQCOutputList->Add(fQC_Tracks_Phi);
    
    fQC_Tracks_V0M                  = new TH2F("fQC_Tracks_V0M",            "Acc. Tracks in Mult.",                             2, -1., 1., kV0M_bin, kV0M_min, kV0M_max);
    fQC_Tracks_V0M                  ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_Tracks_V0M                  ->  GetYaxis()  ->  SetTitle("V0M Multiplicity");
    fQC_Tracks_V0M                  ->  GetZaxis()  ->  SetTitle("Accepted Tracks");
    fQCOutputList->Add(fQC_Tracks_V0M);
    
    fQC_Tracks_TRK                  = new TH2F("fQC_Tracks_TRK",            "Acc. Tracks in Mult.",                             2, -1., 1., kTRK_bin, kTRK_min, kTRK_max);
    fQC_Tracks_TRK                  ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_Tracks_TRK                  ->  GetYaxis()  ->  SetTitle("Tracklets Multiplicity");
    fQC_Tracks_TRK                  ->  GetZaxis()  ->  SetTitle("Accepted Tracks");
    fQCOutputList->Add(fQC_Tracks_TRK);
    
    fQC_Tracks_DCAXY_P              = new TH3F("fQC_Tracks_DCAXY_P",        "Acc. Tracks XY-DCA",                               2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, 1000, -5., 5.);
    fQC_Tracks_DCAXY_P              ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_Tracks_DCAXY_P              ->  GetYaxis()  ->  SetTitle("Track Momentum (GeV/c)");
    fQC_Tracks_DCAXY_P              ->  GetZaxis()  ->  SetTitle("XY-DCA (cm)");
    fQCOutputList->Add(fQC_Tracks_DCAXY_P);
    
    fQC_Tracks_DCAXY_PT             = new TH3F("fQC_Tracks_DCAXY_PT",       "Acc. Tracks XY-DCA",                               2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, 1000, -5., 5.);
    fQC_Tracks_DCAXY_PT             ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_Tracks_DCAXY_PT             ->  GetYaxis()  ->  SetTitle("Track Transverse Momentum (GeV/c)");
    fQC_Tracks_DCAXY_PT             ->  GetZaxis()  ->  SetTitle("XY-DCA (cm)");
    fQCOutputList->Add(fQC_Tracks_DCAXY_PT);
    
    fQC_Tracks_DCAZ_P               = new TH3F("fQC_Tracks_DCAZ_P",          "Acc. Tracks Z-DCA",                               2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, 1000, -5., 5.);
    fQC_Tracks_DCAZ_P               ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_Tracks_DCAZ_P               ->  GetYaxis()  ->  SetTitle("Track Momentum (GeV/c)");
    fQC_Tracks_DCAZ_P               ->  GetZaxis()  ->  SetTitle("Z-DCA (cm)");
    fQCOutputList->Add(fQC_Tracks_DCAZ_P);
    
    fQC_Tracks_DCAZ_PT              = new TH3F("fQC_Tracks_DCAZ_PT",        "Acc. Tracks Z-DCA",                                2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, 1000, -5., 5.);
    fQC_Tracks_DCAZ_PT              ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_Tracks_DCAZ_PT              ->  GetYaxis()  ->  SetTitle("Track Transverse Momentum (GeV/c)");
    fQC_Tracks_DCAZ_PT              ->  GetZaxis()  ->  SetTitle("Z-DCA (cm)");
    fQCOutputList->Add(fQC_Tracks_DCAZ_PT);
    
    fQC_Tracks_TOF_P                = new TH3F("fQC_Tracks_TOF_P",           "Acc. Tracks TOF Sig.",                            2, -1., 1.,  kMOM_bin, kMOM_min, kMOM_max, 120, 0., 1.2);
    fQC_Tracks_TOF_P                ->  GetXaxis()  ->  SetTitle("Kaon Sign");
    fQC_Tracks_TOF_P                ->  GetYaxis()  ->  SetTitle("Kaon Momentum (GeV/c)");
    fQC_Tracks_TOF_P                ->  GetZaxis()  ->  SetTitle("#beta");
    fQCOutputList->Add(fQC_Tracks_TOF_P);
    
    fQC_Tracks_TOF_PT               = new TH3F("fQC_Tracks_TOF_PT",          "Acc. Tracks TOF Sig.",                            2, -1., 1.,  kMOM_bin, kMOM_min, kMOM_max, 120, 0., 1.2);
    fQC_Tracks_TOF_PT               ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_Tracks_TOF_PT               ->  GetYaxis()  ->  SetTitle("Track Transverse Momentum (GeV/c)");
    fQC_Tracks_TOF_PT               ->  GetZaxis()  ->  SetTitle("#beta");
    fQCOutputList->Add(fQC_Tracks_TOF_PT);
    
    fQC_Tracks_TPC_P                = new TH3F("fQC_Tracks_TPC_P",           "Acc. Tracks TPC Sig.",                            2, -1., 1.,  kMOM_bin, kMOM_min, kMOM_max, 100, 0., 500.);
    fQC_Tracks_TPC_P                ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_Tracks_TPC_P                ->  GetYaxis()  ->  SetTitle("Track Momentum (GeV/c)");
    fQC_Tracks_TPC_P                ->  GetZaxis()  ->  SetTitle("dE/dx a.u.");
    fQCOutputList->Add(fQC_Tracks_TPC_P);
    
    fQC_Tracks_TPC_PT               = new TH3F("fQC_Tracks_TPC_PT",          "Acc. Tracks TPC Sig.",                            2, -1., 1.,  kMOM_bin, kMOM_min, kMOM_max, 100, 0., 500.);
    fQC_Tracks_TPC_PT               ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_Tracks_TPC_PT               ->  GetYaxis()  ->  SetTitle("Track Transverse Momentum (GeV/c)");
    fQC_Tracks_TPC_PT               ->  GetZaxis()  ->  SetTitle("dE/dx a.u.");
    fQCOutputList->Add(fQC_Tracks_TPC_PT);
    
    fQC_Tracks_TPC_CLS              = new TH1F("fQC_Tracks_TPC_CLS",        "Acc. Tracks TPC Clusters",                         300, 0., 300. );
    fQC_Tracks_TPC_CLS              ->  GetXaxis()  ->  SetTitle("TPC Clusters");
    fQC_Tracks_TPC_CLS              ->  GetYaxis()  ->  SetTitle("Tracks");
    fQCOutputList->Add(fQC_Tracks_TPC_CLS);
    
    fQC_Tracks_TPC_FRC              = new TH1F("fQC_Tracks_TPC_FRC",        "Acc. Tracks TPC Clusters over Findable",           1000, 0., 1. );
    fQC_Tracks_TPC_FRC              ->  GetXaxis()  ->  SetTitle("TPC Clusters over Findable");
    fQC_Tracks_TPC_FRC              ->  GetYaxis()  ->  SetTitle("Tracks");
    fQCOutputList->Add(fQC_Tracks_TPC_FRC);
    
    fQC_Tracks_TPC_CHI              = new TH1F("fQC_Tracks_TPC_CHI",        "Acc. Tracks TPC Chi^2 over Clusters",              1000, 0., 10. );
    fQC_Tracks_TPC_CHI              ->  GetXaxis()  ->  SetTitle("TPC Chi^2 over Clusters");
    fQC_Tracks_TPC_CHI              ->  GetYaxis()  ->  SetTitle("Tracks");
    fQCOutputList->Add(fQC_Tracks_TPC_CHI);
    
    fQC_Tracks_TPC_CNS              = new TH1F("fQC_Tracks_TPC_CNS",        "Acc. Tracks TPC Chi^2 Constrained Global",         500, 0., 50. );
    fQC_Tracks_TPC_CNS              ->  GetXaxis()  ->  SetTitle("TPC Clusters");
    fQC_Tracks_TPC_CNS              ->  GetYaxis()  ->  SetTitle("Tracks");
    fQCOutputList->Add(fQC_Tracks_TPC_CNS);
    
    fQC_Tracks_ITS_CHI              = new TH1F("fQC_Tracks_ITS_CHI",        "Acc. Tracks ITS Chi^2 over Clusters",              500, 0., 50. );
    fQC_Tracks_ITS_CHI              ->  GetXaxis()  ->  SetTitle("ITS Chi^2 over Clusters");
    fQC_Tracks_ITS_CHI              ->  GetYaxis()  ->  SetTitle("Tracks");
    fQCOutputList->Add(fQC_Tracks_ITS_CHI);
    
    //_____________________________________________________________________________
    //
    //  KAONS
    //_____________________________________________________________________________
    fQC_Kaons_Momentum              = new TH2F("fQC_Kaons_Momentum",        "Acc. Kaons Momentum",                              2, -1., 1., kMOM_bin, kMOM_min, kMOM_max);
    fQC_Kaons_Momentum              ->  GetXaxis()  ->  SetTitle("Kaon Sign");
    fQC_Kaons_Momentum              ->  GetYaxis()  ->  SetTitle("Kaon Momentum (GeV/c)");
    fQC_Kaons_Momentum              ->  GetZaxis()  ->  SetTitle("Accepted Tracks");
    fQCOutputList->Add(fQC_Kaons_Momentum);
    
    fQC_Kaons_TMomentum             = new TH2F("fQC_Kaons_TMomentum",       "Acc. Kaons Transverse Momentum",                   2, -1., 1., kMOM_bin, kMOM_min, kMOM_max);
    fQC_Kaons_TMomentum             ->  GetXaxis()  ->  SetTitle("Kaon Sign");
    fQC_Kaons_TMomentum             ->  GetYaxis()  ->  SetTitle("Kaon Transverse Momentum (GeV/c)");
    fQC_Kaons_TMomentum             ->  GetZaxis()  ->  SetTitle("Accepted Tracks");
    fQCOutputList->Add(fQC_Kaons_TMomentum);
    
    fQC_Kaons_Eta                   = new TH2F("fQC_Kaons_Eta",             "Acc. Kaons #eta",                                  2, -1., 1., kETA_bin, kETA_min, kETA_max);
    fQC_Kaons_Eta                   ->  GetXaxis()  ->  SetTitle("Kaon Sign");
    fQC_Kaons_Eta                   ->  GetYaxis()  ->  SetTitle("Kaon #eta");
    fQC_Kaons_Eta                   ->  GetZaxis()  ->  SetTitle("Accepted Tracks");
    fQCOutputList->Add(fQC_Kaons_Eta);
    
    fQC_Kaons_Phi                   = new TH2F("fQC_Kaons_Phi",             "Acc. Kaons #phi",                                  2, -1., 1., kPhA_bin, kPhA_min, kPhA_max);
    fQC_Kaons_Phi                   ->  GetXaxis()  ->  SetTitle("Kaon Sign");
    fQC_Kaons_Phi                   ->  GetYaxis()  ->  SetTitle("Kaon #varphi");
    fQC_Kaons_Phi                   ->  GetZaxis()  ->  SetTitle("Accepted Tracks");
    fQCOutputList->Add(fQC_Kaons_Phi);
    
    fQC_Kaons_V0M                   = new TH2F("fQC_Kaons_V0M",            "Acc. Kaons in Mult.",                               2, -1., 1., kV0M_bin, kV0M_min, kV0M_max);
    fQC_Kaons_V0M                   ->  GetXaxis()  ->  SetTitle("Kaon Sign");
    fQC_Kaons_V0M                   ->  GetYaxis()  ->  SetTitle("V0M Multiplicity");
    fQC_Kaons_V0M                   ->  GetZaxis()  ->  SetTitle("Accepted Tracks");
    fQCOutputList->Add(fQC_Kaons_V0M);
    
    fQC_Kaons_TRK                   = new TH2F("fQC_Kaons_TRK",            "Acc. Kaons in Mult.",                               2, -1., 1., kTRK_bin, kTRK_min, kTRK_max);
    fQC_Kaons_TRK                   ->  GetXaxis()  ->  SetTitle("Kaon Sign");
    fQC_Kaons_TRK                   ->  GetYaxis()  ->  SetTitle("Tracklets Multiplicity");
    fQC_Kaons_TRK                   ->  GetZaxis()  ->  SetTitle("Accepted Tracks");
    fQCOutputList->Add(fQC_Kaons_TRK);
    
    fQC_Kaons_DCAXY_P               = new TH3F("fQC_Kaons_DCAXY_P",        "Acc. Kaons XY-DCA",                                 2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, 1000, -5., 5.);
    fQC_Kaons_DCAXY_P               ->  GetXaxis()  ->  SetTitle("Kaon Sign");
    fQC_Kaons_DCAXY_P               ->  GetYaxis()  ->  SetTitle("Kaon Momentum (GeV/c)");
    fQC_Kaons_DCAXY_P               ->  GetZaxis()  ->  SetTitle("XY-DCA (cm)");
    fQCOutputList->Add(fQC_Kaons_DCAXY_P);
    
    fQC_Kaons_DCAXY_PT              = new TH3F("fQC_Kaons_DCAXY_PT",       "Acc. Kaons XY-DCA",                                 2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, 1000, -5., 5.);
    fQC_Kaons_DCAXY_PT              ->  GetXaxis()  ->  SetTitle("Kaon Sign");
    fQC_Kaons_DCAXY_PT              ->  GetYaxis()  ->  SetTitle("Kaon Transverse Momentum (GeV/c)");
    fQC_Kaons_DCAXY_PT              ->  GetZaxis()  ->  SetTitle("XY-DCA (cm)");
    fQCOutputList->Add(fQC_Kaons_DCAXY_PT);
    
    fQC_Kaons_DCAZ_P                = new TH3F("fQC_Kaons_DCAZ_P",          "Acc. Kaons Z-DCA",                                 2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, 1000, -5., 5.);
    fQC_Kaons_DCAZ_P                ->  GetXaxis()  ->  SetTitle("Kaon Sign");
    fQC_Kaons_DCAZ_P                ->  GetYaxis()  ->  SetTitle("Kaon Momentum (GeV/c)");
    fQC_Kaons_DCAZ_P                ->  GetZaxis()  ->  SetTitle("Z-DCA (cm)");
    fQCOutputList->Add(fQC_Kaons_DCAZ_P);
    
    fQC_Kaons_DCAZ_PT               = new TH3F("fQC_Kaons_DCAZ_P",          "Acc. Kaons Z-DCA",                                 2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, 1000, -5., 5.);
    fQC_Kaons_DCAZ_PT               ->  GetXaxis()  ->  SetTitle("Kaon Sign");
    fQC_Kaons_DCAZ_PT               ->  GetYaxis()  ->  SetTitle("Kaon Transverse Momentum (GeV/c)");
    fQC_Kaons_DCAZ_PT               ->  GetZaxis()  ->  SetTitle("Z-DCA (cm)");
    fQCOutputList->Add(fQC_Kaons_DCAZ_PT);
    
    fQC_Kaons_TOF_P                 = new TH3F("fQC_Kaons_TOF_P",           "Acc. Kaons TOF Sig.",                              2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, 120, 0., 1.2);
    fQC_Kaons_TOF_P                 ->  GetXaxis()  ->  SetTitle("Kaon Sign");
    fQC_Kaons_TOF_P                 ->  GetYaxis()  ->  SetTitle("Kaon Momentum (GeV/c)");
    fQC_Kaons_TOF_P                 ->  GetZaxis()  ->  SetTitle("#beta");
    fQCOutputList->Add(fQC_Kaons_TOF_P);
    
    fQC_Kaons_TOF_PT                = new TH3F("fQC_Kaons_TOF_PT",          "Acc. Kaons TOF Sig.",                              2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, 120, 0., 1.2);
    fQC_Kaons_TOF_PT                ->  GetXaxis()  ->  SetTitle("Kaon Sign");
    fQC_Kaons_TOF_PT                ->  GetYaxis()  ->  SetTitle("Kaon Transverse Momentum (GeV/c)");
    fQC_Kaons_TOF_PT                ->  GetZaxis()  ->  SetTitle("#beta");
    fQCOutputList->Add(fQC_Kaons_TOF_PT);
    
    fQC_Kaons_TPC_P                 = new TH3F("fQC_Kaons_TPC_P",           "Acc. Kaons TPC Sig.",                              2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, 100, 0., 500.);
    fQC_Kaons_TPC_P                 ->  GetXaxis()  ->  SetTitle("Kaon Sign");
    fQC_Kaons_TPC_P                 ->  GetYaxis()  ->  SetTitle("Kaon Momentum (GeV/c)");
    fQC_Kaons_TPC_P                 ->  GetZaxis()  ->  SetTitle("dE/dx a.u.");
    fQCOutputList->Add(fQC_Kaons_TPC_P);
    
    fQC_Kaons_TPC_PT                = new TH3F("fQC_Kaons_TPC_PT",          "Acc. Kaons TPC Sig.",                              2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, 100, 0., 500.);
    fQC_Kaons_TPC_PT                ->  GetXaxis()  ->  SetTitle("Kaon Sign");
    fQC_Kaons_TPC_PT                ->  GetYaxis()  ->  SetTitle("Kaon Transverse Momentum (GeV/c)");
    fQC_Kaons_TPC_PT                ->  GetZaxis()  ->  SetTitle("dE/dx a.u.");
    fQCOutputList->Add(fQC_Kaons_TPC_PT);
    
    //_____________________________________________________________________________
    //
    //  PID
    //_____________________________________________________________________________
    fQC_PID_TOF_Kaons_P             = new TH3F("fQC_PID_TOF_Kaons_P",   "Tracks TOF Kaons nSigma",                              2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, sPID_bin, sPID_min, sPID_max);
    fQC_PID_TOF_Kaons_P             ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_PID_TOF_Kaons_P             ->  GetYaxis()  ->  SetTitle("Track Momentum (GeV/c)");
    fQC_PID_TOF_Kaons_P             ->  GetZaxis()  ->  SetTitle("n#sigma_{kaons}^{TOF}");
    fQCOutputList->Add(fQC_PID_TOF_Kaons_P);
        
    fQC_PID_TOF_Kaons_PT            = new TH3F("fQC_PID_TOF_Kaons_PT",  "Tracks TOF Kaons nSigma",                              2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, sPID_bin, sPID_min, sPID_max);
    fQC_PID_TOF_Kaons_PT            ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_PID_TOF_Kaons_PT            ->  GetYaxis()  ->  SetTitle("Track Transverse Momentum (GeV/c)");
    fQC_PID_TOF_Kaons_PT            ->  GetZaxis()  ->  SetTitle("n#sigma_{kaons}^{TOF}");
    fQCOutputList->Add(fQC_PID_TOF_Kaons_PT);
    
    fQC_PID_TPC_Kaons_P             = new TH3F("fQC_PID_TPC_Kaons_P",   "Tracks TPC Kaons nSigma",                              2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, sPID_bin, sPID_min, sPID_max);
    fQC_PID_TPC_Kaons_P             ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_PID_TPC_Kaons_P             ->  GetYaxis()  ->  SetTitle("Track Momentum (GeV/c)");
    fQC_PID_TPC_Kaons_P             ->  GetZaxis()  ->  SetTitle("n#sigma_{kaons}^{TPC}");
    fQCOutputList->Add(fQC_PID_TPC_Kaons_P);
    
    fQC_PID_TPC_Kaons_PT            = new TH3F("fQC_PID_TPC_Kaons_PT",  "Tracks TPC Kaons nSigma",                              2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, sPID_bin, sPID_min, sPID_max);
    fQC_PID_TPC_Kaons_PT            ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_PID_TPC_Kaons_PT            ->  GetYaxis()  ->  SetTitle("Track Transverse Momentum (GeV/c)");
    fQC_PID_TPC_Kaons_PT            ->  GetZaxis()  ->  SetTitle("n#sigma_{kaons}^{TPC}");
    fQCOutputList->Add(fQC_PID_TPC_Kaons_PT);
    
    fQC_PID_TPC_TOF_Kaons_PT        = new TH3F("fQC_PID_TPC_TOF_Kaons_PT",  "Tracks TPC Kaons nSigma",                          2, -1., 1., sPID_bin, sPID_min, sPID_max, sPID_bin, sPID_min, sPID_max);
    fQC_PID_TPC_TOF_Kaons_PT        ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_PID_TPC_TOF_Kaons_PT        ->  GetYaxis()  ->  SetTitle("n#sigma_{kaons}^{TOF}");
    fQC_PID_TPC_TOF_Kaons_PT        ->  GetZaxis()  ->  SetTitle("n#sigma_{kaons}^{TPC}");
    fQCOutputList->Add(fQC_PID_TPC_TOF_Kaons_PT);
    
    // KAONS _____________________________________________________________________________
    fQC_PID_TOF_NSig_SEL_Kaons_P    = new TH3F("fQC_PID_TOF_NSig_SEL_Kaons_P",   "Selected tracks TOF Kaons nSigma",                     2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, sPID_bin, sPID_min, sPID_max);
    fQC_PID_TOF_NSig_SEL_Kaons_P    ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_PID_TOF_NSig_SEL_Kaons_P    ->  GetYaxis()  ->  SetTitle("Track Momentum (GeV/c)");
    fQC_PID_TOF_NSig_SEL_Kaons_P    ->  GetZaxis()  ->  SetTitle("n#sigma_{kaons}^{TOF}");
    fQCOutputList->Add(fQC_PID_TOF_NSig_SEL_Kaons_P);
    
    fQC_PID_TOF_NSig_SEL_Kaons_PT   = new TH3F("fQC_PID_TOF_NSig_SEL_Kaons_PT",  "Selected tracks TOF Kaons nSigma",                     2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, sPID_bin, sPID_min, sPID_max);
    fQC_PID_TOF_NSig_SEL_Kaons_PT   ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_PID_TOF_NSig_SEL_Kaons_PT   ->  GetYaxis()  ->  SetTitle("Track Transverse Momentum (GeV/c)");
    fQC_PID_TOF_NSig_SEL_Kaons_PT   ->  GetZaxis()  ->  SetTitle("n#sigma_{kaons}^{TOF}");
    fQCOutputList->Add(fQC_PID_TOF_NSig_SEL_Kaons_PT);
    
    fQC_PID_TPC_NSig_SEL_Kaons_P    = new TH3F("fQC_PID_TPC_NSig_SEL_Kaons_P",   "Selected tracks TPC Kaons nSigma",                     2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, sPID_bin, sPID_min, sPID_max);
    fQC_PID_TPC_NSig_SEL_Kaons_P    ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_PID_TPC_NSig_SEL_Kaons_P    ->  GetYaxis()  ->  SetTitle("Track Momentum (GeV/c)");
    fQC_PID_TPC_NSig_SEL_Kaons_P    ->  GetZaxis()  ->  SetTitle("n#sigma_{kaons}^{TPC}");
    fQCOutputList->Add(fQC_PID_TPC_NSig_SEL_Kaons_P);
    
    fQC_PID_TPC_NSig_SEL_Kaons_PT   = new TH3F("fQC_PID_TPC_NSig_SEL_Kaons_PT",  "Selected tracks TPC Kaons nSigma",                     2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, sPID_bin, sPID_min, sPID_max);
    fQC_PID_TPC_NSig_SEL_Kaons_PT   ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_PID_TPC_NSig_SEL_Kaons_PT   ->  GetYaxis()  ->  SetTitle("Track Transverse Momentum (GeV/c)");
    fQC_PID_TPC_NSig_SEL_Kaons_PT   ->  GetZaxis()  ->  SetTitle("n#sigma_{kaons}^{TPC}");
    fQCOutputList->Add(fQC_PID_TPC_NSig_SEL_Kaons_PT);
    
    fQC_PID_TOF_Sgnl_SEL_Kaons_P    = new TH3F("fQC_PID_TOF_Sgnl_SEL_Kaons_P",   "Selected tracks TOF Kaons nSigma",                     2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, sPID_bin, sPID_min, sPID_max);
    fQC_PID_TOF_Sgnl_SEL_Kaons_P    ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_PID_TOF_Sgnl_SEL_Kaons_P    ->  GetYaxis()  ->  SetTitle("Track Momentum (GeV/c)");
    fQC_PID_TOF_Sgnl_SEL_Kaons_P    ->  GetZaxis()  ->  SetTitle("#beta_{TOF}");
    fQCOutputList->Add(fQC_PID_TOF_Sgnl_SEL_Kaons_P);
    
    fQC_PID_TOF_Sgnl_SEL_Kaons_PT   = new TH3F("fQC_PID_TOF_Sgnl_SEL_Kaons_PT",  "Selected tracks TOF Kaons nSigma",                     2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, sPID_bin, sPID_min, sPID_max);
    fQC_PID_TOF_Sgnl_SEL_Kaons_PT   ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_PID_TOF_Sgnl_SEL_Kaons_PT   ->  GetYaxis()  ->  SetTitle("Track Transverse Momentum (GeV/c)");
    fQC_PID_TOF_Sgnl_SEL_Kaons_PT   ->  GetZaxis()  ->  SetTitle("#beta_{TOF}");
    fQCOutputList->Add(fQC_PID_TOF_Sgnl_SEL_Kaons_PT);
    
    fQC_PID_TPC_Sgnl_SEL_Kaons_P    = new TH3F("fQC_PID_TPC_Sgnl_SEL_Kaons_P",   "Selected tracks TPC Kaons nSigma",                     2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, sPID_bin, sPID_min, sPID_max);
    fQC_PID_TPC_Sgnl_SEL_Kaons_P    ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_PID_TPC_Sgnl_SEL_Kaons_P    ->  GetYaxis()  ->  SetTitle("Track Momentum (GeV/c)");
    fQC_PID_TPC_Sgnl_SEL_Kaons_P    ->  GetZaxis()  ->  SetTitle("dE/dx a.u.");
    fQCOutputList->Add(fQC_PID_TPC_Sgnl_SEL_Kaons_P);
    
    fQC_PID_TPC_Sgnl_SEL_Kaons_PT   = new TH3F("fQC_PID_TPC_Sgnl_SEL_Kaons_PT",  "Selected tracks TPC Kaons nSigma",                     2, -1., 1., kMOM_bin, kMOM_min, kMOM_max, sPID_bin, sPID_min, sPID_max);
    fQC_PID_TPC_Sgnl_SEL_Kaons_PT   ->  GetXaxis()  ->  SetTitle("Track Sign");
    fQC_PID_TPC_Sgnl_SEL_Kaons_PT   ->  GetYaxis()  ->  SetTitle("Track Transverse Momentum (GeV/c)");
    fQC_PID_TPC_Sgnl_SEL_Kaons_PT   ->  GetZaxis()  ->  SetTitle("dE/dx a.u.");
    fQCOutputList->Add(fQC_PID_TPC_Sgnl_SEL_Kaons_PT);
    
    //_____________________________________________________________________________
    //
    //  GENERAL
    //_____________________________________________________________________________
    fQC_Phi_InvMass_Rec             = new TH2F("fQC_Phi_InvMass_Rec",             "Invariant Mass Reconstructed #phi",          1000, 0.99, 1.08, kMOM_bin, kMOM_min, kMOM_max);
    fQC_Phi_InvMass_Rec             ->  GetXaxis()  ->  SetTitle("M_{K^{+}K^{-}}");
    fQC_Phi_InvMass_Rec             ->  GetYaxis()  ->  SetTitle("Track Transverse Momentum (GeV/c)");
    fQC_Phi_InvMass_Rec             ->  GetZaxis()  ->  SetTitle("Counts");
    fQCOutputList->Add(fQC_Phi_InvMass_Rec);
    
    fQC_Phi_InvMass_Gen             = new TH2F("fQC_Phi_InvMass_Gen",             "Invariant Mass Generated #phi",              1000, 0.99, 1.08, kMOM_bin, kMOM_min, kMOM_max);
    fQC_Phi_InvMass_Gen             ->  GetXaxis()  ->  SetTitle("M_{K^{+}K^{-}}");
    fQC_Phi_InvMass_Gen             ->  GetYaxis()  ->  SetTitle("Track Transverse Momentum (GeV/c)");
    fQC_Phi_InvMass_Gen             ->  GetZaxis()  ->  SetTitle("Counts");
    fQCOutputList->Add(fQC_Phi_InvMass_Gen);
    
    fQC_Phi_InvMass_Eff             = new TH2F("fQC_Phi_InvMass_Eff",             "Invariant Mass Efficiency #phi",             1000, 0.99, 1.08, kMOM_bin, kMOM_min, kMOM_max);
    fQC_Phi_InvMass_Eff             ->  GetXaxis()  ->  SetTitle("M_{K^{+}K^{-}}");
    fQC_Phi_InvMass_Eff             ->  GetYaxis()  ->  SetTitle("Track Transverse Momentum (GeV/c)");
    fQC_Phi_InvMass_Eff             ->  GetZaxis()  ->  SetTitle("Counts");
    fQCOutputList->Add(fQC_Phi_InvMass_Eff);
    
    // TODO: ADD VETO
    
    PostData(2, fQCOutputList);
    
    // Where i is the slot that writes to a tree
    OpenFile(3);
    
    // PhiCandidate Tree Set-Up
    fPhiCandidate = new TTree   (Form("PhiCandidate_%s",fRunName.Data()),    "Data Tree for Phi Candidates");
    fPhiCandidate->Branch       ("EventMask",       &fEventMask,        "fEventMask/b");
    fPhiCandidate->Branch       ("Multiplicity",    &fCurrent_V0M,      "fMultiplicity/F");
    fPhiCandidate->Branch       ("Spherocity",      &fCurrent_SPH,      "fSpherocity/F");
    fPhiCandidate->Branch       ("RTransverse",     &fCurrent_RT,       "fRTransverse/F");
    fPhiCandidate->Branch       ("nPhi",            &fnPhi,             "fnPhi/I");
    fPhiCandidate->Branch       ("Px",              &fPhiPx,            "fPhiPx[fnPhi]/F");
    fPhiCandidate->Branch       ("Py",              &fPhiPy,            "fPhiPy[fnPhi]/F");
    fPhiCandidate->Branch       ("Pz",              &fPhiPz,            "fPhiPz[fnPhi]/F");
    fPhiCandidate->Branch       ("InvMass",         &fInvMass,          "fInvMass[fnPhi]/F");
    fPhiCandidate->Branch       ("iKaon",           &fiKaon,            "fiKaon[fnPhi]/I");
    fPhiCandidate->Branch       ("jKaon",           &fjKaon,            "fjKaon[fnPhi]/I");
    if ( kMCbool )  fPhiCandidate->Branch   ("TrueInvMass",          &fTrueInvMass,           "fTrueInvMass[fnPhi]/F");
    
    if ( kPhibool )                 PostData(3, fPhiCandidate);
    
    // KaonCandidate Tree Set-Up
    fKaonCandidate = new TTree (Form("KaonCandidate_%s",fRunName.Data()),    "Data Tree for Kaon Candidates");
    fKaonCandidate->Branch     ("EventMask",        &fEventMask,        "fEventMask/b");
    fKaonCandidate->Branch     ("Multiplicity",     &fCurrent_V0M,      "fMultiplicity/F");
    fKaonCandidate->Branch     ("nKaon",            &fnKaon,            "fnKaon/b");
    fKaonCandidate->Branch     ("Px",               &fKaonPx,           "fKaonPx[fnKaon]/F");
    fKaonCandidate->Branch     ("Py",               &fKaonPy,           "fKaonPy[fnKaon]/F");
    fKaonCandidate->Branch     ("Pz",               &fKaonPz,           "fKaonPz[fnKaon]/F");
    fKaonCandidate->Branch     ("Charge",           &fCharge,           "fCharge[fnKaon]/B");
    fKaonCandidate->Branch     ("TOFSigma",         &fTOFSigma,         "fTOFSigma[fnKaon]/B");
    fKaonCandidate->Branch     ("TPCSigma",         &fTPCSigma,         "fTPCSigma[fnKaon]/B");
    
    if ( kKaonbool )                PostData(4, fKaonCandidate);

    fPhiEfficiency = new TTree  (Form("PhiEfficiency_%s",fRunName.Data()),   "MC Tree for Phi Efficiency");
    fPhiEfficiency->Branch      ("EventMask",       &fEventMask,        "fEventMask/b");
    fPhiEfficiency->Branch      ("TrueEventMask",   &fTrueEventMask,    "fTrueEventMask/b");
    fPhiEfficiency->Branch      ("Multiplicity",    &fCurrent_V0M,      "fMultiplicity/F");
    fPhiEfficiency->Branch      ("Spherocity",      &fCurrent_SPH,      "fSpherocity/F");
    fPhiEfficiency->Branch      ("nPhi",            &fnPhiTru,          "fnPhiTru/I");
    fPhiEfficiency->Branch      ("Px",              &fPhiTruPx,         "fPhiTruPx[fnPhiTru]/F");
    fPhiEfficiency->Branch      ("Py",              &fPhiTruPy,         "fPhiTruPy[fnPhiTru]/F");
    fPhiEfficiency->Branch      ("Pz",              &fPhiTruPz,         "fPhiTruPz[fnPhiTru]/F");
    fPhiEfficiency->Branch      ("Selection",       &fSelection,        "fSelection[fnPhiTru]/b");
    
    if ( kPhibool   &&  kMCbool )   PostData(5, fPhiEfficiency);
    
    fKaonEfficiency = new TTree (Form("KaonEfficiency_%s",fRunName.Data()),  "MC Tree for Kaon Efficiency");
    fKaonEfficiency->Branch     ("EventMask",       &fEventMask,        "fEventMask/b");
    fKaonEfficiency->Branch     ("TrueEventMask",   &fTrueEventMask,    "fTrueEventMask/b");
    fKaonEfficiency->Branch     ("Multiplicity",    &fCurrent_V0M,     "fMultiplicity/F");
    
    if ( kKaonbool  &&  kMCbool )   PostData(6, fKaonEfficiency);
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::UserExec( Option_t* )                      {
    //  Check the Event is available and within requirements
    if ( !fIsEventCandidate() )    return;
    //
    //  Check the Event type
    fIsEventMultiplicityAvailable();
    fIsEventPileUp();
    uCalculateSpherocity();
    uCalculateRT();
    uGuessCollisionSystem();
    //
    // Setting utility variables
    Int_t           nTrack(fAOD->GetNumberOfTracks());
    TLorentzVector  fKaon1, fKaon2, fPhi;
    //
    // Looping over tracks
    for ( Int_t iTrack(0); iTrack < nTrack; iTrack++ )
    {
        // Recovering Track
        fCurrent_Track  =   static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
        
        // Check the track has due requirements
        if ( !fIsTrackCandidate() ) continue;
        
        // Fill the QC for Event Tracks
        fQC_TRK             ();
        
        // Fill the QC for PID
        fQC_PID             ();
        
        // Check the PID is present and within requirements
        if ( !fIsKaonCandidate() ) continue;
        
        // Fill the QC for Kaons Event Tracks
        fQC_TRK_Kaons       ();
        
        // Fill the QC for Kaons PID
        fQC_PID_Sel_Kaons   ();

        // Filling the Kaon Tree
        fKaonPx[fnKaon] =   fCurrent_Track->Px();
        fKaonPy[fnKaon] =   fCurrent_Track->Py();
        fKaonPz[fnKaon] =   fCurrent_Track->Pz();
        fCharge[fnKaon] =   static_cast<Char_t>(fCurrent_Track->GetSign());
        fKaonLabels[fnKaon] =   fCurrent_Track->GetLabel();
        fnKaon++;
    }
    //Coupling Kaons
    for ( Int_t iKaon(0); iKaon < fnKaon; iKaon++)
    {
        // Storing first Kaon kinematics and sign
        fKaon1.SetXYZM(fKaonPx[iKaon],fKaonPy[iKaon],fKaonPz[iKaon],.493677);
        for ( Int_t jKaon(iKaon+1); jKaon < fnKaon; jKaon++)
        {
            if ( fCharge[iKaon] == fCharge[jKaon] ) continue;
            
            // Storing second Kaon Kinematics and combining the two
            fKaon2.SetXYZM(fKaonPx[jKaon],fKaonPy[jKaon],fKaonPz[jKaon],.493677);
            
            // Check the Phi is a good candidate
            fPhi = (fKaon1 + fKaon2);
            if ( !fIsPhiCandidate(fPhi) ) continue;
            
            fPhiPx[fnPhi]       =   fPhi.Px();
            fPhiPy[fnPhi]       =   fPhi.Py();
            fPhiPz[fnPhi]       =   fPhi.Pz();
            fInvMass[fnPhi]     =   (fPhi).Mag();
            fiKaon[fnPhi]       =   iKaon;
            fjKaon[fnPhi]       =   jKaon;
            fTrueInvMass[fnPhi] =   0;
            if ( kMCbool )  {
                auto    fTrue_iKaon =   static_cast<AliAODMCParticle*>(AODMCTrackArray->At(fKaonLabels[iKaon]));
                auto    fTrue_jKaon =   static_cast<AliAODMCParticle*>(AODMCTrackArray->At(fKaonLabels[jKaon]));
                if ( fIsCandidateTruPhi(fTrue_iKaon,fTrue_jKaon) ) {
                    TLorentzVector  fTrue_iKaon_Vector, fTrue_jKaon_Vector;
                    fTrue_iKaon->Momentum(fTrue_iKaon_Vector);
                    fTrue_jKaon->Momentum(fTrue_jKaon_Vector);
                    auto    fTrue_Phi_Vector    =   fTrue_iKaon_Vector + fTrue_jKaon_Vector;
                    fTrueInvMass[fnPhi]         =   fTrue_Phi_Vector.Mag();
                    //
                    fQC_Phi_InvMass_Gen ->  Fill((fTrue_Phi_Vector).Mag(),(fTrue_Phi_Vector).Pt());
                    fQC_Phi_InvMass_Rec ->  Fill((fPhi).Mag(),(fPhi).Pt());
                    //
                    fPhiRecParticles[fnPhiRec]  =   static_cast<AliAODMCParticle*>(AODMCTrackArray->At(fTrue_iKaon->GetMother()));
                    fnPhiRec++;
                }
            }
            fnPhi++;
        }
    }
    // Loop over all primary MC particle
    if ( kMCbool )
    {
        Int_t           nTrack(AODMCTrackArray->GetEntriesFast());
        for ( Int_t iTrack(0); iTrack < nTrack; iTrack++ )
        {
            AliAODMCParticle* fPhiTru = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(iTrack));
            
            if ( !fIsPhi( fPhiTru ) ) continue;
            
            // Kinematics
            fPhiTruPx[fnPhiTru]         =   fPhiTru->Px();
            fPhiTruPy[fnPhiTru]         =   fPhiTru->Py();
            fPhiTruPz[fnPhiTru]         =   fPhiTru->Pz();
            
            // Options
            fSelection[fnPhiTru]        =   0;
            if( fIsPhiGen(fPhiTru) )    fSelection[fnPhiTru] = 1;
            if( fIsPhiRec(fPhiTru) )    fSelection[fnPhiTru] = 2;
            
            fnPhiTru++;
        }
    }
    // Saving output
    fFillTrees();
    fPostData();
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fSetEventCountLabels ( TH1D * fEvCount )   {
    fEvCount->GetXaxis()->SetBinLabel(1,"ALL");
    fEvCount->GetXaxis()->SetBinLabel(2,"fAOD-fMCD");
    fEvCount->GetXaxis()->SetBinLabel(3,"TrackArray");
    fEvCount->GetXaxis()->SetBinLabel(4,"NoTrigger");
    fEvCount->GetXaxis()->SetBinLabel(5,"PIDResponse");
    fEvCount->GetXaxis()->SetBinLabel(6,"IncompleteDAQ");
    fEvCount->GetXaxis()->SetBinLabel(7,"NoSPDVTX");
    fEvCount->GetXaxis()->SetBinLabel(8,"TRK-SPD Mismatch");
    fEvCount->GetXaxis()->SetBinLabel(9,"VTX<Cut");
    fEvCount->GetXaxis()->SetBinLabel(10,"Accepted");
    fEvCount->GetXaxis()->SetBinLabel(11,"HasMult");
    fEvCount->GetXaxis()->SetBinLabel(12,"Pile-Up");
    fEvCount->GetXaxis()->SetBinLabel(13,"Pile-Up in Mult");
    fEvCount->GetXaxis()->SetBinLabel(14,"NoPhiCand");
    fEvCount->GetXaxis()->SetBinLabel(15,"OFPhiCand");
    fEvCount->GetXaxis()->SetBinLabel(16,"NoPhTCand");
    fEvCount->GetXaxis()->SetBinLabel(17,"OFPhTCand");
    fEvCount->GetXaxis()->SetBinLabel(18,"NoKaonCand");
    fEvCount->GetXaxis()->SetBinLabel(19,"OFKaonCand");
    fEvCount->GetXaxis()->SetBinLabel(20,"NoKaoTCand");
    fEvCount->GetXaxis()->SetBinLabel(21,"OFKaoTCand");
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsEventCandidate ()                       {
    //
    //  Counting the Triggered events
    //
    fQC_Event_Enum_FLL->Fill("ALL",1);
    //
    //>-> Starting Mandatory Checks
    //>->____________________________________________________________________//
    //
    // Recovering AOD(MC) Event
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    fMCD = dynamic_cast<AliMCEvent*>(MCEvent());
    //
    // Check the event is there
    if ( !fAOD || (!fMCD && kMCbool) )  {
        fFillEventEnumerate("fAOD-fMCD");
        fPostData();
        return false;
    }
    //
    // Recover and Check the MC tracks
    if ( kMCbool )      AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if ( !AODMCTrackArray && kMCbool )  {
        fFillEventEnumerate("TrackArray");
        fPostData();
        return false;
    }
    //
    if ( ( ( kTriggerMask & ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()) == 0 ) )  {
        fFillEventEnumerate("NoTrigger");
        fStoreTruePhi(kTRU_NOTRG);
        fPostData();
        return false;
    }
    //
    // Define and Fetch PID with Manager
    AliAnalysisManager* fAnalysisManager = AliAnalysisManager::GetAnalysisManager();
    if ( fAnalysisManager )    {
        AliInputEventHandler* fInputHandler = (AliInputEventHandler*)(fAnalysisManager->GetInputEventHandler());
        if (fInputHandler)   fPIDResponse = fInputHandler->GetPIDResponse();
    }
    //
    // Check the PID is present
    if ( !fPIDResponse )    {
        fFillEventEnumerate("PIDResponse");
        fPostData();
        return false;
    }
    //
    // Check the DAQ is present and working
    if ( fAOD->IsIncompleteDAQ() )    {
        fFillEventEnumerate("IncompleteDAQ");
        fPostData();
        return false;
    }
    //
    //  Setting to default Event Variables
    fEventMask      =   0;
    fTrueEventMask  =   0;
    fCurrent_V0M    =   150.;
    fCurrent_Run    =   fAOD->GetRunNumber();
    fnPhi           =   0;
    fnPhiTru        =   0;
    fnPhiRec        =   0;
    fnKaon          =   0;
    fnPhiTru        =   0;
    fnKaonTru       =   0;
    //
    //<-<____________________________________________________________________//
    //
    //>-> Vertex Analysis Cuts
    //>->____________________________________________________________________//
    //
    // Check the Primary Vertex is there and is correct
    if ( !fIsPrimaryVertexCandidate() )   return false;
    //
    //<-<____________________________________________________________________//
    //
    //  Event Counting
    fFillEventEnumerate("Accepted");
    //
    return true;
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsPrimaryVertexCandidate ( )              {
    // Recovering Primary Vertex from General methods and SPD
    auto    PrimaryVertexSPD    = fAOD->GetPrimaryVertexSPD();
    auto    PrimaryVertexTRK    = fAOD->GetPrimaryVertex();
            fPrimaryVertex      = PrimaryVertexTRK;
    
    // Requires the vertex is reconstructed by the SPD
    if ( !PrimaryVertexSPD  ||  PrimaryVertexSPD->GetNContributors() < 1 )
    {
        fFillEventEnumerate("NoSPDVTX");
        fStoreTruePhi(kTRU_NOSPDVTX);
        fPostData();
        return false;
    }

    // In lack of the general Method reconstructed Vertex, take the SPD reconstruction
    if ( !PrimaryVertexTRK  ||  PrimaryVertexTRK->GetNContributors() < 1 )
    {
        fPrimaryVertex = PrimaryVertexSPD;
    }
    
    // If Both are present and reconstructed discard event if the two mismatch for more than 0.5cm
    else
    {
        auto VertexZSPD = PrimaryVertexSPD->GetZ();
        auto VertexZTRK = PrimaryVertexTRK->GetZ();
        if ( std::fabs(VertexZSPD-VertexZTRK) > 0.5 )
        {
            fFillEventEnumerate("TRK-SPD Mismatch");
            fStoreTruePhi(kTRU_TRKSPDMM);
            fPostData();
            return false;
        }
    }
    
    // Fill the Vertex Z position histogram
    fQC_Event_Vertex_Fll->Fill(fPrimaryVertex->GetZ());
    
    if ( std::fabs(fPrimaryVertex->GetZ()) > kVertexCut )
    {
        fFillEventEnumerate("VTX<Cut");
        fStoreTruePhi(kTRU_VTXCUT);
        fPostData();
        return false;
    }
    
    // Fill the Vertex Z position histogram
    fQC_Event_Vertex_Cut->Fill(fPrimaryVertex->GetZ());
    
    fPostData();
    return  true;
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsEventMultiplicityAvailable ()           {
    // Recovering Multiplicity information
    AliMultSelection   *fMultSelectio2  = (AliMultSelection*) fAOD->FindListObject("MultSelection");
    if ( !fMultSelection )  fCurrent_V0M    =   102.;
    // if pPb should be V0A
                    fCurrent_V0M    =   fMultSelectio2->GetMultiplicityPercentile("V0M");
    if ( fIs_p_Pb ) fCurrent_V0M    =   fMultSelectio2->GetMultiplicityPercentile("V0A");
    fCurrent_TRK    =   (dynamic_cast<AliAODHeader*>(fAOD->GetHeader()))->GetRefMultiplicityComb08();
    //
    //  INELGT0 Check
    if ( fCurrent_TRK != -1 && fCurrent_TRK != -2 ) {
        for ( Int_t iTRK = 0; iTRK < fAOD->GetMultiplicity()->GetNumberOfTracklets(); iTRK++ )  {
            if ( TMath::Abs(fAOD->GetMultiplicity()->GetEta(iTRK)) < 1.0 )  {
                fSetEventMask(0);
                break;
            }
        }
    }
    //
    if ( !fCheckMask(0) )       fCurrent_V0M = 104;
    //
    if ( fCurrent_V0M == -200 ) fCurrent_V0M = 154;
    if ( fCurrent_V0M == -201 ) fCurrent_V0M = 156;
    if ( fCurrent_V0M == -202 ) fCurrent_V0M = 158;
    if ( fCurrent_V0M == -203 ) fCurrent_V0M = 160;
    if ( fCurrent_V0M == -204 ) fCurrent_V0M = 162;
    if ( ( fCurrent_V0M >= 0 ) && ( fCurrent_V0M <= 100 ) )   {
        if ( fAOD->IsPileupFromSPD() )              fCurrent_V0M += 300;
        else if ( fAOD->IsPileupFromSPDInMultBins() )    fCurrent_V0M += 400;
    }
    //
    // Fill the QC on Multiplicity
    fQC_Event_Enum_V0M->Fill(fCurrent_V0M);
    fQC_Event_Enum_TRK->Fill(fCurrent_TRK);
    fQC_Event_Enum_V0T->Fill(fCurrent_V0M,fCurrent_TRK);
    //
    return true;
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsEventPileUp()                           {
    bool fResult    =   false;
    if ( fAOD->IsPileupFromSPD() )  {
        fQC_Event_Enum_FLL->Fill("Pile-Up",1);
        fSetEventMask(1);
        fResult     =   true;
    }
    if ( fAOD->IsPileupFromSPDInMultBins() )    {
        fQC_Event_Enum_FLL->Fill("Pile-Up in Mult",1);
        fSetEventMask(2);
        fResult     =   true;
    }
    return fResult;
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsSPDClusterVsTrackletBG()                {
    bool fResult    =   false;
    return fResult;
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fCheckINELgt0( AliAODMCParticle* fCurrent_Particle )   {
    if ( fCheckTrueMask(kTRU_INELGT0) )   return true;
    if ( !fCurrent_Particle->IsPrimary() ) return false;
    if ( !(fCurrent_Particle->Charge()   !=  0) ) return false;
    if ( !(TMath::Abs(fCurrent_Particle->Eta()) <= 1) ) return false;
    fQC_Event_Enum_FLL->Fill(12);
    fSetTrueEventMask(kTRU_INELGT0);
    return true;
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fFillEventEnumerate ( Int_t iIndex )       {
    for ( Int_t iFill = 1; iFill <= iIndex; iFill++ )   {
        fQC_Event_Enum_FLL->Fill(iFill);
    }
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fFillEventEnumerate ( TString iBinName )   {
    for ( Int_t iBin = 1; iBin <= fQC_Event_Enum_FLL->GetNbinsX(); iBin++ )    {
        if ( strcmp(iBinName.Data(),fQC_Event_Enum_FLL->GetXaxis()->GetBinLabel(iBin)) == 0 )  {
            fFillEventEnumerate(iBin-1);
        }
    }
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsTrackCandidate ( )                      {
    // Check the track is there and has proper kinematics
    if ( !fCurrent_Track )                              return false;
    fCurrent_Track_Eta      = fCurrent_Track->Eta();
    fCurrent_Track_Charge   = fCurrent_Track->GetSign();
    fCurrent_Track_TransMom = fCurrent_Track->Pt();
    fCurrent_Track_Momentum = fCurrent_Track->P();
    fCurrent_Track_Phi      = fCurrent_Track->Phi()*360/(2*TMath::Pi());
    fCurrent_Track          ->GetImpactParameters(fCurrent_Track_DCAXY,fCurrent_Track_DCAZ);
    
    //  Kinematics Cuts
    if (  std::fabs(fCurrent_Track_TransMom) < 0.15 ||  std::fabs(fCurrent_Track_Eta) > 0.80 || std::fabs(fCurrent_Track_Charge) != 1 )  return false;
    
    //  If a Filterbit is set, use it
    if ( kFilterBit >= 0 ) return fCurrent_Track->TestFilterBit(BIT(kFilterBit));
    /*
    if ( kFilterBit >= 0 ) {
        if ( !fCurrent_Track->TestFilterBit(BIT(kFilterBit)) )   return false;
        else return true;
    }
     */

    //  ------------------------ CUSTOM IMPLEMENTATION OF FILTERBIT 5 COMMON FOR ALL SYSTEMS
    
    //  ------------------------ GENERAL
    
    //  -   //  KINK DAUGHTERS
    if ( fPrimaryVertex->GetType() == AliAODVertex::kKink )                                                     return false;
    
    //  -   //  DCA-XY
    if ( fabs(fCurrent_Track_DCAXY)   > kNSigmaPtDepXYDCA*(0.0026+0.0050/(pow(fCurrent_Track_TransMom,1.01))) ) return false;
    
    //  -   //  DCA-Z
    if ( fabs(fCurrent_Track_DCAZ)   > kDCAzCut )                                                               return false;
        
    //  ------------------------ TPC
    
    //  -   //  TPC CROSSED ROWS
    if ( fCurrent_Track->GetTPCCrossedRows() < kMinTPCclusters )                                                return false;
    
    //  -   //  TPC CROSSED OVER FINDABLE
    if ( ( ( fCurrent_Track->GetTPCCrossedRows() / fCurrent_Track->GetTPCNclsF() ) < kTPCClsOverFndbl ) && ( fCurrent_Track->GetTPCNclsF() != 0 ) ) return false;
    
    //  -   //  CHI2 PER CLUSTER
    if ( ( ( 1.*fCurrent_Track->GetTPCchi2() ) / ( 1.*fCurrent_Track->GetTPCNcls() ) > kChi2TPCcluster ) && ( fCurrent_Track->GetTPCNcls() != 0 ) ) return false;
    
    //  -   //  CHI2 GLOBAL
    if ( ( fCurrent_Track->GetChi2TPCConstrainedVsGlobal() > kChi2TPCGlobal ) )                                 return false;

    //  -   //  REFIT
    if ( fCurrent_Track->IsOn(0x4) == false )                                                                   return false;
    
    //  ------------------------ ITS
    
    //  -   //  REFIT
    if ( fCurrent_Track->IsOn(0x40) == false )                                                                  return false;
        
    //  -   //  HIT ON SPD
    if ( !fCurrent_Track->HasPointOnITSLayer(0) && !fCurrent_Track->HasPointOnITSLayer(1) )                     return false;
        
    //  -   //  CHI2 PER CLUSTER
    if ( ( ( 1.*fCurrent_Track->GetITSchi2() ) / ( 1.*fCurrent_Track->GetITSNcls() ) > kChi2ITScluster ) && ( fCurrent_Track->GetITSNcls() != 0 ) ) return false;
    
    if ( fIs_p_p )  {
        //  ------------------------ CUSTOM IMPLEMENTATION FOR pp
        
        //  Nothing here
        
    }   else if ( fIs_p_Pb || fIs_Pb_Pb ) {
        //  ------------------------ CUSTOM IMPLEMENTATION FOR pPb AND PbPb
        
        //  ------------------------ TPC
        
        //  -   //  ACCEPT EDGES
        Float_t fDeadZoneWidth=2;
        Float_t fCutGeoNcrNclLength=130.;
        Float_t fCutGeoNcrNclGeom1Pt=1.5;
        Float_t fCutGeoNcrNclFractionNcr=0;
        Float_t fCutGeoNcrNclFractionNcl=0;
        //
        Double_t xyz[3],pxpypz[3],cv[21];
        fCurrent_Track->GetXYZ(xyz);
        pxpypz[0]=fCurrent_Track->Px();
        pxpypz[1]=fCurrent_Track->Py();
        pxpypz[2]=fCurrent_Track->Pz();
        fCurrent_Track->GetCovarianceXYZPxPyPz(cv);
        Short_t sign = (Short_t)fCurrent_Track->Charge();
        //
        AliExternalTrackParam*              fCurrentTrackExternalParameters = new AliExternalTrackParam (xyz,pxpypz,cv,sign);
        Float_t kTPCActiveLenght    =       uTrackLengthInActiveTPC( fCurrentTrackExternalParameters, fDeadZoneWidth, 220. );
        if ( fCutGeoNcrNclLength > 0 ) {
            fCutGeoNcrNclLength -=  TMath::Power(TMath::Abs(1./(fCurrent_Track_Momentum)),fCutGeoNcrNclGeom1Pt);
            if ( kTPCActiveLenght                       < fCutGeoNcrNclLength )                             return false;
            if ( fCurrent_Track->GetTPCCrossedRows()    < fCutGeoNcrNclLength*fCutGeoNcrNclFractionNcr )    return false;
            if ( fCurrent_Track->GetTPCncls()           < fCutGeoNcrNclLength*fCutGeoNcrNclFractionNcl )    return false;
        }
    }
    
    return  true;
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fQC_TRK( )                                 {
    fQC_Tracks_Momentum ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Momentum );
    fQC_Tracks_TMomentum->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_TransMom );
    fQC_Tracks_Eta      ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Eta      );
    fQC_Tracks_Phi      ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Phi      );
    fQC_Tracks_V0M      ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_V0M            );
    fQC_Tracks_TRK      ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_TRK            );
    fQC_Tracks_DCAXY_P  ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Momentum,    fCurrent_Track_DCAXY);
    fQC_Tracks_DCAXY_PT ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_TransMom,    fCurrent_Track_DCAXY);
    fQC_Tracks_DCAZ_P   ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Momentum,    fCurrent_Track_DCAZ);
    fQC_Tracks_DCAZ_PT  ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_TransMom,    fCurrent_Track_DCAZ);
    fQC_Tracks_TPC_CLS  ->Fill(fCurrent_Track->GetTPCCrossedRows() );
    if ( fCurrent_Track->GetTPCNclsF()  != 0 )  fQC_Tracks_TPC_FRC  ->Fill( ( fCurrent_Track->GetTPCCrossedRows() / fCurrent_Track->GetTPCNclsF() ) );
    if ( fCurrent_Track->GetTPCNcls()   != 0 )  fQC_Tracks_TPC_CHI  ->Fill( ( 1.*fCurrent_Track->GetTPCchi2() ) / ( 1.*fCurrent_Track->GetTPCNcls() ) );
    fQC_Tracks_TPC_CNS  ->Fill(fCurrent_Track->GetChi2TPCConstrainedVsGlobal() );
    if ( fCurrent_Track->GetITSNcls()   != 0 )  fQC_Tracks_ITS_CHI  ->Fill( ( 1.*fCurrent_Track->GetITSchi2() ) / ( 1.*fCurrent_Track->GetITSNcls() ) );
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fQC_TRK_Kaons( )                           {
    fQC_Kaons_Momentum  ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Momentum );
    fQC_Kaons_TMomentum ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_TransMom );
    fQC_Kaons_Eta       ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Eta      );
    fQC_Kaons_Phi       ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Phi      );
    fQC_Kaons_V0M       ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_V0M            );
    fQC_Kaons_TRK       ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_TRK            );
    //fQC_Kaons_DCAXY_P   ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Momentum,    fCurrent_Track_DCAXY);
    //fQC_Kaons_DCAXY_PT  ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_TransMom,    fCurrent_Track_DCAXY);
    //fQC_Kaons_DCAZ_P    ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Momentum,    fCurrent_Track_DCAZ);
    //fQC_Kaons_DCAZ_PT   ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_TransMom,    fCurrent_Track_DCAZ);
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fQC_PID ( )                                {
    fIsTPCAvailable     = (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, fCurrent_Track) == AliPIDResponse::kDetPidOk);
    fIsTOFAvailable     = (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, fCurrent_Track) == AliPIDResponse::kDetPidOk);
    fBetaFromTOFSignal  = fTOFBeta();
    fTPCSignal          = fCurrent_Track->GetTPCsignal();
    
    // PID Signal histograms
    if ( fIsTPCAvailable )    {
        fQC_Tracks_TPC_P    ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Momentum,    fTPCSignal);
        fQC_Tracks_TPC_PT   ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_TransMom,    fTPCSignal);
    }
    if ( fIsTOFAvailable )    {
        fQC_Tracks_TOF_P    ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Momentum,    fBetaFromTOFSignal);
        fQC_Tracks_TOF_PT   ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_TransMom,    fBetaFromTOFSignal);
    }
    
    // Sigma PID QC for all particles
    fQC_PID_Kaons();
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fQC_PID_Kaons( )                           {
    auto fSigma_TOF     = (fPIDResponse->NumberOfSigmasTOF(fCurrent_Track,AliPID::kKaon));
    auto fSigma_TPC     = (fPIDResponse->NumberOfSigmasTPC(fCurrent_Track,AliPID::kKaon));
    
    // PID Sigmas histograms
    if ( fIsTPCAvailable )    {
        fQC_PID_TPC_Kaons_P ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Momentum,             fSigma_TPC);
        fQC_PID_TPC_Kaons_PT->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_TransMom,             fSigma_TPC);
    }
    if ( fIsTOFAvailable )    {
        fQC_PID_TOF_Kaons_P ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Momentum,             fSigma_TOF);
        fQC_PID_TOF_Kaons_PT->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_TransMom,             fSigma_TOF);
    }
    if ( fIsTPCAvailable && fIsTOFAvailable )
    {
        fQC_PID_TPC_TOF_Kaons_PT    ->Fill(fCurrent_Track_Charge*0.5,   fSigma_TOF,     fSigma_TPC);
    }
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fQC_PID_Sel_Kaons( )                       {
    auto fSigma_TOF     = (fPIDResponse->NumberOfSigmasTOF(fCurrent_Track,AliPID::kKaon));
    auto fSigma_TPC     = (fPIDResponse->NumberOfSigmasTPC(fCurrent_Track,AliPID::kKaon));
    
    // PID Signal histograms
    if ( fIsTPCAvailable )    {
        fQC_PID_TPC_NSig_SEL_Kaons_P    ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Momentum,    fSigma_TPC);
        fQC_PID_TPC_NSig_SEL_Kaons_PT   ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_TransMom,    fSigma_TPC);
        fQC_PID_TPC_Sgnl_SEL_Kaons_P    ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Momentum,    fTPCSignal);
        fQC_PID_TPC_Sgnl_SEL_Kaons_PT   ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_TransMom,    fTPCSignal);
        fQC_Kaons_TPC_P                 ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Momentum,    fTPCSignal);
        fQC_Kaons_TPC_PT                ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_TransMom,    fTPCSignal);
    }
    if ( fIsTOFAvailable )    {
        fQC_PID_TOF_NSig_SEL_Kaons_P    ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Momentum,    fSigma_TOF);
        fQC_PID_TOF_NSig_SEL_Kaons_PT   ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_TransMom,    fSigma_TOF);
        fQC_PID_TOF_Sgnl_SEL_Kaons_P    ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Momentum,    fBetaFromTOFSignal);
        fQC_PID_TOF_Sgnl_SEL_Kaons_PT   ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_TransMom,    fBetaFromTOFSignal);
        fQC_Kaons_TOF_P                 ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_Momentum,    fBetaFromTOFSignal);
        fQC_Kaons_TOF_PT                ->Fill(fCurrent_Track_Charge*0.5,   fCurrent_Track_TransMom,    fBetaFromTOFSignal);
    }
    if ( fIsTPCAvailable && fIsTOFAvailable )   {
        // TODO: TPC TOF VETO and SELECTED nsigma vs nsgima
    }
}

//_____________________________________________________________________________

Double_t    AliAnalysisTaskPhiCount::fTOFBeta( )                                {
    if ( !fIsTOFAvailable ) return -1.;
    Float_t fTrackTime(fCurrent_Track->GetTOFsignal()*1e-12);
    Float_t fTrackLngt(fCurrent_Track->GetIntegratedLength()*1e-2);
    if ( fTrackTime <= 0. || fTrackLngt <= 0. ) return -1.;
    return  fTrackLngt/(fTrackTime*TMath::C());
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsKaonCandidate ( )                       {
    auto fSigma_TOF     = fPIDResponse->NumberOfSigmasTOF(fCurrent_Track,AliPID::kKaon);
    auto fSigma_TPC     = fPIDResponse->NumberOfSigmasTPC(fCurrent_Track,AliPID::kKaon);
    
    if ( std::fabs(fSigma_TPC) <= 10. && std::fabs(fSigma_TOF) <= 10. )   {
        fTOFSigma[fnKaon]   =   (Char_t)(fSigma_TOF*10.);
        fTPCSigma[fnKaon]   =   (Char_t)(fSigma_TPC*10.);
    }
    if ( fCurrent_Run >= 115881 && fCurrent_Run <= 139513 )   {
        // Custom 2010
        if ( !fIsTPCAvailable || (fIsTOFAvailable && std::fabs(fSigma_TOF) > kSgTOF_Veto) )         return false;
        if ( std::fabs(fSigma_TPC) > 7. )                                                           return false;
        if ( fCurrent_Track_Momentum >= 0.28 )   {
            if ( fIsTOFAvailable && std::fabs(fSigma_TPC) > kSgTPC_TOFVt )                          return false;
            if (!fIsTOFAvailable && std::fabs(fSigma_TPC) > kSgTPC_Alone )                          return false;
        }
    }   else    {
        // Standard
        if ( !fIsTPCAvailable || ( fIsTOFAvailable && std::fabs(fSigma_TOF) >= kSgTOF_Veto ) )      return false;
        if (  fIsTOFAvailable && std::fabs(fSigma_TPC) > kSgTPC_TOFVt )                             return false;
        if ( !fIsTOFAvailable && std::fabs(fSigma_TPC) > kSgTPC_Alone )                             return false;
        return true;
    }
    return true;
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsPhiCandidate ( TLorentzVector fPhi )    {
    if ( fPhi.Mag() < 0.75 || fPhi.Mag() > 1.25 ) return false;
    return true;
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsCandidateTruPhi ( AliAODMCParticle* piKaon, AliAODMCParticle* pjKaon )  {
    if ( !piKaon || !pjKaon ) return false;
    if ( piKaon->GetMother() == -1 ) return false;
    if ( pjKaon->GetMother() == -1 ) return false;
    auto    fiKaonMother    =   static_cast<AliAODMCParticle*>(AODMCTrackArray->At(piKaon->GetMother()));
    auto    fjKaonMother    =   static_cast<AliAODMCParticle*>(AODMCTrackArray->At(pjKaon->GetMother()));
    if  ( fiKaonMother->GetPdgCode() != fjKaonMother->GetPdgCode() || fjKaonMother->GetPdgCode() != 333 ) return false;
    return fiKaonMother == fjKaonMother;
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsPhiGen ( AliAODMCParticle* particle )   {
    if ( particle->GetNDaughters() != 2 ) return false;
    auto const Dau1 =   static_cast<AliAODMCParticle*>  (AODMCTrackArray->At(particle->GetDaughterFirst()));
    auto const Dau2 =   static_cast<AliAODMCParticle*>  (AODMCTrackArray->At(particle->GetDaughterLast()));
    if ( !Dau1 || !Dau2 ) return false;
    return  ( ( Dau1->GetPdgCode() == -Dau2->GetPdgCode() ) && ( abs(Dau1->GetPdgCode()) == 321 ) );
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsPhiRec ( AliAODMCParticle* particle )   {
    // looping over all recordable phis
    for ( int iPhi = 0; iPhi < fnPhiRec; iPhi++ )
    {
        if ( fPhiRecParticles[iPhi] == particle ) return true;
    }
    return false;
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fIsPhi ( AliAODMCParticle* particle )      {
        if ( !particle ) return false;
        if ( particle->GetPdgCode() == 333 ) return true;
        else return false;
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::Terminate( Option_t* )                     {
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fSetEventMask( Int_t iMaskBit )            {
    auto    fResult =   fEventMask  |   BIT( iMaskBit );
    fEventMask      =   fResult;
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fSetTrueEventMask( Int_t iMaskBit )        {
    auto    fResult     =   fTrueEventMask  |   ( iMaskBit );
    fTrueEventMask      =   fResult;
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fCheckMask( Int_t iMaskBit )               {
    return  ( ( BIT(iMaskBit) & fEventMask )        ==  BIT(iMaskBit) );
}

//_____________________________________________________________________________

bool        AliAnalysisTaskPhiCount::fCheckTrueMask( Int_t iMaskBit )           {
    return  ( ( BIT(iMaskBit) & fTrueEventMask )    ==  BIT(iMaskBit) );
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fPostData( )                               {
    // Post-data for TLists
    PostData(1, fAnalysisOutputList);
    PostData(2, fQCOutputList);
    // Post-data for TTrees
    if ( kPhibool )                 PostData(3, fPhiCandidate);
    if ( kKaonbool )                PostData(4, fKaonCandidate);
    if ( kPhibool   &&  kMCbool )   PostData(5, fPhiEfficiency);
    if ( kKaonbool  &&  kMCbool )   PostData(6, fKaonEfficiency);
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fFillTrees( )                              {
    if ( fnPhi      != 0 )  {
        fPhiCandidate   ->Fill();
    }
    if ( fnPhi      == 0 )  {
        fQC_Event_Enum_FLL->Fill("NoPhiCand",1);
    }
    if ( fnPhi      >= 153 ){
        fQC_Event_Enum_FLL->Fill("OFPhiCand",1);
    }
    if ( fnPhiTru   != 0 )  {
        fPhiEfficiency  ->Fill();
    }
    if ( fnPhiTru   == 0 )  {
        fQC_Event_Enum_FLL->Fill("NoPhTCand",1);
    }
    if ( fnPhiTru   >= 153 ){
        fQC_Event_Enum_FLL->Fill("OFPhTCand",1);
    }
    if ( fnKaon     != 0 )  {
        //fKaonCandidate   ->Fill();
    }
    if ( fnKaon     == 0 )  {
        fQC_Event_Enum_FLL->Fill("NoKaonCand",1);
    }
    if ( fnKaon     >= 153 ){
        fQC_Event_Enum_FLL->Fill("OFKaonCand",1);
    }
    if ( fnKaonTru  != 0 )  {
        //fKaonEfficiency   ->Fill();
    }
    if ( fnKaonTru  == 0 )  {
        fQC_Event_Enum_FLL->Fill("NoKaoTCand",1);
    }
    if ( fnKaonTru  >= 153 ){
        fQC_Event_Enum_FLL->Fill("OFKaoTCand",1);
    }
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fSetZero( )                                {
    fnPhi = 0;
    fnPhiTru = 0;
    fnKaon = 0;
    fnKaonTru = 0;
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::fStoreTruePhi ( Int_t iMaskBit )           {
    // Loop over all primary MC particle
    if ( !kMCbool ) return;
    AliMCVertex*   PrimaryVertexMCT = (AliMCVertex*)fMCD->GetPrimaryVertex();
    if ( !PrimaryVertexMCT ) return;
    if ( std::fabs(PrimaryVertexMCT->GetZ()) < kVertexCut ) fSetTrueEventMask(kTRU_HAST10VTX);
    Int_t           nTrack(AODMCTrackArray->GetEntriesFast());
    for ( Int_t iTrack(0); iTrack < nTrack; iTrack++ )
    {
        AliAODMCParticle* fPhiTru = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(iTrack));
        
        if ( !fIsPhi( fPhiTru ) ) continue;
        
        fCheckINELgt0( fPhiTru );
        
        // Kinematics
        fPhiTruPx[fnPhiTru]        =   fPhiTru->Px();
        fPhiTruPy[fnPhiTru]        =   fPhiTru->Py();
        fPhiTruPz[fnPhiTru]        =   fPhiTru->Pz();
            
        // Options
        fSelection[fnPhiTru]       =   0;
        if( fIsPhiGen(fPhiTru) )    fSelection[fnPhiTru] = 1;
        
        fnPhiTru++;
    }
    fSetTrueEventMask(iMaskBit);
    fPhiEfficiency  ->Fill();
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::uCalculateSpherocity ( )                   {
    //
    fCurrent_SPH                =   -1.;
    if ( !kComputeSpherocity ) return;
    auto    nTracks             =   fAOD->GetNumberOfTracks();
    //
    //  Requires at least 10 good tracks, so 10 tracks overall
    if ( nTracks < 10 ) return;
    else                fCurrent_SPH     =   -2.;
    //
    auto    nTrackPassed        =   0;
    auto    fTotalTransMom      =   0.;
    auto    fTransverseMomenta  =   new Double_t[1000];
    auto    fParticlePhi        =   new Double_t[1000];
    for ( Int_t iTrack = 0; iTrack < nTracks; iTrack++ )    {
        fCurrent_Track  =   static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
        //
        //  Selecting Good tracks
        if ( !fCurrent_Track->TestFilterBit(BIT(0)) )   continue;
        if ( !fCurrent_Track->IsOn(0x40) )              continue;
        if ( !fCurrent_Track->IsOn(0x4) )               continue;
        if ( TMath::Abs(fCurrent_Track->Eta()) > 0.8 )  continue;
        if ( fCurrent_Track->Pt() < 0.15 )              continue;
        //
        //  Storing necessary Info
        fTotalTransMom                     +=   kSpherocityPTWeight? fCurrent_Track->Pt() : 1;
        fTransverseMomenta[nTrackPassed]    =   kSpherocityPTWeight? fCurrent_Track->Pt() : 1;
        fParticlePhi[nTrackPassed]          =   fCurrent_Track->Phi();
        nTrackPassed++;
    }
    //
    //  Requires at least 10 good tracks
    if ( nTrackPassed < 10 ) return;
    else                fCurrent_SPH     =   1.;
    //
    auto    fNumerator  =   0.;
    for ( Int_t iAngle = 0; iAngle < 36000; iAngle++ )   {
        auto    fCurrent_Angle  =   (TMath::Pi()*iAngle)/18000.;
        auto    fXProjection    =   TMath::Cos(fCurrent_Angle);
        auto    fYProjection    =   TMath::Sin(fCurrent_Angle);
        for ( Int_t iGoodTrk = 0; iGoodTrk < nTrackPassed; iGoodTrk++ )  {
            auto    fXTransMom  =   fTransverseMomenta[iGoodTrk] * TMath::Cos( fParticlePhi[iGoodTrk] );
            auto    fYTransMom  =   fTransverseMomenta[iGoodTrk] * TMath::Sin( fParticlePhi[iGoodTrk] );
            fNumerator         +=   TMath::Abs(fXProjection*fYTransMom+fYProjection*fXTransMom);
        }
        auto    fCurrentSpherocity  =   TMath::Power( (fNumerator/fTotalTransMom),2 );
        if ( fCurrentSpherocity < fCurrent_SPH ) fCurrent_SPH =   fCurrentSpherocity;
    }
    fQC_Event_Spherocity->Fill(fCurrent_SPH);
    //
    return;
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::uCalculateRT ( )                           {
    //
    fCurrent_RT                 =   -1.;
    if ( !kComputeRT ) return;
    auto    nTracks             =   fAOD->GetNumberOfTracks();
    //
    //  Find a suitable trigger particle
    auto    fTriggerPT  =   -1.;
    auto    fTriggerPhi =   -1.;
    for ( Int_t iTrack = 0; iTrack < nTracks; iTrack++ )    {
        fCurrent_Track  =   static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
        //
        //  Selecting Good tracks
        if ( !fCurrent_Track->TestFilterBit(BIT(0)) )                   continue;
        if ( !fCurrent_Track->IsOn(0x40) )                              continue; // ITS Refit
        if ( !fCurrent_Track->IsOn(0x4) )                               continue; // TPC Refit
        if ( TMath::Abs(fCurrent_Track->Eta()) > 0.8 )                  continue;
        if ( fCurrent_Track->Pt() < 5. || fCurrent_Track->Pt() > 40. )  continue;
        //
        //  Storing necessary Info, choosing the highest pT particle
        if ( fTriggerPT > fCurrent_Track->Pt() )                        continue;
        fTriggerPT  = fCurrent_Track->Pt();
        fTriggerPhi = fCurrent_Track->Phi();
    }
    if ( fTriggerPT == -1. ) return;
    else    fCurrent_RT         =   -2.;
    //
    //  Calculate RT
    auto    nRForward       =   0;
    auto    nRBackward      =   0;
    auto    nRTransverse    =   0;
    for ( Int_t iTrack = 0; iTrack < nTracks; iTrack++ )    {
        fCurrent_Track  =   static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
        //
        //  Selecting Good tracks
        if ( !fCurrent_Track->TestFilterBit(BIT(0)) )   continue;
        if ( !fCurrent_Track->IsOn(0x40) )              continue;
        if ( !fCurrent_Track->IsOn(0x4) )               continue;
        if ( TMath::Abs(fCurrent_Track->Eta()) > 0.8 )  continue;
        if ( fCurrent_Track->Pt() < 0.15 )              continue;
        //
        //  Storing necessary Info
        auto    kDeltaPhi   =   (fTriggerPhi - fCurrent_Track->Phi())*180/(TMath::Pi());
        kDeltaPhi < 0 ? kDeltaPhi += 360 : kDeltaPhi += 0;
        if      ( ( kDeltaPhi > 300. ) || ( kDeltaPhi < 60.  ) )    nRForward++;
        else if ( ( kDeltaPhi > 120. ) && ( kDeltaPhi < 240. ) )    nRBackward++;
        else                                                        nRTransverse++;
    }
    fQC_Event_FullRF    ->  Fill(nRForward);
    fQC_Event_FullRB    ->  Fill(nRBackward);
    fQC_Event_FullRT    ->  Fill(nRTransverse);
    fCurrent_RT         =   nRTransverse;
    //
    return;
}

//_____________________________________________________________________________

Float_t     AliAnalysisTaskPhiCount::uTrackLengthInActiveTPC ( AliExternalTrackParam* fCurrentTrackExternalParameters, Double_t deltaY, Double_t deltaZ ) {
    const   Double_t    kInnerRadius    =   85;
    const   Double_t    kOuterRadius    =   245;
    const   Double_t    kMagField       =   fAOD->GetMagneticField();
    const   Double_t    kB2C            =   0.299792458e-3;             //  Mistery variable...
    Double_t            fResult         =   0;
    Double_t            kXYZ[3];
    Double_t            kPxPyPz[3];
    //
    //  Recover XYZ coordinates and momentum
    if ( fCurrentTrackExternalParameters->GetXYZAt  ( kInnerRadius, kMagField, kXYZ ) ) {
        fCurrentTrackExternalParameters->GetPxPyPzAt( kInnerRadius, kMagField, kPxPyPz );
    } else {
        fCurrentTrackExternalParameters->GetXYZ     ( kXYZ );
        fCurrentTrackExternalParameters->GetPxPyPz  ( kPxPyPz );
    }
    //
    Float_t     fCurrent_DCA        =   -fCurrentTrackExternalParameters->GetD(0,0,kMagField);
    Float_t     fCurrent_Radius     =   TMath::Abs(1./( fCurrentTrackExternalParameters->GetC(kMagField) ));
    Float_t     fCurrent_Sign       =   fCurrentTrackExternalParameters->GetSign()*TMath::Sign(1.,kMagField)*(-1.);
    Float_t     fCurrent_R0         =   TMath::Sqrt(kXYZ[0]*kXYZ[0]+kXYZ[1]*kXYZ[1]);
    Float_t     fCurrent_PhiR0      =   TMath::ATan2(kXYZ[1],kXYZ[0]);
    Float_t     fCurrent_dPhiR0     =   -TMath::ASin((fCurrent_Track_DCAXY*fCurrent_Track_DCAXY-2*fCurrent_Track_DCAXY*fCurrent_Radius*fCurrent_Sign+fCurrent_R0*fCurrent_R0)/(2*fCurrent_R0*(fCurrent_Track_DCAXY-fCurrent_Radius*fCurrent_Sign)));
    Float_t     fCurrent_Phi0       =   fCurrent_PhiR0 - fCurrent_dPhiR0;
    //
    for ( Double_t uRadiusScan = kInnerRadius; uRadiusScan <= kOuterRadius; uRadiusScan++ )  {
        Float_t fSinPhi =   (fCurrent_Track_DCAXY*fCurrent_Track_DCAXY-2*fCurrent_Track_DCAXY*fCurrent_Radius*fCurrent_Sign+uRadiusScan*uRadiusScan)/(2*uRadiusScan*(fCurrent_Track_DCAXY-fCurrent_Radius*fCurrent_Sign));
        if ( TMath::Abs( fSinPhi ) >= 1 ) continue;
        Float_t     fdPhi       =   -TMath::ASin(fSinPhi);
        Float_t     fPhi        =   fCurrent_Phi0 + fdPhi;
        Int_t       nSector     =   TMath::Nint(9*fPhi/(TMath::Pi()));
        Float_t     fdPhiEdge   =   fPhi-(nSector*TMath::Pi()/9);
        Float_t     fX          =   fCurrent_Radius*TMath::Cos(fPhi) - kXYZ[0];
        Float_t     fY          =   fCurrent_Radius*TMath::Sin(fPhi) - kXYZ[1];
        Float_t     fdeltaPhi   =   2*TMath::ASin(0.5*TMath::Sqrt(fX*fX+fY*fY)/fCurrent_Radius);
        Float_t     fZ          =   kXYZ[2]+fdeltaPhi*fCurrent_Radius*fCurrentTrackExternalParameters->GetTgl();
        if ( TMath::Abs( fdPhiEdge*uRadiusScan ) > deltaY && TMath::Abs( fZ ) < deltaZ ) {
            fResult++;
        }
    }
    return fResult;
}

//_____________________________________________________________________________

void        AliAnalysisTaskPhiCount::uGuessCollisionSystem ( )                  {
    fIs_p_p     =   kTRUE;
    fIs_p_Pb    =   kFALSE;
    fIs_Pb_Pb   =   kFALSE;
    if ( ( fCurrent_Run >= 195344 && fCurrent_Run <= 195677 ) || ( fCurrent_Run >= 265309 && fCurrent_Run <= 267166 ) ) {
        fIs_p_p     =   kFALSE;
        fIs_p_Pb    =   kTRUE;
        fIs_Pb_Pb   =   kFALSE;
        fQC_Event_Enum_V0M          ->  GetXaxis()  ->  SetTitle("V0A Multiplicity");
        fQC_Event_Enum_V0T          ->  GetXaxis()  ->  SetTitle("V0A Multiplicity");
        fQC_Tracks_V0M              ->  GetYaxis()  ->  SetTitle("V0A Multiplicity");
    }
}
