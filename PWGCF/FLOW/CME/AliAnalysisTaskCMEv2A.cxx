#include <iostream>
#include <cstdlib>
#include <sys/time.h>

// ROOT classes
#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TExMap.h"
#include "TRandom3.h"
#include "TMath.h"

// Alice analysis base class
#include "AliAnalysisTaskSE.h"

// Alice analysis additional classes
#include "AliAnalysisManager.h"
//#include "AliInputEventHandler.h"

// Alice AOD classes
#include "AliAODInputHandler.h"
//#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"

// Alice classes
#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliAnalysisUtils.h"

// Alice MC classes
#include "AliMCEvent.h"
//#include "AliMCEventHandler.h"
#include "AliAODMCParticle.h"

// Alice "V" classes
#include "AliVParticle.h"
//#include "AliVEvent.h"
//#include "AliVVertex.h"
//#include "AliVVZERO.h"

// Alice PID classes
//#include "AliAODPid.h"
//#include "AliAODpidUtil.h"
#include "AliPID.h"
//#include "AliPIDCombined.h"
#include "AliPIDResponse.h"


// This class
#include "AliAnalysisTaskCMEv2A.h"

ClassImp(AliAnalysisTaskCMEv2A); // import class inheriting from TObject

using std::cout;
using std::endl;

const float pi = 3.141592653589793;

// function prototype...
float GetDPhiStar(float phi1, float pt1, float charge1, float phi2, float pt2, float charge2, float radius, float bSign);


//------------------------------------------------------------------
AliAnalysisTaskCMEv2A::AliAnalysisTaskCMEv2A() : AliAnalysisTaskSE(),
    debug(0),
    doMC(false),
    trigger(AliVEvent::kMB),
    dopupcut(true),
    doeffcorr(true),
    centhandle(1),
    fbit(128),
    zvtxcut(10.0),
    centcut(5.0),
    nclscut(60),
    dcacutz(5.0),
    dcacutxy(5.0),
    dodcacuts(false),
    outeta(0.8),
    ineta(0.5),
    excleta(0.0),
    ptmin(0.2),
    ptmax(5.0),
    doacuts(false),
    nspid(2.0),
    cbinlo(3),
    cbinhi(4),
    donested(true),
    dopaircut(true),
    centlo(20.0),
    centhi(60.0),
    fOutputList(0x0),     
    fHistDCACOVerrors(0x0),
    fHistPt(0x0),
    fHistPhi(0x0),
    fHistEta(0x0),
    fHistCharge(0x0),
    fHistTPCncls(0x0),
    fHistDedx(0x0),
    fHistDCAxy(0x0),
    fHistDCAz(0x0),
    fHistDCAxyAfter(0x0),
    fHistDCAzAfter(0x0),
    fHistPosPt(0x0),
    fHistPosPhi(0x0),
    fHistPosEta(0x0),
    fHistNegPt(0x0),
    fHistNegPhi(0x0),
    fHistNegEta(0x0),
    fHistPtPion(0x0),
    fHistPtPionH(0x0),
    fHistPtPionL(0x0),
    fHistNsigmaPion(0x0),
    fHistPtProt(0x0),
    fHistPtProtH(0x0),
    fHistPtProtL(0x0),
    fHistNsigmaProt(0x0),
    fHistPtMC(0x0),
    fHistPhiMC(0x0),
    fHistEtaMC(0x0),
    fHistChargeMC(0x0),
    fHistTPCnclsMC(0x0),
    fHistDedxMC(0x0),
    fHistDCAxyMC(0x0),
    fHistDCAzMC(0x0),
    fHistPlaneV0h2(0x0),
    fHistPlaneV0Ah2(0x0),
    fHistPlaneV0Ch2(0x0),
    fHistPlaneV0ACDCh2(0x0),
    fHistZVtx(0x0),
    fHistZVtxD(0x0),
    fHistZVtxMC(0x0),
    fHistZVtxDiffMC(0x0),
    fHistCentTRK(0x0),
    fHistCentV0M(0x0),
    fHistCentDIFF(0x0),
    fHistCentDIAG(0x0),
    fHistCtrkDIAG(0x0),
    fHistVtxRDIAG(0x0),
    fHistVtxSDIAG(0x0),
    fHistVtxRDIBG(0x0),
    fHistVtxSDIBG(0x0),
    fHistCentTRKAVEkMB(0x0),
    fHistCentV0MAVEkMB(0x0),
    fHistCentTRKAVEkCentral(0x0),
    fHistCentV0MAVEkCentral(0x0),
    fHistCentTRKAVEkSemiCentral(0x0),
    fHistCentV0MAVEkSemiCentral(0x0),
    fHistCentTRKAVEkA3(0x0),
    fHistCentV0MAVEkA3(0x0),
    fHistCentTRKAVEkSel(0x0),
    fHistCentV0MAVEkSel(0x0),
    
    h_MCq22_cent(0x0),
    h_MCq22_centP(0x0),
    h_MCq22_centN(0x0),
    h_MCAq22_cent(0x0),
    h_MCAq22_centP(0x0),
    h_MCAq22_centN(0x0),
    h_MCq22gap0_cent(0x0),
    h_MCq22gap0_centP(0x0),
    h_MCq22gap0_centN(0x0),
    h_MCq22gap1_cent(0x0),
    h_MCq22gap1_centP(0x0),
    h_MCq22gap1_centN(0x0),
    
    hMC_AT_X_deta(0x0),
    hMC_AT_X_detaP(0x0),
    hMC_AT_X_detaN(0x0),
    hMC_diffq22_X_deta(0x0),
    hMC_diffq22_X_detaP(0x0),
    hMC_diffq22_X_detaN(0x0),
    hMC_ATdiffq22_X_deta(0x0),
    hMC_ATdiffq22_X_detaP(0x0),
    hMC_ATdiffq22_X_detaN(0x0),
    hMC_diffq32_X_deta(0x0),
    hMC_diffq32_X_detaP(0x0),
    hMC_diffq32_X_detaN(0x0),
    hMC_ATdiffq32_X_deta(0x0),
    hMC_ATdiffq32_X_detaP(0x0),
    hMC_ATdiffq32_X_detaN(0x0),
    hMC_diffq42_X_deta(0x0),
    hMC_diffq42_X_detaP(0x0),
    hMC_diffq42_X_detaN(0x0),
    hMC_ATdiffq42_X_deta(0x0),
    hMC_ATdiffq42_X_detaP(0x0),
    hMC_ATdiffq42_X_detaN(0x0),
    
    fProfMeanChargePt(0x0),
    fProfMeanChargeEta(0x0),
    
    h_X22_cent(0x0),
    h_X22_centP(0x0),
    h_X22_centN(0x0),
    h_Y22_cent(0x0),
    h_Y22_centP(0x0),
    h_Y22_centN(0x0),
    h_X22_lo_cent(0x0),
    h_X22_lo_centP(0x0),
    h_X22_lo_centN(0x0),
    h_Y22_lo_cent(0x0),
    h_Y22_lo_centP(0x0),
    h_Y22_lo_centN(0x0),
    h_X22_li_cent(0x0),
    h_X22_li_centP(0x0),
    h_X22_li_centN(0x0),
    h_Y22_li_cent(0x0),
    h_Y22_li_centP(0x0),
    h_Y22_li_centN(0x0),
    h_X22_ri_cent(0x0),
    h_X22_ri_centP(0x0),
    h_X22_ri_centN(0x0),
    h_Y22_ri_cent(0x0),
    h_Y22_ri_centP(0x0),
    h_Y22_ri_centN(0x0),
    h_X22_ro_cent(0x0),
    h_X22_ro_centP(0x0),
    h_X22_ro_centN(0x0),
    h_Y22_ro_cent(0x0),
    h_Y22_ro_centP(0x0),
    h_Y22_ro_centN(0x0),
    
    h_q22_cent(0x0),
    h_q22_centP(0x0),
    h_q22_centN(0x0),
    h_q22_centPP(0x0),
    h_q22_centNN(0x0),
    h_q23_cent(0x0),
    h_q23_centP(0x0),
    h_q23_centN(0x0),
    h_q24_cent(0x0),
    h_q24_centP(0x0),
    h_q24_centN(0x0),
    h_q22gap0_cent(0x0),
    h_q22gap0_centP(0x0),
    h_q22gap0_centN(0x0),
    h_q22gap1_cent(0x0),
    h_q22gap1_centP(0x0),
    h_q22gap1_centN(0x0),
    
    h_q32_cent(0x0),
    h_q32_centP(0x0),
    h_q32_centN(0x0),
    h_q32_centPP(0x0),
    h_q32_centNN(0x0),
    h_q32gap0_cent(0x0),
    h_q32gap0_centP(0x0),
    h_q32gap0_centN(0x0),
    h_q32gap1_cent(0x0),
    h_q32gap1_centP(0x0),
    h_q32gap1_centN(0x0),
    
    h_q42_cent(0x0),
    h_q42_centP(0x0),
    h_q42_centN(0x0),
    h_q42_centPP(0x0),
    h_q42_centNN(0x0),
    h_q42gap0_cent(0x0),
    h_q42gap0_centP(0x0),
    h_q42gap0_centN(0x0),
    h_q42gap1_cent(0x0),
    h_q42gap1_centP(0x0),
    h_q42gap1_centN(0x0),
    
    //TProfile *h_q22qasym_cent[10];
    //TProfile *h_q22qasymP_cent[10];
    //TProfile *h_q22qasymN_cent[10];
    //TProfile *h_q22qasym_gap0_cent[10];
    //TProfile *h_q22qasymP_gap0_cent[10];
    //TProfile *h_q22qasymN_gap0_cent[10];
    //TProfile *h_q22qasym_gap1_cent[10];
    //TProfile *h_q22qasymP_gap1_cent[10];
    //TProfile *h_q22qasymN_gap1_cent[10];
    //TProfile *h_q24qasym_cent[10];
    //TProfile *h_q24qasymP_cent[10];
    //TProfile *h_q24qasymN_cent[10];
    //TProfile *h_q32qasym_cent[10];
    //TProfile *h_q32qasymP_cent[10];
    //TProfile *h_q32qasymN_cent[10];
    //TProfile *h_q32qasym_gap0_cent[10];
    //TProfile *h_q32qasymP_gap0_cent[10];
    //TProfile *h_q32qasymN_gap0_cent[10];
    //TProfile *h_q32qasym_gap1_cent[10];
    //TProfile *h_q32qasymP_gap1_cent[10];
    //TProfile *h_q32qasymN_gap1_cent[10];
    //TProfile *h_q34qasym_cent[10];
    //TProfile *h_q34qasymP_cent[10];
    //TProfile *h_q34qasymN_cent[10];
    //TProfile *h_q42qasym_cent[10];
    //TProfile *h_q42qasymP_cent[10];
    //TProfile *h_q42qasymN_cent[10];
    //TProfile *h_q42qasym_gap0_cent[10];
    //TProfile *h_q42qasymP_gap0_cent[10];
    //TProfile *h_q42qasymN_gap0_cent[10];
    //TProfile *h_q42qasym_gap1_cent[10];
    //TProfile *h_q42qasymP_gap1_cent[10];
    //TProfile *h_q42qasymN_gap1_cent[10];
    //TProfile *h_q44qasym_cent[10];
    //TProfile *h_q44qasymP_cent[10];
    //TProfile *h_q44qasymN_cent[10];
    
    h_eta_pos_F1(0x0),
    h_eta_pos_F3(0x0),
    h_eta_neg_F1(0x0),
    h_eta_neg_F3(0x0),
    h_cut_eta_pos_F1(0x0),
    h_cut_eta_pos_F3(0x0),
    h_cut_eta_neg_F1(0x0),
    h_cut_eta_neg_F3(0x0),
    
    h_AT_eta(0x0),
    h_AT_etaF1(0x0),
    h_AT_etaF3(0x0),
    h_AT_etaMC(0x0),
    h2_AT_eta(0x0),
    h2_AT_etaF1(0x0),
    h2_AT_etaF3(0x0),
    h2_AT_etaMC(0x0),
    
    h_AT_cut_eta(0x0),
    h_AT_cut_etaF1(0x0),
    h_AT_cut_etaF3(0x0),
    h_AT_cut_etaMC(0x0),
    h2_AT_cut_eta(0x0),
    h2_AT_cut_etaF1(0x0),
    h2_AT_cut_etaF3(0x0),
    h2_AT_cut_etaMC(0x0),
    
    h_A_cent(0x0),
    h_rA_cent(0x0),
    h_r1A_cent(0x0),
    h_r2A_cent(0x0),
    h_A_centF1(0x0),
    h_A_centF3(0x0),
    h_A_centMC(0x0),
    h2_A_cent(0x0),
    h2_rA_cent(0x0),
    h2_r1A_cent(0x0),
    h2_r2A_cent(0x0),
    h2_A_centF1(0x0),
    h2_A_centF3(0x0),
    h2_A_centMC(0x0),
    
    h_Aq22_cent(0x0),
    h_Aq22_centP(0x0),
    h_Aq22_centN(0x0),
    h_Aq22_centPP(0x0),
    h_Aq22_centNN(0x0),
    h_Aq22_SL_cent(0x0),
    h_Aq22_SL_centP(0x0),
    h_Aq22_SL_centN(0x0),
    h_Aq22_SLL_cent(0x0),
    h_Aq22_SLL_centP(0x0),
    h_Aq22_SLL_centN(0x0),
    h_Aq22_SLR_cent(0x0),
    h_Aq22_SLR_centP(0x0),
    h_Aq22_SLR_centN(0x0),
    h_Aq32_cent(0x0),
    h_Aq32_centP(0x0),
    h_Aq32_centN(0x0),
    h_Aq32_centPP(0x0),
    h_Aq32_centNN(0x0),
    h_Aq32_SL_cent(0x0),
    h_Aq32_SL_centP(0x0),
    h_Aq32_SL_centN(0x0),
    h_Aq32_SLL_cent(0x0),
    h_Aq32_SLL_centP(0x0),
    h_Aq32_SLL_centN(0x0),
    h_Aq32_SLR_cent(0x0),
    h_Aq32_SLR_centP(0x0),
    h_Aq32_SLR_centN(0x0),
    h_Aq42_cent(0x0),
    h_Aq42_centP(0x0),
    h_Aq42_centN(0x0),
    h_Aq42_centPP(0x0),
    h_Aq42_centNN(0x0),
    h_Aq42_SL_cent(0x0),
    h_Aq42_SL_centP(0x0),
    h_Aq42_SL_centN(0x0),
    h_Aq42_SLL_cent(0x0),
    h_Aq42_SLL_centP(0x0),
    h_Aq42_SLL_centN(0x0),
    h_Aq42_SLR_cent(0x0),
    h_Aq42_SLR_centP(0x0),
    h_Aq42_SLR_centN(0x0),
    
    h_diffq22_pt(0x0),
    h_diffq22_ptP(0x0),
    h_diffq22_ptN(0x0),
    h_diffAq22_pt(0x0),
    h_diffAq22_ptP(0x0),
    h_diffAq22_ptN(0x0),
    h_diffq22_eta(0x0),
    h_diffq22_etaP(0x0),
    h_diffq22_etaN(0x0),
    h_diffAq22_eta(0x0),
    h_diffAq22_etaP(0x0),
    h_diffAq22_etaN(0x0),
    h_diffq32_pt(0x0),
    h_diffq32_ptP(0x0),
    h_diffq32_ptN(0x0),
    h_diffAq32_pt(0x0),
    h_diffAq32_ptP(0x0),
    h_diffAq32_ptN(0x0),
    h_diffq32_eta(0x0),
    h_diffq32_etaP(0x0),
    h_diffq32_etaN(0x0),
    h_diffAq32_eta(0x0),
    h_diffAq32_etaP(0x0),
    h_diffAq32_etaN(0x0),
    h_diffq42_pt(0x0),
    h_diffq42_ptP(0x0),
    h_diffq42_ptN(0x0),
    h_diffAq42_pt(0x0),
    h_diffAq42_ptP(0x0),
    h_diffAq42_ptN(0x0),
    h_diffq42_eta(0x0),
    h_diffq42_etaP(0x0),
    h_diffq42_etaN(0x0),
    h_diffAq42_eta(0x0),
    h_diffAq42_etaP(0x0),
    h_diffAq42_etaN(0x0),
    
    h_AT_X_deta(0x0),
    h_AT_X_detaP(0x0),
    h_AT_X_detaN(0x0),
    h_diffq22_X_deta(0x0),
    h_diffq22_X_detaP(0x0),
    h_diffq22_X_detaN(0x0),
    h_ATdiffq22_X_deta(0x0),
    h_ATdiffq22_X_detaP(0x0),
    h_ATdiffq22_X_detaN(0x0),
    h_diffq32_X_deta(0x0),
    h_diffq32_X_detaP(0x0),
    h_diffq32_X_detaN(0x0),
    h_ATdiffq32_X_deta(0x0),
    h_ATdiffq32_X_detaP(0x0),
    h_ATdiffq32_X_detaN(0x0),
    h_diffq42_X_deta(0x0),
    h_diffq42_X_detaP(0x0),
    h_diffq42_X_detaN(0x0),
    h_ATdiffq42_X_deta(0x0),
    h_ATdiffq42_X_detaP(0x0),
    h_ATdiffq42_X_detaN(0x0),
    
    h_AT_X_deta_F1(0x0),
    h_AT_X_detaP_F1(0x0),
    h_AT_X_detaN_F1(0x0),
    h_diffq22_X_deta_F1(0x0),
    h_diffq22_X_detaP_F1(0x0),
    h_diffq22_X_detaN_F1(0x0),
    h_ATdiffq22_X_deta_F1(0x0),
    h_ATdiffq22_X_detaP_F1(0x0),
    h_ATdiffq22_X_detaN_F1(0x0),
    h_diffq32_X_deta_F1(0x0),
    h_diffq32_X_detaP_F1(0x0),
    h_diffq32_X_detaN_F1(0x0),
    h_ATdiffq32_X_deta_F1(0x0),
    h_ATdiffq32_X_detaP_F1(0x0),
    h_ATdiffq32_X_detaN_F1(0x0),
    h_diffq42_X_deta_F1(0x0),
    h_diffq42_X_detaP_F1(0x0),
    h_diffq42_X_detaN_F1(0x0),
    h_ATdiffq42_X_deta_F1(0x0),
    h_ATdiffq42_X_detaP_F1(0x0),
    h_ATdiffq42_X_detaN_F1(0x0),
    
    h_AT_X_deta_F3(0x0),
    h_AT_X_detaP_F3(0x0),
    h_AT_X_detaN_F3(0x0),
    h_diffq22_X_deta_F3(0x0),
    h_diffq22_X_detaP_F3(0x0),
    h_diffq22_X_detaN_F3(0x0),
    h_ATdiffq22_X_deta_F3(0x0),
    h_ATdiffq22_X_detaP_F3(0x0),
    h_ATdiffq22_X_detaN_F3(0x0),
    h_diffq32_X_deta_F3(0x0),
    h_diffq32_X_detaP_F3(0x0),
    h_diffq32_X_detaN_F3(0x0),
    h_ATdiffq32_X_deta_F3(0x0),
    h_ATdiffq32_X_detaP_F3(0x0),
    h_ATdiffq32_X_detaN_F3(0x0),
    h_diffq42_X_deta_F3(0x0),
    h_diffq42_X_detaP_F3(0x0),
    h_diffq42_X_detaN_F3(0x0),
    h_ATdiffq42_X_deta_F3(0x0),
    h_ATdiffq42_X_detaP_F3(0x0),
    h_ATdiffq42_X_detaN_F3(0x0),
    
    h_ATXdiffq22_X_detaP(0x0),
    h_ATXdiffq22_X_detaN(0x0),
    h_ATYdiffq22_X_detaP(0x0),
    h_ATYdiffq22_X_detaN(0x0),
    h_ATZdiffq22_X_detaP(0x0),
    h_ATZdiffq22_X_detaN(0x0),
    
    h_ATXdiffq32_X_detaP(0x0),
    h_ATXdiffq32_X_detaN(0x0),
    h_ATYdiffq32_X_detaP(0x0),
    h_ATYdiffq32_X_detaN(0x0),
    h_ATZdiffq32_X_detaP(0x0),
    h_ATZdiffq32_X_detaN(0x0),
    
    h_ATXdiffq42_X_detaP(0x0),
    h_ATXdiffq42_X_detaN(0x0),
    h_ATYdiffq42_X_detaP(0x0),
    h_ATYdiffq42_X_detaN(0x0),
    h_ATZdiffq42_X_detaP(0x0),
    h_ATZdiffq42_X_detaN(0x0),
    
    h_AT_S_deta(0x0),
    h_AT_S_detaP(0x0),
    h_AT_S_detaN(0x0),
    h_diffq22_S_deta(0x0),
    h_diffq22_S_detaP(0x0),
    h_diffq22_S_detaN(0x0),
    h_ATdiffq22_S_deta(0x0),
    h_ATdiffq22_S_detaP(0x0),
    h_ATdiffq22_S_detaN(0x0),
    h_diffq32_S_deta(0x0),
    h_diffq32_S_detaP(0x0),
    h_diffq32_S_detaN(0x0),
    h_ATdiffq32_S_deta(0x0),
    h_ATdiffq32_S_detaP(0x0),
    h_ATdiffq32_S_detaN(0x0),
    h_diffq42_S_deta(0x0),
    h_diffq42_S_detaP(0x0),
    h_diffq42_S_detaN(0x0),
    h_ATdiffq42_S_deta(0x0),
    h_ATdiffq42_S_detaP(0x0),
    h_ATdiffq42_S_detaN(0x0),
    
    h_AT_X_dpt(0x0),
    h_AT_X_dptP(0x0),
    h_AT_X_dptN(0x0),
    h_diffq22_X_dpt(0x0),
    h_diffq22_X_dptP(0x0),
    h_diffq22_X_dptN(0x0),
    h_ATdiffq22_X_dpt(0x0),
    h_ATdiffq22_X_dptP(0x0),
    h_ATdiffq22_X_dptN(0x0),
    h_diffq32_X_dpt(0x0),
    h_diffq32_X_dptP(0x0),
    h_diffq32_X_dptN(0x0),
    h_ATdiffq32_X_dpt(0x0),
    h_ATdiffq32_X_dptP(0x0),
    h_ATdiffq32_X_dptN(0x0),
    h_diffq42_X_dpt(0x0),
    h_diffq42_X_dptP(0x0),
    h_diffq42_X_dptN(0x0),
    h_ATdiffq42_X_dpt(0x0),
    h_ATdiffq42_X_dptP(0x0),
    h_ATdiffq42_X_dptN(0x0),
    
    h_a1q22_cent(0x0),
    h_a1q22_centP(0x0),
    h_a1q22_centN(0x0),
    
    h_a2q22_cent(0x0),
    h_a2q22_centP(0x0),
    h_a2q22_centN(0x0),
    
    h_rAq22_X1_cent(0x0),
    h_rAq22_X1_centP(0x0),
    h_rAq22_X1_centN(0x0),
    h_rAq22_X2_cent(0x0),
    h_rAq22_X2_centP(0x0),
    h_rAq22_X2_centN(0x0),
    h_rAq22_X3_cent(0x0),
    h_rAq22_X3_centP(0x0),
    h_rAq22_X3_centN(0x0),
    h_rAq22_X4_cent(0x0),
    h_rAq22_X4_centP(0x0),
    h_rAq22_X4_centN(0x0)
{
  // Default class constructor

  cout<<"Default class constructor called.  Prepare for fun analysis time!"<<endl;
  this->SetParameters();

}



//--------------------------------------------------------------------------------------
AliAnalysisTaskCMEv2A::AliAnalysisTaskCMEv2A(const char *name) : AliAnalysisTaskSE(name),
    debug(0),
    doMC(false),
    trigger(AliVEvent::kMB),
    dopupcut(true),
    doeffcorr(true),
    centhandle(1),
    fbit(128),
    zvtxcut(10.0),
    centcut(5.0),
    nclscut(60),
    dcacutz(5.0),
    dcacutxy(5.0),
    dodcacuts(false),
    outeta(0.8),
    ineta(0.5),
    excleta(0.0),
    ptmin(0.2),
    ptmax(5.0),
    doacuts(false),
    nspid(2.0),
    cbinlo(3),
    cbinhi(4),
    donested(true),
    dopaircut(true),
    centlo(20.0),
    centhi(60.0),
    fOutputList(0x0),     
    fHistDCACOVerrors(0x0),
    fHistPt(0x0),
    fHistPhi(0x0),
    fHistEta(0x0),
    fHistCharge(0x0),
    fHistTPCncls(0x0),
    fHistDedx(0x0),
    fHistDCAxy(0x0),
    fHistDCAz(0x0),
    fHistDCAxyAfter(0x0),
    fHistDCAzAfter(0x0),
    fHistPosPt(0x0),
    fHistPosPhi(0x0),
    fHistPosEta(0x0),
    fHistNegPt(0x0),
    fHistNegPhi(0x0),
    fHistNegEta(0x0),
    fHistPtPion(0x0),
    fHistPtPionH(0x0),
    fHistPtPionL(0x0),
    fHistNsigmaPion(0x0),
    fHistPtProt(0x0),
    fHistPtProtH(0x0),
    fHistPtProtL(0x0),
    fHistNsigmaProt(0x0),
    fHistPtMC(0x0),
    fHistPhiMC(0x0),
    fHistEtaMC(0x0),
    fHistChargeMC(0x0),
    fHistTPCnclsMC(0x0),
    fHistDedxMC(0x0),
    fHistDCAxyMC(0x0),
    fHistDCAzMC(0x0),
    fHistPlaneV0h2(0x0),
    fHistPlaneV0Ah2(0x0),
    fHistPlaneV0Ch2(0x0),
    fHistPlaneV0ACDCh2(0x0),
    fHistZVtx(0x0),
    fHistZVtxD(0x0),
    fHistZVtxMC(0x0),
    fHistZVtxDiffMC(0x0),
    fHistCentTRK(0x0),
    fHistCentV0M(0x0),
    fHistCentDIFF(0x0),
    fHistCentDIAG(0x0),
    fHistCtrkDIAG(0x0),
    fHistVtxRDIAG(0x0),
    fHistVtxSDIAG(0x0),
    fHistVtxRDIBG(0x0),
    fHistVtxSDIBG(0x0),
    fHistCentTRKAVEkMB(0x0),
    fHistCentV0MAVEkMB(0x0),
    fHistCentTRKAVEkCentral(0x0),
    fHistCentV0MAVEkCentral(0x0),
    fHistCentTRKAVEkSemiCentral(0x0),
    fHistCentV0MAVEkSemiCentral(0x0),
    fHistCentTRKAVEkA3(0x0),
    fHistCentV0MAVEkA3(0x0),
    fHistCentTRKAVEkSel(0x0),
    fHistCentV0MAVEkSel(0x0),
    
    h_MCq22_cent(0x0),
    h_MCq22_centP(0x0),
    h_MCq22_centN(0x0),
    h_MCAq22_cent(0x0),
    h_MCAq22_centP(0x0),
    h_MCAq22_centN(0x0),
    h_MCq22gap0_cent(0x0),
    h_MCq22gap0_centP(0x0),
    h_MCq22gap0_centN(0x0),
    h_MCq22gap1_cent(0x0),
    h_MCq22gap1_centP(0x0),
    h_MCq22gap1_centN(0x0),
    
    hMC_AT_X_deta(0x0),
    hMC_AT_X_detaP(0x0),
    hMC_AT_X_detaN(0x0),
    hMC_diffq22_X_deta(0x0),
    hMC_diffq22_X_detaP(0x0),
    hMC_diffq22_X_detaN(0x0),
    hMC_ATdiffq22_X_deta(0x0),
    hMC_ATdiffq22_X_detaP(0x0),
    hMC_ATdiffq22_X_detaN(0x0),
    hMC_diffq32_X_deta(0x0),
    hMC_diffq32_X_detaP(0x0),
    hMC_diffq32_X_detaN(0x0),
    hMC_ATdiffq32_X_deta(0x0),
    hMC_ATdiffq32_X_detaP(0x0),
    hMC_ATdiffq32_X_detaN(0x0),
    hMC_diffq42_X_deta(0x0),
    hMC_diffq42_X_detaP(0x0),
    hMC_diffq42_X_detaN(0x0),
    hMC_ATdiffq42_X_deta(0x0),
    hMC_ATdiffq42_X_detaP(0x0),
    hMC_ATdiffq42_X_detaN(0x0),
    
    fProfMeanChargePt(0x0),
    fProfMeanChargeEta(0x0),
    
    h_X22_cent(0x0),
    h_X22_centP(0x0),
    h_X22_centN(0x0),
    h_Y22_cent(0x0),
    h_Y22_centP(0x0),
    h_Y22_centN(0x0),
    h_X22_lo_cent(0x0),
    h_X22_lo_centP(0x0),
    h_X22_lo_centN(0x0),
    h_Y22_lo_cent(0x0),
    h_Y22_lo_centP(0x0),
    h_Y22_lo_centN(0x0),
    h_X22_li_cent(0x0),
    h_X22_li_centP(0x0),
    h_X22_li_centN(0x0),
    h_Y22_li_cent(0x0),
    h_Y22_li_centP(0x0),
    h_Y22_li_centN(0x0),
    h_X22_ri_cent(0x0),
    h_X22_ri_centP(0x0),
    h_X22_ri_centN(0x0),
    h_Y22_ri_cent(0x0),
    h_Y22_ri_centP(0x0),
    h_Y22_ri_centN(0x0),
    h_X22_ro_cent(0x0),
    h_X22_ro_centP(0x0),
    h_X22_ro_centN(0x0),
    h_Y22_ro_cent(0x0),
    h_Y22_ro_centP(0x0),
    h_Y22_ro_centN(0x0),
    
    h_q22_cent(0x0),
    h_q22_centP(0x0),
    h_q22_centN(0x0),
    h_q22_centPP(0x0),
    h_q22_centNN(0x0),
    h_q23_cent(0x0),
    h_q23_centP(0x0),
    h_q23_centN(0x0),
    h_q24_cent(0x0),
    h_q24_centP(0x0),
    h_q24_centN(0x0),
    h_q22gap0_cent(0x0),
    h_q22gap0_centP(0x0),
    h_q22gap0_centN(0x0),
    h_q22gap1_cent(0x0),
    h_q22gap1_centP(0x0),
    h_q22gap1_centN(0x0),
    
    h_q32_cent(0x0),
    h_q32_centP(0x0),
    h_q32_centN(0x0),
    h_q32_centPP(0x0),
    h_q32_centNN(0x0),
    h_q32gap0_cent(0x0),
    h_q32gap0_centP(0x0),
    h_q32gap0_centN(0x0),
    h_q32gap1_cent(0x0),
    h_q32gap1_centP(0x0),
    h_q32gap1_centN(0x0),
    
    h_q42_cent(0x0),
    h_q42_centP(0x0),
    h_q42_centN(0x0),
    h_q42_centPP(0x0),
    h_q42_centNN(0x0),
    h_q42gap0_cent(0x0),
    h_q42gap0_centP(0x0),
    h_q42gap0_centN(0x0),
    h_q42gap1_cent(0x0),
    h_q42gap1_centP(0x0),
    h_q42gap1_centN(0x0),
    
    //TProfile *h_q22qasym_cent[10];
    //TProfile *h_q22qasymP_cent[10];
    //TProfile *h_q22qasymN_cent[10];
    //TProfile *h_q22qasym_gap0_cent[10];
    //TProfile *h_q22qasymP_gap0_cent[10];
    //TProfile *h_q22qasymN_gap0_cent[10];
    //TProfile *h_q22qasym_gap1_cent[10];
    //TProfile *h_q22qasymP_gap1_cent[10];
    //TProfile *h_q22qasymN_gap1_cent[10];
    //TProfile *h_q24qasym_cent[10];
    //TProfile *h_q24qasymP_cent[10];
    //TProfile *h_q24qasymN_cent[10];
    //TProfile *h_q32qasym_cent[10];
    //TProfile *h_q32qasymP_cent[10];
    //TProfile *h_q32qasymN_cent[10];
    //TProfile *h_q32qasym_gap0_cent[10];
    //TProfile *h_q32qasymP_gap0_cent[10];
    //TProfile *h_q32qasymN_gap0_cent[10];
    //TProfile *h_q32qasym_gap1_cent[10];
    //TProfile *h_q32qasymP_gap1_cent[10];
    //TProfile *h_q32qasymN_gap1_cent[10];
    //TProfile *h_q34qasym_cent[10];
    //TProfile *h_q34qasymP_cent[10];
    //TProfile *h_q34qasymN_cent[10];
    //TProfile *h_q42qasym_cent[10];
    //TProfile *h_q42qasymP_cent[10];
    //TProfile *h_q42qasymN_cent[10];
    //TProfile *h_q42qasym_gap0_cent[10];
    //TProfile *h_q42qasymP_gap0_cent[10];
    //TProfile *h_q42qasymN_gap0_cent[10];
    //TProfile *h_q42qasym_gap1_cent[10];
    //TProfile *h_q42qasymP_gap1_cent[10];
    //TProfile *h_q42qasymN_gap1_cent[10];
    //TProfile *h_q44qasym_cent[10];
    //TProfile *h_q44qasymP_cent[10];
    //TProfile *h_q44qasymN_cent[10];
    
    h_eta_pos_F1(0x0),
    h_eta_pos_F3(0x0),
    h_eta_neg_F1(0x0),
    h_eta_neg_F3(0x0),
    h_cut_eta_pos_F1(0x0),
    h_cut_eta_pos_F3(0x0),
    h_cut_eta_neg_F1(0x0),
    h_cut_eta_neg_F3(0x0),
    
    h_AT_eta(0x0),
    h_AT_etaF1(0x0),
    h_AT_etaF3(0x0),
    h_AT_etaMC(0x0),
    h2_AT_eta(0x0),
    h2_AT_etaF1(0x0),
    h2_AT_etaF3(0x0),
    h2_AT_etaMC(0x0),
    
    h_AT_cut_eta(0x0),
    h_AT_cut_etaF1(0x0),
    h_AT_cut_etaF3(0x0),
    h_AT_cut_etaMC(0x0),
    h2_AT_cut_eta(0x0),
    h2_AT_cut_etaF1(0x0),
    h2_AT_cut_etaF3(0x0),
    h2_AT_cut_etaMC(0x0),
    
    h_A_cent(0x0),
    h_rA_cent(0x0),
    h_r1A_cent(0x0),
    h_r2A_cent(0x0),
    h_A_centF1(0x0),
    h_A_centF3(0x0),
    h_A_centMC(0x0),
    h2_A_cent(0x0),
    h2_rA_cent(0x0),
    h2_r1A_cent(0x0),
    h2_r2A_cent(0x0),
    h2_A_centF1(0x0),
    h2_A_centF3(0x0),
    h2_A_centMC(0x0),
    
    h_Aq22_cent(0x0),
    h_Aq22_centP(0x0),
    h_Aq22_centN(0x0),
    h_Aq22_centPP(0x0),
    h_Aq22_centNN(0x0),
    h_Aq22_SL_cent(0x0),
    h_Aq22_SL_centP(0x0),
    h_Aq22_SL_centN(0x0),
    h_Aq22_SLL_cent(0x0),
    h_Aq22_SLL_centP(0x0),
    h_Aq22_SLL_centN(0x0),
    h_Aq22_SLR_cent(0x0),
    h_Aq22_SLR_centP(0x0),
    h_Aq22_SLR_centN(0x0),
    h_Aq32_cent(0x0),
    h_Aq32_centP(0x0),
    h_Aq32_centN(0x0),
    h_Aq32_centPP(0x0),
    h_Aq32_centNN(0x0),
    h_Aq32_SL_cent(0x0),
    h_Aq32_SL_centP(0x0),
    h_Aq32_SL_centN(0x0),
    h_Aq32_SLL_cent(0x0),
    h_Aq32_SLL_centP(0x0),
    h_Aq32_SLL_centN(0x0),
    h_Aq32_SLR_cent(0x0),
    h_Aq32_SLR_centP(0x0),
    h_Aq32_SLR_centN(0x0),
    h_Aq42_cent(0x0),
    h_Aq42_centP(0x0),
    h_Aq42_centN(0x0),
    h_Aq42_centPP(0x0),
    h_Aq42_centNN(0x0),
    h_Aq42_SL_cent(0x0),
    h_Aq42_SL_centP(0x0),
    h_Aq42_SL_centN(0x0),
    h_Aq42_SLL_cent(0x0),
    h_Aq42_SLL_centP(0x0),
    h_Aq42_SLL_centN(0x0),
    h_Aq42_SLR_cent(0x0),
    h_Aq42_SLR_centP(0x0),
    h_Aq42_SLR_centN(0x0),
    
    h_diffq22_pt(0x0),
    h_diffq22_ptP(0x0),
    h_diffq22_ptN(0x0),
    h_diffAq22_pt(0x0),
    h_diffAq22_ptP(0x0),
    h_diffAq22_ptN(0x0),
    h_diffq22_eta(0x0),
    h_diffq22_etaP(0x0),
    h_diffq22_etaN(0x0),
    h_diffAq22_eta(0x0),
    h_diffAq22_etaP(0x0),
    h_diffAq22_etaN(0x0),
    h_diffq32_pt(0x0),
    h_diffq32_ptP(0x0),
    h_diffq32_ptN(0x0),
    h_diffAq32_pt(0x0),
    h_diffAq32_ptP(0x0),
    h_diffAq32_ptN(0x0),
    h_diffq32_eta(0x0),
    h_diffq32_etaP(0x0),
    h_diffq32_etaN(0x0),
    h_diffAq32_eta(0x0),
    h_diffAq32_etaP(0x0),
    h_diffAq32_etaN(0x0),
    h_diffq42_pt(0x0),
    h_diffq42_ptP(0x0),
    h_diffq42_ptN(0x0),
    h_diffAq42_pt(0x0),
    h_diffAq42_ptP(0x0),
    h_diffAq42_ptN(0x0),
    h_diffq42_eta(0x0),
    h_diffq42_etaP(0x0),
    h_diffq42_etaN(0x0),
    h_diffAq42_eta(0x0),
    h_diffAq42_etaP(0x0),
    h_diffAq42_etaN(0x0),
    
    h_AT_X_deta(0x0),
    h_AT_X_detaP(0x0),
    h_AT_X_detaN(0x0),
    h_diffq22_X_deta(0x0),
    h_diffq22_X_detaP(0x0),
    h_diffq22_X_detaN(0x0),
    h_ATdiffq22_X_deta(0x0),
    h_ATdiffq22_X_detaP(0x0),
    h_ATdiffq22_X_detaN(0x0),
    h_diffq32_X_deta(0x0),
    h_diffq32_X_detaP(0x0),
    h_diffq32_X_detaN(0x0),
    h_ATdiffq32_X_deta(0x0),
    h_ATdiffq32_X_detaP(0x0),
    h_ATdiffq32_X_detaN(0x0),
    h_diffq42_X_deta(0x0),
    h_diffq42_X_detaP(0x0),
    h_diffq42_X_detaN(0x0),
    h_ATdiffq42_X_deta(0x0),
    h_ATdiffq42_X_detaP(0x0),
    h_ATdiffq42_X_detaN(0x0),
    
    h_AT_X_deta_F1(0x0),
    h_AT_X_detaP_F1(0x0),
    h_AT_X_detaN_F1(0x0),
    h_diffq22_X_deta_F1(0x0),
    h_diffq22_X_detaP_F1(0x0),
    h_diffq22_X_detaN_F1(0x0),
    h_ATdiffq22_X_deta_F1(0x0),
    h_ATdiffq22_X_detaP_F1(0x0),
    h_ATdiffq22_X_detaN_F1(0x0),
    h_diffq32_X_deta_F1(0x0),
    h_diffq32_X_detaP_F1(0x0),
    h_diffq32_X_detaN_F1(0x0),
    h_ATdiffq32_X_deta_F1(0x0),
    h_ATdiffq32_X_detaP_F1(0x0),
    h_ATdiffq32_X_detaN_F1(0x0),
    h_diffq42_X_deta_F1(0x0),
    h_diffq42_X_detaP_F1(0x0),
    h_diffq42_X_detaN_F1(0x0),
    h_ATdiffq42_X_deta_F1(0x0),
    h_ATdiffq42_X_detaP_F1(0x0),
    h_ATdiffq42_X_detaN_F1(0x0),
    
    h_AT_X_deta_F3(0x0),
    h_AT_X_detaP_F3(0x0),
    h_AT_X_detaN_F3(0x0),
    h_diffq22_X_deta_F3(0x0),
    h_diffq22_X_detaP_F3(0x0),
    h_diffq22_X_detaN_F3(0x0),
    h_ATdiffq22_X_deta_F3(0x0),
    h_ATdiffq22_X_detaP_F3(0x0),
    h_ATdiffq22_X_detaN_F3(0x0),
    h_diffq32_X_deta_F3(0x0),
    h_diffq32_X_detaP_F3(0x0),
    h_diffq32_X_detaN_F3(0x0),
    h_ATdiffq32_X_deta_F3(0x0),
    h_ATdiffq32_X_detaP_F3(0x0),
    h_ATdiffq32_X_detaN_F3(0x0),
    h_diffq42_X_deta_F3(0x0),
    h_diffq42_X_detaP_F3(0x0),
    h_diffq42_X_detaN_F3(0x0),
    h_ATdiffq42_X_deta_F3(0x0),
    h_ATdiffq42_X_detaP_F3(0x0),
    h_ATdiffq42_X_detaN_F3(0x0),
    
    h_ATXdiffq22_X_detaP(0x0),
    h_ATXdiffq22_X_detaN(0x0),
    h_ATYdiffq22_X_detaP(0x0),
    h_ATYdiffq22_X_detaN(0x0),
    h_ATZdiffq22_X_detaP(0x0),
    h_ATZdiffq22_X_detaN(0x0),
    
    h_ATXdiffq32_X_detaP(0x0),
    h_ATXdiffq32_X_detaN(0x0),
    h_ATYdiffq32_X_detaP(0x0),
    h_ATYdiffq32_X_detaN(0x0),
    h_ATZdiffq32_X_detaP(0x0),
    h_ATZdiffq32_X_detaN(0x0),
    
    h_ATXdiffq42_X_detaP(0x0),
    h_ATXdiffq42_X_detaN(0x0),
    h_ATYdiffq42_X_detaP(0x0),
    h_ATYdiffq42_X_detaN(0x0),
    h_ATZdiffq42_X_detaP(0x0),
    h_ATZdiffq42_X_detaN(0x0),
    
    h_AT_S_deta(0x0),
    h_AT_S_detaP(0x0),
    h_AT_S_detaN(0x0),
    h_diffq22_S_deta(0x0),
    h_diffq22_S_detaP(0x0),
    h_diffq22_S_detaN(0x0),
    h_ATdiffq22_S_deta(0x0),
    h_ATdiffq22_S_detaP(0x0),
    h_ATdiffq22_S_detaN(0x0),
    h_diffq32_S_deta(0x0),
    h_diffq32_S_detaP(0x0),
    h_diffq32_S_detaN(0x0),
    h_ATdiffq32_S_deta(0x0),
    h_ATdiffq32_S_detaP(0x0),
    h_ATdiffq32_S_detaN(0x0),
    h_diffq42_S_deta(0x0),
    h_diffq42_S_detaP(0x0),
    h_diffq42_S_detaN(0x0),
    h_ATdiffq42_S_deta(0x0),
    h_ATdiffq42_S_detaP(0x0),
    h_ATdiffq42_S_detaN(0x0),
    
    h_AT_X_dpt(0x0),
    h_AT_X_dptP(0x0),
    h_AT_X_dptN(0x0),
    h_diffq22_X_dpt(0x0),
    h_diffq22_X_dptP(0x0),
    h_diffq22_X_dptN(0x0),
    h_ATdiffq22_X_dpt(0x0),
    h_ATdiffq22_X_dptP(0x0),
    h_ATdiffq22_X_dptN(0x0),
    h_diffq32_X_dpt(0x0),
    h_diffq32_X_dptP(0x0),
    h_diffq32_X_dptN(0x0),
    h_ATdiffq32_X_dpt(0x0),
    h_ATdiffq32_X_dptP(0x0),
    h_ATdiffq32_X_dptN(0x0),
    h_diffq42_X_dpt(0x0),
    h_diffq42_X_dptP(0x0),
    h_diffq42_X_dptN(0x0),
    h_ATdiffq42_X_dpt(0x0),
    h_ATdiffq42_X_dptP(0x0),
    h_ATdiffq42_X_dptN(0x0),
    
    h_a1q22_cent(0x0),
    h_a1q22_centP(0x0),
    h_a1q22_centN(0x0),
    
    h_a2q22_cent(0x0),
    h_a2q22_centP(0x0),
    h_a2q22_centN(0x0),
    
    h_rAq22_X1_cent(0x0),
    h_rAq22_X1_centP(0x0),
    h_rAq22_X1_centN(0x0),
    h_rAq22_X2_cent(0x0),
    h_rAq22_X2_centP(0x0),
    h_rAq22_X2_centN(0x0),
    h_rAq22_X3_cent(0x0),
    h_rAq22_X3_centP(0x0),
    h_rAq22_X3_centN(0x0),
    h_rAq22_X4_cent(0x0),
    h_rAq22_X4_centP(0x0),
    h_rAq22_X4_centN(0x0)
{
  // Class Constructor with name

  cout<<"User defined class constructor called.  Prepare for fun analysis time!"<<endl;
  this->SetParameters();

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0,TChain::Class());

  // Output slot #0 is reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  //DefineOutput(0,TTree::Class());
  DefineOutput(1,TList::Class());

  // user can change as needed

}



//---------------------------------------------------
AliAnalysisTaskCMEv2A::~AliAnalysisTaskCMEv2A()
{
  // Default class destructor

  cout<<"Default class destructor called.  Analysis fun time has ended"<<endl;

}



//---------------------------------------------------
void AliAnalysisTaskCMEv2A::UserCreateOutputObjects()
{

  // Create output objects, like histograms - called once

  cout<<"UserCreateOutputObjects called, now making histrograms and things"<<endl;



  // -------------------------- //
  // --- create output list --- //
  // -------------------------- //

  fOutputList = new TList();
  fOutputList->SetName(GetName());
  fOutputList->SetOwner(kTRUE);



  // ------------------------- //
  // --- create histograms --- //
  // ------------------------- //

  fHistDCACOVerrors = new TH1F("fHistDCACOVerrors","",100,0.0,1.0);
  fHistPt = new TH1F("fHistPt","",50,0.0,5.0);
  fHistPhi = new TH1F("fHistPhi","",63,0.0,6.3);
  fHistEta = new TH1F("fHistEta","",300,-1.5,1.5);
  fHistCharge = new TH1F("fHistCharge","",3,-1.5,1.5);
  fHistTPCncls = new TH1F("fHistTPCncls","",161,-0.5,160.5);
  fHistDedx = new TH1F("fHistDedx","",1000,-20,980);
  fHistDCAxy = new TH1F("fHistDCAxy","",350,-3.5,3.5);
  fHistDCAz = new TH1F("fHistDCAz","",350,-3.5,3.5);
  fHistDCAxyAfter = new TH1F("fHistDCAxyAfter","",350,-3.5,3.5);
  fHistDCAzAfter = new TH1F("fHistDCAzAfter","",350,-3.5,3.5);
  fOutputList->Add(fHistDCACOVerrors);
  fOutputList->Add(fHistPt);
  fOutputList->Add(fHistPhi);
  fOutputList->Add(fHistEta);
  fOutputList->Add(fHistCharge);
  fOutputList->Add(fHistTPCncls);
  fOutputList->Add(fHistDedx);
  fOutputList->Add(fHistDCAxy);
  fOutputList->Add(fHistDCAz);
  fOutputList->Add(fHistDCAxyAfter);
  fOutputList->Add(fHistDCAzAfter);

  fHistPosPt = new TH1F("fHistPosPt","",50,0.0,5.0);
  fHistPosPhi = new TH1F("fHistPosPhi","",63,0.0,6.3);
  fHistPosEta = new TH1F("fHistPosEta","",300,-1.5,1.5);
  fOutputList->Add(fHistPosPt);
  fOutputList->Add(fHistPosPhi);
  fOutputList->Add(fHistPosEta);

  fHistNegPt = new TH1F("fHistNegPt","",50,0.0,5.0);
  fHistNegPhi = new TH1F("fHistNegPhi","",63,0.0,6.3);
  fHistNegEta = new TH1F("fHistNegEta","",300,-1.5,1.5);
  fOutputList->Add(fHistNegPt);
  fOutputList->Add(fHistNegPhi);
  fOutputList->Add(fHistNegEta);



  fHistPtPion = new TH1F("fHistPtPion","",50,0.0,5.0);
  fHistPtPionH = new TH1F("fHistPtPionH","",50,0.0,5.0);
  fHistPtPionL = new TH1F("fHistPtPionL","",50,0.0,5.0);
  fHistNsigmaPion = new TH1F("fHistNsigmaPion","",100,-5.0,5.0);
  fHistPtProt = new TH1F("fHistPtProt","",50,0.0,5.0);
  fHistPtProtH = new TH1F("fHistPtProtH","",50,0.0,5.0);
  fHistPtProtL = new TH1F("fHistPtProtL","",50,0.0,5.0);
  fHistNsigmaProt = new TH1F("fHistNsigmaProt","",100,-5.0,5.0);
  fOutputList->Add(fHistPtPion);
  fOutputList->Add(fHistPtPionH);
  fOutputList->Add(fHistPtPionL);
  fOutputList->Add(fHistNsigmaPion);
  fOutputList->Add(fHistPtProt);
  fOutputList->Add(fHistPtProtH);
  fOutputList->Add(fHistPtProtL);
  fOutputList->Add(fHistNsigmaProt);


  fHistPtMC = new TH1F("fHistPtMC","",50,0.0,5.0);
  fHistPhiMC = new TH1F("fHistPhiMC","",63,0.0,6.3);
  fHistEtaMC = new TH1F("fHistEtaMC","",110,-1.1,1.1);
  fHistChargeMC = new TH1F("fHistChargeMC","",3,-1.5,1.5);
  fHistTPCnclsMC = new TH1F("fHistTPCnclsMC","",161,-0.5,160.5);
  fHistDedxMC = new TH1F("fHistDedxMC","",1000,-20,980);
  fHistDCAxyMC = new TH1F("fHistDCAxyMC","",350,-3.5,3.5);
  fHistDCAzMC = new TH1F("fHistDCAzMC","",350,-3.5,3.5);
  fOutputList->Add(fHistPtMC);
  fOutputList->Add(fHistPhiMC);
  fOutputList->Add(fHistEtaMC);
  fOutputList->Add(fHistChargeMC);
  fOutputList->Add(fHistTPCnclsMC);
  fOutputList->Add(fHistDedxMC);
  fOutputList->Add(fHistDCAxyMC);
  fOutputList->Add(fHistDCAzMC);




  fHistPlaneV0h2 = new TH1F("fHistPlaneV0h2","",640,-3.2,3.2);
  fHistPlaneV0Ah2 = new TH1F("fHistPlaneV0Ah2","",640,-3.2,3.2);
  fHistPlaneV0Ch2 = new TH1F("fHistPlaneV0Ch2","",640,-3.2,3.2);
  fHistPlaneV0ACDCh2 = new TH1F("fHistPlaneV0ACDCh2","",640,-3.2,3.2);
  fOutputList->Add(fHistPlaneV0h2);
  fOutputList->Add(fHistPlaneV0Ah2);
  fOutputList->Add(fHistPlaneV0Ch2);
  fOutputList->Add(fHistPlaneV0ACDCh2);


  fHistZVtx = new TH1F("fHistZVtx","",80,-20,20);
  fHistZVtxD = new TH1F("fHistZVtxD","",100,-5,5);
  fHistCentTRK = new TH1F("fHistCentTRK","",100,0,100);
  fHistCentV0M = new TH1F("fHistCentV0M","",100,0,100);
  fHistCentDIFF = new TH1F("fHistCentDIFF","",162,-40.5,40.5);
  fHistCentDIAG = new TH1F("fHistCentDIAG","",6,-3.5,2.5);
  fOutputList->Add(fHistZVtx);
  fOutputList->Add(fHistZVtxD);
  fOutputList->Add(fHistCentTRK);
  fOutputList->Add(fHistCentV0M);
  fOutputList->Add(fHistCentDIFF);
  fOutputList->Add(fHistCentDIAG);

  fHistZVtxMC = new TH1F("fHistZVtxMC","",80,-20,20);
  fHistZVtxDiffMC = new TH1F("fHistZVtxDiffMC","",100,-5,5);
  fOutputList->Add(fHistZVtxMC);
  fOutputList->Add(fHistZVtxDiffMC);

  fHistCtrkDIAG = new TH1F("fHistCtrkDIAG","",1001,-0.5,1000.5);
  fHistVtxRDIAG = new TH1F("fHistVtxRDIAG","",100,-5,5);
  fHistVtxSDIAG = new TH1F("fHistVtxSDIAG","",100,-5,5);
  fHistVtxRDIBG = new TH1F("fHistVtxRDIBG","",100,-5,5);
  fHistVtxSDIBG = new TH1F("fHistVtxSDIBG","",100,-5,5);
  fOutputList->Add(fHistCtrkDIAG);
  fOutputList->Add(fHistVtxRDIAG);
  fOutputList->Add(fHistVtxSDIAG);
  fOutputList->Add(fHistVtxRDIBG);
  fOutputList->Add(fHistVtxSDIBG);


  fHistCentTRKAVEkMB = new TH1F("fHistCentTRKAVEkMB","",100,0,100);
  fHistCentV0MAVEkMB = new TH1F("fHistCentV0MAVEkMB","",100,0,100);
  fHistCentTRKAVEkCentral = new TH1F("fHistCentTRKAVEkCentral","",100,0,100);
  fHistCentV0MAVEkCentral = new TH1F("fHistCentV0MAVEkCentral","",100,0,100);
  fHistCentTRKAVEkSemiCentral = new TH1F("fHistCentTRKAVEkSemiCentral","",100,0,100);
  fHistCentV0MAVEkSemiCentral = new TH1F("fHistCentV0MAVEkSemiCentral","",100,0,100);
  fHistCentTRKAVEkA3 = new TH1F("fHistCentTRKAVEkA3","",100,0,100);
  fHistCentV0MAVEkA3 = new TH1F("fHistCentV0MAVEkA3","",100,0,100);
  fHistCentTRKAVEkSel = new TH1F("fHistCentTRKAVEkSel","",100,0,100);
  fHistCentV0MAVEkSel = new TH1F("fHistCentV0MAVEkSel","",100,0,100);
  fOutputList->Add(fHistCentTRKAVEkMB);
  fOutputList->Add(fHistCentV0MAVEkMB);
  fOutputList->Add(fHistCentTRKAVEkCentral);
  fOutputList->Add(fHistCentV0MAVEkCentral);
  fOutputList->Add(fHistCentTRKAVEkSemiCentral);
  fOutputList->Add(fHistCentV0MAVEkSemiCentral);
  fOutputList->Add(fHistCentTRKAVEkA3);
  fOutputList->Add(fHistCentV0MAVEkA3);
  fOutputList->Add(fHistCentTRKAVEkSel);
  fOutputList->Add(fHistCentV0MAVEkSel);


  // -------------------------- //
  // --- flow and asymmetry --- //
  // -------------------------- //

  float tpmax = 1e10;
  float tpmin = -1e10;

  fProfMeanChargePt = new TProfile("fProfMeanChargePt","",50,0,5.0,tpmin,tpmax,"");
  fProfMeanChargeEta = new TProfile("fProfMeanChargeEta","",300,-1.5,1.5,tpmin,tpmax,"");

  h_eta_pos_F1 = new TH1F("h_eta_pos_F1","",300,-1.5,1.5);
  h_eta_neg_F1 = new TH1F("h_eta_neg_F1","",300,-1.5,1.5);
  h_eta_pos_F3 = new TH1F("h_eta_pos_F3","",300,-1.5,1.5);
  h_eta_neg_F3 = new TH1F("h_eta_neg_F3","",300,-1.5,1.5);
  fOutputList->Add(h_eta_pos_F1);
  fOutputList->Add(h_eta_neg_F1);
  fOutputList->Add(h_eta_pos_F3);
  fOutputList->Add(h_eta_neg_F3);
  h_cut_eta_pos_F1 = new TH1F("h_cut_eta_pos_F1","",300,-1.5,1.5);
  h_cut_eta_neg_F1 = new TH1F("h_cut_eta_neg_F1","",300,-1.5,1.5);
  h_cut_eta_pos_F3 = new TH1F("h_cut_eta_pos_F3","",300,-1.5,1.5);
  h_cut_eta_neg_F3 = new TH1F("h_cut_eta_neg_F3","",300,-1.5,1.5);
  fOutputList->Add(h_cut_eta_pos_F1);
  fOutputList->Add(h_cut_eta_neg_F1);
  fOutputList->Add(h_cut_eta_pos_F3);
  fOutputList->Add(h_cut_eta_neg_F3);

  // --- charge asymmetry (average charge in event)
  h_AT_eta = new TProfile("h_AT_eta","",150,-1.5,1.5,tpmin,tpmax,"");
  h_AT_etaF1 = new TProfile("h_AT_etaF1","",150,-1.5,1.5,tpmin,tpmax,"");
  h_AT_etaF3 = new TProfile("h_AT_etaF3","",150,-1.5,1.5,tpmin,tpmax,"");
  h_AT_etaMC = new TProfile("h_AT_etaMC","",150,-1.5,1.5,tpmin,tpmax,"");
  fOutputList->Add(h_AT_eta);
  fOutputList->Add(h_AT_etaF1);
  fOutputList->Add(h_AT_etaF3);
  fOutputList->Add(h_AT_etaMC);

  // --- charge asymmetry (average charge in event)
  h2_AT_eta = new TH2F("h2_AT_eta","",150,-1.5,1.5,100,-1.0,1.0);
  h2_AT_etaF1 = new TH2F("h2_AT_etaF1","",150,-1.5,1.5,100,-1.0,1.0);
  h2_AT_etaF3 = new TH2F("h2_AT_etaF3","",150,-1.5,1.5,100,-1.0,1.0);
  h2_AT_etaMC = new TH2F("h2_AT_etaMC","",150,-1.5,1.5,100,-1.0,1.0);
  fOutputList->Add(h2_AT_eta);
  fOutputList->Add(h2_AT_etaF1);
  fOutputList->Add(h2_AT_etaF3);
  fOutputList->Add(h2_AT_etaMC);

  // --- charge asymmetry (average charge in event)
  h_AT_cut_eta = new TProfile("h_AT_cut_eta","",150,-1.5,1.5,tpmin,tpmax,"");
  h_AT_cut_etaF1 = new TProfile("h_AT_cut_etaF1","",150,-1.5,1.5,tpmin,tpmax,"");
  h_AT_cut_etaF3 = new TProfile("h_AT_cut_etaF3","",150,-1.5,1.5,tpmin,tpmax,"");
  h_AT_cut_etaMC = new TProfile("h_AT_cut_etaMC","",150,-1.5,1.5,tpmin,tpmax,"");
  fOutputList->Add(h_AT_cut_eta);
  fOutputList->Add(h_AT_cut_etaF1);
  fOutputList->Add(h_AT_cut_etaF3);
  fOutputList->Add(h_AT_cut_etaMC);

  // --- charge asymmetry (average charge in event)
  h2_AT_cut_eta = new TH2F("h2_AT_cut_eta","",150,-1.5,1.5,100,-1.0,1.0);
  h2_AT_cut_etaF1 = new TH2F("h2_AT_cut_etaF1","",150,-1.5,1.5,100,-1.0,1.0);
  h2_AT_cut_etaF3 = new TH2F("h2_AT_cut_etaF3","",150,-1.5,1.5,100,-1.0,1.0);
  h2_AT_cut_etaMC = new TH2F("h2_AT_cut_etaMC","",150,-1.5,1.5,100,-1.0,1.0);
  fOutputList->Add(h2_AT_cut_eta);
  fOutputList->Add(h2_AT_cut_etaF1);
  fOutputList->Add(h2_AT_cut_etaF3);
  fOutputList->Add(h2_AT_cut_etaMC);

  // --- charge asymmetry (average charge in event)
  h_A_cent = new TProfile("h_A_cent","",100,0,100,tpmin,tpmax,"");
  h_rA_cent = new TProfile("h_rA_cent","",100,0,100,tpmin,tpmax,"");
  h_r1A_cent = new TProfile("h_r1A_cent","",100,0,100,tpmin,tpmax,"");
  h_r2A_cent = new TProfile("h_r2A_cent","",100,0,100,tpmin,tpmax,"");
  h_A_centF1 = new TProfile("h_A_centF1","",100,0,100,tpmin,tpmax,"");
  h_A_centF3 = new TProfile("h_A_centF3","",100,0,100,tpmin,tpmax,"");
  h_A_centMC = new TProfile("h_A_centMC","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_A_cent);
  fOutputList->Add(h_rA_cent);
  fOutputList->Add(h_r1A_cent);
  fOutputList->Add(h_r2A_cent);
  fOutputList->Add(h_A_centF1);
  fOutputList->Add(h_A_centF3);
  fOutputList->Add(h_A_centMC);

  // --- charge asymmetry (average charge in event)
  h2_A_cent = new TH2F("h2_A_cent","",100,0,100,100,-1.0,1.0);
  h2_rA_cent = new TH2F("h2_rA_cent","",100,0,100,100,-1.0,1.0);
  h2_r1A_cent = new TH2F("h2_r1A_cent","",100,0,100,100,-1.0,1.0);
  h2_r2A_cent = new TH2F("h2_r2A_cent","",100,0,100,100,-1.0,1.0);
  h2_A_centF1 = new TH2F("h2_A_centF1","",100,0,100,100,-1.0,1.0);
  h2_A_centF3 = new TH2F("h2_A_centF3","",100,0,100,100,-1.0,1.0);
  h2_A_centMC = new TH2F("h2_A_centMC","",100,0,100,100,-1.0,1.0);
  fOutputList->Add(h2_A_cent);
  fOutputList->Add(h2_rA_cent);
  fOutputList->Add(h2_r1A_cent);
  fOutputList->Add(h2_r2A_cent);
  fOutputList->Add(h2_A_centF1);
  fOutputList->Add(h2_A_centF3);
  fOutputList->Add(h2_A_centMC);

  // --- MONTE CARLO second harmonic vs centrality
  h_MCq22_cent = new TProfile("h_MCq22_cent","",100,0,100,tpmin,tpmax,"");
  h_MCq22_centP = new TProfile("h_MCq22_centP","",100,0,100,tpmin,tpmax,"");
  h_MCq22_centN = new TProfile("h_MCq22_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_MCq22_cent);
  fOutputList->Add(h_MCq22_centP);
  fOutputList->Add(h_MCq22_centN);
  h_MCAq22_cent = new TProfile("h_MCAq22_cent","",100,0,100,tpmin,tpmax,"");
  h_MCAq22_centP = new TProfile("h_MCAq22_centP","",100,0,100,tpmin,tpmax,"");
  h_MCAq22_centN = new TProfile("h_MCAq22_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_MCAq22_cent);
  fOutputList->Add(h_MCAq22_centP);
  fOutputList->Add(h_MCAq22_centN);
  h_MCq22gap0_cent = new TProfile("h_MCq22gap0_cent","",100,0,100,tpmin,tpmax,"");
  h_MCq22gap0_centP = new TProfile("h_MCq22gap0_centP","",100,0,100,tpmin,tpmax,"");
  h_MCq22gap0_centN = new TProfile("h_MCq22gap0_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_MCq22gap0_cent);
  fOutputList->Add(h_MCq22gap0_centP);
  fOutputList->Add(h_MCq22gap0_centN);
  h_MCq22gap1_cent = new TProfile("h_MCq22gap1_cent","",100,0,100,tpmin,tpmax,"");
  h_MCq22gap1_centP = new TProfile("h_MCq22gap1_centP","",100,0,100,tpmin,tpmax,"");
  h_MCq22gap1_centN = new TProfile("h_MCq22gap1_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_MCq22gap1_cent);
  fOutputList->Add(h_MCq22gap1_centP);
  fOutputList->Add(h_MCq22gap1_centN);

  hMC_AT_X_deta = new TProfile("hMC_AT_X_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  hMC_AT_X_detaP = new TProfile("hMC_AT_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  hMC_AT_X_detaN = new TProfile("hMC_AT_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(hMC_AT_X_deta);
  fOutputList->Add(hMC_AT_X_detaP);
  fOutputList->Add(hMC_AT_X_detaN);
  //hMC_diffq22_X_deta = new TProfile("hMC_diffq22_X_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  hMC_diffq22_X_detaP = new TProfile("hMC_diffq22_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  hMC_diffq22_X_detaN = new TProfile("hMC_diffq22_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  //hMC_ATdiffq22_X_deta = new TProfile("hMC_ATdiffq22_X_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  hMC_ATdiffq22_X_detaP = new TProfile("hMC_ATdiffq22_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  hMC_ATdiffq22_X_detaN = new TProfile("hMC_ATdiffq22_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(hMC_diffq22_X_detaP);
  fOutputList->Add(hMC_diffq22_X_detaN);
  fOutputList->Add(hMC_ATdiffq22_X_detaP);
  fOutputList->Add(hMC_ATdiffq22_X_detaN);
  //hMC_diffq32_X_deta = new TProfile("hMC_diffq32_X_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  hMC_diffq32_X_detaP = new TProfile("hMC_diffq32_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  hMC_diffq32_X_detaN = new TProfile("hMC_diffq32_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  //hMC_ATdiffq32_X_deta = new TProfile("hMC_ATdiffq32_X_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  hMC_ATdiffq32_X_detaP = new TProfile("hMC_ATdiffq32_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  hMC_ATdiffq32_X_detaN = new TProfile("hMC_ATdiffq32_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(hMC_diffq32_X_detaP);
  fOutputList->Add(hMC_diffq32_X_detaN);
  fOutputList->Add(hMC_ATdiffq32_X_detaP);
  fOutputList->Add(hMC_ATdiffq32_X_detaN);
  //hMC_diffq42_X_deta = new TProfile("hMC_diffq42_X_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  hMC_diffq42_X_detaP = new TProfile("hMC_diffq42_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  hMC_diffq42_X_detaN = new TProfile("hMC_diffq42_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  //hMC_ATdiffq42_X_deta = new TProfile("hMC_ATdiffq42_X_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  hMC_ATdiffq42_X_detaP = new TProfile("hMC_ATdiffq42_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  hMC_ATdiffq42_X_detaN = new TProfile("hMC_ATdiffq42_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(hMC_diffq42_X_detaP);
  fOutputList->Add(hMC_diffq42_X_detaN);
  fOutputList->Add(hMC_ATdiffq42_X_detaP);
  fOutputList->Add(hMC_ATdiffq42_X_detaN);


  // --- q-vector components to check recentering
  h_X22_cent = new TProfile("h_X22_cent","",100,0,100,tpmin,tpmax,"");
  h_X22_centP = new TProfile("h_X22_centP","",100,0,100,tpmin,tpmax,"");
  h_X22_centN = new TProfile("h_X22_centN","",100,0,100,tpmin,tpmax,"");
  h_Y22_cent = new TProfile("h_Y22_cent","",100,0,100,tpmin,tpmax,"");
  h_Y22_centP = new TProfile("h_Y22_centP","",100,0,100,tpmin,tpmax,"");
  h_Y22_centN = new TProfile("h_Y22_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_X22_cent);
  fOutputList->Add(h_X22_centP);
  fOutputList->Add(h_X22_centN);
  fOutputList->Add(h_Y22_cent);
  fOutputList->Add(h_Y22_centP);
  fOutputList->Add(h_Y22_centN);
  // --- outer left
  h_X22_lo_cent = new TProfile("h_X22_lo_cent","",100,0,100,tpmin,tpmax,"");
  h_X22_lo_centP = new TProfile("h_X22_lo_centP","",100,0,100,tpmin,tpmax,"");
  h_X22_lo_centN = new TProfile("h_X22_lo_centN","",100,0,100,tpmin,tpmax,"");
  h_Y22_lo_cent = new TProfile("h_Y22_lo_cent","",100,0,100,tpmin,tpmax,"");
  h_Y22_lo_centP = new TProfile("h_Y22_lo_centP","",100,0,100,tpmin,tpmax,"");
  h_Y22_lo_centN = new TProfile("h_Y22_lo_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_X22_lo_cent);
  fOutputList->Add(h_X22_lo_centP);
  fOutputList->Add(h_X22_lo_centN);
  fOutputList->Add(h_Y22_lo_cent);
  fOutputList->Add(h_Y22_lo_centP);
  fOutputList->Add(h_Y22_lo_centN);
  // --- inner left
  h_X22_li_cent = new TProfile("h_X22_li_cent","",100,0,100,tpmin,tpmax,"");
  h_X22_li_centP = new TProfile("h_X22_li_centP","",100,0,100,tpmin,tpmax,"");
  h_X22_li_centN = new TProfile("h_X22_li_centN","",100,0,100,tpmin,tpmax,"");
  h_Y22_li_cent = new TProfile("h_Y22_li_cent","",100,0,100,tpmin,tpmax,"");
  h_Y22_li_centP = new TProfile("h_Y22_li_centP","",100,0,100,tpmin,tpmax,"");
  h_Y22_li_centN = new TProfile("h_Y22_li_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_X22_li_cent);
  fOutputList->Add(h_X22_li_centP);
  fOutputList->Add(h_X22_li_centN);
  fOutputList->Add(h_Y22_li_cent);
  fOutputList->Add(h_Y22_li_centP);
  fOutputList->Add(h_Y22_li_centN);
  // --- inner right
  h_X22_ri_cent = new TProfile("h_X22_ri_cent","",100,0,100,tpmin,tpmax,"");
  h_X22_ri_centP = new TProfile("h_X22_ri_centP","",100,0,100,tpmin,tpmax,"");
  h_X22_ri_centN = new TProfile("h_X22_ri_centN","",100,0,100,tpmin,tpmax,"");
  h_Y22_ri_cent = new TProfile("h_Y22_ri_cent","",100,0,100,tpmin,tpmax,"");
  h_Y22_ri_centP = new TProfile("h_Y22_ri_centP","",100,0,100,tpmin,tpmax,"");
  h_Y22_ri_centN = new TProfile("h_Y22_ri_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_X22_ri_cent);
  fOutputList->Add(h_X22_ri_centP);
  fOutputList->Add(h_X22_ri_centN);
  fOutputList->Add(h_Y22_ri_cent);
  fOutputList->Add(h_Y22_ri_centP);
  fOutputList->Add(h_Y22_ri_centN);
  // --- outer right
  h_X22_ro_cent = new TProfile("h_X22_ro_cent","",100,0,100,tpmin,tpmax,"");
  h_X22_ro_centP = new TProfile("h_X22_ro_centP","",100,0,100,tpmin,tpmax,"");
  h_X22_ro_centN = new TProfile("h_X22_ro_centN","",100,0,100,tpmin,tpmax,"");
  h_Y22_ro_cent = new TProfile("h_Y22_ro_cent","",100,0,100,tpmin,tpmax,"");
  h_Y22_ro_centP = new TProfile("h_Y22_ro_centP","",100,0,100,tpmin,tpmax,"");
  h_Y22_ro_centN = new TProfile("h_Y22_ro_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_X22_ro_cent);
  fOutputList->Add(h_X22_ro_centP);
  fOutputList->Add(h_X22_ro_centN);
  fOutputList->Add(h_Y22_ro_cent);
  fOutputList->Add(h_Y22_ro_centP);
  fOutputList->Add(h_Y22_ro_centN);

  // --- second harmonic vs centrality
  h_q22_cent = new TProfile("h_q22_cent","",100,0,100,tpmin,tpmax,"");
  h_q22_centP = new TProfile("h_q22_centP","",100,0,100,tpmin,tpmax,"");
  h_q22_centN = new TProfile("h_q22_centN","",100,0,100,tpmin,tpmax,"");
  h_q22_centPP = new TProfile("h_q22_centPP","",100,0,100,tpmin,tpmax,"");
  h_q22_centNN = new TProfile("h_q22_centNN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_q22_cent);
  fOutputList->Add(h_q22_centP);
  fOutputList->Add(h_q22_centN);
  fOutputList->Add(h_q22_centPP);
  fOutputList->Add(h_q22_centNN);
  h_q23_cent = new TProfile("h_q23_cent","",100,0,100,tpmin,tpmax,"");
  h_q23_centP = new TProfile("h_q23_centP","",100,0,100,tpmin,tpmax,"");
  h_q23_centN = new TProfile("h_q23_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_q23_cent);
  fOutputList->Add(h_q23_centP);
  fOutputList->Add(h_q23_centN);
  h_q24_cent = new TProfile("h_q24_cent","",100,0,100,tpmin,tpmax,"");
  h_q24_centP = new TProfile("h_q24_centP","",100,0,100,tpmin,tpmax,"");
  h_q24_centN = new TProfile("h_q24_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_q24_cent);
  fOutputList->Add(h_q24_centP);
  fOutputList->Add(h_q24_centN);
  h_q22gap0_cent = new TProfile("h_q22gap0_cent","",100,0,100,tpmin,tpmax,"");
  h_q22gap0_centP = new TProfile("h_q22gap0_centP","",100,0,100,tpmin,tpmax,"");
  h_q22gap0_centN = new TProfile("h_q22gap0_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_q22gap0_cent);
  fOutputList->Add(h_q22gap0_centP);
  fOutputList->Add(h_q22gap0_centN);
  h_q22gap1_cent = new TProfile("h_q22gap1_cent","",100,0,100,tpmin,tpmax,"");
  h_q22gap1_centP = new TProfile("h_q22gap1_centP","",100,0,100,tpmin,tpmax,"");
  h_q22gap1_centN = new TProfile("h_q22gap1_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_q22gap1_cent);
  fOutputList->Add(h_q22gap1_centP);
  fOutputList->Add(h_q22gap1_centN);

  // --- third harmonic vs centrality
  h_q32_cent = new TProfile("h_q32_cent","",100,0,100,tpmin,tpmax,"");
  h_q32_centP = new TProfile("h_q32_centP","",100,0,100,tpmin,tpmax,"");
  h_q32_centN = new TProfile("h_q32_centN","",100,0,100,tpmin,tpmax,"");
  h_q32_centPP = new TProfile("h_q32_centPP","",100,0,100,tpmin,tpmax,"");
  h_q32_centNN = new TProfile("h_q32_centNN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_q32_cent);
  fOutputList->Add(h_q32_centP);
  fOutputList->Add(h_q32_centN);
  fOutputList->Add(h_q32_centPP);
  fOutputList->Add(h_q32_centNN);
  h_q32gap0_cent = new TProfile("h_q32gap0_cent","",100,0,100,tpmin,tpmax,"");
  h_q32gap0_centP = new TProfile("h_q32gap0_centP","",100,0,100,tpmin,tpmax,"");
  h_q32gap0_centN = new TProfile("h_q32gap0_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_q32gap0_cent);
  fOutputList->Add(h_q32gap0_centP);
  fOutputList->Add(h_q32gap0_centN);
  h_q32gap1_cent = new TProfile("h_q32gap1_cent","",100,0,100,tpmin,tpmax,"");
  h_q32gap1_centP = new TProfile("h_q32gap1_centP","",100,0,100,tpmin,tpmax,"");
  h_q32gap1_centN = new TProfile("h_q32gap1_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_q32gap1_cent);
  fOutputList->Add(h_q32gap1_centP);
  fOutputList->Add(h_q32gap1_centN);

  // --- fourth harmonic vs centrality
  h_q42_cent = new TProfile("h_q42_cent","",100,0,100,tpmin,tpmax,"");
  h_q42_centP = new TProfile("h_q42_centP","",100,0,100,tpmin,tpmax,"");
  h_q42_centN = new TProfile("h_q42_centN","",100,0,100,tpmin,tpmax,"");
  h_q42_centPP = new TProfile("h_q42_centPP","",100,0,100,tpmin,tpmax,"");
  h_q42_centNN = new TProfile("h_q42_centNN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_q42_cent);
  fOutputList->Add(h_q42_centP);
  fOutputList->Add(h_q42_centN);
  fOutputList->Add(h_q42_centPP);
  fOutputList->Add(h_q42_centNN);
  h_q42gap0_cent = new TProfile("h_q42gap0_cent","",100,0,100,tpmin,tpmax,"");
  h_q42gap0_centP = new TProfile("h_q42gap0_centP","",100,0,100,tpmin,tpmax,"");
  h_q42gap0_centN = new TProfile("h_q42gap0_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_q42gap0_cent);
  fOutputList->Add(h_q42gap0_centP);
  fOutputList->Add(h_q42gap0_centN);
  h_q42gap1_cent = new TProfile("h_q42gap1_cent","",100,0,100,tpmin,tpmax,"");
  h_q42gap1_centP = new TProfile("h_q42gap1_centP","",100,0,100,tpmin,tpmax,"");
  h_q42gap1_centN = new TProfile("h_q42gap1_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_q42gap1_cent);
  fOutputList->Add(h_q42gap1_centP);
  fOutputList->Add(h_q42gap1_centN);

  // harmonics vs charge asymmetry in centrality bins
  for(int i=0; i<10; i++)
    {
      // --- second harmonic
      h_q22qasym_cent[i] = new TProfile(Form("h_q22qasym_cent_%d",i),Form("h_q22qasym_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q22qasymP_cent[i] = new TProfile(Form("h_q22qasymP_cent_%d",i),Form("h_q22qasymP_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q22qasymN_cent[i] = new TProfile(Form("h_q22qasymN_cent_%d",i),Form("h_q22qasymN_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q22qasym_gap0_cent[i] = new TProfile(Form("h_q22qasym_gap0_cent_%d",i),Form("h_q22qasym_gap0_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q22qasymP_gap0_cent[i] = new TProfile(Form("h_q22qasymP_gap0_cent_%d",i),Form("h_q22qasymP_gap0_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q22qasymN_gap0_cent[i] = new TProfile(Form("h_q22qasymN_gap0_cent_%d",i),Form("h_q22qasymN_gap0_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q22qasym_gap1_cent[i] = new TProfile(Form("h_q22qasym_gap1_cent_%d",i),Form("h_q22qasym_gap1_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q22qasymP_gap1_cent[i] = new TProfile(Form("h_q22qasymP_gap1_cent_%d",i),Form("h_q22qasymP_gap1_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q22qasymN_gap1_cent[i] = new TProfile(Form("h_q22qasymN_gap1_cent_%d",i),Form("h_q22qasymN_gap1_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q24qasym_cent[i] = new TProfile(Form("h_q24qasym_cent_%d",i),Form("h_q24qasym_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q24qasymP_cent[i] = new TProfile(Form("h_q24qasymP_cent_%d",i),Form("h_q24qasymP_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q24qasymN_cent[i] = new TProfile(Form("h_q24qasymN_cent_%d",i),Form("h_q24qasymN_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");

      // --- third harmonic
      h_q32qasym_cent[i] = new TProfile(Form("h_q32qasym_cent_%d",i),Form("h_q32qasym_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q32qasymP_cent[i] = new TProfile(Form("h_q32qasymP_cent_%d",i),Form("h_q32qasymP_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q32qasymN_cent[i] = new TProfile(Form("h_q32qasymN_cent_%d",i),Form("h_q32qasymN_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q32qasym_gap0_cent[i] = new TProfile(Form("h_q32qasym_gap0_cent_%d",i),Form("h_q32qasym_gap0_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q32qasymP_gap0_cent[i] = new TProfile(Form("h_q32qasymP_gap0_cent_%d",i),Form("h_q32qasymP_gap0_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q32qasymN_gap0_cent[i] = new TProfile(Form("h_q32qasymN_gap0_cent_%d",i),Form("h_q32qasymN_gap0_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q32qasym_gap1_cent[i] = new TProfile(Form("h_q32qasym_gap1_cent_%d",i),Form("h_q32qasym_gap1_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q32qasymP_gap1_cent[i] = new TProfile(Form("h_q32qasymP_gap1_cent_%d",i),Form("h_q32qasymP_gap1_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q32qasymN_gap1_cent[i] = new TProfile(Form("h_q32qasymN_gap1_cent_%d",i),Form("h_q32qasymN_gap1_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q34qasym_cent[i] = new TProfile(Form("h_q34qasym_cent_%d",i),Form("h_q34qasym_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q34qasymP_cent[i] = new TProfile(Form("h_q34qasymP_cent_%d",i),Form("h_q34qasymP_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q34qasymN_cent[i] = new TProfile(Form("h_q34qasymN_cent_%d",i),Form("h_q34qasymN_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");

      // --- fourth harmonic
      h_q42qasym_cent[i] = new TProfile(Form("h_q42qasym_cent_%d",i),Form("h_q42qasym_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q42qasymP_cent[i] = new TProfile(Form("h_q42qasymP_cent_%d",i),Form("h_q42qasymP_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q42qasymN_cent[i] = new TProfile(Form("h_q42qasymN_cent_%d",i),Form("h_q42qasymN_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q42qasym_gap0_cent[i] = new TProfile(Form("h_q42qasym_gap0_cent_%d",i),Form("h_q42qasym_gap0_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q42qasymP_gap0_cent[i] = new TProfile(Form("h_q42qasymP_gap0_cent_%d",i),Form("h_q42qasymP_gap0_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q42qasymN_gap0_cent[i] = new TProfile(Form("h_q42qasymN_gap0_cent_%d",i),Form("h_q42qasymN_gap0_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q42qasym_gap1_cent[i] = new TProfile(Form("h_q42qasym_gap1_cent_%d",i),Form("h_q42qasym_gap1_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q42qasymP_gap1_cent[i] = new TProfile(Form("h_q42qasymP_gap1_cent_%d",i),Form("h_q42qasymP_gap1_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q42qasymN_gap1_cent[i] = new TProfile(Form("h_q42qasymN_gap1_cent_%d",i),Form("h_q42qasymN_gap1_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q44qasym_cent[i] = new TProfile(Form("h_q44qasym_cent_%d",i),Form("h_q44qasym_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q44qasymP_cent[i] = new TProfile(Form("h_q44qasymP_cent_%d",i),Form("h_q44qasymP_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
      h_q44qasymN_cent[i] = new TProfile(Form("h_q44qasymN_cent_%d",i),Form("h_q44qasymN_cent_%d",i),220,-1.1,1.1,tpmin,tpmax,"");
    }


  // harmonics vs charge asymmetry in centrality bins
  for(int i=cbinlo; i<cbinhi; i++)
    {
      // --- second harmonic
      fOutputList->Add(h_q22qasym_cent[i]);
      fOutputList->Add(h_q22qasymP_cent[i]);
      fOutputList->Add(h_q22qasymN_cent[i]);
      fOutputList->Add(h_q22qasym_gap0_cent[i]);
      fOutputList->Add(h_q22qasymP_gap0_cent[i]);
      fOutputList->Add(h_q22qasymN_gap0_cent[i]);
      fOutputList->Add(h_q22qasym_gap1_cent[i]);
      fOutputList->Add(h_q22qasymP_gap1_cent[i]);
      fOutputList->Add(h_q22qasymN_gap1_cent[i]);
      fOutputList->Add(h_q24qasym_cent[i]);
      fOutputList->Add(h_q24qasymP_cent[i]);
      fOutputList->Add(h_q24qasymN_cent[i]);

      // --- third harmonic
      fOutputList->Add(h_q32qasym_cent[i]);
      fOutputList->Add(h_q32qasymP_cent[i]);
      fOutputList->Add(h_q32qasymN_cent[i]);
      fOutputList->Add(h_q32qasym_gap0_cent[i]);
      fOutputList->Add(h_q32qasymP_gap0_cent[i]);
      fOutputList->Add(h_q32qasymN_gap0_cent[i]);
      fOutputList->Add(h_q32qasym_gap1_cent[i]);
      fOutputList->Add(h_q32qasymP_gap1_cent[i]);
      fOutputList->Add(h_q32qasymN_gap1_cent[i]);
      fOutputList->Add(h_q34qasym_cent[i]);
      fOutputList->Add(h_q34qasymP_cent[i]);
      fOutputList->Add(h_q34qasymN_cent[i]);

      // --- fourth harmonic
      fOutputList->Add(h_q42qasym_cent[i]);
      fOutputList->Add(h_q42qasymP_cent[i]);
      fOutputList->Add(h_q42qasymN_cent[i]);
      fOutputList->Add(h_q42qasym_gap0_cent[i]);
      fOutputList->Add(h_q42qasymP_gap0_cent[i]);
      fOutputList->Add(h_q42qasymN_gap0_cent[i]);
      fOutputList->Add(h_q42qasym_gap1_cent[i]);
      fOutputList->Add(h_q42qasymP_gap1_cent[i]);
      fOutputList->Add(h_q42qasymN_gap1_cent[i]);
      fOutputList->Add(h_q44qasym_cent[i]);
      fOutputList->Add(h_q44qasymP_cent[i]);
      fOutputList->Add(h_q44qasymN_cent[i]);
    }



  // --------------------------------------- //
  // --- three-particle-like correlators --- //
  // --------------------------------------- //


  // --- second harmonic
  h_Aq22_cent = new TProfile("h_Aq22_cent","",100,0,100,tpmin,tpmax,"");
  h_Aq22_centP = new TProfile("h_Aq22_centP","",100,0,100,tpmin,tpmax,"");
  h_Aq22_centN = new TProfile("h_Aq22_centN","",100,0,100,tpmin,tpmax,"");
  h_Aq22_centPP = new TProfile("h_Aq22_centPP","",100,0,100,tpmin,tpmax,"");
  h_Aq22_centNN = new TProfile("h_Aq22_centNN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_Aq22_cent);
  fOutputList->Add(h_Aq22_centP);
  fOutputList->Add(h_Aq22_centN);
  h_Aq22_SL_cent = new TProfile("h_Aq22_SL_cent","",100,0,100,tpmin,tpmax,"");
  h_Aq22_SL_centP = new TProfile("h_Aq22_SL_centP","",100,0,100,tpmin,tpmax,"");
  h_Aq22_SL_centN = new TProfile("h_Aq22_SL_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_Aq22_SL_cent);
  fOutputList->Add(h_Aq22_SL_centP);
  fOutputList->Add(h_Aq22_SL_centN);
  h_Aq22_SLL_cent = new TProfile("h_Aq22_SLL_cent","",100,0,100,tpmin,tpmax,"");
  h_Aq22_SLL_centP = new TProfile("h_Aq22_SLL_centP","",100,0,100,tpmin,tpmax,"");
  h_Aq22_SLL_centN = new TProfile("h_Aq22_SLL_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_Aq22_SLL_cent);
  fOutputList->Add(h_Aq22_SLL_centP);
  fOutputList->Add(h_Aq22_SLL_centN);
  h_Aq22_SLR_cent = new TProfile("h_Aq22_SLR_cent","",100,0,100,tpmin,tpmax,"");
  h_Aq22_SLR_centP = new TProfile("h_Aq22_SLR_centP","",100,0,100,tpmin,tpmax,"");
  h_Aq22_SLR_centN = new TProfile("h_Aq22_SLR_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_Aq22_SLR_cent);
  fOutputList->Add(h_Aq22_SLR_centP);
  fOutputList->Add(h_Aq22_SLR_centN);

  // --- third harmonic
  h_Aq32_cent = new TProfile("h_Aq32_cent","",100,0,100,tpmin,tpmax,"");
  h_Aq32_centP = new TProfile("h_Aq32_centP","",100,0,100,tpmin,tpmax,"");
  h_Aq32_centN = new TProfile("h_Aq32_centN","",100,0,100,tpmin,tpmax,"");
  h_Aq32_centPP = new TProfile("h_Aq32_centPP","",100,0,100,tpmin,tpmax,"");
  h_Aq32_centNN = new TProfile("h_Aq32_centNN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_Aq32_cent);
  fOutputList->Add(h_Aq32_centP);
  fOutputList->Add(h_Aq32_centN);
  fOutputList->Add(h_Aq32_centPP);
  fOutputList->Add(h_Aq32_centNN);
  h_Aq32_SL_cent = new TProfile("h_Aq32_SL_cent","",100,0,100,tpmin,tpmax,"");
  h_Aq32_SL_centP = new TProfile("h_Aq32_SL_centP","",100,0,100,tpmin,tpmax,"");
  h_Aq32_SL_centN = new TProfile("h_Aq32_SL_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_Aq32_SL_cent);
  fOutputList->Add(h_Aq32_SL_centP);
  fOutputList->Add(h_Aq32_SL_centN);
  h_Aq32_SLL_cent = new TProfile("h_Aq32_SLL_cent","",100,0,100,tpmin,tpmax,"");
  h_Aq32_SLL_centP = new TProfile("h_Aq32_SLL_centP","",100,0,100,tpmin,tpmax,"");
  h_Aq32_SLL_centN = new TProfile("h_Aq32_SLL_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_Aq32_SLL_cent);
  fOutputList->Add(h_Aq32_SLL_centP);
  fOutputList->Add(h_Aq32_SLL_centN);
  h_Aq32_SLR_cent = new TProfile("h_Aq32_SLR_cent","",100,0,100,tpmin,tpmax,"");
  h_Aq32_SLR_centP = new TProfile("h_Aq32_SLR_centP","",100,0,100,tpmin,tpmax,"");
  h_Aq32_SLR_centN = new TProfile("h_Aq32_SLR_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_Aq32_SLR_cent);
  fOutputList->Add(h_Aq32_SLR_centP);
  fOutputList->Add(h_Aq32_SLR_centN);

  // --- fourth harmonic
  h_Aq42_cent = new TProfile("h_Aq42_cent","",100,0,100,tpmin,tpmax,"");
  h_Aq42_centP = new TProfile("h_Aq42_centP","",100,0,100,tpmin,tpmax,"");
  h_Aq42_centN = new TProfile("h_Aq42_centN","",100,0,100,tpmin,tpmax,"");
  h_Aq42_centPP = new TProfile("h_Aq42_centPP","",100,0,100,tpmin,tpmax,"");
  h_Aq42_centNN = new TProfile("h_Aq42_centNN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_Aq42_cent);
  fOutputList->Add(h_Aq42_centP);
  fOutputList->Add(h_Aq42_centN);
  fOutputList->Add(h_Aq42_centPP);
  fOutputList->Add(h_Aq42_centNN);
  h_Aq42_SL_cent = new TProfile("h_Aq42_SL_cent","",100,0,100,tpmin,tpmax,"");
  h_Aq42_SL_centP = new TProfile("h_Aq42_SL_centP","",100,0,100,tpmin,tpmax,"");
  h_Aq42_SL_centN = new TProfile("h_Aq42_SL_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_Aq42_SL_cent);
  fOutputList->Add(h_Aq42_SL_centP);
  fOutputList->Add(h_Aq42_SL_centN);
  h_Aq42_SLL_cent = new TProfile("h_Aq42_SLL_cent","",100,0,100,tpmin,tpmax,"");
  h_Aq42_SLL_centP = new TProfile("h_Aq42_SLL_centP","",100,0,100,tpmin,tpmax,"");
  h_Aq42_SLL_centN = new TProfile("h_Aq42_SLL_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_Aq42_SLL_cent);
  fOutputList->Add(h_Aq42_SLL_centP);
  fOutputList->Add(h_Aq42_SLL_centN);
  h_Aq42_SLR_cent = new TProfile("h_Aq42_SLR_cent","",100,0,100,tpmin,tpmax,"");
  h_Aq42_SLR_centP = new TProfile("h_Aq42_SLR_centP","",100,0,100,tpmin,tpmax,"");
  h_Aq42_SLR_centN = new TProfile("h_Aq42_SLR_centN","",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_Aq42_SLR_cent);
  fOutputList->Add(h_Aq42_SLR_centP);
  fOutputList->Add(h_Aq42_SLR_centN);




  // ------------------------------------
  // --- differential cumulant histograms
  // ------------------------------------

  h_diffq22_pt = new TProfile("h_diffq22_pt","",50,0,5,tpmin,tpmax,"");
  h_diffq22_ptP = new TProfile("h_diffq22_ptP","",50,0,5,tpmin,tpmax,"");
  h_diffq22_ptN = new TProfile("h_diffq22_ptN","",50,0,5,tpmin,tpmax,"");
  h_diffAq22_pt = new TProfile("h_diffAq22_pt","",50,0,5,tpmin,tpmax,"");
  h_diffAq22_ptP = new TProfile("h_diffAq22_ptP","",50,0,5,tpmin,tpmax,"");
  h_diffAq22_ptN = new TProfile("h_diffAq22_ptN","",50,0,5,tpmin,tpmax,"");
  fOutputList->Add(h_diffq22_pt);
  fOutputList->Add(h_diffq22_ptP);
  fOutputList->Add(h_diffq22_ptN);
  fOutputList->Add(h_diffAq22_pt);
  fOutputList->Add(h_diffAq22_ptP);
  fOutputList->Add(h_diffAq22_ptN);
  h_diffq22_eta = new TProfile("h_diffq22_eta","",32,-0.8,0.8,tpmin,tpmax,"");
  h_diffq22_etaP = new TProfile("h_diffq22_etaP","",32,-0.8,0.8,tpmin,tpmax,"");
  h_diffq22_etaN = new TProfile("h_diffq22_etaN","",32,-0.8,0.8,tpmin,tpmax,"");
  h_diffAq22_eta = new TProfile("h_diffAq22_eta","",32,-0.8,0.8,tpmin,tpmax,"");
  h_diffAq22_etaP = new TProfile("h_diffAq22_etaP","",32,-0.8,0.8,tpmin,tpmax,"");
  h_diffAq22_etaN = new TProfile("h_diffAq22_etaN","",32,-0.8,0.8,tpmin,tpmax,"");
  fOutputList->Add(h_diffq22_eta);
  fOutputList->Add(h_diffq22_etaP);
  fOutputList->Add(h_diffq22_etaN);
  fOutputList->Add(h_diffAq22_eta);
  fOutputList->Add(h_diffAq22_etaP);
  fOutputList->Add(h_diffAq22_etaN);

  h_diffq32_pt = new TProfile("h_diffq32_pt","",50,0,5,tpmin,tpmax,"");
  h_diffq32_ptP = new TProfile("h_diffq32_ptP","",50,0,5,tpmin,tpmax,"");
  h_diffq32_ptN = new TProfile("h_diffq32_ptN","",50,0,5,tpmin,tpmax,"");
  h_diffAq32_pt = new TProfile("h_diffAq32_pt","",50,0,5,tpmin,tpmax,"");
  h_diffAq32_ptP = new TProfile("h_diffAq32_ptP","",50,0,5,tpmin,tpmax,"");
  h_diffAq32_ptN = new TProfile("h_diffAq32_ptN","",50,0,5,tpmin,tpmax,"");
  fOutputList->Add(h_diffq32_pt);
  fOutputList->Add(h_diffq32_ptP);
  fOutputList->Add(h_diffq32_ptN);
  fOutputList->Add(h_diffAq32_pt);
  fOutputList->Add(h_diffAq32_ptP);
  fOutputList->Add(h_diffAq32_ptN);
  h_diffq32_eta = new TProfile("h_diffq32_eta","",32,-0.8,0.8,tpmin,tpmax,"");
  h_diffq32_etaP = new TProfile("h_diffq32_etaP","",32,-0.8,0.8,tpmin,tpmax,"");
  h_diffq32_etaN = new TProfile("h_diffq32_etaN","",32,-0.8,0.8,tpmin,tpmax,"");
  h_diffAq32_eta = new TProfile("h_diffAq32_eta","",32,-0.8,0.8,tpmin,tpmax,"");
  h_diffAq32_etaP = new TProfile("h_diffAq32_etaP","",32,-0.8,0.8,tpmin,tpmax,"");
  h_diffAq32_etaN = new TProfile("h_diffAq32_etaN","",32,-0.8,0.8,tpmin,tpmax,"");
  fOutputList->Add(h_diffq32_eta);
  fOutputList->Add(h_diffq32_etaP);
  fOutputList->Add(h_diffq32_etaN);
  fOutputList->Add(h_diffAq32_eta);
  fOutputList->Add(h_diffAq32_etaP);
  fOutputList->Add(h_diffAq32_etaN);

  h_diffq42_pt = new TProfile("h_diffq42_pt","",50,0,5,tpmin,tpmax,"");
  h_diffq42_ptP = new TProfile("h_diffq42_ptP","",50,0,5,tpmin,tpmax,"");
  h_diffq42_ptN = new TProfile("h_diffq42_ptN","",50,0,5,tpmin,tpmax,"");
  h_diffAq42_pt = new TProfile("h_diffAq42_pt","",50,0,5,tpmin,tpmax,"");
  h_diffAq42_ptP = new TProfile("h_diffAq42_ptP","",50,0,5,tpmin,tpmax,"");
  h_diffAq42_ptN = new TProfile("h_diffAq42_ptN","",50,0,5,tpmin,tpmax,"");
  fOutputList->Add(h_diffq42_pt);
  fOutputList->Add(h_diffq42_ptP);
  fOutputList->Add(h_diffq42_ptN);
  fOutputList->Add(h_diffAq42_pt);
  fOutputList->Add(h_diffAq42_ptP);
  fOutputList->Add(h_diffAq42_ptN);
  h_diffq42_eta = new TProfile("h_diffq42_eta","",32,-0.8,0.8,tpmin,tpmax,"");
  h_diffq42_etaP = new TProfile("h_diffq42_etaP","",32,-0.8,0.8,tpmin,tpmax,"");
  h_diffq42_etaN = new TProfile("h_diffq42_etaN","",32,-0.8,0.8,tpmin,tpmax,"");
  h_diffAq42_eta = new TProfile("h_diffAq42_eta","",32,-0.8,0.8,tpmin,tpmax,"");
  h_diffAq42_etaP = new TProfile("h_diffAq42_etaP","",32,-0.8,0.8,tpmin,tpmax,"");
  h_diffAq42_etaN = new TProfile("h_diffAq42_etaN","",32,-0.8,0.8,tpmin,tpmax,"");
  fOutputList->Add(h_diffq42_eta);
  fOutputList->Add(h_diffq42_etaP);
  fOutputList->Add(h_diffq42_etaN);
  fOutputList->Add(h_diffAq42_eta);
  fOutputList->Add(h_diffAq42_etaP);
  fOutputList->Add(h_diffAq42_etaN);



  h_AT_X_deta = new TProfile("h_AT_X_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  h_AT_X_detaP = new TProfile("h_AT_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_AT_X_detaN = new TProfile("h_AT_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_AT_X_deta);
  fOutputList->Add(h_AT_X_detaP);
  fOutputList->Add(h_AT_X_detaN);
  //h_diffq22_X_deta = new TProfile("h_diffq22_X_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq22_X_detaP = new TProfile("h_diffq22_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq22_X_detaN = new TProfile("h_diffq22_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  //h_ATdiffq22_X_deta = new TProfile("h_ATdiffq22_X_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq22_X_detaP = new TProfile("h_ATdiffq22_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq22_X_detaN = new TProfile("h_ATdiffq22_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_diffq22_X_detaP);
  fOutputList->Add(h_diffq22_X_detaN);
  fOutputList->Add(h_ATdiffq22_X_detaP);
  fOutputList->Add(h_ATdiffq22_X_detaN);
  //h_diffq32_X_deta = new TProfile("h_diffq32_X_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq32_X_detaP = new TProfile("h_diffq32_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq32_X_detaN = new TProfile("h_diffq32_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  //h_ATdiffq32_X_deta = new TProfile("h_ATdiffq32_X_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq32_X_detaP = new TProfile("h_ATdiffq32_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq32_X_detaN = new TProfile("h_ATdiffq32_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_diffq32_X_detaP);
  fOutputList->Add(h_diffq32_X_detaN);
  fOutputList->Add(h_ATdiffq32_X_detaP);
  fOutputList->Add(h_ATdiffq32_X_detaN);
  //h_diffq42_X_deta = new TProfile("h_diffq42_X_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq42_X_detaP = new TProfile("h_diffq42_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq42_X_detaN = new TProfile("h_diffq42_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  //h_ATdiffq42_X_deta = new TProfile("h_ATdiffq42_X_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq42_X_detaP = new TProfile("h_ATdiffq42_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq42_X_detaN = new TProfile("h_ATdiffq42_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_diffq42_X_detaP);
  fOutputList->Add(h_diffq42_X_detaN);
  fOutputList->Add(h_ATdiffq42_X_detaP);
  fOutputList->Add(h_ATdiffq42_X_detaN);


  // --- now field selection
  // --- field 1
  h_AT_X_deta_F1 = new TProfile("h_AT_X_deta_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  h_AT_X_detaP_F1 = new TProfile("h_AT_X_detaP_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  h_AT_X_detaN_F1 = new TProfile("h_AT_X_detaN_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_AT_X_deta_F1);
  fOutputList->Add(h_AT_X_detaP_F1);
  fOutputList->Add(h_AT_X_detaN_F1);
  //h_diffq22_X_deta_F1 = new TProfile("h_diffq22_X_deta_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq22_X_detaP_F1 = new TProfile("h_diffq22_X_detaP_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq22_X_detaN_F1 = new TProfile("h_diffq22_X_detaN_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  //h_ATdiffq22_X_deta_F1 = new TProfile("h_ATdiffq22_X_deta_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq22_X_detaP_F1 = new TProfile("h_ATdiffq22_X_detaP_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq22_X_detaN_F1 = new TProfile("h_ATdiffq22_X_detaN_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_diffq22_X_detaP_F1);
  fOutputList->Add(h_diffq22_X_detaN_F1);
  fOutputList->Add(h_ATdiffq22_X_detaP_F1);
  fOutputList->Add(h_ATdiffq22_X_detaN_F1);
  //h_diffq32_X_deta_F1 = new TProfile("h_diffq32_X_deta_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq32_X_detaP_F1 = new TProfile("h_diffq32_X_detaP_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq32_X_detaN_F1 = new TProfile("h_diffq32_X_detaN_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  //h_ATdiffq32_X_deta_F1 = new TProfile("h_ATdiffq32_X_deta_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq32_X_detaP_F1 = new TProfile("h_ATdiffq32_X_detaP_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq32_X_detaN_F1 = new TProfile("h_ATdiffq32_X_detaN_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_diffq32_X_detaP_F1);
  fOutputList->Add(h_diffq32_X_detaN_F1);
  fOutputList->Add(h_ATdiffq32_X_detaP_F1);
  fOutputList->Add(h_ATdiffq32_X_detaN_F1);
  //h_diffq42_X_deta_F1 = new TProfile("h_diffq42_X_deta_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq42_X_detaP_F1 = new TProfile("h_diffq42_X_detaP_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq42_X_detaN_F1 = new TProfile("h_diffq42_X_detaN_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  //h_ATdiffq42_X_deta_F1 = new TProfile("h_ATdiffq42_X_deta_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq42_X_detaP_F1 = new TProfile("h_ATdiffq42_X_detaP_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq42_X_detaN_F1 = new TProfile("h_ATdiffq42_X_detaN_F1","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_diffq42_X_detaP_F1);
  fOutputList->Add(h_diffq42_X_detaN_F1);
  fOutputList->Add(h_ATdiffq42_X_detaP_F1);
  fOutputList->Add(h_ATdiffq42_X_detaN_F1);
  // --- field 3
  h_AT_X_deta_F3 = new TProfile("h_AT_X_deta_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  h_AT_X_detaP_F3 = new TProfile("h_AT_X_detaP_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  h_AT_X_detaN_F3 = new TProfile("h_AT_X_detaN_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_AT_X_deta_F3);
  fOutputList->Add(h_AT_X_detaP_F3);
  fOutputList->Add(h_AT_X_detaN_F3);
  //h_diffq22_X_deta_F3 = new TProfile("h_diffq22_X_deta_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq22_X_detaP_F3 = new TProfile("h_diffq22_X_detaP_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq22_X_detaN_F3 = new TProfile("h_diffq22_X_detaN_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  //h_ATdiffq22_X_deta_F3 = new TProfile("h_ATdiffq22_X_deta_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq22_X_detaP_F3 = new TProfile("h_ATdiffq22_X_detaP_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq22_X_detaN_F3 = new TProfile("h_ATdiffq22_X_detaN_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_diffq22_X_detaP_F3);
  fOutputList->Add(h_diffq22_X_detaN_F3);
  fOutputList->Add(h_ATdiffq22_X_detaP_F3);
  fOutputList->Add(h_ATdiffq22_X_detaN_F3);
  //h_diffq32_X_deta_F3 = new TProfile("h_diffq32_X_deta_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq32_X_detaP_F3 = new TProfile("h_diffq32_X_detaP_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq32_X_detaN_F3 = new TProfile("h_diffq32_X_detaN_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  //h_ATdiffq32_X_deta_F3 = new TProfile("h_ATdiffq32_X_deta_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq32_X_detaP_F3 = new TProfile("h_ATdiffq32_X_detaP_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq32_X_detaN_F3 = new TProfile("h_ATdiffq32_X_detaN_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_diffq32_X_detaP_F3);
  fOutputList->Add(h_diffq32_X_detaN_F3);
  fOutputList->Add(h_ATdiffq32_X_detaP_F3);
  fOutputList->Add(h_ATdiffq32_X_detaN_F3);
  //h_diffq42_X_deta_F3 = new TProfile("h_diffq42_X_deta_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq42_X_detaP_F3 = new TProfile("h_diffq42_X_detaP_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq42_X_detaN_F3 = new TProfile("h_diffq42_X_detaN_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  //h_ATdiffq42_X_deta_F3 = new TProfile("h_ATdiffq42_X_deta_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq42_X_detaP_F3 = new TProfile("h_ATdiffq42_X_detaP_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq42_X_detaN_F3 = new TProfile("h_ATdiffq42_X_detaN_F3","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_diffq42_X_detaP_F3);
  fOutputList->Add(h_diffq42_X_detaN_F3);
  fOutputList->Add(h_ATdiffq42_X_detaP_F3);
  fOutputList->Add(h_ATdiffq42_X_detaN_F3);





  h_ATXdiffq22_X_detaP = new TProfile("h_ATXdiffq22_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATXdiffq22_X_detaN = new TProfile("h_ATXdiffq22_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_ATXdiffq22_X_detaP);
  fOutputList->Add(h_ATXdiffq22_X_detaN);
  h_ATYdiffq22_X_detaP = new TProfile("h_ATYdiffq22_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATYdiffq22_X_detaN = new TProfile("h_ATYdiffq22_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_ATYdiffq22_X_detaP);
  fOutputList->Add(h_ATYdiffq22_X_detaN);
  h_ATZdiffq22_X_detaP = new TProfile("h_ATZdiffq22_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATZdiffq22_X_detaN = new TProfile("h_ATZdiffq22_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_ATZdiffq22_X_detaP);
  fOutputList->Add(h_ATZdiffq22_X_detaN);

  h_ATXdiffq32_X_detaP = new TProfile("h_ATXdiffq32_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATXdiffq32_X_detaN = new TProfile("h_ATXdiffq32_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_ATXdiffq32_X_detaP);
  fOutputList->Add(h_ATXdiffq32_X_detaN);
  h_ATYdiffq32_X_detaP = new TProfile("h_ATYdiffq32_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATYdiffq32_X_detaN = new TProfile("h_ATYdiffq32_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_ATYdiffq32_X_detaP);
  fOutputList->Add(h_ATYdiffq32_X_detaN);
  h_ATZdiffq32_X_detaP = new TProfile("h_ATZdiffq32_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATZdiffq32_X_detaN = new TProfile("h_ATZdiffq32_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_ATZdiffq32_X_detaP);
  fOutputList->Add(h_ATZdiffq32_X_detaN);

  h_ATXdiffq42_X_detaP = new TProfile("h_ATXdiffq42_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATXdiffq42_X_detaN = new TProfile("h_ATXdiffq42_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_ATXdiffq42_X_detaP);
  fOutputList->Add(h_ATXdiffq42_X_detaN);
  h_ATYdiffq42_X_detaP = new TProfile("h_ATYdiffq42_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATYdiffq42_X_detaN = new TProfile("h_ATYdiffq42_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_ATYdiffq42_X_detaP);
  fOutputList->Add(h_ATYdiffq42_X_detaN);
  h_ATZdiffq42_X_detaP = new TProfile("h_ATZdiffq42_X_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATZdiffq42_X_detaN = new TProfile("h_ATZdiffq42_X_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_ATZdiffq42_X_detaP);
  fOutputList->Add(h_ATZdiffq42_X_detaN);



  h_AT_X_dpt = new TProfile("h_AT_X_dpt","",32,-1.6,1.6,tpmin,tpmax,"");
  h_AT_X_dptP = new TProfile("h_AT_X_dptP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_AT_X_dptN = new TProfile("h_AT_X_dptN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_AT_X_dpt);
  fOutputList->Add(h_AT_X_dptP);
  fOutputList->Add(h_AT_X_dptN);
  //h_diffq22_X_dpt = new TProfile("h_diffq22_X_dpt","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq22_X_dptP = new TProfile("h_diffq22_X_dptP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq22_X_dptN = new TProfile("h_diffq22_X_dptN","",32,-1.6,1.6,tpmin,tpmax,"");
  //h_ATdiffq22_X_dpt = new TProfile("h_ATdiffq22_X_dpt","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq22_X_dptP = new TProfile("h_ATdiffq22_X_dptP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq22_X_dptN = new TProfile("h_ATdiffq22_X_dptN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_diffq22_X_dptP);
  fOutputList->Add(h_diffq22_X_dptN);
  fOutputList->Add(h_ATdiffq22_X_dptP);
  fOutputList->Add(h_ATdiffq22_X_dptN);
  //h_diffq32_X_dpt = new TProfile("h_diffq32_X_dpt","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq32_X_dptP = new TProfile("h_diffq32_X_dptP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq32_X_dptN = new TProfile("h_diffq32_X_dptN","",32,-1.6,1.6,tpmin,tpmax,"");
  //h_ATdiffq32_X_dpt = new TProfile("h_ATdiffq32_X_dpt","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq32_X_dptP = new TProfile("h_ATdiffq32_X_dptP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq32_X_dptN = new TProfile("h_ATdiffq32_X_dptN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_diffq32_X_dptP);
  fOutputList->Add(h_diffq32_X_dptN);
  fOutputList->Add(h_ATdiffq32_X_dptP);
  fOutputList->Add(h_ATdiffq32_X_dptN);
  //h_diffq42_X_dpt = new TProfile("h_diffq42_X_dpt","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq42_X_dptP = new TProfile("h_diffq42_X_dptP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq42_X_dptN = new TProfile("h_diffq42_X_dptN","",32,-1.6,1.6,tpmin,tpmax,"");
  //h_ATdiffq42_X_dpt = new TProfile("h_ATdiffq42_X_dpt","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq42_X_dptP = new TProfile("h_ATdiffq42_X_dptP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq42_X_dptN = new TProfile("h_ATdiffq42_X_dptN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_diffq42_X_dptP);
  fOutputList->Add(h_diffq42_X_dptN);
  fOutputList->Add(h_ATdiffq42_X_dptP);
  fOutputList->Add(h_ATdiffq42_X_dptN);



  // --- now Xsub...

  h_AT_S_deta = new TProfile("h_AT_S_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  h_AT_S_detaP = new TProfile("h_AT_S_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_AT_S_detaN = new TProfile("h_AT_S_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_AT_S_deta);
  fOutputList->Add(h_AT_S_detaP);
  fOutputList->Add(h_AT_S_detaN);
  //h_diffq22_S_deta = new TProfile("h_diffq22_S_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq22_S_detaP = new TProfile("h_diffq22_S_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq22_S_detaN = new TProfile("h_diffq22_S_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  //h_ATdiffq22_S_deta = new TProfile("h_ATdiffq22_S_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq22_S_detaP = new TProfile("h_ATdiffq22_S_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq22_S_detaN = new TProfile("h_ATdiffq22_S_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_diffq22_S_detaP);
  fOutputList->Add(h_diffq22_S_detaN);
  fOutputList->Add(h_ATdiffq22_S_detaP);
  fOutputList->Add(h_ATdiffq22_S_detaN);
  //h_diffq32_S_deta = new TProfile("h_diffq32_S_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq32_S_detaP = new TProfile("h_diffq32_S_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq32_S_detaN = new TProfile("h_diffq32_S_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  //h_ATdiffq32_S_deta = new TProfile("h_ATdiffq32_S_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq32_S_detaP = new TProfile("h_ATdiffq32_S_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq32_S_detaN = new TProfile("h_ATdiffq32_S_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_diffq32_S_detaP);
  fOutputList->Add(h_diffq32_S_detaN);
  fOutputList->Add(h_ATdiffq32_S_detaP);
  fOutputList->Add(h_ATdiffq32_S_detaN);
  //h_diffq42_S_deta = new TProfile("h_diffq42_S_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq42_S_detaP = new TProfile("h_diffq42_S_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_diffq42_S_detaN = new TProfile("h_diffq42_S_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  //h_ATdiffq42_S_deta = new TProfile("h_ATdiffq42_S_deta","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq42_S_detaP = new TProfile("h_ATdiffq42_S_detaP","",32,-1.6,1.6,tpmin,tpmax,"");
  h_ATdiffq42_S_detaN = new TProfile("h_ATdiffq42_S_detaN","",32,-1.6,1.6,tpmin,tpmax,"");
  fOutputList->Add(h_diffq42_S_detaP);
  fOutputList->Add(h_diffq42_S_detaN);
  fOutputList->Add(h_ATdiffq42_S_detaP);
  fOutputList->Add(h_ATdiffq42_S_detaN);

  // ---

  // --- random fun

  // ---

  h_a1q22_cent = new TProfile("h_a1q22_cent","q22 vs centrality",100,0,100,tpmin,tpmax,"");
  h_a1q22_centP = new TProfile("h_a1q22_centP","q22(pos) vs centrality",100,0,100,tpmin,tpmax,"");
  h_a1q22_centN = new TProfile("h_a1q22_centN","q22(neg) vs centrality",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_a1q22_cent);
  fOutputList->Add(h_a1q22_centP);
  fOutputList->Add(h_a1q22_centN);

  h_a2q22_cent = new TProfile("h_a2q22_cent","q22 vs centrality",100,0,100,tpmin,tpmax,"");
  h_a2q22_centP = new TProfile("h_a2q22_centP","q22(pos) vs centrality",100,0,100,tpmin,tpmax,"");
  h_a2q22_centN = new TProfile("h_a2q22_centN","q22(neg) vs centrality",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_a2q22_cent);
  fOutputList->Add(h_a2q22_centP);
  fOutputList->Add(h_a2q22_centN);

  h_rAq22_X1_cent = new TProfile("h_rAq22_X1_cent","q22 vs centrality",100,0,100,tpmin,tpmax,"");
  h_rAq22_X1_centP = new TProfile("h_rAq22_X1_centP","q22(pos) vs centrality",100,0,100,tpmin,tpmax,"");
  h_rAq22_X1_centN = new TProfile("h_rAq22_X1_centN","q22(neg) vs centrality",100,0,100,tpmin,tpmax,"");
  h_rAq22_X2_cent = new TProfile("h_rAq22_X2_cent","q22 vs centrality",100,0,100,tpmin,tpmax,"");
  h_rAq22_X2_centP = new TProfile("h_rAq22_X2_centP","q22(pos) vs centrality",100,0,100,tpmin,tpmax,"");
  h_rAq22_X2_centN = new TProfile("h_rAq22_X2_centN","q22(neg) vs centrality",100,0,100,tpmin,tpmax,"");
  h_rAq22_X3_cent = new TProfile("h_rAq22_X3_cent","q22 vs centrality",100,0,100,tpmin,tpmax,"");
  h_rAq22_X3_centP = new TProfile("h_rAq22_X3_centP","q22(pos) vs centrality",100,0,100,tpmin,tpmax,"");
  h_rAq22_X3_centN = new TProfile("h_rAq22_X3_centN","q22(neg) vs centrality",100,0,100,tpmin,tpmax,"");
  h_rAq22_X4_cent = new TProfile("h_rAq22_X4_cent","q22 vs centrality",100,0,100,tpmin,tpmax,"");
  h_rAq22_X4_centP = new TProfile("h_rAq22_X4_centP","q22(pos) vs centrality",100,0,100,tpmin,tpmax,"");
  h_rAq22_X4_centN = new TProfile("h_rAq22_X4_centN","q22(neg) vs centrality",100,0,100,tpmin,tpmax,"");
  fOutputList->Add(h_rAq22_X1_cent);
  fOutputList->Add(h_rAq22_X1_centP);
  fOutputList->Add(h_rAq22_X1_centN);
  fOutputList->Add(h_rAq22_X2_cent);
  fOutputList->Add(h_rAq22_X2_centP);
  fOutputList->Add(h_rAq22_X2_centN);
  fOutputList->Add(h_rAq22_X3_cent);
  fOutputList->Add(h_rAq22_X3_centP);
  fOutputList->Add(h_rAq22_X3_centN);
  fOutputList->Add(h_rAq22_X4_cent);
  fOutputList->Add(h_rAq22_X4_centP);
  fOutputList->Add(h_rAq22_X4_centN);



  // ---------------------------- //
  // --- done with histograms --- //
  // ---------------------------- //


  // ------------------------------------ //
  // --- send data to the output list --- //
  // ------------------------------------ // 
  PostData(1,fOutputList);
}



//----------------------------------------------
void AliAnalysisTaskCMEv2A::UserExec(Option_t *) 
{

  // Main analyis loop, called for each event

  if(debug>0) cout<<"Processing event with debug = "<<debug<<endl;

  // Analysis manager
  AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
  if(!manager)
    {
      cout<<"FATAL: could not get Analysis Manager."<<endl;
      return;
    }

  // Input handler
  AliAODInputHandler *handler = (AliAODInputHandler *)manager->GetInputEventHandler();
  if(!handler)
    {
      cout<<"FATAL: could not get Input Handler."<<endl;
      return;
    }


  // AOD Event object from tree
  AliAODEvent *fAOD = (AliAODEvent *)InputEvent();
  if(!fAOD)
    {
      if(debug>-1) cout<<"ERROR: AOD event object not available. Discarding event..."<<endl;
      return;
    }


  int runnumber = fAOD->GetRunNumber();
  float mag = fAOD->GetMagneticField();

  if(debug>0)
    {
      cout<<"runnumber is "<<runnumber<<endl;
      cout<<"magnetic field is "<<mag<<endl;
    }

  // MC Event object from tree
  AliMCEvent *fMC = MCEvent();
  if(!fMC&&doMC)
    {
      if(debug>-1) cout<<"ERROR: MC event object not available. Discarding event..."<<endl;
      return;
    }

  // get pid
  AliPIDResponse *fPID = handler->GetPIDResponse();
  if(!fPID && !doMC)    // use PID only if no MC 
    {
      if(debug>-1) cout<<"ERROR: PIDResponse object not available. Discarding event..."<<endl;
      return;
    }

  // // Eventplane object (from AOD)
  // AliEventplane *fEventplane = fAOD->GetEventplane();
  // if(!fEventplane)
  //   {
  //     if(debug>-1) cout<<"ERROR: Eventplane object not available. Discarding event..."<<endl;
  //     return;
  //   }

  // float psi_V0_h2 = fEventplane->GetEventplane("V0",fAOD,2);
  // float psi_V0A_h2 = fEventplane->GetEventplane("V0A",fAOD,2);
  // float psi_V0C_h2 = fEventplane->GetEventplane("V0C",fAOD,2);

  // fHistPlaneV0h2->Fill(psi_V0_h2);
  // fHistPlaneV0Ah2->Fill(psi_V0A_h2);
  // fHistPlaneV0Ch2->Fill(psi_V0C_h2);
  // fHistPlaneV0ACDCh2->Fill(psi_V0A_h2-psi_V0C_h2);



  // Centrality object (from AOD)
  AliCentrality *fCentrality = fAOD->GetCentrality();
  if(!fCentrality)
    {
      if(debug>-1) cout<<"ERROR: Centrality object not available. Discarding event..."<<endl;
      return;
    }

  float centTRK = fCentrality->GetCentralityPercentile("TRK");
  float centV0M = fCentrality->GetCentralityPercentile("V0M");
  float centSPD = fCentrality->GetCentralityPercentile("CL1");//outer SPD?
  float cent = centTRK;
  if(centhandle==2) cent = centV0M;
  if(centhandle==3) cent = centSPD;

  int icent = int(cent)/10;

  int centstatus = 0;
  if(centTRK<0.0||centV0M<0.0) centstatus = -1;
  if(centTRK<0.0&&centV0M<0.0) centstatus = -2;
  if(centTRK>0.0||centV0M>0.0) centstatus = 1;
  if(centTRK>0.0&&centV0M>0.0) centstatus = 2;
  if(centTRK==0.0&&centV0M==0.0) centstatus = 0;

  if(debug>0)
    {
      cout<<"centTRK "<<centTRK<<endl;
      cout<<"centV0M "<<centV0M<<endl;
      cout<<"centSPD "<<centSPD<<endl;
      cout<<"centrality selection is "<<centhandle<<endl;
      cout<<"cent is "<<cent<<endl;
    }
  fHistCentTRK->Fill(centTRK);
  fHistCentV0M->Fill(centV0M);
  fHistCentDIFF->Fill(centTRK-centV0M);

  //ULong64_t mask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  ULong64_t mask = handler->IsEventSelected();
  ULong64_t amb = AliVEvent::kMB;
  ULong64_t acn = AliVEvent::kCentral;
  ULong64_t asc = AliVEvent::kSemiCentral;

  if(debug>0) cout<<"trigger selection is "<<trigger<<endl;
  if(debug>0) cout<<"trigger mask is "<<mask<<endl;

  if(mask&amb)
    {
      fHistCentTRKAVEkMB->Fill(centTRK);
      fHistCentV0MAVEkMB->Fill(centV0M);
    }
  if(mask&(acn&(amb|asc)))
    {
      fHistCentTRKAVEkCentral->Fill(centTRK);
      fHistCentV0MAVEkCentral->Fill(centV0M);
    }
  if(mask&(asc|amb))
    {
      fHistCentTRKAVEkSemiCentral->Fill(centTRK);
      fHistCentV0MAVEkSemiCentral->Fill(centV0M);
    }
  if(mask&(amb|acn|asc))
    {
      fHistCentTRKAVEkA3->Fill(centTRK);
      fHistCentV0MAVEkA3->Fill(centV0M);
    }


  // Ionut
  Bool_t isMB = (fAOD->GetTriggerMask() & (ULong64_t(1)<<1));
  Bool_t isCentral = (fAOD->GetTriggerMask() & (ULong64_t(1)<<4));
  Bool_t isSemiCentral = (fAOD->GetTriggerMask() & (ULong64_t(1)<<7));

  cout<<"other trigger mask is "<<fAOD->GetTriggerMask()<<endl;

  //if(mask&trigger)
  //if(isSemiCentral||isMB)
  if(isCentral||isSemiCentral||isMB)
    {
      fHistCentTRKAVEkSel->Fill(centTRK);
      fHistCentV0MAVEkSel->Fill(centV0M);
    }
  else
    {
      if(debug>0) cout<<"wrong trigger, rejecting event"<<endl;
      return;
    }

  if(fabs(centTRK-centV0M)>centcut)
    {
      if(debug>0) cout<<"centrality difference outside cut, rejecting event"<<endl;
      return;
    }



  // AOD vertex objects
  AliAODVertex *fVtx = fAOD->GetPrimaryVertex();
  AliAODVertex *fVtxSPD = fAOD->GetPrimaryVertexSPD();
  if(!fVtx)
    {
      if(debug>-1) cout<<"ERROR: Vertex object not available. Discarding event..."<<endl;
      return;
    }
  float zvtxV0 = fVtx->GetZ();
  float zvtxSPD = fVtxSPD->GetZ();
  if(debug>0) cout<<"zvtxV0 is "<<zvtxV0<<endl;
  if(debug>0) cout<<"zvtxSPD is "<<zvtxSPD<<endl;
  if(centstatus==2)
    {
      fHistZVtx->Fill(zvtxV0);
      fHistZVtxD->Fill(zvtxV0-zvtxSPD);
    }
  if(fabs(zvtxV0)>zvtxcut)
    {
      if(debug>0) cout<<"vertex outside cut, rejecting event"<<endl;
      return;
    }
  float eventX = fVtx->GetX();
  float eventY = fVtx->GetY();
  float eventZ = fVtx->GetZ();



  // // get V0 object
  // AliAODVZERO *fV0 = fAOD->GetVZEROData();
  // for(int i=0; i<64; i++)
  //   {
  //     float phiV0 = pi/4.0*(0.5+i%8);
  //     float multV0 = fV0->GetMultiplicity(i);
  //   }


  int d_ntrk = fAOD->GetNumberOfTracks();
  int d_ntrkMC = 0;
  if(fMC) d_ntrkMC = fMC->GetNumberOfTracks();

  if(centstatus==-2&&d_ntrk>0) centstatus = -3;
  if(debug>0) cout<<"there are "<<d_ntrk<<" tracks in this event"<<endl;
  if(debug>0&&fMC) cout<<"there are "<<d_ntrkMC<<" Monte Carlo tracks in this event"<<endl;
  if(debug>0) cout<<"centrality diagnostic is "<<centstatus<<endl;

  fHistCentDIAG->Fill(centstatus);
  if(centstatus==-3)
    {
      fHistCtrkDIAG->Fill(d_ntrk);
      fHistVtxRDIAG->Fill(zvtxV0);
      fHistVtxSDIAG->Fill(zvtxSPD);
    }
  if(centstatus==-2)
    {
      fHistVtxRDIBG->Fill(zvtxV0);
      fHistVtxSDIBG->Fill(zvtxSPD);
    }

  // --- AliAnalysisUtils for pileup cuts
  if(dopupcut)
    {
      AliAnalysisUtils *fUtils = new AliAnalysisUtils();
      if(!fUtils)
  	{
  	  if(debug>-1) cout<<"ERROR: cannot find AliAnalysisUtils..."<<endl;
  	  return;
  	}
      fUtils->SetUseOutOfBunchPileUp(true);
      bool pileup = fUtils->IsPileUpEvent(fAOD);
      //bool pileup = fAOD->IsPileUpFromSPD(3,0.8,3.0,2.0,5.0); // default parameters in AliAODEvent.h
      if(pileup)
  	{
  	  if(debug>0) cout<<"Rejecting event for pileup (AliAnalysisUtils)"<<endl;
  	  return;
  	}
      if(fUtils) delete fUtils;
      //
      // ---
      //
      /*
      const AliAODVertex* vtxSPD = fAOD->GetPrimaryVertexSPD();
      if(!vtxSPD||vtxSPD->GetNContributors()<=0)
	{
	  if(debug>0) cout<<"rejecting pileup event (no vtxSPD or zero contributors) "<<vtxSPD<<endl;
	  return;
	}
      
      const AliAODVertex* vtxTPC = 0;
      int nVertices = fAOD->GetNumberOfVertices();
      if(debug>0) cout<<"number of vertices is "<<nVertices<<endl;
      for(int iVertices = 0; iVertices < nVertices; iVertices++)
	{
	  const AliAODVertex* vertex = fAOD->GetVertex(iVertices);
	  if(debug>0) cout<<"vertex type is "<<vertex->GetType()<<endl;
	  if(vertex->GetType()!=AliAODVertex::kMainTPC) continue;
	  vtxTPC = vertex;
	}
      if(!vtxTPC||vtxTPC->GetNContributors()<=0)
	{
	  if(debug>0) cout<<"rejecting pileup event (no vtxTPC or zero contributors)"<<vtxTPC<<endl;
	  return;
	}
      
      float diffZ = vtxSPD->GetZ() - vtxTPC->GetZ();
      if(fabs(diffZ)>2.0)
	{
	  if(debug>0) cout<<"rejecting pileup event with vtxTPC "<<vtxTPC->GetZ()<<" vtxSPD "<<vtxSPD->GetZ()<<endl;
	  return;
	}
      */
      // 
      // ---
      //
    }


  int ntrk = 0;
  int ntrkpos = 0;
  int ntrkneg = 0;
  int ntrkL = 0;
  int ntrkposL = 0;
  int ntrknegL = 0;
  int ntrkR = 0;
  int ntrkposR = 0;
  int ntrknegR = 0;
  int ntrkA1 = 0;
  int ntrkposA1 = 0;
  int ntrknegA1 = 0;
  int ntrkA2 = 0;
  int ntrkposA2 = 0;
  int ntrknegA2 = 0;
  int ntrkMC = 0;
  int ntrkposMC = 0;
  int ntrknegMC = 0;

  int cutntrk = 0;
  int cutntrkpos = 0;
  int cutntrkneg = 0;
  int cutntrkL = 0;
  int cutntrkposL = 0;
  int cutntrknegL = 0;
  int cutntrkR = 0;
  int cutntrkposR = 0;
  int cutntrknegR = 0;
  int cutntrkMC = 0;
  int cutntrkposMC = 0;
  int cutntrknegMC = 0;

  const int hmax = 9; // number of harmonics including 0 for multiplicity counting
  float tpcXplo[hmax], tpcYplo[hmax], tpcXpro[hmax], tpcYpro[hmax];
  float tpcXnlo[hmax], tpcYnlo[hmax], tpcXnro[hmax], tpcYnro[hmax];
  float tpcXpli[hmax], tpcYpli[hmax], tpcXpri[hmax], tpcYpri[hmax];
  float tpcXnli[hmax], tpcYnli[hmax], tpcXnri[hmax], tpcYnri[hmax];
  float tpcXpl[hmax], tpcYpl[hmax], tpcXpr[hmax], tpcYpr[hmax];
  float tpcXnl[hmax], tpcYnl[hmax], tpcXnr[hmax], tpcYnr[hmax];
  for(int i=0; i<hmax;i++)
    {
      tpcXplo[i] = 0.0;
      tpcYplo[i] = 0.0;
      tpcXpro[i] = 0.0;
      tpcYpro[i] = 0.0;
      tpcXnlo[i] = 0.0;
      tpcYnlo[i] = 0.0;
      tpcXnro[i] = 0.0;
      tpcYnro[i] = 0.0;
      //
      tpcXpli[i] = 0.0;
      tpcYpli[i] = 0.0;
      tpcXpri[i] = 0.0;
      tpcYpri[i] = 0.0;
      tpcXnli[i] = 0.0;
      tpcYnli[i] = 0.0;
      tpcXnri[i] = 0.0;
      tpcYnri[i] = 0.0;
      //
      tpcXpl[i] = 0.0;
      tpcYpl[i] = 0.0;
      tpcXpr[i] = 0.0;
      tpcYpr[i] = 0.0;
      tpcXnl[i] = 0.0;
      tpcYnl[i] = 0.0;
      tpcXnr[i] = 0.0;
      tpcYnr[i] = 0.0;
    }

  float MCtpcXplo[hmax], MCtpcYplo[hmax], MCtpcXpro[hmax], MCtpcYpro[hmax];
  float MCtpcXnlo[hmax], MCtpcYnlo[hmax], MCtpcXnro[hmax], MCtpcYnro[hmax];
  float MCtpcXpli[hmax], MCtpcYpli[hmax], MCtpcXpri[hmax], MCtpcYpri[hmax];
  float MCtpcXnli[hmax], MCtpcYnli[hmax], MCtpcXnri[hmax], MCtpcYnri[hmax];
  float MCtpcXpl[hmax], MCtpcYpl[hmax], MCtpcXpr[hmax], MCtpcYpr[hmax];
  float MCtpcXnl[hmax], MCtpcYnl[hmax], MCtpcXnr[hmax], MCtpcYnr[hmax];
  for(int i=0; i<hmax;i++)
    {
      MCtpcXplo[i] = 0.0;
      MCtpcYplo[i] = 0.0;
      MCtpcXpro[i] = 0.0;
      MCtpcYpro[i] = 0.0;
      MCtpcXnlo[i] = 0.0;
      MCtpcYnlo[i] = 0.0;
      MCtpcXnro[i] = 0.0;
      MCtpcYnro[i] = 0.0;
      //
      MCtpcXpli[i] = 0.0;
      MCtpcYpli[i] = 0.0;
      MCtpcXpri[i] = 0.0;
      MCtpcYpri[i] = 0.0;
      MCtpcXnli[i] = 0.0;
      MCtpcYnli[i] = 0.0;
      MCtpcXnri[i] = 0.0;
      MCtpcYnri[i] = 0.0;
      //
      MCtpcXpl[i] = 0.0;
      MCtpcYpl[i] = 0.0;
      MCtpcXpr[i] = 0.0;
      MCtpcYpr[i] = 0.0;
      MCtpcXnl[i] = 0.0;
      MCtpcYnl[i] = 0.0;
      MCtpcXnr[i] = 0.0;
      MCtpcYnr[i] = 0.0;
    }

  // ---------------------------------- //
  // --- Now looping over MC tracks --- //
  // ---------------------------------- //
  if(fMC)
    {
      // variables for 3rd particle
      float MCpt3[d_ntrkMC];
      float MCeta3[d_ntrkMC];
      //float MCphi3[d_ntrkMC];//new...
      int MCcharge3[d_ntrkMC];
      // initial variables to out of range values
      // for cases when the track loop skips certain tracks
      for(int i=0; i<d_ntrk; i++)
	{
	  MCpt3[i] = -99;
	  MCeta3[i] = -99;
	  //MCphi3[i] = -99;
	  MCcharge3[i] = 0;
	}
      // --- first MC track loop
      for(int itrkMC = 0; itrkMC<d_ntrkMC; itrkMC++)
	{
	  //AliVParticle *trackMC = fMC->GetTrack(itrkMC);
	  AliAODMCParticle *trackMC = (AliAODMCParticle *)fMC->GetTrack(itrkMC);
	  if(!trackMC)
	    {
	      if(debug>0) cout<<"ERROR: Could not receive track "<<itrkMC<<" (mc loop)"<<endl;
	      continue;
	    }
	  
	  float pt = trackMC->Pt();
	  fHistPtMC->Fill(pt);
	  
	  float phi = trackMC->Phi();
	  float eta = trackMC->Eta();
	  int charge = trackMC->Charge();
	  if(charge==0) continue; // new
	  charge /= 3; // new
	  // int ncls = trackMC->GetMCTPCNcls();
	  // float dedx = trackMC->GetMCTPCsignal();
	  // float dcaxy = trackMC->DCA();
	  // float dcaz = trackMC->ZAtDCA();
	  bool pos = charge>0.0;
	  bool neg = charge<0.0;
	  
	  fHistPhiMC->Fill(phi);
	  fHistEtaMC->Fill(eta);
	  fHistChargeMC->Fill(charge);
	  // fHistMCTPCnclsMC->Fill(ncls);
	  // fHistDedxMC->Fill(dedx);
	  // fHistDCAxyMC->Fill(dcaxy);
	  // fHistDCAzMC->Fill(dcaz);

	  // apply kinematic cuts
	  if(fabs(eta)>outeta) continue; 
	  if(pt<ptmin||pt>ptmax) continue;

	  //if(charge!=-3&&charge!=+3) continue;// x3 by convention
	  if(!trackMC->IsPrimary()) continue;
	  if(!trackMC->IsPhysicalPrimary()) continue;
	  //if(abs(p0->GetPdgCode())==11) continue; //electrons
	  //if(abs(p0->GetPdgCode())==13) continue; //electrons

	  ntrkMC++;
	  if(pos) ntrkposMC++;
	  if(neg) ntrknegMC++;

	  if(pt>ptmin&&pt<ptmax&&fabs(eta)<outeta&&fabs(eta)>excleta)
	    {
	      cutntrkMC++;
	      if(pos) cutntrkposMC++;
	      if(neg) cutntrknegMC++;
	    }

	  if(pt>ptmin&&pt<ptmax)
	    {
	      // --- outer
	      if(eta>-outeta&&eta<-ineta)
		{
		  for(int i=1; i<hmax;i++)
		    {
		      if(charge>0)
			{
			  MCtpcXplo[i] += cos(i*phi);
			  MCtpcYplo[i] += sin(i*phi);
			  MCtpcXpl[i] += cos(i*phi);
			  MCtpcYpl[i] += sin(i*phi);
			}
		      if(charge<0)
			{
			  MCtpcXnlo[i] += cos(i*phi);
			  MCtpcYnlo[i] += sin(i*phi);
			  MCtpcXnl[i] += cos(i*phi);
			  MCtpcYnl[i] += sin(i*phi);
			}
		    }
		  if(charge>0)
		    {
		      MCtpcXplo[0] += 1.0;
		      MCtpcYplo[0] += pt;
		      MCtpcXpl[0] += 1.0;
		      MCtpcYpl[0] += pt;
		    }	
		  if(charge<0)
		    {
		      MCtpcXnlo[0] += 1.0;
		      MCtpcYnlo[0] += pt;
		      MCtpcXnl[0] += 1.0;
		      MCtpcYnl[0] += pt;
		    }
		} // end left half of MCTPC
	      //----------------------------
	      if(eta>ineta&&eta<outeta)
		{
		  for(int i=1; i<hmax; i++)
		    {
		      if(charge>0)
			{
			  MCtpcXpro[i] += cos(i*phi);
			  MCtpcYpro[i] += sin(i*phi);
			  MCtpcXpr[i] += cos(i*phi);
			  MCtpcYpr[i] += sin(i*phi);
			}
		      if(charge<0)
			{
			  MCtpcXnro[i] += cos(i*phi);
			  MCtpcYnro[i] += sin(i*phi);
			  MCtpcXnr[i] += cos(i*phi);
			  MCtpcYnr[i] += sin(i*phi);
			}
		    }
		  if(charge>0)
		    {
		      MCtpcXpro[0] += 1.0;
		      MCtpcYpro[0] += pt;
		      MCtpcXpr[0] += 1.0;
		      MCtpcYpr[0] += pt;
		    }	
		  if(charge<0)
		    {
		      MCtpcXnro[0] += 1.0;
		      MCtpcYnro[0] += pt;
		      MCtpcXnr[0] += 1.0;
		      MCtpcYnr[0] += pt;
		    }
		} // end right half of MCTPC
	      // --- inner
	      if(eta>-ineta&&eta<-excleta)
		{
		  for(int i=1; i<hmax;i++)
		    {
		      if(charge>0)
			{
			  MCtpcXpli[i] += cos(i*phi);
			  MCtpcYpli[i] += sin(i*phi);
			  MCtpcXpl[i] += cos(i*phi);
			  MCtpcYpl[i] += sin(i*phi);
			}
		      if(charge<0)
			{
			  MCtpcXnli[i] += cos(i*phi);
			  MCtpcYnli[i] += sin(i*phi);
			  MCtpcXnl[i] += cos(i*phi);
			  MCtpcYnl[i] += sin(i*phi);
			}
		    }
		  if(charge>0)
		    {
		      MCtpcXpli[0] += 1.0;
		      MCtpcYpli[0] += pt;
		      MCtpcXpl[0] += 1.0;
		      MCtpcYpl[0] += pt;
		    }	
		  if(charge<0)
		    {
		      MCtpcXnli[0] += 1.0;
		      MCtpcYnli[0] += pt;
		      MCtpcXnl[0] += 1.0;
		      MCtpcYnl[0] += pt;
		    }
		} // end left half of MCTPC
	      //----------------------------
	      if(eta>excleta&&eta<ineta)
		{
		  for(int i=1; i<hmax; i++)
		    {
		      if(charge>0)
			{
			  MCtpcXpri[i] += cos(i*phi);
			  MCtpcYpri[i] += sin(i*phi);
			  MCtpcXpr[i] += cos(i*phi);
			  MCtpcYpr[i] += sin(i*phi);
			}
		      if(charge<0)
			{
			  MCtpcXnri[i] += cos(i*phi);
			  MCtpcYnri[i] += sin(i*phi);
			  MCtpcXnr[i] += cos(i*phi);
			  MCtpcYnr[i] += sin(i*phi);
			}
		    }
		  if(charge>0)
		    {
		      MCtpcXpri[0] += 1.0;
		      MCtpcYpri[0] += pt;
		      MCtpcXpr[0] += 1.0;
		      MCtpcYpr[0] += pt;
		    }	
		  if(charge<0)
		    {
		      MCtpcXnri[0] += 1.0;
		      MCtpcYnri[0] += pt;
		      MCtpcXnr[0] += 1.0;
		      MCtpcYnr[0] += pt;
		    }
		} // end right half of MCTPC
	    } // end pt selection
	  
	  h_AT_etaMC->Fill(cent,charge);
	  h2_AT_etaMC->Fill(cent,charge);

	  MCpt3[itrkMC] = pt;
	  MCeta3[itrkMC] = eta;
	  //MCphi3[itrkMC] = phi;
	  MCcharge3[itrkMC] = charge;
	  
	} // end of first MC track loop 
      
      // if(debug>0)
      // 	{
      // 	  cout<<"number of tracks with wrong filter bit "<<badbit<<endl;
      // 	  cout<<"difference = ntrk - nbad = "<<d_ntrk-badbit<<endl;
      // 	}
      float MCtpcX[9], MCtpcY[9], MCtpcQQ[9];//, qq[9];
      float MCtpcXp[9], MCtpcYp[9], MCtpcQQp[9];//, qqp[9]; // pos
      float MCtpcXn[9], MCtpcYn[9], MCtpcQQn[9];//, qqn[9]; // neg
      
      float MCqasymm = float(ntrkposMC-ntrknegMC)/float(ntrkMC);
      if(doacuts)
	{
	  MCqasymm = float(cutntrkposMC-cutntrknegMC)/float(cutntrkMC);
	}
      h_A_centMC->Fill(cent,MCqasymm);
      h2_A_centMC->Fill(cent,MCqasymm);
      
      
      for(int i=0; i<6; i++)
	{
	  MCtpcX[i]=MCtpcXpl[i]+MCtpcXnl[i]+MCtpcXpr[i]+MCtpcXnr[i];
	  MCtpcY[i]=MCtpcYpl[i]+MCtpcYnl[i]+MCtpcYpr[i]+MCtpcYnr[i];
	  MCtpcQQ[i]=MCtpcX[i]*MCtpcX[i]+MCtpcY[i]*MCtpcY[i];
	  //qq[i]=sqrt(MCtpcQQ[i]/MCtpcX[0]);
	  // pos
	  MCtpcXp[i]=MCtpcXpl[i]+MCtpcXpr[i];
	  MCtpcYp[i]=MCtpcYpl[i]+MCtpcYpr[i];
	  MCtpcQQp[i]=MCtpcXp[i]*MCtpcXp[i]+MCtpcYp[i]*MCtpcYp[i];
	  //qqp[i]=sqrt(MCtpcQQp[i]/MCtpcXp[0]);
	  // neg
	  MCtpcXn[i]=MCtpcXnl[i]+MCtpcXnr[i];
	  MCtpcYn[i]=MCtpcYnl[i]+MCtpcYnr[i];
	  MCtpcQQn[i]=MCtpcXn[i]*MCtpcXn[i]+MCtpcYn[i]*MCtpcYn[i];
	  //qqn[i]=sqrt(MCtpcQQn[i]/MCtpcXn[0]);
	}
      
      float M = MCtpcX[0];
      float W_2 = M*(M-1);
      float Mp = MCtpcXp[0];
      float Wp_2 = Mp*(Mp-1);
      float Mn = MCtpcXn[0];
      float Wn_2 = Mn*(Mn-1);
      
      float MCtpcXl2 = MCtpcXnl[2]+MCtpcXpl[2];
      float MCtpcYl2 = MCtpcYnl[2]+MCtpcYpl[2];
      float MCtpcXr2 = MCtpcXnr[2]+MCtpcXpr[2];
      float MCtpcYr2 = MCtpcYnr[2]+MCtpcYpr[2];
      float MCtpcXl0 = MCtpcXnl[0]+MCtpcXpl[0];
      float MCtpcXr0 = MCtpcXnr[0]+MCtpcXpr[0];
      
      float MCtpcXl2o = MCtpcXnlo[2]+MCtpcXplo[2];
      float MCtpcYl2o = MCtpcYnlo[2]+MCtpcYplo[2];
      float MCtpcXr2o = MCtpcXnro[2]+MCtpcXpro[2];
      float MCtpcYr2o = MCtpcYnro[2]+MCtpcYpro[2];
      float MCtpcXl0o = MCtpcXnlo[0]+MCtpcXplo[0];
      float MCtpcXr0o = MCtpcXnro[0]+MCtpcXpro[0];
      
      float MCq22ev = (MCtpcQQ[2]-M)/W_2;
      float MCq22Pev = (MCtpcQQp[2]-Mp)/Wp_2;
      float MCq22Nev = (MCtpcQQn[2]-Mn)/Wn_2;
      float MCq22gap0ev = (MCtpcXl2*MCtpcXr2+MCtpcYl2*MCtpcYr2)/(MCtpcXl0*MCtpcXr0);
      float MCq22gap0Pev = (MCtpcXpl[2]*MCtpcXpr[2]+MCtpcYpl[2]*MCtpcYpr[2])/(MCtpcXpl[0]*MCtpcXpr[0]);
      float MCq22gap0Nev = (MCtpcXnl[2]*MCtpcXnr[2]+MCtpcYnl[2]*MCtpcYnr[2])/(MCtpcXnl[0]*MCtpcXnr[0]);
      float MCq22gap1ev = (MCtpcXl2o*MCtpcXr2o+MCtpcYl2o*MCtpcYr2o)/(MCtpcXl0o*MCtpcXr0o);
      float MCq22gap1Pev = (MCtpcXplo[2]*MCtpcXpro[2]+MCtpcYplo[2]*MCtpcYpro[2])/(MCtpcXplo[0]*MCtpcXpro[0]);
      float MCq22gap1Nev = (MCtpcXnlo[2]*MCtpcXnro[2]+MCtpcYnlo[2]*MCtpcYnro[2])/(MCtpcXnlo[0]*MCtpcXnro[0]);
      h_MCq22_cent->Fill(cent,MCq22ev);
      h_MCq22_centP->Fill(cent,MCq22Pev);
      h_MCq22_centN->Fill(cent,MCq22Nev);
      h_MCAq22_cent->Fill(cent,MCq22ev*MCqasymm);
      h_MCAq22_centP->Fill(cent,MCq22Pev*MCqasymm);
      h_MCAq22_centN->Fill(cent,MCq22Nev*MCqasymm);
      h_MCq22gap0_cent->Fill(cent,MCq22gap0ev);
      h_MCq22gap0_centP->Fill(cent,MCq22gap0Pev);
      h_MCq22gap0_centN->Fill(cent,MCq22gap0Nev);
      h_MCq22gap1_cent->Fill(cent,MCq22gap1ev);
      h_MCq22gap1_centP->Fill(cent,MCq22gap1Pev);
      h_MCq22gap1_centN->Fill(cent,MCq22gap1Nev);
      
      
      
      if(debug>0) cout<<"there are "<<ntrkMC<<" Monte Carlo tracks within the cuts in this event"<<endl;
      
      // --- fsecond MC track loop
      for(int itrkMC = 0; itrkMC<d_ntrkMC; itrkMC++)
	{
	  //AliVParticle *trackMC = fMC->GetTrack(itrkMC);
	  AliAODMCParticle *trackMC = (AliAODMCParticle *)fMC->GetTrack(itrkMC);
	  if(!trackMC)
	    {
	      if(debug>0) cout<<"ERROR: Could not receive track "<<itrkMC<<" (mc loop)"<<endl;
	      continue;
	    }
	  
	  float pt = trackMC->Pt();
	  
	  float phi = trackMC->Phi();
	  float eta = trackMC->Eta();
	  int charge = trackMC->Charge();
	  if(charge==0) continue; // new
	  charge /= 3; // new
	  bool pos = charge>0.0;
	  bool neg = charge<0.0;
	  
	  
	  // apply kinematic cuts
	  if(fabs(eta)>outeta) continue; 
	  if(pt<ptmin||pt>ptmax) continue;

	  //if(charge!=-3&&charge!=+3) continue;// x3 by convention
	  if(!trackMC->IsPrimary()) continue;
	  if(!trackMC->IsPhysicalPrimary()) continue;

	  //if(abs(p0->GetPdgCode())==11) continue; //electrons
	  //if(abs(p0->GetPdgCode())==13) continue; //electrons


	  float MCtrkXpl[hmax], MCtrkYpl[hmax], MCtrkXpr[hmax], MCtrkYpr[hmax];
	  float MCtrkXnl[hmax], MCtrkYnl[hmax], MCtrkXnr[hmax], MCtrkYnr[hmax];
	  for(int i=0; i<hmax;i++)
	    {
	      MCtrkXpl[i] = 0.0;
	      MCtrkYpl[i] = 0.0;
	      MCtrkXpr[i] = 0.0;
	      MCtrkYpr[i] = 0.0;
	      MCtrkXnl[i] = 0.0;
	      MCtrkYnl[i] = 0.0;
	      MCtrkXnr[i] = 0.0;
	      MCtrkYnr[i] = 0.0;
	    }
	  // very similar to above, except += changed to =
	  if(pt>ptmin&&pt<ptmax)
	    {
	      if(eta>-outeta&&eta<-excleta)
		{
		  // assign values to array elements 1..8
		  for(int i=1; i<hmax;i++)
		    {
		      if(charge>0)
			{
			  MCtrkXpl[i] = cos(i*phi);
			  MCtrkYpl[i] = sin(i*phi);
			}
		      if(charge<0)
			{
			  MCtrkXnl[i] = cos(i*phi);
			  MCtrkYnl[i] = sin(i*phi);
			}
		    }
		  // assign values to array element 0
		  if(charge>0)
		    {
		      MCtrkXpl[0] = 1.0;
		      MCtrkYpl[0] = pt;
		    }	
		  if(charge<0)
		    {
		      MCtrkXnl[0] = 1.0;
		      MCtrkYnl[0] = pt;
		    }
		}
	      //----------------------------
	      if(eta>excleta&&eta<outeta) // right half of TPC
		{
		  for(int i=1; i<hmax; i++)
		    {
		      if(charge>0)
			{
			  MCtrkXpr[i] = cos(i*phi);
			  MCtrkYpr[i] = sin(i*phi);
			}
		      if(charge<0)
			{
			  MCtrkXnr[i] = cos(i*phi);
			  MCtrkYnr[i] = sin(i*phi);
			}
		    }
		  if(charge>0)
		    {
		      MCtrkXpr[0] = 1.0;
		      MCtrkYpr[0] = pt;
		    }	
		  if(charge<0)
		    {
		      MCtrkXnr[0] = 1.0;
		      MCtrkYnr[0] = pt;
		    }
		} // end right half of TPC
	    } // end pT selection
	  
	  // COME BACK HERE
	  float MCtrkX[9];//,  MCtrkY[9];
	  float MCtrkXp[9], MCtrkYp[9]; // pos
	  float MCtrkXn[9], MCtrkYn[9]; // neg
	  for(int i=0; i<6; i++)
	    {
	      MCtrkX[i]=MCtrkXpl[i]+MCtrkXnl[i]+MCtrkXpr[i]+MCtrkXnr[i];
	      //MCtrkY[i]=MCtrkYpl[i]+MCtrkYnl[i]+MCtrkYpr[i]+MCtrkYnr[i];
	      // pos
	      MCtrkXp[i]=MCtrkXpl[i]+MCtrkXpr[i];
	      MCtrkYp[i]=MCtrkYpl[i]+MCtrkYpr[i];
	      // neg
	      MCtrkXn[i]=MCtrkXnl[i]+MCtrkXnr[i];
	      MCtrkYn[i]=MCtrkYnl[i]+MCtrkYnr[i];
	    }
	  
	  // sanity checks to prevent division by zero and other problems
	  if(MCtrkX[0]<1) continue;
	  if(MCtpcX[0]<2) continue;
	  if(MCtpcXp[0]<2) continue;
	  if(MCtpcXn[0]<2) continue;
	  
	  // --- calculate differential cumulants from Q-vector components
	  // --- reference flow is from whole TPC
	  //float MCdiffq22ev = ((MCtrkX[2]*MCtpcX[2])+(MCtrkY[2]*MCtpcY[2])-1)/(MCtpcX[0]-1);
	  float MCdiffq22Pev = ((MCtrkXp[2]*MCtpcX[2])+(MCtrkYp[2]*MCtpcY[2])-1)/(MCtpcX[0]-1);
	  float MCdiffq22Nev = ((MCtrkXn[2]*MCtpcX[2])+(MCtrkYn[2]*MCtpcY[2])-1)/(MCtpcX[0]-1);
	  // --- third harmonic
	  //float MCdiffq32ev = ((MCtrkX[3]*MCtpcX[3])+(MCtrkY[3]*MCtpcY[3])-1)/(MCtpcX[0]-1);
	  float MCdiffq32Pev = ((MCtrkXp[3]*MCtpcX[3])+(MCtrkYp[3]*MCtpcY[3])-1)/(MCtpcX[0]-1);
	  float MCdiffq32Nev = ((MCtrkXn[3]*MCtpcX[3])+(MCtrkYn[3]*MCtpcY[3])-1)/(MCtpcX[0]-1);
	  // --- fourth harmonic
	  //float MCdiffq42ev = ((MCtrkX[4]*MCtpcX[4])+(MCtrkY[4]*MCtpcY[4])-1)/(MCtpcX[0]-1);
	  float MCdiffq42Pev = ((MCtrkXp[4]*MCtpcX[4])+(MCtrkYp[4]*MCtpcY[4])-1)/(MCtpcX[0]-1);
	  float MCdiffq42Nev = ((MCtrkXn[4]*MCtpcX[4])+(MCtrkYn[4]*MCtpcY[4])-1)/(MCtpcX[0]-1);


	  if(donested)
	    {
	      for(int i=0; i<d_ntrkMC; i++)
		{
		  // ----------------------------------------------------------------------------------------
		  // --- want delta eta and delta pt dependence of differential cumulant weighted with charge
		  // --- need nested track loop to get eta1 and eta3
		  // ----------------------------------------------------------------------------------------
		  // make sure eta is in range, -99 should be autorejected by this
		  //if(MCeta3[i]==-99) continue;
		  if(MCeta3[i]<-outeta||MCeta3[i]>outeta)
		    {
		      //cout<<"MCeta3 out of range!!! "<<MCeta3[i]<<endl;
		      continue; 
		    }
		  if(eta<-outeta||eta>outeta)
		    {
		      cout<<"eta out of range!!! "<<eta<<endl;
		      continue; 
		    }
		  if(MCpt3[i]<ptmin||MCpt3[i]>ptmax)
		    {
		      //cout<<"MCpt3 out of range!!! "<<MCpt3[i]<<endl;
		      continue;
		    }
		  if(pt<ptmin||pt>ptmax)
		    {
		      cout<<"pt out of range!!!"<<endl;
		      continue;
		    }
		  // charge==0 should be rejected by the above cut???
		  if(MCcharge3[i]!=1&&MCcharge3[i]!=-1)
		    {
		      //cout<<"WTF!"<<endl;
		      continue;
		    }
		  // VERY IMPORTANT remove auto-correlations with 3rd particle
		  if(i==itrkMC) continue;
		  
		  float DETA = eta-MCeta3[i];
		  //float DPHI = phi-MCphi3[i];
		  //float DPT = pt-MCpt3[i];

		  // ---

		  hMC_AT_X_deta->Fill(DETA,MCcharge3[i]);
		  if(pos) hMC_AT_X_detaP->Fill(DETA,MCcharge3[i]);
		  if(neg) hMC_AT_X_detaN->Fill(DETA,MCcharge3[i]);
		  if(pos) hMC_diffq22_X_detaP->Fill(DETA,MCdiffq22Pev);
		  if(neg) hMC_diffq22_X_detaN->Fill(DETA,MCdiffq22Nev);
		  if(pos) hMC_diffq32_X_detaP->Fill(DETA,MCdiffq32Pev);
		  if(neg) hMC_diffq32_X_detaN->Fill(DETA,MCdiffq32Nev);
		  if(pos) hMC_diffq42_X_detaP->Fill(DETA,MCdiffq42Pev);
		  if(neg) hMC_diffq42_X_detaN->Fill(DETA,MCdiffq42Nev);
		  if(pos) hMC_ATdiffq22_X_detaP->Fill(DETA,MCdiffq22Pev*MCcharge3[i]);
		  if(neg) hMC_ATdiffq22_X_detaN->Fill(DETA,MCdiffq22Nev*MCcharge3[i]);
		  if(pos) hMC_ATdiffq32_X_detaP->Fill(DETA,MCdiffq32Pev*MCcharge3[i]);
		  if(neg) hMC_ATdiffq32_X_detaN->Fill(DETA,MCdiffq32Nev*MCcharge3[i]);
		  if(pos) hMC_ATdiffq42_X_detaP->Fill(DETA,MCdiffq42Pev*MCcharge3[i]);
		  if(neg) hMC_ATdiffq42_X_detaN->Fill(DETA,MCdiffq42Nev*MCcharge3[i]);

		} // nested track loop

	    } // check on donested

	} // end of second MC track loop 
      
    } // check on existence of MC
  
  
  
  // ------------------------------- //
  // --- Now looping over tracks --- //
  // ------------------------------- //
  
  // --- track loop for mapping matrix
  TExMap *trackMap = new TExMap();
  for(int itrk=0; itrk<d_ntrk; itrk++)
    {
      AliAODTrack *track = (AliAODTrack *)fAOD->GetTrack(itrk);
      if(!track)
  	{
  	  if(debug>0) cout<<"ERROR: Could not retrieve AODtrack "<<itrk<<endl;
  	  continue;
  	}
      int gid = track->GetID();
      if(track->TestFilterBit(fbit)) trackMap->Add(gid,itrk);
    }



  // variables for 3rd particle
  float pt3[d_ntrk];
  float eta3[d_ntrk];
  float phi3[d_ntrk];//new...
  int charge3[d_ntrk];
  // initial variables to out of range values
  // for cases when the track loop skips certain tracks
  for(int i=0; i<d_ntrk; i++)
    {
      pt3[i] = -99;
      eta3[i] = -99;
      phi3[i] = -99;
      charge3[i] = 0;
    }
  
  int badbit = 0; // counter for number of tracks with wrong filter bit
  int dcacoverrors = 0;
  // --- main track loop
  for(int itrk = 0; itrk<d_ntrk; itrk++)
    {
      AliAODTrack *track = (AliAODTrack *)fAOD->GetTrack(itrk);
      if(!track)
	{
	  if(debug>0) cout<<"ERROR: Could not retrieve AODtrack "<<itrk<<endl;
	  continue;
	}

      if(!track->TestFilterBit(fbit))
	{
	  badbit++; // count tracks with wrong filter bit
	  if(debug>15) cout<<"wrong filter bit for track "<<itrk<<endl;
	  continue;
	}
      
      int gid = track->GetID();
      AliAODTrack *PIDtrack;
      if(gid>=0) PIDtrack = track;
      else PIDtrack = (AliAODTrack *)fAOD->GetTrack(trackMap->GetValue(-1-gid));


      // if(debug>15)
      // 	{
      // 	  cout<<"track index is "<<itrk<<endl;
      // 	  cout<<"track ID is "<<gid<<endl;
      // 	}

      // --- get track variables      
      float pt = track->Pt();
      float phi = track->Phi();
      float eta = track->Eta();
      int charge = track->Charge();
      int ncls = track->GetTPCNcls();
      float dedx = track->GetTPCsignal();

      float Tdcaxy = track->DCA();
      float Tdcaz = track->ZAtDCA();


      float dcaxy = -999;
      float dcax = -999;
      float dcay = -999;
      float dcaz = -999;

      double r[3];
      bool dcaflag = track->GetXYZ(r);
      //cout<<"fbit is "<<fbit<<" GetXYZ is "<<dcaflag<<endl;
      //printf("fbit is %d GetXYZ is %d PropagateToDCA is %d Tdcaz is %f pt is %f \n\n",fbit,dcaflag,proptodca,Tdcaz,pt);

      double DCA[2]; // dca
      double COV[3]; // covariance
      bool proptodca = track->PropagateToDCA(fVtx,mag,100.0,DCA,COV);
      if(!proptodca) {dcacoverrors++;}
      if(!proptodca&&debug>0) {cout<<"No DCACOV for you!"<<endl;}

      if(dcaflag)
      	{
      	  dcaxy = r[0];
      	  dcaz = r[1];
      	  if(debug>5) cout<<"GetXYZ is true, filter bit is "<<fbit<<endl;
      	}
      else
      	{
      	  dcax = r[0] - eventX;
      	  dcay = r[1] - eventY;
      	  dcaz = r[2] - eventZ;
      	  // --- need the sign convention for dcaxy...
      	  dcaxy = sqrt(dcax*dcax+dcay*dcay);
      	  //if((float)dcaxy!=(float)fabs((float)DCA[0])&&debug>4) cout<<"hmm... "<<dcaxy<<" "<<DCA[0]<<endl;
      	  // --- set dcaxy to value from PropagateToDCA to get correct sign
      	  dcaxy = DCA[0];
      	  // --- dcaz on the other hand is unambiguous
      	  //if(dcaz!=(float)DCA[1]&&debug>4) cout<<"hmm... "<<dcaz<<" "<<DCA[1]<<endl;
      	  if(debug>5) cout<<"GetXYZ is false, filter bit is "<<fbit<<endl;
      	}

      if(debug>5)
      	{
      	  cout<<"r[0] is "<<r[0]<<" ";
      	  cout<<"r[1] is "<<r[1]<<" ";
      	  cout<<"r[2] is "<<r[2]<<endl;
      	  cout<<"eventX is "<<eventX<<" ";
      	  cout<<"eventY is "<<eventY<<" ";
      	  cout<<"eventZ is "<<eventZ<<endl;
      	  cout<<"Tdcaxy is "<<Tdcaxy<<" and dcaxy is "<<dcaxy<<" and DCACOV xy is "<<DCA[0]<<endl;
      	  cout<<"Tdcaz is "<<Tdcaz<<" and dcaz is "<<dcaz<<" and DCACOV z is "<<DCA[1]<<endl;
      	}

      if((fbit==128||fbit==272)&&!dcaflag)
      	{
      	  dcaxy = Tdcaxy;
      	  dcaz = Tdcaz;
      	  if(debug>4)
      	    {
      	      cout<<"GetXYZ false but should be true, flipping values"<<endl;
      	      cout<<"Tdcaxy is "<<Tdcaxy<<" and dcaxy is "<<dcaxy<<" and DCACOV xy is "<<DCA[0]<<endl;
      	      cout<<"Tdcaz is "<<Tdcaz<<" and dcaz is "<<dcaz<<" and DCACOV z is "<<DCA[1]<<endl;
      	    }
      	}

      // --- some diagnostic histograms for tracks
      fHistPt->Fill(pt);
      fHistPhi->Fill(phi);
      fHistEta->Fill(eta);
      fHistCharge->Fill(charge);
      fHistTPCncls->Fill(ncls);
      fHistDedx->Fill(dedx);
      fHistDCAxy->Fill(dcaxy);
      fHistDCAz->Fill(dcaz);

      if(charge>0)
	{
	  fHistPosPt->Fill(pt);
	  fHistPosPhi->Fill(phi);
	  fHistPosEta->Fill(eta);
	}
      if(charge<0)
	{
	  fHistNegPt->Fill(pt);
	  fHistNegPhi->Fill(phi);
	  fHistNegEta->Fill(eta);
	}
      fProfMeanChargePt->Fill(pt,charge);
      fProfMeanChargeEta->Fill(eta,charge);

      // --- some extra verbose stuff
      if(debug>8&&(pt<ptmin||pt>ptmax)) cout<<"pt = "<<pt<<" out of range for Q-vectors"<<endl;
      if(debug>9) cout<<"number of tpc clusters is "<<ncls<<endl;

      // --- track cut on number of TPC clusters
      if(ncls<nclscut)
	{
	  if(debug>5) cout<<"number of tpc clusters is too low "<<ncls<<endl;
	  continue;
	}
      
      // --- track cut on dca
      if(fabs(dcaxy)>dcacutxy||fabs(dcaz)>dcacutz)
	{
	  if(debug>7)
	    {
	      cout<<"dca out of range..."<<endl;
	      cout<<"dcaxy cut is "<<dcacutxy<<endl;
	      cout<<"dcaz cut is "<<dcacutz<<endl;
	      cout<<"dcaxy is "<<dcaxy<<endl;
	      cout<<"dcaz is "<<dcaz<<endl;
	    }
	  if(dodcacuts) continue; // problems depending on filter bit...
	}
      // Prabhat
      // fHistPhi->Fill(phi);
      // fHistEta->Fill(eta);
      fHistDCAxyAfter->Fill(dcaxy);
      fHistDCAzAfter->Fill(dcaz);

      float nsigmapion = 0.;
      float nsigmakaon = 0.;
      float nsigmaprot = 0.;
      float nsigmaelec = 0.;

      if(!doMC){
	nsigmapion = fPID->NumberOfSigmasTPC(PIDtrack,(AliPID::EParticleType)AliPID::kPion);
	nsigmakaon = fPID->NumberOfSigmasTPC(PIDtrack,(AliPID::EParticleType)AliPID::kKaon);
	nsigmaprot = fPID->NumberOfSigmasTPC(PIDtrack,(AliPID::EParticleType)AliPID::kProton);
	nsigmaelec = fPID->NumberOfSigmasTPC(PIDtrack,(AliPID::EParticleType)AliPID::kElectron);
      }

      bool isPion = fabs(nsigmapion) <= nspid;
      bool isKaon = fabs(nsigmakaon) <= nspid;
      bool isProt = fabs(nsigmaprot) <= nspid;
      bool isElec = fabs(nsigmaelec) <= nspid;

      bool isPionH = isPion && (!isKaon) && (!isProt);
      bool isProtH = (!isPion) && (!isKaon) && isProt;
      bool isPionL = isPion && (!isKaon) && (!isProt) && (!isElec);
      bool isProtL = (!isPion) && (!isKaon) && isProt && (!isElec);

      if(isPion) fHistPtPion->Fill(pt);
      if(isPionH) fHistPtPionH->Fill(pt);
      if(isPionL) fHistPtPionL->Fill(pt);
      fHistNsigmaPion->Fill(nsigmapion);

      if(isProt) fHistPtProt->Fill(pt);
      if(isProtH) fHistPtProtH->Fill(pt);
      if(isProtL) fHistPtProtL->Fill(pt);
      fHistNsigmaProt->Fill(nsigmapion);

      // --- count track numbers
      ntrk++;
      if(charge>0) ntrkpos++;
      if(charge<0) ntrkneg++;
      if(eta<0)
	{
	  ntrkL++;
	  if(charge>0) ntrkposL++;
	  if(charge<0) ntrknegL++;
	}
      if(eta>0)
	{
	  ntrkR++;
	  if(charge>0) ntrkposR++;
	  if(charge<0) ntrknegR++;
	}
      h_AT_eta->Fill(eta,charge);
      if(mag<-4.0) h_AT_etaF1->Fill(eta,charge);
      if(mag>4.0) h_AT_etaF3->Fill(eta,charge);
      h2_AT_eta->Fill(eta,charge);
      if(mag<-4.0) h2_AT_etaF1->Fill(eta,charge);
      if(mag>4.0) h2_AT_etaF3->Fill(eta,charge);
      //
      if(mag<-4.0&&charge>0) h_eta_pos_F1->Fill(eta);
      if(mag<-4.0&&charge<0) h_eta_neg_F1->Fill(eta);
      if(mag>4.0&&charge>0) h_eta_pos_F3->Fill(eta);
      if(mag>4.0&&charge<0) h_eta_neg_F3->Fill(eta);

      if(doeffcorr)
      	{
      	  if(gRandom) delete gRandom;
      	  gRandom = new TRandom3(0);
      	  gRandom->SetSeed(0);
      	  double rand = gRandom->Rndm();
      	  int ptbin = int(pt*10);
      	  if(ptbin>49) ptbin = 49;
      	  float eff = effTPC[ptbin];
      	  if(fbit==272) eff = effHYB[ptbin];
      	  if(rand>eff)
      	    {
      	      if(debug>5) cout<<"doing random throw for efficiency correction: pt is "<<pt<<" eff is "<<eff<<" rand is "<<rand<<endl;
      	      continue;
      	    }
      	}

      if(pt>ptmin&&pt<ptmax&&fabs(eta)<outeta&&fabs(eta)>excleta)
	{
	  cutntrk++;
	  if(charge>0) cutntrkpos++;
	  if(charge<0) cutntrkneg++;
	  if(eta<0)
	    {
	      cutntrkL++;
	      if(charge>0) cutntrkposL++;
	      if(charge<0) cutntrknegL++;
	    }
	  if(eta>0)
	    {
	      cutntrkR++;
	      if(charge>0) cutntrkposR++;
	      if(charge<0) cutntrknegR++;
	    }
	  h_AT_cut_eta->Fill(eta,charge);
	  if(mag<-4.0) h_AT_cut_etaF1->Fill(eta,charge);
	  if(mag>4.0) h_AT_cut_etaF3->Fill(eta,charge);
	  h2_AT_cut_eta->Fill(eta,charge);
	  if(mag<-4.0) h2_AT_cut_etaF1->Fill(eta,charge);
	  if(mag>4.0) h2_AT_cut_etaF3->Fill(eta,charge);
	  //
	  if(cent>=centlo&&cent<centhi)
	    {
	      if(mag<-4.0&&charge>0) h_cut_eta_pos_F1->Fill(eta);
	      if(mag<-4.0&&charge<0) h_cut_eta_neg_F1->Fill(eta);
	      if(mag>4.0&&charge>0) h_cut_eta_pos_F3->Fill(eta);
	      if(mag>4.0&&charge<0) h_cut_eta_neg_F3->Fill(eta);
	    }
	}
      
      // ------------------------------------ //
      // --- now fill Q-vector components --- //
      // ------------------------------------ //

      if(pt>ptmin&&pt<ptmax)
	{
	  // --- outer
	  if(eta>-outeta&&eta<-ineta)
	    {
	      for(int i=1; i<hmax;i++)
		{
		  if(charge>0)
		    {
		      tpcXplo[i] += cos(i*phi);
		      tpcYplo[i] += sin(i*phi);
		      tpcXpl[i] += cos(i*phi);
		      tpcYpl[i] += sin(i*phi);
		    }
		  if(charge<0)
		    {
		      tpcXnlo[i] += cos(i*phi);
		      tpcYnlo[i] += sin(i*phi);
		      tpcXnl[i] += cos(i*phi);
		      tpcYnl[i] += sin(i*phi);
		    }
		}
	      if(charge>0)
		{
		  tpcXplo[0] += 1.0;
		  tpcYplo[0] += pt;
		  tpcXpl[0] += 1.0;
		  tpcYpl[0] += pt;
		}	
	      if(charge<0)
		{
		  tpcXnlo[0] += 1.0;
		  tpcYnlo[0] += pt;
		  tpcXnl[0] += 1.0;
		  tpcYnl[0] += pt;
		}
	    } // end left half of TPC
	  //----------------------------
	  if(eta>ineta&&eta<outeta)
	    {
	      for(int i=1; i<hmax; i++)
		{
		  if(charge>0)
		    {
		      tpcXpro[i] += cos(i*phi);
		      tpcYpro[i] += sin(i*phi);
		      tpcXpr[i] += cos(i*phi);
		      tpcYpr[i] += sin(i*phi);
		    }
		  if(charge<0)
		    {
		      tpcXnro[i] += cos(i*phi);
		      tpcYnro[i] += sin(i*phi);
		      tpcXnr[i] += cos(i*phi);
		      tpcYnr[i] += sin(i*phi);
		    }
		}
	      if(charge>0)
		{
		  tpcXpro[0] += 1.0;
		  tpcYpro[0] += pt;
		  tpcXpr[0] += 1.0;
		  tpcYpr[0] += pt;
		}	
	      if(charge<0)
		{
		  tpcXnro[0] += 1.0;
		  tpcYnro[0] += pt;
		  tpcXnr[0] += 1.0;
		  tpcYnr[0] += pt;
		}
	    } // end right half of TPC
	  // --- inner
	  if(eta>-ineta&&eta<-excleta)
	    {
	      for(int i=1; i<hmax;i++)
		{
		  if(charge>0)
		    {
		      tpcXpli[i] += cos(i*phi);
		      tpcYpli[i] += sin(i*phi);
		      tpcXpl[i] += cos(i*phi);
		      tpcYpl[i] += sin(i*phi);
		    }
		  if(charge<0)
		    {
		      tpcXnli[i] += cos(i*phi);
		      tpcYnli[i] += sin(i*phi);
		      tpcXnl[i] += cos(i*phi);
		      tpcYnl[i] += sin(i*phi);
		    }
		}
	      if(charge>0)
		{
		  tpcXpli[0] += 1.0;
		  tpcYpli[0] += pt;
		  tpcXpl[0] += 1.0;
		  tpcYpl[0] += pt;
		}	
	      if(charge<0)
		{
		  tpcXnli[0] += 1.0;
		  tpcYnli[0] += pt;
		  tpcXnl[0] += 1.0;
		  tpcYnl[0] += pt;
		}
	    } // end left half of TPC
	  //----------------------------
	  if(eta>excleta&&eta<ineta)
	    {
	      for(int i=1; i<hmax; i++)
		{
		  if(charge>0)
		    {
		      tpcXpri[i] += cos(i*phi);
		      tpcYpri[i] += sin(i*phi);
		      tpcXpr[i] += cos(i*phi);
		      tpcYpr[i] += sin(i*phi);
		    }
		  if(charge<0)
		    {
		      tpcXnri[i] += cos(i*phi);
		      tpcYnri[i] += sin(i*phi);
		      tpcXnr[i] += cos(i*phi);
		      tpcYnr[i] += sin(i*phi);
		    }
		}
	      if(charge>0)
		{
		  tpcXpri[0] += 1.0;
		  tpcYpri[0] += pt;
		  tpcXpr[0] += 1.0;
		  tpcYpr[0] += pt;
		}	
	      if(charge<0)
		{
		  tpcXnri[0] += 1.0;
		  tpcYnri[0] += pt;
		  tpcXnr[0] += 1.0;
		  tpcYnr[0] += pt;
		}
	    } // end right half of TPC
	} // end pt selection

      // set values for 3rd particle variables
      pt3[itrk] = pt;
      eta3[itrk] = eta;
      phi3[itrk] = phi;
      charge3[itrk] = (int)charge;


      // ------------------------------------ //
      // --- that's it for the track loop --- //
      // ------------------------------------ //

    } // end of first track loop 

  if(debug>0)
    {
      cout<<"filter bit is "<<fbit<<endl;
      cout<<"number of tracks with wrong filter bit "<<badbit<<endl;
      cout<<"difference = ntrk - nbad = "<<d_ntrk-badbit<<endl;
      cout<<"number of dcacov errors "<<dcacoverrors<<endl;
    }
  float fuckyou = float(dcacoverrors)/float(d_ntrk);
  fHistDCACOVerrors->Fill(fuckyou);

  float tpcX[9], tpcY[9], tpcQQ[9];//, qq[9];
  float tpcXp[9], tpcYp[9], tpcQQp[9], tpcQQpp[9];//, qqp[9]; // pos
  float tpcXn[9], tpcYn[9], tpcQQn[9], tpcQQnn[9];//, qqn[9]; // neg
  
  for(int i=0; i<6; i++)
    {
      tpcX[i]=tpcXpl[i]+tpcXnl[i]+tpcXpr[i]+tpcXnr[i];
      tpcY[i]=tpcYpl[i]+tpcYnl[i]+tpcYpr[i]+tpcYnr[i];
      tpcQQ[i]=tpcX[i]*tpcX[i]+tpcY[i]*tpcY[i];
      //qq[i]=sqrt(tpcQQ[i]/tpcX[0]);
      // pos
      tpcXp[i]=tpcXpl[i]+tpcXpr[i];
      tpcYp[i]=tpcYpl[i]+tpcYpr[i];
      tpcQQp[i]=tpcXp[i]*tpcX[i]+tpcYp[i]*tpcY[i];
      tpcQQpp[i]=tpcXp[i]*tpcXp[i]+tpcYp[i]*tpcYp[i];
      //qqp[i]=sqrt(tpcQQp[i]/tpcXp[0]);
      // neg
      tpcXn[i]=tpcXnl[i]+tpcXnr[i];
      tpcYn[i]=tpcYnl[i]+tpcYnr[i];
      tpcQQn[i]=tpcXn[i]*tpcX[i]+tpcYn[i]*tpcY[i];
      tpcQQnn[i]=tpcXn[i]*tpcXn[i]+tpcYn[i]*tpcYn[i];
      //qqn[i]=sqrt(tpcQQn[i]/tpcXn[0]);
    }
  




  float M = tpcX[0];
  float W_2 = M*(M-1);
  float Mp = tpcXp[0];
  float Wp_2 = Mp*(Mp-1);
  float Mn = tpcXn[0];
  float Wn_2 = Mn*(Mn-1);

  // prevent division by zero
  if(M<2) return;
  if(Mp<2) return;
  if(Mn<2) return;

  if(ntrkL<1) return;
  if(ntrkR<1) return;

  float qasymm = (float)(ntrkpos-ntrkneg)/ntrk;
  float qasymmL = (float)(ntrkposL-ntrknegL)/ntrkL;
  float qasymmR = (float)(ntrkposR-ntrknegR)/ntrkR;

  if(doacuts)
    {
      qasymm = (float)(cutntrkpos-cutntrkneg)/cutntrk;
      qasymmL = (float)(cutntrkposL-cutntrknegL)/cutntrkL;
      qasymmR = (float)(cutntrkposR-cutntrknegR)/cutntrkR;
    }

  h_A_cent->Fill(cent,qasymm);
  if(mag<-4.0) h_A_centF1->Fill(cent,qasymm);
  if(mag>4.0) h_A_centF3->Fill(cent,qasymm);
  h2_A_cent->Fill(cent,qasymm);
  if(mag<-4.0) h2_A_centF1->Fill(cent,qasymm);
  if(mag>4.0) h2_A_centF3->Fill(cent,qasymm);

  float Xl2 = tpcXnl[2]+tpcXpl[2];
  float Yl2 = tpcYnl[2]+tpcYpl[2];
  float Xr2 = tpcXnr[2]+tpcXpr[2];
  float Yr2 = tpcYnr[2]+tpcYpr[2];
  float Xl3 = tpcXnl[3]+tpcXpl[3];
  float Yl3 = tpcYnl[3]+tpcYpl[3];
  float Xr3 = tpcXnr[3]+tpcXpr[3];
  float Yr3 = tpcYnr[3]+tpcYpr[3];
  float Xl4 = tpcXnl[4]+tpcXpl[4];
  float Yl4 = tpcYnl[4]+tpcYpl[4];
  float Xr4 = tpcXnr[4]+tpcXpr[4];
  float Yr4 = tpcYnr[4]+tpcYpr[4];
  float Xl0 = tpcXnl[0]+tpcXpl[0];
  float Xr0 = tpcXnr[0]+tpcXpr[0];
  
  float Xl2o = tpcXnlo[2]+tpcXplo[2];
  float Yl2o = tpcYnlo[2]+tpcYplo[2];
  float Xr2o = tpcXnro[2]+tpcXpro[2];
  float Yr2o = tpcYnro[2]+tpcYpro[2];
  float Xl3o = tpcXnlo[3]+tpcXplo[3];
  float Yl3o = tpcYnlo[3]+tpcYplo[3];
  float Xr3o = tpcXnro[3]+tpcXpro[3];
  float Yr3o = tpcYnro[3]+tpcYpro[3];
  float Xl4o = tpcXnlo[4]+tpcXplo[4];
  float Yl4o = tpcYnlo[4]+tpcYplo[4];
  float Xr4o = tpcXnro[4]+tpcXpro[4];
  float Yr4o = tpcYnro[4]+tpcYpro[4];
  float Xl0o = tpcXnlo[0]+tpcXplo[0];
  float Xr0o = tpcXnro[0]+tpcXpro[0];

  // --- q vector components for recentering if needed
  h_X22_cent->Fill(cent,(tpcX[2]/tpcX[0]));
  h_X22_centP->Fill(cent,(tpcXp[2]/tpcXp[0]));
  h_X22_centN->Fill(cent,(tpcXn[2]/tpcXn[0]));
  h_Y22_cent->Fill(cent,(tpcY[2]/tpcX[0]));
  h_Y22_centP->Fill(cent,(tpcYp[2]/tpcXp[0]));
  h_Y22_centN->Fill(cent,(tpcYn[2]/tpcXn[0]));
  // --- outer left
  h_X22_lo_cent->Fill(cent,(tpcX[2]/tpcX[0]));
  h_X22_lo_centP->Fill(cent,(tpcXplo[2]/tpcXplo[0]));
  h_X22_lo_centN->Fill(cent,(tpcXnlo[2]/tpcXnlo[0]));
  h_Y22_lo_cent->Fill(cent,(tpcY[2]/tpcX[0]));
  h_Y22_lo_centP->Fill(cent,(tpcYplo[2]/tpcXplo[0]));
  h_Y22_lo_centN->Fill(cent,(tpcYnlo[2]/tpcXnlo[0]));
  // --- inner left
  h_X22_li_cent->Fill(cent,(tpcX[2]/tpcX[0]));
  h_X22_li_centP->Fill(cent,(tpcXpli[2]/tpcXpli[0]));
  h_X22_li_centN->Fill(cent,(tpcXnli[2]/tpcXnli[0]));
  h_Y22_li_cent->Fill(cent,(tpcY[2]/tpcX[0]));
  h_Y22_li_centP->Fill(cent,(tpcYpli[2]/tpcXpli[0]));
  h_Y22_li_centN->Fill(cent,(tpcYnli[2]/tpcXnli[0]));
  // --- inner right
  h_X22_ri_cent->Fill(cent,(tpcX[2]/tpcX[0]));
  h_X22_ri_centP->Fill(cent,(tpcXpri[2]/tpcXpri[0]));
  h_X22_ri_centN->Fill(cent,(tpcXnri[2]/tpcXnri[0]));
  h_Y22_ri_cent->Fill(cent,(tpcY[2]/tpcX[0]));
  h_Y22_ri_centP->Fill(cent,(tpcYpri[2]/tpcXpri[0]));
  h_Y22_ri_centN->Fill(cent,(tpcYnri[2]/tpcXnri[0]));
  // --- outer right
  h_X22_ro_cent->Fill(cent,(tpcX[2]/tpcX[0]));
  h_X22_ro_centP->Fill(cent,(tpcXpro[2]/tpcXpro[0]));
  h_X22_ro_centN->Fill(cent,(tpcXnro[2]/tpcXnro[0]));
  h_Y22_ro_cent->Fill(cent,(tpcY[2]/tpcX[0]));
  h_Y22_ro_centP->Fill(cent,(tpcYpro[2]/tpcXpro[0]));
  h_Y22_ro_centN->Fill(cent,(tpcYnro[2]/tpcXnro[0]));
  
  float q22ev = (tpcQQ[2]-M)/W_2;
  float q22Pev = (tpcQQp[2]-Mp)/(Mp*M-Mp);
  float q22Nev = (tpcQQn[2]-Mn)/(Mn*M-Mn);
  float q22PPev = (tpcQQpp[2]-Mp)/Wp_2;
  float q22NNev = (tpcQQnn[2]-Mn)/Wn_2;
  float q22gap0ev = (Xl2*Xr2+Yl2*Yr2)/(Xl0*Xr0);
  float q22gap0Pev = (tpcXpl[2]*tpcXpr[2]+tpcYpl[2]*tpcYpr[2])/(tpcXpl[0]*tpcXpr[0]);
  float q22gap0Nev = (tpcXnl[2]*tpcXnr[2]+tpcYnl[2]*tpcYnr[2])/(tpcXnl[0]*tpcXnr[0]);
  float q22gap1ev = (Xl2o*Xr2o+Yl2o*Yr2o)/(Xl0o*Xr0o);
  float q22gap1Pev = (tpcXplo[2]*tpcXpro[2]+tpcYplo[2]*tpcYpro[2])/(tpcXplo[0]*tpcXpro[0]);
  float q22gap1Nev = (tpcXnlo[2]*tpcXnro[2]+tpcYnlo[2]*tpcYnro[2])/(tpcXnlo[0]*tpcXnro[0]);
  float q22SLev = (Xl2*Xr2+Yl2*Yr2)/(Xl0*Xr0);
  float q22SLPev = (tpcXpl[2]*Xr2+tpcYpl[2]*Yr2)/(tpcXpl[0]*Xr0);
  float q22SLNev = (tpcXnl[2]*Xr2+tpcYnl[2]*Yr2)/(tpcXnl[0]*Xr0);
  h_q22_cent->Fill(cent,q22ev);
  h_q22_centP->Fill(cent,q22Pev);
  h_q22_centN->Fill(cent,q22Nev);
  h_q22_centPP->Fill(cent,q22PPev);
  h_q22_centNN->Fill(cent,q22NNev);
  h_q22gap0_cent->Fill(cent,q22gap0ev);
  h_q22gap0_centP->Fill(cent,q22gap0Pev);
  h_q22gap0_centN->Fill(cent,q22gap0Nev);
  h_q22gap1_cent->Fill(cent,q22gap1ev);
  h_q22gap1_centP->Fill(cent,q22gap1Pev);
  h_q22gap1_centN->Fill(cent,q22gap1Nev);
  
  float q32ev = (tpcQQ[3]-M)/W_2;
  float q32Pev = (tpcQQp[3]-Mp)/(Mp*M-Mp);
  float q32Nev = (tpcQQn[3]-Mn)/(Mn*M-Mn);
  float q32PPev = (tpcQQpp[3]-Mp)/Wp_2;
  float q32NNev = (tpcQQnn[3]-Mn)/Wn_2;
  float q32gap0ev = (Xl3*Xr3+Yl3*Yr3)/(Xl0*Xr0);
  float q32gap0Pev = (tpcXpl[3]*tpcXpr[3]+tpcYpl[3]*tpcYpr[3])/(tpcXpl[0]*tpcXpr[0]);
  float q32gap0Nev = (tpcXnl[3]*tpcXnr[3]+tpcYnl[3]*tpcYnr[3])/(tpcXnl[0]*tpcXnr[0]);
  float q32gap1ev = (Xl3o*Xr3o+Yl3o*Yr3o)/(Xl0o*Xr0o);
  float q32gap1Pev = (tpcXplo[3]*tpcXpro[3]+tpcYplo[3]*tpcYpro[3])/(tpcXplo[0]*tpcXpro[0]);
  float q32gap1Nev = (tpcXnlo[3]*tpcXnro[3]+tpcYnlo[3]*tpcYnro[3])/(tpcXnlo[0]*tpcXnro[0]);
  float q32SLev = (Xl3*Xr3+Yl3*Yr3)/(Xl0*Xr0);
  float q32SLPev = (tpcXpl[3]*Xr3+tpcYpl[3]*Yr3)/(tpcXpl[0]*Xr0);
  float q32SLNev = (tpcXnl[3]*Xr3+tpcYnl[3]*Yr3)/(tpcXnl[0]*Xr0);

  float q42ev = (tpcQQ[4]-M)/W_2;
  float q42Pev = (tpcQQp[4]-Mp)/(Mp*M-Mp);
  float q42Nev = (tpcQQn[4]-Mn)/(Mn*M-Mn);
  float q42PPev = (tpcQQp[4]-Mp)/Wp_2;
  float q42NNev = (tpcQQn[4]-Mn)/Wn_2;
  float q42gap0ev = (Xl4*Xr4+Yl4*Yr4)/(Xl0*Xr0);
  float q42gap0Pev = (tpcXpl[4]*tpcXpr[4]+tpcYpl[4]*tpcYpr[4])/(tpcXpl[0]*tpcXpr[0]);
  float q42gap0Nev = (tpcXnl[4]*tpcXnr[4]+tpcYnl[4]*tpcYnr[4])/(tpcXnl[0]*tpcXnr[0]);
  float q42gap1ev = (Xl4o*Xr4o+Yl4o*Yr4o)/(Xl0o*Xr0o);
  float q42gap1Pev = (tpcXplo[4]*tpcXpro[4]+tpcYplo[4]*tpcYpro[4])/(tpcXplo[0]*tpcXpro[0]);
  float q42gap1Nev = (tpcXnlo[4]*tpcXnro[4]+tpcYnlo[4]*tpcYnro[4])/(tpcXnlo[0]*tpcXnro[0]);
  float q42SLev = (Xl4*Xr4+Yl4*Yr4)/(Xl0*Xr0);
  float q42SLPev = (tpcXpl[4]*Xr4+tpcYpl[4]*Yr4)/(tpcXpl[0]*Xr0);
  float q42SLNev = (tpcXnl[4]*Xr4+tpcYnl[4]*Yr4)/(tpcXnl[0]*Xr0);


  // 3p correlator  
  float q23ev = calc3(tpcX[1],tpcY[1],tpcX[2],tpcY[2],M);
  float q23Pev = calc3(tpcXp[1],tpcYp[1],tpcXp[2],tpcYp[2],Mp);
  float q23Nev = calc3(tpcXn[1],tpcYn[1],tpcXn[2],tpcYn[2],Mn);
  // 4p correlator
  float q24ev = calc4(tpcX[2],tpcY[2],tpcX[4],tpcY[4],M);
  float q24Pev = calc4(tpcXp[2],tpcYp[2],tpcXp[4],tpcYp[4],Mp);
  float q24Nev = calc4(tpcXn[2],tpcYn[2],tpcXn[4],tpcYn[4],Mn);



  h_q22_cent->Fill(cent,q22ev);
  h_q22_centP->Fill(cent,q22Pev);
  h_q22_centN->Fill(cent,q22Nev);
  h_q22_centPP->Fill(cent,q22PPev);
  h_q22_centNN->Fill(cent,q22NNev);
  h_q23_cent->Fill(cent,q23ev);
  h_q23_centP->Fill(cent,q23Pev);
  h_q23_centN->Fill(cent,q23Nev);
  h_q24_cent->Fill(cent,q24ev);
  h_q24_centP->Fill(cent,q24Pev);
  h_q24_centN->Fill(cent,q24Nev);
  h_q22gap0_cent->Fill(cent,q22gap0ev);
  h_q22gap0_centP->Fill(cent,q22gap0Pev);
  h_q22gap0_centN->Fill(cent,q22gap0Nev);
  h_q22gap1_cent->Fill(cent,q22gap1ev);
  h_q22gap1_centP->Fill(cent,q22gap1Pev);
  h_q22gap1_centN->Fill(cent,q22gap1Nev);

  h_q32_cent->Fill(cent,q32ev);
  h_q32_centP->Fill(cent,q32Pev);
  h_q32_centN->Fill(cent,q32Nev);
  h_q32gap0_cent->Fill(cent,q32gap0ev);
  h_q32gap0_centP->Fill(cent,q32gap0Pev);
  h_q32gap0_centN->Fill(cent,q32gap0Nev);
  h_q32gap1_cent->Fill(cent,q32gap1ev);
  h_q32gap1_centP->Fill(cent,q32gap1Pev);
  h_q32gap1_centN->Fill(cent,q32gap1Nev);

  h_q42_cent->Fill(cent,q42ev);
  h_q42_centP->Fill(cent,q42Pev);
  h_q42_centN->Fill(cent,q42Nev);
  h_q42gap0_cent->Fill(cent,q42gap0ev);
  h_q42gap0_centP->Fill(cent,q42gap0Pev);
  h_q42gap0_centN->Fill(cent,q42gap0Nev);
  h_q42gap1_cent->Fill(cent,q42gap1ev);
  h_q42gap1_centP->Fill(cent,q42gap1Pev);
  h_q42gap1_centN->Fill(cent,q42gap1Nev);

  if(fabs(qasymm)<10.0&&fabs(qasymmL)<10.0&&fabs(qasymmR)<10.0)
    {
      h_q22qasym_cent[icent]->Fill(qasymm,q22ev);
      h_q22qasymP_cent[icent]->Fill(qasymm,q22Pev);
      h_q22qasymN_cent[icent]->Fill(qasymm,q22Nev);
      h_q22qasym_gap0_cent[icent]->Fill(qasymm,q22gap0ev);
      h_q22qasymP_gap0_cent[icent]->Fill(qasymm,q22gap0Pev);
      h_q22qasymN_gap0_cent[icent]->Fill(qasymm,q22gap0Nev);
      h_q22qasym_gap1_cent[icent]->Fill(qasymm,q22gap1ev);
      h_q22qasymP_gap1_cent[icent]->Fill(qasymm,q22gap1Pev);
      h_q22qasymN_gap1_cent[icent]->Fill(qasymm,q22gap1Nev);
      
      h_q32qasym_cent[icent]->Fill(qasymm,q32ev);
      h_q32qasymP_cent[icent]->Fill(qasymm,q32Pev);
      h_q32qasymN_cent[icent]->Fill(qasymm,q32Nev);
      h_q32qasym_gap0_cent[icent]->Fill(qasymm,q32gap0ev);
      h_q32qasymP_gap0_cent[icent]->Fill(qasymm,q32gap0Pev);
      h_q32qasymN_gap0_cent[icent]->Fill(qasymm,q32gap0Nev);
      h_q32qasym_gap1_cent[icent]->Fill(qasymm,q32gap1ev);
      h_q32qasymP_gap1_cent[icent]->Fill(qasymm,q32gap1Pev);
      h_q32qasymN_gap1_cent[icent]->Fill(qasymm,q32gap1Nev);

      h_q42qasym_cent[icent]->Fill(qasymm,q42ev);
      h_q42qasymP_cent[icent]->Fill(qasymm,q42Pev);
      h_q42qasymN_cent[icent]->Fill(qasymm,q42Nev);
      h_q42qasym_gap0_cent[icent]->Fill(qasymm,q42gap0ev);
      h_q42qasymP_gap0_cent[icent]->Fill(qasymm,q42gap0Pev);
      h_q42qasymN_gap0_cent[icent]->Fill(qasymm,q42gap0Nev);
      h_q42qasym_gap1_cent[icent]->Fill(qasymm,q42gap1ev);
      h_q42qasymP_gap1_cent[icent]->Fill(qasymm,q42gap1Pev);
      h_q42qasymN_gap1_cent[icent]->Fill(qasymm,q42gap1Nev);
    }

  h_Aq22_cent->Fill(cent,q22ev*qasymm);
  h_Aq22_centP->Fill(cent,q22Pev*qasymm);
  h_Aq22_centN->Fill(cent,q22Nev*qasymm);
  h_Aq22_centPP->Fill(cent,q22PPev*qasymm);
  h_Aq22_centNN->Fill(cent,q22NNev*qasymm);
  h_Aq22_SL_cent->Fill(cent,q22SLev*qasymm);
  h_Aq22_SL_centP->Fill(cent,q22SLPev*qasymm);
  h_Aq22_SL_centN->Fill(cent,q22SLNev*qasymm);
  h_Aq22_SLL_cent->Fill(cent,q22SLev*qasymmL);
  h_Aq22_SLL_centP->Fill(cent,q22SLPev*qasymmL);
  h_Aq22_SLL_centN->Fill(cent,q22SLNev*qasymmL);
  h_Aq22_SLR_cent->Fill(cent,q22SLev*qasymmR);
  h_Aq22_SLR_centP->Fill(cent,q22SLPev*qasymmR);
  h_Aq22_SLR_centN->Fill(cent,q22SLNev*qasymmR);
  
  h_Aq32_cent->Fill(cent,q32ev*qasymm);
  h_Aq32_centP->Fill(cent,q32Pev*qasymm);
  h_Aq32_centN->Fill(cent,q32Nev*qasymm);
  h_Aq32_centPP->Fill(cent,q32PPev*qasymm);
  h_Aq32_centNN->Fill(cent,q32NNev*qasymm);
  h_Aq32_SL_cent->Fill(cent,q32SLev*qasymm);
  h_Aq32_SL_centP->Fill(cent,q32SLPev*qasymm);
  h_Aq32_SL_centN->Fill(cent,q32SLNev*qasymm);
  h_Aq32_SLL_cent->Fill(cent,q32SLev*qasymmL);
  h_Aq32_SLL_centP->Fill(cent,q32SLPev*qasymmL);
  h_Aq32_SLL_centN->Fill(cent,q32SLNev*qasymmL);
  h_Aq32_SLR_cent->Fill(cent,q32SLev*qasymmR);
  h_Aq32_SLR_centP->Fill(cent,q32SLPev*qasymmR);
  h_Aq32_SLR_centN->Fill(cent,q32SLNev*qasymmR);
  
  h_Aq42_cent->Fill(cent,q42ev*qasymm);
  h_Aq42_centP->Fill(cent,q42Pev*qasymm);
  h_Aq42_centN->Fill(cent,q42Nev*qasymm);
  h_Aq42_centPP->Fill(cent,q42PPev*qasymm);
  h_Aq42_centNN->Fill(cent,q42NNev*qasymm);
  h_Aq42_SL_cent->Fill(cent,q42SLev*qasymm);
  h_Aq42_SL_centP->Fill(cent,q42SLPev*qasymm);
  h_Aq42_SL_centN->Fill(cent,q42SLNev*qasymm);
  h_Aq42_SLL_cent->Fill(cent,q42SLev*qasymmL);
  h_Aq42_SLL_centP->Fill(cent,q42SLPev*qasymmL);
  h_Aq42_SLL_centN->Fill(cent,q42SLNev*qasymmL);
  h_Aq42_SLR_cent->Fill(cent,q42SLev*qasymmR);
  h_Aq42_SLR_centP->Fill(cent,q42SLPev*qasymmR);
  h_Aq42_SLR_centN->Fill(cent,q42SLNev*qasymmR);
  


  float a1Xp[hmax], a1Yp[hmax], a2Xp[hmax], a2Yp[hmax];
  float a1Xn[hmax], a1Yn[hmax], a2Xn[hmax], a2Yn[hmax];
  
  for(int i=0; i<hmax;i++)
    {
      a1Xp[i] = 0.0;
      a1Yp[i] = 0.0;
      a1Xn[i] = 0.0;
      a1Yn[i] = 0.0;
      a2Xn[i] = 0.0;
      a2Yn[i] = 0.0;
      a2Xp[i] = 0.0;
      a2Yp[i] = 0.0;
    }
  
  // --- now second track loop

  for(int itrk = 0; itrk<d_ntrk; itrk++)
    {
      AliAODTrack *track = (AliAODTrack *)fAOD->GetTrack(itrk);
      if(!track)
	{
	  if(debug>0) cout<<"ERROR: Could not retrieve AODtrack "<<itrk<<endl;
	  continue;
	}

      if(!track->TestFilterBit(fbit))
	{
	  if(debug>15) cout<<"second track loop: wrong filter bit for track "<<itrk<<endl;
	  continue;
	}
      
      // int gid = track->GetID();
      // AliAODTrack *PIDtrack;
      // if(gid>=0) PIDtrack = track;
      // else PIDtrack = fAOD->GetTrack(trackMap->GetValue(-1-gid));

      // --- get track variables      
      float pt = track->Pt();
      float phi = track->Phi();
      float eta = track->Eta();
      int charge = track->Charge();
      int ncls = track->GetTPCNcls();
      //float dedx = track->GetTPCsignal();

      float Tdcaxy = track->DCA();
      float Tdcaz = track->ZAtDCA();


      float dcaxy = -999;
      float dcax = -999;
      float dcay = -999;
      float dcaz = -999;

      double r[3];
      bool dcaflag = track->GetXYZ(r);

      double DCA[2]; // dca
      double COV[3]; // covariance
      bool proptodca = track->PropagateToDCA(fVtx,mag,100.0,DCA,COV);
      if(!proptodca&&debug>15) cout<<"Second loop no DCACOV for you!"<<endl;

      if(dcaflag)
      	{
      	  dcaxy = r[0];
      	  dcaz = r[1];
      	  //if(debug>6) cout<<"GetXYZ is true, filter bit is "<<fbit<<endl;
      	}
      else
      	{
      	  dcax = r[0] - eventX;
      	  dcay = r[1] - eventY;
      	  dcaz = r[2] - eventZ;
      	  // --- need the sign convention for dcaxy...
      	  dcaxy = sqrt(dcax*dcax+dcay*dcay);
      	  //if((float)dcaxy!=(float)fabs((float)DCA[0])&&debug>6) cout<<"hmm... "<<dcaxy<<" "<<DCA[0]<<endl;
      	  // --- set dcaxy to value from PropagateToDCA to get correct sign
      	  dcaxy = DCA[0];
      	  // --- dcaz on the other hand is unambiguous
      	  //if(dcaz!=(float)DCA[1]&&debug>6) cout<<"hmm... "<<dcaz<<" "<<DCA[1]<<endl;
      	  //if(debug>6) cout<<"GetXYZ is false, filter bit is "<<fbit<<endl;
      	}

      // if(debug>6)
      // 	{
      // 	  cout<<"r[0] is "<<r[0]<<" ";
      // 	  cout<<"r[1] is "<<r[1]<<" ";
      // 	  cout<<"r[2] is "<<r[2]<<endl;
      // 	  cout<<"eventX is "<<eventX<<" ";
      // 	  cout<<"eventY is "<<eventY<<" ";
      // 	  cout<<"eventZ is "<<eventZ<<endl;
      // 	  cout<<"Tdcaxy is "<<Tdcaxy<<" and dcaxy is "<<dcaxy<<" and DCACOV xy is "<<DCA[0]<<endl;
      // 	  cout<<"Tdcaz is "<<Tdcaz<<" and dcaz is "<<dcaz<<" and DCACOV z is "<<DCA[1]<<endl;
      // 	}

      if((fbit==128||fbit==272)&&!dcaflag)
      	{
      	  dcaxy = Tdcaxy;
      	  dcaz = Tdcaz;
      	  // if(debug>4)
      	  //   {
      	  //     cout<<"GetXYZ false but should be true, flipping values"<<endl;
      	  //     cout<<"Tdcaxy is "<<Tdcaxy<<" and dcaxy is "<<dcaxy<<" and DCACOV xy is "<<DCA[0]<<endl;
      	  //     cout<<"Tdcaz is "<<Tdcaz<<" and dcaz is "<<dcaz<<" and DCACOV z is "<<DCA[1]<<endl;
      	  //   }
      	}

      // --- track cuts
      if(ncls<nclscut) continue;
      if(pt<ptmin||pt>ptmax) continue;
      if(fabs(eta)>outeta||fabs(eta)<excleta) continue;
      if(dodcacuts&&(fabs(dcaxy)>dcacutxy||fabs(dcaz)>dcacutz)) continue;// problems depending on filter bit

      bool pos = charge>0.0;
      bool neg = charge<0.0;


      if(gRandom) delete gRandom;
      gRandom = new TRandom3(0);
      gRandom->SetSeed(0);
      double rand = gRandom->Rndm();
      bool sub1 = false;
      bool sub2 = false;
      if(rand<0.5) sub1 = true;
      else sub2 = true;
      
      if(sub1)
	{
	  ntrkA1++;	
	  if(charge>0) ntrkposA1++;
	  if(charge<0) ntrknegA1++;
	}
      if(sub2)
	{
	  ntrkA2++;	
	  if(charge>0) ntrkposA2++;
	  if(charge<0) ntrknegA2++;
	}

      if(sub1) // random selection
	{
	  for(int i=1; i<hmax;i++)
	    {
	      if(charge>0)
		{
		  a1Xp[i] += cos(i*phi);
		  a1Yp[i] += sin(i*phi);
		}
	      if(charge<0)
		{
		  a1Xn[i] += cos(i*phi);
		  a1Yn[i] += sin(i*phi);
		}
	    }
	  // assign values to array element 0
	  if(charge>0)
	    {
	      a1Xp[0] += 1.0;
	      a1Yp[0] += pt;
	    }	
	  if(charge<0)
	    {
	      a1Xn[0] += 1.0;
	      a1Yn[0] += pt;
	    }
	}
      
      if(sub2) // random selection
	{
	  for(int i=1; i<hmax;i++)
	    {
	      if(charge>0)
		{
		  a2Xp[i] += cos(i*phi);
		  a2Yp[i] += sin(i*phi);
		}
	      if(charge<0)
		{
		  a2Xn[i] += cos(i*phi);
		  a2Yn[i] += sin(i*phi);
		}
	    }
	  // assign values to array element 0
	  if(charge>0)
	    {
	      a2Xp[0] += 1.0;
	      a2Yp[0] += pt;
	    }	
	  if(charge<0)
	    {
	      a2Xn[0] += 1.0;
	      a2Yn[0] += pt;
	    }
	}
      


      // -------------------------------------------------------
      // --- now calculate differential cumulants for each track
      // -------------------------------------------------------	      
      float trkXpl[hmax], trkYpl[hmax], trkXpr[hmax], trkYpr[hmax];
      float trkXnl[hmax], trkYnl[hmax], trkXnr[hmax], trkYnr[hmax];
      for(int i=0; i<hmax;i++)
	{
	  trkXpl[i] = 0.0;
	  trkYpl[i] = 0.0;
	  trkXpr[i] = 0.0;
	  trkYpr[i] = 0.0;
	  trkXnl[i] = 0.0;
	  trkYnl[i] = 0.0;
	  trkXnr[i] = 0.0;
	  trkYnr[i] = 0.0;
	}
      // very similar to above, except += changed to =
      if(pt>ptmin&&pt<ptmax)
	{
	  if(eta>-outeta&&eta<-excleta)
	    {
	      // assign values to array elements 1..8
	      for(int i=1; i<hmax;i++)
		{
		  if(charge>0)
		    {
		      trkXpl[i] = cos(i*phi);
		      trkYpl[i] = sin(i*phi);
		    }
		  if(charge<0)
		    {
		      trkXnl[i] = cos(i*phi);
		      trkYnl[i] = sin(i*phi);
		    }
		}
	      // assign values to array element 0
	      if(charge>0)
		{
		  trkXpl[0] = 1.0;
		  trkYpl[0] = pt;
		}	
	      if(charge<0)
		{
		  trkXnl[0] = 1.0;
		  trkYnl[0] = pt;
		}
	    }
	  //----------------------------
	  if(eta>excleta&&eta<outeta) // right half of TPC
	    {
	      for(int i=1; i<hmax; i++)
		{
		  if(charge>0)
		    {
		      trkXpr[i] = cos(i*phi);
		      trkYpr[i] = sin(i*phi);
		    }
		  if(charge<0)
		    {
		      trkXnr[i] = cos(i*phi);
		      trkYnr[i] = sin(i*phi);
		    }
		}
	      if(charge>0)
		{
		  trkXpr[0] = 1.0;
		  trkYpr[0] = pt;
		}	
	      if(charge<0)
		{
		  trkXnr[0] = 1.0;
		  trkYnr[0] = pt;
		}
	    } // end right half of TPC
	} // end pT selection
      
      // COME BACK HERE
      float trkX[9],  trkY[9];
      float trkXp[9], trkYp[9]; // pos
      float trkXn[9], trkYn[9]; // neg
      for(int i=0; i<6; i++)
	{
	  trkX[i]=trkXpl[i]+trkXnl[i]+trkXpr[i]+trkXnr[i];
	  trkY[i]=trkYpl[i]+trkYnl[i]+trkYpr[i]+trkYnr[i];
	  // pos
	  trkXp[i]=trkXpl[i]+trkXpr[i];
	  trkYp[i]=trkYpl[i]+trkYpr[i];
	  // neg
	  trkXn[i]=trkXnl[i]+trkXnr[i];
	  trkYn[i]=trkYnl[i]+trkYnr[i];
	}
      
      // sanity checks to prevent division by zero and other problems
      if(trkX[0]<1) continue;
      if(tpcX[0]<2) continue;
      if(tpcXp[0]<2) continue;
      if(tpcXn[0]<2) continue;
      
      // --- calculate differential cumulants from Q-vector components
      // --- reference flow is from whole TPC
      float diffq22ev = ((trkX[2]*tpcX[2])+(trkY[2]*tpcY[2])-1)/(tpcX[0]-1);
      float diffq22Pev = ((trkXp[2]*tpcX[2])+(trkYp[2]*tpcY[2])-1)/(tpcX[0]-1);
      float diffq22Nev = ((trkXn[2]*tpcX[2])+(trkYn[2]*tpcY[2])-1)/(tpcX[0]-1);
      // --- third harmonic
      float diffq32ev = ((trkX[3]*tpcX[3])+(trkY[3]*tpcY[3])-1)/(tpcX[0]-1);
      float diffq32Pev = ((trkXp[3]*tpcX[3])+(trkYp[3]*tpcY[3])-1)/(tpcX[0]-1);
      float diffq32Nev = ((trkXn[3]*tpcX[3])+(trkYn[3]*tpcY[3])-1)/(tpcX[0]-1);
      // --- fourth harmonic
      float diffq42ev = ((trkX[4]*tpcX[4])+(trkY[4]*tpcY[4])-1)/(tpcX[0]-1);
      float diffq42Pev = ((trkXp[4]*tpcX[4])+(trkYp[4]*tpcY[4])-1)/(tpcX[0]-1);
      float diffq42Nev = ((trkXn[4]*tpcX[4])+(trkYn[4]*tpcY[4])-1)/(tpcX[0]-1);
      
      // --- differential cumulants test with differing RPs
      float Xdiffq22Pev = ((trkXp[2]*tpcX[2])+(trkYp[2]*tpcY[2])-1)/(tpcX[0]-1);
      float Xdiffq22Nev = ((trkXn[2]*tpcX[2])+(trkYn[2]*tpcY[2])-1)/(tpcX[0]-1);
      float Ydiffq22Pev = ((trkXp[2]*tpcXp[2])+(trkYp[2]*tpcYp[2])-1)/(tpcXp[0]-1);
      float Ydiffq22Nev = ((trkXn[2]*tpcXn[2])+(trkYn[2]*tpcYn[2])-1)/(tpcXn[0]-1);
      float Zdiffq22Pev = ((trkXp[2]*tpcXn[2])+(trkYp[2]*tpcYn[2]))/(tpcXn[0]);
      float Zdiffq22Nev = ((trkXn[2]*tpcXp[2])+(trkYn[2]*tpcYp[2]))/(tpcXp[0]);

      float Xdiffq32Pev = ((trkXp[3]*tpcX[3])+(trkYp[3]*tpcY[3])-1)/(tpcX[0]-1);
      float Xdiffq32Nev = ((trkXn[3]*tpcX[3])+(trkYn[3]*tpcY[3])-1)/(tpcX[0]-1);
      float Ydiffq32Pev = ((trkXp[3]*tpcXp[3])+(trkYp[3]*tpcYp[3])-1)/(tpcXp[0]-1);
      float Ydiffq32Nev = ((trkXn[3]*tpcXn[3])+(trkYn[3]*tpcYn[3])-1)/(tpcXn[0]-1);
      float Zdiffq32Pev = ((trkXp[3]*tpcXn[3])+(trkYp[3]*tpcYn[3]))/(tpcXn[0]);
      float Zdiffq32Nev = ((trkXn[3]*tpcXp[3])+(trkYn[3]*tpcYp[3]))/(tpcXp[0]);

      float Xdiffq42Pev = ((trkXp[4]*tpcX[4])+(trkYp[4]*tpcY[4])-1)/(tpcX[0]-1);
      float Xdiffq42Nev = ((trkXn[4]*tpcX[4])+(trkYn[4]*tpcY[4])-1)/(tpcX[0]-1);
      float Ydiffq42Pev = ((trkXp[4]*tpcXp[4])+(trkYp[4]*tpcYp[4])-1)/(tpcXp[0]-1);
      float Ydiffq42Nev = ((trkXn[4]*tpcXn[4])+(trkYn[4]*tpcYn[4])-1)/(tpcXn[0]-1);
      float Zdiffq42Pev = ((trkXp[4]*tpcXn[4])+(trkYp[4]*tpcYn[4]))/(tpcXn[0]);
      float Zdiffq42Nev = ((trkXn[4]*tpcXp[4])+(trkYn[4]*tpcYp[4]))/(tpcXp[0]);

      // --- histograms to check different reference flow
      // if(pos) h_Xdiffq22_ptP->Fill(pt,Xdiffq22Pev);
      // if(neg) h_Xdiffq22_ptN->Fill(pt,Xdiffq22Nev);
      // if(pos) h_Ydiffq22_ptP->Fill(pt,Ydiffq22Pev);
      // if(neg) h_Ydiffq22_ptN->Fill(pt,Ydiffq22Nev);
      // if(pos) h_Zdiffq22_ptP->Fill(pt,Zdiffq22Pev);
      // if(neg) h_Zdiffq22_ptN->Fill(pt,Zdiffq22Nev);
	  
      if(cent>=centlo&&cent<centhi)
	{
	  // --- differential cumulant		  
	  h_diffq22_pt->Fill(pt,diffq22ev);
	  h_diffq22_eta->Fill(eta,diffq22ev);
	  // --- differential cumulant weighted with qasymm
	  h_diffAq22_pt->Fill(pt,diffq22ev*qasymm);
	  h_diffAq22_eta->Fill(eta,diffq22ev*qasymm);
	  if(pos)
	    {
	      h_diffq22_ptP->Fill(pt,diffq22Pev);
	      h_diffq22_etaP->Fill(eta,diffq22Pev);
	      h_diffAq22_ptP->Fill(pt,diffq22Pev*qasymm);
	      h_diffAq22_etaP->Fill(eta,diffq22Pev*qasymm);
	    }
	  if(neg)
	    {
	      h_diffq22_ptN->Fill(pt,diffq22Nev);
	      h_diffq22_etaN->Fill(eta,diffq22Nev);
	      h_diffAq22_ptN->Fill(pt,diffq22Nev*qasymm);
	      h_diffAq22_etaN->Fill(eta,diffq22Nev*qasymm);
	    }
	  // --- same for third harmonic
	  h_diffq32_pt->Fill(pt,diffq32ev);
	  h_diffq32_eta->Fill(eta,diffq32ev);
	  h_diffAq32_pt->Fill(pt,diffq32ev*qasymm);
	  h_diffAq32_eta->Fill(eta,diffq32ev*qasymm);
	  if(pos)
	    {
	      h_diffq32_ptP->Fill(pt,diffq32Pev);
	      h_diffq32_etaP->Fill(eta,diffq32Pev);
	      h_diffAq32_ptP->Fill(pt,diffq32Pev*qasymm);
	      h_diffAq32_etaP->Fill(eta,diffq32Pev*qasymm);
	    }
	  if(neg)
	    {
	      h_diffq32_ptN->Fill(pt,diffq32Nev);
	      h_diffq32_etaN->Fill(eta,diffq32Nev);
	      h_diffAq32_ptN->Fill(pt,diffq32Nev*qasymm);
	      h_diffAq32_etaN->Fill(eta,diffq32Nev*qasymm);
	    }
	  // --- same for fourth harmonic
	  h_diffq42_pt->Fill(pt,diffq42ev);
	  h_diffq42_eta->Fill(eta,diffq42ev);
	  h_diffAq42_pt->Fill(pt,diffq42ev*qasymm);
	  h_diffAq42_eta->Fill(eta,diffq42ev*qasymm);
	  if(pos)
	    {
	      h_diffq42_ptP->Fill(pt,diffq42Pev);
	      h_diffq42_etaP->Fill(eta,diffq42Pev);
	      h_diffAq42_ptP->Fill(pt,diffq42Pev*qasymm);
	      h_diffAq42_etaP->Fill(eta,diffq42Pev*qasymm);
	    }
	  if(neg)
	    {
	      h_diffq42_ptN->Fill(pt,diffq42Nev);
	      h_diffq42_etaN->Fill(eta,diffq42Nev);
	      h_diffAq42_ptN->Fill(pt,diffq42Nev*qasymm);
	      h_diffAq42_etaN->Fill(eta,diffq42Nev*qasymm);
	    }

	  // selection to do nested track loop
	  if(donested)
	    {
	      for(int i=0; i<d_ntrk; i++)
		{
		  // ----------------------------------------------------------------------------------------
		  // --- want delta eta and delta pt dependence of differential cumulant weighted with charge
		  // --- need nested track loop to get eta1 and eta3
		  // ----------------------------------------------------------------------------------------
		  // make sure eta is in range, -99 should be autorejected by this
		  //if(eta3[i]==-99) continue;
		  if(eta3[i]<-outeta||eta3[i]>outeta)
		    {
		      //cout<<"eta3 out of range!!! "<<eta3[i]<<endl;
		      continue; 
		    }
		  if(eta<-outeta||eta>outeta)
		    {
		      cout<<"eta out of range!!! "<<eta<<endl;
		      continue; 
		    }
		  // charge==0 should be rejected by the above cut
		  if(charge3[i]==0) {cout<<"WTF!"<<endl;continue;}
		  // VERY IMPORTANT remove auto-correlations with 3rd particle
		  if(i==itrk) continue;
		  
		  float DETA = eta-eta3[i];
		  float DPHI = phi-phi3[i];
		  float DPT = pt-pt3[i];

		  // magnetic field sign
		  //bSign = (eventMain->GetMagneticField() > 0) ? 1 : -1;
		  float bSign = 1;
		  if(mag<0) bSign = -1;

		  // double deta = firstEta - secondEta[j];
		  // double dphi = firstPhi - secondPhi[j];
		  float deta = DETA;
		  float dphi = DPHI;
		  float fHBTCutValue = 0.02;
		  // optimization
		  if(dopaircut&&(fabs(deta)<fHBTCutValue*2.5*3)) //fHBTCutValue = 0.02 [default for dphicorrelations]
		    {
		      if(debug>4)
		  	{
		  	  cout<<"inside deta cut, evaluating"<<endl;
		  	  cout<<"deta is "<<deta<<" and dphi is "<<dphi<<endl;
		  	}
		      // phi in rad
		      // float phi1rad = firstPhi;
		      // float phi2rad = secondPhi[j];
   
		      // check first boundaries to see if is worth to loop and find the minimum
		      // float dphistar1 = GetDPhiStar(phi1rad, firstPt, charge1, phi2rad, secondPt[j], charge2, 0.8, bSign);
		      // float dphistar2 = GetDPhiStar(phi1rad, firstPt, charge1, phi2rad, secondPt[j], charge2, 2.5, bSign);
		      float dphistar1 = GetDPhiStar(phi, pt, charge, phi3[i], pt3[i], charge3[i], 0.8, bSign);
		      float dphistar2 = GetDPhiStar(phi, pt, charge, phi3[i], pt3[i], charge3[i], 2.5, bSign);
   
		      const float kLimit = fHBTCutValue * 3;
   
		      float dphistarminabs = 1e5;
		      //float dphistarmin = 1e5;
   
		      if(fabs(dphistar1) < kLimit || fabs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0 )
		  	{
		  	  if(debug>4)
		  	    {
		  	      cout<<"inside detadphi cut, evaluating"<<endl;
		  	      cout<<"deta is "<<deta<<" and dphi is "<<dphi<<endl;
		  	    }
		  	  for(double rad=0.8; rad<2.51; rad+=0.01)
		  	    {
		  	      //float dphistar = GetDPhiStar(phi1rad, firstPt, charge1, phi2rad, secondPt[j], charge2, rad, bSign);
		  	      float dphistar = GetDPhiStar(phi, pt, charge, phi3[i], pt3[i], charge3[i], 0.8, bSign);
		  	      float dphistarabs = fabs(dphistar);
		  	      if(dphistarabs < dphistarminabs)
		  		{
		  		  //dphistarmin = dphistar;
		  		  dphistarminabs = dphistarabs;
		  		}
		  	    }
		  	  if(dphistarminabs < fHBTCutValue && fabs(deta) < fHBTCutValue)
		  	    {
		  	      if(debug>4)
		  		printf("HBT: Removed track pair %d %d with [[%f %f]] %f %f %f | %f %f %d %f %f %d %f",
		  		       itrk, i, deta, dphi, dphistarminabs, dphistar1, dphistar2, phi, pt, charge, phi3[i], pt3[i], charge3[i], bSign);
		  	      continue;
		  	    }
		  	} // loop over angles
		    } // selection on pairs close together




		  
		  h_AT_X_deta->Fill(DETA,charge3[i]);
		  if(pos) h_AT_X_detaP->Fill(DETA,charge3[i]);
		  if(neg) h_AT_X_detaN->Fill(DETA,charge3[i]);
		  h_AT_X_dpt->Fill(DPT,charge3[i]);
		  if(pos) h_AT_X_dptP->Fill(DPT,charge3[i]);
		  if(neg) h_AT_X_dptN->Fill(DPT,charge3[i]);
		  // --- particles 1, 2, and 3 in whole TPC
		  if(pos) h_diffq22_X_detaP->Fill(DETA,diffq22Pev);
		  if(neg) h_diffq22_X_detaN->Fill(DETA,diffq22Nev);
		  if(pos) h_diffq32_X_detaP->Fill(DETA,diffq32Pev);
		  if(neg) h_diffq32_X_detaN->Fill(DETA,diffq32Nev);
		  if(pos) h_diffq42_X_detaP->Fill(DETA,diffq42Pev);
		  if(neg) h_diffq42_X_detaN->Fill(DETA,diffq42Nev);
		  if(pos) h_ATdiffq22_X_detaP->Fill(DETA,diffq22Pev*charge3[i]);
		  if(neg) h_ATdiffq22_X_detaN->Fill(DETA,diffq22Nev*charge3[i]);
		  if(pos) h_ATdiffq32_X_detaP->Fill(DETA,diffq32Pev*charge3[i]);
		  if(neg) h_ATdiffq32_X_detaN->Fill(DETA,diffq32Nev*charge3[i]);
		  if(pos) h_ATdiffq42_X_detaP->Fill(DETA,diffq42Pev*charge3[i]);
		  if(neg) h_ATdiffq42_X_detaN->Fill(DETA,diffq42Nev*charge3[i]);
		  if(pos) h_diffq22_X_dptP->Fill(DPT,diffq22Pev);
		  if(neg) h_diffq22_X_dptN->Fill(DPT,diffq22Nev);
		  if(pos) h_diffq32_X_dptP->Fill(DPT,diffq32Pev);
		  if(neg) h_diffq32_X_dptN->Fill(DPT,diffq32Nev);
		  if(pos) h_diffq42_X_dptP->Fill(DPT,diffq42Pev);
		  if(neg) h_diffq42_X_dptN->Fill(DPT,diffq42Nev);
		  if(pos) h_ATdiffq22_X_dptP->Fill(DPT,diffq22Pev*charge3[i]);
		  if(neg) h_ATdiffq22_X_dptN->Fill(DPT,diffq22Nev*charge3[i]);
		  if(pos) h_ATdiffq32_X_dptP->Fill(DPT,diffq32Pev*charge3[i]);
		  if(neg) h_ATdiffq32_X_dptN->Fill(DPT,diffq32Nev*charge3[i]);
		  if(pos) h_ATdiffq42_X_dptP->Fill(DPT,diffq42Pev*charge3[i]);
		  if(neg) h_ATdiffq42_X_dptN->Fill(DPT,diffq42Nev*charge3[i]);

		  if(pos) h_ATXdiffq22_X_detaP->Fill(DETA,Xdiffq22Pev*charge3[i]);
		  if(neg) h_ATXdiffq22_X_detaN->Fill(DETA,Xdiffq22Nev*charge3[i]);
		  if(pos) h_ATYdiffq22_X_detaP->Fill(DETA,Ydiffq22Pev*charge3[i]);
		  if(neg) h_ATYdiffq22_X_detaN->Fill(DETA,Ydiffq22Nev*charge3[i]);
		  if(pos) h_ATZdiffq22_X_detaP->Fill(DETA,Zdiffq22Pev*charge3[i]);
		  if(neg) h_ATZdiffq22_X_detaN->Fill(DETA,Zdiffq22Nev*charge3[i]);

		  if(pos) h_ATXdiffq32_X_detaP->Fill(DETA,Xdiffq32Pev*charge3[i]);
		  if(neg) h_ATXdiffq32_X_detaN->Fill(DETA,Xdiffq32Nev*charge3[i]);
		  if(pos) h_ATYdiffq32_X_detaP->Fill(DETA,Ydiffq32Pev*charge3[i]);
		  if(neg) h_ATYdiffq32_X_detaN->Fill(DETA,Ydiffq32Nev*charge3[i]);
		  if(pos) h_ATZdiffq32_X_detaP->Fill(DETA,Zdiffq32Pev*charge3[i]);
		  if(neg) h_ATZdiffq32_X_detaN->Fill(DETA,Zdiffq32Nev*charge3[i]);

		  if(pos) h_ATXdiffq42_X_detaP->Fill(DETA,Xdiffq42Pev*charge3[i]);
		  if(neg) h_ATXdiffq42_X_detaN->Fill(DETA,Xdiffq42Nev*charge3[i]);
		  if(pos) h_ATYdiffq42_X_detaP->Fill(DETA,Ydiffq42Pev*charge3[i]);
		  if(neg) h_ATYdiffq42_X_detaN->Fill(DETA,Ydiffq42Nev*charge3[i]);
		  if(pos) h_ATZdiffq42_X_detaP->Fill(DETA,Zdiffq42Pev*charge3[i]);
		  if(neg) h_ATZdiffq42_X_detaN->Fill(DETA,Zdiffq42Nev*charge3[i]);

		  
		  float meanchargeeta = 0;
		  int index = int((eta3[i]+1.0)*50);
		  if(index<0||index>100)
		    {
		      cout<<"out of index!!!"<<endl;
		      continue; // safety, should be redundant
		    }
		  //cout<<"eta is "<<eta3[i]<<" index is "<<index<<endl;
		  if(mag<-4.0)
		    {
		      if(fbit==1) meanchargeeta = MCEF1_fb1[index];
		      if(fbit==128) meanchargeeta = MCEF1_fb128[index];
		      if(fbit==272) meanchargeeta = MCEF1_fb272[index];
		      //if(fbit==1&&dodcacuts) meanchargeeta = MCEF1_dca[index];
		    }
		  if(mag>4.0)
		    {
		      if(fbit==1) meanchargeeta = MCEF3_fb1[index];
		      if(fbit==128) meanchargeeta = MCEF3_fb128[index];
		      if(fbit==272) meanchargeeta = MCEF3_fb272[index];
		      //if(fbit==1&&dodcacuts) meanchargeeta = MCEF3_dca[index];
		    }
		  float charge3subeta = charge3[i] - meanchargeeta;
		  
		  h_AT_S_deta->Fill(DETA,charge3subeta);
		  if(pos) h_AT_S_detaP->Fill(DETA,charge3subeta);
		  if(neg) h_AT_S_detaN->Fill(DETA,charge3subeta);
		  // --- particles 1, 2, and 3 in whole TPC
		  if(pos) h_diffq22_S_detaP->Fill(DETA,diffq22Pev);
		  if(neg) h_diffq22_S_detaN->Fill(DETA,diffq22Nev);
		  if(pos) h_diffq32_S_detaP->Fill(DETA,diffq32Pev);
		  if(neg) h_diffq32_S_detaN->Fill(DETA,diffq32Nev);
		  if(pos) h_diffq42_S_detaP->Fill(DETA,diffq42Pev);
		  if(neg) h_diffq42_S_detaN->Fill(DETA,diffq42Nev);
		  if(pos) h_ATdiffq22_S_detaP->Fill(DETA,diffq22Pev*charge3subeta);
		  if(neg) h_ATdiffq22_S_detaN->Fill(DETA,diffq22Nev*charge3subeta);
		  if(pos) h_ATdiffq32_S_detaP->Fill(DETA,diffq32Pev*charge3subeta);
		  if(neg) h_ATdiffq32_S_detaN->Fill(DETA,diffq32Nev*charge3subeta);
		  if(pos) h_ATdiffq42_S_detaP->Fill(DETA,diffq42Pev*charge3subeta);
		  if(neg) h_ATdiffq42_S_detaN->Fill(DETA,diffq42Nev*charge3subeta);

		  // --- now field selection for delta eta stuff
		  if(mag<-4.0)
		    {
		      h_AT_X_deta_F1->Fill(DETA,charge3[i]);
		      if(pos) h_AT_X_detaP_F1->Fill(DETA,charge3[i]);
		      if(neg) h_AT_X_detaN_F1->Fill(DETA,charge3[i]);
		      if(pos) h_diffq22_X_detaP_F1->Fill(DETA,diffq22Pev);
		      if(neg) h_diffq22_X_detaN_F1->Fill(DETA,diffq22Nev);
		      if(pos) h_diffq32_X_detaP_F1->Fill(DETA,diffq32Pev);
		      if(neg) h_diffq32_X_detaN_F1->Fill(DETA,diffq32Nev);
		      if(pos) h_diffq42_X_detaP_F1->Fill(DETA,diffq42Pev);
		      if(neg) h_diffq42_X_detaN_F1->Fill(DETA,diffq42Nev);
		      if(pos) h_ATdiffq22_X_detaP_F1->Fill(DETA,diffq22Pev*charge3[i]);
		      if(neg) h_ATdiffq22_X_detaN_F1->Fill(DETA,diffq22Nev*charge3[i]);
		      if(pos) h_ATdiffq32_X_detaP_F1->Fill(DETA,diffq32Pev*charge3[i]);
		      if(neg) h_ATdiffq32_X_detaN_F1->Fill(DETA,diffq32Nev*charge3[i]);
		      if(pos) h_ATdiffq42_X_detaP_F1->Fill(DETA,diffq42Pev*charge3[i]);
		      if(neg) h_ATdiffq42_X_detaN_F1->Fill(DETA,diffq42Nev*charge3[i]);
		    }
		  if(mag>4.0)
		    {
		      h_AT_X_deta_F3->Fill(DETA,charge3[i]);
		      if(pos) h_AT_X_detaP_F3->Fill(DETA,charge3[i]);
		      if(neg) h_AT_X_detaN_F3->Fill(DETA,charge3[i]);
		      if(pos) h_diffq22_X_detaP_F3->Fill(DETA,diffq22Pev);
		      if(neg) h_diffq22_X_detaN_F3->Fill(DETA,diffq22Nev);
		      if(pos) h_diffq32_X_detaP_F3->Fill(DETA,diffq32Pev);
		      if(neg) h_diffq32_X_detaN_F3->Fill(DETA,diffq32Nev);
		      if(pos) h_diffq42_X_detaP_F3->Fill(DETA,diffq42Pev);
		      if(neg) h_diffq42_X_detaN_F3->Fill(DETA,diffq42Nev);
		      if(pos) h_ATdiffq22_X_detaP_F3->Fill(DETA,diffq22Pev*charge3[i]);
		      if(neg) h_ATdiffq22_X_detaN_F3->Fill(DETA,diffq22Nev*charge3[i]);
		      if(pos) h_ATdiffq32_X_detaP_F3->Fill(DETA,diffq32Pev*charge3[i]);
		      if(neg) h_ATdiffq32_X_detaN_F3->Fill(DETA,diffq32Nev*charge3[i]);
		      if(pos) h_ATdiffq42_X_detaP_F3->Fill(DETA,diffq42Pev*charge3[i]);
		      if(neg) h_ATdiffq42_X_detaN_F3->Fill(DETA,diffq42Nev*charge3[i]);
		    }

		} // end nested track loop
	      
	    } // end if statement to check for nested track loop
	  
	} // centrality selection for differential cumulants
      
    } // end of second track loop

  float qasymA1 = (float)(ntrkposA1-ntrknegA1)/(ntrkA1);
  float qasymA2 = (float)(ntrkposA2-ntrknegA2)/(ntrkA2);

  h_rA_cent->Fill(cent,qasymA1);
  h_r1A_cent->Fill(cent,qasymA1);
  h_r2A_cent->Fill(cent,qasymA2);

  h2_rA_cent->Fill(cent,qasymA1);
  h2_r1A_cent->Fill(cent,qasymA1);
  h2_r2A_cent->Fill(cent,qasymA2);

  // random subevents 1
  float a1XX[9],  a1YY[9],  a1QQ[9];//,  a1qq[9];
  float a1XXp[9], a1YYp[9], a1QQp[9];//, a1qqp[9]; // pos
  float a1XXn[9], a1YYn[9], a1QQn[9];//, a1qqn[9]; // neg
  for(int i=0; i<6; i++)
    {
      a1XX[i]=a1Xp[i]+a1Xn[i];
      a1YY[i]=a1Yp[i]+a1Yn[i];
      a1QQ[i]=a1XX[i]*a1XX[i]+a1YY[i]*a1YY[i];
      //a1qq[i]=sqrt(a1QQ[i]/a1XX[0]);
      // pos
      a1XXp[i]=a1Xp[i];
      a1YYp[i]=a1Yp[i];
      a1QQp[i]=a1XXp[i]*a1XXp[i]+a1YYp[i]*a1YYp[i];
      //a1qqp[i]=sqrt(a1QQp[i]/a1XXp[0]);
      // neg
      a1XXn[i]=a1Xn[i];
      a1YYn[i]=a1Yn[i];
      a1QQn[i]=a1XXn[i]*a1XXn[i]+a1YYn[i]*a1YYn[i];
      //a1qqn[i]=sqrt(a1QQn[i]/a1XXn[0]);
    }
  // random subevents 2
  float a2XX[9],  a2YY[9],  a2QQ[9];//,  a2qq[9];
  float a2XXp[9], a2YYp[9], a2QQp[9];//, a2qqp[9]; // pos
  float a2XXn[9], a2YYn[9], a2QQn[9];//, a2qqn[9]; // neg
  for(int i=0; i<6; i++)
    {
      a2XX[i]=a2Xp[i]+a2Xn[i];
      a2YY[i]=a2Yp[i]+a2Yn[i];
      a2QQ[i]=a2XX[i]*a2XX[i]+a2YY[i]*a2YY[i];
      //a2qq[i]=sqrt(a2QQ[i]/a2XX[0]);
      // pos
      a2XXp[i]=a2Xp[i];
      a2YYp[i]=a2Yp[i];
      a2QQp[i]=a2XXp[i]*a2XXp[i]+a2YYp[i]*a2YYp[i];
      //a2qqp[i]=sqrt(a2QQp[i]/a2XXp[0]);
      // neg
      a2XXn[i]=a2Xn[i];
      a2YYn[i]=a2Yn[i];
      a2QQn[i]=a2XXn[i]*a2XXn[i]+a2YYn[i]*a2YYn[i];
      //a2qqn[i]=sqrt(a1QQn[i]/a1XXn[0]);
    }
  
  float a1M = a1XX[0];
  float a1W_2 = a1M*(a1M-(a1M/M));
  float a1Mp = a1XXp[0];
  float a1Wp_2 = a1Mp*(a1Mp-(a1Mp/Mp));
  float a1Mn = a1XXn[0];
  float a1Wn_2 = a1Mn*(a1Mn-(a1Mn/Mn));
  
  float a2M = a2XX[0];
  float a2W_2 = a2M*(a2M-(a2M/M));
  float a2Mp = a2XXp[0];
  float a2Wp_2 = a2Mp*(a2Mp-(a2Mp/Mp));
  float a2Mn = a2XXn[0];
  float a2Wn_2 = a2Mn*(a2Mn-(a2Mn/Mn));
  
  float a1q22ev = (a1QQ[2]-a1M)/a1W_2;
  float a1q22Pev = (a1QQp[2]-a1Mp)/a1Wp_2; // pos
  float a1q22Nev = (a1QQn[2]-a1Mn)/a1Wn_2; // neg
  
  float a2q22ev = (a2QQ[2]-a2M)/a2W_2;
  float a2q22Pev = (a2QQp[2]-a2Mp)/a2Wp_2; // pos
  float a2q22Nev = (a2QQn[2]-a2Mn)/a2Wn_2; // neg
  
  h_a1q22_cent->Fill(cent,a1q22ev);
  h_a1q22_centP->Fill(cent,a1q22Pev);
  h_a1q22_centN->Fill(cent,a1q22Nev);
  
  h_a2q22_cent->Fill(cent,a2q22ev);
  h_a2q22_centP->Fill(cent,a2q22Pev);
  h_a2q22_centN->Fill(cent,a2q22Nev);
  
  h_rAq22_X1_cent->Fill(cent,a1q22ev*qasymA1);
  h_rAq22_X1_centP->Fill(cent,a1q22Pev*qasymA1);
  h_rAq22_X1_centN->Fill(cent,a1q22Nev*qasymA1);
  h_rAq22_X2_cent->Fill(cent,a2q22ev*qasymA2);
  h_rAq22_X2_centP->Fill(cent,a2q22Pev*qasymA2);
  h_rAq22_X2_centN->Fill(cent,a2q22Nev*qasymA2);
  h_rAq22_X3_cent->Fill(cent,a1q22ev*qasymA2);
  h_rAq22_X3_centP->Fill(cent,a1q22Pev*qasymA2);
  h_rAq22_X3_centN->Fill(cent,a1q22Nev*qasymA2);
  h_rAq22_X4_cent->Fill(cent,a2q22ev*qasymA1);
  h_rAq22_X4_centP->Fill(cent,a2q22Pev*qasymA1);
  h_rAq22_X4_centN->Fill(cent,a2q22Nev*qasymA1);

  
  
  // ------------------------------------ //
  // --- send data to the output list --- //
  // ------------------------------------ //
  
  PostData(1,fOutputList);
  
  // ------------------------- //
  // --- end of event loop --- //
  // ------------------------- //

}      

void AliAnalysisTaskCMEv2A::SetParameters()
{

  cout<<"Now setting parameters!"<<endl;

  // --- these are the TPC efficiencies
  effTPC[0] = 0;//pt=0.05
  effTPC[1] = 0;//pt=0.15
  effTPC[2] = 0.895961;//pt=0.25
  effTPC[3] = 0.92504;//pt=0.35
  effTPC[4] = 0.944326;//pt=0.45
  effTPC[5] = 0.948261;//pt=0.55
  effTPC[6] = 0.948522;//pt=0.65
  effTPC[7] = 0.946941;//pt=0.75
  effTPC[8] = 0.948537;//pt=0.85
  effTPC[9] = 0.95194;//pt=0.95
  effTPC[10] = 0.95228;//pt=1.05
  effTPC[11] = 0.954775;//pt=1.15
  effTPC[12] = 0.956055;//pt=1.25
  effTPC[13] = 0.955758;//pt=1.35
  effTPC[14] = 0.951268;//pt=1.45
  effTPC[15] = 0.952081;//pt=1.55
  effTPC[16] = 0.947991;//pt=1.65
  effTPC[17] = 0.943755;//pt=1.75
  effTPC[18] = 0.938106;//pt=1.85
  effTPC[19] = 0.932583;//pt=1.95
  effTPC[20] = 0.927206;//pt=2.05
  effTPC[21] = 0.921111;//pt=2.15
  effTPC[22] = 0.913127;//pt=2.25
  effTPC[23] = 0.910508;//pt=2.35
  effTPC[24] = 0.905838;//pt=2.45
  effTPC[25] = 0.90243;//pt=2.55
  effTPC[26] = 0.897893;//pt=2.65
  effTPC[27] = 0.895459;//pt=2.75
  effTPC[28] = 0.89308;//pt=2.85
  effTPC[29] = 0.895013;//pt=2.95
  effTPC[30] = 0.890872;//pt=3.05
  effTPC[31] = 0.892046;//pt=3.15
  effTPC[32] = 0.890133;//pt=3.25
  effTPC[33] = 0.894786;//pt=3.35
  effTPC[34] = 0.891886;//pt=3.45
  effTPC[35] = 0.893649;//pt=3.55
  effTPC[36] = 0.891418;//pt=3.65
  effTPC[37] = 0.897408;//pt=3.75
  effTPC[38] = 0.890273;//pt=3.85
  effTPC[39] = 0.896526;//pt=3.95
  effTPC[40] = 0.892236;//pt=4.05
  effTPC[41] = 0.895075;//pt=4.15
  effTPC[42] = 0.895159;//pt=4.25
  effTPC[43] = 0.889562;//pt=4.35
  effTPC[44] = 0.890777;//pt=4.45
  effTPC[45] = 0.892925;//pt=4.55
  effTPC[46] = 0.890358;//pt=4.65
  effTPC[47] = 0.897014;//pt=4.75
  effTPC[48] = 0.886445;//pt=4.85
  effTPC[49] = 0.890928;//pt=4.95
  
  float max = effTPC[0];
  int maxbin = 0;
  for(int i=0; i<49; i++)
    {
      if(effTPC[i] > max)
	{
	  max = effTPC[i];
	  maxbin = i;
	}
    }
  cout<<"max for TPC is "<<max<<endl;
  cout<<"maxbin for TPC is "<<maxbin<<endl;
  for(int i=0; i<50; i++)
    {
      effTPC[i] /= max;
    }
 
  effHYB[0] = 0;//pt=0.05
  effHYB[1] = 0;//pt=0.15
  effHYB[2] = 0.797902;//pt=0.25
  effHYB[3] = 0.842964;//pt=0.35
  effHYB[4] = 0.872516;//pt=0.45
  effHYB[5] = 0.883225;//pt=0.55
  effHYB[6] = 0.886472;//pt=0.65
  effHYB[7] = 0.886011;//pt=0.75
  effHYB[8] = 0.883906;//pt=0.85
  effHYB[9] = 0.882904;//pt=0.95
  effHYB[10] = 0.881162;//pt=1.05
  effHYB[11] = 0.883611;//pt=1.15
  effHYB[12] = 0.884113;//pt=1.25
  effHYB[13] = 0.883676;//pt=1.35
  effHYB[14] = 0.875165;//pt=1.45
  effHYB[15] = 0.87379;//pt=1.55
  effHYB[16] = 0.871452;//pt=1.65
  effHYB[17] = 0.868833;//pt=1.75
  effHYB[18] = 0.863982;//pt=1.85
  effHYB[19] = 0.858165;//pt=1.95
  effHYB[20] = 0.850305;//pt=2.05
  effHYB[21] = 0.845429;//pt=2.15
  effHYB[22] = 0.839418;//pt=2.25
  effHYB[23] = 0.836466;//pt=2.35
  effHYB[24] = 0.833551;//pt=2.45
  effHYB[25] = 0.831225;//pt=2.55
  effHYB[26] = 0.826823;//pt=2.65
  effHYB[27] = 0.822582;//pt=2.75
  effHYB[28] = 0.82375;//pt=2.85
  effHYB[29] = 0.824866;//pt=2.95
  effHYB[30] = 0.822234;//pt=3.05
  effHYB[31] = 0.825215;//pt=3.15
  effHYB[32] = 0.82225;//pt=3.25
  effHYB[33] = 0.825158;//pt=3.35
  effHYB[34] = 0.822447;//pt=3.45
  effHYB[35] = 0.825239;//pt=3.55
  effHYB[36] = 0.825619;//pt=3.65
  effHYB[37] = 0.831917;//pt=3.75
  effHYB[38] = 0.82425;//pt=3.85
  effHYB[39] = 0.830073;//pt=3.95
  effHYB[40] = 0.82808;//pt=4.05
  effHYB[41] = 0.829005;//pt=4.15
  effHYB[42] = 0.831138;//pt=4.25
  effHYB[43] = 0.827304;//pt=4.35
  effHYB[44] = 0.830745;//pt=4.45
  effHYB[45] = 0.832802;//pt=4.55
  effHYB[46] = 0.828073;//pt=4.65
  effHYB[47] = 0.832246;//pt=4.75
  effHYB[48] = 0.82791;//pt=4.85
  effHYB[49] = 0.831579;//pt=4.95

  max = effHYB[0];
  maxbin = 0;
  for(int i=0; i<49; i++)
    {
      if(effHYB[i] > max)
	{
	  max = effHYB[i];
	  maxbin = i;
	}
    }
  cout<<"max for HYB is "<<max<<endl;
  cout<<"maxbin for HYB is "<<maxbin<<endl;
  for(int i=0; i<50; i++)
    {
      effHYB[i] /= max;
    }
 
  // ------------------------------------
  // --- now for mean charge subtractions
  // ------------------------------------

  // --- filter bit 1 (global tracks)
  MCEF1_fb1[0] = 0.00591299;	MCEF3_fb1[0] = 0.00383748;//eta is -0.99
  MCEF1_fb1[1] = 0.00588992;	MCEF3_fb1[1] = 0.003667;//eta is -0.97
  MCEF1_fb1[2] = 0.00570083;	MCEF3_fb1[2] = 0.00361253;//eta is -0.95
  MCEF1_fb1[3] = 0.00560151;	MCEF3_fb1[3] = 0.00333053;//eta is -0.93
  MCEF1_fb1[4] = 0.00576251;	MCEF3_fb1[4] = 0.00340385;//eta is -0.91
  MCEF1_fb1[5] = 0.00554351;	MCEF3_fb1[5] = 0.00292592;//eta is -0.89
  MCEF1_fb1[6] = 0.0056127;	MCEF3_fb1[6] = 0.00332485;//eta is -0.87
  MCEF1_fb1[7] = 0.00541225;	MCEF3_fb1[7] = 0.00331356;//eta is -0.85
  MCEF1_fb1[8] = 0.00525774;	MCEF3_fb1[8] = 0.00370341;//eta is -0.83
  MCEF1_fb1[9] = 0.00535911;	MCEF3_fb1[9] = 0.00380439;//eta is -0.81
  MCEF1_fb1[10] = 0.00532029;	MCEF3_fb1[10] = 0.00418114;//eta is -0.79
  MCEF1_fb1[11] = 0.00542025;	MCEF3_fb1[11] = 0.00387263;//eta is -0.77
  MCEF1_fb1[12] = 0.00520206;	MCEF3_fb1[12] = 0.00414251;//eta is -0.75
  MCEF1_fb1[13] = 0.00513343;	MCEF3_fb1[13] = 0.00400127;//eta is -0.73
  MCEF1_fb1[14] = 0.00543818;	MCEF3_fb1[14] = 0.00416512;//eta is -0.71
  MCEF1_fb1[15] = 0.00553948;	MCEF3_fb1[15] = 0.00406533;//eta is -0.69
  MCEF1_fb1[16] = 0.00556544;	MCEF3_fb1[16] = 0.00428908;//eta is -0.67
  MCEF1_fb1[17] = 0.00558886;	MCEF3_fb1[17] = 0.00435018;//eta is -0.65
  MCEF1_fb1[18] = 0.0057369;	MCEF3_fb1[18] = 0.00445061;//eta is -0.63
  MCEF1_fb1[19] = 0.00616797;	MCEF3_fb1[19] = 0.00443836;//eta is -0.61
  MCEF1_fb1[20] = 0.00604426;	MCEF3_fb1[20] = 0.00466612;//eta is -0.59
  MCEF1_fb1[21] = 0.00581836;	MCEF3_fb1[21] = 0.00473329;//eta is -0.57
  MCEF1_fb1[22] = 0.00573374;	MCEF3_fb1[22] = 0.00429524;//eta is -0.55
  MCEF1_fb1[23] = 0.00563837;	MCEF3_fb1[23] = 0.0047092;//eta is -0.53
  MCEF1_fb1[24] = 0.00628575;	MCEF3_fb1[24] = 0.0050631;//eta is -0.51
  MCEF1_fb1[25] = 0.00606812;	MCEF3_fb1[25] = 0.00480595;//eta is -0.49
  MCEF1_fb1[26] = 0.00614111;	MCEF3_fb1[26] = 0.00474149;//eta is -0.47
  MCEF1_fb1[27] = 0.0060786;	MCEF3_fb1[27] = 0.00478341;//eta is -0.45
  MCEF1_fb1[28] = 0.00621897;	MCEF3_fb1[28] = 0.00499422;//eta is -0.43
  MCEF1_fb1[29] = 0.00633416;	MCEF3_fb1[29] = 0.00498704;//eta is -0.41
  MCEF1_fb1[30] = 0.00638646;	MCEF3_fb1[30] = 0.00503349;//eta is -0.39
  MCEF1_fb1[31] = 0.00658823;	MCEF3_fb1[31] = 0.00485434;//eta is -0.37
  MCEF1_fb1[32] = 0.0064771;	MCEF3_fb1[32] = 0.00482841;//eta is -0.35
  MCEF1_fb1[33] = 0.00694801;	MCEF3_fb1[33] = 0.00471101;//eta is -0.33
  MCEF1_fb1[34] = 0.00636795;	MCEF3_fb1[34] = 0.00491829;//eta is -0.31
  MCEF1_fb1[35] = 0.00650891;	MCEF3_fb1[35] = 0.00498886;//eta is -0.29
  MCEF1_fb1[36] = 0.00653225;	MCEF3_fb1[36] = 0.00481843;//eta is -0.27
  MCEF1_fb1[37] = 0.00657161;	MCEF3_fb1[37] = 0.00511179;//eta is -0.25
  MCEF1_fb1[38] = 0.00640258;	MCEF3_fb1[38] = 0.00511465;//eta is -0.23
  MCEF1_fb1[39] = 0.00626697;	MCEF3_fb1[39] = 0.00563881;//eta is -0.21
  MCEF1_fb1[40] = 0.00691174;	MCEF3_fb1[40] = 0.00534971;//eta is -0.19
  MCEF1_fb1[41] = 0.00620562;	MCEF3_fb1[41] = 0.00557349;//eta is -0.17
  MCEF1_fb1[42] = 0.00649612;	MCEF3_fb1[42] = 0.00498584;//eta is -0.15
  MCEF1_fb1[43] = 0.0065798;	MCEF3_fb1[43] = 0.00542058;//eta is -0.13
  MCEF1_fb1[44] = 0.00642303;	MCEF3_fb1[44] = 0.00613666;//eta is -0.11
  MCEF1_fb1[45] = 0.00653218;	MCEF3_fb1[45] = 0.00546797;//eta is -0.09
  MCEF1_fb1[46] = 0.0064394;	MCEF3_fb1[46] = 0.00555092;//eta is -0.07
  MCEF1_fb1[47] = 0.00622478;	MCEF3_fb1[47] = 0.00629308;//eta is -0.05
  MCEF1_fb1[48] = 0.00639885;	MCEF3_fb1[48] = 0.0063706;//eta is -0.03
  MCEF1_fb1[49] = 0.00666821;	MCEF3_fb1[49] = 0.00631535;//eta is -0.01
  MCEF1_fb1[50] = 0.00653648;	MCEF3_fb1[50] = 0.00710942;//eta is 0.01
  MCEF1_fb1[51] = 0.00638104;	MCEF3_fb1[51] = 0.00709238;//eta is 0.03
  MCEF1_fb1[52] = 0.0066302;	MCEF3_fb1[52] = 0.00752929;//eta is 0.05
  MCEF1_fb1[53] = 0.00660631;	MCEF3_fb1[53] = 0.00756394;//eta is 0.07
  MCEF1_fb1[54] = 0.0068361;	MCEF3_fb1[54] = 0.00763778;//eta is 0.09
  MCEF1_fb1[55] = 0.00680126;	MCEF3_fb1[55] = 0.00766098;//eta is 0.11
  MCEF1_fb1[56] = 0.00673586;	MCEF3_fb1[56] = 0.0077562;//eta is 0.13
  MCEF1_fb1[57] = 0.0066831;	MCEF3_fb1[57] = 0.00749726;//eta is 0.15
  MCEF1_fb1[58] = 0.00632194;	MCEF3_fb1[58] = 0.00814879;//eta is 0.17
  MCEF1_fb1[59] = 0.00615406;	MCEF3_fb1[59] = 0.00773741;//eta is 0.19
  MCEF1_fb1[60] = 0.00640291;	MCEF3_fb1[60] = 0.00709483;//eta is 0.21
  MCEF1_fb1[61] = 0.00626011;	MCEF3_fb1[61] = 0.00748346;//eta is 0.23
  MCEF1_fb1[62] = 0.00588407;	MCEF3_fb1[62] = 0.00746276;//eta is 0.25
  MCEF1_fb1[63] = 0.00615514;	MCEF3_fb1[63] = 0.0073329;//eta is 0.27
  MCEF1_fb1[64] = 0.00581758;	MCEF3_fb1[64] = 0.00749753;//eta is 0.29
  MCEF1_fb1[65] = 0.00552913;	MCEF3_fb1[65] = 0.00758203;//eta is 0.31
  MCEF1_fb1[66] = 0.00618476;	MCEF3_fb1[66] = 0.00736015;//eta is 0.33
  MCEF1_fb1[67] = 0.00571961;	MCEF3_fb1[67] = 0.00717504;//eta is 0.35
  MCEF1_fb1[68] = 0.00596315;	MCEF3_fb1[68] = 0.00731385;//eta is 0.37
  MCEF1_fb1[69] = 0.00617447;	MCEF3_fb1[69] = 0.00715388;//eta is 0.39
  MCEF1_fb1[70] = 0.00548738;	MCEF3_fb1[70] = 0.00691641;//eta is 0.41
  MCEF1_fb1[71] = 0.00526914;	MCEF3_fb1[71] = 0.00712093;//eta is 0.43
  MCEF1_fb1[72] = 0.00549674;	MCEF3_fb1[72] = 0.00689862;//eta is 0.45
  MCEF1_fb1[73] = 0.00556518;	MCEF3_fb1[73] = 0.00701435;//eta is 0.47
  MCEF1_fb1[74] = 0.005665;	MCEF3_fb1[74] = 0.00676622;//eta is 0.49
  MCEF1_fb1[75] = 0.00552991;	MCEF3_fb1[75] = 0.00679297;//eta is 0.51
  MCEF1_fb1[76] = 0.00569837;	MCEF3_fb1[76] = 0.00662394;//eta is 0.53
  MCEF1_fb1[77] = 0.00562372;	MCEF3_fb1[77] = 0.00638628;//eta is 0.55
  MCEF1_fb1[78] = 0.00543001;	MCEF3_fb1[78] = 0.00657229;//eta is 0.57
  MCEF1_fb1[79] = 0.00524306;	MCEF3_fb1[79] = 0.00635786;//eta is 0.59
  MCEF1_fb1[80] = 0.00535293;	MCEF3_fb1[80] = 0.00641706;//eta is 0.61
  MCEF1_fb1[81] = 0.00514873;	MCEF3_fb1[81] = 0.00636306;//eta is 0.63
  MCEF1_fb1[82] = 0.00510695;	MCEF3_fb1[82] = 0.00638945;//eta is 0.65
  MCEF1_fb1[83] = 0.00484419;	MCEF3_fb1[83] = 0.0067127;//eta is 0.67
  MCEF1_fb1[84] = 0.00484068;	MCEF3_fb1[84] = 0.00621438;//eta is 0.69
  MCEF1_fb1[85] = 0.00480704;	MCEF3_fb1[85] = 0.00606007;//eta is 0.71
  MCEF1_fb1[86] = 0.00478109;	MCEF3_fb1[86] = 0.00596953;//eta is 0.73
  MCEF1_fb1[87] = 0.00488645;	MCEF3_fb1[87] = 0.00594765;//eta is 0.75
  MCEF1_fb1[88] = 0.00459445;	MCEF3_fb1[88] = 0.00579464;//eta is 0.77
  MCEF1_fb1[89] = 0.00463147;	MCEF3_fb1[89] = 0.0056488;//eta is 0.79
  MCEF1_fb1[90] = 0.00464974;	MCEF3_fb1[90] = 0.00567165;//eta is 0.81
  MCEF1_fb1[91] = 0.00458123;	MCEF3_fb1[91] = 0.00532498;//eta is 0.83
  MCEF1_fb1[92] = 0.0047586;	MCEF3_fb1[92] = 0.00472428;//eta is 0.85
  MCEF1_fb1[93] = 0.00492789;	MCEF3_fb1[93] = 0.00514747;//eta is 0.87
  MCEF1_fb1[94] = 0.00464254;	MCEF3_fb1[94] = 0.00544923;//eta is 0.89
  MCEF1_fb1[95] = 0.00453585;	MCEF3_fb1[95] = 0.00511096;//eta is 0.91
  MCEF1_fb1[96] = 0.00427511;	MCEF3_fb1[96] = 0.00556982;//eta is 0.93
  MCEF1_fb1[97] = 0.00463089;	MCEF3_fb1[97] = 0.00584684;//eta is 0.95
  MCEF1_fb1[98] = 0.00446822;	MCEF3_fb1[98] = 0.00615053;//eta is 0.97
  MCEF1_fb1[99] = 0.00433282;	MCEF3_fb1[99] = 0.0063621;//eta is 0.99

  // --- filter bit 128 (TPC only tracks)
  MCEF1_fb128[0] = 0.00591299;	MCEF3_fb128[0] = 0.00383748;//eta is -0.99
  MCEF1_fb128[1] = 0.00588992;	MCEF3_fb128[1] = 0.003667;//eta is -0.97
  MCEF1_fb128[2] = 0.00570083;	MCEF3_fb128[2] = 0.00361253;//eta is -0.95
  MCEF1_fb128[3] = 0.00560151;	MCEF3_fb128[3] = 0.00333053;//eta is -0.93
  MCEF1_fb128[4] = 0.00576251;	MCEF3_fb128[4] = 0.00340385;//eta is -0.91
  MCEF1_fb128[5] = 0.00554351;	MCEF3_fb128[5] = 0.00292592;//eta is -0.89
  MCEF1_fb128[6] = 0.0056127;	MCEF3_fb128[6] = 0.00332485;//eta is -0.87
  MCEF1_fb128[7] = 0.00541225;	MCEF3_fb128[7] = 0.00331356;//eta is -0.85
  MCEF1_fb128[8] = 0.00525774;	MCEF3_fb128[8] = 0.00370341;//eta is -0.83
  MCEF1_fb128[9] = 0.00535911;	MCEF3_fb128[9] = 0.00380439;//eta is -0.81
  MCEF1_fb128[10] = 0.00532029;	MCEF3_fb128[10] = 0.00418114;//eta is -0.79
  MCEF1_fb128[11] = 0.00542025;	MCEF3_fb128[11] = 0.00387263;//eta is -0.77
  MCEF1_fb128[12] = 0.00520206;	MCEF3_fb128[12] = 0.00414251;//eta is -0.75
  MCEF1_fb128[13] = 0.00513343;	MCEF3_fb128[13] = 0.00400127;//eta is -0.73
  MCEF1_fb128[14] = 0.00543818;	MCEF3_fb128[14] = 0.00416512;//eta is -0.71
  MCEF1_fb128[15] = 0.00553948;	MCEF3_fb128[15] = 0.00406533;//eta is -0.69
  MCEF1_fb128[16] = 0.00556544;	MCEF3_fb128[16] = 0.00428908;//eta is -0.67
  MCEF1_fb128[17] = 0.00558886;	MCEF3_fb128[17] = 0.00435018;//eta is -0.65
  MCEF1_fb128[18] = 0.0057369;	MCEF3_fb128[18] = 0.00445061;//eta is -0.63
  MCEF1_fb128[19] = 0.00616797;	MCEF3_fb128[19] = 0.00443836;//eta is -0.61
  MCEF1_fb128[20] = 0.00604426;	MCEF3_fb128[20] = 0.00466612;//eta is -0.59
  MCEF1_fb128[21] = 0.00581836;	MCEF3_fb128[21] = 0.00473329;//eta is -0.57
  MCEF1_fb128[22] = 0.00573374;	MCEF3_fb128[22] = 0.00429524;//eta is -0.55
  MCEF1_fb128[23] = 0.00563837;	MCEF3_fb128[23] = 0.0047092;//eta is -0.53
  MCEF1_fb128[24] = 0.00628575;	MCEF3_fb128[24] = 0.0050631;//eta is -0.51
  MCEF1_fb128[25] = 0.00606812;	MCEF3_fb128[25] = 0.00480595;//eta is -0.49
  MCEF1_fb128[26] = 0.00614111;	MCEF3_fb128[26] = 0.00474149;//eta is -0.47
  MCEF1_fb128[27] = 0.0060786;	MCEF3_fb128[27] = 0.00478341;//eta is -0.45
  MCEF1_fb128[28] = 0.00621897;	MCEF3_fb128[28] = 0.00499422;//eta is -0.43
  MCEF1_fb128[29] = 0.00633416;	MCEF3_fb128[29] = 0.00498704;//eta is -0.41
  MCEF1_fb128[30] = 0.00638646;	MCEF3_fb128[30] = 0.00503349;//eta is -0.39
  MCEF1_fb128[31] = 0.00658823;	MCEF3_fb128[31] = 0.00485434;//eta is -0.37
  MCEF1_fb128[32] = 0.0064771;	MCEF3_fb128[32] = 0.00482841;//eta is -0.35
  MCEF1_fb128[33] = 0.00694801;	MCEF3_fb128[33] = 0.00471101;//eta is -0.33
  MCEF1_fb128[34] = 0.00636795;	MCEF3_fb128[34] = 0.00491829;//eta is -0.31
  MCEF1_fb128[35] = 0.00650891;	MCEF3_fb128[35] = 0.00498886;//eta is -0.29
  MCEF1_fb128[36] = 0.00653225;	MCEF3_fb128[36] = 0.00481843;//eta is -0.27
  MCEF1_fb128[37] = 0.00657161;	MCEF3_fb128[37] = 0.00511179;//eta is -0.25
  MCEF1_fb128[38] = 0.00640258;	MCEF3_fb128[38] = 0.00511465;//eta is -0.23
  MCEF1_fb128[39] = 0.00626697;	MCEF3_fb128[39] = 0.00563881;//eta is -0.21
  MCEF1_fb128[40] = 0.00691174;	MCEF3_fb128[40] = 0.00534971;//eta is -0.19
  MCEF1_fb128[41] = 0.00620562;	MCEF3_fb128[41] = 0.00557349;//eta is -0.17
  MCEF1_fb128[42] = 0.00649612;	MCEF3_fb128[42] = 0.00498584;//eta is -0.15
  MCEF1_fb128[43] = 0.0065798;	MCEF3_fb128[43] = 0.00542058;//eta is -0.13
  MCEF1_fb128[44] = 0.00642303;	MCEF3_fb128[44] = 0.00613666;//eta is -0.11
  MCEF1_fb128[45] = 0.00653218;	MCEF3_fb128[45] = 0.00546797;//eta is -0.09
  MCEF1_fb128[46] = 0.0064394;	MCEF3_fb128[46] = 0.00555092;//eta is -0.07
  MCEF1_fb128[47] = 0.00622478;	MCEF3_fb128[47] = 0.00629308;//eta is -0.05
  MCEF1_fb128[48] = 0.00639885;	MCEF3_fb128[48] = 0.0063706;//eta is -0.03
  MCEF1_fb128[49] = 0.00666821;	MCEF3_fb128[49] = 0.00631535;//eta is -0.01
  MCEF1_fb128[50] = 0.00653648;	MCEF3_fb128[50] = 0.00710942;//eta is 0.01
  MCEF1_fb128[51] = 0.00638104;	MCEF3_fb128[51] = 0.00709238;//eta is 0.03
  MCEF1_fb128[52] = 0.0066302;	MCEF3_fb128[52] = 0.00752929;//eta is 0.05
  MCEF1_fb128[53] = 0.00660631;	MCEF3_fb128[53] = 0.00756394;//eta is 0.07
  MCEF1_fb128[54] = 0.0068361;	MCEF3_fb128[54] = 0.00763778;//eta is 0.09
  MCEF1_fb128[55] = 0.00680126;	MCEF3_fb128[55] = 0.00766098;//eta is 0.11
  MCEF1_fb128[56] = 0.00673586;	MCEF3_fb128[56] = 0.0077562;//eta is 0.13
  MCEF1_fb128[57] = 0.0066831;	MCEF3_fb128[57] = 0.00749726;//eta is 0.15
  MCEF1_fb128[58] = 0.00632194;	MCEF3_fb128[58] = 0.00814879;//eta is 0.17
  MCEF1_fb128[59] = 0.00615406;	MCEF3_fb128[59] = 0.00773741;//eta is 0.19
  MCEF1_fb128[60] = 0.00640291;	MCEF3_fb128[60] = 0.00709483;//eta is 0.21
  MCEF1_fb128[61] = 0.00626011;	MCEF3_fb128[61] = 0.00748346;//eta is 0.23
  MCEF1_fb128[62] = 0.00588407;	MCEF3_fb128[62] = 0.00746276;//eta is 0.25
  MCEF1_fb128[63] = 0.00615514;	MCEF3_fb128[63] = 0.0073329;//eta is 0.27
  MCEF1_fb128[64] = 0.00581758;	MCEF3_fb128[64] = 0.00749753;//eta is 0.29
  MCEF1_fb128[65] = 0.00552913;	MCEF3_fb128[65] = 0.00758203;//eta is 0.31
  MCEF1_fb128[66] = 0.00618476;	MCEF3_fb128[66] = 0.00736015;//eta is 0.33
  MCEF1_fb128[67] = 0.00571961;	MCEF3_fb128[67] = 0.00717504;//eta is 0.35
  MCEF1_fb128[68] = 0.00596315;	MCEF3_fb128[68] = 0.00731385;//eta is 0.37
  MCEF1_fb128[69] = 0.00617447;	MCEF3_fb128[69] = 0.00715388;//eta is 0.39
  MCEF1_fb128[70] = 0.00548738;	MCEF3_fb128[70] = 0.00691641;//eta is 0.41
  MCEF1_fb128[71] = 0.00526914;	MCEF3_fb128[71] = 0.00712093;//eta is 0.43
  MCEF1_fb128[72] = 0.00549674;	MCEF3_fb128[72] = 0.00689862;//eta is 0.45
  MCEF1_fb128[73] = 0.00556518;	MCEF3_fb128[73] = 0.00701435;//eta is 0.47
  MCEF1_fb128[74] = 0.005665;	MCEF3_fb128[74] = 0.00676622;//eta is 0.49
  MCEF1_fb128[75] = 0.00552991;	MCEF3_fb128[75] = 0.00679297;//eta is 0.51
  MCEF1_fb128[76] = 0.00569837;	MCEF3_fb128[76] = 0.00662394;//eta is 0.53
  MCEF1_fb128[77] = 0.00562372;	MCEF3_fb128[77] = 0.00638628;//eta is 0.55
  MCEF1_fb128[78] = 0.00543001;	MCEF3_fb128[78] = 0.00657229;//eta is 0.57
  MCEF1_fb128[79] = 0.00524306;	MCEF3_fb128[79] = 0.00635786;//eta is 0.59
  MCEF1_fb128[80] = 0.00535293;	MCEF3_fb128[80] = 0.00641706;//eta is 0.61
  MCEF1_fb128[81] = 0.00514873;	MCEF3_fb128[81] = 0.00636306;//eta is 0.63
  MCEF1_fb128[82] = 0.00510695;	MCEF3_fb128[82] = 0.00638945;//eta is 0.65
  MCEF1_fb128[83] = 0.00484419;	MCEF3_fb128[83] = 0.0067127;//eta is 0.67
  MCEF1_fb128[84] = 0.00484068;	MCEF3_fb128[84] = 0.00621438;//eta is 0.69
  MCEF1_fb128[85] = 0.00480704;	MCEF3_fb128[85] = 0.00606007;//eta is 0.71
  MCEF1_fb128[86] = 0.00478109;	MCEF3_fb128[86] = 0.00596953;//eta is 0.73
  MCEF1_fb128[87] = 0.00488645;	MCEF3_fb128[87] = 0.00594765;//eta is 0.75
  MCEF1_fb128[88] = 0.00459445;	MCEF3_fb128[88] = 0.00579464;//eta is 0.77
  MCEF1_fb128[89] = 0.00463147;	MCEF3_fb128[89] = 0.0056488;//eta is 0.79
  MCEF1_fb128[90] = 0.00464974;	MCEF3_fb128[90] = 0.00567165;//eta is 0.81
  MCEF1_fb128[91] = 0.00458123;	MCEF3_fb128[91] = 0.00532498;//eta is 0.83
  MCEF1_fb128[92] = 0.0047586;	MCEF3_fb128[92] = 0.00472428;//eta is 0.85
  MCEF1_fb128[93] = 0.00492789;	MCEF3_fb128[93] = 0.00514747;//eta is 0.87
  MCEF1_fb128[94] = 0.00464254;	MCEF3_fb128[94] = 0.00544923;//eta is 0.89
  MCEF1_fb128[95] = 0.00453585;	MCEF3_fb128[95] = 0.00511096;//eta is 0.91
  MCEF1_fb128[96] = 0.00427511;	MCEF3_fb128[96] = 0.00556982;//eta is 0.93
  MCEF1_fb128[97] = 0.00463089;	MCEF3_fb128[97] = 0.00584684;//eta is 0.95
  MCEF1_fb128[98] = 0.00446822;	MCEF3_fb128[98] = 0.00615053;//eta is 0.97
  MCEF1_fb128[99] = 0.00433282;	MCEF3_fb128[99] = 0.0063621;//eta is 0.99

  // --- filter bit 272 (hybrid tracks)
  MCEF1_fb272[0] = 0.00591299;	MCEF3_fb272[0] = 0.00383748;//eta is -0.99
  MCEF1_fb272[1] = 0.00588992;	MCEF3_fb272[1] = 0.003667;//eta is -0.97
  MCEF1_fb272[2] = 0.00570083;	MCEF3_fb272[2] = 0.00361253;//eta is -0.95
  MCEF1_fb272[3] = 0.00560151;	MCEF3_fb272[3] = 0.00333053;//eta is -0.93
  MCEF1_fb272[4] = 0.00576251;	MCEF3_fb272[4] = 0.00340385;//eta is -0.91
  MCEF1_fb272[5] = 0.00554351;	MCEF3_fb272[5] = 0.00292592;//eta is -0.89
  MCEF1_fb272[6] = 0.0056127;	MCEF3_fb272[6] = 0.00332485;//eta is -0.87
  MCEF1_fb272[7] = 0.00541225;	MCEF3_fb272[7] = 0.00331356;//eta is -0.85
  MCEF1_fb272[8] = 0.00525774;	MCEF3_fb272[8] = 0.00370341;//eta is -0.83
  MCEF1_fb272[9] = 0.00535911;	MCEF3_fb272[9] = 0.00380439;//eta is -0.81
  MCEF1_fb272[10] = 0.00532029;	MCEF3_fb272[10] = 0.00418114;//eta is -0.79
  MCEF1_fb272[11] = 0.00542025;	MCEF3_fb272[11] = 0.00387263;//eta is -0.77
  MCEF1_fb272[12] = 0.00520206;	MCEF3_fb272[12] = 0.00414251;//eta is -0.75
  MCEF1_fb272[13] = 0.00513343;	MCEF3_fb272[13] = 0.00400127;//eta is -0.73
  MCEF1_fb272[14] = 0.00543818;	MCEF3_fb272[14] = 0.00416512;//eta is -0.71
  MCEF1_fb272[15] = 0.00553948;	MCEF3_fb272[15] = 0.00406533;//eta is -0.69
  MCEF1_fb272[16] = 0.00556544;	MCEF3_fb272[16] = 0.00428908;//eta is -0.67
  MCEF1_fb272[17] = 0.00558886;	MCEF3_fb272[17] = 0.00435018;//eta is -0.65
  MCEF1_fb272[18] = 0.0057369;	MCEF3_fb272[18] = 0.00445061;//eta is -0.63
  MCEF1_fb272[19] = 0.00616797;	MCEF3_fb272[19] = 0.00443836;//eta is -0.61
  MCEF1_fb272[20] = 0.00604426;	MCEF3_fb272[20] = 0.00466612;//eta is -0.59
  MCEF1_fb272[21] = 0.00581836;	MCEF3_fb272[21] = 0.00473329;//eta is -0.57
  MCEF1_fb272[22] = 0.00573374;	MCEF3_fb272[22] = 0.00429524;//eta is -0.55
  MCEF1_fb272[23] = 0.00563837;	MCEF3_fb272[23] = 0.0047092;//eta is -0.53
  MCEF1_fb272[24] = 0.00628575;	MCEF3_fb272[24] = 0.0050631;//eta is -0.51
  MCEF1_fb272[25] = 0.00606812;	MCEF3_fb272[25] = 0.00480595;//eta is -0.49
  MCEF1_fb272[26] = 0.00614111;	MCEF3_fb272[26] = 0.00474149;//eta is -0.47
  MCEF1_fb272[27] = 0.0060786;	MCEF3_fb272[27] = 0.00478341;//eta is -0.45
  MCEF1_fb272[28] = 0.00621897;	MCEF3_fb272[28] = 0.00499422;//eta is -0.43
  MCEF1_fb272[29] = 0.00633416;	MCEF3_fb272[29] = 0.00498704;//eta is -0.41
  MCEF1_fb272[30] = 0.00638646;	MCEF3_fb272[30] = 0.00503349;//eta is -0.39
  MCEF1_fb272[31] = 0.00658823;	MCEF3_fb272[31] = 0.00485434;//eta is -0.37
  MCEF1_fb272[32] = 0.0064771;	MCEF3_fb272[32] = 0.00482841;//eta is -0.35
  MCEF1_fb272[33] = 0.00694801;	MCEF3_fb272[33] = 0.00471101;//eta is -0.33
  MCEF1_fb272[34] = 0.00636795;	MCEF3_fb272[34] = 0.00491829;//eta is -0.31
  MCEF1_fb272[35] = 0.00650891;	MCEF3_fb272[35] = 0.00498886;//eta is -0.29
  MCEF1_fb272[36] = 0.00653225;	MCEF3_fb272[36] = 0.00481843;//eta is -0.27
  MCEF1_fb272[37] = 0.00657161;	MCEF3_fb272[37] = 0.00511179;//eta is -0.25
  MCEF1_fb272[38] = 0.00640258;	MCEF3_fb272[38] = 0.00511465;//eta is -0.23
  MCEF1_fb272[39] = 0.00626697;	MCEF3_fb272[39] = 0.00563881;//eta is -0.21
  MCEF1_fb272[40] = 0.00691174;	MCEF3_fb272[40] = 0.00534971;//eta is -0.19
  MCEF1_fb272[41] = 0.00620562;	MCEF3_fb272[41] = 0.00557349;//eta is -0.17
  MCEF1_fb272[42] = 0.00649612;	MCEF3_fb272[42] = 0.00498584;//eta is -0.15
  MCEF1_fb272[43] = 0.0065798;	MCEF3_fb272[43] = 0.00542058;//eta is -0.13
  MCEF1_fb272[44] = 0.00642303;	MCEF3_fb272[44] = 0.00613666;//eta is -0.11
  MCEF1_fb272[45] = 0.00653218;	MCEF3_fb272[45] = 0.00546797;//eta is -0.09
  MCEF1_fb272[46] = 0.0064394;	MCEF3_fb272[46] = 0.00555092;//eta is -0.07
  MCEF1_fb272[47] = 0.00622478;	MCEF3_fb272[47] = 0.00629308;//eta is -0.05
  MCEF1_fb272[48] = 0.00639885;	MCEF3_fb272[48] = 0.0063706;//eta is -0.03
  MCEF1_fb272[49] = 0.00666821;	MCEF3_fb272[49] = 0.00631535;//eta is -0.01
  MCEF1_fb272[50] = 0.00653648;	MCEF3_fb272[50] = 0.00710942;//eta is 0.01
  MCEF1_fb272[51] = 0.00638104;	MCEF3_fb272[51] = 0.00709238;//eta is 0.03
  MCEF1_fb272[52] = 0.0066302;	MCEF3_fb272[52] = 0.00752929;//eta is 0.05
  MCEF1_fb272[53] = 0.00660631;	MCEF3_fb272[53] = 0.00756394;//eta is 0.07
  MCEF1_fb272[54] = 0.0068361;	MCEF3_fb272[54] = 0.00763778;//eta is 0.09
  MCEF1_fb272[55] = 0.00680126;	MCEF3_fb272[55] = 0.00766098;//eta is 0.11
  MCEF1_fb272[56] = 0.00673586;	MCEF3_fb272[56] = 0.0077562;//eta is 0.13
  MCEF1_fb272[57] = 0.0066831;	MCEF3_fb272[57] = 0.00749726;//eta is 0.15
  MCEF1_fb272[58] = 0.00632194;	MCEF3_fb272[58] = 0.00814879;//eta is 0.17
  MCEF1_fb272[59] = 0.00615406;	MCEF3_fb272[59] = 0.00773741;//eta is 0.19
  MCEF1_fb272[60] = 0.00640291;	MCEF3_fb272[60] = 0.00709483;//eta is 0.21
  MCEF1_fb272[61] = 0.00626011;	MCEF3_fb272[61] = 0.00748346;//eta is 0.23
  MCEF1_fb272[62] = 0.00588407;	MCEF3_fb272[62] = 0.00746276;//eta is 0.25
  MCEF1_fb272[63] = 0.00615514;	MCEF3_fb272[63] = 0.0073329;//eta is 0.27
  MCEF1_fb272[64] = 0.00581758;	MCEF3_fb272[64] = 0.00749753;//eta is 0.29
  MCEF1_fb272[65] = 0.00552913;	MCEF3_fb272[65] = 0.00758203;//eta is 0.31
  MCEF1_fb272[66] = 0.00618476;	MCEF3_fb272[66] = 0.00736015;//eta is 0.33
  MCEF1_fb272[67] = 0.00571961;	MCEF3_fb272[67] = 0.00717504;//eta is 0.35
  MCEF1_fb272[68] = 0.00596315;	MCEF3_fb272[68] = 0.00731385;//eta is 0.37
  MCEF1_fb272[69] = 0.00617447;	MCEF3_fb272[69] = 0.00715388;//eta is 0.39
  MCEF1_fb272[70] = 0.00548738;	MCEF3_fb272[70] = 0.00691641;//eta is 0.41
  MCEF1_fb272[71] = 0.00526914;	MCEF3_fb272[71] = 0.00712093;//eta is 0.43
  MCEF1_fb272[72] = 0.00549674;	MCEF3_fb272[72] = 0.00689862;//eta is 0.45
  MCEF1_fb272[73] = 0.00556518;	MCEF3_fb272[73] = 0.00701435;//eta is 0.47
  MCEF1_fb272[74] = 0.005665;	MCEF3_fb272[74] = 0.00676622;//eta is 0.49
  MCEF1_fb272[75] = 0.00552991;	MCEF3_fb272[75] = 0.00679297;//eta is 0.51
  MCEF1_fb272[76] = 0.00569837;	MCEF3_fb272[76] = 0.00662394;//eta is 0.53
  MCEF1_fb272[77] = 0.00562372;	MCEF3_fb272[77] = 0.00638628;//eta is 0.55
  MCEF1_fb272[78] = 0.00543001;	MCEF3_fb272[78] = 0.00657229;//eta is 0.57
  MCEF1_fb272[79] = 0.00524306;	MCEF3_fb272[79] = 0.00635786;//eta is 0.59
  MCEF1_fb272[80] = 0.00535293;	MCEF3_fb272[80] = 0.00641706;//eta is 0.61
  MCEF1_fb272[81] = 0.00514873;	MCEF3_fb272[81] = 0.00636306;//eta is 0.63
  MCEF1_fb272[82] = 0.00510695;	MCEF3_fb272[82] = 0.00638945;//eta is 0.65
  MCEF1_fb272[83] = 0.00484419;	MCEF3_fb272[83] = 0.0067127;//eta is 0.67
  MCEF1_fb272[84] = 0.00484068;	MCEF3_fb272[84] = 0.00621438;//eta is 0.69
  MCEF1_fb272[85] = 0.00480704;	MCEF3_fb272[85] = 0.00606007;//eta is 0.71
  MCEF1_fb272[86] = 0.00478109;	MCEF3_fb272[86] = 0.00596953;//eta is 0.73
  MCEF1_fb272[87] = 0.00488645;	MCEF3_fb272[87] = 0.00594765;//eta is 0.75
  MCEF1_fb272[88] = 0.00459445;	MCEF3_fb272[88] = 0.00579464;//eta is 0.77
  MCEF1_fb272[89] = 0.00463147;	MCEF3_fb272[89] = 0.0056488;//eta is 0.79
  MCEF1_fb272[90] = 0.00464974;	MCEF3_fb272[90] = 0.00567165;//eta is 0.81
  MCEF1_fb272[91] = 0.00458123;	MCEF3_fb272[91] = 0.00532498;//eta is 0.83
  MCEF1_fb272[92] = 0.0047586;	MCEF3_fb272[92] = 0.00472428;//eta is 0.85
  MCEF1_fb272[93] = 0.00492789;	MCEF3_fb272[93] = 0.00514747;//eta is 0.87
  MCEF1_fb272[94] = 0.00464254;	MCEF3_fb272[94] = 0.00544923;//eta is 0.89
  MCEF1_fb272[95] = 0.00453585;	MCEF3_fb272[95] = 0.00511096;//eta is 0.91
  MCEF1_fb272[96] = 0.00427511;	MCEF3_fb272[96] = 0.00556982;//eta is 0.93
  MCEF1_fb272[97] = 0.00463089;	MCEF3_fb272[97] = 0.00584684;//eta is 0.95
  MCEF1_fb272[98] = 0.00446822;	MCEF3_fb272[98] = 0.00615053;//eta is 0.97
  MCEF1_fb272[99] = 0.00433282;	MCEF3_fb272[99] = 0.0063621;//eta is 0.99

  // --- filter bit 1 with tight DCA cut...
  MCEF1_dca[0] = 0.00591299;	MCEF3_dca[0] = 0.00383748;//eta is -0.99
  MCEF1_dca[1] = 0.00588992;	MCEF3_dca[1] = 0.003667;//eta is -0.97
  MCEF1_dca[2] = 0.00570083;	MCEF3_dca[2] = 0.00361253;//eta is -0.95
  MCEF1_dca[3] = 0.00560151;	MCEF3_dca[3] = 0.00333053;//eta is -0.93
  MCEF1_dca[4] = 0.00576251;	MCEF3_dca[4] = 0.00340385;//eta is -0.91
  MCEF1_dca[5] = 0.00554351;	MCEF3_dca[5] = 0.00292592;//eta is -0.89
  MCEF1_dca[6] = 0.0056127;	MCEF3_dca[6] = 0.00332485;//eta is -0.87
  MCEF1_dca[7] = 0.00541225;	MCEF3_dca[7] = 0.00331356;//eta is -0.85
  MCEF1_dca[8] = 0.00525774;	MCEF3_dca[8] = 0.00370341;//eta is -0.83
  MCEF1_dca[9] = 0.00535911;	MCEF3_dca[9] = 0.00380439;//eta is -0.81
  MCEF1_dca[10] = 0.00532029;	MCEF3_dca[10] = 0.00418114;//eta is -0.79
  MCEF1_dca[11] = 0.00542025;	MCEF3_dca[11] = 0.00387263;//eta is -0.77
  MCEF1_dca[12] = 0.00520206;	MCEF3_dca[12] = 0.00414251;//eta is -0.75
  MCEF1_dca[13] = 0.00513343;	MCEF3_dca[13] = 0.00400127;//eta is -0.73
  MCEF1_dca[14] = 0.00543818;	MCEF3_dca[14] = 0.00416512;//eta is -0.71
  MCEF1_dca[15] = 0.00553948;	MCEF3_dca[15] = 0.00406533;//eta is -0.69
  MCEF1_dca[16] = 0.00556544;	MCEF3_dca[16] = 0.00428908;//eta is -0.67
  MCEF1_dca[17] = 0.00558886;	MCEF3_dca[17] = 0.00435018;//eta is -0.65
  MCEF1_dca[18] = 0.0057369;	MCEF3_dca[18] = 0.00445061;//eta is -0.63
  MCEF1_dca[19] = 0.00616797;	MCEF3_dca[19] = 0.00443836;//eta is -0.61
  MCEF1_dca[20] = 0.00604426;	MCEF3_dca[20] = 0.00466612;//eta is -0.59
  MCEF1_dca[21] = 0.00581836;	MCEF3_dca[21] = 0.00473329;//eta is -0.57
  MCEF1_dca[22] = 0.00573374;	MCEF3_dca[22] = 0.00429524;//eta is -0.55
  MCEF1_dca[23] = 0.00563837;	MCEF3_dca[23] = 0.0047092;//eta is -0.53
  MCEF1_dca[24] = 0.00628575;	MCEF3_dca[24] = 0.0050631;//eta is -0.51
  MCEF1_dca[25] = 0.00606812;	MCEF3_dca[25] = 0.00480595;//eta is -0.49
  MCEF1_dca[26] = 0.00614111;	MCEF3_dca[26] = 0.00474149;//eta is -0.47
  MCEF1_dca[27] = 0.0060786;	MCEF3_dca[27] = 0.00478341;//eta is -0.45
  MCEF1_dca[28] = 0.00621897;	MCEF3_dca[28] = 0.00499422;//eta is -0.43
  MCEF1_dca[29] = 0.00633416;	MCEF3_dca[29] = 0.00498704;//eta is -0.41
  MCEF1_dca[30] = 0.00638646;	MCEF3_dca[30] = 0.00503349;//eta is -0.39
  MCEF1_dca[31] = 0.00658823;	MCEF3_dca[31] = 0.00485434;//eta is -0.37
  MCEF1_dca[32] = 0.0064771;	MCEF3_dca[32] = 0.00482841;//eta is -0.35
  MCEF1_dca[33] = 0.00694801;	MCEF3_dca[33] = 0.00471101;//eta is -0.33
  MCEF1_dca[34] = 0.00636795;	MCEF3_dca[34] = 0.00491829;//eta is -0.31
  MCEF1_dca[35] = 0.00650891;	MCEF3_dca[35] = 0.00498886;//eta is -0.29
  MCEF1_dca[36] = 0.00653225;	MCEF3_dca[36] = 0.00481843;//eta is -0.27
  MCEF1_dca[37] = 0.00657161;	MCEF3_dca[37] = 0.00511179;//eta is -0.25
  MCEF1_dca[38] = 0.00640258;	MCEF3_dca[38] = 0.00511465;//eta is -0.23
  MCEF1_dca[39] = 0.00626697;	MCEF3_dca[39] = 0.00563881;//eta is -0.21
  MCEF1_dca[40] = 0.00691174;	MCEF3_dca[40] = 0.00534971;//eta is -0.19
  MCEF1_dca[41] = 0.00620562;	MCEF3_dca[41] = 0.00557349;//eta is -0.17
  MCEF1_dca[42] = 0.00649612;	MCEF3_dca[42] = 0.00498584;//eta is -0.15
  MCEF1_dca[43] = 0.0065798;	MCEF3_dca[43] = 0.00542058;//eta is -0.13
  MCEF1_dca[44] = 0.00642303;	MCEF3_dca[44] = 0.00613666;//eta is -0.11
  MCEF1_dca[45] = 0.00653218;	MCEF3_dca[45] = 0.00546797;//eta is -0.09
  MCEF1_dca[46] = 0.0064394;	MCEF3_dca[46] = 0.00555092;//eta is -0.07
  MCEF1_dca[47] = 0.00622478;	MCEF3_dca[47] = 0.00629308;//eta is -0.05
  MCEF1_dca[48] = 0.00639885;	MCEF3_dca[48] = 0.0063706;//eta is -0.03
  MCEF1_dca[49] = 0.00666821;	MCEF3_dca[49] = 0.00631535;//eta is -0.01
  MCEF1_dca[50] = 0.00653648;	MCEF3_dca[50] = 0.00710942;//eta is 0.01
  MCEF1_dca[51] = 0.00638104;	MCEF3_dca[51] = 0.00709238;//eta is 0.03
  MCEF1_dca[52] = 0.0066302;	MCEF3_dca[52] = 0.00752929;//eta is 0.05
  MCEF1_dca[53] = 0.00660631;	MCEF3_dca[53] = 0.00756394;//eta is 0.07
  MCEF1_dca[54] = 0.0068361;	MCEF3_dca[54] = 0.00763778;//eta is 0.09
  MCEF1_dca[55] = 0.00680126;	MCEF3_dca[55] = 0.00766098;//eta is 0.11
  MCEF1_dca[56] = 0.00673586;	MCEF3_dca[56] = 0.0077562;//eta is 0.13
  MCEF1_dca[57] = 0.0066831;	MCEF3_dca[57] = 0.00749726;//eta is 0.15
  MCEF1_dca[58] = 0.00632194;	MCEF3_dca[58] = 0.00814879;//eta is 0.17
  MCEF1_dca[59] = 0.00615406;	MCEF3_dca[59] = 0.00773741;//eta is 0.19
  MCEF1_dca[60] = 0.00640291;	MCEF3_dca[60] = 0.00709483;//eta is 0.21
  MCEF1_dca[61] = 0.00626011;	MCEF3_dca[61] = 0.00748346;//eta is 0.23
  MCEF1_dca[62] = 0.00588407;	MCEF3_dca[62] = 0.00746276;//eta is 0.25
  MCEF1_dca[63] = 0.00615514;	MCEF3_dca[63] = 0.0073329;//eta is 0.27
  MCEF1_dca[64] = 0.00581758;	MCEF3_dca[64] = 0.00749753;//eta is 0.29
  MCEF1_dca[65] = 0.00552913;	MCEF3_dca[65] = 0.00758203;//eta is 0.31
  MCEF1_dca[66] = 0.00618476;	MCEF3_dca[66] = 0.00736015;//eta is 0.33
  MCEF1_dca[67] = 0.00571961;	MCEF3_dca[67] = 0.00717504;//eta is 0.35
  MCEF1_dca[68] = 0.00596315;	MCEF3_dca[68] = 0.00731385;//eta is 0.37
  MCEF1_dca[69] = 0.00617447;	MCEF3_dca[69] = 0.00715388;//eta is 0.39
  MCEF1_dca[70] = 0.00548738;	MCEF3_dca[70] = 0.00691641;//eta is 0.41
  MCEF1_dca[71] = 0.00526914;	MCEF3_dca[71] = 0.00712093;//eta is 0.43
  MCEF1_dca[72] = 0.00549674;	MCEF3_dca[72] = 0.00689862;//eta is 0.45
  MCEF1_dca[73] = 0.00556518;	MCEF3_dca[73] = 0.00701435;//eta is 0.47
  MCEF1_dca[74] = 0.005665;	MCEF3_dca[74] = 0.00676622;//eta is 0.49
  MCEF1_dca[75] = 0.00552991;	MCEF3_dca[75] = 0.00679297;//eta is 0.51
  MCEF1_dca[76] = 0.00569837;	MCEF3_dca[76] = 0.00662394;//eta is 0.53
  MCEF1_dca[77] = 0.00562372;	MCEF3_dca[77] = 0.00638628;//eta is 0.55
  MCEF1_dca[78] = 0.00543001;	MCEF3_dca[78] = 0.00657229;//eta is 0.57
  MCEF1_dca[79] = 0.00524306;	MCEF3_dca[79] = 0.00635786;//eta is 0.59
  MCEF1_dca[80] = 0.00535293;	MCEF3_dca[80] = 0.00641706;//eta is 0.61
  MCEF1_dca[81] = 0.00514873;	MCEF3_dca[81] = 0.00636306;//eta is 0.63
  MCEF1_dca[82] = 0.00510695;	MCEF3_dca[82] = 0.00638945;//eta is 0.65
  MCEF1_dca[83] = 0.00484419;	MCEF3_dca[83] = 0.0067127;//eta is 0.67
  MCEF1_dca[84] = 0.00484068;	MCEF3_dca[84] = 0.00621438;//eta is 0.69
  MCEF1_dca[85] = 0.00480704;	MCEF3_dca[85] = 0.00606007;//eta is 0.71
  MCEF1_dca[86] = 0.00478109;	MCEF3_dca[86] = 0.00596953;//eta is 0.73
  MCEF1_dca[87] = 0.00488645;	MCEF3_dca[87] = 0.00594765;//eta is 0.75
  MCEF1_dca[88] = 0.00459445;	MCEF3_dca[88] = 0.00579464;//eta is 0.77
  MCEF1_dca[89] = 0.00463147;	MCEF3_dca[89] = 0.0056488;//eta is 0.79
  MCEF1_dca[90] = 0.00464974;	MCEF3_dca[90] = 0.00567165;//eta is 0.81
  MCEF1_dca[91] = 0.00458123;	MCEF3_dca[91] = 0.00532498;//eta is 0.83
  MCEF1_dca[92] = 0.0047586;	MCEF3_dca[92] = 0.00472428;//eta is 0.85
  MCEF1_dca[93] = 0.00492789;	MCEF3_dca[93] = 0.00514747;//eta is 0.87
  MCEF1_dca[94] = 0.00464254;	MCEF3_dca[94] = 0.00544923;//eta is 0.89
  MCEF1_dca[95] = 0.00453585;	MCEF3_dca[95] = 0.00511096;//eta is 0.91
  MCEF1_dca[96] = 0.00427511;	MCEF3_dca[96] = 0.00556982;//eta is 0.93
  MCEF1_dca[97] = 0.00463089;	MCEF3_dca[97] = 0.00584684;//eta is 0.95
  MCEF1_dca[98] = 0.00446822;	MCEF3_dca[98] = 0.00615053;//eta is 0.97
  MCEF1_dca[99] = 0.00433282;	MCEF3_dca[99] = 0.0063621;//eta is 0.99



}



float AliAnalysisTaskCMEv2A::calc3(float xn, float yn, float x2n, float y2n, float M)
{

  float Qn2 = xn*xn+yn*yn;
  float Q2n2 = x2n*x2n+y2n*y2n;
  float Qn2d = xn*xn-yn*yn;

  float W_3 = M*(M-1)*(M-2);

  float first = ((x2n*Qn2d)-(2*Qn2))/W_3;
  float second = (Q2n2-(2*M))/W_3;

  return first-second;

}



float AliAnalysisTaskCMEv2A::calc4(float xn, float yn, float x2n, float y2n, float M)
{

  float Qn2 = xn*xn+yn*yn;
  float Q2n2 = x2n*x2n+y2n*y2n;
  float Qn4 = Qn2*Qn2;
  float Qn2d = xn*xn-yn*yn;

  float W_4 = M*(M-1)*(M-2)*(M-3);

  float first = (Qn4+Q2n2-(2*(x2n*Qn2d)))/W_4;
  float second = 2*((2*(M-2)*Qn2)-(M*(M-3)))/W_4;

  return first-second;

}


// function definition
float GetDPhiStar(float phi1, float pt1, float charge1, float phi2, float pt2, float charge2, float radius, float bSign)
{ 
  //
  // calculates dphistar
  //
  float dphistar = phi1 - phi2 - charge1 * bSign * asin(0.075 * radius / pt1) + charge2 * bSign * asin(0.075 * radius / pt2);

  //static const double kPi = pi;
  float kPi = pi;

  // circularity
  //   if (dphistar > 2 * kPi)
  //     dphistar -= 2 * kPi;
  //   if (dphistar < -2 * kPi)
  //     dphistar += 2 * kPi;

  if (dphistar > kPi)
    dphistar = kPi * 2 - dphistar;
  if (dphistar < -kPi)
    dphistar = -kPi * 2 - dphistar;
  if (dphistar > kPi) // might look funny but is needed
    dphistar = kPi * 2 - dphistar;

  return dphistar;

}


