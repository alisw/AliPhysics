#include "MiniV0.h"

using Lifetimes::MiniV0;

const int MiniV0::fgkV0cosPA_n = 50000;
const float MiniV0::fgkV0cosPA_f = 0.95f;
const float MiniV0::fgkV0cosPA_l = 1.f;
const float MiniV0::fgkV0cosPA_w =
    (MiniV0::fgkV0cosPA_l - MiniV0::fgkV0cosPA_f) /
    MiniV0::fgkV0cosPA_n;

const int MiniV0::fgkV0chi2_n = 100;
const float MiniV0::fgkV0chi2_f = 0.f;
const float MiniV0::fgkV0chi2_l = 10.f;
const float MiniV0::fgkV0chi2_w =
    (MiniV0::fgkV0chi2_l - MiniV0::fgkV0chi2_f) / MiniV0::fgkV0chi2_n;

const int MiniV0::fgkDCAProng2PV_n = 250;
const float MiniV0::fgkDCAProng2PV_f = 0.f;
const float MiniV0::fgkDCAProng2PV_l = 0.25f;
const float MiniV0::fgkDCAProng2PV_w =
    (MiniV0::fgkDCAProng2PV_l - MiniV0::fgkDCAProng2PV_f) /
    MiniV0::fgkDCAProng2PV_n;

const int MiniV0::fgkDCAProngs_n = 250;
const float MiniV0::fgkDCAProngs_f = 0.f;
const float MiniV0::fgkDCAProngs_l = 2.f;
const float MiniV0::fgkDCAProngs_w =
    (MiniV0::fgkDCAProngs_l - MiniV0::fgkDCAProngs_f) /
    MiniV0::fgkDCAProngs_n;

const int MiniV0::fgkArmAlpha_n = 250;
const float MiniV0::fgkArmAlpha_f = -1.f;
const float MiniV0::fgkArmAlpha_l = 1.f;
const float MiniV0::fgkArmAlpha_w =
    (MiniV0::fgkArmAlpha_l - MiniV0::fgkArmAlpha_f) /
    MiniV0::fgkArmAlpha_n;

const int MiniV0::fgkArmPt_n = 254;
const float MiniV0::fgkArmPt_f = 0.f;
const float MiniV0::fgkArmPt_l = 0.254f;
const float MiniV0::fgkArmPt_w =
    (MiniV0::fgkArmPt_l - MiniV0::fgkArmPt_f) / MiniV0::fgkArmPt_n;

const int MiniV0::fgkXedOverFindable_n = 250;
const float MiniV0::fgkXedOverFindable_f = 0.f;
const float MiniV0::fgkXedOverFindable_l = 1.f;
const float MiniV0::fgkXedOverFindable_w =
    (MiniV0::fgkXedOverFindable_l - MiniV0::fgkXedOverFindable_f) /
    MiniV0::fgkXedOverFindable_n;

const int MiniV0::fgkChi2xCluster_n = 250;
const float MiniV0::fgkChi2xCluster_f = 0.f;
const float MiniV0::fgkChi2xCluster_l = 10.f;
const float MiniV0::fgkChi2xCluster_w =
    (MiniV0::fgkChi2xCluster_l - MiniV0::fgkChi2xCluster_f) /
    MiniV0::fgkChi2xCluster_n;

const int MiniV0::fgkTPCsigma_n = 12;
const float MiniV0::fgkTPCsigma_f = 0.f;
const float MiniV0::fgkTPCsigma_l = 6.f;
const float MiniV0::fgkTPCsigma_w =
    (MiniV0::fgkTPCsigma_l - MiniV0::fgkTPCsigma_f) /
    MiniV0::fgkTPCsigma_n;

const int MiniV0::fgkEta_n = 200;
const float MiniV0::fgkEta_f = -1.f;
const float MiniV0::fgkEta_l = 1.f;
const float MiniV0::fgkEta_w =
    (MiniV0::fgkEta_l - MiniV0::fgkEta_f) / MiniV0::fgkEta_n;