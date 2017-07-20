/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

///////////////////////////////////////////////////////////////////////////
//                Dielectron Variables Manager class                     //
//                                                                       //
/*

*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliDielectronVarManager.h"

ClassImp(AliDielectronVarManager)

const char* AliDielectronVarManager::fgkParticleNames[AliDielectronVarManager::kNMaxValues][3] = {
  {"Px",                     "#it{p}_{x}",                                         "(GeV/#it{c})"},
  {"Py",                     "#it{p}_{y}",                                         "(GeV/#it{c})"},
  {"Pz",                     "#it{p}_{z}",                                         "(GeV/#it{c})"},
  {"Pt",                     "#it{p}_{T}",                                         "(GeV/#it{c})"},
  {"PtMC",                   "#it{p}_{T} (MC)",                                    "(GeV/#it{c})"},
  {"PtSq",                   "#it{p}_{T}^{2}",                                     "((GeV/#it{c})^{2})"},
  {"P",                      "#it{p}",                                             "(GeV/#it{c})"},
  {"PMC",                    "#it{p} (MC)",                                        "(GeV/#it{c})"},
  {"Xv",                     "x_{vtx}",                                            "(cm)"},
  {"Yv",                     "y_{vtx}",                                            "(cm)"},
  {"Zv",                     "z_{vtx}",                                            "(cm)"},
  {"OneOverPt",              "1/^{}#it{p}_{T}",                                    "(GeV/#it{c})^{-1}"},
  {"Phi",                    "#phi",                                               "(rad.)"},
  {"PhiMC",                  "#phi (MC)",                                          "(rad.)"},
  {"Theta",                  "#theta",                                             "(rad.)"},
  {"Eta",                    "#eta",                                               ""},
  {"EtaMC",                  "#eta (MC)",                                          ""},
  {"Y",                      "#it{y}",                                             ""},
  {"E",                      "E",                                                  "(GeV)"},
  {"M",                      "m_{inv}",                                            "(GeV/#it{c^{2}})"},
  {"MCorr",                      "m_{inv}, corrected",                                            "(GeV/#it{c^{2}})"},
  {"MMC",                    "m_{inv} (MC)",                                       "(GeV/#it{c^{2}})"},
  {"Charge",                 "q",                                                  "(e)"},
  {"NclsITS",                "N_{cls}^{ITS}",                                      ""},
  {"ITSFakeFlag",            "ITS Fake Flag",                                      ""},
  {"ITSchi2PerCluster",      "#chi^{2}/^{}N_{cls}^{ITS}",                          ""},
  {"NclsTPC",                "N_{cls}^{TPC}",                                      ""},
  {"NclsSTPC",               "N_{shar.cls}^{TPC}",                                 ""},
  {"NclsSFracTPC",           "N_{shar.cls}^{TPC}/^{}N_{cls}^{TPC}",                ""},
  {"NclsSITS",               "N_{shar.cls}^{ITS}",                                 ""},
  {"NclsSFracITS",           "N_{shar.cls}^{ITS}/^{}N_{cls}^{ITS}",                ""},
  {"NclsSMapITS",            "ITS_shared_cluster_map",                             ""},
  {"NclsTPCiter1",           "N_{1st.iter.cls}^{TPC}",                             ""},
  {"NFclsTPC",               "N_{find.cls}^{TPC}",                                 ""},
  {"NFclsTPCrobust",         "N_{cross.rows}^{TPC}",                               ""},
  {"NFclsTPCrobustFraction", "N_{cls}^{TPC}/^{}N_{find.cls}^{TPC}",                ""},//TODO: check
  {"NFclsTPCfracCrossedRows","N_{cross.rows}^{TPC}/^{}N_{find.cls}^{TPC}",         ""},
  {"Chi2TPCConstrainedVsGlobal",             "Chi^{2}_{TPC-constrained vs. Global}",                                             ""},
  {"TPCsignalN",             "N_{d#it{E}/d#it{x} points}^{TPC}",                                       ""},
  {"TPCsignalNfrac",         "N_{d#it{E}/d#it{x} points}^{TPC}/^{}N_{cls}^{TPC}",                      ""},
  {"TPCchi2PerCluster",      "#chi^{2}/^{}N_{cls}^{TPC}",                          ""},
  {"TPCclsDiff",             "N_{d#it{E}/d#it{x} points}^{TPC} - N_{cls}^{TPC}",   ""},
  {"TPCclsSegments",         "N_{segments}^{TPC}",                                 ""},
  {"TPCclsIRO",              "N_{clsIRO}^{TPC}",                                   ""},
  {"TPCclsORO",              "N_{clsORO}^{TPC}",                                   ""},
  {"TPCactiveLength",        "active length^{TPC}",                                ""},
  {"TPCGeomLength",          "geom. length_{pTdependent}^{TPC}",                   ""},



  {"TrackStatus",            "TrackStatus",                                        ""},
  {"FilterBit",              "AOD filter bit",                                     ""},

  {"NclsTRD",                "N_{cls}^{TRD}",                                      ""},
  {"TRDntracklets",          "N_{trkl}^{TRD}",                                     ""},
  {"TRDpidQuality",          "N_{PID.trkl}^{TRD}",                                 ""},
  {"TRDchi2",                "#chi^{2}",					   ""},
  {"TRDchi2PerTracklet",     "#chi^{2}/^{}N_{trkl}^{TRD}",                         ""},
  {"TRDpidProb_Electrons",   "P(PID_{e}^{TRD})",                                   ""},
  {"TRDpidProb_Pions",       "P(PID_{#pi}^{TRD})",                                 ""},
  {"TRDpidProb2D_Electrons", "P(2dim. PID_{e}^{TRD})",                             ""},
  {"TRDpidProb2D_Pions",     "P(2dim. PID_{#pi}^{TRD})",                           ""},
  {"TRDpidProb2D_Protons",   "P(2dim. PID_{p}^{TRD})",                             ""},
  {"TRDpidProb3D_Electrons", "P(3dim. PID_{e}^{TRD})",                             ""},
  {"TRDpidProb3D_Pions",     "P(3dim. PID_{#pi}^{TRD})",                           ""},
  {"TRDpidProb3D_Protons",   "P(3dim. PID_{p}^{TRD})",                             ""},
  {"TRDpidProb7D_Electrons", "P(7dim. PID_{e}^{TRD})",	                           ""},
  {"TRDpidProb7D_Pions",     "P(7dim. PID_{#pi}^{TRD})",                           ""},
  {"TRDpidProb7D_Protons",   "P(7dim. PID_{p}^{TRD})",                             ""},
  {"TRDphi",                 "#phi^{TRD}",                                         ""},
  {"TRDpidEffLeg",           "#varepsilon^{TRD legs}(PID)",                        ""},
  {"TRD_signal",             "TRD d#it{E}/d#it{x} + TR",                           "(a.u.)"},
  {"TRDeta",                 "#eta^{TRD}",                                         ""},
  {"InTRDacceptance",        "in TRD acceptance",                                  ""},

  {"ImpactParXY",            "#it{dca}_{xy}",                                      "(cm)"},
  {"ImpactParZ",             "#it{dca}_{z}",                                       "(cm)"},
  {"ImpactParXYsigma",       "#it{dca}_{xy}",                                      "(#sigma)"},
  {"ImpactParZsigma",        "#it{dca}_{z}",                                       "(#sigma)"},
  {"TrackLength",            "#it{l}_{track}",                                     "(cm)"},
  {"DistPrimToSecVtxXYMC",   "distance prim vtx to sec vtx in XY",                 "(cm)"},
  {"DistPrimToSecVtxZMC",    "distance prim vtx to sec vtx inZ",                   "(cm)"},

  {"PdgCode",                "PDG code",                                           ""},
  {"PdgCodeMother",          "mothers PDG code",                                   ""},
  {"PdgCodeGrandMother",     "grand mothers PDG code",                             ""},
  {"HasCocktailMother",      "mother from AliGenCocktail",                         ""},
  {"HasCocktailGrandMother", "grand mother from AliGenCocktail",                   ""},

  {"NumberOfDaughters",      "N_{daughters}",                                      ""},
  {"HaveSameMother",         "HaveSameMother",                                     ""},
  {"IsJpsiPrimary",          "IsJpsiPrimary",                                      ""},
  {"NumberOfJPsisIncl",      "N_{incl. J/^{}#psi}",                                ""},
  {"NumberOfJPsisPrompt",    "N_{prompt J/^{}#psi}",                               ""},
  {"NumberOfJPsisNPrompt",   "N_{non prompt J/^{}#psi}",                           ""},

  {"ITS_signal",             "ITS d#it{E}/d#it{x}",                                "(keV/^{}300#mum)"},//TODO: check units
  {"SSD1_signal",            "SSD 1st layer d#it{E}/d#it{x}",                      "(keV/^{}300#mum)"},
  {"SSD2_signal",            "SSD 2nd layer d#it{E}/d#it{x}",                      "(keV/^{}300#mum)"},
  {"SDD1_signal",            "SDD 1st layer d#it{E}/d#it{x}",                      "(keV/^{}300#mum)"},
  {"SDD2_signal",            "SDD 2nd layer d#it{E}/d#it{x}",                      "(keV/^{}300#mum)"},
  {"ITS_clusterMap",         "ITS_clusterMap",                                     ""},
  {"ITSLayerFirstCls",       "1st ITS cluster",                                    "(layer)"},

  {"ITS_nSigma_ElectronsRaw","n#sigma_{e}^{ITS} (raw)",                            ""},
  {"ITS_nSigma_Electrons",   "n#sigma_{e}^{ITS}",                                  ""},
  {"ITS_nSigma_Pions",       "n#sigma_{#pi}^{ITS}",                                ""},
  {"ITS_nSigma_Muons",       "n#sigma_{#mu}^{ITS}",                                ""},
  {"ITS_nSigma_Kaons",       "n#sigma_{K}^{ITS}",                                  ""},
  {"ITS_nSigma_Protons",     "n#sigma_{p}^{ITS}",                                  ""},

  {"P_InnerParam",           "#it{p}_{inner wall}^{TPC}",                          "(GeV/#it{c})"},
  {"P_OuterParam",           "#it{p}_{outer wall}^{TPC}",                          "(GeV/#it{c})"},
  {"Y_signed_InnerParam",    "sign. y_{inner wall}^{TPC}",                         "(GeV/#it{c})"},
  {"TPC_signal",             "TPC d#it{E}/d#it{x}",                                "(a.u.)"},
  {"TOF_signal",             "TOF signal",                                         "(ps)"},
  {"TOF_beta",               "TOF #beta",                                          "(#it{c})"},
  {"TOF_PIDbit",             "TOF PID bit",                                        ""},
  {"TOF_mismProb",           "TOF mismatch probability"                            ""},

  {"TPC_nSigma_ElectronsRaw","n#sigma_{e}^{TPC} (raw)",                            ""},
  {"TPC_nSigma_Electrons",   "n#sigma_{e}^{TPC}",                                  ""},
  {"TPC_nSigma_Pions",       "n#sigma_{#pi}^{TPC}",                                ""},
  {"TPC_nSigma_Muons",       "n#sigma_{#mu}^{TPC}",                                ""},
  {"TPC_nSigma_Kaons",       "n#sigma_{K}^{TPC}",                                  ""},
  {"TPC_nSigma_Protons",     "n#sigma_{p}^{TPC}",                                  ""},

  {"TOF_nSigma_ElectronsRaw","n#sigma_{e}^{TOF} (raw)",                            ""},
  {"TOF_nSigma_Electrons",   "n#sigma_{e}^{TOF}",                                  ""},
  {"TOF_nSigma_Pions",       "n#sigma_{#pi}^{TOF}",                                ""},
  {"TOF_nSigma_Muons",       "n#sigma_{#mu}^{TOF}",                                ""},
  {"TOF_nSigma_Kaons",       "n#sigma_{K}^{TOF}",                                  ""},
  {"TOF_nSigma_Protons",     "n#sigma_{p}^{TOF}",                                  ""},

  {"EMCAL_nSigma_Electrons", "n#sigma_{e}^{EMCAL}",                                ""},
  {"EMCAL_EoverP",           "E^{EMCAL}/^{}#it{p}",                                "(#it{c})"},
  {"EMCAL_E",                "E^{EMCAL}",                                          "(GeV)"},
  {"EMCAL_NCells",           "N_{cell}^{EMCAL}",                                   ""},
  {"EMCAL_M02",              "M02 EMCAL showershape param.",                       ""},
  {"EMCAL_M20",              "M20 EMCAL showershape param.",                       ""},
  {"EMCAL_Dispersion",       "EMCAL dispersion param.",                            ""},

  {"EffLeg",                 "A#times#varepsilon",                                    ""},
  {"OneOverEffLeg",          "(A#times#varepsilon)^{-1}",                             ""},
  {"V0Index0",               "V0Index0",                                           ""},
  {"KinkIndex0",             "KinkIndex0",                                         ""},

  {"TRDonlineA",             "trans. offset from nominal pr. vtx.",           "(cm?)"},
  {"TRDonlineLayerMask",     "TRD online Layer Mask",                              ""},
  {"TRDonlineFirstLayer",    "1st tracklet layer assoc. to online trk.",     "(layer)"},
  {"TRDonlinePID",           "TRD online PID",                              "(0-256)"},
  {"TRDonlinePt",            "TRD online Pt",                               "(GeV/c)"},
  {"TRDonlineStack",         "Stack number of matched GTU track",                  ""},
  {"TRDonlineSector",        "Sector number of matched GTU track",                 ""},
  {"TRDonlineTrackInTime",   "Time-out Flag for GTU track",                "(boolean)"},
  {"TRDolineFlagsTiming",    "TRD online timing flags",                            ""},
  {"TRDonlineLabel",         "TRD online track label",                              ""},
  {"TRDonlineNTracklets",    "N. of tracklets assoc. to TRD online trk.",           ""},


  //
  {"Chi2NDF",                "#chi^{2}/^{}ndf",                                    ""},
  {"DecayLength",            "l_{decay}",                                          "(cm)"},
  {"R",                      "d(#vec{x}_{vtx},#vec{x}_{part.})",                   "(cm)"},
  {"OpeningAngle",           "#varphi",                                            "(rad.)"},
  {"OpeningAngleCorr",       "#varphi, corrected",                                            "(rad.)"},
  {"OpeningAngleXY",         "#varphi in XY",                                      "(rad.)"},
  {"OpeningAngleRZ",         "#varphi in RZ",                                      "(rad.)"},
  {"TriangularConversionCut","#Phi_{v}^{pair} - 21M",                              ""},
  {"CosPointingAngle",       "cos(#theta)",                                        "(rad.)"},
  {"ArmAlpha",               "#alpha^{arm}",                                       ""},
  {"ArmPt",                  "#it{p}_{T}^{arm}",                                   "(GeV/#it{c})"},
  {"ThetaHE",                "cos(#theta_{HE})",                                   ""},
  {"PhiHE",                  "#phi_{HE}",                                          "(rad.)"},
  {"ThetaSqHE",              "cos^{2}(#theta_{HE})",                               ""},
  {"Cos2PhiHE",              "cos(2#phi_{HE})",                                    ""},
  {"CosTilPhiHE",            "cos(#phi_{HE})",                                     ""},
  {"ThetaCS",                "cos(#theta_{CS})",                                   ""},
  {"PhiCS",                  "#phi_{CS}",                                          "(rad.)"},
  {"ThetaSqCS",              "cos^{2}(#theta_{CS})",                               ""},
  {"PsiPair",                "#Psi^{pair}",                                        "(rad.)"},
  {"PhivPair",               "#Phi_{v}^{pair}",                                    "(rad.)"},
  {"DeltaCotTheta",          "#Delta cot(#Theta)",                                 ""},
  {"PairPlaneAngle1A",       "#Phi_{1A}",                                          "(rad.)"},
  {"PairPlaneAngle2A",       "#Phi_{2A}",                                          "(rad.)"},
  {"PairPlaneAngle3A",       "#Phi_{3A}",                                          "(rad.)"},
  {"PairPlaneAngle4A",       "#Phi_{4A}",                                          "(rad.)"},
  {"PairPlaneAngle1C",       "#Phi_{1C}",                                          "(rad.)"},
  {"PairPlaneAngle2C",       "#Phi_{2C}",                                          "(rad.)"},
  {"PairPlaneAngle3C",       "#Phi_{3C}",                                          "(rad.)"},
  {"PairPlaneAngle4C",       "#Phi_{4C}",                                          "(rad.)"},
  {"PairPlaneAngle1AC",      "#Phi_{1AC}",                                         "(rad.)"},
  {"PairPlaneAngle2AC",      "#Phi_{2AC}",                                         "(rad.)"},
  {"PairPlaneAngle3AC",      "#Phi_{3AC}",                                         "(rad.)"},
  {"PairPlaneAngle4AC",      "#Phi_{4AC}",                                         "(rad.)"},
  {"PairPlaneAngle1Ran",     "#Phi_{1Ran}",                                        "(rad.)"},
  {"PairPlaneAngle2Ran",     "#Phi_{2Ran}",                                        "(rad.)"},
  {"PairPlaneAngle3Ran",     "#Phi_{3Ran}",                                        "(rad.)"},
  {"PairPlaneAngle4Ran",     "#Phi_{4Ran}",                                        "(rad.)"},
  {"RandomRP",               "#Phi_{RanRP}",                                       "(rad.)"},
  {"DeltaPhiRandomRP",       "#Delta #Phi_{RanRP}",                                ""},

  {"PairPlaneMagInPro",      "cos(#Phi)"                                           ""},

  {"Cos2PhiCS",              "cos(2#phi_{CS})",                                    ""},
  {"CosTilPhiCS",            "cos(#phi_{CS})",                                     ""},
  {"CosPhiH2",               "cos(2#phi)",                                         ""},
  {"SinPhiH2",               "sin(2#phi)",                                         ""},
  {"DeltaPhiV0ArpH2",        "#phi^{pair}-#Psi^{V0A}",                             ""},
  {"DeltaPhiV0CrpH2",        "#phi^{pair}-#Psi^{V0C}",                             ""},
  {"DeltaPhiV0ACrpH2",       "#phi^{pair}-#Psi^{V0AC}",                            ""},
  {"V0ArpH2FlowV2",          "cos(2(#phi^{pair}-#Psi^{V0A}))",                     ""},
  {"V0CrpH2FlowV2",          "cos(2(#phi^{pair}-#Psi^{V0C}))",                     ""},
  {"V0ACrpH2FlowV2",         "cos(2(#phi^{pair}-#Psi^{V0AC}))",                    ""},
  {"DeltaPhiv0ArpH2",        "#phi^{pair}-#Psi^{V0A}",                             ""},
  {"DeltaPhiv0CrpH2",        "#phi^{pair}-#Psi^{V0C}",                             ""},
  {"DeltaPhiv0ACrpH2",       "#phi^{pair}-#Psi^{V0AC}",                            ""},
  {"DeltaPhiTPCrpH2",        "#phi^{pair}-#Psi^{TPC}",                             ""},
  {"v0ArpH2FlowV2",          "cos(2(#phi^{pair}-#Psi^{V0A}))",                     ""},
  {"v0CrpH2FlowV2",          "cos(2(#phi^{pair}-#Psi^{V0C}))",                     ""},
  {"v0ACrpH2FlowV2",         "cos(2(#phi^{pair}-#Psi^{V0AC}))",                    ""},
  {"TPCrpH2FlowV2",          "cos(2(#phi^{pair}-#Psi^{TPC}))",                     ""},
  {"TPCrpH2FlowV2Sin",       "sin(2(#phi^{pair}-#Psi^{TPC}))",                     ""},
  {"LegDistance",            "d^{legs}",                                           "(cm)"},
  {"LegDistanceXY",          "d^{legs}_{xy}",                                      "(cm)"},
  {"DeltaEta",               "#Delta #eta",                                        ""},
  {"DeltaPhi",               "#Delta #phi",                                        ""},
  {"DeltaPhiChargeOrdered",  "#Delta #phi_{charge ordered}",                       ""},
  {"Merr",                   "m_{inv} error",                                      "(GeV/#it{c}^{2})"},
  {"DCA",                    "#it{dca}",                                           "(cm)"},
  {"PairType",               "PairType",                                           ""},
  {"PseudoProperTime",       "#tau",                                               "(#mus)"}, //TODO: check unit
  {"PseudoProperTimeErr",    "#tau error",                                         "(#mus)"},
  {"PseudoProperTimeResolution", "(#tau-#tau^{MC truth})",                         "(#mus)"},
  {"PseudoProperTimePull",   "#frac{(#tau-#tau^{MC truth})}{#tau error}",          ""},
  {"TRDpidEffPair",          "#varepsilon^{TRD pair}(PID)",                           ""},
  {"MomAsymDau1",            "#it{p}^{leg1}/#it{p}^{pair}",                        ""},
  {"MomAsymDau2",            "#it{p}^{leg2}/#it{p}^{pair}",                        ""},
  {"EffPair",                "A#times#varepsilon",                                 ""},
  {"OneOverEffPair",         "(A#times#varepsilon)^{-1}",                          ""},
  {"OneOverEffPairSq",       "(A#times#varepsilon)^{-2}",                          ""},
  {"Rndm",                   "P",                                                  ""},
  {"Leg1DCAsigXY",          "#it{dca}_{xy}^{leg1}",                               "(#sigma)"},
  {"Leg1DCAabsXY",          "#it{dca}_{xy}^{leg1}",                               "(cm)"},
  {"Leg1resXY",             "#Delta #it{dca}_{xy}^{leg1}",                        ""},
  //pair dca variables
  {"PairDCAsigXY",          "#it{dca}_{xy}",                                      "(#sigma)"},
  {"PairDCAsigZ",           "#it{dca}_{z}",                                       "(#sigma)"},
  {"PairDCAabsXY",          "#it{dca}_{xy}",                                      "(cm)"},
  {"PairDCAabsZ",           "#it{dca}_{z}",                                       "(cm)"},
  {"PairLinDCAsigXY",       "#it{dca}_{xy} (lin.)",                               "(#sigma)"},
  {"PairLinDCAsigZ",        "#it{dca}_{z} (lin.)",                                "(#sigma)"},
  {"PairLinDCAabsXY",       "#it{dca}_{xy} (lin.)",                               "(cm)"},
  {"PairLinDCAabsZ",        "#it{dca}_{z} (lin.)",                                "(cm)"},
  {"Pairs",                  "pairs/event",                                        ""},
  //
  {"Xprim",                  "x_{prim.vtx}",                                       "(cm)"},
  {"Yprim",                  "y_{prim.vtx}",                                       "(cm)"},
  {"Zprim",                  "z_{prim.vtx}",                                       "(cm)"},
  {"XprimMCtruth",           "x_{MCtruth.prim.vtx}",                               "(cm)"},
  {"YprimMCtruth",           "y_{MCtruth.prim.vtx}",                               "(cm)"},
  {"ZprimMCtruth",           "z_{MCtruth.prim.vtx}",                               "(cm)"},
  {"XRes",                   "#Delta x_{prim.vtx}",                                "(cm)"},
  {"YRes",                   "#Delta y_{prim.vtx}",                                "(cm)"},
  {"ZRes",                   "#Delta z_{prim.vtx}",                                "(cm)"},
  {"PhiMaxPt",               "#phi(#it{p}_{T}^{lead})",                            "(rad.)"},
  {"MaxPt",                  "#it{p}_{T}^{lead}",                                  "(GeV/#it{c})"},

  {"v0ArpH2",                "#Psi^{V0A}",                                         ""},
  {"v0CrpH2",                "#Psi^{V0C}",                                         ""},
  {"v0ACrpH2",               "#Psi^{V0AC}",                                        ""},
  {"v0AxH2",                 "Q_{x}^{V0A}",                                        ""},
  {"v0AyH2",                 "Q_{y}^{V0A}",                                        ""},
  {"v0CxH2",                 "Q_{x}^{V0C}",                                        ""},
  {"v0CyH2",                 "Q_{y}^{V0C}",                                        ""},
  {"v0ACxH2",                "Q_{x}^{V0AC}",                                       ""},
  {"v0ACyH2",                "Q_{y}^{V0AC}",                                       ""},
  {"v0AmagH2",               "|#vec{Q}^{V0A}|",                                     ""},
  {"v0CmagH2",               "|#vec{Q}^{V0C}|",                                     ""},
  {"v0ACmagH2",              "|#vec{Q}^{V0AC}|",                                    ""},
  {"v0A0rpH2",               "#Psi^{V0A}_{ring 0}",                                ""},
  {"v0A3rpH2",               "#Psi^{V0A}_{ring 3}",                                ""},
  {"v0C0rpH2",               "#Psi^{V0C}_{ring 0}",                                ""},
  {"v0C3rpH2",               "#Psi^{V0C}_{ring 3}",                                ""},
  {"v0ATPCDiffH2",           "cos(2(#Psi^{V0A}-#Psi^{TPC}))",                      ""},
  {"v0CTPCDiffH2",           "cos(2(#Psi^{V0C}-#Psi^{TPC}))",                      ""},
  {"v0Av0CDiffH2",           "cos(2(#Psi^{V0A}-#Psi^{V0C}))",                      ""},
  {"v0Av0C0DiffH2",          "cos(2(#Psi^{V0A}-#Psi^{V0C}_{ring 0}))",             ""},
  {"v0Av0C3DiffH2",          "cos(2(#Psi^{V0A}-#Psi^{V0C}_{ring 3}))",             ""},
  {"v0Cv0A0DiffH2",          "cos(2(#Psi^{V0C}-#Psi^{V0A}_{ring 0}))",             ""},
  {"v0Cv0A3DiffH2",          "cos(2(#Psi^{V0C}-#Psi^{V0A}_{ring 3}))",             ""},
  {"v0A0v0A3DiffH2",         "cos(2(#Psi^{V0A}_{ring 0}-#Psi^{V0A}_{ring 3}))",    ""},
  {"v0C0v0C3DiffH2",         "cos(2(#Psi^{V0C}_{ring 0}-#Psi^{V0C}_{ring 3}))",    ""},

  {"MultV0A",                "multiplicity V0A",                                   ""},
  {"MultV0C",                "multiplicity V0C",                                   ""},
  {"MultV0",                 "multiplicity V0",                                    ""},
  {"EqMultV0A",              "equalized multiplicity V0A",                         ""},
  {"EqMultV0C",              "equalized multiplicity V0C",                         ""},
  {"EqMultV0",               "equalized multiplicity V0",                          ""},
  {"AdcV0A",                 "AdcV0A",                                             ""},
  {"AdcV0C",                 "AdcV0C",                                             ""},
  {"AdcV0",                  "AdcV0",                                              ""},
  {"VZERO_ch0",   "VZERO_ch0", ""},
  {"VZERO_ch1",   "VZERO_ch1", ""},
  {"VZERO_ch2",   "VZERO_ch2", ""},
  {"VZERO_ch3",   "VZERO_ch3", ""},
  {"VZERO_ch4",   "VZERO_ch4", ""},
  {"VZERO_ch5",   "VZERO_ch5", ""},
  {"VZERO_ch6",   "VZERO_ch6", ""},
  {"VZERO_ch7",   "VZERO_ch7", ""},
  {"VZERO_ch8",   "VZERO_ch8", ""},
  {"VZERO_ch9",   "VZERO_ch9", ""},
  {"VZERO_ch10",   "VZERO_ch10", ""},
  {"VZERO_ch11",   "VZERO_ch11", ""},
  {"VZERO_ch12",   "VZERO_ch12", ""},
  {"VZERO_ch13",   "VZERO_ch13", ""},
  {"VZERO_ch14",   "VZERO_ch14", ""},
  {"VZERO_ch15",   "VZERO_ch15", ""},
  {"VZERO_ch16",   "VZERO_ch16", ""},
  {"VZERO_ch17",   "VZERO_ch17", ""},
  {"VZERO_ch18",   "VZERO_ch18", ""},
  {"VZERO_ch19",   "VZERO_ch19", ""},
  {"VZERO_ch20",   "VZERO_ch20", ""},
  {"VZERO_ch21",   "VZERO_ch21", ""},
  {"VZERO_ch22",   "VZERO_ch22", ""},
  {"VZERO_ch23",   "VZERO_ch23", ""},
  {"VZERO_ch24",   "VZERO_ch24", ""},
  {"VZERO_ch25",   "VZERO_ch25", ""},
  {"VZERO_ch26",   "VZERO_ch26", ""},
  {"VZERO_ch27",   "VZERO_ch27", ""},
  {"VZERO_ch28",   "VZERO_ch28", ""},
  {"VZERO_ch29",   "VZERO_ch29", ""},
  {"VZERO_ch30",   "VZERO_ch30", ""},
  {"VZERO_ch31",   "VZERO_ch31", ""},
  {"VZERO_ch32",   "VZERO_ch32", ""},
  {"VZERO_ch33",   "VZERO_ch33", ""},
  {"VZERO_ch34",   "VZERO_ch34", ""},
  {"VZERO_ch35",   "VZERO_ch35", ""},
  {"VZERO_ch36",   "VZERO_ch36", ""},
  {"VZERO_ch37",   "VZERO_ch37", ""},
  {"VZERO_ch38",   "VZERO_ch38", ""},
  {"VZERO_ch39",   "VZERO_ch39", ""},
  {"VZERO_ch40",   "VZERO_ch40", ""},
  {"VZERO_ch41",   "VZERO_ch41", ""},
  {"VZERO_ch42",   "VZERO_ch42", ""},
  {"VZERO_ch43",   "VZERO_ch43", ""},
  {"VZERO_ch44",   "VZERO_ch44", ""},
  {"VZERO_ch45",   "VZERO_ch45", ""},
  {"VZERO_ch46",   "VZERO_ch46", ""},
  {"VZERO_ch47",   "VZERO_ch47", ""},
  {"VZERO_ch48",   "VZERO_ch48", ""},
  {"VZERO_ch49",   "VZERO_ch49", ""},
  {"VZERO_ch50",   "VZERO_ch50", ""},
  {"VZERO_ch51",   "VZERO_ch51", ""},
  {"VZERO_ch52",   "VZERO_ch52", ""},
  {"VZERO_ch53",   "VZERO_ch53", ""},
  {"VZERO_ch54",   "VZERO_ch54", ""},
  {"VZERO_ch55",   "VZERO_ch55", ""},
  {"VZERO_ch56",   "VZERO_ch56", ""},
  {"VZERO_ch57",   "VZERO_ch57", ""},
  {"VZERO_ch58",   "VZERO_ch58", ""},
  {"VZERO_ch59",   "VZERO_ch59", ""},
  {"VZERO_ch60",   "VZERO_ch60", ""},
  {"VZERO_ch61",   "VZERO_ch61", ""},
  {"VZERO_ch62",   "VZERO_ch62", ""},
  {"VZERO_ch63",   "VZERO_ch63", ""},
  {"V0AxH2",                 "Q_{x}^{V0A}",                                        ""},
  {"V0AyH2",                 "Q_{y}^{V0A}",                                        ""},
  {"V0ArpH2",                "#Psi^{V0A}",                                         ""},
  {"V0CxH2",                 "Q_{x}^{V0C}",                                        ""},
  {"V0CyH2",                 "Q_{y}^{V0C}",                                        ""},
  {"V0CrpH2",                "#Psi^{V0C}",                                         ""},
  {"V0ACxH2",                "Q_{x}^{V0AC}",                                       ""},
  {"V0ACyH2",                "Q_{y}^{V0AC}",                                       ""},
  {"V0ACrpH2",               "#Psi^{V0AC}",                                        ""},

  {"V0ArpResH2",                   "V0ArpResH2", ""},
  {"V0CrpResH2",                   "V0CrpResH2", ""},
  {"V0ACrpResH2",                   "V0ACrpResH2", ""},
  {"V0XaXcH2",               "Q_{x}^{V0A}#timesQ_{x}^{V0C}",                       ""},
  {"V0XaYaH2",               "Q_{x}^{V0A}#timesQ_{y}^{V0A}",                       ""},
  {"V0XaYcH2",               "Q_{x}^{V0A}#timesQ_{y}^{V0C}",                       ""},
  {"V0YaXcH2",               "Q_{y}^{V0A}#timesQ_{x}^{V0C}",                       ""},
  {"V0YaYcH2",               "Q_{y}^{V0A}#timesQ_{Y}^{V0C}",                       ""},
  {"V0XcYcH2",               "Q_{X}^{V0C}#timesQ_{Y}^{V0C}",                       ""},
  {"V0ATPCDiffH2",           "cos(2(#Psi^{V0A}-#Psi^{TPC}))",                      ""},
  {"V0CTPCDiffH2",           "cos(2(#Psi^{V0C}-#Psi^{TPC}))",                      ""},
  {"V0AV0CDiffH2",           "cos(2(#Psi^{V0A}-#Psi^{V0C}))",                      ""},
  {"TPCxH2",                 "Q_{x}^{TPC}",                                        ""},
  {"TPCyH2",                 "Q_{y}^{TPC}",                                        ""},
  {"TPCmagH2",               "|#vec{Q}^{TPC}|",                                    ""},
  {"TPCrpH2",                "#Psi^{TPC}",                                         ""},
  {"CosTPCrpH2",             "cos(2#Psi^{TPC})",                                   ""},
  {"SinTPCrpH2",             "sin(2#Psi^{TPC})",                                   ""},
  {"TPCsub1xH2",             "Q_{x}^{TPCsub1}",                                    ""},
  {"TPCsub1yH2",             "Q_{y}^{TPCsub1}",                                    ""},
  {"TPCsub1rpH2",            "#Psi^{TPCsub1}",                                     ""},
  {"TPCsub2xH2",             "Q_{x}^{TPCsub2}",                                    ""},
  {"TPCsub2yH2",             "Q_{y}^{TPCsub2}",                                    ""},
  {"TPCsub2rpH2",            "#Psi^{TPCsub2}",                                     ""},
  {"TPCsub12DiffH2",         "cos(2(#Psi^{TPCsub1}-#Psi^{TPCsub2}))",              ""},
  {"TPCsub12DiffH2Sin",      "sin(2(#Psi^{TPCsub1}-#Psi^{TPCsub2}))",              ""},

  {"TPCxH2uc",               "Q_{x}^{TPC} (uncorr.)",                              ""},
  {"TPCyH2uc",               "Q_{y}^{TPC} (uncorr.)",                              ""},
  {"TPCmagH2uc",             "|#vec{Q}^{TPC}| (uncorr.)",                         ""},
  {"TPCrpH2uc",              "#Psi^{TPC} (uncorr.)",                               ""},
  {"TPCsub1xH2uc",           "Q_{x}^{TPCsub1} (uncorr.)",                          ""},
  {"TPCsub1yH2uc",           "Q_{y}^{TPCsub1} (uncorr.)",                          ""},
  {"TPCsub1rpH2uc",          "#Psi^{TPCsub1} (uncorr.)",                           ""},
  {"TPCsub2xH2uc",           "Q_{x}^{TPCsub2} (uncorr.)",                          ""},
  {"TPCsub2yH2uc",           "Q_{y}^{TPCsub2} (uncorr.)",                          ""},
  {"TPCsub2rpH2uc",          "#Psi^{TPCsub2} (uncorr.)",                           ""},
  {"TPCsub12DiffH2uc",       "cos(2(#Psi^{TPCsub1}-#Psi^{TPCsub2})) (uncorr.)",    ""},

  {"ZDCArpH1",               "#Psi^{ZDCA}"                                         ""},
  {"ZDCCrpH1",               "#Psi^{ZDCA}"                                         ""},
  {"ZDCACrpH1",              "#Psi^{ZDCA}"                                         ""},
  {"kZDCrpResH1",            ""                                                    ""},
  {"kv0ZDCrpRes",            ""                                                    ""},
  // Beginning of Eventplane variables from Qn Framework
  // Eventplanes for 2nd harmonic from QnCorrections framework est. 2016
  {"QnTPCrpH2",                "#Psi^{TPC}",                                ""},
  {"QnTPCxH2",                 "Q_{x}^{TPC}",                               ""},
  {"QnTPCyH2",                 "Q_{y}^{TPC}",                               ""},
  {"QnTPCaSiderpH2",                "#Psi^{TPC_{A-Side}}",                                ""},
  {"QnTPCaSidexH2",                 "Q_{x}^{TPC_{A-Side}}",                               ""},
  {"QnTPCaSideyH2",                 "Q_{y}^{TPC_{A-Side}}",                               ""},
  {"QnTPCcSiderpH2",                "#Psi^{TPC_{C-Side}}",                                ""},
  {"QnTPCcSidexH2",                 "Q_{x}^{TPC_{C-Side}}",                               ""},
  {"QnTPCcSideyH2",                 "Q_{y}^{TPC_{C-Side}}",                               ""},
  {"QnV0ArpH2",                "#Psi^{V0A}",                                ""},
  {"QnV0AxH2",                 "Q_{x}^{V0A}",                               ""},
  {"QnV0AyH2",                 "Q_{y}^{V0A}",                               ""},
  {"QnV0CrpH2",                "#Psi^{V0C}",                                ""},
  {"QnV0CxH2",                 "Q_{x}^{V0C}",                               ""},
  {"QnV0CyH2",                 "Q_{y}^{V0C}",                               ""},
  {"QnV0rpH2",                "#Psi^{V0}",                                 ""},
  {"QnV0xH2",                 "Q_{x}^{V0}",                                ""},
  {"QnV0yH2",                 "Q_{y}^{V0}",                                ""},
  {"QnSPDrpH2",                "#Psi^{SPD}",                                ""},
  {"QnSPDxH2",                 "Q_{x}^{SPD}",                               ""},
  {"QnSPDyH2",                 "Q_{y}^{SPD}",                               ""},
  {"QnFMDArpH2",               "#Psi^{FMDA}",                               ""},
  {"QnFMDAxH2",                "Q_{x}^{FMDA}",                              ""},
  {"QnFMDAyH2",                "Q_{y}^{FMDA}",                              ""},
  {"QnFMDCrpH2",               "#Psi^{FMDC}",                               ""},
  {"QnFMDCxH2",                "Q_{x}^{FMDC}",                              ""},
  {"QnFMDCyH2",                "Q_{y}^{FMDC}",                              ""},

  {"QnTPC_V0ADiffH2",          "cos(2(#Psi^{TPC}-#Psi^{V0A}))",             ""},
  {"QnTPC_V0CDiffH2",          "cos(2(#Psi^{TPC}-#Psi^{V0C}))",             ""},
  {"QnTPC_SPDDiffH2",          "cos(2(#Psi^{TPC}-#Psi^{SPD}))",             ""},
  {"QnTPC_FMDADiffH2",         "cos(2(#Psi^{TPC}-#Psi^{FMDA}))",            ""},
  {"QnTPC_FMDCDiffH2",         "cos(2(#Psi^{TPC}-#Psi^{FMDC}))",            ""},
  {"QnTPCaSide_V0DiffH2",      "cos(2(#Psi^{TPC_{A-Side}}-#Psi^{V0}))",     ""},
  {"QnTPCaSide_V0ADiffH2",     "cos(2(#Psi^{TPC_{A-Side}}-#Psi^{V0A}))",    ""},
  {"QnTPCaSide_V0CDiffH2",     "cos(2(#Psi^{TPC_{A-Side}}-#Psi^{V0C}))",    ""},
  {"QnTPCaSide_TPCcSideDiffH2", "cos(2(#Psi^{TPC_{A-Side}}-#Psi^{TPC_{C-Side}}))",             ""},
  {"QnTPCcSide_V0DiffH2",      "cos(2(#Psi^{TPC_{C-Side}}-#Psi^{V0}))",     ""},
  {"QnTPCcSide_V0ADiffH2",     "cos(2(#Psi^{TPC_{C-Side}}-#Psi^{V0A}))",    ""},
  {"QnTPCcSide_V0CDiffH2",     "cos(2(#Psi^{TPC_{C-Side}}-#Psi^{V0C}))",    ""},
  {"QnV0A_V0CDiffH2",          "cos(2(#Psi^{V0A}-#Psi^{V0C}))",             ""},
  {"QnV0A_SPDDiffH2",          "cos(2(#Psi^{V0A}-#Psi^{SPD}))",             ""},
  {"QnV0A_FMDADiffH2",         "cos(2(#Psi^{V0A}-#Psi^{FMDA}))",            ""},
  {"QnV0A_FMDCDiffH2",         "cos(2(#Psi^{V0A}-#Psi^{FMDC}))",            ""},
  {"QnV0C_SPDDiffH2",          "cos(2(#Psi^{V0C}-#Psi^{SPD}))",             ""},
  {"QnV0C_FMDADiffH2",         "cos(2(#Psi^{V0C}-#Psi^{FMDA}))",            ""},
  {"QnV0C_FMDCDiffH2",         "cos(2(#Psi^{V0C}-#Psi^{FMDC}))",            ""},
  {"QnSPD_FMDADiffH2",         "cos(2(#Psi^{SPD}-#Psi^{FMDA}))",            ""},
  {"QnSPD_FMDCDiffH2",         "cos(2(#Psi^{SPD}-#Psi^{FMDC}))",            ""},
  {"QnFMDA_FMDCDiffH2",        "cos(2(#Psi^{FMDA}-#Psi^{FMDC}))",           ""},
// XX YY XY YX relations
  {"QnTPCx_V0AxCorrH2",         "TPC^{Qn}_{x} * V0A^{Qn}_{x}",                 ""},
  {"QnTPCx_V0AyCorrH2",         "TPC^{Qn}_{x} * V0A^{Qn}_{y}",                 ""},
  {"QnTPCy_V0AxCorrH2",         "TPC^{Qn}_{y} * V0A^{Qn}_{x}",                 ""},
  {"QnTPCy_V0AyCorrH2",         "TPC^{Qn}_{y} * V0A^{Qn}_{y}",                 ""},
  {"QnTPCx_V0CxCorrH2",         "TPC^{Qn}_{x} * V0C^{Qn}_{x}",                 ""},
  {"QnTPCx_V0CyCorrH2",         "TPC^{Qn}_{x} * V0C^{Qn}_{y}",                 ""},
  {"QnTPCy_V0CxCorrH2",         "TPC^{Qn}_{y} * V0C^{Qn}_{x}",                 ""},
  {"QnTPCy_V0CyCorrH2",         "TPC^{Qn}_{y} * V0C^{Qn}_{y}",                 ""},
  {"QnTPCx_SPDxCorrH2",         "TPC^{Qn}_{x} * SPD^{Qn}_{x}",                 ""},
  {"QnTPCx_SPDyCorrH2",         "TPC^{Qn}_{x} * SPD^{Qn}_{y}",                 ""},
  {"QnTPCy_SPDxCorrH2",         "TPC^{Qn}_{y} * SPD^{Qn}_{x}",                 ""},
  {"QnTPCy_SPDyCorrH2",         "TPC^{Qn}_{y} * SPD^{Qn}_{y}",                 ""},
  {"QnTPCx_FMDAxCorrH2",        "TPC^{Qn}_{x} * FMDA^{Qn}_{x}",                ""},
  {"QnTPCx_FMDAyCorrH2",        "TPC^{Qn}_{x} * FMDA^{Qn}_{y}",                ""},
  {"QnTPCy_FMDAxCorrH2",        "TPC^{Qn}_{y} * FMDA^{Qn}_{x}",                ""},
  {"QnTPCy_FMDAyCorrH2",        "TPC^{Qn}_{y} * FMDA^{Qn}_{y}",                ""},
  {"QnTPCx_FMDCxCorrH2",        "TPC^{Qn}_{x} * FMDC^{Qn}_{x}",                ""},
  {"QnTPCx_FMDCyCorrH2",        "TPC^{Qn}_{x} * FMDC^{Qn}_{y}",                ""},
  {"QnTPCy_FMDCxCorrH2",        "TPC^{Qn}_{y} * FMDC^{Qn}_{x}",                ""},
  {"QnTPCy_FMDCyCorrH2",        "TPC^{Qn}_{y} * FMDC^{Qn}_{y}",                ""},
  {"QnV0Ax_V0CxCorrH2",         "V0A^{Qn}_{x} * V0C^{Qn}_{x}",                 ""},
  {"QnV0Ax_V0CyCorrH2",         "V0A^{Qn}_{x} * V0C^{Qn}_{y}",                 ""},
  {"QnV0Ay_V0CxCorrH2",         "V0A^{Qn}_{y} * V0C^{Qn}_{x}",                 ""},
  {"QnV0Ay_V0CyCorrH2",         "V0A^{Qn}_{y} * V0C^{Qn}_{y}",                 ""},
  {"QnV0Ax_SPDxCorrH2",         "V0A^{Qn}_{x} * SPD^{Qn}_{x}",                 ""},
  {"QnV0Ax_SPDyCorrH2",         "V0A^{Qn}_{x} * SPD^{Qn}_{y}",                 ""},
  {"QnV0Ay_SPDxCorrH2",         "V0A^{Qn}_{y} * SPD^{Qn}_{x}",                 ""},
  {"QnV0Ay_SPDyCorrH2",         "V0A^{Qn}_{y} * SPD^{Qn}_{y}",                 ""},
  {"QnV0Ax_FMDAxCorrH2",        "V0A^{Qn}_{x} * FMDA^{Qn}_{x}",                ""},
  {"QnV0Ax_FMDAyCorrH2",        "V0A^{Qn}_{x} * FMDA^{Qn}_{y}",                ""},
  {"QnV0Ay_FMDAxCorrH2",        "V0A^{Qn}_{y} * FMDA^{Qn}_{x}",                ""},
  {"QnV0Ay_FMDAyCorrH2",        "V0A^{Qn}_{y} * FMDA^{Qn}_{y}",                ""},
  {"QnV0Ax_FMDCxCorrH2",        "V0A^{Qn}_{x} * FMDC^{Qn}_{x}",                ""},
  {"QnV0Ax_FMDCyCorrH2",        "V0A^{Qn}_{x} * FMDC^{Qn}_{y}",                ""},
  {"QnV0Ay_FMDCxCorrH2",        "V0A^{Qn}_{y} * FMDC^{Qn}_{x}",                ""},
  {"QnV0Ay_FMDCyCorrH2",        "V0A^{Qn}_{y} * FMDC^{Qn}_{y}",                ""},
  {"QnV0Cx_FMDAxCorrH2",        "V0C^{Qn}_{x} * FMDA^{Qn}_{x}",                ""},
  {"QnV0Cx_FMDAyCorrH2",        "V0C^{Qn}_{x} * FMDA^{Qn}_{y}",                ""},
  {"QnV0Cy_FMDAxCorrH2",        "V0C^{Qn}_{y} * FMDA^{Qn}_{x}",                ""},
  {"QnV0Cy_FMDAyCorrH2",        "V0C^{Qn}_{y} * FMDA^{Qn}_{y}",                ""},
  {"QnV0Cx_FMDCxCorrH2",        "V0C^{Qn}_{x} * FMDC^{Qn}_{x}",                ""},
  {"QnV0Cx_FMDCyCorrH2",        "V0C^{Qn}_{x} * FMDC^{Qn}_{y}",                ""},
  {"QnV0Cy_FMDCxCorrH2",        "V0C^{Qn}_{y} * FMDC^{Qn}_{x}",                ""},
  {"QnV0Cy_FMDCyCorrH2",        "V0C^{Qn}_{y} * FMDC^{Qn}_{y}",                ""},
  {"QnV0Cx_SPDxCorrH2",         "V0C^{Qn}_{x} * SPD^{Qn}_{x}",                 ""},
  {"QnV0Cx_SPDyCorrH2",         "V0C^{Qn}_{x} * SPD^{Qn}_{y}",                 ""},
  {"QnV0Cy_SPDxCorrH2",         "V0C^{Qn}_{y} * SPD^{Qn}_{x}",                 ""},
  {"QnV0Cy_SPDyCorrH2",         "V0C^{Qn}_{y} * SPD^{Qn}_{y}",                 ""},
  {"QnSPDx_FMDAxCorrH2",        "SPD^{Qn}_{x} * FMDA^{Qn}_{x}",                ""},
  {"QnSPDx_FMDAyCorrH2",        "SPD^{Qn}_{x} * FMDA^{Qn}_{y}",                ""},
  {"QnSPDy_FMDAxCorrH2",        "SPD^{Qn}_{y} * FMDA^{Qn}_{x}",                ""},
  {"QnSPDy_FMDAyCorrH2",        "SPD^{Qn}_{y} * FMDA^{Qn}_{y}",                ""},
  {"QnSPDx_FMDCxCorrH2",        "SPD^{Qn}_{x} * FMDC^{Qn}_{x}",                ""},
  {"QnSPDx_FMDCyCorrH2",        "SPD^{Qn}_{x} * FMDC^{Qn}_{y}",                ""},
  {"QnSPDy_FMDCxCorrH2",        "SPD^{Qn}_{y} * FMDC^{Qn}_{x}",                ""},
  {"QnSPDy_FMDCyCorrH2",        "SPD^{Qn}_{y} * FMDC^{Qn}_{y}",                ""},
  {"QnFMDAx_FMDCxCorrH2",       "FMDA^{Qn}_{x} * FMDC^{Qn}_{x}",               ""},
  {"QnFMDAx_FMDCyCorrH2",       "FMDA^{Qn}_{x} * FMDC^{Qn}_{y}",               ""},
  {"QnFMDAy_FMDCxCorrH2",       "FMDA^{Qn}_{y} * FMDC^{Qn}_{x}",               ""},
  {"QnFMDAy_FMDCyCorrH2",       "FMDA^{Qn}_{y} * FMDC^{Qn}_{y}",               ""},
//
  {"QnDeltaPhiTPCrpH2",         "#phi^{pair}-#Psi^{TPC}",                     ""},
  {"QnDeltaPhiTrackTPCrpH2",    "#phi^{Track}-#Psi^{TPC}",                     ""},
  {"QnDeltaPhiV0ArpH2",         "#phi^{pair}-#Psi^{V0A}",                     ""},
  {"QnDeltaPhiV0CrpH2",         "#phi^{pair}-#Psi^{V0C}",                     ""},
  {"QnDeltaPhiV0rpH2",          "#phi^{pair}-#Psi^{V0}",                      ""},
  {"QnDeltaPhiSPDrpH2",         "#phi^{pair}-#Psi^{SPD}",                     ""},
  {"QnTPCrpH2FlowV2",           "cos(2(#phi^{pair}-#Psi^{TPC}))",             ""},
  {"QnV0ArpH2FlowV2",           "cos(2(#phi^{pair}-#Psi^{V0A}))",             ""},
  {"QnV0CrpH2FlowV2",           "cos(2(#phi^{pair}-#Psi^{V0C}))",             ""},
  {"QnV0rpH2FlowV2",            "cos(2(#phi^{pair}-#Psi^{V0}))",              ""},
  {"QnSPDrpH2FlowV2",           "cos(2(#phi^{pair}-#Psi^{SPD}))",             ""},
  // End of Eventplane variables from Qn Framework

  {"NTrk",                   "N_{trk}",                                            ""},
  {"Tracks",                 "tracks/per event",                                             ""},
  {"NVtxContrib",            "N_{vtx. contrib.}",                                  ""},
  {"NVtxContribTPC",         "N_{vtx. contrib.}^{TPC}",                            ""},
  {"Nacc",                   "N_{acc} #cbar_{#||{#eta}<0.9}",                      ""},
  {"MatchEffITSTPC",        "N_{trk}^{TPC}/N_{trk}^{ITS} #cbar_{#||{#eta}<0.9}",  ""},
  {"MatchEffITSTPCinPlane", "N_{trk}^{TPC}/N_{trk}^{ITS} #cbar_{#||{#eta}<0.9, #||{#phi^{Track}-#Psi^{TPC}} < #frac{#pi}{4} || #frac{3#pi}{4} < #||{#phi^{Track}-#Psi^{TPC}} < #pi}",  ""},
  {"MatchEffITSTPCoutPlane", "N_{trk}^{TPC}/N_{trk}^{ITS} #cbar_{#||{#eta}<0.9, #frac{#pi}{4} < #||{#phi^{Track}-#Psi^{TPC}} < #frac{3#pi}{4}}",  ""},
  {"NaccTrcklts",            "N_{acc. trkl} #cbar_{#||{#eta}<1.6}",                ""},
  {"NaccTrcklts09",          "N_{acc. trkl} #cbar_{#||{#eta}<0.9}",                ""},
  {"NaccTrcklts10",          "N_{acc. trkl} #cbar_{#||{#eta}<1.0}",                ""},
  {"NaccTrcklts0916",        "N_{acc. trkl} #cbar_{-1.6<#eta<-0.9}^{0.9<#eta<1.6}",""},
  {"NaccTrckltsCorr",        "N_{acc. trkl}^{corr.} #cbar_{#||{#eta}<1.6}",        ""},
  {"NaccTrcklts10Corr",      "N_{acc. trkl}^{corr.} #cbar_{#||{#eta}<1.0}",        ""},
  {"NaccTrcklts0916Corr",    "N_{acc. trkl}^{corr.} #cbar_{-1.6<#eta<-0.9}^{0.9<#eta<1.6}",""},

  {"NaccTrckltsEsd05",       "N_{acc. trkl} #cbar_{#||{#eta}<0.5} (SPD)",          ""},
  {"NaccTrckltsEsd10",       "N_{acc. trkl} #cbar_{#||{#eta}<1.0} (SPD)",          ""},
  {"NaccTrckltsEsd16",       "N_{acc. trkl} #cbar_{#||{#eta}<1.6} (SPD)",          ""},
  {"NaccTrckltsEsd05Corr",   "N_{acc. trkl}^{corr.} #cbar_{#||{#eta}<0.5} (SPD)",  ""},
  {"NaccTrckltsEsd10Corr",   "N_{acc. trkl}^{corr.} #cbar_{#||{#eta}<1.0} (SPD)",  ""},
  {"NaccTrckltsEsd16Corr",   "N_{acc. trkl}^{corr.} #cbar_{#||{#eta}<1.6} (SPD)",  ""},
  {"NaccItsTpcEsd05",        "N_{acc. trkl} #cbar_{#||{#eta}<0.5} (ITS+TPC)",      ""},
  {"NaccItsTpcEsd10",        "N_{acc. trkl} #cbar_{#||{#eta}<1.0} (ITS+TPC)",      ""},
  {"NaccItsTpcEsd16",        "N_{acc. trkl} #cbar_{#||{#eta}<1.6} (ITS+TPC)",      ""},
  {"NaccItsTpcEsd05Corr",    "N_{acc. trkl}^{corr.} #cbar_{#||{#eta}<0.5} (ITS+TPC)", ""},
  {"NaccItsTpcEsd10Corr",    "N_{acc. trkl}^{corr.} #cbar_{#||{#eta}<1.0} (ITS+TPC)", ""},
  {"NaccItsTpcEsd16Corr",    "N_{acc. trkl}^{corr.} #cbar_{#||{#eta}<1.6} (ITS+TPC)", ""},
  {"NaccItsPureEsd05",       "N_{acc. trkl} #cbar_{#||{#eta}<0.5} (ITS)",          ""},
  {"NaccItsPureEsd10",       "N_{acc. trkl} #cbar_{#||{#eta}<1.0} (ITS)",          ""},
  {"NaccItsPureEsd16",       "N_{acc. trkl} #cbar_{#||{#eta}<1.6} (ITS)",          ""},
  {"NaccItsPureEsd05Corr",   "N_{acc. trkl}^{corr.} #cbar_{#||{#eta}<0.5} (ITS)",  ""},
  {"NaccItsPureEsd10Corr",   "N_{acc. trkl}^{corr.} #cbar_{#||{#eta}<1.0} (ITS)",  ""},
  {"NaccItsPureEsd16Corr",   "N_{acc. trkl}^{corr.} #cbar_{#||{#eta}<1.6} (ITS)",  ""},

  {"RefMult",                "N_{trk}^{ref}",                                      ""},
  {"RefMultTPConly",         "N_{trk}^{TPConly}",                                  ""},
  {"NTPCtrkswITSout",         "N_{trk}^{TPCwITSout}",                               ""},
  {"NTPCclsEvent",    "N_{cls}^{TPC}",                                      ""},
  {"RefMultOvRefMultTPConly",  "N_{trk}/N_{trk}^{TPConly}",                        ""},
  {"Nch",                    "N_{ch} #cbar_{#||{#eta}<1.6}",                       ""},
  {"NchJpsiExcl",           "N_{ch} #cbar_{#||{#eta}<1.6} without J/#psi daughters ",  ""},
  {"Nch05",                  "N_{ch} #cbar_{#||{#eta}<0.5}",                       ""},
  {"Nch05JpsiExcl",           "N_{ch} #cbar_{#||{#eta}<0.5} without J/#psi daughters ",  ""},
  {"Nch10",                  "N_{ch} #cbar_{#||{#eta}<1.0}",                       ""},
  {"Nch10JpsiExcl",           "N_{ch} #cbar_{#||{#eta}<1.0} without J/#psi daughters ",  ""},
  {"Centrality",             "centrality_{V0M}",                                  "(%)"},
  {"CentralityV0A",          "centrality_{V0A}",                                  "(%)"},
  {"CentralityV0C",          "centrality_{V0C}",                                  "(%)"},
  {"CentralityZNA",          "centrality_{ZNA}",                                  "(%)"},
  {"CentralitySPD",          "centrality_{SPD}",                                  "(%)"},
  {"CentralityNew",             "centrality_{V0M}",                               "(%)"},
  {"CentralityCL0",              "centrality_{CL0}",                              "(%)"},
  {"CentralityCL1",              "centrality_{CL1}",                              "(%)"},
  {"CentralitySPDClusters",      "centrality_{SPDclusters}",                      "(%)"},
  {"CentralitySPDTracklets",     "centrality_{SPDTracklets}",                     "(%)"},
  {"CentralityCL0plus05",        "centrality_{V0M AP +0.5%}",                     "(%)"},
  {"CentralityCL0minus05",       "centrality_{V0M AP -0.5%}",                     "(%)"},
  {"CentralityCL0plus10",        "centrality_{V0M AP +1.0%}",                     "(%)"},
  {"CentralityCL0minus10",       "centrality_{V0M AP -1.0%}",                     "(%)"},



  {"TriggerInclONL",         "online trigger bit (inclusive)",                     ""},
  {"TriggerInclOFF",         "offline trigger bit (inclusive)",                    ""},
  {"TriggerExclOFF",         "offline trigger bit (exclusive)",                    ""},
  {"Nevents",                "N_{evt}",                                            ""},
  {"RunNumber",              "run",                                                ""},
  {"MixingBin",              "mixing bin",                                         ""},
  {"LegSource",              "Leg source",                                         ""}
};

AliPIDResponse* AliDielectronVarManager::fgPIDResponse      = 0x0;
AliVEvent*      AliDielectronVarManager::fgEvent            = 0x0;
AliEventplane*  AliDielectronVarManager::fgTPCEventPlane    = 0x0;
AliKFVertex*    AliDielectronVarManager::fgKFVertex         = 0x0;
TProfile*       AliDielectronVarManager::fgMultEstimatorAvg[7][9] = {{0x0}};
TH3D*           AliDielectronVarManager::fgTRDpidEff[10][4] = {{0x0}};
TObject*        AliDielectronVarManager::fgLegEffMap           = 0x0;
TObject*        AliDielectronVarManager::fgPairEffMap          = 0x0;
TBits*          AliDielectronVarManager::fgFillMap          = 0x0;
Double_t        AliDielectronVarManager::fgTRDpidEffCentRanges[10][4] = {{0.0}};
TString         AliDielectronVarManager::fgVZEROCalibrationFile = "";
TString         AliDielectronVarManager::fgVZERORecenteringFile = "";
TString         AliDielectronVarManager::fgZDCRecenteringFile = "";
TProfile2D*     AliDielectronVarManager::fgVZEROCalib[64] = {0x0};
TProfile2D*     AliDielectronVarManager::fgVZERORecentering[2][2] = {{0x0,0x0},{0x0,0x0}};
TProfile3D*     AliDielectronVarManager::fgZDCRecentering[3][2] = {{0x0,0x0},{0x0,0x0},{0x0,0x0}};
AliDielectronQnEPcorrection* AliDielectronVarManager::fgQnEPacRemoval = 0x0;
Bool_t          AliDielectronVarManager::fgEventPlaneACremoval = kFALSE;
TString         AliDielectronVarManager::fgQnVectorNorm = "";
Int_t           AliDielectronVarManager::fgCurrentRun = -1;
Double_t        AliDielectronVarManager::fgData[AliDielectronVarManager::kNMaxValues] = {0.};
//________________________________________________________________
AliDielectronVarManager::AliDielectronVarManager() :
  TNamed("AliDielectronVarManager","AliDielectronVarManager")
{
  //
  // Default constructor
  //
  for(Int_t i=0; i<6; ++i)
    for(Int_t j=0; j<9; ++j)
      fgMultEstimatorAvg[i][j] = 0x0;
  for(Int_t i=0; i<10; ++i)
    for(Int_t j=0; j<4; ++j)
      fgTRDpidEff[i][j] = 0x0;
  for(Int_t i=0; i<64; ++i) fgVZEROCalib[i] = 0x0;
  for(Int_t i=0; i<2; ++i) {
    for(Int_t j=0; j<2; ++j) fgVZERORecentering[i][j] = 0x0;
  }
  for(Int_t i=0; i<3; ++i)
    for(Int_t j=0; j<2; ++j) fgZDCRecentering[i][j] = 0x0;

  gRandom->SetSeed();
}

//________________________________________________________________
AliDielectronVarManager::AliDielectronVarManager(const char* name, const char* title) :
  TNamed(name,title)
{
  //
  // Named constructor
  //
  for(Int_t i=0; i<6; ++i)
    for(Int_t j=0; j<9; ++j)
      fgMultEstimatorAvg[i][j] = 0x0;
  for(Int_t i=0; i<10; ++i)
    for(Int_t j=0; j<4; ++j)
      fgTRDpidEff[i][j] = 0x0;
  for(Int_t i=0; i<64; ++i) fgVZEROCalib[i] = 0x0;
  for(Int_t i=0; i<2; ++i)
    for(Int_t j=0; j<2; ++j)
      fgVZERORecentering[i][j] = 0x0;
  for(Int_t i=0; i<3; ++i)
    for(Int_t j=0; j<2; ++j) fgZDCRecentering[i][j] = 0x0;


  gRandom->SetSeed();
}

//________________________________________________________________
AliDielectronVarManager::~AliDielectronVarManager()
{
  //
  // Default destructor
  //
  for(Int_t i=0; i<6; ++i)
    for(Int_t j=0; j<9; ++j)
      if(fgMultEstimatorAvg[i][j]) delete fgMultEstimatorAvg[i][j];
  for(Int_t i=0; i<10; ++i)
    for(Int_t j=0; j<4; ++j)
      if(fgTRDpidEff[i][j]) delete fgTRDpidEff[i][j];
  for(Int_t i=0; i<64; ++i)
    if(fgVZEROCalib[i]) delete fgVZEROCalib[i];
  for(Int_t i=0; i<2; ++i)
    for(Int_t j=0; j<2; ++j)
      if(fgVZERORecentering[i][j]) delete fgVZERORecentering[i][j];
  for(Int_t i=0; i<3; ++i)
    for(Int_t j=0; j<2; ++j)
      if(fgZDCRecentering[i][j]) delete fgZDCRecentering[i][j];

}

//________________________________________________________________
UInt_t AliDielectronVarManager::GetValueType(const char* valname) {
  //
  // Get value type by value name
  //

  TString name(valname);
  for(UInt_t i=0; i<AliDielectronVarManager::kNMaxValues; i++) {
    if(!name.CompareTo(fgkParticleNames[i][0])) return i;
  }
  return -1;
}
