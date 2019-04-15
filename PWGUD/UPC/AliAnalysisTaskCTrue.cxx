/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* AliAnalysisTaskCTrue
 *
 * simple task to study CTRUE events to estimate number of UPC
 * events lost due to vetoes.
 *
 * This analysis task is greatly inspired from the
 * AliAnalysisTaskTransTask.
 * Basically it is a merge of that analysis task with the macros
 * that had been passed onto me.
 * The only thing that changes a bit is that this analysis task should
 * ultimately lead to a simultaneous analysis of all the possible classes
 * for the CTRUE events, be they B, E, or A and C.
 * The basic idea is to have a switch and go to fill the respective
 * histograms. This way it is easier for book keeping.
 */


// c++ headers
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
// #include <pair>
#include <algorithm>
#include <stdlib.h>

// root headers
#include <TMath.h>
#include "TH1I.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH1.h"
#include <TFile.h>
#include <TF2.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TLegend.h>
// #include <TGraphErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraph2D.h>
#include <TStopwatch.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>
#include <TLatex.h>
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TParticle.h"
#include "TObjString.h"
#include "TList.h"
#include "TChain.h"


// aliroot headers
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMuonTrackCuts.h"
#include "AliAODVertex.h"
#include "AliAnalysisTaskCTrue.h"

class AliAnalysisTaskCTrue;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskCTrue) // classimp: necessary for root

AliAnalysisTaskCTrue::AliAnalysisTaskCTrue()
    : AliAnalysisTaskSE(),
      fAOD(0),
      fOutputList(0),
      fAnaTree(0),
      fRunNum(0),
      fTracklets(0),
      fCtrue(0),
      fL0inputs(0),
      fL1inputs(0),
      fZem1Energy(0),
      fZem2Energy(0),
      fZNCEnergy(0),
      fZNAEnergy(0),
      fZPCEnergy(0),
      fZPAEnergy(0),
      fZNATime(0),
      fZNCTime(0),
      fZNATDC{0, 0, 0, 0},
      fZNCTDC{0, 0, 0, 0},
      fZPATDC{0, 0, 0, 0},
      fZPCTDC{0, 0, 0, 0},
      fV0ADecision(-10),
      fV0CDecision(-10),
      fADADecision(-10),
      fADCDecision(-10),
      fIR1Map(0),
      fIR2Map(0),
      fBCrossNum(0),
      fCounterH(0),
      fCounterTryingH(0),
      // trigger inputs
      inputId_0VBA(0),
      inputId_0VBC(0),
      inputId_0UBA(0),
      inputId_0UBC(0),
      inputId_0SH1(0),
      inputId_0STG(0),
      inputId_1ZED(0),
      inputId_0MUL(0),
      inputId_0OM2(0),
      inputId_0VOM(0),
      TotalWeightForEfficiency(0),
      fVectorEfficiency(0),
      // plots for CTRUE
      fEntriesAgainstRunNumberProperlyH(0),
      fCTrueBEventsPerRunNumberH(0),
      fCTrueBEventsPerRunNumberConditionsH(0),
      fTrackletsPerRunNumberH(0),
      fVBAforRunNumberH(0),
      fVBATrackletsForRunNumberH(0),
      fVDAforRunNumberH(0),
      fVDATrackletsForRunNumberH(0),
      fVBAandVDAforRunNumberH(0),
      fVBAandVDATrackletsForRunNumberH(0),
      fVBCforRunNumberH(0),
      fVBCTrackletsForRunNumberH(0),
      fVDCforRunNumberH(0),
      fVDCTrackletsForRunNumberH(0),
      fVBCandVDCforRunNumberH(0),
      fVBCandVDCTrackletsForRunNumberH(0),
      fRecurringVetoAH(0),
      fRecurringVetoAtrackletsH(0),
      fRecurringVetoCH(0),
      fRecurringVetoCtrackletsH(0),
      fCTRUEBfV0CDecisionVSfADCDecisionH(0),
      fCTRUEBfV0ADecisionVSfADADecisionH(0)
      // // result TGraphErrors
      // fOnlyTrackletsGE(0),
      // fVbaGE(0),
      // fVbaTrackletsGE(0),
      // fVdaGE(0),
      // fVdaTrackletsGE(0),
      // fVbaAndVdaGE(0),
      // fVbaAndVdaTrackletsGE(0),
      // fVbcGE(0),
      // fVbcTrackletsGE(0),
      // fVdcGE(0),
      // fVdcTrackletsGE(0),
      // fVbcAndVdcGE(0),
      // fVbcAndVdcTrackletsGE(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCTrue::AliAnalysisTaskCTrue(const char* name)
    : AliAnalysisTaskSE(name),
      fAOD(0),
      fOutputList(0),
      fAnaTree(0),
      fRunNum(0),
      fTracklets(0),
      fCtrue(0),
      fL0inputs(0),
      fL1inputs(0),
      fZem1Energy(0),
      fZem2Energy(0),
      fZNCEnergy(0),
      fZNAEnergy(0),
      fZPCEnergy(0),
      fZPAEnergy(0),
      fZNATime(0),
      fZNCTime(0),
      fZNATDC{0, 0, 0, 0},
      fZNCTDC{0, 0, 0, 0},
      fZPATDC{0, 0, 0, 0},
      fZPCTDC{0, 0, 0, 0},
      fV0ADecision(-10),
      fV0CDecision(-10),
      fADADecision(-10),
      fADCDecision(-10),
      fIR1Map(0),
      fIR2Map(0),
      fBCrossNum(0),
      fCounterH(0),
      fCounterTryingH(0),
      // trigger inputs
      inputId_0VBA(0),
      inputId_0VBC(0),
      inputId_0UBA(0),
      inputId_0UBC(0),
      inputId_0SH1(0),
      inputId_0STG(0),
      inputId_1ZED(0),
      inputId_0MUL(0),
      inputId_0OM2(0),
      inputId_0VOM(0),
      TotalWeightForEfficiency(0),
      fVectorEfficiency(0),
      // plots for CTRUE
      fEntriesAgainstRunNumberProperlyH(0),
      fCTrueBEventsPerRunNumberH(0),
      fCTrueBEventsPerRunNumberConditionsH(0),
      fTrackletsPerRunNumberH(0),
      fVBAforRunNumberH(0),
      fVBATrackletsForRunNumberH(0),
      fVDAforRunNumberH(0),
      fVDATrackletsForRunNumberH(0),
      fVBAandVDAforRunNumberH(0),
      fVBAandVDATrackletsForRunNumberH(0),
      fVBCforRunNumberH(0),
      fVBCTrackletsForRunNumberH(0),
      fVDCforRunNumberH(0),
      fVDCTrackletsForRunNumberH(0),
      fVBCandVDCforRunNumberH(0),
      fVBCandVDCTrackletsForRunNumberH(0),
      fRecurringVetoAH(0),
      fRecurringVetoAtrackletsH(0),
      fRecurringVetoCH(0),
      fRecurringVetoCtrackletsH(0),
      fCTRUEBfV0CDecisionVSfADCDecisionH(0),
      fCTRUEBfV0ADecisionVSfADADecisionH(0)
      // // result TGraphErrors
      // fOnlyTrackletsGE(0),
      // fVbaGE(0),
      // fVbaTrackletsGE(0),
      // fVdaGE(0),
      // fVdaTrackletsGE(0),
      // fVbaAndVdaGE(0),
      // fVbaAndVdaTrackletsGE(0),
      // fVbcGE(0),
      // fVbcTrackletsGE(0),
      // fVdcGE(0),
      // fVdcTrackletsGE(0),
      // fVbcAndVdcGE(0),
      // fVbcAndVdcTrackletsGE(0)
{
    // fMapGoodRunsToMuAndWeight.insert( std::make_pair( 0, std::make_pair(0,0) ) );
    Set2018PbPb();
    FillGoodRunMapInfo(fMapGoodRunsToMuAndWeight);
    fVectorEfficiency = new std::vector<Double_t>[3];

    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TTree::Class());
    DefineOutput(2, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskCTrue::~AliAnalysisTaskCTrue()
{
    // destructor
    if(fOutputList) {delete fOutputList;}
    if(fAnaTree)    {delete fAnaTree;}
}
//_____________________________________________________________________________
void AliAnalysisTaskCTrue::FillGoodRunMapInfo(AliRunWithMuAndWeight &fMapGoodRunsToMuAndWeight)
{
  fMapGoodRunsToMuAndWeight.clear();
  const Int_t nGoodRunsLHC18q = 128;
  Int_t GoodRunsLHC18q[nGoodRunsLHC18q] = {
      295585, 295586, 295587, 295588, 295589, 295612, 295615, 295665, 295666,
      295667, 295668, 295671, 295673, 295675, 295676, 295677, 295714, 295716,
      295717, 295718, 295719, 295723, 295725, 295753, 295754, 295755, 295758,
      295759, 295762, 295763, 295786, 295788, 295791, 295816, 295818, 295819,
      295822, 295825, 295826, 295829, 295831, 295854, 295855, 295856, 295859,
      295860, 295861, 295863, 295881, 295908, 295909, 295910, 295913, 295936,
      295937, 295941, 295942, 295943, 295945, 295947, 296061, 296062, 296063,
      296065, 296066, 296068, 296123, 296128, 296132, 296133, 296134, 296135,
      296142, 296143, 296191, 296192, 296194, 296195, 296196, 296197, 296198,
      296241, 296242, 296243, 296244, 296246, 296247, 296269, 296270, 296273,
      296279, 296280, 296303, 296304, 296307, 296309, 296312, 296377, 296378,
      296379, 296380, 296381, 296383, 296414, 296419, 296420, 296423, 296424,
      296433, 296472, 296509, 296510, 296511, 296514, 296516, 296547, 296548,
      296549, 296550, 296551, 296552, 296553, 296615, 296616, 296618, 296619,
      296622, 296623
  };
  Double_t muLHC18q[nGoodRunsLHC18q] = {
      0.00133815,  0.00119445,  0.00104683,  0.000956058, 0.000804214,
      0.00110139,  0.000901009, 0.00115179,  0.000967337, 0.0005478,
      0.000432373, 0.000888828, 0.000789307, 0.000734837, 0.000683085,
      0.000611429, 0.00146325,  0.00131991,  0.0012197,   0.00112291,
      0.00101555,  0.000819151, 0.000648845, 0.00143418,  0.00143705,
      0.00139216,  0.00102991,  0.000822833, 0.0007025,   0.000645264,
      0.0014374,   0.00126414,  0.000874272, 0.00118069,  0.00117413,
      0.00117183,  0.00103607,  0.000921924, 0.000838568, 0.000740449,
      0.000674052, 0.00108059,  0.0010723,   0.00107055,  0.00102866,
      0.000945887, 0.000830439, 0.000745441, 0.000635393, 0.00108034,
      0.00108607,  0.00108486,  0.00107213,  0.00190871,  0.00183262,
      0.00170843,  0.00159405,  0.00146588,  0.00131433,  0.00116036,
      0.00202716,  0.00202178,  0.00188615,  0.00170681,  0.00160002,
      0.00148213,  0.00110643,  0.00110174,  0.00110478,  0.00110096,
      0.00108669,  0.000872939, 0.000656138, 0.000575556, 0.00108822,
      0.00108259,  0.00109788,  0.00103322,  0.000883841, 0.000702431,
      0.000588796, 0.00105827,  0.00110793,  0.001102,    0.00106385,
      0.000790883, 0.000671304, 0.00111481,  0.00110992,  0.00106459,
      0.000737928, 0.000657009, 0.00111748,  0.0011062,   0.00100514,
      0.000784927, 0.000608221, 0.00111039,  0.00104972,  0.000787295,
      0.000614063, 0.000485616, 0.000338809, 0.00111327,  0.000969593,
      0.00080144,  0.000665887, 0.000582766, 0.000467699, 0.000499902,
      0.00111818,  0.00110917,  0.00111334,  0.00110977,  0.00108202,
      0.00111723,  0.00111212,  0.00110996,  0.00102283,  0.000778352,
      0.000681623, 0.000623996, 0.00151161,  0.00150392,  0.00150577,
      0.00150665,  0.00123625,  0.00102727
  };
  Double_t WeightLHC18q[nGoodRunsLHC18q] = {
      0.000382169, 0.00114974,  0.000527755, 0.000654247, 0.00135902,  0.00214732,
      0.000270681, 0.00166607,  0.00158705,  0.000468779, 0.000615827, 0.000208677,
      0.000199886, 8.60933e-05, 0.000204084, 0.000167014, 0.000217069, 0.00157793,
      0.000184305, 0.000161771, 0.000185254, 0.00235891,  0.00414343,  0.00179014,
      0.00523221,  0.00473166,  0.00881799,  0.00248334,  0.00128146,  0.00483142,
      0.00348662,  0.0142747,   0.00399741,  0.00561353,  0.000702419, 0.0127938,
      0.0104932,   0.00119163,  0.00736677,  0.00435634,  0.00363853,  0.00820298,
      0.00813363,  0.00821918,  0.00489486,  0.00388531,  0.00497056,  0.00373657,
      0.00332877,  0.0142636,   0.00366746,  0.0149605,   0.0145732,   0.00685847,
      0.00188732,  0.00780298,  0.00895397,  0.00779501,  0.00948346,  0.0125277,
      0.00603172,  0.00841124,  0.0125618,   0.011493,    0.00359157,  0.00922165,
      0.00235857,  0.00207436,  0.00981862,  0.00806622,  0.0182107,   0.0109069,
      0.00833837,  0.00248847,  0.0235232,   0.00231616,  0.0134406,   0.00343533,
      0.0112833,   0.0096573,   0.00379312,  0.00390879,  0.00439739,  0.00724271,
      0.038691,    0.00848318,  0.00543936,  0.0177398,   0.00698462,  0.0334586,
      0.00187607,  0.0074507,   0.00931917,  0.0282783,   0.0135473,   0.00972059,
      0.00983623,  0.0281006,   0.0249078,   0.00969403,  0.0133598,   0.0066676,
      0.00700083,  0.0225363,   0.0127199,   0.00653244,  0.00738829,  0.00178683,
      0.0207996,   0.0039916,   0.0143404,   0.0420479,   0.0120012,   0.00226385,
      0.00283461,  0.00500594,  0.00636363,  0.0224793,   0.0184395,   0.00934622,
      0.00223976,  0.00327839,  0.00724002,  0.00249328,  0.0105512,   0.00717872,
      0.00326346,  0.00999924
  };
  const Int_t nGoodRunsLHC18r = 97;
  Int_t GoodRunsLHC18r[nGoodRunsLHC18r] = {
      296690, 296691, 296694, 296749, 296750, 296781, 296784, 296785, 296786,
      296787, 296791, 296793, 296794, 296799, 296836, 296838, 296839, 296848,
      296849, 296850, 296851, 296852, 296890, 296894, 296899, 296900, 296903,
      296930, 296931, 296932, 296934, 296935, 296938, 296941, 296966, 296967,
      296968, 296969, 296971, 296975, 296976, 296979, 297029, 297031, 297035,
      297085, 297117, 297118, 297119, 297123, 297124, 297128, 297129, 297132,
      297133, 297193, 297194, 297196, 297218, 297219, 297221, 297222, 297278,
      297310, 297312, 297315, 297317, 297363, 297366, 297367, 297372, 297379,
      297380, 297405, 297408, 297413, 297414, 297415, 297441, 297442, 297446,
      297450, 297451, 297452, 297479, 297481, 297483, 297512, 297537, 297540,
      297541, 297542, 297544, 297558, 297588, 297590, 297595
  };
  Double_t muLHC18r[nGoodRunsLHC18r] = {
      0.00150804,  0.00150801,  0.00149768,  0.00109303,  0.00109267,
      0.000981497, 0.000982827, 0.0009818,   0.000979027, 0.000981129,
      0.00098248,  0.000983228, 0.000981465, 0.000958229, 0.000982716,
      0.00098012,  0.000979896, 0.000983688, 0.000979531, 0.000882725,
      0.000762617, 0.000700615, 0.000974207, 0.000977741, 0.000977755,
      0.000974729, 0.00091077,  0.000981661, 0.000978104, 0.000979356,
      0.000980038, 0.000980667, 0.000980108, 0.000968508, 0.000984226,
      0.000983014, 0.00120283,  0.00120154,  0.00120146,  0.00119561,
      0.00105211,  0.000688807, 0.000978908, 0.000977728, 0.000975686,
      0.000986131, 0.000981894, 0.000981222, 0.000981072, 0.000979822,
      0.000978465, 0.000979583, 0.000979019, 0.000907481, 0.000784105,
      0.000980856, 0.000980109, 0.000964012, 0.000977384, 0.000969542,
      0.00097847,  0.000948652, 0.000979739, 0.00111646,  0.00111092,
      0.00110773,  0.00110729,  0.00098122,  0.000977469, 0.000977749,
      0.000976518, 0.000975277, 0.000879141, 0.000978345, 0.00114814,
      0.000979787, 0.000978775, 0.000958658, 0.000980799, 0.000975778,
      0.000977995, 0.000976664, 0.000933499, 0.000854953, 0.000980073,
      0.000978981, 0.000971931, 0.00098025,  0.000983023, 0.000976553,
      0.000978715, 0.000979646, 0.000978238, 0.00083193,  0.000982532,
      0.000979778, 0.000976327
  };
  Double_t WeightLHC18r[nGoodRunsLHC18r] = {
      0.021228,   0.00201425, 0.0160721,  0.0287813,  0.0254943,  0.00253203,
      0.00929424, 0.00590839, 0.00233344, 0.0100592,  0.0023444,  0.00432022,
      0.00971419, 0.00859668, 0.00467966, 0.0016817,  0.00919513, 0.00678453,
      0.0368429,  0.00878416, 0.00275611, 0.00296214, 0.0250074,  0.0144777,
      0.0066115,  0.00861758, 0.00321703, 0.00451202, 0.00163823, 0.00371339,
      0.00814181, 0.0139112,  0.00516374, 0.00907643, 0.0104466,  0.00249321,
      0.00987657, 0.00584283, 0.00213909, 0.0231388,  0.00346098, 0.00340511,
      0.0225513,  0.018831,   0.00404415, 0.00302694, 0.00727273, 0.0075655,
      0.00831988, 0.0104444,  0.00197968, 0.00746397, 0.00876137, 0.00872452,
      0.0036229,  0.0236602,  0.0323994,  0.00658033, 0.0198793,  0.0326065,
      0.00876722, 0.00531796, 0.00186332, 0.00207408, 0.00743601, 0.0246112,
      0.0133604,  0.00591977, 0.006592,   0.00976672, 0.00998204, 0.0234177,
      0.00479541, 0.0018597,  0.0126967,  0.00928581, 0.00685629, 0.0211784,
      0.01574,    0.0061538,  0.0251848,  0.00604442, 0.00412684, 0.00356461,
      0.023992,   0.0331809,  0.00618034, 0.00490881, 0.00560265, 0.00194602,
      0.0124206,  0.00480562, 0.0226116,  0.0014836,  0.0163819,  0.00974661,
      0.00513985
  };
  Int_t sizeOfLHC18q = 0;
  Int_t sizeOfLHC18r = 0;
  for ( Int_t iLoopQ = 0; iLoopQ < 128; iLoopQ++ ) {
        fMapGoodRunsToMuAndWeight[GoodRunsLHC18q[iLoopQ]] = std::make_pair( muLHC18q[iLoopQ],
                                                                            WeightLHC18q[iLoopQ]
                                                                            );
        sizeOfLHC18q++;
  }
  // for ( Int_t GoodRunNumberLHC18r : GoodRunsLHC18r ) {
  for ( Int_t iLoopR = 0; iLoopR < 97; iLoopR++ ) {
        fMapGoodRunsToMuAndWeight[GoodRunsLHC18r[iLoopR]] = std::make_pair( muLHC18r[iLoopR],
                                                                            WeightLHC18r[iLoopR]
                                                                            );
        sizeOfLHC18r++;
  }
  Int_t counterRuns = 0;
  for ( auto itr = fMapGoodRunsToMuAndWeight.begin(); itr != fMapGoodRunsToMuAndWeight.end(); ++itr) {
    counterRuns++;
    cout << "Run number:  " << itr->first
         << "  Mu:        " << (itr->second).first
         << "  Weight:    " << (itr->second).second;
    cout << endl;
  }
  cout << "Total runs:" << counterRuns << endl;
  if ( counterRuns != ( nGoodRunsLHC18q + nGoodRunsLHC18r ) ) cout << "OPS!!!!!!" << endl;
}
//_____________________________________________________________________________
void AliAnalysisTaskCTrue::UserCreateOutputObjects()
{
    // define list of histos
    fOutputList = new TList();
    fOutputList ->SetOwner(kTRUE);

    // define hist and add it to the list
    fCounterH = new TH1I("fCounterH", "fCounterH", 100, 1, 101);
    fOutputList->Add(fCounterH);

    /* - Here, the analysis tree gets defined and its branches chosen too.
       -
     */
    fAnaTree = new TTree("fAnaTree", "fAnaTree");
    fAnaTree ->Branch("fRunNum",      &fRunNum,      "fRunNum/I"     );
    fAnaTree ->Branch("fTracklets",   &fTracklets,   "fTracklets/I"  );
    fAnaTree ->Branch("fCtrue",       &fCtrue,       "fCtrue/I"      );
    fAnaTree ->Branch("fC1zed",       &fC1zed,       "fC1zed/I"      );
    fAnaTree ->Branch("fL0inputs",    &fL0inputs,    "L0inputs/i"    );
    fAnaTree ->Branch("fL1inputs",    &fL1inputs,    "L1inputs/i"    );
    fAnaTree ->Branch("fZem1Energy",  &fZem1Energy,  "fZem1Energy/D" );
    fAnaTree ->Branch("fZem2Energy",  &fZem2Energy,  "fZem2Energy/D" );
    fAnaTree ->Branch("fZNCEnergy",   &fZNCEnergy,   "fZNCEnergy/D"  );
    fAnaTree ->Branch("fZNAEnergy",   &fZNAEnergy,   "fZNAEnergy/D"  );
    fAnaTree ->Branch("fZPCEnergy",   &fZPCEnergy,   "fZPCEnergy/D"  );
    fAnaTree ->Branch("fZPAEnergy",   &fZPAEnergy,   "fZPAEnergy/D"  );
    fAnaTree ->Branch("fZNATDC",      &fZNATDC[0],   "fZNATDC[4]/D"  );
    fAnaTree ->Branch("fZNCTDC",      &fZNCTDC[0],   "fZNCTDC[4]/D"  );
    fAnaTree ->Branch("fZPATDC",      &fZPATDC[0],   "fZPATDC[4]/D"  );
    fAnaTree ->Branch("fZPCTDC",      &fZPCTDC[0],   "fZPCTDC[4]/D"  );
    fAnaTree ->Branch("fZNATime",     &fZNATime,     "fZNATime/D"    );
    fAnaTree ->Branch("fZNCTime",     &fZNCTime,     "fZNCTime/D"    );
    fAnaTree ->Branch("fV0ADecision", &fV0ADecision, "fV0ADecision/I");
    fAnaTree ->Branch("fV0CDecision", &fV0CDecision, "fV0CDecision/I");
    fAnaTree ->Branch("fADADecision", &fADADecision, "fADADecision/I");
    fAnaTree ->Branch("fADCDecision", &fADCDecision, "fADCDecision/I");
    fAnaTree ->Branch("fIR1Map",      &fIR1Map                       );
    fAnaTree ->Branch("fIR2Map",      &fIR2Map                       );
    fAnaTree ->Branch("fBCrossNum",   &fBCrossNum,   "fBCrossNum/s"  );



    // fEntriesAgainstRunNumberProperlyH = new TH1F("fEntriesAgainstRunNumberProperlyH", "fEntriesAgainstRunNumberProperlyH", 10000, 290000, 300000);
    fEntriesAgainstRunNumberProperlyH = new TH1F("fEntriesAgainstRunNumberProperlyH", "fEntriesAgainstRunNumberProperlyH", 3, 0, 3);
    fEntriesAgainstRunNumberProperlyH->SetStats(0);
    fEntriesAgainstRunNumberProperlyH->SetFillColor(38);
    // fEntriesAgainstRunNumberProperlyH->SetCanExtend(TH1::kAllAxes);
    fEntriesAgainstRunNumberProperlyH->LabelsDeflate();
    fOutputList->Add(fEntriesAgainstRunNumberProperlyH);

    /*_________________________________________________________________________
    ***************************************************************************
    *                                                                         *
    *                                                                         *
    *                             CTRUE == 1 (B)                              *
    *                                                                         *
    *                                                                         *
    ***************************************************************************
    __________________________________________________________________________*/

    fCTrueBEventsPerRunNumberH = new TH1F("fCTrueBEventsPerRunNumberH", "fCTrueBEventsPerRunNumberH", 3, 0, 3);
    fCTrueBEventsPerRunNumberH->SetStats(0);
    fCTrueBEventsPerRunNumberH->SetFillColor(38);
    // fCTrueBEventsPerRunNumberH->SetCanExtend(TH1::kAllAxes);
    fCTrueBEventsPerRunNumberH->LabelsDeflate();
    fOutputList->Add(fCTrueBEventsPerRunNumberH);

    fCTrueBEventsPerRunNumberConditionsH = new TH1F("fCTrueBEventsPerRunNumberConditionsH", "fCTrueBEventsPerRunNumberConditionsH", 3, 0, 3);
    fCTrueBEventsPerRunNumberConditionsH->SetStats(0);
    fCTrueBEventsPerRunNumberConditionsH->SetFillColor(38);
    // fCTrueBEventsPerRunNumberConditionsH->SetCanExtend(TH1::kAllAxes);
    fCTrueBEventsPerRunNumberConditionsH->LabelsDeflate();
    fOutputList->Add(fCTrueBEventsPerRunNumberConditionsH);

    fTrackletsPerRunNumberH = new TH1F("fTrackletsPerRunNumberH", "fTrackletsPerRunNumberH", 3, 0, 3);
    fTrackletsPerRunNumberH->SetStats(0);
    fTrackletsPerRunNumberH->SetFillColor(38);
    // fTrackletsPerRunNumberH->SetCanExtend(TH1::kAllAxes);
    fTrackletsPerRunNumberH->LabelsDeflate();
    fOutputList->Add(fTrackletsPerRunNumberH);

    //_______________________________
    // VOA plots

    fRecurringVetoAH = new TH1F("fRecurringVetoAH", "fRecurringVetoAH", 3, 0, 3);
    fRecurringVetoAH->SetStats(0);
    fRecurringVetoAH->SetFillColor(38);
    // fRecurringVetoAH->SetCanExtend(TH1::kAllAxes);
    fRecurringVetoAH->LabelsDeflate();
    fOutputList->Add(fRecurringVetoAH);

    fRecurringVetoAtrackletsH = new TH1F("fRecurringVetoAtrackletsH", "fRecurringVetoAtrackletsH", 3, 0, 3);
    fRecurringVetoAtrackletsH->SetStats(0);
    fRecurringVetoAtrackletsH->SetFillColor(38);
    // fRecurringVetoAtrackletsH->SetCanExtend(TH1::kAllAxes);
    fRecurringVetoAtrackletsH->LabelsDeflate();
    fOutputList->Add(fRecurringVetoAtrackletsH);

    fVBAforRunNumberH = new TH1F("fVBAforRunNumberH", "fVBAforRunNumberH", 3, 0, 3);
    fVBAforRunNumberH->SetStats(0);
    fVBAforRunNumberH->SetFillColor(38);
    // fVBAforRunNumberH->SetCanExtend(TH1::kAllAxes);
    fVBAforRunNumberH->LabelsDeflate();
    fOutputList->Add(fVBAforRunNumberH);

    fVBATrackletsForRunNumberH = new TH1F("fVBATrackletsForRunNumberH", "fVBATrackletsForRunNumberH", 3, 0, 3);
    fVBATrackletsForRunNumberH->SetStats(0);
    fVBATrackletsForRunNumberH->SetFillColor(38);
    // fVBATrackletsForRunNumberH->SetCanExtend(TH1::kAllAxes);
    fVBATrackletsForRunNumberH->LabelsDeflate();
    fOutputList->Add(fVBATrackletsForRunNumberH);

    fVDAforRunNumberH = new TH1F("fVDAforRunNumberH", "fVDAforRunNumberH", 3, 0, 3);
    fVDAforRunNumberH->SetStats(0);
    fVDAforRunNumberH->SetFillColor(38);
    // fVDAforRunNumberH->SetCanExtend(TH1::kAllAxes);
    fVDAforRunNumberH->LabelsDeflate();
    fOutputList->Add(fVDAforRunNumberH);

    fVDATrackletsForRunNumberH = new TH1F("fVDATrackletsForRunNumberH", "fVDATrackletsForRunNumberH", 3, 0, 3);
    fVDATrackletsForRunNumberH->SetStats(0);
    fVDATrackletsForRunNumberH->SetFillColor(38);
    // fVDATrackletsForRunNumberH->SetCanExtend(TH1::kAllAxes);
    fVDATrackletsForRunNumberH->LabelsDeflate();
    fOutputList->Add(fVDATrackletsForRunNumberH);

    fVBAandVDAforRunNumberH = new TH1F("fVBAandVDAforRunNumberH", "fVBAandVDAforRunNumberH", 3, 0, 3);
    fVBAandVDAforRunNumberH->SetStats(0);
    fVBAandVDAforRunNumberH->SetFillColor(38);
    // fVBAandVDAforRunNumberH->SetCanExtend(TH1::kAllAxes);
    fVBAandVDAforRunNumberH->LabelsDeflate();
    fOutputList->Add(fVBAandVDAforRunNumberH);

    fVBAandVDATrackletsForRunNumberH = new TH1F("fVBAandVDATrackletsForRunNumberH", "fVBAandVDATrackletsForRunNumberH", 3, 0, 3);
    fVBAandVDATrackletsForRunNumberH->SetStats(0);
    fVBAandVDATrackletsForRunNumberH->SetFillColor(38);
    // fVBAandVDATrackletsForRunNumberH->SetCanExtend(TH1::kAllAxes);
    fVBAandVDATrackletsForRunNumberH->LabelsDeflate();
    fOutputList->Add(fVBAandVDATrackletsForRunNumberH);

    //_______________________________
    // VOC plots

    fRecurringVetoCH = new TH1F("fRecurringVetoCH", "fRecurringVetoCH", 3, 0, 3);
    fRecurringVetoCH->SetStats(0);
    fRecurringVetoCH->SetFillColor(38);
    // fRecurringVetoCH->SetCanExtend(TH1::kAllAxes);
    fRecurringVetoCH->LabelsDeflate();
    fOutputList->Add(fRecurringVetoCH);

    fRecurringVetoCtrackletsH = new TH1F("fRecurringVetoCtrackletsH", "fRecurringVetoCtrackletsH", 3, 0, 3);
    fRecurringVetoCtrackletsH->SetStats(0);
    fRecurringVetoCtrackletsH->SetFillColor(38);
    // fRecurringVetoCtrackletsH->SetCanExtend(TH1::kAllAxes);
    fRecurringVetoCtrackletsH->LabelsDeflate();
    fOutputList->Add(fRecurringVetoCtrackletsH);

    fVBCforRunNumberH = new TH1F("fVBCforRunNumberH", "fVBCforRunNumberH", 3, 0, 3);
    fVBCforRunNumberH->SetStats(0);
    fVBCforRunNumberH->SetFillColor(38);
    // fVBCforRunNumberH->SetCanExtend(TH1::kAllAxes);
    fVBCforRunNumberH->LabelsDeflate();
    fOutputList->Add(fVBCforRunNumberH);

    fVBCTrackletsForRunNumberH = new TH1F("fVBCTrackletsForRunNumberH", "fVBCTrackletsForRunNumberH", 3, 0, 3);
    fVBCTrackletsForRunNumberH->SetStats(0);
    fVBCTrackletsForRunNumberH->SetFillColor(38);
    // fVBCTrackletsForRunNumberH->SetCanExtend(TH1::kAllAxes);
    fVBCTrackletsForRunNumberH->LabelsDeflate();
    fOutputList->Add(fVBCTrackletsForRunNumberH);

    fVDCforRunNumberH = new TH1F("fVDCforRunNumberH", "fVDCforRunNumberH", 3, 0, 3);
    fVDCforRunNumberH->SetStats(0);
    fVDCforRunNumberH->SetFillColor(38);
    // fVDCforRunNumberH->SetCanExtend(TH1::kAllAxes);
    fVDCforRunNumberH->LabelsDeflate();
    fOutputList->Add(fVDCforRunNumberH);

    fVDCTrackletsForRunNumberH = new TH1F("fVDCTrackletsForRunNumberH", "fVDCTrackletsForRunNumberH", 3, 0, 3);
    fVDCTrackletsForRunNumberH->SetStats(0);
    fVDCTrackletsForRunNumberH->SetFillColor(38);
    // fVDCTrackletsForRunNumberH->SetCanExtend(TH1::kAllAxes);
    fVDCTrackletsForRunNumberH->LabelsDeflate();
    fOutputList->Add(fVDCTrackletsForRunNumberH);

    fVBCandVDCforRunNumberH = new TH1F("fVBCandVDCforRunNumberH", "fVBCandVDCforRunNumberH", 3, 0, 3);
    fVBCandVDCforRunNumberH->SetStats(0);
    fVBCandVDCforRunNumberH->SetFillColor(38);
    // fVBCandVDCforRunNumberH->SetCanExtend(TH1::kAllAxes);
    fVBCandVDCforRunNumberH->LabelsDeflate();
    fOutputList->Add(fVBCandVDCforRunNumberH);

    fVBCandVDCTrackletsForRunNumberH = new TH1F("fVBCandVDCTrackletsForRunNumberH", "fVBCandVDCTrackletsForRunNumberH", 3, 0, 3);
    fVBCandVDCTrackletsForRunNumberH->SetStats(0);
    fVBCandVDCTrackletsForRunNumberH->SetFillColor(38);
    // fVBCandVDCTrackletsForRunNumberH->SetCanExtend(TH1::kAllAxes);
    fVBCandVDCTrackletsForRunNumberH->LabelsDeflate();
    fOutputList->Add(fVBCandVDCTrackletsForRunNumberH);

    fCTRUEBfV0CDecisionVSfADCDecisionH = new TH2F(
            "fCTRUEBfV0CDecisionVSfADCDecisionH",
            "fCTRUEBfV0CDecisionVSfADCDecisionH",
            2000, -10, 10,
            2000, -10, 10
            );
    fOutputList->Add(fCTRUEBfV0CDecisionVSfADCDecisionH);

    fCTRUEBfV0ADecisionVSfADADecisionH = new TH2F(
            "fCTRUEBfV0ADecisionVSfADADecisionH",
            "fCTRUEBfV0ADecisionVSfADADecisionH",
            2000, -10, 10,
            2000, -10, 10
            );
    fOutputList->Add(fCTRUEBfV0ADecisionVSfADADecisionH);


    //_______________________________
    // TGRAPHERRORS results for the efficiency computation.
    // fVbaGE = new TGraphErrors();
    // fVbaGE->SetName("fVbaGE");
    // fOutputList->Add(fVbaGE);
    // fVbaGE->SetPoint(0,0,0);
    //
    // fVbaTrackletsGE = new TGraphErrors();
    // fVbaTrackletsGE->SetName("fVbaTrackletsGE");
    // fOutputList->Add(fVbaTrackletsGE);
    // fVbaTrackletsGE->SetPoint(0,0,0);
    //
    //
    // fVdaGE = new TGraphErrors();
    // fVdaGE->SetName("fVdaGE");
    // fOutputList->Add(fVdaGE);
    // fVdaGE->SetPoint(0,0,0);
    //
    // fVdaTrackletsGE = new TGraphErrors();
    // fVdaTrackletsGE->SetName("fVdaTrackletsGE");
    // fOutputList->Add(fVdaTrackletsGE);
    // fVdaTrackletsGE->SetPoint(0,0,0);
    //
    // fVbaAndVdaGE = new TGraphErrors();
    // fVbaAndVdaGE->SetName("fVbaAndVdaGE");
    // fOutputList->Add(fVbaAndVdaGE);
    // fVbaAndVdaGE->SetPoint(0,0,0);
    //
    // fVbaAndVdaTrackletsGE = new TGraphErrors();
    // fVbaAndVdaTrackletsGE->SetName("fVbaAndVdaTrackletsGE");
    // fOutputList->Add(fVbaAndVdaTrackletsGE);
    // fVbaAndVdaTrackletsGE->SetPoint(0,0,0);


    // post output
    PostData(1, fAnaTree);
    PostData(2, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskCTrue::UserExec(Option_t *)
{
  /* - We are inside the UserExec() function, so we
     - invoke the filling of the fCounterH as to remind
     - us of how many events we got.
     -
   */
  fCounterH->Fill(1);

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) {
    PostData(2, fOutputList);
    return;
  }
  /* - Found the input event from the AOD.
     - It is worth remembering how many times
     - it happened...
     -
   */
  fCounterH->Fill(2);

  /* - Trigger selection:
     - here we verify which trigger was the event selected upon.
     - The useful triggers are only those needed for the CTRUE
     - analysis. Hence, if we cannot find them we return...
     -
   */
  fCtrue = -1;
  TString trigger = fAOD->GetFiredTriggerClasses();
  if (trigger.Contains("CTRUE-B")) fCtrue = 1;
  if (trigger.Contains("CTRUE-A")) fCtrue = 2;
  if (trigger.Contains("CTRUE-C")) fCtrue = 3;
  if (trigger.Contains("CTRUE-E")) fCtrue = 4;
  fC1zed = -1;
  if (trigger.Contains("C1ZED-B")) fC1zed = 1;
  if (trigger.Contains("C1ZED-A")) fC1zed = 2;
  if (trigger.Contains("C1ZED-C")) fC1zed = 3;
  if (trigger.Contains("C1ZED-E")) fC1zed = 4;
  if (fCtrue == -1 && fC1zed == -1) {
    PostData(2, fOutputList);
    return;
  }
  fCounterH->Fill(3);

  /* - Requiring a sensible amount of tracks.
     - This means at least one...
     -
   */
  if (fAOD->GetNumberOfTracks() > 0) {
    PostData(2, fOutputList);
    return;
  }
  fCounterH->Fill(4);

  /* - EVENT INFO:
     - fRunNum:    the number of the actual considered run Int_t;
     - fTracklets: the number of SPD tracklets;
     - fBCrossNum: the number of the bunch crossing;
     -
   */
  fRunNum    = fAOD->GetRunNumber();
  fTracklets = fAOD->GetTracklets()->GetNumberOfTracklets();
  fBCrossNum = fAOD->GetBunchCrossNumber();

  // trigger inputs
  fL0inputs  = fAOD->GetHeader()->GetL0TriggerInputs();
  fL1inputs  = fAOD->GetHeader()->GetL1TriggerInputs();

  //Past-future protection maps
  fIR1Map    = fAOD->GetHeader()->GetIRInt1InteractionMap();
  fIR2Map    = fAOD->GetHeader()->GetIRInt2InteractionMap();

  /* - ZDC: we try to find the ZDC object data in the nano-AOD. If we cannot,
     - we return, because there would be no way to actually select the events
     - otherwise! We are here, so we could even check if there is a discrepancy
     - between good events with and without ZDC's information. Or at least, this
     - is my impression of it (filling fCounterH). ZDC information:
     - fZem1Energy:
     - fZem2Energy:
     - fZNAEnergy:
     - fZNCEnergy:
     - fZPAEnergy:
     - fZPCEnergy:
     - fZNATime:
     - fZNCTime:
     - fZNATDC[i]:
     - fZNCTDC[i]:
     - fZPATDC[i]:
     - fZPCTDC[i]:
   */
  AliAODZDC *dataZDC = dynamic_cast<AliAODZDC*>(fAOD->GetZDCData());
  if(!dataZDC) {
    PostData(2, fOutputList);
    return;
  }
  fCounterH->Fill(5);

  fZem1Energy = dataZDC->GetZEM1Energy();
  fZem2Energy = dataZDC->GetZEM2Energy();
  fZNAEnergy  = dataZDC->GetZNATowerEnergy()[0];
  fZNCEnergy  = dataZDC->GetZNCTowerEnergy()[0];
  fZPAEnergy  = dataZDC->GetZPATowerEnergy()[0];
  fZPCEnergy  = dataZDC->GetZPCTowerEnergy()[0];
  /* - NEW: after UPC meeting 5/3/2019
     - On ZDC timing. Usually we use time information from TDCs corresponding to
     - the common PMT (reads all four ZN sectors) on both sides. Each AOD event
     - contains information on up to four consecutive timing hits from these
     - TDCs within +/-12 bcs around the trigger bunch crossing. These hits are
     - stored in fZNATDCm and fZNCTDCm arrays:
     - https://github.com/alisw/AliRoot/blob/master/STEER/AOD/AliAODZDC.h#L153
     - and can be accessed as in:
     -
     - AliAODZDC* aodZDC = aod->GetZDCData();
     - for (Int_t i=0;i<4;i++) fZNATDC[i] = aodZDC->GetZNATDCm(i);
     - for (Int_t i=0;i<4;i++) fZNCTDC[i] = aodZDC->GetZNCTDCm(i);
     -
     - These hits may come from hadronic or EMD processes in neighbouring bcs.
     - In Pb-Pb we usually have 0-2 hits within +/-12 bcs mainly due to EMD.
     - Unused timing slots in these arrays are filled with large negative value
     - (-999). In order to check if there was a timing hit in the trigger bc,
     - you have to check if at least one timing hit out of four is within +/-2
     - ns around 0.
     -
     - Regarding these getters GetZNATime() and GetZNCTime(), defined here:
     - https://github.com/alisw/AliRoot/blob/master/STEER/AOD/AliAODZDC.h#L51
     - They are outdated because, as mentioned here, they return timing
     - information from the first slot in those arrays (fZNATDCm[0], fZNCTDCm[0]):
     - https://github.com/alisw/AliRoot/blob/master/STEER/AOD/AliAODZDC.h#L145
     - The first hit often corresponds to previous bunch crossings (e.g. EMD),
     - while interesting hit around 0 may be stored in the next slots.
     -
   */
  fZNATime    = dataZDC->GetZNATime();
  fZNCTime    = dataZDC->GetZNCTime();

  for (Int_t i=0;i<4;i++) fZNATDC[i] = dataZDC->GetZNATDCm(i);
  for (Int_t i=0;i<4;i++) fZNCTDC[i] = dataZDC->GetZNCTDCm(i);
  for (Int_t i=0;i<4;i++) fZPATDC[i] = dataZDC->GetZPATDCm(i);
  for (Int_t i=0;i<4;i++) fZPCTDC[i] = dataZDC->GetZPCTDCm(i);

  /* - V0: we try to find the V0 object data in the nano-AOD. If we cannot,
     - we return, because there would be no way to actually select the events
     - otherwise! We are here, so we could even check if there is a discrepancy
     - between good events with and without V0's information. Or at least, this
     - is my impression of it (filling fCounterH). V0 information:
     - fV0ADecision;
     - fV0CDecision.
     -
  */
  AliVVZERO *dataVZERO = dynamic_cast<AliVVZERO*>(fAOD->GetVZEROData());
  if(!dataVZERO) {PostData(2, fOutputList); return;}
  fCounterH->Fill(6);

  fV0ADecision = dataVZERO->GetV0ADecision();
  fV0CDecision = dataVZERO->GetV0CDecision();

  /* - AD: we try to find the AD object data in the nano-AOD. If we cannot,
     - we return, because there would be no way to actually select the events
     - otherwise! We are here, so we could even check if there is a discrepancy
     - between good events with and without AD's information. Or at least, this
     - is my impression of it (filling fCounterH). AD information:
     - fADADecision: small detector in ALICE, ADA and ADC at large distances;
     - fADCDecision: again, maybe check whether it is cells or boolean, same as V0.
  */
  AliVAD *dataAD = dynamic_cast<AliVAD*>(fAOD->GetADData());
  if(dataAD) {
    fCounterH->Fill(7);
    fADADecision = dataAD->GetADADecision();
    fADCDecision = dataAD->GetADCDecision();
  }

  /* - Now we have gathered all the information needed for the offline
     - analysis and we can fill the analysis tree.
     - This can be kept for a cross check with the offline macro.
     - Until this point this analysis task is exactly the same as
     - the AliAnalysisTaskTransTask.
     - But now we can implement the full analysis, so that the
     - LEGO train can do everything for us!
     -
   */
  fAnaTree->Fill();


  //_______________________________
  /* - FULL ANALYSIS
     -
     - The analysis is roughly divided in 4 parts:
     - 1) configuration B
     - 2) configuration C
     - 3) configuration E
     - 4) configuration A
     -
     - Until now the macros I was given were only able to
     - analyse type B only. In the latter stages of this analysis
     - task I will try to incorporate even the other three cases.
     - The separation is done at the level of histograms already.
     - What will happen is that I will fill everything inside a
     - if "loop".
     -
     - We first remember that if CTRUE:
     - 1 == B
     - 2 == A
     - 3 == C
     - 4 == E
     - and we will fill the respective histograms...
     -
     - The efficiency retrieval is instead done at the
     - Terminate() level!! So jump up to that point to
     - check the comments... This is needed because this type
     - of analysis require the full statistics available per
     - run!
     -
     - NB: this is the first prototype for the analysis task.
     - If perhaps the Terminate() method doesn't work
     - somehow, another way is feasible. Given that the
     - runs are usually stuck together, I just have to check
     - if the following is a different run and then compute the efficiency
     - for the previous at this point...
     - But this is a messier method and not totally safe.
     -
   */

  // glabal bools
  Bool_t f0VBA = 0;
  Bool_t f0VBC = 0;
  Bool_t f0UBA = 0;
  Bool_t f0UBC = 0;
  Bool_t f1ZED = 0;
  Bool_t f0MUL = 0;
  Bool_t f0STG = 0;

  /* - This type of instruction forces the label to receive a string
     - type of name. So it gets automatically changed to a label of
     - our custom. Every time it gets filled with a weight equal to 1,
     - meaning that it works effectively as a good histogram.
     -
   */
  fEntriesAgainstRunNumberProperlyH->Fill( Form("%d", fRunNum) , 1 );

  // CTRUE == 1 (B)
  if ( fCtrue == 1 ) {
      fCTrueBEventsPerRunNumberH->Fill( Form("%d", fRunNum) , 1 );

      /* - First steps towards a OCDB integration...
         - What I would like to do is to retrieve the
         - L0 inforamtion from the OCDB during runtime.
         - This way even if the Trigger Bit String changes
         - between runs it would be possible to keep track
         - of it and have more realistic results.
         - Until now this had been fixed once for all runs.
         - This is not a bad assumption, considering
         - that 2018 data have been taken right before the
         - end of Run2, meaning that there should not be that
         - much variation between triggers in two consecutive
         - runs!
         -
         -
         -
         - gROOT->Reset();
         - AliCDBManager *man = AliCDBManager::Instance();
         - man->Init();
         - // man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
         - man->SetDefaultStorage("alien://folder=/alice/data/2009/OCDB/");
         - man->SetRun(104065);
         - AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
         - AliTriggerConfiguration* rc = dynamic_cast<AliTriggerConfiguration*>(entry->GetObject());
         - rc->Print();
        */

      // get L0 trigger flags
      f0VBA = fL0inputs & 1 << (inputId_0VBA-1);
      f0VBC = fL0inputs & 1 << (inputId_0VBC-1);
      f0UBA = fL0inputs & 1 << (inputId_0UBA-1);
      f0UBC = fL0inputs & 1 << (inputId_0UBC-1);


      /* - Truth table for fV0CDecision vs fADCDecision
         - (and the same for the other side).
         - This shoudl show the various cells. So True-True and
         - the many configurations.
         -
       */
      fCTRUEBfV0CDecisionVSfADCDecisionH->Fill( fV0CDecision, fADCDecision );
      fCTRUEBfV0ADecisionVSfADADecisionH->Fill( fV0ADecision, fADADecision );

      /* - The raw signal goes down right before it counts. 20ns maybe
         - with -500V. The trigger comes with a discriminator and this tells
         - you when you go more than the threshold value.
         - fV0CDecision tells you when you went over threshold!
         -
         - ADC instead has a reasonably wide timing window, and it
         - measures the charge collected. It measures the integral
         - of the voltage in time, so the overall charge that was
         - deposited and THIS is what is digitized.
         -
         - We can then check offline ONLY if the number in ADC is
         - greater than some random number!!
         -
         -
         - Try to do a plot fV0CDecision vs fADCDecision!
         - Like a truth table...
         -
       */

      if ( !f0UBA && !f0UBC && fADCDecision == 0 && !f0VBC && fV0CDecision == 0 ) {

          if ( !f0VBA && fV0ADecision == 0 ){
            fCTrueBEventsPerRunNumberConditionsH                 ->Fill( Form("%d", fRunNum) , 1 );
            if ( fTracklets > 0 ) fTrackletsPerRunNumberH        ->Fill( Form("%d", fRunNum) , 1 );
          }

          fRecurringVetoAH                                       ->Fill( Form("%d", fRunNum) , 1 );
          if ( fTracklets < 1 ) fRecurringVetoAtrackletsH        ->Fill( Form("%d", fRunNum) , 1 );
          if ( f0VBA  ) {
              fVBAforRunNumberH                                  ->Fill( Form("%d", fRunNum) , 1 );
              if( fTracklets < 1 ) fVBATrackletsForRunNumberH    ->Fill( Form("%d", fRunNum) , 1 );
          }
          if ( fV0ADecision != 0 ) {
              fVDAforRunNumberH                                  ->Fill( Form("%d", fRunNum) , 1 );
              if( fTracklets < 1 ) fVDATrackletsForRunNumberH    ->Fill( Form("%d", fRunNum) , 1 );
          }
          if ( f0VBA || fV0ADecision != 0 ) {
            fVBAandVDAforRunNumberH                              ->Fill( Form("%d", fRunNum) , 1 );
            if( fTracklets < 1 ) fVBAandVDATrackletsForRunNumberH->Fill( Form("%d", fRunNum) , 1 );

          }
      }


      if ( !f0UBA && !f0UBC && fADCDecision == 0 && !f0VBA && fV0ADecision == 0 ) {
          fRecurringVetoCH                                       ->Fill( Form("%d", fRunNum) , 1 );
          if ( fTracklets < 1 ) fRecurringVetoCtrackletsH        ->Fill( Form("%d", fRunNum) , 1 );
          if ( f0VBC  ) {
              fVBCforRunNumberH                                  ->Fill( Form("%d", fRunNum) , 1 );
              if( fTracklets < 1 ) fVBCTrackletsForRunNumberH    ->Fill( Form("%d", fRunNum) , 1 );
          }
          if ( fV0CDecision != 0 ) {
              fVDCforRunNumberH                                  ->Fill( Form("%d", fRunNum) , 1 );
              if( fTracklets < 1 ) fVDCTrackletsForRunNumberH    ->Fill( Form("%d", fRunNum) , 1 );
          }
          if ( f0VBC || fV0CDecision != 0 ) {
            fVBCandVDCforRunNumberH                              ->Fill( Form("%d", fRunNum) , 1 );
            if( fTracklets < 1 ) fVBCandVDCTrackletsForRunNumberH->Fill( Form("%d", fRunNum) , 1 );
          }
      }

  }

  PostData(1, fAnaTree);
  PostData(2, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskCTrue::Set2018PbPb()
{

  // --------------------------------
  // defines variables with the index of the
  // different L0/L1 trigger elements according to
  // Run 296786 (LHC18r), checked with 295666 (LHC18q)
  // --------------------------------

  inputId_0VBA = 1; // 0VBA: >=1 V0A cell fired in BB timing gate
  inputId_0VBC = 2; // 0VBC: >=1 V0C cell fired in BB timing gate
  inputId_0UBA = 4; // 0UBA: >=1 ADA cell fired in BB timing gate
  inputId_0UBC = 5; // 0UBC: >=1 ADC cell fired in BB timing gate
  inputId_0VOM = 7;
  inputId_0STG = 16;
  inputId_0MUL = 19;
  inputId_1ZED = 15;

}
//_____________________________________________________________________________
void AliAnalysisTaskCTrue::Set2015PbPb()
{

  // --------------------------------
  // defines variables with the index of the
  // different L0/L1 trigger elements according to
  // Run 245064
  // --------------------------------

  inputId_0VBA = 1;  // 0VBA: >=1 V0A cell fired in BB timing gate
  inputId_0VBC = 2;  // 0VBC: >=1 V0C cell fired in BB timing gate
  inputId_0UBA = 7;  // 0UBA: >=1 ADA cell fired in BB timing gate
  inputId_0UBC = 8;  // 0UBC: >=1 ADC cell fired in BB timing gate
  inputId_0SH1 = 13; // 0SH1: >=2 outer FO hits
  inputId_1ZED = 15;

}
//_____________________________________________________________________________
void AliAnalysisTaskCTrue::SetXeXe()
{

  // --------------------------------
  // defines variables with the index of the
  // different L0/L1 trigger elements according to
  // Run 280235
  // --------------------------------

  inputId_0VBA = 1; // 0VBA: >=1 V0A cell fired in BB timing gate
  inputId_0VBC = 2; // 0VBC: >=1 V0C cell fired in BB timing gate
  inputId_0UBA = 6; // 0UBA: >=1 ADA cell fired in BB timing gate
  inputId_0UBC = 7; // 0UBC: >=1 ADC cell fired in BB timing gate
  inputId_0SH1 = 8; // 0SH1: >=2 outer and >=2 inner FO hits
  inputId_1ZED = 15;

  inputId_0OM2 = 9;
}
//_____________________________________________________________________________
// Double_t AliAnalysisTaskCTrue::BinomialError(Double_t cut, Double_t CTRUE)
// // use Root formula for binomial error
// {
//   Double_t efficiency = cut/CTRUE;
//   // Double_t e1 = TMath::Sqrt(cut);   // seem unused. I commented them out!!
//   // Double_t e2 = TMath::Sqrt(ctrue); // seem unused. I commented them out!!
//   /* - If I did the math right this should be equivalent to the new method:
//      -
//    */
//   // return TMath::Sqrt(TMath::Abs( ((1.0-2.0*efficiency)*cut + efficiency*efficiency*ctrue)/(ctrue*ctrue) ));
//   return TMath::Sqrt(TMath::Abs( efficiency*(1-efficiency)/CTRUE ));
// }
// //_____________________________________________________________________________
// Double_t AliAnalysisTaskCTrue::FitPolinomial(Double_t *x, Double_t *p)
// // fit model of a pol1 passing through 0
// {
//   return p[0]+p[1]*x[0];
// }
// //_____________________________________________________________________________
// void AliAnalysisTaskCTrue::ComputeEfficiency( Int_t n,
//                                               Double_t *weight,
//                                               Double_t *mu,
// 	                                            Double_t p0,
//                                               Double_t p0e,
//                                               Double_t p1,
//                                               Double_t p1e,
// 	                                            Double_t *eff
//                                               // std::vector<Double_t>* fVectorEfficiencyC
//                                               )
// {
//   // Double_t w_tot = 0.0;
//   // for(Int_t i=0;i<n;i++) {
//   //   Double_t w = weight[i];
//   //   w_tot += w;
//   //   Double_t p = p0+p1*mu[i];
//   //   eff[0] += (w*TMath::Exp(-p));
//   //   Double_t pm = p0+(p1-p1e)*mu[i];
//   //   eff[1] += (w*TMath::Exp(-pm));
//   //   Double_t pp = p0+(p1+p1e)*mu[i];
//   //   eff[2] += (w*TMath::Exp(-pp));
//   //   //  cout << i << " w " << w << " p " << p << " eff " <<  (w*TMath::Exp(-p)) << endl;
//   // }
//   // if (w_tot < 1e-12) return;
//   // for(Int_t i=0;i<3;i++) eff[i] = eff[i]/w_tot;
//   /*
//   cout << " efficiency  = " << eff[0]
//        << " + " << (eff[1]-eff[0])
//        << " - " << (eff[0]-eff[2])
//        << endl;
//   */
//   for(Int_t i=0;i<n;i++) {
//     TotalWeightForEfficiency += weight[i];
//     Double_t p  = p0+      p1*mu[i];
//     Double_t pm = p0+(p1-p1e)*mu[i];
//     Double_t pp = p0+(p1+p1e)*mu[i];
//
//     fVectorEfficiency[0].push_back( weight[i]*TMath::Exp(-p ) );
//     fVectorEfficiency[1].push_back( weight[i]*TMath::Exp(-pm) );
//     fVectorEfficiency[2].push_back( weight[i]*TMath::Exp(-pp) );
//   }
// }
// //_____________________________________________________________________________
// Double_t* AliAnalysisTaskCTrue::NormalizeEfficiency()
// {
//   Double_t effInTheCycle[3];
//   effInTheCycle[0] = effInTheCycle[1] = effInTheCycle[2] = 0;
//   for ( Int_t iLoop = 0; iLoop < 3; iLoop++){
//     effInTheCycle[iLoop] = 0;
//     for ( auto it = fVectorEfficiency[0].begin(); it != fVectorEfficiency[0].end(); it ++ ) {
//       effInTheCycle[iLoop] += (*it);
//     }
//     effInTheCycle[iLoop] /= TotalWeightForEfficiency;
//   }
//   Double_t* returningEff = new Double_t[3];
//   returningEff[0] = effInTheCycle[0];
//   returningEff[1] = effInTheCycle[1] - effInTheCycle[0];
//   returningEff[1] = effInTheCycle[0] - effInTheCycle[2];
//
//   return returningEff;
// }
//_____________________________________________________________________________
// void  AliAnalysisTaskCTrue::DoPlot( TH1F*         signal,
//                                     TH1F*         background,
//                                     TGraphErrors* graphToBeFilled
//                                     )
// {
//   // decide if constraint the origin of the line fit
//   // to pass through the origin
//   Bool_t fix_p0 = kFALSE;
//
//   Int_t nbinsx = signal->GetXaxis()->GetNbins();
//   const char* binLabel;
//   Double_t x, y, yError, weight = 0;
//
//
//   /* - Firstly, we get the label of the bin.
//      - This label is actually the run number!!
//      - Then we look for this "key" inside the map.
//      - Hence, we retrieve this information.
//      - And now we can finish the analysis.
//      -
//    */
//   for (Int_t iLoop = 1; iLoop < nbinsx+1; iLoop++) {
//     binLabel = signal->GetXaxis()->GetBinLabel(iLoop);
//     // cout << "binLabel:   " << binLabel << atoi(binLabel) << endl;
//     if ( atoi(binLabel) == 0 ) continue;
//     if ( signal    ->GetBinContent(iLoop) < 0.5 ) continue;
//     if ( background->GetBinContent(iLoop) < 0.5 ) continue;
//
//     auto iterator = fMapGoodRunsToMuAndWeight.find( atoi(binLabel) );
//     if ( iterator == fMapGoodRunsToMuAndWeight.end() ) {
//         cout << "Missing run number: not a problem but be careful" << endl;
//     } else {
//         x      = (iterator->second).first;
//         cout << "x      = " << x      << endl;
//         weight = (iterator->second).second;
//         cout << "weight = " << weight << endl;
//     }
//
//     y      = (Double_t)signal->GetBinContent(iLoop) / (Double_t)background->GetBinContent(iLoop);
//     cout << "y      = " << y      << endl;
//     yError = BinomialError( (Double_t)signal->GetBinContent(iLoop), (Double_t)background->GetBinContent(iLoop) );
//     cout << "yError = " << yError << endl;
//
//     fVbcGE->Draw();
//     Int_t nOfPoints = fVbcGE->GetN();
//     cout << "nOfPoints   " << nOfPoints << endl;
//     fVbcGE->SetPoint     ( nOfPoints, x, y     );
//     fVbcGE->SetPointError( nOfPoints, 0, yError);
//
//   }
//   // // prepare for fractions to store weights to total efficiency
//   // Double_t *weight = new Double_t [nRuns];
//   // for (Int_t i=0;i<nRuns;i++) weight[i]=0.0;
//   //
//   // // prepare graphs
//   // Double_t ymin =10;
//   // Double_t ymax =-10;
//   // TGraphErrors *gr = new TGraphErrors();
//   // for(Int_t i=0;i<nRuns;i++) {
//   //   Double_t ye = berr(signal[i],empty[i]);
//   //   if (ymin>y-ye) ymin = y-ye;
//   //   if (ymax<y+ye) ymax = y+ye;
//   //   // if (y < 0.01) continue;
//   //   // fill graph
//   //   Int_t n = gr->GetN();
//   //   gr->SetPoint(n,x,y);
//   //   gr->SetPointError(n,0,ye);
//   //   weight[n]=w_all[i];
//   // }
//   //
//   //
//   // //plot graph
//   // // No stats and no titles
//   // gStyle->SetOptTitle(1);
//   // gStyle->SetOptStat(0);
//   // gStyle->SetOptFit(1);
//   //
//   //
//   // TCanvas *c = new TCanvas(title,title,1200,600);
//   // c->Divide(1,1);
//   // c->cd(1);
//   // gr->SetMarkerStyle(20);  gr->SetMarkerColor(kRed); gr->SetLineColor(kRed);
//   // // make some white space in the canvas
//   // ymin *= 0.5;
//   // if (ymin<0) ymin = 0.001;
//   // ymax *= 1.5;
//   // // define the histo
//   // Double_t mu_max = 0.0025;
//   // TH1F* frame1 = gPad->DrawFrame(0.,ymin,mu_max,ymax);
//   // frame1->GetYaxis()->SetTitleOffset(1.3);
//   // frame1->GetYaxis()->SetLabelSize(0.025);
//   // frame1->SetTitle(Form("%s;#mu;Prob.",title));
//   // TF1 *pol = new TF1("pol",fit_p1, 0, 1,2);
//   // if (gr->GetN()>0) {
//   //   gr->Draw("p,same");
//   //   pol->SetParNames("p0","p1");
//   //   pol->SetParameter(0,0.0);
//   //   pol->SetParameter(1,1.0);
//   //   if (fix_p0) pol->FixParameter(0,0.0);
//   //   gr->Fit("pol");
//   // } else { cout << " gr is empty " << endl;}
//   // gPad->Update();
//   //
//   // // compute efficiency
//   // Double_t eff[3]; // mean, minus, plus
//   // eff[0]=eff[1]=eff[2]=0.0;
//   // DoEff(gr->GetN(),weight,gr->GetX(),pol->GetParameter(0), pol->GetParError(0),
// 	// pol->GetParameter(1), pol->GetParError(1), eff);
//   //
//   // // print eff in the canvas
//   // TLatex* l = new TLatex();
//   // l->SetTextFont(42);
//   // l->SetTextSize(0.045);
//   // l->DrawLatex(mu_max*0.05,ymax*0.9,Form("#varepsilon = %.5f + %.5f - %.5f",eff[0],(eff[1]-eff[0]),(eff[0]-eff[2])));
//   //
//   // // save
//   // c->Print(Form("Canvas_%s.pdf",title));
//   // // clean up
//   // delete [] weight;
// }
//_____________________________________________________________________________
void AliAnalysisTaskCTrue::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
  // define hist and add it to the list
  // fCounterTryingH = new TH1I("fCounterTryingH", "fCounterTryingH", 100, 1, 101);
  // fOutputList->Add(fCounterTryingH);
  //
  // DoPlot( fVBCforRunNumberH,                fRecurringVetoCH,          fVbcGE );
  // DoPlot( fVDCforRunNumberH,                fRecurringVetoCH,          fVdcGE );
  // DoPlot( fVBCTrackletsForRunNumberH,       fRecurringVetoCtrackletsH, fVbcTrackletsGE );
  // DoPlot( fVDCTrackletsForRunNumberH,       fRecurringVetoCtrackletsH, fVdcTrackletsGE );
  // DoPlot( fVBCandVDCforRunNumberH,          fRecurringVetoCH,          fVbcAndVdcGE );
  // DoPlot( fVBCandVDCTrackletsForRunNumberH, fRecurringVetoCtrackletsH, fVbcAndVdcTrackletsGE );
  // DoPlot( fVBAforRunNumberH,                fRecurringVetoAH,          fVbaGE );
  // DoPlot( fVDAforRunNumberH,                fRecurringVetoAH,          fVdaGE );
  // DoPlot( fVBATrackletsForRunNumberH,       fRecurringVetoAtrackletsH, fVbaTrackletsGE );
  // DoPlot( fVDATrackletsForRunNumberH,       fRecurringVetoAtrackletsH, fVdaTrackletsGE );
  // DoPlot( fVBAandVDAforRunNumberH,          fRecurringVetoAH,          fVbaAndVdaGE );
  // DoPlot( fVBAandVDATrackletsForRunNumberH, fRecurringVetoAtrackletsH, fVbaAndVdaTrackletsGE );

}
//_____________________________________________________________________________
