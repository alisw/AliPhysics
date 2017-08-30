/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/***********************************
 * ZDC Event Plane                 *
 *                                 *
 * author: Jacopo Margutti         *
 * email:  jacopo.margutti@cern.ch *
 ***********************************/

#define AliAnalysisTaskZDCEP_cxx

#include "Riostream.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "TChain.h"

#include "TFile.h"
#include "TList.h"
#include "TGraph.h"
#include "TParticle.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TArrow.h"
#include "TPaveLabel.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "AliFlowEventSimple.h"
#include "AliFlowVector.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowAnalysisCRC.h"
#include "TRandom.h"
#include "TF1.h"
#include "TNtuple.h"
#include "THnSparse.h"
#include "AliAnalysisTaskZDCEP.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODZDC.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "AliFlowEvent.h"

class TH1;
class TH2;
class TGraph;
class TPave;
class TLatex;
class TMarker;
class TRandom3;
class TObjArray;
class TList;
class TCanvas;
class TSystem;
class TROOT;
class TVector2;
class AliFlowVector;

AliAnalysisTaskZDCEP::AliAnalysisTaskZDCEP () :
AliAnalysisTaskSE (),
fOutputList(0x0),
fHistList(0x0),
fZDCGainAlpha(0.395),
fZDCCalibList(0x0),
fTowerEqList(0x0),
fCachedRunNum(0),
fAnalysisUtils(0x0),
fMultSelection(0x0),
fbFlagIsPosMagField(kFALSE),
fFlowEvent(NULL)
{
  for(Int_t k=0; k<4; k++) {
    fZDCQHist[k] = NULL;
    fZDCVtxHist[k] = NULL;
    fZDCEcomTotHist[k] = NULL;
    for(Int_t c=0; c<10; c++) {
      fZDCVtxCenHist[c][k] = NULL;
    }
    fZDCVtxFitHist[k] = NULL;
    for(Int_t i=0; i<3; i++) {
      fZDCVtxFitCenProjHist[k][i] = NULL;
    }
  }
  for(Int_t c=0; c<10; c++) {
    for(Int_t k=0; k<8; k++) {
      fZDCVtxCenHistMagPol[c][k] = NULL;
    }
  }
  for(Int_t i=0; i<10; i++) {
    for(Int_t z=0; z<10; z++) {
      for(Int_t k=0; k<4; k++) {
        fZDCQVecVtxCenEZDC3D[i][z][k] = NULL;
      }
    }
  }
  for(Int_t c=0; c<2; c++) {
    for(Int_t i=0; i<5; i++) {
      fTowerGainEq[c][i] =  NULL;
    }
  }
  for (Int_t i=0; i<10; i++) {
    fCRCZDCQVecDummyEZDCBins[i] = NULL;
  }
  for (Int_t i=0; i<2; i++) fZDCFlowVect[i] = NULL;
  
  Int_t dRun15o[] = {244917, 244918, 244975, 244980, 244982, 244983, 245064, 245066, 245068, 246390, 246391, 246392, 246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683, 245145, 245146, 245151, 245152, 245231, 245232, 245259, 245343, 245345, 245346, 245347, 245349, 245353, 245396, 245397, 245401, 245407, 245409, 245441, 245446, 245450, 245454, 245496, 245497, 245501, 245504, 245505, 245507, 245535, 245540, 245542, 245543, 245544, 245545, 245554};
  Double_t dVtxPosX15o[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,7.619407e-02, 7.612905e-02, 7.609009e-02, 7.610981e-02, 7.608885e-02, 7.609981e-02, 7.559263e-02, 7.563009e-02, 7.551201e-02, 7.570994e-02, 7.571927e-02, 7.575639e-02, 7.571133e-02, 7.570653e-02, 7.528412e-02, 7.535235e-02, 7.539954e-02, 7.535435e-02, 7.541641e-02, 7.543658e-02, 7.527343e-02, 7.526024e-02, 7.528295e-02, 7.533821e-02, 7.540461e-02, 7.538317e-02, 7.531677e-02, 7.539861e-02, 7.537667e-02, 7.659318e-02, 7.656796e-02, 7.662898e-02, 7.664257e-02, 7.597872e-02, 7.597437e-02, 7.599091e-02, 7.601310e-02, 7.000359e-02, 6.999659e-02, 6.992559e-02, 6.996793e-02, 7.028519e-02, 7.032696e-02, 7.033503e-02, 6.952509e-02, 6.956378e-02, 6.952446e-02, 6.959759e-02, 6.956048e-02, 6.933134e-02, 6.932882e-02, 6.939338e-02, 6.950613e-02, 6.943631e-02, 6.946196e-02, 6.950454e-02, 7.030973e-02, 7.030203e-02, 7.032272e-02, 7.030936e-02, 7.038967e-02, 7.035136e-02, 7.024752e-02, 6.942316e-02, 6.940115e-02, 6.936367e-02, 6.860689e-02, 6.881501e-02, 6.886743e-02, 6.932714e-02, 6.970325e-02, 6.966504e-02, 6.957355e-02, 6.932303e-02, 6.938184e-02, 6.944933e-02, 6.952461e-02, 6.964167e-02, 6.793435e-02, 6.802185e-02, 6.801235e-02, 6.804823e-02, 6.842972e-02, 6.839652e-02, 6.851932e-02, 6.976507e-02, 6.989692e-02, 6.994544e-02, 6.994261e-02, 6.997887e-02, 7.001687e-02, 6.934462e-02, 6.958349e-02, 6.907266e-02, 6.905944e-02, 6.895395e-02, 7.006562e-02, 7.008493e-02, 7.012736e-02, 6.964645e-02, 6.960466e-02, 6.962255e-02, 6.979086e-02, 6.985343e-02, 6.983755e-02, 6.957177e-02, 6.875991e-02, 6.871756e-02, 6.871021e-02, 6.871769e-02, 6.869493e-02, 6.874049e-02, 6.860300e-02};
  Double_t dVtxPosY15o[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,3.361709e-01, 3.361818e-01, 3.362205e-01, 3.363199e-01, 3.363092e-01, 3.362369e-01, 3.374328e-01, 3.374148e-01, 3.375140e-01, 3.361514e-01, 3.361743e-01, 3.362329e-01, 3.361395e-01, 3.361633e-01, 3.367675e-01, 3.366963e-01, 3.366845e-01, 3.366490e-01, 3.366937e-01, 3.366825e-01, 3.373764e-01, 3.373762e-01, 3.373721e-01, 3.373705e-01, 3.373943e-01, 3.373675e-01, 3.374071e-01, 3.373368e-01, 3.373442e-01, 3.375773e-01, 3.375333e-01, 3.377335e-01, 3.378285e-01, 3.362674e-01, 3.362492e-01, 3.362604e-01, 3.363473e-01, 3.295003e-01, 3.295046e-01, 3.295761e-01, 3.296100e-01, 3.291527e-01, 3.292071e-01, 3.290824e-01, 3.299371e-01, 3.300008e-01, 3.300078e-01, 3.300391e-01, 3.300740e-01, 3.300345e-01, 3.300776e-01, 3.301195e-01, 3.289427e-01, 3.289736e-01, 3.296084e-01, 3.297025e-01, 3.297724e-01, 3.298166e-01, 3.298278e-01, 3.298682e-01, 3.297381e-01, 3.296875e-01, 3.297720e-01, 3.298361e-01, 3.298561e-01, 3.299325e-01, 3.300111e-01, 3.301161e-01, 3.302630e-01, 3.289954e-01, 3.292915e-01, 3.293319e-01, 3.294174e-01, 3.314355e-01, 3.314431e-01, 3.316189e-01, 3.318682e-01, 3.323906e-01, 3.315020e-01, 3.312268e-01, 3.310778e-01, 3.310524e-01, 3.314478e-01, 3.312986e-01, 3.311297e-01, 3.324064e-01, 3.322524e-01, 3.322019e-01, 3.321221e-01, 3.321050e-01, 3.319118e-01, 3.317922e-01, 3.314658e-01, 3.315735e-01, 3.316331e-01, 3.316525e-01, 3.308030e-01, 3.308038e-01, 3.306947e-01, 3.305741e-01, 3.316492e-01, 3.316117e-01, 3.314973e-01, 3.314110e-01, 3.313450e-01, 3.313649e-01, 3.325841e-01, 3.324226e-01, 3.323649e-01, 3.323381e-01, 3.322566e-01, 3.322077e-01, 3.320860e-01};
  Double_t dVtxPosZ15o[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,5.559279e-01, 3.535446e-01, 4.846955e-01, 4.525585e-01, 3.684501e-01, 2.485494e-01, 2.372653e-01, 1.707859e-01, 3.314213e-01, 1.709195e-01, 2.209753e-01, 3.125757e-01, 3.422085e-01, 3.868156e-01, 4.859695e-01, 4.780697e-01, 4.400149e-01, 4.014992e-01, 3.049883e-01, 3.708501e-01, 3.883566e-01, 3.940632e-01, 4.197670e-01, 3.938399e-01, 3.814413e-01, 3.335539e-01, 3.181929e-01, 2.300734e-01, 2.722395e-01, 5.241033e-01, 3.225908e-01, 1.925791e-01, 1.892765e-01, 3.384066e-01, 2.026459e-01, 2.495699e-01, 3.569992e-01, 3.891381e-01, 4.603724e-01, 3.696685e-01, 3.002207e-01, 2.929533e-01, 3.095468e-01, 3.517200e-01, 2.784445e-01, 3.866626e-01, 3.058719e-01, 3.336752e-01, 3.226473e-01, 3.222815e-01, 3.428469e-01, 3.728514e-01, 2.858642e-01, 2.832485e-01, 3.378933e-01, 3.547548e-01, 3.799414e-01, 4.043543e-01, 4.314049e-01, 4.141138e-01, 3.888746e-01, 4.103586e-01, 3.871045e-01, 4.614473e-01, 4.023404e-01, 4.203531e-01, 4.401272e-01, 6.450558e-01, 6.819582e-01, 2.588529e-01, 3.693471e-01, 3.990708e-01, 3.813842e-01, 3.471682e-01, 3.356156e-01, 2.550150e-01, 3.830723e-01, 4.293259e-01, 4.723797e-01, 4.684324e-01, 4.609304e-01, 4.554974e-01, 4.523016e-01, 3.769890e-01, 4.485548e-01, 5.024484e-01, 5.200088e-01, 5.261731e-01, 5.392851e-01, 5.399264e-01, 5.155504e-01, 4.267668e-01, 5.348764e-01, 4.526746e-01, 4.045626e-01, 4.261759e-01, 5.889205e-01, 6.364843e-01, 5.896163e-01, 3.768637e-01, 4.440771e-01, 4.687029e-01, 4.794467e-01, 4.313422e-01, 3.954777e-01, 3.983129e-01, 3.608064e-01, 2.627038e-01, 3.665826e-01, 4.275667e-01, 3.335445e-01, 3.250815e-01, 3.022907e-01};
  for(Int_t r=0; r<fnRun; r++) {
    fRunList[r] = dRun15o[r];
  }
  fAvVtxPosX=TArrayD(fnRun,dVtxPosX15o);
  fAvVtxPosY=TArrayD(fnRun,dVtxPosY15o);
  fAvVtxPosZ=TArrayD(fnRun,dVtxPosZ15o);
}

//=====================================================================

AliAnalysisTaskZDCEP::AliAnalysisTaskZDCEP(const  char* name)
: AliAnalysisTaskSE(name),
fOutputList(0x0),
fHistList(0x0),
fZDCGainAlpha(0.395),
fZDCCalibList(0x0),
fTowerEqList(0x0),
fCachedRunNum(0),
fAnalysisUtils(0x0),
fMultSelection(0x0),
fbFlagIsPosMagField(kFALSE),
fFlowEvent(NULL)
{
  for(Int_t k=0; k<4; k++) {
    fZDCQHist[k] = NULL;
    fZDCVtxHist[k] = NULL;
    fZDCEcomTotHist[k] = NULL;
    for(Int_t c=0; c<10; c++) {
      fZDCVtxCenHist[c][k] = NULL;
    }
    fZDCVtxFitHist[k] = NULL;
    for(Int_t i=0; i<3; i++) {
      fZDCVtxFitCenProjHist[k][i] = NULL;
    }
  }
  for(Int_t c=0; c<10; c++) {
    for(Int_t k=0; k<8; k++) {
      fZDCVtxCenHistMagPol[c][k] = NULL;
    }
  }
  for(Int_t i=0; i<10; i++) {
    for(Int_t z=0; z<10; z++) {
      for(Int_t k=0; k<4; k++) {
        fZDCQVecVtxCenEZDC3D[i][z][k] = NULL;
      }
    }
  }
  for(Int_t c=0; c<2; c++) {
    for(Int_t i=0; i<5; i++) {
      fTowerGainEq[c][i] =  NULL;
    }
  }
  for (Int_t i=0; i<10; i++) {
    fCRCZDCQVecDummyEZDCBins[i] = NULL;
  }
  for (Int_t i=0; i<2; i++) fZDCFlowVect[i] = NULL;
  
  Int_t dRun15o[] = {244917, 244918, 244975, 244980, 244982, 244983, 245064, 245066, 245068, 246390, 246391, 246392, 246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683, 245145, 245146, 245151, 245152, 245231, 245232, 245259, 245343, 245345, 245346, 245347, 245349, 245353, 245396, 245397, 245401, 245407, 245409, 245441, 245446, 245450, 245454, 245496, 245497, 245501, 245504, 245505, 245507, 245535, 245540, 245542, 245543, 245544, 245545, 245554};
  Double_t dVtxPosX15o[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,7.619407e-02, 7.612905e-02, 7.609009e-02, 7.610981e-02, 7.608885e-02, 7.609981e-02, 7.559263e-02, 7.563009e-02, 7.551201e-02, 7.570994e-02, 7.571927e-02, 7.575639e-02, 7.571133e-02, 7.570653e-02, 7.528412e-02, 7.535235e-02, 7.539954e-02, 7.535435e-02, 7.541641e-02, 7.543658e-02, 7.527343e-02, 7.526024e-02, 7.528295e-02, 7.533821e-02, 7.540461e-02, 7.538317e-02, 7.531677e-02, 7.539861e-02, 7.537667e-02, 7.659318e-02, 7.656796e-02, 7.662898e-02, 7.664257e-02, 7.597872e-02, 7.597437e-02, 7.599091e-02, 7.601310e-02, 7.000359e-02, 6.999659e-02, 6.992559e-02, 6.996793e-02, 7.028519e-02, 7.032696e-02, 7.033503e-02, 6.952509e-02, 6.956378e-02, 6.952446e-02, 6.959759e-02, 6.956048e-02, 6.933134e-02, 6.932882e-02, 6.939338e-02, 6.950613e-02, 6.943631e-02, 6.946196e-02, 6.950454e-02, 7.030973e-02, 7.030203e-02, 7.032272e-02, 7.030936e-02, 7.038967e-02, 7.035136e-02, 7.024752e-02, 6.942316e-02, 6.940115e-02, 6.936367e-02, 6.860689e-02, 6.881501e-02, 6.886743e-02, 6.932714e-02, 6.970325e-02, 6.966504e-02, 6.957355e-02, 6.932303e-02, 6.938184e-02, 6.944933e-02, 6.952461e-02, 6.964167e-02, 6.793435e-02, 6.802185e-02, 6.801235e-02, 6.804823e-02, 6.842972e-02, 6.839652e-02, 6.851932e-02, 6.976507e-02, 6.989692e-02, 6.994544e-02, 6.994261e-02, 6.997887e-02, 7.001687e-02, 6.934462e-02, 6.958349e-02, 6.907266e-02, 6.905944e-02, 6.895395e-02, 7.006562e-02, 7.008493e-02, 7.012736e-02, 6.964645e-02, 6.960466e-02, 6.962255e-02, 6.979086e-02, 6.985343e-02, 6.983755e-02, 6.957177e-02, 6.875991e-02, 6.871756e-02, 6.871021e-02, 6.871769e-02, 6.869493e-02, 6.874049e-02, 6.860300e-02};
  Double_t dVtxPosY15o[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,3.361709e-01, 3.361818e-01, 3.362205e-01, 3.363199e-01, 3.363092e-01, 3.362369e-01, 3.374328e-01, 3.374148e-01, 3.375140e-01, 3.361514e-01, 3.361743e-01, 3.362329e-01, 3.361395e-01, 3.361633e-01, 3.367675e-01, 3.366963e-01, 3.366845e-01, 3.366490e-01, 3.366937e-01, 3.366825e-01, 3.373764e-01, 3.373762e-01, 3.373721e-01, 3.373705e-01, 3.373943e-01, 3.373675e-01, 3.374071e-01, 3.373368e-01, 3.373442e-01, 3.375773e-01, 3.375333e-01, 3.377335e-01, 3.378285e-01, 3.362674e-01, 3.362492e-01, 3.362604e-01, 3.363473e-01, 3.295003e-01, 3.295046e-01, 3.295761e-01, 3.296100e-01, 3.291527e-01, 3.292071e-01, 3.290824e-01, 3.299371e-01, 3.300008e-01, 3.300078e-01, 3.300391e-01, 3.300740e-01, 3.300345e-01, 3.300776e-01, 3.301195e-01, 3.289427e-01, 3.289736e-01, 3.296084e-01, 3.297025e-01, 3.297724e-01, 3.298166e-01, 3.298278e-01, 3.298682e-01, 3.297381e-01, 3.296875e-01, 3.297720e-01, 3.298361e-01, 3.298561e-01, 3.299325e-01, 3.300111e-01, 3.301161e-01, 3.302630e-01, 3.289954e-01, 3.292915e-01, 3.293319e-01, 3.294174e-01, 3.314355e-01, 3.314431e-01, 3.316189e-01, 3.318682e-01, 3.323906e-01, 3.315020e-01, 3.312268e-01, 3.310778e-01, 3.310524e-01, 3.314478e-01, 3.312986e-01, 3.311297e-01, 3.324064e-01, 3.322524e-01, 3.322019e-01, 3.321221e-01, 3.321050e-01, 3.319118e-01, 3.317922e-01, 3.314658e-01, 3.315735e-01, 3.316331e-01, 3.316525e-01, 3.308030e-01, 3.308038e-01, 3.306947e-01, 3.305741e-01, 3.316492e-01, 3.316117e-01, 3.314973e-01, 3.314110e-01, 3.313450e-01, 3.313649e-01, 3.325841e-01, 3.324226e-01, 3.323649e-01, 3.323381e-01, 3.322566e-01, 3.322077e-01, 3.320860e-01};
  Double_t dVtxPosZ15o[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,5.559279e-01, 3.535446e-01, 4.846955e-01, 4.525585e-01, 3.684501e-01, 2.485494e-01, 2.372653e-01, 1.707859e-01, 3.314213e-01, 1.709195e-01, 2.209753e-01, 3.125757e-01, 3.422085e-01, 3.868156e-01, 4.859695e-01, 4.780697e-01, 4.400149e-01, 4.014992e-01, 3.049883e-01, 3.708501e-01, 3.883566e-01, 3.940632e-01, 4.197670e-01, 3.938399e-01, 3.814413e-01, 3.335539e-01, 3.181929e-01, 2.300734e-01, 2.722395e-01, 5.241033e-01, 3.225908e-01, 1.925791e-01, 1.892765e-01, 3.384066e-01, 2.026459e-01, 2.495699e-01, 3.569992e-01, 3.891381e-01, 4.603724e-01, 3.696685e-01, 3.002207e-01, 2.929533e-01, 3.095468e-01, 3.517200e-01, 2.784445e-01, 3.866626e-01, 3.058719e-01, 3.336752e-01, 3.226473e-01, 3.222815e-01, 3.428469e-01, 3.728514e-01, 2.858642e-01, 2.832485e-01, 3.378933e-01, 3.547548e-01, 3.799414e-01, 4.043543e-01, 4.314049e-01, 4.141138e-01, 3.888746e-01, 4.103586e-01, 3.871045e-01, 4.614473e-01, 4.023404e-01, 4.203531e-01, 4.401272e-01, 6.450558e-01, 6.819582e-01, 2.588529e-01, 3.693471e-01, 3.990708e-01, 3.813842e-01, 3.471682e-01, 3.356156e-01, 2.550150e-01, 3.830723e-01, 4.293259e-01, 4.723797e-01, 4.684324e-01, 4.609304e-01, 4.554974e-01, 4.523016e-01, 3.769890e-01, 4.485548e-01, 5.024484e-01, 5.200088e-01, 5.261731e-01, 5.392851e-01, 5.399264e-01, 5.155504e-01, 4.267668e-01, 5.348764e-01, 4.526746e-01, 4.045626e-01, 4.261759e-01, 5.889205e-01, 6.364843e-01, 5.896163e-01, 3.768637e-01, 4.440771e-01, 4.687029e-01, 4.794467e-01, 4.313422e-01, 3.954777e-01, 3.983129e-01, 3.608064e-01, 2.627038e-01, 3.665826e-01, 4.275667e-01, 3.335445e-01, 3.250815e-01, 3.022907e-01};
  for(Int_t r=0; r<fnRun; r++) {
    fRunList[r] = dRun15o[r];
  }
  fAvVtxPosX=TArrayD(fnRun,dVtxPosX15o);
  fAvVtxPosY=TArrayD(fnRun,dVtxPosY15o);
  fAvVtxPosZ=TArrayD(fnRun,dVtxPosZ15o);
  
  DefineInput(0,TChain::Class());
  DefineOutput(1,AliFlowEventSimple::Class());
}

//=====================================================================

AliAnalysisTaskZDCEP::~AliAnalysisTaskZDCEP()
{
  // Destructor
  delete fOutputList;
  delete fHistList;
  delete fFlowEvent;
  if(fAnalysisUtils) delete fAnalysisUtils;
  if(fMultSelection) delete fMultSelection;
}

//=====================================================================

void AliAnalysisTaskZDCEP::UserCreateOutputObjects ()
{
  //  create  a new  TList  that  OWNS  its  objects
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  
  for(Int_t c=0;c<2;c++) {
    fZDCFlowVect[c] = new AliFlowVector();
    fOutputList->Add(fZDCFlowVect[c]);
  }
  
  fHistList = new TList();
  fHistList->SetOwner(kTRUE);
  
  for(Int_t k=0; k<4; k++) {
    fZDCQHist[k] = new TProfile();
    fHistList->Add(fZDCQHist[k]);
    fZDCVtxHist[k] = new TProfile3D();
    fHistList->Add(fZDCVtxHist[k]);
    fZDCEcomTotHist[k] = new TProfile2D();
    fHistList->Add(fZDCEcomTotHist[k]);
    for(Int_t c=0; c<10; c++) {
      fZDCVtxCenHist[c][k] = new TProfile3D();
      fHistList->Add(fZDCVtxCenHist[c][k]);
    }
    fZDCVtxFitHist[k] = new TH3D();
    fHistList->Add(fZDCVtxFitHist[k]);
    for(Int_t i=0; i<3; i++) {
      fZDCVtxFitCenProjHist[k][i] = new TH1D();
      fHistList->Add(fZDCVtxFitCenProjHist[k][i]);
    }
  }
  for(Int_t c=0; c<10; c++) {
    for(Int_t k=0; k<8; k++) {
      fZDCVtxCenHistMagPol[c][k] = new TProfile3D();
      fHistList->Add(fZDCVtxCenHistMagPol[c][k]);
    }
  }
  for(Int_t i=0; i<10; i++) {
    for(Int_t z=0; z<10; z++) {
      for(Int_t k=0; k<4; k++) {
        fZDCQVecVtxCenEZDC3D[i][z][k] = new TH3D();
        fHistList->Add(fZDCQVecVtxCenEZDC3D[i][z][k]);
      }
    }
  }
  for(Int_t c=0; c<2; c++) {
    for(Int_t i=0; i<5; i++) {
      fTowerGainEq[c][i] = new TH1D();
      fHistList->Add(fTowerGainEq[c][i]);
    }
  }
  
  Double_t DummyEZDCBins[10][11] = {{-3.000000e+02, -4.008000e+01, -2.658000e+01, -1.686000e+01, -8.520000e+00, -7.200000e-01, 7.080000e+00, 1.542000e+01, 2.520000e+01, 3.888000e+01, 3.000000e+02},{-3.000000e+02, -3.690000e+01, -2.436000e+01, -1.530000e+01, -7.560000e+00, -3.000000e-01, 6.960000e+00, 1.476000e+01, 2.388000e+01, 3.666000e+01, 3.000000e+02},{-3.000000e+02, -3.522000e+01, -2.316000e+01, -1.446000e+01, -7.020000e+00, -6.000000e-02, 6.900000e+00, 1.434000e+01, 2.310000e+01, 3.534000e+01, 3.000000e+02},{-3.000000e+02, -3.528000e+01, -2.322000e+01, -1.452000e+01, -7.080000e+00, -1.200000e-01, 6.840000e+00, 1.434000e+01, 2.310000e+01, 3.528000e+01, 3.000000e+02},{-3.000000e+02, -3.666000e+01, -2.412000e+01, -1.506000e+01, -7.320000e+00, -6.000000e-02, 7.200000e+00, 1.500000e+01, 2.412000e+01, 3.684000e+01, 3.000000e+02},{-3.000000e+02, -3.936000e+01, -2.580000e+01, -1.602000e+01, -7.680000e+00, 1.200000e-01, 7.920000e+00, 1.632000e+01, 2.616000e+01, 3.990000e+01, 3.000000e+02},{-3.000000e+02, -4.416000e+01, -2.880000e+01, -1.776000e+01, -8.280000e+00, 5.400000e-01, 9.420000e+00, 1.890000e+01, 3.000000e+01, 4.554000e+01, 3.000000e+02},{-3.000000e+02, -5.262000e+01, -3.384000e+01, -2.028000e+01, -8.700000e+00, 2.100000e+00, 1.296000e+01, 2.454000e+01, 3.816000e+01, 5.712000e+01, 3.000000e+02},{-3.000000e+02, -6.588000e+01, -4.122000e+01, -2.340000e+01, -8.160000e+00, 6.060000e+00, 2.028000e+01, 3.552000e+01, 5.340000e+01, 7.830000e+01, 3.000000e+02},{-3.000000e+02, -8.844000e+01, -5.556000e+01, -3.186000e+01, -1.158000e+01, 7.380000e+00, 2.634000e+01, 4.662000e+01, 7.038000e+01, 1.034400e+02, 3.000000e+02}};
  for (Int_t i=0; i<10; i++) {
    fCRCZDCQVecDummyEZDCBins[i] = new TH1D(Form("fCRCZDCQVecDummyEZDCBins[%d]",i),Form("fCRCZDCQVecDummyEZDCBins[%d]",i),10,DummyEZDCBins[i]);
    fHistList->Add(fCRCZDCQVecDummyEZDCBins[i]);
  }
  
  fAnalysisUtils = new AliAnalysisUtils;
  fAnalysisUtils->SetUseMVPlpSelection(kTRUE);
  fAnalysisUtils->SetUseOutOfBunchPileUp(kTRUE);
  
  fFlowEvent = new AliFlowEvent(1);
}

//=====================================================================

void AliAnalysisTaskZDCEP::UserExec(Option_t *)
{
  //  get  an  event  from  the  analysis  manager
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!aod) return;
  Float_t Centrality = 0.;
  fMultSelection = (AliMultSelection*)aod->FindListObject("MultSelection");
  if(!fMultSelection) {
    AliWarning("WARNING: AliMultSelection object not found ! \n");
  } else {
    Centrality = fMultSelection->GetMultiplicityPercentile("V0M");
  }
  Int_t RunNum = aod->GetRunNumber();
  Int_t RunBin=-1, bin=0;
  for(Int_t c=0;c<fnRun;c++) {
    if(fRunList[c]==RunNum) RunBin=bin;
    else bin++;
  }
  if(RunBin==-1) return;
  Int_t fCenBin = GetCenBin(Centrality);
  
  if(RunNum!=fCachedRunNum) {
    fbFlagIsPosMagField = kFALSE;
    Int_t dRun15hPos[] = {246390, 246391, 246392, 246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424};
    for (Int_t i=0; i<40; i++) {
      if(RunNum==dRun15hPos[i]) fbFlagIsPosMagField = kTRUE;
    }
  }
  
  // get primary vertex position
  Double_t fVtxPos[3]={0.,0.,0.};
  fVtxPos[0] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetX();
  fVtxPos[1] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetY();
  fVtxPos[2] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetZ();
  Double_t fVtxPosCor[3] = {fVtxPos[0]-fAvVtxPosX[RunBin],fVtxPos[1]-fAvVtxPosY[RunBin],fVtxPos[2]-fAvVtxPosZ[RunBin]};
  
  // zdc selection
  AliAODZDC *aodZDC = aod->GetZDCData();
  
  const Double_t * towZNCraw = aodZDC->GetZNCTowerEnergy();
  const Double_t * towZNAraw = aodZDC->GetZNATowerEnergy();
  
  // Get centroid from ZDCs *******************************************************
  
  Double_t Enucl = (RunNum < 209122 ? 1380. : 2511.);
  Double_t xyZNC[2]={0.,0.}, xyZNA[2]={0.,0.};
  Double_t towZNC[5]={0.}, towZNA[5]={0.};
  
  Double_t ZNCcalib=1., ZNAcalib=1.;
  
  // equalize gain of all towers
  if(RunNum!=fCachedRunNum) {
    for(Int_t i=0; i<5; i++) {
      fTowerGainEq[0][i] = (TH1D*)(fTowerEqList->FindObject(Form("fZNCTower[%d][%d]",RunNum,i)));
      fTowerGainEq[1][i] = (TH1D*)(fTowerEqList->FindObject(Form("fZNATower[%d][%d]",RunNum,i)));
    }
  }
  for(Int_t i=0; i<5; i++) {
    if(fTowerGainEq[0][i]) towZNC[i] = towZNCraw[i]*fTowerGainEq[0][i]->GetBinContent(fTowerGainEq[0][i]->FindBin(Centrality));
    if(fTowerGainEq[1][i]) towZNA[i] = towZNAraw[i]*fTowerGainEq[1][i]->GetBinContent(fTowerGainEq[1][i]->FindBin(Centrality));
  }
  
  if(RunNum>=245829) towZNA[2] = 0.;
  Double_t zncEnergy=0., znaEnergy=0.;
  for(Int_t i=0; i<5; i++){
    zncEnergy += towZNC[i];
    znaEnergy += towZNA[i];
  }
  if(RunNum>=245829) znaEnergy *= 8./7.;
  Double_t fZNCen = towZNC[0]/Enucl;
  Double_t fZNAen = towZNA[0]/Enucl;
  
  const Double_t x[4] = {-1.75, 1.75, -1.75, 1.75};
  const Double_t y[4] = {-1.75, -1.75, 1.75, 1.75};
  Double_t numXZNC=0., numYZNC=0., denZNC=0., cZNC, wZNC, EZNC, SumEZNC=0.;
  Double_t numXZNA=0., numYZNA=0., denZNA=0., cZNA, wZNA, EZNA, SumEZNA=0., BadChOr;
  Bool_t fAllChONZNC=kTRUE, fAllChONZNA=kTRUE;
  
  for(Int_t i=0; i<4; i++){
    // get energy
    EZNC = towZNC[i+1];
    SumEZNC += EZNC;
    
    // build centroid
    wZNC = TMath::Power(EZNC, fZDCGainAlpha);
    numXZNC += x[i]*wZNC;
    numYZNC += y[i]*wZNC;
    denZNC += wZNC;
    
    // get energy
    if(i==1) {
      EZNA = towZNA[0]-towZNA[1]-towZNA[3]-towZNA[4];
    } else {
      EZNA = towZNA[i+1];
    }
    SumEZNA += EZNA;
    
    // build centroid
    wZNA = TMath::Power(EZNA, fZDCGainAlpha);
    numXZNA += x[i]*wZNA;
    numYZNA += y[i]*wZNA;
    denZNA += wZNA;
  }
  if(denZNC>0.){
    Double_t nSpecnC = SumEZNC/Enucl;
    cZNC = 1.89358-0.71262/(nSpecnC+0.71789);
    xyZNC[0] = cZNC*numXZNC/denZNC;
    xyZNC[1] = cZNC*numYZNC/denZNC;
    denZNC *= cZNC;
  }
  else{
    xyZNC[0] = xyZNC[1] = 0.;
  }
  if(denZNA>0.){
    Double_t nSpecnA = SumEZNA/Enucl;
    cZNA = 1.89358-0.71262/(nSpecnA+0.71789);
    xyZNA[0] = cZNA*numXZNA/denZNA;
    xyZNA[1] = cZNA*numYZNA/denZNA;
    denZNA *= cZNA;
  }
  else{
    xyZNA[0] = xyZNA[1] = 0.;
  }
  
  fZDCFlowVect[0]->Set(xyZNC[0],xyZNC[1]);
  fZDCFlowVect[1]->Set(xyZNA[0],xyZNA[1]);
  fZDCFlowVect[0]->SetMult(denZNC);
  fZDCFlowVect[1]->SetMult(denZNA);
  
  // RE-CENTER ZDC Q-VECTORS ***************************************************
  
  Int_t qb[4] = {0};
  if(fbFlagIsPosMagField) { qb[0]=0; qb[1]=1; qb[2]=4; qb[3]=5; }
  else                    { qb[0]=2; qb[1]=3; qb[2]=6; qb[3]=7; }
  
  // get re-centered QM*
//  Double_t QMCrec = denZNC;
//  Double_t QMArec = denZNA;
//  if(fAvEZDCCRbRPro && fAvEZDCARbRPro) {
//    Int_t runbin = fAvEZDCCRbRPro->GetXaxis()->FindBin(Form("%d",RunNum));
//    Int_t cenbin = fAvEZDCCRbRPro->GetYaxis()->FindBin(Centrality);
//    QMCrec -= fAvEZDCCRbRPro->GetBinContent(runbin,cenbin);
//    QMArec -= fAvEZDCARbRPro->GetBinContent(runbin,cenbin);
//  }
  
  Bool_t IsGoodEvent = kTRUE;
  
  if(fZDCCalibList) {
    if(RunNum!=fCachedRunNum) {
      // get histos of run
      fZDCQHist[0] = (TProfile*)(fZDCCalibList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fCRCZDCQVecC[%d][%d]",RunNum,0)));
      fZDCQHist[1] = (TProfile*)(fZDCCalibList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fCRCZDCQVecC[%d][%d]",RunNum,1)));
      fZDCQHist[2] = (TProfile*)(fZDCCalibList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fCRCZDCQVecA[%d][%d]",RunNum,0)));
      fZDCQHist[3] = (TProfile*)(fZDCCalibList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fCRCZDCQVecA[%d][%d]",RunNum,1)));
      
      for(Int_t k=0; k<4; k++) {
        fZDCVtxHist[k] = (TProfile3D*)(fZDCCalibList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fCRCZDCQVecVtxPos540[%d][%d]",RunNum,k)));
        
        fZDCEcomTotHist[k] = (TProfile2D*)(fZDCCalibList->FindObject(Form("fCRCZDCQVecEComTot[%d]",k)));
      }
      
      for(Int_t c=0; c<10; c++) {
        for(Int_t k=0; k<4; k++) {
          fZDCVtxCenHist[c][k] = (TProfile3D*)(fZDCCalibList->FindObject(Form("fCRCZDCQVecVtxPosCen[%d][%d]",c,k)));
        }
      }
      
      for(Int_t k=0; k<4; k++) {
        fZDCVtxFitHist[k] = (TH3D*)fZDCCalibList->FindObject(Form("TH3SlopeRunCenVtx[%d]",k));
        if(fZDCVtxFitHist[k]) {
          fZDCVtxFitHist[k]->Sumw2(kFALSE);
          Int_t runbin = fZDCVtxFitHist[k]->GetXaxis()->FindBin(Form("%d",RunNum));
          fZDCVtxFitHist[k]->GetXaxis()->SetRange(runbin,runbin);
          for(Int_t i=0; i<3; i++) {
            fZDCVtxFitHist[k]->GetZaxis()->SetRange(i+1,i+1);
            fZDCVtxFitCenProjHist[k][i] = (TH1D*)fZDCVtxFitHist[k]->Project3D("y")->Clone(Form("proj[%d][%d]",k,i));
          }
        }
      }
      
      for(Int_t c=0; c<10; c++) {
        for(Int_t k=0; k<8; k++) {
          fZDCVtxCenHistMagPol[c][k] = (TProfile3D*)(fZDCCalibList->FindObject(Form("fZDCVtxCenHistMagPol[%d][%d]",c,k)));
        }
      }
      
//      for(Int_t i=0; i<10; i++) {
//        for(Int_t z=0; z<10; z++) {
//          for(Int_t k=0; k<4; k++) {
//            fZDCQVecVtxCenEZDC3D[i][z][k] = (TH3D*)(fZDCCalibList->FindObject(Form("ZDCQVecVtxCenEZDC3D[%d][%d][%d]",i,z,qb[k])));
//          }
//        }
//      }
      
    }
    
    // ZDCN-C
    Double_t QCRe = fZDCFlowVect[0]->X();
    Double_t QCIm = fZDCFlowVect[0]->Y();
    Double_t QMC  = fZDCFlowVect[0]->GetMult();
    // ZDCN-A
    Double_t QARe = fZDCFlowVect[1]->X();
    Double_t QAIm = fZDCFlowVect[1]->Y();
    Double_t QMA  = fZDCFlowVect[1]->GetMult();
    
    Double_t QCReR=QCRe, QCImR=QCIm, QAReR=QARe, QAImR=QAIm;
    
    // STEP #1: re-center vs centrality (1%) vs run number
    
    if (fZDCQHist[0]) {
      Double_t AvQCRe = fZDCQHist[0]->GetBinContent(fZDCQHist[0]->FindBin(Centrality));
      Double_t SDQCRe = fZDCQHist[0]->GetBinError(fZDCQHist[0]->FindBin(Centrality));
      Double_t AvQCIm = fZDCQHist[1]->GetBinContent(fZDCQHist[1]->FindBin(Centrality));
      Double_t SDQCIm = fZDCQHist[1]->GetBinError(fZDCQHist[1]->FindBin(Centrality));
      
      Double_t AvQARe = fZDCQHist[2]->GetBinContent(fZDCQHist[2]->FindBin(Centrality));
      Double_t SDQARe = fZDCQHist[2]->GetBinError(fZDCQHist[2]->FindBin(Centrality));
      Double_t AvQAIm = fZDCQHist[3]->GetBinContent(fZDCQHist[3]->FindBin(Centrality));
      Double_t SDQAIm = fZDCQHist[3]->GetBinError(fZDCQHist[3]->FindBin(Centrality));
      
      if(AvQCRe && AvQCIm && QMC>0. && sqrt(QCRe*QCRe+QCIm*QCIm)>1.E-6) {
        QCReR = QCRe-AvQCRe;
        QCImR = QCIm-AvQCIm;
        fZDCFlowVect[0]->Set(QCReR,QCImR);
      }
      
      if(AvQARe && AvQAIm && QMA>0. && sqrt(QARe*QARe+QAIm*QAIm)>1.E-6) {
        QAReR = QARe-AvQARe;
        QAImR = QAIm-AvQAIm;
        fZDCFlowVect[1]->Set(QAReR,QAImR);
      }
    }
    
    // STEP #2: re-center vs primary vtx vs centrality (10%)
    
    if (fZDCVtxCenHist[fCenBin][0]) {
      Bool_t pass = kTRUE;
      if(fVtxPosCor[0] < fZDCVtxCenHist[fCenBin][0]->GetXaxis()->GetXmin() || fVtxPosCor[0] > fZDCVtxCenHist[fCenBin][0]->GetXaxis()->GetXmax()) pass = kFALSE;
      if(fVtxPosCor[1] < fZDCVtxCenHist[fCenBin][0]->GetYaxis()->GetXmin() || fVtxPosCor[1] > fZDCVtxCenHist[fCenBin][0]->GetYaxis()->GetXmax()) pass = kFALSE;
      if(fVtxPosCor[2] < fZDCVtxCenHist[fCenBin][0]->GetZaxis()->GetXmin() || fVtxPosCor[2] > fZDCVtxCenHist[fCenBin][0]->GetZaxis()->GetXmax()) pass = kFALSE;
      if(!pass) {
        IsGoodEvent = kFALSE;
      } else {
        QCReR -= fZDCVtxCenHist[fCenBin][0]->GetBinContent(fZDCVtxCenHist[fCenBin][0]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
        QCImR -= fZDCVtxCenHist[fCenBin][1]->GetBinContent(fZDCVtxCenHist[fCenBin][1]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
        QAReR -= fZDCVtxCenHist[fCenBin][2]->GetBinContent(fZDCVtxCenHist[fCenBin][2]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
        QAImR -= fZDCVtxCenHist[fCenBin][3]->GetBinContent(fZDCVtxCenHist[fCenBin][3]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
        fZDCFlowVect[0]->Set(QCReR,QCImR);
        fZDCFlowVect[1]->Set(QAReR,QAImR);
      }
    }
    
    // STEP #3: re-center vs primary vtx vs run number
    
    if (fZDCVtxHist[0]) {
      QCReR -= fZDCVtxHist[0]->GetBinContent(fZDCVtxHist[0]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
      QCImR -= fZDCVtxHist[1]->GetBinContent(fZDCVtxHist[1]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
      fZDCFlowVect[0]->Set(QCReR,QCImR);
      QAReR -= fZDCVtxHist[2]->GetBinContent(fZDCVtxHist[2]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
      QAImR -= fZDCVtxHist[3]->GetBinContent(fZDCVtxHist[3]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
      fZDCFlowVect[1]->Set(QAReR,QAImR);
    }
    
    // STEP #4: re-center vs centrality vs energy in the common tower
    
    if(fZDCEcomTotHist[0]) {
      QCReR -= fZDCEcomTotHist[0]->GetBinContent(fZDCEcomTotHist[0]->FindBin(Centrality,fZNCen));
      QCImR -= fZDCEcomTotHist[1]->GetBinContent(fZDCEcomTotHist[1]->FindBin(Centrality,fZNCen));
      fZDCFlowVect[0]->Set(QCReR,QCImR);
      QAReR -= fZDCEcomTotHist[2]->GetBinContent(fZDCEcomTotHist[2]->FindBin(Centrality,fZNAen));
      QAImR -= fZDCEcomTotHist[3]->GetBinContent(fZDCEcomTotHist[3]->FindBin(Centrality,fZNAen));
      fZDCFlowVect[1]->Set(QAReR,QAImR);
    }
    
    // STEP #5: re-center vs vtx vs cen vs run number (through fits)
    
    if(fZDCVtxFitHist[0]) {
      for (Int_t i=0; i<3; i++) {
        QCReR -= fVtxPosCor[i]*fZDCVtxFitCenProjHist[0][i]->Interpolate(Centrality);
        QCImR -= fVtxPosCor[i]*fZDCVtxFitCenProjHist[1][i]->Interpolate(Centrality);
        QAReR -= fVtxPosCor[i]*fZDCVtxFitCenProjHist[2][i]->Interpolate(Centrality);
        QAImR -= fVtxPosCor[i]*fZDCVtxFitCenProjHist[3][i]->Interpolate(Centrality);
      }
      fZDCFlowVect[0]->Set(QCReR,QCImR);
      fZDCFlowVect[1]->Set(QAReR,QAImR);
    }
    
    // second iteration (2D)
    
    if (fZDCVtxCenHistMagPol[fCenBin][0]) {
      if(fbFlagIsPosMagField) {
        QCReR -= fZDCVtxCenHistMagPol[fCenBin][0]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][0]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
        QCImR -= fZDCVtxCenHistMagPol[fCenBin][1]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][1]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
        QAReR -= fZDCVtxCenHistMagPol[fCenBin][4]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][4]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
        QAImR -= fZDCVtxCenHistMagPol[fCenBin][5]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][5]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
      } else {
        QCReR -= fZDCVtxCenHistMagPol[fCenBin][2]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][2]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
        QCImR -= fZDCVtxCenHistMagPol[fCenBin][3]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][3]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
        QAReR -= fZDCVtxCenHistMagPol[fCenBin][6]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][6]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
        QAImR -= fZDCVtxCenHistMagPol[fCenBin][7]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][7]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
      }
      fZDCFlowVect[0]->Set(QCReR,QCImR);
      fZDCFlowVect[1]->Set(QAReR,QAImR);
    }
    
    // STEP #6: re-center vs centrality vs total energy vs vtx
    
//    Int_t EZDCCBin = fCRCZDCQVecDummyEZDCBins[fCenBin]->GetXaxis()->FindBin(QMCrec)-1;
//    Int_t EZDCABin = fCRCZDCQVecDummyEZDCBins[fCenBin]->GetXaxis()->FindBin(QMArec)-1;
//    
//    if(fZDCQVecVtxCenEZDC3D[0][0][0]) {
//      printf("doing step 6 \n");
//      Bool_t pass2=kTRUE;
//      // exclude events with vtx outside of range
//      if(fVtxPosCor[0] < fZDCQVecVtxCenEZDC3D[0][0][0]->GetXaxis()->GetXmin() || fVtxPosCor[0] > fZDCQVecVtxCenEZDC3D[0][0][0]->GetXaxis()->GetXmax()) pass2 = kFALSE;
//      if(fVtxPosCor[1] < fZDCQVecVtxCenEZDC3D[0][0][0]->GetYaxis()->GetXmin() || fVtxPosCor[1] > fZDCQVecVtxCenEZDC3D[0][0][0]->GetYaxis()->GetXmax()) pass2 = kFALSE;
//      if(fVtxPosCor[2] < fZDCQVecVtxCenEZDC3D[0][0][0]->GetZaxis()->GetXmin() || fVtxPosCor[2] > fZDCQVecVtxCenEZDC3D[0][0][0]->GetZaxis()->GetXmax()) pass2 = kFALSE;
//      // exclude events with very low or very high total energy
//      if(fabs(QMCrec)>100. || fabs(QMArec)>100.) pass2 = kFALSE;
//      // get EZDC bin
//      Int_t EZDCCBin = fCRCZDCQVecDummyEZDCBins[fCenBin]->GetXaxis()->FindBin(QMCrec)-1;
//      Int_t EZDCABin = fCRCZDCQVecDummyEZDCBins[fCenBin]->GetXaxis()->FindBin(QMArec)-1;
//      if(EZDCCBin<0) EZDCCBin=0;
//      if(EZDCCBin>9) EZDCCBin=9;
//      if(EZDCABin<0) EZDCABin=0;
//      if(EZDCABin>9) EZDCABin=9;
//      if(pass2) {
//        // check if possible to interpolate
//        Bool_t bInterp = kTRUE;
//        Int_t bx = fZDCQVecVtxCenEZDC3D[0][0][0]->GetXaxis()->FindBin(fVtxPosCor[0]);
//        Int_t by = fZDCQVecVtxCenEZDC3D[0][0][0]->GetYaxis()->FindBin(fVtxPosCor[1]);
//        Int_t bz = fZDCQVecVtxCenEZDC3D[0][0][0]->GetZaxis()->FindBin(fVtxPosCor[2]);
//        if(bx==1 || bx==fZDCQVecVtxCenEZDC3D[0][0][0]->GetXaxis()->GetNbins()) bInterp = kFALSE;
//        if(by==1 || by==fZDCQVecVtxCenEZDC3D[0][0][0]->GetYaxis()->GetNbins()) bInterp = kFALSE;
//        if(bz==1 || bz==fZDCQVecVtxCenEZDC3D[0][0][0]->GetZaxis()->GetNbins()) bInterp = kFALSE;
//        if(bInterp) {
//          QCReR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCCBin][0]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
//          QCImR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCCBin][1]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
//          QAReR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCABin][2]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
//          QAImR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCABin][3]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
//        } else {
//          QCReR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCCBin][0]->GetBinContent(bx,by,bz);
//          QCImR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCCBin][1]->GetBinContent(bx,by,bz);
//          QAReR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCABin][2]->GetBinContent(bx,by,bz);
//          QAImR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCABin][3]->GetBinContent(bx,by,bz);
//        }
//      } else {
//        IsGoodEvent = kFALSE;
//      }
//    }
    
  } else {
    printf("WARNING: no list provided for ZDC Q-vector re-centering ! \n");
  }
  
  Double_t xyZNCfinal[2]={fZDCFlowVect[0]->X(),fZDCFlowVect[0]->Y()};
  Double_t xyZNAfinal[2]={-fZDCFlowVect[1]->X(),fZDCFlowVect[1]->Y()}; // this is not a bug: QAReR --> -QAReR
  if(!IsGoodEvent) {
    xyZNCfinal[0]=0.; xyZNCfinal[1]=0.;
    xyZNAfinal[0]=0.; xyZNAfinal[1]=0.;
  }
  fFlowEvent->SetZDC2Qsub(xyZNCfinal,denZNC,xyZNAfinal,denZNA);
  
  // save run number
  fCachedRunNum = RunNum;
  
  PostData(1, fFlowEvent);
}

//=====================================================================

void AliAnalysisTaskZDCEP::GetZDCQVectors(Double_t QAX, Double_t QAY, Double_t QCX, Double_t QCY)
{
  QAX = fZDCFlowVect[1]->X();
  QAY = fZDCFlowVect[1]->Y();
  QCX = fZDCFlowVect[0]->X();
  QCY = fZDCFlowVect[0]->Y();
}

//=====================================================================

Int_t AliAnalysisTaskZDCEP::GetCenBin(Double_t Centrality)
{
  Int_t CenBin=-1;
  if (Centrality>0. && Centrality<5.) CenBin=0;
  if (Centrality>5. && Centrality<10.) CenBin=1;
  if (Centrality>10. && Centrality<20.) CenBin=2;
  if (Centrality>20. && Centrality<30.) CenBin=3;
  if (Centrality>30. && Centrality<40.) CenBin=4;
  if (Centrality>40. && Centrality<50.) CenBin=5;
  if (Centrality>50. && Centrality<60.) CenBin=6;
  if (Centrality>60. && Centrality<70.) CenBin=7;
  if (Centrality>70. && Centrality<80.) CenBin=8;
  if (Centrality>80. && Centrality<90.) CenBin=9;
  if (Centrality>90. && Centrality<100.) CenBin=10;
  return CenBin;
} // end of AliAnalysisTaskZDCEP::GetCenBin(Double_t Centrality)

//=====================================================================

void AliAnalysisTaskZDCEP::Terminate(Option_t */*option*/)
{
}


