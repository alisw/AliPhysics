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

/*************************************************
 * Qvec event                                    *
 *                                               *
 * author: Shi Qiu                               *
 *         (s.qiu@nikhef.nl)                     *
 *************************************************/

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
#include "TProfile3D.h"
#include "TMath.h"
#include "TArrow.h"
#include "TPaveLabel.h"
#include "TCanvas.h"
#include "AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple.h"
#include "AliLog.h"
#include "TRandom.h"
#include "TF1.h"
#include "TNtuple.h"
#include "THnSparse.h"
#include <complex>
#include <cmath>

ClassImp(AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple)

AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple::AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple()
{
  fRunNum = 0;
  fCentrality = 0;
  fVtxPosX = 0;
  fVtxPosY = 0;
  fVtxPosZ = 0;
  
  // period, orbit number, bunch cross, time stamp
  fOrbitNumber = 0;
  
  // VZ eta < 0
  fVZCRe = 0;
  fVZCIm = 0;
  fVZCM = 0;
  // VZ eta > 0
  fVZARe = 0;
  fVZAIm = 0;
  fVZAM = 0;
  
  // ZNC each tow energy
  fTowZNCraw0 = 0;
  fTowZNCraw1 = 0;
  fTowZNCraw2 = 0;
  fTowZNCraw3 = 0;
  fTowZNCraw4 = 0;
  
  // ZNA each tow energy
  fTowZNAraw0 = 0;
  fTowZNAraw1 = 0;
  fTowZNAraw2 = 0;
  fTowZNAraw3 = 0;
  fTowZNAraw4 = 0;
  
  // TPC 
  fRPReTPC = 0; // Qvec TPC cos(phi) 
  fRPImTPC = 0; // Qvec TPC sin(phi) 
  fRPMultTPC = 0; 
  
  fRP2ReTPC = 0; // Qvec TPC cos(2phi)
  fRP2ImTPC = 0; // Qvec TPC sin(2phi)
  
  fPOIPosReTPC = 0; // Cos(phi_POIPos)
  fPOIPosMult = 0;
  fPOIPosImTPC = 0;

  fPOINegReTPC = 0; // Cos(phi_POINeg)
  fPOINegMult = 0;
  fPOINegImTPC = 0;

  fPOIPos2ReTPC = 0; // Cos(2phi_POIPos)
  fPOIPos2ImTPC = 0;

  fPOINeg2ReTPC = 0; // Cos(2phi_POINeg)
  fPOINeg2ImTPC = 0;

  f2pCorrelatorCos2PsiDiff2PsiV0CRP = 0; //<cos(2psi1-2phi_V0C)>
  f2pCorrelatorCos2PsiDiff2PsiV0ARP = 0; //<cos(2psi1-2phi_V0A)>

  f2pCorrelatorCos2PsiDiff2PsiZDCCRP = 0; //<cos(2psi1-2phi_ZDCC)>
  f2pCorrelatorCos2PsiDiff2PsiZDCARP = 0; //<cos(2psi1-2phi_ZDCA)>
  f2pCorrelatorCos2PsiDiff2PsiZDCCARP = 0; //<cos(2psi1-2phi_ZDCCA)>

	    
  fNITCosPsidiff2PsiV0COS = 0; // <<cos(psi1-2phi_V0C)>> 
  fNITSinPsidiff2PsiV0COS = 0; // <<sin(psi1-2phi_V0C)>>
  fNITCosPsidiff2PsiV0AOS = 0; // <<cos(psi1-2phi_V0A)>>
  fNITSinPsidiff2PsiV0AOS = 0; // <<sin(psi1-2phi_V0A)>>

  fNITCosPsidiff2PsiZDCCOS = 0; // <<cos(psi1-2phi_ZDCC)>>
  fNITSinPsidiff2PsiZDCCOS = 0; // <<sin(psi1-2phi_ZDCC)>>
  fNITCosPsidiff2PsiZDCAOS = 0; // <<cos(psi1-2phi_ZDCA)>>
  fNITSinPsidiff2PsiZDCAOS = 0; // <<sin(psi1-2phi_ZDCA)>>
  fNITCosPsidiff2PsiZDCCAOS = 0; // <<cos(psi1-2phi_ZDCCA)>>
  fNITSinPsidiff2PsiZDCCAOS = 0; // <<sin(psi1-2phi_ZDCCA)>>


  fNITCosPsidiff2PsiV0CPOIPos = 0; // <<cos(psi1-2phi_V0C)>> 
  fNITSinPsidiff2PsiV0CPOIPos = 0; // <<sin(psi1-2phi_V0C)>>
  fNITCosPsidiff2PsiV0APOIPos = 0; // <<cos(psi1-2phi_V0A)>>
  fNITSinPsidiff2PsiV0APOIPos = 0; // <<sin(psi1-2phi_V0A)>>

  fNITCosPsidiff2PsiZDCCPOIPos = 0; // <<cos(psi1-2phi_ZDCC)>>
  fNITSinPsidiff2PsiZDCCPOIPos = 0; // <<sin(psi1-2phi_ZDCC)>>
  fNITCosPsidiff2PsiZDCAPOIPos = 0; // <<cos(psi1-2phi_ZDCA)>>
  fNITSinPsidiff2PsiZDCAPOIPos = 0; // <<sin(psi1-2phi_ZDCA)>>
  fNITCosPsidiff2PsiZDCCAPOIPos = 0; // <<cos(psi1-2phi_ZDCCA)>>
  fNITSinPsidiff2PsiZDCCAPOIPos = 0; // <<sin(psi1-2phi_ZDCCA)>>

  fNITCosPsidiff2PsiV0CPOINeg = 0; // <<cos(psi1-2phi_V0C)>> 
  fNITSinPsidiff2PsiV0CPOINeg = 0; // <<sin(psi1-2phi_V0C)>>
  fNITCosPsidiff2PsiV0APOINeg = 0; // <<cos(psi1-2phi_V0A)>>
  fNITSinPsidiff2PsiV0APOINeg = 0; // <<sin(psi1-2phi_V0A)>>

  fNITCosPsidiff2PsiZDCCPOINeg = 0; // <<cos(psi1-2phi_ZDCC)>>
  fNITSinPsidiff2PsiZDCCPOINeg = 0; // <<sin(psi1-2phi_ZDCC)>>
  fNITCosPsidiff2PsiZDCAPOINeg = 0; // <<cos(psi1-2phi_ZDCA)>>
  fNITSinPsidiff2PsiZDCAPOINeg = 0; // <<sin(psi1-2phi_ZDCA)>>
  fNITCosPsidiff2PsiZDCCAPOINeg = 0; // <<cos(psi1-2phi_ZDCCA)>>
  fNITSinPsidiff2PsiZDCCAPOINeg = 0; // <<sin(psi1-2phi_ZDCCA)>>


  f2pCorrelatorCosPsiDiff = 0; // <cos(dPsi1-dPsi2)>
  f2pCorrelatorCos2PsiDiff = 0; // <cos(2(dPsi1-dPsi2))>
  f2pCorrelatorRPMult = 0;

  f2pCorrelatorCosPsiSumPOIOS = 0; // <cos(dPsi1+dPsi2)>
  f2pCorrelatorPOIOSMult = 0;
  f2pCorrelatorSinPsiSumPOIOS = 0; // <sin(dPsi1+dPsi2)>
  f2pCorrelatorCosPsiDiffPOIOS = 0; // <cos(dPsi1-dPsi2)>
  f2pCorrelatorCos2PsiDiffPOIOS = 0; // <cos(2(dPsi1-dPsi2))>

  f2pCorrelatorCosPsiSumPOIPP = 0; // <cos(dPsi1+dPsi2)>
  f2pCorrelatorPOIPPMult = 0;
  f2pCorrelatorSinPsiSumPOIPP = 0; // <sin(dPsi1+dPsi2)>
  f2pCorrelatorCosPsiDiffPOIPP = 0; // <cos(dPsi1-dPsi2)>
  f2pCorrelatorCos2PsiDiffPOIPP = 0; // <cos(2(dPsi1-dPsi2))>

  f2pCorrelatorCosPsiSumPOINN = 0; // <cos(dPsi1+dPsi2)>
  f2pCorrelatorPOINNMult = 0;
  f2pCorrelatorSinPsiSumPOINN = 0; // <sin(dPsi1+dPsi2)>
  f2pCorrelatorCosPsiDiffPOINN = 0; // <cos(dPsi1-dPsi2)>
  f2pCorrelatorCos2PsiDiffPOINN = 0; // <cos(2(dPsi1-dPsi2))>
}

AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple::~AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple() {}
