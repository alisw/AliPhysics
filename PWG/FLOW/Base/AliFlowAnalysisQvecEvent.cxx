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
#include "AliFlowAnalysisQvecEvent.h"
#include "AliLog.h"
#include "TRandom.h"
#include "TF1.h"
#include "TNtuple.h"
#include "THnSparse.h"
#include <complex>
#include <cmath>

ClassImp(AliFlowAnalysisQvecEvent)

AliFlowAnalysisQvecEvent::AliFlowAnalysisQvecEvent()
{
  fRunNum = 0;
  fCentrality = 0;
  fVtxPosX = 0;
  fVtxPosY = 0;
  fVtxPosZ = 0;
  // ZDC-C (eta < -8.8)
  fZCRe = 0;
  fZCIm = 0;
  fZCM = 0;
  // ZDC-A (eta > 8.8)
  fZARe = 0;
  fZAIm = 0;
  fZAM = 0;
  // VZ eta < 0
  fVZCRe = 0;
  fVZCIm = 0;
  fVZCM = 0;
  // VZ eta > 0
  fVZARe = 0;
  fVZAIm = 0;
  fVZAM = 0;
  
  // TPC (|eta|<0.8)
  fTPCRePosChPosEta = 0; // w * cos(theta+) eta+
  fTPCImPosChPosEta = 0; // w * sin(theta+) eta+
  fTPCMPosChPosEta = 0;   // w eta+
  fTPCRePosChNegEta = 0; // w * cos(theta+) eta-
  fTPCImPosChNegEta = 0; // w * sin(theta+) eta-
  fTPCMPosChNegEta = 0;    // w eta-
  fTPCReNegChPosEta = 0; // w * cos(theta-) eta+
  fTPCImNegChPosEta = 0; // w * sin(theta-) eta+
  fTPCMNegChPosEta = 0;    // w eta+
  fTPCReNegChNegEta = 0; // w * cos(theta-) eta-
  fTPCImNegChNegEta = 0; // w * sin(theta-) eta-
  fTPCMNegChNegEta = 0;    // w eta-
  
  fTPC2RePosChPosEta = 0; // w * cos(2theta+) eta+
  fTPC2ImPosChPosEta = 0; // w * sin(2theta+) eta+
  fTPC2Re2PosChPosEta = 0; // w^2 * cos(2theta+) eta+
  fTPC2Im2PosChPosEta = 0; // w^2 * sin(2theta+) eta+
  fTPC2MPosChPosEta = 0;   // w^2 eta+
  fTPC2RePosChNegEta = 0; // w * cos(2theta+) eta-
  fTPC2ImPosChNegEta = 0; // w * sin(2theta+) eta-
  fTPC2Re2PosChNegEta = 0; // w^2 * cos(2theta+) eta-
  fTPC2Im2PosChNegEta = 0; // w^2 * sin(2theta+) eta-
  fTPC2MPosChNegEta = 0;    // w^2 eta-
  fTPC2ReNegChPosEta = 0; // w * cos(2theta-) eta+
  fTPC2ImNegChPosEta = 0; // w * sin(2theta-) eta+
  fTPC2Re2NegChPosEta = 0; // w^2 * cos(2theta-) eta+
  fTPC2Im2NegChPosEta = 0; // w^2 * sin(2theta-) eta+
  fTPC2MNegChPosEta = 0;    // w^2 eta+
  fTPC2ReNegChNegEta = 0; // w * cos(2theta-) eta-
  fTPC2ImNegChNegEta = 0; // w * sin(2theta-) eta-
  fTPC2Re2NegChNegEta = 0; // w^2 * cos(2theta-) eta-
  fTPC2Im2NegChNegEta = 0; // w^2 * sin(2theta-) eta-
  fTPC2MNegChNegEta = 0;    // w^2 eta-
	
}

AliFlowAnalysisQvecEvent::~AliFlowAnalysisQvecEvent() {}
