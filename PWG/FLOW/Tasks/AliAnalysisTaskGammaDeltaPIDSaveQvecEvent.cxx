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
#include "AliAnalysisTaskGammaDeltaPIDSaveQvecEvent.h"
#include "AliLog.h"
#include "TRandom.h"
#include "TF1.h"
#include "TNtuple.h"
#include "THnSparse.h"
#include <complex>
#include <cmath>

ClassImp(AliAnalysisTaskGammaDeltaPIDSaveQvecEvent)

AliAnalysisTaskGammaDeltaPIDSaveQvecEvent::AliAnalysisTaskGammaDeltaPIDSaveQvecEvent()
{
  fRunNum = 0;
  fCentrality = 0;
  fVtxPosX = 0;
  fVtxPosY = 0;
  fVtxPosZ = 0;
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
  
  // TPC Pion (|eta|<0.8)
  fTPCPionRePosChPosEta = 0; // w * cos(theta+) eta+
  fTPCPionImPosChPosEta = 0; // w * sin(theta+) eta+
  fTPCPionMPosChPosEta = 0;   // w eta+
  fTPCPionRePosChNegEta = 0; // w * cos(theta+) eta-
  fTPCPionImPosChNegEta = 0; // w * sin(theta+) eta-
  fTPCPionMPosChNegEta = 0;    // w eta-
  fTPCPionReNegChPosEta = 0; // w * cos(theta-) eta+
  fTPCPionImNegChPosEta = 0; // w * sin(theta-) eta+
  fTPCPionMNegChPosEta = 0;    // w eta+
  fTPCPionReNegChNegEta = 0; // w * cos(theta-) eta-
  fTPCPionImNegChNegEta = 0; // w * sin(theta-) eta-
  fTPCPionMNegChNegEta = 0;    // w eta-
  
  fTPCPion2RePosChPosEta = 0; // w * cos(2theta+) eta+
  fTPCPion2ImPosChPosEta = 0; // w * sin(2theta+) eta+
  fTPCPion2Re2PosChPosEta = 0; // w^2 * cos(2theta+) eta+
  fTPCPion2Im2PosChPosEta = 0; // w^2 * sin(2theta+) eta+
  fTPCPion2MPosChPosEta = 0;   // w^2 eta+
  fTPCPion2RePosChNegEta = 0; // w * cos(2theta+) eta-
  fTPCPion2ImPosChNegEta = 0; // w * sin(2theta+) eta-
  fTPCPion2Re2PosChNegEta = 0; // w^2 * cos(2theta+) eta-
  fTPCPion2Im2PosChNegEta = 0; // w^2 * sin(2theta+) eta-
  fTPCPion2MPosChNegEta = 0;    // w^2 eta-
  fTPCPion2ReNegChPosEta = 0; // w * cos(2theta-) eta+
  fTPCPion2ImNegChPosEta = 0; // w * sin(2theta-) eta+
  fTPCPion2Re2NegChPosEta = 0; // w^2 * cos(2theta-) eta+
  fTPCPion2Im2NegChPosEta = 0; // w^2 * sin(2theta-) eta+
  fTPCPion2MNegChPosEta = 0;    // w^2 eta+
  fTPCPion2ReNegChNegEta = 0; // w * cos(2theta-) eta-
  fTPCPion2ImNegChNegEta = 0; // w * sin(2theta-) eta-
  fTPCPion2Re2NegChNegEta = 0; // w^2 * cos(2theta-) eta-
  fTPCPion2Im2NegChNegEta = 0; // w^2 * sin(2theta-) eta-
  fTPCPion2MNegChNegEta = 0;    // w^2 eta-
  
  // TPC Kaon (|eta|<0.8)
  fTPCKaonRePosChPosEta = 0; // w * cos(theta+) eta+
  fTPCKaonImPosChPosEta = 0; // w * sin(theta+) eta+
  fTPCKaonMPosChPosEta = 0;   // w eta+
  fTPCKaonRePosChNegEta = 0; // w * cos(theta+) eta-
  fTPCKaonImPosChNegEta = 0; // w * sin(theta+) eta-
  fTPCKaonMPosChNegEta = 0;    // w eta-
  fTPCKaonReNegChPosEta = 0; // w * cos(theta-) eta+
  fTPCKaonImNegChPosEta = 0; // w * sin(theta-) eta+
  fTPCKaonMNegChPosEta = 0;    // w eta+
  fTPCKaonReNegChNegEta = 0; // w * cos(theta-) eta-
  fTPCKaonImNegChNegEta = 0; // w * sin(theta-) eta-
  fTPCKaonMNegChNegEta = 0;    // w eta-
  
  fTPCKaon2RePosChPosEta = 0; // w * cos(2theta+) eta+
  fTPCKaon2ImPosChPosEta = 0; // w * sin(2theta+) eta+
  fTPCKaon2Re2PosChPosEta = 0; // w^2 * cos(2theta+) eta+
  fTPCKaon2Im2PosChPosEta = 0; // w^2 * sin(2theta+) eta+
  fTPCKaon2MPosChPosEta = 0;   // w^2 eta+
  fTPCKaon2RePosChNegEta = 0; // w * cos(2theta+) eta-
  fTPCKaon2ImPosChNegEta = 0; // w * sin(2theta+) eta-
  fTPCKaon2Re2PosChNegEta = 0; // w^2 * cos(2theta+) eta-
  fTPCKaon2Im2PosChNegEta = 0; // w^2 * sin(2theta+) eta-
  fTPCKaon2MPosChNegEta = 0;    // w^2 eta-
  fTPCKaon2ReNegChPosEta = 0; // w * cos(2theta-) eta+
  fTPCKaon2ImNegChPosEta = 0; // w * sin(2theta-) eta+
  fTPCKaon2Re2NegChPosEta = 0; // w^2 * cos(2theta-) eta+
  fTPCKaon2Im2NegChPosEta = 0; // w^2 * sin(2theta-) eta+
  fTPCKaon2MNegChPosEta = 0;    // w^2 eta+
  fTPCKaon2ReNegChNegEta = 0; // w * cos(2theta-) eta-
  fTPCKaon2ImNegChNegEta = 0; // w * sin(2theta-) eta-
  fTPCKaon2Re2NegChNegEta = 0; // w^2 * cos(2theta-) eta-
  fTPCKaon2Im2NegChNegEta = 0; // w^2 * sin(2theta-) eta-
  fTPCKaon2MNegChNegEta = 0;    // w^2 eta-
  
  // TPC Proton (|eta|<0.8)
  fTPCProtonRePosChPosEta = 0; // w * cos(theta+) eta+
  fTPCProtonImPosChPosEta = 0; // w * sin(theta+) eta+
  fTPCProtonMPosChPosEta = 0;   // w eta+
  fTPCProtonRePosChNegEta = 0; // w * cos(theta+) eta-
  fTPCProtonImPosChNegEta = 0; // w * sin(theta+) eta-
  fTPCProtonMPosChNegEta = 0;    // w eta-
  fTPCProtonReNegChPosEta = 0; // w * cos(theta-) eta+
  fTPCProtonImNegChPosEta = 0; // w * sin(theta-) eta+
  fTPCProtonMNegChPosEta = 0;    // w eta+
  fTPCProtonReNegChNegEta = 0; // w * cos(theta-) eta-
  fTPCProtonImNegChNegEta = 0; // w * sin(theta-) eta-
  fTPCProtonMNegChNegEta = 0;    // w eta-
  
  fTPCProton2RePosChPosEta = 0; // w * cos(2theta+) eta+
  fTPCProton2ImPosChPosEta = 0; // w * sin(2theta+) eta+
  fTPCProton2Re2PosChPosEta = 0; // w^2 * cos(2theta+) eta+
  fTPCProton2Im2PosChPosEta = 0; // w^2 * sin(2theta+) eta+
  fTPCProton2MPosChPosEta = 0;   // w^2 eta+
  fTPCProton2RePosChNegEta = 0; // w * cos(2theta+) eta-
  fTPCProton2ImPosChNegEta = 0; // w * sin(2theta+) eta-
  fTPCProton2Re2PosChNegEta = 0; // w^2 * cos(2theta+) eta-
  fTPCProton2Im2PosChNegEta = 0; // w^2 * sin(2theta+) eta-
  fTPCProton2MPosChNegEta = 0;    // w^2 eta-
  fTPCProton2ReNegChPosEta = 0; // w * cos(2theta-) eta+
  fTPCProton2ImNegChPosEta = 0; // w * sin(2theta-) eta+
  fTPCProton2Re2NegChPosEta = 0; // w^2 * cos(2theta-) eta+
  fTPCProton2Im2NegChPosEta = 0; // w^2 * sin(2theta-) eta+
  fTPCProton2MNegChPosEta = 0;    // w^2 eta+
  fTPCProton2ReNegChNegEta = 0; // w * cos(2theta-) eta-
  fTPCProton2ImNegChNegEta = 0; // w * sin(2theta-) eta-
  fTPCProton2Re2NegChNegEta = 0; // w^2 * cos(2theta-) eta-
  fTPCProton2Im2NegChNegEta = 0; // w^2 * sin(2theta-) eta-
  fTPCProton2MNegChNegEta = 0;    // w^2 eta-
	
}

AliAnalysisTaskGammaDeltaPIDSaveQvecEvent::~AliAnalysisTaskGammaDeltaPIDSaveQvecEvent() {}
