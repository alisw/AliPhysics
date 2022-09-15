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
  fRawPeriod = 0;
  fRawOrbitNumber24 = 0;
  fOrbitNumber = 0;
  fBunchCrossNumber = 0;
  fTimeStamp = 0;
  
  // VZ eta < 0
  fVZCRe = 0;
  fVZCIm = 0;
  fVZCM = 0;
  // VZ eta > 0
  fVZARe = 0;
  fVZAIm = 0;
  fVZAM = 0;
  
  // V0 tow
  fTowV0Craw0 = 0;
  fTowV0Craw1 = 0;
  fTowV0Craw2 = 0;
  fTowV0Craw3 = 0;
  fTowV0Craw4 = 0;
  fTowV0Craw5 = 0;
  fTowV0Craw6 = 0;
  fTowV0Craw7 = 0;
  fTowV0Craw8 = 0;
  fTowV0Craw9 = 0;
  fTowV0Craw10 = 0;
  fTowV0Craw11 = 0;
  fTowV0Craw12 = 0;
  fTowV0Craw13 = 0;
  fTowV0Craw14 = 0;
  fTowV0Craw15 = 0;
  fTowV0Craw16 = 0;
  fTowV0Craw17 = 0;
  fTowV0Craw18 = 0;
  fTowV0Craw19 = 0;
  fTowV0Craw20 = 0;
  fTowV0Craw21 = 0;
  fTowV0Craw22 = 0;
  fTowV0Craw23 = 0;
  fTowV0Craw24 = 0;
  fTowV0Craw25 = 0;
  fTowV0Craw26 = 0;
  fTowV0Craw27 = 0;
  fTowV0Craw28 = 0;
  fTowV0Craw29 = 0;
  fTowV0Craw30 = 0;
  fTowV0Craw31 = 0;

  fTowV0Araw0 = 0;
  fTowV0Araw1 = 0;
  fTowV0Araw2 = 0;
  fTowV0Araw3 = 0;
  fTowV0Araw4 = 0;
  fTowV0Araw5 = 0;
  fTowV0Araw6 = 0;
  fTowV0Araw7 = 0;
  fTowV0Araw8 = 0;
  fTowV0Araw9 = 0;
  fTowV0Araw10 = 0;
  fTowV0Araw11 = 0;
  fTowV0Araw12 = 0;
  fTowV0Araw13 = 0;
  fTowV0Araw14 = 0;
  fTowV0Araw15 = 0;
  fTowV0Araw16 = 0;
  fTowV0Araw17 = 0;
  fTowV0Araw18 = 0;
  fTowV0Araw19 = 0;
  fTowV0Araw20 = 0;
  fTowV0Araw21 = 0;
  fTowV0Araw22 = 0;
  fTowV0Araw23 = 0;
  fTowV0Araw24 = 0;
  fTowV0Araw25 = 0;
  fTowV0Araw26 = 0;
  fTowV0Araw27 = 0;
  fTowV0Araw28 = 0;
  fTowV0Araw29 = 0;
  fTowV0Araw30 = 0;
  fTowV0Araw31 = 0;
  
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
  
  // TPC (Pos Eta 0.1<|eta|<0.8, Neg Eta -0.8<|eta|<-0.1, SubPos Eta 0<|eta|<0.1, SubNeg Eta -0.1<|eta|<0)
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
  
  fTPC4Re2PosChPosEta = 0; // w^2*cos(4phi+) eta+
  fTPC4Im2PosChPosEta = 0; // w^2*sin(4phi+) eta+
  fTPC2Re3PosChPosEta = 0; // w^3*cos(2phi+) eta+
  fTPC2Im3PosChPosEta = 0; // w^3*sin(2phi+) eta+
  fTPC0MPosChPosEta = 0;   // w^0 eta+
  fTPC3MPosChPosEta = 0;   // w^3 eta+
  fTPC4MPosChPosEta = 0;   // w^4 eta+
  fTPC4Re2PosChNegEta = 0; // w^2*cos(4phi+) eta-
  fTPC4Im2PosChNegEta = 0; // w^2*sin(4phi+) eta-
  fTPC2Re3PosChNegEta = 0; // w^3*cos(2phi+) eta-
  fTPC2Im3PosChNegEta = 0; // w^3*sin(2phi+) eta-
  fTPC0MPosChNegEta = 0;   // w^0 eta-
  fTPC3MPosChNegEta = 0;   // w^3 eta-
  fTPC4MPosChNegEta = 0;   // w^4 eta-
  fTPC4Re2NegChPosEta = 0; // w^2*cos(4phi-) eta+
  fTPC4Im2NegChPosEta = 0; // w^2*sin(4phi-) eta+
  fTPC2Re3NegChPosEta = 0; // w^3*cos(2phi-) eta+
  fTPC2Im3NegChPosEta = 0; // w^3*sin(2phi-) eta+
  fTPC0MNegChPosEta = 0;   // w^0 eta+
  fTPC3MNegChPosEta = 0;   // w^3 eta+
  fTPC4MNegChPosEta = 0;   // w^4 eta+
  fTPC4Re2NegChNegEta = 0; // w^2*cos(4phi-) eta-
  fTPC4Im2NegChNegEta = 0; // w^2*sin(4phi-) eta-
  fTPC2Re3NegChNegEta = 0; // w^3*cos(2phi-) eta-
  fTPC2Im3NegChNegEta = 0; // w^3*sin(2phi-) eta-
  fTPC0MNegChNegEta = 0;   // w^0 eta-
  fTPC3MNegChNegEta = 0;   // w^3 eta-
  fTPC4MNegChNegEta = 0;   // w^4 eta-
  
  fTPCRePosChSubPosEta = 0; // w * cos(theta+) eta+
  fTPCImPosChSubPosEta = 0; // w * sin(theta+) eta+
  fTPCMPosChSubPosEta = 0;   // w eta+
  fTPCRePosChSubNegEta = 0; // w * cos(theta+) eta-
  fTPCImPosChSubNegEta = 0; // w * sin(theta+) eta-
  fTPCMPosChSubNegEta = 0;    // w eta-
  fTPCReNegChSubPosEta = 0; // w * cos(theta-) eta+
  fTPCImNegChSubPosEta = 0; // w * sin(theta-) eta+
  fTPCMNegChSubPosEta = 0;    // w eta+
  fTPCReNegChSubNegEta = 0; // w * cos(theta-) eta-
  fTPCImNegChSubNegEta = 0; // w * sin(theta-) eta-
  fTPCMNegChSubNegEta = 0;    // w eta-
  
  fTPC2RePosChSubPosEta = 0; // w * cos(2theta+) eta+
  fTPC2ImPosChSubPosEta = 0; // w * sin(2theta+) eta+
  fTPC2Re2PosChSubPosEta = 0; // w^2 * cos(2theta+) eta+
  fTPC2Im2PosChSubPosEta = 0; // w^2 * sin(2theta+) eta+
  fTPC2MPosChSubPosEta = 0;   // w^2 eta+
  fTPC2RePosChSubNegEta = 0; // w * cos(2theta+) eta-
  fTPC2ImPosChSubNegEta = 0; // w * sin(2theta+) eta-
  fTPC2Re2PosChSubNegEta = 0; // w^2 * cos(2theta+) eta-
  fTPC2Im2PosChSubNegEta = 0; // w^2 * sin(2theta+) eta-
  fTPC2MPosChSubNegEta = 0;    // w^2 eta-
  fTPC2ReNegChSubPosEta = 0; // w * cos(2theta-) eta+
  fTPC2ImNegChSubPosEta = 0; // w * sin(2theta-) eta+
  fTPC2Re2NegChSubPosEta = 0; // w^2 * cos(2theta-) eta+
  fTPC2Im2NegChSubPosEta = 0; // w^2 * sin(2theta-) eta+
  fTPC2MNegChSubPosEta = 0;    // w^2 eta+
  fTPC2ReNegChSubNegEta = 0; // w * cos(2theta-) eta-
  fTPC2ImNegChSubNegEta = 0; // w * sin(2theta-) eta-
  fTPC2Re2NegChSubNegEta = 0; // w^2 * cos(2theta-) eta-
  fTPC2Im2NegChSubNegEta = 0; // w^2 * sin(2theta-) eta-
  fTPC2MNegChSubNegEta = 0;    // w^2 eta-
  
  fTPC4Re2PosChSubPosEta = 0; // w^2*cos(4phi+) eta+
  fTPC4Im2PosChSubPosEta = 0; // w^2*sin(4phi+) eta+
  fTPC2Re3PosChSubPosEta = 0; // w^3*cos(2phi+) eta+
  fTPC2Im3PosChSubPosEta = 0; // w^3*sin(2phi+) eta+
  fTPC0MPosChSubPosEta = 0;   // w^0 eta+
  fTPC3MPosChSubPosEta = 0;   // w^3 eta+
  fTPC4MPosChSubPosEta = 0;   // w^4 eta+
  fTPC4Re2PosChSubNegEta = 0; // w^2*cos(4phi+) eta-
  fTPC4Im2PosChSubNegEta = 0; // w^2*sin(4phi+) eta-
  fTPC2Re3PosChSubNegEta = 0; // w^3*cos(2phi+) eta-
  fTPC2Im3PosChSubNegEta = 0; // w^3*sin(2phi+) eta-
  fTPC0MPosChSubNegEta = 0;   // w^0 eta-
  fTPC3MPosChSubNegEta = 0;   // w^3 eta-
  fTPC4MPosChSubNegEta = 0;   // w^4 eta-
  fTPC4Re2NegChSubPosEta = 0; // w^2*cos(4phi-) eta+
  fTPC4Im2NegChSubPosEta = 0; // w^2*sin(4phi-) eta+
  fTPC2Re3NegChSubPosEta = 0; // w^3*cos(2phi-) eta+
  fTPC2Im3NegChSubPosEta = 0; // w^3*sin(2phi-) eta+
  fTPC0MNegChSubPosEta = 0;   // w^0 eta+
  fTPC3MNegChSubPosEta = 0;   // w^3 eta+
  fTPC4MNegChSubPosEta = 0;   // w^4 eta+
  fTPC4Re2NegChSubNegEta = 0; // w^2*cos(4phi-) eta-
  fTPC4Im2NegChSubNegEta = 0; // w^2*sin(4phi-) eta-
  fTPC2Re3NegChSubNegEta = 0; // w^3*cos(2phi-) eta-
  fTPC2Im3NegChSubNegEta = 0; // w^3*sin(2phi-) eta-
  fTPC0MNegChSubNegEta = 0;   // w^0 eta-
  fTPC3MNegChSubNegEta = 0;   // w^3 eta-
  fTPC4MNegChSubNegEta = 0;   // w^4 eta-
}

AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple::~AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple() {}
