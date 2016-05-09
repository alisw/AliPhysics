#include "TRandom.h"
#include "TList.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"

#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliAnalysisUtils.h"
#include "TParticle.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEventHandler.h"

#include "AliAODMCHeader.h"
#include "AliJetContainer.h"
#include "AliPicoTrack.h"
#include "TMath.h"
#include "AliRDHFJetsCuts.h"
#include "AliAnalysisTaskHFJetIPQA.h"
#include "AliExternalTrackParam.h"
#include "AliVertexerTracks.h"
#include "AliHFJetsTagging.h"
#include "AliKFParticle.h"
#include <vector>
#include <algorithm>
using std::min;
using std::cout;
using std::endl;
using std::vector;
using std::pair;
ClassImp(AliAnalysisTaskHFJetIPQA)
// ######################################################################################## CONSTRUCTORS
AliAnalysisTaskHFJetIPQA::AliAnalysisTaskHFJetIPQA(): AliAnalysisTaskEmcalJet("AliAnalysisTaskHFJetIPQA", kTRUE),
fJetCutsHF(new AliRDHFJetsCuts()),
fHFJetUtils(0x0),
fh1dEventRejectionRDHFCuts(0x0),
fh1dVertexZ(0x0),
fh1dVertexZAccepted(0x0),
fh1dVertexR(0x0),
fh1dVertexRAccepted(0x0),
fh1dTracksAccepeted(0x0),
fh1dTracksImpParXY(0x0),
fh1dTracksImpParXYZ(0x0),
fh1dTracksImpParXYSignificance(0x0),
fh1dTracksImpParXYZSignificance(0x0),
fh1dTracksImpParXYTruth(0x0),
fh1dTracksImpParXYZTruth(0x0),
fh1dTracksImpParXYResidualTruth(0x0),
fh1dTracksImpParXYZResidualTruth(0x0),
fh1dTracksImpParXY_McCorr(0x0),
fh1dTracksImpParXYZ_McCorr(0x0),
fh1dTracksImpParXYSignificance_McCorr(0x0),
fh1dTracksImpParXYZSignificance_McCorr(0x0),
fh2dVertexChi2NDFNESDTracks(0x0),
fh1dJetGenPt(0x0),
fh1dJetGenPtUnidentified(0x0),
fh1dJetGenPtudsg(0x0),
fh1dJetGenPtc(0x0),
fh1dJetGenPtb(0x0),
fh1dJetRecPt(0x0),
fh1dJetRecPtAccepted(0x0),
fh1dJetRecEtaPhiAccepted(0x0),
fh1dJetRecPtUnidentified(0x0),
fh1dJetRecPtudsg(0x0),
fh1dJetRecPtc(0x0),
fh1dJetRecPtb(0x0),
fh1dJetRecPtUnidentifiedAccepted(0x0),
fh1dJetRecPtudsgAccepted(0x0),
fh1dJetRecPtcAccepted(0x0),
fh1dJetRecPtbAccepted(0x0),
fh2dJetGenPtVsJetRecPt(0x0),
fh2dJetSignedImpParXY(0x0),
fh2dJetSignedImpParXYUnidentified(0x0),
fh2dJetSignedImpParXYudsg(0x0),
fh2dJetSignedImpParXYb(0x0),
fh2dJetSignedImpParXYbNonBDecay(0x0),
fh2dJetSignedImpParXYc(0x0),
fh2dJetSignedImpParXYSignificance(0x0),
fh2dJetSignedImpParXYSignificanceUnidentified(0x0),
fh2dJetSignedImpParXYSignificanceudsg(0x0),
fh2dJetSignedImpParXYSignificanceb(0x0),
fh2dJetSignedImpParXYSignificancebNonBDecay(0x0),
fh2dJetSignedImpParXYSignificancec(0x0),
fh2dJetSignedImpParXYZ(0x0),
fh2dJetSignedImpParXYZUnidentified(0x0),
fh2dJetSignedImpParXYZudsg(0x0),
fh2dJetSignedImpParXYZb(0x0),
fh2dJetSignedImpParXYZbNonBDecay(0x0),
fh2dJetSignedImpParXYZc(0x0),
fh2dJetSignedImpParXYZSignificance(0x0),
fh2dJetSignedImpParXYZSignificanceUnidentified(0x0),
fh2dJetSignedImpParXYZSignificanceudsg(0x0),
fh2dJetSignedImpParXYZSignificanceb(0x0),
fh2dJetSignedImpParXYZSignificancebNonBDecay(0x0),
fh2dJetSignedImpParXYZSignificancec(0x0),
fh2dJetSignedImpParXY_McCorr(0x0),
fh2dJetSignedImpParXYUnidentified_McCorr(0x0),
fh2dJetSignedImpParXYudsg_McCorr(0x0),
fh2dJetSignedImpParXYb_McCorr(0x0),
fh2dJetSignedImpParXYb_McCorrNonBDecay(0x0),
fh2dJetSignedImpParXYc_McCorr(0x0),
fh2dJetSignedImpParXYSignificance_McCorr(0x0),
fh2dJetSignedImpParXYSignificanceUnidentified_McCorr(0x0),
fh2dJetSignedImpParXYSignificanceudsg_McCorr(0x0),
fh2dJetSignedImpParXYSignificanceb_McCorr(0x0),
fh2dJetSignedImpParXYSignificanceb_McCorrNonBDecay(0x0),
fh2dJetSignedImpParXYSignificancec_McCorr(0x0),
fh2dJetSignedImpParXYZ_McCorr(0x0),
fh2dJetSignedImpParXYZUnidentified_McCorr(0x0),
fh2dJetSignedImpParXYZudsg_McCorr(0x0),
fh2dJetSignedImpParXYZb_McCorr(0x0),
fh2dJetSignedImpParXYZb_McCorr_McCorrNonBDecay(0x0),
fh2dJetSignedImpParXYZc_McCorr(0x0),
fh2dJetSignedImpParXYZSignificance_McCorr(0x0),
fh2dJetSignedImpParXYZSignificanceUnidentified_McCorr(0x0),
fh2dJetSignedImpParXYZSignificanceudsg_McCorr(0x0),
fh2dJetSignedImpParXYZSignificanceb_McCorr(0x0),
fh2dJetSignedImpParXYZSignificanceb_McCorrNonBDecay(0x0),
fh2dJetSignedImpParXYZSignificancec_McCorr(0x0),
fh2dJetSignedImpParXYFirst(0x0),
fh2dJetSignedImpParXYUnidentifiedFirst(0x0),
fh2dJetSignedImpParXYudsgFirst(0x0),
fh2dJetSignedImpParXYbFirst(0x0),
fh2dJetSignedImpParXYbFirstNonBDecay(0x0),
fh2dJetSignedImpParXYcFirst(0x0),
fh2dJetSignedImpParXYSignificanceFirst(0x0),
fh2dJetSignedImpParXYSignificanceUnidentifiedFirst(0x0),
fh2dJetSignedImpParXYSignificanceudsgFirst(0x0),
fh2dJetSignedImpParXYSignificancebFirst(0x0),
fh2dJetSignedImpParXYSignificancebFirstNonBDecay(0x0),
fh2dJetSignedImpParXYSignificancecFirst(0x0),
fh2dJetSignedImpParXYZFirst(0x0),
fh2dJetSignedImpParXYZUnidentifiedFirst(0x0),
fh2dJetSignedImpParXYZudsgFirst(0x0),
fh2dJetSignedImpParXYZbFirst(0x0),
fh2dJetSignedImpParXYZbFirstNonBDecay(0x0),
fh2dJetSignedImpParXYZcFirst(0x0),
fh2dJetSignedImpParXYZSignificanceFirst(0x0),
fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst(0x0),
fh2dJetSignedImpParXYZSignificanceudsgFirst(0x0),
fh2dJetSignedImpParXYZSignificancebFirst(0x0),
fh2dJetSignedImpParXYZSignificancebFirstNonBDecay(0x0),
fh2dJetSignedImpParXYZSignificancecFirst(0x0),
fh2dJetSignedImpParXYSecond(0x0),
fh2dJetSignedImpParXYUnidentifiedSecond(0x0),
fh2dJetSignedImpParXYudsgSecond(0x0),
fh2dJetSignedImpParXYbSecond(0x0),
fh2dJetSignedImpParXYbSecondNonBDecay(0x0),
fh2dJetSignedImpParXYcSecond(0x0),
fh2dJetSignedImpParXYSignificanceSecond(0x0),
fh2dJetSignedImpParXYSignificanceUnidentifiedSecond(0x0),
fh2dJetSignedImpParXYSignificanceudsgSecond(0x0),
fh2dJetSignedImpParXYSignificancebSecond(0x0),
fh2dJetSignedImpParXYSignificancebSecondNonBDecay(0x0),
fh2dJetSignedImpParXYSignificancecSecond(0x0),
fh2dJetSignedImpParXYZSecond(0x0),
fh2dJetSignedImpParXYZUnidentifiedSecond(0x0),
fh2dJetSignedImpParXYZudsgSecond(0x0),
fh2dJetSignedImpParXYZbSecond(0x0),
fh2dJetSignedImpParXYZbSecondNonBDecay(0x0),
fh2dJetSignedImpParXYZcSecond(0x0),
fh2dJetSignedImpParXYZSignificanceSecond(0x0),
fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond(0x0),
fh2dJetSignedImpParXYZSignificanceudsgSecond(0x0),
fh2dJetSignedImpParXYZSignificancebSecond(0x0),
fh2dJetSignedImpParXYZSignificancebSecondNonBDecay(0x0),
fh2dJetSignedImpParXYZSignificancecSecond(0x0),
fh2dJetSignedImpParXYThird(0x0),
fh2dJetSignedImpParXYUnidentifiedThird(0x0),
fh2dJetSignedImpParXYudsgThird(0x0),
fh2dJetSignedImpParXYbThird(0x0),
fh2dJetSignedImpParXYbThirdNonBDecay(0x0),
fh2dJetSignedImpParXYcThird(0x0),
fh2dJetSignedImpParXYSignificanceThird(0x0),
fh2dJetSignedImpParXYSignificanceUnidentifiedThird(0x0),
fh2dJetSignedImpParXYSignificanceudsgThird(0x0),
fh2dJetSignedImpParXYSignificancebThird(0x0),
fh2dJetSignedImpParXYSignificancebThirdNonBDecay(0x0),
fh2dJetSignedImpParXYSignificancecThird(0x0),
fh2dJetSignedImpParXYZThird(0x0),
fh2dJetSignedImpParXYZUnidentifiedThird(0x0),
fh2dJetSignedImpParXYZudsgThird(0x0),
fh2dJetSignedImpParXYZbThird(0x0),
fh2dJetSignedImpParXYZbThirdNonBDecay(0x0),
fh2dJetSignedImpParXYZcThird(0x0),
fh2dJetSignedImpParXYZSignificanceThird(0x0),
fh2dJetSignedImpParXYZSignificanceUnidentifiedThird(0x0),
fh2dJetSignedImpParXYZSignificanceudsgThird(0x0),
fh2dJetSignedImpParXYZSignificancebThird(0x0),
fh2dJetSignedImpParXYZSignificancebThirdNonBDecay(0x0),
fh2dJetSignedImpParXYZSignificancecThird(0x0),
fh2dJetSignedImpParXYFirst_McCorr(0x0),
fh2dJetSignedImpParXYUnidentifiedFirst_McCorr(0x0),
fh2dJetSignedImpParXYudsgFirst_McCorr(0x0),
fh2dJetSignedImpParXYbFirst_McCorr(0x0),
fh2dJetSignedImpParXYbFirst_McCorrNonBDecay(0x0),
fh2dJetSignedImpParXYcFirst_McCorr(0x0),
fh2dJetSignedImpParXYSignificanceFirst_McCorr(0x0),
fh2dJetSignedImpParXYSignificanceUnidentifiedFirst_McCorr(0x0),
fh2dJetSignedImpParXYSignificanceudsgFirst_McCorr(0x0),
fh2dJetSignedImpParXYSignificancebFirst_McCorr(0x0),
fh2dJetSignedImpParXYSignificancebFirst_McCorrNonBDecay(0x0),
fh2dJetSignedImpParXYSignificancecFirst_McCorr(0x0),
fh2dJetSignedImpParXYZFirst_McCorr(0x0),
fh2dJetSignedImpParXYZUnidentifiedFirst_McCorr(0x0),
fh2dJetSignedImpParXYZudsgFirst_McCorr(0x0),
fh2dJetSignedImpParXYZbFirst_McCorr(0x0),
fh2dJetSignedImpParXYZbFirst_McCorrNonBDecay(0x0),
fh2dJetSignedImpParXYZcFirst_McCorr(0x0),
fh2dJetSignedImpParXYZSignificanceFirst_McCorr(0x0),
fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst_McCorr(0x0),
fh2dJetSignedImpParXYZSignificanceudsgFirst_McCorr(0x0),
fh2dJetSignedImpParXYZSignificancebFirst_McCorr(0x0),
fh2dJetSignedImpParXYZSignificancebFirst_McCorrNonBDecay(0x0),
fh2dJetSignedImpParXYZSignificancecFirst_McCorr(0x0),
fh2dJetSignedImpParXYSecond_McCorr(0x0),
fh2dJetSignedImpParXYUnidentifiedSecond_McCorr(0x0),
fh2dJetSignedImpParXYudsgSecond_McCorr(0x0),
fh2dJetSignedImpParXYbSecond_McCorr(0x0),
fh2dJetSignedImpParXYbSecond_McCorrNonBDecay(0x0),
fh2dJetSignedImpParXYcSecond_McCorr(0x0),
fh2dJetSignedImpParXYSignificanceSecond_McCorr(0x0),
fh2dJetSignedImpParXYSignificanceUnidentifiedSecond_McCorr(0x0),
fh2dJetSignedImpParXYSignificanceudsgSecond_McCorr(0x0),
fh2dJetSignedImpParXYSignificancebSecond_McCorr(0x0),
fh2dJetSignedImpParXYSignificancebSecond_McCorrNonBDecay(0x0),
fh2dJetSignedImpParXYSignificancecSecond_McCorr(0x0),
fh2dJetSignedImpParXYZSecond_McCorr(0x0),
fh2dJetSignedImpParXYZUnidentifiedSecond_McCorr(0x0),
fh2dJetSignedImpParXYZudsgSecond_McCorr(0x0),
fh2dJetSignedImpParXYZbSecond_McCorrNonBDecay(0x0),
fh2dJetSignedImpParXYZcSecond_McCorr(0x0),
fh2dJetSignedImpParXYZSignificanceSecond_McCorr(0x0),
fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond_McCorr(0x0),
fh2dJetSignedImpParXYZSignificanceudsgSecond_McCorr(0x0),
fh2dJetSignedImpParXYZSignificancebSecond_McCorr(0x0),
fh2dJetSignedImpParXYZSignificancebSecond_McCorrNonBDecay(0x0),
fh2dJetSignedImpParXYZSignificancecSecond_McCorr(0x0),
fh2dJetSignedImpParXYThird_McCorr(0x0),
fh2dJetSignedImpParXYUnidentifiedThird_McCorr(0x0),
fh2dJetSignedImpParXYudsgThird_McCorr(0x0),
fh2dJetSignedImpParXYbThird_McCorr(0x0),
fh2dJetSignedImpParXYbThird_McCorrNonBDecay(0x0),
fh2dJetSignedImpParXYcThird_McCorr(0x0),
fh2dJetSignedImpParXYSignificanceThird_McCorr(0x0),
fh2dJetSignedImpParXYSignificanceUnidentifiedThird_McCorr(0x0),
fh2dJetSignedImpParXYSignificanceudsgThird_McCorr(0x0),
fh2dJetSignedImpParXYSignificancebThird_McCorr(0x0),
fh2dJetSignedImpParXYSignificancebThird_McCorrNonBDecay(0x0),
fh2dJetSignedImpParXYSignificancecThird_McCorr(0x0),
fh2dJetSignedImpParXYZThird_McCorr(0x0),
fh2dJetSignedImpParXYZUnidentifiedThird_McCorr(0x0),
fh2dJetSignedImpParXYZudsgThird_McCorr(0x0),
fh2dJetSignedImpParXYZbThird_McCorr(0x0),
fh2dJetSignedImpParXYZbThird_McCorrNonBDecay(0x0),
fh2dJetSignedImpParXYZcThird_McCorr(0x0),
fh2dJetSignedImpParXYZSignificanceThird_McCorr(0x0),
fh2dJetSignedImpParXYZSignificanceUnidentifiedThird_McCorr(0x0),
fh2dJetSignedImpParXYZSignificanceudsgThird_McCorr(0x0),
fh2dJetSignedImpParXYZSignificancebThird_McCorr(0x0),
fh2dJetSignedImpParXYZSignificancebThird_McCorrNonBDecay(0x0),
fh2dJetSignedImpParXYZSignificancecThird_McCorr(0x0),
fh1dTracksImpParXY_electron(0x0),
fh1dTracksImpParXYZ_electron(0x0),
fh1dTracksImpParXYSignificance_electron(0x0),
fh1dTracksImpParXYZSignificance_electron(0x0),
//electrons
fh2dJetSignedImpParXY_electron(0x0),
fh2dJetSignedImpParXYUnidentified_electron(0x0),
fh2dJetSignedImpParXYudsg_electron(0x0),
fh2dJetSignedImpParXYb_electron(0x0),
fh2dJetSignedImpParXYc_electron(0x0),
fh2dJetSignedImpParXYSignificance_electron(0x0),
fh2dJetSignedImpParXYSignificanceUnidentified_electron(0x0),
fh2dJetSignedImpParXYSignificanceudsg_electron(0x0),
fh2dJetSignedImpParXYSignificanceb_electron(0x0),
fh2dJetSignedImpParXYSignificancec_electron(0x0),
fh2dJetSignedImpParXYZ_electron(0x0),
fh2dJetSignedImpParXYZUnidentified_electron(0x0),
fh2dJetSignedImpParXYZudsg_electron(0x0),
fh2dJetSignedImpParXYZb_electron(0x0),
fh2dJetSignedImpParXYZc_electron(0x0),
fh2dJetSignedImpParXYZSignificance_electron(0x0),
fh2dJetSignedImpParXYZSignificanceUnidentified_electron(0x0),
fh2dJetSignedImpParXYZSignificanceudsg_electron(0x0),
fh2dJetSignedImpParXYZSignificanceb_electron(0x0),
fh2dJetSignedImpParXYZSignificancec_electron(0x0),
fh2dJetSignedImpParXYFirst_electron(0x0),
fh2dJetSignedImpParXYUnidentifiedFirst_electron(0x0),
fh2dJetSignedImpParXYudsgFirst_electron(0x0),
fh2dJetSignedImpParXYbFirst_electron(0x0),
fh2dJetSignedImpParXYcFirst_electron(0x0),
fh2dJetSignedImpParXYSignificanceFirst_electron(0x0),
fh2dJetSignedImpParXYSignificanceUnidentifiedFirst_electron(0x0),
fh2dJetSignedImpParXYSignificanceudsgFirst_electron(0x0),
fh2dJetSignedImpParXYSignificancebFirst_electron(0x0),
fh2dJetSignedImpParXYSignificancecFirst_electron(0x0),
fh2dJetSignedImpParXYZFirst_electron(0x0),
fh2dJetSignedImpParXYZUnidentifiedFirst_electron(0x0),
fh2dJetSignedImpParXYZudsgFirst_electron(0x0),
fh2dJetSignedImpParXYZbFirst_electron(0x0),
fh2dJetSignedImpParXYZcFirst_electron(0x0),
fh2dJetSignedImpParXYZSignificanceFirst_electron(0x0),
fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst_electron(0x0),
fh2dJetSignedImpParXYZSignificanceudsgFirst_electron(0x0),
fh2dJetSignedImpParXYZSignificancebFirst_electron(0x0),
fh2dJetSignedImpParXYZSignificancecFirst_electron(0x0),
fh2dJetSignedImpParXYSecond_electron(0x0),
fh2dJetSignedImpParXYUnidentifiedSecond_electron(0x0),
fh2dJetSignedImpParXYudsgSecond_electron(0x0),
fh2dJetSignedImpParXYbSecond_electron(0x0),
fh2dJetSignedImpParXYcSecond_electron(0x0),
fh2dJetSignedImpParXYSignificanceSecond_electron(0x0),
fh2dJetSignedImpParXYSignificanceUnidentifiedSecond_electron(0x0),
fh2dJetSignedImpParXYSignificanceudsgSecond_electron(0x0),
fh2dJetSignedImpParXYSignificancebSecond_electron(0x0),
fh2dJetSignedImpParXYSignificancecSecond_electron(0x0),
fh2dJetSignedImpParXYZSecond_electron(0x0),
fh2dJetSignedImpParXYZUnidentifiedSecond_electron(0x0),
fh2dJetSignedImpParXYZudsgSecond_electron(0x0),
fh2dJetSignedImpParXYZbSecond_electron(0x0),
fh2dJetSignedImpParXYZcSecond_electron(0x0),
fh2dJetSignedImpParXYZSignificanceSecond_electron(0x0),
fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond_electron(0x0),
fh2dJetSignedImpParXYZSignificanceudsgSecond_electron(0x0),
fh2dJetSignedImpParXYZSignificancebSecond_electron(0x0),
fh2dJetSignedImpParXYZSignificancecSecond_electron(0x0),
fh2dJetSignedImpParXYThird_electron(0x0),
fh2dJetSignedImpParXYUnidentifiedThird_electron(0x0),
fh2dJetSignedImpParXYudsgThird_electron(0x0),
fh2dJetSignedImpParXYbThird_electron(0x0),
fh2dJetSignedImpParXYcThird_electron(0x0),
fh2dJetSignedImpParXYSignificanceThird_electron(0x0),
fh2dJetSignedImpParXYSignificanceUnidentifiedThird_electron(0x0),
fh2dJetSignedImpParXYSignificanceudsgThird_electron(0x0),
fh2dJetSignedImpParXYSignificancebThird_electron(0x0),
fh2dJetSignedImpParXYSignificancecThird_electron(0x0),
fh2dJetSignedImpParXYZThird_electron(0x0),
fh2dJetSignedImpParXYZUnidentifiedThird_electron(0x0),
fh2dJetSignedImpParXYZudsgThird_electron(0x0),
fh2dJetSignedImpParXYZbThird_electron(0x0),
fh2dJetSignedImpParXYZcThird_electron(0x0),
fh2dJetSignedImpParXYZSignificanceThird_electron(0x0),
fh2dJetSignedImpParXYZSignificanceUnidentifiedThird_electron(0x0),
fh2dJetSignedImpParXYZSignificanceudsgThird_electron(0x0),
fh2dJetSignedImpParXYZSignificancebThird_electron(0x0),
fh2dJetSignedImpParXYZSignificancecThird_electron(0x0),
fMCArray(0x0),
fMCEvent(0x0),
fESDTrackCut(0x0),
fUtils(new AliAnalysisUtils())
{
	for(int i =0 ; i<498;++i){
		for(int j =0 ; j<16;++j){
			fBackgroundFactorLinus[j][i]=1; //set default to 1
		}}
	fEnableV0GammaRejection = kFALSE;
}
// ######################################################################################## CONSTRUCTORS
AliAnalysisTaskHFJetIPQA::AliAnalysisTaskHFJetIPQA(const char *name): AliAnalysisTaskEmcalJet(name, kTRUE),fJetCutsHF(new AliRDHFJetsCuts()),
		fHFJetUtils(0x0),
		fh1dEventRejectionRDHFCuts(0x0),
		fh1dVertexZ(0x0),
		fh1dVertexZAccepted(0x0),
		fh1dVertexR(0x0),
		fh1dVertexRAccepted(0x0),
		fh1dTracksAccepeted(0x0),
		fh1dTracksImpParXY(0x0),
		fh1dTracksImpParXYZ(0x0),
		fh1dTracksImpParXYSignificance(0x0),
		fh1dTracksImpParXYZSignificance(0x0),
		fh1dTracksImpParXYTruth(0x0),
		fh1dTracksImpParXYZTruth(0x0),
		fh1dTracksImpParXYResidualTruth(0x0),
		fh1dTracksImpParXYZResidualTruth(0x0),
		fh1dTracksImpParXY_McCorr(0x0),
		fh1dTracksImpParXYZ_McCorr(0x0),
		fh1dTracksImpParXYSignificance_McCorr(0x0),
		fh1dTracksImpParXYZSignificance_McCorr(0x0),
		fh2dVertexChi2NDFNESDTracks(0x0),
		fh1dJetGenPt(0x0),
		fh1dJetGenPtUnidentified(0x0),
		fh1dJetGenPtudsg(0x0),
		fh1dJetGenPtc(0x0),
		fh1dJetGenPtb(0x0),
		fh1dJetRecPt(0x0),
		fh1dJetRecPtAccepted(0x0),
		fh1dJetRecEtaPhiAccepted(0x0),
		fh1dJetRecPtUnidentified(0x0),
		fh1dJetRecPtudsg(0x0),
		fh1dJetRecPtc(0x0),
		fh1dJetRecPtb(0x0),
		fh1dJetRecPtUnidentifiedAccepted(0x0),
		fh1dJetRecPtudsgAccepted(0x0),
		fh1dJetRecPtcAccepted(0x0),
		fh1dJetRecPtbAccepted(0x0),
		fh2dJetGenPtVsJetRecPt(0x0),
		fh2dJetSignedImpParXY(0x0),
		fh2dJetSignedImpParXYUnidentified(0x0),
		fh2dJetSignedImpParXYudsg(0x0),
		fh2dJetSignedImpParXYb(0x0),
		fh2dJetSignedImpParXYbNonBDecay(0x0),
		fh2dJetSignedImpParXYc(0x0),
		fh2dJetSignedImpParXYSignificance(0x0),
		fh2dJetSignedImpParXYSignificanceUnidentified(0x0),
		fh2dJetSignedImpParXYSignificanceudsg(0x0),
		fh2dJetSignedImpParXYSignificanceb(0x0),
		fh2dJetSignedImpParXYSignificancebNonBDecay(0x0),
		fh2dJetSignedImpParXYSignificancec(0x0),
		fh2dJetSignedImpParXYZ(0x0),
		fh2dJetSignedImpParXYZUnidentified(0x0),
		fh2dJetSignedImpParXYZudsg(0x0),
		fh2dJetSignedImpParXYZb(0x0),
		fh2dJetSignedImpParXYZbNonBDecay(0x0),
		fh2dJetSignedImpParXYZc(0x0),
		fh2dJetSignedImpParXYZSignificance(0x0),
		fh2dJetSignedImpParXYZSignificanceUnidentified(0x0),
		fh2dJetSignedImpParXYZSignificanceudsg(0x0),
		fh2dJetSignedImpParXYZSignificanceb(0x0),
		fh2dJetSignedImpParXYZSignificancebNonBDecay(0x0),
		fh2dJetSignedImpParXYZSignificancec(0x0),
		fh2dJetSignedImpParXY_McCorr(0x0),
		fh2dJetSignedImpParXYUnidentified_McCorr(0x0),
		fh2dJetSignedImpParXYudsg_McCorr(0x0),
		fh2dJetSignedImpParXYb_McCorr(0x0),
		fh2dJetSignedImpParXYb_McCorrNonBDecay(0x0),
		fh2dJetSignedImpParXYc_McCorr(0x0),
		fh2dJetSignedImpParXYSignificance_McCorr(0x0),
		fh2dJetSignedImpParXYSignificanceUnidentified_McCorr(0x0),
		fh2dJetSignedImpParXYSignificanceudsg_McCorr(0x0),
		fh2dJetSignedImpParXYSignificanceb_McCorr(0x0),
		fh2dJetSignedImpParXYSignificanceb_McCorrNonBDecay(0x0),
		fh2dJetSignedImpParXYSignificancec_McCorr(0x0),
		fh2dJetSignedImpParXYZ_McCorr(0x0),
		fh2dJetSignedImpParXYZUnidentified_McCorr(0x0),
		fh2dJetSignedImpParXYZudsg_McCorr(0x0),
		fh2dJetSignedImpParXYZb_McCorr(0x0),
		fh2dJetSignedImpParXYZb_McCorr_McCorrNonBDecay(0x0),
		fh2dJetSignedImpParXYZc_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificance_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificanceUnidentified_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificanceudsg_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificanceb_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificanceb_McCorrNonBDecay(0x0),
		fh2dJetSignedImpParXYZSignificancec_McCorr(0x0),
		fh2dJetSignedImpParXYFirst(0x0),
		fh2dJetSignedImpParXYUnidentifiedFirst(0x0),
		fh2dJetSignedImpParXYudsgFirst(0x0),
		fh2dJetSignedImpParXYbFirst(0x0),
		fh2dJetSignedImpParXYbFirstNonBDecay(0x0),
		fh2dJetSignedImpParXYcFirst(0x0),
		fh2dJetSignedImpParXYSignificanceFirst(0x0),
		fh2dJetSignedImpParXYSignificanceUnidentifiedFirst(0x0),
		fh2dJetSignedImpParXYSignificanceudsgFirst(0x0),
		fh2dJetSignedImpParXYSignificancebFirst(0x0),
		fh2dJetSignedImpParXYSignificancebFirstNonBDecay(0x0),
		fh2dJetSignedImpParXYSignificancecFirst(0x0),
		fh2dJetSignedImpParXYZFirst(0x0),
		fh2dJetSignedImpParXYZUnidentifiedFirst(0x0),
		fh2dJetSignedImpParXYZudsgFirst(0x0),
		fh2dJetSignedImpParXYZbFirst(0x0),
		fh2dJetSignedImpParXYZbFirstNonBDecay(0x0),
		fh2dJetSignedImpParXYZcFirst(0x0),
		fh2dJetSignedImpParXYZSignificanceFirst(0x0),
		fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst(0x0),
		fh2dJetSignedImpParXYZSignificanceudsgFirst(0x0),
		fh2dJetSignedImpParXYZSignificancebFirst(0x0),
		fh2dJetSignedImpParXYZSignificancebFirstNonBDecay(0x0),
		fh2dJetSignedImpParXYZSignificancecFirst(0x0),
		fh2dJetSignedImpParXYSecond(0x0),
		fh2dJetSignedImpParXYUnidentifiedSecond(0x0),
		fh2dJetSignedImpParXYudsgSecond(0x0),
		fh2dJetSignedImpParXYbSecond(0x0),
		fh2dJetSignedImpParXYbSecondNonBDecay(0x0),
		fh2dJetSignedImpParXYcSecond(0x0),
		fh2dJetSignedImpParXYSignificanceSecond(0x0),
		fh2dJetSignedImpParXYSignificanceUnidentifiedSecond(0x0),
		fh2dJetSignedImpParXYSignificanceudsgSecond(0x0),
		fh2dJetSignedImpParXYSignificancebSecond(0x0),
		fh2dJetSignedImpParXYSignificancebSecondNonBDecay(0x0),
		fh2dJetSignedImpParXYSignificancecSecond(0x0),
		fh2dJetSignedImpParXYZSecond(0x0),
		fh2dJetSignedImpParXYZUnidentifiedSecond(0x0),
		fh2dJetSignedImpParXYZudsgSecond(0x0),
		fh2dJetSignedImpParXYZbSecond(0x0),
		fh2dJetSignedImpParXYZbSecondNonBDecay(0x0),
		fh2dJetSignedImpParXYZcSecond(0x0),
		fh2dJetSignedImpParXYZSignificanceSecond(0x0),
		fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond(0x0),
		fh2dJetSignedImpParXYZSignificanceudsgSecond(0x0),
		fh2dJetSignedImpParXYZSignificancebSecond(0x0),
		fh2dJetSignedImpParXYZSignificancebSecondNonBDecay(0x0),
		fh2dJetSignedImpParXYZSignificancecSecond(0x0),
		fh2dJetSignedImpParXYThird(0x0),
		fh2dJetSignedImpParXYUnidentifiedThird(0x0),
		fh2dJetSignedImpParXYudsgThird(0x0),
		fh2dJetSignedImpParXYbThird(0x0),
		fh2dJetSignedImpParXYbThirdNonBDecay(0x0),
		fh2dJetSignedImpParXYcThird(0x0),
		fh2dJetSignedImpParXYSignificanceThird(0x0),
		fh2dJetSignedImpParXYSignificanceUnidentifiedThird(0x0),
		fh2dJetSignedImpParXYSignificanceudsgThird(0x0),
		fh2dJetSignedImpParXYSignificancebThird(0x0),
		fh2dJetSignedImpParXYSignificancebThirdNonBDecay(0x0),
		fh2dJetSignedImpParXYSignificancecThird(0x0),
		fh2dJetSignedImpParXYZThird(0x0),
		fh2dJetSignedImpParXYZUnidentifiedThird(0x0),
		fh2dJetSignedImpParXYZudsgThird(0x0),
		fh2dJetSignedImpParXYZbThird(0x0),
		fh2dJetSignedImpParXYZbThirdNonBDecay(0x0),
		fh2dJetSignedImpParXYZcThird(0x0),
		fh2dJetSignedImpParXYZSignificanceThird(0x0),
		fh2dJetSignedImpParXYZSignificanceUnidentifiedThird(0x0),
		fh2dJetSignedImpParXYZSignificanceudsgThird(0x0),
		fh2dJetSignedImpParXYZSignificancebThird(0x0),
		fh2dJetSignedImpParXYZSignificancebThirdNonBDecay(0x0),
		fh2dJetSignedImpParXYZSignificancecThird(0x0),
		fh2dJetSignedImpParXYFirst_McCorr(0x0),
		fh2dJetSignedImpParXYUnidentifiedFirst_McCorr(0x0),
		fh2dJetSignedImpParXYudsgFirst_McCorr(0x0),
		fh2dJetSignedImpParXYbFirst_McCorr(0x0),
		fh2dJetSignedImpParXYbFirst_McCorrNonBDecay(0x0),
		fh2dJetSignedImpParXYcFirst_McCorr(0x0),
		fh2dJetSignedImpParXYSignificanceFirst_McCorr(0x0),
		fh2dJetSignedImpParXYSignificanceUnidentifiedFirst_McCorr(0x0),
		fh2dJetSignedImpParXYSignificanceudsgFirst_McCorr(0x0),
		fh2dJetSignedImpParXYSignificancebFirst_McCorr(0x0),
		fh2dJetSignedImpParXYSignificancebFirst_McCorrNonBDecay(0x0),
		fh2dJetSignedImpParXYSignificancecFirst_McCorr(0x0),
		fh2dJetSignedImpParXYZFirst_McCorr(0x0),
		fh2dJetSignedImpParXYZUnidentifiedFirst_McCorr(0x0),
		fh2dJetSignedImpParXYZudsgFirst_McCorr(0x0),
		fh2dJetSignedImpParXYZbFirst_McCorr(0x0),
		fh2dJetSignedImpParXYZbFirst_McCorrNonBDecay(0x0),
		fh2dJetSignedImpParXYZcFirst_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificanceFirst_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificanceudsgFirst_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificancebFirst_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificancebFirst_McCorrNonBDecay(0x0),
		fh2dJetSignedImpParXYZSignificancecFirst_McCorr(0x0),
		fh2dJetSignedImpParXYSecond_McCorr(0x0),
		fh2dJetSignedImpParXYUnidentifiedSecond_McCorr(0x0),
		fh2dJetSignedImpParXYudsgSecond_McCorr(0x0),
		fh2dJetSignedImpParXYbSecond_McCorr(0x0),
		fh2dJetSignedImpParXYbSecond_McCorrNonBDecay(0x0),
		fh2dJetSignedImpParXYcSecond_McCorr(0x0),
		fh2dJetSignedImpParXYSignificanceSecond_McCorr(0x0),
		fh2dJetSignedImpParXYSignificanceUnidentifiedSecond_McCorr(0x0),
		fh2dJetSignedImpParXYSignificanceudsgSecond_McCorr(0x0),
		fh2dJetSignedImpParXYSignificancebSecond_McCorr(0x0),
		fh2dJetSignedImpParXYSignificancebSecond_McCorrNonBDecay(0x0),
		fh2dJetSignedImpParXYSignificancecSecond_McCorr(0x0),
		fh2dJetSignedImpParXYZSecond_McCorr(0x0),
		fh2dJetSignedImpParXYZUnidentifiedSecond_McCorr(0x0),
		fh2dJetSignedImpParXYZudsgSecond_McCorr(0x0),
		fh2dJetSignedImpParXYZbSecond_McCorrNonBDecay(0x0),
		fh2dJetSignedImpParXYZcSecond_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificanceSecond_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificanceudsgSecond_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificancebSecond_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificancebSecond_McCorrNonBDecay(0x0),
		fh2dJetSignedImpParXYZSignificancecSecond_McCorr(0x0),
		fh2dJetSignedImpParXYThird_McCorr(0x0),
		fh2dJetSignedImpParXYUnidentifiedThird_McCorr(0x0),
		fh2dJetSignedImpParXYudsgThird_McCorr(0x0),
		fh2dJetSignedImpParXYbThird_McCorr(0x0),
		fh2dJetSignedImpParXYbThird_McCorrNonBDecay(0x0),
		fh2dJetSignedImpParXYcThird_McCorr(0x0),
		fh2dJetSignedImpParXYSignificanceThird_McCorr(0x0),
		fh2dJetSignedImpParXYSignificanceUnidentifiedThird_McCorr(0x0),
		fh2dJetSignedImpParXYSignificanceudsgThird_McCorr(0x0),
		fh2dJetSignedImpParXYSignificancebThird_McCorr(0x0),
		fh2dJetSignedImpParXYSignificancebThird_McCorrNonBDecay(0x0),
		fh2dJetSignedImpParXYSignificancecThird_McCorr(0x0),
		fh2dJetSignedImpParXYZThird_McCorr(0x0),
		fh2dJetSignedImpParXYZUnidentifiedThird_McCorr(0x0),
		fh2dJetSignedImpParXYZudsgThird_McCorr(0x0),
		fh2dJetSignedImpParXYZbThird_McCorr(0x0),
		fh2dJetSignedImpParXYZbThird_McCorrNonBDecay(0x0),
		fh2dJetSignedImpParXYZcThird_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificanceThird_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificanceUnidentifiedThird_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificanceudsgThird_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificancebThird_McCorr(0x0),
		fh2dJetSignedImpParXYZSignificancebThird_McCorrNonBDecay(0x0),
		fh2dJetSignedImpParXYZSignificancecThird_McCorr(0x0),
		fh1dTracksImpParXY_electron(0x0),
		fh1dTracksImpParXYZ_electron(0x0),
		fh1dTracksImpParXYSignificance_electron(0x0),
		fh1dTracksImpParXYZSignificance_electron(0x0),
		//electrons
		fh2dJetSignedImpParXY_electron(0x0),
		fh2dJetSignedImpParXYUnidentified_electron(0x0),
		fh2dJetSignedImpParXYudsg_electron(0x0),
		fh2dJetSignedImpParXYb_electron(0x0),
		fh2dJetSignedImpParXYc_electron(0x0),
		fh2dJetSignedImpParXYSignificance_electron(0x0),
		fh2dJetSignedImpParXYSignificanceUnidentified_electron(0x0),
		fh2dJetSignedImpParXYSignificanceudsg_electron(0x0),
		fh2dJetSignedImpParXYSignificanceb_electron(0x0),
		fh2dJetSignedImpParXYSignificancec_electron(0x0),
		fh2dJetSignedImpParXYZ_electron(0x0),
		fh2dJetSignedImpParXYZUnidentified_electron(0x0),
		fh2dJetSignedImpParXYZudsg_electron(0x0),
		fh2dJetSignedImpParXYZb_electron(0x0),
		fh2dJetSignedImpParXYZc_electron(0x0),
		fh2dJetSignedImpParXYZSignificance_electron(0x0),
		fh2dJetSignedImpParXYZSignificanceUnidentified_electron(0x0),
		fh2dJetSignedImpParXYZSignificanceudsg_electron(0x0),
		fh2dJetSignedImpParXYZSignificanceb_electron(0x0),
		fh2dJetSignedImpParXYZSignificancec_electron(0x0),
		fh2dJetSignedImpParXYFirst_electron(0x0),
		fh2dJetSignedImpParXYUnidentifiedFirst_electron(0x0),
		fh2dJetSignedImpParXYudsgFirst_electron(0x0),
		fh2dJetSignedImpParXYbFirst_electron(0x0),
		fh2dJetSignedImpParXYcFirst_electron(0x0),
		fh2dJetSignedImpParXYSignificanceFirst_electron(0x0),
		fh2dJetSignedImpParXYSignificanceUnidentifiedFirst_electron(0x0),
		fh2dJetSignedImpParXYSignificanceudsgFirst_electron(0x0),
		fh2dJetSignedImpParXYSignificancebFirst_electron(0x0),
		fh2dJetSignedImpParXYSignificancecFirst_electron(0x0),
		fh2dJetSignedImpParXYZFirst_electron(0x0),
		fh2dJetSignedImpParXYZUnidentifiedFirst_electron(0x0),
		fh2dJetSignedImpParXYZudsgFirst_electron(0x0),
		fh2dJetSignedImpParXYZbFirst_electron(0x0),
		fh2dJetSignedImpParXYZcFirst_electron(0x0),
		fh2dJetSignedImpParXYZSignificanceFirst_electron(0x0),
		fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst_electron(0x0),
		fh2dJetSignedImpParXYZSignificanceudsgFirst_electron(0x0),
		fh2dJetSignedImpParXYZSignificancebFirst_electron(0x0),
		fh2dJetSignedImpParXYZSignificancecFirst_electron(0x0),
		fh2dJetSignedImpParXYSecond_electron(0x0),
		fh2dJetSignedImpParXYUnidentifiedSecond_electron(0x0),
		fh2dJetSignedImpParXYudsgSecond_electron(0x0),
		fh2dJetSignedImpParXYbSecond_electron(0x0),
		fh2dJetSignedImpParXYcSecond_electron(0x0),
		fh2dJetSignedImpParXYSignificanceSecond_electron(0x0),
		fh2dJetSignedImpParXYSignificanceUnidentifiedSecond_electron(0x0),
		fh2dJetSignedImpParXYSignificanceudsgSecond_electron(0x0),
		fh2dJetSignedImpParXYSignificancebSecond_electron(0x0),
		fh2dJetSignedImpParXYSignificancecSecond_electron(0x0),
		fh2dJetSignedImpParXYZSecond_electron(0x0),
		fh2dJetSignedImpParXYZUnidentifiedSecond_electron(0x0),
		fh2dJetSignedImpParXYZudsgSecond_electron(0x0),
		fh2dJetSignedImpParXYZbSecond_electron(0x0),
		fh2dJetSignedImpParXYZcSecond_electron(0x0),
		fh2dJetSignedImpParXYZSignificanceSecond_electron(0x0),
		fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond_electron(0x0),
		fh2dJetSignedImpParXYZSignificanceudsgSecond_electron(0x0),
		fh2dJetSignedImpParXYZSignificancebSecond_electron(0x0),
		fh2dJetSignedImpParXYZSignificancecSecond_electron(0x0),
		fh2dJetSignedImpParXYThird_electron(0x0),
		fh2dJetSignedImpParXYUnidentifiedThird_electron(0x0),
		fh2dJetSignedImpParXYudsgThird_electron(0x0),
		fh2dJetSignedImpParXYbThird_electron(0x0),
		fh2dJetSignedImpParXYcThird_electron(0x0),
		fh2dJetSignedImpParXYSignificanceThird_electron(0x0),
		fh2dJetSignedImpParXYSignificanceUnidentifiedThird_electron(0x0),
		fh2dJetSignedImpParXYSignificanceudsgThird_electron(0x0),
		fh2dJetSignedImpParXYSignificancebThird_electron(0x0),
		fh2dJetSignedImpParXYSignificancecThird_electron(0x0),
		fh2dJetSignedImpParXYZThird_electron(0x0),
		fh2dJetSignedImpParXYZUnidentifiedThird_electron(0x0),
		fh2dJetSignedImpParXYZudsgThird_electron(0x0),
		fh2dJetSignedImpParXYZbThird_electron(0x0),
		fh2dJetSignedImpParXYZcThird_electron(0x0),
		fh2dJetSignedImpParXYZSignificanceThird_electron(0x0),
		fh2dJetSignedImpParXYZSignificanceUnidentifiedThird_electron(0x0),
		fh2dJetSignedImpParXYZSignificanceudsgThird_electron(0x0),
		fh2dJetSignedImpParXYZSignificancebThird_electron(0x0),
		fh2dJetSignedImpParXYZSignificancecThird_electron(0x0),
		fMCArray(0x0),
		fMCEvent(0x0),
		fESDTrackCut(0x0)
		, fUtils(new AliAnalysisUtils())
{

	for(int i =0 ; i<498;++i){
		for(int j =0 ; j<16;++j){
			fBackgroundFactorLinus[j][i]=1; //set default to 1
		}}

	fEnableV0GammaRejection = kFALSE;

}
// ########################################################################################  Main Loop
Bool_t AliAnalysisTaskHFJetIPQA::Run()
{

	if(fESD) fIsEsd = kTRUE;
	AliAODEvent* aev = NULL;
	AliESDEvent* eev = NULL;
	fMCArray     	 = NULL;
	fMCEvent		 = NULL;

	if (!fESD) 	aev = dynamic_cast<AliAODEvent*>(InputEvent());
	else	 	eev = dynamic_cast<AliESDEvent*>(InputEvent());

	if(fIsPythia){
		if(!fESD){
			fMCArray= dynamic_cast<TClonesArray*>(aev->FindListObject(AliAODMCParticle::StdBranchName()));
		}
		else{
			fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
		}
	}
	//Main loop over AOD tracks with filterbit 4  || or ESD tracks with filter
	double weight =1;
	int nTracksInEvent = 0 ;

	if(fESD)nTracksInEvent  = eev ->GetNumberOfTracks() ;
	else nTracksInEvent  = aev ->GetNumberOfTracks() ;

	AliAODTrack * trackAOD =  NULL;
	AliESDtrack * trackESD =  NULL;

	for(long itrack= 0; itrack<nTracksInEvent;++itrack)
	{
		if (!fESD){
			trackAOD = (AliAODTrack*)aev->GetTrack(itrack);
			if (!trackAOD) continue;
			fh1dTracksAccepeted->SetBinContent(1,fh1dTracksAccepeted->GetBinContent(1)+1);

			if(!((trackAOD)->TestFilterBit(1 << 4))) continue; // to avoid double counting in stats hist
		}
		else {
			trackESD = (AliESDtrack*)eev->GetTrack(itrack);
			if (!trackESD) continue;
			fh1dTracksAccepeted->SetBinContent(1,fh1dTracksAccepeted->GetBinContent(1)+1);
		}

		if(!fESD && !IsTrackAccepted(trackAOD)) {
			fh1dTracksAccepeted->SetBinContent(3,fh1dTracksAccepeted->GetBinContent(3)+1);
			continue;
		}
		else if(fESD && !IsTrackAccepted(trackESD)) {
			fh1dTracksAccepeted->SetBinContent(3,fh1dTracksAccepeted->GetBinContent(3)+1);
			continue;
		}

		fh1dTracksAccepeted->SetBinContent(2,fh1dTracksAccepeted->GetBinContent(2)+1);
		//Calculate impact parameters and fill histograms
		double dca[2] = {-99999,-99999};
		double cov[3] = {-99999,-99999,-99999};

		bool hasIPSuccess =kFALSE;
		if (!fESD && CalculateTrackImpactParameter(trackAOD,dca,cov))hasIPSuccess =kTRUE;
		if (fESD && CalculateTrackImpactParameter(trackESD,dca,cov))hasIPSuccess =kTRUE;

		if(hasIPSuccess)
		{
			weight =1;
			fh1dTracksImpParXY->Fill(GetValImpactParameter(kXY,dca,cov));
			fh1dTracksImpParXYZ->Fill(GetValImpactParameter(kXYZ,dca,cov));
			fh1dTracksImpParXYSignificance->Fill(GetValImpactParameter(kXYSig,dca,cov));
			fh1dTracksImpParXYZSignificance->Fill(GetValImpactParameter(kXYZSig,dca,cov));
			if(fIsPythia){
				bool iselectron = false;
				bool isfromBMeson = false;
				if(!fESD)weight = GetMonteCarloCorrectionFactor(trackAOD,iselectron,isfromBMeson);
				else if(fESD)weight = GetMonteCarloCorrectionFactor(trackESD,iselectron,isfromBMeson);
				fh1dTracksImpParXY_McCorr->Fill(GetValImpactParameter(kXY,dca,cov),weight);
				fh1dTracksImpParXYZ_McCorr->Fill(GetValImpactParameter(kXYZ,dca,cov),weight);
				fh1dTracksImpParXYSignificance_McCorr->Fill(GetValImpactParameter(kXYSig,dca,cov),weight);
				fh1dTracksImpParXYZSignificance_McCorr->Fill(GetValImpactParameter(kXYZSig,dca,cov),weight);
				if(iselectron){
					fh1dTracksImpParXY_electron->Fill(GetValImpactParameter(kXY,dca,cov),weight);
					fh1dTracksImpParXYZ_electron->Fill(GetValImpactParameter(kXYZ,dca,cov),weight);
					fh1dTracksImpParXYSignificance_electron->Fill(GetValImpactParameter(kXYSig,dca,cov),weight);
					fh1dTracksImpParXYZSignificance_electron->Fill(GetValImpactParameter(kXYZSig,dca,cov),weight);
				}
				double dcaMC[2] = {-99999,-99999};
				double covMC[3] = {-99999,-99999,-99999};
				Bool_t  hasSuccess = kFALSE;

				if(!fESD && CalculateTrackImpactParameterTruth(trackAOD,dcaMC,covMC))hasSuccess = kTRUE;
				if(fESD && CalculateTrackImpactParameterTruth(trackESD,dcaMC,covMC))hasSuccess = kTRUE;

				if( hasSuccess){
					fh1dTracksImpParXYTruth->Fill(GetValImpactParameter(kXY,dcaMC,covMC));
					fh1dTracksImpParXYZTruth->Fill(GetValImpactParameter(kXYZ,dcaMC,covMC));
					// Fill residual plots
					double residualxy = TMath::Abs(GetValImpactParameter(kXY,dca,cov)) - TMath::Abs(GetValImpactParameter(kXY,dcaMC,covMC));
					residualxy /= TMath::Sqrt(cov[0]);
					fh1dTracksImpParXYResidualTruth->Fill(residualxy);
					double residualxyz = TMath::Abs(GetValImpactParameter(kXYZ,dca,cov)) - TMath::Abs(GetValImpactParameter(kXYZ,dcaMC,covMC));
					residualxyz /= 	GetValImpactParameter(kXYZSigmaOnly,dca,cov);
					fh1dTracksImpParXYZResidualTruth->Fill(residualxyz);
				}
			}
		}
	}

	// Main part jet analysis
	//preparation
	AliJetContainer * jetconrec = 0x0;
	jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(0));
	AliJetContainer * jetcongen = 0x0;
	AliEmcalJet * jetgen  = 0x0;
	AliAODMCParticle* partonAOD = NULL;
	AliMCParticle* partonESD = NULL;

	if(fIsPythia)
	{
		jetcongen = static_cast<AliJetContainer*>(fJetCollArray.At(1));
		if(!MatchJetsGeometricDefault()) cout << "Error running jet matching!" << endl;
		jetcongen->ResetCurrentID();
		// Fill gen. level jet histograms
		while ((jetgen = jetcongen->GetNextJet()))
		{
			if (!jetgen) continue;
			Int_t jetflavour =0;
			Int_t partonpdg=0;
			AliAODMCParticle* parton = NULL;
			if(!fESD){
				partonAOD = fHFJetUtils->IsMCJetParton(fMCArray, jetgen, 0.4);
			}
			else{
				partonESD = fHFJetUtils->IsMCJetParton(MCEvent(), jetgen, 0.4);
			}
			if(!(partonAOD) && !(partonESD)) jetflavour =0;
			else
			{
				if(!fESD)partonpdg = abs(partonAOD->PdgCode());
				else if(fESD)partonpdg = abs(partonESD->PdgCode());
				if(partonpdg==1||partonpdg==2||partonpdg==3||partonpdg==21 )jetflavour=1;
				else if(partonpdg==4)jetflavour=2;
				else if(partonpdg==5)jetflavour=3;
			}

			fh1dJetGenPt->Fill(jetgen->Pt());
			if(jetflavour ==0)
				fh1dJetGenPtUnidentified->Fill(jetgen->Pt());
			else if(jetflavour ==1)
				fh1dJetGenPtudsg->Fill(jetgen->Pt());
			else if(jetflavour ==2)
				fh1dJetGenPtc->Fill(jetgen->Pt());
			else if(jetflavour ==3)
				fh1dJetGenPtb->Fill(jetgen->Pt());
		}
		jetcongen->ResetCurrentID();
		jetconrec->ResetCurrentID();
	}

	// loop rec level jets
	AliEmcalJet * jetrec  = 0x0;
	AliEmcalJet * jetmatched  = 0x0;
	jetconrec->ResetCurrentID();
	double jetpt=0;
	double jetptmc=0;

	while ((jetrec = jetconrec->GetNextJet()))
	{
		//	Printf("%s:%i",__FUNCTION__,__LINE__);
		jetpt= jetrec->Pt();
		//	Printf("%s:%i",__FUNCTION__,__LINE__);

		if(!(jetconrec->GetRhoParameter() == 0x0))
		{
			jetpt = jetpt - jetconrec->GetRhoVal() * jetrec->Area();
		}

		if(fIsPythia){
			if (jetrec->MatchedJet()) {
				double genpt = jetrec->MatchedJet()->Pt();
				if(!(jetcongen->GetRhoParameter() == 0x0)){
					genpt = genpt - jetcongen->GetRhoVal() * jetrec->MatchedJet()->Area();
				}
				fh2dJetGenPtVsJetRecPt->Fill(genpt,jetpt);
			}
		}

		// make inclusive signed imp. parameter constituent histograms
		AliAODTrack* trackAOD = 0x0;
		AliESDtrack* trackESD = 0x0;
		Int_t ntracks = (Int_t)jetrec->GetNumberOfTracks();

		double dca[2] = {-99999,-99999};
		double cov[3] = {-99999,-99999,-99999};
		double sign=0;
		//Printf("%s:%i",__FUNCTION__,__LINE__);

		Int_t jetflavour =0;
		Int_t partonpdg=0;
		if(fIsPythia){
			jetmatched = 0x0;
			jetmatched =jetrec->MatchedJet();

			partonAOD = NULL;
			partonESD = NULL;

			if(jetmatched){
				if(!fESD){
					partonAOD = fHFJetUtils->IsMCJetParton(fMCArray, jetmatched, 0.4);
				}
				else{
					partonESD = fHFJetUtils->IsMCJetParton(MCEvent(), jetmatched, 0.4);
				}

				if((!partonAOD)&&(!partonESD)) jetflavour =0;
				else{
					if(!fESD)partonpdg = abs(partonAOD->PdgCode());
					else if(fESD)partonpdg = abs(partonESD->PdgCode());
					if(partonpdg==1||partonpdg==2||partonpdg==3||partonpdg==21 )jetflavour=1;
					else if(partonpdg==4)jetflavour=2;
					else if(partonpdg==5)jetflavour=3;
				}
			}
		}
		//	Printf("%s:%i",__FUNCTION__,__LINE__);

			fh1dJetRecPt->Fill(jetrec->Pt());
			if(fIsPythia){
				if(jetflavour==0) fh1dJetRecPtUnidentified->Fill(jetpt);
				else if(jetflavour==1)fh1dJetRecPtudsg->Fill(jetpt);
				else if(jetflavour==2)fh1dJetRecPtc->Fill(jetpt);
				else if(jetflavour==3)fh1dJetRecPtb->Fill(jetpt);
			}
			if(!(fJetCutsHF->IsJetSelected(jetrec))) continue;
			//	Printf("%s:%i",__FUNCTION__,__LINE__);

			fh1dJetRecEtaPhiAccepted->Fill(jetrec->Eta(),jetrec->Phi());
			fh1dJetRecPtAccepted->Fill(jetpt);
			if(fIsPythia){
				if(jetflavour==0) fh1dJetRecPtUnidentifiedAccepted->Fill(jetpt);
				else if(jetflavour==1)fh1dJetRecPtudsgAccepted->Fill(jetpt);
				else if(jetflavour==2)fh1dJetRecPtcAccepted->Fill(jetpt);
				else if(jetflavour==3)fh1dJetRecPtbAccepted->Fill(jetpt);
			}
			std::vector<myvaluetuple> sImpParXY,sImpParXYZ,sImpParXYSig,sImpParXYZSig;
			for(Int_t itrack = 0; itrack < ntracks; ++itrack)
			{
				double dcatrackjet =999;
				double lineardecaylenth = 999;

				if (!fESD)trackAOD = (((AliAODTrack*)((jetconrec->GetParticleContainer())->GetParticle(jetrec->TrackAt(itrack)))));
				else 	if (fESD)trackESD = (((AliESDtrack*)((jetconrec->GetParticleContainer())->GetParticle(jetrec->TrackAt(itrack)))));
				//	Printf("%s:%i",__FUNCTION__,__LINE__);

				if(!fESD && !trackAOD) 	continue;
				else if(fESD && !trackESD) 	continue;
				//		Printf("%s:%i",__FUNCTION__,__LINE__);

				if (!fESD  && !IsTrackAccepted(trackAOD)) continue;
				if (fESD  && !IsTrackAccepted(trackESD)) continue;
				//		Printf("%s:%i",__FUNCTION__,__LINE__);

				Bool_t hasSIP =kFALSE;
				//		Printf("%s:%i",__FUNCTION__,__LINE__);

				if(!fESD && CalculateJetSignedTrackImpactParameter(trackAOD,jetrec,dca,cov,sign,dcatrackjet,lineardecaylenth))hasSIP =kTRUE;
				if(fESD && CalculateJetSignedTrackImpactParameter(trackESD,jetrec,dca,cov,sign,dcatrackjet,lineardecaylenth))hasSIP =kTRUE;
				////		Printf("%s:%i",__FUNCTION__,__LINE__);


				if(hasSIP)
				{
					//Select only tracks with dca rphi <1 cm and dca z < 2 cm
					// linear decay length < 10 cm
					// dca track to jet < 0.07 cm
					if(abs(dca[0])>1.) continue;
					if(abs(dca[1])>2.) continue;
					if(lineardecaylenth > 10.) continue;
					if (dcatrackjet > 0.07) continue;

					double cursImParXY =TMath::Abs(GetValImpactParameter(kXY,dca,cov))*sign;
					double cursImParXYZ =TMath::Abs(GetValImpactParameter(kXYZ,dca,cov))*sign;
					double cursImParXYSig =TMath::Abs(GetValImpactParameter(kXYSig,dca,cov))*sign;
					double cursImParXYZSig =TMath::Abs(GetValImpactParameter(kXYZSig,dca,cov))*sign;
					bool iselectron =false;
					bool isfromBMeson = false;

					fh2dJetSignedImpParXY->Fill(jetpt,cursImParXY);
					fh2dJetSignedImpParXYZ->Fill(jetpt,cursImParXYZ);
					fh2dJetSignedImpParXYSignificance->Fill(jetpt,cursImParXYSig);
					fh2dJetSignedImpParXYZSignificance->Fill(jetpt,cursImParXYZSig);

					if(fIsPythia)
					{
						weight =1;
						if (!fESD){
							weight = GetMonteCarloCorrectionFactor(trackAOD,iselectron,isfromBMeson);
							//		Printf("%s:%i",__FUNCTION__,__LINE__);

						}
						else if (fESD){
							weight = GetMonteCarloCorrectionFactor(trackESD,iselectron,isfromBMeson);
							//		Printf("%s:%i",__FUNCTION__,__LINE__);

						}
						fh2dJetSignedImpParXY_McCorr->Fill(jetpt,cursImParXY,weight);
						fh2dJetSignedImpParXYZ_McCorr->Fill(jetpt,cursImParXYZ,weight);
						fh2dJetSignedImpParXYSignificance_McCorr->Fill(jetpt,cursImParXYSig,weight);
						fh2dJetSignedImpParXYZSignificance_McCorr->Fill(jetpt,cursImParXYZSig,weight);
						if(iselectron){
							fh2dJetSignedImpParXY_electron->Fill(jetpt,cursImParXY,weight);
							fh2dJetSignedImpParXYZ_electron->Fill(jetpt,cursImParXYZ,weight);
							fh2dJetSignedImpParXYSignificance_electron->Fill(jetpt,cursImParXYSig,weight);
							fh2dJetSignedImpParXYZSignificance_electron->Fill(jetpt,cursImParXYZSig,weight);
						}

						if(jetflavour ==0){
							fh2dJetSignedImpParXYUnidentified->Fill(jetpt,cursImParXY);
							fh2dJetSignedImpParXYZUnidentified->Fill(jetpt,cursImParXYZ);
							fh2dJetSignedImpParXYSignificanceUnidentified->Fill(jetpt,cursImParXYSig);
							fh2dJetSignedImpParXYZSignificanceUnidentified->Fill(jetpt,cursImParXYZSig);
							fh2dJetSignedImpParXYUnidentified_McCorr->Fill(jetpt,cursImParXY,weight);
							fh2dJetSignedImpParXYZUnidentified_McCorr->Fill(jetpt,cursImParXYZ,weight);
							fh2dJetSignedImpParXYSignificanceUnidentified_McCorr->Fill(jetpt,cursImParXYSig,weight);
							fh2dJetSignedImpParXYZSignificanceUnidentified_McCorr->Fill(jetpt,cursImParXYZSig,weight);
							if(iselectron){
								fh2dJetSignedImpParXYUnidentified_electron->Fill(jetpt,cursImParXY,weight);
								fh2dJetSignedImpParXYZUnidentified_electron->Fill(jetpt,cursImParXYZ,weight);
								fh2dJetSignedImpParXYSignificanceUnidentified_electron->Fill(jetpt,cursImParXYSig,weight);
								fh2dJetSignedImpParXYZSignificanceUnidentified_electron->Fill(jetpt,cursImParXYZSig,weight);
							}
						}
						else if(jetflavour ==1){
							fh2dJetSignedImpParXYudsg->Fill(jetpt,cursImParXY);
							fh2dJetSignedImpParXYZudsg->Fill(jetpt,cursImParXYZ);
							fh2dJetSignedImpParXYSignificanceudsg->Fill(jetpt,cursImParXYSig);
							fh2dJetSignedImpParXYZSignificanceudsg->Fill(jetpt,cursImParXYZSig);
							fh2dJetSignedImpParXYudsg_McCorr->Fill(jetpt,cursImParXY,weight);
							fh2dJetSignedImpParXYZudsg_McCorr->Fill(jetpt,cursImParXYZ,weight);
							fh2dJetSignedImpParXYSignificanceudsg_McCorr->Fill(jetpt,cursImParXYSig,weight);
							fh2dJetSignedImpParXYZSignificanceudsg_McCorr->Fill(jetpt,cursImParXYZSig,weight);
							if(iselectron){
								fh2dJetSignedImpParXYudsg_electron->Fill(jetpt,cursImParXY,weight);
								fh2dJetSignedImpParXYZudsg_electron->Fill(jetpt,cursImParXYZ,weight);
								fh2dJetSignedImpParXYSignificanceudsg_electron->Fill(jetpt,cursImParXYSig,weight);
								fh2dJetSignedImpParXYZSignificanceudsg_electron->Fill(jetpt,cursImParXYZSig,weight);
							}
						}
						else if(jetflavour ==2){
							fh2dJetSignedImpParXYc->Fill(jetpt,cursImParXY);
							fh2dJetSignedImpParXYZc->Fill(jetpt,cursImParXYZ);
							fh2dJetSignedImpParXYSignificancec->Fill(jetpt,cursImParXYSig);
							fh2dJetSignedImpParXYZSignificancec->Fill(jetpt,cursImParXYZSig);
							fh2dJetSignedImpParXYc_McCorr->Fill(jetpt,cursImParXY,weight);
							fh2dJetSignedImpParXYZc_McCorr->Fill(jetpt,cursImParXYZ,weight);
							fh2dJetSignedImpParXYSignificancec_McCorr->Fill(jetpt,cursImParXYSig,weight);
							fh2dJetSignedImpParXYZSignificancec_McCorr->Fill(jetpt,cursImParXYZSig,weight);
							if(iselectron){
								fh2dJetSignedImpParXYc_electron->Fill(jetpt,cursImParXY,weight);
								fh2dJetSignedImpParXYZc_electron->Fill(jetpt,cursImParXYZ,weight);
								fh2dJetSignedImpParXYSignificancec_electron->Fill(jetpt,cursImParXYSig,weight);
								fh2dJetSignedImpParXYZSignificancec_electron->Fill(jetpt,cursImParXYZSig,weight);
							}
						}
						else if(jetflavour ==3){

							fh2dJetSignedImpParXYb->Fill(jetpt,cursImParXY);
							fh2dJetSignedImpParXYZb->Fill(jetpt,cursImParXYZ);
							fh2dJetSignedImpParXYSignificanceb->Fill(jetpt,cursImParXYSig);
							fh2dJetSignedImpParXYZSignificanceb->Fill(jetpt,cursImParXYZSig);

							if(isfromBMeson){
								fh2dJetSignedImpParXYb_McCorr->Fill(jetpt,cursImParXY,weight);
								fh2dJetSignedImpParXYZb_McCorr->Fill(jetpt,cursImParXYZ,weight);
								fh2dJetSignedImpParXYSignificanceb_McCorr->Fill(jetpt,cursImParXYSig,weight);
								fh2dJetSignedImpParXYZSignificanceb_McCorr->Fill(jetpt,cursImParXYZSig,weight);
							}
							else{
								fh2dJetSignedImpParXYb_McCorrNonBDecay->Fill(jetpt,cursImParXY,weight);
								fh2dJetSignedImpParXYZb_McCorr_McCorrNonBDecay->Fill(jetpt,cursImParXYZ,weight);
								fh2dJetSignedImpParXYSignificanceb_McCorrNonBDecay->Fill(jetpt,cursImParXYSig,weight);
								fh2dJetSignedImpParXYZSignificanceb_McCorrNonBDecay->Fill(jetpt,cursImParXYZSig,weight);
							}

							if(iselectron){
								fh2dJetSignedImpParXYb_electron->Fill(jetpt,cursImParXY,weight);
								fh2dJetSignedImpParXYZb_electron->Fill(jetpt,cursImParXYZ,weight);
								fh2dJetSignedImpParXYSignificanceb_electron->Fill(jetpt,cursImParXYSig,weight);
								fh2dJetSignedImpParXYZSignificanceb_electron->Fill(jetpt,cursImParXYZSig,weight);
							}
						}
					}
					myvaluetuple a(cursImParXY, weight,iselectron,isfromBMeson);
					sImpParXY.push_back(a);
					myvaluetuple b(cursImParXYZ, weight,iselectron,isfromBMeson);
					sImpParXYZ.push_back(b);
					myvaluetuple c(cursImParXYSig, weight,iselectron,isfromBMeson);
					sImpParXYSig.push_back(c);
					myvaluetuple d(cursImParXYZSig, weight,iselectron,isfromBMeson);
					sImpParXYZSig.push_back(d);
				}
			}
			// end of track loop
			std::sort(sImpParXY.begin(),sImpParXY.end(), AliAnalysisTaskHFJetIPQA::mysort);
			std::sort(sImpParXYZ.begin(),sImpParXYZ.end(), AliAnalysisTaskHFJetIPQA::mysort);
			std::sort(sImpParXYSig.begin(),sImpParXYSig.end(), AliAnalysisTaskHFJetIPQA::mysort);
			std::sort(sImpParXYZSig.begin(),sImpParXYZSig.end(), AliAnalysisTaskHFJetIPQA::mysort);

			std::reverse(sImpParXY.begin(),sImpParXY.end());
			std::reverse(sImpParXYZ.begin(),sImpParXYZ.end());
			std::reverse(sImpParXYSig.begin(),sImpParXYSig.end());
			std::reverse(sImpParXYZSig.begin(),sImpParXYZSig.end());


			//Ordered n=1,2,3 sip
			if (sImpParXY.size()>0){
				fh2dJetSignedImpParXYFirst->Fill(jetpt,sImpParXY.at(0).first);
				fh2dJetSignedImpParXYZFirst->Fill(jetpt,sImpParXYZ.at(0).first);
				fh2dJetSignedImpParXYSignificanceFirst->Fill(jetpt,sImpParXYSig.at(0).first);
				fh2dJetSignedImpParXYZSignificanceFirst->Fill(jetpt,sImpParXYZSig.at(0).first);


				if(fIsPythia){
					fh2dJetSignedImpParXYFirst_McCorr->Fill(jetpt,sImpParXY.at(0).first,sImpParXY.at(0).second);
					fh2dJetSignedImpParXYZFirst_McCorr->Fill(jetpt,sImpParXYZ.at(0).first,sImpParXYZ.at(0).second);
					fh2dJetSignedImpParXYSignificanceFirst_McCorr->Fill(jetpt,sImpParXYSig.at(0).first,sImpParXYSig.at(0).second);
					fh2dJetSignedImpParXYZSignificanceFirst_McCorr->Fill(jetpt,sImpParXYZSig.at(0).first,sImpParXYZSig.at(0).second);

					if(jetflavour ==0){
						fh2dJetSignedImpParXYUnidentifiedFirst->Fill(jetpt,sImpParXY.at(0).first);
						fh2dJetSignedImpParXYZUnidentifiedFirst->Fill(jetpt,sImpParXYZ.at(0).first);
						fh2dJetSignedImpParXYSignificanceUnidentifiedFirst->Fill(jetpt,sImpParXYSig.at(0).first);
						fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst->Fill(jetpt,sImpParXYZSig.at(0).first);

						fh2dJetSignedImpParXYUnidentifiedFirst_McCorr->Fill(jetpt,sImpParXY.at(0).first,sImpParXY.at(0).second);
						fh2dJetSignedImpParXYZUnidentifiedFirst_McCorr->Fill(jetpt,sImpParXYZ.at(0).first,sImpParXYZ.at(0).second);
						fh2dJetSignedImpParXYSignificanceUnidentifiedFirst_McCorr->Fill(jetpt,sImpParXYSig.at(0).first,sImpParXYSig.at(0).second);
						fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst_McCorr->Fill(jetpt,sImpParXYZSig.at(0).first,sImpParXYZSig.at(0).second);

						//electron
						if(sImpParXY.at(0).is_electron)		fh2dJetSignedImpParXYUnidentifiedFirst_electron->Fill(jetpt,sImpParXY.at(0).first,sImpParXY.at(0).second);
						if(sImpParXYZ.at(0).is_electron)	fh2dJetSignedImpParXYZUnidentifiedFirst_electron->Fill(jetpt,sImpParXYZ.at(0).first,sImpParXYZ.at(0).second);
						if(sImpParXYSig.at(0).is_electron)	fh2dJetSignedImpParXYSignificanceUnidentifiedFirst_electron->Fill(jetpt,sImpParXYSig.at(0).first,sImpParXYSig.at(0).second);
						if(sImpParXYZSig.at(0).is_electron)	fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst_electron->Fill(jetpt,sImpParXYZSig.at(0).first,sImpParXYZSig.at(0).second);
					}
					else if(jetflavour ==1){
						fh2dJetSignedImpParXYudsgFirst->Fill(jetpt,sImpParXY.at(0).first);
						fh2dJetSignedImpParXYZudsgFirst->Fill(jetpt,sImpParXYZ.at(0).first);
						fh2dJetSignedImpParXYSignificanceudsgFirst->Fill(jetpt,sImpParXYSig.at(0).first);
						fh2dJetSignedImpParXYZSignificanceudsgFirst->Fill(jetpt,sImpParXYZSig.at(0).first);

						fh2dJetSignedImpParXYudsgFirst_McCorr->Fill(jetpt,sImpParXY.at(0).first,sImpParXY.at(0).second);
						fh2dJetSignedImpParXYZudsgFirst_McCorr->Fill(jetpt,sImpParXYZ.at(0).first, sImpParXYZ.at(0).second );
						fh2dJetSignedImpParXYSignificanceudsgFirst_McCorr->Fill(jetpt,sImpParXYSig.at(0).first,sImpParXYSig.at(0).second);
						fh2dJetSignedImpParXYZSignificanceudsgFirst_McCorr->Fill(jetpt,sImpParXYZSig.at(0).first,sImpParXYZSig.at(0).second);

						if(sImpParXY.at(0).is_electron)		fh2dJetSignedImpParXYudsgFirst_electron->Fill(jetpt,sImpParXY.at(0).first,sImpParXY.at(0).second);
						if(sImpParXYZ.at(0).is_electron)	fh2dJetSignedImpParXYZudsgFirst_electron->Fill(jetpt,sImpParXYZ.at(0).first,sImpParXYZ.at(0).second);
						if(sImpParXYSig.at(0).is_electron)	fh2dJetSignedImpParXYSignificanceudsgFirst_electron->Fill(jetpt,sImpParXYSig.at(0).first,sImpParXYSig.at(0).second);
						if(sImpParXYZSig.at(0).is_electron)	fh2dJetSignedImpParXYZSignificanceudsgFirst_electron->Fill(jetpt,sImpParXYZSig.at(0).first,sImpParXYZSig.at(0).second);
					}
					else if(jetflavour ==2){
						fh2dJetSignedImpParXYcFirst->Fill(jetpt,sImpParXY.at(0).first);
						fh2dJetSignedImpParXYZcFirst->Fill(jetpt,sImpParXYZ.at(0).first);
						fh2dJetSignedImpParXYSignificancecFirst->Fill(jetpt,sImpParXYSig.at(0).first);
						fh2dJetSignedImpParXYZSignificancecFirst->Fill(jetpt,sImpParXYZSig.at(0).first);

						fh2dJetSignedImpParXYcFirst_McCorr->Fill(jetpt,sImpParXY.at(0).first,sImpParXY.at(0).second);
						fh2dJetSignedImpParXYZcFirst_McCorr->Fill(jetpt,sImpParXYZ.at(0).first,sImpParXYZ.at(0).second);
						fh2dJetSignedImpParXYSignificancecFirst_McCorr->Fill(jetpt,sImpParXYSig.at(0).first,sImpParXYSig.at(0).second);
						fh2dJetSignedImpParXYZSignificancecFirst_McCorr->Fill(jetpt,sImpParXYZSig.at(0).first,sImpParXYZSig.at(0).second);


						if(sImpParXY.at(0).is_electron)	fh2dJetSignedImpParXYcFirst_electron->Fill(jetpt,sImpParXY.at(0).first,sImpParXY.at(0).second);
						if(sImpParXYZ.at(0).is_electron)fh2dJetSignedImpParXYZcFirst_electron->Fill(jetpt,sImpParXYZ.at(0).first,sImpParXYZ.at(0).second);
						if(sImpParXYSig.at(0).is_electron)	fh2dJetSignedImpParXYSignificancecFirst_electron->Fill(jetpt,sImpParXYSig.at(0).first,sImpParXYSig.at(0).second);
						if(sImpParXYZSig.at(0).is_electron)	fh2dJetSignedImpParXYZSignificancecFirst_electron->Fill(jetpt,sImpParXYZSig.at(0).first,sImpParXYZSig.at(0).second);
					}
					else if(jetflavour ==3){
						fh2dJetSignedImpParXYbFirst->Fill(jetpt,sImpParXY.at(0).first);
						fh2dJetSignedImpParXYZbFirst->Fill(jetpt,sImpParXYZ.at(0).first);
						fh2dJetSignedImpParXYSignificancebFirst->Fill(jetpt,sImpParXYSig.at(0).first);
						fh2dJetSignedImpParXYZSignificancebFirst->Fill(jetpt,sImpParXYZSig.at(0).first);

						if(sImpParXY.at(0).is_fromB) fh2dJetSignedImpParXYbFirst_McCorr->Fill(jetpt,sImpParXY.at(0).first,sImpParXY.at(0).second);
						if(sImpParXYZ.at(0).is_fromB) fh2dJetSignedImpParXYZbFirst_McCorr->Fill(jetpt,sImpParXYZ.at(0).first,sImpParXYZ.at(0).second);
						if(sImpParXYSig.at(0).is_fromB) fh2dJetSignedImpParXYSignificancebFirst_McCorr->Fill(jetpt,sImpParXYSig.at(0).first,sImpParXYSig.at(0).second);
						if(sImpParXYZSig.at(0).is_fromB) fh2dJetSignedImpParXYZSignificancebFirst_McCorr->Fill(jetpt,sImpParXYZSig.at(0).first,sImpParXYZSig.at(0).second);

						if(!sImpParXY.at(0).is_fromB)fh2dJetSignedImpParXYbFirst_McCorrNonBDecay->Fill(jetpt,sImpParXY.at(0).first,sImpParXY.at(0).second);
						if(!sImpParXYZ.at(0).is_fromB)fh2dJetSignedImpParXYZbFirst_McCorrNonBDecay->Fill(jetpt,sImpParXYZ.at(0).first,sImpParXYZ.at(0).second);
						if(!sImpParXYSig.at(0).is_fromB)fh2dJetSignedImpParXYSignificancebFirst_McCorrNonBDecay->Fill(jetpt,sImpParXYSig.at(0).first,sImpParXYSig.at(0).second);
						if(!sImpParXYZSig.at(0).is_fromB)fh2dJetSignedImpParXYZSignificancebFirst_McCorrNonBDecay->Fill(jetpt,sImpParXYZSig.at(0).first,sImpParXYZSig.at(0).second);

						if(sImpParXY.at(0).is_electron)	fh2dJetSignedImpParXYbFirst_electron->Fill(jetpt,sImpParXY.at(0).first,sImpParXY.at(0).second);
						if(sImpParXYZ.at(0).is_electron)fh2dJetSignedImpParXYZbFirst_electron->Fill(jetpt,sImpParXYZ.at(0).first,sImpParXYZ.at(0).second);
						if(sImpParXYSig.at(0).is_electron) fh2dJetSignedImpParXYSignificancebFirst_electron->Fill(jetpt,sImpParXYSig.at(0).first,sImpParXYSig.at(0).second);
						if(sImpParXYZSig.at(0).is_electron) fh2dJetSignedImpParXYZSignificancebFirst_electron->Fill(jetpt,sImpParXYZSig.at(0).first,sImpParXYZSig.at(0).second);
					}
				}
			}

			//Second largest
			if (sImpParXY.size()>1)
			{

				fh2dJetSignedImpParXYSecond->Fill(jetpt,sImpParXY.at(1).first);
				fh2dJetSignedImpParXYZSecond->Fill(jetpt,sImpParXYZ.at(1).first);
				fh2dJetSignedImpParXYSignificanceSecond->Fill(jetpt,sImpParXYSig.at(1).first);
				fh2dJetSignedImpParXYZSignificanceSecond->Fill(jetpt,sImpParXYZSig.at(1).first);


				if(fIsPythia){
					fh2dJetSignedImpParXYSecond_McCorr->Fill(jetpt,sImpParXY.at(1).first,sImpParXY.at(1).second);
					fh2dJetSignedImpParXYZSecond_McCorr->Fill(jetpt,sImpParXYZ.at(1).first,sImpParXYZ.at(1).second);
					fh2dJetSignedImpParXYSignificanceSecond_McCorr->Fill(jetpt,sImpParXYSig.at(1).first,sImpParXYSig.at(1).second);
					fh2dJetSignedImpParXYZSignificanceSecond_McCorr->Fill(jetpt,sImpParXYZSig.at(1).first,sImpParXYZSig.at(1).second);

					if(jetflavour ==0){
						fh2dJetSignedImpParXYUnidentifiedSecond->Fill(jetpt,sImpParXY.at(1).first);
						fh2dJetSignedImpParXYZUnidentifiedSecond->Fill(jetpt,sImpParXYZ.at(1).first);
						fh2dJetSignedImpParXYSignificanceUnidentifiedSecond->Fill(jetpt,sImpParXYSig.at(1).first);
						fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond->Fill(jetpt,sImpParXYZSig.at(1).first);

						fh2dJetSignedImpParXYUnidentifiedSecond_McCorr->Fill(jetpt,sImpParXY.at(1).first,sImpParXY.at(1).second);
						fh2dJetSignedImpParXYZUnidentifiedSecond_McCorr->Fill(jetpt,sImpParXYZ.at(1).first,sImpParXYZ.at(1).second);
						fh2dJetSignedImpParXYSignificanceUnidentifiedSecond_McCorr->Fill(jetpt,sImpParXYSig.at(1).first,sImpParXYSig.at(1).second);
						fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond_McCorr->Fill(jetpt,sImpParXYZSig.at(1).first,sImpParXYZSig.at(1).second);

						if(sImpParXY.at(1).is_electron)		fh2dJetSignedImpParXYUnidentifiedSecond_electron->Fill(jetpt,sImpParXY.at(1).first,sImpParXY.at(1).second);
						if(sImpParXYZ.at(1).is_electron)	fh2dJetSignedImpParXYZUnidentifiedSecond_electron->Fill(jetpt,sImpParXYZ.at(1).first,sImpParXYZ.at(1).second);
						if(sImpParXYSig.at(1).is_electron)	fh2dJetSignedImpParXYSignificanceUnidentifiedSecond_electron->Fill(jetpt,sImpParXYSig.at(1).first,sImpParXYSig.at(1).second);
						if(sImpParXYZSig.at(1).is_electron)	fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond_electron->Fill(jetpt,sImpParXYZSig.at(1).first,sImpParXYZSig.at(1).second);
					}
					else if(jetflavour ==1){
						fh2dJetSignedImpParXYudsgSecond->Fill(jetpt,sImpParXY.at(1).first);
						fh2dJetSignedImpParXYZudsgSecond->Fill(jetpt,sImpParXYZ.at(1).first);
						fh2dJetSignedImpParXYSignificanceudsgSecond->Fill(jetpt,sImpParXYSig.at(1).first);
						fh2dJetSignedImpParXYZSignificanceudsgSecond->Fill(jetpt,sImpParXYZSig.at(1).first);

						fh2dJetSignedImpParXYudsgSecond_McCorr->Fill(jetpt,sImpParXY.at(1).first,sImpParXY.at(1).second);
						fh2dJetSignedImpParXYZudsgSecond_McCorr->Fill(jetpt,sImpParXYZ.at(1).first,sImpParXYZ.at(1).second);
						fh2dJetSignedImpParXYSignificanceudsgSecond_McCorr->Fill(jetpt,sImpParXYSig.at(1).first,sImpParXYSig.at(1).second);
						fh2dJetSignedImpParXYZSignificanceudsgSecond_McCorr->Fill(jetpt,sImpParXYZSig.at(1).first,sImpParXYZSig.at(1).second);

						if(sImpParXY.at(1).is_electron)		fh2dJetSignedImpParXYudsgSecond_electron->Fill(jetpt,sImpParXY.at(1).first,sImpParXY.at(1).second);
						if(sImpParXYZ.at(1).is_electron)	fh2dJetSignedImpParXYZudsgSecond_electron->Fill(jetpt,sImpParXYZ.at(1).first,sImpParXYZ.at(1).second);
						if(sImpParXYSig.at(1).is_electron)	fh2dJetSignedImpParXYSignificanceudsgSecond_electron->Fill(jetpt,sImpParXYSig.at(1).first,sImpParXYSig.at(1).second);
						if(sImpParXYZSig.at(1).is_electron)	fh2dJetSignedImpParXYZSignificanceudsgSecond_electron->Fill(jetpt,sImpParXYZSig.at(1).first,sImpParXYZSig.at(1).second);
					}
					else if(jetflavour ==2){
						fh2dJetSignedImpParXYcSecond->Fill(jetpt,sImpParXY.at(1).first);
						fh2dJetSignedImpParXYZcSecond->Fill(jetpt,sImpParXYZ.at(1).first);
						fh2dJetSignedImpParXYSignificancecSecond->Fill(jetpt,sImpParXYSig.at(1).first);
						fh2dJetSignedImpParXYZSignificancecSecond->Fill(jetpt,sImpParXYZSig.at(1).first);

						fh2dJetSignedImpParXYcSecond_McCorr->Fill(jetpt,sImpParXY.at(1).first,sImpParXY.at(1).second);
						fh2dJetSignedImpParXYZcSecond_McCorr->Fill(jetpt,sImpParXYZ.at(1).first,sImpParXYZ.at(1).second);
						fh2dJetSignedImpParXYSignificancecSecond_McCorr->Fill(jetpt,sImpParXYSig.at(1).first,sImpParXYSig.at(1).second);
						fh2dJetSignedImpParXYZSignificancecSecond_McCorr->Fill(jetpt,sImpParXYZSig.at(1).first,sImpParXYZSig.at(1).second);



						if(sImpParXY.at(1).is_electron)			fh2dJetSignedImpParXYcSecond_electron->Fill(jetpt,sImpParXY.at(1).first,sImpParXY.at(1).second);
						if(sImpParXYZ.at(1).is_electron)		fh2dJetSignedImpParXYZcSecond_electron->Fill(jetpt,sImpParXYZ.at(1).first,sImpParXYZ.at(1).second);
						if(sImpParXYSig.at(1).is_electron)		fh2dJetSignedImpParXYSignificancecSecond_electron->Fill(jetpt,sImpParXYSig.at(1).first,sImpParXYSig.at(1).second);
						if(sImpParXYZSig.at(1).is_electron)		fh2dJetSignedImpParXYZSignificancecSecond_electron->Fill(jetpt,sImpParXYZSig.at(1).first,sImpParXYZSig.at(1).second);
					}
					else if(jetflavour ==3)
					{
						fh2dJetSignedImpParXYbSecond->Fill(jetpt,sImpParXY.at(1).first);
						fh2dJetSignedImpParXYZbSecond->Fill(jetpt,sImpParXYZ.at(1).first);
						fh2dJetSignedImpParXYSignificancebSecond->Fill(jetpt,sImpParXYSig.at(1).first);
						fh2dJetSignedImpParXYZSignificancebSecond->Fill(jetpt,sImpParXYZSig.at(1).first);

						if(sImpParXY.at(1).is_fromB)fh2dJetSignedImpParXYbSecond_McCorr->Fill(jetpt,sImpParXY.at(1).first,sImpParXY.at(1).second);
						if(sImpParXYZ.at(1).is_fromB)fh2dJetSignedImpParXYZbSecond_McCorr->Fill(jetpt,sImpParXYZ.at(1).first,sImpParXYZ.at(1).second);
						if(sImpParXYSig.at(1).is_fromB)fh2dJetSignedImpParXYSignificancebSecond_McCorr->Fill(jetpt,sImpParXYSig.at(1).first,sImpParXYSig.at(1).second);
						if(sImpParXYZSig.at(1).is_fromB)fh2dJetSignedImpParXYZSignificancebSecond_McCorr->Fill(jetpt,sImpParXYZSig.at(1).first,sImpParXYZSig.at(1).second);

						if(!sImpParXY.at(1).is_fromB)	fh2dJetSignedImpParXYbSecond_McCorrNonBDecay->Fill(jetpt,sImpParXY.at(1).first,sImpParXY.at(1).second);
						if(!sImpParXYZ.at(1).is_fromB)	fh2dJetSignedImpParXYZbSecond_McCorrNonBDecay->Fill(jetpt,sImpParXYZ.at(1).first,sImpParXYZ.at(1).second);
						if(!sImpParXYSig.at(1).is_fromB)	fh2dJetSignedImpParXYSignificancebSecond_McCorrNonBDecay->Fill(jetpt,sImpParXYSig.at(1).first,sImpParXYSig.at(1).second);
						if(!sImpParXYZSig.at(1).is_fromB)	fh2dJetSignedImpParXYZSignificancebSecond_McCorrNonBDecay->Fill(jetpt,sImpParXYZSig.at(1).first,sImpParXYZSig.at(1).second);

						if(sImpParXY.at(1).is_electron)		fh2dJetSignedImpParXYbSecond_electron->Fill(jetpt,sImpParXY.at(1).first,sImpParXY.at(1).second);
						if(sImpParXYZ.at(1).is_electron)	fh2dJetSignedImpParXYZbSecond_electron->Fill(jetpt,sImpParXYZ.at(1).first,sImpParXYZ.at(1).second);
						if(sImpParXYSig.at(1).is_electron)	fh2dJetSignedImpParXYSignificancebSecond_electron->Fill(jetpt,sImpParXYSig.at(1).first,sImpParXYSig.at(1).second);
						if(sImpParXYZSig.at(1).is_electron)	fh2dJetSignedImpParXYZSignificancebSecond_electron->Fill(jetpt,sImpParXYZSig.at(1).first,sImpParXYZSig.at(1).second);
					}
				}
			}
			//Third largest

			if (sImpParXY.size()>2)
			{

				//Printf("%f %f %f",sImpParXY.at(0).second,sImpParXY.at(1).second,sImpParXY.at(2).second);
				fh2dJetSignedImpParXYThird->Fill(jetpt,sImpParXY.at(2).first);
				fh2dJetSignedImpParXYZThird->Fill(jetpt,sImpParXYZ.at(2).first);
				fh2dJetSignedImpParXYSignificanceThird->Fill(jetpt,sImpParXYSig.at(2).first);
				fh2dJetSignedImpParXYZSignificanceThird->Fill(jetpt,sImpParXYZSig.at(2).first);


				if(fIsPythia){
					fh2dJetSignedImpParXYThird_McCorr->Fill(jetpt,sImpParXY.at(2).first,sImpParXY.at(2).second);
					fh2dJetSignedImpParXYZThird_McCorr->Fill(jetpt,sImpParXYZ.at(2).first,sImpParXYZ.at(2).second);
					fh2dJetSignedImpParXYSignificanceThird_McCorr->Fill(jetpt,sImpParXYSig.at(2).first,sImpParXYSig.at(2).second);
					fh2dJetSignedImpParXYZSignificanceThird_McCorr->Fill(jetpt,sImpParXYZSig.at(2).first,sImpParXYZSig.at(2).second);





					if(jetflavour ==0){
						fh2dJetSignedImpParXYUnidentifiedThird->Fill(jetpt,sImpParXY.at(2).first);
						fh2dJetSignedImpParXYZUnidentifiedThird->Fill(jetpt,sImpParXYZ.at(2).first);
						fh2dJetSignedImpParXYSignificanceUnidentifiedThird->Fill(jetpt,sImpParXYSig.at(2).first);
						fh2dJetSignedImpParXYZSignificanceUnidentifiedThird->Fill(jetpt,sImpParXYZSig.at(2).first);

						fh2dJetSignedImpParXYUnidentifiedThird_McCorr->Fill(jetpt,sImpParXY.at(2).first,sImpParXY.at(2).second);
						fh2dJetSignedImpParXYZUnidentifiedThird_McCorr->Fill(jetpt,sImpParXYZ.at(2).first,sImpParXYZ.at(2).second);
						fh2dJetSignedImpParXYSignificanceUnidentifiedThird_McCorr->Fill(jetpt,sImpParXYSig.at(2).first,sImpParXYSig.at(2).second);
						fh2dJetSignedImpParXYZSignificanceUnidentifiedThird_McCorr->Fill(jetpt,sImpParXYZSig.at(2).first,sImpParXYZSig.at(2).second);


						if(sImpParXY.at(2).is_electron)			fh2dJetSignedImpParXYUnidentifiedThird_electron->Fill(jetpt,sImpParXY.at(1).first,sImpParXY.at(2).second);
						if(sImpParXYZ.at(2).is_electron)		fh2dJetSignedImpParXYZUnidentifiedThird_electron->Fill(jetpt,sImpParXYZ.at(2).first,sImpParXYZ.at(2).second);
						if(sImpParXYSig.at(2).is_electron)		fh2dJetSignedImpParXYSignificanceUnidentifiedThird_electron->Fill(jetpt,sImpParXYSig.at(2).first,sImpParXYSig.at(2).second);
						if(sImpParXYZSig.at(2).is_electron)		fh2dJetSignedImpParXYZSignificanceUnidentifiedThird_electron->Fill(jetpt,sImpParXYZSig.at(2).first,sImpParXYZSig.at(2).second);
					}
					else if(jetflavour ==1){
						fh2dJetSignedImpParXYudsgThird->Fill(jetpt,sImpParXY.at(2).first);
						fh2dJetSignedImpParXYZudsgThird->Fill(jetpt,sImpParXYZ.at(2).first);
						fh2dJetSignedImpParXYSignificanceudsgThird->Fill(jetpt,sImpParXYSig.at(2).first);
						fh2dJetSignedImpParXYZSignificanceudsgThird->Fill(jetpt,sImpParXYZSig.at(2).first);

						fh2dJetSignedImpParXYudsgThird_McCorr->Fill(jetpt,sImpParXY.at(2).first,sImpParXY.at(2).second);
						fh2dJetSignedImpParXYZudsgThird_McCorr->Fill(jetpt,sImpParXYZ.at(2).first,sImpParXYZ.at(2).second);
						fh2dJetSignedImpParXYSignificanceudsgThird_McCorr->Fill(jetpt,sImpParXYSig.at(2).first,sImpParXYSig.at(2).second);
						fh2dJetSignedImpParXYZSignificanceudsgThird_McCorr->Fill(jetpt,sImpParXYZSig.at(2).first,sImpParXYZSig.at(2).second);

						if(sImpParXY.at(2).is_electron)		fh2dJetSignedImpParXYudsgThird_electron->Fill(jetpt,sImpParXY.at(1).first,sImpParXY.at(2).second);
						if(sImpParXYZ.at(2).is_electron)	fh2dJetSignedImpParXYZudsgThird_electron->Fill(jetpt,sImpParXYZ.at(2).first,sImpParXYZ.at(2).second);
						if(sImpParXYSig.at(2).is_electron)	fh2dJetSignedImpParXYSignificanceudsgThird_electron->Fill(jetpt,sImpParXYSig.at(2).first,sImpParXYSig.at(2).second);
						if(sImpParXYZSig.at(2).is_electron)	fh2dJetSignedImpParXYZSignificanceudsgThird_electron->Fill(jetpt,sImpParXYZSig.at(2).first,sImpParXYZSig.at(2).second);
					}
					else if(jetflavour ==2){
						fh2dJetSignedImpParXYcThird->Fill(jetpt,sImpParXY.at(2).first);
						fh2dJetSignedImpParXYZcThird->Fill(jetpt,sImpParXYZ.at(2).first);
						fh2dJetSignedImpParXYSignificancecThird->Fill(jetpt,sImpParXYSig.at(2).first);
						fh2dJetSignedImpParXYZSignificancecThird->Fill(jetpt,sImpParXYZSig.at(2).first);

						fh2dJetSignedImpParXYcThird_McCorr->Fill(jetpt,sImpParXY.at(2).first,sImpParXY.at(2).second);
						fh2dJetSignedImpParXYZcThird_McCorr->Fill(jetpt,sImpParXYZ.at(2).first,sImpParXYZ.at(2).second);
						fh2dJetSignedImpParXYSignificancecThird_McCorr->Fill(jetpt,sImpParXYSig.at(2).first,sImpParXYSig.at(2).second);
						fh2dJetSignedImpParXYZSignificancecThird_McCorr->Fill(jetpt,sImpParXYZSig.at(2).first,sImpParXYZSig.at(2).second);


						if(sImpParXY.at(2).is_electron)			fh2dJetSignedImpParXYcThird_electron->Fill(jetpt,sImpParXY.at(1).first,sImpParXY.at(2).second);
						if(sImpParXYZ.at(2).is_electron)		fh2dJetSignedImpParXYZcThird_electron->Fill(jetpt,sImpParXYZ.at(2).first,sImpParXYZ.at(2).second);
						if(sImpParXYSig.at(2).is_electron)		fh2dJetSignedImpParXYSignificancecThird_electron->Fill(jetpt,sImpParXYSig.at(2).first,sImpParXYSig.at(2).second);
						if(sImpParXYZSig.at(2).is_electron)  	fh2dJetSignedImpParXYZSignificancecThird_electron->Fill(jetpt,sImpParXYZSig.at(2).first,sImpParXYZSig.at(2).second);
					}
					else if(jetflavour ==3){
						fh2dJetSignedImpParXYbThird->Fill(jetpt,sImpParXY.at(2).first);
						fh2dJetSignedImpParXYZbThird->Fill(jetpt,sImpParXYZ.at(2).first);
						fh2dJetSignedImpParXYSignificancebThird->Fill(jetpt,sImpParXYSig.at(2).first);
						fh2dJetSignedImpParXYZSignificancebFirst->Fill(jetpt,sImpParXYZSig.at(2).first);

						if(sImpParXY.at(2).is_fromB) fh2dJetSignedImpParXYbThird_McCorr->Fill(jetpt,sImpParXY.at(2).first,sImpParXY.at(2).second);
						if(sImpParXYZ.at(2).is_fromB)fh2dJetSignedImpParXYZbThird_McCorr->Fill(jetpt,sImpParXYZ.at(2).first,sImpParXYZ.at(2).second);
						if(sImpParXYSig.at(2).is_fromB) fh2dJetSignedImpParXYSignificancebThird_McCorr->Fill(jetpt,sImpParXYSig.at(2).first,sImpParXYSig.at(2).second);
						if(sImpParXYZSig.at(2).is_fromB)fh2dJetSignedImpParXYZSignificancebFirst_McCorr->Fill(jetpt,sImpParXYZSig.at(2).first,sImpParXYZSig.at(2).second);
						if(!sImpParXY.at(2).is_fromB) fh2dJetSignedImpParXYbThird_McCorrNonBDecay->Fill(jetpt,sImpParXY.at(2).first,sImpParXY.at(2).second);
						if(!sImpParXYZ.at(2).is_fromB) fh2dJetSignedImpParXYZbThird_McCorrNonBDecay->Fill(jetpt,sImpParXYZ.at(2).first,sImpParXYZ.at(2).second);
						if(!sImpParXYSig.at(2).is_fromB) fh2dJetSignedImpParXYSignificancebThird_McCorrNonBDecay->Fill(jetpt,sImpParXYSig.at(2).first,sImpParXYSig.at(2).second);
						if(!sImpParXYZSig.at(2).is_fromB) fh2dJetSignedImpParXYZSignificancebFirst_McCorrNonBDecay->Fill(jetpt,sImpParXYZSig.at(2).first,sImpParXYZSig.at(2).second);
						if(sImpParXY.at(2).is_electron) 	fh2dJetSignedImpParXYbThird_electron->Fill(jetpt,sImpParXY.at(2).first,sImpParXY.at(2).second);
						if(sImpParXYZ.at(2).is_electron)	fh2dJetSignedImpParXYZbThird_electron->Fill(jetpt,sImpParXYZ.at(2).first,sImpParXYZ.at(2).second);
						if(sImpParXYSig.at(2).is_electron)	fh2dJetSignedImpParXYSignificancebThird_electron->Fill(jetpt,sImpParXYSig.at(2).first,sImpParXYSig.at(2).second);
						if(sImpParXYZSig.at(2).is_electron)	fh2dJetSignedImpParXYZSignificancebThird_electron->Fill(jetpt,sImpParXYZSig.at(2).first,sImpParXYZSig.at(2).second);
					}
				}
			}

			sImpParXY.clear();
			sImpParXYZ.clear();
			sImpParXYSig.clear();
			sImpParXYZSig.clear();
	}
	return kTRUE;
}


Bool_t AliAnalysisTaskHFJetIPQA::IsSelected(AliVEvent *event, Int_t &WhyRejected,ULong_t &RejectionBits){
	WhyRejected =0;
	Int_t fMinVtxContr=1;
	Int_t fMinVtxType=3;
	Bool_t accept=kTRUE;
	Int_t fMinContrPileup = 5;
	Float_t fMinDzPileup = 0.8;
	Double_t fMaxVtxZ = 10;
	RejectionBits=0;
	//Physics Selection Cut

	    Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
	    if(!isSelected) {
	      if(accept) WhyRejected=7;
	      RejectionBits+=1<<kPhysicsSelection;
	      accept=kFALSE;
	    }

	    // vertex requirements
	    const AliVVertex *vertex = event->GetPrimaryVertex();


	    if(!vertex || vertex->GetNContributors() ==0){
	      accept=kFALSE;
	      if(!vertex)
	    	   RejectionBits+=1<<kNoVertex;
	      else {
	    	    RejectionBits+=1<<kNoContributors;
	      }
	    }else{
	    	const AliVVertex* trkVtx = dynamic_cast<const AliVVertex*>(event->GetPrimaryVertex());
	    	const AliVVertex* spdVtx = dynamic_cast<const AliVVertex*>(event->GetPrimaryVertexSPD());
	    	TString vtxTtl = trkVtx->GetTitle();
	    	if(!vtxTtl.Contains("VertexerTracks"))
	    	{
	    		  accept=kFALSE;
	    		  RejectionBits+=1<<kNoVertexTracks;

	    	}
	    	if(trkVtx->GetNContributors()<2)	{
	    	    		 accept=kFALSE;
	    	    		 RejectionBits+=1<<kTooFewVtxContrib;
	    	    	}
	    	Double_t cov[6] = { 0 };
	    	spdVtx->GetCovarianceMatrix(cov);
	    	Double_t zRes = TMath::Sqrt(cov[5]);
	    	if(spdVtx->IsFromVertexerZ() && (zRes > 0.25)) {
	    		  accept=kFALSE;
	    		  RejectionBits+=1<<kVertexZResolution;
	    	}

	    	if((TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ()) > 0.5)) {
	  		  accept=kFALSE;
	  		  RejectionBits+=1<<kDeltaVertexZ;
	    	}
	    	if(spdVtx->GetNContributors() <1) {
	    		 accept=kFALSE;
	    		RejectionBits+=1<<kVertexZContrib;
	    	}

	    	 if(TMath::Abs(trkVtx->GetZ())>=fMaxVtxZ) {
	    		        RejectionBits+=1<<kZVtxOutFid;
	    		        if(accept) WhyRejected=6;
	    		        accept=kFALSE;
	    	 }
	    }
//Pileup
	        Int_t cutc=(Int_t)fMinContrPileup;
	        Double_t cutz=(Double_t)fMinDzPileup;
	        if(event->IsPileupFromSPD(5, 0.8, 3.0, 2.0, 5.0)) {
	          if(accept) WhyRejected=1;
	          RejectionBits+=1<<kPileupSPD;
	          accept=kFALSE;
	        }
//Special out-of bunch pileup cuts
	    	// SPD Cluster vs Tracklet plot to estimate pileup effect
	    	Int_t nClustersLayer0 = event->GetNumberOfITSClusters(0);
	    	Int_t nClustersLayer1 = event->GetNumberOfITSClusters(1);
	    	Int_t nTracklets = event->GetMultiplicity()->GetNumberOfTracklets();
			if(nClustersLayer0 + nClustersLayer1 > 65 + 4 * nTracklets){
				 accept=kFALSE;
				RejectionBits+=1<<kSPDClusterCut;
			}
			if(fUtils->IsPileUpMV(event)){

				 accept=kFALSE;
				RejectionBits+=1<<kMVPileup;
			}




return accept;
}

// ######################################################################################## Event Selection
Bool_t AliAnalysisTaskHFJetIPQA::IsEventSelected()	{
	AliAODEvent* aev = NULL;

	Int_t WhyRejected =0;
	ULong_t RejectionBits=0;


	if(!fESD)
	{
		aev = dynamic_cast<AliAODEvent*>(InputEvent());
		if(aev && aev->GetPrimaryVertex() && aev->GetPrimaryVertex()->GetNContributors()>0){
			fh1dVertexZ->Fill(aev->GetPrimaryVertex()->GetZ());
			double vtxx = aev->GetPrimaryVertex()->GetX();
			double vtxy = aev->GetPrimaryVertex()->GetY();
			fh1dVertexR->Fill(vtxx,vtxy);
		}else return kFALSE;

/*
 *
 * fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(1,"Accepted");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(2,"Rejected");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(3,"DueToPhysicsSelection");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(4,"DueCentralitySelection");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(5,"ZVertexOutsideFiducialRegion");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(6,"IsEventRejectedDueToNotRecoVertex");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(7,"IsEventRejectedDueToVertexContributors");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(8,"DueToTrigger");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(9,"DueToSPDTrackletClusterCut");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(10,"DueToMVPileup");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(11,"NoVertexTracks");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(12,"NoContributors");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(13,"DeltaVertexZSPDTracks");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(14,"ZVertexResolution");

 */
		if(!(IsSelected(aev,WhyRejected,RejectionBits)))
		{

			fh1dEventRejectionRDHFCuts->SetBinContent(2,fh1dEventRejectionRDHFCuts->GetBinContent(2)+1);
			if(RejectionBits&(1<<kPhysicsSelection))fh1dEventRejectionRDHFCuts->SetBinContent(3,fh1dEventRejectionRDHFCuts->GetBinContent(3)+1);
			else if(RejectionBits&(1<<kNoVertex))fh1dEventRejectionRDHFCuts->SetBinContent(6,fh1dEventRejectionRDHFCuts->GetBinContent(6)+1);
			else if(RejectionBits&(1<<kNoVertexTracks))fh1dEventRejectionRDHFCuts->SetBinContent(11,fh1dEventRejectionRDHFCuts->GetBinContent(11)+1);
			else if(RejectionBits&(1<<kNoContributors))fh1dEventRejectionRDHFCuts->SetBinContent(12,fh1dEventRejectionRDHFCuts->GetBinContent(12)+1);
			else if(RejectionBits&(1<<kTooFewVtxContrib))fh1dEventRejectionRDHFCuts->SetBinContent(7,fh1dEventRejectionRDHFCuts->GetBinContent(7)+1);
			else if(RejectionBits&(1<<kVertexZContrib))fh1dEventRejectionRDHFCuts->SetBinContent(15,fh1dEventRejectionRDHFCuts->GetBinContent(15)+1);
			else if(RejectionBits&(1<<kVertexZResolution))fh1dEventRejectionRDHFCuts->SetBinContent(14,fh1dEventRejectionRDHFCuts->GetBinContent(14)+1);
			else if(RejectionBits&(1<<kDeltaVertexZ))fh1dEventRejectionRDHFCuts->SetBinContent(13,fh1dEventRejectionRDHFCuts->GetBinContent(13)+1);
			else if(RejectionBits&(1<<kZVtxOutFid))fh1dEventRejectionRDHFCuts->SetBinContent(5,fh1dEventRejectionRDHFCuts->GetBinContent(5)+1);
			else if(RejectionBits&(1<<kOutsideCentrality))fh1dEventRejectionRDHFCuts->SetBinContent(4,fh1dEventRejectionRDHFCuts->GetBinContent(4)+1);
			else if(RejectionBits&(1<<kNotSelTrigger))fh1dEventRejectionRDHFCuts->SetBinContent(8,fh1dEventRejectionRDHFCuts->GetBinContent(8)+1);
			else if(RejectionBits&(1<<kSPDClusterCut))fh1dEventRejectionRDHFCuts->SetBinContent(9,fh1dEventRejectionRDHFCuts->GetBinContent(9)+1);
			else if(RejectionBits&(1<<kMVPileup))fh1dEventRejectionRDHFCuts->SetBinContent(10,fh1dEventRejectionRDHFCuts->GetBinContent(10)+1);


			return kFALSE;
		}else {
			fh1dEventRejectionRDHFCuts->SetBinContent(1,fh1dEventRejectionRDHFCuts->GetBinContent(1)+1);
			fh1dVertexZAccepted->Fill(aev->GetPrimaryVertex()->GetZ());
			double vtxx = aev->GetPrimaryVertex()->GetX();
			double vtxy = aev->GetPrimaryVertex()->GetY();
			fh1dVertexRAccepted->Fill(vtxx,vtxy);
			fh2dVertexChi2NDFNESDTracks->Fill(aev->GetPrimaryVertex()->GetChi2perNDF(),aev->GetNumberOfESDTracks());
			return kTRUE;
		}
	}
	AliESDEvent* eev = NULL;
	if(fESD)
	{
		eev = dynamic_cast<AliESDEvent*>(InputEvent());
		if(eev && eev->GetPrimaryVertex() && eev->GetPrimaryVertex()->GetNContributors()>0){
			fh1dVertexZ->Fill(eev->GetPrimaryVertex()->GetZ());
			double vtxx = eev->GetPrimaryVertex()->GetX();
			double vtxy = eev->GetPrimaryVertex()->GetY();
			fh1dVertexR->Fill(vtxx,vtxy);
		}else return kFALSE;

		if(!(IsSelected(eev,WhyRejected,RejectionBits)))
		{
			fh1dEventRejectionRDHFCuts->SetBinContent(2,fh1dEventRejectionRDHFCuts->GetBinContent(2)+1);
					if(RejectionBits&(1<<kPhysicsSelection))fh1dEventRejectionRDHFCuts->SetBinContent(3,fh1dEventRejectionRDHFCuts->GetBinContent(3)+1);
					else if(RejectionBits&(1<<kNoVertex))fh1dEventRejectionRDHFCuts->SetBinContent(6,fh1dEventRejectionRDHFCuts->GetBinContent(6)+1);
					else if(RejectionBits&(1<<kNoVertexTracks))fh1dEventRejectionRDHFCuts->SetBinContent(11,fh1dEventRejectionRDHFCuts->GetBinContent(11)+1);
					else if(RejectionBits&(1<<kNoContributors))fh1dEventRejectionRDHFCuts->SetBinContent(12,fh1dEventRejectionRDHFCuts->GetBinContent(12)+1);
					else if(RejectionBits&(1<<kTooFewVtxContrib))fh1dEventRejectionRDHFCuts->SetBinContent(7,fh1dEventRejectionRDHFCuts->GetBinContent(7)+1);
					else if(RejectionBits&(1<<kVertexZContrib))fh1dEventRejectionRDHFCuts->SetBinContent(15,fh1dEventRejectionRDHFCuts->GetBinContent(15)+1);
					else if(RejectionBits&(1<<kVertexZResolution))fh1dEventRejectionRDHFCuts->SetBinContent(14,fh1dEventRejectionRDHFCuts->GetBinContent(14)+1);
					else if(RejectionBits&(1<<kDeltaVertexZ))fh1dEventRejectionRDHFCuts->SetBinContent(13,fh1dEventRejectionRDHFCuts->GetBinContent(13)+1);
					else if(RejectionBits&(1<<kZVtxOutFid))fh1dEventRejectionRDHFCuts->SetBinContent(5,fh1dEventRejectionRDHFCuts->GetBinContent(5)+1);
					else if(RejectionBits&(1<<kOutsideCentrality))fh1dEventRejectionRDHFCuts->SetBinContent(4,fh1dEventRejectionRDHFCuts->GetBinContent(4)+1);
					else if(RejectionBits&(1<<kNotSelTrigger))fh1dEventRejectionRDHFCuts->SetBinContent(8,fh1dEventRejectionRDHFCuts->GetBinContent(8)+1);
					else if(RejectionBits&(1<<kSPDClusterCut))fh1dEventRejectionRDHFCuts->SetBinContent(9,fh1dEventRejectionRDHFCuts->GetBinContent(9)+1);
					else if(RejectionBits&(1<<kMVPileup))fh1dEventRejectionRDHFCuts->SetBinContent(10,fh1dEventRejectionRDHFCuts->GetBinContent(10)+1);
					return kFALSE;
		}else {
			fh1dEventRejectionRDHFCuts->SetBinContent(1,fh1dEventRejectionRDHFCuts->GetBinContent(1)+1);
			fh1dVertexZAccepted->Fill(eev->GetPrimaryVertex()->GetZ());
			double vtxx = eev->GetPrimaryVertex()->GetX();
			double vtxy = eev->GetPrimaryVertex()->GetY();
			fh1dVertexRAccepted->Fill(vtxx,vtxy);
			fh2dVertexChi2NDFNESDTracks->Fill(eev->GetPrimaryVertex()->GetChi2perNDF(),eev->GetNumberOfTracks());
			return kTRUE;
		}
	}

	return kFALSE;
}
// ######################################################################################## Init histograms
void AliAnalysisTaskHFJetIPQA::UserCreateOutputObjects(){
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (mgr) {
		AliVEventHandler *evhand = mgr->GetInputEventHandler();
		if (evhand) {
			if (evhand->InheritsFrom("AliESDInputHandler")) {
				fIsEsd = kTRUE;
				fESD=kTRUE;
			}
			else {
				fIsEsd = kFALSE;
				fESD=kFALSE;
			}
		}
		else {
			AliError("Event handler not found!");
		}
	}
	else {
		AliError("Analysis manager not found!");
	}


	const Int_t nBins2dSignificance =500;
	const Int_t nBins3dSignificance =500;
	const Int_t nBins2d=250;
	const Int_t nBins3d =250;
	fHFJetUtils = new AliHFJetsTagging("fHFJetUtils");
	fNparents = 7;

	if (!fOutput) fOutput = new AliEmcalList();
	fOutput->SetOwner(kTRUE);

	//Event selection histograms
	fh1dEventRejectionRDHFCuts = new TH1D("fh1dEventRejectionRDHFCuts;reason;count","Rejection reasons",15,0,15);
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(1,"Accepted");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(2,"Rejected");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(3,"DueToPhysicsSelection");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(4,"DueCentralitySelection");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(5,"ZVertexOutsideFiducialRegion");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(6,"IsEventRejectedDueToNotRecoVertex");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(7,"IsEventRejectedDueToVertexContributors");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(8,"DueToTrigger");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(9,"DueToSPDTrackletClusterCut");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(10,"DueToMVPileup");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(11,"NoVertexTracks");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(12,"NoContributorsVertexTracks");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(13,"DeltaVertexZSPDTracks");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(14,"ZVertexResolution");
	fh1dEventRejectionRDHFCuts->GetXaxis()->SetBinLabel(15,"VertexZContributors");



	//Vertex Z before and after
	fh1dVertexZ = new TH1D("fh1dVertexZ","Vertex Z before Event selection;primary vertex z (cm);count",500,-30,30);
	fh1dVertexZAccepted = new TH1D("fh1dVertexZAccepted","Vertex Z after Event selection;primary vertex z (cm);count",500,-30,30);
	//Vertex R before and after
	fh1dVertexR = new TH2D("fh1dVertexR","Vertex R before Event selection;primary vertex xy (cm);x;y",500,-0.5,0.5,500,-0.5,0.5);
	fh1dVertexRAccepted = new TH2D("fh1dVertexRAccepted","Vertex R after Event selection;primary vertex xy (cm);x;y",500,-0.5,0.5,500,-0.5,0.5);

	// Vertex Chi2/NDF vs TPC track multiplicity(ESD tracks)
	fh2dVertexChi2NDFNESDTracks = new TH2D("fh1dVertexChi2NDFNESDTracks","Vertex Chi2/NDF vs # tracks ESD;vertex #chi^{2}/NDF;# tracks esd",200,0,10,500,0,500);
	// AOD tracks accepted
	fh1dTracksAccepeted = new TH1D("fh1dTracksAccepeted","# tracks before/after cuts;;",3,0,3);
	fh1dTracksAccepeted->GetXaxis()->SetBinLabel(1,"total");
	fh1dTracksAccepeted->GetXaxis()->SetBinLabel(2,"accepted");
	fh1dTracksAccepeted->GetXaxis()->SetBinLabel(3,"rejected");
	// Tracks impact parameter histograms
	fh1dTracksImpParXY = new TH1D("fh1dTracksImpParXY","radial imp. parameter ;impact parameter xy (cm);a.u.",nBins2d,-1.,1.);
	fh1dTracksImpParXYZ = new TH1D("fh1dTracksImpParXYZ","3d imp. parameter ;impact parameter 3d (cm);a.u.",nBins3d,0,1.);
	fh1dTracksImpParXYSignificance = new TH1D("fh1dTracksImpParXYSignificance","radial imp. parameter ;impact parameter xy significance;a.u.",nBins2dSignificance,-100,100.);
	fh1dTracksImpParXYZSignificance = new TH1D("fh1dTracksImpParXYZSignificance","3d imp. parameter ;impact parameter 3d significance;a.u.",nBins3dSignificance/2,0.,100.);

	fh1dJetRecEtaPhiAccepted = new TH2D("fh1dJetRecEtaPhiAccepted","detector level jet;#eta;phi",200,-0.5,0.5,200,0.,TMath::TwoPi());

	if (fIsPythia){
		fh1dTracksImpParXY_McCorr = new TH1D("fh1dTracksImpParXY_McCorr","radial imp. parameter (after correction);impact parameter xy (cm);a.u.",nBins2d,-1.,1.);
		fh1dTracksImpParXYZ_McCorr = new TH1D("fh1dTracksImpParXYZ_McCorr","3d imp. parameter (after correction);impact parameter 3d (cm);a.u.",nBins3d,0,1.);
		fh1dTracksImpParXYSignificance_McCorr = new TH1D("fh1dTracksImpParXYSignificance_McCorr","radial imp. parameter (after correction);impact parameter xy significance;a.u.",nBins2dSignificance,-100,100.);
		fh1dTracksImpParXYZSignificance_McCorr = new TH1D("fh1dTracksImpParXYZSignificance_McCorr","3d imp. parameter (after correction);impact parameter 3d significance;a.u.",nBins3dSignificance/2,0.,100.);

		fh1dTracksImpParXY_electron = new TH1D("fh1dTracksImpParXY_electron","radial imp. parameter (electrons);impact parameter xy (cm);a.u.",nBins2d,-1.,1.);
		fh1dTracksImpParXYZ_electron = new TH1D("fh1dTracksImpParXYZ_electron","3d imp. parameter (electrons);impact parameter 3d (cm);a.u.",nBins3d,0,1.);
		fh1dTracksImpParXYSignificance_electron = new TH1D("fh1dTracksImpParXYSignificance_electron","radial imp. parameter (electrons);impact parameter xy significance;a.u.",nBins2dSignificance,-100,100.);
		fh1dTracksImpParXYZSignificance_electron = new TH1D("fh1dTracksImpParXYZSignificance_electron","3d imp. parameter (electrons);impact parameter 3d significance;a.u.",nBins3dSignificance/2,0.,100.);


		fh1dTracksImpParXYTruth = new TH1D("fh1dTracksImpParXYTruth","True: radial imp. parameter ;impact parameter xy (cm);a.u.",nBins2d,-1.,1.);
		fh1dTracksImpParXYZTruth = new TH1D("fh1dTracksImpParXYZTruth","True: 3d imp. parameter ;impact parameter 3d (cm);a.u.",nBins3d,0,1.);
		fh1dTracksImpParXYResidualTruth  = new TH1D ("fh1dTracksImpParXYResidualTruth","Residual radial imp. parameter; #frac{|DCA_{xy}| - |DCA^{Truth}_{xy}|}{#sigma_{xy}} (N#sigma);a.u.",1000,-5,5);
		fh1dTracksImpParXYZResidualTruth  = new TH1D ("fh1dTracksImpParXYZResidualTruth","Residual 3d imp. parameter; #frac{|DCA_{xyz}| - |DCA^{Truth}_{xyz}|}{#sigma_{xyz}} (N#sigma);a.u.",1000,-5,5);
		fh1dJetGenPt = new TH1D("fh1dJetGenPt","generator level jets;pt (GeV/c); count",500,0,250);
		fh1dJetGenPtUnidentified = new TH1D("fh1dJetGenPtUnidentified","generator level jets (no flavour assigned);pt (GeV/c); count",500,0,250);
		fh1dJetGenPtudsg = new TH1D("fh1dJetGenPtudsg","generator level udsg jets;pt (GeV/c); count",500,0,250);
		fh1dJetGenPtc = new TH1D("fh1dJetGenPtc","generator level c jets;pt (GeV/c); count",500,0,250);
		fh1dJetGenPtb = new TH1D("fh1dJetGenPtb","generator level b jets;pt (GeV/c); count",500,0,250);

		fh2dJetGenPtVsJetRecPt = new TH2D("fh2dJetGenPtVsJetRecPt","detector momentum response;gen pt;rec pt",500,0,250,500,0,250);

		fh2dJetSignedImpParXYUnidentified = new TH2D("fh2dJetSignedImpParXYUnidentified","fh2dJetSignedImpParXYZSignificanceUnidentified;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZUnidentified = new TH2D("fh2dJetSignedImpParXYZUnidentified","fh2dJetSignedImpParXYZUnidentified;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceUnidentified = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentified","fh2dJetSignedImpParXYSignificanceUnidentified;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceUnidentified = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentified","fh2dJetSignedImpParXYZSignificanceUnidentified;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYudsg = new TH2D("fh2dJetSignedImpParXYudsg","fh2dJetSignedImpParXYudsg;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZudsg  = new TH2D("fh2dJetSignedImpParXYZudsg","fh2dJetSignedImpParXYZudsg;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceudsg  = new TH2D("fh2dJetSignedImpParXYSignificanceudsg","fh2dJetSignedImpParXYSignificanceudsg;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceudsg  = new TH2D("fh2dJetSignedImpParXYZSignificanceudsg","fh2dJetSignedImpParXYZSignificanceudsg;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYc= new TH2D("fh2dJetSignedImpParXYc","fh2dJetSignedImpParXYc;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZc  = new TH2D("fh2dJetSignedImpParXYZc","fh2dJetSignedImpParXYZc;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancec  = new TH2D("fh2dJetSignedImpParXYSignificancec","fh2dJetSignedImpParXYSignificancec;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancec  = new TH2D("fh2dJetSignedImpParXYZSignificancec","fh2dJetSignedImpParXYZSignificancec;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYb= new TH2D("fh2dJetSignedImpParXYb","fh2dJetSignedImpParXYb;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZb  = new TH2D("fh2dJetSignedImpParXYZb","fh2dJetSignedImpParXYZb;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceb  = new TH2D("fh2dJetSignedImpParXYSignificanceb","fh2dJetSignedImpParXYSignificanceb;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceb  = new TH2D("fh2dJetSignedImpParXYZSignificanceb","fh2dJetSignedImpParXYZSignificanceb;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
		//N=1
		fh2dJetSignedImpParXYUnidentifiedFirst= new TH2D("fh2dJetSignedImpParXYUnidentifiedFirst","fh2dJetSignedImpParXYUnidentifiedFirst;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZUnidentifiedFirst = new TH2D("fh2dJetSignedImpParXYZUnidentifiedFirst","fh2dJetSignedImpParXYZUnidentifiedFirst;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceUnidentifiedFirst = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentifiedFirst","fh2dJetSignedImpParXYSignificanceUnidentifiedFirst;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst","fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYudsgFirst = new TH2D("fh2dJetSignedImpParXYudsgFirst","fh2dJetSignedImpParXYudsgFirst;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZudsgFirst  = new TH2D("fh2dJetSignedImpParXYZudsgFirst","fh2dJetSignedImpParXYZudsgFirst;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceudsgFirst  = new TH2D("fh2dJetSignedImpParXYSignificanceudsgFirst","fh2dJetSignedImpParXYSignificanceudsgFirst;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceudsgFirst  = new TH2D("fh2dJetSignedImpParXYZSignificanceudsgFirst","fh2dJetSignedImpParXYZSignificanceudsgFirst;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYcFirst= new TH2D("fh2dJetSignedImpParXYcFirst","fh2dJetSignedImpParXYcFirst;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZcFirst  = new TH2D("fh2dJetSignedImpParXYZcFirst","fh2dJetSignedImpParXYZcFirst;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancecFirst  = new TH2D("fh2dJetSignedImpParXYSignificancecFirst","fh2dJetSignedImpParXYSignificancecFirst;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancecFirst  = new TH2D("fh2dJetSignedImpParXYZSignificancecFirst","fh2dJetSignedImpParXYZSignificancecFirst;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYbFirst= new TH2D("fh2dJetSignedImpParXYbFirst","fh2dJetSignedImpParXYbFirst;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZbFirst  = new TH2D("fh2dJetSignedImpParXYZbFirst","fh2dJetSignedImpParXYZbFirst;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancebFirst  = new TH2D("fh2dJetSignedImpParXYSignificancebFirst","fh2dJetSignedImpParXYSignificancebFirst;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancebFirst  = new TH2D("fh2dJetSignedImpParXYZSignificancebFirst","fh2dJetSignedImpParXYZSignificancebFirst;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		//N=2
		fh2dJetSignedImpParXYUnidentifiedSecond = new TH2D("fh2dJetSignedImpParXYUnidentifiedSecond","fh2dJetSignedImpParXYUnidentifiedSecond;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZUnidentifiedSecond = new TH2D("fh2dJetSignedImpParXYZUnidentifiedSecond","fh2dJetSignedImpParXYZUnidentifiedSecond;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceUnidentifiedSecond = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentifiedSecond","fh2dJetSignedImpParXYSignificanceUnidentifiedSecond;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond","fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYudsgSecond = new TH2D("fh2dJetSignedImpParXYudsgSecond","fh2dJetSignedImpParXYudsgSecond;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZudsgSecond  = new TH2D("fh2dJetSignedImpParXYZudsgSecond","fh2dJetSignedImpParXYZudsgSecond;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceudsgSecond  = new TH2D("fh2dJetSignedImpParXYSignificanceudsgSecond","fh2dJetSignedImpParXYSignificanceudsgSecond;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceudsgSecond  = new TH2D("fh2dJetSignedImpParXYZSignificanceudsgSecond","fh2dJetSignedImpParXYZSignificanceudsgSecond;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYcSecond= new TH2D("fh2dJetSignedImpParXYcSecond","fh2dJetSignedImpParXYcSecond;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZcSecond  = new TH2D("fh2dJetSignedImpParXYZcSecond","fh2dJetSignedImpParXYZcSecond;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancecSecond  = new TH2D("fh2dJetSignedImpParXYSignificancecSecond","fh2dJetSignedImpParXYSignificancecSecond;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancecSecond  = new TH2D("fh2dJetSignedImpParXYZSignificancecSecond","fh2dJetSignedImpParXYZSignificancecSecond;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYbSecond= new TH2D("fh2dJetSignedImpParXYbSecond","fh2dJetSignedImpParXYbSecond;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZbSecond  = new TH2D("fh2dJetSignedImpParXYZbSecond","fh2dJetSignedImpParXYZbSecond;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancebSecond  = new TH2D("fh2dJetSignedImpParXYSignificancebSecond","fh2dJetSignedImpParXYSignificancebSecond;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancebSecond  = new TH2D("fh2dJetSignedImpParXYZSignificancebSecond","fh2dJetSignedImpParXYZSignificancebSecond;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
		//N=3
		fh2dJetSignedImpParXYUnidentifiedThird = new TH2D("fh2dJetSignedImpParXYUnidentifiedThird","fh2dJetSignedImpParXYUnidentifiedThird;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZUnidentifiedThird = new TH2D("fh2dJetSignedImpParXYZUnidentifiedThird","fh2dJetSignedImpParXYZUnidentifiedThird;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceUnidentifiedThird = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentifiedThird","fh2dJetSignedImpParXYSignificanceUnidentifiedThird;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceUnidentifiedThird = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentifiedThird","fh2dJetSignedImpParXYZSignificanceUnidentifiedThird;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYudsgThird = new TH2D("fh2dJetSignedImpParXYudsgThird","fh2dJetSignedImpParXYudsgThird;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZudsgThird  = new TH2D("fh2dJetSignedImpParXYZudsgThird","fh2dJetSignedImpParXYZudsgThird;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceudsgThird  = new TH2D("fh2dJetSignedImpParXYSignificanceudsgThird","fh2dJetSignedImpParXYSignificanceudsgThird;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceudsgThird  = new TH2D("fh2dJetSignedImpParXYZSignificanceudsgThird","fh2dJetSignedImpParXYZSignificanceudsgThird;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYcThird= new TH2D("fh2dJetSignedImpParXYcThird","fh2dJetSignedImpParXYcThird;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZcThird  = new TH2D("fh2dJetSignedImpParXYZcThird","fh2dJetSignedImpParXYZcThird;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancecThird  = new TH2D("fh2dJetSignedImpParXYSignificancecThird","fh2dJetSignedImpParXYSignificancecThird;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancecThird  = new TH2D("fh2dJetSignedImpParXYZSignificancecThird","fh2dJetSignedImpParXYZSignificancecThird;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYbThird= new TH2D("fh2dJetSignedImpParXYbThird","fh2dJetSignedImpParXYbThird;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZbThird  = new TH2D("fh2dJetSignedImpParXYZbThird","fh2dJetSignedImpParXYZbThird;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancebThird  = new TH2D("fh2dJetSignedImpParXYSignificancebThird","fh2dJetSignedImpParXYSignificancebThird;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancebThird  = new TH2D("fh2dJetSignedImpParXYZSignificancebThird","fh2dJetSignedImpParXYZSignificancebThird;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
		// Histograms with Monte Carlo correction factors
		fh2dJetSignedImpParXYUnidentified_McCorr = new TH2D("fh2dJetSignedImpParXYUnidentified_McCorr","fh2dJetSignedImpParXYZSignificanceUnidentified (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZUnidentified_McCorr = new TH2D("fh2dJetSignedImpParXYZUnidentified_McCorr","fh2dJetSignedImpParXYZUnidentified (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,1000,-1,1);
		fh2dJetSignedImpParXYSignificanceUnidentified_McCorr = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentified_McCorr","fh2dJetSignedImpParXYSignificanceUnidentified (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceUnidentified_McCorr = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentified_McCorr","fh2dJetSignedImpParXYZSignificanceUnidentified (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYudsg_McCorr = new TH2D("fh2dJetSignedImpParXYudsg_McCorr","fh2dJetSignedImpParXYudsg (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZudsg_McCorr  = new TH2D("fh2dJetSignedImpParXYZudsg_McCorr","fh2dJetSignedImpParXYZudsg (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceudsg_McCorr  = new TH2D("fh2dJetSignedImpParXYSignificanceudsg_McCorr","fh2dJetSignedImpParXYSignificanceudsg (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceudsg_McCorr  = new TH2D("fh2dJetSignedImpParXYZSignificanceudsg_McCorr","fh2dJetSignedImpParXYZSignificanceudsg (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYc_McCorr= new TH2D("fh2dJetSignedImpParXYc_McCorr","fh2dJetSignedImpParXYc (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZc_McCorr  = new TH2D("fh2dJetSignedImpParXYZc_McCorr","fh2dJetSignedImpParXYZc (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancec_McCorr  = new TH2D("fh2dJetSignedImpParXYSignificancec_McCorr","fh2dJetSignedImpParXYSignificancec (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancec_McCorr  = new TH2D("fh2dJetSignedImpParXYZSignificancec_McCorr","fh2dJetSignedImpParXYZSignificancec (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYb_McCorr= new TH2D("fh2dJetSignedImpParXYb_McCorr","fh2dJetSignedImpParXYb (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYb_McCorrNonBDecay= new TH2D("fh2dJetSignedImpParXYb_McCorrNonBDecay","fh2dJetSignedImpParXYb (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZb_McCorr  = new TH2D("fh2dJetSignedImpParXYZb_McCorr","fh2dJetSignedImpParXYZb (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYZb_McCorr_McCorrNonBDecay  = new TH2D("fh2dJetSignedImpParXYZb_McCorr_McCorrNonBDecay","fh2dJetSignedImpParXYZb (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceb_McCorr  = new TH2D("fh2dJetSignedImpParXYSignificanceb_McCorr","fh2dJetSignedImpParXYSignificanceb (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYSignificanceb_McCorrNonBDecay  = new TH2D("fh2dJetSignedImpParXYSignificanceb_McCorrNonBDecay","fh2dJetSignedImpParXYSignificanceb (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceb_McCorr  = new TH2D("fh2dJetSignedImpParXYZSignificanceb_McCorr","fh2dJetSignedImpParXYZSignificanceb (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceb_McCorrNonBDecay  = new TH2D("fh2dJetSignedImpParXYZSignificanceb_McCorrNonBDecay","fh2dJetSignedImpParXYZSignificanceb (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		//N=1
		fh2dJetSignedImpParXYUnidentifiedFirst_McCorr= new TH2D("fh2dJetSignedImpParXYUnidentifiedFirst_McCorr","fh2dJetSignedImpParXYUnidentifiedFirst (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZUnidentifiedFirst_McCorr = new TH2D("fh2dJetSignedImpParXYZUnidentifiedFirst_McCorr","fh2dJetSignedImpParXYZUnidentifiedFirst (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,1000,-1,1);
		fh2dJetSignedImpParXYSignificanceUnidentifiedFirst_McCorr = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentifiedFirst_McCorr","fh2dJetSignedImpParXYSignificanceUnidentifiedFirst (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst_McCorr = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst_McCorr","fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYudsgFirst_McCorr = new TH2D("fh2dJetSignedImpParXYudsgFirst_McCorr","fh2dJetSignedImpParXYudsgFirst (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZudsgFirst_McCorr  = new TH2D("fh2dJetSignedImpParXYZudsgFirst_McCorr","fh2dJetSignedImpParXYZudsgFirst (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceudsgFirst_McCorr  = new TH2D("fh2dJetSignedImpParXYSignificanceudsgFirst_McCorr","fh2dJetSignedImpParXYSignificanceudsgFirst (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceudsgFirst_McCorr  = new TH2D("fh2dJetSignedImpParXYZSignificanceudsgFirst_McCorr","fh2dJetSignedImpParXYZSignificanceudsgFirst (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYcFirst_McCorr= new TH2D("fh2dJetSignedImpParXYcFirst_McCorr","fh2dJetSignedImpParXYcFirst (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZcFirst_McCorr  = new TH2D("fh2dJetSignedImpParXYZcFirst_McCorr","fh2dJetSignedImpParXYZcFirst (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancecFirst_McCorr  = new TH2D("fh2dJetSignedImpParXYSignificancecFirst_McCorr","fh2dJetSignedImpParXYSignificancecFirst (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancecFirst_McCorr  = new TH2D("fh2dJetSignedImpParXYZSignificancecFirst_McCorr","fh2dJetSignedImpParXYZSignificancecFirst (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYbFirst_McCorr= new TH2D("fh2dJetSignedImpParXYbFirst_McCorr","fh2dJetSignedImpParXYbFirst (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYbFirst_McCorrNonBDecay= new TH2D("fh2dJetSignedImpParXYbFirst_McCorrNonBDecay","fh2dJetSignedImpParXYbFirst (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZbFirst_McCorr  = new TH2D("fh2dJetSignedImpParXYZbFirst_McCorr","fh2dJetSignedImpParXYZbFirst (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYZbFirst_McCorrNonBDecay  = new TH2D("fh2dJetSignedImpParXYZbFirst_McCorrNonBDecay","fh2dJetSignedImpParXYZbFirst (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancebFirst_McCorr  = new TH2D("fh2dJetSignedImpParXYSignificancebFirst_McCorr","fh2dJetSignedImpParXYSignificancebFirst (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYSignificancebFirst_McCorrNonBDecay  = new TH2D("fh2dJetSignedImpParXYSignificancebFirst_McCorrNonBDecay","fh2dJetSignedImpParXYSignificancebFirst (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancebFirst_McCorr  = new TH2D("fh2dJetSignedImpParXYZSignificancebFirst_McCorr","fh2dJetSignedImpParXYZSignificancebFirst (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancebFirst_McCorrNonBDecay  = new TH2D("fh2dJetSignedImpParXYZSignificancebFirst_McCorrNonBDecay","fh2dJetSignedImpParXYZSignificancebFirst (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		//N=2
		fh2dJetSignedImpParXYUnidentifiedSecond_McCorr = new TH2D("fh2dJetSignedImpParXYUnidentifiedSecond_McCorr","fh2dJetSignedImpParXYUnidentifiedSecond (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZUnidentifiedSecond_McCorr = new TH2D("fh2dJetSignedImpParXYZUnidentifiedSecond_McCorr","fh2dJetSignedImpParXYZUnidentifiedSecond (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceUnidentifiedSecond_McCorr = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentifiedSecond_McCorr","fh2dJetSignedImpParXYSignificanceUnidentifiedSecond (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond_McCorr = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond_McCorr","fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYudsgSecond_McCorr = new TH2D("fh2dJetSignedImpParXYudsgSecond_McCorr","fh2dJetSignedImpParXYudsgSecond (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZudsgSecond_McCorr = new TH2D("fh2dJetSignedImpParXYZudsgSecond_McCorr","fh2dJetSignedImpParXYZudsgSecond (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceudsgSecond_McCorr  = new TH2D("fh2dJetSignedImpParXYSignificanceudsgSecond_McCorr","fh2dJetSignedImpParXYSignificanceudsgSecond (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceudsgSecond_McCorr  = new TH2D("fh2dJetSignedImpParXYZSignificanceudsgSecond_McCorr","fh2dJetSignedImpParXYZSignificanceudsgSecond (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYcSecond_McCorr= new TH2D("fh2dJetSignedImpParXYcSecond_McCorr","fh2dJetSignedImpParXYcSecond (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZcSecond_McCorr  = new TH2D("fh2dJetSignedImpParXYZcSecond_McCorr","fh2dJetSignedImpParXYZcSecond (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancecSecond_McCorr  = new TH2D("fh2dJetSignedImpParXYSignificancecSecond_McCorr","fh2dJetSignedImpParXYSignificancecSecond (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancecSecond_McCorr  = new TH2D("fh2dJetSignedImpParXYZSignificancecSecond_McCorr","fh2dJetSignedImpParXYZSignificancecSecond (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);


		fh2dJetSignedImpParXYbSecond_McCorr= new TH2D("fh2dJetSignedImpParXYbSecond_McCorr","fh2dJetSignedImpParXYbSecond (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYbSecond_McCorrNonBDecay= new TH2D("fh2dJetSignedImpParXYbSecond_McCorrNonBDecay","fh2dJetSignedImpParXYbSecond (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZbSecond_McCorr  = new TH2D("fh2dJetSignedImpParXYZbSecond_McCorr","fh2dJetSignedImpParXYZbSecond (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYZbSecond_McCorrNonBDecay  = new TH2D("fh2dJetSignedImpParXYZbSecond_McCorrNonBDecay","fh2dJetSignedImpParXYZbSecond (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancebSecond_McCorr  = new TH2D("fh2dJetSignedImpParXYSignificancebSecond_McCorr","fh2dJetSignedImpParXYSignificancebSecond (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYSignificancebSecond_McCorrNonBDecay  = new TH2D("fh2dJetSignedImpParXYSignificancebSecond_McCorrNonBDecay","fh2dJetSignedImpParXYSignificancebSecond (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancebSecond_McCorr  = new TH2D("fh2dJetSignedImpParXYZSignificancebSecond_McCorr","fh2dJetSignedImpParXYZSignificancebSecond (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancebSecond_McCorrNonBDecay  = new TH2D("fh2dJetSignedImpParXYZSignificancebSecond_McCorrNonBDecay","fh2dJetSignedImpParXYZSignificancebSecond (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);


		//N=3
		fh2dJetSignedImpParXYUnidentifiedThird_McCorr = new TH2D("fh2dJetSignedImpParXYUnidentifiedThird_McCorr","fh2dJetSignedImpParXYUnidentifiedThird (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZUnidentifiedThird_McCorr = new TH2D("fh2dJetSignedImpParXYZUnidentifiedThird_McCorr","fh2dJetSignedImpParXYZUnidentifiedThird (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceUnidentifiedThird_McCorr = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentifiedThird_McCorr","fh2dJetSignedImpParXYSignificanceUnidentifiedThird (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceUnidentifiedThird_McCorr = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentifiedThird_McCorr","fh2dJetSignedImpParXYZSignificanceUnidentifiedThird (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYudsgThird_McCorr = new TH2D("fh2dJetSignedImpParXYudsgThird_McCorr","fh2dJetSignedImpParXYudsgThird (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZudsgThird_McCorr  = new TH2D("fh2dJetSignedImpParXYZudsgThird_McCorr","fh2dJetSignedImpParXYZudsgThird (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceudsgThird_McCorr  = new TH2D("fh2dJetSignedImpParXYSignificanceudsgThird_McCorr","fh2dJetSignedImpParXYSignificanceudsgThird (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceudsgThird_McCorr  = new TH2D("fh2dJetSignedImpParXYZSignificanceudsgThird_McCorr","fh2dJetSignedImpParXYZSignificanceudsgThird (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYcThird_McCorr= new TH2D("fh2dJetSignedImpParXYcThird_McCorr","fh2dJetSignedImpParXYcThird (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZcThird_McCorr  = new TH2D("fh2dJetSignedImpParXYZcThird_McCorr","fh2dJetSignedImpParXYZcThird (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancecThird_McCorr  = new TH2D("fh2dJetSignedImpParXYSignificancecThird_McCorr","fh2dJetSignedImpParXYSignificancecThird (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancecThird_McCorr  = new TH2D("fh2dJetSignedImpParXYZSignificancecThird_McCorr","fh2dJetSignedImpParXYZSignificancecThird (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);


		fh2dJetSignedImpParXYbThird_McCorr= new TH2D("fh2dJetSignedImpParXYbThird_McCorr","fh2dJetSignedImpParXYbThird (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYbThird_McCorrNonBDecay= new TH2D("fh2dJetSignedImpParXYbThird_McCorrNonBDecay","fh2dJetSignedImpParXYbThird (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);

		fh2dJetSignedImpParXYZbThird_McCorr = new TH2D("fh2dJetSignedImpParXYZbThird_McCorr","fh2dJetSignedImpParXYZbThird (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYZbThird_McCorrNonBDecay = new TH2D("fh2dJetSignedImpParXYZbThird_McCorrNonBDecay","fh2dJetSignedImpParXYZbThird (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);

		fh2dJetSignedImpParXYSignificancebThird_McCorr  = new TH2D("fh2dJetSignedImpParXYSignificancebThird_McCorr","fh2dJetSignedImpParXYSignificancebThird (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYSignificancebThird_McCorrNonBDecay  = new TH2D("fh2dJetSignedImpParXYSignificancebThird_McCorrNonBDecay","fh2dJetSignedImpParXYSignificancebThird (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancebThird_McCorr  = new TH2D("fh2dJetSignedImpParXYZSignificancebThird_McCorr","fh2dJetSignedImpParXYZSignificancebThird (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancebThird_McCorrNonBDecay  = new TH2D("fh2dJetSignedImpParXYZSignificancebThird_McCorrNonBDecay","fh2dJetSignedImpParXYZSignificancebThird (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);





		//electrons
		// Histograms with Monte Carlo correction factors
		fh2dJetSignedImpParXYUnidentified_electron = new TH2D("fh2dJetSignedImpParXYUnidentified_electron","fh2dJetSignedImpParXYZSignificanceUnidentified (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZUnidentified_electron = new TH2D("fh2dJetSignedImpParXYZUnidentified_electron","fh2dJetSignedImpParXYZUnidentified (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceUnidentified_electron = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentified_electron","fh2dJetSignedImpParXYSignificanceUnidentified (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceUnidentified_electron = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentified_electron","fh2dJetSignedImpParXYZSignificanceUnidentified (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYudsg_electron = new TH2D("fh2dJetSignedImpParXYudsg_electron","fh2dJetSignedImpParXYudsg (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZudsg_electron  = new TH2D("fh2dJetSignedImpParXYZudsg_electron","fh2dJetSignedImpParXYZudsg (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceudsg_electron  = new TH2D("fh2dJetSignedImpParXYSignificanceudsg_electron","fh2dJetSignedImpParXYSignificanceudsg (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceudsg_electron  = new TH2D("fh2dJetSignedImpParXYZSignificanceudsg_electron","fh2dJetSignedImpParXYZSignificanceudsg (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYc_electron= new TH2D("fh2dJetSignedImpParXYc_electron","fh2dJetSignedImpParXYc (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZc_electron  = new TH2D("fh2dJetSignedImpParXYZc_electron","fh2dJetSignedImpParXYZc (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancec_electron  = new TH2D("fh2dJetSignedImpParXYSignificancec_electron","fh2dJetSignedImpParXYSignificancec (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancec_electron  = new TH2D("fh2dJetSignedImpParXYZSignificancec_electron","fh2dJetSignedImpParXYZSignificancec (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYb_electron= new TH2D("fh2dJetSignedImpParXYb_electron","fh2dJetSignedImpParXYb (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZb_electron  = new TH2D("fh2dJetSignedImpParXYZb_electron","fh2dJetSignedImpParXYZb (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceb_electron  = new TH2D("fh2dJetSignedImpParXYSignificanceb_electron","fh2dJetSignedImpParXYSignificanceb (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceb_electron  = new TH2D("fh2dJetSignedImpParXYZSignificanceb_electron","fh2dJetSignedImpParXYZSignificanceb (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);


		//N=1
		fh2dJetSignedImpParXYUnidentifiedFirst_electron= new TH2D("fh2dJetSignedImpParXYUnidentifiedFirst_electron","fh2dJetSignedImpParXYUnidentifiedFirst (electrons);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZUnidentifiedFirst_electron = new TH2D("fh2dJetSignedImpParXYZUnidentifiedFirst_electron","fh2dJetSignedImpParXYZUnidentifiedFirst (electrons);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceUnidentifiedFirst_electron = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentifiedFirst_electron","fh2dJetSignedImpParXYSignificanceUnidentifiedFirst (electrons);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst_electron = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst_electron","fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst (electrons);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYudsgFirst_electron = new TH2D("fh2dJetSignedImpParXYudsgFirst_electron","fh2dJetSignedImpParXYudsgFirst (electrons);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZudsgFirst_electron  = new TH2D("fh2dJetSignedImpParXYZudsgFirst_electron","fh2dJetSignedImpParXYZudsgFirst (electrons);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceudsgFirst_electron  = new TH2D("fh2dJetSignedImpParXYSignificanceudsgFirst_electron","fh2dJetSignedImpParXYSignificanceudsgFirst (electrons);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceudsgFirst_electron  = new TH2D("fh2dJetSignedImpParXYZSignificanceudsgFirst_electron","fh2dJetSignedImpParXYZSignificanceudsgFirst (electrons);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYcFirst_electron= new TH2D("fh2dJetSignedImpParXYcFirst_electron","fh2dJetSignedImpParXYcFirst (electrons);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZcFirst_electron  = new TH2D("fh2dJetSignedImpParXYZcFirst_electron","fh2dJetSignedImpParXYZcFirst (electrons);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancecFirst_electron  = new TH2D("fh2dJetSignedImpParXYSignificancecFirst_electron","fh2dJetSignedImpParXYSignificancecFirst (electrons);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancecFirst_electron  = new TH2D("fh2dJetSignedImpParXYZSignificancecFirst_electron","fh2dJetSignedImpParXYZSignificancecFirst (electrons);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYbFirst_electron= new TH2D("fh2dJetSignedImpParXYbFirst_electron","fh2dJetSignedImpParXYbFirst (electrons);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZbFirst_electron  = new TH2D("fh2dJetSignedImpParXYZbFirst_electron","fh2dJetSignedImpParXYZbFirst (electrons);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancebFirst_electron  = new TH2D("fh2dJetSignedImpParXYSignificancebFirst_electron","fh2dJetSignedImpParXYSignificancebFirst (electrons);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancebFirst_electron  = new TH2D("fh2dJetSignedImpParXYZSignificancebFirst_electron","fh2dJetSignedImpParXYZSignificancebFirst (electrons);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		//N=2
		fh2dJetSignedImpParXYUnidentifiedSecond_electron = new TH2D("fh2dJetSignedImpParXYUnidentifiedSecond_electron","fh2dJetSignedImpParXYUnidentifiedSecond (electrons);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZUnidentifiedSecond_electron = new TH2D("fh2dJetSignedImpParXYZUnidentifiedSecond_electron","fh2dJetSignedImpParXYZUnidentifiedSecond (electrons);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceUnidentifiedSecond_electron = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentifiedSecond_electron","fh2dJetSignedImpParXYSignificanceUnidentifiedSecond (electrons);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond_electron = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond_electron","fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond (electrons);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYudsgSecond_electron = new TH2D("fh2dJetSignedImpParXYudsgSecond_electron","fh2dJetSignedImpParXYudsgSecond (electrons);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZudsgSecond_electron = new TH2D("fh2dJetSignedImpParXYZudsgSecond_electron","fh2dJetSignedImpParXYZudsgSecond (electrons);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceudsgSecond_electron  = new TH2D("fh2dJetSignedImpParXYSignificanceudsgSecond_electron","fh2dJetSignedImpParXYSignificanceudsgSecond (electrons);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceudsgSecond_electron  = new TH2D("fh2dJetSignedImpParXYZSignificanceudsgSecond_electron","fh2dJetSignedImpParXYZSignificanceudsgSecond (electrons);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYcSecond_electron= new TH2D("fh2dJetSignedImpParXYcSecond_electron","fh2dJetSignedImpParXYcSecond (electrons);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZcSecond_electron  = new TH2D("fh2dJetSignedImpParXYZcSecond_electron","fh2dJetSignedImpParXYZcSecond (electrons);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancecSecond_electron  = new TH2D("fh2dJetSignedImpParXYSignificancecSecond_electron","fh2dJetSignedImpParXYSignificancecSecond (electrons);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancecSecond_electron  = new TH2D("fh2dJetSignedImpParXYZSignificancecSecond_electron","fh2dJetSignedImpParXYZSignificancecSecond (electrons);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYbSecond_electron= new TH2D("fh2dJetSignedImpParXYbSecond_electron","fh2dJetSignedImpParXYbSecond (electrons);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZbSecond_electron  = new TH2D("fh2dJetSignedImpParXYZbSecond_electron","fh2dJetSignedImpParXYZbSecond (electrons);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancebSecond_electron  = new TH2D("fh2dJetSignedImpParXYSignificancebSecond_electron","fh2dJetSignedImpParXYSignificancebSecond (electrons);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancebSecond_electron  = new TH2D("fh2dJetSignedImpParXYZSignificancebSecond_electron","fh2dJetSignedImpParXYZSignificancebSecond (electrons);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
		//N=3
		fh2dJetSignedImpParXYUnidentifiedThird_electron = new TH2D("fh2dJetSignedImpParXYUnidentifiedThird_electron","fh2dJetSignedImpParXYUnidentifiedThird (electrons);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZUnidentifiedThird_electron = new TH2D("fh2dJetSignedImpParXYZUnidentifiedThird_electron","fh2dJetSignedImpParXYZUnidentifiedThird (electrons);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceUnidentifiedThird_electron = new TH2D("fh2dJetSignedImpParXYSignificanceUnidentifiedThird_electron","fh2dJetSignedImpParXYSignificanceUnidentifiedThird (electrons);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceUnidentifiedThird_electron = new TH2D("fh2dJetSignedImpParXYZSignificanceUnidentifiedThird_electron","fh2dJetSignedImpParXYZSignificanceUnidentifiedThird (electrons);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYudsgThird_electron = new TH2D("fh2dJetSignedImpParXYudsgThird_electron","fh2dJetSignedImpParXYudsgThird (electrons);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZudsgThird_electron  = new TH2D("fh2dJetSignedImpParXYZudsgThird_electron","fh2dJetSignedImpParXYZudsgThird (electrons);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceudsgThird_electron  = new TH2D("fh2dJetSignedImpParXYSignificanceudsgThird_electron","fh2dJetSignedImpParXYSignificanceudsgThird (electrons);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceudsgThird_electron  = new TH2D("fh2dJetSignedImpParXYZSignificanceudsgThird_electron","fh2dJetSignedImpParXYZSignificanceudsgThird (electrons);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYcThird_electron= new TH2D("fh2dJetSignedImpParXYcThird_electron","fh2dJetSignedImpParXYcThird (electrons);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZcThird_electron  = new TH2D("fh2dJetSignedImpParXYZcThird_electron","fh2dJetSignedImpParXYZcThird (electrons);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancecThird_electron  = new TH2D("fh2dJetSignedImpParXYSignificancecThird_electron","fh2dJetSignedImpParXYSignificancecThird (electrons);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancecThird_electron  = new TH2D("fh2dJetSignedImpParXYZSignificancecThird_electron","fh2dJetSignedImpParXYZSignificancecThird (electrons);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		fh2dJetSignedImpParXYbThird_electron= new TH2D("fh2dJetSignedImpParXYbThird_electron","fh2dJetSignedImpParXYbThird (electrons);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZbThird_electron = new TH2D("fh2dJetSignedImpParXYZbThird_electron","fh2dJetSignedImpParXYZbThird (electrons);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificancebThird_electron  = new TH2D("fh2dJetSignedImpParXYSignificancebThird_electron","fh2dJetSignedImpParXYSignificancebThird (electrons);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificancebThird_electron  = new TH2D("fh2dJetSignedImpParXYZSignificancebThird_electron","fh2dJetSignedImpParXYZSignificancebThird (electrons);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

	}
	// Jet histograms
	fh1dJetRecPt = new TH1D("fh1dJetRecPt","detector level jets;pt (GeV/c); count",500,0,250);
	fh1dJetRecPtAccepted = new TH1D("fh1dJetRecPtAccepted","accepted detector level jets;pt (GeV/c); count",500,0,250);

	if(fIsPythia){
		fh1dJetRecPtUnidentified = new TH1D("fh1dJetRecPtUnidentified","detector level jets;pt (GeV/c); count",500,0,250);
		fh1dJetRecPtudsg = new TH1D("fh1dJetRecPtudsg","detector level jets;pt (GeV/c); count",500,0,250);
		fh1dJetRecPtc = new TH1D("fh1dJetRecPtc","detector level jets;pt (GeV/c); count",500,0,250);
		fh1dJetRecPtb = new TH1D("fh1dJetRecPtb","detector level jets;pt (GeV/c); count",500,0,250);
		fh1dJetRecPtUnidentifiedAccepted = new TH1D("fh1dJetRecPtUnidentifiedAccepted","detector level jets;pt (GeV/c); count",500,0,250);
		fh1dJetRecPtudsgAccepted = new TH1D("fh1dJetRecPtudsgAccepted","detector level jets;pt (GeV/c); count",500,0,250);
		fh1dJetRecPtcAccepted= new TH1D("fh1dJetRecPtcAccepted","detector level jets;pt (GeV/c); count",500,0,250);
		fh1dJetRecPtbAccepted = new TH1D("fh1dJetRecPtbAccepted","detector level jets;pt (GeV/c); count",500,0,250);
	}


	fh2dJetSignedImpParXY = new TH2D("fh2dJetSignedImpParXY","fh2dJetSignedImpParXY;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
	fh2dJetSignedImpParXYZ = new TH2D("fh2dJetSignedImpParXYZ","fh2dJetSignedImpParXYZ;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
	fh2dJetSignedImpParXYSignificance = new TH2D("fh2dJetSignedImpParXYSignificance","fh2dJetSignedImpParXYSignificance;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
	fh2dJetSignedImpParXYZSignificance = new TH2D("fh2dJetSignedImpParXYZSignificance","fh2dJetSignedImpParXYZSignificance;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
	//N=1
	fh2dJetSignedImpParXYFirst = new TH2D("fh2dJetSignedImpParXYFirst","fh2dJetSignedImpParXYFirst;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
	fh2dJetSignedImpParXYZFirst = new TH2D("fh2dJetSignedImpParXYZFirst","fh2dJetSignedImpParXYZFirst;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
	fh2dJetSignedImpParXYSignificanceFirst = new TH2D("fh2dJetSignedImpParXYSignificanceFirst","fh2dJetSignedImpParXYSignificanceFirst;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
	fh2dJetSignedImpParXYZSignificanceFirst = new TH2D("fh2dJetSignedImpParXYZSignificanceFirst","fh2dJetSignedImpParXYZSignificanceFirst;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
	//N=2
	fh2dJetSignedImpParXYSecond = new TH2D("fh2dJetSignedImpParXYSecond","fh2dJetSignedImpParXYSecond;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
	fh2dJetSignedImpParXYZSecond = new TH2D("fh2dJetSignedImpParXYZSecond","fh2dJetSignedImpParXYZSecond;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
	fh2dJetSignedImpParXYSignificanceSecond = new TH2D("fh2dJetSignedImpParXYSignificanceSecond","fh2dJetSignedImpParXYSignificanceSecond;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
	fh2dJetSignedImpParXYZSignificanceSecond = new TH2D("fh2dJetSignedImpParXYZSignificanceSecond","fh2dJetSignedImpParXYZSignificanceThird;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
	//N=3
	fh2dJetSignedImpParXYThird = new TH2D("fh2dJetSignedImpParXYThird","fh2dJetSignedImpParXYThird;pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
	fh2dJetSignedImpParXYZThird = new TH2D("fh2dJetSignedImpParXYZThird","fh2dJetSignedImpParXYZThird;pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
	fh2dJetSignedImpParXYSignificanceThird = new TH2D("fh2dJetSignedImpParXYSignificanceThird","fh2dJetSignedImpParXYSignificanceThird;pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
	fh2dJetSignedImpParXYZSignificanceThird = new TH2D("fh2dJetSignedImpParXYZSignificanceThird","fh2dJetSignedImpParXYZSignificanceThird;pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

	if(fIsPythia){
		fh2dJetSignedImpParXY_McCorr = new TH2D("fh2dJetSignedImpParXY_McCorr","fh2dJetSignedImpParXY (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZ_McCorr = new TH2D("fh2dJetSignedImpParXYZ_McCorr","fh2dJetSignedImpParXYZ (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificance_McCorr = new TH2D("fh2dJetSignedImpParXYSignificance_McCorr","fh2dJetSignedImpParXYSignificance (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificance_McCorr = new TH2D("fh2dJetSignedImpParXYZSignificance_McCorr","fh2dJetSignedImpParXYZSignificance (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
		//N=1
		fh2dJetSignedImpParXYFirst_McCorr = new TH2D("fh2dJetSignedImpParXYFirst_McCorr","fh2dJetSignedImpParXYFirst (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZFirst_McCorr = new TH2D("fh2dJetSignedImpParXYZFirst_McCorr","fh2dJetSignedImpParXYZFirst (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceFirst_McCorr = new TH2D("fh2dJetSignedImpParXYSignificanceFirst_McCorr","fh2dJetSignedImpParXYSignificanceFirst (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceFirst_McCorr = new TH2D("fh2dJetSignedImpParXYZSignificanceFirst_McCorr","fh2dJetSignedImpParXYZSignificanceFirst (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
		//N=2
		fh2dJetSignedImpParXYSecond_McCorr = new TH2D("fh2dJetSignedImpParXYSecond_McCorr","fh2dJetSignedImpParXYSecond (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZSecond_McCorr = new TH2D("fh2dJetSignedImpParXYZSecond_McCorr","fh2dJetSignedImpParXYZSecond (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceSecond_McCorr = new TH2D("fh2dJetSignedImpParXYSignificanceSecond_McCorr","fh2dJetSignedImpParXYSignificanceSecond (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceSecond_McCorr = new TH2D("fh2dJetSignedImpParXYZSignificanceSecond_McCorr","fh2dJetSignedImpParXYZSignificanceThird (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
		//N=3
		fh2dJetSignedImpParXYThird_McCorr = new TH2D("fh2dJetSignedImpParXYThird_McCorr","fh2dJetSignedImpParXYThird (after correction);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZThird_McCorr = new TH2D("fh2dJetSignedImpParXYZThird_McCorr","fh2dJetSignedImpParXYZThird (after correction);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceThird_McCorr = new TH2D("fh2dJetSignedImpParXYSignificanceThird_McCorr","fh2dJetSignedImpParXYSignificanceThird (after correction);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceThird_McCorr = new TH2D("fh2dJetSignedImpParXYZSignificanceThird_McCorr","fh2dJetSignedImpParXYZSignificanceThird (after correction);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		// electron
		fh2dJetSignedImpParXY_electron = new TH2D("fh2dJetSignedImpParXY_electron","fh2dJetSignedImpParXY (electron);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZ_electron = new TH2D("fh2dJetSignedImpParXYZ_electron","fh2dJetSignedImpParXYZ (electron);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificance_electron = new TH2D("fh2dJetSignedImpParXYSignificance_electron","fh2dJetSignedImpParXYSignificance (electrons);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificance_electron = new TH2D("fh2dJetSignedImpParXYZSignificance_electron","fh2dJetSignedImpParXYZSignificance (electrons);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

		//N=1
		fh2dJetSignedImpParXYFirst_electron = new TH2D("fh2dJetSignedImpParXYFirst_electron","fh2dJetSignedImpParXYFirst (electrons);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZFirst_electron = new TH2D("fh2dJetSignedImpParXYZFirst_electron","fh2dJetSignedImpParXYZFirst (electrons);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceFirst_electron = new TH2D("fh2dJetSignedImpParXYSignificanceFirst_electron","fh2dJetSignedImpParXYSignificanceFirst (electrons);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceFirst_electron = new TH2D("fh2dJetSignedImpParXYZSignificanceFirst_electron","fh2dJetSignedImpParXYZSignificanceFirst (electrons);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
		//N=2
		fh2dJetSignedImpParXYSecond_electron = new TH2D("fh2dJetSignedImpParXYSecond_electron","fh2dJetSignedImpParXYSecond (electrons);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZSecond_electron = new TH2D("fh2dJetSignedImpParXYZSecond_electron","fh2dJetSignedImpParXYZSecond (electrons);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceSecond_electron = new TH2D("fh2dJetSignedImpParXYSignificanceSecond_electron","fh2dJetSignedImpParXYSignificanceSecond (electrons);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceSecond_electron = new TH2D("fh2dJetSignedImpParXYZSignificanceSecond_electron","fh2dJetSignedImpParXYZSignificanceThird (electrons);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);
		//N=3
		fh2dJetSignedImpParXYThird_electron = new TH2D("fh2dJetSignedImpParXYThird_electron","fh2dJetSignedImpParXYThird (electrons);pt_jet (GeV/c); radial imp. parameter (cm);a.u.",500,0.,250,nBins2d,-1,1);
		fh2dJetSignedImpParXYZThird_electron = new TH2D("fh2dJetSignedImpParXYZThird_electron","fh2dJetSignedImpParXYZThird (electrons);pt_jet (GeV/c); 3d imp. parameter (cm);a.u.",500,0.,250,nBins3d,-1,1);
		fh2dJetSignedImpParXYSignificanceThird_electron = new TH2D("fh2dJetSignedImpParXYSignificanceThird_electron","fh2dJetSignedImpParXYSignificanceThird (electrons);pt_jet (GeV/c); radial imp. parameter significance;a.u.",500,0.,250,nBins2dSignificance,-100,100);
		fh2dJetSignedImpParXYZSignificanceThird_electron = new TH2D("fh2dJetSignedImpParXYZSignificanceThird_electron","fh2dJetSignedImpParXYZSignificanceThird (electrons);pt_jet (GeV/c); 3d imp. parameter significance;a.u.",500,0.,250,nBins3dSignificance,-100,100);

	}

	//Add to output list
	fOutput->Add(fh1dEventRejectionRDHFCuts);
	fOutput->Add(fh1dVertexZ);
	fOutput->Add(fh1dVertexZAccepted);
	fOutput->Add(fh1dVertexR);
	fOutput->Add(fh1dVertexRAccepted);
	fOutput->Add(fh1dTracksAccepeted);
	fOutput->Add(fh1dTracksImpParXY);
	fOutput->Add(fh1dTracksImpParXYZ);
	fOutput->Add(fh1dTracksImpParXYSignificance);
	fOutput->Add(fh1dTracksImpParXYZSignificance);
	fOutput->Add(fh2dVertexChi2NDFNESDTracks);
	if(fIsPythia){

		fOutput->Add(fh1dTracksImpParXYTruth);
		fOutput->Add(fh1dTracksImpParXYZTruth);
		fOutput->Add(fh1dTracksImpParXYResidualTruth);
		fOutput->Add(fh1dTracksImpParXYZResidualTruth);
		fOutput->Add(fh1dJetGenPt);
		fOutput->Add(fh1dJetGenPtUnidentified);
		fOutput->Add(fh1dJetGenPtudsg);
		fOutput->Add(fh1dJetGenPtc);
		fOutput->Add(fh1dJetGenPtb);
		fOutput->Add(fh2dJetGenPtVsJetRecPt);


		fOutput->Add(fh2dJetSignedImpParXYUnidentified);
		fOutput->Add(fh2dJetSignedImpParXYZUnidentified);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentified);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentified);
		fOutput->Add(fh2dJetSignedImpParXYudsg);
		fOutput->Add(fh2dJetSignedImpParXYZudsg);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceudsg);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsg);
		fOutput->Add(fh2dJetSignedImpParXYc);
		fOutput->Add(fh2dJetSignedImpParXYZc);
		fOutput->Add(fh2dJetSignedImpParXYSignificancec);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancec);
		fOutput->Add(fh2dJetSignedImpParXYb);
		fOutput->Add(fh2dJetSignedImpParXYZb);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceb);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceb);



		fOutput->Add(fh2dJetSignedImpParXYUnidentifiedFirst);
		fOutput->Add(fh2dJetSignedImpParXYZUnidentifiedFirst);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentifiedFirst);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst);
		fOutput->Add(fh2dJetSignedImpParXYudsgFirst);
		fOutput->Add(fh2dJetSignedImpParXYZudsgFirst);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceudsgFirst);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsgFirst);
		fOutput->Add(fh2dJetSignedImpParXYcFirst);
		fOutput->Add(fh2dJetSignedImpParXYZcFirst);
		fOutput->Add(fh2dJetSignedImpParXYSignificancecFirst);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancecFirst);
		fOutput->Add(fh2dJetSignedImpParXYbFirst);
		fOutput->Add(fh2dJetSignedImpParXYZbFirst);
		fOutput->Add(fh2dJetSignedImpParXYSignificancebFirst);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancebFirst);

		fOutput->Add(fh2dJetSignedImpParXYUnidentifiedSecond);
		fOutput->Add(fh2dJetSignedImpParXYZUnidentifiedSecond);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentifiedSecond);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond);
		fOutput->Add(fh2dJetSignedImpParXYudsgSecond);
		fOutput->Add(fh2dJetSignedImpParXYZudsgSecond);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceudsgSecond);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsgSecond);
		fOutput->Add(fh2dJetSignedImpParXYcSecond);
		fOutput->Add(fh2dJetSignedImpParXYZcSecond);
		fOutput->Add(fh2dJetSignedImpParXYSignificancecSecond);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancecSecond);
		fOutput->Add(fh2dJetSignedImpParXYbSecond);
		fOutput->Add(fh2dJetSignedImpParXYZbSecond);
		fOutput->Add(fh2dJetSignedImpParXYSignificancebSecond);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancebSecond);

		fOutput->Add(fh2dJetSignedImpParXYUnidentifiedThird);
		fOutput->Add(fh2dJetSignedImpParXYZUnidentifiedThird);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentifiedThird);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentifiedThird);
		fOutput->Add(fh2dJetSignedImpParXYudsgThird);
		fOutput->Add(fh2dJetSignedImpParXYZudsgThird);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceudsgThird);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsgThird);
		fOutput->Add(fh2dJetSignedImpParXYcThird);
		fOutput->Add(fh2dJetSignedImpParXYZcThird);
		fOutput->Add(fh2dJetSignedImpParXYSignificancecThird);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancecThird);
		fOutput->Add(fh2dJetSignedImpParXYbThird);
		fOutput->Add(fh2dJetSignedImpParXYZbThird);
		fOutput->Add(fh2dJetSignedImpParXYSignificancebThird);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancebThird);

		// Add corrected Monte Carlo plots

		fOutput->Add(fh1dTracksImpParXY_McCorr);
		fOutput->Add(fh1dTracksImpParXYZ_McCorr);
		fOutput->Add(fh1dTracksImpParXYSignificance_McCorr);
		fOutput->Add(fh1dTracksImpParXYZSignificance_McCorr);

		fOutput->Add(fh2dJetSignedImpParXYUnidentified_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZUnidentified_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentified_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentified_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYudsg_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZudsg_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceudsg_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsg_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYc_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZc_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificancec_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancec_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYb_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZb_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceb_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceb_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYb_McCorrNonBDecay);
		fOutput->Add(fh2dJetSignedImpParXYZb_McCorr_McCorrNonBDecay);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceb_McCorrNonBDecay);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceb_McCorrNonBDecay);

		fOutput->Add(fh2dJetSignedImpParXYUnidentifiedFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZUnidentifiedFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentifiedFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYudsgFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZudsgFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceudsgFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsgFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYcFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZcFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificancecFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancecFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYbFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZbFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificancebFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancebFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYbFirst_McCorrNonBDecay);
		fOutput->Add(fh2dJetSignedImpParXYZbFirst_McCorrNonBDecay);
		fOutput->Add(fh2dJetSignedImpParXYSignificancebFirst_McCorrNonBDecay);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancebFirst_McCorrNonBDecay);


		fOutput->Add(fh2dJetSignedImpParXYUnidentifiedSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZUnidentifiedSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentifiedSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYudsgSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZudsgSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceudsgSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsgSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYcSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZcSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificancecSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancecSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYbSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZbSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificancebSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancebSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYbSecond_McCorrNonBDecay);
		fOutput->Add(fh2dJetSignedImpParXYZbSecond_McCorrNonBDecay);
		fOutput->Add(fh2dJetSignedImpParXYSignificancebSecond_McCorrNonBDecay);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancebSecond_McCorrNonBDecay);

		fOutput->Add(fh2dJetSignedImpParXYUnidentifiedThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZUnidentifiedThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentifiedThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentifiedThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYudsgThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZudsgThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceudsgThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsgThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYcThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZcThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificancecThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancecThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYbThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZbThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificancebThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancebThird_McCorr);

		fOutput->Add(fh2dJetSignedImpParXYbThird_McCorrNonBDecay);
		fOutput->Add(fh2dJetSignedImpParXYZbThird_McCorrNonBDecay);
		fOutput->Add(fh2dJetSignedImpParXYSignificancebThird_McCorrNonBDecay);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancebThird_McCorrNonBDecay);

		fOutput->Add(fh2dJetSignedImpParXY_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZ_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificance_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificance_McCorr);

		fOutput->Add(fh2dJetSignedImpParXYFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceFirst_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceFirst_McCorr);

		fOutput->Add(fh2dJetSignedImpParXYSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceSecond_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceSecond_McCorr);

		fOutput->Add(fh2dJetSignedImpParXYThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYZThird_McCorr);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceThird_McCorr);
		//electrons---------

		fOutput->Add(fh1dTracksImpParXY_electron);


		fOutput->Add(fh1dTracksImpParXYZ_electron);
		fOutput->Add(fh1dTracksImpParXYSignificance_electron);
		fOutput->Add(fh1dTracksImpParXYZSignificance_electron);

		fOutput->Add(fh2dJetSignedImpParXYUnidentified_electron);
		fOutput->Add(fh2dJetSignedImpParXYZUnidentified_electron);

		fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentified_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentified_electron);

		fOutput->Add(fh2dJetSignedImpParXYudsg_electron);
		fOutput->Add(fh2dJetSignedImpParXYZudsg_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceudsg_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsg_electron);
		fOutput->Add(fh2dJetSignedImpParXYc_electron);
		fOutput->Add(fh2dJetSignedImpParXYZc_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificancec_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancec_electron);
		fOutput->Add(fh2dJetSignedImpParXYb_electron);
		fOutput->Add(fh2dJetSignedImpParXYZb_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceb_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceb_electron);



		fOutput->Add(fh2dJetSignedImpParXYUnidentifiedFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYZUnidentifiedFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentifiedFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentifiedFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYudsgFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYZudsgFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceudsgFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsgFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYcFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYZcFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificancecFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancecFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYbFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYZbFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificancebFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancebFirst_electron);

		fOutput->Add(fh2dJetSignedImpParXYUnidentifiedSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYZUnidentifiedSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentifiedSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentifiedSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYudsgSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYZudsgSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceudsgSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsgSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYcSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYZcSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificancecSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancecSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYbSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYZbSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificancebSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancebSecond_electron);

		fOutput->Add(fh2dJetSignedImpParXYUnidentifiedThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYZUnidentifiedThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceUnidentifiedThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceUnidentifiedThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYudsgThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYZudsgThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceudsgThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceudsgThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYcThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYZcThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificancecThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancecThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYbThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYZbThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificancebThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificancebThird_electron);

		fOutput->Add(fh2dJetSignedImpParXY_electron);
		fOutput->Add(fh2dJetSignedImpParXYZ_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificance_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificance_electron);
		fOutput->Add(fh2dJetSignedImpParXYFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYZFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceFirst_electron);
		fOutput->Add(fh2dJetSignedImpParXYSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYZSignificanceSecond_electron);
		fOutput->Add(fh2dJetSignedImpParXYThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYZThird_electron);
		fOutput->Add(fh2dJetSignedImpParXYSignificanceThird_electron);


	}
	fOutput->Add(fh1dJetRecPt);
	fOutput->Add(fh1dJetRecPtAccepted);
	fOutput->Add(fh1dJetRecEtaPhiAccepted);

	if(fIsPythia){

		fOutput->Add(fh1dJetRecPtUnidentified);
		fOutput->Add(fh1dJetRecPtudsg);
		fOutput->Add(fh1dJetRecPtc);
		fOutput->Add(fh1dJetRecPtb);
		fOutput->Add(fh1dJetRecPtUnidentifiedAccepted);
		fOutput->Add(fh1dJetRecPtudsgAccepted);
		fOutput->Add(fh1dJetRecPtcAccepted);
		fOutput->Add(fh1dJetRecPtbAccepted);
	}

	fOutput->Add(fh2dJetSignedImpParXY);
	fOutput->Add(fh2dJetSignedImpParXYZ);
	fOutput->Add(fh2dJetSignedImpParXYSignificance);
	fOutput->Add(fh2dJetSignedImpParXYZSignificance);
	fOutput->Add(fh2dJetSignedImpParXYFirst);
	fOutput->Add(fh2dJetSignedImpParXYZFirst);
	fOutput->Add(fh2dJetSignedImpParXYSignificanceFirst);
	fOutput->Add(fh2dJetSignedImpParXYZSignificanceFirst);
	fOutput->Add(fh2dJetSignedImpParXYSecond);
	fOutput->Add(fh2dJetSignedImpParXYZSecond);
	fOutput->Add(fh2dJetSignedImpParXYSignificanceSecond);
	fOutput->Add(fh2dJetSignedImpParXYZSignificanceSecond);
	fOutput->Add(fh2dJetSignedImpParXYThird);
	fOutput->Add(fh2dJetSignedImpParXYZThird);
	fOutput->Add(fh2dJetSignedImpParXYSignificanceThird);
	fOutput->Add(fh2dJetSignedImpParXYZSignificanceThird);

	TIter next(fOutput);
	while (TObject *obj = next.Next()){
		if(obj->IsA() == TH1D::Class() || obj->IsA() == TH2D::Class()){
			((TH1*)obj)->Sumw2();
		}
	}
	PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}
// ######################################################################################## Calculate impact parameters
Bool_t AliAnalysisTaskHFJetIPQA::CalculateTrackImpactParameter(AliAODTrack * track,double *impar, double * cov)
{
	AliAODEvent* aev = dynamic_cast<AliAODEvent*>(InputEvent());
	AliAODVertex *vtxAOD = aev->GetPrimaryVertex();
	AliAODVertex *vtxAODNew=vtxAOD;
	if(!vtxAOD) return kFALSE;
	TString title=vtxAOD->GetTitle();
	if(!title.Contains("VertexerTracks")) return kFALSE;
	AliESDVertex *vtxESDNew =0x0;
	Bool_t recalculate = kFALSE;
	if( vtxAOD->GetNContributors() < 30){
		recalculate=kTRUE;
		AliVertexerTracks *vertexer = new AliVertexerTracks(aev->GetMagneticField());
		Int_t ndg = 1;
		vertexer->SetITSMode();
		vertexer->SetMinClusters(3);
		vertexer->SetConstraintOff();
		if(title.Contains("WithConstraint")) {
			Float_t diamondcovxy[3];
			aev->GetDiamondCovXY(diamondcovxy);
			Double_t pos[3]={aev->GetDiamondX(),aev->GetDiamondY(),0.};
			Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
			AliESDVertex *diamond = new AliESDVertex(pos,cov,1.,1);
			vertexer->SetVtxStart(diamond);
			delete diamond; diamond=NULL;
		}
		Int_t skipped[1] = {-1};
		Int_t id = (Int_t)track->GetID();
		if(id<0) return kFALSE;
		skipped[0] = id;
		vertexer->SetSkipTracks(1,skipped);
		vtxESDNew = vertexer->FindPrimaryVertex(aev);
		delete vertexer; vertexer=NULL;
		if(!vtxESDNew) return kFALSE;
		if(vtxESDNew->GetNContributors()<=0) {
			delete vtxESDNew; vtxESDNew=NULL;
			return kFALSE;
		}
		// convert to AliAODVertex
		Double_t pos[3],cova[6],chi2perNDF;
		vtxESDNew->GetXYZ(pos); // position
		vtxESDNew->GetCovMatrix(cova); //covariance matrix
		chi2perNDF = vtxESDNew->GetChi2toNDF();
		delete vtxESDNew; vtxESDNew=NULL;
		vtxAODNew = new AliAODVertex(pos,cova,chi2perNDF);
	}
	// Calculate Impact Parameters
	AliExternalTrackParam etp; etp.CopyFromVTrack(track);
	if(etp.PropagateToDCA(vtxAODNew,aev->GetMagneticField(),3.,impar,cov))
	{
		if(recalculate)
			delete vtxAODNew;
		return kTRUE;
	}
	else{
		if(recalculate)
			delete vtxAODNew;
		return kFALSE;

	}
}
// ######################################################################################## Calculate impact parameters
Bool_t AliAnalysisTaskHFJetIPQA::CalculateTrackImpactParameter(AliESDtrack * track,double *impar, double * cov)
{
	//
	// Copied from AliHFEextraCuts::GetHFEImpactParameters  impact parameter (with recalculated primary vertex)
	//
	////		Printf("%s:%i",__FUNCTION__,__LINE__);


	const Double_t kBeampiperadius=3.;
	AliVVertex *vtxESDSkip = NULL;
	AliESDEvent * eev = (AliESDEvent*)InputEvent();
	////		Printf("%s:%i",__FUNCTION__,__LINE__);


	if(!eev){
		AliDebug(1, "No Input event available\n");
		return kFALSE;
	}
	////		Printf("%s:%i",__FUNCTION__,__LINE__);


	TString type = track->IsA()->GetName();
	Double_t dcaD[2]={-999.,-999.},
			covD[3]={-999.,-999.,-999.};
	Bool_t isRecalcVertex(kFALSE);
	////		Printf("%s:%i",__FUNCTION__,__LINE__);


	if(!TString(track->IsA()->GetName()).CompareTo("AliESDtrack")){
		//case of ESD tracks
		////		Printf("%s:%i",__FUNCTION__,__LINE__);


		AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(InputEvent());
		if(!esdevent) {
			AliDebug(1, "No esd event available\n");
			return kFALSE;
		}
		AliVVertex *vtxESDSkip = (	 AliVVertex *)esdevent->GetPrimaryVertex();
		if(!vtxESDSkip) return kFALSE;
		//case ESD track: take copy constructor

		//		Printf("%s:%i",__FUNCTION__,__LINE__);

		const AliESDtrack *tmptrack = dynamic_cast<const AliESDtrack *>(track);
		if(tmptrack){
			//		Printf("%s:%i",__FUNCTION__,__LINE__);

			if( vtxESDSkip->GetNContributors() < 30){ // if vertex contributor is smaller than 30, recalculate the primary vertex

				AliVertexerTracks vertexer(eev->GetMagneticField());
				vertexer.SetITSMode();
				vertexer.SetMinClusters(4);
				Int_t skipped[2];
				skipped[0] = track->GetID();
				vertexer.SetSkipTracks(1,skipped);
				//diamond constraint
				vertexer.SetConstraintOn();
				Float_t diamondcovxy[3];
				esdevent->GetDiamondCovXY(diamondcovxy);
				Double_t pos[3]={eev->GetDiamondX(),eev->GetDiamondY(),0.};
				Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
				AliESDVertex diamond(pos,cov,1.,1);
				vertexer.SetVtxStart(&diamond);
				vtxESDSkip = vertexer.FindPrimaryVertex(eev);
				isRecalcVertex = kTRUE;
			}
			//		Printf("%s:%i",__FUNCTION__,__LINE__);

			if(vtxESDSkip){
				AliESDtrack esdtrack(*tmptrack);
				if(esdtrack.PropagateToDCA(vtxESDSkip, eev->GetMagneticField(), kBeampiperadius, impar, cov)){
					if(isRecalcVertex) delete vtxESDSkip;
					return kTRUE;
				}
				else{
					delete vtxESDSkip;
				}
				if(isRecalcVertex) delete vtxESDSkip;
				return kFALSE;
			}
		}
	}
	//		Printf("%s:%i",__FUNCTION__,__LINE__);

	return kFALSE;
}

// ######################################################################################## Calculate impact parameters based on MC event vertex and MC particle information (no special mass treatment)
Bool_t AliAnalysisTaskHFJetIPQA::CalculateTrackImpactParameterTruth(AliAODTrack * track,double *impar, double * cov)
{
	AliAODMCParticle *pMC = 0x0;
	AliAODEvent* aev = dynamic_cast<AliAODEvent*>(InputEvent());
	AliAODMCHeader* mcheader = dynamic_cast<AliAODMCHeader*>(aev->FindListObject(AliAODMCHeader::StdBranchName()));
	if (!mcheader) return kFALSE;

	TClonesArray * fMCparticles = dynamic_cast<TClonesArray*>(aev->FindListObject(AliAODMCParticle::StdBranchName()));
	if (!fMCparticles) return kFALSE;
	if(track->GetLabel()>-1)
		pMC = dynamic_cast<AliAODMCParticle*>(fMCparticles->At(track->GetLabel()));
	if (!pMC) return kFALSE;

	Double_t pos[3]={0,0,0};
	mcheader->GetVertex(pos);

	Double_t cova[6]={0,0,0,0,0,0};
	Double_t chi2perNDF =0;
	AliAODVertex *vtxAODNew = new AliAODVertex(pos,cova,chi2perNDF);
	double xpart[3] = {0,0,0};
	pMC->XvYvZv(xpart);
	double ppart[3] = {0,0,0};
	pMC->PxPyPz(ppart);
	double cv[21] ;
	for (int i=0;i<21;++i)cv[i] =0.;
	AliExternalTrackParam trackparam(xpart,ppart,cv,(TMath::Sign((Short_t)1,(Short_t)pMC->Charge())));
	if(trackparam.PropagateToDCA(vtxAODNew,aev->GetMagneticField(),3.,impar,cov))
	{
		delete vtxAODNew;
		return kTRUE;
	}
	else{
		delete vtxAODNew;
		return kFALSE;

	}

	return kFALSE;
}
// ######################################################################################## Calculate impact parameters based on MC event vertex and MC particle information (no special mass treatment)
Bool_t AliAnalysisTaskHFJetIPQA::CalculateTrackImpactParameterTruth(AliESDtrack * track,double *impar, double * cov)
{
	AliMCParticle *pMC = 0x0;
	AliESDEvent* eev = dynamic_cast<AliESDEvent*>(InputEvent());
	////		Printf("%s:%i",__FUNCTION__,__LINE__);



	AliMCEvent * MCparticles = dynamic_cast<AliMCEvent*>(MCEvent());
	if (!MCparticles) return kFALSE;
	pMC = dynamic_cast<AliMCParticle*>(MCEvent()->GetTrack(abs(track->GetLabel())));
	if (!pMC) return kFALSE;
	////		Printf("%s:%i",__FUNCTION__,__LINE__);


	Double_t pos[3]={0,0,0};
	MCparticles->GetPrimaryVertex()->GetXYZ(pos);
	////		Printf("%s:%i",__FUNCTION__,__LINE__);


	Double_t cova[6]={0,0,0,0,0,0};
	Double_t chi2perNDF =0;
	const AliVVertex *vtxAODNew =  MCparticles->GetPrimaryVertex();
	double xpart[3] = {0,0,0};
	pMC->XvYvZv(xpart);
	double ppart[3] = {0,0,0};
	pMC->PxPyPz(ppart);
	double cv[21] ;
	for (int i=0;i<21;++i)cv[i] =0.;
	AliExternalTrackParam trackparam(xpart,ppart,cv,(TMath::Sign((Short_t)1,(Short_t)pMC->Charge())));
	if(trackparam.PropagateToDCA(vtxAODNew,eev->GetMagneticField(),3.,impar,cov))
	{
		//delete vtxAODNew;
		return kTRUE;
	}
	else{
		//delete vtxAODNew;
		return kFALSE;

	}
	////		Printf("%s:%i",__FUNCTION__,__LINE__);


	return kFALSE;
}

// ######################################################################################## Calculate signed  impact parameters

Bool_t AliAnalysisTaskHFJetIPQA::CalculateJetSignedTrackImpactParameter(AliAODTrack * track,AliEmcalJet * jet ,double *impar, double * cov, double &sign, double &dcajetrack, double &lineardecaylength){
	AliAODEvent* aev = dynamic_cast<AliAODEvent*>(InputEvent());
	AliAODVertex *vtxAOD = aev->GetPrimaryVertex();
	if(!vtxAOD) return kFALSE;
	TString title=vtxAOD->GetTitle();
	if(!title.Contains("VertexerTracks")) return kFALSE;
	AliVertexerTracks *vertexer = new AliVertexerTracks(aev->GetMagneticField());
	Int_t ndg = 1;
	vertexer->SetITSMode();
	vertexer->SetMinClusters(3);
	vertexer->SetConstraintOff();
	if(title.Contains("WithConstraint")) {
		Float_t diamondcovxy[3];
		aev->GetDiamondCovXY(diamondcovxy);
		Double_t pos[3]={aev->GetDiamondX(),aev->GetDiamondY(),0.};
		Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
		AliESDVertex *diamond = new AliESDVertex(pos,cov,1.,1);
		vertexer->SetVtxStart(diamond);
		delete diamond; diamond=NULL;
	}
	Int_t skipped[1] = {-1};
	Int_t id = (Int_t)track->GetID();
	if(id<0) return kFALSE;
	skipped[0] = id;
	vertexer->SetSkipTracks(1,skipped);
	AliESDVertex *vtxESDNew = vertexer->FindPrimaryVertex(aev);
	delete vertexer; vertexer=NULL;
	if(!vtxESDNew) return kFALSE;
	if(vtxESDNew->GetNContributors()<=0) {
		delete vtxESDNew; vtxESDNew=NULL;
		return kFALSE;
	}
	// convert to AliAODVertex
	Double_t pos[3],cova[6],chi2perNDF;
	vtxESDNew->GetXYZ(pos); // position
	vtxESDNew->GetCovMatrix(cova); //covariance matrix
	chi2perNDF = vtxESDNew->GetChi2toNDF();
	delete vtxESDNew; vtxESDNew=NULL;
	AliAODVertex *vtxAODNew = new AliAODVertex(pos,cova,chi2perNDF);

	// Calculate Impact Parameters
	AliExternalTrackParam etp; etp.CopyFromVTrack(track);
	if(etp.PropagateToDCA(vtxAODNew,aev->GetMagneticField(),3.,impar,cov))
	{
		//Calculate Sign
		Double_t posdcatrack[3]= {0.,0.,0.};
		etp.GetXYZ(posdcatrack);
		Double_t ipvector3[3] = { posdcatrack[0] - pos[0], posdcatrack[1] - pos[1], posdcatrack[2] - pos[2] };
		sign =TMath::Sign(1.,ipvector3[0]*jet->Px() +ipvector3[1]*jet->Py()+ipvector3[2]*jet->Pz() );
		// Calculate decay legnth and track jet DCA against new vertex
		Double_t bpos[3] = { 0,0,0 };
		vtxAODNew->GetXYZ(bpos);
		Double_t bpxpypz[3] = { jet->Px(), jet->Py(), jet->Pz() };
		Double_t bcv[21] = { 0 };
		AliExternalTrackParam bjetparam(bpos, bpxpypz, bcv, (Short_t)0);
		Double_t xa = 0., xb = 0.;
		bjetparam.GetDCA(&etp, aev->GetMagneticField(), xa, xb);
		Double_t xyz[3] = { 0., 0., 0. };
		Double_t xyzb[3] = { 0., 0., 0. };
		bjetparam.GetXYZAt(xa, aev->GetMagneticField(), xyz);
		etp.GetXYZAt(xb, aev->GetMagneticField(), xyzb);
		double  bdecaylength =
				TMath::Sqrt(
						(bpos[0] - xyz[0]) * (bpos[0] - xyz[0]) +
						(bpos[1] - xyz[1]) * (bpos[1] - xyz[1]) +
						(bpos[2] - xyz[2]) * (bpos[2] - xyz[2]));
		dcajetrack =
				TMath::Sqrt(
						(xyzb[0] - xyz[0]) * (xyzb[0] - xyz[0]) +
						(xyzb[1] - xyz[1]) * (xyzb[1] - xyz[1]) +
						(xyzb[2] - xyz[2]) * (xyzb[2]- xyz[2]));
		if(bdecaylength>0) lineardecaylength=bdecaylength;
		delete vtxAODNew;
		return kTRUE;
	}
	else{
		delete vtxAODNew;
		return kFALSE;

	}
}

// ######################################################################################## Calculate signed  impact parameters

Bool_t AliAnalysisTaskHFJetIPQA::CalculateJetSignedTrackImpactParameter(AliESDtrack * track,AliEmcalJet * jet ,double *impar, double * cov, double &sign, double &dcajetrack, double &lineardecaylength){
	AliESDEvent* eev = dynamic_cast<AliESDEvent*>(InputEvent());
	//
	// Copied from AliHFEextraCuts::GetHFEImpactParameters  impact parameter (with recalculated primary vertex)
	//
	AliVVertex *vtxESDSkip = NULL;

	if(!eev){
		AliDebug(1, "No Input event available\n");
		return kFALSE;
	}
	TString type = track->IsA()->GetName();
	const Double_t kBeampiperadius=3.;
	Double_t dcaD[2]={-999.,-999.},
			covD[3]={-999.,-999.,-999.};

	Double_t pos[3] = {0.,0.,0.};
	Bool_t isRecalcVertex(kFALSE);

	if(!TString(track->IsA()->GetName()).CompareTo("AliESDtrack")){
		//case of ESD tracks
		AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(eev);
		if(!esdevent) {
			AliDebug(1, "No esd event available\n");
			return kFALSE;
		}
		const AliVVertex *vtxESDSkip = esdevent->GetPrimaryVertex();
		if(!vtxESDSkip) return kFALSE;
		//case ESD track: take copy constructor
		const AliESDtrack *tmptrack = dynamic_cast<const AliESDtrack *>(track);
		if(tmptrack){

			if( vtxESDSkip->GetNContributors() < 30){ // if vertex contributor is smaller than 30, recalculate the primary vertex

				AliVertexerTracks vertexer(eev->GetMagneticField());
				vertexer.SetITSMode();
				vertexer.SetMinClusters(4);
				Int_t skipped[2];
				skipped[0] = track->GetID();
				vertexer.SetSkipTracks(1,skipped);
				//diamond constraint
				vertexer.SetConstraintOn();
				Float_t diamondcovxy[3];
				esdevent->GetDiamondCovXY(diamondcovxy);
				Double_t pos[3]={eev->GetDiamondX(),eev->GetDiamondY(),0.};
				Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
				AliESDVertex diamond(pos,cov,1.,1);
				vertexer.SetVtxStart(&diamond);
				vtxESDSkip = vertexer.FindPrimaryVertex(eev);
				isRecalcVertex = kTRUE;
			}

			if(vtxESDSkip){
				AliESDtrack esdtrack(*tmptrack);
				if(esdtrack.PropagateToDCA(vtxESDSkip, eev->GetMagneticField(), kBeampiperadius, impar, cov)){
					//*****//
					//Calculate Sign
					Double_t posdcatrack[3]= {0.,0.,0.};
					esdtrack.GetXYZ(posdcatrack);
					vtxESDSkip->GetXYZ(pos);
					Double_t ipvector3[3] = { posdcatrack[0] - pos[0], posdcatrack[1] - pos[1], posdcatrack[2] - pos[2] };
					sign =TMath::Sign(1.,ipvector3[0]*jet->Px() +ipvector3[1]*jet->Py()+ipvector3[2]*jet->Pz() );
					// Calculate decay legnth and track jet DCA against new vertex
					Double_t bpos[3] = { 0,0,0 };
					vtxESDSkip->GetXYZ(bpos);
					Double_t bpxpypz[3] = { jet->Px(), jet->Py(), jet->Pz() };
					Double_t bcv[21] = { 0 };
					AliExternalTrackParam bjetparam(bpos, bpxpypz, bcv, (Short_t)0);
					Double_t xa = 0., xb = 0.;
					bjetparam.GetDCA(&esdtrack, eev->GetMagneticField(), xa, xb);
					Double_t xyz[3] = { 0., 0., 0. };
					Double_t xyzb[3] = { 0., 0., 0. };
					bjetparam.GetXYZAt(xa, eev->GetMagneticField(), xyz);
					esdtrack.GetXYZAt(xb, eev->GetMagneticField(), xyzb);
					double  bdecaylength =
							TMath::Sqrt(
									(bpos[0] - xyz[0]) * (bpos[0] - xyz[0]) +
									(bpos[1] - xyz[1]) * (bpos[1] - xyz[1]) +
									(bpos[2] - xyz[2]) * (bpos[2] - xyz[2]));
					dcajetrack =
							TMath::Sqrt(
									(xyzb[0] - xyz[0]) * (xyzb[0] - xyz[0]) +
									(xyzb[1] - xyz[1]) * (xyzb[1] - xyz[1]) +
									(xyzb[2] - xyz[2]) * (xyzb[2]- xyz[2]));
					if(bdecaylength>0) lineardecaylength=bdecaylength;
					if(isRecalcVertex) delete vtxESDSkip;
					return kTRUE;
				}
				else{
					delete vtxESDSkip;
				}
				if(isRecalcVertex) delete vtxESDSkip;
				return kFALSE;
			}
		}
	}
	return kFALSE;
}
// ######################################################################################## Post-process ImpPar
Double_t AliAnalysisTaskHFJetIPQA::GetValImpactParameter(TTypeImpPar type,double *impar, double * cov)
{
	double result =-999999;
	double dFdx = 0;
	double dFdy = 0;

	switch(type){
	case kXY:
		result = impar[0];
		break;
	case kXYSig:
		result = impar[0]/TMath::Sqrt(cov[0]);
		break;
	case kXYZ:
		result = TMath::Sqrt(impar[0]*impar[0]+impar[1]*impar[1]);
		break;
	case kXYZSig:
		result =  TMath::Sqrt(impar[0]*impar[0]+impar[1]*impar[1]);
		dFdx = 2*impar[0]/result;
		dFdy = 2*impar[1]/result;
		result /=TMath::Sqrt(cov[0]*dFdx*dFdx + cov[2]*dFdy*dFdy + 2* cov[1] *dFdx*dFdy);
		break;
	case kXYZSigmaOnly:
		dFdx = 2*impar[0]/result;
		dFdy = 2*impar[1]/result;
		result =TMath::Sqrt(cov[0]*dFdx*dFdx + cov[2]*dFdy*dFdy + 2* cov[1] *dFdx*dFdy);
		break;

	default:
		break;
	}
	return result;
}
// ########################################################################################Track Selection
Bool_t AliAnalysisTaskHFJetIPQA::IsTrackAccepted(AliESDtrack* track){
	if(!fESDTrackCut->AcceptTrack(track)) return kFALSE;
	if(track->Pt() < 1.)return kFALSE;
	ULong_t status = track->GetStatus();
	if(!(status & AliESDtrack::kTPCrefit))return kFALSE;
	if(!(status & AliESDtrack::kITSrefit))return kFALSE;
	if(!(track->HasPointOnITSLayer(0)&& track->HasPointOnITSLayer(1))) return kFALSE;
	Float_t cRatioTPC = track->GetTPCNclsF() > 0. ?
			static_cast<Float_t>(track->GetTPCNcls()) / static_cast<Float_t>(track->GetTPCNclsF()) :
			1.;
	if(cRatioTPC < 0.6) return kFALSE;
	if(track->GetNcls(1)<100) return kFALSE;
	if(track->GetNcls(0)<3) return kFALSE;
	if(track->GetTPCchi2()/track->GetNcls(1)>=4) return kFALSE;
	if(fEnableV0GammaRejection)if(IsV0PhotonFromBeamPipeDaughter(track)) return kFALSE;
	return kTRUE;
}
// ########################################################################################Track Selection
Bool_t AliAnalysisTaskHFJetIPQA::IsTrackAccepted(AliAODTrack* track){
	if(!((track)->TestFilterBit(1 << 4))) return kFALSE;
	if(track->Pt() < 1.)return kFALSE;
	ULong_t status = track->GetStatus();
	if(!(status & AliAODTrack::kTPCrefit))return kFALSE;
	if(!(status & AliAODTrack::kITSrefit))return kFALSE;
	if(!(track->HasPointOnITSLayer(0)&& track->HasPointOnITSLayer(1))) return kFALSE;
	Float_t cRatioTPC = track->GetTPCNclsF() > 0. ?
			static_cast<Float_t>(track->GetTPCNcls()) / static_cast<Float_t>(track->GetTPCNclsF()) :
			1.;
	if(cRatioTPC < 0.6) return kFALSE;
	if(track->GetNcls(1)<100) return kFALSE;
	if(track->GetNcls(0)<3) return kFALSE;
	if(track->GetTPCchi2()/track->GetNcls(1)>=4) return kFALSE;
	if(fEnableV0GammaRejection)if(IsV0PhotonFromBeamPipeDaughter(track)) return kFALSE;
	return kTRUE;
}
// ######################################################################################## Jet matching 1/4
Bool_t AliAnalysisTaskHFJetIPQA::MatchJetsGeometricDefault()
{
	AliJetContainer *jets1 = static_cast<AliJetContainer*>(fJetCollArray.At(0));
	AliJetContainer *jets2 = static_cast<AliJetContainer*>(fJetCollArray.At(1));
	double matchingpar1 =0.25;
	double matchingpar2 =0.25;
	if (!jets1 || !jets1->GetArray() || !jets2 || !jets2->GetArray()) return kFALSE;
	DoJetLoop();
	AliEmcalJet* jet1 = 0;
	jets1->ResetCurrentID();
	while ((jet1 = jets1->GetNextJet())) {
		AliEmcalJet *jet2 = jet1->ClosestJet();
		if (!jet2) continue;
		if (jet2->ClosestJet() != jet1) continue;
		if (jet1->ClosestJetDistance() > matchingpar1 || jet2->ClosestJetDistance() > matchingpar2) continue;
		// Matched jet found
		jet1->SetMatchedToClosest(1);
		jet2->SetMatchedToClosest(1);
	}
	return kTRUE;
}
// ######################################################################################## Jet matching 2/4
void AliAnalysisTaskHFJetIPQA::DoJetLoop()
{
	// Do the jet loop.
	double minjetpt =1.;
	AliJetContainer *jets1 = static_cast<AliJetContainer*>(fJetCollArray.At(0));
	AliJetContainer *jets2 = static_cast<AliJetContainer*>(fJetCollArray.At(1));
	if (!jets1 || !jets1->GetArray() || !jets2 || !jets2->GetArray()) return;
	AliEmcalJet* jet1 = 0;
	AliEmcalJet* jet2 = 0;
	jets2->ResetCurrentID();
	while ((jet2 = jets2->GetNextJet())) jet2->ResetMatching();
	jets1->ResetCurrentID();
	while ((jet1 = jets1->GetNextJet())) {
		jet1->ResetMatching();
		if (jet1->MCPt() < minjetpt) continue;
		jets2->ResetCurrentID();
		while ((jet2 = jets2->GetNextJet())) {
			SetMatchingLevel(jet1, jet2, 1);
		} // jet2 loop
	} // jet1 loop
}
// ######################################################################################## Jet matching 3/4
Double_t AliAnalysisTaskHFJetIPQA::GetWeightFactorLinus( AliAODMCParticle * mcpart,bool &isTrackFromPrompt){
	if(!mcpart) {
		return 1;
	}
	//******Check the following cases*****
	// 1. is from primary D Meson
	// 2. is from primary K0s or Lambda
	// 3. is primary charged pion kaon or proton/antiproton
	// 4. is from primary neutral pion /rho / eta / etap / phi / omega meson
	// check if particle is a primary charged pion or proton
	Int_t pPdgcode = abs(mcpart->GetPdgCode());
	Int_t pMotherLabel = 	mcpart->GetMother();
	Bool_t found = kFALSE;
	Double_t pTWeight =0;
	Int_t foundPdg =-1;
	Int_t correctionidx =-1;
	Int_t maPdgcode = 0;
	Int_t maPdgcodeOld = 0;
	AliAODMCParticle * motherOld =0x0;
	bool isprim = false;
	AliAODMCParticle * mcpartMother = 0x0;
	if(pMotherLabel<2) isprim =true;
	else mcpartMother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(pMotherLabel));
	if (!mcpartMother )isprim =true;
	if(!isprim){
		maPdgcode = abs(mcpartMother->GetPdgCode());
		if(maPdgcode>0 && maPdgcode<6) isprim =true;
		if(mcpartMother->GetStatus()>21) isprim =true;
		if(maPdgcode==21) isprim =true;
		if(IsPromptBMeson(mcpartMother)) {
			isTrackFromPrompt = kTRUE;
			if(GetBMesonWeight(mcpartMother,maPdgcode,pTWeight,correctionidx)){
				found = kTRUE;
				foundPdg =maPdgcode;
			}
		}
	}
	if(!isTrackFromPrompt){
		if(pPdgcode == bProton && 	isprim ){
			//Is Primary proton
			found = kTRUE;
			foundPdg =bProton;
			correctionidx = bIdxProton;
			pTWeight =mcpart->Pt();
		}
		else if(pPdgcode == bPi&& 	isprim){
			//Is Primary charged pion
			found = kTRUE;
			foundPdg =bPi;
			correctionidx = bIdxPi;
			pTWeight =mcpart->Pt();
		}
		else if(pPdgcode == bKaon&& 	isprim ){
			//Is Primary charged kaon
			found = kTRUE;
			pTWeight =mcpart->Pt();
			foundPdg =bKaon;
			correctionidx = bIdxKaon;
		}
		else if(!isprim && ParticleIsPossibleSource(maPdgcode)){
			if (IsSelectionParticle(mcpartMother,maPdgcode,pTWeight,correctionidx)) {
				found = kTRUE;
				foundPdg =maPdgcode;
			}
		}
		else {
			while (mcpartMother->GetMother() >0 && !isTrackFromPrompt){
				maPdgcodeOld = maPdgcode;
				motherOld =mcpartMother;
				isprim=false;
				pMotherLabel = 	mcpartMother->GetMother();
				mcpartMother = 0x0;
				mcpartMother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(pMotherLabel)));
				if (!mcpartMother )isprim =true;
				maPdgcode = abs(mcpartMother->GetPdgCode());

				if(IsPromptBMeson(mcpartMother)) {
					isTrackFromPrompt = kTRUE;
					if(GetBMesonWeight(mcpartMother,maPdgcode,pTWeight,correctionidx)){
						found = kTRUE;
						foundPdg =maPdgcode;
						break;
					}}
				if(!isprim){
					maPdgcode = abs(mcpartMother->GetPdgCode());
					if(maPdgcode>0 && maPdgcode<6) isprim =true;
					if(mcpartMother->GetStatus()>21) isprim =true;
					if(maPdgcode==21) isprim =true;
				}
				if(!isprim && ParticleIsPossibleSource(maPdgcode)){
					if (IsSelectionParticle(mcpartMother,maPdgcode,pTWeight,correctionidx)) {
						found = kTRUE;
						foundPdg =maPdgcode;
						break;
					}
				}
				else if (!isprim) continue;
				else break;
			}
		}

	}
	if (!found) {
		return 1;
	}
	Double_t factor = 1;
	//calculate pt of array entry
	// 0.1 -25.GeV 0.05 per bin 498 bins
	double wpos = ((pTWeight - 0.15)/ 0.05);
	double  fractpart, intpart;
	fractpart = modf (wpos , &intpart);
	if (fractpart > 0) intpart = intpart + 1;
	int  bin = floor(intpart);
	if (bin > 497) bin = 497;			// above weight definition
	if (pTWeight < 0.1+ 1E-5) bin = 0; //below weight definition
	factor = fBackgroundFactorLinus[correctionidx][bin];
	if (factor <= 0) return 1;
	else
		return factor;
}

// ######################################################################################## Jet matching 3/4
Double_t AliAnalysisTaskHFJetIPQA::GetWeightFactorLinus( AliMCParticle * mcpart,bool &isTrackFromPrompt){
	if(!mcpart) {
		return 1;
	}
	//******Check the following cases*****
	// 1. is from primary D Meson
	// 2. is from primary K0s or Lambda
	// 3. is primary charged pion kaon or proton/antiproton
	// 4. is from primary neutral pion /rho / eta / etap / phi / omega meson
	// check if particle is a primary charged pion or proton
	Int_t pPdgcode = abs(mcpart->PdgCode());
	Int_t pMotherLabel = mcpart->GetMother();
	Bool_t found = kFALSE;
	Double_t pTWeight =0;
	Int_t foundPdg =-1;
	Int_t correctionidx =-1;
	Int_t maPdgcode = 0;
	Int_t maPdgcodeOld = 0;
	AliMCParticle * motherOld =0x0;
	bool isprim = false;
	AliMCParticle * mcpartMother = 0x0;
	if(pMotherLabel<2) isprim =true;
	else mcpartMother = ((AliMCParticle*)MCEvent()->GetTrack(abs(pMotherLabel)));
	if (!mcpartMother )isprim =true;
	if(!isprim && mcpartMother){
		maPdgcode = abs(mcpartMother->PdgCode());
		if(maPdgcode>0 && maPdgcode<6) isprim =true;
		if(maPdgcode==21) isprim =true;
		if(IsPromptBMeson(mcpartMother)) {
			isTrackFromPrompt = kTRUE;
			if(GetBMesonWeight(mcpartMother,maPdgcode,pTWeight,correctionidx)){
				found = kTRUE;
				foundPdg =maPdgcode;
			}
		}
	}
	if(!isTrackFromPrompt){
		if(pPdgcode == bProton && 	isprim ){
			//Is Primary proton
			found = kTRUE;
			foundPdg =bProton;
			correctionidx = bIdxProton;
			pTWeight =mcpart->Pt();
		}
		else if(pPdgcode == bPi&& 	isprim){
			//Is Primary charged pion
			found = kTRUE;
			foundPdg =bPi;
			correctionidx = bIdxPi;
			pTWeight =mcpart->Pt();
		}
		else if(pPdgcode == bKaon&& 	isprim ){
			//Is Primary charged kaon
			found = kTRUE;
			pTWeight =mcpart->Pt();
			foundPdg =bKaon;
			correctionidx = bIdxKaon;
		}
		else if(!isprim && ParticleIsPossibleSource(maPdgcode)){
			if (IsSelectionParticle(mcpartMother,maPdgcode,pTWeight,correctionidx)) {
				found = kTRUE;
				foundPdg =maPdgcode;
			}
		}
		else {
			while (mcpartMother->GetMother() >0 && !isTrackFromPrompt){
				maPdgcodeOld = maPdgcode;
				motherOld =mcpartMother;
				isprim=false;
				pMotherLabel = 	mcpartMother->GetMother();
				mcpartMother = 0x0;
				mcpartMother =  ((AliMCParticle*)fMCEvent->GetTrack(abs(pMotherLabel)));
				if (!mcpartMother )isprim =true;
				maPdgcode = abs(mcpartMother->PdgCode());
				if(IsPromptBMeson(mcpartMother)) {
					isTrackFromPrompt = kTRUE;
					if(GetBMesonWeight(mcpartMother,maPdgcode,pTWeight,correctionidx)){
						found = kTRUE;
						foundPdg =maPdgcode;
						break;
					}}
				if(!isprim){
					maPdgcode = abs(mcpartMother->PdgCode());
					if(maPdgcode>0 && maPdgcode<6) isprim =true;
					if(maPdgcode==21) isprim =true;
				}
				if(!isprim && ParticleIsPossibleSource(maPdgcode)){
					if (IsSelectionParticle(mcpartMother,maPdgcode,pTWeight,correctionidx)) {
						found = kTRUE;
						foundPdg =maPdgcode;
						break;
					}
				}
				else if (!isprim) continue;
				else break;
			}
		}

	}
	if (!found) {

		return 1;
	}
	// Do the weighting
	Double_t factor = 1;
	//calculate pt of array entry
	// 0.1 -25.GeV 0.05 per bin 498 bins
	double wpos = ((pTWeight - 0.15)/ 0.05);
	double  fractpart, intpart;
	fractpart = modf (wpos , &intpart);
	if (fractpart > 0) intpart = intpart + 1;
	int  bin = floor(intpart);
	if (bin > 497) bin = 497;			// above weight definition
	if (pTWeight < 0.1+ 1E-5) bin = 0; //below weight definition
	factor = fBackgroundFactorLinus[correctionidx][bin];
	if (factor <= 0) return 1;
	else
		return factor;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsSecondaryFromWeakDecay( AliAODMCParticle * particle ) {
	// If a particle is not a physical primary, check if it comes from weak decay
	Int_t mfl = 0;
	Int_t indexMoth = particle->GetMother();
	if(indexMoth < 0) return kFALSE; // if index mother < 0 and not a physical primary, is a non-stable product or one of the beams
	AliAODMCParticle* moth = dynamic_cast<AliAODMCParticle *>(fMCArray->At(indexMoth));
	Int_t pcodemoth = TMath::Abs(particle->PdgCode());
	Int_t codemoth = TMath::Abs(moth->PdgCode());
	// mass of the flavour
	mfl = Int_t (codemoth / TMath::Power(10, Int_t(TMath::Log10(codemoth))));
	if(mfl == 4|| mfl ==5) return kTRUE;
	if(pcodemoth==bPi0){
		if(codemoth == bK0s) return kTRUE;
		if(codemoth == bK0l) return kTRUE;
		if(TMath::Abs(codemoth) == bKaon) return kTRUE;
		if(TMath::Abs(codemoth) == bLambda) return kTRUE;
		if(codemoth == 221 || codemoth == 223 || codemoth == 333 || codemoth == 331 || codemoth == 113 || codemoth == bRhoPlus) return kTRUE;
	}
	else if (pcodemoth==bPhi){
		if(codemoth == 111 || codemoth == 221 || codemoth == 223 || codemoth == 331 || codemoth == 113|| codemoth == bRhoPlus) return kTRUE;
	}
	else if (pcodemoth==bOmega){
		if(codemoth == 111 || codemoth == 221 || codemoth == 333 || codemoth == 331 || codemoth == 113|| codemoth == bRhoPlus) return kTRUE;
	}
	else if (pcodemoth==bEtaPrime){
		if(codemoth == 111 || codemoth == 221 || codemoth == 223 || codemoth == 333 || codemoth == 113|| codemoth == bRhoPlus) return kTRUE;
	}
	else if (pcodemoth==bRho){
		if(codemoth== 111 || codemoth == 221 || codemoth == 223 || codemoth == 333 || codemoth == 331|| codemoth == bRhoPlus) return kTRUE;
	}
	else if (pcodemoth==bEta){
		if(codemoth == 111 || codemoth == 223 || codemoth == 333 || codemoth == 331 || codemoth == 113|| codemoth == bRhoPlus) return kTRUE;    }
	return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsSecondaryFromWeakDecay( AliMCParticle * particle ) {
	// If a particle is not a physical primary, check if it comes from weak decay
	Int_t mfl = 0;
	Int_t indexMoth = particle->GetMother();
	if(indexMoth < 0) return kFALSE; // if index mother < 0 and not a physical primary, is a non-stable product or one of the beams
	AliMCParticle* moth = ((AliMCParticle*)fMCEvent->GetTrack(abs(indexMoth)));
	Int_t pcodemoth = TMath::Abs(particle->PdgCode());
	Int_t codemoth = TMath::Abs(moth->PdgCode());
	// mass of the flavour
	mfl = Int_t (codemoth / TMath::Power(10, Int_t(TMath::Log10(codemoth))));
	if(mfl == 4|| mfl ==5) return kTRUE;
	if(pcodemoth==bPi0){
		if(codemoth == bK0s) return kTRUE;
		if(codemoth == bK0l) return kTRUE;
		if(TMath::Abs(codemoth) == bKaon) return kTRUE;
		if(TMath::Abs(codemoth) == bLambda) return kTRUE;
		if(codemoth == 221 || codemoth == 223 || codemoth == 333 || codemoth == 331 || codemoth == 113 || codemoth == bRhoPlus) return kTRUE;
	}
	else if (pcodemoth==bPhi){
		if(codemoth == 111 || codemoth == 221 || codemoth == 223 || codemoth == 331 || codemoth == 113|| codemoth == bRhoPlus) return kTRUE;
	}
	else if (pcodemoth==bOmega){
		if(codemoth == 111 || codemoth == 221 || codemoth == 333 || codemoth == 331 || codemoth == 113|| codemoth == bRhoPlus) return kTRUE;
	}
	else if (pcodemoth==bEtaPrime){
		if(codemoth == 111 || codemoth == 221 || codemoth == 223 || codemoth == 333 || codemoth == 113|| codemoth == bRhoPlus) return kTRUE;
	}
	else if (pcodemoth==bRho){
		if(codemoth== 111 || codemoth == 221 || codemoth == 223 || codemoth == 333 || codemoth == 331|| codemoth == bRhoPlus) return kTRUE;
	}
	else if (pcodemoth==bEta){
		if(codemoth == 111 || codemoth == 223 || codemoth == 333 || codemoth == 331 || codemoth == 113|| codemoth == bRhoPlus) return kTRUE;    }
	return kFALSE;
}

// ########################################################################################

Bool_t AliAnalysisTaskHFJetIPQA::GetBMesonWeight( AliAODMCParticle * mcpart ,int &pdg,double &pT,int &idx  )
{
	pT = mcpart->Pt();
	switch(pdg){
	case bBPlus:
		idx = bIdxBPlus;
		return kTRUE;
	case bB0:
		idx = bIdxB0;
		return kTRUE;
	case bLambdaB:
		idx = bIdxLambdaB;
		return kTRUE;
		break;
	case bBStarPlus:
		idx = bIdxBStarPlus;
		return kTRUE;
		break;
	default:
		break;
	}
	return kFALSE;
}
// ########################################################################################

Bool_t AliAnalysisTaskHFJetIPQA::GetBMesonWeight( AliMCParticle * mcpart ,int &pdg,double &pT,int &idx  )
{
	pT = mcpart->Pt();
	switch(pdg){
	case bBPlus:
		idx = bIdxBPlus;
		return kTRUE;
	case bB0:
		idx = bIdxB0;
		return kTRUE;
	case bLambdaB:
		idx = bIdxLambdaB;
		return kTRUE;
		break;
	case bBStarPlus:
		idx = bIdxBStarPlus;
		return kTRUE;
		break;
	default:
		break;
	}
	return kFALSE;
}

// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticle( AliAODMCParticle *  mcpart ,int &pdg,double &pT,int &idx ){
	pT = mcpart->Pt();
	bool isprim =false;
	int mpdg = 0;
	if(mcpart->GetMother()>-1)mpdg = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(mcpart->GetMother())))->GetPdgCode();
	if(mpdg>0 && mpdg<6) isprim =true;
	if(mpdg==21) isprim =true;
	if(dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(mcpart->GetMother())))->GetStatus()>11) isprim =true;
	switch(pdg){
	case bPi0:
		if(!IsSecondaryFromWeakDecay(mcpart) ){
			idx = bIdxPi0;
			return kTRUE;
		}
		break;
	case bEta:
		if(!IsSecondaryFromWeakDecay(mcpart) ){
			idx = bIdxEta;
			return kTRUE;
		}
		break;
	case bEtaPrime:
		if(!IsSecondaryFromWeakDecay(mcpart) ){
			idx = bIdxEtaPrime;
			return kTRUE;
		}
		break;
	case bOmega:
		idx = bIdxOmega;
		if(!IsSecondaryFromWeakDecay(mcpart)){
			return kTRUE;
		}
		break;
	case bPhi:
		idx = bIdxPhi;
		if(!IsSecondaryFromWeakDecay(mcpart)){
			return kTRUE;
		}
		break;
	case bRho:
		idx = bIdxRho;
		if(!IsSecondaryFromWeakDecay(mcpart)){
			return kTRUE;
		}
		break;
	case bD0:
		idx = bIdxD0;
		if(IsPromptDMeson(mcpart)){
			return kTRUE;
		}
		break;
	case bDPlus:
		idx = bIdxDPlus;
		if(IsPromptDMeson(mcpart)){
			return kTRUE;
		}
		break;
	case bIdxDSPlus:
		idx = bDSPlus;
		if(IsPromptDMeson(mcpart)){
			return kTRUE;
		}
		break;
	case bDStarPlus:
		idx = bIdxDStarPlus;
		if(IsPromptDMeson(mcpart)){
			return kTRUE;
		}
		break;
	case bLambdaC:
		idx = bIdxLambdaC;
		if(IsPromptDMeson(mcpart)){
			return kTRUE;
		}
		break;
	case bLambda:
		idx = bIdxLambda;
		if(mcpart->IsPhysicalPrimary())
		{
			return kTRUE;
		}
		break;
	case bK0s:
		idx = bIdxK0s;

		if(mcpart->IsPhysicalPrimary())
		{
			return kTRUE;
		}
		break;
	case bPi:
		if(isprim){
			idx =bIdxPi;
			return kTRUE;
		}
		break;
	case bKaon:
		if(isprim){
			idx =bIdxKaon;
			return kTRUE;
		}
		break;
	default:
		break;
	}
	return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsSelectionParticle( AliMCParticle *  mcpart ,int &pdg,double &pT,int &idx ){
	pT = mcpart->Pt();
	bool isprim =false;
	int mpdg = 0;
	if(mcpart->GetMother()>-1)mpdg = ((AliMCParticle*)fMCEvent->GetTrack(abs(mcpart->GetMother())))->PdgCode();
	if(mpdg>0 && mpdg<6) isprim =true;
	if(mpdg==21) isprim =true;
	switch(pdg){
	case bPi0:
		if(!IsSecondaryFromWeakDecay(mcpart) ){
			idx = bIdxPi0;
			return kTRUE;
		}
		break;
	case bEta:
		if(!IsSecondaryFromWeakDecay(mcpart) ){
			idx = bIdxEta;
			return kTRUE;
		}
		break;
	case bEtaPrime:
		if(!IsSecondaryFromWeakDecay(mcpart) ){
			idx = bIdxEtaPrime;
			return kTRUE;
		}
		break;
	case bOmega:
		idx = bIdxOmega;
		if(!IsSecondaryFromWeakDecay(mcpart)){
			return kTRUE;
		}
		break;
	case bPhi:
		idx = bIdxPhi;
		if(!IsSecondaryFromWeakDecay(mcpart)){
			return kTRUE;
		}

		break;
	case bRho:
		idx = bIdxRho;
		if(!IsSecondaryFromWeakDecay(mcpart)){
			return kTRUE;
		}
		break;
	case bD0:
		idx = bIdxD0;
		if(IsPromptDMeson(mcpart)){
			return kTRUE;
		}
		break;
	case bDPlus:
		idx = bIdxDPlus;
		if(IsPromptDMeson(mcpart)){
			return kTRUE;
		}
		break;
	case bIdxDSPlus:
		idx = bDSPlus;
		if(IsPromptDMeson(mcpart)){
			return kTRUE;
		}
		break;
	case bDStarPlus:
		idx = bIdxDStarPlus;
		if(IsPromptDMeson(mcpart)){
			return kTRUE;
		}
		break;
	case bLambdaC:
		idx = bIdxLambdaC;
		if(IsPromptDMeson(mcpart)){
			return kTRUE;
		}
		break;
	case bLambda:
		idx = bIdxLambda;
		if(MCEvent()->IsPhysicalPrimary(mcpart->GetLabel()))
		{
			return kTRUE;
		}
		break;
	case bK0s:
		idx = bIdxK0s;
		if(MCEvent()->IsPhysicalPrimary(mcpart->GetLabel()))
		{
			return kTRUE;
		}
		break;
	case bPi:
		if(isprim){
			idx =bIdxPi;
			return kTRUE;
		}
		break;
	case bKaon:
		if(isprim){
			idx =bIdxKaon;
			return kTRUE;
		}
		break;
	default:
		break;
	}
	return kFALSE;
}
// ########################################################################################
bool AliAnalysisTaskHFJetIPQA::IsPromptBMeson(AliAODMCParticle * part )
{
	if(!part) return false;
	int pdg = TMath::Abs(part->GetPdgCode());
	if ((pdg >= 500 && pdg < 600 )||(pdg >= 5000 && pdg < 6000 ))
	{
		Int_t imo =  part->GetMother();
		AliAODMCParticle* pm = dynamic_cast<AliAODMCParticle *>(fMCArray->At(imo));
		Int_t mpdg = TMath::Abs(pm->GetPdgCode());
		if (!(mpdg >5000 && mpdg <6000) && !(mpdg >500 && mpdg <600))
			return true;
	}
	return false;
}
// ########################################################################################
bool AliAnalysisTaskHFJetIPQA::IsPromptBMeson(AliMCParticle * part )
{
	if(!part) return false;
	int pdg = TMath::Abs(part->PdgCode());
	if ((pdg >= 500 && pdg < 600 )||(pdg >= 5000 && pdg < 6000 ))
	{
		Int_t imo =  part->GetMother();
		AliMCParticle* pm = dynamic_cast<AliMCParticle *>(MCEvent()->GetTrack(imo));
		Int_t mpdg = TMath::Abs(pm->PdgCode());
		if (!(mpdg >5000 && mpdg <6000) && !(mpdg >500 && mpdg <600))
			return true;
	}
	return false;
}
// ########################################################################################
bool AliAnalysisTaskHFJetIPQA::IsPromptDMeson(AliAODMCParticle * part )
{
	if(!part) return false;
	int pdg = TMath::Abs(part->GetPdgCode());
	if ((pdg >= 400 && pdg < 500 )||(pdg >= 4000 && pdg < 5000 ))
	{
		Int_t imo =  part->GetMother();
		AliAODMCParticle* pm = dynamic_cast<AliAODMCParticle *>(fMCArray->At(imo));
		Int_t mpdg = TMath::Abs(pm->GetPdgCode());
		if (!(mpdg >4000 && mpdg <6000) && !(mpdg >400 && mpdg <600))
			return true;
	}
	return false;
}
// ########################################################################################
bool AliAnalysisTaskHFJetIPQA::IsPromptDMeson(AliMCParticle * part )
{
	if(!part) return false;
	int pdg = TMath::Abs(part->PdgCode());
	if ((pdg >= 400 && pdg < 500 )||(pdg >= 4000 && pdg < 5000 ))
	{
		Int_t imo =  part->GetMother();
		AliMCParticle* pm = ((AliMCParticle*)fMCEvent->GetTrack(abs(imo)));
		Int_t mpdg = TMath::Abs(pm->PdgCode());
		if (!(mpdg >4000 && mpdg <6000) && !(mpdg >400 && mpdg <600))
			return true;
	}
	return false;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::ParticleIsPossibleSource(int pdg){
	int pos[19] = {bPi0,bEta,bEtaPrime,bPhi,bRho,bOmega,bK0s,bLambda,bOmegaBaryon,bXiBaryon,bD0,bDPlus,bDStarPlus,bDSPlus,bLambdaB,bLambdaC,bBPlus,bB0,bBStarPlus};
	for (int i =0 ;i<19 ;++i){
		if (abs(pdg)==pos[i] ) return kTRUE;
	}
	return kFALSE;
}
// ######################################################################################## Jet matching 3/4
void AliAnalysisTaskHFJetIPQA::SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, int matching)
{
	Double_t d1 = -1;
	Double_t d2 = -1;

	switch (matching) {
	case 1:
		GetGeometricalMatchingLevel(jet1,jet2,d1);
		d2 = d1;
		break;
	default:
		break;
	}
	if (d1 >= 0) {

		if (d1 < jet1->ClosestJetDistance()) {
			jet1->SetSecondClosestJet(jet1->ClosestJet(), jet1->ClosestJetDistance());
			jet1->SetClosestJet(jet2, d1);
		}
		else if (d1 < jet1->SecondClosestJetDistance()) {
			jet1->SetSecondClosestJet(jet2, d1);
		}
	}
	if (d2 >= 0) {
		if (d2 < jet2->ClosestJetDistance()) {
			jet2->SetSecondClosestJet(jet2->ClosestJet(), jet2->ClosestJetDistance());
			jet2->SetClosestJet(jet1, d2);
		}
		else if (d2 < jet2->SecondClosestJetDistance()) {
			jet2->SetSecondClosestJet(jet1, d2);
		}
	}
}

// ######################################################################################## Jet matching 4/4
void AliAnalysisTaskHFJetIPQA::GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const
{
	Double_t deta = jet2->Eta() - jet1->Eta();
	Double_t dphi = jet2->Phi() - jet1->Phi();
	dphi = TVector2::Phi_mpi_pi(dphi);
	d = TMath::Sqrt(deta * deta + dphi * dphi);
}
// ######################################################################################## Monte Carlo correction factors
Double_t AliAnalysisTaskHFJetIPQA::GetMonteCarloCorrectionFactor(AliAODTrack* track,bool &ise,bool &fromB){
	AliAODMCParticle *pMC = 0x0;
	AliAODEvent* aev = dynamic_cast<AliAODEvent*>(InputEvent());
	pMC = GetMCTrack(track);
	if(pMC) 	if(TMath::Abs(pMC->PdgCode())==11)ise = true;
	double val = GetWeightFactorLinus(pMC,fromB);
	if(val > 0 )return val;
	return 1.;
}
// ######################################################################################## Monte Carlo correction factors
Double_t AliAnalysisTaskHFJetIPQA::GetMonteCarloCorrectionFactor(AliESDtrack* track,bool &ise,bool &fromB){
	AliMCParticle *pMC = 0x0;
	AliESDEvent* eev = dynamic_cast<AliESDEvent*>(InputEvent());
	pMC = ((AliMCParticle*)fMCEvent->GetTrack(abs(track->GetLabel())));
	if(!pMC) return 1;
	if(pMC) 	if(TMath::Abs(pMC->PdgCode())==11)ise = true;
	double val = GetWeightFactorLinus(pMC,fromB);
	if(val > 0 )return val;
	return 1.;
}
/*
// ######################################################################################## Monte Carlo correction factors
Double_t AliAnalysisTaskHFJetIPQA::GetMonteCarloCorrectionFactor(AliAODTrack* track){
	// Todo: Best handle on missing mc information ... currently scale with 1
	AliAODMCParticle *pMC = 0x0;
	AliAODEvent* aev = dynamic_cast<AliAODEvent*>(InputEvent());
	AliAODMCHeader* mcheader = dynamic_cast<AliAODMCHeader*>(aev->FindListObject(AliAODMCHeader::StdBranchName()));
	if (!mcheader) return 1;
	TClonesArray * fMCparticles = dynamic_cast<TClonesArray*>(aev->FindListObject(AliAODMCParticle::StdBranchName()));
	if (!fMCparticles) return 1;
	if(track->GetLabel()>-1)
		pMC = dynamic_cast<AliAODMCParticle*>(fMCparticles->At(track->GetLabel()));
	if (!pMC) return 1;

	vector <int> chain;
	vector <int> chain_pt;

	// follow till primary particle
	AliAODMCParticle * mother =0x0;

	chain.push_back(pMC->GetPdgCode());
	chain_pt.push_back(pMC->Pt());
	// Check if mother is the projectile
	if(pMC->GetMother()>1){
		int label = pMC->GetMother();
		mother = dynamic_cast<AliAODMCParticle*>(fMCparticles->At(label));
		while(!(mother == 0 )){
			if(!mother) break;
			chain.push_back(mother->GetPdgCode());
			chain_pt.push_back(mother->Pt());

			if(mother->GetMother()>1){
				label = mother->GetMother();
				mother =0x0;
				mother = dynamic_cast<AliAODMCParticle*>(fMCparticles->At(label));
				continue;
			} else break;
		}
	}
	// Print decay history (debug only... to be removed)
	// todo store map to have momentum for correction
	std::reverse(chain.begin(),chain.end());
	std::reverse(chain_pt.begin(),chain_pt.end());
	Printf("╭┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈╮");
	for (int i=0; i<chain.size();++i){
		if(i==0) printf("┊");
		printf("%i",chain.at(i));
		if(i<chain.size()-1) printf("->");
		else printf("┊\n");
	}
	Printf("╰┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈╯");

	return 1.;
}*/
// ######################################################################################## Monte Carlo correction factors
Int_t AliAnalysisTaskHFJetIPQA::GetElecSource(const AliAODMCParticle *  const mcpart,Double_t &mpt) const
{
	//
	// Function for AliAODMCParticle
	//

	if (!fMCArray) return -1;
	////////////////
	Int_t origin = -1;
	Bool_t isFinalOpenCharm = kFALSE;

	Int_t iLabel = mcpart->GetMother();
	if ((iLabel<0) || (iLabel>=fMCArray->GetEntriesFast())){
		AliDebug(1, "label is out of range, return\n");
		return -1;
	}

	AliAODMCParticle *mctrack = NULL; // will change all the time
	Int_t tmpMomLabel=0;
	if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))) return -1;
	AliAODMCParticle *partMother = mctrack;
	AliAODMCParticle *partMotherCopy = mctrack;
	Int_t maPdgcode = mctrack->GetPdgCode();
	mpt = partMother->Pt();

	Int_t grmaPdgcode;
	Int_t ggrmaPdgcode;

	// if the mother is charmed hadron

	if(TMath::Abs(maPdgcode)==443){
		//
		// J/spi
		//
		Int_t jLabel = partMother->GetMother();
		if ((jLabel>=0) && (jLabel<fMCArray->GetEntriesFast())){
			if((mctrack = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(jLabel))))){
				Int_t grandMaPDG = mctrack->GetPdgCode();
				mpt = mctrack->Pt();
				if((int(TMath::Abs(grandMaPDG)/100.)%10) == kBeauty || (int(TMath::Abs(grandMaPDG)/1000.)%10) == kBeauty) {
					return kB2Jpsi;
				}
			}
		}
		return kJpsi;
	}
	else if ( (int(TMath::Abs(maPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(maPdgcode)/1000.)%10) == kCharm ) {
		//
		// charm
		//
		for (Int_t i=0; i<fNparents; i++){
			if (TMath::Abs(maPdgcode)==fParentSelect[0][i]){
				isFinalOpenCharm = kTRUE;
			}
		}
		if (!isFinalOpenCharm) {
			return -1;
		}

		// iterate until you find B hadron as a mother or become top ancester
		for (Int_t i=1; i<fgkMaxIter; i++){

			Int_t jLabel = partMother->GetMother();
			if (jLabel == -1){
				return kDirectCharm;
			}
			if ((jLabel<0) || (jLabel>=fMCArray->GetEntriesFast())){
				AliDebug(1, "Stack label is negative, return\n");
				return -1;
			}

			// if there is an ancester
			if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(jLabel))))) {
				return -1;
			}
			Int_t grandMaPDG = mctrack->GetPdgCode();
			for (Int_t j=0; j<fNparents; j++){
				if (TMath::Abs(grandMaPDG)==fParentSelect[1][j]){
					mpt = mctrack->Pt();
					return kBeautyCharm;
				}
			}
			partMother = mctrack;
		} // end of iteration

	} // end of if
	else if ( (int(TMath::Abs(maPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(maPdgcode)/1000.)%10) == kBeauty ) {
		//
		// beauty
		//
		for (Int_t i=0; i<fNparents; i++){
			if (TMath::Abs(maPdgcode)==fParentSelect[1][i]){
				mpt = partMotherCopy->Pt();
				return kDirectBeauty;
			}
		}
	} // end of if
	else if ( TMath::Abs(maPdgcode) == 22 ) {
		//
		//conversion
		//
		tmpMomLabel = partMotherCopy->GetMother();
		if(tmpMomLabel==-1) return kGamma;
		if((tmpMomLabel<0) || (tmpMomLabel>=fMCArray->GetEntriesFast())) {
			return -1;
		}
		if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(tmpMomLabel))))) {
			return -1;
		}
		partMother = mctrack;
		maPdgcode = partMother->GetPdgCode();

		// check if the ligth meson is the decay product of heavy mesons
		tmpMomLabel = partMother->GetMother();
		if((tmpMomLabel>=0) && (tmpMomLabel<fMCArray->GetEntriesFast())) {//grandgrandmother
			if((mctrack = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(tmpMomLabel))))) {
				partMother = mctrack;
				grmaPdgcode = partMother->GetPdgCode();
				mpt = partMother->Pt();

				if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kBeauty ) {
					return kGammaB2M;
				}
				if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kCharm ) {
					return kGammaD2M;
				}

				tmpMomLabel = partMother->GetMother();
				if((tmpMomLabel>=0) && (tmpMomLabel<fMCArray->GetEntriesFast())) {//grandgrandgrandmother
					if((mctrack = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(tmpMomLabel))))) {
						partMother = mctrack;
						ggrmaPdgcode = partMother->GetPdgCode();
						mpt = partMother->Pt();

						if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == kBeauty ) {
							return kGammaB2M;
						}
						if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == kCharm ) {
							return kGammaD2M;
						}
					}
				}//grandgrandgrandmother

				if ( TMath::Abs(maPdgcode) == 111 ) {
					if(grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
					else if(grmaPdgcode == 310) return kGammaK0s2P;
					else if(grmaPdgcode == 130) return kGammaK0l2P;
					else if(TMath::Abs(grmaPdgcode) == 321) return kGammaK2P;
					else if(TMath::Abs(grmaPdgcode) == 3122) return kGammaLamda2P;
					else if(grmaPdgcode == 3222) return kGammaSigma2P;
					return kGammaPi0;
				}
				else if ( TMath::Abs(maPdgcode) == 221 ) {
					if(grmaPdgcode == 111 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
					return kGammaEta;
				}
				else if ( TMath::Abs(maPdgcode) == 223 ) {
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
					return kGammaOmega;
				}
				else if ( TMath::Abs(maPdgcode) == 333 ) {
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
					return kGammaPhi;
				}
				else if ( TMath::Abs(maPdgcode) == 331 ) {
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 113) return kGammaM2M;
					return kGammaEtaPrime;
				}
				else if ( TMath::Abs(maPdgcode) == 113 ) {
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331) return kGammaM2M;
					return kGammaRho0;
				}
				else origin = kElse;//grandgrandmother but nothing we identify
			}//mctrack grandgrandmother
		}
		else {
			// grandmother is primary
			if ( TMath::Abs(maPdgcode) == 111 ) {
				return kGammaPi0;
			}
			else if ( TMath::Abs(maPdgcode) == 221 ) {
				return kGammaEta;
			}
			else if ( TMath::Abs(maPdgcode) == 223 ) {
				return kGammaOmega;
			}
			else if ( TMath::Abs(maPdgcode) == 333 ) {
				return kGammaPhi;
			}
			else if ( TMath::Abs(maPdgcode) == 331 ) {
				return kGammaEtaPrime;
			}
			else if ( TMath::Abs(maPdgcode) == 113 ) {
				return kGammaRho0;
			}
			else origin = kElse;//grandmother is primary but nothing we identify
		}

		return origin;

	}
	else {
		//
		// check if the ligth meson is the decay product of heavy mesons
		//
		tmpMomLabel = partMotherCopy->GetMother();
		if((tmpMomLabel>=0) && (tmpMomLabel<fMCArray->GetEntriesFast())) {//grandmother
			if((mctrack = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(tmpMomLabel))))) {
				partMother = mctrack;
				grmaPdgcode = partMother->GetPdgCode();

				if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kBeauty ) {
					return kB2M;
				}
				if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kCharm ) {
					return kD2M;
				}

				tmpMomLabel = partMother->GetMother();
				if((tmpMomLabel>=0) && (tmpMomLabel<fMCArray->GetEntriesFast())) {//grandgrandmother
					if((mctrack = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(tmpMomLabel))))) {
						partMother = mctrack;
						ggrmaPdgcode = partMother->GetPdgCode();

						if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == kBeauty ) {
							return kB2M;
						}
						if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == kCharm ) {
							return kD2M;
						}
					}
				}//grandgrandmother

				if ( TMath::Abs(maPdgcode) == 111 ) {
					if(grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
					else if(grmaPdgcode == 310) return kK0s2P;
					else if(grmaPdgcode == 130) return kK0l2P;
					else if(TMath::Abs(grmaPdgcode) == 321) return kK2P;
					else if(TMath::Abs(grmaPdgcode) == 3122) return kLamda2P;
					else if(grmaPdgcode == 3222) return kSigma2P;
					return kPi0;
				}
				else if ( TMath::Abs(maPdgcode) == 221 ) {
					if(grmaPdgcode == 111 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
					return kEta;
				}
				else if ( TMath::Abs(maPdgcode) == 223 ) {
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
					return kOmega;
				}
				else if ( TMath::Abs(maPdgcode) == 333 ) {
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
					return kPhi;
				}
				else if ( TMath::Abs(maPdgcode) == 331 ) {
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 113) return kM2M;
					return kEtaPrime;
				}
				else if ( TMath::Abs(maPdgcode) == 113 ) {
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331) return kM2M;
					return kRho0;
				}
				else if ( TMath::Abs(maPdgcode) == 321 ) {
					return kKe3;
				}
				else if ( TMath::Abs(maPdgcode) == 130 ) {
					return kK0L;
				}
				else origin = kElse;//grandmother but nothing we identidy
			}//mctrack grandmother
		}
		else {
			// no grandmother
			if ( TMath::Abs(maPdgcode) == 111 ) {
				return kPi0;
			}
			else if ( TMath::Abs(maPdgcode) == 221 ) {
				return kEta;
			}
			else if ( TMath::Abs(maPdgcode) == 223 ) {
				return kOmega;
			}
			else if ( TMath::Abs(maPdgcode) == 333 ) {
				return kPhi;
			}
			else if ( TMath::Abs(maPdgcode) == 331 ) {
				return kEtaPrime;
			}
			else if ( TMath::Abs(maPdgcode) == 113 ) {
				return kRho0;
			}
			else if ( TMath::Abs(maPdgcode) == 321 ) {
				return kKe3;
			}
			else if ( TMath::Abs(maPdgcode) == 130 ) {
				return kK0L;
			}
			else origin = kElse;//mother but nothing we identify
		}
	}//mother is something different from J/psi,charm,beauty or gamma
	return origin;
}

//#########
//__________________________________________
Double_t AliAnalysisTaskHFJetIPQA::GetWeightFactor(const AliAODMCParticle * const mcpart, const Int_t iBgLevel)
{
	//
	// Get weighting factor for the realistic background estimation, for three possible background yield levels, indicated by the argument "iLevel": the best estimate (0), the lower uncertainty level (1), and the upper uncertainty level (2)
	//
	Double_t weightElecBg = 0.; // make 0 again
	Double_t mesonPt = 0.;
	Double_t mesonMotherPt = 0.;
	Double_t bgcategory = 0.;
	Bool_t condition = kTRUE;
	Int_t mArr = -1;
	Double_t mpt=0;
	if (!mcpart) return 0;
	Int_t mesonID = GetElecSource(mcpart,mpt);
	//Printf("mesonID = %i pt = %f",mesonID,mpt);
	//return 0.;
	if(mesonID==kGammaPi0 || mesonID==kPi0) mArr=0;                //pion
	else if(mesonID==kGammaEta || mesonID==kEta) mArr=1;           //eta
	else if(mesonID==kGammaOmega || mesonID==kOmega) mArr=2;       //omega
	else if(mesonID==kGammaPhi || mesonID==kPhi) mArr=3;           //phi
	else if(mesonID==kGammaEtaPrime || mesonID==kEtaPrime) mArr=4; //etaprime
	else if(mesonID==kGammaRho0 || mesonID==kRho0) mArr=5;         //rho
	else if(mesonID==kKe3 || mesonID==kGammaK2P|| mesonID==kK2P) mArr=6; //ke3 or K->pi->e
	else if(mesonID==kK0L || mesonID==kGammaK0s2P|| mesonID==kK0s2P) mArr=7; //K0L->e+X or k0s->pi->e
	else if(mesonID==kGammaLamda2P|| mesonID==kLamda2P) mArr=8;    //lambda->pi->e

	Double_t datamc[30]={-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999,-999};
	Double_t xr[3]={-999,-999,-999};
	datamc[0] = mesonID;
	datamc[17] = mcpart->Pt(); //electron pt
	datamc[18] = mcpart->Eta(); //electron eta

	mcpart->XvYvZv(xr);
	datamc[9] = TMath::Sqrt(xr[0]*xr[0]+xr[1]*xr[1]);
	datamc[10] = xr[2];
	datamc[19] = (mcpart->IsPrimary()) ? 1 : 0;
	datamc[20] = (mcpart->IsPhysicalPrimary()) ? 1 : 0;

	datamc[24] = mcpart->Label();

	if(!(mArr<0)){

		AliAODMCParticle *mctrackmother = NULL; // will change all the time

		if((mesonID>=kGammaPi0 && mesonID<=kGammaRho0) || mesonID==kGammaK2P || mesonID==kGammaK0s2P || mesonID==kGammaLamda2P) { // conversion electron
			Int_t iLabel = mcpart->GetMother(); //gamma label
			if(!(mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))) return 0;
			iLabel = mctrackmother->GetMother(); //gamma's mother's label
			if(!(mctrackmother= dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))) return 0;
			mesonPt = mctrackmother->Pt(); //meson pt
			datamc[28] = mctrackmother->Eta(); //meson eta
			bgcategory = 1.;
			if(TMath::Abs((mctrackmother->Y()))>0.8) condition = kFALSE;
			if(!mctrackmother->IsPrimary()) condition = kFALSE;
			datamc[1] = bgcategory;
			datamc[2] = mesonPt;
			mctrackmother->XvYvZv(xr);
			datamc[11] = TMath::Sqrt(xr[0]*xr[0]+xr[1]*xr[1]);
			datamc[12] = xr[2];

			datamc[21] = (mctrackmother->IsPrimary()) ? 1 : 0;
			datamc[22] = (mctrackmother->IsPhysicalPrimary()) ? 1 : 0;

			//bgcategory 2, 3 is not defined for AOD

			iLabel=mctrackmother->GetMother(); // gamma's mother's mother
			datamc[26] = iLabel;
			if((mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))){
				datamc[3]=mctrackmother->PdgCode();
				mesonMotherPt=mctrackmother->Pt();
				datamc[4]=mesonMotherPt;
				datamc[29] = mctrackmother->Eta(); //meson mother's eta
				if(TMath::Abs(mctrackmother->PdgCode())==310){
					datamc[21] = (mctrackmother->IsPrimary()) ? 1 : 0;
					datamc[22] = (mctrackmother->IsPhysicalPrimary()) ? 1 : 0;
					iLabel=mctrackmother->GetMother(); // gamma's mother's mother's mother
					if((mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))){
						datamc[5]=mctrackmother->PdgCode();
						iLabel=mctrackmother->GetMother(); // gamma's mother's mother's mother
						if((mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))){
							datamc[6]=mctrackmother->PdgCode();
						}
					}
				}
			}
		}
		else{ // nonHFE except for the conversion electron
			Int_t iLabel = mcpart->GetMother(); //meson label
			if(!(mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))) return 0;
			mesonPt = mctrackmother->Pt(); //meson pt
			datamc[28] = mctrackmother->Eta(); //meson eta
			if(mesonID==kEta) bgcategory = -1.41; // to consider new branching ratio for the eta Dalitz decay
			else 	bgcategory = -1.;
			if(TMath::Abs((mctrackmother->Y()))>0.8) condition = kFALSE;
			if(!mctrackmother->IsPrimary()) condition = kFALSE;
			datamc[1] = bgcategory;
			datamc[2] = mesonPt;
			datamc[23] = mctrackmother->PdgCode();
			mctrackmother->XvYvZv(xr);
			datamc[11] = TMath::Sqrt(xr[0]*xr[0]+xr[1]*xr[1]);
			datamc[12] = xr[2];

			datamc[21] = (mctrackmother->IsPrimary()) ? 1 : 0;
			datamc[22] = (mctrackmother->IsPhysicalPrimary()) ? 1 : 0;

			//bgcategory 2, 3 is not defined for AOD

			iLabel=mctrackmother->GetMother(); // mesons' mother
			datamc[26] = iLabel;
			if((mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))){
				datamc[3]=mctrackmother->PdgCode();
				mesonMotherPt=mctrackmother->Pt();
				datamc[4]=mesonMotherPt;
				datamc[29] = mctrackmother->Eta(); //meson mother's eta
				if(TMath::Abs(mctrackmother->PdgCode())==310){
					datamc[21] = (mctrackmother->IsPrimary()) ? 1 : 0;
					datamc[22] = (mctrackmother->IsPhysicalPrimary()) ? 1 : 0;
					iLabel=mctrackmother->GetMother(); // mesons' mother's mother
					if((mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))){
						datamc[5]=mctrackmother->PdgCode();
						iLabel=mctrackmother->GetMother(); // meson's mother's mother's mother
						if((mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))){
							datamc[6]=mctrackmother->PdgCode();
						}
					}
				}
			}
		}

		int kBgPtBins = 45;
		weightElecBg=fBackgroundFactor[mArr][kBgPtBins-1];

		if(mArr<=5 || mesonID==kKe3 || mesonID==kK0L){
			for(int ii=0; ii<kBgPtBins; ii++){
				if((mesonPt > fBackgroundFactorBins[ii]) && (mesonPt < fBackgroundFactorBins[ii+1])){
					weightElecBg = fBackgroundFactor[mArr][ii];

					break;
				}
			}
		}
		else{
			for(int ii=0; ii<kBgPtBins; ii++){
				if((mesonMotherPt > fBackgroundFactorBins[ii]) && (mesonMotherPt < fBackgroundFactorBins[ii+1])){
					weightElecBg = fBackgroundFactor[mArr][ii];
					break;
				}
			}
		}
	}

	if(!condition) weightElecBg=0.;


	Double_t returnval = bgcategory*weightElecBg;

	return returnval;
}

Float_t AliAnalysisTaskHFJetIPQA::GetRapidity(const TParticle *part){
	//
	// return rapidity
	//
	Float_t rapidity;
	if(!((part->Energy() - part->Pz())*(part->Energy() + part->Pz())>0)) rapidity=-999;
	else rapidity = 0.5*(TMath::Log((part->Energy()+part->Pz()) / (part->Energy()-part->Pz())));
	return rapidity;
}


bool AliAnalysisTaskHFJetIPQA::mysort(const myvaluetuple& i, const myvaluetuple& j)
{
	if(i.first <= j.first)
		return false;
	else
		return true;
}
//

AliAODMCParticle* AliAnalysisTaskHFJetIPQA::GetMCTrack( const AliAODTrack* _track)
{
	//
	// return MC track
	//
	if(!fMCArray) { AliError("No fMCArray"); return NULL;}
	Int_t nStack = fMCArray->GetEntriesFast();
	Int_t label  = TMath::Abs(_track->GetLabel()); // negative label indicate poor matching quality
	if(label > nStack) return NULL;
	AliAODMCParticle *mctrack = (AliAODMCParticle*)fMCArray->At(label);
	return mctrack;
}

//
Bool_t AliAnalysisTaskHFJetIPQA::IsV0PhotonFromBeamPipeDaughter(const AliAODTrack* track)
{
	if(!track)return kFALSE;
	AliAODv0* v0aod = 0x0;
	int posid = -1;
	int negid = -1;
	int trackid = -1;
	Double_t P[3];


	for(int i = 0; i < InputEvent()->GetNumberOfV0s(); ++i) {
		P[0]=0.;
		P[1]=0.;
		P[2]=0.;
		v0aod = ((AliAODEvent*)InputEvent())->GetV0(i);
		if (!v0aod->GetOnFlyStatus()) continue;
		posid = v0aod->GetPosID();
		negid = v0aod->GetNegID();
		trackid = track->GetID();
		if(posid == trackid || negid == trackid) {
			P[0] = v0aod->DecayVertexV0X();
			P[1] = v0aod->DecayVertexV0Y();
			P[2] = v0aod->DecayVertexV0Z();
			double Radius = sqrt(P[0] * P[0] + P[1] * P[1]);
			if(Radius < 800.) {
				//Try To construct gamma from daughters

				AliVTrack* posTrack = ((AliVTrack*)(InputEvent()->GetTrack(posid))) ;
				AliVTrack* negTrack = ((AliVTrack*)(InputEvent()->GetTrack(negid))) ;

				if(!posTrack) continue;
				if(!negTrack) continue;
				AliKFParticle::SetField((((AliAODEvent*)(InputEvent()))->GetMagneticField()));
				AliKFParticle pos(*posTrack,-11);
				AliKFParticle neg(*negTrack,11);
				AliKFParticle partGamma(pos,neg,true);

				// Caluclate Armenteros Podolanski cut
				double pPlus[3]  ={	pos.GetPx() ,	pos.GetPy() ,	pos.GetPz() };
				double pMinus[3] ={	neg.GetPx() ,	neg.GetPy() ,	neg.GetPz() };
				double PV0[3] ={	partGamma.GetPx() ,	partGamma.GetPy() ,	partGamma.GetPz() };

				double PLPlus = (PV0[0]*pPlus[0] +PV0[1]*pPlus[1]+PV0[2]*pPlus[2])/fabs(partGamma.GetP());
				double PLMinus = (PV0[0]*pMinus[0] +PV0[1]*pMinus[1]+PV0[2]*pMinus[2])/fabs(partGamma.GetP());

				/*
				AliAODMCParticle * pTruePos = GetMCTrack(((AliAODTrack*)(InputEvent()->GetTrack(posid))));
				AliAODMCParticle * pMother =0x0;
				if (pTruePos)pMother =  dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(pTruePos->GetMother())));
				int pdg = 0;
				if (pMother) pdg = pMother->PdgCode();
				 */
				double alpha = PLPlus-PLMinus/(PLPlus+PLMinus);

				TVector3 pPlusV3(pPlus[0],pPlus[1],pPlus[2]);
				TVector3 pMinusV3(pMinus[0],pMinus[1],pMinus[2]);

				TVector3 pV0V3(PV0[0],PV0[1],PV0[2]);
				TVector3 v = pPlusV3.Cross(pV0V3);
				double pPerp = v.Mag();
				double qT =pPerp /pV0V3.Mag();
				// Calculate psi Pair value
				double xipair = acos((pPlusV3 * pMinusV3)/(pPlusV3.Mag() * pMinusV3.Mag()));

				double theataPlus = ((AliVTrack*)(InputEvent()->GetTrack(posid)))->Theta();
				double theataMinus = ((AliVTrack*)(InputEvent()->GetTrack(negid)))->Theta();
				double delta_theta =theataMinus-theataPlus;

				double psi_pair = asin(	delta_theta/xipair);


				Printf("qt %f/%f",v0aod->PtArmV0(),qT );



				//if(pdg==22)Printf("%i:  X^2/NDF %f  qT %f Psipair %f Mother PDG %i",i,partGamma.GetChi2()/partGamma.GetNDF(),qT,psi_pair,pdg );

				if((partGamma.Chi2()/partGamma.NDF() <30.) && (qT < 0.05) && (fabs(psi_pair) < 0.05) ) {
					Printf("Rejected as Gamma");
					return kTRUE;
				}
			}
			return kFALSE;
		}
	}
	return kFALSE;
}
// ########################################################################################
Bool_t AliAnalysisTaskHFJetIPQA::IsV0PhotonFromBeamPipeDaughter(const AliESDtrack* track)
{
	//Still need inclusion of Track + PID cuts


	if(!track)return kFALSE;
	AliESDv0* v0esd = 0x0;
	int posid = -1;
	int negid = -1;
	int trackid = -1;
	Double_t P[3];
	for(int i = 0; i < InputEvent()->GetNumberOfV0s(); ++i) {
		P[0]=0.;
		P[1]=0.;
		P[2]=0.;
		v0esd = ((AliESDEvent*)InputEvent())->GetV0(i);
		if (!v0esd->GetOnFlyStatus())  continue;

		v0esd->XvYvZv(P);

		double Radius = sqrt(P[0] * P[0] + P[1] * P[1]);
		if(Radius < 800.) {
			//Try To construct gamma from daughters

			AliKFParticle::SetField((((AliESDEvent*)(InputEvent()))->GetMagneticField()));
			AliKFParticle pos(*(v0esd->GetParamP()),-11);
			AliKFParticle neg(*(v0esd->GetParamN()),11);
			AliKFParticle partGamma(pos,neg,true);
			// Caluclate Armenteros Podolanski cut
			int NIndex =v0esd->GetNindex();
			int PIndex =v0esd->GetPindex();
			int PLabel = ((AliESDEvent*)InputEvent())->GetTrack(PIndex)->GetLabel();
			int NLabel = ((AliESDEvent*)InputEvent())->GetTrack(NIndex)->GetLabel();
			//if(!(PLabel==abs(track->GetLabel()))&&(!(NLabel==abs(track->GetLabel())))) continue;
			//DEBUG
			/*AliMCParticle * pTruePos = (AliMCParticle*)fMCEvent->GetTrack(abs(PLabel));
			AliMCParticle * pMother =0x0;
			if (pTruePos)pMother =  (AliMCParticle*)fMCEvent->GetTrack(pTruePos->GetMother()); //dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(pTruePos->GetMother())));
			int pdg = 0;
			if (pMother) pdg = pMother->PdgCode();
			AliMCParticle * pTrueNeg = (AliMCParticle*)fMCEvent->GetTrack(abs(NLabel));
			AliMCParticle * nMother =0x0;
			if (pTrueNeg)nMother =  (AliMCParticle*)fMCEvent->GetTrack(pTrueNeg->GetMother()); //dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(pTruePos->GetMother())));
			int pdgn = 0;
			if (nMother) pdgn = nMother->PdgCode();
			*/
			//DEBUG
			double alpha =0;
			double qT = GetArmenteros(v0esd,11,-11,alpha);
			double psi_pair = 	GetPsiPair(v0esd);
		/*	if(pdg==22)Printf("pP X^2/NDF %f  qT %f Psipair %f Mother PDG %i",partGamma.GetChi2()/partGamma.GetNDF(),qT,psi_pair,pdg );
			if(pdgn==22)Printf("nP X^2/NDF %f  qT %f Psipair %f Mother PDG %i",partGamma.GetChi2()/partGamma.GetNDF(),qT,psi_pair,pdgn );
*/
			if((partGamma.GetChi2()/partGamma.GetNDF() <30.) && (qT < 0.05) && (fabs(psi_pair) < 0.05) ) {
				Printf("Rejected as Gamma");
				return kTRUE;

			}
			return kFALSE;
		}
	}
	return kFALSE;
}
//###############################################################################################################
Double_t AliAnalysisTaskHFJetIPQA::GetArmenteros(AliESDv0 * v0, int pidneg,int pidpos ,double &alpha){
	AliKFParticle pos(*(v0->GetParamP()),pidpos);
	AliKFParticle neg(*(v0->GetParamN()),pidneg);
	AliKFParticle partGamma(pos,neg,true);
	Double_t armenteros[2]={0,0};
	partGamma.GetArmenterosPodolanski(pos,neg, armenteros);
	alpha= armenteros[1];
	return armenteros[0];
}
//###############################################################################################################
Double_t AliAnalysisTaskHFJetIPQA::GetPsiPair(AliESDv0 * v0){
	AliExternalTrackParam nt(*v0->GetParamN());
	AliExternalTrackParam pt(*v0->GetParamP());
	Float_t magField = fInputEvent->GetMagneticField();
	Double_t xyz[3] = {0.,0.,0.};
	v0->GetXYZ(xyz[0],xyz[1],xyz[2]);
	Double_t mn[3] = {0,0,0};
	Double_t mp[3] = {0,0,0};
	v0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter;
	v0->GetPPxPyPz(mp[0],mp[1],mp[2]);
	Double_t deltat = 1.;
	deltat = TMath::ATan(mp[2]/(TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1])+1.e-13)) - TMath::ATan(mn[2]/(TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1])+1.e-13));//difference of angles of the two daughter tracks with z-axis
	Double_t radiussum = TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]) + 50;//radius to which tracks shall be propagated
	Double_t momPosProp[3] = {0,0,0};
	Double_t momNegProp[3] = {0,0,0};
	Double_t psiPair = 4.;
	if(nt.PropagateTo(radiussum,magField) == 0) return psiPair; //propagate tracks to the outside -> Better Purity and Efficiency
	if(pt.PropagateTo(radiussum,magField) == 0) return psiPair; //propagate tracks to the outside -> Better Purity and Efficiency
	pt.GetPxPyPz(momPosProp);//Get momentum vectors of tracks after propagation
	nt.GetPxPyPz(momNegProp);
	Double_t pEle =
			TMath::Sqrt(momNegProp[0]*momNegProp[0]+momNegProp[1]*momNegProp[1]+momNegProp[2]*momNegProp[2]);//absolute momentum value of negative daughter
	Double_t pPos =
			TMath::Sqrt(momPosProp[0]*momPosProp[0]+momPosProp[1]*momPosProp[1]+momPosProp[2]*momPosProp[2]);//absolute momentum value of positive daughter
	Double_t scalarproduct =
			momPosProp[0]*momNegProp[0]+momPosProp[1]*momNegProp[1]+momPosProp[2]*momNegProp[2];//scalar product of propagated positive and negative daughters' momenta
	Double_t chipair = TMath::ACos(scalarproduct/(pEle*pPos));//Angle between propagated daughter tracks
	psiPair = TMath::ASin(deltat/chipair);
	return psiPair;
}

