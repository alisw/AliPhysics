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


/* AliAnalysisTaskEPCorrAA.cxx
Author : minwoo.kim@cern.ch

 */
#include "AliAnalysisTaskEPCorrAA.h"

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TMath.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliCentrality.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliHeader.h"
#include "AliPID.h"
#include "AliPIDResponse.h"

#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"

//#include "AliEPSelectionTask.h"
#include "AliEventplane.h"
#include "AliVZEROEPSelectionTask.h"

#include "AliExternalTrackParam.h"

#include "AliEventPoolManager.h"

#include "AliVEvent.h"

using namespace std;

ClassImp(AliAnalysisTaskEPCorrAA)
ClassImp(AliCorrReducedTrackAA)

/*
//________________________________________________________________________
AliAnalysisTaskEPCorrAA::AliAnalysisTaskEPCorrAA() // All data members should be initialised here
	:AliAnalysisTaskSE(),
	fOutput(0),
//	fTrackCuts(0),
	fPoolMgr(0x0),
	fMyprimRecoTracks(0x0),
	fTracksMixing(0x0),
	fMinNumTrack(50000),
	fPoolSize(100), // from 1000 (150610)
	fMinNEventsToMix(5),
	fNsigmaCut(3.0),
	fHistZvertex(0),
	fPIDResponse(0),
	fHistCent(0) // The last in the above list should not have a comma after it
{
	// Dummy constructor ALWAYS needed for I/O.
	for(int ic=0; ic<kCentBin;ic++) {
//		fHistNsigmaTPCpT[ic]=0x0;
//		fHistNsigmaTOFpT[ic]=0x0;
//		fHistNsigmaITSpT[ic]=0x0;
		for(int iptt=0; iptt<kpTBin;iptt++) {
//			fHistNsigmaITSTPC[ic][iptt]=0x0;
//			fHistNsigmaTPCTOF[ic][iptt]=0x0;
//			fHistNsigmaITSTOF[ic][iptt]=0x0;
			for(int ipta=0; ipta<kpTBin; ipta++) {
				for(int iz=0; iz<kZvertBin;iz++) {
					for(int ipid=0; ipid<kPID;ipid++) {
						fHistdEtadPhiSame[ic][iz][iptt][ipta][ipid]=0x0;
						fHistdEtadPhiMixed[ic][iz][iptt][ipta][ipid]=0x0;
					}
					fHistdEtadPhiSameH[ic][iz][iptt][ipta]=0x0;
					fHistdEtadPhiMixedH[ic][iz][iptt][ipta]=0x0;
					fHistdEtadPhiSameMone[ic][iz][iptt][ipta]=0x0;
					fHistdEtadPhiMixedMone[ic][iz][iptt][ipta]=0x0;
					fHistdEtadPhiSameMtwo[ic][iz][iptt][ipta]=0x0;
					fHistdEtadPhiMixedMtwo[ic][iz][iptt][ipta]=0x0;
					fHistdEtadPhiSameL[ic][iz][iptt][ipta]=0x0;
					fHistdEtadPhiMixedL[ic][iz][iptt][ipta]=0x0;
				}
			}
		}
		for(int iz=0; iz<kZvertBin;iz++) {
			fHistPtSame[ic][iz]=0x0;
			fHistPtSameH[ic][iz]=0x0;
			fHistPtSameMone[ic][iz]=0x0;
			fHistPtSameMtwo[ic][iz]=0x0;
			fHistPtSameL[ic][iz]=0x0;
			fHistPtMixed[ic][iz]=0x0;
			fHistNevtSame[ic][iz]=0x0;
			fHistNevtMixed[ic][iz]=0x0;
		}
		fHistCentBin[ic]=0x0;
	}
	for(int iz=0; iz<kZvertBin;iz++) {
		fHistZvertexBin[iz]=0x0;
	}
    for(int ipid=0; ipid<kPID; ipid++) {
        fHistPt[ipid]=0x0;
        fHistEta[ipid]=0x0;
        fHistPhi[ipid]=0x0;
    }


}
*/

//________________________________________________________________________
AliAnalysisTaskEPCorrAA::AliAnalysisTaskEPCorrAA(const char *name) // All data members should be initialised here
:AliAnalysisTaskSE(name),
	fOutput(0),
//	fTrackCuts(0),
	fPoolMgr(0x0),
	fMyprimRecoTracks(0x0),
	fTracksMixing(0x0),
	fMinNumTrack(50000),
	fPoolSize(100), // from 1000 (150610)
	fMinNEventsToMix(5),
	fNsigmaCut(3.0),
//	fHistPt(0), 
//	fHistEta(0), 
//	fHistPhi(0), 
	fHistZvertex(0),
	fPIDResponse(0),
	fHistCent(0) // The last in the above list should not have a comma after it
{
	// Constructor
	// Define input and output slots here (never in the dummy constructor)
	// Input slot #0 works with a TChain - it is connected to the default input container
	// Output slot #1 writes into a TH1 container
	for(int ic=0; ic<kCentBin;ic++) {
//		fHistNsigmaTPCpT[ic]=0x0;
//		fHistNsigmaTOFpT[ic]=0x0;
//		fHistNsigmaITSpT[ic]=0x0;
		for(int iptt=0; iptt<kpTBin;iptt++) {
//			fHistNsigmaITSTPC[ic][iptt]=0x0;
//			fHistNsigmaTPCTOF[ic][iptt]=0x0;
//			fHistNsigmaITSTOF[ic][iptt]=0x0;
			for(int ipta=0; ipta<kpTBin; ipta++) {
				for(int iz=0; iz<kZvertBin;iz++) {
					for(int ipid=0; ipid<kPID;ipid++) {
						fHistdEtadPhiSame[ic][iz][iptt][ipta][ipid]=0x0;
						fHistdEtadPhiMixed[ic][iz][iptt][ipta][ipid]=0x0;
					}
					fHistdEtadPhiSameH[ic][iz][iptt][ipta]=0x0;
					fHistdEtadPhiMixedH[ic][iz][iptt][ipta]=0x0;
					fHistdEtadPhiSameMone[ic][iz][iptt][ipta]=0x0;
					fHistdEtadPhiMixedMone[ic][iz][iptt][ipta]=0x0;
					fHistdEtadPhiSameMtwo[ic][iz][iptt][ipta]=0x0;
					fHistdEtadPhiMixedMtwo[ic][iz][iptt][ipta]=0x0;
					fHistdEtadPhiSameL[ic][iz][iptt][ipta]=0x0;
					fHistdEtadPhiMixedL[ic][iz][iptt][ipta]=0x0;

				}
			}
		}
		for(int iz=0; iz<kZvertBin;iz++) {
			fHistPtSame[ic][iz]=0x0;
			fHistPtSameH[ic][iz]=0x0;
			fHistPtSameMone[ic][iz]=0x0;
			fHistPtSameMtwo[ic][iz]=0x0;
			fHistPtSameL[ic][iz]=0x0;
			fHistPtMixed[ic][iz]=0x0;
			fHistNevtSame[ic][iz]=0x0;
			fHistNevtMixed[ic][iz]=0x0;
		}
		fHistCentBin[ic]=0x0;
	}
	for(int iz=0; iz<kZvertBin;iz++) {
		fHistZvertexBin[iz]=0x0;
	}
    for(int ipid=0; ipid<kPID; ipid++) {
        fHistPt[ipid]=0x0;
        fHistEta[ipid]=0x0;
        fHistPhi[ipid]=0x0;
    }

	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());                                            // for output list
}

//________________________________________________________________________
AliAnalysisTaskEPCorrAA::~AliAnalysisTaskEPCorrAA()
{
	// Destructor. Clean-up the output list, but not the histograms that are put inside
	// (the list is owner and will clean-up these histograms). Protect in PROOF case.
	if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
		delete fOutput;
	}
//	if(fTrackCuts) {delete fTrackCuts; fTrackCuts=0x0;}
	if(fPoolMgr) {delete fPoolMgr; fPoolMgr = 0x0;}
    if(fPIDResponse) {delete fPIDResponse; fPIDResponse=0x0;}
    if(fMyprimRecoTracks) {delete fMyprimRecoTracks; fMyprimRecoTracks=0x0;}
    if(fTracksMixing) {delete fTracksMixing; fTracksMixing=0x0;}

}

//________________________________________________________________________
void AliAnalysisTaskEPCorrAA::UserCreateOutputObjects()
{
	// Create histograms
	// Called once (on the worker node)

	fOutput = new TList();
	fOutput->SetOwner();  // IMPORTANT!
/*
	fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE, 0); // clusterCut==0 is set by default
	fTrackCuts->SetMaxDCAToVertexZ(3.2);
	fTrackCuts->SetMaxDCAToVertexXY(2.4);
	fTrackCuts->SetDCAToVertex2D(kTRUE); //  Tracks are accepted if sqrt((DCAXY / fCutMaxDCAToVertexXY)^2 + (DCAZ / fCutMaxDCAToVertexZ)^2) < 1 AND sqrt((DCAXY / fCutMinDCAToVertexXY)^2 + (DCAZ / fCutMinDCAToVertexZ)^2) > 1 (from AliESDtrackCut.h)
*/

//	fTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
//	fTrackCuts->SetMinNClustersTPC(70);

	// === Primary Track Selection ===
	//
	// The definition of a primary track is taken from the ALICE Twiki
	// page https://twiki.cern.ch/twiki/bin/view/ALICE/SelectionOfPrimaryTracksForPpDataAnalysis
	// using the following parameters for a standard dN/dPt analysis:
	//  track quality cuts:
	//          esdTrackCuts->SetMinNClustersTPC(70);
	//          esdTrackCuts->SetMaxChi2PerClusterTPC(4);
	//          esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
	//          esdTrackCuts->SetRequireTPCRefit(kTRUE);
	//          esdTrackCuts->SetRequireITSRefit(kTRUE);
	//          esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
	//                                                                                     AliESDtrackCuts::kAny);
	//  dca:
	//          if(selPrimaries) {
	//                  // 7*(0.0026+0.0050/pt^1.01)
	//                  esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
	//          }
	//          esdTrackCuts->SetMaxDCAToVertexZ(2);
	//          esdTrackCuts->SetDCAToVertex2D(kFALSE);
	//          esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
	//
	// The Primary Track Selection is implemented here by creating an
	// AliESDtrackCuts object, with kTRUE argument to choose primary tracks.
	//
	// By default, it is set to the above conditions which are suitable for
	// a standard inclusive dN/dPt analysis. For others, such as identified
	// dN/dPt or strangeness as well as others, follow the above link for
	// the specific changes to include in the selection.

	// To change cuts after selecting some default set, one can use 
	// esdtrackcuts->SetMinNClustersTPC(70) for example

    char hname[10000]; char htit[10000];

	// Create histograms
    for(int ipid=0; ipid<kPID; ipid++) {
        Int_t ptbins = 100;
        Float_t ptlow = 0.1, ptup = 20.1;
        sprintf(hname, "fHistPt%d", ipid);
        sprintf(htit, "p_{T} distribution for reconstructed PID : %d", ipid);
        fHistPt[ipid] = new TH1F(hname, htit, ptbins, ptlow, ptup);
        fHistPt[ipid]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistPt[ipid]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
        fHistPt[ipid]->SetMarkerStyle(kFullCircle);

        Int_t etabins = 40;
        Float_t etalow = -2.0, etaup = 2.0;
        sprintf(hname, "fHistEta%d", ipid);
        sprintf(htit, "#eta distribution for reconstructed PID : %d", ipid);
        fHistEta[ipid] = new TH1F(hname, htit ,etabins, etalow, etaup);
        fHistEta[ipid]->GetXaxis()->SetTitle("#eta");
        fHistEta[ipid]->GetYaxis()->SetTitle("counts");

        Int_t phibins = 160;
        Float_t philow = -1.0, phiup = 7.0;
        sprintf(hname, "fHistPhi%d", ipid);
        sprintf(htit, "#phi distribution for reconstructed PID : %d", ipid);
        fHistPhi[ipid] = new TH1F(hname, htit,phibins, philow, phiup);
        fHistPhi[ipid]->GetXaxis()->SetTitle("#phi");
        fHistPhi[ipid]->GetYaxis()->SetTitle("counts");

    }

	// For events
	int zbins = 200;
	float zlow = -10., zup = 10;
	fHistZvertex = new TH1F("fHistZvertex","Z vertex of events", zbins , zlow, zup);
	fHistZvertex->GetXaxis()->SetTitle("z-vertex (cm)");
	fHistZvertex->GetYaxis()->SetTitle("counts");
/*
	//	Float_t zBinArray[] = {-8, -4, 0, 4, 8}; // 4 bins
	float zBinArray[] = {-8, -6, -4, -2, 0, 2, 4, 6, 8}; // 8 bins
	fHistZvertexBin = new TH1F("fHistZvertexBin","Z vertex of events (zBin)",kZvertBin, zBinArray);
	fHistZvertexBin->GetXaxis()->SetTitle("z-vertex (cm)");
	fHistZvertexBin->GetYaxis()->SetTitle("counts");
*/
	Int_t centbins = 100;
	Float_t centlow = 0., centup = 100.; // by each 1%
	fHistCent = new TH1F("fHistCent","Centrality of events",centbins, centlow, centup);
	fHistCent->GetXaxis()->SetTitle("Centrality(%)");
	fHistCent->GetYaxis()->SetTitle("counts");
/*
	float cBinArray[] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}; // 10 bins
	//	float cBinArray[] = {0, 5, 10, 20, 30, 40, 50, 60, 90}; // 8 bins
	fHistCentBin = new TH1F("fHistCentBin","Centrality of events (cBin)",kCentBin, cBinArray);
	fHistCentBin->GetXaxis()->SetTitle("Centrality(%)");
	fHistCentBin->GetYaxis()->SetTitle("counts");
*/


	// Histograms for 2PC
	//	Int_t npTbins = 12, nPID = 4;
	//cout << "kZvertBin = " << kZvertBin << endl;
	// previous setting until 141017
/*
	float maxEtaRange = 1.0; // reliable eta range of TPC is -0.9 < eta < +0.9
	float maxPhiRange = 3.5; // -pi ~ +pi
	int EtaNBins = 80; // 1 bin size = 0.05 from -2.0 to 2.0
	int PhiNBins = 140; // 1 bin size = 0.10 from -7.0 to 7.0
*/
/*
	// new setting after 141017 (to set bin center at (0,0) exactly
	float maxEtaRange = 1.95 / 2.0; // reliable eta range of TPC is -0.9 < eta < +0.9
	float maxPhiRange = 6.95 / 2.0; // -pi ~ +pi
	int EtaNBins = 39; // 1 bin size = 0.10 from -1.95 to 1.95
	int PhiNBins = 139; // 1 bin size = 0.10 from -6.95 to 6.95
*/
	// new setting to reduce the memory (with preper size of bins) 150609
	float lowerEta = -1.75; float upperEta = 1.75;
	float lowerPhi = -1.75;	float upperPhi = 4.95;
	int EtaNBins = 35; // 1 bin size = 0.10 from -1.95 to 1.95
	int PhiNBins = 67; // 1 bin size = 0.10 from -6.95 to 6.95


	for(int ic=0; ic<kCentBin;ic++) {
		for(int iz=0; iz<kZvertBin;iz++) {
			for(int iptt=0; iptt<kpTBin;iptt++) {
				for(int ipta=0; ipta<kpTBin;ipta++) {
                    for(int ipid=0; ipid<kPID; ipid++) {
                        sprintf(hname, "fHistdEtadPhiSame%02d%02d%02d%02d%02d", ic, iz, iptt, ipta, ipid);
                        sprintf(htit, "Same cBin : %d zBin : %d pTbinT : %d pTbinA : %d PID : %d", ic, iz, iptt, ipta, ipid);
                        fHistdEtadPhiSame[ic][iz][iptt][ipta][ipid] = new TH2F(hname, htit, EtaNBins, lowerEta, upperEta, PhiNBins, lowerPhi, upperPhi);
                        sprintf(hname, "fHistdEtadPhiMixed%02d%02d%02d%02d%02d", ic, iz, iptt, ipta, ipid);
                        sprintf(htit, "Mixed cBin : %d zBin : %d pTbinT : %d pTbinA : %d PID : %d", ic, iz, iptt, ipta, ipid);
                        fHistdEtadPhiMixed[ic][iz][iptt][ipta][ipid] = new TH2F(hname, htit, EtaNBins, lowerEta, upperEta, PhiNBins, lowerPhi, upperPhi);
                    } // end of ipid

					sprintf(hname, "fHistdEtadPhiSameH%02d%02d%02d%02d", ic, iz, iptt, ipta);
					sprintf(htit, "Same H cBin : %d zBin : %d pTbinT : %d pTbinA : %d", ic, iz, iptt, ipta);
					fHistdEtadPhiSameH[ic][iz][iptt][ipta] = new TH2F(hname, htit, EtaNBins, lowerEta, upperEta, PhiNBins, lowerPhi, upperPhi);
					sprintf(hname, "fHistdEtadPhiMixedH%02d%02d%02d%02d", ic, iz, iptt, ipta);
					sprintf(htit, "Mixed H cBin : %d zBin : %d pTbinT : %d pTbinA : %d", ic, iz, iptt, ipta);
					fHistdEtadPhiMixedH[ic][iz][iptt][ipta] = new TH2F(hname, htit, EtaNBins, lowerEta, upperEta, PhiNBins, lowerPhi, upperPhi);
					sprintf(hname, "fHistdEtadPhiSameMone%02d%02d%02d%02d", ic, iz, iptt, ipta);
					sprintf(htit, "Same M (1) cBin : %d zBin : %d pTbinT : %d pTbinA : %d", ic, iz, iptt, ipta);
					fHistdEtadPhiSameMone[ic][iz][iptt][ipta] = new TH2F(hname, htit, EtaNBins, lowerEta, upperEta, PhiNBins, lowerPhi, upperPhi);
					sprintf(hname, "fHistdEtadPhiMixedMone%02d%02d%02d%02d", ic, iz, iptt, ipta);
					sprintf(htit, "Mixed M (1) cBin : %d zBin : %d pTbinT : %d pTbinA : %d", ic, iz, iptt, ipta);
					fHistdEtadPhiMixedMone[ic][iz][iptt][ipta] = new TH2F(hname, htit, EtaNBins, lowerEta, upperEta, PhiNBins, lowerPhi, upperPhi);
					sprintf(hname, "fHistdEtadPhiSameMtwo%02d%02d%02d%02d", ic, iz, iptt, ipta);
					sprintf(htit, "Same M (2) cBin : %d zBin : %d pTbinT : %d pTbinA : %d", ic, iz, iptt, ipta);
					fHistdEtadPhiSameMtwo[ic][iz][iptt][ipta] = new TH2F(hname, htit, EtaNBins, lowerEta, upperEta, PhiNBins, lowerPhi, upperPhi);
					sprintf(hname, "fHistdEtadPhiMixedMtwo%02d%02d%02d%02d", ic, iz, iptt, ipta);
					sprintf(htit, "Mixed M (2) cBin : %d zBin : %d pTbinT : %d pTbinA : %d", ic, iz, iptt, ipta);
					fHistdEtadPhiMixedMtwo[ic][iz][iptt][ipta] = new TH2F(hname, htit, EtaNBins, lowerEta, upperEta, PhiNBins, lowerPhi, upperPhi);

					sprintf(hname, "fHistdEtadPhiSameL%02d%02d%02d%02d", ic, iz, iptt, ipta);
					sprintf(htit, "Same L cBin : %d zBin : %d pTbinT : %d pTbinA : %d", ic, iz, iptt, ipta);
					fHistdEtadPhiSameL[ic][iz][iptt][ipta] = new TH2F(hname, htit, EtaNBins, lowerEta, upperEta, PhiNBins, lowerPhi, upperPhi);
					sprintf(hname, "fHistdEtadPhiMixedL%02d%02d%02d%02d", ic, iz, iptt, ipta);
					sprintf(htit, "Mixed L cBin : %d zBin : %d pTbinT : %d pTbinA : %d", ic, iz, iptt, ipta);
					fHistdEtadPhiMixedL[ic][iz][iptt][ipta] = new TH2F(hname, htit, EtaNBins, lowerEta, upperEta, PhiNBins, lowerPhi, upperPhi);

				} // end of ipta
			} // end of iptt

			Int_t ptbins = 100;
			Float_t ptlow = 0.1, ptup = 20.1;

			sprintf(hname, "fHistPtSame%02d%02d", ic, iz);
			sprintf(htit, "pT, same cBin : %d zBin : %d", ic, iz);
			fHistPtSame[ic][iz] = new TH1F(hname, htit, ptbins, ptlow, ptup);
			fHistPtSame[ic][iz]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
			fHistPtSame[ic][iz]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
			fHistPtSame[ic][iz]->SetMarkerStyle(kFullCircle);

			sprintf(hname, "fHistPtSameH%02d%02d", ic, iz);
			sprintf(htit, "pT, same H cBin : %d zBin : %d", ic, iz);
			fHistPtSameH[ic][iz] = new TH1F(hname, htit, ptbins, ptlow, ptup);
			fHistPtSameH[ic][iz]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
			fHistPtSameH[ic][iz]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
			fHistPtSameH[ic][iz]->SetMarkerStyle(kFullCircle);

			sprintf(hname, "fHistPtSameMone%02d%02d", ic, iz);
			sprintf(htit, "pT, same M cBin : %d zBin : %d", ic, iz);
			fHistPtSameMone[ic][iz] = new TH1F(hname, htit, ptbins, ptlow, ptup);
			fHistPtSameMone[ic][iz]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
			fHistPtSameMone[ic][iz]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
			fHistPtSameMone[ic][iz]->SetMarkerStyle(kFullCircle);

			sprintf(hname, "fHistPtSameMtwo%02d%02d", ic, iz);
			sprintf(htit, "pT, same M cBin : %d zBin : %d", ic, iz);
			fHistPtSameMtwo[ic][iz] = new TH1F(hname, htit, ptbins, ptlow, ptup);
			fHistPtSameMtwo[ic][iz]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
			fHistPtSameMtwo[ic][iz]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
			fHistPtSameMtwo[ic][iz]->SetMarkerStyle(kFullCircle);

			sprintf(hname, "fHistPtSameL%02d%02d", ic, iz);
			sprintf(htit, "pT, same L cBin : %d zBin : %d", ic, iz);
			fHistPtSameL[ic][iz] = new TH1F(hname, htit, ptbins, ptlow, ptup);
			fHistPtSameL[ic][iz]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
			fHistPtSameL[ic][iz]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
			fHistPtSameL[ic][iz]->SetMarkerStyle(kFullCircle);

			sprintf(hname, "fHistPtMixed%02d%02d", ic, iz);
			sprintf(htit, "pT, mixed cBin : %d zBin : %d", ic, iz);
			fHistPtMixed[ic][iz] = new TH1F(hname, htit, ptbins, ptlow, ptup);
			fHistPtMixed[ic][iz]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
			fHistPtMixed[ic][iz]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
			fHistPtMixed[ic][iz]->SetMarkerStyle(kFullCircle);

			// Count for normalization
			int czmixbin = 80; // cBin 10 * zBin 8
			Float_t czlow = -0.5, czup = 79.5;

			sprintf(hname, "fHistNevtSame%02d%02d", ic, iz);
			sprintf(htit, "Nevt,same %02d%02d (zBin*10+cBin)", ic, iz);
			fHistNevtSame[ic][iz] = new TH1F(hname, htit, czmixbin, czlow, czup); 
			fHistNevtSame[ic][iz]->GetXaxis()->SetTitle("zBin*10 + cBin");
			fHistNevtSame[ic][iz]->GetXaxis()->SetTitle("N_{evt, same}");

			sprintf(hname, "fHistNevtMixed%02d%02d", ic, iz);
			sprintf(htit, "Nevt,mixed %02d%02d (zBin*10+cBin)", ic, iz);
			fHistNevtMixed[ic][iz] = new TH1F(hname, htit, czmixbin, czlow, czup); 
			fHistNevtMixed[ic][iz]->GetXaxis()->SetTitle("zBin*10 + cBin");
			fHistNevtMixed[ic][iz]->GetXaxis()->SetTitle("N_{evt, mixed}");

		} // end of iz

		fHistCentBin[ic] = new TH1F(Form("fHistCentBin%d", ic),Form("Centrality of events (cBin=%d)", ic),kCentBin, -0.5,9.5);
		fHistCentBin[ic]->GetXaxis()->SetTitle("CentBin");
		fHistCentBin[ic]->GetYaxis()->SetTitle("counts");

		fHistEPQv[ic] = new TH1F(Form("fHistEPQv%d", ic),Form("Event Plane angle (cBin=%d) with V0A and V0C", ic), 100, -5,5);
		fHistEPQv[ic]->GetXaxis()->SetTitle("#phi_{EP}");
		fHistEPQv[ic]->GetYaxis()->SetTitle("counts");

		fHistEPQvVzeroA[ic] = new TH1F(Form("fHistEPQvVzeroA%d", ic),Form("Event Plane angle (cBin=%d) with V0A", ic), 100, -5,5);
		fHistEPQvVzeroA[ic]->GetXaxis()->SetTitle("#phi_{EP}");
		fHistEPQvVzeroA[ic]->GetYaxis()->SetTitle("counts");

		fHistEPQvVzeroC[ic] = new TH1F(Form("fHistEPQvVzeroC%d", ic),Form("Event Plane angle (cBin=%d) with V0C", ic), 100, -5,5);
		fHistEPQvVzeroC[ic]->GetXaxis()->SetTitle("#phi_{EP}");
		fHistEPQvVzeroC[ic]->GetYaxis()->SetTitle("counts");

		fHistEPQvTPC[ic] = new TH1F(Form("fHistEPQvTPC%d", ic),Form("Event Plane angle (cBin=%d) with TPC", ic), 100, -5,5);
		fHistEPQvTPC[ic]->GetXaxis()->SetTitle("#phi_{EP}");
		fHistEPQvTPC[ic]->GetYaxis()->SetTitle("counts");


		fHistEPQvAandC[ic] = new TH2F(Form("fHistEPQvAandC%d", ic),Form("Event Plane angle (cBin=%d) A vs. C", ic), 100, -5,5, 100, -5, 5);
		fHistEPQvAandC[ic]->GetXaxis()->SetTitle("EP by V0A");
		fHistEPQvAandC[ic]->GetYaxis()->SetTitle("EP by V0C");


		fHistEPQvAandAC[ic] = new TH2F(Form("fHistEPQvAandAC%d", ic),Form("Event Plane angle (cBin=%d) A vs. A-C", ic), 100, -5,5, 100, -5, 5);
		fHistEPQvAandAC[ic]->GetXaxis()->SetTitle("EP by V0A");
		fHistEPQvAandAC[ic]->GetYaxis()->SetTitle("EP by V0A and V0C");

		fHistEPQvCandAC[ic] = new TH2F(Form("fHistEPQvCandAC%d", ic),Form("Event Plane angle (cBin=%d) C vs. A-C", ic), 100, -5,5, 100, -5, 5);
		fHistEPQvCandAC[ic]->GetXaxis()->SetTitle("EP by V0C");
		fHistEPQvCandAC[ic]->GetYaxis()->SetTitle("EP by V0C and V0C");


	} // end of ic

	for(int iz=0; iz<kZvertBin;iz++) {
		fHistZvertexBin[iz] = new TH1F(Form("fHistZvertexBin%d",iz),Form("Z vertex of events(zBin) %d",iz),kZvertBin, -0.5, 7.5);
		fHistZvertexBin[iz]->GetXaxis()->SetTitle("zvertexBin");
		fHistZvertexBin[iz]->GetYaxis()->SetTitle("counts");
	}

    // Add histograms into Output
    for(int ipid=0; ipid<kPID; ipid++) {
        fOutput->Add(fHistPt[ipid]);
        fOutput->Add(fHistEta[ipid]);
        fOutput->Add(fHistPhi[ipid]);
    }


	fOutput->Add(fHistZvertex);

	for(int iz=0; iz<kZvertBin;iz++) {
		fOutput->Add(fHistZvertexBin[iz]);
	}

	fOutput->Add(fHistCent);
	for(int ic=0; ic<kCentBin;ic++) {
		fOutput->Add(fHistCentBin[ic]);
	}

	for(int ic=0; ic<kCentBin;ic++) {
		fOutput->Add(fHistEPQv[ic]);
		fOutput->Add(fHistEPQvVzeroA[ic]);
		fOutput->Add(fHistEPQvVzeroC[ic]);
		fOutput->Add(fHistEPQvTPC[ic]);
		fOutput->Add(fHistEPQvAandC[ic]);
		fOutput->Add(fHistEPQvAandAC[ic]);
		fOutput->Add(fHistEPQvCandAC[ic]);

	}


	for(int ic=0; ic<kCentBin;ic++) {
		for(int iz=0; iz<kZvertBin;iz++) {
			fOutput->Add(fHistNevtSame[ic][iz]);
			fOutput->Add(fHistNevtMixed[ic][iz]);
			fOutput->Add(fHistPtSame[ic][iz]);
			fOutput->Add(fHistPtSameH[ic][iz]);
			fOutput->Add(fHistPtSameMone[ic][iz]);
			fOutput->Add(fHistPtSameMtwo[ic][iz]);
			fOutput->Add(fHistPtSameL[ic][iz]);
			fOutput->Add(fHistPtMixed[ic][iz]);

		}
	}
	for(int ic=0; ic<kCentBin;ic++) {
		for(int iz=0; iz<kZvertBin;iz++) {
			for(int iptt=0; iptt<kpTBin;iptt++) {
				for(int ipta=0; ipta<kpTBin;ipta++) {
					//					cout << "ic = " << ic << " | iz = " << iz << " | iptt = " << iptt << " | ipta = " << ipta << endl;
                    for(int ipid=0; ipid<kPID; ipid++) {
						if(iptt >= ipta) fOutput->Add(fHistdEtadPhiSame[ic][iz][iptt][ipta][ipid]);
						if(iptt >= ipta) fOutput->Add(fHistdEtadPhiMixed[ic][iz][iptt][ipta][ipid]);
					}
					if(iptt >= ipta) fOutput->Add(fHistdEtadPhiSameH[ic][iz][iptt][ipta]);
					if(iptt >= ipta) fOutput->Add(fHistdEtadPhiMixedH[ic][iz][iptt][ipta]);
					if(iptt >= ipta) fOutput->Add(fHistdEtadPhiSameMone[ic][iz][iptt][ipta]);
					if(iptt >= ipta) fOutput->Add(fHistdEtadPhiMixedMone[ic][iz][iptt][ipta]);
					if(iptt >= ipta) fOutput->Add(fHistdEtadPhiSameMtwo[ic][iz][iptt][ipta]);
					if(iptt >= ipta) fOutput->Add(fHistdEtadPhiMixedMtwo[ic][iz][iptt][ipta]);
					if(iptt >= ipta) fOutput->Add(fHistdEtadPhiSameL[ic][iz][iptt][ipta]);
					if(iptt >= ipta) fOutput->Add(fHistdEtadPhiMixedL[ic][iz][iptt][ipta]);


				} // end of ipta
			} // end of iptt
		} // end of iz
	} // end of ic


    // call PIDResponse for PID
    AliAnalysisManager *anamanager = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*)anamanager->GetInputEventHandler();
    fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();


	// call PoolManager setup for mixing
	SetupForMixing();

	PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
	//cout << "DEBUG0" << endl;

}



//________________________________________________________________________
void AliAnalysisTaskEPCorrAA::UserExec(Option_t *) 
{
	// Main loop
	// Called for each event

	// Create pointer to reconstructed event
	AliVEvent *event = InputEvent();
	if (!event) { Printf("ERROR: Could not retrieve event"); return; }

	// If the task accesses MC info, this can be done as in the commented block below:
	/*
	// Create pointer to reconstructed event
	AliMCEvent *mcEvent = MCEvent();
	if (!mcEvent) { Printf("ERROR: Could not retrieve MC event"); return; }
	Printf("MC particles: %d", mcEvent->GetNumberOfTracks());

	// set up a stack for use in check for primary/stable particles
	AliStack* stack = mcEvent->Stack();
	if( !stack ) { Printf( "Stack not available"); return; }
	 */  
	//cout << "DEBUG1" << endl;

	// create pointer to event
	AliAODEvent* aodevent = dynamic_cast<AliAODEvent*>(event);
	if (!aodevent) {
		AliError("Cannot get the AOD event");
		return;
	}  
	//if (!TGeoGlobalMagField::Instance()->GetField()) aodevent->InitMagneticField(); 

	AliAODHeader* aodheader = (AliAODHeader*)aodevent->GetHeader();
	//cout << "DEBUG2" << endl;


	// === Physics Selection Task ===
	// 
	// To perform a physics selection here, a bitwise operation is used against
	// the UInt_t mask which is extracted in the following way:
	//
	  UInt_t mask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();   
	//
	// This can be tested to produce the following
	//
	Bool_t bTargetEvents = ( (mask==AliVEvent::kMB) || (mask==AliVEvent::kCentral) || (mask==AliVEvent::kSemiCentral) ) ? 1 : 0; // check if minimum bias trigger class fired
	


	//  Bool_t bHighMult = (mask == AliVEvent::kHighMult) ? 1 : 0; // check if high multiplicity trigger class fired
	//
	// For more complicated trigger selections, one can directly test both
	// trigger classes and fired trigger inputs for a particular event, for e.g.
	//
	//  Bool_t bCSH1 = (esd->IsTriggerClassFired("CSH1-B-NOPF-ALLNOTRD")) ? 1 : 0;
	//  Bool_t b0SH1 = (esdheader->IsTriggerInputFired("0SH1")) ? 1 : 0;
	//
	// These booleans can then be used to fill different histograms for specific
	// conditions, or summed to make one cut for all events that fill histos.

	if (bTargetEvents==0) return;


	// Do some fast cuts first
	// check for good reconstructed vertex
	if(!(aodevent->GetPrimaryVertex())) return;
	// if vertex is from spd vertexZ, require more stringent cut
//	if (esdevent->GetPrimaryVertex()->IsFromVertexerZ()) {
//		if (esdevent->GetPrimaryVertex()->GetDispersion()>0.02 ||  esdevent->GetPrimaryVertex()->GetZRes()>0.25 ) return; // bad vertex from VertexerZ
//	}

	Float_t bSign = 0;
	bSign = (event->GetMagneticField() > 0) ? 1 : -1;

	Float_t cent = -99;
	cent = aodevent->GetCentrality()->GetCentralityPercentile("V0M");
	fHistCent->Fill(cent);

//	float cBinArray[] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}; // 10 bins
//	float cBinArray[] = {0, 5, 10, 20, 30, 40, 50, 60, 90}; // 8 bins
	float cBinArray[] = {0, 10, 20, 50, 60, 90}; // 5 bins for LHC11h (non flat cent) 
//	float cBinArray[] = {0, 10, 20, 35, 50, 60, 90}; // 6 bins for LHC11h (non flat cent + closer look in mid-cent) 150519

//	float cBinArray[] = {0, 10, 30, 60, 90}; // 4 bins for LHC10h 

	//	float cBinArray[] = {0, 5, 10, 20, 30, 40, 50, 70, 90}; // 8 bins

	// Get cBin from percentage
	int cBin = -1;
	for(int icbin = 0; icbin<kCentBin; icbin++) {
		if(0 > cent) continue;
		float cBinMin = cBinArray[icbin];
		float cBinMax = cBinArray[icbin+1];
		if((cBinMin <= cent) && (cent < cBinMax)) cBin = icbin;	
	}
	if(!(cBin == -1)) fHistCentBin[cBin]->Fill(cBin); // to fill from the 1st bin


	// fill event plane angle
	double EPangle = -10; double EPangleVzeroA = -10; double EPangleVzeroC = -10; double EPangleTPC = -10;

	AliEventplane *evp = aodevent->GetEventplane();
	if(evp) {
		EPangle = evp->GetEventplane("V0",aodevent, 2); // calculated by whold V0 for second harmonics
		EPangleVzeroA = evp->GetEventplane("V0A",aodevent, 2);
	 	EPangleVzeroC = evp->GetEventplane("V0C",aodevent, 2);
		EPangleTPC = evp->GetEventplane("Q"); // default fEventplaneQ

		if(EPangle<0) EPangle = EPangle + TMath::Pi(); // V0 EP angle varies from -pi/2 to +pi/2, so need to shift
		if(EPangleVzeroA<0) EPangleVzeroA = EPangleVzeroA + TMath::Pi(); 
		if(EPangleVzeroC<0) EPangleVzeroC = EPangleVzeroC + TMath::Pi(); 
		if(EPangleTPC<0) EPangleTPC = EPangleTPC + TMath::Pi(); 
	}
	if(!(cBin == -1)) {
		fHistEPQv[cBin]->Fill(EPangle); 
		fHistEPQvVzeroA[cBin]->Fill(EPangleVzeroA); 
		fHistEPQvVzeroC[cBin]->Fill(EPangleVzeroC); 
		fHistEPQvTPC[cBin]->Fill(EPangleTPC); 

		fHistEPQvAandC[cBin]->Fill(EPangleVzeroA, EPangleVzeroC); 
		fHistEPQvAandAC[cBin]->Fill(EPangleVzeroA, EPangle); 
		fHistEPQvCandAC[cBin]->Fill(EPangleVzeroC, EPangle); 


//	cout << "EP angle = " << EPangle << endl;
	}


	//cout << "DEBUG3" << endl;
	// Fill zVertex and zBin
//	float zBinArray[] = {-8, -6, -4, -2, 0, 2, 4, 6, 8}; // 8 bins
	float zBinArray[] = {-7, -5, -3, -1, 1, 3, 5, 7}; // 7 bins

	//	float zBinArray[] = {-8, -4, 0, 4, 8}; // 4 bins
	float zvertex = -999;
	zvertex = aodevent->GetPrimaryVertex()->GetZ();

	// select the events within |vtx-z| < 7 cm 
//	if( TMath::Abs(zvertex) > 7.0 ) return;
	if( TMath::Abs(zvertex) > 3.0 ) return; // for systematic study


	int zBin = -1;
	for(int izbin = 0; izbin<kZvertBin; izbin++) {
		if(zvertex == -999.) continue; // exclude what have initial value
		float zBinMin = zBinArray[izbin];
		float zBinMax = zBinArray[izbin+1];
		if((zBinMin <= zvertex) && (zvertex < zBinMax)) zBin = izbin; // to fill from the 1st bin
	}
	fHistZvertex->Fill(zvertex);
	if(!(zBin == -1)) fHistZvertexBin[zBin]->Fill(zBin); // to fill from the 1st bin

//	if(!(zBin == -1 || cBin == -1)) fHistNevtSame->AddBinContent(zBin*10 + cBin + 1); // 8 zBins and 10 cBin will fill 1~80th bins
	if(!(zBin == -1 || cBin == -1)) fHistNevtSame[cBin][zBin]->Fill(zBin*10 + cBin); // 8 zBins and 10 cBin will be filled

	//	cout << "zvertex = " << zvertex << " zBin = " << zBin << endl;



//	const double nsigmaCut = 3.0;






	//cout << "DEBUG4" << endl;
	// ===== Track loop for reconstructed event in SAME events =====
	Int_t ntracks = aodevent->GetNumberOfTracks();
	//cout << "ntracks = " << ntracks << endl;
	for(Int_t i = 0; i < ntracks; i++) {

		//cout << "DEBUG5" << endl;
		AliAODTrack* aodtrack = (AliAODTrack*)aodevent->GetTrack(i); // pointer to reconstructed to track          
		if(!aodtrack) { 
			AliError(Form("ERROR: Could not retrieve aodtrack %d",i)); 
			continue; 
		}

		//cout << "DEBUG6" << endl;
		// Some MC checks, if MC is used
		if(aodtrack->GetLabel() < 0) continue; // get rid of "ghost" tracks

		// ... and the thorough checking of ESD cuts after.
		// if this is not a primary track, skip to the next one
//		if(!fTrackCuts->AcceptTrack(aodtrack)) continue;
        if(aodtrack->GetType() != AliAODTrack::kPrimary) continue;
//        if(!aodtrack->TestFilterBit(kHybridTrackCut)) continue;
        if(!aodtrack->TestFilterBit(kTPCOnlyTrackCut)) continue;

		float trackpT = aodtrack->Pt();
		float trackEta = aodtrack->Eta(); 
		float trackPhi = aodtrack->Phi();

//		float trackZv = aodtrack->GetZ();
		float trackZv = aodtrack->Charge();

		if( TMath::Abs(trackEta)>0.8 ) continue;

		//		if(i == 4) cout << "trackEta of the 4th track in same event loop = " << trackEta << endl;

		// find corresponding pT bin for many uses
//		float pTBinArray[] = { 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0 }; // 10 bins
//		float pTBinArray[] = { 0.5, 1.0, 2.0, 3.0, 4.0, 8.0, 15.0 }; // 6 bins
		float pTBinArray[] = { 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 15.0 }; // 7 bins
//		float pTBinArray[] = { 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 6.0, 8.0, 15.0 }; // 10 bins

		int pTBinT = -1; // pT bin of triggered particles
		for(int ipTbin = 0; ipTbin<kpTBin; ipTbin++) {
			if(0 > trackpT) continue;
			float pTBinMin = pTBinArray[ipTbin];
			float pTBinMax = pTBinArray[ipTbin+1];
			if((pTBinMin <= trackpT) && (trackpT < pTBinMax)) pTBinT = ipTbin; // to fill from the 1st bin
		}
		//		cout << "trackpT = " << trackpT << " | pTBinT = " << pTBinT << endl;
/*
        int trigID = -1; // 0(hadron) 1(pion) 2(kaon) 3(proton)
        trigID = GetParticleID(aodtrack);
        int trigIDBin = -1;
        if(trigID==1 || trigID==3 || trigID==5 || trigID==7) trigIDBin = 1; // recon. pion
        if(trigID==2 || trigID==3 || trigID==6 || trigID==7) trigIDBin = 2; // recon. kaon
        if(trigID==4 || trigID==5 || trigID==6 || trigID==7) trigIDBin = 3; // recon. proton
*/

        fHistPt[0]->Fill(trackpT); // fill inclusive spectra
        fHistEta[0]->Fill(trackEta);
        fHistPhi[0]->Fill(trackPhi);
/*
        if(!(trigIDBin == -1)) {
          fHistPt[trigIDBin]->Fill(trackpT);
            fHistEta[trigIDBin]->Fill(trackEta);
            fHistPhi[trigIDBin]->Fill(trackPhi);
        }
*/

		// Classify by the angle between trigger and EP
/*
		if(EPangleVzeroA == -10) continue; // exclude default value
		double diffTrigEP = -10;
		diffTrigEP = trackPhi - EPangleVzeroA;
*/
		if(EPangle == -10) continue; // exclude default value
		double diffTrigEP = -10;
		diffTrigEP = trackPhi - EPangle;



		if(diffTrigEP < 0) diffTrigEP = diffTrigEP + 2*TMath::Pi();

		if(!(zBin == -1 || cBin == -1)) {
			fHistPtSame[cBin][zBin]->Fill(trackpT);

			// H : near EP axis, M : intermidiate (one & two seperately to distinguish the shifts), L : perpendicular to EP
			if( (0 < diffTrigEP) && (diffTrigEP < (1./6.)*TMath::Pi()) ) fHistPtSameH[cBin][zBin]->Fill(trackpT);
			if( ((5./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (7./6.)*TMath::Pi()) ) fHistPtSameH[cBin][zBin]->Fill(trackpT);
			if( ((11./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (12./6.)*TMath::Pi()) ) fHistPtSameH[cBin][zBin]->Fill(trackpT);

			if( ((1./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (2./6.)*TMath::Pi()) ) fHistPtSameMone[cBin][zBin]->Fill(trackpT);
			if( ((7./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (8./6.)*TMath::Pi()) ) fHistPtSameMone[cBin][zBin]->Fill(trackpT);

			if( ((4./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (5./6.)*TMath::Pi()) ) fHistPtSameMtwo[cBin][zBin]->Fill(trackpT);
			if( ((10./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (11./6.)*TMath::Pi()) ) fHistPtSameMtwo[cBin][zBin]->Fill(trackpT);

			if( ((2./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (4./6.)*TMath::Pi()) ) fHistPtSameL[cBin][zBin]->Fill(trackpT);
			if( ((8./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (10./6.)*TMath::Pi()) ) fHistPtSameL[cBin][zBin]->Fill(trackpT);

		}

		// ===== another track loop to use dEta-dPhi analysis in same event =====
//		for(int j = i+1; j < ntracks; j++) {
		for(int j = 0; j < ntracks; j++) {

			if(i==j) continue;

			AliAODTrack* aodtrack2 = (AliAODTrack*)aodevent->GetTrack(j); // pointer to reconstructed to track          
			if(!aodtrack2) { 
				AliError(Form("ERROR: Could not retrieve aodtrack %d",j)); 
				continue; 
			}

//			if(!fTrackCuts->AcceptTrack(aodtrack2)) continue;
            if(aodtrack2->GetType() != AliAODTrack::kPrimary) continue;
//            if(!aodtrack2->TestFilterBit(kHybridTrackCut)) continue;
  		    if(!aodtrack2->TestFilterBit(kTPCOnlyTrackCut)) continue;

			float trackpT2 = aodtrack2->Pt();
			float trackEta2 = aodtrack2->Eta(); 
			float trackPhi2 = aodtrack2->Phi();

//			float trackZv2 = aodtrack2->GetZ();
			float trackZv2 = aodtrack2->Charge();

			int pTBinA = -1; // pT bin of associated particles
			for(int ipTbin = 0; ipTbin<kpTBin; ipTbin++) {
				if(0 > trackpT2) continue;
				float pTBinMin = pTBinArray[ipTbin];
				float pTBinMax = pTBinArray[ipTbin+1];
				if((pTBinMin <= trackpT2) && (trackpT2 < pTBinMax)) pTBinA = ipTbin; // to fill from the 1st bin
			}

			if( cBin == -1 || zBin == -1 || pTBinT == -1 || pTBinA == -1 ) continue; // exclude what have initial values



			if( trackpT < trackpT2 ) continue; // select pT_trigg >= pT_assoc

			if( TMath::Abs(trackEta2)>0.8 ) continue;

/*
			// instead of excluding tracks, switch the information of pTt and pTa (dEta, dPhi, bins and diffTrigEP) 
			if( trackpT < trackpT2 ) {
				double bufferpTBin, bufferEta, bufferPhi;
				bufferpTBin = pTBinT;
				pTBinT = pTBinA;
				pTBinA = bufferpTBin;

				bufferEta = trackEta;
				trackEta = trackEta2;
				trackEta2 = bufferEta;

				bufferPhi = trackPhi;
				trackPhi = trackPhi2;
				trackPhi2 = bufferPhi;

				diffTrigEP = trackPhi - EPangle; // reset phiS with 'second' track after the switch
			}
*/

			float deltaEta, deltaPhi; // calculate the difference of two tracks
			deltaEta = trackEta - trackEta2;
			deltaPhi = trackPhi - trackPhi2;

			if(deltaPhi < -0.5*TMath::Pi()) deltaPhi += TMath::TwoPi();
			if(deltaPhi > 1.5*TMath::Pi()) deltaPhi -= TMath::TwoPi();


///* INACTIVATED 140519
			float dPhiStarSame = -99; // for two-track efficiency cut
			dPhiStarSame = CalculatedPhiStar(deltaPhi, deltaEta, trackZv, trackZv2, trackpT, trackpT2, bSign);

//			cout << "dPhiStarSame = " << dPhiStarSame << endl;

			// two-track efficiency cut
//			if(TMath::Abs(dPhiStarSame)<0.02) continue;
//			if(TMath::Abs(deltaEta)<0.02) continue;
            if(TMath::Abs(dPhiStarSame)<0.02 && TMath::Abs(deltaEta)<0.02 ) continue; // should apply 'AND' here! 150901

//*/

/*
            int assocID = -1; // 0(hadron) 1(pion) 2(kaon) 3(proton)
            assocID = GetParticleID(aodtrack2);

            int assocIDBin = -1;
            if(assocID==1 || assocID==3 || assocID==5 || assocID==7) assocIDBin = 1; // recon. pion
            if(assocID==2 || assocID==3 || assocID==6 || assocID==7) assocIDBin = 2; // recon. kaon
 	        if(assocID==4 || assocID==5 || assocID==6 || assocID==7) assocIDBin = 3; // recon. proton
*/

//			cout << "assocID = " << assocID << endl;

            fHistdEtadPhiSame[cBin][zBin][pTBinT][pTBinA][0]->Fill(deltaEta, deltaPhi);
/*
            if(!(assocIDBin == -1)) {
                fHistdEtadPhiSame[cBin][zBin][pTBinT][pTBinA][assocIDBin]->Fill(deltaEta, deltaPhi);
            //          cout << "deltaPhi = " << deltaPhi << " cBin = " << cBin << " zBin = " << zBin << " pTBinT = " << pTBinT << " pTBinA = " << pTBinA << endl;
            }
*/

			
			// H : near EP axis, M : intermidiate, L : perpendicular to EP
			if( (0 < diffTrigEP) && (diffTrigEP < (1./6.)*TMath::Pi()) ) fHistdEtadPhiSameH[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEta, deltaPhi);
			if( ((5./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (7./6.)*TMath::Pi()) ) fHistdEtadPhiSameH[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEta, deltaPhi);
			if( ((11./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (12./6.)*TMath::Pi()) ) fHistdEtadPhiSameH[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEta, deltaPhi);

			if( ((1./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (2./6.)*TMath::Pi()) ) fHistdEtadPhiSameMone[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEta, deltaPhi);
			if( ((7./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (8./6.)*TMath::Pi()) ) fHistdEtadPhiSameMone[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEta, deltaPhi);
			if( ((4./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (5./6.)*TMath::Pi()) ) fHistdEtadPhiSameMtwo[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEta, deltaPhi);
			if( ((10./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (11./6.)*TMath::Pi()) ) fHistdEtadPhiSameMtwo[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEta, deltaPhi);

			if( ((2./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (4./6.)*TMath::Pi()) ) fHistdEtadPhiSameL[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEta, deltaPhi);
			if( ((8./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (10./6.)*TMath::Pi()) ) fHistdEtadPhiSameL[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEta, deltaPhi);


			//			cout << "deltaPhi = " << deltaPhi << " cBin = " << cBin << " zBin = " << zBin << " pTBinT = " << pTBinT << " pTBinA = " << pTBinA << endl;
		} // end of j track loop
	} // end of i track loop

	//	cout << "===== Prepare event-mixing =====" << endl;

	bool kTrackCut = kTRUE; // currently basic track cut is apllied 150429

	fMyprimRecoTracks = AcceptTracksReduced(aodevent, kTrackCut);

	if(fMyprimRecoTracks==0x0) { AliInfo(" ==== fMyprimRecoTracks: Zero track pointer"); return; }

	if(fMyprimRecoTracks->GetEntriesFast() == 0) {
		AliInfo(Form("========== Rejecting event because it has no tracks: %f %d", cent, fMyprimRecoTracks->GetEntriesFast()));
		return;
	}

	bool doMixing; doMixing = kTRUE;
//	if(doMixing) DoMixing(cent, zvertex, EPangleVzeroA, fMyprimRecoTracks);
	if(doMixing) DoMixing(cent, zvertex, EPangle, fMyprimRecoTracks, bSign);

	// NEW HISTO should be filled before this point, as PostData puts the
	// information for this iteration of the UserExec in the container

	PostData(1, fOutput);
	//	cout << "DEBUG12" << endl;

}








//________________________________________________________________________
void AliAnalysisTaskEPCorrAA::SetupForMixing()
{
	const Int_t trackDepth  = fMinNumTrack;
	const Int_t poolsize        = fPoolSize;
//	Double_t centralityBins[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90 }; // 10 bins
//	Double_t centralityBins[] = { 0, 5, 10, 20, 30, 40, 50, 60, 90 }; // 8 bins
	Double_t centralityBins[] = { 0, 10, 20, 50, 60, 90}; // 5 bins for LHC11h (non flat cent) 
//	Double_t centralityBins[] = { 0, 10, 20, 35, 50, 60, 90}; // 6 bins for LHC11h (non flat cent + closer look in mid-cent) 

//	Double_t zBinArray[] = {-8, -6, -4, -2, 0, 2, 4, 6, 8}; // 8 bins
	Double_t zBinArray[] = {-7, -5, -3, -1, 1, 3, 5, 7}; // 8 bins

	//    const Int_t nCentralityBins = length(centralityBins)-1;
	const Int_t nCentralityBins = ( sizeof(centralityBins)/sizeof(*centralityBins) ) -1;
	const Int_t nzBins = ( sizeof(zBinArray)/sizeof(*zBinArray) ) -1;


	fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins/*fCentAxis->GetNbins()+1*/, centralityBins/*const_cast<Double_t*>(fCentAxis->GetXbins()->GetArray())*/, nzBins /*8 bins + 1*/, zBinArray /*const_cast<Double_t*>(fZvtxAxis->GetXbins()->GetArray())*/);

	//	fPoolMgr->SetTargetValues(trackDepth, 0.1, 5);

	fPoolMgr->SetDebug(0);

//	cout << "nCentralityBins = " << nCentralityBins << endl;

}

//________________________________________________________________________
void AliAnalysisTaskEPCorrAA::FillMixedHistos(TObjArray* partNew, TObjArray* partMix, int cBin, int zBin, double EPangle, float bSignHistos, double weight)
{
//	TObjArray *particle = new TObjArray;
//	TObjArray *particleMixed = new TObjArray;

	partNew->SetOwner(kTRUE);
	partMix->SetOwner(kTRUE);

//	particle = (TObjArray*)partNew->Clone();
//	particleMixed = (TObjArray*)partMix->Clone();

    TObjArray *particle = (TObjArray*)partNew->Clone();
    TObjArray *particleMixed = (TObjArray*)partMix->Clone();

	int ntrackFirst = particle->GetEntriesFast();
	double useWeight = weight;

	//	for(int targetindex = 0; targetindex < 5; targetindex++) { // pick a few tracks in the first event
	for(int targetindex = 0; targetindex < ntrackFirst; targetindex++) { // pick a few tracks in the first event

		//	int targetindex = 3; // targeted index of track in the 1st event for mixing

		if(!(particle->At(targetindex))) return;

//		if(!fTrackCuts->AcceptTrack( particle->At(targetindex) )) continue;


		float firstpT = -999; float firstEta = -999; float firstPhi = -999; float firstZv = -999;
		firstpT = ((AliCorrReducedTrackAA*)particle->At(targetindex))->Pt();
		firstEta = ((AliCorrReducedTrackAA*)particle->At(targetindex))->Eta();
		firstPhi = ((AliCorrReducedTrackAA*)particle->At(targetindex))->Phi();
//		firstZv = ((AliCorrReducedTrackAA*)particle->At(targetindex))->Zv(); // not implemented in Reduced track....
        firstZv = ((AliCorrReducedTrackAA*)particle->At(targetindex))->Charge(); 

//        if(((AliAODTrack*)particle->At(targetindex))->TestFilterBit(kTPCOnlyTrackCut)) cout << "first track is TPCOnly" << endl;;
//		cout << ((AliAODTrack*)particle->At(targetindex))->TestFilterBit(kTPCOnlyTrackCut) << endl;

		if( TMath::Abs(firstEta)>0.8 ) continue;


		// find pT bin of first track
//		float pTBinArray[] = { 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0 }; // 10 bins
//		float pTBinArray[] = { 0.5, 1.0, 2.0, 3.0, 4.0, 8.0, 15.0 }; // 6 bins
		float pTBinArray[] = { 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 15.0 }; // 7 bins
//		float pTBinArray[] = { 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 6.0, 8.0, 15.0 }; // 10 bins

		int pTBinT = -1; // pT bin of triggered particles
		for(int ipTbin = 0; ipTbin<kpTBin; ipTbin++) {
			if(0 > firstpT) continue;
			float pTBinMin = pTBinArray[ipTbin];
			float pTBinMax = pTBinArray[ipTbin+1];
			if((pTBinMin <= firstpT) && (firstpT < pTBinMax)) pTBinT = ipTbin; // to fill from the 1st bin
		}

		// fill trigger pT spectra in mixed event
		if(!(zBin == -1 || cBin == -1)) fHistPtMixed[cBin][zBin]->Fill(firstpT);

		int ntrackSecond = particleMixed->GetEntriesFast();

		//	cout << "firstEta = " << firstEta << endl;
		//	cout << "firstPhi = " << firstPhi << endl;

		for(int is = 0; is < ntrackSecond; is++) 
		{ 
//			AliCorrReducedTrackAA* second = dynamic_cast<AliCorrReducedTrackAA*>(particleMixed->At(is));
			//		if(is==targetindex) cout << "second(init) = " << second << endl;
//			if (!second) { AliInfo("AcceptTracks: Could not receive second track"); continue; }
//			if(!fTrackCuts->AcceptTrack(second)) continue;

            if(!(particleMixed->At(is))) continue;

			float secondpT = -1; float secondEta = -99; float secondPhi = -99; float secondZv = -99; int secondPID = -1; // init
/*
			secondpT = second->Pt();
			secondEta = second->Eta();
			secondPhi = second->Phi();
//			secondZv = second->Zv();
			secondZv = second->Charge();
            secondPID = second->GetMyPartID();
*/
            secondpT = ((AliCorrReducedTrackAA*)particleMixed->At(is))->Pt();
            secondEta = ((AliCorrReducedTrackAA*)particleMixed->At(is))->Eta();
            secondPhi = ((AliCorrReducedTrackAA*)particleMixed->At(is))->Phi();
            secondZv = ((AliCorrReducedTrackAA*)particleMixed->At(is))->Charge();
            secondPID = ((AliCorrReducedTrackAA*)particleMixed->At(is))->GetMyPartID();



			int pTBinA = -1; // pT bin of associated particles
			for(int ipTbin = 0; ipTbin<kpTBin; ipTbin++) {
				if(0 > secondpT) continue;
				float pTBinMin = pTBinArray[ipTbin];
				float pTBinMax = pTBinArray[ipTbin+1];
				if((pTBinMin <= secondpT) && (secondpT < pTBinMax)) pTBinA = ipTbin; // to fill from the 1st bin
			}

			if( cBin == -1 || zBin == -1 || pTBinT == -1 || pTBinA == -1 ) continue; // exclude what have initial values


			if( firstpT < secondpT ) continue; // select pT_trigg >= pT_assoc

			if( TMath::Abs(secondEta)>0.8 ) continue;

/*
			// instead of excluding tracks, switch the information of pTt and pTa (dEta, dPhi, bins and diffTrigEP) 
			if( firstpT < secondpT ) {
				double bufferpTBin, bufferEta, bufferPhi;
				bufferpTBin = pTBinT;
				pTBinT = pTBinA;
				pTBinA = bufferpTBin;

				bufferEta = firstEta;
				firstEta = secondEta;
				secondEta = bufferEta;

				bufferPhi = firstPhi;
				firstPhi = secondPhi;
				secondPhi = bufferPhi;

//				diffTrigEP = firstPhi - EPangle; // reset phiS with 'second' track after the switch
			}
*/


			float deltaEtaMixed, deltaPhiMixed; // delta Eta - delat Phi in mixed events
			deltaEtaMixed = firstEta - secondEta;
			deltaPhiMixed = firstPhi - secondPhi;

			if(deltaPhiMixed < -0.5*TMath::Pi()) deltaPhiMixed += TMath::TwoPi();
			if(deltaPhiMixed > 1.5*TMath::Pi()) deltaPhiMixed -= TMath::TwoPi();


///*	INACTIVATED 140519
			float dPhiStarMixed = -99; // for two-track efficiency cut
/*
			cout << "firstEta = " << firstEta << " | firstPhi  = " << firstPhi << endl;
			cout << "secondEta = " << secondEta << " | secondPhi  = " << secondPhi << endl;
			cout << "deltaPhiMixed = " << deltaPhiMixed << " | deltaEtaMixed = " << deltaEtaMixed << endl;
			cout << "firstpT = " << firstpT << " | secondpT = " << secondpT << endl;
			cout << "firstZv = " << firstZv << " | secondZv = " << secondZv << " | bSign = " << bSign << endl;
*/
			dPhiStarMixed = CalculatedPhiStar(deltaPhiMixed, deltaEtaMixed, firstZv, secondZv, firstpT, secondpT, bSignHistos);

//			cout << "dPhiStarMixed = " << dPhiStarMixed << endl;

			// two-track efficiency cut
//			if(TMath::Abs(dPhiStarMixed)<0.02) continue;
//			if(TMath::Abs(deltaEtaMixed)<0.02) continue;
            if(TMath::Abs(dPhiStarMixed)<0.02 && TMath::Abs(deltaEtaMixed)<0.02) continue;

//*/

            fHistdEtadPhiMixed[cBin][zBin][pTBinT][pTBinA][0]->Fill(deltaEtaMixed, deltaPhiMixed, useWeight);

/*
            if(!(secondPID == -1)) {
                fHistdEtadPhiMixed[cBin][zBin][pTBinT][pTBinA][secondPID]->Fill(dEta, dPhi, useWeight);
            }
*/			

			// Classify by the angle between trigger and EP
			if(EPangle == -10) continue; // exclude default value
			double diffTrigEP = -10;
			diffTrigEP = firstPhi - EPangle;
			if(diffTrigEP < 0) diffTrigEP = diffTrigEP + 2*TMath::Pi();
			
			// H : near EP axis, M : intermidiate, L : perpendicular to EP
			if( (0 < diffTrigEP) && (diffTrigEP < (1./6.)*TMath::Pi()) ) fHistdEtadPhiMixedH[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEtaMixed, deltaPhiMixed, useWeight);
			if( ((5./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (7./6.)*TMath::Pi()) ) fHistdEtadPhiMixedH[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEtaMixed, deltaPhiMixed, useWeight);
			if( ((11./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (12./6.)*TMath::Pi()) ) fHistdEtadPhiMixedH[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEtaMixed, deltaPhiMixed, useWeight);

			if( ((1./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (2./6.)*TMath::Pi()) ) fHistdEtadPhiMixedMone[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEtaMixed, deltaPhiMixed, useWeight);
			if( ((7./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (8./6.)*TMath::Pi()) ) fHistdEtadPhiMixedMone[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEtaMixed, deltaPhiMixed, useWeight);
			if( ((4./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (5./6.)*TMath::Pi()) ) fHistdEtadPhiMixedMtwo[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEtaMixed, deltaPhiMixed, useWeight);
			if( ((10./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (11./6.)*TMath::Pi()) ) fHistdEtadPhiMixedMtwo[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEtaMixed, deltaPhiMixed, useWeight);

			if( ((2./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (4./6.)*TMath::Pi()) ) fHistdEtadPhiMixedL[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEtaMixed, deltaPhiMixed, useWeight);
			if( ((8./6.)*TMath::Pi() < diffTrigEP ) && (diffTrigEP < (10./6.)*TMath::Pi()) ) fHistdEtadPhiMixedL[cBin][zBin][pTBinT][pTBinA]->Fill(deltaEtaMixed, deltaPhiMixed, useWeight);





			//			cout << "|dEta| = " << TMath::Abs(dEta) << " cBin = " << cBin << " zBin = " << zBin << " pTBinT = " << pTBinT << " pTBinA = " << pTBinA << endl;
			//		if(is==targetindex) cout << "secondEta = " << secondEta << endl;
			//		if(is==targetindex) cout << "secondPhi = " << secondPhi << endl;
		} // end of is
	} // end of itargetindex

    delete particle;
    delete particleMixed;

}

//________________________________________________________________________
void AliAnalysisTaskEPCorrAA::DoMixing(double cent, double zvertex, double EPangle, TObjArray* trigTracks, float bSignDoMixing)
{
	AliEventPool* pool = fPoolMgr->GetEventPool(cent, zvertex);
	if (!pool) {
		AliInfo(Form("No pool found for centrality = %f, zVtx = %f", cent, zvertex));
	}
	else
	{
		//		cout << "GetCurrentNEvents() = " << pool->GetCurrentNEvents() << endl;
		if (pool->IsReady() || pool->NTracksInPool() > fMinNumTrack / 10 || pool->GetCurrentNEvents() >= fMinNEventsToMix)
//		if (pool->GetCurrentNEvents()==5) // previous set as 5 (141017)
//		if ( (pool->GetCurrentNEvents() >= 5) && ( pool->GetCurrentNEvents() < 10 ) ) // loop over 5~10 events (till 150413)
//		if ( (pool->GetCurrentNEvents() >= 5) && ( pool->GetCurrentNEvents() < 7 ) ) 
//		if ( (pool->GetCurrentNEvents() == 7) ) 
		{
			pool->PrintInfo();
			Int_t nMix = pool->GetCurrentNEvents();

			// Get cBin and zBin to pass those to FillMixedHistos
//			float cBinArray[] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}; // 10 bins
//			float cBinArray[] = {0, 5, 10, 20, 30, 40, 50, 60, 90}; // 8 bins
//			float cBinArray[] = {0, 10, 30, 60, 90}; // 4 bins (from 141223)
			float cBinArray[] = {0, 10, 20, 50, 60, 90}; // 5 bins for LHC11h
//			float cBinArray[] = {0, 10, 20, 35, 50, 60, 90}; // 6 bins for LHC11h 150519

			// Get cBin from percentage
			int cBin = -1;
			for(int icbin = 0; icbin<kCentBin; icbin++) {
				if(0 > cent) continue;
				float cBinMin = cBinArray[icbin];
				float cBinMax = cBinArray[icbin+1];
				if((cBinMin <= cent) && (cent < cBinMax)) cBin = icbin;	
			}

			// Fill zVertex and zBin
//			float zBinArray[] = {-8, -6, -4, -2, 0, 2, 4, 6, 8}; // 8 bins
			float zBinArray[] = {-7, -5, -3, -1, 1, 3, 5, 7}; // 8 bins

			int zBin = -1;
			for(int izbin = 0; izbin<kZvertBin; izbin++) {
				if(zvertex == -999.) continue; // exclude what have initial value
				float zBinMin = zBinArray[izbin];
				float zBinMax = zBinArray[izbin+1];
				if((zBinMin <= zvertex) && (zvertex < zBinMax)) zBin = izbin; // to fill from the 1st bin
			}

			//			cout << "cBin = " << cBin << " | zBin = " << zBin << endl;

			for (Int_t jMix=0; jMix<nMix; jMix++)
			{
				//				cout << "jMix = " << jMix << endl;
				fTracksMixing = pool->GetEvent(jMix);
				if (!fTracksMixing) continue;
/*
				cout << "Apply trackcut on the 2nd track in mixed" << endl;
				fTracksMixing = AcceptTracksReduced(fTracksMixingEvt, kTRUE);
				cout << "End : Apply trackcut on the 2nd track in mixed" << endl;
*/
//				fMyprimRecoTracks = AcceptTracks(fTracksMixing, kTrackCut);

				// count NevtMixed for nomalization
				if(!(zBin == -1 || cBin == -1)) fHistNevtMixed[cBin][zBin]->Fill(zBin*10 + cBin); // 8 zBins and 10 cBin be filled

				// Fill sub-information in jMix-th event to fTracksMixing to proceed to next step!
				FillMixedHistos( trigTracks, fTracksMixing, cBin, zBin, EPangle, bSignDoMixing, 1./nMix );
			} // end of jMix
			pool->UpdatePool(trigTracks);
			//			fPoolMgr->SetupForMixing(); // hmm.. this kind of init can be allowed?
			//			pool = 0x0; // initialize pool to prevent double counted events. Can be done by GetCurrentEvents == 5

		} // end of pool condition cut (end of if)
		else { pool->UpdatePool(trigTracks); }
	}//_____ pool NULL check (end of else)
} // end of DoMixing


/*
double AliAnalysisTaskEPCorrAA::DeltaPhi(double phi1, double phi2) {
	// dphi

	double res = phi1-phi2;

//	double res = TMath::ATan2(TMath::Sin(phi1-phi2), TMath::Cos(phi1-phi2));
//	return res>-(TMath::Pi())*9./20. ? res : 2*(TMath::Pi()) +res ;
	if (res < -0.5*TMath::Pi()) res += TMath::TwoPi();
	if (res > 1.5*TMath::Pi()) res -= TMath::TwoPi();

//	double res = 1.;
	return res;

}
*/




//________________________________________________________________________
TObjArray* AliAnalysisTaskEPCorrAA::AcceptTracksReduced(AliAODEvent *aodevent, bool useCuts) {
	Int_t nTracks = aodevent->GetNumberOfTracks();

	TObjArray* tracksAccepted = new TObjArray;
	tracksAccepted->SetOwner(kTRUE);

	for (Int_t itrk=0; itrk < nTracks; itrk++)
	{
		AliAODTrack* track = dynamic_cast<AliAODTrack*>(aodevent->GetTrack(itrk));
		if (!track) { AliInfo("AcceptTracks: Could not receive track"); continue; }

		// You can put some track cuts here
//		if(useCuts) { if(!fTrackCuts->AcceptTrack(track)) continue; }
        if(useCuts) {
            if(track->GetType() != AliAODTrack::kPrimary) continue;
//            if(!track->TestFilterBit(kHybridTrackCut)) continue;
        	if(!track->TestFilterBit(kTPCOnlyTrackCut)) continue;

        }

//        double nsigmaCut = 3.0;
/*
        int myPartID = -1;
        int myPartIDinit = -1;

        myPartIDinit = GetParticleID(track);
        if(myPartIDinit==1 || myPartIDinit==3 || myPartIDinit==5 || myPartIDinit==7) myPartID = 1; // recon. pion
        if(myPartIDinit==2 || myPartIDinit==3 || myPartIDinit==6 || myPartIDinit==7) myPartID = 2; // recon. kaon
        if(myPartIDinit==4 || myPartIDinit==5 || myPartIDinit==6 || myPartIDinit==7) myPartID = 3; // recon. proton
*/
		int myPartID = 0;

		tracksAccepted->Add(new AliCorrReducedTrackAA(myPartID,track->Eta(),track->Phi(),track->Pt(),track->Zv(),track->Charge()));
		//		tracksAccepted->Add(track);

	} // end of track loop

	return tracksAccepted;

//	delete tracksAccepted;

}

//___________________________________________________________
int AliAnalysisTaskEPCorrAA::GetParticleID(AliAODTrack* track)
{
//  cout << "GetParticleID CALLED!!!" << endl;


    int trackPID = 0;

/* // disabled at 150629 to exclude PID part

    float nsigmaTPC[kPID]; float nsigmaTOF[kPID];
    // float nSigmaITS[kPID];
    int TPCok, TOFok, ITSok = -1;


    if (!fPIDResponse) { AliInfo("GetParticleID: no fPIDResponse to proceed"); return -1; }

//  if(!(fPIDResponse==0)) {cout << "fPIDResponse exists!" << endl; cout << "fPIDResponse = " << fPIDResponse << endl; }

    TPCok = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,(AliVTrack*)track );
    TOFok = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,(AliVTrack*)track );
//  ITSok = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,track );

//  cout << "TPCok = " << TPCok << endl;

    if(!(TPCok==-1 || TOFok==-1)) {

//		cout << "fPIDResponse = " << fPIDResponse << endl;

        nsigmaTPC[0] = fPIDResponse->NumberOfSigmasTPC((AliVTrack*)track, AliPID::kPion);
        nsigmaTPC[1] = fPIDResponse->NumberOfSigmasTPC((AliVTrack*)track, AliPID::kKaon);
        nsigmaTPC[2] = fPIDResponse->NumberOfSigmasTPC((AliVTrack*)track, AliPID::kProton);

        nsigmaTOF[0] = fPIDResponse->NumberOfSigmasTOF((AliVTrack*)track, AliPID::kPion);
        nsigmaTOF[1] = fPIDResponse->NumberOfSigmasTOF((AliVTrack*)track, AliPID::kKaon);
        nsigmaTOF[2] = fPIDResponse->NumberOfSigmasTOF((AliVTrack*)track, AliPID::kProton);


//      cout << "nsigmaTPC[0] = " <<  nsigmaTPC[0] << " nsigmaTPC[1] = " << nsigmaTPC[1] << " nsigmaTPC[2] = " << nsigmaTPC[2] << endl;

//  nsigmaITS[0] = fPIDResponse->NumberOfSigmasITS(track, AliPID::kPion);
//  nsigmaITS[1] = fPIDResponse->NumberOfSigmasITS(track, AliPID::kKaon);
//  nsigmaITS[2] = fPIDResponse->NumberOfSigmasITS(track, AliPID::kProton);

//		cout << "nsigmaTPC[0] = " << nsigmaTPC[0] << endl;
//		cout << "nsigmaTOF[0] = " << nsigmaTOF[0] << endl;


        if(TMath::Sqrt( nsigmaTPC[0]*nsigmaTPC[0] + nsigmaTOF[0]*nsigmaTOF[0] ) < fNsigmaCut) trackPID = trackPID + 1;
        if(TMath::Sqrt( nsigmaTPC[1]*nsigmaTPC[1] + nsigmaTOF[1]*nsigmaTOF[1] ) < fNsigmaCut) trackPID = trackPID + 2;
        if(TMath::Sqrt( nsigmaTPC[2]*nsigmaTPC[2] + nsigmaTOF[2]*nsigmaTOF[2] ) < fNsigmaCut) trackPID = trackPID + 4;

    } // end of detectorOK

//	cout << "fNsigmaCut = " << fNsigmaCut << endl;
//	cout << "trackPID = " << trackPID << endl;
*/
    return trackPID;

/*
if (!trk) { AliInfo(" ==== zero track pointer."); return fPartUndefined; }

Bool_t* pidFlag; pidFlag = new Bool_t[AliPID::kSPECIES];
CalculateNSigmas((AliAODTrack*)trk, centbin, pidFlag, fillQA);
delete[] pidFlag;
Int_t mypid = -1;
mypid = FindNSigma((AliAODTrack*)trk);
//     Printf(" >>>>>>>>>>>>>>>>> mypid = %d",mypid);
return mypid;
*/
}

float AliAnalysisTaskEPCorrAA::CalculatedPhiStar(float dPhi, float dEta, float Zv, float Zv2, float pT, float pT2, float bSigntmp) {
	float dPhiStar = 999; // init
	float dPhiStartemp = 0; // init
	float Bze = 0.075;
	float radius = 0; // init
	float twoTrackEfficiencyCutValue = 0.02;

	const float kLimit = twoTrackEfficiencyCutValue * 3;
	static const Double_t kPi = TMath::Pi();

	if(TMath::Abs(dEta) < twoTrackEfficiencyCutValue * 2.5 * 3) {

		float dphistarInner = dPhi + bSigntmp*TMath::ASin((Zv*Bze*0.8)/(2*pT)) - bSigntmp*TMath::ASin((Zv2*Bze*0.8)/(2*pT2));
		float dphistarOuter = dPhi + bSigntmp*TMath::ASin((Zv*Bze*2.5)/(2*pT)) - bSigntmp*TMath::ASin((Zv2*Bze*2.5)/(2*pT2));

		if (dphistarInner > kPi) dphistarInner = kPi * 2 - dphistarInner;
		if (dphistarInner < -kPi) dphistarInner = -kPi * 2 - dphistarInner;
		if (dphistarInner > kPi) dphistarInner = kPi * 2 - dphistarInner; // might look funny but is needed

		if (dphistarOuter > kPi) dphistarOuter = kPi * 2 - dphistarOuter;
		if (dphistarOuter < -kPi) dphistarOuter = -kPi * 2 - dphistarOuter;
		if (dphistarOuter > kPi) dphistarOuter = kPi * 2 - dphistarOuter; // might look funny but is needed

		if(TMath::Abs(dphistarInner) < kLimit || TMath::Abs(dphistarOuter) < kLimit || dphistarInner * dphistarOuter < 0) {

			// find the smallest deltaPhistar
			for(int ir = 80; ir < 251; ir++) { // radius inside TPC vary for 0.80 ~ 2.50 (m)
				radius = ir*0.01;
				dPhiStartemp = dPhi + bSigntmp*TMath::ASin((Zv*Bze*radius)/(2*pT)) - bSigntmp*TMath::ASin((Zv2*Bze*radius)/(2*pT2));
				dPhiStartemp = TMath::Abs(dPhiStartemp);
				if(ir==80) {dPhiStar = dPhiStartemp;} // init once
				if(dPhiStartemp < dPhiStar) { dPhiStar = dPhiStartemp; }
				//		cout << "ir = " << ir << " | dPhiStar = " << dPhiStar << endl; 
			}

		}

		if (dPhiStar > kPi) dPhiStar = kPi * 2 - dPhiStar;
		if (dPhiStar < -kPi) dPhiStar = -kPi * 2 - dPhiStar;
		if (dPhiStar > kPi) dPhiStar = kPi * 2 - dPhiStar; // might look funny but is needed

	}

	return dPhiStar; // don't neet to concern in case dPhiStar has default value 999. It will pass 0.02 cut anyway

}





//________________________________________________________________________
void AliAnalysisTaskEPCorrAA::Terminate(Option_t *) 
{
	cout << "DEBUG : You are inside Terminate()" << endl;
	// Draw result to screen, or perform fitting, normalizations
	// Called once at the end of the query

    AliAnalysisTaskSE::Terminate();

	fOutput = dynamic_cast<TList*> (GetOutputData(1));
	if(!fOutput) { Printf("ERROR: could not retrieve TList fOutput"); return; }

	// Get the physics selection histograms with the selection statistics
	//AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	//AliESDInputHandler *inputH = dynamic_cast<AliESDInputHandler*>(mgr->GetInputEventHandler());
	//TH2F *histStat = (TH2F*)inputH->GetStatistics();

	// NEW HISTO should be retrieved from the TList container in the above way,
	// so it is available to draw on a canvas such as below

    if(fPoolMgr) delete fPoolMgr;

    return;
//	cout << "DEBUG : End of Treminate" << endl;
}
