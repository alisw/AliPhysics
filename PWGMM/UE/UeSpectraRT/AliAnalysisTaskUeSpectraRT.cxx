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

//Authors: Sergio Iga: sigabuit@cern.ch for the MB dNdpT part
//	   Valentina Zaccolo: vzaccolo@cern.ch for the UE part
//         Aditya Nath Mishra, Aditya.Nath.Mishra@cern.ch for high RT and pT leading part

#include "AliAnalysisTaskUeSpectraRT.h"

// ROOT includes
#include <TList.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TProfile.h>
#include <TParticle.h>
#include <TFile.h>
#include <TVector3.h>

// AliRoot includes
#include <AliAnalysisManager.h>
#include <AliAnalysisFilter.h>
#include <AliInputEventHandler.h>
#include <AliAODInputHandler.h>
#include <AliAODEvent.h>
#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliESDVertex.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrackCuts.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>
#include "AliPhysicsSelection.h"

#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>

#include <TTreeStream.h>

#include <AliHeader.h>

#include "AliCentrality.h"
#include <AliESDv0.h>
#include <AliKFVertex.h>
#include <AliAODVertex.h>

#include <AliAODTrack.h>
#include <AliAODPid.h>
#include <AliAODMCHeader.h>
#include <AliAODHeader.h>


//#include "AliPPVsMultUtils.h"
#include "AliAnalysisUtils.h"

#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliMultSelectionTask.h"


//#include "AliMCSpectraWeights.h"

// STL includes
#include <iostream>
using namespace std;

ClassImp(AliAnalysisTaskUeSpectraRT)
	//_____________________________________________________________________________
	AliAnalysisTaskUeSpectraRT::AliAnalysisTaskUeSpectraRT():
	AliAnalysisTaskSE(),
	fMCEvent(0x0),
	fESD(0x0),
	fAOD(0x0),
	fAnalysisMC(kFALSE),
	fAnalysisCorr(kTRUE),
	fTrackFilterDCA(0x0),
	fTrackFilterMatchEff(0x0),
	fAnalysisType("ESD"),
	ftrigBit(0x0),
	fRandom(0x0),
	fPileUpRej(kFALSE),
// 	fPPVsMultUtils(0x0),
	fVtxCut(10.0),
	fEtaCut(0.9),
	fLeadMin(6.0),
	fTriggeredEventMB(-999),
	fVtxStatus(-999),
	fZvtx(-999),
	fZvtx_SPD(-9999),
	fRun(-999),
	fEventId(-999),
	fdcaxy(-999),
	fdcaz(-999),
	fisPS(kFALSE),
	fisTracklet(kFALSE),
	fisMCvtxInZcut(kFALSE),
	fListOfObjects(0),
	fEvents(0x0),
	fVtxBeforeCuts(0x0),
	fVtxAfterCuts(0x0),
	fn1(0x0),
	fNchargedTrue(0x0),
        isINEL0Rec(kFALSE),
	isINEL0True(kFALSE),
	INEL0(0x0),
	INEL0Gen(0x0),
	TrigINEL0(0x0),
	TrueINEL0(0x0),
	ptvsdcaData(0x0),
	ptvsdcaPrim(0x0),
	ptvsdcaDecs(0x0),
	ptvsdcaMatl(0x0),
	ptvsdcacentralData(0x0),
	ptvsdcacentralPrim(0x0),
	ptvsdcacentralDecs(0x0),
	ptvsdcacentralMatl(0x0),
	effcomputationGen(0x0),
        effcomputationGen1(0x0),
        effcomputationGen2(0x0),
        effcomputationGen3(0x0),
        effcomputationGen4(0x0),
        effcomputationGen5(0x0),
	sigLossTrueINEL0(0x0),
	sigLossTrigINEL0(0x0),
        sigLossTrueINEL01(0x0),
        sigLossTrigINEL01(0x0),
        sigLossTrueINEL02(0x0),
        sigLossTrigINEL02(0x0),
        sigLossTrueINEL03(0x0),
        sigLossTrigINEL03(0x0),
        sigLossTrueINEL04(0x0),
        sigLossTrigINEL04(0x0),
        sigLossTrueINEL05(0x0),
        sigLossTrigINEL05(0x0),
        sigLossTrueINEL06(0x0),
        sigLossTrigINEL06(0x0),
        sigLossTrueINEL07(0x0),
        sigLossTrigINEL07(0x0),
        sigLossTrueINEL08(0x0),
        sigLossTrigINEL08(0x0),
        sigLossTrueINEL09(0x0),
        sigLossTrigINEL09(0x0),
        sigLossTrueINEL010(0x0),
        sigLossTrigINEL010(0x0),
        sigLossTrueINEL011(0x0),
        sigLossTrigINEL011(0x0),
	effcomputationGenPi(0x0),
	effcomputationGenK(0x0),
	effcomputationGenP(0x0),
	effcomputationGenSm(0x0),
	effcomputationGenSp(0x0),
	effcomputationGenO(0x0),
	effcomputationGenXi(0x0),
	effcomputationGenL(0x0),
	effcomputationGenRest(0x0),
        effcomputationRec(0x0),
        effcomputationRec1(0x0),
        effcomputationRec2(0x0),
        effcomputationRec3(0x0),
        effcomputationRec4(0x0),
        effcomputationRec5(0x0),
        effcomputationRec6(0x0),
        effcomputationRec7(0x0),
        effcomputationRec8(0x0),
        effcomputationRec9(0x0),
        effcomputationRec10(0x0),
        effcomputationRec11(0x0),
	effcomputationRecPi(0x0),
	effcomputationRecK(0x0),
	effcomputationRecP(0x0),
	effcomputationRecSm(0x0),
	effcomputationRecSp(0x0),
	effcomputationRecO(0x0),
	effcomputationRecXi(0x0),
	effcomputationRecL(0x0),
	effcomputationRecRest(0x0),
	fPS_MC(0x0),
	nchtrue(0x0),
	fVtxPS_MC(0x0),
	fPS(0x0),
	fVtxPS(0x0),
	Zvtx(0x0),
	Phi(0x0),
	Eta(0x0),
	EtaPhi(0x0),
	fhRTResponse(0x0),
	fAveMultiInTrans(4.939),
	fAveRecMultiInTrans(4.895),
	fAveGenMultiInTrans(7.392),
	//fMCSpectraWeights(0),
	fhRTData(0x0),
        fhRTTrue(0x0),
        fhRTReco(0x0),
	ptiME(0x0),
	pTGen(0x0),
	pTGenTrans(0x0),
	pTGenTow(0x0),
        pTGenAw(0x0),
        ptliGen1(0x0),
        ptliGen2(0x0),
        ptliGen3(0x0),
        ptliGen4(0x0),
        ptliGen5(0x0),
        ptliGen6(0x0),
        ptliGen7(0x0),
        ptliGen8(0x0),
        ptliGen9(0x0),
        ptliGen10(0x0),
        ptliGen11(0x0),
	pTGenTrans1(0x0),
        pTGenTrans2(0x0),
        pTGenTrans3(0x0),
        pTGenTrans4(0x0),
        pTGenTrans5(0x0),
	pTGenTrans6(0x0),
	pTGenTrans7(0x0),
	pTGenTrans8(0x0),
	pTGenTrans9(0x0),
	pTGenTrans10(0x0),
	pTGenTrans11(0x0),
        pTGenTow1(0x0),
        pTGenTow2(0x0),
        pTGenTow3(0x0),
        pTGenTow4(0x0),
        pTGenTow5(0x0),
	pTGenTow6(0x0),
	pTGenTow7(0x0),
	pTGenTow8(0x0),
	pTGenTow9(0x0),
	pTGenTow10(0x0),
	pTGenTow11(0x0),
        pTGenAw1(0x0),
        pTGenAw2(0x0),
        pTGenAw3(0x0),
        pTGenAw4(0x0),
        pTGenAw5(0x0),
	pTGenAw6(0x0),
	pTGenAw7(0x0),
	pTGenAw8(0x0),
	pTGenAw9(0x0),
	pTGenAw10(0x0),
        pTGenAw11(0x0),
        ptli1(0x0),
        ptli2(0x0),
        ptli3(0x0),
        ptli4(0x0),
        ptli5(0x0),
        ptli6(0x0),
        ptli7(0x0),
        ptli8(0x0),
        ptli9(0x0),
        ptli10(0x0),
        ptli11(0x0),
        ptliMB(0x0),
        hRefMult08(0x0),
        hV0Mmult(0x0)

{
	for ( int i = 0 ; i < 18 ; i++ )
	{
	  secondaries[i] = 0;
	  primariesTrackFilter[i] = 0;
	  pti[i] = 0;
          pti1Trans[i] = 0;
          pti2Trans[i] = 0;
          pti3Trans[i] = 0;
          pti4Trans[i] = 0;
          pti5Trans[i] = 0;
	  pti6Trans[i] = 0;
	  pti7Trans[i] = 0;
	  pti8Trans[i] = 0;
	  pti9Trans[i] = 0;
	  pti10Trans[i] = 0;
	  pti11Trans[i] = 0;
	  ptiMBTrans[i] =0;
	  pti1Tow[i] = 0;
          pti2Tow[i] = 0;
          pti3Tow[i] = 0;
          pti4Tow[i] = 0;
          pti5Tow[i] = 0;
	  pti6Tow[i] = 0;
	  pti7Tow[i] = 0;
	  pti8Tow[i] = 0;
	  pti9Tow[i] = 0;
	  pti10Tow[i] = 0;
	  pti11Tow[i] = 0;
	  ptiMBTow[i] =0;
          pti1Aw[i] = 0;
          pti2Aw[i] = 0;
          pti3Aw[i] = 0;
          pti4Aw[i] = 0;
          pti5Aw[i] = 0;
	  pti6Aw[i] = 0;
	  pti7Aw[i] = 0;
	  pti8Aw[i] = 0;
	  pti9Aw[i] = 0;
	  pti10Aw[i] = 0;
	  pti11Aw[i] = 0;
	  ptiMBAw[i] =0;
	  ptiNoLead[i] = 0;
	  fTrackFilter[i] = 0;

	}
}

//_____________________________________________________________________________
AliAnalysisTaskUeSpectraRT::AliAnalysisTaskUeSpectraRT(const char *name):
	AliAnalysisTaskSE(name),
	fMCEvent(0x0),
	fESD(0x0),
	fAOD(0x0),
	fAnalysisMC(kFALSE),
	fAnalysisCorr(kTRUE),
	fTrackFilterDCA(0x0),
	fTrackFilterMatchEff(0x0),
	fAnalysisType("ESD"),
	ftrigBit(0x0),
	fRandom(0x0),
	fPileUpRej(kFALSE),
// 	fPPVsMultUtils(0x0),
	fVtxCut(10.0),
	fEtaCut(0.9),
	fLeadMin(6.0),
	fTriggeredEventMB(-999),
	fVtxStatus(-999),
	fZvtx(-999),
	fZvtx_SPD(-9999),
	fRun(-999),
	fEventId(-999),
	fdcaxy(-999),
	fdcaz(-999),
        fisPS(kFALSE),
	fisTracklet(kFALSE),
	fisMCvtxInZcut(kFALSE),
	fListOfObjects(0),
	fEvents(0x0),
	fVtxBeforeCuts(0x0),
	fVtxAfterCuts(0x0),
	fn1(0x0),
	fNchargedTrue(0x0),
        isINEL0Rec(kFALSE),
	isINEL0True(kFALSE),
	INEL0(0x0),
	INEL0Gen(0x0),
	TrigINEL0(0x0),
	TrueINEL0(0x0),
	ptvsdcaData(0x0),
	ptvsdcaPrim(0x0),
	ptvsdcaDecs(0x0),
	ptvsdcaMatl(0x0),
	ptvsdcacentralData(0x0),
	ptvsdcacentralPrim(0x0),
	ptvsdcacentralDecs(0x0),
	ptvsdcacentralMatl(0x0),
	effcomputationGen(0x0),
        effcomputationGen1(0x0),
        effcomputationGen2(0x0),
        effcomputationGen3(0x0),
        effcomputationGen4(0x0),
        effcomputationGen5(0x0),
	sigLossTrueINEL0(0x0),
	sigLossTrigINEL0(0x0),
        sigLossTrueINEL01(0x0),
        sigLossTrigINEL01(0x0),
        sigLossTrueINEL02(0x0),
        sigLossTrigINEL02(0x0),
        sigLossTrueINEL03(0x0),
        sigLossTrigINEL03(0x0),
        sigLossTrueINEL04(0x0),
        sigLossTrigINEL04(0x0),
        sigLossTrueINEL05(0x0),
        sigLossTrigINEL05(0x0),
	sigLossTrueINEL06(0x0),
        sigLossTrigINEL06(0x0),
        sigLossTrueINEL07(0x0),
        sigLossTrigINEL07(0x0),
        sigLossTrueINEL08(0x0),
        sigLossTrigINEL08(0x0),
        sigLossTrueINEL09(0x0),
        sigLossTrigINEL09(0x0),
        sigLossTrueINEL010(0x0),
        sigLossTrigINEL010(0x0),
        sigLossTrueINEL011(0x0),
        sigLossTrigINEL011(0x0),
	effcomputationGenPi(0x0),
	effcomputationGenK(0x0),
	effcomputationGenP(0x0),
	effcomputationGenSm(0x0),
	effcomputationGenSp(0x0),
	effcomputationGenO(0x0),
	effcomputationGenXi(0x0),
	effcomputationGenL(0x0),
	effcomputationGenRest(0x0),
        effcomputationRec(0x0),
        effcomputationRec1(0x0),
        effcomputationRec2(0x0),
        effcomputationRec3(0x0),
        effcomputationRec4(0x0),
        effcomputationRec5(0x0),
	effcomputationRec6(0x0),
        effcomputationRec7(0x0),
        effcomputationRec8(0x0),
        effcomputationRec9(0x0),
        effcomputationRec10(0x0),
        effcomputationRec11(0x0),
	effcomputationRecPi(0x0),
	effcomputationRecK(0x0),
	effcomputationRecP(0x0),
	effcomputationRecSm(0x0),
	effcomputationRecSp(0x0),
	effcomputationRecO(0x0),
	effcomputationRecXi(0x0),
	effcomputationRecL(0x0),
	effcomputationRecRest(0x0),
	fPS_MC(0x0),
	nchtrue(0x0),
	fVtxPS_MC(0x0),
	fPS(0x0),
	fVtxPS(0x0),
	Zvtx(0x0),
        Phi(0x0),
        Eta(0x0),
        EtaPhi(0x0),
	fhRTResponse(0x0),
	fAveMultiInTrans(4.939),
	fAveRecMultiInTrans(4.895),
	fAveGenMultiInTrans(7.614),
	//fMCSpectraWeights(0),
	fhRTData(0x0),
        fhRTTrue(0x0),
        fhRTReco(0x0),
	ptiME(0x0),
        pTGen(0x0),
        pTGenTrans(0x0),
        pTGenTow(0x0),
        pTGenAw(0x0),
	ptliGen1(0x0),
        ptliGen2(0x0),
        ptliGen3(0x0),
        ptliGen4(0x0),
        ptliGen5(0x0),
        ptliGen6(0x0),
        ptliGen7(0x0),
        ptliGen8(0x0),
        ptliGen9(0x0),
        ptliGen10(0x0),
        ptliGen11(0x0),
        pTGenTrans1(0x0),
        pTGenTrans2(0x0),
        pTGenTrans3(0x0),
        pTGenTrans4(0x0),
        pTGenTrans5(0x0),
	pTGenTrans6(0x0),
	pTGenTrans7(0x0),
	pTGenTrans8(0x0),
	pTGenTrans9(0x0),
	pTGenTrans10(0x0),
	pTGenTrans11(0x0),
        pTGenTow1(0x0),
        pTGenTow2(0x0),
        pTGenTow3(0x0),
        pTGenTow4(0x0),
        pTGenTow5(0x0),
	pTGenTow6(0x0),
	pTGenTow7(0x0),
	pTGenTow8(0x0),
	pTGenTow9(0x0),
	pTGenTow10(0x0),
	pTGenTow11(0x0),
        pTGenAw1(0x0),
        pTGenAw2(0x0),
        pTGenAw3(0x0),
        pTGenAw4(0x0),
        pTGenAw5(0x0),
	pTGenAw6(0x0),
	pTGenAw7(0x0),
	pTGenAw8(0x0),
	pTGenAw9(0x0),
	pTGenAw10(0x0),
	pTGenAw11(0x0),
	ptli1(0x0),
        ptli2(0x0),
        ptli3(0x0),
        ptli4(0x0),
        ptli5(0x0),
        ptli6(0x0),
        ptli7(0x0),
        ptli8(0x0),
        ptli9(0x0),
        ptli10(0x0),
        ptli11(0x0),
        ptliMB(0x0),
	hRefMult08(0x0),
        hV0Mmult(0x0)

{
	for ( int i = 0 ; i < 18 ; i++ )
	{
	  secondaries[i] = 0;
	  primariesTrackFilter[i] = 0;
	  pti[i] = 0;
          pti1Trans[i] = 0;
          pti2Trans[i] = 0;
          pti3Trans[i] = 0;
          pti4Trans[i] = 0;
          pti5Trans[i] = 0;
	  pti6Trans[i] = 0;
	  pti7Trans[i] = 0;
	  pti8Trans[i] = 0;
	  pti9Trans[i] = 0;
	  pti10Trans[i] = 0;
	  pti11Trans[i] = 0;
	  ptiMBTrans[i] =0;
	  pti1Tow[i] = 0;
          pti2Tow[i] = 0;
          pti3Tow[i] = 0;
          pti4Tow[i] = 0;
          pti5Tow[i] = 0;
	  pti6Tow[i] = 0;
	  pti7Tow[i] = 0;
	  pti8Tow[i] = 0;
	  pti9Tow[i] = 0;
	  pti10Tow[i] = 0;
	  pti11Tow[i] = 0;
	  ptiMBTow[i] =0;
          pti1Aw[i] = 0;
          pti2Aw[i] = 0;
          pti3Aw[i] = 0;
          pti4Aw[i] = 0;
          pti5Aw[i] = 0;
	  pti6Aw[i] = 0;
	  pti7Aw[i] = 0;
	  pti8Aw[i] = 0;
	  pti9Aw[i] = 0;
	  pti10Aw[i] = 0;
	  pti11Aw[i] = 0;
	  ptiMBAw[i] =0;
	  ptiNoLead[i] = 0;
	  fTrackFilter[i] = 0;
	}
	DefineOutput(1, TList::Class());

}

//_____________________________________________________________________________
AliAnalysisTaskUeSpectraRT::~AliAnalysisTaskUeSpectraRT()
{
  //
  // Destructor
  //

  if (fListOfObjects) {
    delete fListOfObjects;
    fListOfObjects = 0x0;
  }

/*  if (fPPVsMultUtils)
  {
    delete fPPVsMultUtils;
    fPPVsMultUtils = 0x0;
  }
 */
}

//______________________________________________________________________________
void AliAnalysisTaskUeSpectraRT::UserCreateOutputObjects()
{
  fRandom = new TRandom(0); // 0 means random seed

  //OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();

  //
  // Histograms
  //
  fEvents = new TH1I("fEvents","Number of analyzed events; Events; Counts", 6, 0, 6);
  fListOfObjects->Add(fEvents);

  fVtxBeforeCuts = new TH1F("fVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 12, -30, 30);
  fListOfObjects->Add(fVtxBeforeCuts);

  fVtxAfterCuts = new TH1F("fVtxAfterCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 12, -30, 30);
  fListOfObjects->Add(fVtxAfterCuts);


  fn1=new TH1F("fn1","fn1",10,0,10);
  fListOfObjects->Add(fn1);

  const Int_t nPtBins = 79;

  Double_t xBins[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
			       0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
			       1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
			       2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
			       4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
			       11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
			       26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 55.0,
			       60.0, 65.0, 70.0, 75.0, 80.0, 90.0, 100.0,130.0,160.0,200.0};



  const int nBinsDCAxy = 121;
  double binsDCAxy[] = {-3.025,-2.975,-2.925,-2.875,-2.825,-2.775,-2.725,-2.675,-2.625,-2.575,-2.525,-2.475,-2.425,-2.375,-2.325,-2.275,-2.225, -2.175,-2.125,-2.075,-2.025,-1.975,-1.925,-1.875,-1.825,-1.775,-1.725,-1.675,-1.625,-1.575,-1.525,-1.475,-1.425,-1.375,-1.325,-1.275,-1.225,-1.175,-1.125,-1.075,-1.025,-0.975,-0.925,-0.875,-0.825,-0.775,-0.725,-0.675,-0.625,-0.575,-0.525,-0.475,-0.425,-0.375,-0.325,-0.275,-0.225,-0.175,-0.125,-0.075,-0.025,0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,0.825,0.875,0.925,0.975,1.025,1.075,1.125,1.175,1.225,1.275,1.325,1.375,1.425,1.475,1.525,1.575,1.625,1.675,1.725,1.775,1.825,1.875,1.925,1.975,2.025,2.075,2.125,2.175,2.225,2.275,2.325,2.375,2.425,2.475,2.525,2.575,2.625,2.675,2.725,2.775,2.825,2.875,2.925,2.975,3.025};

  const int nBinsDCAxyCentral = 461;

  double binsDCAxyCentral[]= {-0.2305,-0.2295,-0.2285,-0.2275,-0.2265,-0.2255,-0.2245,-0.2235,-0.2225,-0.2215,-0.2205,-0.2195,-0.2185,-0.2175,-0.2165,-0.2155,-0.2145,-0.2135,-0.2125,-0.2115,-0.2105,-0.2095,-0.2085,-0.2075,-0.2065,-0.2055,-0.2045,-0.2035,-0.2025,-0.2015,-0.2005,-0.1995,-0.1985,-0.1975,-0.1965,-0.1955,-0.1945,-0.1935,-0.1925,-0.1915,-0.1905,-0.1895,-0.1885,-0.1875,-0.1865,-0.1855,-0.1845,-0.1835,-0.1825,-0.1815,-0.1805,-0.1795,-0.1785,-0.1775,-0.1765,-0.1755,-0.1745,-0.1735,-0.1725,-0.1715,-0.1705,-0.1695,-0.1685,-0.1675,-0.1665,-0.1655,-0.1645,-0.1635,-0.1625,-0.1615,-0.1605,-0.1595,-0.1585,-0.1575,-0.1565,-0.1555,-0.1545,-0.1535,-0.1525,-0.1515,-0.1505,-0.1495,-0.1485,-0.1475,-0.1465,-0.1455,-0.1445,-0.1435,-0.1425,-0.1415,-0.1405,-0.1395,-0.1385,-0.1375,-0.1365,-0.1355,-0.1345,-0.1335,-0.1325,-0.1315,-0.1305,-0.1295,-0.1285,-0.1275,-0.1265,-0.1255,-0.124,-0.1235,-0.1225,-0.1215,-0.1205,-0.1195,-0.1185,-0.1175,-0.1165,-0.1155,-0.1145,-0.1135,-0.1125,-0.1115,-0.1105,-0.1095,-0.1085,-0.1075,-0.1065,-0.1055,-0.1045,-0.1035,-0.1025,-0.1015,-0.1005,-0.0995,-0.0985,-0.0975,-0.0965,-0.0955,-0.0945,-0.0935,-0.0925,-0.0915,-0.0905,-0.0895,-0.0885,-0.0875,-0.0865,-0.0855,-0.0845,-0.0835,-0.0825,-0.0815,-0.0805,-0.0795,-0.0785,-0.0775,-0.0765,-0.0755,-0.0745,-0.0735,-0.0725,-0.0715,-0.0705,-0.0695,-0.0685,-0.0675,-0.0665,-0.0655,-0.0645,-0.0635,-0.0625,-0.0615,-0.0605,-0.0595,-0.0585,-0.0575,-0.0565,-0.0555,-0.0545,-0.0535,-0.0525,-0.0515,-0.0505,-0.0495,-0.0485,-0.0475,-0.0465,-0.0455,-0.0445,-0.0435,-0.0425,-0.0415,-0.0405,-0.0395,-0.0385,-0.0375,-0.0365,-0.0355,-0.0345,-0.0335,-0.0325,-0.0315,-0.0305,-0.0295,-0.0285,-0.0275,-0.0265,-0.0255,-0.0245,-0.0235,-0.0225,-0.0215,-0.0205,-0.0195,-0.0185,-0.0175,-0.0165,-0.0155,-0.0145,-0.0135,-0.0125,-0.0115,-0.0105,-0.0095,-0.0085,-0.0075,-0.0065,-0.0055,-0.0045,-0.0035,-0.0025,-0.0015,-0.0005,0.0005,0.0015,0.0025,0.0035,0.0045,0.0055,0.0065,0.0075,0.0085,0.0095,0.0105,0.0115,0.0125,0.0135,0.0145,0.0155,0.0165,0.0175,0.0185,0.0195,0.0205,0.0215,0.0225,0.0235,0.0245,0.0255,0.0265,0.0275,0.0285,0.0295,0.0305,0.0315,0.0325,0.0335,0.0345,0.0355,0.0365,0.0375,0.0385,0.0395,0.0405,0.0415,0.0425,0.0435,0.0445,0.0455,0.0465,0.0475,0.0485,0.0495,0.0505,0.0515,0.0525,0.0535,0.0545,0.0555,0.0565,0.0575,0.0585,0.0595,0.0605,0.0615,0.0625,0.0635,0.0645,0.0655,0.0665,0.0675,0.0685,0.0695,0.0705,0.0715,0.0725,0.0735,0.0745,0.0755,0.0765,0.0775,0.0785,0.0795,0.0805,0.0815,0.0825,0.0835,0.0845,0.0855,0.0865,0.0875,0.0885,0.0895,0.0905,0.0915,0.0925,0.0935,0.0945,0.0955,0.0965,0.0975,0.0985,0.0995,0.1005,0.1015,0.1025,0.1035,0.1045,0.1055,0.1065,0.1075,0.1085,0.1095,0.1105,0.1115,0.1125,0.1135,0.1145,0.1155,0.1165,0.1175,0.1185,0.1195,0.1205,0.1215,0.1225,0.1235,0.1245,0.1255,0.1265,0.1275,0.1285,0.1295,0.1305,0.1315,0.1325,0.1335,0.1345,0.1355,0.1365,0.1375,0.1385,0.1395,0.1405,0.1415,0.1425,0.1435,0.1445,0.1455,0.1465,0.1475,0.1485,0.1495,0.1505,0.1515,0.1525,0.1535,0.1545,0.1555,0.1565,0.1575,0.1585,0.1595,0.1605,0.1615,0.1625,0.1635,0.1645,0.1655,0.1665,0.1675,0.1685,0.1695,0.1705,0.1715,0.1725,0.1735,0.1745,0.1755,0.1765,0.1775,0.1785,0.1795,0.1805,0.1815,0.1825,0.1835,0.1845,0.1855,0.1865,0.1875,0.1885,0.1895,0.1905,0.1915,0.1925,0.1935,0.1945,0.1955,0.1965,0.1975,0.1985,0.1995,0.2005,0.2015,0.2025,0.2035,0.2045,0.2055,0.2065,0.2075,0.2085,0.2095,0.2105,0.2115,0.2125,0.2135,0.2145,0.2155,0.2165,0.2175,0.2185,0.2195,0.2205,0.2215,0.2225,0.2235,0.2245,0.2255,0.2265,0.2275,0.2285,0.2295,0.2305};

  const Int_t nRTBins = 90;
  Double_t RTEdge_Data[nRTBins + 1] = {0};
  for(Int_t i = 0; i <= nRTBins; i++){
   RTEdge_Data[i] = i/fAveMultiInTrans-1/(2*fAveMultiInTrans);
   }

  Double_t RTEdge_Rec[nRTBins + 1] = {0};
  for(Int_t i = 0; i <= nRTBins; i++){
   RTEdge_Rec[i] = i/fAveRecMultiInTrans-1/(2*fAveRecMultiInTrans);
   }

  Double_t RTEdge_Gen[nRTBins + 1] = {0};
  for(Int_t i = 0; i <= nRTBins; i++){
   RTEdge_Gen[i] = i/fAveGenMultiInTrans-1/(2*fAveGenMultiInTrans);
   }

  const Int_t multBins = 200;
  const Double_t xmultBins[201] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
				  20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,
				  40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,
				  60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,
				  80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,
				  100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,
				  120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,
				  140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,
				  160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,
				  180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,
				  200};
  Int_t nphibins = 200;
  Float_t phimin = 0.;
  Float_t phimax = TMath::TwoPi();
  Int_t netabins = 300;
  Float_t etamin = -1.5;
  Float_t etamax = 1.5;

  AliESDtrackCuts* esdTrackCutsRun2[18] = {0};

  for ( int iTc = 0 ; iTc < 18 ; iTc++ )
  {
	// standar parameters ------------------- //
    double maxdcaz = 2.;
    double minratiocrossrowstpcover = 0.8;
    double maxfraclusterstpcshared = 0.4;
    double maxchi2perclustertpc = 4.0;
    double maxchi2perclusterits = 36.;
    double geowidth = 3.;
    double geolenght = 130.;
    double maxchi2tpcglobal = 36.;
	// ------------------------------------- //

	// variations of the track cuts -------- //
    if ( iTc == 1) maxdcaz = 1.0;
    if ( iTc == 2) maxdcaz = 5.0;
    if ( iTc == 5) minratiocrossrowstpcover = 0.7;
    if ( iTc == 6) minratiocrossrowstpcover = 0.9;
    if ( iTc == 7) maxfraclusterstpcshared = 0.2;
    if ( iTc == 8) maxfraclusterstpcshared = 1.0;
    if ( iTc == 9) maxchi2perclustertpc = 3.0;
    if ( iTc == 10) maxchi2perclustertpc = 5.0;
    if ( iTc == 11) maxchi2perclusterits = 25.0;
    if ( iTc == 12) maxchi2perclusterits = 49.0;
    if ( iTc == 14) geowidth = 2.0;
    if ( iTc == 15) geowidth = 4.0;
    if ( iTc == 16) geolenght = 120.0;
    if ( iTc == 17) geolenght = 140.0;

	// variations of the track cuts -------- //


    fTrackFilter[iTc] = new AliAnalysisFilter(Form("fTrackFilter%d",iTc));
    esdTrackCutsRun2[iTc] = new AliESDtrackCuts(Form("esdTrackCutsRun2%d",iTc));

    // TPC

    esdTrackCutsRun2[iTc]->SetCutGeoNcrNcl(geowidth,geolenght,1.5,0.85,0.7);
    esdTrackCutsRun2[iTc]->SetRequireTPCRefit(kTRUE);
    esdTrackCutsRun2[iTc]->SetMinRatioCrossedRowsOverFindableClustersTPC(minratiocrossrowstpcover);
    esdTrackCutsRun2[iTc]->SetMaxChi2PerClusterTPC(maxchi2perclustertpc);
    esdTrackCutsRun2[iTc]->SetMaxFractionSharedTPCClusters(maxfraclusterstpcshared);
    //esdTrackCutsRun2[iTc]->SetMaxChi2TPCConstrainedGlobal(maxchi2tpcglobal); TODO VZ: check this cut

    // ITS
    esdTrackCutsRun2[iTc]->SetRequireITSRefit(kTRUE);
    if ( iTc != 13 )
    esdTrackCutsRun2[iTc]->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);

    esdTrackCutsRun2[iTc]->SetMaxChi2PerClusterITS(maxchi2perclusterits);

    // primary selection
    esdTrackCutsRun2[iTc]->SetDCAToVertex2D(kFALSE);
    esdTrackCutsRun2[iTc]->SetRequireSigmaToVertex(kFALSE);
    esdTrackCutsRun2[iTc]->SetMaxDCAToVertexZ(maxdcaz);
    esdTrackCutsRun2[iTc]->SetAcceptKinkDaughters(kFALSE);

    if ( iTc == 3 )
      // esdTrackCutsRun2[iTc]->SetMaxDCAToVertexXYPtDep("4*(0.0026+0.0050/pt^1.01)");
    	esdTrackCutsRun2[iTc]->SetMaxDCAToVertexXYPtDep("6.5*(0.0026+0.0050/pt^1.01)");
    else if ( iTc == 4 )
      // esdTrackCutsRun2[iTc]->SetMaxDCAToVertexXYPtDep("10*(0.0026+0.0050/pt^1.01)");
    	esdTrackCutsRun2[iTc]->SetMaxDCAToVertexXYPtDep("7.5*(0.0026+0.0050/pt^1.01)");
    else
      esdTrackCutsRun2[iTc]->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); // (7*(------))

    fTrackFilter[iTc]->AddCuts(esdTrackCutsRun2[iTc]);
  }

  // ------------------------------------------------------------------ //

  // track cuts for Feed-Down correction ------------------------------ //

  AliESDtrackCuts* esdTrackCutsDCA = 0;

  double maxdcaz = 2.;
  double minratiocrossrowstpcover = 0.8;
  double maxfraclusterstpcshared = 0.4;
  double maxchi2perclustertpc = 4.0;
  double maxchi2perclusterits = 36.;
  double geowidth = 3.;
  double geolenght = 130.;

  fTrackFilterDCA = new AliAnalysisFilter("fTrackFilterDCA");
  esdTrackCutsDCA = new AliESDtrackCuts("esdTrackCutsDCA");

  // TPC

  esdTrackCutsDCA->SetRequireTPCRefit(kTRUE);
  esdTrackCutsDCA->SetMinRatioCrossedRowsOverFindableClustersTPC(minratiocrossrowstpcover);
  esdTrackCutsDCA->SetMaxChi2PerClusterTPC(maxchi2perclustertpc);
  esdTrackCutsDCA->SetMaxFractionSharedTPCClusters(maxfraclusterstpcshared);

  // ITS
  esdTrackCutsDCA->SetRequireITSRefit(kTRUE);
  esdTrackCutsDCA->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCutsDCA->SetMaxChi2PerClusterITS(maxchi2perclusterits);

  // primary selection
  esdTrackCutsDCA->SetDCAToVertex2D(kFALSE);
  esdTrackCutsDCA->SetRequireSigmaToVertex(kFALSE);
  esdTrackCutsDCA->SetMaxDCAToVertexZ(maxdcaz);
  esdTrackCutsDCA->SetAcceptKinkDaughters(kFALSE);
  esdTrackCutsDCA->SetCutGeoNcrNcl(geowidth,geolenght,1.5,0.85,0.7); // Affects more or less 10% !!
  fTrackFilterDCA->AddCuts(esdTrackCutsDCA);






  //___________________________________________________//
  // track cuts for matching efficiency _______________//

  AliESDtrackCuts* esdTrackCutsMatchEff = 0;

  double maxdcazME = 3.;
  double maxdcaxy = 3.;

  fTrackFilterMatchEff = new AliAnalysisFilter("fTrackFilterMatchEff");
  esdTrackCutsMatchEff = new AliESDtrackCuts("esdTrackCutsMatchEff");

  // TPC
  esdTrackCutsMatchEff->SetRequireTPCRefit(kTRUE);
  esdTrackCutsMatchEff->SetMinRatioCrossedRowsOverFindableClustersTPC(minratiocrossrowstpcover);
  esdTrackCutsMatchEff->SetMaxChi2PerClusterTPC(maxchi2perclustertpc);
  esdTrackCutsMatchEff->SetMaxFractionSharedTPCClusters(maxfraclusterstpcshared);

  // primary selection
  esdTrackCutsMatchEff->SetAcceptKinkDaughters(kFALSE);
  esdTrackCutsMatchEff->SetCutGeoNcrNcl(geowidth,geolenght,1.5,0.85,0.7);
  esdTrackCutsMatchEff->SetMaxDCAToVertexZ(maxdcazME);
  esdTrackCutsMatchEff->SetMaxDCAToVertexXY(maxdcaxy);
  fTrackFilterMatchEff->AddCuts(esdTrackCutsMatchEff);



  // ----------------------------------------------------------------------------------------------------- //

  if ( fAnalysisMC )
  {
	// generated pT
    effcomputationGen = new TH1F("effcomputationGen","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationGen);

    effcomputationGen1 = new TH1F("effcomputationGen1","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationGen1);

    effcomputationGen2 = new TH1F("effcomputationGen2","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationGen2);

    effcomputationGen3 = new TH1F("effcomputationGen3","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationGen3);

    effcomputationGen4 = new TH1F("effcomputationGen4","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationGen4);

    effcomputationGen5 = new TH1F("effcomputationGen5","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationGen5);

	// MC truth histos
    pTGen = new TH1F("pTGen","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGen);

    pTGenTrans = new TH1F("pTGenTrans","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTrans);

    pTGenTow = new TH1F("pTGenTow","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTow);

    pTGenAw = new TH1F("pTGenAw","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenAw);

    //ptleading
    ptliGen1 = new TH1F("ptl1","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptliGen1);
      ptliGen2 = new TH1F("ptl2","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptliGen2);
      ptliGen3 = new TH1F("ptl3","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptliGen3);
      ptliGen4 = new TH1F("ptl4","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptliGen4);
      ptliGen5 = new TH1F("ptl5","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptliGen5);
      ptliGen6 = new TH1F("ptl6","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptliGen6);
      ptliGen7 = new TH1F("ptl7","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptliGen7);
      ptliGen8 = new TH1F("ptl8","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptliGen8);
      ptliGen9 = new TH1F("ptl9","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptliGen9);
      ptliGen10 = new TH1F("ptl10","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptliGen10);
      ptliGen11 = new TH1F("ptl11","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptliGen11);

    //Transverse
    pTGenTrans1 = new TH1F("pTGenTrans1","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTrans1);

    pTGenTrans2 = new TH1F("pTGenTrans2","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTrans2);

    pTGenTrans3 = new TH1F("pTGenTrans3","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTrans3);

    pTGenTrans4 = new TH1F("pTGenTrans4","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTrans4);

    pTGenTrans5 = new TH1F("pTGenTrans5","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTrans5);

    pTGenTrans6 = new TH1F("pTGenTrans6","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTrans6);

    pTGenTrans7 = new TH1F("pTGenTrans7","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTrans7);

    pTGenTrans8 = new TH1F("pTGenTrans8","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTrans8);

    pTGenTrans9 = new TH1F("pTGenTrans9","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTrans9);

    pTGenTrans10 = new TH1F("pTGenTrans0","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTrans10);

    pTGenTrans11 = new TH1F("pTGenTrans11","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTrans11);

    //Toward
    pTGenTow1 = new TH1F("pTGenTow1","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTow1);

    pTGenTow2 = new TH1F("pTGenTow2","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTow2);

    pTGenTow3 = new TH1F("pTGenTow3","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTow3);

    pTGenTow4 = new TH1F("pTGenTow4","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTow4);

    pTGenTow5 = new TH1F("pTGenTow5","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTow5);

    pTGenTow6 = new TH1F("pTGenTow6","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTow6);

    pTGenTow7 = new TH1F("pTGenTow7","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTow7);

    pTGenTow8 = new TH1F("pTGenTow8","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTow8);

    pTGenTow9 = new TH1F("pTGenTow9","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTow9);

    pTGenTow10 = new TH1F("pTGenTow10","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTow10);

    pTGenTow11 = new TH1F("pTGenTow11","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenTow11);

    //Away
    pTGenAw1 = new TH1F("pTGenAw1","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenAw1);

    pTGenAw2 = new TH1F("pTGenAw2","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenAw2);

    pTGenAw3 = new TH1F("pTGenAw3","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenAw3);

    pTGenAw4 = new TH1F("pTGenAw4","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenAw4);

    pTGenAw5 = new TH1F("pTGenAw5","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenAw5);

    pTGenAw6 = new TH1F("pTGenAw6","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenAw6);

    pTGenAw7 = new TH1F("pTGenAw7","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenAw7);

    pTGenAw8 = new TH1F("pTGenAw8","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenAw8);

    pTGenAw9 = new TH1F("pTGenAw9","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenAw9);

    pTGenAw10 = new TH1F("pTGenAw10","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenAw10);

    pTGenAw11 = new TH1F("pTGenAw11","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(pTGenAw11);



	// histos for signal loss

    sigLossTrueINEL0 = new TH1F("sigLossTrueINEL0","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrueINEL0);

    sigLossTrigINEL0 = new TH1F("sigLossTrigINEL0","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrigINEL0);

    sigLossTrueINEL01 = new TH1F("sigLossTrueINEL01","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrueINEL01);

    sigLossTrigINEL01 = new TH1F("sigLossTrigINEL01","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrigINEL01);

    sigLossTrueINEL02 = new TH1F("sigLossTrueINEL02","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrueINEL02);

    sigLossTrigINEL02 = new TH1F("sigLossTrigINEL02","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrigINEL02);

    sigLossTrueINEL03 = new TH1F("sigLossTrueINEL03","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrueINEL03);

    sigLossTrigINEL03 = new TH1F("sigLossTrigINEL03","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrigINEL03);

    sigLossTrueINEL04 = new TH1F("sigLossTrueINEL04","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrueINEL04);

    sigLossTrigINEL04 = new TH1F("sigLossTrigINEL04","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrigINEL04);

    sigLossTrueINEL05 = new TH1F("sigLossTrueINEL05","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrueINEL05);

    sigLossTrigINEL05 = new TH1F("sigLossTrigINEL05","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrigINEL05);

    sigLossTrueINEL06 = new TH1F("sigLossTrueINEL06","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrueINEL06);

    sigLossTrigINEL06 = new TH1F("sigLossTrigINEL06","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrigINEL06);

    sigLossTrueINEL07 = new TH1F("sigLossTrueINEL07","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrueINEL07);

    sigLossTrigINEL07 = new TH1F("sigLossTrigINEL07","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrigINEL07);

    sigLossTrueINEL08 = new TH1F("sigLossTrueINEL08","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrueINEL08);

    sigLossTrigINEL08 = new TH1F("sigLossTrigINEL08","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrigINEL08);

    sigLossTrueINEL09 = new TH1F("sigLossTrueINEL09","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrueINEL09);

    sigLossTrigINEL09 = new TH1F("sigLossTrigINEL09","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrigINEL09);

    sigLossTrueINEL010 = new TH1F("sigLossTrueINEL010","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrueINEL010);

    sigLossTrigINEL010 = new TH1F("sigLossTrigINEL010","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrigINEL010);

    sigLossTrueINEL011 = new TH1F("sigLossTrueINEL011","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrueINEL011);

    sigLossTrigINEL011 = new TH1F("sigLossTrigINEL011","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(sigLossTrigINEL011);

	// histos for event loss (trigger eff)

    TrigINEL0 = new TH1F("TrigINEL0","",15,0,15);
    fListOfObjects->Add(TrigINEL0);

    TrueINEL0 = new TH1F("TrueINEL0","",15,0,15);
    fListOfObjects->Add(TrueINEL0);


	// event counter
    INEL0Gen = new TH1F("INEL0Gen","", 15, 0, 15);
    fListOfObjects->Add(INEL0Gen);

	// histos for Vtx correction MC ( not used, only to see the difference with data)

    nchtrue = new TH1F("nchtrue","N_{ch}^{True}N_{ch}^{True}",multBins,xmultBins);
    fListOfObjects->Add(nchtrue);

    fPS_MC = new TH1F("fPS_MC","",1,0,1);
    fListOfObjects->Add(fPS_MC);

    fVtxPS_MC = new TH1F("fVtxPS_MC","",1,0,1);
    fListOfObjects->Add(fVtxPS_MC);

    // HISTOS FOR FEED_DOWN CORRECTION --------------------------------------

    ptvsdcaPrim = new TH2F("ptvsdcaPrim","pt vs dca Primaries",nPtBins,xBins,nBinsDCAxy,binsDCAxy);
    fListOfObjects->Add(ptvsdcaPrim);

    ptvsdcaDecs = new TH2F("ptvsdcaDecs","pt vs dca Decays",nPtBins,xBins,nBinsDCAxy,binsDCAxy);
    fListOfObjects->Add(ptvsdcaDecs);

    ptvsdcaMatl = new TH2F("ptvsdcaMatl","pt vs dca Material",nPtBins,xBins,nBinsDCAxy,binsDCAxy);
    fListOfObjects->Add(ptvsdcaMatl);


    ptvsdcacentralPrim = new TH2F("ptvsdcacentralPrim","pt vs dca central Primaries",nPtBins,xBins,nBinsDCAxyCentral,binsDCAxyCentral);
    fListOfObjects->Add(ptvsdcacentralPrim);

    ptvsdcacentralDecs = new TH2F("ptvsdcacentralDecs","pt vs dca central Decays",nPtBins,xBins,nBinsDCAxyCentral,binsDCAxyCentral);
    fListOfObjects->Add(ptvsdcacentralDecs);

    ptvsdcacentralMatl = new TH2F("ptvsdcacentralMatl","pt vs dca central Material",nPtBins,xBins,nBinsDCAxyCentral,binsDCAxyCentral);
    fListOfObjects->Add(ptvsdcacentralMatl);

    // -------------------------------------------------------------------


	// histos for particle composition

    //pions
    effcomputationGenPi = new TH1F("effcomputationGenPi","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationGenPi);
    //kaons
    effcomputationGenK = new TH1F("effcomputationGenK","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationGenK);
    //protons
    effcomputationGenP = new TH1F("effcomputationGenP","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationGenP);
    //Sigma-
    effcomputationGenSm = new TH1F("effcomputationGenSm","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationGenSm);
    //sigma+
    effcomputationGenSp = new TH1F("effcomputationGenSp","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationGenSp);
    //omega-
    effcomputationGenO = new TH1F("effcomputationGenO","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationGenO);
    //xi-
    effcomputationGenXi = new TH1F("effcomputationGenXi","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationGenXi);
    //Lambda
    effcomputationGenL = new TH1F("effcomputationGenL","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationGenL);
    //Others
    effcomputationGenRest = new TH1F("effcomputationGenRest","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationGenRest);




    effcomputationRec = new TH1F("effcomputationRec","p_{T}^{rec};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRec);

    effcomputationRec1 = new TH1F("effcomputationRec1","p_{T}^{rec};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRec1);

    effcomputationRec2 = new TH1F("effcomputationRec2","p_{T}^{rec};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRec2);

    effcomputationRec3 = new TH1F("effcomputationRec3","p_{T}^{rec};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRec3);

    effcomputationRec4 = new TH1F("effcomputationRec4","p_{T}^{rec};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRec4);

    effcomputationRec5 = new TH1F("effcomputationRec5","p_{T}^{rec};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRec5);

    effcomputationRec6 = new TH1F("effcomputationRec6","p_{T}^{rec};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRec6);

    effcomputationRec7 = new TH1F("effcomputationRec7","p_{T}^{rec};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRec7);

    effcomputationRec8 = new TH1F("effcomputationRec8","p_{T}^{rec};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRec8);

    effcomputationRec9 = new TH1F("effcomputationRec9","p_{T}^{rec};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRec9);

    effcomputationRec10 = new TH1F("effcomputationRec10","p_{T}^{rec};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRec10);

    effcomputationRec11 = new TH1F("effcomputationRec11","p_{T}^{rec};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRec11);

    //pions
    effcomputationRecPi = new TH1F("effcomputationRecPi","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRecPi);
    //kaons
    effcomputationRecK = new TH1F("effcomputationRecK","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRecK);
    //protons
    effcomputationRecP = new TH1F("effcomputationRecP","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRecP);
    //Sigma-
    effcomputationRecSm = new TH1F("effcomputationRecSm","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRecSm);
    //sigma+
    effcomputationRecSp = new TH1F("effcomputationRecSp","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRecSp);
    //omega-
    effcomputationRecO = new TH1F("effcomputationRecO","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRecO);
    //xi-
    effcomputationRecXi = new TH1F("effcomputationRecXi","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRecXi);
    //Lambda
    effcomputationRecL = new TH1F("effcomputationRecL","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRecL);
    //Others
    effcomputationRecRest = new TH1F("effcomputationRecRest","p_{T}^{gen};#it{p}_{T} (GeV/c)",nPtBins,xBins);
    fListOfObjects->Add(effcomputationRecRest);


    fhRTTrue = new TH1F("fhRTTrue","Event counts vs. R_{T}",nRTBins,RTEdge_Gen);
    fListOfObjects->Add(fhRTTrue);


    fhRTReco = new TH1F("fhRTReco","Event counts vs. R_{T}",nRTBins,RTEdge_Rec);
    fListOfObjects->Add(fhRTReco);

    fhRTResponse = new TH2F("hRTResponse","Event R_{T} response matrix", nRTBins, RTEdge_Gen, nRTBins, RTEdge_Rec);
    fListOfObjects->Add(fhRTResponse);

    for ( int i = 0 ; i < 18 ; i++ ) // rec pT, secondaries, and "primaries" (all) that passes track fileters, i == 0  standard
    {

      secondaries[i] = new TH1F(Form("secondaries%d",i),"p_{T}^{rec};#it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(secondaries[i]);

      primariesTrackFilter[i] = new TH1F(Form("primariesTrackFilter%d",i),"p_{T}^{rec};#it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(primariesTrackFilter[i]);
    }

      primariesTrackFilterME= new TH1F("primariesTrackFilterME","p_{T}^{rec};#it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(primariesTrackFilterME);

  }
  else // for the data
  {
    // histos for FEED_DOWN CORRECTION ------------------------------------

    ptvsdcaData = new TH2F("ptvsdcaData","pt vs dca Data",nPtBins,xBins,nBinsDCAxy,binsDCAxy);
    fListOfObjects->Add(ptvsdcaData);


    ptvsdcacentralData = new TH2F("ptvsdcacentralData","pt vs dca central Data",nPtBins,xBins,nBinsDCAxyCentral,binsDCAxyCentral);
    fListOfObjects->Add(ptvsdcacentralData);


    // ----------------------------------------------------------------


    // event counter

    INEL0 = new TH1F("INEL0","", 15, 0, 15);
    fListOfObjects->Add(INEL0);
	// histos for vtx corrections, data driven

    fPS = new TH1F("fPS","", 15, 0, 15);
    fListOfObjects->Add(fPS);

    fVtxPS = new TH1F("fVtxPS","", 15, 0, 15);
    fListOfObjects->Add(fVtxPS);

    fhRTData = new TH1F("fhRTData","Event counts vs. R_{T}",nRTBins,RTEdge_Rec);
    fListOfObjects->Add(fhRTData);

      ptli1 = new TH1F("ptl1","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptli1);
      ptli2 = new TH1F("ptl2","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptli2);
      ptli3 = new TH1F("ptl3","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptli3);
      ptli4 = new TH1F("ptl4","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptli4);
      ptli5 = new TH1F("ptl5","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptli5);
      ptli6 = new TH1F("ptl6","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptli6);
      ptli7 = new TH1F("ptl7","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptli7);
      ptli8 = new TH1F("ptl8","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptli8);
      ptli9 = new TH1F("ptl9","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptli9);
      ptli10 = new TH1F("ptl10","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptli10);
      ptli11 = new TH1F("ptl11","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptli11);

      ptliMB = new TH1F("ptlMB","p_{T}^{leading}; #it{p}_{T}^{leading} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptliMB);

      hRefMult08 = new TH1F("hRefMult08","Multiplicity (-0.8 < #eta < 0.8);RefMult08;count",multBins,xmultBins);
      fListOfObjects->Add(hRefMult08);

      hV0Mmult = new TH1F("hV0Mmult","V0M ;V0M percentile;count",multBins,xmultBins);
      fListOfObjects->Add(hV0Mmult);



    for ( int i = 0 ; i < 18 ; i++ )
    {
      // Histos for data -----------------------------------------------------------------------------------------------------------------

      pti[i] = new TH1F(Form("ptAll%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti[i]);
      pti1Trans[i] = new TH1F(Form("pt1Trans%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti1Trans[i]);
      pti2Trans[i] = new TH1F(Form("pt2Trans%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti2Trans[i]);
      pti3Trans[i] = new TH1F(Form("pt3Trans%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti3Trans[i]);
      pti4Trans[i] = new TH1F(Form("pt4Trans%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti4Trans[i]);
      pti5Trans[i] = new TH1F(Form("pt5Trans%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti5Trans[i]);
      pti6Trans[i] = new TH1F(Form("pt6Trans%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti6Trans[i]);
      pti7Trans[i] = new TH1F(Form("pt7Trans%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti7Trans[i]);
      pti8Trans[i] = new TH1F(Form("pt8Trans%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti8Trans[i]);
      pti9Trans[i] = new TH1F(Form("pt9Trans%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti9Trans[i]);
      pti10Trans[i] = new TH1F(Form("pt10Trans%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti10Trans[i]);
      pti11Trans[i] = new TH1F(Form("pt11Trans%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti11Trans[i]);

      ptiMBTrans[i] = new TH1F(Form("ptMBTrans%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptiMBTrans[i]);

      pti1Tow[i] = new TH1F(Form("pt1Tow%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti1Tow[i]);
      pti2Tow[i] = new TH1F(Form("pt2Tow%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti2Tow[i]);
      pti3Tow[i] = new TH1F(Form("pt3Tow%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti3Tow[i]);
      pti4Tow[i] = new TH1F(Form("pt4Tow%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti4Tow[i]);
      pti5Tow[i] = new TH1F(Form("pt5Tow%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti5Tow[i]);
      pti6Tow[i] = new TH1F(Form("pt6Tow%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti6Tow[i]);
      pti7Tow[i] = new TH1F(Form("pt7Tow%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti7Tow[i]);
      pti8Tow[i] = new TH1F(Form("pt8Tow%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti8Tow[i]);
      pti9Tow[i] = new TH1F(Form("pt9Tow%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti9Tow[i]);
      pti10Tow[i] = new TH1F(Form("pt10Tow%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti10Tow[i]);
      pti11Tow[i] = new TH1F(Form("pt11Tow%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti11Tow[i]);

      ptiMBTow[i] = new TH1F(Form("ptMBTow%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptiMBTow[i]);

      pti1Aw[i] = new TH1F(Form("pt1Aw%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti1Aw[i]);
      pti2Aw[i] = new TH1F(Form("pt2Aw%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti2Aw[i]);
      pti3Aw[i] = new TH1F(Form("pt3Aw%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti3Aw[i]);
      pti4Aw[i] = new TH1F(Form("pt4Aw%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti4Aw[i]);
      pti5Aw[i] = new TH1F(Form("pt5Aw%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti5Aw[i]);
      pti6Aw[i] = new TH1F(Form("pt6Aw%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti6Aw[i]);
      pti7Aw[i] = new TH1F(Form("pt7Aw%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti7Aw[i]);
      pti8Aw[i] = new TH1F(Form("pt8Aw%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti8Aw[i]);
      pti9Aw[i] = new TH1F(Form("pt9Aw%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti9Aw[i]);
      pti10Aw[i] = new TH1F(Form("pt10Aw%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti10Aw[i]);
      pti11Aw[i] = new TH1F(Form("pt11Aw%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(pti11Aw[i]);

      ptiMBAw[i] = new TH1F(Form("ptMBAw%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptiMBAw[i]);


      ptiNoLead[i] = new TH1F(Form("ptNoLead%d",i),"p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptiNoLead[i]);

    // ----------------------------------------------------------------------------------------------------------------------------------
    }

      ptiME = new TH1F("ptME","p_{T}; #it{p}_{T} (GeV/c)",nPtBins,xBins);
      fListOfObjects->Add(ptiME);

 }

    Zvtx = new TH1F("Zvtx"," ; Vtx_z ; Entries",40,-10,10);
    fListOfObjects->Add(Zvtx);


    Phi = new TH2F ("Phi","#phi",nPtBins, xBins, nphibins,phimin,phimax);
    fListOfObjects->Add(Phi);

    Eta  = new TH2F ("hEta","#eta", nPtBins, xBins, netabins,etamin,etamax);
    fListOfObjects->Add(Eta);

    EtaPhi  = new TH2F ("EtaPhi","eta vs. phi",netabins,etamin,etamax, nphibins,phimin,phimax);
    fListOfObjects->Add(EtaPhi);

  // Post output data.
  PostData(1,fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskUeSpectraRT::UserExec(Option_t *)
{
	// Main loop
	//
	// First we make sure that we have valid input(s)!
	//


	AliVEvent *event = InputEvent();
	if (!event)
	{
		Error("UserExec", "Could not retrieve event");
		return;
	}

	if (fAnalysisMC)
	{
		fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
		if (!fMCEvent)
		{
			cout<<"Could not retrieve MC event"<<endl;
			return;
		}

		AliHeader *header = fMCEvent->Header();
		if(!header) {AliDebug( AliLog::kError , "Header not avaible" ); return; }

	}


	fESD = dynamic_cast<AliESDEvent*>(event);
	if(!fESD)
	{
		Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
		this->Dump();
		return;
	}

	/************ BEGINNING OF EVENT SELECTION *******************/

	// Get trigger decision
	 fTriggeredEventMB = 0; //init
	 if (!fAnalysisMC){  // for data check event selection as well
		if((((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()) & ftrigBit )
	 	{
	 		fTriggeredEventMB = 1;  //event triggered as minimum bias
	 	}
	}
	else {
		if (ftrigBit) {
			fTriggeredEventMB = 1;
		}
	}
	Bool_t SPDvsClustersBG = kFALSE;

	AliAnalysisUtils *AnalysisUtils = new AliAnalysisUtils();
	if (!AnalysisUtils)
	{
		cout<<"------- No AnalysisUtils Object Found --------"<<AnalysisUtils<<endl;
		return;
	}
	else
		SPDvsClustersBG = AnalysisUtils->IsSPDClusterVsTrackletBG(fESD); // We want NO BG


	Bool_t isNotPileUp = !fESD->IsPileupFromSPD(5,0.8);
	Bool_t IncompleteDAQ = fESD->IsIncompleteDAQ(); // we want is not incomplete DAQ


	// vertex
	const AliESDVertex * vertex    =    fESD->GetPrimaryVertex(); // tracks vertex, if not -> spd vertex, if not TPC vertex

  	Bool_t isVtxGood = vertex->GetStatus() && selectVertex2015pp( fESD ,kTRUE,kFALSE,kTRUE); // requires Tracks and spd vertex, and Zconsistency of 5mm

 	double vertex_z = vertex->GetZ();
  	Bool_t isVtxInZCut = (TMath::Abs(vertex_z)   <= fVtxCut); // Zvtx in +- 10



   	// Implement INEL>0
   	const AliMultiplicity* mult = fESD->GetMultiplicity();
   	Bool_t isINEL0 = kFALSE;
   	for (Int_t i = 0; i < mult->GetNumberOfTracklets(); ++i) {
   		if (TMath::Abs(mult->GetEta(i)) < 1.) isINEL0 = kTRUE;
   	}

  	/********** IS PHYSICS SELECTION FLAG ****************************/

  	fisPS = fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG && isNotPileUp;


	// recontructed INEL > 0 is PS + vtx + Zvtx inside +-10 ------
	isINEL0Rec = kFALSE;
        if ( isINEL0 && fisPS && isVtxGood && isVtxInZCut) isINEL0Rec = kTRUE;

	if (fisPS) // VZ hack
		fVtxBeforeCuts->Fill(vertex_z);
	if (isINEL0Rec)
		fVtxAfterCuts->Fill(vertex_z);

	fn1->Fill(0); // all events

	if(fTriggeredEventMB)
		fn1->Fill(1); // triggered events

	if ( fTriggeredEventMB && !IncompleteDAQ )
		fn1->Fill(2); // trigger + IsIncompleteDAQ

	if (fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG)
		fn1->Fill(3); // trigger + IsIncompleteDAQ + BG rejection

	if(fPileUpRej)
	{
		if(fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG && isNotPileUp) //
			fn1->Fill(4); // trigger + IsIncompleteDAQ + BG rejection + PileUp
	}
	if (fisPS)
		fn1->Fill(5); //PS: trigger + IsIncompleteDAQ + BG rejection + PileUp + 1 Tracklet in eta +-1

	if (fisPS && isVtxGood)
		fn1->Fill(6); //PS + GetPrimaryVertex

	if (isINEL0Rec)
		fn1->Fill(7); //PS + GetPrimaryVertex + isVtxInZCut

	// -------------------------------------- multiplcity estimators section ------------------------------------------ //

	ftrackmult08 = -999;
	fv0mpercentile = -999;

	ftrackmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 0.8);     //tracklets
	hRefMult08->Fill(ftrackmult08);

	fMultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection"); // Esto es para 13 TeV
	if (!fMultSelection)
	  cout<<"------- No AliMultSelection Object Found --------"<<fMultSelection<<endl;
	else
	  fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");
	hV0Mmult->Fill(fv0mpercentile);
	cout<<"------- V0M mult ==  "<<fv0mpercentile<<"--------"<<endl;
	// ------------------------------------------ end of mult estimators -------------------------------------------------//


	//________ calculate RT and binned spectra ___


        if(isINEL0Rec && !fAnalysisMC) {
		 AnalyseDataRT(fESD);
	}
	if (fisPS && !fAnalysisMC && fAnalysisCorr){
		 CorrectionsDataRT(fESD, isVtxGood);
	}



	//______________________________________________

	if (fAnalysisMC) // analysis for generated MC
	{
		const AliVVertex *vertexMC = (AliVVertex*) fMCEvent->GetPrimaryVertex();
		fisMCvtxInZcut     = (TMath::Abs(vertexMC->GetZ()) <= fVtxCut);   // ZMCvtx in +- 10
		isINEL0True = isMCEventTrueINEL0(fMCEvent);

		// for trigger efficiency PS / INEL > 0 true
		if ( fisPS ) TrigINEL0->Fill(0);
		if (isINEL0True) TrueINEL0->Fill(0);
		// MC (not used) next two histos for missing vtx correction
		if ( fisPS ) fPS_MC->Fill(0);
	        if (fisPS && isVtxGood) fVtxPS_MC->Fill(0);
	        if (fAnalysisCorr) CorrectionsMCRT(fMCEvent, fESD); // analysis for MC RT binned
		AnalyzeMC(fMCEvent); // analysis for MC MB
        }
	else
	{
	  	if (isINEL0Rec) INEL0->Fill(0);
	  	// Two histos for missing vtx correction
	  	if ( fisPS ) fPS->Fill(0);
	  	if ( fisPS && isVtxGood ) fVtxPS->Fill(0);
 	}


	if(isINEL0Rec) // analysis for data MB
	{
		Zvtx->Fill(vertex_z);
		AnalyzeESD(fESD);
		AnalyzeESDforDCA(fESD);
	}
	//________ post output data __________________________________________
	PostData(1, fListOfObjects);
}

//________________________________________________________________________
void AliAnalysisTaskUeSpectraRT::CorrectionsMCRT(AliMCEvent* fMCEvent, AliESDEvent* esdEvent)
{

	// generated MC stuff
        TObjArray* fCTSTracksMC = new TObjArray();
        Double_t TrackRT  = 0.;
        Double_t nRecTracks = 0;
	Int_t counting =0;
	AliVParticle* LeadingMC = 0;
	TObjArray *regionsMinMaxMC = 0x0;

	for ( int iT = 0 ; iT < fMCEvent->GetNumberOfTracks(); iT++ ) // loop over TRUE MC
        {
		AliVParticle *mcParticle = fMCEvent->GetTrack(iT);
		if (!mcParticle)
          	{
             		cout<<"no mcParticle"<<endl;
             	continue;
          	}
          	if(!(fMCEvent->IsPhysicalPrimary(iT))) {
          		continue;
          	}

          	Double_t eta = mcParticle->Eta();
          	Double_t pt = mcParticle->Pt();

        	Int_t partPDG = TMath::Abs(mcParticle->PdgCode());
	 	if (TMath::Abs(eta) > fEtaCut) {
	        	continue;
	        }
          	if ( pt < 0.15 ){
            		continue;
          	}

		if ( !( TMath::Abs( mcParticle->Charge() ) == 3 ) ) continue;
		fCTSTracksMC->Add(mcParticle);
		if (!mcParticle) continue;
	}



      	// leading object
       	TObjArray* LeadingTrackMC = 0x0;
       	LeadingTrackMC = FindLeadingObjects(fCTSTracksMC);
 	if (LeadingTrackMC){
		LeadingMC = (AliVParticle*)LeadingTrackMC->At(0);
 		// Sorting
		TObjArray *regionSortedParticlesMC = 0x0;
                regionSortedParticlesMC = SortRegions((AliVParticle*)LeadingTrackMC->At(0), fCTSTracksMC);
		// Taking transverse regions
                regionsMinMaxMC = GetMinMaxRegion((TList*)regionSortedParticlesMC->At(2),(TList*)regionSortedParticlesMC->At(3));

		// Taking toward region
                TObjArray *towardMC = GetRegionAwTow((TList*)regionSortedParticlesMC->At(0));

                // Taking away region
                TObjArray *awayMC = GetRegionAwTow((TList*)regionSortedParticlesMC->At(1));

		TList *listMax = (TList*)regionsMinMaxMC->At(0);
                TList *listMin = (TList*)regionsMinMaxMC->At(1);
		TList *listTow = (TList*) towardMC->At(0);
		TList *listAw = (TList*) awayMC->At(0);

                TList *toward = new TList();
                TList *away = new TList();
		TList *transverse = new TList();

                // fill transverse list
                for(Int_t j = 0; j < listMax->GetEntries(); j++){
        		AliVParticle* particle1 = (AliVParticle*)listMax->At(j);
                        transverse->Add(particle1);
                }

                for(Int_t j = 0; j < listMin->GetEntries(); j++){
                        AliVParticle* particle2 = (AliVParticle*)listMin->At(j);
                        transverse->Add(particle2);
                }


                // fill toward list
                for(Int_t j = 0; j < listTow->GetEntries(); j++){
                        AliVParticle* particleTow = (AliVParticle*)listTow->At(j);
                        toward->Add(particleTow);
                }

                // fill away list
                for(Int_t j = 0; j < listAw->GetEntries(); j++){
                        AliVParticle* particleAw = (AliVParticle*)listAw->At(j);
                        away->Add(particleAw);
                }

        	Float_t LeadingPt = LeadingMC->Pt();
        	// select only plateau region
                if(LeadingPt > fLeadMin && LeadingPt < 300.){
			nRecTracks = transverse->GetEntries();
                	TrackRT = nRecTracks/fAveGenMultiInTrans;

			fhRTTrue->Fill(TrackRT);
			if (TrackRT<1.) INEL0Gen->Fill(1);
                        if (TrackRT>=1. && TrackRT<2.) INEL0Gen->Fill(2);
                        if (TrackRT>=2. && TrackRT<3.) INEL0Gen->Fill(3);
                        if (TrackRT>=3. && TrackRT<4.) INEL0Gen->Fill(4);
			if (TrackRT>=4. && TrackRT<5.) INEL0Gen->Fill(5);
			if (TrackRT>=5. && TrackRT<6.) INEL0Gen->Fill(6);
			if (TrackRT>=6. && TrackRT<7.) INEL0Gen->Fill(7);
			if (TrackRT>=7. && TrackRT<8.) INEL0Gen->Fill(8);
			if (TrackRT>=8. && TrackRT<9.) INEL0Gen->Fill(9);
			if (TrackRT>=9. && TrackRT<10.) INEL0Gen->Fill(10);
                        if (TrackRT>=10.) INEL0Gen->Fill(11);

			if (TrackRT<1.) ptliGen1->Fill(LeadingPt);
			if (TrackRT>=1. && TrackRT<2.) ptliGen2->Fill(LeadingPt);
			if (TrackRT>=2. && TrackRT<3.) ptliGen3->Fill(LeadingPt);
			if (TrackRT>=3. && TrackRT<4.) ptliGen4->Fill(LeadingPt);
			if (TrackRT>=4. && TrackRT<5.) ptliGen5->Fill(LeadingPt);
			if (TrackRT>=5. && TrackRT<6.) ptliGen6->Fill(LeadingPt);
			if (TrackRT>=6. && TrackRT<7.) ptliGen7->Fill(LeadingPt);
			if (TrackRT>=7. && TrackRT<8.) ptliGen8->Fill(LeadingPt);
			if (TrackRT>=8. && TrackRT<9.) ptliGen9->Fill(LeadingPt);
			if (TrackRT>=9. && TrackRT<10.) ptliGen10->Fill(LeadingPt);
			if (TrackRT>=10.) ptliGen11->Fill(LeadingPt);

			for(Int_t k = 0; k < nRecTracks; k++)
                        {
				AliVParticle* particle = (AliVParticle*)transverse->At(k);
				Float_t TrackPt  = particle->Pt();
			 	if ((TMath::Abs( particle->Charge() ) == 3))
                                {
                                	if (TrackRT<1.) pTGenTrans1->Fill(TrackPt);
                                        if (TrackRT>=1. && TrackRT<2.) pTGenTrans2->Fill(TrackPt);
                                        if (TrackRT>=2. && TrackRT<3.) pTGenTrans3->Fill(TrackPt);
                                        if (TrackRT>=3. && TrackRT<4.) pTGenTrans4->Fill(TrackPt);
					if (TrackRT>=4. && TrackRT<5.) pTGenTrans5->Fill(TrackPt);
					if (TrackRT>=5. && TrackRT<6.) pTGenTrans6->Fill(TrackPt);
					if (TrackRT>=6. && TrackRT<7.) pTGenTrans7->Fill(TrackPt);
					if (TrackRT>=7. && TrackRT<8.) pTGenTrans8->Fill(TrackPt);
					if (TrackRT>=8. && TrackRT<9.) pTGenTrans9->Fill(TrackPt);
					if (TrackRT>=9. && TrackRT<10.) pTGenTrans10->Fill(TrackPt);
                                        if (TrackRT>=10.) pTGenTrans11->Fill(TrackPt);

					pTGenTrans->Fill(TrackPt); // all particles in transverse
                                }
				if (isINEL0Rec)
          			{

					if ((TMath::Abs( particle->Charge() ) == 3))
					{
              					if (TrackRT<1.) sigLossTrigINEL01->Fill(TrackPt); // this is the same quantity to be used for inclusive tracking efficiency!
                        			if (TrackRT>=1. && TrackRT<2.) sigLossTrigINEL02->Fill(TrackPt);
						if (TrackRT>=2. && TrackRT<3.) sigLossTrigINEL03->Fill(TrackPt);
                        			if (TrackRT>=3. && TrackRT<4.) sigLossTrigINEL04->Fill(TrackPt);
						if (TrackRT>=4. && TrackRT<5.) sigLossTrigINEL05->Fill(TrackPt);
						if (TrackRT>=5. && TrackRT<6.) sigLossTrigINEL06->Fill(TrackPt);
						if (TrackRT>=6. && TrackRT<7.) sigLossTrigINEL07->Fill(TrackPt);
						if (TrackRT>=7. && TrackRT<8.) sigLossTrigINEL08->Fill(TrackPt);
						if (TrackRT>=8. && TrackRT<9.) sigLossTrigINEL09->Fill(TrackPt);
						if (TrackRT>=9. && TrackRT<10.) sigLossTrigINEL010->Fill(TrackPt);
						if (TrackRT>=10.) sigLossTrigINEL011->Fill(TrackPt);
					}
				}
				// redundant if no particle composition
				//if ( !( TMath::Abs( mcParticle->GetPDG()->Charge() ) == 3 ) ) continue;
				if ( isINEL0True && fisMCvtxInZcut) {
					if (TrackRT<1.) sigLossTrueINEL01->Fill(TrackPt);
                                	if (TrackRT>=1. && TrackRT<2.) sigLossTrueINEL02->Fill(TrackPt);
                                	if (TrackRT>=2. && TrackRT<3.) sigLossTrueINEL03->Fill(TrackPt);
                                	if (TrackRT>=3. && TrackRT<4.) sigLossTrueINEL04->Fill(TrackPt);
					if (TrackRT>=4. && TrackRT<5.) sigLossTrueINEL05->Fill(TrackPt);
					if (TrackRT>=5. && TrackRT<6.) sigLossTrueINEL06->Fill(TrackPt);
					if (TrackRT>=6. && TrackRT<7.) sigLossTrueINEL07->Fill(TrackPt);
					if (TrackRT>=7. && TrackRT<8.) sigLossTrueINEL08->Fill(TrackPt);
					if (TrackRT>=8. && TrackRT<9.) sigLossTrueINEL09->Fill(TrackPt);
					if (TrackRT>=9. && TrackRT<10.) sigLossTrueINEL010->Fill(TrackPt);
                                	if (TrackRT>=10.) sigLossTrueINEL011->Fill(TrackPt);
				}
				if ( fisPS ) {
                                        if (TrackRT<1.) TrigINEL0->Fill(1);
                                        if (TrackRT>=1. && TrackRT<2.) TrigINEL0->Fill(2);
                                        if (TrackRT>=2. && TrackRT<3.) TrigINEL0->Fill(3);
                                        if (TrackRT>=3. && TrackRT<4.) TrigINEL0->Fill(4);
					if (TrackRT>=4. && TrackRT<5.) TrigINEL0->Fill(5);
					if (TrackRT>=5. && TrackRT<6.) TrigINEL0->Fill(6);
					if (TrackRT>=6. && TrackRT<7.) TrigINEL0->Fill(7);
					if (TrackRT>=7. && TrackRT<8.) TrigINEL0->Fill(8);
					if (TrackRT>=8. && TrackRT<9.) TrigINEL0->Fill(9);
					if (TrackRT>=9. && TrackRT<10.) TrigINEL0->Fill(10);
					if (TrackRT>=10.) TrigINEL0->Fill(11);
				}
				if (isINEL0True) {
                                        if (TrackRT<1.) TrueINEL0->Fill(1);
                                        if (TrackRT>=1. && TrackRT<2.) TrueINEL0->Fill(2);
                                        if (TrackRT>=2. && TrackRT<3.) TrueINEL0->Fill(3);
                                        if (TrackRT>=3. && TrackRT<4.) TrueINEL0->Fill(4);
					if (TrackRT>=4. && TrackRT<5.) TrueINEL0->Fill(5);
					if (TrackRT>=5. && TrackRT<6.) TrueINEL0->Fill(6);
					if (TrackRT>=6. && TrackRT<7.) TrueINEL0->Fill(7);
					if (TrackRT>=7. && TrackRT<8.) TrueINEL0->Fill(8);
					if (TrackRT>=8. && TrackRT<9.) TrueINEL0->Fill(9);
					if (TrackRT>=9. && TrackRT<10.) TrueINEL0->Fill(10);
                                        if (TrackRT>=10.) TrueINEL0->Fill(11);
				}
			}

		        for(Int_t k = 0; k < toward->GetEntries(); k++)
                        {
                                AliVParticle* particle = (AliVParticle*)toward->At(k);
                                Float_t TrackPt  = particle->Pt();
				if ((TMath::Abs( particle->Charge() ) == 3))
                                {
                                	if (TrackRT<1.) pTGenTow1->Fill(TrackPt);
                                	if (TrackRT>=1. && TrackRT<2.) pTGenTow2->Fill(TrackPt);
                                	if (TrackRT>=2. && TrackRT<3.) pTGenTow3->Fill(TrackPt);
                                	if (TrackRT>=3. && TrackRT<4.) pTGenTow4->Fill(TrackPt);
					if (TrackRT>=4. && TrackRT<5.) pTGenTow5->Fill(TrackPt);
					if (TrackRT>=5. && TrackRT<6.) pTGenTow6->Fill(TrackPt);
					if (TrackRT>=6. && TrackRT<7.) pTGenTow7->Fill(TrackPt);
					if (TrackRT>=7. && TrackRT<8.) pTGenTow8->Fill(TrackPt);
					if (TrackRT>=8. && TrackRT<9.) pTGenTow9->Fill(TrackPt);
					if (TrackRT>=9. && TrackRT<10.) pTGenTow10->Fill(TrackPt);
					if (TrackRT>=10.) pTGenTow11->Fill(TrackPt);

                                	pTGenTow->Fill(TrackPt); // all particles in toward
                        	}
			}


                        for(Int_t k = 0; k < away->GetEntries(); k++)
                        {
                                AliVParticle* particle = (AliVParticle*)away->At(k);
                                Float_t TrackPt  = particle->Pt();
                                if ((TMath::Abs( particle->Charge() ) == 3))
                                {
                                        if (TrackRT<1.) pTGenAw1->Fill(TrackPt);
                                        if (TrackRT>=1. && TrackRT<2.) pTGenAw2->Fill(TrackPt);
                                        if (TrackRT>=2. && TrackRT<3.) pTGenAw3->Fill(TrackPt);
                                        if (TrackRT>=3. && TrackRT<4.) pTGenAw4->Fill(TrackPt);
					if (TrackRT>=4. && TrackRT<5.) pTGenAw5->Fill(TrackPt);
					if (TrackRT>=5. && TrackRT<6.) pTGenAw6->Fill(TrackPt);
					if (TrackRT>=6. && TrackRT<7.) pTGenAw7->Fill(TrackPt);
					if (TrackRT>=7. && TrackRT<8.) pTGenAw8->Fill(TrackPt);
					if (TrackRT>=8. && TrackRT<9.) pTGenAw9->Fill(TrackPt);
					if (TrackRT>=9. && TrackRT<10.) pTGenAw10->Fill(TrackPt);
					if (TrackRT>=10.) pTGenAw11->Fill(TrackPt);

					pTGenAw->Fill(TrackPt); // all particles in away
                                }
                        }
		}
      	}


	//______________________________________________________________
	// reconstructed MC stuff

	fRun  = esdEvent->GetRunNumber();
        fEventId = 0;
        if(esdEvent->GetHeader()) fEventId = GetEventIdAsLong(esdEvent->GetHeader());

	const Int_t nESDTracks = fESD->GetNumberOfTracks();

        TObjArray* fCTSTracks = new TObjArray();
        TrackRT  = 0.;
        nRecTracks = 0;
	AliVParticle* LeadingReco = 0;
	TObjArray *regionsMinMaxRECO = 0x0;

        for(Int_t iT = 0; iT < nESDTracks; iT++)
        {
        	AliESDtrack* esdTrack = esdEvent->GetTrack(iT);
        	if(!esdTrack) continue;
        	Double_t eta = esdTrack->Eta();
        	Double_t pt  = esdTrack->Pt();

        	if ( TMath::Abs(eta) > fEtaCut) continue;
        	if ( !TMath::Abs(pt) > 0.15) continue;

        	Int_t mcLabel = TMath::Abs(esdTrack->GetLabel());
        	AliVParticle *mcParticle = fMCEvent->GetTrack(mcLabel);

        	if(!mcParticle) {
       	        	printf("ERROR: mcParticle not available-------\n"); \
              		continue;
            	}

        	eta = mcParticle->Eta(); // generated eta and pT used intead of recontructed
        	pt = mcParticle->Pt();
        	if ( TMath::Abs(eta) > fEtaCut) continue;
        	if ( !TMath::Abs(pt) > 0.15) continue;

        	Int_t partPDG = TMath::Abs(mcParticle->PdgCode());

        	if ((TMath::Abs(mcParticle->Charge()) == 3) ) // only for charged particles
        	{

        		for ( int i = 0 ; i < 1 ; i++ ) // run over all the track cuts, i == 0 standard
          		{
            			UInt_t selectDebug = 0;
        			if (fTrackFilter[i])
        			{
          				selectDebug = fTrackFilter[i]->IsSelected(esdTrack);
          				if (!selectDebug)
          				continue;
        			}

         		// fill tracks array
        		if (!mcParticle) continue;
			fCTSTracks->Add(mcParticle);
         		}
    		}
   	}

	// leading object
        TObjArray *LeadingTrackReco = FindLeadingObjects(fCTSTracks);
    	if (LeadingTrackReco) {
                LeadingReco = (AliVParticle*)LeadingTrackReco->At(0);
                // Sorting
                TObjArray *regionSortedParticlesRECO = SortRegions((AliVParticle*)LeadingTrackReco->At(0), fCTSTracks);
                // Taking transverse regions
                regionsMinMaxRECO = GetMinMaxRegion((TList*)regionSortedParticlesRECO->At(2),(TList*)regionSortedParticlesRECO->At(3));


                TList *listMax = (TList*)regionsMinMaxRECO->At(0);
                TList *listMin = (TList*)regionsMinMaxRECO->At(1);
                TList *transverse = new TList();

                // fill transverse list
                for(Int_t j = 0; j < listMax->GetEntries(); j++){
                   	AliVParticle* particle1 = (AliVParticle*)listMax->At(j);
                	transverse->Add(particle1);
                }

                for(Int_t j = 0; j < listMin->GetEntries(); j++){
                        AliVParticle* particle2 = (AliVParticle*)listMin->At(j);
                        transverse->Add(particle2);
                }

        	Float_t LeadingPt = LeadingReco->Pt();
        	// select only plateau region
                if(LeadingPt > fLeadMin && LeadingPt < 300.){
			nRecTracks = transverse->GetEntries();
			TrackRT = nRecTracks/fAveRecMultiInTrans;

			fhRTReco->Fill(TrackRT);
			for(Int_t k = 0; k < nRecTracks; k++){
				AliVParticle* particle = (AliVParticle*)transverse->At(k);
				if (particle->IsPhysicalPrimary()) {
                        	        Float_t TrackPt  = particle->Pt();
	                                if (isINEL0Rec)
                                	{
						if (TrackRT<1.) effcomputationRec1->Fill(TrackPt); //  VZ: this is the same quantity to be used for inclusive tracking efficiency!
                        	        	if (TrackRT>=1. && TrackRT<2.) effcomputationRec2->Fill(TrackPt);
                        	        	if (TrackRT>=2. && TrackRT<3.) effcomputationRec3->Fill(TrackPt);
                        	        	if (TrackRT>=3. && TrackRT<4.) effcomputationRec4->Fill(TrackPt);
						if (TrackRT>=4. && TrackRT<5.) effcomputationRec5->Fill(TrackPt);
						if (TrackRT>=5. && TrackRT<6.) effcomputationRec6->Fill(TrackPt);
						if (TrackRT>=6. && TrackRT<7.) effcomputationRec7->Fill(TrackPt);
						if (TrackRT>=7. && TrackRT<8.) effcomputationRec8->Fill(TrackPt);
						if (TrackRT>=8. && TrackRT<9.) effcomputationRec9->Fill(TrackPt);
						if (TrackRT>=9. && TrackRT<10.) effcomputationRec10->Fill(TrackPt);
                        	        	if (TrackRT>=10.) effcomputationRec11->Fill(TrackPt);
					}
				}
			}
		}
	}

	if (LeadingMC && LeadingReco) FillRTResponseMatrix(LeadingMC, LeadingReco, (TList*)regionsMinMaxMC->At(0), (TList*)regionsMinMaxRECO->At(0), (TList*)regionsMinMaxMC->At(1), (TList*)regionsMinMaxRECO->At(1));
}




//________________________________________________________________________
void AliAnalysisTaskUeSpectraRT::AnalyzeMC(AliMCEvent* fMCEvent)
{
	fEvents->Fill(0);

	fNchargedTrue = 0;
        // VZ
	for ( int iT = 0 ; iT < fMCEvent->GetNumberOfTracks(); iT++ ) // loop over TRUE MC
	{
		// VZ
          	TParticle *mcParticle = fMCEvent->GetTrack(iT)->Particle();

          	if (!mcParticle)
	 	{
	     		cout<<"no mcParticle"<<endl;
	     		continue;
	  	}
          	// VZ
     	  	if(!(fMCEvent->IsPhysicalPrimary(iT))) {
          	continue;
          	}

	  	Double_t eta = mcParticle->Eta();
	  	Double_t pt = mcParticle->Pt();
	 	Int_t partPDG = TMath::Abs(mcParticle->GetPdgCode());
	  	if (TMath::Abs(eta) > fEtaCut) {
            		continue;
          	}
	  	if ((TMath::Abs( mcParticle->GetPDG()->Charge()) == 3 )) fNchargedTrue++;
	  	if ( !TMath::Abs(pt) > 0.15 ){
            		continue;
          	}

		// online reweighting for particle composition
		 Double_t weight = 1.;
		//TODO VZ: this is still not working
		//Double_t weight = fMCSpectraWeights->GetMCSpectraWeight(mcParticle ,fMCEvent);
    		// printf("Got weight factor %lf for pid: %d\n", weight, mcParticle->GetPdgCode());

		if ((TMath::Abs( mcParticle->GetPDG()->Charge() ) == 3)) pTGen->Fill(pt);
		if (isINEL0Rec)
	  	{
	    		if ((TMath::Abs( mcParticle->GetPDG()->Charge() ) == 3))
	    		{
	      			effcomputationGen->Fill(pt, weight); //inclusive
	      			// VZ: this is not needed if I use data driven correction
	      			/*if (partPDG == 211)  effcomputationGenPi->Fill(pt, weight); //pions
	      			else if (partPDG == 321)  effcomputationGenK->Fill(pt, weight);//kaons
	      			else if (partPDG == 2212) effcomputationGenP->Fill(pt, weight); //protons
	      			else if (partPDG == 3112) effcomputationGenSm->Fill(pt, weight);//sigma-
	      			else if (partPDG == 3222) effcomputationGenSp->Fill(pt, weight);//sigma+
	      			else if (partPDG == 3334) effcomputationGenO->Fill(pt, weight);//omega-
              			else if (partPDG == 3312) effcomputationGenXi->Fill(pt, weight);//Xi-
              			else effcomputationGenRest->Fill(pt, weight);//rest
	    			*/
			}
	    		// else if (partPDG == 3122) effcomputationGenL->Fill(pt, weight); //Lambda
	    	}

	   	if ( !( TMath::Abs( mcParticle->GetPDG()->Charge() ) == 3 ) ) continue;

	   	if ( isINEL0Rec ) sigLossTrigINEL0->Fill(pt);
	   	if ( isINEL0True && fisMCvtxInZcut) sigLossTrueINEL0->Fill(pt);
	}

	if (isINEL0Rec) nchtrue->Fill(fNchargedTrue);

	fEvents->Fill(1);

}

//________________________________________________________________________
void AliAnalysisTaskUeSpectraRT::AnalyseDataRT(AliESDEvent* esdEvent)
{

        fRun  = esdEvent->GetRunNumber();
        fEventId = 0;
        if(esdEvent->GetHeader()) fEventId = GetEventIdAsLong(esdEvent->GetHeader());

	const Int_t nESDTracks = fESD->GetNumberOfTracks();

        TObjArray* fCTSTracks = new TObjArray();
        Double_t TrackRT  = 0.;
        Double_t nRecTracks = 0;

	for(Int_t iT = 0; iT < nESDTracks; iT++)
        {
		AliVParticle* part = fESD->GetTrack(iT);
                Double_t eta  = part->Eta();
                Double_t pt   = part->Pt();

                if( TMath::Abs(eta) > fEtaCut )
                continue;
                if( !TMath::Abs(pt) > 0.15 )
                continue;

                // VZ: Selecting only default track filter for the moment (too heavy to loop over all 18 masks)
                for ( int i = 0 ; i < 1 ; i++ )
                {
			UInt_t selectDebug = 0;
                        if (fTrackFilter[i])
                        {
                        	selectDebug = fTrackFilter[i]->IsSelected(part);
                                if (!selectDebug)
                                {
                                	continue;
                                }
                       	}
                        // fill tracks array
                       	fCTSTracks->Add(part);
                        if (!part) continue;
                }
	}

       // leading object
        TObjArray *LeadingTrackReco = FindLeadingObjects(fCTSTracks);
	if (LeadingTrackReco) {
		AliVParticle* LeadingReco = 0;
                LeadingReco = (AliVParticle*)LeadingTrackReco->At(0);
                // Sorting
                TObjArray *regionSortedParticlesRECO = SortRegions((AliVParticle*)LeadingTrackReco->At(0), fCTSTracks);
                // Taking transverse regions
                TObjArray *regionsMinMaxRECO = GetMinMaxRegion((TList*)regionSortedParticlesRECO->At(2),(TList*)regionSortedParticlesRECO->At(3));

	    	// Taking toward region
	    	TObjArray *towardRECO = GetRegionAwTow((TList*)regionSortedParticlesRECO->At(0));

                // Taking away region
                TObjArray *awayRECO = GetRegionAwTow((TList*)regionSortedParticlesRECO->At(1));

		INEL0->Fill(12); // filling histogram for events to normalise NoLead
		// for spectra without leading particle
		for(Int_t i=0; i<4; i++) {
   			for(Int_t j=0; j<((TList*)regionSortedParticlesRECO->At(i))->GetEntries(); j++){
				AliVParticle* allNoLead = ((AliVParticle*)((TList*)regionSortedParticlesRECO->At(i))->At(j));
				if(allNoLead != LeadingReco) {
					Float_t TrackPtNoLead  = allNoLead->Pt();
					ptiNoLead[i]->Fill(TrackPtNoLead);
				}
			}
		}


                TList *listMax = (TList*)regionsMinMaxRECO->At(0);
                TList *listMin = (TList*)regionsMinMaxRECO->At(1);
                TList *listTow = (TList*) towardRECO->At(0);
                TList *listAw = (TList*) awayRECO->At(0);

		TList *transverse = new TList();
		TList * toward = new TList();
		TList * away = new TList();

                // fill transverse list
                for(Int_t j = 0; j < listMax->GetEntries(); j++){
                	AliVParticle* particle1 = (AliVParticle*)listMax->At(j);
                        transverse->Add(particle1);
                }

                for(Int_t j = 0; j < listMin->GetEntries(); j++){
                        AliVParticle* particle2 = (AliVParticle*)listMin->At(j);
                        transverse->Add(particle2);
                }

		// fill toward list
  		for(Int_t j = 0; j < listTow->GetEntries(); j++){
   			AliVParticle* particleTow = (AliVParticle*)listTow->At(j);
   			toward->Add(particleTow);
   		}

                // fill away list
                for(Int_t j = 0; j < listAw->GetEntries(); j++){
                        AliVParticle* particleAw = (AliVParticle*)listAw->At(j);
                        away->Add(particleAw);
                }


		Float_t LeadingPt = LeadingReco->Pt();
		Int_t nTotTracks = 0;

		// select only plateau region
                if(LeadingPt > fLeadMin && LeadingPt < 300.){
			nRecTracks = transverse->GetEntries();
			TrackRT = nRecTracks/fAveMultiInTrans;

			fhRTData->Fill(TrackRT);

			if (TrackRT<1.) INEL0->Fill(1);
                        if (TrackRT>=1. && TrackRT<2.) INEL0->Fill(2);
                        if (TrackRT>=2. && TrackRT<3.) INEL0->Fill(3);
                        if (TrackRT>=3. && TrackRT<4.) INEL0->Fill(4);
			if (TrackRT>=4. && TrackRT<5.) INEL0->Fill(5);
			if (TrackRT>=5. && TrackRT<6.) INEL0->Fill(6);
			if (TrackRT>=6. && TrackRT<7.) INEL0->Fill(7);
			if (TrackRT>=7. && TrackRT<8.) INEL0->Fill(8);
			if (TrackRT>=8. && TrackRT<9.) INEL0->Fill(9);
			if (TrackRT>=9. && TrackRT<10.) INEL0->Fill(10);
			if (TrackRT>=10. ) INEL0->Fill(11);

			if (TrackRT<1.) ptli1->Fill(LeadingPt);
			if (TrackRT>=1. && TrackRT<2.) ptli2->Fill(LeadingPt);
			if (TrackRT>=2. && TrackRT<3.) ptli3->Fill(LeadingPt);
			if (TrackRT>=3. && TrackRT<4.) ptli4->Fill(LeadingPt);
			if (TrackRT>=4. && TrackRT<5.) ptli5->Fill(LeadingPt);
			if (TrackRT>=5. && TrackRT<6.) ptli6->Fill(LeadingPt);
			if (TrackRT>=6. && TrackRT<7.) ptli7->Fill(LeadingPt);
			if (TrackRT>=7. && TrackRT<8.) ptli8->Fill(LeadingPt);
			if (TrackRT>=8. && TrackRT<9.) ptli9->Fill(LeadingPt);
			if (TrackRT>=9. && TrackRT<10.) ptli10->Fill(LeadingPt);
			if (TrackRT>=10.) ptli11->Fill(LeadingPt);
			ptliMB->Fill(LeadingPt);

			for(Int_t k = 0; k < nRecTracks; k++)
			  {
			    AliVParticle* particle = (AliVParticle*)transverse->At(k);
			    Float_t TrackPt  = particle->Pt();
			    if (TrackRT<1.) pti1Trans[0]->Fill(TrackPt);
			    if (TrackRT>=1. && TrackRT<2.) pti2Trans[0]->Fill(TrackPt);
			    if (TrackRT>=2. && TrackRT<3.) pti3Trans[0]->Fill(TrackPt);
			    if (TrackRT>=3. && TrackRT<4.) pti4Trans[0]->Fill(TrackPt);
			    if (TrackRT>=4. && TrackRT<5.) pti5Trans[0]->Fill(TrackPt);
			    if (TrackRT>=5. && TrackRT<6.) pti6Trans[0]->Fill(TrackPt);
			    if (TrackRT>=6. && TrackRT<7.) pti7Trans[0]->Fill(TrackPt);
			    if (TrackRT>=7. && TrackRT<8.) pti8Trans[0]->Fill(TrackPt);
			    if (TrackRT>=8. && TrackRT<9.) pti9Trans[0]->Fill(TrackPt);
			    if (TrackRT>=9. && TrackRT<10.) pti10Trans[0]->Fill(TrackPt);
			    if (TrackRT>=10.) pti11Trans[0]->Fill(TrackPt);

			    ptiMBTrans[0]->Fill(TrackPt); // all particles in transverse
			  }

			for(Int_t k = 0; k < toward->GetEntries(); k++)
			  {
			    AliVParticle* particle = (AliVParticle*)toward->At(k);
			    Float_t TrackPt  = particle->Pt();
			    if (TrackRT<1.) pti1Tow[0]->Fill(TrackPt);
			    if (TrackRT>=1. && TrackRT<2.) pti2Tow[0]->Fill(TrackPt);
			    if (TrackRT>=2. && TrackRT<3.) pti3Tow[0]->Fill(TrackPt);
			    if (TrackRT>=3. && TrackRT<4.) pti4Tow[0]->Fill(TrackPt);
			    if (TrackRT>=4. && TrackRT<5.) pti5Tow[0]->Fill(TrackPt);
			    if (TrackRT>=5. && TrackRT<6.) pti6Tow[0]->Fill(TrackPt);
			    if (TrackRT>=6. && TrackRT<7.) pti7Tow[0]->Fill(TrackPt);
			    if (TrackRT>=7. && TrackRT<8.) pti8Tow[0]->Fill(TrackPt);
			    if (TrackRT>=8. && TrackRT<9.) pti9Tow[0]->Fill(TrackPt);
			    if (TrackRT>=9. && TrackRT<10.) pti10Tow[0]->Fill(TrackPt);
			    if (TrackRT>=10.) pti11Tow[0]->Fill(TrackPt);

			    ptiMBTow[0]->Fill(TrackPt); // all particles in toward
			  }

                        for(Int_t k = 0; k < away->GetEntries(); k++)
			  {
                                AliVParticle* particle = (AliVParticle*)away->At(k);
                                Float_t TrackPt  = particle->Pt();
                                if (TrackRT<1.) pti1Aw[0]->Fill(TrackPt);
                                if (TrackRT>=1. && TrackRT<2.) pti2Aw[0]->Fill(TrackPt);
                                if (TrackRT>=2. && TrackRT<3.) pti3Aw[0]->Fill(TrackPt);
                                if (TrackRT>=3. && TrackRT<4.) pti4Aw[0]->Fill(TrackPt);
				if (TrackRT>=4. && TrackRT<5.) pti5Aw[0]->Fill(TrackPt);
				if (TrackRT>=5. && TrackRT<6.) pti6Aw[0]->Fill(TrackPt);
				if (TrackRT>=6. && TrackRT<7.) pti7Aw[0]->Fill(TrackPt);
				if (TrackRT>=7. && TrackRT<8.) pti8Aw[0]->Fill(TrackPt);
				if (TrackRT>=8. && TrackRT<9.) pti9Aw[0]->Fill(TrackPt);
				if (TrackRT>=9. && TrackRT<10.) pti10Aw[0]->Fill(TrackPt);
				if (TrackRT>=10.) pti11Aw[0]->Fill(TrackPt);

                                ptiMBAw[0]->Fill(TrackPt); // all particles in away
                        }

		}
	}

}

//________________________________________________________________________
void AliAnalysisTaskUeSpectraRT::CorrectionsDataRT(AliESDEvent* esdEvent, Bool_t isVtxGood)
{

	fRun  = esdEvent->GetRunNumber();
        fEventId = 0;
        if(esdEvent->GetHeader()) fEventId = GetEventIdAsLong(esdEvent->GetHeader());

        const Int_t nESDTracks = fESD->GetNumberOfTracks();

        TObjArray* fCTSTracks = new TObjArray();
        Double_t TrackRT  = 0.;
        Double_t nRecTracks = 0;
        for(Int_t iT = 0; iT < nESDTracks; iT++)
        {
                AliVParticle* part = fESD->GetTrack(iT);
                Double_t eta  = part->Eta();
                Double_t pt   = part->Pt();

                if( TMath::Abs(eta) > fEtaCut )
                continue;
		if( !TMath::Abs(pt) > 0.15 )
                continue;

                // VZ: Selecting only default track filter for the moment (too heavy to loop over all 18 masks)
                for ( int i = 0 ; i < 1 ; i++ )
                {
                        UInt_t selectDebug = 0;
                        if (fTrackFilter[i])
                        {
                                selectDebug = fTrackFilter[i]->IsSelected(part);
                                if (!selectDebug)
                                {
                                        continue;
                                }
                        }
                        // fill tracks array
                        fCTSTracks->Add(part);
                        if (!part) continue;
		}
	}
        // leading object
        TObjArray *LeadingTrackReco = FindLeadingObjects(fCTSTracks);
        if (LeadingTrackReco){
		AliVParticle* LeadingReco = 0;
                LeadingReco = (AliVParticle*)LeadingTrackReco->At(0);

                // Sorting
                TObjArray *regionSortedParticlesRECO = SortRegions((AliVParticle*)LeadingTrackReco->At(0), fCTSTracks);

                // Taking transverse regions
                TObjArray *regionsMinMaxRECO = GetMinMaxRegion((TList*)regionSortedParticlesRECO->At(2),(TList*)regionSortedParticlesRECO->At(3));

                TList *listMax = (TList*)regionsMinMaxRECO->At(0);
                TList *listMin = (TList*)regionsMinMaxRECO->At(1);
                TList *transverse = new TList();

                // fill transverse list
                for(Int_t j = 0; j < listMax->GetEntries(); j++){
                	AliVParticle* particle1 = (AliVParticle*)listMax->At(j);
                        transverse->Add(particle1);
                }

                for(Int_t j = 0; j < listMin->GetEntries(); j++){
                	AliVParticle* particle2 = (AliVParticle*)listMin->At(j);
                        transverse->Add(particle2);
                }

                Float_t LeadingPt = LeadingReco->Pt();

                // select only plateau region
                if(LeadingPt > fLeadMin && LeadingPt < 300.){
                	nRecTracks = transverse->GetEntries();
                        TrackRT = nRecTracks/fAveMultiInTrans;

			if (TrackRT<1.) fPS->Fill(1);
                        if (TrackRT>=1. && TrackRT<2.) fPS->Fill(2);
                        if (TrackRT>=2. && TrackRT<3.) fPS->Fill(3);
                        if (TrackRT>=3. && TrackRT<4.) fPS->Fill(4);
			if (TrackRT>=4. && TrackRT<5.) fPS->Fill(5);
			if (TrackRT>=5. && TrackRT<6.) fPS->Fill(6);
			if (TrackRT>=6. && TrackRT<7.) fPS->Fill(7);
			if (TrackRT>=7. && TrackRT<8.) fPS->Fill(8);
			if (TrackRT>=8. && TrackRT<9.) fPS->Fill(9);
			if (TrackRT>=9. && TrackRT<10.) fPS->Fill(10);
			if (TrackRT>=10.) fPS->Fill(11);

			fPS->Fill(12); // filling the inclusive case
                 }
	}

	if (isVtxGood){

		const Int_t nESDTracks_vtx = fESD->GetNumberOfTracks();

	        TObjArray* fCTSTracks_vtx = new TObjArray();
	        Double_t TrackRT_vtx  = 0.;
	        Double_t nRecTracks_vtx = 0;
	        for(Int_t iT = 0; iT < nESDTracks_vtx; iT++)
	        {
	                AliVParticle* part_vtx = fESD->GetTrack(iT);
	                Double_t eta_vtx  = part_vtx->Eta();
	                Double_t pt_vtx   = part_vtx->Pt();

	                if( TMath::Abs(eta_vtx) > fEtaCut )
	                continue;
	                if( !TMath::Abs(pt_vtx) > 0.15 )
	                continue;

	                // VZ: Selecting only default track filter for the moment (too heavy to loop over all 18 masks)
	                for ( int i = 0 ; i < 1 ; i++ )
	                {
	                        UInt_t selectDebug_vtx = 0;
	                        if (fTrackFilter[i])
	                        {
	                                selectDebug_vtx = fTrackFilter[i]->IsSelected(part_vtx);
	                                if (!selectDebug_vtx)
	                                {
	                                        continue;
	                                }
	                        }
	                        // fill tracks array
	                        fCTSTracks_vtx->Add(part_vtx);
	                        if (!part_vtx) continue;
			}
		}
	        // leading object
	        TObjArray *LeadingTrackReco_vtx = FindLeadingObjects(fCTSTracks_vtx);
	        if(LeadingTrackReco_vtx) {
			AliVParticle* LeadingReco_vtx = 0;
	        	LeadingReco_vtx = (AliVParticle*)LeadingTrackReco_vtx->At(0);

	                // Sorting
	                TObjArray *regionSortedParticlesRECO_vtx = SortRegions((AliVParticle*)LeadingTrackReco_vtx->At(0), fCTSTracks_vtx);

	                // Taking transverse regions
	                TObjArray *regionsMinMaxRECO_vtx = GetMinMaxRegion((TList*)regionSortedParticlesRECO_vtx->At(2),(TList*)regionSortedParticlesRECO_vtx->At(3));

	                TList *listMax_vtx = (TList*)regionsMinMaxRECO_vtx->At(0);
	                TList *listMin_vtx = (TList*)regionsMinMaxRECO_vtx->At(1);
	                TList *transverse_vtx = new TList();

	                // fill transverse list
	                for(Int_t j = 0; j < listMax_vtx->GetEntries(); j++){
	                	AliVParticle* particle1_vtx = (AliVParticle*)listMax_vtx->At(j);
	                        transverse_vtx->Add(particle1_vtx);
	                }

	                for(Int_t j = 0; j < listMin_vtx->GetEntries(); j++){
	                	AliVParticle* particle2_vtx = (AliVParticle*)listMin_vtx->At(j);
	                        transverse_vtx->Add(particle2_vtx);
	                }

	                Float_t LeadingPt_vtx = LeadingReco_vtx->Pt();
	                // select only plateau region
	                if(LeadingPt_vtx > fLeadMin && LeadingPt_vtx < 300.){
				nRecTracks_vtx = transverse_vtx->GetEntries();
	                        TrackRT_vtx = nRecTracks_vtx/fAveMultiInTrans;
	                        if (TrackRT_vtx<1.) fVtxPS->Fill(1);
	                        if (TrackRT_vtx>=1. && TrackRT_vtx<2.) fVtxPS->Fill(2);
	                        if (TrackRT_vtx>=2. && TrackRT_vtx<3.) fVtxPS->Fill(3);
	                        if (TrackRT_vtx>=3. && TrackRT_vtx<4.) fVtxPS->Fill(4);
				if (TrackRT_vtx>=4. && TrackRT_vtx<5.) fVtxPS->Fill(5);
				if (TrackRT_vtx>=5. && TrackRT_vtx<6.) fVtxPS->Fill(6);
				if (TrackRT_vtx>=6. && TrackRT_vtx<7.) fVtxPS->Fill(7);
				if (TrackRT_vtx>=7. && TrackRT_vtx<8.) fVtxPS->Fill(8);
				if (TrackRT_vtx>=8. && TrackRT_vtx<9.) fVtxPS->Fill(9);
				if (TrackRT_vtx>=9. && TrackRT_vtx<10.) fVtxPS->Fill(10);
				if (TrackRT_vtx>=10.) fVtxPS->Fill(11);

				fVtxPS->Fill(12); // filling inclusive case

                 	}
		}

	}

}


//________________________________________________________________________
void AliAnalysisTaskUeSpectraRT::AnalyzeESD(AliESDEvent* esdEvent)
{


	fRun  = esdEvent->GetRunNumber();
	fEventId = 0;
	if(esdEvent->GetHeader())
		fEventId = GetEventIdAsLong(esdEvent->GetHeader());
	fEvents->Fill(2);

	const Int_t nESDTracks = esdEvent->GetNumberOfTracks();

	if (fAnalysisMC) // for reconstructed MC
	{
	  for ( int iT = 0 ; iT < nESDTracks ; iT++ )
	  {
	    AliESDtrack* esdTrack = esdEvent->GetTrack(iT);
	    if(!esdTrack) continue;
	    Double_t eta = esdTrack->Eta();
	    Double_t pt  = esdTrack->Pt();
	    Double_t phi = esdTrack->Phi();
	    if ( TMath::Abs(eta) > fEtaCut) continue;
	    if ( !TMath::Abs(pt) > 0.15) continue;

	    Int_t mcLabel = TMath::Abs(esdTrack->GetLabel());
            TParticle *mcParticle = fMCEvent->GetTrack(mcLabel)->Particle();

	    if(!mcParticle) {
              printf("ERROR: mcParticle not available-------\n"); \
              continue;
            }

	    eta = mcParticle->Eta(); // generated eta and pT used intead of recontructed
	    pt = mcParticle->Pt();
	    phi = mcParticle->Phi();
	    if ( TMath::Abs(eta) > fEtaCut) continue;
	    if ( !TMath::Abs(pt) > 0.15) continue;

	    Int_t partPDG = TMath::Abs(mcParticle->GetPdgCode());

	    if ((TMath::Abs(mcParticle->GetPDG()->Charge()) == 3) ) // only for charged particles
	    {

		for ( int i = 0 ; i < 18 ; i++ ) // run over all the track cuts, i == 0 standard
		{
	        UInt_t selectDebug = 0;
		if (fTrackFilter[i])
		{
		  selectDebug = fTrackFilter[i]->IsSelected(esdTrack);
		  if (!selectDebug)
		  continue;
	   	}

		primariesTrackFilter[i]->Fill(pt); // "primaries" filtered with track filter
		Phi->Fill(pt, phi);
		Eta->Fill(pt, eta);
		EtaPhi->Fill(eta, phi);


 		// for online particle composition
 		// TODO VZ: this is still not working
 		Double_t weight = 1.;
 		/*if (fMCSpectraWeights) {
			Double_t weight = fMCSpectraWeights->GetMCSpectraWeight(mcParticle ,fMCEvent);
        		//printf("I got the fMCSpectraWeights");
		}
		else {
			//printf("Problem in getting fMCSpectraWeights");
			Double_t weight = 1.;
		}
		*/
                //printf("Got weight factor %lf for pid: %d\n", weight, mcParticle->GetPdgCode());

		if (!(fMCEvent->IsPhysicalPrimary(mcLabel))) // secondary particles
		{
		  secondaries[i]->Fill(pt);
		}
		else
		{
		  if (i==0) effcomputationRec->Fill(pt, weight); // rec pT for efficiency true primary particles
		// this is not needed if I use the data driven correction
		/*  if (i == 0) // for particle composition, only standar track cuts requires (i == 0)
		  {
		    if (partPDG == 211)  effcomputationRecPi->Fill(pt, weight); //pions
		    else if (partPDG == 321)  effcomputationRecK->Fill(pt, weight);//kaons
		    else if (partPDG == 2212) effcomputationRecP->Fill(pt, weight); //protons
		    else if (partPDG == 3112) effcomputationRecSm->Fill(pt, weight);//sigma-
		    else if (partPDG == 3222) effcomputationRecSp->Fill(pt, weight);//sigma+
		    else if (partPDG == 3334) effcomputationRecO->Fill(pt, weight);//omega-
		    else if (partPDG == 3312) effcomputationRecXi->Fill(pt, weight);//Xi-
		    else effcomputationRecRest->Fill(pt);//Rest
		  }

		*/
 		}
	      }

	    // for matching efficiency
	    UInt_t selectDebug = 0;
	    if (fTrackFilterMatchEff) {
                  selectDebug = fTrackFilterMatchEff->IsSelected(esdTrack);
                  if (!selectDebug)
                  continue;
	    }
	    primariesTrackFilterME->Fill(pt);

	    }
	    /*else // only for neutral particles (searching Lambdas)
	    {
	      UInt_t selectDebug = 0;
	      if (fTrackFilter[0])
	      {
	        selectDebug = fTrackFilter[0]->IsSelected(esdTrack);
		if (!selectDebug)
		continue;
	      }

   	      if (!(fMCEvent->IsPhysicalPrimary(mcLabel))) {
                continue;
 	      }

	      if (partPDG == 3122) effcomputationRecL->Fill(pt); //Lambda

	     }*/
	   }
	}
	else // for real data
	{
		for(Int_t iT = 0; iT < nESDTracks; iT++)
		{
			AliESDtrack* esdTrack = esdEvent->GetTrack(iT);

			Double_t eta  = esdTrack->Eta();
			Double_t pt   = esdTrack->Pt();
			Double_t phi  = esdTrack->Phi();

			if( TMath::Abs(eta) > fEtaCut )
			continue;
			if( !TMath::Abs(pt) > 0.15 )
			continue;

			for ( int i = 0 ; i < 18 ; i++ )
 			{
				UInt_t selectDebug = 0;
				if (fTrackFilter[i])
				{
					selectDebug = fTrackFilter[i]->IsSelected(esdTrack);

					if (!selectDebug)
           				  continue;
				}
				pti[i]->Fill(pt);
				Phi->Fill(pt, phi);
				Eta->Fill(pt, eta);
				EtaPhi->Fill(eta, phi);
			}
			UInt_t selectDebug = 0;
			// for matching efficiency calculation (use TPC only tracks)
			if (fTrackFilterMatchEff)
			{
				selectDebug = fTrackFilterMatchEff->IsSelected(esdTrack);

				if (!selectDebug)
					continue;
			}
			ptiME->Fill(pt);
		}//end of track loop
	}
	fEvents->Fill(3);
}

//__________________________________________________________________________
void AliAnalysisTaskUeSpectraRT::AnalyzeESDforDCA(AliESDEvent* esdEvent)
{

	fEvents->Fill(4);

	fRun  = esdEvent->GetRunNumber();
	fEventId = 0;
	if(esdEvent->GetHeader())
	fEventId = GetEventIdAsLong(esdEvent->GetHeader());


	const Int_t nESDTracks = esdEvent->GetNumberOfTracks();

	if (fAnalysisMC)
	{

	  for ( int iT = 0 ; iT < nESDTracks ; iT++ )
	  {

	  AliESDtrack* esdTrack = esdEvent->GetTrack(iT);
	  if(!esdTrack) continue;
	  //track cuts
	  UInt_t selectDebug = 0;
          if (fTrackFilterDCA)
	  {
	  selectDebug = fTrackFilterDCA->IsSelected(esdTrack);
	  if (!selectDebug) continue;
	  }

	  Double_t eta = esdTrack->Eta();
	  Double_t pt  = esdTrack->Pt();

	  if ( TMath::Abs(eta) > fEtaCut) {
            continue;
          }
	  if ( !TMath::Abs(pt) > 0.15){
            continue;
          }

	  Int_t mcLabel = TMath::Abs(esdTrack->GetLabel());
	  TParticle *mcParticle = fMCEvent->GetTrack(mcLabel)->Particle();
          if(!mcParticle)
          {
            printf("----ERROR: mcParticle not available------------------\n");
 	    continue;
          }

	  eta = mcParticle->Eta();
	  pt = mcParticle->Pt();

          if ( TMath::Abs(eta) > fEtaCut) {
            continue;
          }
	  if ( !TMath::Abs(pt) > 0.15) {
            continue;
          }
	  if (!(TMath::Abs(mcParticle->GetPDG()->Charge())  == 3) ){
            continue;
          }

	  esdTrack->GetImpactParameters(fdcaxy,fdcaz);

	  if (!(fMCEvent->IsPhysicalPrimary(mcLabel)))
	  {

            if (fMCEvent->IsSecondaryFromWeakDecay(mcLabel))
            {
	      ptvsdcaDecs->Fill(pt,fdcaxy);
	      ptvsdcacentralDecs->Fill(pt,fdcaxy);
	    }

	    if (fMCEvent->IsSecondaryFromMaterial(mcLabel))
	    {
	      ptvsdcaMatl->Fill(pt,fdcaxy);
	      ptvsdcacentralMatl->Fill(pt,fdcaxy);
	    }

	    continue;
	  }

          ptvsdcaPrim->Fill(pt,fdcaxy);
	  ptvsdcacentralPrim->Fill(pt,fdcaxy);

	}

	}
	else
	{
	  for(Int_t iT = 0; iT < nESDTracks; iT++)
	  {
	    AliESDtrack* esdTrack = esdEvent->GetTrack(iT);
 	    //only golden track cuts
	    UInt_t selectDebug = 0;
	    if (fTrackFilterDCA) {
	      selectDebug = fTrackFilterDCA->IsSelected(esdTrack);
	      if (!selectDebug) continue;
	    }

	    Double_t eta  = esdTrack->Eta();
	    Double_t pt   = esdTrack->Pt();

	    if( TMath::Abs(eta) > fEtaCut ) continue;
	    if( !TMath::Abs(pt) > 0.15 ) continue;

	    esdTrack->GetImpactParameters(fdcaxy,fdcaz);

	    ptvsdcaData->Fill(pt,fdcaxy);

	    ptvsdcacentralData->Fill(pt,fdcaxy);

	   }//end of track loop
	}

	fEvents->Fill(5);

}

//_____________________________________________________________________________
ULong64_t AliAnalysisTaskUeSpectraRT::GetEventIdAsLong(AliVHeader* header)
{
	// To have a unique id for each event in a run!
	// Modified from AliRawReader.h
	return ((ULong64_t)header->GetBunchCrossNumber()+
			(ULong64_t)header->GetOrbitNumber()*3564+
			(ULong64_t)header->GetPeriodNumber()*16777215*3564);
}

//______________________________________________________________________
Bool_t AliAnalysisTaskUeSpectraRT::selectVertex2015pp(AliESDEvent *esd,
			  Bool_t checkSPDres, //enable check on vtx resolution
			  Bool_t requireSPDandTrk, //ask for both trk and SPD vertex
			  Bool_t checkProximity) //apply cut on relative position of spd and trk verteces
{

  if (!esd) return kFALSE;

  const AliESDVertex * trkVertex = esd->GetPrimaryVertexTracks();
  const AliESDVertex * spdVertex = esd->GetPrimaryVertexSPD();
  Bool_t hasSPD = spdVertex->GetStatus();
  Bool_t hasTrk = trkVertex->GetStatus();

  //Note that AliVertex::GetStatus checks that N_contributors is > 0
  //reject events if both are explicitly requested and none is available
  if (requireSPDandTrk && !(hasSPD && hasTrk)) return kFALSE;

  //reject events if none between the SPD or track verteces are available
  //if no trk vertex, try to fall back to SPD vertex;
  if (!hasTrk) {
    if (!hasSPD) return kFALSE;
    //on demand check the spd vertex resolution and reject if not satisfied
    if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
  } else {
    if (hasSPD) {
      //if enabled check the spd vertex resolution and reject if not satisfied
      //if enabled, check the proximity between the spd vertex and trak vertex, and reject if not satisfied
      if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
      if ((checkProximity && TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())>0.5)) return kFALSE;
    }
  }
return kTRUE;
}

//____________________________________________________________________________
Bool_t AliAnalysisTaskUeSpectraRT::IsGoodSPDvertexRes(const AliESDVertex * spdVertex)
{
  if (!spdVertex) return kFALSE;
  if (spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25)) return kFALSE;
  return kTRUE;
}

//__________________________________________________________________________
Bool_t AliAnalysisTaskUeSpectraRT::isMCEventTrueINEL0(AliMCEvent* fMCEvent)
{
  Bool_t isINEL0 = kFALSE;
  for ( int iT = 0 ; iT < fMCEvent->GetNumberOfTracks(); iT++ )
  {
    TParticle *mcParticle = fMCEvent->GetTrack(iT)->Particle();

    if (!mcParticle)
    {
      cout<<"no mcParticle"<<endl;
      continue;
    }

    if(!fMCEvent->IsPhysicalPrimary(iT))
      continue;

    if(!(mcParticle->Pt()>0.0))
      continue;

    Double_t eta = mcParticle->Eta();
    if ( TMath::Abs(eta) > 1.0 )
      continue;

    if ( !( TMath::Abs( mcParticle->GetPDG()->Charge() ) == 3 ) )
      continue;

    isINEL0 = kTRUE;
    break;

  }


  return isINEL0;

}


//___________________________________________________
TObjArray* AliAnalysisTaskUeSpectraRT::FindLeadingObjects(TObjArray *array)
{
  if(!array) return 0;

  Int_t nTracks = array->GetEntriesFast();
  if(!nTracks) return 0;

  TObjArray* tracks = new TObjArray(nTracks);
  for(Int_t ipart = 0; ipart < nTracks; ++ipart){
   AliVParticle* part = (AliVParticle*)(array->At(ipart));

  if(!part) continue;
   tracks->AddLast(part);
   }
  QSortTracks( *tracks, 0, tracks->GetEntriesFast() );

  nTracks = tracks->GetEntriesFast();
  if( !nTracks ) return 0;
  return tracks;
 }


//__________________________________________________
void AliAnalysisTaskUeSpectraRT::QSortTracks(TObjArray &a, Int_t first, Int_t last)
{
// Sort array of TObjArray of tracks by Pt using a quicksort algorithm.
  static TObject *tmp;
  static int i;           // "static" to save stack space
  int j;

  while (last - first > 1) {
    i = first;
    j = last;
    for (;;) {
      while (++i < last && ((AliVParticle*)a[i])->Pt() > ((AliVParticle*)a[first])->Pt() )
        ;
      while (--j > first && ((AliVParticle*)a[j])->Pt() < ((AliVParticle*)a[first])->Pt() )
        ;
      if (i >= j)
        break;

      tmp  = a[i];
      a[i] = a[j];
      a[j] = tmp;
    }
    if (j == first) {
      ++first;
      continue;
    }
    tmp = a[first];
    a[first] = a[j];
    a[j] = tmp;
    if (j - first < last - (j + 1)) {
      QSortTracks(a, first, j);
      first = j + 1;
    } else {
      QSortTracks(a, j + 1, last);
      last = j;
    }
  }
}


//___________________________________________________
TObjArray* AliAnalysisTaskUeSpectraRT::SortRegions(const AliVParticle* leading, TObjArray *array)
{
  if(!array) return 0;

// Assign particles to towards, away or transverse regions.
// Returns a lists of particles for each region.
  static const Double_t k60rad  = 60.*TMath::Pi()/180.;
  static const Double_t k120rad = 120.*TMath::Pi()/180.;

// Define output lists of particles
  TList *toward = new TList();
  TList *away = new TList();
  TList *transverse1 = new TList();
  TList *transverse2 = new TList();

  TObjArray *regionParticles = new TObjArray;
  regionParticles->SetOwner(kTRUE);

  regionParticles->AddLast(toward);
  regionParticles->AddLast(away);
  regionParticles->AddLast(transverse1);
  regionParticles->AddLast(transverse2);

  if(!leading)
   return regionParticles;

// Switch to vector for leading particle
  TVector3 leadVect(leading->Px(),leading->Py(),leading->Pz());

  Int_t nTracks = array->GetEntriesFast();
  if(!nTracks) return 0;

// Loop over tracks
  for(Int_t ipart = 0; ipart < nTracks; ++ipart){
   AliVParticle* part = (AliVParticle*)(array->At(ipart));
   if(!part) continue;

// Switch to vectors for particles
   TVector3 partVect(part->Px(), part->Py(), part->Pz());

   Int_t region = 0;

   Float_t deltaPhi = leadVect.DeltaPhi(partVect);
   if(deltaPhi <= -TMath::PiOver2()) deltaPhi+=TMath::TwoPi();
   if(deltaPhi > 3*TMath::PiOver2()) deltaPhi-=TMath::TwoPi();

   Double_t fUeDeltaPhiMinCut     = TMath::DegToRad()*60.;
   Double_t fUeDeltaPhiMaxCut     = TMath::DegToRad()*120.;

// transverse regions
   if((deltaPhi<-fUeDeltaPhiMinCut) || (deltaPhi >2*fUeDeltaPhiMaxCut))region = -1; //left
   if((deltaPhi > fUeDeltaPhiMinCut) && (deltaPhi <fUeDeltaPhiMaxCut)) region = 1;   //right

   if(deltaPhi > -fUeDeltaPhiMinCut && deltaPhi < fUeDeltaPhiMinCut) region = 2;    //forward
   if(deltaPhi > fUeDeltaPhiMaxCut && deltaPhi < 2*fUeDeltaPhiMaxCut) region = -2;  //backward

// skip leading particle
   if(leading == part) continue;
   if(part->Pt() >= leading->Pt()) continue;
   if(!region)continue;

   if(region == 1) transverse1->Add(part);
   if(region == -1) transverse2->Add(part);
   if(region == 2) toward->Add(part);
   if(region == -2) away->Add(part);
  }//end loop on tracks

  return regionParticles;
 }

//__________________________________________________
TObjArray* AliAnalysisTaskUeSpectraRT::GetMinMaxRegion(TList *transv1, TList *transv2)
{
// Returns two lists of particles, one for MIN and one for MAX region
  Double_t sumpT1 = 0.;
  Double_t sumpT2 = 0.;

  Int_t particles1 = transv1->GetEntries();
  Int_t particles2 = transv2->GetEntries();

// Loop on transverse region 1
  for(Int_t i=0; i<particles1; i++){
   AliVParticle *part = (AliVParticle*)transv1->At(i);
   sumpT1 +=  part->Pt();
   }

// Loop on transverse region 2
  for(Int_t i=0; i<particles2; i++){
   AliVParticle *part = (AliVParticle*)transv2->At(i);
   sumpT2 +=  part->Pt();
   }

  TObjArray *regionParticles = new TObjArray;
  if(sumpT2 >= sumpT1){
   regionParticles->AddLast(transv1); // MIN
   regionParticles->AddLast(transv2); // MAX
   }
  else{
   regionParticles->AddLast(transv2); // MIN
   regionParticles->AddLast(transv1); // MAX
   }

  return regionParticles;
}

//__________________________________________________
TObjArray* AliAnalysisTaskUeSpectraRT::GetRegionAwTow(TList *region)
{
// Returns list of particles
  Double_t sumpT = 0.;

  Int_t particles = region->GetEntries();

// Loop on region
  for(Int_t i=0; i<particles; i++){
   AliVParticle *part = (AliVParticle*)region->At(i);
   sumpT +=  part->Pt();
   }

  TObjArray *regionParticles = new TObjArray;
  regionParticles->AddLast(region);

  return regionParticles;
}

//__________________________________________
void AliAnalysisTaskUeSpectraRT::FillRTResponseMatrix(const AliVParticle* leadingMC, const AliVParticle* leading, TList* listMaxMC, TList* listMax, TList* listMinMC, TList* listMin){
  Double_t  LeadingPt    = leading->Pt();

  if(LeadingPt > fLeadMin && LeadingPt < 300.){
   Double_t nTracksMC = 0.;
   Double_t nTracks   = 0.;
   Double_t TrackRTMC = 0.;
   Double_t TrackRT   = 0.;

   nTracksMC = listMaxMC->GetEntries() + listMinMC->GetEntries();
   nTracks   = listMax->GetEntries() + listMin->GetEntries();

   TrackRTMC = nTracksMC/fAveGenMultiInTrans;
   TrackRT   = nTracks/fAveRecMultiInTrans;

   fhRTResponse->Fill(TrackRTMC, TrackRT);
  }
}
