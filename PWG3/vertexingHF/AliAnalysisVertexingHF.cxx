/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//----------------------------------------------------------------------------
//    Implementation of the heavy-flavour vertexing analysis class
// Candidates are stored in the AOD as objects deriving from AliAODRecoDecay.
// To be used as a task of AliAnalysisManager by means of the interface
// class AliAnalysisTaskSEVertexingHF. 
// An example of usage in the macro AliAnalysisTaskSEVertexingHFTest.C.
//
//  Contact: andrea.dainese@pd.infn.it
//  Contributors: E.Bruna, G.E.Bruno, A.Dainese, C.Di Gliglio,
//                F.Prino, R.Romita, X.M.Zhang
//----------------------------------------------------------------------------
#include <Riostream.h>
#include <TFile.h>
#include <TDatabasePDG.h>
#include <TString.h>
#include <TList.h>
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVTrack.h"
#include "AliVertexerTracks.h"
#include "AliKFVertex.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"
#include "AliNeutralTrackParam.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoDecayHF4Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsJpsitoee.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCutsLctoV0.h"
#include "AliRDHFCutsD0toKpipipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisVertexingHF.h"
#include "AliMixedEvent.h"
#include "AliESDv0.h"
#include "AliAODv0.h"

ClassImp(AliAnalysisVertexingHF)

//----------------------------------------------------------------------------
AliAnalysisVertexingHF::AliAnalysisVertexingHF():
fInputAOD(kFALSE),
fAODMapSize(0),
fAODMap(0),
fBzkG(0.),
fSecVtxWithKF(kFALSE),
fRecoPrimVtxSkippingTrks(kFALSE),
fRmTrksFromPrimVtx(kFALSE),
fV1(0x0),
fD0toKpi(kTRUE),
fJPSItoEle(kTRUE),
f3Prong(kTRUE),
f4Prong(kTRUE),
fDstar(kTRUE),
fCascades(kTRUE),
fLikeSign(kFALSE),
fLikeSign3prong(kFALSE),
fMixEvent(kFALSE),
fTrackFilter(0x0),
fTrackFilterSoftPi(0x0),
fCutsD0toKpi(0x0),
fCutsJpsitoee(0x0),
fCutsDplustoKpipi(0x0),
fCutsDstoKKpi(0x0),
fCutsLctopKpi(0x0),
fCutsLctoV0(0x0),
fCutsD0toKpipipi(0x0),
fCutsDStartoKpipi(0x0),
fListOfCuts(0x0),
fFindVertexForDstar(kTRUE),
fFindVertexForCascades(kTRUE),
fMassCutBeforeVertexing(kFALSE),
fMassCalc2(0),
fMassCalc3(0),
fMassCalc4(0),
fnTrksTotal(0),
fnSeleTrksTotal(0)
{
  // Default constructor

  Double_t d02[2]={0.,0.};
  Double_t d03[3]={0.,0.,0.};
  Double_t d04[4]={0.,0.,0.,0.};
  fMassCalc2 = new AliAODRecoDecay(0x0,2,0,d02);
  fMassCalc3 = new AliAODRecoDecay(0x0,3,1,d03);
  fMassCalc4 = new AliAODRecoDecay(0x0,4,0,d04);

}
//--------------------------------------------------------------------------
AliAnalysisVertexingHF::AliAnalysisVertexingHF(const AliAnalysisVertexingHF &source) : 
TNamed(source),
fInputAOD(source.fInputAOD),
fAODMapSize(source.fAODMapSize),
fAODMap(source.fAODMap),
fBzkG(source.fBzkG),
fSecVtxWithKF(source.fSecVtxWithKF),
fRecoPrimVtxSkippingTrks(source.fRecoPrimVtxSkippingTrks),
fRmTrksFromPrimVtx(source.fRmTrksFromPrimVtx),
fV1(source.fV1),
fD0toKpi(source.fD0toKpi),
fJPSItoEle(source.fJPSItoEle),
f3Prong(source.f3Prong),
f4Prong(source.f4Prong),
fDstar(source.fDstar),
fCascades(source.fCascades),
fLikeSign(source.fLikeSign),
fLikeSign3prong(source.fLikeSign3prong),
fMixEvent(source.fMixEvent),
fTrackFilter(source.fTrackFilter),
fTrackFilterSoftPi(source.fTrackFilterSoftPi),
fCutsD0toKpi(source.fCutsD0toKpi),
fCutsJpsitoee(source.fCutsJpsitoee),
fCutsDplustoKpipi(source.fCutsDplustoKpipi),
fCutsDstoKKpi(source.fCutsDstoKKpi),
fCutsLctopKpi(source.fCutsLctopKpi),
fCutsLctoV0(source.fCutsLctoV0),
fCutsD0toKpipipi(source.fCutsD0toKpipipi),
fCutsDStartoKpipi(source.fCutsDStartoKpipi),
fListOfCuts(source.fListOfCuts),
fFindVertexForDstar(source.fFindVertexForDstar),
fFindVertexForCascades(source.fFindVertexForCascades),
fMassCutBeforeVertexing(source.fMassCutBeforeVertexing),
fMassCalc2(source.fMassCalc2),
fMassCalc3(source.fMassCalc3),
fMassCalc4(source.fMassCalc4),
fnTrksTotal(0),
fnSeleTrksTotal(0)
{
  //
  // Copy constructor
  //
}
//--------------------------------------------------------------------------
AliAnalysisVertexingHF &AliAnalysisVertexingHF::operator=(const AliAnalysisVertexingHF &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;
  fInputAOD = source.fInputAOD;
  fAODMapSize = source.fAODMapSize;
  fBzkG = source.fBzkG;
  fSecVtxWithKF = source.fSecVtxWithKF;
  fRecoPrimVtxSkippingTrks = source.fRecoPrimVtxSkippingTrks;
  fRmTrksFromPrimVtx = source.fRmTrksFromPrimVtx;
  fV1 = source.fV1;
  fD0toKpi = source.fD0toKpi;
  fJPSItoEle = source.fJPSItoEle;
  f3Prong = source.f3Prong;
  f4Prong = source.f4Prong;
  fDstar = source.fDstar;
  fCascades = source.fCascades;
  fLikeSign = source.fLikeSign;
  fLikeSign3prong = source.fLikeSign3prong;
  fMixEvent = source.fMixEvent;
  fTrackFilter = source.fTrackFilter;
  fTrackFilterSoftPi = source.fTrackFilterSoftPi;
  fCutsD0toKpi = source.fCutsD0toKpi;
  fCutsJpsitoee = source.fCutsJpsitoee;
  fCutsDplustoKpipi = source.fCutsDplustoKpipi;
  fCutsDstoKKpi = source.fCutsDstoKKpi;
  fCutsLctopKpi = source.fCutsLctopKpi;
  fCutsLctoV0 = source.fCutsLctoV0;
  fCutsD0toKpipipi = source.fCutsD0toKpipipi;
  fCutsDStartoKpipi = source.fCutsDStartoKpipi;
  fListOfCuts = source.fListOfCuts;
  fFindVertexForDstar = source.fFindVertexForDstar;
  fFindVertexForCascades = source.fFindVertexForCascades;
  fMassCutBeforeVertexing = source.fMassCutBeforeVertexing;
  fMassCalc2 = source.fMassCalc2;
  fMassCalc3 = source.fMassCalc3;
  fMassCalc4 = source.fMassCalc4;

  return *this;
}
//----------------------------------------------------------------------------
AliAnalysisVertexingHF::~AliAnalysisVertexingHF() {
  // Destructor
  if(fV1) { delete fV1; fV1=0; }
  if(fTrackFilter) { delete fTrackFilter; fTrackFilter=0; }
  if(fTrackFilterSoftPi) { delete fTrackFilterSoftPi; fTrackFilterSoftPi=0; }
  if(fCutsD0toKpi) { delete fCutsD0toKpi; fCutsD0toKpi=0; }
  if(fCutsJpsitoee) { delete fCutsJpsitoee; fCutsJpsitoee=0; }
  if(fCutsDplustoKpipi) { delete fCutsDplustoKpipi; fCutsDplustoKpipi=0; }
  if(fCutsDstoKKpi) { delete fCutsDstoKKpi; fCutsDstoKKpi=0; }
  if(fCutsLctopKpi) { delete fCutsLctopKpi; fCutsLctopKpi=0; }
  if(fCutsLctoV0) { delete fCutsLctoV0; fCutsLctoV0=0; }
  if(fCutsD0toKpipipi) { delete fCutsD0toKpipipi; fCutsD0toKpipipi=0; }
  if(fCutsDStartoKpipi) { delete fCutsDStartoKpipi; fCutsDStartoKpipi=0; }
  if(fAODMap) { delete [] fAODMap; fAODMap=0; }
  if(fMassCalc2) { delete fMassCalc2; fMassCalc2=0; }
  if(fMassCalc3) { delete fMassCalc3; fMassCalc3=0; }
  if(fMassCalc4) { delete fMassCalc4; fMassCalc4=0; }
}
//----------------------------------------------------------------------------
TList *AliAnalysisVertexingHF::FillListOfCuts() {
  // Fill list of analysis cuts

  TList *list = new TList();
  list->SetOwner();
  list->SetName("ListOfCuts");
  
  if(fCutsD0toKpi) {
    AliRDHFCutsD0toKpi *cutsD0toKpi = new AliRDHFCutsD0toKpi(*fCutsD0toKpi);
    list->Add(cutsD0toKpi);
  }
  if(fCutsJpsitoee) {
    AliRDHFCutsJpsitoee *cutsJpsitoee = new AliRDHFCutsJpsitoee(*fCutsJpsitoee);
    list->Add(cutsJpsitoee);
  }
  if(fCutsDplustoKpipi) {
    AliRDHFCutsDplustoKpipi *cutsDplustoKpipi = new AliRDHFCutsDplustoKpipi(*fCutsDplustoKpipi);
    list->Add(cutsDplustoKpipi);
  }
  if(fCutsDstoKKpi) {
    AliRDHFCutsDstoKKpi *cutsDstoKKpi = new AliRDHFCutsDstoKKpi(*fCutsDstoKKpi);
    list->Add(cutsDstoKKpi);
  }
  if(fCutsLctopKpi) {
    AliRDHFCutsLctopKpi *cutsLctopKpi = new AliRDHFCutsLctopKpi(*fCutsLctopKpi);
    list->Add(cutsLctopKpi);
  }
  if(fCutsLctoV0){
    AliRDHFCutsLctoV0 *cutsLctoV0 = new AliRDHFCutsLctoV0(*fCutsLctoV0);
    list->Add(cutsLctoV0);
  }
  if(fCutsD0toKpipipi) {
    AliRDHFCutsD0toKpipipi *cutsD0toKpipipi = new AliRDHFCutsD0toKpipipi(*fCutsD0toKpipipi);
    list->Add(cutsD0toKpipipi);
  }
  if(fCutsDStartoKpipi) {
    AliRDHFCutsDStartoKpipi *cutsDStartoKpipi = new AliRDHFCutsDStartoKpipi(*fCutsDStartoKpipi);
    list->Add(cutsDStartoKpipi);
  }
  
  // keep a pointer to the list
  fListOfCuts = list;

  return list;
}
//----------------------------------------------------------------------------
void AliAnalysisVertexingHF::FindCandidates(AliVEvent *event,
					    TClonesArray *aodVerticesHFTClArr,
					    TClonesArray *aodD0toKpiTClArr,
					    TClonesArray *aodJPSItoEleTClArr,
					    TClonesArray *aodCharm3ProngTClArr,
					    TClonesArray *aodCharm4ProngTClArr,
					    TClonesArray *aodDstarTClArr,
					    TClonesArray *aodCascadesTClArr,
					    TClonesArray *aodLikeSign2ProngTClArr,
					    TClonesArray *aodLikeSign3ProngTClArr)
{
  // Find heavy-flavour vertex candidates
  // Input:  ESD or AOD
  // Output: AOD (additional branches added)

  if(!fMixEvent){
    TString evtype = event->IsA()->GetName();
    fInputAOD = ((evtype=="AliAODEvent") ? kTRUE : kFALSE);
  } // if we do mixing AliVEvent is a AliMixedEvent

  if(fInputAOD) {
    AliDebug(2,"Creating HF candidates from AOD");
  } else {
    AliDebug(2,"Creating HF candidates from ESD");
  }

  if(!aodVerticesHFTClArr) {
    printf("ERROR: no aodVerticesHFTClArr");
    return;
  }
  if((fD0toKpi || fDstar) && !aodD0toKpiTClArr) {
    printf("ERROR: no aodD0toKpiTClArr");
    return;
  }
  if(fJPSItoEle && !aodJPSItoEleTClArr) {
    printf("ERROR: no aodJPSItoEleTClArr");
    return;
  }
  if(f3Prong && !aodCharm3ProngTClArr) {
    printf("ERROR: no aodCharm3ProngTClArr");
    return;
  }
  if(f4Prong && !aodCharm4ProngTClArr) {
    printf("ERROR: no aodCharm4ProngTClArr");
    return;
  }
  if(fDstar && !aodDstarTClArr) {
    printf("ERROR: no aodDstarTClArr");
    return;
  }
  if(fCascades && !aodCascadesTClArr){
    printf("ERROR: no aodCascadesTClArr ");
    return;
  }
  if(fLikeSign && !aodLikeSign2ProngTClArr) {
    printf("ERROR: no aodLikeSign2ProngTClArr");
    return;
  }
  if(fLikeSign3prong && f3Prong && !aodLikeSign3ProngTClArr) {
    printf("ERROR: no aodLikeSign3ProngTClArr");
    return;
  }

  // delete candidates from previous event and create references
  Int_t iVerticesHF=0,iD0toKpi=0,iJPSItoEle=0,i3Prong=0,i4Prong=0,iDstar=0,iCascades=0,iLikeSign2Prong=0,iLikeSign3Prong=0;
  aodVerticesHFTClArr->Delete();
  iVerticesHF = aodVerticesHFTClArr->GetEntriesFast();
  TClonesArray &verticesHFRef = *aodVerticesHFTClArr;
  if(fD0toKpi || fDstar)   {
    aodD0toKpiTClArr->Delete();
    iD0toKpi = aodD0toKpiTClArr->GetEntriesFast();
  }
  if(fJPSItoEle) {
    aodJPSItoEleTClArr->Delete();
    iJPSItoEle = aodJPSItoEleTClArr->GetEntriesFast();
  }
  if(f3Prong) {   
    aodCharm3ProngTClArr->Delete();
    i3Prong = aodCharm3ProngTClArr->GetEntriesFast();
  }
  if(f4Prong) {
    aodCharm4ProngTClArr->Delete();
    i4Prong = aodCharm4ProngTClArr->GetEntriesFast();
  }
  if(fDstar) {
    aodDstarTClArr->Delete();
    iDstar = aodDstarTClArr->GetEntriesFast();
  }
  if(fCascades) {
    aodCascadesTClArr->Delete();
    iCascades = aodCascadesTClArr->GetEntriesFast();
  }
  if(fLikeSign) {                                
    aodLikeSign2ProngTClArr->Delete();                     
    iLikeSign2Prong = aodLikeSign2ProngTClArr->GetEntriesFast(); 
  }  
  if(fLikeSign3prong && f3Prong) {                                
    aodLikeSign3ProngTClArr->Delete();                     
    iLikeSign3Prong = aodLikeSign3ProngTClArr->GetEntriesFast(); 
  }  

  TClonesArray &aodD0toKpiRef        = *aodD0toKpiTClArr;
  TClonesArray &aodJPSItoEleRef      = *aodJPSItoEleTClArr;
  TClonesArray &aodCharm3ProngRef    = *aodCharm3ProngTClArr;
  TClonesArray &aodCharm4ProngRef    = *aodCharm4ProngTClArr;
  TClonesArray &aodDstarRef          = *aodDstarTClArr;
  TClonesArray &aodCascadesRef       = *aodCascadesTClArr;
  TClonesArray &aodLikeSign2ProngRef = *aodLikeSign2ProngTClArr;
  TClonesArray &aodLikeSign3ProngRef = *aodLikeSign3ProngTClArr;


  AliAODRecoDecayHF2Prong *io2Prong  = 0;
  AliAODRecoDecayHF3Prong *io3Prong  = 0;
  AliAODRecoDecayHF4Prong *io4Prong  = 0;
  AliAODRecoCascadeHF     *ioCascade = 0;

  Int_t    iTrkP1,iTrkP2,iTrkN1,iTrkN2,iTrkSoftPi,trkEntries,iv0,nv0;
  Double_t xdummy,ydummy,dcap1n1,dcap1n2,dcap2n1,dcap1p2,dcan1n2,dcap2n2,dcaV0,dcaCasc;
  Bool_t   okD0=kFALSE,okJPSI=kFALSE,ok3Prong=kFALSE,ok4Prong=kFALSE;
  Bool_t   okDstar=kFALSE,okD0fromDstar=kFALSE;
  Bool_t   okCascades=kFALSE;
  AliESDtrack *postrack1 = 0;
  AliESDtrack *postrack2 = 0;
  AliESDtrack *negtrack1 = 0;
  AliESDtrack *negtrack2 = 0;
  AliESDtrack *trackPi   = 0;
  //   AliESDtrack *posV0track = 0;
  //   AliESDtrack *negV0track = 0;
  Float_t dcaMax = fCutsD0toKpi->GetDCACut();
  if(fCutsJpsitoee) dcaMax=TMath::Max(dcaMax,fCutsJpsitoee->GetDCACut());
  if(fCutsDplustoKpipi) dcaMax=TMath::Max(dcaMax,fCutsDplustoKpipi->GetDCACut());
  if(fCutsDstoKKpi) dcaMax=TMath::Max(dcaMax,fCutsDstoKKpi->GetDCACut());
  if(fCutsLctopKpi) dcaMax=TMath::Max(dcaMax,fCutsLctopKpi->GetDCACut());
  if(fCutsD0toKpipipi) dcaMax=TMath::Max(dcaMax,fCutsD0toKpipipi->GetDCACut());
  if(fCutsDStartoKpipi) dcaMax=TMath::Max(dcaMax,fCutsDStartoKpipi->GetDCACut());
  
  AliDebug(2,Form(" dca cut set to %f cm",dcaMax));


  // get Bz
  fBzkG = (Double_t)event->GetMagneticField(); 

  trkEntries = (Int_t)event->GetNumberOfTracks();
  AliDebug(1,Form(" Number of tracks: %d",trkEntries));
  fnTrksTotal += trkEntries;

  nv0 = (Int_t)event->GetNumberOfV0s();
  AliDebug(1,Form(" Number of V0s: %d",nv0));

  if( trkEntries<2 && (trkEntries<1 || nv0<1) ) {
    AliDebug(1,Form(" Not enough tracks: %d",trkEntries));
    return;
  }

  // event selection
  if(!fCutsD0toKpi->IsEventSelected(event)) return;

  // call function that applies sigle-track selection,
  // for displaced tracks and soft pions (both charges) for D*,
  // and retrieves primary vertex
  TObjArray seleTrksArray(trkEntries);
  TObjArray tracksAtVertex(trkEntries);
  UChar_t  *seleFlags = new UChar_t[trkEntries]; // bit 0: displaced, bit 1: softpi
  Int_t     nSeleTrks=0;
  Int_t *evtNumber    = new Int_t[trkEntries];
  SelectTracksAndCopyVertex(event,trkEntries,seleTrksArray,tracksAtVertex,nSeleTrks,seleFlags,evtNumber);

  AliDebug(1,Form(" Selected tracks: %d",nSeleTrks));
  fnSeleTrksTotal += nSeleTrks;
    

  TObjArray *twoTrackArray1    = new TObjArray(2);
  TObjArray *twoTrackArray2    = new TObjArray(2);
  TObjArray *twoTrackArrayV0   = new TObjArray(2);
  TObjArray *twoTrackArrayCasc = new TObjArray(2);
  TObjArray *threeTrackArray   = new TObjArray(3);
  TObjArray *fourTrackArray    = new TObjArray(4);

  Double_t dispersion;
  Bool_t isLikeSign2Prong=kFALSE,isLikeSign3Prong=kFALSE;

  AliAODRecoDecayHF   *rd = 0;
  AliAODRecoCascadeHF *rc = 0;
  AliAODv0            *v0 = 0;
  AliESDv0         *esdV0 = 0;

  Bool_t massCutOK=kTRUE;

  // LOOP ON  POSITIVE  TRACKS
  for(iTrkP1=0; iTrkP1<nSeleTrks; iTrkP1++) {

    //if(iTrkP1%1==0) AliDebug(1,Form("  1st loop on pos: track number %d of %d",iTrkP1,nSeleTrks));  
    //if(iTrkP1%1==0) printf("  1st loop on pos: track number %d of %d\n",iTrkP1,nSeleTrks);  

    // get track from tracks array
    postrack1 = (AliESDtrack*)seleTrksArray.UncheckedAt(iTrkP1);

    if(!TESTBIT(seleFlags[iTrkP1],kBitDispl)) continue;

    // Make cascades with V0+track
    // 
    if(fCascades) {
      // loop on V0's
      for(iv0=0; iv0<nv0; iv0++){

	//AliDebug(1,Form("   loop on v0s for track number %d and v0 number %d",iTrkP1,iv0));	

	// Get the V0 
	if(fInputAOD) {
	  v0 = ((AliAODEvent*)event)->GetV0(iv0);
	} else {
	  esdV0 = ((AliESDEvent*)event)->GetV0(iv0);
	}
	if ( (!v0 || !v0->IsA()->InheritsFrom("AliAODv0") ) && 
	     (!esdV0 || !esdV0->IsA()->InheritsFrom("AliESDv0") ) ) continue;
	

	// Get the tracks that form the V0
	//  ( parameters at primary vertex )
	//   and define an AliExternalTrackParam out of them
	AliExternalTrackParam * posV0track;
	AliExternalTrackParam * negV0track;

	if(fInputAOD){
          AliAODTrack *posVV0track = (AliAODTrack*)(v0->GetDaughter(0));
	  AliAODTrack *negVV0track = (AliAODTrack*)(v0->GetDaughter(1));
	  if( !posVV0track || !negVV0track ) continue;
	  //
	  // Apply some basic V0 daughter criteria
	  //
	  // bachelor must not be a v0-track                                                                  
	  if (posVV0track->GetID() == postrack1->GetID() ||
	      negVV0track->GetID() == postrack1->GetID()) continue;
	  // reject like-sign v0                                                                              
	  if ( posVV0track->Charge() == negVV0track->Charge() ) continue;
	  // avoid ghost TPC tracks                                                                           
	  if(!(posVV0track->GetStatus() & AliESDtrack::kTPCrefit) ||
	     !(negVV0track->GetStatus() & AliESDtrack::kTPCrefit)) continue;
	  // Get AliExternalTrackParam out of the AliAODTracks
	  Double_t xyz[3], pxpypz[3], cv[21]; Short_t sign;
	  posVV0track->PxPyPz(pxpypz); 	                  posVV0track->XvYvZv(xyz);
	  posVV0track->GetCovarianceXYZPxPyPz(cv);	  sign=posVV0track->Charge();
	  posV0track = new AliExternalTrackParam(xyz,pxpypz,cv,sign);
	  negVV0track->PxPyPz(pxpypz); 	                  negVV0track->XvYvZv(xyz);
	  negVV0track->GetCovarianceXYZPxPyPz(cv);	  sign=negVV0track->Charge();
	  negV0track = new AliExternalTrackParam(xyz,pxpypz,cv,sign);
	}  else {
	  AliESDtrack *posVV0track = (AliESDtrack*)(event->GetTrack( esdV0->GetPindex() ));
          AliESDtrack *negVV0track = (AliESDtrack*)(event->GetTrack( esdV0->GetNindex() ));
	  if( !posVV0track || !negVV0track ) continue;
	  //
	  // Apply some basic V0 daughter criteria
	  //
	  // bachelor must not be a v0-track                                                                  
	  if (posVV0track->GetID() == postrack1->GetID() ||
	      negVV0track->GetID() == postrack1->GetID()) continue;
	  // reject like-sign v0                                                                              
	  if ( posVV0track->Charge() == negVV0track->Charge() ) continue;
	  // avoid ghost TPC tracks                                                                           
	  if(!(posVV0track->GetStatus() & AliESDtrack::kTPCrefit) ||
	     !(negVV0track->GetStatus() & AliESDtrack::kTPCrefit)) continue;
	  //  reject kinks (only necessary on AliESDtracks)
	  if (posVV0track->GetKinkIndex(0)>0  || negVV0track->GetKinkIndex(0)>0) continue;
	  // Get AliExternalTrackParam out of the AliESDtracks	
	  posV0track = new AliExternalTrackParam(*posVV0track);
	  negV0track = new AliExternalTrackParam(*negVV0track);

	  // Define the AODv0 from ESDv0 if reading ESDs
	  v0 = TransformESDv0toAODv0(esdV0,twoTrackArrayV0);
	}
	if( !posV0track || !negV0track ){
	  AliDebug(1,Form(" Couldn't get the V0 daughters"));
	  continue;
	}

	// fill in the v0 two-external-track-param array
 	twoTrackArrayV0->AddAt(posV0track,0);
 	twoTrackArrayV0->AddAt(negV0track,1);

 	// Get the V0 dca
	dcaV0 = v0->DcaV0Daughters();

	// Define the V0 (neutral) track
	AliNeutralTrackParam *trackV0;
	if(fInputAOD) {
	  const AliVTrack *trackVV0 = dynamic_cast<const AliVTrack*>(v0);
	  if(trackVV0) trackV0 = new AliNeutralTrackParam(trackVV0);
	} else {  
	  Double_t xyz[3], pxpypz[3];
	  esdV0->XvYvZv(xyz);
	  esdV0->PxPyPz(pxpypz);
	  Double_t cv[21]; for(int i=0; i<21; i++) cv[i]=0;
	  trackV0 = new AliNeutralTrackParam(xyz,pxpypz,cv,0);
	}
	// Fill in the object array to create the cascade
	twoTrackArrayCasc->AddAt(postrack1,0);
	twoTrackArrayCasc->AddAt(trackV0,1);
	// Compute the cascade vertex
	AliAODVertex *vertexCasc = 0;
	if(fFindVertexForCascades) {  
	  // DCA between the two tracks
	  dcaCasc = postrack1->GetDCA(trackV0,fBzkG,xdummy,ydummy);
	  // Vertexing+	  
	  vertexCasc = ReconstructSecondaryVertex(twoTrackArrayCasc,dispersion,kFALSE);
	} else {
	  // assume Cascade decays at the primary vertex
	  Double_t pos[3],cov[6],chi2perNDF;
	  fV1->GetXYZ(pos);
	  fV1->GetCovMatrix(cov);
	  chi2perNDF = fV1->GetChi2toNDF();
	  vertexCasc = new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,2);
	  dcaCasc = 0.;
	}
	if(!vertexCasc) {
	  delete posV0track; posV0track=NULL;
	  delete negV0track; negV0track=NULL;
	  delete trackV0; trackV0=NULL;
	  if(!fInputAOD) {delete v0; v0=NULL;}
	  twoTrackArrayV0->Clear();
	  twoTrackArrayCasc->Clear();
	  continue; 
	}

	// Create and store the Cascade if passed the cuts
	ioCascade = MakeCascade(twoTrackArrayCasc,event,vertexCasc,v0,dcaCasc,okCascades);
	if(okCascades && ioCascade) {
	  //AliDebug(1,Form("Storing a cascade object... "));
	  // add the vertex and the cascade to the AOD
	  AliAODVertex *vCasc = new(verticesHFRef[iVerticesHF++])AliAODVertex(*vertexCasc);
	  rc = new(aodCascadesRef[iCascades++])AliAODRecoCascadeHF(*ioCascade);
	  rc->SetSecondaryVtx(vCasc);
	  vCasc->SetParent(rc);
	  rc->SetPrimaryVtxRef((AliAODVertex*)event->GetPrimaryVertex());
	  if(!fInputAOD) vCasc->AddDaughter(v0); // just to fill ref #0 ??
	  AddRefs(vCasc,rc,event,twoTrackArrayCasc); // add the track (proton)
	  vCasc->AddDaughter(v0); // fill the 2prong V0 
	}

	// Clean up 
	delete posV0track; posV0track=NULL;
	delete negV0track; negV0track=NULL;
	delete trackV0; trackV0=NULL;
	twoTrackArrayV0->Clear();
	twoTrackArrayCasc->Clear();
	if(ioCascade) { delete ioCascade; ioCascade=NULL; }
	if(vertexCasc) { delete vertexCasc; vertexCasc=NULL; }
	if(!fInputAOD) {delete v0; v0=NULL;}

      } // end loop on V0's
    } 
  
    // If there is less than 2 particles continue
    if(trkEntries<2) {
      AliDebug(1,Form(" Not enough tracks: %d",trkEntries));
      continue;
    }

    if(postrack1->Charge()<0 && !fLikeSign) continue;

    // LOOP ON  NEGATIVE  TRACKS
    for(iTrkN1=0; iTrkN1<nSeleTrks; iTrkN1++) {

      //if(iTrkN1%1==0) AliDebug(1,Form("    1st loop on neg: track number %d of %d",iTrkN1,nSeleTrks));  
      //if(iTrkN1%1==0) printf("    1st loop on neg: track number %d of %d\n",iTrkN1,nSeleTrks);  

      if(iTrkN1==iTrkP1) continue;

      // get track from tracks array
      negtrack1 = (AliESDtrack*)seleTrksArray.UncheckedAt(iTrkN1);

      if(negtrack1->Charge()>0 && !fLikeSign) continue;

      if(!TESTBIT(seleFlags[iTrkN1],kBitDispl)) continue;

      if(fMixEvent) {
	if(evtNumber[iTrkP1]==evtNumber[iTrkN1]) continue;
      }

      if(postrack1->Charge()==negtrack1->Charge()) { // like-sign 
	isLikeSign2Prong=kTRUE;
	if(!fLikeSign)    continue;
	if(iTrkN1<iTrkP1) continue; // this is needed to avoid double-counting of like-sign
      } else { // unlike-sign
	isLikeSign2Prong=kFALSE;
	if(postrack1->Charge()<0 || negtrack1->Charge()>0) continue;  // this is needed to avoid double-counting of unlike-sign
	if(fMixEvent) {
	  if(evtNumber[iTrkP1]==evtNumber[iTrkN1]) continue;
	}
       
      }

      // back to primary vertex
      //      postrack1->PropagateToDCA(fV1,fBzkG,kVeryBig);
      //      negtrack1->PropagateToDCA(fV1,fBzkG,kVeryBig);
      SetParametersAtVertex(postrack1,(AliExternalTrackParam*)tracksAtVertex.UncheckedAt(iTrkP1));
      SetParametersAtVertex(negtrack1,(AliExternalTrackParam*)tracksAtVertex.UncheckedAt(iTrkN1));

      // DCA between the two tracks
      dcap1n1 = postrack1->GetDCA(negtrack1,fBzkG,xdummy,ydummy);
      if(dcap1n1>dcaMax) { negtrack1=0; continue; }

      // Vertexing
      twoTrackArray1->AddAt(postrack1,0);
      twoTrackArray1->AddAt(negtrack1,1);
      AliAODVertex *vertexp1n1 = ReconstructSecondaryVertex(twoTrackArray1,dispersion);
      if(!vertexp1n1) {
	twoTrackArray1->Clear();
	negtrack1=0; 
	continue; 
      }

      // 2 prong candidate
      if(fD0toKpi || fJPSItoEle || fDstar || fLikeSign) { 
      
	io2Prong = Make2Prong(twoTrackArray1,event,vertexp1n1,dcap1n1,okD0,okJPSI,okD0fromDstar);
      
	if((fD0toKpi && okD0) || (fJPSItoEle && okJPSI) || (isLikeSign2Prong && (okD0 || okJPSI))) {
	  // add the vertex and the decay to the AOD
	  AliAODVertex *v2Prong = new(verticesHFRef[iVerticesHF++])AliAODVertex(*vertexp1n1);
	  if(!isLikeSign2Prong) {
	    if(okD0) {  
	      rd = new(aodD0toKpiRef[iD0toKpi++])AliAODRecoDecayHF2Prong(*io2Prong);
	      rd->SetSecondaryVtx(v2Prong);
	      v2Prong->SetParent(rd);
	      AddRefs(v2Prong,rd,event,twoTrackArray1);
	    }
	    if(okJPSI) {
	      rd = new(aodJPSItoEleRef[iJPSItoEle++])AliAODRecoDecayHF2Prong(*io2Prong);
	      rd->SetSecondaryVtx(v2Prong);
	      if(!okD0) v2Prong->SetParent(rd); // it cannot have two mothers ...
	      AddRefs(v2Prong,rd,event,twoTrackArray1);
	    }
	  } else { // isLikeSign2Prong
	    rd = new(aodLikeSign2ProngRef[iLikeSign2Prong++])AliAODRecoDecayHF2Prong(*io2Prong);
	    rd->SetSecondaryVtx(v2Prong);
	    v2Prong->SetParent(rd);
	    AddRefs(v2Prong,rd,event,twoTrackArray1);
	  }
	  // Set selection bit for PID
	  if(okD0) SetSelectionBitForPID(fCutsD0toKpi,rd,AliRDHFCuts::kD0toKpiPID);
	}
	// D* candidates
	if(fDstar && okD0fromDstar && !isLikeSign2Prong) {
	  // write references in io2Prong
	  if(fInputAOD) {
	    AddDaughterRefs(vertexp1n1,event,twoTrackArray1);
	  } else {
	    vertexp1n1->AddDaughter(postrack1);
	    vertexp1n1->AddDaughter(negtrack1);
	  }
	  io2Prong->SetSecondaryVtx(vertexp1n1);
	  //printf("--->  %d %d %d %d %d\n",vertexp1n1->GetNDaughters(),iTrkP1,iTrkN1,postrack1->Charge(),negtrack1->Charge());
	  // create a track from the D0
	  AliNeutralTrackParam *trackD0 = new AliNeutralTrackParam(io2Prong);

	  // LOOP ON TRACKS THAT PASSED THE SOFT PION CUTS
	  for(iTrkSoftPi=0; iTrkSoftPi<nSeleTrks; iTrkSoftPi++) {

	    if(iTrkSoftPi==iTrkP1 || iTrkSoftPi==iTrkN1) continue;

	    if(!TESTBIT(seleFlags[iTrkSoftPi],kBitSoftPi)) continue;

	    if(fMixEvent) {
	      if(evtNumber[iTrkP1]==evtNumber[iTrkSoftPi] || 
		 evtNumber[iTrkN1]==evtNumber[iTrkSoftPi] || 
		 evtNumber[iTrkP1]==evtNumber[iTrkN1]) continue;
	    }

	    //if(iTrkSoftPi%1==0) AliDebug(1,Form("    1st loop on pi_s: track number %d of %d",iTrkSoftPi,nSeleTrks));  

	    trackD0->PropagateToDCA(fV1,fBzkG,kVeryBig);
	    if(trackD0->GetSigmaY2()<0. || trackD0->GetSigmaZ2()<0.) continue; // this is insipired by the AliITStrackV2::Invariant() checks

	    // get track from tracks array
	    trackPi = (AliESDtrack*)seleTrksArray.UncheckedAt(iTrkSoftPi);
	    //	    trackPi->PropagateToDCA(fV1,fBzkG,kVeryBig);
	    SetParametersAtVertex(trackPi,(AliExternalTrackParam*)tracksAtVertex.UncheckedAt(iTrkSoftPi));
	    twoTrackArrayCasc->AddAt(trackPi,0);
	    twoTrackArrayCasc->AddAt(trackD0,1);

	    AliAODVertex *vertexCasc = 0;

	    if(fFindVertexForDstar) {
	      // DCA between the two tracks
	      dcaCasc = trackPi->GetDCA(trackD0,fBzkG,xdummy,ydummy);
	      // Vertexing
	      vertexCasc = ReconstructSecondaryVertex(twoTrackArrayCasc,dispersion,kFALSE);
	    } else {
	      // assume Dstar decays at the primary vertex
	      Double_t pos[3],cov[6],chi2perNDF;
	      fV1->GetXYZ(pos);
	      fV1->GetCovMatrix(cov);
	      chi2perNDF = fV1->GetChi2toNDF();
	      vertexCasc = new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,2);
	      dcaCasc = 0.;
	    }
	    if(!vertexCasc) { 
	      twoTrackArrayCasc->Clear();
	      trackPi=0; 
	      continue; 
	    }

            ioCascade = MakeCascade(twoTrackArrayCasc,event,vertexCasc,io2Prong,dcaCasc,okDstar);
            if(okDstar) {
	      // add the D0 to the AOD (if not already done)
	      if(!okD0) {
		AliAODVertex *v2Prong = new(verticesHFRef[iVerticesHF++])AliAODVertex(*vertexp1n1);
	        rd = new(aodD0toKpiRef[iD0toKpi++])AliAODRecoDecayHF2Prong(*io2Prong);
		rd->SetSecondaryVtx(v2Prong);
		v2Prong->SetParent(rd);
		AddRefs(v2Prong,rd,event,twoTrackArray1);
		okD0=kTRUE; // this is done to add it only once
	      }
	      // add the vertex and the cascade to the AOD
	      AliAODVertex *vCasc = new(verticesHFRef[iVerticesHF++])AliAODVertex(*vertexCasc); 
	      rc = new(aodDstarRef[iDstar++])AliAODRecoCascadeHF(*ioCascade);
	      rc->SetSecondaryVtx(vCasc);
	      vCasc->SetParent(rc);
	      if(!fInputAOD) vCasc->AddDaughter(rd); // just to fill ref #0 
	      AddRefs(vCasc,rc,event,twoTrackArrayCasc);
	      vCasc->AddDaughter(rd); // add the D0 (in ref #1)
	      // Set selection bit for PID
	      SetSelectionBitForPID(fCutsDStartoKpipi,rc,AliRDHFCuts::kDstarPID);
	    }
	    twoTrackArrayCasc->Clear();
	    trackPi=0; 
	    if(ioCascade) {delete ioCascade; ioCascade=NULL;}
	    delete vertexCasc; vertexCasc=NULL;
	  } // end loop on soft pi tracks

	  if(trackD0) {delete trackD0; trackD0=NULL;}

	}
	if(io2Prong) {delete io2Prong; io2Prong=NULL;}
      }      

      twoTrackArray1->Clear(); 
      if( (!f3Prong && !f4Prong) || 
	  (isLikeSign2Prong && !f3Prong) ) { 
	negtrack1=0; 
	delete vertexp1n1; 
	continue; 
      }

	
      // 2nd LOOP  ON  POSITIVE  TRACKS 
      for(iTrkP2=iTrkP1+1; iTrkP2<nSeleTrks; iTrkP2++) {

	if(iTrkP2==iTrkP1 || iTrkP2==iTrkN1) continue;

	//if(iTrkP2%1==0) AliDebug(1,Form("    2nd loop on pos: track number %d of %d",iTrkP2,nSeleTrks));  

	// get track from tracks array
	postrack2 = (AliESDtrack*)seleTrksArray.UncheckedAt(iTrkP2);

	if(postrack2->Charge()<0) continue; 

	if(!TESTBIT(seleFlags[iTrkP2],kBitDispl)) continue;

	if(fMixEvent) {
	  if(evtNumber[iTrkP1]==evtNumber[iTrkP2] || 
	     evtNumber[iTrkN1]==evtNumber[iTrkP2] ||
	     evtNumber[iTrkP1]==evtNumber[iTrkN1]) continue;
	}

	if(isLikeSign2Prong) { // like-sign pair -> have to build only like-sign triplet 
	  if(!fLikeSign3prong) continue;
	  if(postrack1->Charge()>0) { // ok: like-sign triplet (+++)
	    isLikeSign3Prong=kTRUE;
	  } else { // not ok
	    continue;
	  }
	} else { // normal triplet (+-+)
	  isLikeSign3Prong=kFALSE; 
	  if(fMixEvent) {
	    if(evtNumber[iTrkP1]==evtNumber[iTrkP2] || 
	       evtNumber[iTrkN1]==evtNumber[iTrkP2] ||
	       evtNumber[iTrkP1]==evtNumber[iTrkN1]) continue;
	  }
	}

	// back to primary vertex
	//	postrack1->PropagateToDCA(fV1,fBzkG,kVeryBig);
	//	postrack2->PropagateToDCA(fV1,fBzkG,kVeryBig);
	//	negtrack1->PropagateToDCA(fV1,fBzkG,kVeryBig);
	SetParametersAtVertex(postrack1,(AliExternalTrackParam*)tracksAtVertex.UncheckedAt(iTrkP1));
	SetParametersAtVertex(negtrack1,(AliExternalTrackParam*)tracksAtVertex.UncheckedAt(iTrkN1));
	SetParametersAtVertex(postrack2,(AliExternalTrackParam*)tracksAtVertex.UncheckedAt(iTrkP2));
      
	//printf("********** %d %d %d\n",postrack1->GetID(),postrack2->GetID(),negtrack1->GetID());

	dcap2n1 = postrack2->GetDCA(negtrack1,fBzkG,xdummy,ydummy);
	if(dcap2n1>dcaMax) { postrack2=0; continue; }
	dcap1p2 = postrack2->GetDCA(postrack1,fBzkG,xdummy,ydummy);
	if(dcap1p2>dcaMax) { postrack2=0; continue; }

	// check invariant mass cuts for D+,Ds,Lc
        massCutOK=kTRUE;
	if(f3Prong) {
	  if(postrack2->Charge()>0) {
	    threeTrackArray->AddAt(postrack1,0);
	    threeTrackArray->AddAt(negtrack1,1);
	    threeTrackArray->AddAt(postrack2,2);
	  } else {
	    threeTrackArray->AddAt(negtrack1,0);
	    threeTrackArray->AddAt(postrack1,1);
	    threeTrackArray->AddAt(postrack2,2);
	  }
	  if(fMassCutBeforeVertexing)
	    massCutOK = SelectInvMassAndPt(threeTrackArray);
	}

	if(f3Prong && !massCutOK) {
	  threeTrackArray->Clear();
	  if(!f4Prong) { 
	    postrack2=0; 
	    continue;
	  } 
	} 
	
	// Vertexing
	twoTrackArray2->AddAt(postrack2,0);
	twoTrackArray2->AddAt(negtrack1,1);
	AliAODVertex *vertexp2n1 = ReconstructSecondaryVertex(twoTrackArray2,dispersion);
	if(!vertexp2n1) { 
	  twoTrackArray2->Clear();
	  postrack2=0; 
	  continue;
	}

	// 3 prong candidates
	if(f3Prong && massCutOK) { 

	  AliAODVertex* secVert3PrAOD = ReconstructSecondaryVertex(threeTrackArray,dispersion);
	  io3Prong = Make3Prong(threeTrackArray,event,secVert3PrAOD,dispersion,vertexp1n1,vertexp2n1,dcap1n1,dcap2n1,dcap1p2,ok3Prong);
	  if(ok3Prong) {
	    AliAODVertex *v3Prong = new(verticesHFRef[iVerticesHF++])AliAODVertex(*secVert3PrAOD);
	    if(!isLikeSign3Prong) {
	      rd = new(aodCharm3ProngRef[i3Prong++])AliAODRecoDecayHF3Prong(*io3Prong);
	      rd->SetSecondaryVtx(v3Prong);
	      v3Prong->SetParent(rd);
	      AddRefs(v3Prong,rd,event,threeTrackArray);
	    } else { // isLikeSign3Prong
	      if(fLikeSign3prong){
		rd = new(aodLikeSign3ProngRef[iLikeSign3Prong++])AliAODRecoDecayHF3Prong(*io3Prong);
		rd->SetSecondaryVtx(v3Prong);
		v3Prong->SetParent(rd);
		AddRefs(v3Prong,rd,event,threeTrackArray);
	      }
	    }
	    // Set selection bit for PID
	    SetSelectionBitForPID(fCutsDplustoKpipi,rd,AliRDHFCuts::kDplusPID);
	    SetSelectionBitForPID(fCutsDstoKKpi,rd,AliRDHFCuts::kDsPID);
	    SetSelectionBitForPID(fCutsLctopKpi,rd,AliRDHFCuts::kLcPID);
	  }
	  if(io3Prong) {delete io3Prong; io3Prong=NULL;} 
	  if(secVert3PrAOD) {delete secVert3PrAOD; secVert3PrAOD=NULL;} 
	}

	// 4 prong candidates
	if(f4Prong 
	   // don't make 4 prong with like-sign pairs and triplets
	   && !isLikeSign2Prong && !isLikeSign3Prong
	   // track-to-track dca cuts already now
	   && dcap1n1 < fCutsD0toKpipipi->GetDCACut()
	   && dcap2n1 < fCutsD0toKpipipi->GetDCACut()) {

	  // back to primary vertex
	  //	  postrack1->PropagateToDCA(fV1,fBzkG,kVeryBig);
	  //	  postrack2->PropagateToDCA(fV1,fBzkG,kVeryBig);
	  //	  negtrack1->PropagateToDCA(fV1,fBzkG,kVeryBig);
	  SetParametersAtVertex(postrack1,(AliExternalTrackParam*)tracksAtVertex.UncheckedAt(iTrkP1));
	  SetParametersAtVertex(negtrack1,(AliExternalTrackParam*)tracksAtVertex.UncheckedAt(iTrkN1));
	  SetParametersAtVertex(postrack2,(AliExternalTrackParam*)tracksAtVertex.UncheckedAt(iTrkP2));

	  // Vertexing for these 3 (can be taken from above?)
          threeTrackArray->AddAt(postrack1,0);
          threeTrackArray->AddAt(negtrack1,1);
	  threeTrackArray->AddAt(postrack2,2);
          AliAODVertex* vertexp1n1p2 = ReconstructSecondaryVertex(threeTrackArray,dispersion);

	  // 3rd LOOP  ON  NEGATIVE  TRACKS (for 4 prong) 
	  for(iTrkN2=iTrkN1+1; iTrkN2<nSeleTrks; iTrkN2++) {

	    if(iTrkN2==iTrkP1 || iTrkN2==iTrkP2 || iTrkN2==iTrkN1) continue;

	    //if(iTrkN2%1==0) AliDebug(1,Form("    3rd loop on neg: track number %d of %d",iTrkN2,nSeleTrks));  

	    // get track from tracks array
	    negtrack2 = (AliESDtrack*)seleTrksArray.UncheckedAt(iTrkN2);

	    if(negtrack2->Charge()>0) continue;

	    if(!TESTBIT(seleFlags[iTrkN2],kBitDispl)) continue;
	    if(fMixEvent){ 
	      if(evtNumber[iTrkP1]==evtNumber[iTrkN2] || 
		 evtNumber[iTrkN1]==evtNumber[iTrkN2] || 
		 evtNumber[iTrkP2]==evtNumber[iTrkN2] ||
		 evtNumber[iTrkP1]==evtNumber[iTrkN1] ||
		 evtNumber[iTrkP1]==evtNumber[iTrkP2] ||
		 evtNumber[iTrkN1]==evtNumber[iTrkP2]) continue;
	    }

	    // back to primary vertex
	    // postrack1->PropagateToDCA(fV1,fBzkG,kVeryBig);
	    // postrack2->PropagateToDCA(fV1,fBzkG,kVeryBig);
	    // negtrack1->PropagateToDCA(fV1,fBzkG,kVeryBig);
	    // negtrack2->PropagateToDCA(fV1,fBzkG,kVeryBig);
	    SetParametersAtVertex(postrack1,(AliExternalTrackParam*)tracksAtVertex.UncheckedAt(iTrkP1));
	    SetParametersAtVertex(negtrack1,(AliExternalTrackParam*)tracksAtVertex.UncheckedAt(iTrkN1));
	    SetParametersAtVertex(postrack2,(AliExternalTrackParam*)tracksAtVertex.UncheckedAt(iTrkP2));
	    SetParametersAtVertex(negtrack2,(AliExternalTrackParam*)tracksAtVertex.UncheckedAt(iTrkN2));

	    dcap1n2 = postrack1->GetDCA(negtrack2,fBzkG,xdummy,ydummy);
	    if(dcap1n2 > fCutsD0toKpipipi->GetDCACut()) { negtrack2=0; continue; }
            dcap2n2 = postrack2->GetDCA(negtrack2,fBzkG,xdummy,ydummy);
            if(dcap2n2 > fCutsD0toKpipipi->GetDCACut()) { negtrack2=0; continue; }


	    fourTrackArray->AddAt(postrack1,0);
	    fourTrackArray->AddAt(negtrack1,1);
	    fourTrackArray->AddAt(postrack2,2);
	    fourTrackArray->AddAt(negtrack2,3);

	    // check invariant mass cuts for D0
	    massCutOK=kTRUE;
	    if(fMassCutBeforeVertexing) 
	      massCutOK = SelectInvMassAndPt(fourTrackArray);
	    
	    if(!massCutOK) {
	      fourTrackArray->Clear(); 
	      negtrack2=0; 
	      continue; 
	    }

	    // Vertexing
	    AliAODVertex* secVert4PrAOD = ReconstructSecondaryVertex(fourTrackArray,dispersion);
	    io4Prong = Make4Prong(fourTrackArray,event,secVert4PrAOD,vertexp1n1,vertexp1n1p2,dcap1n1,dcap1n2,dcap2n1,dcap2n2,ok4Prong);
	    if(ok4Prong) {
	      AliAODVertex *v4Prong = new(verticesHFRef[iVerticesHF++])AliAODVertex(*secVert4PrAOD);
	      rd = new(aodCharm4ProngRef[i4Prong++])AliAODRecoDecayHF4Prong(*io4Prong);
	      rd->SetSecondaryVtx(v4Prong);
	      v4Prong->SetParent(rd);
	      AddRefs(v4Prong,rd,event,fourTrackArray);
	    }

	    if(io4Prong) {delete io4Prong; io4Prong=NULL;} 
	    if(secVert4PrAOD) {delete secVert4PrAOD; secVert4PrAOD=NULL;} 
	    fourTrackArray->Clear();
	    negtrack2 = 0;

	  } // end loop on negative tracks

          threeTrackArray->Clear();
	  delete vertexp1n1p2;

	}

	postrack2 = 0;
	delete vertexp2n1;

      } // end 2nd loop on positive tracks

      twoTrackArray2->Clear();
      
      // 2nd LOOP  ON  NEGATIVE  TRACKS (for 3 prong -+-)
      for(iTrkN2=iTrkN1+1; iTrkN2<nSeleTrks; iTrkN2++) {

	if(iTrkN2==iTrkP1 || iTrkN2==iTrkP2 || iTrkN2==iTrkN1) continue;

	//if(iTrkN2%1==0) AliDebug(1,Form("    2nd loop on neg: track number %d of %d",iTrkN2,nSeleTrks));  

	// get track from tracks array
	negtrack2 = (AliESDtrack*)seleTrksArray.UncheckedAt(iTrkN2);

	if(negtrack2->Charge()>0) continue;

	if(!TESTBIT(seleFlags[iTrkN2],kBitDispl)) continue;

	if(fMixEvent) {
	  if(evtNumber[iTrkP1]==evtNumber[iTrkN2] || 
	     evtNumber[iTrkN1]==evtNumber[iTrkN2] ||
	     evtNumber[iTrkP1]==evtNumber[iTrkN1]) continue;
	}

	if(isLikeSign2Prong) { // like-sign pair -> have to build only like-sign triplet 
	  if(!fLikeSign3prong) continue;
	  if(postrack1->Charge()<0) { // ok: like-sign triplet (---)
	    isLikeSign3Prong=kTRUE;
	  } else { // not ok
	    continue;
	  }
	} else { // normal triplet (-+-)
	  isLikeSign3Prong=kFALSE;
	}

	// back to primary vertex
	// postrack1->PropagateToDCA(fV1,fBzkG,kVeryBig);
	// negtrack1->PropagateToDCA(fV1,fBzkG,kVeryBig);
	// negtrack2->PropagateToDCA(fV1,fBzkG,kVeryBig);
	SetParametersAtVertex(postrack1,(AliExternalTrackParam*)tracksAtVertex.UncheckedAt(iTrkP1));
	SetParametersAtVertex(negtrack1,(AliExternalTrackParam*)tracksAtVertex.UncheckedAt(iTrkN1));
	SetParametersAtVertex(negtrack2,(AliExternalTrackParam*)tracksAtVertex.UncheckedAt(iTrkN2));
	//printf("********** %d %d %d\n",postrack1->GetID(),negtrack1->GetID(),negtrack2->GetID());

	dcap1n2 = postrack1->GetDCA(negtrack2,fBzkG,xdummy,ydummy);
	if(dcap1n2>dcaMax) { negtrack2=0; continue; }
	dcan1n2 = negtrack1->GetDCA(negtrack2,fBzkG,xdummy,ydummy);
	if(dcan1n2>dcaMax) { negtrack2=0; continue; }

	threeTrackArray->AddAt(negtrack1,0);
	threeTrackArray->AddAt(postrack1,1);
	threeTrackArray->AddAt(negtrack2,2);

	// check invariant mass cuts for D+,Ds,Lc
        massCutOK=kTRUE;
	if(fMassCutBeforeVertexing && f3Prong) 
	  massCutOK = SelectInvMassAndPt(threeTrackArray);

	if(!massCutOK) { 
	  threeTrackArray->Clear();
	  negtrack2=0; 
	  continue; 
	}
	
	// Vertexing
	twoTrackArray2->AddAt(postrack1,0);
	twoTrackArray2->AddAt(negtrack2,1);

	AliAODVertex *vertexp1n2 = ReconstructSecondaryVertex(twoTrackArray2,dispersion);
	if(!vertexp1n2) { 
	  twoTrackArray2->Clear();
	  negtrack2=0; 
	  continue; 
	}

	if(f3Prong) { 
	  AliAODVertex* secVert3PrAOD = ReconstructSecondaryVertex(threeTrackArray,dispersion);
	  io3Prong = Make3Prong(threeTrackArray,event,secVert3PrAOD,dispersion,vertexp1n1,vertexp1n2,dcap1n1,dcap1n2,dcan1n2,ok3Prong);
	  if(ok3Prong) {
	    AliAODVertex *v3Prong = new(verticesHFRef[iVerticesHF++])AliAODVertex(*secVert3PrAOD);
	    if(!isLikeSign3Prong) {
	      rd = new(aodCharm3ProngRef[i3Prong++])AliAODRecoDecayHF3Prong(*io3Prong);
	      rd->SetSecondaryVtx(v3Prong);
	      v3Prong->SetParent(rd);
	      AddRefs(v3Prong,rd,event,threeTrackArray);
	    } else { // isLikeSign3Prong
	      if(fLikeSign3prong){
		rd = new(aodLikeSign3ProngRef[iLikeSign3Prong++])AliAODRecoDecayHF3Prong(*io3Prong);
		rd->SetSecondaryVtx(v3Prong);
		v3Prong->SetParent(rd);
		AddRefs(v3Prong,rd,event,threeTrackArray);
	      }
	    }
	    // Set selection bit for PID
	    SetSelectionBitForPID(fCutsDplustoKpipi,rd,AliRDHFCuts::kDplusPID);
	    SetSelectionBitForPID(fCutsDstoKKpi,rd,AliRDHFCuts::kDsPID);
	    SetSelectionBitForPID(fCutsLctopKpi,rd,AliRDHFCuts::kLcPID);
	  }
	  if(io3Prong) {delete io3Prong; io3Prong=NULL;} 
	  if(secVert3PrAOD) {delete secVert3PrAOD; secVert3PrAOD=NULL;}
	}
	threeTrackArray->Clear();
	negtrack2 = 0;
	delete vertexp1n2;

      } // end 2nd loop on negative tracks
      
      twoTrackArray2->Clear();
      
      negtrack1 = 0;
      delete vertexp1n1; 
    } // end 1st loop on negative tracks
    
    postrack1 = 0;
  }  // end 1st loop on positive tracks


  AliDebug(1,Form(" Total HF vertices in event = %d;",
		  (Int_t)aodVerticesHFTClArr->GetEntriesFast()));
  if(fD0toKpi) {
    AliDebug(1,Form(" D0->Kpi in event = %d;",
		    (Int_t)aodD0toKpiTClArr->GetEntriesFast()));
  }
  if(fJPSItoEle) {
    AliDebug(1,Form(" JPSI->ee in event = %d;",
		    (Int_t)aodJPSItoEleTClArr->GetEntriesFast()));
  }
  if(f3Prong) {
    AliDebug(1,Form(" Charm->3Prong in event = %d;",
		    (Int_t)aodCharm3ProngTClArr->GetEntriesFast()));
  }
  if(f4Prong) {
    AliDebug(1,Form(" Charm->4Prong in event = %d;\n",
		    (Int_t)aodCharm4ProngTClArr->GetEntriesFast()));
  }
  if(fDstar) {
    AliDebug(1,Form(" D*->D0pi in event = %d;\n",
		    (Int_t)aodDstarTClArr->GetEntriesFast()));
  }
  if(fCascades){
    AliDebug(1,Form(" cascades -> v0 + track in event = %d;\n",
		    (Int_t)aodCascadesTClArr->GetEntriesFast()));
  }
  if(fLikeSign) {
    AliDebug(1,Form(" Like-sign 2Prong in event = %d;\n",
		    (Int_t)aodLikeSign2ProngTClArr->GetEntriesFast()));
  }
  if(fLikeSign3prong && f3Prong) {
    AliDebug(1,Form(" Like-sign 3Prong in event = %d;\n",
		    (Int_t)aodLikeSign3ProngTClArr->GetEntriesFast()));
  }
    

  twoTrackArray1->Delete();  delete twoTrackArray1;
  twoTrackArray2->Delete();  delete twoTrackArray2;
  twoTrackArrayCasc->Delete();  delete twoTrackArrayCasc;
  twoTrackArrayV0->Delete();  delete twoTrackArrayV0;
  threeTrackArray->Clear(); 
  threeTrackArray->Delete(); delete threeTrackArray;
  fourTrackArray->Delete();  delete fourTrackArray;
  delete [] seleFlags; seleFlags=NULL;
  if(evtNumber) {delete [] evtNumber; evtNumber=NULL;}
  tracksAtVertex.Delete();

  if(fInputAOD) {
    seleTrksArray.Delete(); 
    if(fAODMap) { delete fAODMap; fAODMap=NULL; }
  }
  

  //printf("Trks: total %d  sele %d\n",fnTrksTotal,fnSeleTrksTotal);

  return;
}
//----------------------------------------------------------------------------
void AliAnalysisVertexingHF::AddRefs(AliAODVertex *v,AliAODRecoDecayHF *rd,
				     const AliVEvent *event,
				     const TObjArray *trkArray) const
{
  // Add the AOD tracks as daughters of the vertex (TRef)
  // Also add the references to the primary vertex and to the cuts

  if(fInputAOD) {
    AddDaughterRefs(v,event,trkArray);
    rd->SetPrimaryVtxRef((AliAODVertex*)event->GetPrimaryVertex());
  }

  /*
  rd->SetListOfCutsRef((TList*)fListOfCuts);
  //fListOfCuts->Print();
  cout<<fListOfCuts<<endl;
  TList *l=(TList*)rd->GetListOfCuts();
  cout<<l<<endl;
  if(l) {l->Print(); }else{printf("error\n");}
  */

  return;
}	
//----------------------------------------------------------------------------
void AliAnalysisVertexingHF::AddDaughterRefs(AliAODVertex *v,
					     const AliVEvent *event,
					     const TObjArray *trkArray) const
{
  // Add the AOD tracks as daughters of the vertex (TRef)

  Int_t nDg = v->GetNDaughters();
  TObject *dg = 0;
  if(nDg) dg = v->GetDaughter(0);
  
  if(dg) return; // daughters already added

  Int_t nTrks = trkArray->GetEntriesFast();

  AliExternalTrackParam *track = 0;
  AliAODTrack *aodTrack = 0;
  Int_t id;

  for(Int_t i=0; i<nTrks; i++) {
    track = (AliExternalTrackParam*)trkArray->UncheckedAt(i);
    id = (Int_t)track->GetID();
    //printf("---> %d\n",id);
    if(id<0) continue; // this track is a AliAODRecoDecay
    aodTrack = (AliAODTrack*)event->GetTrack(fAODMap[id]);
    v->AddDaughter(aodTrack);
  }

  return;
}	
//---------------------------------------------------------------------------
void AliAnalysisVertexingHF::FixReferences(AliAODEvent *aod)  
{
  // Checks that the references to the daughter tracks are properly
  // assigned and reassigns them if needed
  //


  TClonesArray *inputArray=(TClonesArray*)aod->GetList()->FindObject("VerticesHF");
  if(!inputArray) return;

  AliAODTrack *track = 0;
  AliAODVertex *vertex = 0;

  Bool_t needtofix=kFALSE;
  for(Int_t iv=0; iv<inputArray->GetEntriesFast(); iv++) {
    vertex = (AliAODVertex*)inputArray->UncheckedAt(iv);
    for(Int_t id=0; id<vertex->GetNDaughters(); id++) {
      track = (AliAODTrack*)vertex->GetDaughter(id);
      if(!track->GetStatus()) needtofix=kTRUE;
    }
    if(needtofix) break;
  }

  if(!needtofix) return;


  printf("Fixing references\n");

  fAODMapSize = 100000;
  fAODMap = new Int_t[fAODMapSize];

  for(Int_t i=0; i<aod->GetNumberOfTracks(); i++) {
    track = aod->GetTrack(i);

    // skip pure ITS SA tracks
    if(track->GetStatus()&AliESDtrack::kITSpureSA) continue;

    // skip tracks without ITS
    if(!(track->GetStatus()&AliESDtrack::kITSin)) continue;

    // TEMPORARY: check that the cov matrix is there
    Double_t covtest[21];
    if(!track->GetCovarianceXYZPxPyPz(covtest)) continue;
    //

    fAODMap[(Int_t)track->GetID()] = i;
  }


  Int_t ids[4]={-1,-1,-1,-1};
  for(Int_t iv=0; iv<inputArray->GetEntriesFast(); iv++) {
    Bool_t cascade=kFALSE;
    vertex = (AliAODVertex*)inputArray->UncheckedAt(iv);
    Int_t id=0;
    Int_t nDgs = vertex->GetNDaughters();
    for(id=0; id<nDgs; id++) {
      track = (AliAODTrack*)vertex->GetDaughter(id);
      if(track->Charge()==0) {cascade=kTRUE; continue;} // cascade
      ids[id]=(Int_t)track->GetID();
      vertex->RemoveDaughter(track);
    }
    if(cascade) continue;
    for(id=0; id<nDgs; id++) {
      track = aod->GetTrack(fAODMap[ids[id]]);
      vertex->AddDaughter(track);
    }
    
  }

  return;
}
//----------------------------------------------------------------------------
AliAODRecoCascadeHF* AliAnalysisVertexingHF::MakeCascade(
				   TObjArray *twoTrackArray,AliVEvent *event,
				   AliAODVertex *secVert,
				   AliAODRecoDecayHF2Prong *rd2Prong,
				   Double_t dca,
				   Bool_t &okDstar) const
{
  // Make the cascade as a 2Prong decay and check if it passes Dstar
  // reconstruction cuts

  okDstar = kFALSE;

  Bool_t dummy1,dummy2,dummy3;

  // We use Make2Prong to construct the AliAODRecoCascadeHF
  // (which inherits from AliAODRecoDecayHF2Prong) 
  AliAODRecoCascadeHF *theCascade = 
    (AliAODRecoCascadeHF*)Make2Prong(twoTrackArray,event,secVert,dca,
				     dummy1,dummy2,dummy3);
  if(!theCascade) return 0x0;

  // charge
  AliESDtrack *trackPi = (AliESDtrack*)twoTrackArray->UncheckedAt(0);
  theCascade->SetCharge(trackPi->Charge());

  //--- selection cuts
  //
  AliAODRecoCascadeHF *tmpCascade = new AliAODRecoCascadeHF(*theCascade);
  if(fInputAOD){
    Int_t idSoftPi=(Int_t)trackPi->GetID();
    AliAODTrack* trackPiAOD=(AliAODTrack*)event->GetTrack(fAODMap[idSoftPi]);
    tmpCascade->GetSecondaryVtx()->AddDaughter(trackPiAOD);
  }else{
    tmpCascade->GetSecondaryVtx()->AddDaughter(trackPi);
  }
  tmpCascade->GetSecondaryVtx()->AddDaughter(rd2Prong);

  AliAODVertex *primVertexAOD=0;
  if(!fRecoPrimVtxSkippingTrks && !fRmTrksFromPrimVtx) {
    // take event primary vertex
    primVertexAOD = PrimaryVertex(); 
    tmpCascade->SetOwnPrimaryVtx(primVertexAOD);
    rd2Prong->SetOwnPrimaryVtx(primVertexAOD);
  }
  // select D*->D0pi
  if(fDstar) {
    okDstar = (Bool_t)fCutsDStartoKpipi->IsSelected(tmpCascade,AliRDHFCuts::kCandidate);
    if(okDstar) theCascade->SetSelectionBit(AliRDHFCuts::kDstarCuts);
  }
  tmpCascade->GetSecondaryVtx()->RemoveDaughters();
  tmpCascade->UnsetOwnPrimaryVtx(); 
  delete tmpCascade; tmpCascade=NULL;
  if(!fRecoPrimVtxSkippingTrks && !fRmTrksFromPrimVtx && !fMixEvent) {
    rd2Prong->UnsetOwnPrimaryVtx();
  }
  if(primVertexAOD) {delete primVertexAOD; primVertexAOD=NULL;}
  //---

  
  return theCascade;
}


//----------------------------------------------------------------------------
AliAODRecoCascadeHF* AliAnalysisVertexingHF::MakeCascade(
				   TObjArray *twoTrackArray,AliVEvent *event,
				   AliAODVertex *secVert,
				   AliAODv0 *v0,
				   Double_t dca,
				   Bool_t &okCascades) const
{
  //
  // Make the cascade as a 2Prong decay and check if it passes 
  // cascades reconstruction cuts
  
  //  AliDebug(2,Form("         building the cascade"));
  okCascades= kFALSE; 
  Bool_t dummy1,dummy2,dummy3;

  // We use Make2Prong to construct the AliAODRecoCascadeHF
  // (which inherits from AliAODRecoDecayHF2Prong) 
  AliAODRecoCascadeHF *theCascade = 
    (AliAODRecoCascadeHF*)Make2Prong(twoTrackArray,event,secVert,dca,
				     dummy1,dummy2,dummy3);
  if(!theCascade) return 0x0;

  // bachelor track and charge
  AliESDtrack *trackBachelor = (AliESDtrack*)twoTrackArray->UncheckedAt(0);
  theCascade->SetCharge(trackBachelor->Charge());

  //--- selection cuts
  //

  AliAODRecoCascadeHF *tmpCascade = new AliAODRecoCascadeHF(*theCascade);  
  if(fInputAOD){
    Int_t idBachelor=(Int_t)trackBachelor->GetID();
    AliAODTrack* trackBachelorAOD=(AliAODTrack*)event->GetTrack(fAODMap[idBachelor]);
    tmpCascade->GetSecondaryVtx()->AddDaughter(trackBachelorAOD);
  }else{
    tmpCascade->GetSecondaryVtx()->AddDaughter(trackBachelor);
  }
  tmpCascade->GetSecondaryVtx()->AddDaughter(v0);
  
  AliAODVertex *primVertexAOD=0;
  if(!fRecoPrimVtxSkippingTrks && !fRmTrksFromPrimVtx) {
    // take event primary vertex
    primVertexAOD = PrimaryVertex(); 
    if(!primVertexAOD) primVertexAOD = (AliAODVertex*)event->GetPrimaryVertex(); 
    tmpCascade->SetOwnPrimaryVtx(primVertexAOD);
  }
  
  // select Cascades
  if(fCascades && fInputAOD){
    okCascades = (Bool_t)fCutsLctoV0->IsSelected(tmpCascade,AliRDHFCuts::kCandidate);
  }
  else { 
    //AliDebug(2,Form("The cascade is contructed from ESDs, no cuts are applied")); 
    okCascades=kTRUE; 
  }// no cuts implemented from ESDs
  tmpCascade->GetSecondaryVtx()->RemoveDaughters();
  tmpCascade->UnsetOwnPrimaryVtx(); 
  delete tmpCascade; tmpCascade=NULL;
  if(primVertexAOD) {delete primVertexAOD; primVertexAOD=NULL;}
  //---
  
  return theCascade;
}

//-----------------------------------------------------------------------------
AliAODRecoDecayHF2Prong *AliAnalysisVertexingHF::Make2Prong(
				   TObjArray *twoTrackArray,AliVEvent *event,
				   AliAODVertex *secVert,Double_t dca,
				   Bool_t &okD0,Bool_t &okJPSI,
				   Bool_t &okD0fromDstar) const
{
  // Make 2Prong candidates and check if they pass D0toKpi or BtoJPSI
  // reconstruction cuts
  // G.E.Bruno (J/psi), A.Dainese (D0->Kpi)

  okD0=kFALSE; okJPSI=kFALSE; okD0fromDstar=kFALSE;

  Double_t px[2],py[2],pz[2],d0[2],d0err[2];
  AliESDtrack *postrack = (AliESDtrack*)twoTrackArray->UncheckedAt(0);
  AliESDtrack *negtrack = (AliESDtrack*)twoTrackArray->UncheckedAt(1);

  // propagate tracks to secondary vertex, to compute inv. mass
  postrack->PropagateToDCA(secVert,fBzkG,kVeryBig);
  negtrack->PropagateToDCA(secVert,fBzkG,kVeryBig);

  Double_t momentum[3];
  postrack->GetPxPyPz(momentum);
  px[0] = momentum[0]; py[0] = momentum[1]; pz[0] = momentum[2]; 
  negtrack->GetPxPyPz(momentum);
  px[1] = momentum[0]; py[1] = momentum[1]; pz[1] = momentum[2]; 


  // invariant mass cut (try to improve coding here..)
  Bool_t okMassCut=kFALSE;
  if(!okMassCut && fD0toKpi)   if(SelectInvMassAndPt(0,2,px,py,pz)) okMassCut=kTRUE;
  if(!okMassCut && fJPSItoEle) if(SelectInvMassAndPt(1,2,px,py,pz)) okMassCut=kTRUE;
  if(!okMassCut && fDstar)     if(SelectInvMassAndPt(3,2,px,py,pz)) okMassCut=kTRUE;
  if(!okMassCut && fCascades)  if(SelectInvMassAndPt(5,2,px,py,pz)) okMassCut=kTRUE;
  if(!okMassCut) {
    //AliDebug(2," candidate didn't pass mass cut");
    return 0x0;    
  }
  // primary vertex to be used by this candidate
  AliAODVertex *primVertexAOD  = PrimaryVertex(twoTrackArray,event);
  if(!primVertexAOD) return 0x0;


  Double_t d0z0[2],covd0z0[3];
  postrack->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
  d0[0] = d0z0[0];
  d0err[0] = TMath::Sqrt(covd0z0[0]);
  negtrack->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
  d0[1] = d0z0[0];
  d0err[1] = TMath::Sqrt(covd0z0[0]);
  
  // create the object AliAODRecoDecayHF2Prong
  AliAODRecoDecayHF2Prong *the2Prong = new AliAODRecoDecayHF2Prong(secVert,px,py,pz,d0,d0err,dca);
  the2Prong->SetOwnPrimaryVtx(primVertexAOD);
  UShort_t id[2]={(UShort_t)postrack->GetID(),(UShort_t)negtrack->GetID()};
  the2Prong->SetProngIDs(2,id);
  delete primVertexAOD; primVertexAOD=NULL;

 
  if(postrack->Charge()!=0 && negtrack->Charge()!=0) { // don't apply these cuts if it's a Dstar 

    // Add daughter references already here
    if(fInputAOD) AddDaughterRefs(secVert,(AliAODEvent*)event,twoTrackArray);

    // select D0->Kpi
    if(fD0toKpi)   {
      okD0 = (Bool_t)fCutsD0toKpi->IsSelected(the2Prong,AliRDHFCuts::kCandidate,(AliAODEvent*)event);
      if(okD0) the2Prong->SetSelectionBit(AliRDHFCuts::kD0toKpiCuts);
    }
    //if(fDebug && fD0toKpi) printf("   %d\n",(Int_t)okD0);
    // select J/psi from B
    if(fJPSItoEle)   {
      okJPSI = (Bool_t)fCutsJpsitoee->IsSelected(the2Prong,AliRDHFCuts::kCandidate);
    }
    //if(fDebug && fJPSItoEle) printf("   %d\n",(Int_t)okJPSI);
    // select D0->Kpi from Dstar
    if(fDstar)   {
      okD0fromDstar = (Bool_t)fCutsDStartoKpipi->IsD0FromDStarSelected(the2Prong->Pt(),the2Prong,AliRDHFCuts::kCandidate);
      if(okD0fromDstar) the2Prong->SetSelectionBit(AliRDHFCuts::kD0fromDstarCuts);
    }
    //if(fDebug && fDstar) printf("   %d\n",(Int_t)okD0fromDstar);
  }


  // remove the primary vertex (was used only for selection)
  if(!fRecoPrimVtxSkippingTrks && !fRmTrksFromPrimVtx && !fMixEvent) {
    the2Prong->UnsetOwnPrimaryVtx();
  }
  
  // get PID info from ESD
  Double_t esdpid0[5]={0.,0.,0.,0.,0.};
  if(postrack->GetStatus()&AliESDtrack::kESDpid) postrack->GetESDpid(esdpid0);
  Double_t esdpid1[5]={0.,0.,0.,0.,0.};
  if(negtrack->GetStatus()&AliESDtrack::kESDpid) negtrack->GetESDpid(esdpid1);
  Double_t esdpid[10];
  for(Int_t i=0;i<5;i++) {
    esdpid[i]   = esdpid0[i];
    esdpid[5+i] = esdpid1[i];
  }
  the2Prong->SetPID(2,esdpid);

  return the2Prong;  
}
//----------------------------------------------------------------------------
AliAODRecoDecayHF3Prong* AliAnalysisVertexingHF::Make3Prong(
                             TObjArray *threeTrackArray,AliVEvent *event,
			     AliAODVertex *secVert,Double_t dispersion,
			     const AliAODVertex *vertexp1n1,const AliAODVertex *vertexp2n1,
			     Double_t dcap1n1,Double_t dcap2n1,Double_t dcap1p2,
			     Bool_t &ok3Prong) const
{
  // Make 3Prong candidates and check if they pass Dplus or Ds or Lambdac
  // reconstruction cuts 
  // E.Bruna, F.Prino


  ok3Prong=kFALSE;
  if(!secVert || !vertexp1n1 || !vertexp2n1) return 0x0; 

  Double_t px[3],py[3],pz[3],d0[3],d0err[3];
  Double_t momentum[3];


  AliESDtrack *postrack1 = (AliESDtrack*)threeTrackArray->UncheckedAt(0);
  AliESDtrack *negtrack  = (AliESDtrack*)threeTrackArray->UncheckedAt(1);
  AliESDtrack *postrack2 = (AliESDtrack*)threeTrackArray->UncheckedAt(2);

  postrack1->PropagateToDCA(secVert,fBzkG,kVeryBig);
  negtrack->PropagateToDCA(secVert,fBzkG,kVeryBig);
  postrack2->PropagateToDCA(secVert,fBzkG,kVeryBig);
  postrack1->GetPxPyPz(momentum);
  px[0] = momentum[0]; py[0] = momentum[1]; pz[0] = momentum[2]; 
  negtrack->GetPxPyPz(momentum);
  px[1] = momentum[0]; py[1] = momentum[1]; pz[1] = momentum[2]; 
  postrack2->GetPxPyPz(momentum);
  px[2] = momentum[0]; py[2] = momentum[1]; pz[2] = momentum[2]; 

  // invariant mass cut for D+, Ds, Lc
  Bool_t okMassCut=kFALSE;
  if(fMassCutBeforeVertexing) okMassCut=kTRUE; // mass cut already done and passed 
  if(!okMassCut && f3Prong) if(SelectInvMassAndPt(2,3,px,py,pz)) okMassCut=kTRUE;
  if(!okMassCut) {
    //AliDebug(2," candidate didn't pass mass cut");
    return 0x0;    
  }

  // primary vertex to be used by this candidate
  AliAODVertex *primVertexAOD  = PrimaryVertex(threeTrackArray,event);
  if(!primVertexAOD) return 0x0;

  Double_t d0z0[2],covd0z0[3];
  postrack1->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
  d0[0]=d0z0[0];
  d0err[0] = TMath::Sqrt(covd0z0[0]);
  negtrack->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
  d0[1]=d0z0[0];
  d0err[1] = TMath::Sqrt(covd0z0[0]);
  postrack2->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
  d0[2]=d0z0[0];
  d0err[2] = TMath::Sqrt(covd0z0[0]);


  // create the object AliAODRecoDecayHF3Prong
  Double_t pos[3]; primVertexAOD->GetXYZ(pos);
  Double_t dca[3]={dcap1n1,dcap2n1,dcap1p2};
  Double_t dist12=TMath::Sqrt((vertexp1n1->GetX()-pos[0])*(vertexp1n1->GetX()-pos[0])+(vertexp1n1->GetY()-pos[1])*(vertexp1n1->GetY()-pos[1])+(vertexp1n1->GetZ()-pos[2])*(vertexp1n1->GetZ()-pos[2]));
  Double_t dist23=TMath::Sqrt((vertexp2n1->GetX()-pos[0])*(vertexp2n1->GetX()-pos[0])+(vertexp2n1->GetY()-pos[1])*(vertexp2n1->GetY()-pos[1])+(vertexp2n1->GetZ()-pos[2])*(vertexp2n1->GetZ()-pos[2]));
  Short_t charge=(Short_t)(postrack1->Charge()+postrack2->Charge()+negtrack->Charge());
  AliAODRecoDecayHF3Prong *the3Prong = new AliAODRecoDecayHF3Prong(secVert,px,py,pz,d0,d0err,dca,dispersion,dist12,dist23,charge);
  the3Prong->SetOwnPrimaryVtx(primVertexAOD);
  UShort_t id[3]={(UShort_t)postrack1->GetID(),(UShort_t)negtrack->GetID(),(UShort_t)postrack2->GetID()};
  the3Prong->SetProngIDs(3,id);

  delete primVertexAOD; primVertexAOD=NULL;

  // Add daughter references already here
  if(fInputAOD) AddDaughterRefs(secVert,(AliAODEvent*)event,threeTrackArray);

  // select D+->Kpipi, Ds->KKpi, Lc->pKpi
  if(f3Prong) {
    ok3Prong = kFALSE;
    
    if(fCutsDplustoKpipi->IsSelected(the3Prong,AliRDHFCuts::kCandidate,(AliAODEvent*)event)) {
      ok3Prong = kTRUE;
      the3Prong->SetSelectionBit(AliRDHFCuts::kDplusCuts);
    }
    if(fCutsDstoKKpi->IsSelected(the3Prong,AliRDHFCuts::kCandidate,(AliAODEvent*)event)) {
      ok3Prong = kTRUE;
      the3Prong->SetSelectionBit(AliRDHFCuts::kDsCuts);
    }
    if(fCutsLctopKpi->IsSelected(the3Prong,AliRDHFCuts::kCandidate,(AliAODEvent*)event)) {
      ok3Prong = kTRUE;
      the3Prong->SetSelectionBit(AliRDHFCuts::kLcCuts);
    } 
  }
  //if(fDebug) printf("ok3Prong: %d\n",(Int_t)ok3Prong);

  if(!fRecoPrimVtxSkippingTrks && !fRmTrksFromPrimVtx && !fMixEvent) {
    the3Prong->UnsetOwnPrimaryVtx();
  }

  // get PID info from ESD
  Double_t esdpid0[5]={0.,0.,0.,0.,0.};
  if(postrack1->GetStatus()&AliESDtrack::kESDpid) postrack1->GetESDpid(esdpid0);
  Double_t esdpid1[5]={0.,0.,0.,0.,0.};
  if(negtrack->GetStatus()&AliESDtrack::kESDpid) negtrack->GetESDpid(esdpid1);
  Double_t esdpid2[5]={0.,0.,0.,0.,0.};
  if(postrack2->GetStatus()&AliESDtrack::kESDpid) postrack2->GetESDpid(esdpid2);
  
  Double_t esdpid[15];
  for(Int_t i=0;i<5;i++) {
    esdpid[i]    = esdpid0[i];
    esdpid[5+i]  = esdpid1[i];
    esdpid[10+i] = esdpid2[i];
  }
  the3Prong->SetPID(3,esdpid);

  return the3Prong;
}
//----------------------------------------------------------------------------
AliAODRecoDecayHF4Prong* AliAnalysisVertexingHF::Make4Prong(
                             TObjArray *fourTrackArray,AliVEvent *event,
                             AliAODVertex *secVert,
                             const AliAODVertex *vertexp1n1,
                             const AliAODVertex *vertexp1n1p2,
                             Double_t dcap1n1,Double_t dcap1n2,
			     Double_t dcap2n1,Double_t dcap2n2,
                             Bool_t &ok4Prong) const
{
  // Make 4Prong candidates and check if they pass D0toKpipipi
  // reconstruction cuts
  // G.E.Bruno, R.Romita

  ok4Prong=kFALSE;
  if(!secVert || !vertexp1n1 || !vertexp1n1p2) return 0x0; 

  Double_t px[4],py[4],pz[4],d0[4],d0err[4];//d0z[3];

  AliESDtrack *postrack1 = (AliESDtrack*)fourTrackArray->UncheckedAt(0);
  AliESDtrack *negtrack1 = (AliESDtrack*)fourTrackArray->UncheckedAt(1);
  AliESDtrack *postrack2 = (AliESDtrack*)fourTrackArray->UncheckedAt(2);
  AliESDtrack *negtrack2 = (AliESDtrack*)fourTrackArray->UncheckedAt(3);

  postrack1->PropagateToDCA(secVert,fBzkG,kVeryBig);
  negtrack1->PropagateToDCA(secVert,fBzkG,kVeryBig);
  postrack2->PropagateToDCA(secVert,fBzkG,kVeryBig);
  negtrack2->PropagateToDCA(secVert,fBzkG,kVeryBig);

  Double_t momentum[3];
  postrack1->GetPxPyPz(momentum);
  px[0] = momentum[0]; py[0] = momentum[1]; pz[0] = momentum[2];
  negtrack1->GetPxPyPz(momentum);
  px[1] = momentum[0]; py[1] = momentum[1]; pz[1] = momentum[2];
  postrack2->GetPxPyPz(momentum);
  px[2] = momentum[0]; py[2] = momentum[1]; pz[2] = momentum[2];
  negtrack2->GetPxPyPz(momentum);
  px[3] = momentum[0]; py[3] = momentum[1]; pz[3] = momentum[2];

  // invariant mass cut for rho or D0 (try to improve coding here..)
  Bool_t okMassCut=kFALSE;
  if(fMassCutBeforeVertexing) okMassCut=kTRUE; // mass cut already done and passed 
  if(!okMassCut && !(fCutsD0toKpipipi->GetUsePID())) {      //no PID, to be implemented with PID
    if(SelectInvMassAndPt(4,4,px,py,pz)) okMassCut=kTRUE;
  }
  if(!okMassCut) {
    //if(fDebug) printf(" candidate didn't pass mass cut\n");
    //printf(" candidate didn't pass mass cut\n");
    return 0x0;
  }

  // primary vertex to be used by this candidate
  AliAODVertex *primVertexAOD  = PrimaryVertex(fourTrackArray,event);
  if(!primVertexAOD) return 0x0;

  Double_t d0z0[2],covd0z0[3];
  postrack1->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
  d0[0]=d0z0[0];
  d0err[0] = TMath::Sqrt(covd0z0[0]);
  negtrack1->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
  d0[1]=d0z0[0];
  d0err[1] = TMath::Sqrt(covd0z0[0]);
  postrack2->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
  d0[2]=d0z0[0];
  d0err[2] = TMath::Sqrt(covd0z0[0]);
  negtrack2->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
  d0[3]=d0z0[0];
  d0err[3] = TMath::Sqrt(covd0z0[0]);


  // create the object AliAODRecoDecayHF4Prong
  Double_t pos[3]; primVertexAOD->GetXYZ(pos);
  Double_t dca[6]={dcap1n1,0.,dcap1n2,dcap2n1,0.,dcap2n2};
  Double_t dist12=TMath::Sqrt((vertexp1n1->GetX()-pos[0])*(vertexp1n1->GetX()-pos[0])+(vertexp1n1->GetY()-pos[1])*(vertexp1n1->GetY()-pos[1])+(vertexp1n1->GetZ()-pos[2])*(vertexp1n1->GetZ()-pos[2]));
  Double_t dist3=TMath::Sqrt((vertexp1n1p2->GetX()-pos[0])*(vertexp1n1p2->GetX()-pos[0])+(vertexp1n1p2->GetY()-pos[1])*(vertexp1n1p2->GetY()-pos[1])+(vertexp1n1p2->GetZ()-pos[2])*(vertexp1n1p2->GetZ()-pos[2]));
  Double_t dist4=TMath::Sqrt((secVert->GetX()-pos[0])*(secVert->GetX()-pos[0])+(secVert->GetY()-pos[1])*(secVert->GetY()-pos[1])+(secVert->GetZ()-pos[2])*(secVert->GetZ()-pos[2]));
  Short_t charge=0;
  AliAODRecoDecayHF4Prong *the4Prong = new AliAODRecoDecayHF4Prong(secVert,px,py,pz,d0,d0err,dca,dist12,dist3,dist4,charge);
  the4Prong->SetOwnPrimaryVtx(primVertexAOD);
  UShort_t id[4]={(UShort_t)postrack1->GetID(),(UShort_t)negtrack1->GetID(),(UShort_t)postrack2->GetID(),(UShort_t)negtrack2->GetID()};
  the4Prong->SetProngIDs(4,id);

  delete primVertexAOD; primVertexAOD=NULL;

  ok4Prong=(Bool_t)fCutsD0toKpipipi->IsSelected(the4Prong,AliRDHFCuts::kCandidate);


  if(!fRecoPrimVtxSkippingTrks && !fRmTrksFromPrimVtx && !fMixEvent) {
    the4Prong->UnsetOwnPrimaryVtx();
  }

 
  // get PID info from ESD
  Double_t esdpid0[5]={0.,0.,0.,0.,0.};
  if(postrack1->GetStatus()&AliESDtrack::kESDpid) postrack1->GetESDpid(esdpid0);
  Double_t esdpid1[5]={0.,0.,0.,0.,0.};
  if(negtrack1->GetStatus()&AliESDtrack::kESDpid) negtrack1->GetESDpid(esdpid1);
  Double_t esdpid2[5]={0.,0.,0.,0.,0.};
  if(postrack2->GetStatus()&AliESDtrack::kESDpid) postrack2->GetESDpid(esdpid2);
  Double_t esdpid3[5]={0.,0.,0.,0.,0.};
  if(negtrack2->GetStatus()&AliESDtrack::kESDpid) negtrack2->GetESDpid(esdpid3);

  Double_t esdpid[20];
  for(Int_t i=0;i<5;i++) {
    esdpid[i]    = esdpid0[i];
    esdpid[5+i]  = esdpid1[i];
    esdpid[10+i] = esdpid2[i];
    esdpid[15+i] = esdpid3[i];
  }
  the4Prong->SetPID(4,esdpid);
  
  return the4Prong;
}
//-----------------------------------------------------------------------------
AliAODVertex* AliAnalysisVertexingHF::PrimaryVertex(const TObjArray *trkArray,
						    AliVEvent *event) const
{
  // Returns primary vertex to be used for this candidate
 
  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;


  if(!fRecoPrimVtxSkippingTrks && !fRmTrksFromPrimVtx) { 
    // primary vertex from the input event
    
    vertexESD = new AliESDVertex(*fV1);

  } else {
    // primary vertex specific to this candidate

    Int_t nTrks = trkArray->GetEntriesFast();
    AliVertexerTracks *vertexer = new AliVertexerTracks(event->GetMagneticField());

    if(fRecoPrimVtxSkippingTrks) { 
      // recalculating the vertex
      
      if(strstr(fV1->GetTitle(),"VertexerTracksWithConstraint")) {
	Float_t diamondcovxy[3];
	event->GetDiamondCovXY(diamondcovxy);
	Double_t pos[3]={event->GetDiamondX(),event->GetDiamondY(),0.};
	Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
	AliESDVertex *diamond = new AliESDVertex(pos,cov,1.,1);
	vertexer->SetVtxStart(diamond);
	delete diamond; diamond=NULL;
	if(strstr(fV1->GetTitle(),"VertexerTracksWithConstraintOnlyFitter")) 
	  vertexer->SetOnlyFitter();
      }
      Int_t skipped[1000];
      Int_t nTrksToSkip=0,id;
      AliExternalTrackParam *t = 0;
      for(Int_t i=0; i<nTrks; i++) {
	t = (AliExternalTrackParam*)trkArray->UncheckedAt(i);
	id = (Int_t)t->GetID();
	if(id<0) continue;
	skipped[nTrksToSkip++] = id;
      }
      // TEMPORARY FIX
      // For AOD, skip also tracks without covariance matrix
      if(fInputAOD) {
	Double_t covtest[21];
	for(Int_t j=0; j<event->GetNumberOfTracks(); j++) {
	  AliVTrack *vtrack = (AliVTrack*)event->GetTrack(j);
	  if(!vtrack->GetCovarianceXYZPxPyPz(covtest)) {
	    id = (Int_t)vtrack->GetID();
	    if(id<0) continue;
	    skipped[nTrksToSkip++] = id;
	  }
	}
      }
      for(Int_t ijk=nTrksToSkip; ijk<1000; ijk++) skipped[ijk]=-1;
      //
      vertexer->SetSkipTracks(nTrksToSkip,skipped);
      vertexESD = (AliESDVertex*)vertexer->FindPrimaryVertex(event); 
      
    } else if(fRmTrksFromPrimVtx && nTrks>0) { 
      // removing the prongs tracks
      
      TObjArray rmArray(nTrks);
      UShort_t *rmId = new UShort_t[nTrks];
      AliESDtrack *esdTrack = 0;
      AliESDtrack *t = 0;
      for(Int_t i=0; i<nTrks; i++) {
	t = (AliESDtrack*)trkArray->UncheckedAt(i);
	esdTrack = new AliESDtrack(*t);
	rmArray.AddLast(esdTrack);
	if(esdTrack->GetID()>=0) {
	  rmId[i]=(UShort_t)esdTrack->GetID();
	} else {
	  rmId[i]=9999;
	}
      }
      Float_t diamondxy[2]={event->GetDiamondX(),event->GetDiamondY()};
      vertexESD = vertexer->RemoveTracksFromVertex(fV1,&rmArray,rmId,diamondxy);
      delete [] rmId; rmId=NULL;
      rmArray.Delete();
      
    }

    if(!vertexESD) return vertexAOD;
    if(vertexESD->GetNContributors()<=0) { 
      //AliDebug(2,"vertexing failed"); 
      delete vertexESD; vertexESD=NULL;
      return vertexAOD;
    }

    delete vertexer; vertexer=NULL;

  }

  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  delete vertexESD; vertexESD=NULL;

  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF);

  return vertexAOD;
}
//-----------------------------------------------------------------------------
void AliAnalysisVertexingHF::PrintStatus() const {
  // Print parameters being used

  //printf("Preselections:\n");
  //   fTrackFilter->Dump();
  if(fSecVtxWithKF) {
    printf("Secondary vertex with Kalman filter package (AliKFParticle)\n");
  } else {
    printf("Secondary vertex with AliVertexerTracks\n");
  }
  if(fRecoPrimVtxSkippingTrks) printf("RecoPrimVtxSkippingTrks\n");
  if(fRmTrksFromPrimVtx) printf("RmTrksFromPrimVtx\n");
  if(fD0toKpi) {
    printf("Reconstruct D0->Kpi candidates with cuts:\n");
    if(fCutsD0toKpi) fCutsD0toKpi->PrintAll();
  }
  if(fDstar) {
    printf("Reconstruct D*->D0pi candidates with cuts:\n");
    if(fFindVertexForDstar) {
      printf("    Reconstruct a secondary vertex for the D*\n");
    } else {
      printf("    Assume the D* comes from the primary vertex\n");
    }
    if(fCutsDStartoKpipi) fCutsDStartoKpipi->PrintAll();
  }
  if(fJPSItoEle) {
    printf("Reconstruct J/psi from B candidates with cuts:\n");
    if(fCutsJpsitoee) fCutsJpsitoee->PrintAll();
  }
  if(f3Prong) {
    printf("Reconstruct 3 prong candidates.\n");
    printf("  D+->Kpipi cuts:\n");
    if(fCutsDplustoKpipi) fCutsDplustoKpipi->PrintAll();
    printf("  Ds->KKpi cuts:\n");
    if(fCutsDstoKKpi) fCutsDstoKKpi->PrintAll();
    printf("  Lc->pKpi cuts:\n");
    if(fCutsLctopKpi) fCutsLctopKpi->PrintAll();
  }
  if(f4Prong) {
    printf("Reconstruct 4 prong candidates.\n");
    printf("  D0->Kpipipi cuts:\n");
    if(fCutsD0toKpipipi) fCutsD0toKpipipi->PrintAll();
  }
  if(fCascades) {
    printf("Reconstruct cascades candidates formed with v0s.\n");
    printf("  Lc -> k0s P & Lc -> L Pi cuts:\n");
    if(fCutsLctoV0) fCutsLctoV0->PrintAll();
  }

  return;
}
//-----------------------------------------------------------------------------
AliAODVertex* AliAnalysisVertexingHF::ReconstructSecondaryVertex(TObjArray *trkArray,
								 Double_t &dispersion,Bool_t useTRefArray) const
{
  // Secondary vertex reconstruction with AliVertexerTracks or AliKFParticle

  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;

  if(!fSecVtxWithKF) { // AliVertexerTracks

    AliVertexerTracks *vertexer = new AliVertexerTracks(fBzkG);
    vertexer->SetVtxStart(fV1);
    vertexESD = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
    delete vertexer; vertexer=NULL;

    if(!vertexESD) return vertexAOD;

    if(vertexESD->GetNContributors()!=trkArray->GetEntriesFast()) { 
      //AliDebug(2,"vertexing failed"); 
      delete vertexESD; vertexESD=NULL;
      return vertexAOD;
    }
    
    Double_t vertRadius2=vertexESD->GetXv()*vertexESD->GetXv()+vertexESD->GetYv()*vertexESD->GetYv();
    if(vertRadius2>8.){
      // vertex outside beam pipe, reject candidate to avoid propagation through material
      delete vertexESD; vertexESD=NULL;
      return vertexAOD;
    }

  } else { // Kalman Filter vertexer (AliKFParticle)

    AliKFParticle::SetField(fBzkG);

    AliKFVertex vertexKF;

    Int_t nTrks = trkArray->GetEntriesFast();
    for(Int_t i=0; i<nTrks; i++) {
      AliESDtrack *esdTrack = (AliESDtrack*)trkArray->At(i);
      AliKFParticle daughterKF(*esdTrack,211);
      vertexKF.AddDaughter(daughterKF);
    }
    vertexESD = new AliESDVertex(vertexKF.Parameters(),
				 vertexKF.CovarianceMatrix(),
				 vertexKF.GetChi2(),
				 vertexKF.GetNContributors());

  }

  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  dispersion = vertexESD->GetDispersion();
  delete vertexESD; vertexESD=NULL;

  Int_t nprongs= (useTRefArray ? 0 : trkArray->GetEntriesFast());
  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,nprongs);

  return vertexAOD;
}
//-----------------------------------------------------------------------------
Bool_t AliAnalysisVertexingHF::SelectInvMassAndPt(TObjArray *trkArray) const {
  // Invariant mass cut on tracks

  Int_t nProngs=trkArray->GetEntriesFast();
  Int_t retval=kFALSE;

  Double_t momentum[3];
  Double_t px3[3],py3[3],pz3[3];
  Double_t px4[4],py4[4],pz4[4];
  AliESDtrack *postrack1=0;
  AliESDtrack *postrack2=0;
  AliESDtrack *negtrack1=0;
  AliESDtrack *negtrack2=0;

  switch(nProngs) {
  case 3:
    // invariant mass cut for D+, Ds, Lc
    postrack1 = (AliESDtrack*)trkArray->UncheckedAt(0);
    negtrack1  = (AliESDtrack*)trkArray->UncheckedAt(1);
    postrack2 = (AliESDtrack*)trkArray->UncheckedAt(2);
    postrack1->GetPxPyPz(momentum);
    px3[0] = momentum[0]; py3[0] = momentum[1]; pz3[0] = momentum[2]; 
    negtrack1->GetPxPyPz(momentum);
    px3[1] = momentum[0]; py3[1] = momentum[1]; pz3[1] = momentum[2]; 
    postrack2->GetPxPyPz(momentum);
    px3[2] = momentum[0]; py3[2] = momentum[1]; pz3[2] = momentum[2]; 
    retval = SelectInvMassAndPt(2,3,px3,py3,pz3);
    break;
  case 4:
    // invariant mass cut for D0->4prong
    postrack1 = (AliESDtrack*)trkArray->UncheckedAt(0);
    negtrack1 = (AliESDtrack*)trkArray->UncheckedAt(1);
    postrack2 = (AliESDtrack*)trkArray->UncheckedAt(2);
    negtrack2 = (AliESDtrack*)trkArray->UncheckedAt(3);
    postrack1->GetPxPyPz(momentum);
    px4[0] = momentum[0]; py4[0] = momentum[1]; pz4[0] = momentum[2];
    negtrack1->GetPxPyPz(momentum);
    px4[1] = momentum[0]; py4[1] = momentum[1]; pz4[1] = momentum[2];
    postrack2->GetPxPyPz(momentum);
    px4[2] = momentum[0]; py4[2] = momentum[1]; pz4[2] = momentum[2];
    negtrack2->GetPxPyPz(momentum);
    px4[3] = momentum[0]; py4[3] = momentum[1]; pz4[3] = momentum[2];
    retval = SelectInvMassAndPt(4,4,px4,py4,pz4);
    break;
  default:
    printf("SelectInvMassAndPt(): wrong decay selection\n");
    break;
  }

  return retval;
}
//-----------------------------------------------------------------------------
Bool_t AliAnalysisVertexingHF::SelectInvMassAndPt(Int_t decay,
					     Int_t nprongs,
					     Double_t *px,
					     Double_t *py,
					     Double_t *pz) const {
  // Check invariant mass cut and pt candidate cut

  UInt_t pdg2[2],pdg3[3],pdg4[4];
  Double_t mPDG,minv;
  Double_t minPt=0;

  Bool_t retval=kFALSE;
  switch (decay) 
    { 
    case 0:                  // D0->Kpi
      fMassCalc2->SetPxPyPzProngs(nprongs,px,py,pz);
      // pt cut
      minPt=fCutsD0toKpi->GetMinPtCandidate();
      if(minPt>0.1) 
	if(fMassCalc2->Pt2() < minPt*minPt) break;
      // mass cut
      pdg2[0]=211; pdg2[1]=321;
      mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();
      minv = fMassCalc2->InvMass(nprongs,pdg2);
      if(TMath::Abs(minv-mPDG)<fCutsD0toKpi->GetMassCut()) retval=kTRUE;
      pdg2[0]=321; pdg2[1]=211;
      minv = fMassCalc2->InvMass(nprongs,pdg2);
      if(TMath::Abs(minv-mPDG)<fCutsD0toKpi->GetMassCut()) retval=kTRUE;
      break;
    case 1:                  // JPSI->ee
      fMassCalc2->SetPxPyPzProngs(nprongs,px,py,pz);
      // pt cut
      minPt=fCutsJpsitoee->GetMinPtCandidate();
      if(minPt>0.1) 
	if(fMassCalc2->Pt2() < minPt*minPt) break;
      // mass cut
      pdg2[0]=11; pdg2[1]=11;
      mPDG=TDatabasePDG::Instance()->GetParticle(443)->Mass();
      minv = fMassCalc2->InvMass(nprongs,pdg2);
      if(TMath::Abs(minv-mPDG)<fCutsJpsitoee->GetMassCut()) retval=kTRUE;
      break;
    case 2:                  // D+->Kpipi
      fMassCalc3->SetPxPyPzProngs(nprongs,px,py,pz);
      // pt cut
      minPt=TMath::Min(fCutsDplustoKpipi->GetMinPtCandidate(),fCutsDstoKKpi->GetMinPtCandidate());
      minPt=TMath::Min(minPt,fCutsLctopKpi->GetMinPtCandidate());
      if(minPt>0.1) 
	if(fMassCalc3->Pt2() < minPt*minPt) break;
      // mass cut
      pdg3[0]=211; pdg3[1]=321; pdg3[2]=211;
      mPDG=TDatabasePDG::Instance()->GetParticle(411)->Mass();
      minv = fMassCalc3->InvMass(nprongs,pdg3);
      if(TMath::Abs(minv-mPDG)<fCutsDplustoKpipi->GetMassCut()) retval=kTRUE;
                            // Ds+->KKpi
      pdg3[0]=321; pdg3[1]=321; pdg3[2]=211;
      mPDG=TDatabasePDG::Instance()->GetParticle(431)->Mass();
      minv = fMassCalc3->InvMass(nprongs,pdg3);
      if(TMath::Abs(minv-mPDG)<fCutsDstoKKpi->GetMassCut()) retval=kTRUE;
      pdg3[0]=211; pdg3[1]=321; pdg3[2]=321;
      minv = fMassCalc3->InvMass(nprongs,pdg3);
      if(TMath::Abs(minv-mPDG)<fCutsDstoKKpi->GetMassCut()) retval=kTRUE;
                            // Lc->pKpi
      pdg3[0]=2212; pdg3[1]=321; pdg3[2]=211;
      mPDG=TDatabasePDG::Instance()->GetParticle(4122)->Mass();
      minv = fMassCalc3->InvMass(nprongs,pdg3);
      if(TMath::Abs(minv-mPDG)<fCutsLctopKpi->GetMassCut()) retval=kTRUE;
      pdg3[0]=211; pdg3[1]=321; pdg3[2]=2212;
      minv = fMassCalc3->InvMass(nprongs,pdg3);
      if(TMath::Abs(minv-mPDG)<fCutsLctopKpi->GetMassCut()) retval=kTRUE;
      break;
    case 3:                  // D*->D0pi
      fMassCalc2->SetPxPyPzProngs(nprongs,px,py,pz);
      // pt cut
      minPt=fCutsDStartoKpipi->GetMinPtCandidate();
      if(minPt>0.1) 
	if(fMassCalc2->Pt2() < minPt*minPt) break;
      // mass cut
      pdg2[0]=211; pdg2[1]=421; // in twoTrackArrayCasc we put the pion first
      mPDG=TDatabasePDG::Instance()->GetParticle(413)->Mass();
      minv = fMassCalc2->InvMass(nprongs,pdg2);
      if(TMath::Abs(minv-mPDG)<fCutsDStartoKpipi->GetMassCut()) retval=kTRUE;
      break;
    case 4:                 // D0->Kpipipi without PID
      fMassCalc4->SetPxPyPzProngs(nprongs,px,py,pz);
      // pt cut
      minPt=fCutsD0toKpipipi->GetMinPtCandidate();
      if(minPt>0.1) 
	if(fMassCalc4->Pt2() < minPt*minPt) break;
      // mass cut
      pdg4[0]=321; pdg4[1]=211; pdg4[2]=211; pdg4[3]=211;
      mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();
      minv = fMassCalc4->InvMass(nprongs,pdg4);
      if(TMath::Abs(minv-mPDG)<fCutsD0toKpipipi->GetMassCut()) retval=kTRUE;
      pdg4[0]=211; pdg4[1]=321; pdg4[2]=211; pdg4[3]=211;
      mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();
      minv = fMassCalc4->InvMass(nprongs,pdg4);
      if(TMath::Abs(minv-mPDG)<fCutsD0toKpipipi->GetMassCut()) retval=kTRUE;
      pdg4[0]=211; pdg4[1]=211; pdg4[2]=321; pdg4[3]=211;
      mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();
      minv = fMassCalc4->InvMass(nprongs,pdg4);
      if(TMath::Abs(minv-mPDG)<fCutsD0toKpipipi->GetMassCut()) retval=kTRUE;
      pdg4[0]=211; pdg4[1]=211; pdg4[2]=211; pdg4[3]=321;
      mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();
      minv = fMassCalc4->InvMass(nprongs,pdg4);
      if(TMath::Abs(minv-mPDG)<fCutsD0toKpipipi->GetMassCut()) retval=kTRUE;
      break;
    case 5:
      fMassCalc2->SetPxPyPzProngs(nprongs,px,py,pz);
      mPDG=TDatabasePDG::Instance()->GetParticle(4122)->Mass();
      pdg2[0]=2212;pdg2[1]=310;
      minv=fMassCalc2->InvMass(2,pdg2);
      if(TMath::Abs(minv-mPDG)<fCutsLctoV0->GetMassCut()) retval=kTRUE;
      pdg2[0]=211;pdg2[1]=3122;
      minv=fMassCalc2->InvMass(2,pdg2);
      if(TMath::Abs(minv-mPDG)<fCutsLctoV0->GetMassCut()) retval=kTRUE;
      break;
    default:
      printf("SelectInvMassAndPt(): wrong decay selection\n");
      break;
    }

  return retval;
}
//-----------------------------------------------------------------------------
void AliAnalysisVertexingHF::SelectTracksAndCopyVertex(const AliVEvent *event,
						       Int_t trkEntries,
						       TObjArray &seleTrksArray,TObjArray &tracksAtVertex,
						       Int_t &nSeleTrks,
						       UChar_t *seleFlags,Int_t *evtNumber)
{
  // Apply single-track preselection.
  // Fill a TObjArray with selected tracks (for displaced vertices or
  // soft pion from D*). Selection flag stored in seleFlags.
  // Create the AliESDVertex object (convert from AliAODVertex if necessary)
  // In case of AOD input, also fill fAODMap for track index<->ID

  const AliVVertex *vprimary = event->GetPrimaryVertex();

  if(fV1) { delete fV1; fV1=NULL; }
  if(fAODMap) { delete fAODMap; fAODMap=NULL; }

  Int_t nindices=0;
  UShort_t *indices = 0;
  Double_t pos[3],cov[6];

  if(!fInputAOD) { // ESD
    fV1 = new AliESDVertex(*((AliESDVertex*)vprimary));
  } else {         // AOD
    vprimary->GetXYZ(pos);
    vprimary->GetCovarianceMatrix(cov);
    fV1 = new AliESDVertex(pos,cov,100.,100,vprimary->GetName());
    indices = new UShort_t[event->GetNumberOfTracks()];
    for(Int_t ijk=0; ijk<event->GetNumberOfTracks(); ijk++) indices[ijk]=0;
    fAODMapSize = 100000;
    fAODMap = new Int_t[fAODMapSize];
  }


  Int_t entries = (Int_t)event->GetNumberOfTracks();
  Bool_t okDisplaced=kFALSE,okSoftPi=kFALSE;
  nSeleTrks=0;
 
  // transfer ITS tracks from event to arrays
  for(Int_t i=0; i<entries; i++) {
    AliVTrack *track;
    track = (AliVTrack*)event->GetTrack(i);

    // skip pure ITS SA tracks
    if(track->GetStatus()&AliESDtrack::kITSpureSA) continue;

    // skip tracks without ITS
    if(!(track->GetStatus()&AliESDtrack::kITSin)) continue;

    // TEMPORARY: check that the cov matrix is there
    Double_t covtest[21];
    if(!track->GetCovarianceXYZPxPyPz(covtest)) continue;
    //

    if(fInputAOD) {
      AliAODTrack *aodt = (AliAODTrack*)track;
      if(aodt->GetUsedForPrimVtxFit()) { 
	indices[nindices]=aodt->GetID(); nindices++; 
      }
      fAODMap[(Int_t)aodt->GetID()] = i;
    }

    AliESDtrack *esdt = 0;

    if(!fInputAOD) {
      esdt = (AliESDtrack*)track;
    } else {
      esdt = new AliESDtrack(track);
    }

    // single track selection
    okDisplaced=kFALSE; okSoftPi=kFALSE;
    if(fMixEvent && i<trkEntries){
      evtNumber[i]=((AliMixedEvent*)event)->EventIndex(i);
      const AliVVertex* eventVtx=((AliMixedEvent*)event)->GetEventVertex(i);
      Double_t vtxPos[3],primPos[3],primCov[6],trasl[3];
      eventVtx->GetXYZ(vtxPos);
      vprimary->GetXYZ(primPos);
      eventVtx->GetCovarianceMatrix(primCov);
      for(Int_t ind=0;ind<3;ind++){
	trasl[ind]=vtxPos[ind]-primPos[ind];
      }
      
      Bool_t isTransl=esdt->Translate(trasl,primCov);
      if(!isTransl) {
	delete esdt;
	esdt = NULL;
	continue;
      }
    }

    if(SingleTrkCuts(esdt,okDisplaced,okSoftPi) && nSeleTrks<trkEntries) {
      esdt->PropagateToDCA(fV1,fBzkG,kVeryBig);
      seleTrksArray.AddLast(esdt);
      tracksAtVertex.AddLast(new AliExternalTrackParam(*esdt));

      seleFlags[nSeleTrks]=0;
      if(okDisplaced) SETBIT(seleFlags[nSeleTrks],kBitDispl);
      if(okSoftPi)    SETBIT(seleFlags[nSeleTrks],kBitSoftPi);
      nSeleTrks++;
    } else {
      if(fInputAOD) delete esdt; 
      esdt = NULL;
      continue;
    } 

  } // end loop on tracks

  // primary vertex from AOD
  if(fInputAOD) {
    delete fV1; fV1=NULL;
    vprimary->GetXYZ(pos);
    vprimary->GetCovarianceMatrix(cov);
    Double_t chi2toNDF = vprimary->GetChi2perNDF();
    Int_t ncontr=nindices;
    if(!strcmp(vprimary->GetTitle(),"VertexerTracksWithContraint")) ncontr += 1;
    Double_t chi2=chi2toNDF*(2.*(Double_t)ncontr-3.); 
    fV1 = new AliESDVertex(pos,cov,chi2,ncontr,vprimary->GetName());
    fV1->SetTitle(vprimary->GetTitle());
    fV1->SetIndices(nindices,indices);
  }
  if(indices) { delete [] indices; indices=NULL; }


  return;
}
//-----------------------------------------------------------------------------
void AliAnalysisVertexingHF::SetSelectionBitForPID(AliRDHFCuts *cuts,AliAODRecoDecayHF *rd,Int_t bit) {
  //
  // Set the selection bit for PID
  //
  if(cuts->GetPidHF()) {
    Bool_t usepid=cuts->GetIsUsePID();
    cuts->SetUsePID(kTRUE);
    if(cuts->IsSelectedPID(rd))
      rd->SetSelectionBit(bit);
    cuts->SetUsePID(usepid);
  }
  return;
}
//-----------------------------------------------------------------------------
Bool_t AliAnalysisVertexingHF::SingleTrkCuts(AliESDtrack *trk,
					     Bool_t &okDisplaced,Bool_t &okSoftPi) const 
{
  // Check if track passes some kinematical cuts  

  // this is needed to store the impact parameters
  trk->RelateToVertex(fV1,fBzkG,kVeryBig);

  UInt_t selectInfo;
  //
  // Track selection, displaced tracks
  selectInfo = 0; 
  if(fTrackFilter) {
    selectInfo = fTrackFilter->IsSelected(trk);
  }
  if(selectInfo) okDisplaced=kTRUE;

  // Track selection, soft pions
  selectInfo = 0; 
  if(fDstar && fTrackFilterSoftPi) {
    selectInfo = fTrackFilterSoftPi->IsSelected(trk);
  }
  if(selectInfo) okSoftPi=kTRUE;

  if(okDisplaced || okSoftPi) return kTRUE;

  return kFALSE;
}


//-----------------------------------------------------------------------------
AliAODv0* AliAnalysisVertexingHF::TransformESDv0toAODv0(AliESDv0 *esdV0, TObjArray *twoTrackArrayV0){
  //
  // Transform ESDv0 to AODv0
  //
  //  this function takes the ESDv0 vertex, computes the DCA variables from the ESDv0
  //  and creates an AODv0 out of them
  //
  Double_t vertex[3]; esdV0->GetXYZ(vertex[0],vertex[1],vertex[2]);
  AliAODVertex *vertexV0 = new AliAODVertex(vertex,esdV0->GetChi2V0(),AliAODVertex::kV0,2);

  // create the v0 neutral track to compute the DCA to the primary vertex
  Double_t xyz[3], pxpypz[3];
  esdV0->XvYvZv(xyz);
  esdV0->PxPyPz(pxpypz);
  Double_t cv[21]; for(int i=0; i<21; i++) cv[i]=0;
  AliNeutralTrackParam *trackesdV0 = new AliNeutralTrackParam(xyz,pxpypz,cv,0);
  if(!trackesdV0) {
    delete vertexV0;
    return 0;
  }
  Double_t d0z0[2],covd0z0[3];
  AliAODVertex *primVertexAOD = PrimaryVertex(); 
  trackesdV0->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
  Double_t dcaV0ToPrimVertex = TMath::Sqrt(covd0z0[0]);
  // get the v0 daughters to compute their DCA to the v0 vertex and get their momentum
  Double_t dcaV0DaughterToPrimVertex[2];  
  AliExternalTrackParam *posV0track = (AliExternalTrackParam*)twoTrackArrayV0->UncheckedAt(0);
  AliExternalTrackParam *negV0track = (AliExternalTrackParam*)twoTrackArrayV0->UncheckedAt(1);
  if( !posV0track || !negV0track) {
    if(trackesdV0) {delete trackesdV0; trackesdV0=NULL;}
    delete vertexV0;
    delete primVertexAOD;
    return 0;
  }
  posV0track->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);
  //  if ( covd0z0[0]<=0.) dcaV0DaughterToPrimVertex[0] = 0;
  //  else 
  dcaV0DaughterToPrimVertex[0] = TMath::Sqrt(covd0z0[0]);
  negV0track->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0,covd0z0);  
  //  if ( covd0z0[0]<=0.)dcaV0DaughterToPrimVertex[1] = 0;
  //  else 
  dcaV0DaughterToPrimVertex[1] = TMath::Sqrt(covd0z0[0]);
  Double_t dcaV0Daughters = esdV0->GetDcaV0Daughters();
  Double_t pmom[3],nmom[3];
  esdV0->GetNPxPyPz(nmom[0],nmom[1],nmom[2]);
  esdV0->GetPPxPyPz(pmom[0],pmom[1],pmom[2]);

  AliAODv0 *aodV0 = new AliAODv0(vertexV0,dcaV0Daughters,dcaV0ToPrimVertex,pmom,nmom,dcaV0DaughterToPrimVertex);
  aodV0->SetOnFlyStatus(esdV0->GetOnFlyStatus());

  delete trackesdV0;
  delete primVertexAOD;

  return aodV0;
}
//-----------------------------------------------------------------------------
void AliAnalysisVertexingHF::SetParametersAtVertex(AliESDtrack* esdt, AliExternalTrackParam* extpar) const{
  // Set the stored track parameters at primary vertex into AliESDtrack
  Double_t xyz[3], pxpypz[3], cv[21]; 
  extpar->PxPyPz(pxpypz); 	                  
  extpar->XvYvZv(xyz);
  extpar->GetCovarianceXYZPxPyPz(cv);
  Short_t sign=extpar->Charge();
  esdt->Set(xyz,pxpypz,cv,sign);
  return;
}
