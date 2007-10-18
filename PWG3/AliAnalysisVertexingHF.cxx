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

//----------------------------------------------------------------------------
//    Implementation of the heavy-flavour vertexing analysis class
// Candidates are store as objects deriving from AliAODRecoDecay.
// An example of usage can be found in the macro AliVertexingHFTest.C.
// Can be used as a task of AliAnalysisManager by means of the interface
// class AliAnalysisTaskVertexingHF. 
//
//  Origin: E.Bruna, G.E.Bruno, A.Dainese, F.Prino, R.Romita
//----------------------------------------------------------------------------
#include <TTree.h>
#include <TDatabasePDG.h>
#include <TString.h>
#include "AliESDEvent.h"
#include "AliVertexerTracks.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoDecayHF4Prong.h"
#include "AliAnalysisVertexingHF.h"

ClassImp(AliAnalysisVertexingHF)

//----------------------------------------------------------------------------
AliAnalysisVertexingHF::AliAnalysisVertexingHF():
fRecoPrimVtxSkippingTrks(kFALSE),
fRmTrksFromPrimVtx(kFALSE),
fV1(0x0),
fDebug(0),
fD0toKpi(kTRUE),fJPSItoEle(kTRUE),f3Prong(kTRUE),f4Prong(kTRUE),
fITSrefit(kFALSE),fBothSPD(kTRUE),fMinITSCls(5),
fMinPtCut(0.),fMind0rphiCut(0.)
{
  // Default constructor

  SetD0toKpiCuts();
  SetBtoJPSICuts();
  SetDplusCuts();
}
//--------------------------------------------------------------------------
AliAnalysisVertexingHF::AliAnalysisVertexingHF(const AliAnalysisVertexingHF &source) : 
TNamed(source),
fRecoPrimVtxSkippingTrks(source.fRecoPrimVtxSkippingTrks),
fRmTrksFromPrimVtx(source.fRmTrksFromPrimVtx),
fV1(source.fV1),
fDebug(source.fDebug),
fD0toKpi(source.fD0toKpi),
fJPSItoEle(source.fJPSItoEle),
f3Prong(source.f3Prong),
f4Prong(source.f4Prong),
fITSrefit(source.fITSrefit),
fBothSPD(source.fBothSPD),
fMinITSCls(source.fMinITSCls),
fMinPtCut(source.fMinPtCut),
fMind0rphiCut(source.fMind0rphiCut)
{
  //
  // Copy constructor
  //
  for(Int_t i=0; i<9; i++)  fD0toKpiCuts[i]=source.fD0toKpiCuts[i];
  for(Int_t i=0; i<9; i++)  fBtoJPSICuts[i]=source.fBtoJPSICuts[i];
  for(Int_t i=0; i<12; i++) fDplusCuts[i]=source.fDplusCuts[i];
}
//--------------------------------------------------------------------------
AliAnalysisVertexingHF &AliAnalysisVertexingHF::operator=(const AliAnalysisVertexingHF &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;
  fRecoPrimVtxSkippingTrks = source.fRecoPrimVtxSkippingTrks;
  fRmTrksFromPrimVtx = source.fRmTrksFromPrimVtx;
  fV1 = source.fV1;
  fDebug = source.fDebug;
  fD0toKpi = source.fD0toKpi;
  fJPSItoEle = source.fJPSItoEle;
  f3Prong = source.f3Prong;
  f4Prong = source.f4Prong;
  fITSrefit = source.fITSrefit;
  fBothSPD = source.fBothSPD;
  fMinITSCls = source.fMinITSCls;
  fMinPtCut = source.fMinPtCut;
  fMind0rphiCut = source.fMind0rphiCut;

  for(Int_t i=0; i<9; i++)  fD0toKpiCuts[i]=source.fD0toKpiCuts[i];
  for(Int_t i=0; i<9; i++)  fBtoJPSICuts[i]=source.fBtoJPSICuts[i];
  for(Int_t i=0; i<12; i++) fDplusCuts[i]=source.fDplusCuts[i];

  return *this;
}
//----------------------------------------------------------------------------
AliAnalysisVertexingHF::~AliAnalysisVertexingHF() {
  // Destructor
  if(fV1) { delete fV1; fV1=0; }
}
//----------------------------------------------------------------------------
void AliAnalysisVertexingHF::FindCandidates(AliESDEvent *esd,TTree treeout[])
{
  // Find heavy-flavour vertex candidates

  
  AliAODRecoDecayHF2Prong *io2Prong = new AliAODRecoDecayHF2Prong;
  AliAODRecoDecayHF3Prong *io3Prong = new AliAODRecoDecayHF3Prong;
  AliAODRecoDecayHF4Prong *io4Prong = new AliAODRecoDecayHF4Prong;
  Int_t itree=0;
  Int_t itreeD0toKpi=-1,itreeJPSItoEle=-1,itree3Prong=-1,itree4Prong=-1;
  Int_t initEntriesD0toKpi=0,initEntriesJPSItoEle=0,initEntries3Prong=0,initEntries4Prong=0;
  if(fD0toKpi) {
    itreeD0toKpi=itree;
    //treeout[itree].Print();
    treeout[itree].SetBranchAddress("D0toKpi",&io2Prong);
    itree++;
    initEntriesD0toKpi = treeout[itreeD0toKpi].GetEntries();
  }
  if(fJPSItoEle) {
    itreeJPSItoEle=itree;
    //treeout[itree].Print();
    treeout[itree].SetBranchAddress("JPSItoEle",&io2Prong);
    itree++;
    initEntriesJPSItoEle = treeout[itreeJPSItoEle].GetEntries();
  }
  if(f3Prong) {
    itree3Prong=itree;
    treeout[itree].SetBranchAddress("Charmto3Prong",&io3Prong);
    itree++;
    initEntries3Prong = treeout[itree3Prong].GetEntries();
  }
  if(f4Prong) {
    itree4Prong=itree;
    treeout[itree].SetBranchAddress("D0to4Prong",&io4Prong);
    itree++;
    initEntries4Prong = treeout[itree4Prong].GetEntries();
  }
  delete io2Prong; io2Prong = NULL;
  delete io3Prong; io3Prong = NULL;
  delete io4Prong; io4Prong = NULL;

  Int_t    iTrkP1,iTrkP2,iTrkN1,iTrkN2,trkEntries;
  Int_t    nTrksP=0,nTrksN=0;
  Double_t xdummy,ydummy,dcap1n1,dcap1n2,dcap2n1,dcap1p2,dcan1n2;
  Bool_t   okD0=kFALSE,okJPSI=kFALSE,ok3Prong=kFALSE,ok4Prong=kFALSE;
  AliESDtrack *postrack1 = 0;
  AliESDtrack *postrack2 = 0;
  AliESDtrack *negtrack1 = 0;
  AliESDtrack *negtrack2 = 0;
  Double_t dcaMax = fD0toKpiCuts[1];
  if(dcaMax<fBtoJPSICuts[1]) dcaMax=fBtoJPSICuts[1];
  if(dcaMax<fDplusCuts[11]) dcaMax=fDplusCuts[11];
  if(fDebug) printf(" dca cut set to %f cm\n",dcaMax);

  // create the AliVertexerTracks object for the secondary vertex
  AliVertexerTracks *vertexer2 = new AliVertexerTracks(esd->GetMagneticField());
  vertexer2->SetDebug(0);

  Int_t ev = (Int_t)esd->GetEventNumberInFile();
  printf("--- Finding candidates in event %d\n",ev);

  Double_t bfieldkG=(Double_t)esd->GetMagneticField(); 

  trkEntries = (Int_t)esd->GetNumberOfTracks();
  printf(" Number of tracks: %d\n",trkEntries);
  if(trkEntries<2) return;

  // retrieve primary vertex from the AliESD
  if(!esd->GetPrimaryVertex()) { 
    printf(" No vertex in AliESD\n");
    return;
  }
  AliESDVertex copy(*(esd->GetPrimaryVertex()));
  SetPrimaryVertex(&copy);
  vertexer2->SetVtxStart(&copy);

  // call function which applies sigle-track selection and
  // separetes positives and negatives
  TObjArray trksP(trkEntries/2); // will become TClonesArray
  TObjArray trksN(trkEntries/2); // will become TClonesArray
  SelectTracks(esd,trksP,nTrksP,
	           trksN,nTrksN);

  printf(" Pos. tracks: %d    Neg. tracks: %d\n",nTrksP,nTrksN);

  TObjArray *twoTrackArray1 = new TObjArray(2);
  TObjArray *twoTrackArray2 = new TObjArray(2);
  TObjArray *threeTrackArray = new TObjArray(3);
  TObjArray *fourTrackArray = new TObjArray(4);

  // LOOP ON  POSITIVE  TRACKS
  for(iTrkP1=0; iTrkP1<nTrksP; iTrkP1++) {
    if(fDebug) if(iTrkP1%1==0) printf("  Processing positive track number %d of %d\n",iTrkP1,nTrksP);  
    // get track from track array
    postrack1 = (AliESDtrack*)trksP.UncheckedAt(iTrkP1);

    // LOOP ON  NEGATIVE  TRACKS
    for(iTrkN1=0; iTrkN1<nTrksN; iTrkN1++) {
      if(fDebug) if(iTrkN1%1==0) printf("    Processing negative track number %d of %d\n",iTrkN1,nTrksN);  
      // get track from tracks array
      negtrack1 = (AliESDtrack*)trksN.UncheckedAt(iTrkN1);
      // DCA between the two tracks
      dcap1n1 = postrack1->GetDCA(negtrack1,bfieldkG,xdummy,ydummy);
      if(dcap1n1>dcaMax) { negtrack1=0; continue; }
      // Vertexing with AliVertexerTracks
      twoTrackArray1->AddAt(postrack1,0);
      twoTrackArray1->AddAt(negtrack1,1);
      AliESDVertex *vertexp1n1 = vertexer2->VertexForSelectedTracks(twoTrackArray1);
      if(vertexp1n1->GetNContributors()!=2) { 
	if(fDebug) printf("two-track vertexing failed\n"); 
	negtrack1=0; continue; 
      }
      if(fD0toKpi || fJPSItoEle) { 
	io2Prong = Make2Prong(twoTrackArray1,esd,vertexp1n1,dcap1n1,okD0,okJPSI);
	if(okD0) treeout[itreeD0toKpi].Fill();
	if(okJPSI) treeout[itreeJPSItoEle].Fill();
        delete io2Prong; io2Prong=NULL; 
      }

      twoTrackArray1->Clear(); 
      if(!f3Prong && !f4Prong)  { 
	negtrack1=0; 
	delete vertexp1n1; 
	continue; 
      }
      
      // 2nd LOOP  ON  POSITIVE  TRACKS 
      for(iTrkP2=iTrkP1+1; iTrkP2<nTrksP; iTrkP2++) {
	// get track from tracks array
	postrack2 = (AliESDtrack*)trksP.UncheckedAt(iTrkP2);
	dcap2n1 = postrack2->GetDCA(negtrack1,bfieldkG,xdummy,ydummy);
	if(dcap2n1>dcaMax) { postrack2=0; continue; }
	dcap1p2 = postrack2->GetDCA(postrack1,bfieldkG,xdummy,ydummy);
	if(dcap1p2>dcaMax) { postrack2=0; continue; }

	// Vertexing with AliVertexerTracks
	twoTrackArray2->AddAt(postrack2,0);
	twoTrackArray2->AddAt(negtrack1,1);
	AliESDVertex *vertexp2n1 = vertexer2->VertexForSelectedTracks(twoTrackArray2);
	if(f3Prong) { 
	  threeTrackArray->AddAt(postrack1,0);
	  threeTrackArray->AddAt(negtrack1,1);
	  threeTrackArray->AddAt(postrack2,2);
	  io3Prong = Make3Prong(threeTrackArray,esd,vertexp1n1,vertexp2n1,dcap1n1,dcap2n1,dcap1p2,ok3Prong);
	  if(ok3Prong) treeout[itree3Prong].Fill();
	  if(io3Prong) delete io3Prong; io3Prong=NULL; 
	}
	if(f4Prong) {
	  // 3rd LOOP  ON  NEGATIVE  TRACKS (for 4 prong) 
	  for(iTrkN2=iTrkN1+1; iTrkN2<nTrksN; iTrkN2++) {
	    // get track from tracks array
	    negtrack2 = (AliESDtrack*)trksN.UncheckedAt(iTrkN2);
	    dcap1n2 = postrack1->GetDCA(negtrack2,bfieldkG,xdummy,ydummy);
	    if(dcap1n2>dcaMax) { negtrack2=0; continue; }
	    // Vertexing with AliVertexerTracks
	    fourTrackArray->AddAt(postrack1,0);
	    fourTrackArray->AddAt(negtrack1,1);
	    fourTrackArray->AddAt(postrack2,2);
	    fourTrackArray->AddAt(negtrack2,3);
	    io4Prong = Make4Prong(fourTrackArray,esd,vertexp1n1,vertexp2n1,dcap1n1,dcap1n2,dcap2n1,ok4Prong);
	    if(ok4Prong) treeout[itree4Prong].Fill();
	    delete io4Prong; io4Prong=NULL; 
            fourTrackArray->Clear();
	    negtrack2 = 0;
	  } // end loop on negative tracks
	}
	postrack2 = 0;
	delete vertexp2n1;
      } // end 2nd loop on positive tracks
      twoTrackArray2->Clear();
      
      // 2nd LOOP  ON  NEGATIVE  TRACKS 
      for(iTrkN2=iTrkN1+1; iTrkN2<nTrksN; iTrkN2++) {
	// get track from tracks array
	negtrack2 = (AliESDtrack*)trksN.UncheckedAt(iTrkN2);
	dcap1n2 = postrack1->GetDCA(negtrack2,bfieldkG,xdummy,ydummy);
	if(dcap1n2>dcaMax) { negtrack2=0; continue; }
	dcan1n2 = negtrack1->GetDCA(negtrack2,bfieldkG,xdummy,ydummy);
	if(dcan1n2>dcaMax) { negtrack2=0; continue; }

	// Vertexing with AliVertexerTracks
	twoTrackArray2->AddAt(postrack1,0);
	twoTrackArray2->AddAt(negtrack2,1);
	AliESDVertex *vertexp1n2 = vertexer2->VertexForSelectedTracks(twoTrackArray2);
	if(f3Prong) { 
	  threeTrackArray->AddAt(negtrack1,0);
	  threeTrackArray->AddAt(postrack1,1);
	  threeTrackArray->AddAt(negtrack2,2);
	  io3Prong = Make3Prong(threeTrackArray,esd,vertexp1n1,vertexp1n2,dcap1n1,dcap1n2,dcan1n2,ok3Prong);
	  if(ok3Prong) treeout[itree3Prong].Fill();
	  if(io3Prong) delete io3Prong; io3Prong=NULL; 
	  }
	negtrack2 = 0;
	delete vertexp1n2;
      } // end 2nd loop on negative tracks
      twoTrackArray2->Clear();

      negtrack1 = 0;
      delete vertexp1n1; 
    } // end 1st loop on negative tracks

    postrack1 = 0;
  }  // end 1st loop on positive tracks
    


  delete vertexer2;
  //printf("delete twoTr 1\n");
  twoTrackArray1->Delete(); delete twoTrackArray1;
  //printf("delete twoTr 2\n");
  twoTrackArray2->Delete(); delete twoTrackArray2;
  //printf("delete threeTr 1\n");
  threeTrackArray->Clear(); 
  threeTrackArray->Delete(); delete threeTrackArray;
  //printf("delete fourTr 1\n");
  fourTrackArray->Delete(); delete fourTrackArray;


  // create a copy of this class to be written to output file
  //AliAnalysisVertexingHF *copy = (AliAnalysisVertexingHF*)this->Clone("AnalysisVertexingHF");

  // print statistics
  if(fD0toKpi) {
    printf(" D0->Kpi: event %d = %d; total = %d;\n",
	   (Int_t)esd->GetEventNumberInFile(),
	   (Int_t)treeout[itreeD0toKpi].GetEntries()-initEntriesD0toKpi,
	   (Int_t)treeout[itreeD0toKpi].GetEntries());
  }
  if(fJPSItoEle) {
    printf(" JPSI->ee: event %d = %d; total = %d;\n",
	   (Int_t)esd->GetEventNumberInFile(),
	   (Int_t)treeout[itreeJPSItoEle].GetEntries()-initEntriesJPSItoEle,
	   (Int_t)treeout[itreeJPSItoEle].GetEntries());
  }
  if(f3Prong) {
    printf(" Charm->3Prong: event %d = %d; total = %d;\n",
	   (Int_t)esd->GetEventNumberInFile(),
	   (Int_t)treeout[itree3Prong].GetEntries()-initEntries3Prong,
	   (Int_t)treeout[itree3Prong].GetEntries());
  }
  if(f4Prong) {
    printf(" Charm->4Prong: event %d = %d; total = %d;\n",
	   (Int_t)esd->GetEventNumberInFile(),
	   (Int_t)treeout[itree4Prong].GetEntries()-initEntries4Prong,
	   (Int_t)treeout[itree4Prong].GetEntries());
  }


  return;
}
//----------------------------------------------------------------------------
AliAODRecoDecayHF2Prong *AliAnalysisVertexingHF::Make2Prong(
				   TObjArray *twoTrackArray1,AliESDEvent *esd,
				   AliESDVertex *secVertexESD,Double_t dca,
				   Bool_t &okD0,Bool_t &okJPSI) const
{
  // Make 2Prong candidates and check if they pass D0toKpi or BtoJPSI
  // reconstruction cuts
  // G.E.Bruno (J/psi), A.Dainese (D0->Kpi)

  okD0=kFALSE; okJPSI=kFALSE;

  Double_t px[2],py[2],pz[2],d0[2],d0err[2];
  Double_t bfieldkG=(Double_t)esd->GetMagneticField(); 

  AliESDtrack *postrack = (AliESDtrack*)twoTrackArray1->UncheckedAt(0);
  AliESDtrack *negtrack = (AliESDtrack*)twoTrackArray1->UncheckedAt(1);

  // propagate tracks to secondary vertex, to compute inv. mass
  postrack->RelateToVertex(secVertexESD,bfieldkG,10.);
  negtrack->RelateToVertex(secVertexESD,bfieldkG,10.);

  Double_t momentum[3];
  postrack->GetPxPyPz(momentum);
  px[0] = momentum[0]; py[0] = momentum[1]; pz[0] = momentum[2]; 
  negtrack->GetPxPyPz(momentum);
  px[1] = momentum[0]; py[1] = momentum[1]; pz[1] = momentum[2]; 


  // invariant mass cut (try to improve coding here..)
  Bool_t okMassCut=kFALSE;
  if(!okMassCut && fD0toKpi) if(SelectInvMass(0,2,px,py,pz)) okMassCut=kTRUE;
  if(!okMassCut && fJPSItoEle) if(SelectInvMass(1,2,px,py,pz)) okMassCut=kTRUE;
  if(!okMassCut) {
    if(fDebug) printf(" candidate didn't pass mass cut\n");
    return 0x0;    
  }


  AliESDVertex *primVertex = fV1;  
  AliESDVertex *ownPrimVertex=0;

  // primary vertex from *other* tracks in the event
  if(fRecoPrimVtxSkippingTrks || fRmTrksFromPrimVtx) {
    ownPrimVertex = OwnPrimaryVertex(2,twoTrackArray1,esd);
    if(!ownPrimVertex) {
      return 0x0;
    } else {
      if(ownPrimVertex->GetNContributors()<2) {
	delete ownPrimVertex;
	return 0x0;
      } else {
	primVertex = ownPrimVertex;
      }
    }
  }

  Float_t d0z0[2],covd0z0[3];
  postrack->RelateToVertex(primVertex,bfieldkG,10.);
  postrack->GetImpactParameters(d0z0,covd0z0);
  d0[0] = d0z0[0];
  d0err[0] = TMath::Sqrt(covd0z0[0]);
  negtrack->RelateToVertex(primVertex,bfieldkG,10.);
  negtrack->GetImpactParameters(d0z0,covd0z0);
  d0[1] = d0z0[0];
  d0err[1] = TMath::Sqrt(covd0z0[0]);

  // create the object AliAODRecoDecayHF2Prong
  Double_t pos[3],cov[6];
  secVertexESD->GetXYZ(pos); // position
  secVertexESD->GetCovMatrix(cov); //covariance matrix
  AliAODVertex *secVertexAOD = new AliAODVertex(pos,cov,secVertexESD->GetChi2toNDF());
  primVertex->GetXYZ(pos); // position
  primVertex->GetCovMatrix(cov); //covariance matrix
  AliAODVertex *primVertexAOD = new AliAODVertex(pos,cov,primVertex->GetChi2toNDF());
  AliAODRecoDecayHF2Prong *the2Prong = new AliAODRecoDecayHF2Prong(secVertexAOD,px,py,pz,d0,d0err,dca);
  the2Prong->SetOwnPrimaryVtx(primVertexAOD);

  // select D0->Kpi
  Int_t checkD0,checkD0bar;
  if(fD0toKpi) okD0 = the2Prong->SelectD0(fD0toKpiCuts,checkD0,checkD0bar);
  //if(fDebug && fD0toKpi) printf("   %d\n",(Int_t)okD0);
  // select J/psi from B
  Int_t checkJPSI;
  if(fJPSItoEle) okJPSI = the2Prong->SelectBtoJPSI(fBtoJPSICuts,checkJPSI);
  //if(fDebug && fJPSItoEle) printf("   %d\n",(Int_t)okJPSI);


  if(okD0 || okJPSI) {
    // get PID info from ESD
    Double_t esdpid0[5];
    postrack->GetESDpid(esdpid0);
    Double_t esdpid1[5];
    negtrack->GetESDpid(esdpid1);
    Double_t esdpid[10];
    for(Int_t i=0;i<5;i++) {
      esdpid[i]   = esdpid0[i];
      esdpid[5+i] = esdpid1[i];
    }
    the2Prong->SetPID(2,esdpid);
  }

  if(ownPrimVertex) delete ownPrimVertex;	

  return the2Prong;  
}
//----------------------------------------------------------------------------
AliAODRecoDecayHF3Prong* AliAnalysisVertexingHF::Make3Prong(
                             TObjArray *threeTrackArray,AliESDEvent *esd,
			     AliESDVertex *vertexp1n1,AliESDVertex *vertexp2n1,
			     Double_t dcap1n1,Double_t dcap2n1,Double_t dcap1p2,
			     Bool_t &ok3Prong) const
{
  // Make 3Prong candidates and check if they pass Dplus or Ds or Lambdac
  // reconstruction cuts 
  // E.Bruna, F.Prino

  ok3Prong=kFALSE;
  Double_t px[3],py[3],pz[3],d0[3],d0err[3];//d0z[3];  
  Float_t d0z0[2],covd0z0[3];

  Double_t bfieldkG=(Double_t)esd->GetMagneticField(); 

  AliESDtrack *postrack1 = (AliESDtrack*)threeTrackArray->UncheckedAt(0);
  AliESDtrack *negtrack = (AliESDtrack*)threeTrackArray->UncheckedAt(1);
  AliESDtrack *postrack2 = (AliESDtrack*)threeTrackArray->UncheckedAt(2);

  AliESDVertex *primVertex = fV1;  

  postrack1->RelateToVertex(primVertex,bfieldkG,10.);
  negtrack->RelateToVertex(primVertex,bfieldkG,10.);
  postrack2->RelateToVertex(primVertex,bfieldkG,10.);

  Double_t momentum[3];
  postrack1->GetPxPyPz(momentum);
  px[0] = momentum[0]; py[0] = momentum[1]; pz[0] = momentum[2]; 
  negtrack->GetPxPyPz(momentum);
  px[1] = momentum[0]; py[1] = momentum[1]; pz[1] = momentum[2]; 
  postrack2->GetPxPyPz(momentum);
  px[2] = momentum[0]; py[2] = momentum[1]; pz[2] = momentum[2]; 

  postrack1->GetImpactParameters(d0z0,covd0z0);
  d0[0]=d0z0[0];
  d0err[0] = TMath::Sqrt(covd0z0[0]);
  negtrack->GetImpactParameters(d0z0,covd0z0);
  d0[1]=d0z0[0];
  d0err[1] = TMath::Sqrt(covd0z0[0]);
  postrack2->GetImpactParameters(d0z0,covd0z0);
  d0[2]=d0z0[0];
  d0err[2] = TMath::Sqrt(covd0z0[0]);


  // invariant mass cut for D+ (try to improve coding here..)
  Bool_t okMassCut=kFALSE;
  if(!okMassCut && f3Prong) if(SelectInvMass(2,3,px,py,pz)) okMassCut=kTRUE;
  if(!okMassCut) {
    if(fDebug) printf(" candidate didn't pass mass cut\n");
    return 0x0;    
  }

  //charge
  Short_t charge=(Short_t)(postrack1->GetSign()*postrack2->GetSign()*negtrack->GetSign());

  AliESDVertex *ownPrimVertex = 0;  
  // primary vertex from *other* tracks in the event
  if(fRecoPrimVtxSkippingTrks || fRmTrksFromPrimVtx) {
    ownPrimVertex = OwnPrimaryVertex(3,threeTrackArray,esd);
    if(!ownPrimVertex) {
      return 0x0;
    } else {
      if(ownPrimVertex->GetNContributors()<2) {
	delete ownPrimVertex;
	return 0x0;
      } else {
	primVertex = ownPrimVertex;
      }
    }
  }

  // create the object AliAODRecoDecayHF3Prong
  AliVertexerTracks *vertexer = new AliVertexerTracks(esd->GetMagneticField());
  AliESDVertex* secVert3Prong=(AliESDVertex*)vertexer->VertexForSelectedTracks(threeTrackArray,kTRUE,kTRUE);
  delete vertexer; vertexer=NULL;
  Double_t pos[3],cov[6],sigmavert;
  secVert3Prong->GetXYZ(pos); // position
  secVert3Prong->GetCovMatrix(cov); //covariance matrix
  sigmavert=secVert3Prong->GetDispersion();

  AliAODVertex *secVert3PrAOD = new AliAODVertex(pos,cov,secVert3Prong->GetChi2toNDF());
  primVertex->GetXYZ(pos); // position
  primVertex->GetCovMatrix(cov); //covariance matrix
  AliAODVertex *primVertexAOD = new AliAODVertex(pos,cov,primVertex->GetChi2toNDF());
  Double_t dca[3]={dcap1n1,dcap2n1,dcap1p2};

  Double_t dist12=TMath::Sqrt((vertexp1n1->GetXv()-pos[0])*(vertexp1n1->GetXv()-pos[0])+(vertexp1n1->GetYv()-pos[1])*(vertexp1n1->GetYv()-pos[1])+(vertexp1n1->GetZv()-pos[2])*(vertexp1n1->GetZv()-pos[2]));
  Double_t dist23=TMath::Sqrt((vertexp2n1->GetXv()-pos[0])*(vertexp2n1->GetXv()-pos[0])+(vertexp2n1->GetYv()-pos[1])*(vertexp2n1->GetYv()-pos[1])+(vertexp2n1->GetZv()-pos[2])*(vertexp2n1->GetZv()-pos[2]));

  AliAODRecoDecayHF3Prong *the3Prong = new AliAODRecoDecayHF3Prong(secVert3PrAOD,px,py,pz,d0,d0err,dca,sigmavert,dist12,dist23,charge);
  the3Prong->SetOwnPrimaryVtx(primVertexAOD);


  // select D+->Kpipi
  if(f3Prong) ok3Prong = the3Prong->SelectDplus(fDplusCuts);
  //if(fDebug) printf("ok3Prong: %d\n",(Int_t)ok3Prong);
  if(ok3Prong) {
    // get PID info from ESD
    Double_t esdpid0[5];
    postrack1->GetESDpid(esdpid0);
    Double_t esdpid1[5];
    negtrack->GetESDpid(esdpid1);
    Double_t esdpid2[5];
    postrack2->GetESDpid(esdpid2);


    Double_t esdpid[15];
    for(Int_t i=0;i<5;i++) {
      esdpid[i]   = esdpid0[i];
      esdpid[5+i] = esdpid1[i];
      esdpid[10+i] = esdpid2[i];
    }
    the3Prong->SetPID(3,esdpid);
  }

  if(ownPrimVertex) delete ownPrimVertex;	

  return the3Prong;
}
//----------------------------------------------------------------------------
AliAODRecoDecayHF4Prong* AliAnalysisVertexingHF::Make4Prong(
                             TObjArray *fourTrackArray,AliESDEvent *esd,
			     AliESDVertex *vertexp1n1,AliESDVertex *vertexp2n1,
			     Double_t dcap1n1,Double_t dcap1n2,Double_t dcap2n1,
			     Bool_t &ok4Prong) const
{
  // Make 4Prong candidates and check if they pass D0toKpipipi
  // reconstruction cuts
  // G.E.Bruno, R.Romita

  ok4Prong=kFALSE;

  Double_t px[4],py[4],pz[4],d0[4],d0err[4];//d0z[3];
  //Float_t d0z0[2],covd0z0[3];

  Double_t bfieldkG=dcap1n1*dcap1n2*dcap2n1; // TO BE CHANGED (done just to removed compilation warning about dca... not used)
  bfieldkG=(Double_t)esd->GetMagneticField();

  //charge
  Short_t charge=0;

  AliESDtrack *postrack1 = (AliESDtrack*)fourTrackArray->UncheckedAt(0);
  AliESDtrack *negtrack1 = (AliESDtrack*)fourTrackArray->UncheckedAt(1);
  AliESDtrack *postrack2 = (AliESDtrack*)fourTrackArray->UncheckedAt(2);
  AliESDtrack *negtrack2 = (AliESDtrack*)fourTrackArray->UncheckedAt(3);

  AliESDVertex *primVertex = fV1;

  postrack1->RelateToVertex(primVertex,bfieldkG,10.);
  negtrack1->RelateToVertex(primVertex,bfieldkG,10.);
  postrack2->RelateToVertex(primVertex,bfieldkG,10.);
  negtrack2->RelateToVertex(primVertex,bfieldkG,10.);

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
  //Bool_t okMassCut=kFALSE;
  //if(!okMassCut) if(SelectInvMass(2,3,px,py,pz)) okMassCut=kTRUE;
  //if(!okMassCut) {
  //  if(fDebug) printf(" candidate didn't pass mass cut\n");
  //  return 0x0;
  //}

  AliESDVertex *ownPrimVertex = 0;
  // primary vertex from *other* tracks in the event
  if(fRecoPrimVtxSkippingTrks || fRmTrksFromPrimVtx) {
    ownPrimVertex = OwnPrimaryVertex(4,fourTrackArray,esd);
    if(!ownPrimVertex) {
      return 0x0;
    } else {
      if(ownPrimVertex->GetNContributors()<2) {
        delete ownPrimVertex;
        return 0x0;
      } else {
        primVertex = ownPrimVertex;
      }
    }
  }

  // create the object AliAODRecoDecayHF4Prong
  AliVertexerTracks *vertexer = new AliVertexerTracks(esd->GetMagneticField());
  AliESDVertex* secVert4Prong=(AliESDVertex*)vertexer->VertexForSelectedTracks(fourTrackArray,kTRUE,kTRUE);
  delete vertexer; vertexer=NULL;
  Double_t pos[3],cov[6],sigmavert;
  secVert4Prong->GetXYZ(pos); // position
  secVert4Prong->GetCovMatrix(cov); //covariance matrix
  sigmavert=secVert4Prong->GetDispersion();

  AliAODVertex *secVert4PrAOD = new AliAODVertex(pos,cov,secVert4Prong->GetChi2toNDF());
  primVertex->GetXYZ(pos); // position
  primVertex->GetCovMatrix(cov); //covariance matrix
  AliAODVertex *primVertexAOD = new AliAODVertex(pos,cov,primVertex->GetChi2toNDF());
  //Double_t dca[6]={dcap1n1,dcap2n1,dcap1p2,0.,0.,0.}; //
  Double_t dca[6]={0.,0.,0.,0.,0.,0.}; //  modify it

  Double_t dist12=TMath::Sqrt((vertexp1n1->GetXv()-pos[0])*(vertexp1n1->GetXv()-pos[0])+(vertexp1n1->GetYv()-pos[1])*(vertexp1n1->GetYv()-pos[1])+(vertexp1n1->GetZv()-pos[2])*(vertexp1n1->GetZv()-pos[2]));
  Double_t dist23=TMath::Sqrt((vertexp2n1->GetXv()-pos[0])*(vertexp2n1->GetXv()-pos[0])+(vertexp2n1->GetYv()-pos[1])*(vertexp2n1->GetYv()-pos[1])+(vertexp2n1->GetZv()-pos[2])*(vertexp2n1->GetZv()-pos[2]));
  Double_t dist14=0.; // to be implemented
  Double_t dist34=0.; // to be implemented

  //AliAODRecoDecayHF4Prong *the4Prong = new AliAODRecoDecayHF4Prong(secVert4PrAOD,px,py,pz,d0,d0err,dca,sigmavert,dist12,dist23,charge);
  AliAODRecoDecayHF4Prong *the4Prong = new AliAODRecoDecayHF4Prong(secVert4PrAOD,px,py,pz,d0,d0err,dca,dist12,dist23,dist14,dist34,charge);
  the4Prong->SetOwnPrimaryVtx(primVertexAOD);


  // use the following two lines once AliAODRecoDecayHF4Prong::SelectD0 is available
  // select D0->Kpipipi
  //Int_t checkD0,checkD0bar;   
  // ok4Prong=the4Prong->SelectD0(fD04pCuts,checkD0,checkD0bar); 
  ok4Prong=kFALSE;  //for the time being ...


  // get PID info from ESD
  Double_t esdpid0[5];
  postrack1->GetESDpid(esdpid0);
  Double_t esdpid1[5];
  negtrack1->GetESDpid(esdpid1);
  Double_t esdpid2[5];
  postrack2->GetESDpid(esdpid2);
  Double_t esdpid3[5];
  negtrack2->GetESDpid(esdpid3);

  Double_t esdpid[20];
  for(Int_t i=0;i<5;i++) {
    esdpid[i]   = esdpid0[i];
    esdpid[5+i] = esdpid1[i];
    esdpid[10+i] = esdpid2[i];
    esdpid[15+i] = esdpid3[i];
  }
  the4Prong->SetPID(4,esdpid);

  if(ownPrimVertex) delete ownPrimVertex;

  return the4Prong;
}
//-----------------------------------------------------------------------------
AliESDVertex* AliAnalysisVertexingHF::OwnPrimaryVertex(Int_t ntrks,
						       TObjArray *trkArray,
						       AliESDEvent *esd) const
{
  // Returns primary vertex specific to this candidate
 
  AliVertexerTracks *vertexer1 = new AliVertexerTracks(esd->GetMagneticField());
  AliESDVertex *ownPrimVertex = 0;

  // recalculating the vertex
  if(fRecoPrimVtxSkippingTrks) { 
    if(strstr(fV1->GetTitle(),"VertexerTracksWithConstraint")) {
      Float_t diamondcovxy[3];
      esd->GetDiamondCovXY(diamondcovxy);
      Double_t pos[3]={esd->GetDiamondX(),esd->GetDiamondY(),0.};
      Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.};
      AliESDVertex *diamond = new AliESDVertex(pos,cov,1.,1);
      vertexer1->SetVtxStart(diamond);
      delete diamond; diamond=NULL;
      if(strstr(fV1->GetTitle(),"VertexerTracksWithConstraintOnlyFitter")) 
	vertexer1->SetOnlyFitter();
    }
    Int_t skipped[10];
    AliESDtrack *t = 0;
    for(Int_t i=0; i<ntrks; i++) {
      t = (AliESDtrack*)trkArray->UncheckedAt(i);
      skipped[i] = (Int_t)t->GetID();
    }
    vertexer1->SetSkipTracks(ntrks,skipped);
    ownPrimVertex = (AliESDVertex*)vertexer1->FindPrimaryVertex(esd); 
  }

  // removing the prongs tracks
  if(fRmTrksFromPrimVtx) { 
    TTree *rmTree = new TTree("trksTree","tree with tracks to be removed");
    AliESDtrack *esdTrack = 0;
    rmTree->Branch("tracks","AliESDtrack",&esdTrack);
    AliESDtrack *t = 0;
    for(Int_t i=0; i<ntrks; i++) {
      t = (AliESDtrack*)trkArray->UncheckedAt(i);
      esdTrack = new AliESDtrack(*t);
      rmTree->Fill();
      delete esdTrack;
    }
    Float_t diamondxy[2]={esd->GetDiamondX(),esd->GetDiamondY()};
    ownPrimVertex = vertexer1->RemoveTracksFromVertex(fV1,rmTree,diamondxy);
    delete rmTree; rmTree=NULL;
  }

  delete vertexer1; vertexer1=NULL;

  return ownPrimVertex;
}
//-----------------------------------------------------------------------------
void AliAnalysisVertexingHF::PrintStatus() const {
  // Print parameters being used

  printf("Preselections:\n");
  printf("    fITSrefit   = %d\n",(Int_t)fITSrefit);
  printf("    fBothSPD   = %d\n",(Int_t)fBothSPD);
  printf("    fMinITSCls   = %d\n",fMinITSCls);
  printf("    fMinPtCut   = %f GeV/c\n",fMinPtCut);
  printf("    fMind0rphiCut   = %f cm\n",fMind0rphiCut);
  if(fRecoPrimVtxSkippingTrks) printf("RecoPrimVtxSkippingTrks\n");
  if(fRmTrksFromPrimVtx) printf("RmTrksFromPrimVtx\n");
  if(fD0toKpi) {
    printf("Reconstruct D0->Kpi candidates with cuts:\n");
    printf("    |M-MD0| [GeV]    < %f\n",fD0toKpiCuts[0]);
    printf("    dca    [cm]  < %f\n",fD0toKpiCuts[1]);
    printf("    cosThetaStar     < %f\n",fD0toKpiCuts[2]);
    printf("    pTK     [GeV/c]    > %f\n",fD0toKpiCuts[3]);
    printf("    pTpi    [GeV/c]    > %f\n",fD0toKpiCuts[4]);
    printf("    |d0K|  [cm]  < %f\n",fD0toKpiCuts[5]);
    printf("    |d0pi| [cm]  < %f\n",fD0toKpiCuts[6]);
    printf("    d0d0  [cm^2] < %f\n",fD0toKpiCuts[7]);
    printf("    cosThetaPoint    > %f\n",fD0toKpiCuts[8]);
  }
  if(fJPSItoEle) {
    printf("Reconstruct J/psi from B candidates with cuts:\n");
    printf("    |M-MJPSI| [GeV]    < %f\n",fBtoJPSICuts[0]);
    printf("    dca    [cm]  < %f\n",fBtoJPSICuts[1]);
    printf("    cosThetaStar     < %f\n",fBtoJPSICuts[2]);
    printf("    pTP     [GeV/c]    > %f\n",fBtoJPSICuts[3]);
    printf("    pTN    [GeV/c]    > %f\n",fBtoJPSICuts[4]);
    printf("    |d0P|  [cm]  < %f\n",fBtoJPSICuts[5]);
    printf("    |d0N| [cm]  < %f\n",fBtoJPSICuts[6]);
    printf("    d0d0  [cm^2] < %f\n",fBtoJPSICuts[7]);
    printf("    cosThetaPoint    > %f\n",fBtoJPSICuts[8]);
  }
  if(f3Prong) {
    printf("Reconstruct 3 prong candidates.\n");
    printf("  D+ cuts:\n");
    printf("    |M-MD+| [GeV]    < %f\n",fDplusCuts[0]);
    printf("    pTK     [GeV/c]    > %f\n",fDplusCuts[1]);
    printf("    pTPi    [GeV/c]    > %f\n",fDplusCuts[2]);
    printf("    |d0K|  [cm]  > %f\n",fDplusCuts[3]);
    printf("    |d0Pi| [cm]  > %f\n",fDplusCuts[4]);
    printf("    dist12    [cm]  < %f\n",fDplusCuts[5]);
    printf("    sigmavert [cm]   < %f\n",fDplusCuts[6]);
    printf("    dist prim-sec [cm] > %f\n",fDplusCuts[7]);
    printf("    pM=Max{pT1,pT2,pT3} [GeV/c] > %f\n",fDplusCuts[8]);
    printf("    cosThetaPoint    > %f\n",fDplusCuts[9]);
    printf("    Sum d0^2 [cm^2]  > %f\n",fDplusCuts[10]);
    printf("    dca cut [cm]  < %f\n",fDplusCuts[11]);
  }

  return;
}
//-----------------------------------------------------------------------------
Bool_t AliAnalysisVertexingHF::SelectInvMass(Int_t decay,
					     Int_t nprongs,
					     Double_t *px,
					     Double_t *py,
					     Double_t *pz) const {
  // Check invariant mass cut

  Short_t dummycharge=0;
  Double_t *dummyd0 = new Double_t[nprongs];
  AliAODRecoDecay *rd = new AliAODRecoDecay(0x0,nprongs,dummycharge,px,py,pz,dummyd0);
  delete [] dummyd0;

  UInt_t pdg2[2],pdg3[3];
  Double_t mPDG,minv;

  Bool_t retval=kFALSE;
  switch (decay) 
    { 
    case 0:                  // D0->Kpi
      pdg2[0]=211; pdg2[1]=321;
      mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();
      minv = rd->InvMass(nprongs,pdg2);
      if(TMath::Abs(minv-mPDG)<fD0toKpiCuts[0]) retval=kTRUE;
      pdg2[0]=321; pdg2[1]=211;
      minv = rd->InvMass(nprongs,pdg2);
      if(TMath::Abs(minv-mPDG)<fD0toKpiCuts[0]) retval=kTRUE;
      break;
    case 1:                  // JPSI->ee
      pdg2[0]=11; pdg2[1]=11;
      mPDG=TDatabasePDG::Instance()->GetParticle(443)->Mass();
      minv = rd->InvMass(nprongs,pdg2);
      if(TMath::Abs(minv-mPDG)<fBtoJPSICuts[0]) retval=kTRUE;
      break;
    case 2:                  // D+->Kpipi
      pdg3[0]=211; pdg3[1]=321; pdg3[2]=211;
      mPDG=TDatabasePDG::Instance()->GetParticle(411)->Mass();
      minv = rd->InvMass(nprongs,pdg3);
      if(TMath::Abs(minv-mPDG)<fDplusCuts[0]) retval=kTRUE;
      break;
    default:
      break;
    }

  delete rd;

  return retval;
}
//-----------------------------------------------------------------------------
void AliAnalysisVertexingHF::SelectTracks(AliESDEvent *esd,
					  TObjArray &trksP,Int_t &nTrksP,
					  TObjArray &trksN,Int_t &nTrksN) const
{
  // Fill two TObjArrays with positive and negative tracks and 
  // apply single-track preselection

  nTrksP=0,nTrksN=0;

  Int_t entries = (Int_t)esd->GetNumberOfTracks();
 
  // transfer ITS tracks from ESD to arrays and to a tree
  for(Int_t i=0; i<entries; i++) {

    AliESDtrack *esdtrack = esd->GetTrack(i);
    UInt_t status = esdtrack->GetStatus();

    // require refit in ITS 
    if(fITSrefit && !(status&AliESDtrack::kITSrefit)) {
      if(fDebug) printf("track %d is not kITSrefit\n",i);
      continue;
    }

    // require minimum # of ITS points    
    if(esdtrack->GetNcls(0)<fMinITSCls)  {
      if(fDebug) printf("track %d has %d ITS cls\n",i,esdtrack->GetNcls(0));
      continue;
    }
    // require points on the 2 pixel layers
    if(fBothSPD) 
      if(!TESTBIT(esdtrack->GetITSClusterMap(),0) || 
	 !TESTBIT(esdtrack->GetITSClusterMap(),1)) continue;

    // single track selection
    if(!SingleTrkCuts(*esdtrack,esd->GetMagneticField())) continue;

    if(esdtrack->GetSign()<0) { // negative track
      trksN.AddLast(esdtrack);
      nTrksN++;
    } else {                 // positive track
      trksP.AddLast(esdtrack);
      nTrksP++;
    }

  } // loop on ESD tracks

  return;
}
//-----------------------------------------------------------------------------
void AliAnalysisVertexingHF::SetD0toKpiCuts(Double_t cut0,Double_t cut1,
				   Double_t cut2,Double_t cut3,Double_t cut4,
				   Double_t cut5,Double_t cut6,
				   Double_t cut7,Double_t cut8) 
{
  // Set the cuts for D0 selection
  fD0toKpiCuts[0] = cut0;
  fD0toKpiCuts[1] = cut1;
  fD0toKpiCuts[2] = cut2;
  fD0toKpiCuts[3] = cut3;
  fD0toKpiCuts[4] = cut4;
  fD0toKpiCuts[5] = cut5;
  fD0toKpiCuts[6] = cut6;
  fD0toKpiCuts[7] = cut7;
  fD0toKpiCuts[8] = cut8;

  return;
}
//-----------------------------------------------------------------------------
void AliAnalysisVertexingHF::SetD0toKpiCuts(const Double_t cuts[9]) 
{
  // Set the cuts for D0 selection

  for(Int_t i=0; i<9; i++) fD0toKpiCuts[i] = cuts[i];

  return;
}
//-----------------------------------------------------------------------------
void AliAnalysisVertexingHF::SetBtoJPSICuts(Double_t cut0,Double_t cut1,
				   Double_t cut2,Double_t cut3,Double_t cut4,
				   Double_t cut5,Double_t cut6,
				   Double_t cut7,Double_t cut8) 
{
  // Set the cuts for J/psi from B selection
  fBtoJPSICuts[0] = cut0;
  fBtoJPSICuts[1] = cut1;
  fBtoJPSICuts[2] = cut2;
  fBtoJPSICuts[3] = cut3;
  fBtoJPSICuts[4] = cut4;
  fBtoJPSICuts[5] = cut5;
  fBtoJPSICuts[6] = cut6;
  fBtoJPSICuts[7] = cut7;
  fBtoJPSICuts[8] = cut8;

  return;
}
//-----------------------------------------------------------------------------
void AliAnalysisVertexingHF::SetBtoJPSICuts(const Double_t cuts[9]) 
{
  // Set the cuts for J/psi from B selection

  for(Int_t i=0; i<9; i++) fBtoJPSICuts[i] = cuts[i];

  return;
}
//-----------------------------------------------------------------------------
void AliAnalysisVertexingHF::SetDplusCuts(Double_t cut0,Double_t cut1,
				   Double_t cut2,Double_t cut3,Double_t cut4,
				   Double_t cut5,Double_t cut6,
				   Double_t cut7,Double_t cut8,
				   Double_t cut9,Double_t cut10,Double_t cut11)
{
  // Set the cuts for Dplus->Kpipi selection
  fDplusCuts[0] = cut0;
  fDplusCuts[1] = cut1;
  fDplusCuts[2] = cut2;
  fDplusCuts[3] = cut3;
  fDplusCuts[4] = cut4;
  fDplusCuts[5] = cut5;
  fDplusCuts[6] = cut6;
  fDplusCuts[7] = cut7;
  fDplusCuts[8] = cut8;
  fDplusCuts[9] = cut9;
  fDplusCuts[10] = cut10;
  fDplusCuts[11] = cut11;

  return;
}
//-----------------------------------------------------------------------------
void AliAnalysisVertexingHF::SetDplusCuts(const Double_t cuts[12]) 
{
  // Set the cuts for Dplus selection

  for(Int_t i=0; i<12; i++) fDplusCuts[i] = cuts[i];

  return;
}
//-----------------------------------------------------------------------------
Bool_t AliAnalysisVertexingHF::SingleTrkCuts(AliESDtrack& trk, 
					     Double_t bfieldkG) const 
{
  // Check if track passes some kinematical cuts  
  // Magnetic field "bfieldkG" (kG)

  if(TMath::Abs(1./trk.GetParameter()[4]) < fMinPtCut) {
    printf("pt %f\n",1./trk.GetParameter()[4]);
    return kFALSE;
  }
  trk.RelateToVertex(fV1,bfieldkG,10.);
  Float_t d0z0[2],covd0z0[3];
  trk.GetImpactParameters(d0z0,covd0z0);
  if(TMath::Abs(d0z0[0]) < fMind0rphiCut) {
    printf("d0rphi %f\n",TMath::Abs(d0z0[0]));
    return kFALSE;
  }

  return kTRUE;
}
//-----------------------------------------------------------------------------



