/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// AliAnalysisTask to study the TPC-ITS track matching
//
// Author: A.Dainese, andrea.dainese@pd.infn.it
/////////////////////////////////////////////////////////////

#include <TTree.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TVector3.h>
#include <TGeoManager.h>

#include "AliLog.h"
#include "AliGeomManager.h"
#include "AliTrackPointArray.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliExternalTrackParam.h"
#include "AliTrackMatchingTPCITSCosmics.h"


ClassImp(AliTrackMatchingTPCITSCosmics)


//________________________________________________________________________
AliTrackMatchingTPCITSCosmics::AliTrackMatchingTPCITSCosmics(const char *name):
AliAnalysisTask(name,"task"),
fOnlySPDFO(kFALSE),
fReadHLTESD(kFALSE),
fGeometryFileName("geometry.root"),
fESD(0),
fList(0),
fHistEvCount(0),
fntTrks(0),
fntPairs(0),
fntITSPairs(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList
  DefineOutput(0,TList::Class());  //My private output
}

//________________________________________________________________________
AliTrackMatchingTPCITSCosmics::~AliTrackMatchingTPCITSCosmics()
{
  // Destructor
  if (fList) {
    delete fList;
    fList = 0;
  }
  if (fHistEvCount) {
    delete fHistEvCount;
    fHistEvCount = 0;
  }
  if (fntTrks) {
    delete fntTrks;
    fntTrks = 0;
  }
  if (fntPairs) {
    delete fntPairs;
    fntPairs = 0;
  }
  if (fntITSPairs) {
    delete fntITSPairs;
    fntITSPairs = 0;
  }
}  

//________________________________________________________________________
void AliTrackMatchingTPCITSCosmics::ConnectInputData(Option_t *) 
{
  // Connect ESD
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if(!tree) {
    printf("ERROR: Could not read chain from input slot 0\n");
  } else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if(!esdH) {
      printf("ERROR: Could not get ESDInputHandler\n");
    } else {
      fESD = esdH->GetEvent();
    }
  }
  
  return;
}

//________________________________________________________________________
void AliTrackMatchingTPCITSCosmics::Init()
{
  // Initialization

  return;
}

//________________________________________________________________________
void AliTrackMatchingTPCITSCosmics::CreateOutputObjects()
{
  // Create the output container
  //

  // load the geometry  
  if(!gGeoManager) {    
    AliGeomManager::LoadGeometry(fGeometryFileName.Data());
    if(!gGeoManager) { 
      printf("AliTrackMatchingTPCITSCosmics::CreateOutputObjects(): no geometry loaded \n");
      return;
    }
  }

  // Several histograms are more conveniently managed in a TList
  fList = new TList();
  fList->SetOwner();

  fHistEvCount = new TH1F("fHistEvCount","0: all, 1: SPDFO, 2: vtx, 3: TPCtrks, 4: ITStrks, 5: HLT, 6: HLT&TPCtrks, 7: HLT&ITStrks",11,-0.5,10.5);
  fList->Add(fHistEvCount);

  fntTrks = new TNtuple("fntTrks","TPC-ITS matching","ptTPC:nClsTPC:nClsITSSA:nClsITSMI:SSD1:SSD2:phi:z:dx:dy:dz:drphi:dphi:dtgl");
  fList->Add(fntTrks);
  fntPairs = new TNtuple("fntPairs","pairs at vertex","nClsITSSAin:nClsITSSAout:nClsITSMIin:nClsITSMIout:dxyITSSA:dzITSSA:dxyITSMI:dzITSMI:ptITSSAin:ptITSSAout:ptITSMIin:ptITSMIout:sigmad0MIin:sigmad0MIout");
  fList->Add(fntPairs);
  fntITSPairs = new TNtuple("fntITSPairs","pairs at vertex","nClsITSSAin:nClsITSSAout:dxyITSSA:dzITSSA:ptITSSAin:ptITSSAout:d0mu");
  fList->Add(fntITSPairs);

  return;
}

//________________________________________________________________________
void AliTrackMatchingTPCITSCosmics::Exec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // write ITS AliTrackPoints for selected tracks to fspTree
  
  // check the geometry  
  if(!gGeoManager) { 
    printf("AliTrackMatchingTPCITSCosmics::Exec(): no geometry loaded \n");
    return;
  }

  if(!fESD) {
    printf("AliTrackMatchingTPCITSCosmics::Exec(): no ESD \n");
    return;
  } 


  fHistEvCount->Fill(0);
  PostData(0,fList);

  // check if event is triggered by HLT
  Bool_t hltTrigg=kFALSE;
  if(fReadHLTESD) {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    AliESDEvent *hltESD = esdH->GetHLTEvent();
    if(!hltESD) {
      printf("AliTrackMatchingTPCITSCosmics::Exec(): no HLT ESD \n");
      return;
    } 
    if(hltESD->IsHLTTriggerFired()) {fHistEvCount->Fill(5);hltTrigg=kTRUE;}
  }

  
  TString triggeredClass = fESD->GetFiredTriggerClasses(); 
  if(triggeredClass.Contains("C0SCO-ABCE-NOPF-CENT")) {
    fHistEvCount->Fill(1);
  } else {
    if(fOnlySPDFO) {  PostData(0,fList); return; }
  }
  
  Int_t ntracks = fESD->GetNumberOfTracks();
  printf("CONTR: %d\n",fESD->GetVertex()->GetNContributors());
  if(fESD->GetVertex()->GetNContributors()>=0) {
    fHistEvCount->Fill(2);
  } 
  
  Int_t nTrksITSSA=0,nTrksTPC=0,nTrksTPCITS=0;
  Int_t idxITSSA[10000],idxTPC[10000];
  Int_t idxITSSAin=-1,idxITSSAout=-1;
  Int_t idxTPCin=-1,idxTPCout=-1;
  for(Int_t itr=0; itr<ntracks; itr++) {
    AliESDtrack *t = fESD->GetTrack(itr);
    Int_t nClsTPC = t->GetNcls(1);
    Int_t nClsITS = t->GetNcls(0);
    if(nClsITS>=2 && nClsTPC==0) { // ITS SA
      idxITSSA[nTrksITSSA]=itr;
      nTrksITSSA++;
      //printf("%d\n",itr);
    }
    if(nClsTPC>=50) {              // TPC
      idxTPC[nTrksTPC]=itr;
      nTrksTPC++;
    }
    if(nClsTPC>=50 && nClsITS>=2) {// TPC+ITS
      nTrksTPCITS++;
      /*printf("  --> TPC+ITS: %d + %d (%d,%d,%d,%d,%d,%d) pt = %f\n",nClsTPC,nClsITS,
	(Int_t)(t->HasPointOnITSLayer(0)),
	(Int_t)(t->HasPointOnITSLayer(1)),
	(Int_t)(t->HasPointOnITSLayer(2)),
	(Int_t)(t->HasPointOnITSLayer(3)),
	(Int_t)(t->HasPointOnITSLayer(4)),
	(Int_t)(t->HasPointOnITSLayer(5)),
	t->Pt());*/
    }
  }
  
  //printf("nTrksTPC %d  nTrksITSSA %d\n",nTrksTPC,nTrksITSSA);
  if(nTrksITSSA>0) fHistEvCount->Fill(4);
  if(hltTrigg && nTrksITSSA>0) fHistEvCount->Fill(7);
  if(nTrksTPC>0)   fHistEvCount->Fill(3);
  if(hltTrigg && nTrksTPC>0)   fHistEvCount->Fill(6);
  if(nTrksTPCITS>0)   fHistEvCount->Fill(8);
  if(hltTrigg && nTrksTPCITS>0)   fHistEvCount->Fill(9);
  
  if(nTrksTPC>2 || nTrksTPC==0 || nTrksITSSA>2 || nTrksITSSA==0) { PostData(0,fList); return; }
  
  for(Int_t itr=0; itr<nTrksTPC; itr++) {
    AliESDtrack *t = fESD->GetTrack(idxTPC[itr]);
    if(t->Py()>0) idxTPCin=idxTPC[itr];
    if(t->Py()<0) idxTPCout=idxTPC[itr];
  }
  
  
  for(Int_t itr=0; itr<nTrksITSSA; itr++) {
    AliESDtrack *t = fESD->GetTrack(idxITSSA[itr]);
    if(t->Py()>0) idxITSSAin=idxITSSA[itr];
    if(t->Py()<0) idxITSSAout=idxITSSA[itr];
  }
  
  Double_t xyzITS[3],xyzTPC[3];

  // analysis for inward track
  if(idxITSSAin>-1 && idxTPCin>-1) { 
    //printf(" %d\n",idxITSSAin);
    AliESDtrack *tITS = fESD->GetTrack(idxITSSAin);
    const AliExternalTrackParam *outerITS = tITS->GetOuterParam();
    xyzITS[0]=-1000.;xyzITS[1]=-1000.;xyzITS[2]=-1000.;
    if(outerITS) outerITS->GetXYZAt(50.,fESD->GetMagneticField(),xyzITS);
    AliESDtrack *tTPC = fESD->GetTrack(idxTPCin);
    AliExternalTrackParam *innerTPC = new AliExternalTrackParam(*(tTPC->GetInnerParam()));
    innerTPC->Rotate(tITS->GetAlpha());
    innerTPC->GetXYZAt(50.,fESD->GetMagneticField(),xyzTPC);
    fntTrks->Fill(innerTPC->Pt(),
		 tTPC->GetNcls(1),
		 tITS->GetNcls(0),
		 tTPC->GetNcls(0),
		 tITS->HasPointOnITSLayer(4),
		 tITS->HasPointOnITSLayer(5),
		 TMath::ATan2(xyzITS[1],xyzITS[0]),
		 xyzITS[2],
		 xyzITS[0]-xyzTPC[0],
		 xyzITS[1]-xyzTPC[1],
		 xyzITS[2]-xyzTPC[2],
		 0.5*(TMath::Sqrt(xyzITS[0]*xyzITS[0]+xyzITS[1]*xyzITS[1])+TMath::Sqrt(xyzTPC[0]*xyzTPC[0]+xyzTPC[1]*xyzTPC[1]))*(TMath::ATan2(xyzITS[1],xyzITS[0])-TMath::ATan2(xyzTPC[1],xyzTPC[0])),
		 TMath::ATan2(tITS->Py(),tITS->Px())-TMath::ATan2(innerTPC->Py(),innerTPC->Px()),
		 tITS->Pz()/tITS->Pt()-innerTPC->Pz()/innerTPC->Pt());
    delete innerTPC; innerTPC=NULL;
  }

  // analysis for outward track
  if(idxITSSAout>-1 && idxTPCout>-1) { 
    //printf(" %d\n",idxITSSAout);
    AliESDtrack *tITS = fESD->GetTrack(idxITSSAout);
    const AliExternalTrackParam *outerITS = tITS->GetOuterParam();
    xyzITS[0]=-1000.;xyzITS[1]=-1000.;xyzITS[2]=-1000.;
    if(outerITS) outerITS->GetXYZAt(50.,fESD->GetMagneticField(),xyzITS);
    AliESDtrack *tTPC = fESD->GetTrack(idxTPCout);
    AliExternalTrackParam *innerTPC = new AliExternalTrackParam(*(tTPC->GetInnerParam()));
    innerTPC->Rotate(tITS->GetAlpha());
    innerTPC->GetXYZAt(50.,fESD->GetMagneticField(),xyzTPC);
    fntTrks->Fill(innerTPC->Pt(),
		 tTPC->GetNcls(1),
		 tITS->GetNcls(0),
		 tTPC->GetNcls(0),
		 tITS->HasPointOnITSLayer(4),
		 tITS->HasPointOnITSLayer(5),
		 TMath::ATan2(xyzITS[1],xyzITS[0]),
		 xyzITS[2],
		 xyzITS[0]-xyzTPC[0],
		 xyzITS[1]-xyzTPC[1],
		 xyzITS[2]-xyzTPC[2],
		 0.5*(TMath::Sqrt(xyzITS[0]*xyzITS[0]+xyzITS[1]*xyzITS[1])+TMath::Sqrt(xyzTPC[0]*xyzTPC[0]+xyzTPC[1]*xyzTPC[1]))*(TMath::ATan2(xyzITS[1],xyzITS[0])-TMath::ATan2(xyzTPC[1],xyzTPC[0])),
		 TMath::ATan2(tITS->Py(),tITS->Px())-TMath::ATan2(innerTPC->Py(),innerTPC->Px()),
		 tITS->Pz()/tITS->Pt()-innerTPC->Pz()/innerTPC->Pt()
		 );
    delete innerTPC; innerTPC=NULL;
  }

  // inward-outward track ITS only
  if(idxITSSAin>-1 && idxITSSAout>-1) {
    AliESDtrack *tITSin = fESD->GetTrack(idxITSSAin);
    AliESDtrack *tITSout = fESD->GetTrack(idxITSSAout);
    if(tITSin->HasPointOnITSLayer(0) && 
       tITSin->HasPointOnITSLayer(1) &&
       tITSout->HasPointOnITSLayer(0) &&
       tITSout->HasPointOnITSLayer(1)) {
      Double_t alpha = TMath::ATan2(tITSin->Py(),tITSin->Px());
      tITSin->Propagate(alpha,0.,fESD->GetMagneticField());
      tITSout->Propagate(alpha,0.,fESD->GetMagneticField());
      Float_t d0[2],z0[2];
      d0[0] = tITSin->GetY();
      z0[0] = tITSin->GetZ();
      d0[1] = tITSout->GetY();
      z0[1] = tITSout->GetZ();
      Float_t dxyITS = -(d0[0]-d0[1]);
      Float_t dzITS = z0[0]-z0[1];
      fntITSPairs->Fill(tITSin->GetNcls(0),tITSout->GetNcls(0),dxyITS,dzITS,tITSin->Pt(),tITSout->Pt(),TMath::Abs(d0[0]));
    }
  }

  // inward-outward track TPC+ITS
  if(idxITSSAin>-1 && idxTPCin>-1 && idxITSSAout>-1 && idxTPCout>-1) {
    AliESDtrack *tTPCin = fESD->GetTrack(idxTPCin);
    AliESDtrack *tITSin = fESD->GetTrack(idxITSSAin);
    AliESDtrack *tTPCout = fESD->GetTrack(idxTPCout);
    AliESDtrack *tITSout = fESD->GetTrack(idxITSSAout);
    if(tTPCin->HasPointOnITSLayer(0) && 
       tTPCin->HasPointOnITSLayer(1) &&
       tITSin->HasPointOnITSLayer(0) && 
       tITSin->HasPointOnITSLayer(1) && 
       tTPCout->HasPointOnITSLayer(0) &&
       tTPCout->HasPointOnITSLayer(1) &&
       tITSout->HasPointOnITSLayer(0) &&
       tITSout->HasPointOnITSLayer(1)) {
      Double_t alpha = TMath::ATan2(tTPCin->Py(),tTPCin->Px());
      tTPCin->Propagate(alpha,0.,fESD->GetMagneticField());
      tTPCout->Propagate(alpha,0.,fESD->GetMagneticField());
      Float_t d0[2],z0[2];
      d0[0] = tTPCin->GetY();
      z0[0] = tTPCin->GetZ();
      d0[1] = tTPCout->GetY();
      z0[1] = tTPCout->GetZ();
      Float_t dxyTPC = -(d0[0]-d0[1]);
      Float_t dzTPC = z0[0]-z0[1];
      alpha = TMath::ATan2(tITSin->Py(),tITSin->Px());
      tITSin->Propagate(alpha,0.,fESD->GetMagneticField());
      tITSout->Propagate(alpha,0.,fESD->GetMagneticField());
      d0[0] = tITSin->GetY();
      z0[0] = tITSin->GetZ();
      d0[1] = tITSout->GetY();
      z0[1] = tITSout->GetZ();
      Float_t dxyITS = -(d0[0]-d0[1]);
      Float_t dzITS = z0[0]-z0[1];
      fntPairs->Fill(tITSin->GetNcls(0),tITSout->GetNcls(0),tTPCin->GetNcls(0),tTPCout->GetNcls(0),dxyITS,dzITS,dxyTPC,dzTPC,tITSin->Pt(),tITSout->Pt(),tTPCin->Pt(),tTPCout->Pt(),tTPCin->GetSigmaY2(),tTPCout->GetSigmaY2());
    }
  }

  PostData(0,fList);

  return;
}

//________________________________________________________________________
void AliTrackMatchingTPCITSCosmics::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  AliDebug(2,"AliTrackMatchingTPCITSCosmics: Terminate() \n");

  fList = dynamic_cast<TList*> (GetOutputData(0));
  if (!fList) {     
    printf("ERROR: fList not available\n");
    return;
  }

  return;
}

