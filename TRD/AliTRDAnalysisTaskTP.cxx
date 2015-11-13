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

/* $Id: AliTRDAnalysisTaskTP.cxx 42548 2010-07-27 08:10:51Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Track point maker for the alignment of TRD                            //
//                                                                        //
//  Author:                                                               //
//     Sebastian Huber (S.Huber@gsi.de)                                   // 
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <TROOT.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include <TLinearFitter.h>

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAlignObjParams.h"
#include "AliTrackPointArray.h"
#include "AliESDfriend.h"
#include "AliESDtrack.h"
#include "AliESDfriendTrack.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TCollection.h"
#include "AliLog.h"
#include "AliGeomManager.h"

#include "AliTRDAnalysisTaskTP.h"

ClassImp(AliTRDAnalysisTaskTP)

AliTRDAnalysisTaskTP::AliTRDAnalysisTaskTP()
  :AliAnalysisTaskSE(),
  fArrHists(0x0),
  fArrTTree(0x0),
  fTree(NULL),
  fESD(0x0),
  fModpop(0x0),
  fBug(0x0),
  fNevents(0x0),
  fNtracks(0x0),
  fNAcceptedTracks(0),
  fArray(0x0),
  fFile(0x0)
{
  //
  // Default constructor
  //

}

//____________________________________________________________
AliTRDAnalysisTaskTP::AliTRDAnalysisTaskTP(const char *name) :
  AliAnalysisTaskSE(name),
  fArrHists(0x0),
  fArrTTree(0x0),
  fTree(0x0),
  fESD(0x0),
  fModpop(0x0),
  fBug(0x0),
  fNevents(0),
  fNtracks(0),
  fNAcceptedTracks(0),
  fArray(0x0),
  fFile(0x0)
{
  //
  // Constructor
  //

  DefineOutput(1, TTree::Class());
  DefineOutput(2, TObjArray::Class());

}

//____________________________________________________________
AliTRDAnalysisTaskTP::~AliTRDAnalysisTaskTP() 
{
  //
  // destructor
  //

}

//____________________________________________________________
void AliTRDAnalysisTaskTP::UserCreateOutputObjects() 
{
  //
  // Create the output objects
  //

  AliAlignObjParams alobj;  // initialize align obj.  
  TString option = GetOption();

  if (!fArrHists) fArrHists=new TObjArray;

  fModpop = new TH2D("modpop","modpop",90,-0.5,89.5,30,-0.5,29.5);
  fModpop->SetXTitle("module nr");
  fModpop->SetYTitle("layer nr");
  fArrHists->Add(fModpop);

  OpenFile(1);
  fTree = new TTree("spTree", "Tree with track space point arrays");
  fTree->Branch("SP","AliTrackPointArray", &fArray);

}

//____________________________________________________________
void AliTRDAnalysisTaskTP::UserExec(Option_t *) 
{
  //
  // Exec function
  //

  //AliESDEvent *fESD = dynamic_cast<AliESDEvent *>(fInputEvent);
  fESD = dynamic_cast<AliESDEvent *>(fInputEvent);
  if(!fESD){
  //cout << "ERROR: fESDs not available " << endl;
  return;
  }

  fESDfriend = dynamic_cast<AliESDfriend*>(fESD->FindListObject("AliESDfriend"));
  if (!fESDfriend) {
  //cout << "ERROR: fESDfriends not available " << endl;
  return;
  }

  TLinearFitter fitter(2, "pol1");
  TLinearFitter fitterz(2, "pol1");

  // track cuts
  int tpc = 0; // require tpc
  int ptu = 0; // require certain pt's (magnetic field and tpc presumably on)

  const Float_t kMaxDelta       = 1;
  const Float_t kMinNcl         = 60;
  const Float_t kMinPtLow       = 0.2;
  const Float_t kMinNclLow      = 100;
  const Float_t kMinPt0         = 2;
  const Float_t kMinPt          = 0;
  UInt_t status = AliESDtrack::kTRDrefit; 
  if (tpc) status |= AliESDtrack::kTPCrefit; 

  const Float_t kMinRadius2  = 2*2;
  const Float_t kMaxRadius2  = 400*400;
  const Float_t kDeadSpace   = 4;
  const Float_t kTan = TMath::Tan(10*TMath::DegToRad());
  Int_t ntracks = fESD->GetNumberOfTracks();
  const AliTrackPointArray *array=0;
  AliTrackPointArray *tmpArray = 0;
  // trdarray contains all trd points in this event, for duplication detection
  AliTrackPointArray *trdarray = new AliTrackPointArray(1000);
  int ntrdarray = 0;

  for (Int_t itrack=0; itrack < ntracks; itrack++) { //track loop

    AliESDtrack * track = fESD->GetTrack(itrack);
    fNtracks++;
    if (!track) continue;

    if (track->GetP() < kMinPt) continue;
    if (track->GetKinkIndex(0)!=0) continue;
    if (tpc) if (track->GetTPCNcls()<kMinNcl) continue;
    if (ptu) if (track->GetP() < kMinPtLow) continue;
    if (ptu) if (track->GetP() < kMinPt0 && track->GetTPCNcls()<kMinNclLow) continue;

    AliESDfriendTrack *friendtrack = track->GetFriendTrack(itrack);
    if (!friendtrack){ 
    continue;
    }

    array = friendtrack->GetTrackPointArray();
    if (!array) continue;
    Int_t npoints = array->GetNPoints();
    if (tmpArray) delete tmpArray;
    tmpArray = new AliTrackPointArray(npoints);
    Int_t current = 0;
    int ntpc = 0;
    int ntrd = 0;
    for (Int_t ipoint=0; ipoint<npoints; ipoint++){  

      AliTrackPoint p;
      array->GetPoint(p, ipoint);

      UShort_t volid = array->GetVolumeID()[ipoint];
      Int_t iModule;
      AliGeomManager::ELayerID layer = AliGeomManager::VolUIDToLayer(volid,iModule);
      if ((layer < AliGeomManager::kFirstLayer) || (layer >= AliGeomManager::kLastLayer)) continue;

      if ((iModule >= AliGeomManager::LayerSize(layer)) || (iModule < 0)) continue;

      Float_t r2 = p.GetX()*p.GetX()+p.GetY()*p.GetY();
      if ( r2<kMinRadius2 || r2 > kMaxRadius2 ) continue;



      if (layer>=AliGeomManager::kTPC1 && layer<=AliGeomManager::kTPC2){
	if (p.GetCov()[0]<0 || p.GetCov()[3]<0 ||  p.GetCov()[5]<0) continue;

	AliTrackPoint& plocal = p.MasterToLocal();
	Double_t ylocal  = plocal.GetY();
	Double_t zlocal  = plocal.GetZ();
	Double_t xlocal  = plocal.GetX();
	Float_t edgey = TMath::Abs(plocal.GetX()*kTan);
	Int_t nclose=0;
	fitter.ClearPoints();
	fitterz.ClearPoints();
	for (Int_t jpoint=ipoint-7; jpoint<=ipoint+7; jpoint++){
	  if (jpoint<0 || jpoint>=npoints) continue;
	  if (ipoint==jpoint) continue;
	  UShort_t volidL = array->GetVolumeID()[jpoint];
	  if (volidL!=volid) continue;
	  AliTrackPoint pc;	
	  array->GetPoint(pc, jpoint);
	  AliTrackPoint &pcl=  pc.MasterToLocal();		  
	  Double_t dx = pcl.GetX()-xlocal;
	  fitter.AddPoint(&dx,pcl.GetY(),1);	  
	  fitterz.AddPoint(&dx,pcl.GetZ(),1);	  
	  nclose++;
	}
	if (nclose<6) continue;

	fitter.Eval();
	fitterz.Eval();
	Double_t fity =fitter.GetParameter(0); 
	Double_t fitz =fitterz.GetParameter(0); 
	if (TMath::Abs(ylocal-fity)>kMaxDelta) continue;
	if (TMath::Abs(zlocal-fitz)>kMaxDelta) continue;
	if (TMath::Abs(fity)>edgey-kDeadSpace) continue;
	ntpc++;
      }

      if (layer>=AliGeomManager::kTRD1 && layer<=AliGeomManager::kTRD6){

	trdarray->AddPoint(ntrdarray++,&p);
	fModpop->Fill(iModule,layer);
	ntrd++;
      }
      tmpArray->AddPoint(current,&p);
      current++;
    }
    if (ntpc < 100) continue;
    if (ntrd < 4) continue;
    if (fArray) delete fArray;
    fArray = new AliTrackPointArray(current);
    for (Int_t ipoint=0; ipoint<current; ipoint++){
      AliTrackPoint p;
      tmpArray->GetPoint(p, ipoint);
      fArray->AddPoint(ipoint,&p);
    }
    fNAcceptedTracks++;
    fTree->Fill();
  }
  delete trdarray;
  fNevents++;
  PostData(1,fTree);
  PostData(2,fArrHists);
}

//____________________________________________________________
void AliTRDAnalysisTaskTP::Terminate(Option_t */*option*/) 
{
  //
  // Terminate
  //

}


