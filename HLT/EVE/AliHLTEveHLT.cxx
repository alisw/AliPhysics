/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Svein Lindal <slindal@fys.uio.no   >                  *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliHLTEvePhos.cxx
/// @author Svein Lindal <slindal@fys.uio.no>
/// @brief  HLT class for the HLT EVE display

#include "AliHLTEveHLT.h"
#include "AliHLTHOMERBlockDesc.h"
#include "AliHLTEveBase.h"
#include "AliEveHOMERManager.h"
#include "TEveManager.h"
#include "TEvePointSet.h"
#include "TEveTrack.h"
#include "TCanvas.h"
#include "AliESDEvent.h"
#include "TEveTrackPropagator.h"
#include "AliEveTrack.h"
#include "TEveVSDStructs.h"
#include "TString.h"
#include "TPCLib/tracking-ca/AliHLTTPCCATrackParam.h"
#include "TPCLib/tracking-ca/AliHLTTPCCATrackConvertor.h"
#include "AliEveMagField.h"
#include "TH1F.h"
#include "TH2F.h"

ClassImp(AliHLTEveHLT)

AliHLTEveHLT::AliHLTEveHLT() : 
  AliHLTEveBase(), 
  fTrueField(kFALSE),
  fUseIpOnFailedITS(kFALSE),
  fUseRkStepper(kFALSE),
  fTrackList(NULL),
  fTrCanvas(NULL),
  fHistPt(NULL), 
  fHistP(NULL), 
  fHistEta(NULL),
  fHistTheta(NULL),
  fHistPhi(NULL),
  fHistnClusters(NULL),
  fHistMult(NULL)
{
  // Constructor.
  CreateHistograms();

}

AliHLTEveHLT::~AliHLTEveHLT()
{
  //Destructor, not implemented
  if(fTrackList)
    delete fTrackList;
  fTrackList = NULL;
}

void AliHLTEveHLT::CreateHistograms(){
  //See header file for documentation
  fHistPt        = new TH1F("fHistPt",       "transverse momentum",    100, 0, 10); // KK   
  fHistP         = new TH1F("fHistP",        "signed momentum",        100,-7,  7);	   
  fHistEta       = new TH1F("fHistEta",      "pseudorapidity",         100,-2,  2);	   
  fHistTheta     = new TH1F("fHistTheta",    "polar angle",            180, 0,180);   
  fHistPhi	 = new TH1F("fHistPhi",      "azimuthal angle",        180, 0,360);   
  fHistnClusters = new TH1F("fHistnClusters","TPC clusters per track", 160, 0,160);
  fHistMult      = new TH1F("fHistMult",     "event track multiplicity",50, 0, 50);    
  
  fHistPt   ->SetXTitle("p_{t} (GeV/c)");   // KK
  fHistP    ->SetXTitle("P*charge (GeV/c)");
  fHistEta  ->SetXTitle("#eta");
  fHistTheta->SetXTitle("#theta (degrees)");
  fHistPhi  ->SetXTitle("#phi (degrees)");

}

void AliHLTEveHLT::ProcessBlock(AliHLTHOMERBlockDesc * block) {
  //See header file for documentation
  if ( ! block->GetDataType().CompareTo("ALIESDV0") ) {
    if(!fTrackList) CreateTrackList();
    ProcessEsdBlock(block, fTrackList);
  } 
  
  else if ( ! block->GetDataType().CompareTo("ROOTTOBJ") ) {
    //processROOTTOBJ( block, gHLTText );
  } 

  else if ( ! block->GetDataType().CompareTo("HLTRDLST") ) {
    //processHLTRDLST( block );
  } 

  else if ( !block->GetDataType().CompareTo("ROOTHIST") ) {      
    if( !fCanvas ) { 
      fCanvas = CreateCanvas("Primary Vertex", "Primary Vertex");
      fCanvas->Divide(3, 2);
    }
    ProcessHistograms( block , fCanvas);
  }
  
}


void AliHLTEveHLT::UpdateElements() {
  //See header file for documentation
  if(fCanvas) fCanvas->Update();
  DrawHistograms();
  if(fTrackList) fTrackList->ElementChanged();
}

void AliHLTEveHLT::ResetElements(){
    //See header file for documentation
  if(fTrackList) fTrackList->DestroyElements();
  fHistoCount = 0;

}

void AliHLTEveHLT::ProcessHistograms(AliHLTHOMERBlockDesc * block, TCanvas * canvas) {
  //See header file for documentation
  if ( ! block->GetClassName().CompareTo("TH1F")) {
    TH1F* histo = reinterpret_cast<TH1F*>(block->GetTObject());
    if( histo ){
      TString name(histo->GetName());
      if( !name.CompareTo("primVertexZ") ){
	canvas->cd(2);
	histo->Draw();
      }else if( !name.CompareTo("primVertexX") ){
	canvas->cd(3);
	histo->Draw();
      }else if( !name.CompareTo("primVertexY") ){
	canvas->cd(4);
	histo->Draw();
      }
    }
  }  else if ( ! block->GetClassName().CompareTo("TH2F")) {
    TH2F *hista = reinterpret_cast<TH2F*>(block->GetTObject());
    if (hista ){
       TString name(hista->GetName());
       if( !name.CompareTo("primVertexXY")) {      
	 canvas->cd(1);
	 hista->Draw();
       }
    }
  }
  canvas->cd();




}

void AliHLTEveHLT::CreateTrackList() {
  //See header file for documentation
  fTrackList = new TEveTrackList("ESD Tracks");
  fTrackList->SetMainColor(6);
  gEve->AddElement(fTrackList);
}


void AliHLTEveHLT::ProcessEsdBlock( AliHLTHOMERBlockDesc * block, TEveTrackList * cont ) {
  //See header file for documentation

  AliESDEvent* esd = (AliESDEvent *) (block->GetTObject());
  esd->GetStdContent();
  
  SetUpTrackPropagator(cont->GetPropagator(),-0.1*esd->GetMagneticField(), 520);

  for (Int_t iter = 0; iter < esd->GetNumberOfTracks(); ++iter) {
    AliEveTrack* track = dynamic_cast<AliEveTrack*>(MakeEsdTrack(esd->GetTrack(iter), cont));
    cont->AddElement(track);
   
    fHistPt->Fill(esd->GetTrack(iter)->Pt());   // KK
    fHistP->Fill(esd->GetTrack(iter)->P()*esd->GetTrack(iter)->Charge());
    fHistEta->Fill(esd->GetTrack(iter)->Eta());
    fHistTheta->Fill(esd->GetTrack(iter)->Theta()*TMath::RadToDeg());
    fHistPhi->Fill(esd->GetTrack(iter)->Phi()*TMath::RadToDeg());
    if(esd->GetTrack(iter)->GetStatus()&AliESDtrack::kTPCin || (esd->GetTrack(iter)->GetStatus()&AliESDtrack::kTPCin && esd->GetTrack(iter)->GetStatus()&AliESDtrack::kITSin)){
       fHistnClusters->Fill(esd->GetTrack(iter)->GetTPCNcls());  
    }
  }
  
  fHistMult->Fill(esd->GetNumberOfTracks()); // KK
  
  
  cont->SetTitle(Form("N=%d", esd->GetNumberOfTracks()) );
  cont->MakeTracks();

}


void AliHLTEveHLT::DrawHistograms(){
  //See header file for documentation

  if (!fTrCanvas) {
    fTrCanvas = CreateCanvas("TPC Tr QA", "TPC Track QA");
    fTrCanvas->Divide(4, 2);
  }

  Int_t icd = 1;
  fTrCanvas->cd(icd++);
  fHistPt->Draw();
  fTrCanvas->cd(icd++);
  fHistP->Draw();
  fTrCanvas->cd(icd++);
  fHistEta->Draw();
  fTrCanvas->cd(icd++);
  fHistTheta->Draw();
  fTrCanvas->cd(icd++);
  fHistPhi->Draw();
  fTrCanvas->cd(icd++);
  fHistnClusters->Draw();
  fTrCanvas->cd(icd++);
  fHistMult->Draw();
  fTrCanvas->cd();

  fTrCanvas->Update();

}

AliEveTrack* AliHLTEveHLT::MakeEsdTrack (AliESDtrack *at, TEveTrackList* cont) {
  //See header file for documentation


  const double kCLight = 0.000299792458;
  double bz = - kCLight*10.*( cont->GetPropagator()->GetMagField(0,0,0).fZ);

  Bool_t innerTaken = kFALSE;
  if ( ! at->IsOn(AliESDtrack::kITSrefit) && fUseIpOnFailedITS)
  {
    //tp = at->GetInnerParam();
    innerTaken = kTRUE;
  }

  // Add inner/outer track parameters as path-marks.

  Double_t     pbuf[3], vbuf[3];

  AliExternalTrackParam trackParam = *at;

  // take parameters constrained to vertex (if they are)

  if( at->GetConstrainedParam() ){
    trackParam = *at->GetConstrainedParam();
  }
  else if( at->GetInnerParam() ){
    trackParam = *(at->GetInnerParam());
  }
  if( at->GetStatus()&AliESDtrack::kTRDin ){
    // transport to TRD in
    trackParam = *at;
    trackParam.PropagateTo( 290.45, -10.*( cont->GetPropagator()->GetMagField(0,0,0).fZ) );
  }

  TEveRecTrack rt;
  {
    rt.fLabel  = at->GetLabel();
    rt.fIndex  = (Int_t) at->GetID();
    rt.fStatus = (Int_t) at->GetStatus();
    rt.fSign   = (Int_t) trackParam.GetSign();  
    trackParam.GetXYZ(vbuf);
    trackParam.GetPxPyPz(pbuf);    
    rt.fV.Set(vbuf);
    rt.fP.Set(pbuf);
    Double_t ep = at->GetP(), mc = at->GetMass();
    rt.fBeta = ep/TMath::Sqrt(ep*ep + mc*mc);
  }

  AliEveTrack* track = new AliEveTrack(&rt, cont->GetPropagator());
  track->SetAttLineAttMarker(cont);
  track->SetName(Form("AliEveTrack %d", at->GetID()));
  track->SetElementTitle(CreateTrackTitle(at));
  track->SetSourceObject(at);


  // Set reference points along the trajectory
  // and the last point

  { 
    TEvePathMark startPoint(TEvePathMark::kReference);
    trackParam.GetXYZ(vbuf);
    trackParam.GetPxPyPz(pbuf);    
    startPoint.fV.Set(vbuf);
    startPoint.fP.Set(pbuf);
    rt.fV.Set(vbuf);
    rt.fP.Set(pbuf);
    Double_t ep = at->GetP(), mc = at->GetMass();
    rt.fBeta = ep/TMath::Sqrt(ep*ep + mc*mc);

    track->AddPathMark( startPoint );    
  }


  if( at->GetTPCPoints(2)>80 ){
  
    //
    // use AliHLTTPCCATrackParam propagator 
    // since AliExternalTrackParam:PropagateTo()
    // has an offset at big distances
    //
    
    AliHLTTPCCATrackParam t;
    AliHLTTPCCATrackConvertor::SetExtParam( t, trackParam );
    
    Double_t x0 = trackParam.GetX();
    Double_t dx = at->GetTPCPoints(2) - x0;
    
    //
    // set a reference at the half of trajectory for better drawing
    //
    
    for( double dxx=dx/2; TMath::Abs(dxx)>=1.; dxx*=.9 ){
      if( !t.TransportToX(x0+dxx, bz, .99 ) ) continue;
      AliHLTTPCCATrackConvertor::GetExtParam( t, trackParam, trackParam.GetAlpha() ); 
      trackParam.GetXYZ(vbuf);
      trackParam.GetPxPyPz(pbuf);
      TEvePathMark midPoint(TEvePathMark::kReference);
      midPoint.fV.Set(vbuf);
      midPoint.fP.Set(pbuf);    
      track->AddPathMark( midPoint );
      break;
    }
    
    //
    // Set a reference at the end of the trajectory
    // and a "decay point", to let the event display know where the track ends
    //
    
    for( ; TMath::Abs(dx)>=1.; dx*=.9 ){
      if( !t.TransportToX(x0+dx, bz, .99 ) ) continue;
      AliHLTTPCCATrackConvertor::GetExtParam( t, trackParam, trackParam.GetAlpha() ); 
      trackParam.GetXYZ(vbuf);
      trackParam.GetPxPyPz(pbuf);
      TEvePathMark endPoint(TEvePathMark::kReference);
      TEvePathMark decPoint(TEvePathMark::kDecay);
      endPoint.fV.Set(vbuf);
      endPoint.fP.Set(pbuf);
      decPoint.fV.Set(vbuf);
      decPoint.fP.Set(pbuf);
      track->AddPathMark( endPoint );
      track->AddPathMark( decPoint );
      break;
    }  
  }

  if (at->IsOn(AliESDtrack::kTPCrefit))
  {
    if ( ! innerTaken)
    {
      AddTrackParamToTrack(track, at->GetInnerParam());
    }
    AddTrackParamToTrack(track, at->GetOuterParam());
  }
  return track;
}

void AliHLTEveHLT::SetUpTrackPropagator(TEveTrackPropagator* trkProp, Float_t magF, Float_t maxR) {
  //See header file for documentation

  if (fTrueField) {
    trkProp->SetMagFieldObj(new AliEveMagField);
  
  } else {
    trkProp->SetMagField(magF);
  }
 
  if (fUseRkStepper) {
    trkProp->SetStepper(TEveTrackPropagator::kRungeKutta);
  }

  trkProp->SetMaxR(maxR);
}


void AliHLTEveHLT::AddTrackParamToTrack(AliEveTrack* track, const AliExternalTrackParam* tp) {
  //See header file for documentation

  if (tp == 0)
    return;

  Double_t pbuf[3], vbuf[3];
  tp->GetXYZ(vbuf);
  tp->GetPxPyPz(pbuf);

  TEvePathMark pm(TEvePathMark::kReference);
  pm.fV.Set(vbuf);
  pm.fP.Set(pbuf);
  track->AddPathMark(pm);
}



TString AliHLTEveHLT::CreateTrackTitle(AliESDtrack* t) {
  // Add additional track parameters as a path-mark to track.

  TString s;

  Int_t label = t->GetLabel(), index = t->GetID();
  TString idx(index == kMinInt ? "<undef>" : Form("%d", index));
  TString lbl(label == kMinInt ? "<undef>" : Form("%d", label));

  Double_t p[3], v[3];
  t->GetXYZ(v);
  t->GetPxPyPz(p);
  Double_t pt    = t->Pt();
  Double_t ptsig = TMath::Sqrt(t->GetSigma1Pt2());
  Double_t ptsq  = pt*pt;
  Double_t ptm   = pt / (1.0 + pt*ptsig);
  Double_t ptM   = pt / (1.0 - pt*ptsig);

  s = Form("Index=%s, Label=%s\nChg=%d, Pdg=%d\n"
	   "pT = %.3f + %.3f - %.3f [%.3f]\n"
           "P  = (%.3f, %.3f, %.3f)\n"
           "V  = (%.3f, %.3f, %.3f)\n",
	   idx.Data(), lbl.Data(), t->Charge(), 0,
	   pt, ptM - pt, pt - ptm, ptsig*ptsq,
           p[0], p[1], p[2],
           v[0], v[1], v[2]);

  Int_t   o;
  s += "Det (in,out,refit,pid):\n";
  o  = AliESDtrack::kITSin;
  s += Form("ITS (%d,%d,%d,%d)  ",  t->IsOn(o), t->IsOn(o<<1), t->IsOn(o<<2), t->IsOn(o<<3));
  o  = AliESDtrack::kTPCin;
  s += Form("TPC(%d,%d,%d,%d)\n",   t->IsOn(o), t->IsOn(o<<1), t->IsOn(o<<2), t->IsOn(o<<3));
  o  = AliESDtrack::kTRDin;
  s += Form("TRD(%d,%d,%d,%d) ",    t->IsOn(o), t->IsOn(o<<1), t->IsOn(o<<2), t->IsOn(o<<3));
  o  = AliESDtrack::kTOFin;
  s += Form("TOF(%d,%d,%d,%d)\n",   t->IsOn(o), t->IsOn(o<<1), t->IsOn(o<<2), t->IsOn(o<<3));
  o  = AliESDtrack::kHMPIDout;
  s += Form("HMPID(out=%d,pid=%d)\n", t->IsOn(o), t->IsOn(o<<1));
  s += Form("ESD pid=%d", t->IsOn(AliESDtrack::kESDpid));

  if (t->IsOn(AliESDtrack::kESDpid))
  {
    Double_t pid[5];
    t->GetESDpid(pid);
    s += Form("\n[%.2f %.2f %.2f %.2f %.2f]", pid[0], pid[1], pid[2], pid[3], pid[4]);
  }

  return s;
}

