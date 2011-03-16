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
#include "AliEveHLTEventManager.h"
#include "AliHLTGlobalTriggerDecision.h"
#include "TEveManager.h"
#include "TEvePointSet.h"
#include "TEveTrack.h"
#include "TEveElement.h"
#include "TCanvas.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "TEveTrackPropagator.h"
#include "AliEveTrack.h"
#include "TEveVSDStructs.h"
#include "TString.h"
#include "TPCLib/tracking-ca/AliHLTTPCCATrackParam.h"
#include "TPCLib/tracking-ca/AliHLTTPCCATrackConvertor.h"
#include "AliEveMagField.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TThread.h"

ClassImp(AliHLTEveHLT)

AliHLTEveHLT::AliHLTEveHLT() : 
  AliHLTEveBase("TPC tracks"), 
  fTrueField(kFALSE),
  fUseIpOnFailedITS(kFALSE),
  fUseRkStepper(kFALSE),
  fTrackList(NULL),
  fTrackLists(NULL),
  fOldTrackList(NULL),
  fPointSetVertex(NULL),
  fTrCanvas(NULL),
  fVertexCanvas(NULL),
  fHistEta(NULL),
  fHistPhi(NULL),
  fHistnClusters(NULL),
  fHistMult(NULL),
  fHistDCAr(NULL),
  fTrCount(0), 
  fVCount(0),
  fNTrackBins(10)
{
  // Constructor.
  CreateHistograms();

}
///___________________________________________________________________
AliHLTEveHLT::~AliHLTEveHLT()
{
  //Destructor, not implemented
  if(fTrackList)
    delete fTrackList;
  fTrackList = NULL;

  if(fTrackLists)
    delete fTrackLists;
  fTrackLists = NULL;


}


///________________________________________________________________________
void AliHLTEveHLT::CreateHistograms(){
  //See header file for documentation

  //fHistPt        = new TH1F("fHistPt",       "transverse momentum",    100, 0, 10); // KK   
  //fHistP         = new TH1F("fHistP",        "signed momentum",        100,-7,  7);	   

  //fHistPt   ->SetXTitle("p_{t} (GeV/c)");   // KK
  //fHistP    ->SetXTitle("P*charge (GeV/c)");

  // fHistTheta     = new TH1F("fHistTheta",    "polar angle",            180, 0,180);   
  // fHistTheta->SetXTitle("#theta (degrees)");
  
  fHistEta       = new TH1F("fHistEta",      "pseudorapidity",         100,-2,  2);	   
  fHistEta  ->SetXTitle("#eta");

  fHistPhi	 = new TH1F("fHistPhi",      "azimuthal angle",        180, 0,360);   
  fHistPhi  ->SetXTitle("#phi (degrees)");

  fHistnClusters = new TH1F("fHistnClusters","TPC clusters per track", 160, 0,160);

  fHistMult      = new TH1F("fHistMult",     "event track multiplicity",150, 0, 15000);    
  
  fHistDCAr = new TH1F("fHistDCAr", "DCA r", 200, -100, 100);

}


///_____________________________________________________________________
void AliHLTEveHLT::ProcessBlock(AliHLTHOMERBlockDesc * block) {
  //See header file for documentation
  if ( ! block->GetDataType().CompareTo("ALIESDV0") ) {
    ProcessEsdBlock(block);
  } 
  
  else if ( ! block->GetDataType().CompareTo("ROOTTOBJ") ) {
    //processROOTTOBJ( block, gHLTText );
  } 

  else if ( ! block->GetDataType().CompareTo("HLTRDLST") ) {
    cout << "ignoring hlt rdlst"<<endl;
    //processHLTRDLST( block );
  } 

  else if ( ! block->GetDataType().CompareTo("GLOBTRIG") ) {
    ProcessGlobalTrigger( block );
  } 

  else if ( !block->GetDataType().CompareTo("ROOTHIST") ) {      
    if( !fCanvas ) { 
      fCanvas = CreateCanvas("PVtx/Tr QA", "PrimV");
      fCanvas->Divide(3, 3);
    }
    ProcessHistograms( block , fCanvas);
  }  
}

///____________________________________________________________________________
void AliHLTEveHLT::UpdateElements() {
  //See header file for documentation

  DrawHistograms();

  if(fTrackLists) {
    for(Int_t il = 0; il < fNTrackBins; il++) {
      TEveTrackList * trackList = dynamic_cast<TEveTrackList*>(fTrackLists->FindChild(Form("Tracks_%d", il)));
      if(trackList) trackList->ElementChanged();
    } 
  }

  if(fPointSetVertex) fPointSetVertex->ResetBBox();
  if(fTrCanvas) fTrCanvas->Update();
  if(fVertexCanvas) fVertexCanvas->Update();
  if(fCanvas) fCanvas->Update();
  
}

///_________________________________________________________________________________
void AliHLTEveHLT::ResetElements(){
  //See header file for documentation
  
  cout << "destroy"<<endl;

  if(fTrackLists) {
    for(Int_t il = 0; il < fNTrackBins; il++) {
      TEveTrackList * trackList = dynamic_cast<TEveTrackList*>(fTrackLists->FindChild(Form("Tracks_%d", il)));
      if(trackList) trackList->DestroyElements();
    } 
  }

  
  fTrCount = 0;
  fVCount = 0;

  if(fPointSetVertex) fPointSetVertex->Reset();
  cout<< "reset done"<<endl;
  fHistoCount = 0;

}

///_____________________________________________________________________________________
void * AliHLTEveHLT::DestroyGarbage(void * arg) {
  AliHLTEveHLT * hlt = reinterpret_cast<AliHLTEveHLT*>(arg);
  if(hlt) hlt->DestroyOldTrackList();
  return (void*)0;
}
///_____________________________________________________________________________________
void AliHLTEveHLT::DestroyOldTrackList() {
  cout << "Destroying the old tracklist's elements"<<endl;
  fOldTrackList->DestroyElements();
  cout << "Destroying the old tracklist itself"<<endl;
  fOldTrackList->Destroy();
}

///_____________________________________________________________________________________
void AliHLTEveHLT::ProcessHistograms(AliHLTHOMERBlockDesc * block, TCanvas * canvas) {
  //See header file for documentation
  if (!fTrCanvas) {
    fTrCanvas = CreateCanvas("ESD", "ESD");
    fTrCanvas->Divide(4, 2);
  }

  if(!fVertexCanvas) {
    fVertexCanvas = CreateCanvas("Vtx", "Vtx");
    fVertexCanvas->Divide(4, 2);
  }



  if ( ! block->GetClassName().CompareTo("TH1F")) {
    TH1F* histo = reinterpret_cast<TH1F*>(block->GetTObject());
    if( histo ){
      TString name(histo->GetName());
      cout << "TH1F " <<name << endl;
      if( !name.CompareTo("primVertexZ") ){
	canvas->cd(2);
	histo->Draw();
      }else if( !name.CompareTo("primVertexX") ){
	canvas->cd(3);
	histo->Draw();
      }else if( !name.CompareTo("primVertexY") ){
	canvas->cd(4);
	histo->Draw();
      } else if ( name.Contains("Track")) {
	AddHistogramToCanvas(histo, fTrCanvas, fTrCount);
      } else {
	AddHistogramToCanvas(histo, fVertexCanvas, fVCount);
      }
    }
  }  else if ( ! block->GetClassName().CompareTo("TH2F")) {
    TH2F *hista = reinterpret_cast<TH2F*>(block->GetTObject());
    if (hista ){
       TString name(hista->GetName());
       cout << "TH2F " << name << endl;
       if( !name.CompareTo("primVertexXY")) {      
	 canvas->cd(1);
	 hista->Draw();
       } else if ( name.Contains("Track")) {
	AddHistogramToCanvas(hista, fTrCanvas, fTrCount);
      } else {
	AddHistogramToCanvas(hista, fVertexCanvas, fVCount);
      }
    }
  }
  canvas->cd();
}

///_____________________________________________________________________________________
void AliHLTEveHLT::AddHistogramToCanvas(TH1* histogram, TCanvas * canvas, Int_t &cdCount) {
  canvas->cd(++cdCount);
  histogram->Draw();
}


///________________________________________________________________________________________
void AliHLTEveHLT::CreateTrackList() {
  //See header file for documentation
  fTrackList = new TEveTrackList("ESD Tracks");
  fTrackList->SetMainColor(kOrange);
  AddElement(fTrackList);
}

///________________________________________________________________________________________
void AliHLTEveHLT::CreateTrackLists() {
  //See header file for documentation
  fTrackLists = new TEveElementList("ESD Track lists");
  AddElement(fTrackLists);
  for(Int_t i = 0; i < 10; i++) {
    TEveTrackList * trackList = new TEveTrackList(Form("Tracks_%d", i));
    trackList->SetMainColor(kOrange-i);
    fTrackLists->AddElement(trackList);
  }
}


///_________________________________________________________________________________________
void AliHLTEveHLT::CreateVertexPointSet() {
  //See header file for documentation
  fPointSetVertex = new TEvePointSet("Primary Vertex");
  fPointSetVertex->SetMainColor(6);
  fPointSetVertex->SetMarkerStyle((Style_t)kFullStar);

  AddElement(fPointSetVertex);
}

///________________________________________________________________________
void AliHLTEveHLT::ProcessGlobalTrigger( AliHLTHOMERBlockDesc * block ) {
  //See header file for documentation
  AliHLTGlobalTriggerDecision * decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(block->GetTObject());
  if(decision) decision->Print();

}


///______________________________________________________________________
void AliHLTEveHLT::ProcessEsdBlock( AliHLTHOMERBlockDesc * block) {
  //See header file for documentation

  AliESDEvent* esd = (AliESDEvent *) (block->GetTObject());
  if (!esd) return;
  
  ProcessEsdEvent(esd);
}



///_________________________________________________________________________________
Color_t AliHLTEveHLT::GetColor(Float_t pt) {
  //See header file
  Color_t baseColor = kOrange;
  
  Float_t binlimit = 0.1;
  for(Int_t i = 0; i< 10; i++) {
   
    if (pt < binlimit)
      return baseColor - i;
    
    binlimit +=0.1;
  }
  
  return baseColor - 9;

}

///_________________________________________________________________________________
Int_t AliHLTEveHLT::GetColorBin(Float_t pt) {
  //See header file
  
  Float_t binlimit = 0.1;
  for(Int_t i = 0; i< 10; i++) {
   
    if (pt < binlimit)
      return i;
    
    binlimit +=0.1;
  }
  
  return 9;
  
}

///____________________________________________________________________
void AliHLTEveHLT::ProcessEsdEvent( AliESDEvent * esd ) {
  //See header file for documentation
  if(!fPointSetVertex) CreateVertexPointSet();
  if(!fTrackLists) CreateTrackLists();

  cout << "ProcessESDEvent() :"<< esd->GetEventNumberInFile()<< "  " << esd->GetNumberOfCaloClusters() << " tracks : " << esd->GetNumberOfTracks() << endl;

  //fEventManager->SetRunNumber(esd->GetRunNumber());

  Double_t vertex[3];
  const AliESDVertex * esdVertex = esd->GetPrimaryVertex();
  
  if(esdVertex) {
    esdVertex->GetXYZ(vertex);
    fPointSetVertex->SetNextPoint(vertex[0], vertex[1], vertex[2]);
  }
  
  TEveTrackList * trackList = dynamic_cast<TEveTrackList*>(fTrackLists->FirstChild());
  if(trackList) SetUpTrackPropagator(trackList->GetPropagator(),-0.1*esd->GetMagneticField(), 520);
  for (Int_t iter = 0; iter < esd->GetNumberOfTracks(); ++iter) {

   AliESDtrack * esdTrack = dynamic_cast<AliESDtrack*>(esd->GetTrack(iter));
   FillTrackList(esdTrack);
   FillHistograms(esdTrack);
   
  }
  fHistMult->Fill(esd->GetNumberOfTracks()); // KK
  
  //BALLE hardcoded size
  for(Int_t il = 0; il < fNTrackBins; il++) {
    trackList = dynamic_cast<TEveTrackList*>(fTrackLists->FindChild(Form("Tracks_%d", il)));
    trackList->MakeTracks();
  } 
}
///__________________________________________________________________________
void AliHLTEveHLT::FillTrackList(AliESDtrack * esdTrack) {
  //See header file for documentation
  TEveTrackList * trackList = dynamic_cast<TEveTrackList*>(fTrackLists->FirstChild());
  if (!trackList) return;
  
  AliEveTrack* track = dynamic_cast<AliEveTrack*>(MakeEsdTrack(esdTrack, trackList));        
  Int_t bin = GetColorBin(esdTrack->Pt());
  trackList = dynamic_cast<TEveTrackList*>(fTrackLists->FindChild(Form("Tracks_%d", bin)));
  if(trackList) {
    track->SetAttLineAttMarker(trackList);
    trackList->AddElement(track);
  } else cout << "BALLE"<<endl;


}


///____________________________________________________________________________________
void AliHLTEveHLT::FillHistograms(AliESDtrack * esdTrack) {

  if(esdTrack->GetTPCNcls() == 0) return;
  
  fHistEta->Fill(esdTrack->Eta());
  // fHistTheta->Fill(esdTrack->Theta()*TMath::RadToDeg());
  fHistPhi->Fill(esdTrack->Phi()*TMath::RadToDeg());
  
  
  Float_t DCAr, DCAz = -99;
  esdTrack->GetImpactParametersTPC(DCAr, DCAz);
  fHistDCAr->Fill(DCAr);
  
  
  if(esdTrack->GetStatus()&AliESDtrack::kTPCin || (esdTrack->GetStatus()&AliESDtrack::kTPCin && esdTrack->GetStatus()&AliESDtrack::kITSin)){
    fHistnClusters->Fill(esdTrack->GetTPCNcls());  
  }
}

void AliHLTEveHLT::DrawHistograms(){
  //See header file for documentation
  if(!fCanvas) {
    fCanvas = CreateCanvas("PVtx/Tr QA", "Primary vertex, Track QA");
    fCanvas->Divide(3, 3);
  }

  fCanvas->cd(5);
  fHistEta->Draw();

  fCanvas->cd(6);
  fHistPhi->Draw();

  fCanvas->cd(7);
  fHistnClusters->Draw();

  fCanvas->cd(8);
  fHistMult->Draw();

  fCanvas->cd(9);
  fHistDCAr->Draw();

  fCanvas->cd();

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
    float x = at->GetX();
    float y = at->GetY();
    if(sqrt(x*x+y*y)>.5 ) trackParam = *(at->GetInnerParam());
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

  bool ok = 1;
  if( at->GetOuterParam() && at->GetTPCPoints(2)>80 ){

    //
    // use AliHLTTPCCATrackParam propagator 
    // since AliExternalTrackParam:PropagateTo()
    // has an offset at big distances
    //
    double rot = at->GetOuterParam()->GetAlpha() - trackParam.GetAlpha();
    double crot = cos(rot), srot = sin(rot);
    double xEnd = at->GetTPCPoints(2)*crot -  at->GetTPCPoints(3)*srot;
  // take parameters constrained to vertex (if they are)
 

    AliHLTTPCCATrackParam t;
    AliHLTTPCCATrackConvertor::SetExtParam( t, trackParam );
    
    Double_t x0 = trackParam.GetX();
    Double_t dx = xEnd - x0;
    
    if( dx<0 ) ok = 0;
    //
    // set a reference at the half of trajectory for better drawing
    //
    
    if( ok ){ 
      double dxx=dx/2; 
      if( TMath::Abs(dxx)>=1. ){

	if( !t.TransportToX(x0+dxx, bz, .999 ) ){
	  ok = 0;
	} else {
	  AliExternalTrackParam tt;
	  AliHLTTPCCATrackConvertor::GetExtParam( t, tt, trackParam.GetAlpha() ); 
	  tt.GetXYZ(vbuf);
	  tt.GetPxPyPz(pbuf);
	  TEvePathMark midPoint(TEvePathMark::kReference);
	  midPoint.fV.Set(vbuf);
	  midPoint.fP.Set(pbuf);    
	  track->AddPathMark( midPoint );     
	}
      }
    }
   
    //
    // Set a reference at the end of the trajectory
    // and a "decay point", to let the event display know where the track ends
    //
    
    for( ; ok; dx*=.9 ){

      if( !t.TransportToX(x0+dx, bz, .999 ) ){
	ok = 0; 
	break;
	if( TMath::Abs(dx)<1. ) break;      
	continue;
      }
      break;
    }

    {
      if( !ok ){ 
	AliHLTTPCCATrackConvertor::SetExtParam( t, trackParam );
      }
      AliExternalTrackParam tt;
      AliHLTTPCCATrackConvertor::GetExtParam( t, tt, trackParam.GetAlpha() ); 
      tt.GetXYZ(vbuf);
      tt.GetPxPyPz(pbuf);
      TEvePathMark endPoint(TEvePathMark::kReference);
      TEvePathMark decPoint(TEvePathMark::kDecay);
      endPoint.fV.Set(vbuf);
      endPoint.fP.Set(pbuf);
      decPoint.fV.Set(vbuf);
      decPoint.fP.Set(pbuf);
      track->AddPathMark( endPoint );
      track->AddPathMark( decPoint );
    }  
  }

  if ( ok &&  at->IsOn(AliESDtrack::kTPCrefit))
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


