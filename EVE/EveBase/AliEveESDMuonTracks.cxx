// $Id$

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

/// \ingroup evemacros
/// \file esd_muon_tracks.C
/// \brief Macro to visualise ESD tracks from MUON spectrometer
/// (both tracker and trigger).
///
/// Use esd_muon_tracks(Bool_t showClusters, Bool_t showDigits) in order to run it
///
/// Needs that alieve_init() is already called
///
/// \author P. Pillot, L. Aphecetche; Subatech

#include <AliEveESDMuonTracks.h>

#include <TStyle.h>
#include <TROOT.h>
#include <TEveManager.h>
#include <TEveUtil.h>
#include <TEveTrack.h>
#include <TEvePointSet.h>
#include <TEveQuadSet.h>
#include <TEveVSDStructs.h>


#include <AliMUONTrack.h>
#include <AliMUONTrackExtrap.h>
#include <AliMUONTrackParam.h>
#include <AliMUONConstants.h>
#include <AliMUONCDB.h>
#include <AliMUONGeometryTransformer.h>
#include <AliMUONTriggerCircuit.h>
#include <AliMpCDB.h>

#include <AliESDMuonTrack.h>
#include <AliEveMagField.h>

#include <AliEveEventManager.h>

//______________________________________________________________________________
void AliEveESDMuonTracks::SetupTrackPropagator(TEveTrackPropagator* trkProp, Bool_t tracker, Bool_t trigger)
{
    // set magnetic field
    if (AliMUONTrackExtrap::IsFieldON())
    {
        trkProp->SetMagFieldObj(new AliEveMagField);
    }
    else
    {
        trkProp->SetMagField(0.0);
    }
    trkProp->SetStepper(TEveTrackPropagator::kRungeKutta);
    
    // set propagation range
    trkProp->SetMaxR(1000);
    if (trigger) trkProp->SetMaxZ(-AliMUONConstants::DefaultChamberZ(13)+10.);
    else trkProp->SetMaxZ(-AliMUONConstants::MuonFilterZBeg());
    
    // go through pathmarks
    trkProp->SetFitDaughters(kFALSE);
    trkProp->SetFitReferences(kTRUE);
    trkProp->SetFitDecay(kFALSE);
    trkProp->SetFitCluster2Ds(kFALSE);
    
    // Render the ref pathmarks
    trkProp->SetRnrReferences(kTRUE);
    trkProp->RefPMAtt().SetMarkerSize(0.5);
    if (trigger) trkProp->RefPMAtt().SetMarkerColor(kGreen);
    else trkProp->RefPMAtt().SetMarkerColor(kAzure);
    
    // Render first vertex
    if (tracker)
    {
        trkProp->SetRnrFV(kTRUE);
        if (trigger) trkProp->RefFVAtt().SetMarkerColor(kGreen);
        else trkProp->RefFVAtt().SetMarkerColor(kAzure);
    }
}

//______________________________________________________________________________
void AliEveESDMuonTracks::AddMuonTracks(AliESDEvent* esd, AliMUONESDInterface* data,TEveTrackList* match, TEveTrackList* nomatch, TEveTrackList* ghost)
{
    // load trigger circuit
    static AliMUONTriggerCircuit* gTriggerCircuit = 0x0;
    if (!gTriggerCircuit)
    {
        AliEveEventManager::Instance()->AssertGeometry();
        AliMUONGeometryTransformer* fMUONGeometryTransformer = new AliMUONGeometryTransformer();
        fMUONGeometryTransformer->LoadGeometryData();
        gTriggerCircuit = new AliMUONTriggerCircuit(fMUONGeometryTransformer);
    }
    
    Int_t nTrack(esd->GetNumberOfMuonTracks());
    TEveRecTrack recTrack;
    TEveTrack* track;
    
    // add ESD tracks to the proper list
    for (Int_t n = 0; n < nTrack; ++n)
    {
        AliESDMuonTrack* emt = esd->GetMuonTrack(n);
        
        // fill general info
        UInt_t trackId = emt->GetUniqueID();
        recTrack.fLabel = emt->GetLabel();
        recTrack.fIndex = (Int_t)trackId;
        
        // fill tracker track specific info
        if ( emt->ContainTrackerData() )
        {
            recTrack.fStatus = emt->GetMatchTrigger();
            recTrack.fSign = emt->Charge();
            recTrack.fV.Set(emt->GetNonBendingCoorAtDCA(),emt->GetBendingCoorAtDCA(),emt->GetZ());
            recTrack.fP.Set(emt->PxAtDCA(),emt->PyAtDCA(),emt->PzAtDCA());
            recTrack.fBeta = ( emt->E() > 0 ) ? emt->P()/emt->E() : 0;
            
            // get proper track list
            TEveTrackList* trackList = nomatch;
            if ( emt->GetMatchTrigger() > 0 ) trackList = match;
            
            // produce eve track
            track = new AliEveTrack(&recTrack,trackList->GetPropagator());
            track->SetName(Form("%cmu",emt->Charge()>0 ? '+':'-'));
            track->SetStdTitle();
            track->SetSourceObject(emt); // WARNING: Change the UniqueID of the object!!
            
            // add path mark
            TIter next(data->FindTrack(trackId)->GetTrackParamAtCluster());
            AliMUONTrackParam* param;
            while ( ( param = static_cast<AliMUONTrackParam*>(next()) ) )
            {
                TEveVector v(param->GetNonBendingCoor(),param->GetBendingCoor(),param->GetZ());
                TEveVector p(param->Px(),param->Py(),param->Pz());
                track->AddPathMark(TEvePathMark(TEvePathMark::kReference,v,p));
            }
            
            // add trigger track if any
            if (emt->ContainTriggerData())
            {
                Double_t x11 = gTriggerCircuit->GetX11Pos(emt->LoCircuit(), emt->LoStripY());
                Double_t y11 = gTriggerCircuit->GetY11Pos(emt->LoCircuit(), emt->LoStripX());
                Double_t z11 = gTriggerCircuit->GetZ11Pos(emt->LoCircuit(), emt->LoStripX());
                Double_t y21 = gTriggerCircuit->GetY21Pos(emt->LoCircuit(), emt->LoStripX()+emt->LoDev()+1);
                Double_t z21 = gTriggerCircuit->GetZ21Pos(emt->LoCircuit(), emt->LoStripX()+emt->LoDev()+1);
                Double_t pz  = -emt->PUncorrected(); // max value
                TEveVector v(x11, y11, z11);
                TEveVector p(pz*x11/z11, pz*(y21-y11)/(z21-z11), pz);
                track->AddPathMark(TEvePathMark(TEvePathMark::kReference,v,p));
            }
            
            // add the track to proper list
            track->SetAttLineAttMarker(trackList);
            trackList->AddElement(track);
        }
        else // fill ghost track specific info
        {
            recTrack.fStatus = 0;
            recTrack.fSign = emt->Charge();
            Double_t z11 = (emt->GetZUncorrected() < -1.) ? emt->GetZUncorrected() : (Double_t)AliMUONConstants::DefaultChamberZ(10);
            recTrack.fV.Set(emt->GetNonBendingCoorUncorrected(),emt->GetBendingCoorUncorrected(),z11);
            recTrack.fP.Set(-TMath::Tan(emt->GetThetaXUncorrected()),-TMath::Tan(emt->GetThetaYUncorrected()),-1.);
            
            // produce eve track
            track = new AliEveTrack(&recTrack,ghost->GetPropagator());
            track->SetName("mu");
            track->SetTitle("Trigger only");
            track->SetSourceObject(emt);
            
            // add the track to proper list
            track->SetAttLineAttMarker(ghost);
            ghost->AddElement(track);
        }
        
    }
    
}

//______________________________________________________________________________
void AliEveESDMuonTracks::Draw(Bool_t showClusters, Bool_t showDigits)
{    
    // load ESD
    AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();
    if (esd->GetNumberOfMuonTracks() == 0 && !gEve->GetKeepEmptyCont()) return;
    
    // load field
    AliEveEventManager::AssertMagField();
    if (!AliMUONESDInterface::GetTracker()) AliMUONESDInterface::ResetTracker(AliMUONCDB::LoadRecoParam());
    
    // load mapping
    AliMpCDB::LoadAll(kFALSE);
    
    // convert ESD objects to MUON objects
    AliMUONESDInterface data;
    data.LoadEvent(*esd);
    
    // track containers
    TEveElementList* trackCont = new TEveElementList("ESD MUON Tracks");
    trackCont->SetTitle(Form("N=%d", esd->GetNumberOfMuonTracks()));
    
    TEveTrackList* match = new TEveTrackList("Matched");
    match->SetRnrPoints(kFALSE);
    match->SetRnrLine(kTRUE);
    match->SetLineColor(kGreen);
    SetupTrackPropagator(match->GetPropagator(), kTRUE, kTRUE);
    trackCont->AddElement(match);
    
    TEveTrackList* nomatch = new TEveTrackList("Not matched");
    nomatch->SetRnrPoints(kFALSE);
    nomatch->SetRnrLine(kTRUE);
    nomatch->SetLineColor(kGreen);
    SetupTrackPropagator(nomatch->GetPropagator(), kTRUE, kFALSE);
    trackCont->AddElement(nomatch);
    
    TEveTrackList* ghost = new TEveTrackList("Ghost");
    ghost->SetRnrPoints(kFALSE);
    ghost->SetRnrLine(kTRUE);
    ghost->SetLineColor(kGreen);
    SetupTrackPropagator(ghost->GetPropagator(), kFALSE, kTRUE);
    trackCont->AddElement(ghost);
    
    // cluster container
    TEvePointSet* clusterList = 0x0;
    if (showClusters && (data.GetNClusters() > 0 || gEve->GetKeepEmptyCont()))
    {
        clusterList = new TEvePointSet(10000);
        clusterList->SetName("ESD MUON Clusters");
        clusterList->SetTitle(Form("N=%d",data.GetNClusters()));
        clusterList->SetPickable(kFALSE);
        clusterList->SetMarkerStyle(5);
        clusterList->SetMarkerColor(kYellow);
        clusterList->SetMarkerSize(2.5);
    }
    
    // digit containers
    TEveElementList* digitCont = 0x0;
    TEveQuadSet* bending = 0x0;
    TEveQuadSet* nonBending = 0x0;
    if (showDigits && (data.GetNDigits() > 0 || gEve->GetKeepEmptyCont()))
    {
        digitCont = new TEveElementList("ESD MUON Digits");
        digitCont->SetTitle(Form("N=%d",data.GetNDigits()));
        
        bending = new TEveQuadSet(TEveQuadSet::kQT_RectangleXY, kFALSE, 32);
        bending->SetName("Bending");
        bending->SetRenderMode(TEveDigitSet::kRM_Fill);
        bending->SetPickable(kFALSE);
        digitCont->AddElement(bending);
        
        nonBending = new TEveQuadSet(TEveQuadSet::kQT_RectangleXY, kFALSE, 32);
        nonBending->SetName("Non bending");
        nonBending->SetRenderMode(TEveDigitSet::kRM_Line);
        nonBending->SetPickable(kFALSE);
        digitCont->AddElement(nonBending);
    }
    
    // add tracks to the proper list and propagate them
    AddMuonTracks(esd, &data, match, nomatch, ghost);
    match->SetTitle(Form("N=%d",match->NumChildren()));
    nomatch->SetTitle(Form("N=%d",nomatch->NumChildren()));
    ghost->SetTitle(Form("N=%d",ghost->NumChildren()));
    match->MakeTracks();
    nomatch->MakeTracks();
    ghost->MakeTracks();
    
    // add cluster to the container
    if (clusterList)
    {
        TEveUtil::LoadMacro("muon_clusters.C+");
        TIter next(data.CreateClusterIterator());
        gROOT->ProcessLine(Form("add_muon_clusters((TIter*)%p, (TEvePointSet*)%p);",&next, clusterList));
    }
    
    // add digits to the containers
    if (digitCont)
    {
        TEveUtil::LoadMacro("muon_digits.C+");
        TIter next(data.CreateDigitIterator());
        gROOT->ProcessLine(Form("add_muon_digits((TIter*)%p, (TEveQuadSet*)%p, (TEveQuadSet*)%p, kFALSE);",
                                &next, bending, nonBending));
        
        // set containers' title
        bending->SetTitle(Form("N=%d",bending->GetPlex()->Size()));
        nonBending->SetTitle(Form("N=%d",nonBending->GetPlex()->Size()));
        
        // automatic scaling
        gStyle->SetPalette(1);
        bending->AssertPalette();
        nonBending->AssertPalette();
    }
    
    // add graphic containers
    gEve->DisableRedraw();
    gEve->AddElement(trackCont);
    if (clusterList) gEve->AddElement(clusterList);
    if (digitCont) gEve->AddElement(digitCont);
    gEve->EnableRedraw();
    gEve->Redraw3D();
}
