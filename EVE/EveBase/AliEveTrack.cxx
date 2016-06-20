// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTrack.h"

#include "AliESDtrack.h"
#include "AliAODTrack.h"

#include <TROOT.h>
#include <TMath.h>

#include <TEveUtil.h>
#include <TEvePointSet.h>
#include <TEveElement.h>
#include <TEveManager.h>

#include <TParticle.h>

#include <AliRunLoader.h>
#include <AliStack.h>
#include <AliEveEventManager.h>

#include "its_hits.C"
#include "tof_hits.C"
#include "tpc_hits.C"
#include "trd_hits.C"

//______________________________________________________________________________
// Full description of AliEveTrack
//

ClassImp(AliEveTrack)

//______________________________________________________________________________
AliEveTrack::AliEveTrack() :
TEveTrack()
{
    // Constructor.
}

//______________________________________________________________________________
AliEveTrack::AliEveTrack(TParticle* t, Int_t label, TEveTrackPropagator* prop) :
TEveTrack(t, label, prop)
{
    // Constructor.
}

//______________________________________________________________________________
AliEveTrack::AliEveTrack(TEveMCTrack*  t, TEveTrackPropagator* prop) :
TEveTrack(t, prop)
{
    // Constructor.
}

//______________________________________________________________________________
AliEveTrack::AliEveTrack(TEveRecTrack* t, TEveTrackPropagator* prop) :
TEveTrack(t, prop)
{
    // Constructor.
}

//______________________________________________________________________________
AliEveTrack::AliEveTrack(AliESDtrack* t, TEveTrackPropagator* prop) :
TEveTrack()
{
    // Constructor.
    
    Double_t buf[3];
    t->GetXYZ(buf);    fV.Set(buf);
    t->GetPxPyPz(buf); fP.Set(buf);
    
    Double_t ep = t->GetP(), mc = t->GetMass();
    fBeta = ep/TMath::Sqrt(ep*ep + mc*mc);
    // fPdg = 0; // ??? Use PID ?
    fCharge = TMath::Nint(t->GetSign());
    
    fLabel = t->GetLabel();
    fIndex = t->GetID();
    // fStatus = (Int_t) t->GetStatus(); // RRRR Uncomment for root-5.26.
    
    SetPropagator(prop);
}

//______________________________________________________________________________
AliEveTrack::AliEveTrack(AliAODTrack* t, TEveTrackPropagator* prop) :
TEveTrack()
{
    // Constructor.
    
    Double_t buf[3];
    
    t->GetXYZ(buf); fV.Set(buf);
    t->PxPyPz(buf); fP.Set(buf);
    
    // fBeta = 0; // Unknown, no mass function
    // fPdg = 0;  // ??? Use PID ?
    fCharge= t->Charge();
    
    fLabel = t->GetLabel();
    fIndex = t->GetID();
    // fStatus = (Int_t) t->GetStatus(); // RRRR Uncomment for root-5.26.
    
    SetPropagator(prop);
}

//______________________________________________________________________________
AliEveTrack::AliEveTrack(const AliEveTrack& t) :
TEveTrack(t)
{
    // Copy constructor.
}

//______________________________________________________________________________
AliEveTrack::~AliEveTrack()
{
    // Destructor.
}

//______________________________________________________________________________
void AliEveTrack::SetStartParams(const AliExternalTrackParam* tp)
{
    // Set the initial vertex / momentum of eve track from 'tp'.
    
    Double_t buf[3];
    
    tp->GetXYZ(buf); fV.Set(buf);
    tp->PxPyPz(buf); fP.Set(buf);
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveTrack::ImportHits()
{
    // Import hits with same label as the track.
    
    TEveUtil::LoadMacro("its_hits.C");
    
    TEvePointSet* h = 0;
    //PH The line below is replaced waiting for a fix in Root
    //PH which permits to use variable siza arguments in CINT
    //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
    //PH  h = its_hits("fX:fY:fZ", Form("ITS.fTrack==%d", label));
    char form[1000];
    sprintf(form,"ITS.fTrack==%d", fLabel);
    h = its_hits("fX:fY:fZ", form, this);
    if (h) h->SetMarkerSize(1);
    
    TEveUtil::LoadMacro("tpc_hits.C");
    sprintf(form,"TPC2.fArray.fTrackID==%d", fLabel);
    h = tpc_hits("TPC2.fArray.fR:TPC2.fArray.fFi:TPC2.fArray.fZ",form, this);
    //PH  h = tpc_hits("TPC2.fArray.fR:TPC2.fArray.fFi:TPC2.fArray.fZ",
    //PH	       Form("TPC2.fArray.fTrackID==%d", label));
    if (h) h->SetMarkerSize(1);
    
    TEveUtil::LoadMacro("trd_hits.C");
    sprintf(form,"TRD.fTrack==%d", fLabel);
    h = trd_hits("fX:fY:fZ", form, this);
    if (h) h->SetMarkerSize(1);
    
    TEveUtil::LoadMacro("tof_hits.C");
    sprintf(form,"TOF.fTrack==%d", fLabel);
    h = tof_hits("fX:fY:fZ", form, this);
    if (h) h->SetMarkerSize(1);
    
    gEve->Redraw3D();
    
    
}

//______________________________________________________________________________
void AliEveTrack::ImportClustersFromLabel()
{
    // Import clusters with same label as the track.
    
    AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();
    TEvePointSet* clusters = new TEvePointSet(64);
    clusters->SetOwnIds(kTRUE);
    
    for (Int_t n=0; n<esd->GetNumberOfTracks(); n++)
    {
        AliESDtrack* at = esd->GetTrack(n);
        if (at->GetLabel() == fLabel) {
            const AliTrackPointArray* pArr = at->GetTrackPointArray();
            if (pArr == 0) {
                Warning("clusters_from_label", "TrackPointArray not stored with ESD track.");
                continue;
            }
            Int_t np =  pArr->GetNPoints();
            const Float_t* x = pArr->GetX();
            const Float_t* y = pArr->GetY();
            const Float_t* z = pArr->GetZ();
            for (Int_t i=0; i<np; ++i) {
                clusters->SetNextPoint(x[i], y[i], z[i]);
                AliTrackPoint *atp = new AliTrackPoint;
                pArr->GetPoint(*atp, i);
                clusters->SetPointId(atp);
            }
        }
    }
    
    if(clusters->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
        Warning("clusters_from_label", "No clusters match label '%d'", fLabel);
        delete clusters;
        return;
    }
    
    clusters->SetMarkerStyle(2);
    clusters->SetMarkerSize(0.5);
    clusters->SetMarkerColor(4);

    char form[1000];
    sprintf(form,"Clusters lab=%d", fLabel);
    clusters->SetName(form);
    char tip[1000];
    sprintf(tip,"N=%d", clusters->Size());
    clusters->SetTitle(tip);
    gEve->AddElement(clusters, this);
    gEve->Redraw3D();
    
    return;
}

//______________________________________________________________________________
void AliEveTrack::ImportClustersFromIndex()
{
    // Import clusters marked with same reconstructed track index as the track.
    static const TEveException kEH("AliEveTrack::ImportClustersFromIndex ");
    
    if (fIndex == kMinInt)
        throw kEH + "index not set.";
    
    ImportClustersFromIndex(fIndex);
}

TEvePointSet* AliEveTrack::ImportClustersFromIndex(Int_t index)
{
    AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();
    
    if (index < 0) {
        Warning("AliEveTrack::ImportClustersFromIndex", "index not set.");
        return 0;
    }
    
    if (index >= esd->GetNumberOfTracks()) {
        Warning("AliEveTrack::ImportClustersFromIndex", "index out of range");
        return 0;
    }
    
    TEvePointSet* clusters = new TEvePointSet(64);
    clusters->SetOwnIds(kTRUE);
    
    AliESDtrack* at = esd->GetTrack(index);
    const AliTrackPointArray* pArr = at->GetTrackPointArray();
    if (pArr == 0) {
        Warning("AliEveTrack::ImportClustersFromIndex", "TrackPointArray not stored with ESD track.");
    }
    
    Int_t np =  pArr->GetNPoints();
    const Float_t* x = pArr->GetX();
    const Float_t* y = pArr->GetY();
    const Float_t* z = pArr->GetZ();
    for (Int_t i=0; i<np; ++i) {
        clusters->SetNextPoint(x[i], y[i], z[i]);
        AliTrackPoint *atp = new AliTrackPoint;
        pArr->GetPoint(*atp, i);
        clusters->SetPointId(atp);    }
    
    
    if(clusters->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE) {
        Warning("AliEveTrack::ImportClustersFromIndex", "No clusters for index '%d'", index);
        delete clusters;
        return 0;
    }
    
    clusters->SetMarkerStyle(2);
    clusters->SetMarkerSize(2);
    clusters->SetMarkerColor(4);
    
    clusters->SetName(Form("Clusters idx=%d", index));
    clusters->SetTitle(Form("N=%d", clusters->Size()));
    
    gEve->AddElement(clusters);
    
    if (AliEveMultiView::Instance())
    {
        AliEveMultiView::Instance()->ImportEventRPhi(clusters);
        AliEveMultiView::Instance()->ImportEventRhoZ(clusters);
    }
    
    gEve->Redraw3D();
    
    return clusters;}

/******************************************************************************/

//______________________________________________________________________________
void AliEveTrack::ImportKine()
{
    // Import kinematics of the track's label recursively.
    // Uses macro "kine_tracks.C".
    
    static const TEveException kEH("AliEveTrack::ImportKine ");
    
    if (fLabel == kMinInt)
        throw kEH + "label not set.";
    
    Int_t label;
    if (fLabel < 0) {
        Warning(kEH, "label negative, taking absolute value.");
        label = -fLabel;
    } else {
        label = fLabel;
    }
    
    TEveUtil::LoadMacro("kine_tracks.C");
    gROOT->ProcessLine(Form("kine_track(%d, kTRUE, kTRUE, kTRUE, kTRUE, (TEveElement*)%p);",
                            label, this));
    
}

//______________________________________________________________________________
void AliEveTrack::ImportKineWithArgs(Bool_t importMother, Bool_t importDaugters,
                                     Bool_t colorPdg,     Bool_t recurse)
{
    // Import kinematics of the track's label. Arguments steer the
    // import process:
    //   importMother     import particle with track's label
    //   importDaugters   import direct daughters of label
    //   colorPdg         color kinematics by PDG code
    //   recurse          recursive import of daughters' daughters
    // Uses macro "kine_tracks.C".
    
    static const TEveException kEH("AliEveTrack::ImportKineWithArgs ");
    
    if (fLabel == kMinInt)
        throw kEH + "label not set.";
    
    Int_t label;
    if (fLabel < 0) {
        Warning(kEH, "label negative, taking absolute value.");
        label = -fLabel;
    } else {
        label = fLabel;
    }
    
    TEveUtil::LoadMacro("kine_tracks.C");
    gROOT->ProcessLine(Form("kine_track(%d, %d, %d, %d, %d, (TEveElement*)%p);",
                            label, importMother, importDaugters, colorPdg, recurse, this));
}

//______________________________________________________________________________
void AliEveTrack::PrintKineStack()
{
    // Print kinematics pertaining to track's label.
    // Uses macro "print_kine_from_label.C".
    
    static const TEveException kEH("AliEveTrack::PrintKineStack ");
    
    if (fLabel == kMinInt)
        throw kEH + "label not set.";
    
    Int_t label;
    if (fLabel < 0) {
        Warning(kEH, "label negative, taking absolute value.");
        label = -fLabel;
    } else {
        label = fLabel;
    }
    
    AliRunLoader* rl = AliEveEventManager::AssertRunLoader();
    rl->LoadKinematics();
    AliStack* stack = rl->Stack();
    
    printf("Number primaries %d, all particles %d, label %d\n",
           stack->GetNprimary(), stack->GetNtrack(), label);
    if (label < 0 || label >= stack->GetNtrack()) {
        printf("  Label exceeds available range.\n");
        return;
    }
    
    TParticle* part = stack->Particle(label);
    if(part != 0) {
        part->Print();
        while(part->GetMother(0) >= 0) {
            part = stack->Particle(part->GetMother(0));
            part->Print();
        }
    }
}

//______________________________________________________________________________
void AliEveTrack::SecSelected(TEveTrack* track)
{
    // Emits "SecSelected(TEveTrack*)" signal.
    // Called from TEveTrackGL on secondary-selection.
    
    Emit("SecSelected(TEveTrack*)", (Long_t)track);
    SecSelectedTrack((AliEveTrack*) track);
}

//______________________________________________________________________________
void AliEveTrack::SecSelectedTrack(AliEveTrack* track)
{
    // Emits "SecSelectedTrack(AliEveTrack*)" signal.
    
    Emit("SecSelectedTrack(AliEveTrack*)", (Long_t)track);
}

//______________________________________________________________________________
AliESDtrack* AliEveTrack::GetESDTrack() const
{
    // Return source object dyn-casted to AliESDtrack.
    
    return dynamic_cast<AliESDtrack*>(GetSourceObject());
}

//______________________________________________________________________________
AliAODTrack* AliEveTrack::GetAODTrack() const
{
    // Return source object dyn-casted to AliAODTrack.
    
    return dynamic_cast<AliAODTrack*>(GetSourceObject());
}
