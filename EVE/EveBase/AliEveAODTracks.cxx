//
//  AliEveAODTracks.cxx
//
//  Created by Jeremi Niedziela on 17/11/15.
//
//

#include "AliEveAODTracks.h"

#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TGListTree.h>


#include <AliAODEvent.h>
#include <AliExternalTrackParam.h>
#include <AliEveTrackCounter.h>
#include <AliEveMagField.h>
#include <AliEveEventManager.h>
#include <AliEveInit.h>

#include <iostream>

using namespace std;

AliEveAODTracks::AliEveAODTracks() :
fUseIPonFailedITSrefit(kFALSE),
fTrueField(kTRUE),
fRKstepper(kFALSE),
fDashNoRefit(true),
fDrawNoRefit(false),
fWidth(1)
{
    // default color scheme by category:
    fColorsByCategory[0] = kMagenta;
    fColorsByCategory[1] = kMagenta+1;
    fColorsByCategory[2] = kMagenta+2;
    fColorsByCategory[3] = kRed;
    fColorsByCategory[4] = kRed+1;
    fColorsByCategory[5] = kRed+2;
    fColorsByCategory[6] = kGreen;
    fColorsByCategory[7] = kGreen+1;
    fColorsByCategory[8] = kGreen+2;
}

AliEveAODTracks::~AliEveAODTracks()
{
    
}


void AliEveAODTracks::SetupPropagator(TEveTrackPropagator* trkProp,Float_t magF, Float_t maxR)
{
    if (fTrueField)
    {
        trkProp->SetMagFieldObj(new AliEveMagField);
    }
    else
    {
        trkProp->SetMagField(magF);
    }
    if (fRKstepper)
    {
        trkProp->SetStepper(TEveTrackPropagator::kRungeKutta);
    }
    trkProp->SetMaxR(maxR);
}

TString AliEveAODTracks::GetTitle(AliAODTrack* t)
{
    TString s("");
    /*
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
             idx.Data(), lbl.Data(), t->Charge(), AliPID::ParticleCode(t->GetPID()),
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
    */
    return s;
}

void AliEveAODTracks::AddParam(AliEveTrack* track, const AliExternalTrackParam* tp)
{
    // Add additional track parameters as a path-mark to track.
    
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



AliEveTrack* AliEveAODTracks::MakeTrack(AliAODTrack *at, TEveTrackList* cont)
{
    // Make a standard track representation and put it into given container.
    
    // Choose which parameters to use a track's starting point.
    // If gkFixFailedITSExtr is TRUE (FALSE by default) and
    // if ITS refit failed, take track parameters at inner TPC radius.
    
    Bool_t innerTaken = kFALSE;
    if ( ! at->IsOn(AliESDtrack::kITSrefit) && fUseIPonFailedITSrefit)
    {
        innerTaken = kTRUE;
    }
    
    AliEveTrack* track = new AliEveTrack(at, cont->GetPropagator());
    track->SetAttLineAttMarker(cont);
    track->SetName(Form("AliEveTrack %d", at->GetID()));
    track->SetElementTitle(GetTitle(at));
    track->SetSourceObject(at);
    
    // Add inner/outer track parameters as path-marks.
    if (at->IsOn(AliESDtrack::kTPCrefit))
    {
        if ( ! innerTaken)
        {
            AddParam(track, at->GetInnerParam());
        }
        AddParam(track, at->GetOuterParam());
    }
    
    return track;
}


TEveElementList* AliEveAODTracks::ByPID()
{
    // Import AOD tracks, separate them into several containers by PID
    
    cout<<"*** AliEveAODTracks::ByPID() ***"<<endl;
    
    TEnv settings;
    AliEveInit::GetConfig(&settings);
    Color_t colors[15];
    // default color scheme by type:
    colors[0] = settings.GetValue("tracks.byType.electron",600);
    colors[1] = settings.GetValue("tracks.byType.muon",416);
    colors[2] = settings.GetValue("tracks.byType.pion",632);
    colors[3] = settings.GetValue("tracks.byType.kaon",400);
    colors[4] = settings.GetValue("tracks.byType.proton",797);
    colors[5] = settings.GetValue("tracks.byType.deuteron",797);
    colors[6] = settings.GetValue("tracks.byType.triton",797);
    colors[7] = settings.GetValue("tracks.byType.he3",797);
    colors[8] = settings.GetValue("tracks.byType.alpha",403);
    colors[9] = settings.GetValue("tracks.byType.photon",0);
    colors[10]= settings.GetValue("tracks.byType.pi0",616);
    colors[11]= settings.GetValue("tracks.byType.neutron",900);
    colors[12]= settings.GetValue("tracks.byType.kaon0",801);
    colors[13]= settings.GetValue("tracks.byType.elecon",920);
    colors[14]= settings.GetValue("tracks.byType.unknown",920);
    
    AliAODEvent *aod = AliEveEventManager::GetMaster()->AssertAOD();
    
    TEveElementList* cont = new TEveElementList("AOD Tracks by PID");
    gEve->AddElement(cont);
    
    const Int_t   nCont = 15;
    const Float_t maxR  = 520;
    const Float_t magF  = 0.1*aod->GetMagneticField();
    
    TEveTrackList *tl[nCont];
    Int_t          tc[nCont];
    Int_t          count = 0;
    
    tl[0] = new TEveTrackList("Electrons");
    tl[1] = new TEveTrackList("Muons");
    tl[2] = new TEveTrackList("Pions");
    tl[3] = new TEveTrackList("Kaons");
    tl[4] = new TEveTrackList("Protons");
    tl[5] = new TEveTrackList("Deuterons");
    tl[6] = new TEveTrackList("Tritons");
    tl[7] = new TEveTrackList("He3");
    tl[8] = new TEveTrackList("Alpha");
    tl[9] = new TEveTrackList("Photons");
    tl[10]= new TEveTrackList("Pi0");
    tl[11]= new TEveTrackList("Neutrons");
    tl[12]= new TEveTrackList("Kaon0");
    tl[13]= new TEveTrackList("EleCon");
    tl[14]= new TEveTrackList("Unknown");

    for (int i=0; i<15; i++) {
        tc[i] = 0;
        SetupPropagator(tl[i]->GetPropagator(), magF, maxR);
        tl[i]->SetMainColor(colors[i]);
        tl[i]->SetLineWidth(fWidth);
        cont->AddElement(tl[i]);
    }
    
    int pid = -1;
    AliAODTrack* at = NULL;
    
    for (Int_t n = 0; n < aod->GetNumberOfTracks(); ++n)
    {
        at = (AliAODTrack*)aod->GetTrack(n);
        
        bool good_cont = (at->IsOn(AliESDtrack::kITSin) && (!at->IsOn(AliESDtrack::kTPCin)));
        
        if(good_cont || fDrawNoRefit)
        {
            pid = at->GetMostProbablePID();
            TEveTrackList* tlist = tl[pid];
            ++tc[pid];
            ++count;
            
            AliEveTrack* track = MakeTrack(at, tlist);
            
            track->SetName(Form("AOD Track idx=%d, pid=%d", at->GetID(), pid));
            tlist->AddElement(track);
        }
    }
    for (Int_t ti = 0; ti < nCont; ++ti)
    {
        TEveTrackList* tlist = tl[ti];
        tlist->SetName(Form("%s [%d]", tlist->GetName(), tlist->NumChildren()));
        tlist->SetTitle(Form("N tracks=%d", tc[ti]));
        tlist->MakeTracks();
    }
    cont->SetTitle(Form("N all tracks = %d", count));
    cont->FindListTreeItem(gEve->GetListTree())->SetOpen(kTRUE);
    
    gEve->Redraw3D();
    
    return cont;
}





