//
//  AliEveAODTracks.cxx
//
//  Created by Jeremi Niedziela on 17/11/15.
//
//

#include "AliEveAODTracks.h"

#include <AliEveMagField.h>
#include <AliEveEventManager.h>
#include <AliEveInit.h>
#include <AliAODEvent.h>

#include <TGListTree.h>

#include <iostream>

using namespace std;

AliEveAODTracks::AliEveAODTracks() : fDrawNoRefit(true)
{}

AliEveAODTracks::~AliEveAODTracks()
{}

TString AliEveAODTracks::GetTitle(AliAODTrack* t)
{
    TString s("");
    
    Int_t label = t->GetLabel(), index = t->GetID();
    TString idx(index == kMinInt ? "<undef>" : Form("%d", index));
    TString lbl(label == kMinInt ? "<undef>" : Form("%d", label));
    
    Double_t p[3], v[3];
    t->GetXYZ(v);
    t->GetPxPyPz(p);
    
    s = Form("Index=%s, Label=%s\nChg=%d, Pdg=%d\n"
             "P  = (%.3f, %.3f, %.3f)\n"
             "V  = (%.3f, %.3f, %.3f)\n",
             idx.Data(), lbl.Data(), t->Charge(), AliPID::ParticleCode(t->GetMostProbablePID()),
             p[0], p[1], p[2],
             v[0], v[1], v[2]);
    
    return s;
}

void AliEveAODTracks::AddParam(AliEveTrack* track, const AliExternalTrackParam* tp)
{
    // Add additional track parameters as a path-mark to track.
    
    if (tp == 0) return;
    
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
    
    AliEveTrack* track = new AliEveTrack(at, cont->GetPropagator());
    track->SetAttLineAttMarker(cont);
    track->SetName(Form("AliEveTrack %d", at->GetID()));
    track->SetElementTitle(GetTitle(at));
    track->SetSourceObject(at);
    
    if (at->IsOn(AliAODTrack::kTPCrefit))
    {
        AddParam(track, at->GetInnerParam());
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

    Width_t width = settings.GetValue("tracks.width",1);
    Color_t colors[15];
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
    

    TEveElementList* cont = new TEveElementList("AOD Tracks by PID");
    gEve->AddElement(cont);
    
    const Int_t   nCont = 15;
    TEveTrackList *tl[nCont];
    Int_t          tc[nCont];
    
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
        tl[i]->GetPropagator()->SetMagFieldObj(new AliEveMagField());
        tl[i]->GetPropagator()->SetMaxR(520);
        tl[i]->SetMainColor(colors[i]);
        tl[i]->SetLineWidth(width);
        cont->AddElement(tl[i]);
    }
    
    int pid = -1;
    int count = 0;
    AliAODEvent *aod = AliEveEventManager::Instance()->AssertAOD();
    AliAODTrack* at = NULL;
    
    for (Int_t n = 0; n < aod->GetNumberOfTracks(); ++n)
    {
        at = (AliAODTrack*)aod->GetTrack(n);
        
        bool good_cont = true;
        string trackSelection = settings.GetValue("tracks.selection","");
        
        if(trackSelection == "ITSin_noTPCin"){
            good_cont = at->IsOn(AliESDtrack::kITSin && !at->IsOn(AliESDtrack::kTPCin));
        }
        else if(trackSelection == "noTISpureSA"){
            good_cont = !at->IsOn(AliESDtrack::kITSpureSA);
        }
        else if(trackSelection == "TPCrefit"){
            good_cont = at->IsOn(AliESDtrack::kTPCrefit);
        }
        else if(trackSelection == "TPCrefit_ITSrefit"){
            good_cont = at->IsOn(AliESDtrack::kTPCrefit) && at->IsOn(AliESDtrack::kITSrefit);
        }
        
        if(good_cont)
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





