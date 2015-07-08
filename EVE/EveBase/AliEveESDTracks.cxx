//
//  AliEveESDTracks.cxx
//  xAliRoot
//
//  Created by Jeremi Niedziela on 11/05/15.
//
//

#include "AliEveESDTracks.h"

#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TGListTree.h>
#include <TEveVSDStructs.h>


//#include <AliPWG0Helper.h>
#include <AliESDEvent.h>
#include <AliESDfriend.h>
#include <AliESDfriendTrack.h>
#include <AliExternalTrackParam.h>
#include <AliEveTrackCounter.h>
#include <AliEveMagField.h>
#include <AliEveEventManager.h>
#include <AliEveInit.h>

#include <iostream>

using namespace std;

AliEveESDTracks::AliEveESDTracks() :
fUseIPonFailedITSrefit(kFALSE),
fTrueField(kTRUE),
fRKstepper(kFALSE),
fAnalCuts(0),
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

AliEveESDTracks::~AliEveESDTracks()
{
    
}


void AliEveESDTracks::SetupPropagator(TEveTrackPropagator* trkProp,Float_t magF, Float_t maxR)
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

TString AliEveESDTracks::GetTitle(AliESDtrack* t)
{
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
    
    return s;
}

void AliEveESDTracks::AddParam(AliEveTrack* track, const AliExternalTrackParam* tp)
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

AliEveTrack* AliEveESDTracks::MakeTPCtrack(AliESDtrack *at, AliESDfriendTrack* aft, TEveTrackList* cont)
{
    // Make a TPC track representation and put it into given container.
    
    AliEveTrack* track = new AliEveTrack(at, cont->GetPropagator());
    track->SetAttLineAttMarker(cont);
    track->SetName(Form("AliEveTrack %d", at->GetID()));
    track->SetElementTitle(GetTitle(at));
    track->SetSourceObject(at);
    
    // Add inner/outer track parameters as start point and pathmark.
    if (at->GetInnerParam()) track->SetStartParams(at->GetInnerParam());
    else return NULL;
    if (aft->GetTPCOut()) AddParam(track, aft->GetTPCOut());
    else return NULL;
    
    return track;
}

AliEveTrack* AliEveESDTracks::MakeITSstandaloneTrack(AliESDtrack *at, AliESDfriendTrack* aft, TEveTrackList* cont)
{
    // Make a ITS standalone track representation and put it into given container.
    
    if ( !(!at->IsOn(AliESDtrack::kTPCin) &&
           at->IsOn(AliESDtrack::kITSout)) ) return NULL; //only ITS standalone
    AliEveTrack* track = new AliEveTrack(at, cont->GetPropagator());
    track->SetAttLineAttMarker(cont);
    track->SetName(Form("AliEveTrack %d", at->GetID()));
    track->SetElementTitle(GetTitle(at));
    track->SetSourceObject(at);
    
    // Add inner/outer track parameters as path-marks.
    if (aft->GetITSOut())
    {
        AddParam(track, aft->GetITSOut());
    }
    else return NULL;
    
    return track;
}


AliEveTrack* AliEveESDTracks::MakeITStrack(AliESDtrack *at, AliESDfriendTrack* aft, TEveTrackList* cont)
{
    // Make a ITS track representation and put it into given container.
    
    if ( (!at->IsOn(AliESDtrack::kTPCin) &&
          at->IsOn(AliESDtrack::kITSout)) ) return NULL; //ignore ITS standalone
    AliEveTrack* track = new AliEveTrack(at, cont->GetPropagator());
    track->SetAttLineAttMarker(cont);
    track->SetName(Form("AliEveTrack %d", at->GetID()));
    track->SetElementTitle(GetTitle(at));
    track->SetSourceObject(at);
    
    // Add inner/outer track parameters as path-marks.
    if (aft->GetITSOut())
    {
        AddParam(track, aft->GetITSOut());
    }
    else return NULL;
    
    return track;
}

AliEveTrack* AliEveESDTracks::MakeTrack(AliESDtrack *at, TEveTrackList* cont)
{
    // Make a standard track representation and put it into given container.
    
    // Choose which parameters to use a track's starting point.
    // If gkFixFailedITSExtr is TRUE (FALSE by default) and
    // if ITS refit failed, take track parameters at inner TPC radius.
    
    const AliExternalTrackParam* tp = at;
    
    Bool_t innerTaken = kFALSE;
    if ( ! at->IsOn(AliESDtrack::kITSrefit) && fUseIPonFailedITSrefit)
    {
        tp = at->GetInnerParam();
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


TEveTrackList* AliEveESDTracks::TPCtracks()
{
    AliESDEvent* esd = AliEveEventManager::GetMaster()->AssertESD();
    AliESDfriend* esd_friend = AliEveEventManager::AssertESDfriend();
    
    TEveTrackList* cont = new TEveTrackList("TPC Tracks");
    cont->SetMainColor(kMagenta);
    
    SetupPropagator(cont->GetPropagator(),
                    0.1*esd->GetMagneticField(), 520);
    
    gEve->AddElement(cont);
    
    Int_t count = 0;
    for (Int_t n = 0; n < esd->GetNumberOfTracks(); ++n)
    {
        ++count;
        if (!esd->GetTrack(n)) continue;
        if (!esd_friend->GetTrack(n)) continue;
        AliEveTrack* track = MakeTPCtrack(esd->GetTrack(n), esd_friend->GetTrack(n), cont);
        if (!track) continue;
        
        cont->AddElement(track);
    }
    cont->SetTitle(Form("N=%d", count));
    cont->MakeTracks();
    
    gEve->Redraw3D();
    
    return cont;
}

TEveTrackList* AliEveESDTracks::ITStracks()
{
    AliESDEvent* esd = AliEveEventManager::GetMaster()->AssertESD();
    AliESDfriend* esd_friend = AliEveEventManager::AssertESDfriend();
    
    TEveTrackList* cont = new TEveTrackList("ITS Tracks");
    cont->SetMainColor(kMagenta+3);
    
    SetupPropagator(cont->GetPropagator(),
                    0.1*esd->GetMagneticField(), 520);
    cont->GetPropagator()->SetMaxR(85.0);
    cont->SetLineWidth(1);
    
    gEve->AddElement(cont);
    
    Int_t count = 0;
    for (Int_t n = 0; n < esd->GetNumberOfTracks(); ++n)
    {
        ++count;
        if (!esd->GetTrack(n)) continue;
        if (!esd_friend->GetTrack(n)) continue;
        AliEveTrack* track = MakeITStrack(esd->GetTrack(n), esd_friend->GetTrack(n), cont);
        if (!track) continue;
        
        cont->AddElement(track);
    }
    cont->SetTitle(Form("N=%d", count));
    cont->MakeTracks();
    
    gEve->Redraw3D();
    
    return cont;
}

TEveTrackList* AliEveESDTracks::ITSstandaloneTracks()
{
    AliESDEvent* esd = AliEveEventManager::GetMaster()->AssertESD();
    AliESDfriend* esd_friend = AliEveEventManager::AssertESDfriend();
    
    TEveTrackList* cont = new TEveTrackList("ITS Standalone Tracks");
    cont->SetMainColor(kBlue);
    
    SetupPropagator(cont->GetPropagator(),
                    0.1*esd->GetMagneticField(), 520);
    cont->GetPropagator()->SetMaxR(85.0);
    cont->SetLineWidth(1);
    
    gEve->AddElement(cont);
    
    Int_t count = 0;
    for (Int_t n = 0; n < esd->GetNumberOfTracks(); ++n)
    {
        ++count;
        if (!esd->GetTrack(n)) continue;
        if (!esd_friend->GetTrack(n)) continue;
        AliEveTrack* track = MakeITSstandaloneTrack(esd->GetTrack(n), esd_friend->GetTrack(n), cont);
        if (!track) continue;
        
        cont->AddElement(track);
    }
    cont->SetTitle(Form("N=%d", count));
    cont->MakeTracks();
    
    gEve->Redraw3D();
    
    return cont;
}

TEveTrackList* AliEveESDTracks::Tracks()
{
    AliESDEvent* esd = AliEveEventManager::GetMaster()->AssertESD();
    
    TEveTrackList* cont = new TEveTrackList("ESD Tracks");
    cont->SetMainColor(6);
    
    SetupPropagator(cont->GetPropagator(),
                    0.1*esd->GetMagneticField(), 520);
    
    gEve->AddElement(cont);
    
    Int_t count = 0;
    for (Int_t n = 0; n < esd->GetNumberOfTracks(); ++n)
    {
        ++count;
        AliEveTrack* track = MakeTrack(esd->GetTrack(n), cont);
        
        cont->AddElement(track);
    }
    cont->SetTitle(Form("N=%d", count));
    cont->MakeTracks();
    
    gEve->Redraw3D();
    
    return cont;
}

TEveTrackList* AliEveESDTracks::MItracks()
{
    AliESDEvent* esd = AliEveEventManager::GetMaster()->AssertESD();
    
    TEveTrackList* cont = new TEveTrackList("ESD Tracks MI");
    cont->SetLineColor(5);
    gEve->AddElement(cont);
    
    Int_t count = 0;
    for (Int_t n = 0; n < esd->GetNumberOfTracks(); ++n)
    {
        ++count;
        AliESDtrack* at = esd->GetTrack(n);
        AliEveTrack* l = new AliEveTrack(at, cont->GetPropagator());
        l->SetName(Form("ESDTrackMI %d", at->GetID()));
        l->SetElementTitle(GetTitle(at));
        l->SetAttLineAttMarker(cont);
        l->SetSourceObject(at);
        
        at->FillPolymarker(l, esd->GetMagneticField(), 0, 250, 5);
        
        l->SetLockPoints(kTRUE);
        cont->AddElement(l);
    }
    cont->SetTitle(Form("N=%d", count));
    
    gEve->Redraw3D();
    
    return cont;
}

TEveTrackList* AliEveESDTracks::TracksFromArray(TCollection* col, AliESDEvent* esd)
{
    // Retrieves AliESDTrack's from collection.
    // See example usage with AliAnalysisTrackCuts in the next function.
    
    if (esd == 0) esd = AliEveEventManager::GetMaster()->AssertESD();
    
    TEveTrackList* cont = new TEveTrackList("ESD Tracks");
    cont->SetMainColor(6);
    
    SetupPropagator(cont->GetPropagator(),0.1*esd->GetMagneticField(), 520);
    gEve->AddElement(cont);
    
    Int_t    count = 0;
    TIter    next(col);
    TObject *obj;
    while ((obj = next()) != 0)
    {
        if (obj->IsA()->InheritsFrom("AliESDtrack") == kFALSE)
        {
            Warning("TracksFromArray", "Object '%s', '%s' is not an AliESDtrack.",
                    obj->GetName(), obj->GetTitle());
            continue;
        }
        
        ++count;
        AliESDtrack* at = reinterpret_cast<AliESDtrack*>(obj);
        
        AliEveTrack* track = MakeTrack(at, cont);
        cont->AddElement(track);
    }
    cont->SetTitle(Form("N=%d", count));
    cont->MakeTracks();
    
    gEve->Redraw3D();
    
    return cont;
}

void AliEveESDTracks::AliAnalCutsDemo()
{
    AliESDEvent* esd = AliEveEventManager::GetMaster()->AssertESD();
    
    AliESDtrackCuts atc;
    atc.SetPtRange(0.1, 5);
    atc.SetRapRange(-1, 1);
    
    TracksFromArray(atc.GetAcceptedTracks(esd), esd);
}

Float_t AliEveESDTracks::GetSigmaToVertex(AliESDtrack* esdTrack)
{
    // Taken from: PWG0/esdTrackCuts/AliESDtrackCuts.cxx
    // Float_t AliESDtrackCuts::GetSigmaToVertex(AliESDtrack* esdTrack)
    
    Float_t b[2];
    Float_t bRes[2];
    Float_t bCov[3];
    esdTrack->GetImpactParameters(b,bCov);
    if (bCov[0] <= 0 || bCov[2] <= 0)
    {
        printf("Estimated b resolution lower or equal zero!\n");
        bCov[0] = bCov[2] = 0;
    }
    bRes[0] = TMath::Sqrt(bCov[0]);
    bRes[1] = TMath::Sqrt(bCov[2]);
    
    // -----------------------------------
    // How to get to a n-sigma cut?
    //
    // The accumulated statistics from 0 to d is
    //
    // ->  Erf(d/Sqrt(2)) for a 1-dim gauss (d = n_sigma)
    // ->  1 - Exp(-d**2) for a 2-dim gauss (d*d = dx*dx + dy*dy != n_sigma)
    //
    // It means that for a 2-dim gauss: n_sigma(d) = Sqrt(2)*ErfInv(1 - Exp((-x**2)/2)
    // Can this be expressed in a different way?
    
    if (bRes[0] == 0 || bRes[1] == 0)
        return -1;
    
    Float_t d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));
    
    // stupid rounding problem screws up everything:
    // if d is too big, TMath::Exp(...) gets 0, and TMath::ErfInverse(1) that should be infinite, gets 0 :(
    if (TMath::Exp(-d * d / 2) < 1e-10)
        return 1000;
    
    d = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
    return d;
}

TEveElementList* AliEveESDTracks::ByCategory()
{
    // Import ESD tracks, separate them into several containers
    // according to primary-vertex cut and ITS&TPC refit status.
    
    AliESDEvent* esd = AliEveEventManager::GetMaster()->AssertESD();
    
    TEveElementList* cont = new TEveElementList("ESD Tracks by category");
    gEve->AddElement(cont);
    
    const Int_t   nCont = 9;
    const Float_t maxR  = 520;
    const Float_t magF  = 0.1*esd->GetMagneticField();
    
    TEveTrackList *tl[nCont];
    Int_t          tc[nCont];
    Int_t          count = 0;
    
    tl[0] = new TEveTrackList("Sigma < 3");
    tl[1] = new TEveTrackList("3 < Sigma < 5");
    tl[2] = new TEveTrackList("5 < Sigma");
    tl[3] = new TEveTrackList("no ITS refit; Sigma < 5");
    tl[4] = new TEveTrackList("no ITS refit; Sigma > 5");
    tl[5] = new TEveTrackList("no TPC refit");
    tl[6] = new TEveTrackList("ITS ncl>=3 & SPD Inner");
    tl[7] = new TEveTrackList("ITS ncl>=3 & b<3 cm");
    tl[8] = new TEveTrackList("ITS others");
    
    for (int i=0; i<9; i++) {
        tc[i] = 0;
        SetupPropagator(tl[i]->GetPropagator(), magF, maxR);
        tl[i]->SetMainColor(fColorsByCategory[i]);
        tl[i]->SetLineWidth(fWidth);
        cont->AddElement(tl[i]);
    }
    
    for (Int_t n = 0; n < esd->GetNumberOfTracks(); ++n)
    {
        AliESDtrack* at = esd->GetTrack(n);
        
        Float_t s  = GetSigmaToVertex(at);
        Int_t   ti;
        if      (s <  3) ti = 0;
        else if (s <= 5) ti = 1;
        else             ti = 2;
        
        Int_t    nclits;
        Double_t dtobeam;
        
        if (at->IsOn(AliESDtrack::kITSin) && ! at->IsOn(AliESDtrack::kTPCin))
        {
            UChar_t itsclmap = at->GetITSClusterMap();
            Bool_t  spdinner = (itsclmap & 3) != 0;
            
            nclits = 0;
            for (Int_t iter = 0; iter < 6; ++iter)
                if (itsclmap & (1 << iter)) nclits ++;
            
            Double_t xyz[3];
            at->GetXYZ(xyz);
            dtobeam = TMath::Hypot(xyz[0], xyz[1]);
            
            if ((nclits >= 3) && (spdinner))
                ti = 6;
            else if ((nclits >= 3) && (dtobeam < 3.0))
                ti = 7;
            else
                ti = 8;
        }
        else if (at->IsOn(AliESDtrack::kTPCin) && ! at->IsOn(AliESDtrack::kTPCrefit))
        {
            ti = 5;
        }
        else if ( ! at->IsOn(AliESDtrack::kITSrefit))
        {
            ti = (ti == 2) ? 4 : 3;
        }
        
        TEveTrackList* tlist = tl[ti];
        ++tc[ti];
        ++count;
        
        AliEveTrack* track = MakeTrack(at, tlist);
        if (ti<6)
            track->SetName(Form("ESD Track idx=%d, sigma=%5.3f", at->GetID(), s));
        else
            track->SetName(Form("ESD Track idx=%d, dxy=%5.3f cl=%i", at->GetID(), dtobeam, nclits));
        tlist->AddElement(track);
    }
    
    for (Int_t ti = 0; ti < nCont; ++ti)
    {
        TEveTrackList* tlist = tl[ti];
        tlist->SetName(Form("%s [%d]", tlist->GetName(), tlist->NumChildren()));
        tlist->SetTitle(Form("N tracks=%d", tc[ti]));
        
        tlist->MakeTracks();
        
        //    Bool_t good_cont = ti <= 1;
        Bool_t good_cont = ((ti == 6) || (ti == 7));
        if (AliEveTrackCounter::IsActive())
        {
            AliEveTrackCounter::fgInstance->RegisterTracks(tlist, good_cont);
        }
        else
        {
            if ( ! good_cont && fDashNoRefit)
                tlist->SetLineStyle(6);
        }
        if(!good_cont)
        {
            tlist->SetRnrLine(fDrawNoRefit);
        }
    }
    cont->SetTitle(Form("N all tracks = %d", count));
    // ??? The following does not always work:
    cont->FindListTreeItem(gEve->GetListTree())->SetOpen(kTRUE);
    
    gEve->Redraw3D();
    
    return cont;
}

TEveElementList* AliEveESDTracks::ByType()
{
    // Import ESD tracks, separate them into several containers
    // according to primary-vertex cut and ITS&TPC refit status.
    
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

    
    AliESDEvent *esd = AliEveEventManager::GetMaster()->AssertESD();
    
    TEveElementList* cont = new TEveElementList("ESD Tracks by PID");
    gEve->AddElement(cont);
    
    const Int_t   nCont = 15;
    const Float_t maxR  = 520;
    const Float_t magF  = 0.1*esd->GetMagneticField();
    
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
    AliESDtrack* at = NULL;
    
    for (Int_t n = 0; n < esd->GetNumberOfTracks(); ++n)
    {
        at = esd->GetTrack(n);
        
        bool good_cont = (at->IsOn(AliESDtrack::kITSin) && (!at->IsOn(AliESDtrack::kTPCin)));
        
        if(good_cont || fDrawNoRefit)
        {
            pid = at->GetPID();
            TEveTrackList* tlist = tl[pid];
            ++tc[pid];
            ++count;
            
            AliEveTrack* track = MakeTrack(at, tlist);
            
            track->SetName(Form("ESD Track idx=%d, pid=%d", at->GetID(), pid));
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
    // ??? The following does not always work:
    cont->FindListTreeItem(gEve->GetListTree())->SetOpen(kTRUE);
    
    gEve->Redraw3D();
    
    return cont;
}

TEveElementList* AliEveESDTracks::ByAnalCuts()
{
    AliESDEvent* esd = AliEveEventManager::GetMaster()->AssertESD();
    
    if (fAnalCuts == 0)
    {
        gSystem->Load("libPWGUDbase");
        gROOT->ProcessLine(".L $ALICE_ROOT/PWGUD/CreateStandardCuts.C");
        /*Int_t mode = AliPWG0Helper::kTPC;
         if (TMath::Abs(esd->GetMagneticField()) > 0.01)
         mode |= AliPWG0Helper::kFieldOn;
         gROOT->ProcessLine(Form("fAnalCuts = CreateTrackCuts(%d, kFALSE)", mode));*/
    }
    
    TEveElementList* cont = new TEveElementList("ESD Tracks by Analysis Cuts");
    gEve->AddElement(cont);
    
    const Int_t   nCont = 2;
    const Float_t maxR  = 520;
    const Float_t magF  = 0.1*esd->GetMagneticField();
    
    TEveTrackList *tl[nCont];
    Int_t          tc[nCont];
    Int_t          count = 0;
    
    tl[0] = new TEveTrackList("Passed");
    tc[0] = 0;
    SetupPropagator(tl[0]->GetPropagator(), magF, maxR);
    tl[0]->SetMainColor(3);
    cont->AddElement(tl[0]);
    
    tl[1] = new TEveTrackList("Rejected");
    tc[1] = 0;
    SetupPropagator(tl[1]->GetPropagator(), magF, maxR);
    tl[1]->SetMainColor(kRed);
    cont->AddElement(tl[1]);
    
    for (Int_t n = 0; n < esd->GetNumberOfTracks(); ++n)
    {
        AliESDtrack* at = esd->GetTrack(n);
        
        Float_t s  = GetSigmaToVertex(at);
        Int_t   ti = (fAnalCuts->AcceptTrack(at)) ? 0 : 1;
        
        TEveTrackList* tlist = tl[ti];
        ++tc[ti];
        ++count;
        
        AliEveTrack* track = MakeTrack(at, tlist);
        track->SetName(Form("ESD Track idx=%d, sigma=%5.3f", at->GetID(), s));
        tlist->AddElement(track);
    }
    
    for (Int_t ti = 0; ti < nCont; ++ti)
    {
        TEveTrackList* tlist = tl[ti];
        tlist->SetName(Form("%s [%d]", tlist->GetName(), tlist->NumChildren()));
        tlist->SetTitle(Form("N tracks=%d", tc[ti]));
        
        tlist->MakeTracks();
        
        Bool_t good_cont = ti < 1;
        if (AliEveTrackCounter::IsActive())
        {
            AliEveTrackCounter::fgInstance->RegisterTracks(tlist, good_cont);
        }
        else
        {
            if ( ! good_cont)
                tlist->SetLineStyle(6);
        }
    }
    cont->SetTitle(Form("N all tracks = %d", count));
    // ??? The following does not always work:
    cont->FindListTreeItem(gEve->GetListTree())->SetOpen(kTRUE);
    
    gEve->Redraw3D();
    
    return cont;
}






