// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include <AliEveESDSPDTracklets.h>

#include <TMath.h>
#include <TEveManager.h>

#include <TEveTrack.h>
#include <TEveTrackPropagator.h>

#include <AliESDEvent.h>
#include <AliESDVertex.h>
#include <AliMultiplicity.h>
#include <AliMagF.h>
#include <AliEveEventManager.h>
#include <AliEveTracklet.h>
#include <AliEveTrackCounter.h>

TEveElementList* AliEveESDSPDTracklets::Draw(Float_t radius, Width_t line_width,
                                   Float_t dPhiWindow, Float_t dThetaWindow,
                                   Float_t dPhiShift05T)
{
    // radius - cylindrical radius to which the tracklets should be extrapolated
    
    AliESDEvent     *esd = AliEveEventManager::Instance()->AssertESD();
    const AliESDVertex    *pv  = esd->GetPrimaryVertexSPD();
    const AliMultiplicity *mul = esd->GetMultiplicity();
    
    AliMagF *field = AliEveEventManager::AssertMagField();
    
    TEveElementList* cont = new TEveElementList("SPD Tracklets");
    gEve->AddElement(cont);
    
    TEveTrackList *tg = new TEveTrackList("Good");
    tg->SetMainColor(kCyan);
    tg->SetLineWidth(line_width);
    cont->AddElement(tg);
    
    TEveTrackPropagator* pg = tg->GetPropagator();
    pg->SetMaxR(radius);
    
    TEveTrackList *tb = new TEveTrackList("Bad");
    tb->SetMainColor(kMagenta);
    tb->SetLineWidth(line_width);
    cont->AddElement(tb);
    
    TEveTrackPropagator* pb = tb->GetPropagator();
    pb->SetMaxR(radius);
    
    const Float_t  Bz = TMath::Abs(field->SolenoidField());
    
    const Double_t dPhiShift     = dPhiShift05T / 5.0 * Bz;
    const Double_t dPhiWindow2   = dPhiWindow * dPhiWindow;
    const Double_t dThetaWindow2 = dThetaWindow * dThetaWindow;
    
    for (Int_t i = 0; i < mul->GetNumberOfTracklets(); ++i)
    {
        Float_t theta  = mul->GetTheta(i);
        Float_t phi    = mul->GetPhi(i);
        Float_t dTheta = mul->GetDeltaTheta(i);
        Float_t dPhi   = mul->GetDeltaPhi(i) - dPhiShift;
        Float_t d      = dPhi*dPhi/dPhiWindow2 + dTheta*dTheta/dThetaWindow2;
        
        TEveTrackList* tl = (d < 1.0f) ? tg : tb;
        
        AliEveTracklet *t = new AliEveTracklet(i, pv, theta, phi, tl->GetPropagator());
        t->SetAttLineAttMarker(tl);
        t->SetElementName(Form("Tracklet %d", i));
        t->SetElementTitle(Form("Id = %d\nEta=%.3f, Theta=%.3f, dTheta=%.3f\nPhi=%.3f dPhi=%.3f",
                                i, mul->GetEta(i), theta, dTheta, phi, dPhi));
        
        tl->AddElement(t);
    }
    
    tg->MakeTracks();
    tg->SetTitle(Form("N=%d", tg->NumChildren()));
    
    tb->MakeTracks();
    tb->SetTitle(Form("N=%d", tb->NumChildren()));
    
    if (AliEveTrackCounter::IsActive())
    {
        AliEveTrackCounter::fgInstance->RegisterTracklets(tg, kTRUE);
        AliEveTrackCounter::fgInstance->RegisterTracklets(tb, kFALSE);
    }
    else
    {
        //==========================================
        tb->SetLineStyle(1);
        tb->SetLineWidth(2);
        //==========================================
        
    }
    
    gEve->Redraw3D();
    
    return cont;
}
