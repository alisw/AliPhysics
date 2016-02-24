// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include <TEveManager.h>
#include <TEveCompound.h>


#include <AliESDEvent.h>
#include <AliESDVertex.h>
#include <AliEveEventManager.h>

#include <iostream>

using namespace std;

void AliEvePrimaryVertex::PrimaryVertex(EVertexType type,EVertexStyle style, Bool_t use_sigma, Float_t fx, Float_t fy, Float_t fz)
{
    AliESDEvent  *esd = AliEveEventManager::Instance()->AssertESD();
    
    const AliESDVertex *pv = NULL;
    Color_t color;
    const char* name;
    
    switch (type)
    {
        case kGlobal:
            pv  = esd->GetPrimaryVertex();
            color = 7;
            name = "Primary Vertex";
            break;
        case kSPD:
            pv = esd->GetPrimaryVertexSPD();
            color = 6;
            name = "Primary Vertex SPD";
            break;
        case kTPC:
            pv = esd->GetPrimaryVertexTPC();
            color = 5;
            name = "Primary Vertex TPC";
            break;
        default:break;
    }
    
    if(!pv){
        cout<<"AliEvePrimaryVertex::PrimaryVertex -- Primary vertex not available"<<endl;
        return;
    }
    if (!pv->GetStatus()){
        cout<<"AliEvePrimaryVertex::PrimaryVertex -- Primary vertex not available"<<endl;
        return;
    }
    
    TEveStraightLineSet* ls;
    
    switch (style)
    {
        case kCross:
            ls = MakeVertexCross(pv, use_sigma, fx, fy, fz);
            break;
        case kEllipse:
            ls = MakeVertexEllipse(pv, use_sigma, fx, fy, fz);
            break;
        case kBox:
            ls = MakeVertexBox(pv, use_sigma, fx, fy, fz);
            break;
        default:break;
    }
    
    ls->ApplyVizTag("REC PVTX Box");
    ls->SetMainColor(color);
    TEveCompound* parent = dynamic_cast<TEveCompound*>(AliEveEventManager::Instance()->FindChild(name));
    if (parent == 0)
    {
        parent = new TEveCompound(name);
        parent->OpenCompound();
        parent->SetMainColor(color);
        AliEveEventManager::Instance()->AddElement(parent);
    }
    parent->AddElement(ls);
    gEve->Redraw3D();
}

TEveStraightLineSet* AliEvePrimaryVertex::MakeVertexCross(const AliESDVertex* v, Bool_t use_sigma, Float_t fx, Float_t fy, Float_t fz)
{
    Double_t x[3], e[3];
    v->GetXYZ(x); v->GetSigmaXYZ(e);
    
    TEveStraightLineSet* ls = new TEveStraightLineSet("Cross");
    TString title;
    if (use_sigma)
    {
        e[0] *= fx; e[1] *= fy; e[2] *= fz;
        title += Form("+- %.1f*sigma_x, %.1f*sigma_y, %.1f*sigma_z", fx, fy, fz);
    }
    else
    {
        e[0] = fx; e[1] = fy; e[2] = fz;
        title += Form("+- %.1f cm x %.1f cm x %.1f cm", fx, fy, fz);
    }
    title += Form("\nx=%.5f, y=%.5f, z=%.5f\nsx=%.5f, sy=%.5f, sz=%.5f",
                  x[0], x[1], x[2], e[0], e[1], e[2]);
    ls->SetTitle(title);
    
    ls->AddLine(e[0], 0,    0,   -e[0], 0,    0);
    ls->AddLine(0,    e[1], 0,    0,   -e[1], 0);
    ls->AddLine(0,    0,    e[2], 0,    0,   -e[2]);
    
    ls->RefMainTrans().SetPos(x);
    return ls;
}

TEveStraightLineSet* AliEvePrimaryVertex::MakeVertexEllipse(const AliESDVertex* v, Bool_t use_sigma, Float_t fx, Float_t fy, Float_t fz)
{
    Double_t x[3], e[3];
    v->GetXYZ(x); v->GetSigmaXYZ(e);
    
    TEveStraightLineSet* ls = new TEveStraightLineSet("Ellipse");
    TString title;
    if (use_sigma)
    {
        e[0] *= fx; e[1] *= fy; e[2] *= fz;
        title += Form("+- %.1f*sigma_x, %.1f*sigma_y, %.1f sigma_z", fx, fy, fz);
    }
    else
    {
        e[0] = fx; e[1] = fy; e[2] = fz;
        title += Form("+- %.1f cm x %.1f cm x %.1f cm", fx, fy, fz);
    }
    title += Form("\nx=%.5f, y=%.5f, z=%.5f\nsx=%.5f, sy=%.5f, sz=%.5f",
                  x[0], x[1], x[2], e[0], e[1], e[2]);
    ls->SetTitle(title);
    
    const Int_t   N = 32;
    const Float_t S = 2*TMath::Pi()/N;
    
    Float_t a = e[0], b = e[1];
    for (Int_t i = 0; i<N; i++)
        ls->AddLine(a*TMath::Cos(i*S)  , b*TMath::Sin(i*S)  , 0,
                    a*TMath::Cos(i*S+S), b*TMath::Sin(i*S+S), 0);
    
    a = e[0]; b = e[2];
    for (Int_t i = 0; i<N; i++)
        ls->AddLine(a*TMath::Cos(i*S)  , 0, b*TMath::Sin(i*S),
                    a*TMath::Cos(i*S+S), 0, b*TMath::Sin(i*S+S));
    
    a = e[1]; b = e[2];
    for (Int_t i = 0; i<N; i++)
        ls->AddLine(0, a*TMath::Cos(i*S)  ,  b*TMath::Sin(i*S),
                    0, a*TMath::Cos(i*S+S),  b*TMath::Sin(i*S+S));
    
    ls->RefMainTrans().SetPos(x);
    return ls;
}

TEveStraightLineSet* AliEvePrimaryVertex::MakeVertexBox(const AliESDVertex* v, Bool_t use_sigma, Float_t fx, Float_t fy, Float_t fz)
{
    Double_t x[3], e[3];
    v->GetXYZ(x); v->GetSigmaXYZ(e);
    
    TEveStraightLineSet* ls = new TEveStraightLineSet("Box");
    TString title;
    if (use_sigma)
    {
        e[0] *= fx; e[1] *= fy; e[2] *= fz;
        title += Form("+- %.1f*sigma_x, %.1f*sigma_y, %.1f*sigma_z", fx, fy, fz);
    }
    else
    {
        e[0] = fx; e[1] = fy; e[2] = fz;
        title += Form("+- %.1f cm x %.1f cm x %.1f cm", fx, fy, fz);
    }
    title += Form("\nx=%.5f, y=%.5f, z=%.5f\nsx=%.5f, sy=%.5f, sz=%.5f",
                  x[0], x[1], x[2], e[0], e[1], e[2]);
    ls->SetTitle(title);
    
    // pos z
    ls->AddLine( e[0],  e[1],  e[2],  e[0], -e[1],  e[2]);
    ls->AddLine( e[0], -e[1],  e[2], -e[0], -e[1],  e[2]);
    ls->AddLine(-e[0], -e[1],  e[2], -e[0],  e[1],  e[2]);
    ls->AddLine(-e[0],  e[1],  e[2],  e[0],  e[1],  e[2]);
    // lines along z
    ls->AddLine( e[0],  e[1],  e[2],  e[0],  e[1], -e[2]);
    ls->AddLine( e[0], -e[1],  e[2],  e[0], -e[1], -e[2]);
    ls->AddLine(-e[0], -e[1],  e[2], -e[0], -e[1], -e[2]);
    ls->AddLine(-e[0],  e[1],  e[2], -e[0],  e[1], -e[2]);
    // neg z
    ls->AddLine( e[0],  e[1], -e[2],  e[0], -e[1], -e[2]);
    ls->AddLine( e[0], -e[1], -e[2], -e[0], -e[1], -e[2]);
    ls->AddLine(-e[0], -e[1], -e[2], -e[0],  e[1], -e[2]);
    ls->AddLine(-e[0],  e[1], -e[2],  e[0],  e[1], -e[2]);
    
    ls->RefMainTrans().SetPos(x);
    return ls;
}
