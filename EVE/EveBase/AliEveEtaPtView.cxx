// Author: Jeremi Niedziela 2015

#include "AliEveEtaPtView.h"
#include <TBrowser.h>
#include <TMath.h>
#include <TEveStraightLineSet.h>
#include <TRandom.h>
#include <TEveBoxSet.h>
#include <TGeoTube.h>
#include <TEveGeoShape.h>

using namespace TMath;

ClassImp(AliEveEtaPtView)

AliEveEtaPtView* AliEveEtaPtView::fInstance = 0;

AliEveEtaPtView* AliEveEtaPtView::Instance()
{
    return fInstance;
}

AliEveEtaPtView::AliEveEtaPtView() :
    f3DView(0),
    fGeom(0)
{
    if (fInstance){
        throw TEveException("AliEveEtaPtView::AliEveEtaPtView already instantiated.");
    }
    fInstance = this;

    f3DView = gEve->SpawnNewViewer("3D Eta-Pt View", "");
    f3DView->AddScene(gEve->GetGlobalScene());
    f3DView->AddScene(gEve->GetEventScene());
}

AliEveEtaPtView::~AliEveEtaPtView()
{
    fGeom->DestroyElements();
    delete fGeom;
}

void AliEveEtaPtView::InitGeom()
{
    Float_t x=0, y=0, z=0;

    TGeoTube *tube = new TGeoTube(0,500,100);
    TGeoTube *tube2 = new TGeoTube(0,400,100);
    TGeoTube *tube3 = new TGeoTube(0,400,100);
    TGeoTube *tube4 = new TGeoTube(0,300,100);
    TGeoTube *tube5 = new TGeoTube(0,300,100);

    TEveGeoShape *eveTube = new TEveGeoShape();
    TEveGeoShape *eveTube2 = new TEveGeoShape();
    TEveGeoShape *eveTube3 = new TEveGeoShape();
    TEveGeoShape *eveTube4 = new TEveGeoShape();
    TEveGeoShape *eveTube5 = new TEveGeoShape();
    
    eveTube->SetShape(tube);
    eveTube2->SetShape(tube2);
    eveTube3->SetShape(tube3);
    eveTube4->SetShape(tube4);
    eveTube5->SetShape(tube5);
    
    eveTube->SetMainColor(kGreen);
    eveTube2->SetMainColor(kRed);
    eveTube3->SetMainColor(kRed);
    eveTube4->SetMainColor(kYellow);
    eveTube5->SetMainColor(kYellow);
    
    eveTube->SetMainTransparency(70);
    eveTube2->SetMainTransparency(70);
    eveTube3->SetMainTransparency(70);
    eveTube4->SetMainTransparency(70);
    eveTube5->SetMainTransparency(70);
    
    
    TEveTrans& t = eveTube->RefMainTrans();
    TEveTrans& t2 = eveTube2->RefMainTrans();
    TEveTrans& t3 = eveTube3->RefMainTrans();
    TEveTrans& t4 = eveTube4->RefMainTrans();
    TEveTrans& t5 = eveTube5->RefMainTrans();

    t.SetPos(x, y, z);
    t2.SetPos(x, y, z-200);
    t3.SetPos(x, y, z+200);
    t4.SetPos(x, y, z-400);
    t5.SetPos(x, y, z+400);

    
    gEve->AddElement(eveTube);
    gEve->AddElement(eveTube2);
    gEve->AddElement(eveTube3);
    gEve->AddElement(eveTube4);
    gEve->AddElement(eveTube5);
    
    gEve->Redraw3D(kTRUE);
    
}

