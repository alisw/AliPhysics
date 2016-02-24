// $Id$
// Main authors: Paraskevi Ganoti: 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include <AliEveESDKinks.h>

#include <AliEveEventManager.h>
#include <AliPID.h>
#include <AliESDEvent.h>

#include <TVector3.h>
#include <TEveVector.h>
#include <TEveTrackPropagator.h>
#include <TEveVSDStructs.h>
#include <TEveManager.h>

#include <iostream>

using namespace std;

void AliEveESDKinks::InitRecTrackDother(TEveRecTrack& rt, const AliExternalTrackParam* tp)
{
    Double_t pbuf[3], vbuf[3];
    
    rt.fSign = tp->GetSign();
    tp->GetXYZ(vbuf);     rt.fV.Set(vbuf);
    tp->GetPxPyPz(pbuf);  rt.fP.Set(pbuf);
    
    rt.fBeta = 1; // ep/TMath::Sqrt(ep*ep + mc*mc);
}

void AliEveESDKinks::InitRecTrackDaughter(TEveRecTrack& rt, const AliExternalTrackParam* tp, TEveVector* svt,TEveVector* spt)
{
    rt.fSign = tp->GetSign();
    rt.fV.Set(*svt);
    rt.fP.Set(*spt);
    
    rt.fBeta = 1; // ep/TMath::Sqrt(ep*ep + mc*mc);
}

AliEveKink* AliEveESDKinks::MakeKink(TEveTrackPropagator* rnrStyleMoth,TEveTrackPropagator* rnrStyleDaugh, AliESDtrack* moth, AliESDtrack* daug, AliESDkink* kink, Int_t i)
{
    TEveRecTrack  rcMoth;
    TEveRecTrack  rcDaug;
    TEveRecKink   rcKink;
    
    const TVector3 p1(kink->GetMotherP());
    rcKink.fPMother.Set(p1);
    const TVector3 p2(kink->GetDaughterP());
    rcKink.fPDaughter.Set(p2);
    
    const TVector3 r1(kink->GetPosition());
    rcKink.fVKink.Set(r1);
    for (Int_t j=0; j<3; ++j) rcKink.fKinkAngle[j]=kink->GetAngle(j);

    Double_t r[3], r2[3];
    const AliExternalTrackParam *trackParam = moth->GetTPCInnerParam();
    if(!trackParam)
    {
        trackParam = moth->GetInnerParam();
    }
    if(trackParam){
        trackParam->GetXYZ(r2);
        rcKink.fVMother.Set(r2);
        daug->GetOuterParam()->GetXYZ(r);  rcKink.fVDaughter.Set(r);
        InitRecTrackDother(rcMoth, trackParam);
        rcMoth.fIndex = moth->GetID();
    }
    else{
        cout<<"AliEveESDKinks::MakeKink : couldn't get TPC inner params"<<endl;
    }
    InitRecTrackDaughter(rcDaug, daug->GetOuterParam(), &rcKink.fVKink, &rcKink.fPDaughter);
    rcDaug.fIndex = daug->GetID();
    AliEveKink* myKink = new AliEveKink(&rcMoth, &rcDaug, &rcKink, rnrStyleMoth,rnrStyleDaugh);
    myKink->SetElementName(Form("ESDkink %d  \n", i));
    myKink->SetESDKinkIndex(i);
    for (Int_t j=0; j<3; ++j) myKink->SetKinkAngle(j, kink->GetAngle(j));
    Double_t daugProbability[10];
    Double_t daugP = 0.0;
    daug->GetESDpid(daugProbability);
    daugP = daug->P();
    
    // ****** Tentative particle type "concentrations"
    Double_t c[5]={0.01, 0.01, 0.85, 0.10, 0.05};
    AliPID::SetPriors(c);
    
    AliPID daugPid(daugProbability);
    
    Int_t   daugMostProbPdg =  0;
    
    switch (daugPid.GetMostProbable()){
        case 0:
            daugMostProbPdg =   11; break;
        case 1:
            daugMostProbPdg =   13; break;
        case 2:
            daugMostProbPdg =  211; break;
        default :
            daugMostProbPdg =  13; break;
    }
    
    Float_t daugMaxProbPid  = daugPid.GetProbability(daugPid.GetMostProbable());
    myKink->SetMaxProbPdgPid(daugMostProbPdg,daugMaxProbPid);//????????????
    return myKink;
}


AliEveKinkList* AliEveESDKinks::Draw()
{
    AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();
    AliEveKinkList* cont = new AliEveKinkList("ESD kink");

    cont->SetMainColor(3); // green
    TEveTrackPropagator* rnrStyleMoth = cont->GetPropagatorMoth();
    rnrStyleMoth->SetMagField( 0.1*esd->GetMagneticField() );
    TEveTrackPropagator* rnrStyleDaugh = cont->GetPropagatorDaugh();
    rnrStyleDaugh->SetMagField( 0.1*esd->GetMagneticField() );
    rnrStyleDaugh->SetMaxR(520);

    gEve->AddElement(cont);
    
    Int_t count = 0;
    //   for (Int_t n=0; n<esd->GetNumberOfKinks(); ++n)
    //   {
    //     AliESDkink *kink = esd->GetKink(n);  //???????????
    //     printf("kink number = %d,  label of mother = %d , label of daug = %d --", n, kink->GetLabel(0), kink->GetLabel(1));
    //   }   To be investigated...
    for (Int_t n=0; n<esd->GetNumberOfTracks(); ++n)
    {
        AliESDtrack* mtrack = esd->GetTrack(n);
        if(mtrack->GetKinkIndex(0)<0)
        {
            AliESDkink *kink = new AliESDkink;
            
            kink=esd->GetKink(TMath::Abs(mtrack->GetKinkIndex(0))-1);
            
            for (Int_t m=0; m<esd->GetNumberOfTracks(); ++m)
            {
                AliESDtrack * dtrack = esd->GetTrack(m);

                if((dtrack->GetKinkIndex(0)>0)&&(dtrack->GetKinkIndex(0)==TMath::Abs(mtrack->GetKinkIndex(0))))
                {
                    AliESDtrack* mothTr = esd->GetTrack(n);
                    AliESDtrack* daugTr = esd->GetTrack(m);
                    AliEveKink* myKink = MakeKink(rnrStyleMoth, rnrStyleDaugh, mothTr, daugTr, kink, (TMath::Abs(mtrack->GetKinkIndex(0))-1));

                    if (myKink)
                    {
                        gEve->AddElement(myKink, cont);
                        ++count;
                    }
                }
            }  // inner track loop
            
        }  //mother kink index <0
    } // Outer track loop
    cont->SetTitle("test");
    
    cont->MakeKinks();
    gEve->Redraw3D();
    return cont;
}

TEvePointSet* AliEveESDKinks::DrawPoints()
{
    TEvePointSet* points = new TEvePointSet("Kink vertex locations");
    
    AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();
    
    for (Int_t n=0; n<esd->GetNumberOfTracks(); ++n)
    {
        AliESDtrack* track = esd->GetTrack(n);
        if(track->GetKinkIndex(0)<0){
            
            AliESDkink *kink = esd->GetKink(TMath::Abs(track->GetKinkIndex(0))-1);
            const TVector3 Position(kink->GetPosition());
            points->SetNextPoint(Position.X(), Position.Y(), Position.Z());
            points->SetPointId(kink);
        }
    }

    points->SetTitle(Form("N=%d", points->Size()));
    points->SetMarkerStyle(4);
    points->SetMarkerSize(1.5);
    points->SetMarkerColor(kOrange+8);
    
    gEve->AddElement(points);
    gEve->Redraw3D();
    
    return points;
}

