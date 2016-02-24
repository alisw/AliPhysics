// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#include <AliEveESDCascades.h>

#include <TMath.h>
#include <TEveManager.h>
#include <TEveVSDStructs.h>
#include <TEveTrackPropagator.h>


#include <AliESDEvent.h>
#include <AliESDVertex.h>


#include <AliEveEventManager.h>

AliEveCascadeList* AliEveESDCascades::Draw()
{    
    AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();
    
    AliESDVertex* primVertex = (AliESDVertex*) esd->GetPrimaryVertex();
    
    AliEveCascadeList* cont = new AliEveCascadeList("ESD cascade");
    cont->SetMainColor(kBlue+2);
    TEveTrackPropagator* rnrStyleBac = cont->GetPropagatorBac();
    TEveTrackPropagator* rnrStyleNeg = cont->GetPropagatorNeg();
    TEveTrackPropagator* rnrStylePos = cont->GetPropagatorPos();
    rnrStyleBac->SetMagField( 0.1*esd->GetMagneticField() );
    rnrStyleNeg->SetMagField( 0.1*esd->GetMagneticField() );
    rnrStylePos->SetMagField( 0.1*esd->GetMagneticField() );
    
    gEve->AddElement(cont);
    
    Int_t count = 0;
    for (Int_t n=0; n<esd->GetNumberOfCascades(); ++n)
    {
        AliESDcascade *cascade = esd->GetCascade(n);
        
        Int_t bacInd = cascade->GetBindex();
        AliESDtrack* bacTr = esd->GetTrack(bacInd);
        
        AliEveCascade* myCascade = MakeCascade(rnrStyleBac,rnrStyleNeg,rnrStylePos, primVertex, bacTr, cascade, n);
        if (myCascade)
        {
            gEve->AddElement(myCascade, cont);
            ++count;
        }
    }
    
    cont->SetTitle("Cascade candidates (reco)");
    
    cont->MakeCascades();
    gEve->Redraw3D();
    
    return cont;
}

void AliEveESDCascades::InitRecTrack(TEveRecTrack& rt, const AliExternalTrackParam* tp)
{
    Double_t      pbuf[3], vbuf[3];
    
    rt.fSign = tp->GetSign();
    tp->GetXYZ(vbuf);     rt.fV.Set(vbuf);
    tp->GetPxPyPz(pbuf);  rt.fP.Set(pbuf);
    // Double_t ep = at->GetP(), mc = at->GetMass();
    rt.fBeta = 1; // ep/TMath::Sqrt(ep*ep + mc*mc);
}

AliEveCascade* AliEveESDCascades::MakeCascade(TEveTrackPropagator* rnrStyleBac,TEveTrackPropagator* rnrStyleNeg,TEveTrackPropagator* rnrStylePos, AliESDVertex* primVtx,AliESDtrack* bac, AliESDcascade* cascade, Int_t i)
{
    TEveRecTrack   rcPos;
    TEveRecTrack   rcNeg;
    TEveRecV0      rcV0;
    
    TEveRecTrack   rcBac;
    TEveRecCascade rcCascade;
    
    Double_t v[3];
    cascade->GetXYZ(v[0], v[1], v[2]);
    rcV0.fVCa.Set(v);
    
    cascade->GetParamN()->GetXYZ(v);  rcV0.fVNeg.Set(v);
    cascade->GetParamP()->GetXYZ(v);  rcV0.fVPos.Set(v);
    
    rcV0.fV0Birth.Set(primVtx->GetX(), primVtx->GetY(), primVtx->GetZ());
    
    Double_t pCascade[3]={0.}, pBac[3]={0.}, pNeg[3]={0.}, pPos[3]={0.}, cv[21]={0.};
    //cascade->GetPxPyPz(pCascade[0], pCascade[1], pCascade[2]);
    cascade->GetBPxPyPz(pBac[0], pBac[1], pBac[2]);
    cascade->GetNPxPyPz(pNeg[0], pNeg[1], pNeg[2]);
    cascade->GetPPxPyPz(pPos[0], pPos[1], pPos[2]);
    
    rcCascade.fPBac.Set(pBac);
    rcV0.fPNeg.Set(pNeg);
    rcV0.fPPos.Set(pPos);

    Double_t pLambda = TMath::Sqrt((pNeg[0]+pPos[0])* (pNeg[0]+pPos[0]) +
                                   (pNeg[1]+pPos[1])* (pNeg[1]+pPos[1]) +
                                   (pNeg[2]+pPos[2])* (pNeg[2]+pPos[2]) );
    Double_t pBach   = TMath::Sqrt( pBac[0]*pBac[0] + pBac[1]*pBac[1] + pBac[2]*pBac[2]);
    
    cascade->GetXYZcascade(v[0], v[1], v[2]);
    rcCascade.fCascadeVCa.Set(v);
    
    rcCascade.fCascadeBirth.Set(primVtx->GetX(), primVtx->GetY(), primVtx->GetZ());
    
    // Simulation data not directly available in AliESDcascade
    // rcCascade.fDLabel = cascade->GetBindex();
    
    // Problem: two following lines are not possible: no GetParamB !!
    // cascade->GetParamB()->GetXYZ(v);  rcCascade.fVBac.Set(v);
    // InitRecTrack(rcBac, cascade->GetParamB());
    // Solution: create an AliExternalTrackParam with null cv...
    AliExternalTrackParam *bParam = new AliExternalTrackParam(v,pBac,cv,cascade->Charge());
    InitRecTrack(rcBac,bParam);
    rcBac.fIndex = cascade->GetBindex();
    
    InitRecTrack(rcNeg, cascade->GetParamN());
    rcNeg.fIndex = cascade->GetNindex();
    InitRecTrack(rcPos, cascade->GetParamP());
    rcPos.fIndex = cascade->GetPindex();
    
    AliEveCascade* myCascade = new AliEveCascade(&rcBac, &rcNeg, &rcPos, &rcV0, &rcCascade, rnrStyleBac, rnrStyleNeg, rnrStylePos);
    myCascade->SetElementName(Form("ESDcascade %d", i));
    
    myCascade->SetElementTitle(Form("Info coming directly from AliESDcascade : \n - Charge : %d \n - Cascade decay position : x = %.4f, y = %.4f, z = %.4f, Transv. radius = %.4f cm, Decay Length = %.4f cm\n\n - Pt(Cascade) : %f GeV/c, Ptot(Cascade): %f GeV/c\n - Lambda : px = %.4f, py = %.4f, pz = %.4f, Ptot : %f GeV/c\n - Bach   : px = %.4f, py = %.4f, pz = %.4f, Ptot : %f GeV/c\n\n - Eta : %f\n - Phi : %f deg \n - Theta : %f deg\n - DCA : %f cm \n - Cos(Ptg Angle) : %f \n\n - Eff. mass (Xi hyp) : %f GeV/c2",
                                    cascade->Charge(),
                                    v[0], v[1], v[2], TMath::Sqrt(v[0]*v[0] +v[1]*v[1]), TMath::Sqrt(v[0]*v[0] +v[1]*v[1] + v[2]*v[2] ),
                                    cascade->Pt(), cascade->P(),
                                    pNeg[0]+pPos[0], pNeg[1]+pPos[1], pNeg[2]+pPos[2], pLambda,
                                    pBac[0], pBac[1], pBac[2], pBach,
                                    cascade->Eta(),
                                    cascade->Phi()   * 180/TMath::Pi(),
                                    cascade->Theta() * 180/TMath::Pi(),
                                    cascade->GetDcaXiDaughters(),
                                    cascade->GetCascadeCosineOfPointingAngle(primVtx->GetX(),
                                                                             primVtx->GetY(),
                                                                             primVtx->GetZ()),
                                    cascade->GetEffMassXi()
                                    )
                               );
    
    
    myCascade->SetESDIndex(i);
    myCascade->SetDaughterDCA(cascade->GetDcaXiDaughters());
    myCascade->SetLambdaP( pNeg[0]+pPos[0], pNeg[1]+pPos[1], pNeg[2]+pPos[2] );
    myCascade->SetBachP(   pBac[0], pBac[1], pBac[2]);
    return myCascade;
}

TEvePointSet* AliEveESDCascades::DrawPoints()
{
    TEvePointSet* points = new TEvePointSet("Cascade vertex locations");
    
    AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();
    
    Int_t NCascades = esd->GetNumberOfCascades();
    
    Double_t x, y, z;
    for (Int_t n = 0; n < NCascades; ++n)
    {
        AliESDcascade* av = esd->GetCascade(n);
        av->GetXYZcascade(x, y, z);
        points->SetNextPoint(x, y, z);
        points->SetPointId(av);
    }
    
    points->SetTitle(Form("N=%d", points->Size()));
    points->SetMarkerStyle(4);
    points->SetMarkerSize(1.5);
    points->SetMarkerColor(kMagenta-9);
    
    gEve->AddElement(points);
    gEve->Redraw3D();
    
    return points;
}
