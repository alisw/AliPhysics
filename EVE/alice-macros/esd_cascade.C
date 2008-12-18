// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

void esd_cascade_init_rectrack(TEveRecTrack& rt, AliExternalTrackParam* tp)
{
  Double_t      pbuf[3], vbuf[3];

  rt.fSign = tp->GetSign();
  tp->GetXYZ(vbuf);     rt.fV.Set(vbuf);
  tp->GetPxPyPz(pbuf);  rt.fP.Set(pbuf);
  // Double_t ep = at->GetP(), mc = at->GetMass();
  rt.fBeta = 1; // ep/TMath::Sqrt(ep*ep + mc*mc);
}

AliEveCascade* esd_make_cascade(TEveTrackPropagator* rnrStyle, AliESDVertex* primVtx,
				AliESDtrack* bac, AliESDcascade* cascade, Int_t i)
{
  TEveRecTrack   rcPos;
  TEveRecTrack   rcNeg;
  TEveRecV0      rcV0;

  TEveRecTrack   rcBac;
  TEveRecCascade rcCascade;

  Double_t pNeg[3], pPos[3], pV0[3];
  cascade->GetNPxPyPz(pNeg[0], pNeg[1], pNeg[2]);
  rcV0.fPPos.Set(pPos);
  cascade->GetPPxPyPz(pPos[0], pPos[1], pPos[2]);
  rcV0.fPNeg.Set(pNeg);

  Double_t v[3];
  cascade->GetXYZ(v[0], v[1], v[2]);
  rcV0.fVCa.Set(v);

  cascade->GetParamN()->GetXYZ(v);  rcV0.fVNeg.Set(v);
  cascade->GetParamP()->GetXYZ(v);  rcV0.fVPos.Set(v);

  rcV0.fV0Birth.Set(primVtx->GetXv(), primVtx->GetYv(), primVtx->GetZv());

  Double_t pBac[3], pCascade[3], cv[21]={0.};
  cascade->GetBPxPyPz(pBac[0], pBac[1], pBac[2]);
  rcCascade.fPBac.Set(pBac);
  cascade->GetPxPyPz(pCascade[0], pCascade[1], pCascade[2]);

  cascade->GetXYZcascade(v[0], v[1], v[2]);
  rcCascade.fCascadeVCa.Set(v);

  rcCascade.fCascadeBirth.Set(primVtx->GetXv(), primVtx->GetYv(), primVtx->GetZv());

  // Simulation data not directly available in AliESDcascade
  // rcCascade.fDLabel = cascade->GetBindex();

  // Problem: two following lines are not possible: no GetParamB !!
  // cascade->GetParamB()->GetXYZ(v);  rcCascade.fVBac.Set(v); 
  // esd_cascade_init_rectrack(rcBac, cascade->GetParamB());
  // Solution: create an AliExternalTrackParam with null cv...
  AliExternalTrackParam *bParam = new AliExternalTrackParam(v,pBac,cv,cascade->Charge());
  esd_cascade_init_rectrack(rcBac,bParam);
  rcBac.fIndex = cascade->GetBindex();

  AliEveCascade* myCascade = new AliEveCascade(&rcBac, &rcV0, &rcCascade, rnrStyle);
  myCascade->SetElementName(Form("ESDcascade %d", i));
  myCascade->SetElementTitle(Form("DCA %f",
                             cascade->GetDcaXiDaughters()));
  myCascade->SetESDIndex(i);
  myCascade->SetDaughterDCA(cascade->GetDcaXiDaughters());

  return myCascade;
}


AliEveCascadeList* esd_cascade()
{
  AliESDEvent* esd = AliEveEventManager::AssertESD();

  AliESDVertex* primVertex = (AliESDVertex*) esd->GetPrimaryVertex();

  AliEveCascadeList* cont = new AliEveCascadeList("ESD cascade");
  cont->SetMainColor(3); // green
  TEveTrackPropagator* rnrStyle = cont->GetPropagator();
  rnrStyle->SetMagField( 0.1*esd->GetMagneticField() );

  gEve->AddElement(cont);

  Int_t count = 0;
  for (Int_t n=0; n<esd->GetNumberOfCascades(); ++n)
  {
    AliESDcascade *cascade = esd->GetCascade(n);

    Int_t bacInd = cascade->GetBindex();
    AliESDtrack* bacTr = esd->GetTrack(bacInd);

    AliEveCascade* myCascade = esd_make_cascade(rnrStyle, primVertex, bacTr, cascade, n);
    if (myCascade)
    {
      gEve->AddElement(myCascade, cont);
      ++count;
    }
  }

  cont->SetTitle("test");

  cont->MakeCascades();
  gEve->Redraw3D();

  return cont;
}
