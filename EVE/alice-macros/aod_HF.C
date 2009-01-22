// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

void esd_v0_init_rectrack(TEveRecTrack& rt, AliExternalTrackParam* tp)
{
  Double_t      pbuf[3], vbuf[3];

  rt.fSign = tp->GetSign();
  tp->GetXYZ(vbuf);     rt.fV.Set(vbuf);
  tp->GetPxPyPz(pbuf);  rt.fP.Set(pbuf);
  // Double_t ep = at->GetP(), mc = at->GetMass();
  rt.fBeta = 1; // ep/TMath::Sqrt(ep*ep + mc*mc);
}

AliEveV0* aod_make_HF(TEveTrackPropagator* rnrStyle, AliESDVertex* primVtx,
		      AliESDtrack* neg, AliESDtrack* pos, AliAODRecoDecayHF* rd, Int_t i)
{
  TEveRecTrack  rcPos;
  TEveRecTrack  rcNeg;
  TEveRecV0     rcV0;

  Double_t p[3];
  p[0]=rd->PxProng(0);  p[1]=rd->PyProng(0);  p[2]=rd->PzProng(0);
  rcV0.fPPos.Set(p);
  p[0]=rd->PxProng(1);  p[1]=rd->PyProng(1);  p[2]=rd->PzProng(1);
  rcV0.fPNeg.Set(p);

  p[0]=rd->Px();  p[1]=rd->Py();  p[2]=rd->Pz();

  Double_t v[3] = {rd->Xv(),rd->Yv(),rd->Zv()}
  printf("vertex %f %f %f\n",v[0],v[1],v[2]);
  rcV0.fVCa.Set(v);

  neg->GetXYZ(v);  rcV0.fVNeg.Set(v);
  pos->GetXYZ(v);  rcV0.fVPos.Set(v);

  rcV0.fV0Birth.Set(primVtx->GetXv(), primVtx->GetYv(), primVtx->GetZv());

  // Simulation data not directly available in AliESDv0
  //rcV0.fDLabel[0] = v0->GetNindex();
  //rcV0.fDLabel[1] = v0->GetPindex();

  esd_v0_init_rectrack(rcNeg, neg);
  //rcNeg.fIndex = v0->GetNindex();
  esd_v0_init_rectrack(rcPos, pos);
  //rcPos.fIndex = v0->GetPindex();

  AliEveV0* myD0 = new AliEveV0(&rcNeg, &rcPos, &rcV0, rnrStyle);
  myD0->SetElementName(Form("D0->Kpi %d", i));
  myD0->SetElementTitle(Form("CosPointingAngle %f",
                             rd->CosPointingAngle()));
  //myV0->SetESDIndex(i);
  //myV0->SetOnFlyStatus(v0->GetOnFlyStatus());
  //myV0->SetDaughterDCA(v0->GetDcaV0Daughters());

  return myD0;
}


AliEveV0List* aod_HF()
{
  AliAODEvent* aod = AliEveEventManager::AssertAOD();
  AliESDEvent* esd = AliEveEventManager::AssertESD();

  AliESDVertex* primVertex = (AliESDVertex*) esd->GetPrimaryVertex();

  AliEveV0List* cont = new AliEveV0List("AOD HF vertices");
  cont->SetMainColor(2);
  TEveTrackPropagator* rnrStyle = cont->GetPropagator();
  rnrStyle->SetMagField( 0.1*aod->GetMagneticField() );

  gEve->AddElement(cont);


  TEvePointSet* pointsD0toKpi = new TEvePointSet("D0->Kpi vertex locations");

  // load D0->Kpi candidates
  TClonesArray *arrayD0toKpi = 
    (TClonesArray*)aod->GetList()->FindObject("D0toKpi"); 
     
  // load 3prong candidates
  TClonesArray *array3Prong = 
    (TClonesArray*)aod->GetList()->FindObject("Charm3Prong"); 


  Int_t countD0 = 0;
  for (Int_t iD0toKpi=0; iD0toKpi<arrayD0toKpi->GetEntriesFast(); iD0toKpi++) {
    AliAODRecoDecayHF2Prong *rd = (AliAODRecoDecayHF2Prong*)arrayD0toKpi->UncheckedAt(iD0toKpi);
    Bool_t unsetvtx=kFALSE;
    if(!rd->GetOwnPrimaryVtx()) {
      rd->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
      unsetvtx=kTRUE;
    }

    AliESDtrack* negTr = esd->GetTrack(rd->GetProngID(0));
    AliESDtrack* posTr = esd->GetTrack(rd->GetProngID(1));

    AliVVertex *secv = rd->GetSecondaryVtx();
    negTr->PropagateToDCA((AliESDVertex*)secv,aod->GetMagneticField(),100.);
    posTr->PropagateToDCA((AliESDVertex*)secv,aod->GetMagneticField(),100.);

    AliEveV0* myD0 = aod_make_HF(rnrStyle,primVertex,negTr,posTr,rd,iD0toKpi);
    if (myD0) {
      gEve->AddElement(myD0,cont);
      countD0++;
    }

    pointsD0toKpi->SetNextPoint(rd->Xv(),rd->Yv(),rd->Zv());
    pointsD0toKpi->SetPointId(rd);

    if(unsetvtx) rd->UnsetOwnPrimaryVtx();
  }

  cont->SetTitle("test");

  cont->MakeV0s();
  gEve->Redraw3D();

  pointsD0toKpi->SetTitle(Form("N=%d",pointsD0toKpi->Size()));
  pointsD0toKpi->SetMarkerStyle(4);
  pointsD0toKpi->SetMarkerSize(1.5);
  pointsD0toKpi->SetMarkerColor(kViolet);

  gEve->AddElement(pointsD0toKpi);
  gEve->Redraw3D();

  return cont;
}
















