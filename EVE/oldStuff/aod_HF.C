// $Id$
// Main author: Davide Caffarri 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TClonesArray.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TEveVSDStructs.h>
#include <TEveTrackPropagator.h>
#include <TEvePointSet.h>
#include <TEveManager.h>
#include <TEveUtil.h>

#include <AliExternalTrackParam.h>
#include <AliVVertex.h>
#include <AliAODVertex.h>
#include <AliAODEvent.h>
#include <AliAODTrack.h>
#include <AliAODMCParticle.h>
#include <AliESDtrack.h>
#include <AliESDEvent.h>
#include <PWGHF/vertexingHF/AliAODRecoDecayHF.h>
#include <PWGHF/vertexingHF/AliAODRecoDecayHF2Prong.h>
#include <PWGHF/vertexingHF/macros/LoadLibraries.C>
#include <AliEveHF.h>
#include <AliEveEventManager.h>
#else
class AliAODRecoDecayHF;
#endif

void aod_hf_init_rectrack(TEveRecTrack& rt, AliExternalTrackParam* tp)
{
  Double_t      pbuf[3], vbuf[3];
  rt.fSign = tp->GetSign();
  tp->GetXYZ(vbuf);     rt.fV.Set(vbuf);
  tp->GetPxPyPz(pbuf);  rt.fP.Set(pbuf);
  // Double_t ep = at->GetP(), mc = at->GetMass();
  rt.fBeta = 1; // ep/TMath::Sqrt(ep*ep + mc*mc);
}


AliEveHF* aod_make_HF(TEveTrackPropagator* rnrStyle, AliAODVertex* primAODVtx,
		      AliESDtrack* neg, AliESDtrack* pos, AliAODRecoDecayHF* rd, Int_t i)
{
  TEveRecTrack  rcPos;
  TEveRecTrack  rcNeg;
  //TEveRecV0     rcV0;

  /*
    Double_t p[3];
    p[0]=rd->PxProng(0);  p[1]=rd->PyProng(0);  p[2]=rd->PzProng(0);
    rcV0.fPPos.Set(p);
    p[0]=rd->PxProng(1);  p[1]=rd->PyProng(1);  p[2]=rd->PzProng(1);
    rcV0.fPNeg.Set(p);

    p[0]=rd->Px();  p[1]=rd->Py();  p[2]=rd->Pz();
  */

  Double_t primVtx[3]={primAODVtx->GetX(), primAODVtx->GetY(), primAODVtx->GetZ()};


  Double_t v[3] = {rd->Xv(),rd->Yv(),rd->Zv()};
  printf("vertex %f %f %f\n",v[0],v[1],v[2]);

  aod_hf_init_rectrack(rcNeg, neg);
  //rcNeg.fIndex = v0->GetNindex();
  aod_hf_init_rectrack(rcPos, pos);
  //rcPos.fIndex = v0->GetPindex();

  AliEveHF* myHF = new AliEveHF(&rcNeg, &rcPos, primVtx ,rd, rnrStyle);
  myHF->SetElementName(Form("D0->Kpi %d", i));
  myHF->SetElementTitle(Form("CosPointingAngle %f", rd->CosPointingAngle()));
  myHF->SetAODIndex(i);
  return myHF;
}

AliEveHFList* aod_HF()
{
  Bool_t useParFiles=kFALSE;
  
//  TEveUtil::LoadMacro("$ALICE_ROOT/PWGHF/vertexingHF/macros/LoadLibraries.C+");
//  LoadLibraries(useParFiles);
  
  AliAODEvent* aod = AliEveEventManager::AssertAOD();
  AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();

  /*
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW");
    gSystem->Load("libPWG3base");
    gSystem->Load("libPWG3vertexingHF");
  */

  // load MC particles
  TClonesArray *mcArray =
    (TClonesArray*) aod->FindListObject(AliAODMCParticle::StdBranchName());
  if (!mcArray) {
    printf("MC particles branch not found!\n");
    return 0;
  }

  AliAODVertex* primVtx_aod = (AliAODVertex*) aod->GetPrimaryVertex();
  // AliESDVertex *primVtx_esd = (AliESDVertex*) esd->GetPrimaryVertex();

  AliEveHFList* cont = new AliEveHFList("AOD HF vertices");
  cont->SetMainColor(2);
  TEveTrackPropagator* rnrStyle = cont->GetPropagator();
  rnrStyle->SetMagField( 0.1*aod->GetMagneticField() );

  gEve->AddElement(cont);

  TEvePointSet* pointsD0toKpi = new TEvePointSet("D0->Kpi vertex locations");

  // load D0->Kpi candidates
  TClonesArray *arrayD0toKpi = (TClonesArray*) aod->FindListObject("D0toKpi");

  // load 3prong candidates
  // TClonesArray *array3Prong =
  // (TClonesArray*)aod->GetList()->FindObject("Charm3Prong");

  Int_t countD0 = 0;
  for (Int_t iD0toKpi=0; iD0toKpi<arrayD0toKpi->GetEntriesFast(); iD0toKpi++)
  {
    AliAODRecoDecayHF2Prong *rd = (AliAODRecoDecayHF2Prong*)arrayD0toKpi->UncheckedAt(iD0toKpi);
    Bool_t unsetvtx=kFALSE;
    if (!rd->GetOwnPrimaryVtx()) {
      rd->SetOwnPrimaryVtx(primVtx_aod);
      unsetvtx=kTRUE;
    }
    // REAL D0 particle. If you want to draw only real D0 un-comment these lines
    //Int_t labD0 = rd->MatchToMC(421,mcArray);
    //if(labD0<0) continue;

    AliAODTrack *negAODTr = dynamic_cast<AliAODTrack *>(rd->GetDaughter(0));
    AliAODTrack *posAODTr = dynamic_cast<AliAODTrack *>(rd->GetDaughter(1));

    AliVVertex  *secv = rd->GetSecondaryVtx();

    AliESDtrack *negTr = new AliESDtrack(negAODTr);
    AliESDtrack *posTr = new AliESDtrack(posAODTr);

    negTr->PropagateToDCA((AliAODVertex*)secv,aod->GetMagneticField(),100.);
    posTr->PropagateToDCA((AliAODVertex*)secv,aod->GetMagneticField(),100.);

    AliEveHF* myD0 = aod_make_HF(rnrStyle,primVtx_aod,negTr,posTr,rd,iD0toKpi);
    if (myD0) {
      gEve->AddElement(myD0,cont);
      countD0++;
    }

    pointsD0toKpi->SetNextPoint(rd->Xv(),rd->Yv(),rd->Zv());
    pointsD0toKpi->SetPointId(rd);

    if(unsetvtx) {
      rd->UnsetOwnPrimaryVtx();
    }
  }

  //cont->SetTitle("test");

  cont->MakeHFs();
  gEve->Redraw3D();

  pointsD0toKpi->SetTitle(Form("N=%d", pointsD0toKpi->Size()));
  pointsD0toKpi->SetMarkerStyle(4);
  pointsD0toKpi->SetMarkerSize(1.5);
  pointsD0toKpi->SetMarkerColor(kViolet);

  gEve->AddElement(pointsD0toKpi);
  gEve->Redraw3D();

  return cont;
}
