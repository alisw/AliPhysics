// Launches the list analyser and loads kinks into the list.
// If you already have a list (or rather a list analyser) with objects and you want to add the kinks to this list, do the following:
// Right-click the list analyser in the eve-browser and select "ExportToCint". Choose e.g. "list" for the name in the window that pops up.
// In the console type ".x ana_list_load_kinks.C(list)"
// For more information please see "ana_list.C" or have a look at the class documentation.


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGLViewer.h>
#include <TEveManager.h>
#include <TEveTrackPropagator.h>

#include <AliTRDarrayADC.h>
#include <AliTRDReconstructor.h>
#include <AliTRDtrackV1.h>
#include <AliESDkink.h>
#include <AliESDEvent.h>
#include <AliESDfriend.h>
#include <AliEveTRDData.h>
#include <AliEveEventManager.h>
#include <AliEveKink.h>
#include <AliEveListAnalyser.h>
#endif

void esd_kink_init_rectrackMother(TEveRecTrack& rt, const AliExternalTrackParam* tp)
{
  Double_t pbuf[3], vbuf[3];

  rt.fSign = tp->GetSign();
  tp->GetXYZ(vbuf);     rt.fV.Set(vbuf);
  tp->GetPxPyPz(pbuf);  rt.fP.Set(pbuf);
  
  rt.fBeta = 1; // ep/TMath::Sqrt(ep*ep + mc*mc);
}

void esd_kink_init_rectrackDaughter(TEveRecTrack& rt, const AliExternalTrackParam* tp, TEveVector* svt,TEveVector* spt)
{
  rt.fSign = tp->GetSign();
  rt.fV.Set(*svt);
  rt.fP.Set(*spt);
  
  rt.fBeta = 1; // ep/TMath::Sqrt(ep*ep + mc*mc);
}

AliEveKink* esd_make_kink(TEveTrackPropagator* rnrStyleMoth,TEveTrackPropagator* rnrStyleDaugh, AliESDtrack* moth, AliESDtrack* daug, AliESDkink* kink, Int_t i)
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

  moth->GetTPCInnerParam()->GetXYZ(r2);  rcKink.fVMother.Set(r2);
  daug->GetOuterParam()->GetXYZ(r);  rcKink.fVDaughter.Set(r);

  esd_kink_init_rectrackMother(rcMoth, (moth->GetTPCInnerParam()));  
  rcMoth.fIndex = moth->GetID();
  
  esd_kink_init_rectrackDaughter(rcDaug, daug->GetOuterParam(), &rcKink.fVKink, &rcKink.fPDaughter);
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

void ana_list_load_kinks(AliEveListAnalyser* objects = 0, TEveElement *cont = 0)
{
  Bool_t noListProvided = kTRUE;
  if (objects)  noListProvided = kFALSE;

  // Link data containers
  AliESDfriend *eventESDfriend = 0x0;
  if(!(eventESDfriend = AliEveEventManager::AssertESDfriend())){
    Warning("ana_list_load_tracks", "AliESDfriend not found");
    return;
  }

  AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();

  AliEveEventManager::Instance()->AssertGeometry();

  if (!objects) objects = new AliEveListAnalyser("Analysis Objects");
  
  // Kinks
  AliEveKinkList* list = new AliEveKinkList("ESD kink");
  list->SetMainColor(3); // green
  TEveTrackPropagator* rnrStyleMoth = list->GetPropagatorMoth();
  rnrStyleMoth->SetMagField( 0.1*esd->GetMagneticField() );
  TEveTrackPropagator* rnrStyleDaugh = list->GetPropagatorDaugh();
  rnrStyleDaugh->SetMagField( 0.1*esd->GetMagneticField() ); 
  rnrStyleDaugh->SetMaxR(520);
  //gEve->AddElement(list);

// Number of elements already in the list
  Int_t nOld = objects->NumChildren();

  Int_t count = 0;

  for (Int_t n = 0; n < esd->GetNumberOfTracks(); ++n)
  { 
    AliESDtrack* mtrack = esd->GetTrack(n);
    if (!mtrack) continue;
    // GetKinkIndex < 0 in the following!!! => -GetKinkIndex > 0
    if(mtrack->GetKinkIndex(0) < 0)
    {
      AliESDkink *kink = new AliESDkink;
    
      kink = esd->GetKink(-mtrack->GetKinkIndex(0) - 1);
      if (!kink)  continue;      

      for (Int_t m = 0; m < esd->GetNumberOfTracks(); ++m)
      { 
        AliESDtrack * dtrack = esd->GetTrack(m);
        if (!dtrack)  continue;

        if((dtrack->GetKinkIndex(0) > 0 ) && (dtrack->GetKinkIndex(0) == -mtrack->GetKinkIndex(0))) 
        {
          AliESDtrack* mothTr = esd->GetTrack(n);
          AliESDtrack* daugTr = esd->GetTrack(m);
          if (!mothTr || !daugTr) continue;

          AliEveKink* myKink = esd_make_kink(rnrStyleMoth, rnrStyleDaugh, mothTr, daugTr, kink, (-mtrack->GetKinkIndex(0)-1));
          if (myKink)
          {
            myKink->SetUserData(kink);
            myKink->SetMarkerColor(3);
            gEve->AddElement(myKink, list);
            objects->AddElement(myKink);
            myKink->SetName(Form("[%4d] kink", count + nOld));
            ++count;
          }
        }
      }  // inner track loop
    }  //mother kink index <0
  } // Outer track loop

  list->MakeKinks();

  objects->SetTitle(Form("Objects %d", objects->NumChildren()));
  objects->StampObjProps();
 
  // If a new list analyser has been created, add it to eve.
  if (noListProvided)  gEve->AddElement(objects, cont);

  gEve->Redraw3D();

  //  TGLViewer *v = gEve->GetDefaultGLViewer();
  //  v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  //  ((TGLOrthoCamera&)v->CurrentCamera()).SetEnableRotate(kTRUE);
  //  v->UpdateScene();

  return;
}
