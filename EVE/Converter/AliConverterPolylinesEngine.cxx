//
// Created by mgrochow on 9/2/15.
//

#include "AliConverterPolylinesEngine.h"
#include <AliEveMUONTrack.h>
#include <AliEveKink.h>
#include <AliExternalFormatConverter.h>


//void AliConverterPolylinesEngine::PopulateEventWithMuonTracks(AliMinimalisticEvent &event) const
//{
//    AliEveMUONTrackList* lt = new AliEveMUONTrackList("ESD-Tracks");
//    TEveRecTrack rt;
//    AliESDMuonTrack *mTrack;
//    AliEveMUONTrack *track = new AliEveMUONTrack(&rt, lt->GetPropagator());
//    for (Int_t muonTrack = 0; muonTrack < fESDEvent->GetNumberOfMuonTracks(); muonTrack++){
//        mTrack = fESDEvent->GetMuonTrack(muonTrack);
//        if (mTrack->GetNHit() == 0)
//            continue;
//
//        rt.fLabel = muonTrack;
//
//        track->MakeESDTrack(mTrack);
//        lt->AddElement(track);
//
//        Int_t charge = mTrack->Charge();
//        Double_t energy = mTrack->E();
//        Int_t parentID = AliMinimalisticTrack::fgkImaginaryParent;
//        Double_t signedPT = mTrack->Pt()*charge;
//        Double_t mass = mTrack->M();
//        Double_t PxPyPz[3]; mTrack->PxPyPz(PxPyPz);
//        Double_t startXYZ[3];
//        startXYZ[0] = mTrack->GetNonBendingCoorAtDCA();
//        startXYZ[1] = mTrack->GetBendingCoorAtDCA();
//        startXYZ[2] = mTrack->GetZ();
//        Double_t theta = mTrack->Theta();
//        Double_t phi = mTrack->Phi();
//    }
//    lt->HackMomentumLimits();
//    delete track;
//    delete lt;
//}
//
//void AliExternalFormatConverter::AddPolyLinesToKinkTrack(
//        Int_t kinkID, AliMinimalisticTrack &mTrack, AliMinimalisticTrack &dTrack) const
//{
//    AliEveKinkList* cont = new AliEveKinkList("ESD kink");
//
//    TEveTrackPropagator* rnrStyleMoth = cont->GetPropagatorMoth();
//    TEveTrackPropagator* rnrStyleDaugh = cont->GetPropagatorDaugh();
//
//    rnrStyleMoth->SetMagField(0.1 * fESDEvent->GetMagneticField());
//    rnrStyleDaugh->SetMagField(0.1 * fESDEvent->GetMagneticField());
//    rnrStyleDaugh->SetMaxR(520);
//
//    AliESDkink *kink = fESDEvent->GetKink(kinkID);
//    AliESDtrack* moth = fESDEvent->GetTrack(kink->GetIndex(0));
//    AliESDtrack* daug = fESDEvent->GetTrack(kink->GetIndex(1));
//
//    TEveRecTrack rcMoth;
//    TEveRecTrack rcDaug;
//    TEveRecKink rcKink;
//
//    rcKink.fPMother.Set(kink->GetMotherP());
//    rcKink.fPDaughter.Set(kink->GetDaughterP());
//    rcKink.fVKink.Set(kink->GetPosition());
//
//    Double_t pbuf[3], vbuf[3];
//
//    rcMoth.fSign = moth->GetTPCInnerParam()->GetSign();
//    moth->GetTPCInnerParam()->GetXYZ(vbuf); rcMoth.fV.Set(vbuf);
//    moth->GetTPCInnerParam()->GetPxPyPz(pbuf); rcMoth.fP.Set(pbuf);
//
//    rcDaug.fSign = daug->GetOuterParam()->GetSign();
//    rcDaug.fV.Set(rcKink.fVKink);
//    rcDaug.fP.Set(rcKink.fPDaughter);
//
//    AliEveKink* myKink = new AliEveKink(&rcMoth, &rcDaug, &rcKink, rnrStyleMoth,rnrStyleDaugh);
//    myKink->SetESDKinkIndex(kinkID);
//    gEve->AddElement(myKink, cont);
//    cont->MakeKinks();
//
//    std::vector<TEveVector4D> mPoints = rnrStyleMoth->GetPoints();
//    std::vector<TEveVector4D> dPoints = rnrStyleDaugh->GetPoints();
//    std::cout << "Przed funkcja wielkosc: " << mPoints.size() << std::endl;
//    InsertPolyPoints(mTrack, mPoints);
//    InsertPolyPoints(dTrack, dPoints);
//    delete myKink;
//    delete cont;
//}
//
//void AliExternalFormatConverter::AddPolyLinesToV0Track(
//        Int_t v0ID, AliMinimalisticTrack &negativeTrack, AliMinimalisticTrack &positiveTrack) const
//{
//    AliEveV0List* cont = new AliEveV0List("ESD v0");
//    TEveTrackPropagator* positivePropagator = cont->GetPropagatorPos();
//    TEveTrackPropagator* negativePropagator = cont->GetPropagatorNeg();
//    positivePropagator->SetMagField( 0.1*fESDEvent->GetMagneticField() );
//    negativePropagator->SetMagField( 0.1*fESDEvent->GetMagneticField() );
//
//    gEve->AddElement(cont);
//    AliESDv0 *v0 = fESDEvent->GetV0(v0ID);
//
//    TEveRecTrack rcPos;
//    TEveRecTrack rcNeg;
//    TEveRecV0 rcV0;
//
//    Double_t pbuf[3], vbuf[3];
//
//    rcNeg.fSign = v0->GetParamN()->GetSign();
//    v0->GetParamN()->GetXYZ(vbuf); rcNeg.fV.Set(vbuf);
//    v0->GetParamN()->GetPxPyPz(pbuf); rcNeg.fP.Set(pbuf);
//    rcPos.fSign = v0->GetParamP()->GetSign();
//    v0->GetParamP()->GetXYZ(vbuf); rcPos.fV.Set(vbuf);
//    v0->GetParamP()->GetPxPyPz(pbuf); rcPos.fP.Set(pbuf);
//
//    rcNeg.fIndex = v0->GetNindex();
//    rcPos.fIndex = v0->GetPindex();
//
//    AliEveV0* myV0 = new AliEveV0(&rcNeg, &rcPos, &rcV0, negativePropagator, positivePropagator);
//
//    myV0->SetESDIndex(v0ID);
//    myV0->SetOnFlyStatus(v0->GetOnFlyStatus());
//    myV0->SetDaughterDCA(v0->GetDcaV0Daughters());
//
//    gEve->AddElement(myV0, cont);
//
//    cont->MakeV0s();
//
//    std::vector<TEveVector4D> negPoints = negativePropagator->GetPoints();
//    std::vector<TEveVector4D> posPoints = positivePropagator->GetPoints();
//    InsertPolyPoints(negativeTrack, negPoints);
//    InsertPolyPoints(positiveTrack, posPoints);
//    delete myV0;
//    delete cont;
//}
//
//void AliExternalFormatConverter::AddPolylinesToCascade(
//        Int_t cascadeID,
//        AliMinimalisticTrack &negativeTrack,
//        AliMinimalisticTrack &positiveTrack,
//        AliMinimalisticTrack &bachelorTrack) const
//{
//    AliESDVertex* primVtx = (AliESDVertex*) fESDEvent->GetPrimaryVertex();
//
//    AliEveCascadeList* cont = new AliEveCascadeList("ESD cascade");
//    TEveTrackPropagator* rnrStyleBac = cont->GetPropagatorBac();
//    TEveTrackPropagator* rnrStyleNeg = cont->GetPropagatorNeg();
//    TEveTrackPropagator* rnrStylePos = cont->GetPropagatorPos();
//    rnrStyleBac->SetMagField( 0.1 * fESDEvent->GetMagneticField() );
//    rnrStyleNeg->SetMagField( 0.1 * fESDEvent->GetMagneticField() );
//    rnrStylePos->SetMagField( 0.1 * fESDEvent->GetMagneticField() );
//
//    gEve->AddElement(cont);
//
//    AliESDcascade *cascade = fESDEvent->GetCascade(cascadeID);
//
//    TEveRecTrack rcPos,rcNeg,rcBac;
//    TEveRecV0 rcV0;
//    TEveRecCascade rcCascade;
//
//    Double_t v[3],pBac[3]={0.}, pNeg[3]={0.}, pPos[3]={0.}, cv[21]={0.},pbuf[3], vbuf[3];
//
//    cascade->GetBPxPyPz(pBac[0], pBac[1], pBac[2]);
//    cascade->GetNPxPyPz(pNeg[0], pNeg[1], pNeg[2]);
//    cascade->GetPPxPyPz(pPos[0], pPos[1], pPos[2]);
//    cascade->GetXYZcascade(v[0], v[1], v[2]);
//
//    AliExternalTrackParam *bacParam = new AliExternalTrackParam(v,pBac,cv,cascade->Charge());
//
//    rcBac.fSign = bacParam->GetSign();
//    rcBac.fIndex = cascade->GetBindex();
//    bacParam->GetXYZ(vbuf); rcBac.fV.Set(vbuf);
//    bacParam->GetPxPyPz(pbuf); rcBac.fP.Set(pbuf);
//
//    rcNeg.fSign = cascade->GetParamN()->GetSign();
//    rcNeg.fIndex = cascade->GetNindex();
//    cascade->GetParamN()->GetXYZ(vbuf); rcNeg.fV.Set(vbuf);
//    cascade->GetParamN()->GetPxPyPz(pbuf); rcNeg.fP.Set(pbuf);
//
//    rcPos.fSign = cascade->GetParamP()->GetSign();
//    rcPos.fIndex = cascade->GetPindex();
//    cascade->GetParamP()->GetXYZ(vbuf); rcPos.fV.Set(vbuf);
//    cascade->GetParamP()->GetPxPyPz(pbuf); rcPos.fP.Set(pbuf);
//
//    AliEveCascade* myCascade = new AliEveCascade(&rcBac, &rcNeg, &rcPos, &rcV0, &rcCascade, rnrStyleBac, rnrStyleNeg, rnrStylePos);
//
//    myCascade->SetESDIndex(cascadeID);
//    myCascade->SetDaughterDCA(cascade->GetDcaXiDaughters());
//    myCascade->SetLambdaP( pNeg[0]+pPos[0], pNeg[1]+pPos[1], pNeg[2]+pPos[2] );
//    myCascade->SetBachP( pBac[0], pBac[1], pBac[2]);
//
//    gEve->AddElement(myCascade, cont);
//    cont->MakeCascades();
//
//    std::vector<TEveVector4D> negPoints = rnrStyleNeg->GetPoints();
//    std::vector<TEveVector4D> posPoints = rnrStylePos->GetPoints();
//    std::vector<TEveVector4D> bacPoints = rnrStyleBac->GetPoints();
//    InsertPolyPoints(negativeTrack, negPoints);
//    InsertPolyPoints(positiveTrack, posPoints);
//    InsertPolyPoints(bachelorTrack, bacPoints);
//    delete myCascade;
//    delete bacParam;
//    delete cont;
//}
//
//void AliExternalFormatConverter::InsertPolyPoints(
//        AliMinimalisticTrack &Track, std::vector<TEveVector4D> &Points) const
//{
//    for(std::vector<TEveVector4D>::iterator iter = Points.begin(); iter != Points.end(); ++iter){
//        Track.AddPolyPoint(*iter);
//    }
//}
//
//void AliExternalFormatConverter::AddPolylinesToMuonTracks(
//        Int_t trackNumber, AliMinimalisticTrack &minimalisticMuonTrack) const
//{
//    // setup CDB
//    AliCDBManager* cdb = AliCDBManager::Instance();
//    cdb->SetDefaultStorage(Form("local://%s/../src/OCDB",gSystem->Getenv("ALICE_ROOT")));
//    if (!cdb->IsDefaultStorageSet()){
//        std::cerr<<"\n\nCould not set CDB.\n\n"<<std::endl;
//        return;
//    }
//    cdb->SetRun(fESDEvent->GetRunNumber());
//
//    fESDEvent->InitMagneticField();
//    if (!AliMUONESDInterface::GetTracker())
//        AliMUONESDInterface::ResetTracker(AliMUONCDB::LoadRecoParam());
//    AliMpCDB::LoadAll(kFALSE);
//    AliMUONESDInterface data;
//    data.LoadEvent(*fESDEvent);
//    AliESDMuonTrack *muonTrack = fESDEvent->GetMuonTrack(trackNumber);
//
//    TEveTrackList* eveTrackList = new TEveTrackList("MUON list");
//    std::string trackType;
//
//    AliMagF *field=NULL;
//
//    if (TGeoGlobalMagField::Instance()->GetField())
//        field = dynamic_cast<AliMagF*>(TGeoGlobalMagField::Instance()->GetField());
//    if(!field)
//        std::cerr<<"\n\n\nCouldn't set magnetic field!!\n\n\n"<<std::endl;
//
//    TEveTrackPropagator *trkProp = eveTrackList->GetPropagator();
//
//    AliEveMagField *aliEveMagField = new AliEveMagField(field);
//    trkProp->SetMagFieldObj(aliEveMagField);
//    trkProp->SetStepper(TEveTrackPropagator::kRungeKutta);
//    trkProp->SetFitDaughters(kFALSE);
//    trkProp->SetFitReferences(kTRUE);
//    trkProp->SetFitDecay(kFALSE);
//    trkProp->SetFitCluster2Ds(kFALSE);
//    trkProp->SetRnrReferences(kTRUE);
//    trkProp->SetMaxR(1000);
//
//    TObjArray *trackParams = nullptr;
//
//    if (muonTrack->ContainTrackerData() )
//    {
//        trackParams = data.FindTrack(muonTrack->GetUniqueID())->GetTrackParamAtCluster();
//        if(muonTrack->GetMatchTrigger() > 0){
//            minimalisticMuonTrack.SetTrackType(AliMinimalisticTrack::kMuonMatched);
//            trkProp->SetMaxZ(-AliMUONConstants::DefaultChamberZ(13)+10.);
//        } else {
//            minimalisticMuonTrack.SetTrackType(AliMinimalisticTrack::kMuonNotMatched);
//            trkProp->SetMaxZ(-AliMUONConstants::MuonFilterZBeg());
//        }
//    } else {
//        minimalisticMuonTrack.SetTrackType(AliMinimalisticTrack::kMuonGhost);;
//        trkProp->SetMaxZ(-AliMUONConstants::DefaultChamberZ(13)+10.);
//    }
//    eveTrackList->SetMainColor(kGreen);
//
//    TEveRecTrack recTrack;
//    TEveTrack* track;
//
//    if (muonTrack->ContainTrackerData()) {
//        track = new AliEveTrack(&recTrack, eveTrackList->GetPropagator());
//        track->SetSourceObject(muonTrack);
//
//        TIter next(trackParams);
//        AliMUONTrackParam* param;
//        while ((param = static_cast<AliMUONTrackParam*>(next()))){
//            TEveVector v(param->GetNonBendingCoor(),param->GetBendingCoor(),param->GetZ());
//            TEveVector p(param->Px(),param->Py(),param->Pz());
//            track->AddPathMark(TEvePathMark(TEvePathMark::kReference, v, p));
//        }
//        if (muonTrack->ContainTriggerData())
//        {
//            AliMUONTriggerCircuit* gTriggerCircuit = nullptr;
//            if (AliGeomManager::GetGeometry() == 0)
//            {
//                gGeoManager = 0;
//
//                AliGeomManager::LoadGeometry();
//                if (!AliGeomManager::GetGeometry())
//                {
//                    std::cout<<"\n\nCan not load geometry.\n\n"<<std::endl;
//                }
//                if (!AliGeomManager::ApplyAlignObjsFromCDB("ITS TPC TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO ACORDE"))
//                {
//                    std::cout<<"\n\nMismatch of alignable volumes. Proceeding.\n\n"<<std::endl;
//                }
//                AliGeomManager::GetGeometry()->DefaultColors();
//            }
//            gGeoManager = AliGeomManager::GetGeometry();
//
//            AliMUONGeometryTransformer* fMUONGeometryTransformer = new AliMUONGeometryTransformer();
//            fMUONGeometryTransformer->LoadGeometryData();
//            gTriggerCircuit = new AliMUONTriggerCircuit(fMUONGeometryTransformer);
//
//            Double_t x11 = gTriggerCircuit->GetX11Pos(muonTrack->LoCircuit(), muonTrack->LoStripY());
//            Double_t y11 = gTriggerCircuit->GetY11Pos(muonTrack->LoCircuit(), muonTrack->LoStripX());
//            Double_t z11 = gTriggerCircuit->GetZ11Pos(muonTrack->LoCircuit(), muonTrack->LoStripX());
//            Double_t y21 = gTriggerCircuit->GetY21Pos(muonTrack->LoCircuit(), muonTrack->LoStripX()+muonTrack->LoDev()+1);
//            Double_t z21 = gTriggerCircuit->GetZ21Pos(muonTrack->LoCircuit(), muonTrack->LoStripX()+muonTrack->LoDev()+1);
//            Double_t pz  = -muonTrack->PUncorrected(); // max value
//            TEveVector v(x11, y11, z11);
//            TEveVector p(pz*x11/z11, pz*(y21-y11)/(z21-z11), pz);
//            track->AddPathMark(TEvePathMark(TEvePathMark::kReference,v,p));
//            delete gTriggerCircuit;
//            delete fMUONGeometryTransformer;
//        }
//
//    } else {// ghost tracks (trigger only)
//        recTrack.fStatus = 0;
//        recTrack.fSign = muonTrack->Charge();
//        Double_t z11 = (muonTrack->GetZUncorrected() < -1.) ? muonTrack->GetZUncorrected() : (Double_t)AliMUONConstants::DefaultChamberZ(10);
//        recTrack.fV.Set(muonTrack->GetNonBendingCoorUncorrected(),muonTrack->GetBendingCoorUncorrected(),z11);
//        recTrack.fP.Set(-TMath::Tan(muonTrack->GetThetaXUncorrected()),-TMath::Tan(muonTrack->GetThetaYUncorrected()),-1.);
//
//        track = new AliEveTrack(&recTrack,eveTrackList->GetPropagator());
//        track->SetSourceObject(muonTrack);
//    }
//    track->SetAttLineAttMarker(eveTrackList);
//    eveTrackList->AddElement(track);
//    eveTrackList->MakeTracks();
//
//    gEve->AddElement(eveTrackList);
//
//    std::vector<TEveVector4D > muonPoints = eveTrackList->GetPropagator()->GetPoints();
//
//    InsertPolyPoints(minimalisticMuonTrack, muonPoints);
//    gEve->Redraw3D();
//    delete track;
//    delete eveTrackList;
//    delete aliEveMagField;
//}
