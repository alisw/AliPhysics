//
// Created by mgrochow on 9/2/15.
//

#include "AliConverterPolylinesEngine.h"
#include <TSystem.h>
#include <TGeoManager.h>
#include <TEveTrackPropagator.h>
#include <TGeoGlobalMagField.h>
#include <TEveTrack.h>
#include <TEveManager.h>
#include <TInterpreter.h>

#include <AliEveTrack.h>
#include <AliEveMUONTrack.h>
#include <AliESDcascade.h>
#include <AliESDkink.h>
#include <AliESDMuonTrack.h>
#include <AliGeomManager.h>
#include <AliEveKink.h>
#include <AliEveCascade.h>
#include <AliCDBManager.h>
#include <AliMUONESDInterface.h>
#include <AliMUONTrack.h>
#include <AliMUONTrackExtrap.h>
#include <AliMUONTrackParam.h>
#include <AliMUONConstants.h>
#include <AliMUONCDB.h>
#include <AliMUONGeometryTransformer.h>
#include <AliMUONTriggerCircuit.h>
#include <AliMpCDB.h>
#include <AliMagF.h>
#include <AliEveMagField.h>
#include <AliEveV0.h>


void AliConverterPolylinesEngine::InitializeEngine(AliESDEvent *event)
{
    fESDEvent = event;
    fESDEvent->InitMagneticField();
    AssertMagField();
    AssertGeometry();
}

void AliConverterPolylinesEngine::AssertGeometry() const
{
    if (AliGeomManager::GetGeometry() == 0)
    {
        gGeoManager = 0;

        AliGeomManager::LoadGeometry();
        if ( ! AliGeomManager::GetGeometry())
        {
            std::cerr<<"\n\nCan not load geometry.\n\n"<<std::endl;
        }
        if ( ! AliGeomManager::ApplyAlignObjsFromCDB("ITS TPC TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO ACORDE"))
        {
            std::cerr<<"\n\nMismatch of alignable volumes. Proceeding.\n\n"<<std::endl;
        }
        AliGeomManager::GetGeometry()->DefaultColors();
    }
    gGeoManager = AliGeomManager::GetGeometry();
}

void AliConverterPolylinesEngine::AssertMagField() const
{
    // setup CDB
    AliCDBManager* cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage(Form("local://%s/../src/OCDB",gSystem->Getenv("ALICE_ROOT")));
    if (!cdb->IsDefaultStorageSet())
    {
        std::cerr<<"\n\nCould not set CDB.\n\n"<<std::endl;
        return;
    }
    cdb->SetRun(fESDEvent->GetRunNumber());
}

AliConverterPolylinesEngine::AliConverterPolylinesEngine()
{
    Char_t **dummy_1;
    Int_t dummy_2 = 0;
    fApp = new TRint("PolyLinesEnging", &dummy_2, dummy_1, 0, 0, kTRUE);
    TEveManager::Create(kFALSE);
}

AliConverterPolylinesEngine::~AliConverterPolylinesEngine()
{
    fApp->Terminate(0);
    if (fApp)
        delete fApp;
    TEveManager::Terminate();
}

void AliConverterPolylinesEngine::AddPolylinesToMuonTracks(
        Int_t trackNumber,
        AliMinimalisticTrack &minimalisticMuonTrack
) const
{
    if (!AliMUONESDInterface::GetTracker())
        AliMUONESDInterface::ResetTracker(AliMUONCDB::LoadRecoParam());
    AliMpCDB::LoadAll(kFALSE);
    AliMUONESDInterface data;
    data.LoadEvent(*fESDEvent);

    TEveTrackList* trackList = new TEveTrackList("MUON Tracks");
    trackList->SetRnrPoints(kFALSE);
    trackList->SetRnrLine(kTRUE);
    trackList->SetLineColor(kGreen);

    AliESDMuonTrack* emt = fESDEvent->GetMuonTrack(trackNumber);
    AliMagF *field=NULL;

    if (TGeoGlobalMagField::Instance()->GetField())
        field = dynamic_cast<AliMagF*>(TGeoGlobalMagField::Instance()->GetField());
    if(!field)
        std::cerr<<"\n\n\nCouldn't set magnetic field!!\n\n\n"<<std::endl;

    AliEveMagField *aliEveMagField = new AliEveMagField(field);

    TEveTrackPropagator *trkProp = trackList->GetPropagator();

    // set magnetic field
    if (AliMUONTrackExtrap::IsFieldON()){trkProp->SetMagFieldObj(aliEveMagField);}
    else{trkProp->SetMagField(0.0);}
    trkProp->SetStepper(TEveTrackPropagator::kRungeKutta);
    trkProp->SetMaxR(1000);
    trkProp->SetFitDaughters(kFALSE);
    trkProp->SetFitReferences(kTRUE);
    trkProp->SetFitDecay(kFALSE);
    trkProp->SetFitCluster2Ds(kFALSE);
    trkProp->SetRnrReferences(kTRUE);
    trkProp->RefPMAtt().SetMarkerSize(0.5);
    trkProp->RefPMAtt().SetMarkerColor(kGreen);

    TEveRecTrack recTrack;
    TEveTrack* track;

    // fill general info
    UInt_t trackId = emt->GetUniqueID();
    recTrack.fLabel = emt->GetLabel();
    recTrack.fIndex = (Int_t)trackId;

    if (emt->ContainTrackerData()) // in tracker
    {
        if(emt->GetMatchTrigger() > 0) { // matched
            trkProp->SetMaxZ(-AliMUONConstants::DefaultChamberZ(13)+10.);
            minimalisticMuonTrack.SetTrackType(kMuonMatched);
        }
        else { // not matched
            trkProp->SetMaxZ(-AliMUONConstants::MuonFilterZBeg());
            minimalisticMuonTrack.SetTrackType(kMuonNotMatched);
        }
        recTrack.fStatus = emt->GetMatchTrigger();
        recTrack.fSign = emt->Charge();
        recTrack.fV.Set(emt->GetNonBendingCoorAtDCA(),emt->GetBendingCoorAtDCA(),emt->GetZ());
        recTrack.fP.Set(emt->PxAtDCA(),emt->PyAtDCA(),emt->PzAtDCA());
        recTrack.fBeta = ( emt->E() > 0 ) ? emt->P()/emt->E() : 0;

        // produce eve track
        track = new AliEveTrack(&recTrack,trackList->GetPropagator());
        track->SetStdTitle();
        track->SetSourceObject(emt); // WARNING: Change the UniqueID of the object!!

        // add path mark
        TIter next(data.FindTrack(trackId)->GetTrackParamAtCluster());
        AliMUONTrackParam* param;
        while ( ( param = static_cast<AliMUONTrackParam*>(next()) ) )
        {
            TEveVector v(param->GetNonBendingCoor(),param->GetBendingCoor(),param->GetZ());
            TEveVector p(param->Px(),param->Py(),param->Pz());
            track->AddPathMark(TEvePathMark(TEvePathMark::kReference,v,p));
        }

        if (emt->ContainTriggerData())        // add trigger track if any
        {
            // load trigger circuit
            AliMUONTriggerCircuit* gTriggerCircuit = 0x0;
            AssertGeometry();
            AliMUONGeometryTransformer* fMUONGeometryTransformer = new AliMUONGeometryTransformer();
            fMUONGeometryTransformer->LoadGeometryData();
            gTriggerCircuit = new AliMUONTriggerCircuit(fMUONGeometryTransformer);
            Double_t x11 = gTriggerCircuit->GetX11Pos(emt->LoCircuit(), emt->LoStripY());
            Double_t y11 = gTriggerCircuit->GetY11Pos(emt->LoCircuit(), emt->LoStripX());
            Double_t z11 = gTriggerCircuit->GetZ11Pos(emt->LoCircuit(), emt->LoStripX());
            Double_t y21 = gTriggerCircuit->GetY21Pos(emt->LoCircuit(), emt->LoStripX()+emt->LoDev()+1);
            Double_t z21 = gTriggerCircuit->GetZ21Pos(emt->LoCircuit(), emt->LoStripX()+emt->LoDev()+1);
            Double_t pz  = -emt->PUncorrected(); // max value
            TEveVector v(x11, y11, z11);
            TEveVector p(pz*x11/z11, pz*(y21-y11)/(z21-z11), pz);
            track->AddPathMark(TEvePathMark(TEvePathMark::kReference,v,p));
        }
    }
    else // ghost
    {
        trkProp->SetMaxZ(-AliMUONConstants::DefaultChamberZ(13)+10.);

        recTrack.fStatus = 0;
        recTrack.fSign = emt->Charge();
        Double_t z11 = (emt->GetZUncorrected() < -1.) ? emt->GetZUncorrected() : (Double_t)AliMUONConstants::DefaultChamberZ(10);
        recTrack.fV.Set(emt->GetNonBendingCoorUncorrected(),emt->GetBendingCoorUncorrected(),z11);
        recTrack.fP.Set(-TMath::Tan(emt->GetThetaXUncorrected()),-TMath::Tan(emt->GetThetaYUncorrected()),-1.);

        // produce eve track
        track = new AliEveTrack(&recTrack,trkProp);
        track->SetSourceObject(emt);
        minimalisticMuonTrack.SetTrackType(kMuonGhost);
    }
    track->SetAttLineAttMarker(trackList);
    trackList->AddElement(track);

    trackList->MakeTracks();
    gEve->AddElement(trackList);

    std::vector<TEveVector4D > muonPoints = trkProp->GetLastPoints();
    InsertPolyPoints(minimalisticMuonTrack, muonPoints);
}

void AliConverterPolylinesEngine::AddPolylinesToMinimalisticTrack(
        Int_t trackID,
        AliMinimalisticTrack &minimalisticTrack
) const
{
    Int_t maxTrackRadius = 520;
    TEveTrackList* tEveTrackList = new TEveTrackList("ESD-Tracks");
    TEveTrackPropagator *pTrackPropagator = tEveTrackList->GetPropagator();

    pTrackPropagator->SetMagField(-fESDEvent->GetMagneticField() / 10.0);
    pTrackPropagator->SetMaxR(maxTrackRadius);
    AliESDtrack *track = fESDEvent->GetTrack(trackID);
    AliEveTrack* evetrack = new AliEveTrack(track, pTrackPropagator);
    evetrack->SetSourceObject(track);
    // Add inner/outer track parameters as path-marks.
    if (track->IsOn(AliESDtrack::kTPCrefit))
    {
        Double_t pbuf[3], vbuf[3];
        if (track->GetInnerParam() != 0) {
            track->GetInnerParam()->GetXYZ(vbuf);
            track->GetInnerParam()->GetPxPyPz(pbuf);
        }
        if (track->GetOuterParam() != 0) {
            track->GetOuterParam()->GetXYZ(vbuf);
            track->GetOuterParam()->GetPxPyPz(pbuf);
        }
        TEvePathMark pm(TEvePathMark::kReference);
        pm.fV.Set(vbuf);
        pm.fP.Set(pbuf);
        evetrack->AddPathMark(pm);
    }
    tEveTrackList->SetChildClass(evetrack->Class());
    tEveTrackList->AddElement(evetrack);
    evetrack->MakeTrack();
    tEveTrackList->MakeTracks();

    std::vector<TEveVector4D> pointsVec = pTrackPropagator->GetLastPoints();
    InsertPolyPoints(minimalisticTrack, pointsVec);
}

void AliConverterPolylinesEngine::AddPolyLinesToKinkTrack(
        Int_t kinkID,
        AliMinimalisticTrack &mTrack,
        AliMinimalisticTrack &dTrack
) const
{
    AliEveKinkList* cont = new AliEveKinkList("ESD kink");

    TEveTrackPropagator* rnrStyleMoth = cont->GetPropagatorMoth();
    TEveTrackPropagator* rnrStyleDaugh = cont->GetPropagatorDaugh();

    rnrStyleMoth->SetMagField(0.1 * fESDEvent->GetMagneticField());
    rnrStyleDaugh->SetMagField(0.1 * fESDEvent->GetMagneticField());
    rnrStyleDaugh->SetMaxR(520);

    AliESDkink *kink = fESDEvent->GetKink(kinkID);
    AliESDtrack* moth = fESDEvent->GetTrack(kink->GetIndex(0));
    AliESDtrack* daug = fESDEvent->GetTrack(kink->GetIndex(1));

    TEveRecTrack rcMoth;
    TEveRecTrack rcDaug;
    TEveRecKink rcKink;

    rcKink.fPMother.Set(kink->GetMotherP());
    rcKink.fPDaughter.Set(kink->GetDaughterP());
    rcKink.fVKink.Set(kink->GetPosition());

    Double_t pbuf[3], vbuf[3];

    if (!moth->GetTPCInnerParam())
        return;
    rcMoth.fSign = moth->GetTPCInnerParam()->GetSign();
    moth->GetTPCInnerParam()->GetXYZ(vbuf); rcMoth.fV.Set(vbuf);
    moth->GetTPCInnerParam()->GetPxPyPz(pbuf); rcMoth.fP.Set(pbuf);
    if (!daug->GetOuterParam())
        return;
    rcDaug.fSign = daug->GetOuterParam()->GetSign();
    rcDaug.fV.Set(rcKink.fVKink);
    rcDaug.fP.Set(rcKink.fPDaughter);

    AliEveKink* myKink = new AliEveKink(&rcMoth, &rcDaug, &rcKink, rnrStyleMoth,rnrStyleDaugh);
    if (!myKink)
        return;
    myKink->SetESDKinkIndex(kinkID);
    gEve->AddElement(myKink, cont);
    cont->MakeKinks();

    std::vector<TEveVector4D> mPoints = rnrStyleMoth->GetLastPoints();
    std::vector<TEveVector4D> dPoints = rnrStyleDaugh->GetLastPoints();
    InsertPolyPoints(mTrack, mPoints);
    InsertPolyPoints(dTrack, dPoints);
    delete myKink;
    delete cont;
}

void AliConverterPolylinesEngine::AddPolyLinesToV0Track(
        Int_t v0ID,
        AliMinimalisticTrack &negativeTrack,
        AliMinimalisticTrack &positiveTrack) const
{
    AliEveV0List* cont = new AliEveV0List("ESD v0");
    TEveTrackPropagator* positivePropagator = cont->GetPropagatorPos();
    TEveTrackPropagator* negativePropagator = cont->GetPropagatorNeg();
    positivePropagator->SetMagField( 0.1*fESDEvent->GetMagneticField() );
    negativePropagator->SetMagField( 0.1*fESDEvent->GetMagneticField() );

    gEve->AddElement(cont);
    AliESDv0 *v0 = fESDEvent->GetV0(v0ID);

    TEveRecTrack rcPos;
    TEveRecTrack rcNeg;
    TEveRecV0 rcV0;

    Double_t pbuf[3], vbuf[3];

    rcNeg.fSign = v0->GetParamN()->GetSign();
    v0->GetParamN()->GetXYZ(vbuf); rcNeg.fV.Set(vbuf);
    v0->GetParamN()->GetPxPyPz(pbuf); rcNeg.fP.Set(pbuf);
    rcPos.fSign = v0->GetParamP()->GetSign();
    v0->GetParamP()->GetXYZ(vbuf); rcPos.fV.Set(vbuf);
    v0->GetParamP()->GetPxPyPz(pbuf); rcPos.fP.Set(pbuf);

    rcNeg.fIndex = v0->GetNindex();
    rcPos.fIndex = v0->GetPindex();

    AliEveV0* myV0 = new AliEveV0(&rcNeg, &rcPos, &rcV0, negativePropagator, positivePropagator);

    myV0->SetESDIndex(v0ID);
    myV0->SetOnFlyStatus(v0->GetOnFlyStatus());
    myV0->SetDaughterDCA(v0->GetDcaV0Daughters());

    gEve->AddElement(myV0, cont);

    cont->MakeV0s();

    std::vector<TEveVector4D> negPoints = negativePropagator->GetLastPoints();
    std::vector<TEveVector4D> posPoints = positivePropagator->GetLastPoints();
    InsertPolyPoints(negativeTrack, negPoints);
    InsertPolyPoints(positiveTrack, posPoints);
    delete myV0;
    delete cont;
}

void AliConverterPolylinesEngine::AddPolylinesToCascade(
        Int_t cascadeID,
        AliMinimalisticTrack &negativeTrack,
        AliMinimalisticTrack &positiveTrack,
        AliMinimalisticTrack &bachelorTrack
) const
{
    AliEveCascadeList* cont = new AliEveCascadeList("ESD cascade");
    TEveTrackPropagator* rnrStyleBac = cont->GetPropagatorBac();
    TEveTrackPropagator* rnrStyleNeg = cont->GetPropagatorNeg();
    TEveTrackPropagator* rnrStylePos = cont->GetPropagatorPos();
    rnrStyleBac->SetMagField(0.1 * fESDEvent->GetMagneticField());
    rnrStyleNeg->SetMagField(0.1 * fESDEvent->GetMagneticField());
    rnrStylePos->SetMagField(0.1 * fESDEvent->GetMagneticField());

    gEve->AddElement(cont);

    AliESDcascade *cascade = fESDEvent->GetCascade(cascadeID);

    TEveRecTrack rcPos,rcNeg,rcBac;
    TEveRecV0 rcV0;
    TEveRecCascade rcCascade;

    Double_t v[3],pBac[3]={0.}, pNeg[3]={0.}, pPos[3]={0.}, cv[21]={0.},pbuf[3], vbuf[3];

    cascade->GetBPxPyPz(pBac[0], pBac[1], pBac[2]);
    cascade->GetNPxPyPz(pNeg[0], pNeg[1], pNeg[2]);
    cascade->GetPPxPyPz(pPos[0], pPos[1], pPos[2]);
    cascade->GetXYZcascade(v[0], v[1], v[2]);

    AliExternalTrackParam *bacParam = new AliExternalTrackParam(v,pBac,cv,cascade->Charge());

    rcBac.fSign = bacParam->GetSign();
    rcBac.fIndex = cascade->GetBindex();
    bacParam->GetXYZ(vbuf); rcBac.fV.Set(vbuf);
    bacParam->GetPxPyPz(pbuf); rcBac.fP.Set(pbuf);

    rcNeg.fSign = cascade->GetParamN()->GetSign(); // Shouldn't it always be -1 ?
    rcNeg.fIndex = cascade->GetNindex();
    cascade->GetParamN()->GetXYZ(vbuf); rcNeg.fV.Set(vbuf);
    cascade->GetParamN()->GetPxPyPz(pbuf); rcNeg.fP.Set(pbuf);

    rcPos.fSign = cascade->GetParamP()->GetSign();  // Shouldn't it always be +1 ?
    rcPos.fIndex = cascade->GetPindex();
    cascade->GetParamP()->GetXYZ(vbuf); rcPos.fV.Set(vbuf);
    cascade->GetParamP()->GetPxPyPz(pbuf); rcPos.fP.Set(pbuf);

    AliEveCascade* myCascade = new AliEveCascade(
            &rcBac, &rcNeg, &rcPos, &rcV0, &rcCascade, rnrStyleBac, rnrStyleNeg, rnrStylePos
    );

    myCascade->SetESDIndex(cascadeID);
    myCascade->SetDaughterDCA(cascade->GetDcaXiDaughters());
    myCascade->SetLambdaP( pNeg[0]+pPos[0], pNeg[1]+pPos[1], pNeg[2]+pPos[2] );
    myCascade->SetBachP( pBac[0], pBac[1], pBac[2]);

    gEve->AddElement(myCascade, cont);
    cont->MakeCascades();

    std::vector<TEveVector4D> negPoints = rnrStyleNeg->GetLastPoints();
    std::vector<TEveVector4D> posPoints = rnrStylePos->GetLastPoints();
    std::vector<TEveVector4D> bacPoints = rnrStyleBac->GetLastPoints();
    InsertPolyPoints(negativeTrack, negPoints);
    InsertPolyPoints(positiveTrack, posPoints);
    InsertPolyPoints(bachelorTrack, bacPoints);
    delete myCascade;
    delete bacParam;
    delete cont;
}

void AliConverterPolylinesEngine::InsertPolyPoints(
        AliMinimalisticTrack &Track, std::vector<TEveVector4D> &Points
) const
{
    for(std::vector<TEveVector4D>::iterator iter = Points.begin(); iter != Points.end(); ++iter){
        Track.AddPolyPoint(*iter);
    }
}
