#include "AliVertexerHyperTriton2Body.h"

#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include <AliPIDResponse.h>
#include <AliTrackerBase.h>

#include <AliMCEvent.h>

#include <iostream>

ClassImp(AliVertexerHyperTriton2Body)

    AliVertexerHyperTriton2Body::AliVertexerHyperTriton2Body()
    : TNamed(), fHe3Cuts{nullptr}, fPiCuts{nullptr}, fLikeSign{false}, fRotation{false}, fLambda{false},
      //________________________________________________
      //Flags for V0 vertexer
      fMC{kFALSE},
      fkDoV0Refit(kFALSE),
      fMaxIterationsWhenMinimizing(27),
      fkPreselectX{true},
      fkXYCase1(kTRUE),
      fkXYCase2(kTRUE),
      fkResetInitialPositions(kFALSE),
      fkDoImprovedDCAV0DauPropagation(kFALSE),
      fkDoMaterialCorrection(kFALSE),
      fMinPtV0(2), //pre-selection
      fMaxPtV0(1000),
      fMinXforXYtest(-3.0),
      fV0VertexerSels{},
      fMassRange{2.9, 3.1},
      fMaxCt{45},
      fMagneticField{0.},
      fPrimaryVertexX{0.},
      fPrimaryVertexY{0.},
      fPrimaryVertexZ{0.},
      fPID{nullptr},
      fSpline{nullptr}

//________________________________________________
{
    SetupStandardVertexing();
}

//________________________________________________________________________
void AliVertexerHyperTriton2Body::SetupStandardVertexing()
//Meant to store standard re-vertexing configuration
{
    //Tell the task to re-run vertexers
    SetDoV0Refit(kTRUE);

    //V0-Related topological selections
    SetV0VertexerDCAFirstToPV(0.05);
    SetV0VertexerDCASecondToPV(0.05);
    SetV0VertexerDCAV0Daughters(1.20);
    SetV0VertexerCosinePA(0.98);
    SetV0VertexerMinRadius(0.9);
    SetV0VertexerMaxRadius(200);
}

//________________________________________________________________________
void AliVertexerHyperTriton2Body::SetupLooseVertexing()
//Meant to store standard re-vertexing configuration
{
    SetDoV0Refit(kTRUE);

    //V0-Related topological selections
    SetV0VertexerDCAFirstToPV(0.05);
    SetV0VertexerDCASecondToPV(0.05);
    SetV0VertexerDCAV0Daughters(1.60);
    SetV0VertexerCosinePA(0.99);
    SetV0VertexerMinRadius(0.9);
    SetV0VertexerMaxRadius(200);
}

//________________________________________________________________________
std::vector<AliESDv0> AliVertexerHyperTriton2Body::Tracks2V0vertices(AliESDEvent *event, AliPIDResponse *pid, AliMCEvent *mcEvent)
{
    //--------------------------------------------------------------------
    //This function reconstructs V0 vertices
    //--------------------------------------------------------------------
    std::vector<AliESDv0> v0s;

    fPID = pid;
    const AliESDVertex *vtxT3D = event->GetPrimaryVertex();

    fPrimaryVertexX = vtxT3D->GetX();
    fPrimaryVertexY = vtxT3D->GetY();
    fPrimaryVertexZ = vtxT3D->GetZ();

    Long_t nentr = event->GetNumberOfTracks();
    fMagneticField = event->GetMagneticField();

    if (nentr < 2)
        return v0s;

    std::vector<int> tracks[2][2];
    if (!mcEvent || fMC == kFALSE)
    {
        SelectTracks(event, tracks);
    }
    else
    {
        SelectTracksMC(event, mcEvent, tracks);
    }

    auto CreateV0 = [&](int nidx, AliESDtrack *ntrk, int pidx, AliESDtrack *ptrk) {
        Double_t lNegMassForTracking = ntrk->GetMassForTracking();
        Double_t lPosMassForTracking = ptrk->GetMassForTracking();

        float ndcaxy, ndcaz;
        ntrk->GetImpactParameters(ndcaxy, ndcaz);
        float pdcaxy, pdcaz;
        ptrk->GetImpactParameters(pdcaxy, pdcaz);
        if (std::abs(ndcaxy) < fV0VertexerSels[1] || std::abs(pdcaxy) < fV0VertexerSels[2])
            return;

        AliExternalTrackParam nt(*ntrk), pt(*ptrk);

        AliExternalTrackParam *ntp = &nt, *ptp = &pt;
        Double_t xn, xp, dca;

        //Improved call: use own function, including XY-pre-opt stage

        //Re-propagate to closest position to the primary vertex if asked to do so
        if (fkResetInitialPositions)
        {
            Double_t dztemp[2], covartemp[3];
            //Safety margin: 250 -> exceedingly large... not sure this makes sense, but ok
            ntp->PropagateToDCA(vtxT3D, fMagneticField, 250, dztemp, covartemp);
            ptp->PropagateToDCA(vtxT3D, fMagneticField, 250, dztemp, covartemp);
        }

        if (fkDoImprovedDCAV0DauPropagation)
        {
            //Improved: use own call
            dca = GetDCAV0Dau(ptp, ntp, xp, xn, fMagneticField, lNegMassForTracking, lPosMassForTracking);
        }
        else
        {
            //Old: use old call
            dca = nt.GetDCA(&pt, fMagneticField, xn, xp);
        }

        if (dca > fV0VertexerSels[3])
            return;

        if ((xn + xp) > 2 * fV0VertexerSels[6] && fkPreselectX)
            return;
        if ((xn + xp) < 2 * fV0VertexerSels[5] && fkPreselectX)
            return;

        if (!fkDoMaterialCorrection)
        {
            nt.PropagateTo(xn, fMagneticField);
            pt.PropagateTo(xp, fMagneticField);
        }
        else
        {
            AliTrackerBase::PropagateTrackTo(ntp, xn, lNegMassForTracking, 3, kFALSE, 0.75, kFALSE, kTRUE);
            AliTrackerBase::PropagateTrackTo(ptp, xp, lPosMassForTracking, 3, kFALSE, 0.75, kFALSE, kTRUE);
        }

        //select maximum eta range (after propagation)
        AliESDv0 vertex(nt, nidx, pt, pidx);

        //Experimental: refit V0 if asked to do so
        if (fkDoV0Refit)
            vertex.Refit();

        //No selection: it was not previously applied, don't  apply now.
        //if (vertex.GetChi2V0() > fChi2max) continue;

        Double_t x = vertex.Xv(), y = vertex.Yv();
        Double_t r2 = x * x + y * y;
        if (r2 < fV0VertexerSels[5] * fV0VertexerSels[5])
            return;
        if (r2 > fV0VertexerSels[6] * fV0VertexerSels[6])
            return;

        AliPID::EParticleType fatParticle = fLambda ? AliPID::kProton : AliPID::kHe3;
        Int_t posCharge = fLambda ? 1 : (std::abs(fPID->NumberOfSigmasTPC(ptrk, fatParticle)) < 5) + 1;
        Int_t negCharge = fLambda ? 1 : (std::abs(fPID->NumberOfSigmasTPC(ntrk, fatParticle)) < 5) + 1;
        Double_t posMass = posCharge > 1 ? AliPID::ParticleMass(fatParticle) : AliPID::ParticleMass(AliPID::kPion);
        Double_t negMass = negCharge > 1 ? AliPID::ParticleMass(fatParticle) : AliPID::ParticleMass(AliPID::kPion);

        Double_t posMom[3], negMom[3];
        LVector_t posVector, negVector, hyperVector;
        vertex.GetNPxPyPz(negMom[0], negMom[1], negMom[2]);
        vertex.GetPPxPyPz(posMom[0], posMom[1], posMom[2]);
        posVector.SetCoordinates(posCharge * posMom[0], posCharge * posMom[1], posCharge * posMom[2], posMass);
        negVector.SetCoordinates(negCharge * negMom[0], negCharge * negMom[1], negCharge * negMom[2], negMass);
        hyperVector = posVector + negVector;

        if (hyperVector.M() < fMassRange[0] || hyperVector.M() > fMassRange[1]) //selection on hypertriton invariant mass
            return;

        Double_t momV0[3] = {hyperVector.Px(), hyperVector.Py(), hyperVector.Pz()};
        Double_t deltaPos[3]; //vector between the reference point and the V0 vertex
        Double_t SPos[3];
        vertex.GetXYZ(SPos[0], SPos[1], SPos[2]);

        deltaPos[0] = SPos[0] - fPrimaryVertexX;
        deltaPos[1] = SPos[1] - fPrimaryVertexY;
        deltaPos[2] = SPos[2] - fPrimaryVertexZ;
        Double_t momV02 = momV0[0] * momV0[0] + momV0[1] * momV0[1] + momV0[2] * momV0[2];
        Double_t deltaPos2 = deltaPos[0] * deltaPos[0] + deltaPos[1] * deltaPos[1] + deltaPos[2] * deltaPos[2];
        Double_t ct = 2.99131 * TMath::Sqrt(deltaPos2 / momV02);
        if (ct > fMaxCt)
            return;

        double cpa = (deltaPos[0] * momV0[0] +
                      deltaPos[1] * momV0[1] +
                      deltaPos[2] * momV0[2]) /
                     TMath::Sqrt(momV02 * deltaPos2);
        //   Float_t cpa = vertex.GetV0CosineOfPointingAngle(fPrimaryVertexX, fPrimaryVertexY, fPrimaryVertexZ);

        //Simple cosine cut (no pt dependence for now)
        if (fSpline)
        {
            double cpaCut = fSpline->Eval(ct);
            if (cpa < cpaCut)
                return;
        }
        else
        {
            if (cpa < fV0VertexerSels[4])
                return;
        }

        vertex.SetDcaV0Daughters(dca);
        vertex.SetV0CosineOfPointingAngle(cpa);
        vertex.ChangeMassHypothesis(kK0Short);

        //pre-select on pT
        double lTransvMom = std::hypot(momV0[0], momV0[1]);
        if (lTransvMom < fMinPtV0)
            return;
        if (lTransvMom > fMaxPtV0)
            return;
        v0s.push_back(vertex);
    };

    if (!fLikeSign)
    {
        for (int index{0}; index < 2; ++index)
        {
            for (auto &nidx : tracks[1][index])
            {
                AliESDtrack *ntrk = event->GetTrack(nidx);
                if (!ntrk)
                    continue;
                if (fRotation && index == 1)
                {
                    double params[5]{ntrk->GetY(), ntrk->GetZ(), -ntrk->GetSnp(), ntrk->GetTgl(), ntrk->GetSigned1Pt()};
                    ntrk->SetParamOnly(ntrk->GetX(), ntrk->GetAlpha(), params);
                }
                for (auto &pidx : tracks[0][index == 1 ? 0 : 1])
                {
                    AliESDtrack *ptrk = event->GetTrack(pidx);
                    if (!ptrk)
                        continue;
                    if (fRotation && index == 0)
                    {
                        double params[5]{ptrk->GetY(), ptrk->GetZ(), -ptrk->GetSnp(), ptrk->GetTgl(), ptrk->GetSigned1Pt()};
                        ptrk->SetParamOnly(ptrk->GetX(), ptrk->GetAlpha(), params);
                    }

                    CreateV0(nidx, ntrk, pidx, ptrk);

                    if (fRotation && index == 0)
                    { /// Restore the params
                        double params[5]{ptrk->GetY(), ptrk->GetZ(), -ptrk->GetSnp(), ptrk->GetTgl(), ptrk->GetSigned1Pt()};
                        ptrk->SetParamOnly(ptrk->GetX(), ptrk->GetAlpha(), params);
                    }
                }
                if (fRotation && index == 1)
                {
                    double params[5]{ntrk->GetY(), ntrk->GetZ(), -ntrk->GetSnp(), ntrk->GetTgl(), ntrk->GetSigned1Pt()};
                    ntrk->SetParamOnly(ntrk->GetX(), ntrk->GetAlpha(), params);
                }
            }
        }
    }
    else
    {
        for (int index{0}; index < 2; ++index)
        {
            for (int charge{0}; charge < 2; ++charge)
            {
                for (auto &nidx : tracks[charge][index])
                {
                    AliESDtrack *ntrk = event->GetTrack(nidx);
                    if (!ntrk)
                        continue;
                    for (auto &pidx : tracks[charge][index > 0 ? 0 : 1])
                    {
                        AliESDtrack *ptrk = event->GetTrack(pidx);
                        if (!ptrk)
                            continue;
                        CreateV0(nidx, ntrk, pidx, ptrk);
                    }
                }
            }
        }
    }
    return v0s;
}

//________________________________________________________________________
Double_t AliVertexerHyperTriton2Body::Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const
{
    //--------------------------------------------------------------------
    // This function calculates locally a 2x2 determinant
    //--------------------------------------------------------------------
    return a00 * a11 - a01 * a10;
}

//________________________________________________________________________
Double_t AliVertexerHyperTriton2Body::Det(Double_t a00, Double_t a01, Double_t a02,
                                          Double_t a10, Double_t a11, Double_t a12,
                                          Double_t a20, Double_t a21, Double_t a22) const
{
    //--------------------------------------------------------------------
    // This function calculates locally a 3x3 determinant
    //--------------------------------------------------------------------
    return a00 * Det(a11, a12, a21, a22) - a01 * Det(a10, a12, a20, a22) + a02 * Det(a10, a11, a20, a21);
}

//________________________________________________________________________
void AliVertexerHyperTriton2Body::Evaluate(const Double_t *h, Double_t t,
                                           Double_t r[3],  //radius vector
                                           Double_t g[3],  //first defivatives
                                           Double_t gg[3]) //second derivatives
{
    //--------------------------------------------------------------------
    // Calculate position of a point on a track and some derivatives
    //--------------------------------------------------------------------
    Double_t phase = h[4] * t + h[2];
    Double_t sn = TMath::Sin(phase), cs = TMath::Cos(phase);

    r[0] = h[5];
    r[1] = h[0];
    if (TMath::Abs(h[4]) > kAlmost0)
    {
        r[0] += (sn - h[6]) / h[4];
        r[1] -= (cs - h[7]) / h[4];
    }
    else
    {
        r[0] += t * cs;
        r[1] -= -t * sn;
    }
    r[2] = h[1] + h[3] * t;

    g[0] = cs;
    g[1] = sn;
    g[2] = h[3];

    gg[0] = -h[4] * sn;
    gg[1] = h[4] * cs;
    gg[2] = 0.;
}

//________________________________________________________________________
void AliVertexerHyperTriton2Body::CheckChargeV0(AliESDv0 *v0)
{
    // This function checks charge of negative and positive daughter tracks.
    // If incorrectly defined (onfly vertexer), swaps out.
    if (v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0)
    {
        //V0 daughter track swapping is required! Note: everything is swapped here... P->N, N->P
        Long_t lCorrectNidx = v0->GetPindex();
        Long_t lCorrectPidx = v0->GetNindex();
        Double32_t lCorrectNmom[3];
        Double32_t lCorrectPmom[3];
        v0->GetPPxPyPz(lCorrectNmom[0], lCorrectNmom[1], lCorrectNmom[2]);
        v0->GetNPxPyPz(lCorrectPmom[0], lCorrectPmom[1], lCorrectPmom[2]);

        AliExternalTrackParam lCorrectParamN(
            v0->GetParamP()->GetX(),
            v0->GetParamP()->GetAlpha(),
            v0->GetParamP()->GetParameter(),
            v0->GetParamP()->GetCovariance());
        AliExternalTrackParam lCorrectParamP(
            v0->GetParamN()->GetX(),
            v0->GetParamN()->GetAlpha(),
            v0->GetParamN()->GetParameter(),
            v0->GetParamN()->GetCovariance());
        lCorrectParamN.SetMostProbablePt(v0->GetParamP()->GetMostProbablePt());
        lCorrectParamP.SetMostProbablePt(v0->GetParamN()->GetMostProbablePt());

        //Get Variables___________________________________________________
        Double_t lDcaV0Daughters = v0->GetDcaV0Daughters();
        Double_t lCosPALocal = v0->GetV0CosineOfPointingAngle();
        Bool_t lOnFlyStatusLocal = v0->GetOnFlyStatus();

        //Create Replacement Object_______________________________________
        AliESDv0 *v0correct = new AliESDv0(lCorrectParamN, lCorrectNidx, lCorrectParamP, lCorrectPidx);
        v0correct->SetDcaV0Daughters(lDcaV0Daughters);
        v0correct->SetV0CosineOfPointingAngle(lCosPALocal);
        v0correct->ChangeMassHypothesis(kK0Short);
        v0correct->SetOnFlyStatus(lOnFlyStatusLocal);

        //Reverse Cluster info..._________________________________________
        v0correct->SetClusters(v0->GetClusters(1), v0->GetClusters(0));

        *v0 = *v0correct;
        //Proper cleanup..._______________________________________________
        v0correct->Delete();
        v0correct = 0x0;

        //Just another cross-check and output_____________________________
        if (v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0)
        {
            AliWarning("Found Swapped Charges, tried to correct but something FAILED!");
        }
        else
        {
            //AliWarning("Found Swapped Charges and fixed.");
        }
        //________________________________________________________________
    }
    else
    {
        //Don't touch it! ---
        //Printf("Ah, nice. Charges are already ordered...");
    }
    return;
}

Double_t AliVertexerHyperTriton2Body::GetDCAV0Dau(AliExternalTrackParam *pt, AliExternalTrackParam *nt, Double_t &xp, Double_t &xn, double b, Double_t lNegMassForTracking, Double_t lPosMassForTracking)
{
    //--------------------------------------------------------------
    // Propagates this track and the argument track to the position of the
    // distance of closest approach.
    // Returns the (weighed !) distance of closest approach.
    //--------------------------------------------------------------

    //if( fkDoPureGeometricMinimization ){
    //Override uncertainties with small values -> pure geometry
    //dx2 = 1e-10;
    //dy2 = 1e-10;
    //dz2 = 1e-10;
    //}

    Double_t p1[8];
    nt->GetHelixParameters(p1, b);
    p1[6] = TMath::Sin(p1[2]);
    p1[7] = TMath::Cos(p1[2]);
    Double_t p2[8];
    pt->GetHelixParameters(p2, b);
    p2[6] = TMath::Sin(p2[2]);
    p2[7] = TMath::Cos(p2[2]);

    //Minimum X: allow for negative X if it means we're still *after* the primary vertex in the track ref frame
    Double_t lMinimumX = fMinXforXYtest; //
    //Maximum X: some very big value, should not be a problem
    Double_t lMaximumX = 300;

    if (fkDoImprovedDCAV0DauPropagation)
    {
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        // V0 preprocessing: analytical estimate of DCAxy position
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        Double_t nhelix[6], phelix[6];
        nt->GetHelixParameters(nhelix, b);
        pt->GetHelixParameters(phelix, b);
        Double_t lNegCenterR[2], lPosCenterR[2];

        //Negative track parameters in XY
        GetHelixCenter(nt, lNegCenterR, b);
        Double_t xNegCenter = lNegCenterR[0];
        Double_t yNegCenter = lNegCenterR[1];
        Double_t NegRadius = TMath::Abs(1. / nhelix[4]);

        //Positive track parameters in XY
        GetHelixCenter(pt, lPosCenterR, b);
        Double_t xPosCenter = lPosCenterR[0];
        Double_t yPosCenter = lPosCenterR[1];
        Double_t PosRadius = TMath::Abs(1. / phelix[4]);

        //Define convenient coordinate system
        //Logical zero: position of negative center
        Double_t ux = xPosCenter - xNegCenter;
        Double_t uy = yPosCenter - yNegCenter;

        //Check center-to-center distance
        Double_t lDist = TMath::Sqrt(
            TMath::Power(xNegCenter - xPosCenter, 2) +
            TMath::Power(yNegCenter - yPosCenter, 2));
        //Normalize ux, uz to unit vector
        ux /= lDist;
        uy /= lDist;

        //Calculate perpendicular vector (normalized)
        Double_t vx = -uy;
        Double_t vy = +ux;

        Double_t lPreprocessDCAxy = 1e+3;    //define outside scope
        Double_t lPreprocessxp = pt->GetX(); //start at current location
        Double_t lPreprocessxn = nt->GetX(); //start at current location

        //============================================================
        //Pre-optimization in the XY plane: cases considered here
        //============================================================
        //
        //  Case 1: Circles do not touch, centers far away
        //          (D > R1 + R2)
        //
        //  Case 2: Circles touch, centers at reasonable distance wrt D
        //          (D < R1 + R2) && (D > |R1-R2|)
        //
        //  Case 3: Circles do not touch, one inside the other
        //          (D < |R1-R2|)
        //
        //  Cases 1 and 2 are treated. Case 3 is not treated (unlikely
        //  to be a problem with unlike-sign charged tracks): brute
        //  force minimization takes place in any case
        //
        //============================================================

        //______________________
        //fast skipper: if XY plane pre-optimization says they're far, they're far! don't insist
        //TODO: check if this is relevant for HyperTriton
        // if ( fkSkipLargeXYDCA ) {
        //     if( lDist > NegRadius + PosRadius + 2*fV0VertexerSels[3] ) return 2000;
        //     if( lDist < TMath::Abs(NegRadius - PosRadius) - 2*fV0VertexerSels[3] ) return 2000;
        // }

        //______________________
        //CASE 1
        if ((lDist > NegRadius + PosRadius) && fkXYCase1)
        {
            //================================================================
            //Case 1: distance bigger than sum of radii ("gamma-like")
            //        re-position tracks along the center-to-center axis
            //Re-position negative track
            Double_t xNegOptPosition = xNegCenter + NegRadius * ux;
            Double_t yNegOptPosition = yNegCenter + NegRadius * uy;
            Double_t csNeg = TMath::Cos(nt->GetAlpha());
            Double_t snNeg = TMath::Sin(nt->GetAlpha());
            Double_t xThisNeg = xNegOptPosition * csNeg + yNegOptPosition * snNeg;

            //Re-position positive track
            Double_t xPosOptPosition = xPosCenter - PosRadius * ux;
            Double_t yPosOptPosition = yPosCenter - PosRadius * uy;
            Double_t csPos = TMath::Cos(pt->GetAlpha());
            Double_t snPos = TMath::Sin(pt->GetAlpha());
            Double_t xThisPos = xPosOptPosition * csPos + yPosOptPosition * snPos;

            if (xThisNeg < lMaximumX && xThisPos < lMaximumX && xThisNeg > lMinimumX && xThisPos > lMinimumX)
            {
                Bool_t lPropagA = kFALSE, lPropagB = kFALSE;
                Double_t lCase1NegR[3];
                lPropagA = nt->GetXYZAt(xThisNeg, b, lCase1NegR);
                Double_t lCase1PosR[3];
                lPropagB = pt->GetXYZAt(xThisPos, b, lCase1PosR);
                if (lPropagA && lPropagB)
                {
                    lPreprocessDCAxy = TMath::Sqrt(
                        TMath::Power(lCase1NegR[0] - lCase1PosR[0], 2) +
                        TMath::Power(lCase1NegR[1] - lCase1PosR[1], 2) +
                        TMath::Power(lCase1NegR[2] - lCase1PosR[2], 2));
                    //Pass coordinates
                    if (lPreprocessDCAxy < 999)
                    {
                        lPreprocessxp = xThisPos;
                        lPreprocessxn = xThisNeg;
                    }
                }
            }
            //================================================================
        }

        //______________________
        //CASE 2
        if ((lDist > TMath::Abs(NegRadius - PosRadius)) && (lDist < NegRadius + PosRadius) && fkXYCase2)
        {
            //================================================================
            //Case 2: distance smaller than sum of radii (cowboy/sailor configs)

            //Calculate coordinate for radical line
            Double_t lRadical = (lDist * lDist - PosRadius * PosRadius + NegRadius * NegRadius) / (2 * lDist);

            //Calculate absolute displacement from center-to-center axis
            Double_t lDisplace = (0.5 / lDist) * TMath::Sqrt(
                                                     (-lDist + PosRadius - NegRadius) *
                                                     (-lDist - PosRadius + NegRadius) *
                                                     (-lDist + PosRadius + NegRadius) *
                                                     (lDist + PosRadius + NegRadius));

            Double_t lCase2aDCA = 1e+3;
            Double_t lCase2bDCA = 1e+3;

            //2 cases: positive and negative displacement
            Double_t xNegOptPosition[2], yNegOptPosition[2], xPosOptPosition[2], yPosOptPosition[2];
            Double_t csNeg, snNeg, csPos, snPos;
            Double_t xThisNeg[2], xThisPos[2];

            csNeg = TMath::Cos(nt->GetAlpha());
            snNeg = TMath::Sin(nt->GetAlpha());
            csPos = TMath::Cos(pt->GetAlpha());
            snPos = TMath::Sin(pt->GetAlpha());

            //Case 2a: Positive displacement along v vector
            //Re-position negative track
            xNegOptPosition[0] = xNegCenter + lRadical * ux + lDisplace * vx;
            yNegOptPosition[0] = yNegCenter + lRadical * uy + lDisplace * vy;
            xThisNeg[0] = xNegOptPosition[0] * csNeg + yNegOptPosition[0] * snNeg;
            //Re-position positive track
            xPosOptPosition[0] = xNegCenter + lRadical * ux + lDisplace * vx;
            yPosOptPosition[0] = yNegCenter + lRadical * uy + lDisplace * vy;
            xThisPos[0] = xPosOptPosition[0] * csPos + yPosOptPosition[0] * snPos;

            //Case 2b: Negative displacement along v vector
            //Re-position negative track
            xNegOptPosition[1] = xNegCenter + lRadical * ux - lDisplace * vx;
            yNegOptPosition[1] = yNegCenter + lRadical * uy - lDisplace * vy;
            xThisNeg[1] = xNegOptPosition[1] * csNeg + yNegOptPosition[1] * snNeg;
            //Re-position positive track
            xPosOptPosition[1] = xNegCenter + lRadical * ux - lDisplace * vx;
            yPosOptPosition[1] = yNegCenter + lRadical * uy - lDisplace * vy;
            xThisPos[1] = xPosOptPosition[1] * csPos + yPosOptPosition[1] * snPos;

            //Test the two cases, please

            //Case 2a
            if (xThisNeg[0] < lMaximumX && xThisPos[0] < lMaximumX && xThisNeg[0] > lMinimumX && xThisPos[0] > lMinimumX)
            {
                Bool_t lPropagA = kFALSE, lPropagB = kFALSE;
                Double_t lCase2aNegR[3];
                lPropagA = nt->GetXYZAt(xThisNeg[0], b, lCase2aNegR);
                Double_t lCase2aPosR[3];
                lPropagB = pt->GetXYZAt(xThisPos[0], b, lCase2aPosR);
                if (lPropagA && lPropagB)
                {
                    lCase2aDCA = TMath::Sqrt(
                        TMath::Power(lCase2aNegR[0] - lCase2aPosR[0], 2) +
                        TMath::Power(lCase2aNegR[1] - lCase2aPosR[1], 2) +
                        TMath::Power(lCase2aNegR[2] - lCase2aPosR[2], 2));
                }
                else
                {
                    for (Int_t ic = 0; ic < 3; ic++)
                        lCase2aNegR[ic] = 0;
                    for (Int_t ic = 0; ic < 3; ic++)
                        lCase2aPosR[ic] = 0;
                    lCase2aDCA = 1e+4;
                }
            }

            //Case 2b
            if (xThisNeg[1] < lMaximumX && xThisPos[1] < lMaximumX && xThisNeg[1] > lMinimumX && xThisPos[1] > lMinimumX)
            {
                Bool_t lPropagA = kFALSE, lPropagB = kFALSE;
                Double_t lCase2bNegR[3];
                lPropagA = nt->GetXYZAt(xThisNeg[1], b, lCase2bNegR);
                Double_t lCase2bPosR[3];
                lPropagB = pt->GetXYZAt(xThisPos[1], b, lCase2bPosR);
                if (lPropagA && lPropagB)
                {
                    lCase2bDCA = TMath::Sqrt(
                        TMath::Power(lCase2bNegR[0] - lCase2bPosR[0], 2) +
                        TMath::Power(lCase2bNegR[1] - lCase2bPosR[1], 2) +
                        TMath::Power(lCase2bNegR[2] - lCase2bPosR[2], 2));
                }
                else
                {
                    for (Int_t ic = 0; ic < 3; ic++)
                        lCase2bNegR[ic] = 0;
                    for (Int_t ic = 0; ic < 3; ic++)
                        lCase2bPosR[ic] = 0;
                    lCase2bDCA = 1e+4;
                }
            }

            //Minor detail: all things being equal, prefer closest X
            Double_t lCase2aSumX = xThisPos[0] + xThisNeg[0];
            Double_t lCase2bSumX = xThisPos[1] + xThisNeg[1];

            Double_t lDCAxySmallestR = lCase2aDCA;
            Double_t lxpSmallestR = xThisPos[0];
            Double_t lxnSmallestR = xThisNeg[0];

            Double_t lDCAxyLargestR = lCase2bDCA;
            Double_t lxpLargestR = xThisPos[1];
            Double_t lxnLargestR = xThisNeg[1];

            if (lCase2bSumX + 1e-6 < lCase2aSumX)
            {
                lDCAxySmallestR = lCase2bDCA;
                lxpSmallestR = xThisPos[1];
                lxnSmallestR = xThisNeg[1];
                lDCAxyLargestR = lCase2aDCA;
                lxpLargestR = xThisPos[0];
                lxnLargestR = xThisNeg[0];
            }

            //Pass conclusion to lPreprocess variables, please
            lPreprocessDCAxy = lDCAxySmallestR;
            lPreprocessxp = lxpSmallestR;
            lPreprocessxn = lxnSmallestR;
            if (lDCAxyLargestR + 1e-6 < lDCAxySmallestR)
            { //beware epsilon: numerical calculations are unstable here
                lPreprocessDCAxy = lDCAxyLargestR;
                lPreprocessxp = lxpLargestR;
                lPreprocessxn = lxnLargestR;
            }
            //Protection against something too crazy, please
            if (lPreprocessDCAxy > 999)
            {
                lPreprocessxp = pt->GetX(); //start at current location
                lPreprocessxn = nt->GetX(); //start at current location
            }
        }
        //End of preprocessing stage!
        //at this point lPreprocessxp, lPreprocessxn are already good starting points: update helixparams
        if (lPreprocessDCAxy < 999)
        { //some improvement... otherwise discard in all cases, please
            if (fkDoMaterialCorrection)
            {
                AliTrackerBase::PropagateTrackTo(nt, lPreprocessxn, lNegMassForTracking, 3, kFALSE, 0.75, kFALSE, kTRUE);
                AliTrackerBase::PropagateTrackTo(pt, lPreprocessxp, lPosMassForTracking, 3, kFALSE, 0.75, kFALSE, kTRUE);
            }
            else
            {
                nt->PropagateTo(lPreprocessxn, b);
                pt->PropagateTo(lPreprocessxp, b);
            }
        }

        //don't redefine!
        nt->GetHelixParameters(p1, b);
        p1[6] = TMath::Sin(p1[2]);
        p1[7] = TMath::Cos(p1[2]);
        pt->GetHelixParameters(p2, b);
        p2[6] = TMath::Sin(p2[2]);
        p2[7] = TMath::Cos(p2[2]);
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    }

    Double_t dy2 = nt->GetSigmaY2() + pt->GetSigmaY2();
    Double_t dz2 = nt->GetSigmaZ2() + pt->GetSigmaZ2();
    Double_t dx2 = dy2;

    Double_t r1[3], g1[3], gg1[3];
    Double_t t1 = 0.;
    Evaluate(p1, t1, r1, g1, gg1);
    Double_t r2[3], g2[3], gg2[3];
    Double_t t2 = 0.;
    Evaluate(p2, t2, r2, g2, gg2);

    Double_t dx = r2[0] - r1[0], dy = r2[1] - r1[1], dz = r2[2] - r1[2];
    Double_t dm = dx * dx / dx2 + dy * dy / dy2 + dz * dz / dz2;

    Int_t max = fMaxIterationsWhenMinimizing;
    while (max--)
    {
        Double_t gt1 = -(dx * g1[0] / dx2 + dy * g1[1] / dy2 + dz * g1[2] / dz2);
        Double_t gt2 = +(dx * g2[0] / dx2 + dy * g2[1] / dy2 + dz * g2[2] / dz2);
        Double_t h11 = (g1[0] * g1[0] - dx * gg1[0]) / dx2 +
                       (g1[1] * g1[1] - dy * gg1[1]) / dy2 +
                       (g1[2] * g1[2] - dz * gg1[2]) / dz2;
        Double_t h22 = (g2[0] * g2[0] + dx * gg2[0]) / dx2 +
                       (g2[1] * g2[1] + dy * gg2[1]) / dy2 +
                       (g2[2] * g2[2] + dz * gg2[2]) / dz2;
        Double_t h12 = -(g1[0] * g2[0] / dx2 + g1[1] * g2[1] / dy2 + g1[2] * g2[2] / dz2);

        Double_t det = h11 * h22 - h12 * h12;

        Double_t dt1, dt2;
        if (TMath::Abs(det) < 1.e-33)
        {
            //(quasi)singular Hessian
            dt1 = -gt1;
            dt2 = -gt2;
        }
        else
        {
            dt1 = -(gt1 * h22 - gt2 * h12) / det;
            dt2 = -(h11 * gt2 - h12 * gt1) / det;
        }

        if ((dt1 * gt1 + dt2 * gt2) > 0)
        {
            dt1 = -dt1;
            dt2 = -dt2;
        }

        //check delta(phase1) ?
        //check delta(phase2) ?

        if (TMath::Abs(dt1) / (TMath::Abs(t1) + 1.e-3) < 1.e-4)
            if (TMath::Abs(dt2) / (TMath::Abs(t2) + 1.e-3) < 1.e-4)
            {
                if ((gt1 * gt1 + gt2 * gt2) > 1.e-4 / dy2 / dy2)
                    AliDebug(1, " stopped at not a stationary point !");
                Double_t lmb = h11 + h22;
                lmb = lmb - TMath::Sqrt(lmb * lmb - 4 * det);
                if (lmb < 0.)
                    AliDebug(1, " stopped at not a minimum !");
                break;
            }

        Double_t dd = dm;
        for (Int_t div = 1;; div *= 2)
        {
            Evaluate(p1, t1 + dt1, r1, g1, gg1);
            Evaluate(p2, t2 + dt2, r2, g2, gg2);
            dx = r2[0] - r1[0];
            dy = r2[1] - r1[1];
            dz = r2[2] - r1[2];
            dd = dx * dx / dx2 + dy * dy / dy2 + dz * dz / dz2;
            if (dd < dm)
                break;
            dt1 *= 0.5;
            dt2 *= 0.5;
            if (div > 512)
            {
                AliDebug(1, " overshoot !");
                break;
            }
        }
        dm = dd;

        t1 += dt1;
        t2 += dt2;
    }

    if (max <= 0)
        AliDebug(1, " too many iterations !");

    Double_t cs = TMath::Cos(nt->GetAlpha());
    Double_t sn = TMath::Sin(nt->GetAlpha());
    xn = r1[0] * cs + r1[1] * sn;

    cs = TMath::Cos(pt->GetAlpha());
    sn = TMath::Sin(pt->GetAlpha());
    xp = r2[0] * cs + r2[1] * sn;

    return TMath::Sqrt(dm * TMath::Sqrt(dy2 * dz2));
}

///________________________________________________________________________
void AliVertexerHyperTriton2Body::GetHelixCenter(const AliExternalTrackParam *track, Double_t center[2], double b)
{
    // Copied from AliV0ReaderV1::GetHelixCenter
    // Get Center of the helix track parametrization

    Int_t charge = track->Charge();

    Double_t helix[6];
    track->GetHelixParameters(helix, b);

    Double_t xpos = helix[5];
    Double_t ypos = helix[0];
    Double_t radius = TMath::Abs(1. / helix[4]);
    Double_t phi = helix[2];
    if (phi < 0)
    {
        phi = phi + 2 * TMath::Pi();
    }
    phi -= TMath::Pi() / 2.;
    Double_t xpoint = radius * TMath::Cos(phi);
    Double_t ypoint = radius * TMath::Sin(phi);
    if (b < 0 && charge > 0)
    {
        xpoint = -xpoint;
        ypoint = -ypoint;
    }
    if (b > 0 && charge < 0)
    {
        xpoint = -xpoint;
        ypoint = -ypoint;
    }
    center[0] = xpos + xpoint;
    center[1] = ypos + ypoint;
    return;
}

void AliVertexerHyperTriton2Body::SelectTracks(AliESDEvent *event, std::vector<int> tracks[2][2])
{

    if (!fPiCuts)
        fPiCuts = SetPionTPCTrackCuts();
    if (!fHe3Cuts)
        fHe3Cuts = SetHe3TPCTrackCuts();

    fMagneticField = event->GetMagneticField();

    AliPID::EParticleType fatParticle = fLambda ? AliPID::kProton : AliPID::kHe3;
    for (int i = 0; i < event->GetNumberOfTracks(); i++)
    {
        AliESDtrack *esdTrack = event->GetTrack(i);

        float d, z;
        esdTrack->GetImpactParameters(d, z);
        if (TMath::Abs(d) < fV0VertexerSels[2])
            continue;
        if (TMath::Abs(d) > fV0VertexerSels[6])
            continue;

        const int index = int(esdTrack->GetSign() < 0.);
        if (std::abs(fPID->NumberOfSigmasTPC(esdTrack, fatParticle)) < 5 &&
            fHe3Cuts->AcceptTrack(esdTrack))
            tracks[index][0].push_back(i);
        else if (fPiCuts->AcceptTrack(esdTrack))
            tracks[index][1].push_back(i);
    }
}

void AliVertexerHyperTriton2Body::SelectTracksMC(AliESDEvent *event, AliMCEvent *mcEvent, std::vector<int> tracks[2][2])
{

    for (int i = 0; i < event->GetNumberOfTracks(); i++)
    {
        AliESDtrack *esdTrack = event->GetTrack(i);

        int label = std::abs(esdTrack->GetLabel());

        if (!fPiCuts)
            fPiCuts = SetPionTPCTrackCuts();
        if (!fHe3Cuts)
            fHe3Cuts = SetHe3TPCTrackCuts();

        fMagneticField = event->GetMagneticField();

        float d, z;
        esdTrack->GetImpactParameters(d, z);
        if (TMath::Abs(d) < fV0VertexerSels[2])
            continue;
        if (TMath::Abs(d) > fV0VertexerSels[6])
            continue;

        AliVParticle *part = mcEvent->GetTrack(label);
        AliMCParticle *mother = mcEvent->MotherOfParticle(label);
        int fatParticle = fLambda ? 3122 : 1010010030;
        if (!mother || !part)
            continue;
        if (std::abs(mother->PdgCode()) != fatParticle)
            continue;
        const int index = int(esdTrack->GetSign() < 0.);
        const int pion = std::abs(part->PdgCode()) == 211;
        if (pion == 1)
        {
            if (fPiCuts->AcceptTrack(esdTrack))
                tracks[index][pion].push_back(i);
        }

        else
        {
            if (fHe3Cuts->AcceptTrack(esdTrack))
                tracks[index][pion].push_back(i);
        }
    }
}
