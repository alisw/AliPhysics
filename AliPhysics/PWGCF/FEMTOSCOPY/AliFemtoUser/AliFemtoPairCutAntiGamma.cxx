/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoPairCutAntiGamma - a pair cut which checks     //
// for some pair qualities that attempt to identify slit/doubly            //
// reconstructed tracks and also selects pairs based on their separation   //
// at the entrance to the TPC                                              //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoPairCutAntiGamma.cxx,v 1.1.2.1 2007/10/19 13:35:33 akisiel Exp $
 *
 * Author: Adam Kisiel, Ohio State, kisiel@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   a cut to remove "shared" and "split" pairs
 *
 ***************************************************************************
 *
 *
 **************************************************************************/

#include "AliFemtoPairCutAntiGamma.h"
#include <string>
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoPairCutAntiGamma)
#endif

//__________________
AliFemtoPairCutAntiGamma::AliFemtoPairCutAntiGamma():
  AliFemtoShareQualityPairCut(),
  fMaxEEMinv(0.0),
  fMaxDTheta(0.0),
  fDTPCMin(0),
  fMinAvgsep(0),
  fDataType(kESD),
  fNanoAODAnalysis(kFALSE)
{
}
//__________________
AliFemtoPairCutAntiGamma::AliFemtoPairCutAntiGamma(const AliFemtoPairCutAntiGamma& c):
  AliFemtoShareQualityPairCut(c),
  fMaxEEMinv(c.fMaxEEMinv),
  fMaxDTheta(c.fMaxDTheta),
  fDTPCMin(c.fDTPCMin),
  fMinAvgsep(c.fMinAvgsep),
  fDataType(c.fDataType),
  fNanoAODAnalysis(c.fNanoAODAnalysis)
{
}

AliFemtoPairCutAntiGamma&
AliFemtoPairCutAntiGamma::operator=(const AliFemtoPairCutAntiGamma& c)
{
    if (this != &c) {
        AliFemtoShareQualityPairCut::operator=(c);
        fMaxEEMinv = c.fMaxEEMinv;
        fMaxDTheta = c.fMaxDTheta;
        fDTPCMin = c.fDTPCMin;
        fMinAvgsep = c.fMinAvgsep;
        fDataType = c.fDataType;
        fNanoAODAnalysis = c.fNanoAODAnalysis;
    }

    return *this;

}
//__________________
AliFemtoPairCutAntiGamma::~AliFemtoPairCutAntiGamma(){
    /* no-op */
}
//__________________
bool AliFemtoPairCutAntiGamma::Pass(const AliFemtoPair* pair)
{
    // Accept pairs based on their TPC entrance separation and
    // quality and sharity
    bool temp = true;

    if(fDataType==kKine)
        return true;

    const AliFemtoTrack *track1 = pair->Track1()->Track(),
                        *track2 = pair->Track2()->Track();

    const AliFemtoThreeVector p1 = track1->P(),
                              p2 = track2->P();

    double me = 0.000511;

    if ((track1->Charge() * track2->Charge()) < 0.0) {

        double theta1 = p1.Theta();
        double theta2 = p2.Theta();
        double dtheta = TMath::Abs(theta1 - theta2);

        double e1 = TMath::Sqrt(me*me + p1.Mag2());
        double e2 = TMath::Sqrt(me*me + p2.Mag2());

        double minv = 2*me*me + 2*(e1*e2 - p1.Dot(p2));
        if ((TMath::Abs(minv) < fMaxEEMinv) && (dtheta < fMaxDTheta)) {
            temp = false;
        }
    }

    if(temp && fNanoAODAnalysis ) return true;
    // check separation at TPC entrance
    bool tempTPCEntrance = true;

    if(fDataType==kESD || fDataType==kAOD)
    {
        auto &epoint1 = track1->NominalTpcEntrancePoint(),
             &epoint2 = track2->NominalTpcEntrancePoint();

        tempTPCEntrance = (epoint1 - epoint2).Mag2() > fDTPCMin * fDTPCMin;
    }

    // check average separation
    bool avgsepCheck = true;

    if(fDataType==kESD || fDataType==kAOD)
    {
        // sums the separation magnitude
        double avgSep = 0.0;
        int count = 0;

        // loop through the 8 points of the 'NominalTpcPoint' methods
        for (int i = 0; i < 8; i++)
        {
            const AliFemtoThreeVector &point_1 = track1->NominalTpcPoint(i);
            const AliFemtoThreeVector &point_2 = track2->NominalTpcPoint(i);

            if (TpcPointIsUnset(point_1) || TpcPointIsUnset(point_2)) {
                break;
            }

            avgSep += (point_1 - point_2).Mag();
            count++;
        }
        avgSep /= count;
        avgsepCheck = avgSep > fMinAvgsep;
    }

    

    if (temp && tempTPCEntrance && avgsepCheck)
    {

        temp = AliFemtoShareQualityPairCut::Pass(pair);
        if (temp) {fNPairsPassed++;}
        else fNPairsFailed++;

        return temp;
    }
    else
    {
        fNPairsFailed++;
        return false;
    }
}

//__________________
AliFemtoString AliFemtoPairCutAntiGamma::Report()
{
    // Prepare a report from the execution
    AliFemtoString report = "AliFemtoPairCutAntiGamma Pair Cut - remove pairs possibly coming from Gamma conversions\n";
    report += Form("Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
    return report;
}

//__________________
TList *AliFemtoPairCutAntiGamma::ListSettings()
{
    // return a list of settings in a writable form
    TList *tListSetttings =  AliFemtoShareQualityPairCut::ListSettings();

    tListSetttings->AddVector(
        new TObjString(Form("AliFemtoPairCutAntiGamma.maxeeminv=%f", fMaxEEMinv)),
        new TObjString(Form("AliFemtoPairCutAntiGamma.maxdtheta=%f", fMaxDTheta)),
        new TObjString(Form("AliFemtoPairCutAntiGamma.dtpcmin=%f", fDTPCMin)),
        new TObjString(Form("AliFemtoPairCutAntiGamma.minavgsep=%f", fMinAvgsep)),
        nullptr
    );

    return tListSetttings;
}

void AliFemtoPairCutAntiGamma::SetMaxEEMinv(Double_t maxeeminv)
{
    fMaxEEMinv = maxeeminv;
}

void AliFemtoPairCutAntiGamma::SetMaxThetaDiff(Double_t maxdtheta)
{
    fMaxDTheta = maxdtheta;
}

void AliFemtoPairCutAntiGamma::SetTPCEntranceSepMinimum(double dtpc)
{
    fDTPCMin = dtpc;
}

void AliFemtoPairCutAntiGamma::SetAvgsepMinimum(double minAvgsep)
{
    fMinAvgsep = minAvgsep;
}

void AliFemtoPairCutAntiGamma::SetDataType(AliFemtoDataType type)
{
    fDataType = type;
}

bool AliFemtoPairCutAntiGamma::TpcPointIsUnset(const AliFemtoThreeVector& v)
{
    return v.x() < -9000. ||
    v.y() < -9000. ||
    v.z() < -9000.;
}
