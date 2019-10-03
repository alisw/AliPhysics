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
 * $Id: AliFemtoPairCutAntiGamma.h,v 1.1.2.1 2007/10/19 13:35:33 akisiel Exp $
 *
 * Author: Adam Kisiel, Ohio State University, kisiel@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   a cut to remove "shared" and "split" pairs
 *
 ***************************************************************************
 *
 *
 **************************************************************************/


#ifndef ALIFEMTOPAIRCUTANTIGAMMA_H
#define ALIFEMTOPAIRCUTANTIGAMMA_H

#include "AliFemtoPairCut.h"
#include "AliFemtoShareQualityPairCut.h"

class AliFemtoPairCutAntiGamma : public AliFemtoShareQualityPairCut{
public:

    AliFemtoPairCutAntiGamma();
    AliFemtoPairCutAntiGamma(const AliFemtoPairCutAntiGamma& c);
    virtual ~AliFemtoPairCutAntiGamma();
    AliFemtoPairCutAntiGamma& operator=(const AliFemtoPairCutAntiGamma& c);

    virtual bool Pass(const AliFemtoPair* pair);
    virtual AliFemtoString Report();
    virtual TList *ListSettings();
    virtual AliFemtoPairCutAntiGamma* Clone() const;
    void SetMaxEEMinv(Double_t maxeeminv);
    void SetMaxThetaDiff(Double_t maxdtheta);
    void SetTPCEntranceSepMinimum(double dtpc);
    void SetAvgsepMinimum(double minAvgsep);
    /* void SetTPCExitSepMinimum(double dtpc); */
    void SetDataType(AliFemtoDataType type);
    void SetNanoAODAnalysis(bool aNanoAOD);

    double GetMinAvgSep() const
      { return fMinAvgsep; }

    double GetMaxEEMinv() const
      { return fMaxEEMinv; }

    double GetMaxDTheta() const
      { return fMaxDTheta; }

    double GetDTPCMin() const
      { return fDTPCMin; }

    int GetDataType() const
      { return fDataType; }

protected:
    Double_t fMaxEEMinv; // Maximum allowed ee Minv
    Double_t fMaxDTheta; // Maximum polar angle difference
    Double_t fDTPCMin;          // Minimum allowed pair nominal separation at the entrance to the TPC
    Double_t fMinAvgsep;        // Minimum allowed pair avergae separation
    AliFemtoDataType fDataType; //Use ESD / AOD / Kinematics.
    bool fNanoAODAnalysis;


private:
    bool TpcPointIsUnset(const AliFemtoThreeVector& v);
#ifdef __ROOT__
    ClassDef(AliFemtoPairCutAntiGamma, 0)
#endif
};

inline AliFemtoPairCutAntiGamma* AliFemtoPairCutAntiGamma::Clone() const
{
  return new AliFemtoPairCutAntiGamma(*this);
}
inline void AliFemtoPairCutAntiGamma::SetNanoAODAnalysis(bool aNanoAOD) {fNanoAODAnalysis = aNanoAOD;}

#endif
