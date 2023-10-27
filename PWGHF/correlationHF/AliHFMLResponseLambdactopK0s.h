#ifndef AliHFMLResponseLambdactopK0s_H
#define AliHFMLResponseLambdactopK0s_H

// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliHFMLResponseLambdactopK0s
// \brief helper class to handle application of ML models for Lc->pK0s analysis trained
// with python libraries
// \authors:
// A. Palasciano, antonio.palasciano@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include "AliHFMLResponse.h"

class AliHFMLResponseLambdactopK0s : public AliHFMLResponse
{
public:
    AliHFMLResponseLambdactopK0s();
    AliHFMLResponseLambdactopK0s(const Char_t *name, const Char_t *title, const std::string configfilepath);
    virtual ~AliHFMLResponseLambdactopK0s();

    AliHFMLResponseLambdactopK0s(const AliHFMLResponseLambdactopK0s &source);
    AliHFMLResponseLambdactopK0s& operator=(const AliHFMLResponseLambdactopK0s& source);
private:
        int fNProngs = 3;
        int fCandType = 0;                                               /// flag for candidate type (bit map above)
        float fInvMass = -999.;                                          /// candidate invariant mass
        float fPt = -999.;                                               /// candidate pt
        float fDecayLength = -999.;                                      /// candidate decay length
        float fDecayLengthXY = -999.;                                    /// candidate decay length in the transverse plane
        float fNormDecayLengthXY = -999.;                                /// candidate normalised decay length in the transverse plane
        float fCosP = -999.;                                             /// candidate cosine of pointing angle
        float fCosPXY = -999.;                                           /// candidate cosine of pointing angle in the transcverse plane
        float fImpParXY = -999.;                                         /// candidate impact parameter in the transverse plane
        float fDCA = -999.;        
        float fImpParProng[4] = {-999., -999., -999., -999.};            /// prong impact parameter
        float fNormd0MeasMinusExp = -999.;                               /// candidate topomatic variable
        float fMassK0s = -999.;                                          /// mass of the pipi in case of pK0s hypothesis
        float fMassL = -999.;                                            /// mass of the ppi in case of piL hypothesis
        float fImpParV0 = -999.;                                         /// impact parameter of V0
        float fDecayLengthV0 = -999.;                                    /// decay length of V0
        float fDCAV0 = -999.;                                            /// DCA V0 prongs
        float fPtV0 = -999.;                                             /// pT of V0
        float fcTauK0s = -999.;                                          /// cTau of K0s
        float fcTauL = -999.;                                            /// cTau of L
        float fCosPV0 = -999.;                                           /// cosine of pointing angle of V0
        float fsignd0 = -999.;                                           /// signed d0 proton (different from standard d0)
        float fArmqTOverAlphaV0 = -999.;                                 /// Armenteros qT/|alpha| of the V0
        bool fUseKFParticle = false;                                     /// flag to enable/disable recalculation of Lc decay vertex with KFParticle package
        float fKFTopoChi2 = -999.;                                       /// KF chi2 to primary vertex with topological constraint

protected:
    virtual void SetMapOfVariables(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int /*masshypo*/);

    /// \cond CLASSIMP
    ClassDef(AliHFMLResponseLambdactopK0s, 1); ///
    /// \endcond
};
#endif
