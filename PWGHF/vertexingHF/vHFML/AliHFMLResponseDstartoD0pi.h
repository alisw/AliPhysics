#ifndef ALIHFMLRESPONSEDSTARTOD0PI_H
#define ALIHFMLRESPONSEDSTARTOD0PI_H

// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliHFMLResponseDstartoD0pi
// \brief helper class to handle application of ML models for D+ analyses trained
// with python libraries
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include "AliHFMLResponse.h"

class AliHFMLResponseDstartoD0pi : public AliHFMLResponse
{
public:
    AliHFMLResponseDstartoD0pi();
    AliHFMLResponseDstartoD0pi(const Char_t *name, const Char_t *title, const std::string configfilepath);
    virtual ~AliHFMLResponseDstartoD0pi();

    AliHFMLResponseDstartoD0pi(const AliHFMLResponseDstartoD0pi &source);
    AliHFMLResponseDstartoD0pi& operator=(const AliHFMLResponseDstartoD0pi& source);

protected:
    virtual void SetMapOfVariables(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int /*masshypo*/);

    /// \cond CLASSIMP
    ClassDef(AliHFMLResponseDstartoD0pi, 1); ///
    /// \endcond
};
#endif
