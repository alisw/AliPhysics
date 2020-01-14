#ifndef ALIHFMLRESPONSEDPLUSTOKPIPI_H
#define ALIHFMLRESPONSEDPLUSTOKPIPI_H

// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliHFMLResponseDplustoKpipi
// \brief helper class to handle application of ML models for D+ analyses trained
// with python libraries
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include "AliHFMLResponse.h"

class AliHFMLResponseDplustoKpipi : public AliHFMLResponse
{
public:
    AliHFMLResponseDplustoKpipi();
    AliHFMLResponseDplustoKpipi(const Char_t *name, const Char_t *title, const std::string configfilepath);
    virtual ~AliHFMLResponseDplustoKpipi();

    AliHFMLResponseDplustoKpipi(const AliHFMLResponseDplustoKpipi &source);
    AliHFMLResponseDplustoKpipi& operator=(const AliHFMLResponseDplustoKpipi& source);

protected:
    virtual void SetMapOfVariables(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int /*masshypo*/);

    /// \cond CLASSIMP
    ClassDef(AliHFMLResponseDplustoKpipi, 1); ///
    /// \endcond
};
#endif
