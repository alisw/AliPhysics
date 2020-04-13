#ifndef ALIHFMLRESPONSEDSTOKKPI_H
#define ALIHFMLRESPONSEDSTOKKPI_H

// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliHFMLResponseDstoKKpi
// \brief helper class to handle application of ML models for Ds analyses trained
// with python libraries
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include "AliHFMLResponse.h"

class AliHFMLResponseDstoKKpi : public AliHFMLResponse
{
public:
    AliHFMLResponseDstoKKpi();
    AliHFMLResponseDstoKKpi(const Char_t *name, const Char_t *title, const std::string configfilepath);
    virtual ~AliHFMLResponseDstoKKpi();

    AliHFMLResponseDstoKKpi(const AliHFMLResponseDstoKKpi &source);
    AliHFMLResponseDstoKKpi& operator=(const AliHFMLResponseDstoKKpi& source);

protected:
    virtual void SetMapOfVariables(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int masshypo);

    /// \cond CLASSIMP
    ClassDef(AliHFMLResponseDstoKKpi, 1); ///
    /// \endcond
};
#endif
