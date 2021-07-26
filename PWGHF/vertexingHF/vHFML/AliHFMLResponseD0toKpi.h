#ifndef ALIHFMLRESPONSED0TOKPI_H
#define ALIHFMLRESPONSED0TOKPI_H

// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliHFMLResponseD0toKpi
// \brief helper class to handle application of ML models for D0 analyses trained
// with python libraries
// \authors:
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include "AliHFMLResponse.h"

class AliHFMLResponseD0toKpi : public AliHFMLResponse
{
public:
    AliHFMLResponseD0toKpi();
    AliHFMLResponseD0toKpi(const Char_t *name, const Char_t *title, const std::string configfilepath);
    virtual ~AliHFMLResponseD0toKpi();

    AliHFMLResponseD0toKpi(const AliHFMLResponseD0toKpi &source);
    AliHFMLResponseD0toKpi& operator=(const AliHFMLResponseD0toKpi& source);

protected:
    virtual void SetMapOfVariables(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int /*masshypo*/);

    /// \cond CLASSIMP
    ClassDef(AliHFMLResponseD0toKpi, 1); ///
    /// \endcond
};
#endif
