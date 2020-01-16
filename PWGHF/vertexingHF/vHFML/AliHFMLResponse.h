#ifndef ALIHFMLRESPONSE_H
#define ALIHFMLRESPONSE_H

// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliHFMLResponse
// \brief helper class to handle application of ML models trained with python libraries
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <map>

#include "AliMLResponse.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODPidHF.h"

class AliHFMLResponse : public AliMLResponse
{
public:
    AliHFMLResponse();
    AliHFMLResponse(const Char_t *name, const Char_t *title, const std::string configfilepath);
    virtual ~AliHFMLResponse();

    AliHFMLResponse(const AliHFMLResponse &source);
    AliHFMLResponse& operator=(const AliHFMLResponse& source);

    /// methods to get ML response
    using AliMLResponse::IsSelected; // exposes function from mother class
    bool IsSelected(double &prob, AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF = nullptr, int masshypo = 0);
    using AliMLResponse::Predict; // exposes function from mother class
    double Predict(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF = nullptr, int masshypo = 0);

    /// method to get variable (feature) from map
    double GetVariable(std::string name = "") {return fVars[name];}

protected:
    /// method used to define map of name <-> variables (features) --> to be implemented for each derived class
    virtual void SetMapOfVariables(AliAODRecoDecayHF * /*cand*/, double /*bfield*/, AliAODPidHF * /*pidHF*/, int /*masshypo*/) { return; }

    std::map<std::string, double> fVars;       /// map of variables (features) that can be used for the ML model application

    /// \cond CLASSIMP
    ClassDef(AliHFMLResponse, 2); ///
    /// \endcond
};
#endif
