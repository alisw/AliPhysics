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

#include "AliHFMLResponse.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliHFMLResponse);
/// \endcond

//________________________________________________________________
AliHFMLResponse::AliHFMLResponse() : AliMLResponse(),
                                     fVars{}
{
    //
    // Default constructor
    //
}

//________________________________________________________________
AliHFMLResponse::AliHFMLResponse(const Char_t *name, const Char_t *title, 
                                 const std::string configfilepath) : AliMLResponse(name, title),
                                                                     fVars{}
{
    //
    // Standard constructor
    //

    SetConfigFilePath(configfilepath);
}

//________________________________________________________________
AliHFMLResponse::~AliHFMLResponse()
{
    //
    // Destructor
    //
}

//--------------------------------------------------------------------------
AliHFMLResponse::AliHFMLResponse(const AliHFMLResponse &source) : AliMLResponse(source),
                                                                  fVars(source.fVars)
{
    //
    // Copy constructor
    //
}

AliHFMLResponse &AliHFMLResponse::operator=(const AliHFMLResponse &source)
{
    //
    // assignment operator
    //
    if (&source == this)
        return *this;

    AliMLResponse::operator=(source);
    fVars = source.fVars;

    return *this;
}

//________________________________________________________________
double AliHFMLResponse::Predict(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int masshypo)
{
    SetMapOfVariables(cand, bfield, pidHF, masshypo);
    if (fVars.empty())
    {
        AliWarning("Map of features empty!");
        return -999.;
    }

    return Predict(cand->Pt(), fVars);
}

//________________________________________________________________
bool AliHFMLResponse::IsSelected(double &prob, AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int masshypo)
{   
    SetMapOfVariables(cand, bfield, pidHF, masshypo);
    if (fVars.empty())
    {
        AliWarning("Map of features empty!");
        return -999.;
    }

    return IsSelected(cand->Pt(), fVars, prob);
}
