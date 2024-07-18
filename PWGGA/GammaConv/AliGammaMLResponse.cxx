// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Orgammanization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliGammaMLResponse
// \brief helper class to handle application of ML models trained with python libraries
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include "AliGammaMLResponse.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliGammaMLResponse);
/// \endcond

//________________________________________________________________
AliGammaMLResponse::AliGammaMLResponse() : AliMLResponse(), fVars{}
{
    //
    // Default constructor
    //
}

//________________________________________________________________
AliGammaMLResponse::AliGammaMLResponse(const Char_t *name, const Char_t *title, 
                                 const std::string configfilepath) : AliMLResponse(name, title),
                                                                     fVars{}
{
    //
    // Standard constructor
    //
    SetConfigFilePath(configfilepath);
}

//________________________________________________________________
AliGammaMLResponse::~AliGammaMLResponse()
{
    //
    // Destructor
    //
}

//--------------------------------------------------------------------------
AliGammaMLResponse::AliGammaMLResponse(const AliGammaMLResponse &source) : AliMLResponse(source),
                                                                  fVars(source.fVars)
{
    //
    // Copy constructor
    //
}

AliGammaMLResponse &AliGammaMLResponse::operator=(const AliGammaMLResponse &source)
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
double AliGammaMLResponse::Predict(AliAODConversionPhoton *cand, AliVEvent* fInputEvent, AliConversionPhotonCuts* fiPhotonCut, AliV0ReaderV1* fV0Reader)
{
    SetMapOfVariables(cand, fInputEvent, fiPhotonCut, fV0Reader);
    if (fVars.empty())
    {
        AliWarning("Map of features empty!");
        return -999.;
    }

    return Predict(cand->Pt(), fVars);
}

//________________________________________________________________
bool AliGammaMLResponse::IsSelected(double &prob, AliAODConversionPhoton *cand, AliVEvent* fInputEvent, AliConversionPhotonCuts* fiPhotonCut, AliV0ReaderV1* fV0Reader)
{   
    SetMapOfVariables(cand, fInputEvent, fiPhotonCut, fV0Reader);
    if (fVars.empty())
    {
        AliWarning("Map of features empty!");
        return false;
    }

    return IsSelected(cand->Pt(), fVars, prob);
}

//________________________________________________________________
bool AliGammaMLResponse::PredictMultiClass(std::vector<double> &outScores, AliAODConversionPhoton *cand, AliVEvent* fInputEvent, AliConversionPhotonCuts* fiPhotonCut, AliV0ReaderV1* fV0Reader)
{
    SetMapOfVariables(cand, fInputEvent, fiPhotonCut, fV0Reader);
    if (fVars.empty())
    {
        AliWarning("Map of features empty!");
        return false;
    }

    return PredictMultiClass(cand->Pt(), fVars, outScores);
}

//________________________________________________________________
bool AliGammaMLResponse::IsSelectedMultiClass(std::vector<double> &outScores, AliAODConversionPhoton *cand, AliVEvent* fInputEvent, AliConversionPhotonCuts* fiPhotonCut, AliV0ReaderV1* fV0Reader)
{   
    SetMapOfVariables(cand, fInputEvent, fiPhotonCut, fV0Reader);
    if (fVars.empty())
    {
        AliWarning("Map of features empty!");
        return false;
    }

    return IsSelectedMultiClass(cand->Pt(), fVars, outScores);
}
