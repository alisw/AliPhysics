#ifndef ALIGammaMLRESPONSE_H
#define ALIGammaMLRESPONSE_H

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

#include <string>
#include <vector>
#include <map>

#include "../../ML/AliMLResponse.h"
#include "AliAODConversionPhoton.h"
#include "AliConversionPhotonCuts.h"
#include "AliV0ReaderV1.h"
#include "AliVTrack.h"

class AliGammaMLResponse : public AliMLResponse
{
public:
    AliGammaMLResponse();
    AliGammaMLResponse(const Char_t *name, const Char_t *title, const std::string configfilepath);
    virtual ~AliGammaMLResponse();

    AliGammaMLResponse(const AliGammaMLResponse &source);
    AliGammaMLResponse& operator=(const AliGammaMLResponse& source);

    /// methods to get ML response
    using AliMLResponse::IsSelected; // exposes function from mother class
    bool IsSelected(double &prob, AliAODConversionPhoton *cand, AliVEvent* fInputEvent, AliConversionPhotonCuts* fiPhotonCut, AliV0ReaderV1* fV0Reader);
    using AliMLResponse::Predict; // exposes function from mother class
    double Predict(AliAODConversionPhoton *cand, AliVEvent* fInputEvent, AliConversionPhotonCuts* fiPhotonCut, AliV0ReaderV1* fV0Reader);
    using AliMLResponse::IsSelectedMultiClass; // exposes function from mother class
    bool IsSelectedMultiClass(std::vector<double> &outScores, AliAODConversionPhoton *cand, AliVEvent* fInputEvent, AliConversionPhotonCuts* fiPhotonCut, AliV0ReaderV1* fV0Reader);
    using AliMLResponse::PredictMultiClass; // exposes function from mother class
    bool PredictMultiClass(std::vector<double> &outScores, AliAODConversionPhoton *cand, AliVEvent* fInputEvent, AliConversionPhotonCuts* fiPhotonCut, AliV0ReaderV1* fV0Reader);
    /// method to get variable (feature) from map
    double GetVariable(std::string name = "") {return fVars[name];}

protected:
    /// method used to define map of name <-> variables (features) --> to be implemented for each derived class
    virtual void SetMapOfVariables(AliAODConversionPhoton *cand, AliVEvent* fInputEvent, AliConversionPhotonCuts* fiPhotonCut, AliV0ReaderV1* fV0Reader) { return; }

    std::map<std::string, double> fVars;       /// map of variables (features) that can be used for the ML model application

    /// \cond CLASSIMP
    ClassDef(AliGammaMLResponse, 2); ///
    /// \endcond
};
#endif