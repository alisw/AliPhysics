/************************************************************************************
 * Copyright (C) 2019, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#ifndef __ALILUNDPLANEHELPER_H__
#define __ALILUNDPLANEHELPER_H__
#include <TObject.h>
#include <exception>
#include <string>
#include <vector>

class AliEmcalJet;
class AliParticleContainer;
class AliClusterContainer;

namespace PWGJE {

namespace EMCALJetTasks {

/**
 * @class AliLundPlaneException
 * @brief Exception handling errors during iterative declustering or access to incomplete lund parameter set
 * @ingroup PWGJEBASE
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since June 5, 2019
 * 
 * Handling all errors occuring during the iterative declustering process. The kinds of error can appear:
 * - Fastjet error during declustering
 * - Access to incomplete lund plane parameter set
 */
class AliLundPlaneException : public std::exception {
public:

  /**
   * @brief Error types for the AliLundPlaneException
   */
  enum ErrorType_t {
    kFastjetError,      ///< Exception coming from FastJet
    kVertexNotSet,      ///< Vertex position missing
    kParamError,        ///< Access to uninitialized parameter
    kUndef              ///< Error source not defined
  };

  /**
   * @brief Dummy constructor
   */
  AliLundPlaneException(): fErrorType(kUndef), fErrorMessage("") {}

  /**
   * @brief Constructor with error message 
   * @param errtyep Type of the error
   * @param error Message content
   */
  AliLundPlaneException(ErrorType_t errtype, const char *error): fErrorType(errtype), fErrorMessage(error) {}

  /**
   * @brief Destructor
   */
  virtual ~AliLundPlaneException() throw() {}

  /**
   * @brief Error message of the exception, define when exception is thrown
   * @return Message content 
   */
  const char *what() const throw() { return fErrorMessage.data(); };

  /**
   * @brief Get type of the error
   * @return Type of the error
   */
  ErrorType_t GetErrorType() const throw() {return fErrorType; }

private:
  ErrorType_t         fErrorType;            ///< Type of the error
  std::string         fErrorMessage;         ///< Error Message
};

/**
 * @class AliLundPlaneParameters
 * @brief Container for iterative declustering
 * @ingroup PWGJEBASE
 * @authoer Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since June 5, 2019
 * 
 * The class containing the paramters of a splitting obtained during
 * iterative declustering. Parameters are:
 * - ln(\f$\Delta\f$R)  -  angle between subjets
 * - ln(\f$p_{t,rel}\f$)
 * - \f$p_{t}\f$ of the lower subjet
 * - Splitting number in the cluster tree
 */
class AliLundPlaneParameters {
public:

  /**
   * @brief Dummy constructor, for I/O purpose, should not be used
   */
  AliLundPlaneParameters() :
    fLnDeltaR(0),
    fLnPtrel(0),
    fPtLower(0),
    fNSplitting(0),
    fInitialized(false)
  {

  }

  /**
   * @brief Construct a new AliLundPlane parameter set and initialize it with values
   * 
   * @param lndeltaR  ln(\f$\Delta\f$R)
   * @param lnprel ln(\f$p_{t,rel}\f$)
   * @param ptlower \f$p_{t}\f$ of the lower subjet
   * @param nsplitting Splitting number in the cluster tree
   */
  AliLundPlaneParameters(Double_t lndeltaR, Double_t lnprel, Double_t ptlower, int nsplitting):
    fLnDeltaR(lndeltaR),
    fLnPtrel(lnprel),
    fPtLower(ptlower),
    fNSplitting(nsplitting),
    fInitialized(true)
  {
  }

  /**
   * @brief Destructor
   */
  ~AliLundPlaneParameters() {}

  /**
   * @brief Get the angle between the subjets (expressed as ln(\f$\Delta\f$R)
   * @return Angle between subjets
   * @throw AliLundPlaneException
   */
  Double_t GetLnDeltaR() const { ProtectAccessUninit(); return fLnDeltaR; }

  /**
   * @brief Get ln(\f$p_{t,rel}\f$)
   * @return ln(\f$p_{t,rel}\f$)
   * @throw AliLundPlaneException
   */
  Double_t GetLnPtrel() const { ProtectAccessUninit(); return fLnPtrel; }
  
  /**
   * @brief Get the \f$p_{t}\f$ of the lower subjet
   * @return \f$p_{t}\f$ of the lower subjet
   * @throw AliLundPlaneException
   */
  Double_t GetPtLower() const { ProtectAccessUninit(); return fPtLower; }

  /**
   * @brief Get hte splitting number in the cluster tree
   * @return Splitting number in the cluster tree
   * @throw AliLundPlaneException
   */
  Int_t GetNSplittings() const { ProtectAccessUninit(); return fNSplitting; } 

  /**
   * @brief Check whether parameters are initialized
   * @return True if the parameters are initialized (with values), false otherwise
   * @throw AliLundPlaneException
   */
  Bool_t IsInitialized() const { return fInitialized; }

private:
  /**
   * @brief Protection against uninitialized values
   * @throw ALiLundPlaneException
   * 
   * Function throwing an exception in case the lund plane 
   * parameter set is uninitialized. To be called in functions
   * accessing any of the parameters
   */
  void ProtectAccessUninit() const {
    if(!fInitialized) throw AliLundPlaneException(AliLundPlaneException::kParamError, "Access to uninitialized splitting");
  }

  Double_t        fLnDeltaR;            ///< ln(\f$\DeltaR\f$)
  Double_t        fLnPtrel;             ///< ln(\f$p_{t,rel}\f$)
  Double_t        fPtLower;             ///< \f$p_{t}\f$ of the lower subjet
  Int_t           fNSplitting;          ///< Splitting number in the cluster tree
  Bool_t          fInitialized;         ///< Init status of the parameters
};

/**
 * @class AliLundPlaneData
 * @brief Container for iterative declustering
 * @ingroup PWGJEBASE
 * @authoer Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since June 5, 2019
 * 
 * Container class handling the parameters for all splittings of a jet. Parameters
 * are defined in AliLundPlaneParameters. Each set of parameters belongs to a new 
 * splitting. Splittings are added with InsertSplitting. GetSplittings returns a 
 * list of all splittings found for the jet.
 */
class AliLundPlaneData {
public:
  /**
   * @brief Constructor
   */
  AliLundPlaneData() {}

  /**
   * @brief Destructor
   */
  ~AliLundPlaneData() {}

  /**
   * @brief Insert new splitting into the splitting container.
   * @param splitting Parameters for the new splitting
   */
  void InsertSplitting(const AliLundPlaneParameters &splitting) { fSplittings.push_back(splitting); }

  /**
   * @brief Get list of all splittings handled by this container
   * @return List of splittings
   */
  const std::vector<AliLundPlaneParameters> &GetSplittings() const { return fSplittings; }

private:
    std::vector<AliLundPlaneParameters>             fSplittings;        ///< List of splittings
};

/**
 * @class AliLundPlaneHelper
 * @brief Helper class for iterative declustering of a jet
 * @ingroup PWGJEBASE
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory 
 */
class AliLundPlaneHelper : public TObject {
public:
    /**
     * @brief Constructor
     */
    AliLundPlaneHelper();

    /**
     * @brief Destructor
     */
    virtual ~AliLundPlaneHelper() {}

    /**
     * @brief Set hard cutof
     * 
     * @param cutoff 
     */
    void SetHardCutoff(Double_t cutoff) { fHardCutoff = cutoff; }

    /**
     * @brief Run iterative declustering of the jet
     * 
     * Reclusters jet with the Cambridge/Aachen algorithm and traverses the cluster tree.
     * For each splitting the lund plane parameters defined in AliLundPlaneParamters are
     * evaluated and added to the splitting container. The declustering is done always on the 
     * harder subjet until the abort condition (no more child subjets or subjet > hard cutoff) 
     * is fulfilled. Needs track/particle/cluster containers in order to have acces to the jet 
     * constituents.
     * 
     * @param jet Jet to be tested
     * @param tracks Container with tracks / particles (optional only in case of neutral jets)
     * @param clusters Container with calorimeter clusters (optional, for det. level full jets)
     * @param vertexpos Position of the primary vertex (optional, mandatory in case a cluster container is set)
     * @return Container with lund plane parameters for each splitting
     */
    AliLundPlaneData Evaluate(const AliEmcalJet &jet, const AliParticleContainer *tracks, const AliClusterContainer *clusters = nullptr, Double_t *vertexpos = nullptr);    

private:
    Double_t                fHardCutoff;      ///< Hard cutoff (abort condition for iterative declustering)

    ClassDef(AliLundPlaneHelper, 1);
};

}

}
#endif
