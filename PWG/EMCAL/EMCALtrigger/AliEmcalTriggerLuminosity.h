/************************************************************************************
 * Copyright (C) 2021, Copyright Holders of the ALICE Collaboration                 *
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
#include <TObject.h>
#include <exception>
#include <map>
#include <sstream>
#include <string>
#include <vector>

class TH1;
class TH2;
class TList;

namespace PWG {

namespace EMCAL {

/**
 * @class AliEmcalTriggerLuminosity
 * @brief Calculator for the integrated fuminosity of EMCAL triggers based on normalization task output
 * @ingroup EMCALTRGFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Sept 30, 2021 
 */
class AliEmcalTriggerLuminosity : public TObject {
public:

  /**
   * @class UninitException
   * @brief Handling errors in evaluation due to uninitialized data
   * 
   * Handled fields can be
   * - Collision system (relevant for livetime correction)
   * - Year (relevant for reference cross section)
   * - Raw event counter
   * - Cluster counters (relevant for livetime correction)
   * Uninitialzied counters are handled specifically to the dataset
   * while all other fields are handled in the Evaluate function.
   */
  class UninitException : public std::exception {
    public:
      /**
       * @brief Constructor
       * @param fields List of required information which were uninitialized
       */
      UninitException(const std::vector<std::string> &fields) : std::exception(), fFields(fields), fMessage("") {
        buildErrorMessage();
      }

      /**
       * @brief Destructor
       */
      virtual ~UninitException() throw() {}

      /**
       * @brief Display error message
       * @return Error message
       */
      virtual const char *what() const throw() { return fMessage.c_str(); }

      /**
       * @brief Access to list of uninitialized fields
       * @return List of required information which were uninitialized
       */
      const std::vector<std::string> &getFields() const { return fFields; }

    private:
      std::vector<std::string> fFields;   ///< Fields which were required but missing
      std::string fMessage;               ///< Buffer for error message

      /**
       * @brief Internal helper building the error message
       * 
       * Displaying all fields which were required but not set. Handled
       * as helper function to move this part to the implementation file
       * in order to be able to make use of c++11 code which cannot be handled
       * by stupid ROOT5.
       */
      void buildErrorMessage();
  };

  /**
   * @class InputDataException
   * @brief Handling errors due to invalid input files (task output not from normalizatinon task)
   */
  class InputDataException : public std::exception {
    public:
      /**
       * @brief Constructor
       * @param filename Name of the file raising the exception
       * @param directory Name of the directory raising the exception
       */
      InputDataException(const char *filename, const char *directory): std::exception(), fFilename(filename), fMessage("") {
        std::stringstream msgbuilder;
        msgbuilder << "No trigger norm data found in file " << fFilename << ", directory " << fDirname;
        fMessage = msgbuilder.str();
      }

      /**
       * @brief Destructor
       */
      virtual ~InputDataException() throw() {}

      /**
       * @brief Get name of the file raising the exception
       * @return Filename
       */
      const std::string &getFilename() const { return fFilename; }

      /**
       * @brief Get the name of the directory raising the exception
       * @return Directory name  
       */
      const std::string &getDirname() const { return fDirname; }

      /**
       * @brief Display error message
       * @return Error message
       */
      virtual const char *what() const throw() { return fMessage.c_str(); }

    private:
      std::string fFilename;      ///< Name of the file read
      std::string fDirname;
      std::string fMessage;       ///< Buffer for error message
  };

  /**
   * @class TriggerNotFoundException
   * @brief Handling errors requesting non-existing trigger classes
   */
  class TriggerNotFoundException : public std::exception {
    public:
      /**
       * @brief Constructor
       * @param trigger Trigger class raising the exception
       */
      TriggerNotFoundException(const char *trigger) : std::exception(), fTriggerClass(trigger), fMessage() {
        std::stringstream msgbuilder;
        msgbuilder << "Luminosity not found for trigger " << fTriggerClass;
        fMessage = msgbuilder.str();
      }

      /**
       * @brief Destructor
       */
      virtual ~TriggerNotFoundException() throw() {}

      /**
       * @brief Display error message
       * @return Error message
       */
      virtual const char *what() const throw() { return fMessage.c_str(); }

      /**
       * @brief Access to trigger class raising the exception
       * @return Name of the trigger class
       */
      const std::string &getTriggerClass() const { return fTriggerClass; }

    private:
      std::string fTriggerClass;    ///< Trigger class raising the exception
      std::string fMessage;         ///< Buffer for error message
  };

  /**
   * @class AmbiguityException
   * @brief Handling of ambiguous information in certain fields (year, collision system)
   */
  class AmbiguityException : public std::exception { 
    public:
      /**
       * @brief Constructor
       * @param field Field with ambiguous information
       */
      AmbiguityException(const char *field) : std::exception(), fField(field), fMessage() {
        std::stringstream msgbuilder;
        msgbuilder << "Ambiguity in " <<  field;
        fMessage = msgbuilder.str();
      }

      /**
       * @brief Destructor
       */
      virtual ~AmbiguityException() throw() {}

      /**
       * @brief Display error message
       * @return Error message
       */
      virtual const char *what() const throw() { return fMessage.c_str(); }

    private:
      std::string fField;           ///< Fields (year, collision system) with ambiguous information
      std::string fMessage;         ///< Buffer for error message
  };
  
  /**
   * @enum CollisionType_t
   * @brief Collision systems supported by the normalization task 
   */
  enum CollisionType_t {
    kPP13TeV,         ///< pp, \f$ \sqrt{s}\f$ = 13 TeV
    kPP5TeV,          ///< pp, \f$ \sqrt{s}\f$ = 5.02 TeV
    kPPB5TeV,         ///< p-Pb, \f$ \sqrt{s_{NN}}\f$ = 5.02 TeV
    kPPB8TeV,         ///< p-Pb, \f$ \sqrt{s_{NN}}\f$ = 8.16 TeV
    kPBPB5TeV,        ///< Pb-Pb, \f$ \sqrt{s_{NN}}\f$ = 5.02 TeV
    kUnknown          ///< Unknown collision system
  };

  /**
   * @brief Construct a new Ali Emcal Trigger Luminosity object
   */
  AliEmcalTriggerLuminosity();

  /**
   * @brief Constructor, initializing from a file with normalization task output
   * @param filename Name of the file with the normalization task output
   * @param dirname  Directory of the normalization task output inside the file
   * @throw AmbiguityException Year or collision system contain ambiguous information
   * @throw InputDataException Input data is not from the normalization task
   */
  AliEmcalTriggerLuminosity(const char *filename, const char *dirname  = "EmcalTriggerNormtask");

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalTriggerLuminosity();

  /**
   * @brief Load normalization histograms from an input list
   * @param histlist List with the normalization histograms
   * @throw AmbiguityException Year or collision system contain ambiguous information
   */
  void LoadHistograms(const TList &histlist);

  /**
   * @brief Initialize with normalization task output from an input file
   * @param filename File with the normmalization task output
   * @param dirname Directory in file where to find the normalization task output
   * @throw AmbiguityException Year or collision system contain ambiguous information
   * @throw InputDataException Input data is not from the normalization task
   */
  void InitFromFile(const char *filename, const char *dirname);

  /**
   * @brief Evaluate all luminosities
   * @throw UninitException Several required information was not initialized
   * 
   * Evaluating luminosities based on year and collision system of the 
   * data used for the normalization task output. Luminosities are calculated
   * for all trigger classes in the normalization histogram
   */
  void Evaluate();

  /**
   * @brief Get the integrated luminosity evaluated for a certain trigger 
   * @param trigger Trigger for which to determine the luminosity
   * @return Integrated luminosity
   * @throw TriggerNotFoundException in case of access to unsupported trigger
   * 
   * Access to luminosities calculated in Evaluate.
   */
  double GetLuminosityForTrigger(const char *trigger);

  /**
   * @brief Set the collision system of the data used in the input file
   * @param ycoltypeear Collision system of the data used in the input file
   * 
   * Only necessary for old normalization task results. In new normalization task results
   * this is taken automatically from corresponding histograms
   */
  void SetCollisionType(CollisionType_t coltype) { fCollisionType = coltype; }

  /**
   * @brief Set the year of the data used in the input file
   * @param year Year of the data used in the input file
   * 
   * Only necessary for old normalization task results. In new normalization task results
   * this is taken automatically from corresponding histograms
   */
  void SetYear(int year) { fYear = year; }

protected:

  /**
   * @brief Luminosity evaluation method for pp, $\sqrt{s} = 13 TeV\f$, 2016-2018
   */
  void evaluatePP13TeV();

  /**
   * @brief Luminosity evaluation method for p-Pb, $\sqrt{s_{NN}} = 8.16 TeV\f$, 2016
   */
  void evaluatePPB8TeV();

  /**
   * @brief Determine year of the input data from year counter histogram
   * @param yearhist Year counter histogram
   * @return  Year of the dataset
   */
  int determineYear(const TH1 * const yearhist) const;

  /**
   * @brief Determine collision system of the input data from collision system counter histogram
   * @param collisionhist Collision system counter from normalization task
   * @return Collision system of the dataset
   */
  CollisionType_t determineCollisionType(const TH1 *const collisionhist) const;

  /**
   * @brief Converter for collision type to string
   * @param coltype CollisionType_t representation 
   * @return Collision type label
   */
  std::string getCollisionLabel(CollisionType_t coltype) const;

  /**
   * @brief Converter for collision type from string
   * @param collabel String representation of collision type
   * @return CollisionType_t representation
   */
  CollisionType_t getCollisionType(const std::string &collabel) const; 
  
  /**
   * @brief Get numnber of events recorded in a certain trigger cluster for a trigger class
   * @param clustercounter Input cluster counter histogram
   * @param clustername Name of the trigger cluster
   * @return Number of events
   */
  double getTriggerClusterCounts(TH1 * clustercounter, const std::string &clustername) const;

private:
  AliEmcalTriggerLuminosity(const AliEmcalTriggerLuminosity &);
  AliEmcalTriggerLuminosity &operator=(const AliEmcalTriggerLuminosity &);

  CollisionType_t fCollisionType;                       ///< Collision system of the input file
  int fYear;                                            ///< Year of the input file
  TH2 *fLuminosityHist;                                 ///< Raw luminosity histogram
  std::map<std::string, TH1 *> fClusterCounters;        ///< Cluster counter histogram for various trigger classes
  std::map<std::string, double> fLuminosities;          ///< Evaluated integrated luminosities for all trigger classes

  ClassDef(AliEmcalTriggerLuminosity, 1);
};

}

}