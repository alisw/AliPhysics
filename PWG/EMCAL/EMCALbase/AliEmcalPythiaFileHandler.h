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
#ifndef ALIEMCALPYTHIAFILEHANDLER_H
#define ALIEMCALPYTHIAFILEHANDLER_H

#include <exception>
#include <string>

namespace PWG {

namespace EMCAL {

/**
 * @struct AliEmcalPythiaCrossSectionInfo
 * @brief Container for cross section and Number of trials information in PYTHIA cross section file
 */
struct AliEmcalPythiaCrossSectionData {
  double fCrossSection;       ///< Cross section
  double fNTrials;            ///< Number of trials
};

/**
 * @class     AliEmcalPythiaFileHandler
 * @brief     Hander for access to PYTHIA cross section file
 * @ingroup   EMCALCOREFW
 * @author    Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since     Feb 19, 2021
 * 
 * Singleton instance caching the cross section and the number of trials 
 * from PYTHIA cross section files. Cached values can be obtained by all 
 * tasks in the train. Current files must be indicated via the function 
 * UpdateFile, which will fetch new values in case the file name differs 
 * from the file connected to the current cache.
 */
class AliEmcalPythiaFileHandler {
public:

  /**
   * @class FileNotFoundException
   * @brief Error handling for access to non-existing files
   */
  class FileNotFoundException : public std::exception {
    public:
      /**
       * @brief Constructor
       * @param filename Name of the file raising the exception
       */
      FileNotFoundException(const char *filename) : std::exception(), fFilename(filename), fMessage() {
        fMessage = "File " + fFilename + " not available or zombie";
      }
      /**
       * @brief Destructor
       */
      virtual ~FileNotFoundException() throw() {}

      /**
       * @brief Providing error message
       * @return Error message of the exception
       */
      virtual const char *what() const throw() {
        return fMessage.c_str();
      }

      /**
       * @brief Get the name of the file raising the exception
       * @return Name of the file raising the exception
       */
      const std::string &GetFilename() const throw() { return fFilename; }
      
    private:
      std::string fFilename;    ///< File raising the exception
      std::string fMessage;     ///< Cache for error message
  };

  /**
   * @class FileContentException
   * @brief Error handling for files with corrupted content
   */
  class FileContentException : public std::exception {
    public:
     /**
      * @brief Constructor
      * @param filename Name of the file raising the exception
      */
      FileContentException(const char *filename) : std::exception(), fFilename(filename), fMessage() {
        fMessage = "File " + fFilename + " does not contain cross section information";
      }
      
      /**
       * @brief Destructor
       */
      virtual ~FileContentException() throw() {}

      /*
       * @brief Providing error message
       * @return Error message of the exception
       */
      virtual const char *what() const throw() {
        return fMessage.c_str();
      }

      /**
       * @brief Get the name of the file raising the exception
       * @return Name of the file raising the exception
       */
      const std::string &GetFilename() const throw() { return fFilename; }

    private:
      std::string fFilename;  ///< File raising the exception
      std::string fMessage;   ///< Cache for error message
  };

  /**
   * @class UninitializedException
   * @brief Access to data where initialization failed
   */
  class UninitializedException : public std::exception {
    public:
      /**
       * @enum ValueType_t
       * @brief Specification of content that failed during access  
       */
      enum ValueType_t {
        kCrossSection,    ///< Cross section information
        kNtrials          ///< Number of trials information
      };

      /**
       * @brief Constructor
       * @param valtype Access type raising the exception
       */
      UninitializedException(ValueType_t valtype) : std::exception(), fValueType(valtype), fMessage() {
        std::string valname;
        switch(valtype){
          case ValueType_t::kCrossSection: valname = "cross section"; break;
          case ValueType_t::kNtrials: valname = "number of trials"; break;
        };
        fMessage = "Access to uninitialized value: " + valname;
      }

      /**
       * @brief Destructor
       */
      virtual ~UninitializedException() throw() {}

      /*
       * @brief Providing error message
       * @return Error message of the exception
       */
      virtual const char *what() const throw() {
        return fMessage.c_str();
      }

      /**
       * @brief Get the type of information raising the exception
       * @return Type of object raising the exception
       */
      ValueType_t getValueType() const throw() { return fValueType; }

    private:
      ValueType_t fValueType;   ///< Type of value raising the exception
      std::string fMessage;     /// Cache for error message
 };

  /**
   * @brief Get current instance of hte file handler
   * @return Current instance
   */ 
  static AliEmcalPythiaFileHandler *Instance();

  /**
   * @brief Destructor
   */
  ~AliEmcalPythiaFileHandler();

  /**
   * @brief Get cross section and number 
   * @param filename Full path to the file analysed
   * @return Cross section and number of trials in cross section file
   * @throw FileNotFoundException if requested file is not found
   * @throw FileContentException if requested file does not contain cross section information
   * @throw UninitializedExeption if the requested file is the same but the data is not properly initialized
   * 
   * Updating the content of the current cache based on the request of the
   * associated analysis file. The filename must be the full path to the
   * file analysed (either and AliESDs.root or AliAOD.root file), as the 
   * associated file is determined from the full path.
   */
  AliEmcalPythiaCrossSectionData GetCrossSectionAndNTrials(const char *filename);

private:

  /**
   *  @brief Constructor, private as singleton object
   */
  AliEmcalPythiaFileHandler();
  AliEmcalPythiaFileHandler(const AliEmcalPythiaFileHandler &);
  AliEmcalPythiaFileHandler &operator=(const AliEmcalPythiaFileHandler &);

  /**
   * @brief Updating internal cache for the new file requested
   * @param filename Name of the ROOT file currently analysed
   * @throw FileNotFoundException if requested file is not found
   * @throw FileContentException if requested file does not contain cross section information
   * 
   * Reading the connected cross section file connected to the 
   * file name required and updating the values stored in the 
   * internal cache. The filename must contain the full path
   * of the file currently analysed. The connected file is
   * determined based on the path of the file analysed and is
   * expected to be in the same directory. Following files are
   * used:
   * - pyxsec.root for AliESDs.root
   * - pyxsec_hists.root for AliAOD.root
   */
  void UpdateCache(const char *filename);

  /**
   * @brief Read cross section information from file pyxsec.root (ESD case)
   * @param pyxsecfile Name of the cross section file
   * @throw FileNotFoundException if requested file is not found
   * @throw FileContentException if requested file does not contain cross section information
   */
  void UpdateFromXsecFile(const char *pyxsecfile);

  /**
   * @brief Read cross section informatino from file pyxsec_hists.root (AOD case)
   * @param pyxsechistfile Name of the cross section file
   * @throw FileNotFoundException if requested file is not found
   * @throw FileContentException if requested file does not contain cross section information
   */
  void UpdateFromXsecHistFile(const char *pyxsechistfile);

  AliEmcalPythiaCrossSectionData fCrossSectionNrials;   ///< Cross section and number of trials in current file
  std::string   fCurrentFile;                           ///< Current file used for cross section handling;
  bool          fInitialized;                           ///< Check whehther handler has been initialized

  static AliEmcalPythiaFileHandler *fgInstance;  ///< File handler singleton instance

};

}

}

#endif