/************************************************************************************
 * Copyright (C) 2020, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef ALIEMCALTRIGGERALIAS_H
#define ALIEMCALTRIGGERALIAS_H

#include <iosfwd>
#include <TObject.h>
#include <TList.h>

// operator<< has to be forward declared carefully to stay in the global namespace so that it works with CINT.
// For generally how to keep the operator in the global namespace, See: https://stackoverflow.com/a/38801633
namespace PWG { namespace EMCAL { class AliEmcalTriggerAlias; } }
std::ostream & operator<< (std::ostream &in, const PWG::EMCAL::AliEmcalTriggerAlias &alias);

namespace PWG {

namespace EMCAL {

/**
 * @class AliEmcalTriggerAlias
 * @brief Class for trigger aliases
 * @ingroup EMCALTRGFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Jan. 20, 2020
 * 
 * A trigger alias is a collection of trigger classes connected to the same trigger
 * condition. In case the trigger condition is fulfilled, any of the trigger aliases
 * is accepted. In order to check for a specific class name check
 * 
 * ~~~.{cxx}
 * AliEmcalTriggerAlias myalias("EGA;EG1");
 * bool accepted = myalias.HasTriggerClass("EG1");  //true
 * bool notaccepted = myalias.HasTriggerClass("EG2");  //false
 * ~~~.{cxx}
 */
class AliEmcalTriggerAlias : public TObject {
public:

  /**
   * @brief Default I/O constructor
   */
  AliEmcalTriggerAlias();

  /**
   * @brief Constructor, creates a new trigger alias from a list of trigger names
   * @param triggernames List of the trigger classes handled by the alias, separated by ";"
   */
  AliEmcalTriggerAlias(const char *triggernames);


  /**
   * @brief Constructor, creates a new trigger alias from a list of trigger names
   * @param triggernames List of trigger classes handled
   */
  AliEmcalTriggerAlias(const TList &triggernames);

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalTriggerAlias() {}

  /**
   * @brief check whether trigger name is handled by the trigger alias
   * @param triggername Name of the trigger class to be checked
   * @return True if the trigger class is handled, false otherwise
   */
  Bool_t HasTriggerClass(const char *triggername) const;

  /**
   * @brief Set the trigger classes handled by the trigger alias
   * @param triggernames List of the trigger classes handled by the alias, separated by ";"
   */
  void SetTriggerClasses(const char *triggernames);

  /**
   * @brief Set the trigger classes handled by the trigger alias
   * @param triggernames List of trigger classes handled
   */
  void SetTriggerClasses(const TList &triggernames);

  /**
   * @brief Print trigger alias on a stream
   * @param stream Stream where the trigger alias is printed on
   * 
   * Helper function used in the streaming operator
   * to print the trigger classes handled by the alias
   */
  void PrintStream(std::ostream &stream) const;

private:

  /**
   * @brief Decode trigger classes handled by the alias
   * @param triggernames String representation of the trigger classes
   * 
   * Helper function used in the constructor or the setter
   * decoding the trigger classes handled by the trigger alias from
   * a string. Trigger classes are separated by ";"
   */
  void DecodeTriggerClasses(const char *triggernames);

  TList    fTriggerClasses;       ///< List of trigger classes handled by Alias

  ClassDef(AliEmcalTriggerAlias,1)
};

}

}

#endif