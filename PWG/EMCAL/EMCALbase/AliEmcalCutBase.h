/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef ALIEMCALCUTSBASE_H
#define ALIEMCALCUTSBASE_H

#include <TNamed.h>
#include <AliEmcalTrackSelResultPtr.h>

namespace PWG {

namespace EMCAL {

/**
 * @class AliEmcalCutBase
 * @brief Interface for a cut class returning selection status and user information
 * @ingroup EMCALCOREFW
 * @author: Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Dec 11, 2017
 * 
 * This class extends the functionality of AliVCuts in the way that the cut implementation
 * can add user information to the selection results: The method IsSelected of class AliVCuts
 * returns a bool, representing whether an object was selected or not. For some selections,
 * i.e. the hybrid track selection, where certain track types are defined, this is oversimplified.
 * Instead one wants to return a selection status (track selected or not) and a user information
 * which in case of hybrid tracks stores the hybrid track type. This is done using a smart pointer
 * approach implemented in AliEmcalTrackSelResultPtr. 
 * 
 * User classes inheriting from AliEmcalCutBase must implement the function IsSelected(const TObject *),
 * this time returning an AliEmcalTrackSelResultPtr with
 * - Selection Status
 * - Pointer to track processed
 * - User information (optional)
 */
class AliEmcalCutBase : public TNamed {
public:
  AliEmcalCutBase() {}
  AliEmcalCutBase(const char *name, const char *title);
  virtual ~AliEmcalCutBase() {}

  virtual AliEmcalTrackSelResultPtr IsSelected(TObject *o) = 0;

private:

  ClassDef(AliEmcalCutBase, 1);
};

}

}
#endif