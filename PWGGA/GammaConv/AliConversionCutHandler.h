/**************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/
#ifndef ALICONVERSIONCUTHANDLER_H
#define ALICONVERSIONCUTHANDLER_H

#include <Rtypes.h>
class TString;

class AliConversionCutHandler {
public:
	AliConversionCutHandler(Int_t nMax=10);
	AliConversionCutHandler(const AliConversionCutHandler &ref);
	AliConversionCutHandler &operator=(const AliConversionCutHandler &ref);
	virtual ~AliConversionCutHandler();

	void AddCut(TString eventCut, TString photonCut, TString mesonCut);
	void AddCut(TString eventCut, TString photonCut, TString mesonCut, TString clusterCut);

	Bool_t AreValid() const {return fValidCuts;}
	Int_t GetNCuts() const { if(fValidCuts) return fNCuts; else return 0; }
	TString GetEventCut(Int_t i) const;
	TString GetPhotonCut(Int_t i) const;
	TString GetMesonCut(Int_t i) const;
	TString GetClusterCut(Int_t i) const;

private:
	Bool_t fValidCuts;              ///< Test if all cuts are valid
	Int_t fNCuts;                   ///< Number of cuts handled so far
	Int_t fNMaxCuts;                ///< Max. number of cuts that can be handled
    TString* fEventCutArray;        ///< Array of event cuts
    TString* fPhotonCutArray;       ///< Array of photon cuts
    TString* fMesonCutArray;        ///< Array of meson cuts
    TString* fClusterCutArray;      ///< Array of calo cluster cuts
};

#endif /* ALICONVERSIONCUTHANDLER_H */
