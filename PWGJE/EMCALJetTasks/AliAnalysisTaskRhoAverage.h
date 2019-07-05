/************************************************************************************
 * Copyright (C) 2012, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef ALIANALYSISTASKRHOAVERAGE_H
#define ALIANALYSISTASKRHOAVERAGE_H

// $Id$

#include "AliAnalysisTaskRhoBase.h"

/**
 * @class AliAnalysisTaskRhoAverage
 * @brief Calculation of rho, method: median all particle pt / multiplicity density.
 * @ingroup PWGJEBASE
 * @author Salvatore Aiola
 * @since May 24, 2012
 */
class AliAnalysisTaskRhoAverage : public AliAnalysisTaskRhoBase {

 public:
  /**
   * @brief Default constructor.
   */
  AliAnalysisTaskRhoAverage();

  /**
   * @brief Constructor.
   * @param name Name of the rho task
   * @param histo If true QA/Debug histograms are created
   */
  AliAnalysisTaskRhoAverage(const char *name, Bool_t histo=kFALSE);
  
  /**
   * @brief Destructor
   */
  virtual ~AliAnalysisTaskRhoAverage() {}

  void             SetRhoType(Int_t t)             { fRhoType       = t    ; }
  void             SetExcludeLeadPart(UInt_t n)    { fNExclLeadPart = n    ; }
  void             SetUseMedian(Bool_t b=kTRUE)    { fUseMedian     = b    ; }
  
 protected:
  void             ExecOnce();

  /**
   * @brief Run the analysis.
   * @return Always true
   */
  Bool_t           Run();

  Int_t            fRhoType       ; ///< rho type: 0 = charged+neutral, 1 = charged, 2 = neutral
  UInt_t           fNExclLeadPart ; ///< number of leading particles to be excluded from the median calculation
  Bool_t           fUseMedian     ; ///< whether or not use the median to calculate rho (mean is used if false)
  Double_t         fTotalArea     ; //!<! total area

  AliAnalysisTaskRhoAverage(const AliAnalysisTaskRhoAverage&);             // not implemented
  AliAnalysisTaskRhoAverage& operator=(const AliAnalysisTaskRhoAverage&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoAverage, 4); // Rho task
};
#endif
