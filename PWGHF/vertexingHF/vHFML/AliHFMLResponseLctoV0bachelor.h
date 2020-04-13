#ifndef ALIHFMLRESPONSELCTOV0BACHELOR_H
#define ALIHFMLRESPONSELCTOV0BACHELOR_H

// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliHFMLResponseLctoV0bachelor
// \brief helper class to handle application of ML models for Ds analyses trained
// with python libraries
// \authors:
// L. Vermunt, luuk.vermunt@cern.ch
/////////////////////////////////////////////////////////////////////////////////////////

#include "AliHFMLResponse.h"

class AliHFMLResponseLctoV0bachelor : public AliHFMLResponse
{
public:
  AliHFMLResponseLctoV0bachelor();
  AliHFMLResponseLctoV0bachelor(const Char_t *name, const Char_t *title, const std::string configfilepath);
  virtual ~AliHFMLResponseLctoV0bachelor();
  
  AliHFMLResponseLctoV0bachelor(const AliHFMLResponseLctoV0bachelor &source);
  AliHFMLResponseLctoV0bachelor& operator=(const AliHFMLResponseLctoV0bachelor& source);
  
protected:
  virtual void SetMapOfVariables(AliAODRecoDecayHF *cand, double bfield, AliAODPidHF *pidHF, int /*masshypo*/);
  
  /// \cond CLASSIMP
  ClassDef(AliHFMLResponseLctoV0bachelor, 1); ///
  /// \endcond
};
#endif
