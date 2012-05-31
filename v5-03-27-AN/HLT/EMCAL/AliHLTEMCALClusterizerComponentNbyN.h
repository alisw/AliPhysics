
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Anders Knospe <anders.knospe@cern.ch>                         *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#ifndef ALIHLTEMCALCLUSTERIZERCOMPONENTNBYN_H
#define ALIHLTEMCALCLUSTERIZERCOMPONENTNBYN_H

#include "AliHLTEMCALClusterizerComponent.h"


class AliHLTEMCALClusterizerComponentNbyN : public AliHLTEMCALClusterizerComponent
{

public:

    /** Constructor */
    AliHLTEMCALClusterizerComponentNbyN();

    /** Destructor */
    virtual ~AliHLTEMCALClusterizerComponentNbyN();

    /** interface function, see @ref AliHLTComponent for description */
    const char* GetComponentID();

    /** interface function, see @ref AliHLTComponent for description */
    AliHLTComponent* Spawn();

protected:

    /** interface function, see @ref AliHLTComponent for description */
    int ScanConfigurationArgument ( int argc, const char** argv );

    /** interface function, see @ref AliHLTComponent for description */
    virtual int DoInit(int argc, const char** argv);

      /** Copy constructor,  not implemented */
    AliHLTEMCALClusterizerComponentNbyN(const AliHLTEMCALClusterizerComponent& );;
  
  /** Assignment operator, not implemented */
    AliHLTEMCALClusterizerComponentNbyN & operator = (const AliHLTEMCALClusterizerComponentNbyN);

    
    ClassDef(AliHLTEMCALClusterizerComponent, 0);
};

#endif // ALIHLTEMCALCLUSTERIZERCOMPONENTNBYN_H
