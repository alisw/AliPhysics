
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Oystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#ifndef ALIHLTPHOSCLUSTERIZERCOMPONENTNBYN_H
#define ALIHLTPHOSCLUSTERIZERCOMPONENTNBYN_H

#include "AliHLTPHOSClusterizerComponent.h"


class AliHLTPHOSClusterizerComponentNbyN : public AliHLTPHOSClusterizerComponent
{

public:

    /** Constructor */
    AliHLTPHOSClusterizerComponentNbyN();

    /** Destructor */
    virtual ~AliHLTPHOSClusterizerComponentNbyN();

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
    AliHLTPHOSClusterizerComponentNbyN(const AliHLTPHOSClusterizerComponent& );;
  
  /** Assignment operator, not implemented */
    AliHLTPHOSClusterizerComponentNbyN & operator = (const AliHLTPHOSClusterizerComponentNbyN);

    
    ClassDef(AliHLTPHOSClusterizerComponent, 0);
};

#endif // ALIHLTPHOSCLUSTERIZERCOMPONENTNBYN_H
