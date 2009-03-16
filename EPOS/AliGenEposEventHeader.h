/*
 * AliGenEposEventHeader.h
 * 
 * Header for EPOS generated event.
 *
 *      Author: Piotr Ostrowski
 */

#ifndef ALIGENEPOSEVENTHEADER_H_
#define ALIGENEPOSEVENTHEADER_H_

#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"

class AliGenEposEventHeader : public AliGenEventHeader, public AliCollisionGeometry
{
public:
	AliGenEposEventHeader(const char* name);
	AliGenEposEventHeader();
	virtual ~AliGenEposEventHeader() {}


protected:

private:
	ClassDef(AliGenEposEventHeader,1)
};


#endif /* ALIGENEPOSEVENTHEADER_H_ */
