#ifndef ALIGCONE_H
#define ALIGCONE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGTube.h"

class AliGCone: public AliGTube {

    protected:
        Float_t fRmax2;  /* outside radius at the high z limit */
        Float_t fRmin2;  /* inside radius at the high z limit  */

    public:
        AliGCone(); /* Default Constructor */
        AliGCone( Text_t *name, Text_t *title, Float_t dz, Float_t rmin1, Float_t rmax1, Float_t rmin2, Float_t rmax2 );
        AliGCone( Text_t *name, Text_t *title, Float_t dz, Float_t rmax1, Float_t rmax2=0 );
        AliGCone(AliGCone *cone);
	virtual ~AliGCone(); /* Destructor */

        Float_t GetRmin2() {return fRmin2;}
        Float_t GetRmax2() {return fRmax2;}

        virtual void  DrawShape(Option_t *option); // *MENU*
        virtual void  Draw(Option_t *option);
        virtual void  SetPoints(Float_t *buff);

    ClassDef(AliGCone,1) // Simple cone class
};
#endif
