#ifndef TMEMC_H
#define TMEMC_H

/******************************************************************************
*                                                                             *
* TMEMC                                                                       *
*                                                                             *
* EMC Module Class                                                            *
*                                                                             *
* Author: Timo Alho                                                           *
*                                                                             *
******************************************************************************/

#include "TMCal.h"

class TMEMC: public TMCal {

 RQ_OBJECT("TMEMC")

 public:
 
 TMEMC(const TGWindow *p, UInt_t w, UInt_t h);
 
 ClassDef(TMEMC, 1); // EmCal Module Base Class

};

#endif
