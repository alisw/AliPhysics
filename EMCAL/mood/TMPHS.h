#ifndef TMPHS_H
#define TMPHS_H

/******************************************************************************
*                                                                             *
* TMPHS                                                                       *
*                                                                             *
* PHS Module Class                                                            *
*                                                                             *
* Author: Timo Alho                                                           *
*                                                                             *
******************************************************************************/

#include "TMCal.h"


class TMPHS: public TMCal {

  RQ_OBJECT("TMPHS")

  public:

  TMPHS(const TGWindow *p, UInt_t w, UInt_t h);
 
  ClassDef(TMPHS, 1);
  
};

#endif
