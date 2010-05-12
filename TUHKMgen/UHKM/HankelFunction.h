//////////////////////////////////////////////////////////////////////////////////       
//                                                                              //
//        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna     //
//      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru  //
//                           November. 2, 2005                                  //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

#ifndef HANKELFUNCTION_INCLUDED
#define HANKELFUNCTION_INCLUDED

#include <Rtypes.h>

Double_t HankelK0(Double_t x);
Double_t HankelK1(Double_t x);
// compute modified Hankel function of the second,...,order
Double_t HankelKn(Int_t n, Double_t x);

#endif
