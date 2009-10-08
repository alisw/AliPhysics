//-*- Mode: C++ -*-
// $Id: AliHLTPHOSConstants.h 34622 2009-09-04 13:22:01Z odjuvsla $

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2006                                       *
 *                                                                        * 
 * Author: Per Thomas Hille perthi@fys.uio.no for the ALICE DCS Project.  *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                * 
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#ifndef ALIHLTEMCALCONSTANTS_H
#define ALIHLTEMCALCONSTANTS_H

namespace EmcalHLTConst
{
  //  const int NZROWSRCU     =   56;            /**<Number of rows per RCU*/       
  // const int NXCOLUMNSRCU  =   16; 
  const int NZROWSMOD      =  48;            /**<Number of rows per module*/       
  const int NXCOLUMNSMOD   =  24;            /**<Number of columns per module*/ 

  const int NRCUSPERSECTOR = 4;

  /*
  const int MIN_DDL_NUMBER = 4608;
  const int MAX_DDL_NUMBER = 4631;
  const int N_DDLS =  MAX_DDL_NUMBER - MIN_DDL_NUMBER +1;
  */

#ifndef __CINT__
  const unsigned char PFVECTORDIR[] = "/HLT/PHOS/PFVectors";
#endif
  const int PFDEFAULTNSAMPLES = 70;
  const int PFDEFAULTSTARTINDEX = 0;
  const double DEFAULTTAU = 0.2;                /**<Assume that the signal rise time of the altrp pulses is 2 us (nominal value of the electronics)*/
  const int DEFAULTFS = 10;                  /**<Assume that the signal is samples with 10 MHZ samle rate*/
  const int NMODULES    =      13;           /**<Number of modules of the EMCAL detector*/
  //  const int NRCUS       =      2;            /**<Number of RCUs per Module*/
  const int NRCUSPERMODULE =  2 ;            /**<Number of RCUs per Module*/
  const int NFEECS         =  9;             /**<Number of Frontend cards per branch*/
}


#endif


