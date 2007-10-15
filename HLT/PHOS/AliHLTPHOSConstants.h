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

#ifndef ALIHLTPHOSCONSTANTS_H
#define ALIHLTPHOSCONSTANTS_H

namespace PhosHLTConst
{
  const int MAX_HOSTS = 20;
  const int DEFAULT_EVENT_PORT = 42001;
  const int MAX_BIN_VALUE = 1023;
  const unsigned int HIGH_GAIN    =   1;
  const unsigned int LOW_GAIN     =   0;

  const unsigned int ALTRO_MAX_SAMPLES = 1008;                           /**<The maximum number of samples of the ALTRO*/
  //  const unsigned int ALTRO_MAX_TRALER_SIZE = 7;  
  //  const unsigned int  DDL_BLOCK_SIZE = 5;

  const unsigned int N_ZROWS_RCU     =   28;                    /**<Number of rows per module*/       
  const unsigned int N_XCOLUMNS_RCU  =   32; 
  const unsigned int N_ZROWS_MOD      =  56;                    /**<Number of rows per module*/       
  const unsigned int N_XCOLUMNS_MOD   =  64;                 /**<Number of columns per module*/ 
  const unsigned int N_GAINS         =   2;                             /**<Number of gains per ALTRO channel*/
  const unsigned int N_DATATYPES     =   10;    

  const unsigned int  PF_MAX_PATH_LENGTH = 256;

#ifndef __CINT__
  const unsigned char PF_VECTOR_DIR[] = "/HLT/PHOS/PFVectors";
#endif

    const unsigned int PF_DEFAULT_N_SAMPLES = 70;
  const unsigned int PF_DEFAULT_STARTINDEX = 0;

  const unsigned int DEFAULT_TAU = 2;    /**<Assume that the signal rise time of the altrp pulses is 2 us (nominal value of the electronics)*/
  const unsigned int  DEFAULT_FS = 10;   /**<Assume that the signal is samples with 10 MHZ samle rate*/

  const unsigned int MODULE_0     = 0;
  const unsigned int MODULE_1     = 1;
  const unsigned int MODULE_2     = 2;
  const unsigned int MODULE_3     = 3;
  const unsigned int MODULE_4     = 4;

  const unsigned int CSPS_PER_FEE    = 32;
  const unsigned int RCU_0       = 0;
  const unsigned int RCU_1       = 1;
  const unsigned int RCU_2       = 2;
  const unsigned int RCU_3       = 3;

  const unsigned int Z_0         = 0;
  const unsigned int Z_1         = 1;
  const unsigned int X_0         = 0;
  const unsigned int X_1         = 1;

  const unsigned int N_MODULES    =      5;                             /**<Number of modules of the PHOS detector*/
  const unsigned int N_RCUS       =      4;                             /**<Number of RCUs per Module*/
  const unsigned int N_RCUS_PER_MODULE =  4 ;                            /**<Number of RCUs per Module*/
  const unsigned int N_RCUS_PER_TOTAL =  N_MODULES*N_RCUS_PER_MODULE;   /**<Total number of RCUs for PHOS*/
  const unsigned int N_FEECS         =  14;                             /**<Number of Frontend cards per branch*/
  const unsigned int N_ALTROS        =   4;                             /**<Number of ALTROs per frontend card*/
  const unsigned int N_ALTROCHANNELS =  16;
  const unsigned int N_BRANCHES      =   2;      
}


#endif


