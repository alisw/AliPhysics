#ifndef ALIPHOSCOMMONDEFS_H
#define ALIPHOSCOMMONDEFS_H

//Hardware constants
#define N_MODULES          5                             /**<Number of modules of the PHOS detector*/
#define N_RCUS             4                             /**<Number of RCUs per Module*/
#define N_RCUS_PER_MODULE  4                             /**<Number of RCUs per Module*/
#define N_RCUS_PER_TOTAL   N_MODULES*N_RCUS_PER_MODULE   /**<Total number of RCUs for PHOS*/
#define N_BRANCHES         2                             /**<Number of branches per RCU*/
#define N_FEECS           14                             /**<Number of Frontend cards per branch*/
#define N_ALTROS           4                             /**<Number of ALTROs per frontend card*/
#define N_ALTROCHANNELS   16                             /**<Number of readout channles per ALTRO*/
#define ALTRO_MAX_SAMPLES 1008                           /**<The maximum number of samples of the ALTRO*/

//Geometry constants
//#define N_MODULES          5                           /**<Number of modules of the PHOS detector*/

//#define N_ROWS_MOD         64                          /**<Number of rows per module*/       
//#define N_COLUMNS_MOD      56                          /**<Number of columns per module*/ 

#define N_ROWS_MOD         56                            /**<Number of rows per module*/       
#define N_COLUMNS_MOD      64                            /**<Number of columns per module*/ 

#define N_ROWS_RCU         28                            /**<Number of rows per module*/       
#define N_COLUMNS_RCU      32   
#define N_ZROWS_RCU        N_ROWS_RCU                                 /**<Number of rows per module*/       
#define N_XCOLUMNS_RCU     N_COLUMNS_RCU 

#define N_ZROWS_MOD        N_ROWS_MOD                            /**<Number of rows per module*/       
#define N_XCOLUMNS_MOD     N_COLUMNS_MOD                            /**<Number of columns per module*/ 

//#define N_ROWS_RCU         32                            /**<Number of rows per module*/       
//#define N_COLUMNS_RCU      28                            /**<Number of columns per module*/ 



//#define N_ZROWS_RCU        N_ROWS_RCU                                 /**<Number of rows per module*/       
//#define N_XCOLUMNS_RCU     N_COLUMNS_RCU 

#define N_GAINS            2                             /**<Number of gains per ALTRO channel*/

#define N_DATATYPES        10 

//peakfinder constatnts
#define PF_MAX_PATH_LENGTH 256
#define PF_VECTOR_DIR "/HLT/PHOS/PFVectors"
#define PF_DEFAULT_N_SAMPLES 70
#define PF_DEFAULT_STARTINDEX 0

//analysis constatnts
#define LOW_GAIN 1
#define HIGH_GAIN 0


//general altro signal constatnts
#define DEFAULT_TAU 2    /**<Assume that the signal rise time of the altrp pulses is 2 us (nominal value of the electronics)*/
#define DEFAULT_FS  10   /**<Assume that the signal is samples with 10 MHZ samle rate*/

#endif
