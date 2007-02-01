#ifndef ALIPHOSCOMMONDEFS_H
#define ALIPHOSCOMMONDEFS_H

//Hardware constants
#define N_RCUS             4   /**<Number of RCUs per Module*/
#define N_BRANCHES         2   /**<Number of branches per RCU*/
#define N_FEECS           14   /**<Number of Frontend cards per branch*/
#define N_ALTROS           4   /**<Number of ALTROs per frontend card*/
#define N_ALTROCHANNELS   16   /**<Number of readout channles per ALTRO*/

//Geometry constants
#define N_MODULES          5   /**<Number of modules of the PHOS detector*/
#define N_ROWS_MOD         64  /**<Number of rows per module*/       
#define N_COLUMNS_MOD      56  /**<Number of columns per module*/ 
#define N_ROWS_RCU         32  /**<Number of rows per module*/       
#define N_COLUMNS_RCU      28  /**<Number of columns per module*/ 
#define N_GAINS            2   /**<Number of gains per ALTRO channel*/

#endif
