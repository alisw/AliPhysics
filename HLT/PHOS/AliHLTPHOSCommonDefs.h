#ifndef ALIPHOSCOMMONDEFS_H
#define ALIPHOSCOMMONDEFS_H

//Hardware constants
#define PHOS_RCUS             4   /**<Number of RCUs per Module*/
#define PHOS_BRANCHES         2   /**<Number of branches per RCU*/
#define PHOS_FEECS           14   /**<Number of Frontend cards per branch*/
#define PHOS_ALTROS           4   /**<Number of ALTROs per frontend card*/
#define PHOS_ALTROCHANNELS   16   /**<Number of readout channles per ALTRO*/

//Geometry constants
#define PHOS_MODULES          5   /**<Number of modules of the PHOS detector*/
#define PHOS_ROWS            64   /**<Number of rows per module*/       
#define PHOS_COLUMNS         56   /**<Number of columns per module*/ 
#define PHOS_GAINS            2   /**<Number of gains per ALTRO channel*/

#endif
