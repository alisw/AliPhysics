/********************************************************************************
 *
 * equipment.h: DATE equipment data format
 *
 * Revision History:
 *    V01.1  AV        12-Aug-98 Changed order
 *    V01.0  MG & AV   09-Mar-98 Initial version
 *
 *******************************************************************************/

//We use this copy of equipment.h file
//to avoid nececcity compilation of aliroot with DATE 

#ifndef __AliPHOSequipment_h__
#define __AliPHOSequipment_h__

struct equipmentHeaderStruct {
  short headerExtLen;   /* length in bytes of an optional extension */
  short type;           /* the equipment type identifier */
  char reserved;        /* reserved byte */
  char rawByteAlign;    /* length (in bytes) of the word read from hw */
  short equipmentId;    /* equipment identifier */
  long rawDataLen;      /* length (in bytes) of the data block (no header) */
};

struct AliPHOSequipmentStruct {
  struct equipmentHeaderStruct equipmentHeader;  /* the equipment header */
  unsigned short rawData[1];     /* the equipment raw data */
};

#endif
