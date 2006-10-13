/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
Revision 1.8  2006/08/22 13:30:17  arcelli
removal of effective c++ warnings (C.Zampolli)

Revision 1.7  2006/08/10 14:46:54  decaro
TOF raw data format: updated version

Revision 1.6.1  2006/06/28 A. De Caro, R. Preghenella:
        Update TOF raw data format
        according to the final version
        (see the ALICE internal note in preparation
         'ALICE TOF raw data format')
        Added the methods for the correspoonding numbering
         between the equipment IDs and the volume IDs:
           Equip2VolNPlate(...)
           Equip2VolNStrip(...)
           Equip2VolNPad(...)

Revision 0.02  2005/07/28 A. De Caro:
        Update format TOF raw data
               (temporary solution) 
        Correction of few wrong corrispondences
               between 'software' and 'hardware' numberings

Revision 0.01  2005/07/22 A. De Caro
        Implement methods Next()
	                  GetSector(),
	                  GetPlate(),
	                  GetStrip(),
	                  GetPadZ(),
	                  GetPadX()
*/

////////////////////////////////////////////////////////////////////////
//                                                                    //
//     This class provides access to TOF raw data in DDL files.       //
//                                                                    //
//      It loops over all TOF raw data given by the AliRawReader.     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliRawReader.h"

#include "AliTOFGeometry.h"
#include "AliTOFGeometryV5.h"
#include "AliTOFRawStream.h"


/******************************************
GENERAL DATA FORMAT                    
******************************************/

//filler
#define FILLER 0x70000000

//word type mask/position
#define WORD_TYPE_MASK 0xf0000000
#define WORD_TYPE_POSITION 28

//global header word required bit pattern
#define GLOBAL_HEADER 0x40000000

//global trailer word required bit pattern
#define GLOBAL_TRAILER 0x50000000

//error word required bit pattern
#define ERROR 0x30000000

//header slot ID mask/position
#define HEADER_SLOT_ID_MASK 0x0000000f
#define HEADER_SLOT_ID_POSITION 0

//word types
#define GLOBAL_HEADER_TYPE 4
#define GLOBAL_TRAILER_TYPE 5
#define ERROR_TYPE 6
#define FILLER_TYPE 7
#define TRM_CHAIN0_HEADER_TYPE 0
#define TRM_CHAIN0_TRAILER_TYPE 1
#define TRM_CHAIN1_HEADER_TYPE 2
#define TRM_CHAIN1_TRAILER_TYPE 3

//slot types
#define DRM_ID_NUMBER 1
#define LTM_ID_NUMBER 2


/******************************************
DRM DATA FORMAT                        
******************************************/

//DRM global header word required bit pattern
#define DRM_GLOBAL_HEADER 0x40000001

//DRM event words mask/position
#define DRM_EVENT_WORDS_MASK 0x001ffff0
#define DRM_EVENT_WORDS_POSITION 4

//DRM DRM ID mask/position
#define DRM_DRM_ID_MASK 0x0fe00000
#define DRM_DRM_ID_POSITION 21

//DRM status header 1 word required bit pattern
#define DRM_STATUS_HEADER_1 0x40000001

//DRM slot ID mask/position
#define DRM_SLOT_ID_MASK 0x00007ff0
#define DRM_SLOT_ID_POSITION 4

//DRM C-bit mask/position
#define DRM_C_BIT_MASK 0x00008000
#define DRM_C_BIT_POSITION 15

//DRM status header 2 word required bit pattern
#define DRM_STATUS_HEADER_2 0x40000001

//DRM enable ID mask/position
#define DRM_ENABLE_ID_MASK 0x00007ff0
#define DRM_ENABLE_ID_POSITION 4

//DRM fault ID mask/position
#define DRM_FAULT_ID_MASK 0x07ff0000
#define DRM_FAULT_ID_POSITION 16

//DRM status header 3 word required bit pattern
#define DRM_STATUS_HEADER_3 0x40000001

//DRM TTC event counter mask/position
#define DRM_TTC_EVENT_COUNTER_MASK 0x0ffffff0
#define DRM_TTC_EVENT_COUNTER_POSITION 4

//DRM event CRC mask/position
//#define DRM_EVENT_CRC_MASK 0x001ffff0
#define DRM_EVENT_CRC_MASK 0x000ffff0
#define DRM_EVENT_CRC_POSITION 4

//DRM global trailer word required bit pattern
#define DRM_GLOBAL_TRAILER 0x50000001

//DRM local event counter mask/position
#define DRM_LOCAL_EVENT_COUNTER_MASK 0x0000fff0
#define DRM_LOCAL_EVENT_COUNTER_POSITION 4


/******************************************
TRM DATA FORMAT                        
******************************************/

//TRM global header word required bit pattern
#define TRM_GLOBAL_HEADER 0x40000000

//TRM slot ID mask/position
#define TRM_SLOT_ID_MASK 0x0000000f
#define TRM_SLOT_ID_POSITION 0

//TRM event words mask/position
#define TRM_EVENT_WORDS_MASK 0x0001fff0
#define TRM_EVENT_WORDS_POSITION 4

//TRM ACQ-bits mask/position
#define TRM_ACQ_BITS_MASK 0x00060000
#define TRM_ACQ_BITS_POSITION 17

//TRM L-bit mask/position
#define TRM_L_BIT_MASK 0x00080000
#define TRM_L_BIT_POSITION 19

//TRM chain-0 header word required bit pattern
#define TRM_CHAIN_0_HEADER 0x00000000

//TRM chain-1 header word required bit pattern
#define TRM_CHAIN_1_HEADER 0x20000000

//TRM bunch ID mask/position
#define TRM_BUNCH_ID_MASK 0x0000fff0
#define TRM_BUNCH_ID_POSITION 4

//TRM PB24 temp mask/position
#define TRM_PB24_TEMP_MASK 0x00ff0000
#define TRM_PB24_TEMP_POSITION 16

//TRM PB24 ID mask/position
#define TRM_PB24_ID_MASK 0x07000000
#define TRM_PB24_ID_POSITION 24

//TRM TS-bit mask/position
#define TRM_TS_BIT_MASK 0x08000000
#define TRM_TS_BIT_POSITION 27

//TRM chain-0 trailer word required bit pattern
#define TRM_CHAIN_0_TRAILER 0x10000000

//TRM chain-1 trailer word required bit pattern
#define TRM_CHAIN_1_TRAILER 0x30000000

//TRM status mask/position
#define TRM_STATUS_MASK 0x0000000f
#define TRM_STATUS_POSITION 0


//TDC digit

//TRM TDC digit word required bit pattern
#define TRM_TDC_DIGIT 0x8000000

//TRM digit time mask/position
#define TRM_DIGIT_TIME_MASK 0x00001fff
#define TRM_DIGIT_TIME_POSITION 0

//TRM long digit time mask/position
#define TRM_LONG_DIGIT_TIME_MASK 0x001fffff
#define TRM_LONG_DIGIT_TIME_POSITION 0

//TRM TOT width mask/position
#define TRM_TOT_WIDTH_MASK 0x001fe000
#define TRM_TOT_WIDTH_POSITION 13

//TRM chan mask/position
#define TRM_CHAN_MASK 0x00e00000
#define TRM_CHAN_POSITION 21

//TRM TDC ID mask/position
#define TRM_TDC_ID_MASK 0x0f000000
#define TRM_TDC_ID_POSITION 24

//TRM E-bit mask/position
#define TRM_E_BIT_MASK 0x10000000
#define TRM_E_BIT_POSITION 28

//TRM PS-bits mask/position
#define TRM_PS_BITS_MASK 0x60000000
#define TRM_PS_BITS_POSITION 29


//TRM errors

//TRM TDC error word required bit pattern
#define TRM_TDC_ERROR 0x6000000

//TRM TDC diagnostic error word required bit pattern
#define TRM_TDC_DIAGNOSTIC_ERROR 0x6f00000

//TRM TDC error flags mask/position
#define TRM_TDC_ERROR_FLAGS_MASK 0x00007fff
#define TRM_TDC_ERROR_FLAGS_POSITION 0

//TRM TDC error TDC ID mask/position
#define TRM_TDC_ERROR_TDC_ID_MASK 0x0f00000
#define TRM_TDC_ERROR_TDC_ID_POSITION 24

//TRM TDC fault chip flag ID mask/position
#define TRM_TDC_ERROR_FAULT_CHIP_FLAG_ID_MASK 0x00007fff
#define TRM_TDC_ERROR_FAULT_CHIP_FLAG_ID_POSITION 0

//TRM TDC error C-bit mask/position
#define TRM_TDC_ERROR_C_BIT_MASK 0x00008000
#define TRM_TDC_ERROR_C_BIT_POSITION 15

//TRM TDC JTAG error code mask/position
#define TRM_TDC_ERROR_JTAG_ERROR_CODE_MASK 0x000007ff
#define TRM_TDC_ERROR_JTAG_ERROR_CODE_POSITION 0

//TRM TDC disgnostic error TDC ID mask/position
#define TRM_TDC_DIAGNOSTIC_ERROR_TDC_ID_MASK 0x00007800
#define TRM_TDC_DIAGNOSTIC_ERROR_TDC_ID_POSITION 11 

//TRM global trailer word required bit pattern
//#define TRM_GLOBAL_TRAILER 0x50000000
#define TRM_GLOBAL_TRAILER 0x5000000f

//TRM event CRC mask/position
#define TRM_EVENT_CRC_MASK 0x0000fff0
#define TRM_EVENT_CRC_POSITION 4

//TRM event counter mask/position
#define TRM_EVENT_COUNTER_MASK 0x0fff0000
#define TRM_EVENT_COUNTER_POSITION 16


/******************************************
LTM DATA FORMAT                        
******************************************/

//LTM global header word required bit pattern
#define LTM_GLOBAL_HEADER 0x40000002

//LTM event words mask/position
#define LTM_EVENT_WORDS_MASK 0x0001fff0
#define LTM_EVENT_WORDS_POSITION 4

//LTM C-bit mask/position
#define LTM_C_BIT_MASK 0x00020000
#define LTM_C_BIT_POSITION 17

//LTM fault mask/position
#define LTM_FAULT_MASK 0x00fc0000
#define LTM_FAULT_POSITION 18

//PDL data 

//PDL value 1 mask/position
#define LTM_PDL_VALUE_1_MASK 0x000000ff
#define LTM_PDL_VALUE_1_POSITION 0

//PDL value 2 mask/position
#define LTM_PDL_VALUE_2_MASK 0x0000ff00
#define LTM_PDL_VALUE_2_POSITION 8

//PDL value 3 mask/position
#define LTM_PDL_VALUE_3_MASK 0x00ff0000
#define LTM_PDL_VALUE_3_POSITION 16

//PDL value 4 mask/position
#define LTM_PDL_VALUE_4_MASK 0xff000000
#define LTM_PDL_VALUE_4_POSITION 24

//ADC data 

//ADC value 1 mask/position
#define LTM_ADC_VALUE_1_MASK 0x000003ff
#define LTM_ADC_VALUE_1_POSITION 0

//ADC value 2 mask/position
#define LTM_ADC_VALUE_2_MASK 0x000ffc00
#define LTM_ADC_VALUE_2_POSITION 10

//ADC value 3 mask/position
#define LTM_ADC_VALUE_3_MASK 0x3ff00000
#define LTM_ADC_VALUE_3_POSITION 20

//LTM global trailer word required bit pattern
#define LTM_GLOBAL_TRAILER 0x50000002

//LTM event CRC mask/position
#define LTM_EVENT_CRC_MASK 0x0000fff0
#define LTM_EVENT_CRC_POSITION 4

//LTM event number mask/position
#define LTM_EVENT_NUMBER_MASK 0x0fff0000
#define LTM_EVENT_NUMBER_POSITION 16


ClassImp(AliTOFRawStream)


//_____________________________________________________________________________
AliTOFRawStream::AliTOFRawStream(AliRawReader* rawReader):
  fRawReader(0x0),
  fDDL(-1),
  fTRM(-1),
  fTDC(-1),
  fTRMchain(-1),
  fTDCchannel(-1),
  fTof(-1),
  fToT(-1),
  fErrorFlag(-1),

  fSector(-1),
  fPlate(-1),
  fStrip(-1),
  fPadX(-1),
  fPadZ(-1),
  fTOFGeometry(new AliTOFGeometryV5()),
  fWordType(-1),
  fSlotID(-1),
  fACQ(-1),
  fPSbit(-1),
  fTime(-1),
  fTDCerrorFlag(-1),
  fInsideDRM(kFALSE),
  fInsideTRM(kFALSE),
  fInsideLTM(kFALSE),
  fInsideTRMchain0(kFALSE),
  fInsideTRMchain1(kFALSE),
  fLeadingOrphane(kFALSE)
{
  //
  // create an object to read TOF raw digits
  //

  fRawReader = rawReader;
  /*
  fDDL = -1;
  fTRM = -1;
  fTDC = -1;
  fTRMchain = -1;
  fTDCchannel = -1;
  fTof = -1;
  fToT = -1;
  fErrorFlag = -1;

  fSector = -1;
  fPlate = -1;
  fStrip = -1;
  fPadX = -1;
  fPadZ = -1;
  fTOFGeometry = new AliTOFGeometryV5();
  */
  fRawReader->Select("TOF");
  /*
  fWordType = -1;
  fSlotID = -1;
  fACQ = -1;
  fPSbit = -1;
  fTime = -1;
  fTDCerrorFlag = -1;
  fInsideDRM = kFALSE;
  fInsideTRM = kFALSE;
  fInsideLTM = kFALSE;
  fInsideTRMchain0 = kFALSE;
  fInsideTRMchain1 = kFALSE;
  fLeadingOrphane = kFALSE;
  */
}

//_____________________________________________________________________________
AliTOFRawStream::AliTOFRawStream():
  fRawReader(0x0),
  fDDL(-1),
  fTRM(-1),
  fTDC(-1),
  fTRMchain(-1),
  fTDCchannel(-1),
  fTof(-1),
  fToT(-1),
  fErrorFlag(-1),

  fSector(-1),
  fPlate(-1),
  fStrip(-1),
  fPadX(-1),
  fPadZ(-1),
  fTOFGeometry(new AliTOFGeometryV5()),
  fWordType(-1),
  fSlotID(-1),
  fACQ(-1),
  fPSbit(-1),
  fTime(-1),
  fTDCerrorFlag(-1),
  fInsideDRM(kFALSE),
  fInsideTRM(kFALSE),
  fInsideLTM(kFALSE),
  fInsideTRMchain0(kFALSE),
  fInsideTRMchain1(kFALSE),
  fLeadingOrphane(kFALSE)
{
  //
  // default ctr
  //
  /*
  fRawReader = 0x0;
  fDDL = -1;
  fTRM = -1;
  fTDC = -1;
  fTRMchain = -1;
  fTDCchannel = -1;
  fTof = -1;
  fToT = -1;
  fErrorFlag = -1;

  fSector = -1;
  fPlate = -1;
  fStrip = -1;
  fPadX = -1;
  fPadZ = -1;
  fTOFGeometry = new AliTOFGeometryV5();
  fWordType = -1;
  fSlotID = -1;
  fACQ = -1;
  fPSbit = -1;
  fTime = -1;
  fTDCerrorFlag = -1;
  fInsideDRM = kFALSE;
  fInsideTRM = kFALSE;
  fInsideLTM = kFALSE;
  fInsideTRMchain0 = kFALSE;
  fInsideTRMchain1 = kFALSE;
  fLeadingOrphane = kFALSE;
  */
}

//_____________________________________________________________________________
AliTOFRawStream::AliTOFRawStream(const AliTOFRawStream& stream) :
  TObject(stream),
  fRawReader(0x0),
  fDDL(-1),
  fTRM(-1),
  fTDC(-1),
  fTRMchain(-1),
  fTDCchannel(-1),
  fTof(-1),
  fToT(-1),
  fErrorFlag(-1),

  fSector(-1),
  fPlate(-1),
  fStrip(-1),
  fPadX(-1),
  fPadZ(-1),
  fTOFGeometry(new AliTOFGeometryV5()),
  fWordType(-1),
  fSlotID(-1),
  fACQ(-1),
  fPSbit(-1),
  fTime(-1),
  fTDCerrorFlag(-1),
  fInsideDRM(kFALSE),
  fInsideTRM(kFALSE),
  fInsideLTM(kFALSE),
  fInsideTRMchain0(kFALSE),
  fInsideTRMchain1(kFALSE),
  fLeadingOrphane(kFALSE)
{
  //
  // copy constructor
  //

  fRawReader = stream.fRawReader;
  fDDL = stream.fDDL;
  fTRM = stream.fTRM;
  fTDC = stream.fTDC;
  fTRMchain = stream.fTRMchain;
  fTDCchannel = stream.fTDCchannel;
  fTof = stream.fTof;
  fToT = stream.fToT;
  fErrorFlag = stream.fErrorFlag;

  fSector = stream.fSector;
  fPlate = stream.fPlate;
  fStrip = stream.fStrip;
  fPadX = stream.fPadX;
  fPadZ = stream.fPadZ;

  fTOFGeometry = stream.fTOFGeometry;

  fWordType = stream.fWordType;
  fSlotID = stream.fSlotID;
  fACQ = stream.fACQ;
  fPSbit = stream.fPSbit;
  fTime = stream.fTime;
  fTDCerrorFlag = stream.fTDCerrorFlag;
  fInsideDRM = stream.fInsideDRM;
  fInsideTRM = stream.fInsideTRM;
  fInsideLTM = stream.fInsideLTM;
  fInsideTRMchain0 = stream.fInsideTRMchain0;
  fInsideTRMchain1 = stream.fInsideTRMchain1;
  fLeadingOrphane = stream.fLeadingOrphane;

}

//_____________________________________________________________________________
AliTOFRawStream& AliTOFRawStream::operator = (const AliTOFRawStream& stream)
{
  //
  // assignment operator
  //

  fRawReader = stream.fRawReader;
  fDDL = stream.fDDL;
  fTRM = stream.fTRM;
  fTDC = stream.fTDC;
  fTRMchain = stream.fTRMchain;
  fTDCchannel = stream.fTDCchannel;
  fTof = stream.fTof;
  fToT = stream.fToT;
  fErrorFlag = stream.fErrorFlag;

  fSector = stream.fSector;
  fPlate = stream.fPlate;
  fStrip = stream.fStrip;
  fPadX = stream.fPadX;
  fPadZ = stream.fPadZ;

  fTOFGeometry = stream.fTOFGeometry;

  fWordType = stream.fWordType;
  fSlotID = stream.fSlotID;
  fACQ = stream.fACQ;
  fPSbit = stream.fPSbit;
  fTime = stream.fTime;
  fTDCerrorFlag = stream.fTDCerrorFlag;
  fInsideDRM = stream.fInsideDRM;
  fInsideTRM = stream.fInsideTRM;
  fInsideLTM = stream.fInsideLTM;
  fInsideTRMchain0 = stream.fInsideTRMchain0;
  fInsideTRMchain1 = stream.fInsideTRMchain1;
  fLeadingOrphane = stream.fLeadingOrphane;

  return *this;

}

//_____________________________________________________________________________
AliTOFRawStream::~AliTOFRawStream()
{
// destructor

  fTOFGeometry = 0;

}


//_____________________________________________________________________________
Bool_t AliTOFRawStream::Next()
{
  //
  // Read next 32-bit word in TOF raw data files
  // returns kFALSE if there is no word left
  //

  UInt_t data;

  if (!fRawReader->ReadNextInt(data)) return kFALSE;

  if (fSector!=-1 && fPlate!=-1 && fStrip!=-1 && fPadZ!=-1 && fPadX!=-1) {
    fSector = -1;
    fPlate = -1;
    fStrip = -1;
    fPadZ = -1;
    fPadX = -1;
  }


  fDDL  = fRawReader->GetDDLID();

  // orphane digits
  AliTOFtdcDigit orphaneLeadingDigit={0,0,0,0,0,0,0};

  fWordType = GetField(data,WORD_TYPE_MASK,WORD_TYPE_POSITION);

  switch (fWordType) { // switch word type

  case GLOBAL_HEADER_TYPE: // global header
    fSlotID = GetField(data, HEADER_SLOT_ID_MASK, HEADER_SLOT_ID_POSITION);
    fTRM = fSlotID;


    switch (fSlotID) { // switch global header slot ID

    case DRM_ID_NUMBER: //DRM global header
      if (fInsideDRM) { // unexpected DRM global headers -> exit
	break;
      }
      fInsideDRM = kTRUE; // DRM global header accepted
      break;

    case LTM_ID_NUMBER: // LTM global header
      if (fInsideLTM) { // unexpected LTM global headers -> exit
	break;
      }
      fInsideLTM = kTRUE; // LTM global header accepted
      break;

    case  3: //TRM header
    case  4: //TRM header
    case  5: //TRM header
    case  6: //TRM header
    case  7: //TRM header
    case  8: //TRM header
    case  9: //TRM header
    case 10: //TRM header
    case 11: //TRM header
    case 12: //TRM header
      if (fInsideTRM) { // unexpected TRM global headers -> exit
	break;
      }
      fInsideTRM = kTRUE; // TRM global header accepted
      fACQ =  GetField(data,TRM_ACQ_BITS_MASK,TRM_ACQ_BITS_POSITION);
      break;

    default: // unexpected global header slot ID
      break;

    } //end switch global header slot id

    break;


  case GLOBAL_TRAILER_TYPE: // global trailer
    fSlotID = GetField(data,HEADER_SLOT_ID_MASK,HEADER_SLOT_ID_POSITION);
    

    switch (fSlotID) { // switch global trailer slot ID

    case DRM_ID_NUMBER: // DRM global trailer
      if (!fInsideDRM) { // unexpected DRM global trailers -> exit
	break;
      }
      fInsideDRM = kFALSE; // DRM global trailer accepted
      fInsideTRM = kFALSE;
      fInsideLTM = kFALSE;
      fInsideTRMchain0 = kFALSE;
      fInsideTRMchain1 = kFALSE;
      fLeadingOrphane = kFALSE;
      fSector = -1;
      fPlate = -1;
      fStrip = -1;
      fPadZ = -1;
      fPadX = -1;
      fDDL = -1;
      fTRM = -1;
      fTDC = -1;
      fTRMchain = -1;
      fTDCchannel = -1;
      fTof = -1;
      fToT = -1;
      fErrorFlag = -1;
      fACQ = -1;
      fPSbit = -1;
      fTime = -1;
      fTDCerrorFlag = -1;
      break;
    case LTM_ID_NUMBER: // LTM global trailer
      if (!fInsideLTM) { // unexpected LTM global trailer -> exit
	break;
      }
      fInsideLTM = kFALSE; // LTM global trailer accepted
      break;
    case 15: //TRM global trailer
      if (!fInsideTRM) { // unexpected TRM global trailers -> exit
	break;
      }
      fInsideTRM = kFALSE; // TRM global trailer accepted
      break;
    default: // unexpected global trailer slot ID
      break;
    } //end switch global trailer slot id


    break;


  case ERROR_TYPE: // TDC error
    fTDC = GetField(data,TRM_TDC_ERROR_TDC_ID_MASK,TRM_TDC_ERROR_TDC_ID_POSITION);
    fTDCerrorFlag = GetField(data,TRM_TDC_ERROR_FLAGS_MASK,TRM_TDC_ERROR_FLAGS_POSITION);
    break;


  case FILLER_TYPE: // filler
    break;


  default: // other word types

    if (fInsideTRM) { // inside TRM

      switch (fWordType) { // switch word type inside TRM
      case TRM_CHAIN0_HEADER_TYPE: // TRM chain0 header
	if (fInsideTRMchain0) { // unexpected TRM chain0 header
	  break;
	}
	fInsideTRMchain0 = kTRUE;
	fTRMchain = 0;
	break;
      case TRM_CHAIN0_TRAILER_TYPE: // TRM chain0 trailer
	if (!fInsideTRMchain0) { // unexpected TRM chain0 trailer
	  break;
	}
	fInsideTRMchain0 = kFALSE;
	fTRMchain = -1;
	break;
      case TRM_CHAIN1_HEADER_TYPE: // TRM chain1 header
	if (fInsideTRMchain1) { // unexpected TRM chain1 header
	  break;
	}
	fInsideTRMchain1 = kTRUE;
	fTRMchain = 1;
	break;
      case TRM_CHAIN1_TRAILER_TYPE: // TRM chain1 trailer
	if (!fInsideTRMchain1) { // unexpected TRM chain1 trailer
	  break;
	}
	fInsideTRMchain1 = kFALSE;
	fTRMchain = -1;
	break;
      } // end switch word type inside TRM

    } // end if (fInsideTRM)

      
    if (
	((fInsideTRMchain0&&!fInsideTRMchain1) || (!fInsideTRMchain0&&fInsideTRMchain1)) 
	&& fWordType!=TRM_CHAIN0_HEADER_TYPE && fWordType!=TRM_CHAIN0_TRAILER_TYPE
	&& fWordType!=TRM_CHAIN1_HEADER_TYPE && fWordType!=TRM_CHAIN1_TRAILER_TYPE
	){ // inside TRM chains
      fPSbit = GetField(data,TRM_PS_BITS_MASK,TRM_PS_BITS_POSITION);
      fTDC = GetField(data,TRM_TDC_ID_MASK,TRM_TDC_ID_POSITION);
      fTDCchannel = GetField(data,TRM_CHAN_MASK,TRM_CHAN_POSITION);

      switch (fPSbit) { // switch fPSbit bits inside TRM chains
      case 0: // packing ok, digit time and tot
	fToT = GetField(data,TRM_TOT_WIDTH_MASK,TRM_TOT_WIDTH_POSITION);
	fTime = GetField(data,TRM_DIGIT_TIME_MASK,TRM_DIGIT_TIME_POSITION);
	fTof = fTime;
	SetSector();
	SetPlate();
	SetStrip();
	SetPadZ();
	SetPadX();
	break;

      case 1: // leading edge digit, long digit time, no TOT
	fToT = -1;
	fTime = GetField(data,TRM_LONG_DIGIT_TIME_MASK,TRM_LONG_DIGIT_TIME_POSITION);
	fTof = fTime;
	SetSector();
	SetPlate();
	SetStrip();
	SetPadZ();
	SetPadX();
	// always set it as orphane leading
	fLeadingOrphane=1;
	orphaneLeadingDigit.fSlotID = fSlotID;
	orphaneLeadingDigit.fChain = fTRMchain;
	orphaneLeadingDigit.fPS = fPSbit;
	orphaneLeadingDigit.fTDC = fTDC;
	orphaneLeadingDigit.fChannel = fTDCchannel;
	orphaneLeadingDigit.fTOT = fToT;
	orphaneLeadingDigit.fTime = fTime;
	break;

      case 2: // trailing edge digit, long digit time, no TOT
	fToT = -1;
	fTime = GetField(data,TRM_LONG_DIGIT_TIME_MASK,TRM_LONG_DIGIT_TIME_POSITION);
	fTof = fTime;
	SetSector();
	SetPlate();
	SetStrip();
	SetPadZ();
	SetPadX();
	if (fACQ!=3) // check if packing is disabled
	  break;
	if (!fLeadingOrphane) // check for a orphane leading edge
	  break;
	if (orphaneLeadingDigit.fSlotID != fSlotID ||
	    orphaneLeadingDigit.fChain != fTRMchain ||
	    orphaneLeadingDigit.fTDC != fTDC ||
	    orphaneLeadingDigit.fChannel != fTDCchannel) // check leading edge compatibility (fSlotID, fTRMchain, fTDC, fTDCchannel)
	  break;
	fLeadingOrphane = 0; // orphane leading is no longer orphane
	SetSector();
	SetPlate();
	SetStrip();
	SetPadZ();
	SetPadX();
	break;
      case 3: // TOT overflow
	fToT = GetField(data,TRM_TOT_WIDTH_MASK,TRM_TOT_WIDTH_POSITION);
	fTime = GetField(data,TRM_DIGIT_TIME_MASK,TRM_DIGIT_TIME_POSITION);
	fTof = fTime;
	SetSector();
	SetPlate();
	SetStrip();
	SetPadZ();
	SetPadX();
	break;
      } // end switch fPSbit bits inside TRM chains


    } // end if is inside TRM chains

  } // end switch on fWordType

  return kTRUE;

}
//_____________________________________________________________________________

void AliTOFRawStream::SetSector()
{
  //
  // Evaluate the TOF sector number -> [ 0;17]
  // corresponding to the TOF equipment IDs:
  //                                  fDDL        -> [ 0;71]
  //                                  fTRM        -> [ 3;12]
  //                                  fTRMchain   -> [0;  1]
  //                                  fTDC        -> [ 0;14]
  //                                  fTDCchannel -> [ 0; 7]
  //

  Int_t iSector = -1;

  if (!(fDDL==-1)) iSector = Int_t((Float_t)(fDDL)/AliTOFGeometry::NDDL());

  fSector = iSector;

}
//_____________________________________________________________________________


void AliTOFRawStream::SetPlate()
{
  //
  // Evaluate the TOF plate number ->[ 0; 4]
  // corresponding to the TOF equipment IDs:
  //                                  fDDL        -> [ 0;71]
  //                                  fTRM        -> [ 3;12]
  //                                  fTRMchain   -> [0;  1]
  //                                  fTDC        -> [ 0;14]
  //                                  fTDCchannel -> [ 0; 7]
  //

  Int_t iPlate = -1;
  if (!(fDDL==-1 || fTRM==-1 || fTDC==-1
	|| fSector==-1))
    iPlate = Equip2VolNplate(GetDDLnumberPerSector(fDDL), fTRM, fTDC);

  fPlate = iPlate;

}
//_____________________________________________________________________________

void AliTOFRawStream::SetStrip()
{
  //
  // Evaluate the TOF strip number per module -> [ 0; 14/18]
  // corresponding to the TOF equipment IDs:
  //                                  fDDL        -> [ 0;71]
  //                                  fTRM        -> [ 3;12]
  //                                  fTRMchain   -> [0;  1]
  //                                  fTDC        -> [ 0;14]
  //                                  fTDCchannel -> [ 0; 7]
  //

  Int_t iStrip = -1;

  if (!(fDDL==-1 || fTRM==-1 || fTDC==-1
	|| fSector==-1 || fPlate==-1))
    iStrip = Equip2VolNstrip(GetDDLnumberPerSector(fDDL), fTRM, fTDC);

  fStrip = iStrip;

}
//_____________________________________________________________________________

void AliTOFRawStream::SetPadZ()
{
  //
  // Evaluate the TOF padRow number per strip -> [ 0; 1]
  // corresponding to the TOF equipment IDs:
  //                                  fDDL        -> [ 0;71]
  //                                  fTRM        -> [ 3;12]
  //                                  fTRMchain   -> [0;  1]
  //                                  fTDC        -> [ 0;14]
  //                                  fTDCchannel -> [ 0; 7]
  //

  Int_t iPadZ = -1;

  if (!(fDDL==-1 || fTRM==-1 || fTRMchain==-1 || fTDC==-1 || fTDCchannel==-1
	|| fSector==-1 || fPlate==-1 || fStrip==-1))
    {
      Int_t iPadAlongTheStrip = Equip2VolNpad(GetDDLnumberPerSector(fDDL), fTRMchain, fTDC, fTDCchannel);
      if (iPadAlongTheStrip!=-1)
	iPadZ  = iPadAlongTheStrip%AliTOFGeometry::NpadZ();
    }

  fPadZ = iPadZ;

}
//_____________________________________________________________________________

void AliTOFRawStream::SetPadX()
{
  //
  // Evaluate the TOF pad number per strip padRow -> [ 0;47]
  // corresponding to the TOF equipment IDs:
  //                                  fDDL        -> [ 0;71]
  //                                  fTRM        -> [ 3;12]
  //                                  fTRMchain   -> [0;  1]
  //                                  fTDC        -> [ 0;14]
  //                                  fTDCchannel -> [ 0; 7]
  //

  Int_t iPadX = -1;

  if (!(fDDL==-1 || fTRM==-1 || fTRMchain==-1 || fTDC==-1 || fTDCchannel==-1
	|| fSector==-1 || fPlate==-1 || fStrip==-1))
    {
      Int_t iPadAlongTheStrip = Equip2VolNpad(GetDDLnumberPerSector(fDDL), fTRMchain, fTDC, fTDCchannel);
      if (iPadAlongTheStrip!=-1)
	iPadX  = (Int_t)(iPadAlongTheStrip/(Float_t(AliTOFGeometry::NpadZ())));
    }

  fPadX = iPadX;

}

//----------------------------------------------------------------------------
Int_t AliTOFRawStream::GetField(UInt_t word, Int_t fieldMask, Int_t fieldPosition) const
{
  // 
  // 
  // 

  return ((word & fieldMask) >> fieldPosition);
}

//----------------------------------------------------------------------------
Int_t AliTOFRawStream::Equip2VolNplate(Int_t iDDL, Int_t nTRM, Int_t nTDC) const 
{
  //
  // Returns the TOF plate number [0;4]
  // corresponding to the TOF equipment ID numbers:
  //                          iDDL -> DDL number per sector [0;3]
  //                          nTRM -> TRM number [3;12]
  //                          nTDC -> TDC number [0;14]
  //

  Int_t iPlate = -1;
  if (iDDL==0) {

    if (nTRM>=4 && nTRM<7) {
      iPlate = 0;
    } else if (nTRM==7) {
      if (nTDC<12) iPlate = 0;
      else iPlate = 1;
    } else if (nTRM>=8 && nTRM<11) {
      iPlate = 1;
    } else if (nTRM==11) {
      if (nTDC<9) iPlate = 1;
      else iPlate = 2;
    }else if (nTRM==12) {
      iPlate = 2;
    } 

  } else if (iDDL==1) {

    if (nTRM==3) {
      if (nTDC<3) iPlate = 0;
    } else if (nTRM>=4 && nTRM<7) {
      iPlate = 0;
    } else if (nTRM==7) {
      if (nTDC<6) iPlate = 1;
      else iPlate = 0;
    } else if (nTRM>=8 && nTRM<11) {
      iPlate = 1;
    } else if (nTRM==11) {
      if (nTDC<9) iPlate = 2;
      else iPlate = 1;
    } else if (nTRM==12) {
      iPlate = 2;
    } 

  } else if (iDDL==2) {

    if (nTRM>=4 && nTRM<7) {
      iPlate = 4;
    } else if (nTRM==7) {
      if (nTDC<12) iPlate = 4;
      else iPlate = 3;
    } else if (nTRM>=8 && nTRM<11) {
      iPlate = 3;
    } else if (nTRM==11) {
      if (nTDC<9) iPlate = 3;
      else iPlate = 2;
    }else if (nTRM==12) {
      iPlate = 2;
    } 

  }  else if (iDDL==3) {

    if (nTRM==3) {
      if (nTDC<3) iPlate = 4;
    } else if (nTRM>=4 && nTRM<7) {
      iPlate = 4;
    } else if (nTRM==7) {
      if (nTDC<6) iPlate = 3;
      else iPlate = 4;
    } else if (nTRM>=8 && nTRM<11) {
      iPlate = 3;
    } else if (nTRM==11) {
      if (nTDC<9) iPlate = 2;
      else iPlate = 3;
    } else if (nTRM==12) {
      iPlate = 2;
    } 

  }

  return iPlate;

}

//----------------------------------------------------------------------------
Int_t AliTOFRawStream::Equip2VolNstrip(Int_t iDDL, Int_t nTRM, Int_t nTDC) const 
{
  //
  // Returns the TOF strip number per module:
  //                                [0;14], in the central plates,
  //                                [0;18], in the intermediate and external plates
  // corresponding to the TOF equipment ID numbers:
  //                                iDDL -> DDL number per sector [0;3]
  //                                nTRM -> TRM number [3;12]
  //                                nTDC -> TDC number [0;14]
  //

  Int_t iStrip = -1;

  if (iDDL==0) {

    if (nTRM== 4) iStrip =  (Int_t)(nTDC/3.);
    else if (nTRM== 5) iStrip =  5 + (Int_t)(nTDC/3.);
    else if (nTRM== 6) iStrip = 10 + (Int_t)(nTDC/3.);
    else if (nTRM== 7) {
      if (nTDC<12) iStrip =  15 + (Int_t)(nTDC/3.);
      else iStrip = (Int_t)(nTDC/3.) -  4;
    }
    else if (nTRM== 8) iStrip =  1 + (Int_t)(nTDC/3.);
    else if (nTRM== 9) iStrip =  6 + (Int_t)(nTDC/3.);
    else if (nTRM==10) iStrip = 11 + (Int_t)(nTDC/3.);
    else if (nTRM==11) {
      if (nTDC<9) iStrip = 16 + (Int_t)(nTDC/3.);
      else iStrip = (Int_t)(nTDC/3.) -  3;
    }
    else if (nTRM==12) iStrip =  2 + (Int_t)(nTDC/3.);

  } else if (iDDL==1) {

    if (nTRM==3 && nTDC<3) iStrip = (Int_t)(nTDC/3.);
    else if (nTRM== 4) iStrip =  5 - (Int_t)(nTDC/3.);
    else if (nTRM== 5) iStrip = 10 - (Int_t)(nTDC/3.);
    else if (nTRM== 6) iStrip = 15 - (Int_t)(nTDC/3.);
    else if (nTRM== 7) {
      if (nTDC<6) iStrip =  1 - (Int_t)(nTDC/3.);
      else iStrip = 20 - (Int_t)(nTDC/3.);
    }
    else if (nTRM== 8) iStrip =  6 - (Int_t)(nTDC/3.);
    else if (nTRM== 9) iStrip = 11 - (Int_t)(nTDC/3.);
    else if (nTRM==10) iStrip = 16 - (Int_t)(nTDC/3.);
    else if (nTRM==11) {
      if (nTDC<9) iStrip =  2 - (Int_t)(nTDC/3.);
      else iStrip = 21 - (Int_t)(nTDC/3.);
    }
    else if (nTRM==12) iStrip =  7 - (Int_t)(nTDC/3.);

  } else if (iDDL==2) {

    if (nTRM== 4) iStrip =  18 - (Int_t)(nTDC/3.);
    else if (nTRM== 5) iStrip = 18 - ( 5 + (Int_t)(nTDC/3.));
    else if (nTRM== 6) iStrip = 18 - (10 + (Int_t)(nTDC/3.));
    else if (nTRM== 7) {
      if (nTDC<12) iStrip =  18 - (15 + (Int_t)(nTDC/3.));
      else iStrip = 18 - ((Int_t)(nTDC/3.) -  4);
    }
    else if (nTRM== 8) iStrip = 18 - ( 1 + (Int_t)(nTDC/3.));
    else if (nTRM== 9) iStrip = 18 - ( 6 + (Int_t)(nTDC/3.));
    else if (nTRM==10) iStrip = 18 - (11 + (Int_t)(nTDC/3.));
    else if (nTRM==11) {
      if (nTDC<9) iStrip = 18 - (16 + (Int_t)(nTDC/3.));
      else iStrip = 14 - ((Int_t)(nTDC/3.) -  3);
    }
    else if (nTRM==12) iStrip = 14 - ( 2 + (Int_t)(nTDC/3.));

  } else if (iDDL==3) {

    if (nTRM==3 && nTDC<3) iStrip = 18 - (Int_t)(nTDC/3.);
    else if (nTRM== 4) iStrip = 18 - ( 5 - (Int_t)(nTDC/3.));
    else if (nTRM== 5) iStrip = 18 - (10 - (Int_t)(nTDC/3.));
    else if (nTRM== 6) iStrip = 18 - (15 - (Int_t)(nTDC/3.));
    else if (nTRM== 7) {
      if (nTDC<6) iStrip =  18 - (1 - (Int_t)(nTDC/3.));
      else iStrip = 18 - (20 - (Int_t)(nTDC/3.));
    }
    else if (nTRM== 8) iStrip = 18 - ( 6 - (Int_t)(nTDC/3.));
    else if (nTRM== 9) iStrip = 18 - (11 - (Int_t)(nTDC/3.));
    else if (nTRM==10) iStrip = 18 - (16 - (Int_t)(nTDC/3.));
    else if (nTRM==11) {
      if (nTDC<9) iStrip = 14 - ( 2 - (Int_t)(nTDC/3.));
      else iStrip = 18 - (21 - (Int_t)(nTDC/3.));
    }
    else if (nTRM==12) iStrip = 14 - ( 7 - (Int_t)(nTDC/3.));

  } 

  return iStrip;

}

//----------------------------------------------------------------------------
Int_t AliTOFRawStream::Equip2VolNpad(Int_t iDDL, Int_t iChain, Int_t nTDC,
				Int_t iCH) const 
{
  //
  // Returns the TOF pad number per strip [0;95]
  // corresponding to the TOF equipment ID numbers:
  //                          iDDL -> DDL number per sector [0;3]
  //                        iChain -> TRM chain number [0;1]
  //                          nTDC -> TDC number [0;14]
  //                           iCH -> TDC channel number [0;7]
  //

  Int_t iPadAlongTheStrip = -1;

  Int_t iTDClocal = nTDC%3 + (1-iChain)*3;

  if (iDDL==0 || iDDL==3) iTDClocal = 5 - iTDClocal;
  else if (iDDL==1 || iDDL==2) iTDClocal = 6 + (5 - iTDClocal);

  Int_t iCHlocal = iCH;
  if (iDDL==0 || iDDL==3) iCHlocal = 7 - iCH;

  iPadAlongTheStrip = iTDClocal*AliTOFGeometry::NCh() + iCHlocal;

  if (((iDDL==1 || iDDL==2) && iPadAlongTheStrip< AliTOFGeometry::NpadX()) ||
      ((iDDL==0 || iDDL==3) && iPadAlongTheStrip>=AliTOFGeometry::NpadX()))
    AliError("Problems with the padX number!");

  return iPadAlongTheStrip;

}

//----------------------------------------------------------------------------
Int_t AliTOFRawStream::GetSectorNumber(Int_t nDDL) const
{
  //
  // Returns the sector number [0;17]
  // corresponing to the assigned DRM/DDL number [0;71]
  //

  Int_t iSector = Int_t((Float_t)(nDDL)/AliTOFGeometry::NDDL());

  return iSector;

}
//----------------------------------------------------------------------------
Int_t AliTOFRawStream::GetDDLnumberPerSector(Int_t nDDL) const
{
  //
  // Return the DRM/DDL number per sector [0;3]
  // corresponing to the assigned DRM/DDL number [0;71]
  //

  Int_t iDDL = nDDL%AliTOFGeometry::NDDL();

  return iDDL;

}

//----------------------------------------------------------------------------
void AliTOFRawStream::EquipmentId2VolumeId(Int_t nDDL, Int_t nTRM, Int_t iChain,
					Int_t nTDC, Int_t iCH,
					Int_t *volume) const
{
  //
  // To convert:
  //            nDDL   (variable in [0;71]) -> number of the DDL file 
  //            nTRM   (variable in [3;12]) -> number of the TRM slot
  //            iChain (variable in [0; 1]) -> number of the TRM chain
  //            nTDC   (variable in [0;14]) -> number of the TDC
  //            iCH    (variable in [0; 7]) -> number of the TDC channel
  //
  // in:
  //      sector number, i.e. volume[0] (variable in [0,17])
  //      plate  number, i.e. volume[1] (variable in [0, 5])
  //      strip  number, i.e. volume[2] (variable in [0,14/18])
  //      padX   number, i.e. volume[3] (variable in [0,47])
  //      padZ   number, i.e. volume[4] (variable in [0, 1])
  //

  Int_t iDDL = GetDDLnumberPerSector(nDDL);

  Int_t iSector = GetSectorNumber(nDDL);

  Int_t iPlate = Equip2VolNplate(iDDL, nTRM, nTDC);
  if (iPlate==-1) AliError("Problems with the plate number!");

  Int_t iStrip = Equip2VolNstrip(iDDL, nTRM, nTDC);
  if (iStrip==-1) AliError("Problems with the strip number!");

  Int_t iPadAlongTheStrip  = Equip2VolNpad(iDDL, iChain, nTDC, iCH);
  if (iPadAlongTheStrip==-1)
    AliError("Problems with the pad number along the strip!");

  Int_t iPadX  = (Int_t)(iPadAlongTheStrip/(Float_t(AliTOFGeometry::NpadZ())));
  Int_t iPadZ  = iPadAlongTheStrip%AliTOFGeometry::NpadZ();

  volume[0] = iSector;
  volume[1] = iPlate;
  volume[2] = iStrip;
  volume[3] = iPadX;
  volume[4] = iPadZ;

}
