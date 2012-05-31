// $Id$

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

//
// Mapping class fro mapping
// from hardware address to geometrical address
// Authors: Per Thomas Hille, Oystein Djuvsland
//

#include "AliHLTPHOSMapper.h"
#include "unistd.h"
#include <iostream>
#include "AliHLTPHOSCoordinate.h"
#include "AliPHOSRecoParam.h"
#include "AliAltroMapping.h"
#include "AliCaloAltroMapping.h"
#include "TObjArray.h"


AliHLTPHOSMapper::AliHLTPHOSMapper():
        AliHLTCaloMapper(0,"PHOS")
        ,fIsInitializedMapping(false)
        ,fDDLMapInitialised(false)
        ,fModuleId(-1)
        ,fDDLId(-1)
{
    sprintf(fFilepath, "./");
}


AliHLTPHOSMapper::~AliHLTPHOSMapper()
{
    delete []  fHw2geomapPtr;
    fHw2geomapPtr = 0;
}


Bool_t
AliHLTPHOSMapper::InitAltroMapping(const unsigned long specification)
{
    // Loads mapping between Altro addresses and geometrical addresses from file

    //  HLTError("Initialising ALTRO map");

    if (!fDDLMapInitialised) InitDDLSpecificationMapping();

    fDDLId = GetDDLFromSpec(specification);
    Int_t modId = GetModuleFromSpec(specification);

    const TObjArray* maps = AliPHOSRecoParam::GetMappings();
    if (!maps)
    {
        HLTError("Cannot retrieve ALTRO mappings!!");
        fIsInitializedMapping = false;
        return false;
    }

    AliCaloAltroMapping *map = dynamic_cast<AliCaloAltroMapping*>(maps->At(modId*fCaloConstants->GetNRCUSPERMODULE()));
    
    if(!map)
    {
        HLTError("Cannot retrieve ALTRO mappings!!");
        fIsInitializedMapping = false;
        return false;
    }

    if ( modId != fModuleId )
    {
        fModuleId = modId;

//        int nChannels = 0;
        int maxaddr = fCaloConstants->GetMAXHWADDRESSES();
        //int tmpHwaddr = 0;
        //int tmpZRow = 0;
        //int tmpXCol = 0;
        //int tmpGain = 0;
//        int res = 0;
        if (fHw2geomapPtr)
        {
            delete [] fHw2geomapPtr;
        }
        fHw2geomapPtr = new fAltromap[maxaddr];
	Int_t hwAdds[maxaddr];
	Int_t nCh = 0;
        for (int x = 0; x < fCaloConstants->GetNXCOLUMNSRCU(); x++)
        {
            for (int z = 0; z < fCaloConstants->GetNZROWSRCU(); z++)
            {
	       for(int g = 0; g < fCaloConstants->GetNGAINS(); g++)
	       {
		  hwAdds[nCh] = map->GetHWAddress(x, z, g);
		  nCh++;
	       }
            }
        }
        for ( int i=0; i < nCh; i ++ )
        {
	   Int_t add = hwAdds[i];
            if (map->GetSector(add) < 2)
            {
                fHw2geomapPtr[add].fXCol = map->GetPadRow(add);
                fHw2geomapPtr[add].fZRow = map->GetPad(add);
                fHw2geomapPtr[add].fGain = map->GetSector(add);
            }
        }
        fIsInitializedMapping = true;
    }
    else
    {
        fIsInitializedMapping = false;
    }

    return fIsInitializedMapping;
}


void
AliHLTPHOSMapper::InitDDLSpecificationMapping()
{
    fSpecificationMapPtr = new fDDLSpecificationMap[fCaloConstants->GetNMODULES()*fCaloConstants->GetNRCUSPERMODULE()];
    //  HLTError("NUMBER OF DDLs: %d, map ptr: %d", fCaloConstants->GetNMODULES()*fCaloConstants->GetNRCUSPERMODULE(), fSpecificationMapPtr);
    for (Int_t ddl = 0; ddl < fCaloConstants->GetNMODULES()*fCaloConstants->GetNRCUSPERMODULE(); ddl++)
    {

        fSpecificationMapPtr[ddl].fModId = ddl/fCaloConstants->GetNRCUSPERMODULE();

        if (ddl%4 == 0)
        {
            fSpecificationMapPtr[ddl].fRcuX = 0;
            fSpecificationMapPtr[ddl].fRcuZ = 0;
        }

        else if (ddl%4 == 1)
        {
            fSpecificationMapPtr[ddl].fRcuX = 1;
            fSpecificationMapPtr[ddl].fRcuZ = 0;
        }

        else if ( ddl%4 == 2)
        {
            fSpecificationMapPtr[ddl].fRcuX = 2;
            fSpecificationMapPtr[ddl].fRcuZ = 0;
        }
        else
        {
            fSpecificationMapPtr[ddl].fRcuX = 3;
            fSpecificationMapPtr[ddl].fRcuZ = 0;
        }

        fSpecificationMapPtr[ddl].fRcuZOffset = fCaloConstants->GetNZROWSRCU()*(fSpecificationMapPtr[ddl].fRcuZ);
        fSpecificationMapPtr[ddl].fRcuXOffset = fCaloConstants->GetNXCOLUMNSRCU()*(fSpecificationMapPtr[ddl].fRcuX);

    }
    fDDLMapInitialised = true;
}


bool
AliHLTPHOSMapper::GetIsInitializedMapping()
{
    return  fIsInitializedMapping;
}

Int_t
AliHLTPHOSMapper::GetChannelID(Int_t hwAddress)
{
    if (!fDDLMapInitialised) InitDDLSpecificationMapping();

    //  HLTError("HW add: %d -> x: %d, z: %d, gain: %d", fHw2geomapPtr[hwAddress].fXCol + fSpecificationMapPtr[fDDLId].fRcuXOffset,
    //	   fHw2geomapPtr[hwAddress].fZRow + fSpecificationMapPtr[fDDLId].fRcuZOffset,
    //	   fHw2geomapPtr[hwAddress].fGain);
    return ((fHw2geomapPtr[hwAddress].fXCol + fSpecificationMapPtr[fDDLId].fRcuXOffset) |

            ((fHw2geomapPtr[hwAddress].fZRow + fSpecificationMapPtr[fDDLId].fRcuZOffset) << 6) |
            (fHw2geomapPtr[hwAddress].fGain << 12) |
            fSpecificationMapPtr[fDDLId].fModId << 13);
}

Int_t
AliHLTPHOSMapper::GetChannelID(AliHLTUInt32_t specification, Int_t hwAddress)
{

    if (!fDDLMapInitialised) InitDDLSpecificationMapping();

    Short_t index = 0;

    if (specification == 0x00001) index = 0;
    else if (specification == 0x00002) index = 1;
    else if (specification == 0x00004) index = 2;
    else if (specification == 0x00008) index = 3;

    else if (specification == 0x00010) index = 4;
    else if (specification == 0x00020) index = 5;
    else if (specification == 0x00040) index = 6;
    else if (specification == 0x00080) index = 7;

    else if (specification == 0x00100) index = 8;
    else if (specification == 0x00200) index = 9;
    else if (specification == 0x00400) index = 10;
    else if (specification == 0x00800) index = 11;

    else if (specification == 0x01000) index = 12;
    else if (specification == 0x02000) index = 13;
    else if (specification == 0x04000) index = 14;
    else if (specification == 0x08000) index = 15;

    else if (specification == 0x10000) index = 16;
    else if (specification == 0x20000) index = 17;
    else if (specification == 0x40000) index = 18;
    else if (specification == 0x80000) index = 19;

    else HLTError("Specification 0x%X not consistent with single DDL in PHOS", specification);

    //  HLTError("Channel ID: 0x%X Coordinates: x = %d, z = %d, gain = %d", ((fHw2geomapPtr[hwAddress].fXCol + fSpecificationMapPtr[index].fRcuXOffset) |((fHw2geomapPtr[hwAddress].fZRow + fSpecificationMapPtr[index].fRcuZOffset) << 6) | (fHw2geomapPtr[hwAddress].fGain << 12) | fSpecificationMapPtr[index].fModId << 13),
    //  	   fHw2geomapPtr[hwAddress].fXCol,
    //     fHw2geomapPtr[hwAddress].fZRow,
    //  	   fHw2geomapPtr[hwAddress].fGain);

    /*    HLTError("HW add: %d -> x: %d, z: %d, gain: %d", hwAddress, fHw2geomapPtr[hwAddress].fXCol + fSpecificationMapPtr[index].fRcuXOffset,
      	   fHw2geomapPtr[hwAddress].fZRow + fSpecificationMapPtr[index].fRcuZOffset,
      	   fHw2geomapPtr[hwAddress].fGain);*/
    //  HLTError("RCU X offset: %d", fSpecificationMapPtr[index].fRcuXOffset);
    //  HLTError("RCU Z offset: %d", fSpecificationMapPtr[index].fRcuZOffset);
    return ((fHw2geomapPtr[hwAddress].fXCol + fSpecificationMapPtr[index].fRcuXOffset) |
            ((fHw2geomapPtr[hwAddress].fZRow + fSpecificationMapPtr[index].fRcuZOffset) << 6) |
            (fHw2geomapPtr[hwAddress].fGain << 12) |
            fSpecificationMapPtr[index].fModId << 13);
}



// void
// AliHLTPHOSMapper::GetChannelCoord(const UShort_t channelId, UShort_t* channelCoord)
// {
//   channelCoord[0] = channelId&0x3f;
//   channelCoord[1] = (channelId >> 6)&0x3f;
//   channelCoord[2] = (channelId >> 12)&0x1;
//   channelCoord[3] = (channelId >> 13)&0x1f;
//   //  printf("Channel ID: 0x%X Coordinates: x = %d, z = %d, gain = %d\n", channelId, channelCoord[0], channelCoord[1], channelCoord[2]);
// }
//
//
//
// void
// AliHLTPHOSMapper::ChannelId2Coordinate(const UShort_t channelId,    AliHLTPHOSCoordinate &channelCoord)
// {
//   channelCoord.fX = channelId&0x3f;
//   channelCoord.fZ = (channelId >> 6)&0x3f;
//   channelCoord.fGain = (channelId >> 12)&0x1;
//   channelCoord.fModuleId  = (channelId >> 13)&0x1f;
//   //  printf("Channel ID: 0x%X Coordinates: x = %d, z = %d, gain = %d\n", channelId, channelCoord[0], channelCoord[1], channelCoord[2]);
// }
//
//
//
// void
// AliHLTPHOSMapper::GetLocalCoord(const UShort_t channelId, Float_t* channelCoord)
// {
//   channelCoord[0] = (static_cast<Float_t>(channelId&0x3f) - NXCOLUMNSMOD/2)* fCellStep;
//   channelCoord[1] = (static_cast<Float_t>((channelId >> 6)&0x3f) - NZROWSMOD/2) * fCellStep;
//   //  printf("Local coordinates: x = %f, z = %f\n", channelCoord[0], channelCoord[1]);
// }

Int_t
AliHLTPHOSMapper::GetDDLFromSpec(AliHLTUInt32_t specification)
{
    Int_t index = -1;
    if (specification == 0x00001) index = 0;
    else if (specification == 0x00002) index = 1;
    else if (specification == 0x00004) index = 2;
    else if (specification == 0x00008) index = 3;

    else if (specification == 0x00010) index = 4;
    else if (specification == 0x00020) index = 5;
    else if (specification == 0x00040) index = 6;
    else if (specification == 0x00080) index = 7;

    else if (specification == 0x00100) index = 8;
    else if (specification == 0x00200) index = 9;
    else if (specification == 0x00400) index = 10;
    else if (specification == 0x00800) index = 11;

    else if (specification == 0x01000) index = 12;
    else if (specification == 0x02000) index = 13;
    else if (specification == 0x04000) index = 14;
    else if (specification == 0x08000) index = 15;

    else if (specification == 0x10000) index = 16;
    else if (specification == 0x20000) index = 17;
    else if (specification == 0x40000) index = 18;
    else if (specification == 0x80000) index = 19;

    else HLTError("Specification 0x%X not consistent with single DDL in PHOS", specification);

    return index;
}

Int_t
AliHLTPHOSMapper::GetModuleFromSpec(AliHLTUInt32_t specification)
{
    Int_t module = -1;

    if (specification & 0xf) module = 0;
    else if ((specification >> 4) & 0xf) module = 1;
    else if ((specification >> 8) & 0xf) module = 2;
    else if ((specification >> 12) & 0xf) module = 3;
    else if ((specification >> 16) & 0xf) module = 4;

    else HLTError("Specification 0x%X not consistent with single module in PHOS", specification);

    return module;
}
