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

//-------------------------------------------------------------------------
//
// AliTPCmapper
// Authors: Christian.Lippmann@cern.ch, J.Wiechula@gsi.de
// Class to map detector coordinates (row, pad, sector, ...) to
// hardware coordinates (RCU, Branch, FEC, Altro, channel, Equipment ID, ...)
// 
// There are two different ways to number padrows:
// 1) local padrow: for each ROC, 0 ... 62 for an IROC, 0 ... 95 for an OROC,
// 2) global padrow: for each sector, from 0 ... 158.
// If the global numbering is used, it is denoted by the variable name
// globalpadrow in this class.
//
// There are two different ways to number sectors:
// 1) Sectors contain one IROC and one OROC and are counted from 0 to 17 on
//    each of the two sides (A=0 and C=1),
// 2) ROCs are numbered from 0 to 71 where the ROCs 0 ... 35 are IROCS and
//    ROCs 36 ... 71 are OROCs. A ROC is often named "sector" in aliroot,
//    which can be very confusing!
//
//-------------------------------------------------------------------------

#include <TMath.h>
#include <TSystem.h>
#include <TString.h>

#include "AliTPCmapper.h"
#include "AliTPCAltroMapping.h"
#include "AliTPCROC.h"
#include "AliLog.h"
#include "AliDAQ.h"

ClassImp(AliTPCmapper)
//______________________________________________________________
AliTPCmapper::AliTPCmapper() :
  fNside(0),
  fNsector(0),
  fNrcu(0),
  fNbranch(0),
  fNaltro(0),
  fNchannel(0),
  fNpadrow(0),
  fNpadrowIROC(0),
  fNpadrowOROC(0),
  fTpcDdlOffset(0)
{
  //
  // Constructor
  //
  for ( Int_t i = 0; i < 6; i++ )  fMapping[i]=0;
}

//______________________________________________________________
AliTPCmapper::AliTPCmapper(const char * dirname) :
  fNside(0),
  fNsector(0),
  fNrcu(0),
  fNbranch(0),
  fNaltro(0),
  fNchannel(0),
  fNpadrow(0),
  fNpadrowIROC(0),
  fNpadrowOROC(0),
  fTpcDdlOffset(0)
{
  //
  // Constructor
  //
  // dirname - specify the directory with the ascii Altro mapping files
  //
  Init(dirname);
}

//______________________________________________________________
AliTPCmapper::~AliTPCmapper()
{
  // Destructor

  for ( Int_t i = 0; i < fNrcu; i++ ) {
    delete fMapping[i];
    fMapping[i] = 0;
  }
}


//_____________________________________________________________________________
AliTPCmapper::AliTPCmapper(const AliTPCmapper& mapper) :
  TObject(mapper),
  fNside(mapper.fNside),
  fNsector(mapper.fNsector),
  fNrcu(mapper.fNrcu),
  fNbranch(mapper.fNbranch),
  fNaltro(mapper.fNaltro),
  fNchannel(mapper.fNchannel),
  fNpadrow(mapper.fNpadrow),
  fNpadrowIROC(mapper.fNpadrowIROC),
  fNpadrowOROC(mapper.fNpadrowOROC),
  fTpcDdlOffset(mapper.fTpcDdlOffset)
{
  // Copy Constructor
  for ( Int_t i = 0; i < 6; i++ )  fMapping[i]=0;
  for ( Int_t i = 0; i < fNrcu; i++ ) fMapping[i] = mapper.fMapping[i];
}

//_____________________________________________________________________________
AliTPCmapper& AliTPCmapper::operator = (const AliTPCmapper& mapper)
{
  // Assignment operator

  if(&mapper == this) return *this;
  ((TObject *)this)->operator=(mapper);

  for ( Int_t i = 0; i < fNrcu; i++ ) fMapping[i] = mapper.fMapping[i];

  fNside = mapper.fNside;
  fNsector = mapper.fNsector;
  fNrcu = mapper.fNrcu;
  fNbranch = mapper.fNbranch;
  fNaltro = mapper.fNaltro;
  fNchannel = mapper.fNchannel;
  fNpadrow = mapper.fNpadrow;
  fNpadrowIROC = mapper.fNpadrowIROC;
  fNpadrowOROC = mapper.fNpadrowOROC;
  fTpcDdlOffset = mapper.fTpcDdlOffset;

  return *this;
}

//______________________________________________________________
void AliTPCmapper::Init(const char *dirname)
{
  // Initialize all
  fNside    = 2;
  fNsector  = 18;
  fNrcu     = 6;
  fNbranch  = 2;
  fNaltro   = 8;
  fNchannel = 16;

  // Load and read mapping files. AliTPCAltroMapping contains the mapping for
  // each patch (rcu).
  TString path;
  if (dirname==0){
    path  =gSystem->Getenv("ALICE_ROOT");
    path += "/TPC/mapping/Patch";
  }else{
    path  = dirname;
    path +="Patch";
  }

  TString path2;
  for(Int_t i = 0; i < fNrcu; i++) {
    path2 = path;
    path2 += i;
    path2 += ".data";
    fMapping[i] = new AliTPCAltroMapping(path2.Data());
  }

  // Get instance of AliTPCROC object
  AliTPCROC *fROC = AliTPCROC::Instance();
  fNpadrowIROC = fROC->GetNRows(0);
  fNpadrowOROC = fROC->GetNRows(36);
  fNpadrow     = fNpadrowIROC+fNpadrowOROC;

  AliDAQ daq;
  fTpcDdlOffset = daq.DdlIDOffset("TPC");

}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetHWAddress(Int_t roc, Int_t padrow, Int_t pad) const
{
  // Get the hardware address from pad coordinates for a given ROC
  Int_t patch = GetPatch(roc, padrow, pad);
  if ( (patch >= fNrcu) || (patch < 0) ) return -1;
  return fMapping[patch]->GetHWAddress(padrow, pad, roc);
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetHWAddressSector(Int_t globalpadrow, Int_t pad) const
{
  // Get the hardware address from pad coordinates
  Int_t patch = 0;
  Int_t hwAddress=-1;
  if ( globalpadrow < fNpadrowIROC   ) {
    patch = GetPatch(0,  globalpadrow, pad);
    if (patch>-1)
      hwAddress = fMapping[patch]->GetHWAddress(globalpadrow, pad, 0);
  } else if ( globalpadrow < fNpadrow ) {
    patch = GetPatch(36, globalpadrow - fNpadrowIROC, pad);
    if (patch>-1)
      hwAddress = fMapping[patch]->GetHWAddress(globalpadrow - fNpadrowIROC, pad, 36);
  } else {
    AliWarning(Form("Padrow outside range (globalpadrow %d) !", globalpadrow));
    hwAddress = -1;
  }
  return hwAddress;
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetRcu(Int_t roc, Int_t padrow, Int_t pad) const
{
  // Get the patch (rcu) index from the pad coordinates. The Roc index is
  // needed as well to determine if it is IROC or OROC. 
  return GetPatch(roc, padrow, pad);
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetPatch(Int_t roc, Int_t padrow, Int_t pad) const
{
  // Get the patch (rcu) index from the pad coordinates. The Roc index is
  // needed as well to determine if it is IROC or OROC. 

  if ( (padrow < 0) || (pad < 0) || (roc < 0) ) {
    AliWarning(Form("Pad coordinates outside range (padrow %d, pad %d, roc %d) !", padrow, pad, roc));
    return -1;
  }

  if ( roc < 36 ) {
    // IROC (ROCs 0 ... 35)
    Int_t padsInRow = GetNpads(padrow);
    if ( (padsInRow < 0) || (pad >= padsInRow) ) {
      AliWarning(Form("Pad index outside range (padrow %d, pad %d, roc %d) !", padrow, pad, roc));
      return -1;
    }
    if ( padrow < 30 ) {                  return 0;
    } else if ( padrow == 30 ) {          // padrow 30 is shared between rcus 0 and 1
      if ( (pad < 37) || (pad > 48) )     return 1;
      else                                return 0;
    } else if ( padrow < fNpadrowIROC ) { return 1;
    } else {
      AliWarning(Form("Padrow outside range (padrow %d, roc %d) !", padrow, roc));
      return -1;
    }
  } else if ( roc < 72 ) {
    // OROC (ROCs 36 ... 71)
    Int_t padsInRow = GetNpads(fNpadrowIROC+padrow);
    if ( (padsInRow < 0) || (pad >= padsInRow) ) {
      AliWarning(Form("Pad index outside range (padrow %d, pad %d, roc %d) !", padrow, pad, roc));
      return -1;
    }
    if ( padrow < 27 ) {                  return 2;
    } else if ( padrow == 27 ) {          // padrow 27 is shared between rcus 2 and 3
      if ( (pad >= 43) && (pad <= 46) )   return 3;
      else                                return 2;
    } else if ( padrow < 54 ) {           return 3;
    } else if ( padrow < 76 ) {           return 4;
    } else if ( padrow == 76) {           // padrow 76 is shared between rcus 4 and 5
      if ( (pad >= 33) && (pad <= 88) )   return 5;
      else                                return 4;
    } else if ( padrow < fNpadrowOROC ) { return 5;
    } else {
      AliWarning(Form("Padrow outside range (padrow %d, roc %d) !", padrow, roc));
      return -1;
    }
  }
  return -1;
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetRcuSector(Int_t globalpadrow, Int_t pad) const
{
  // Get the patch (rcu) index from the pad coordinates for a sector
  return GetPatchSector(globalpadrow, pad);
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetPatchSector(Int_t globalpadrow, Int_t pad) const
{
  // Get the patch (rcu) index from the pad coordinates for a sector
  if ( globalpadrow >= fNpadrow ) {
    AliWarning(Form("Padrow outside range (globalpadrow %d) !", globalpadrow));
    return -1;
  }
  if ( globalpadrow < fNpadrowIROC ) return GetPatch(0,  globalpadrow, pad);
  else                               return GetPatch(36, globalpadrow-fNpadrowIROC, pad);
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetPadRow(Int_t patch, Int_t hwAddress) const
{
  // Get Pad Row (for a ROC) from the hardware address
  if ( (patch >= fNrcu) || (patch < 0) ) {
    AliWarning(Form("Patch index outside range (patch %d) !", patch));
    return -1;
  }
  return fMapping[patch]->GetPadRow(hwAddress);
}


//_____________________________________________________________________________
  Int_t AliTPCmapper::GetGlobalPadRow(Int_t patch, Int_t hwAddress) const
{
  // Get Pad Row (for full sector) from the hardware address
  if ( patch < 2 ) return GetPadRow(patch, hwAddress);
  else             return GetPadRow(patch, hwAddress) + fNpadrowIROC;
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetPad(Int_t patch, Int_t hwAddress) const
{
  // Get Pad index from the hardware address
  if ( (patch >= fNrcu) || (patch < 0) ) {
    AliWarning(Form("Patch index outside range (patch %d) !", patch));
    return -1;
  }
  return fMapping[patch]->GetPad(hwAddress);
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetPadRow(Int_t patch, Int_t branch, Int_t fec, Int_t chip,
			      Int_t channel) const
{
  // Get pad row (for a ROC) from hardware coordinates
  if ( (patch >= fNrcu) || (branch >= fNbranch) || (chip >= fNaltro) || (channel >= fNchannel)
       || (patch < 0) || (branch < 0) || (chip < 0) || (channel < 0) ) {
    AliWarning(Form("Coordinates outside range (patch %d, branch %d, fec %d, chip %d, channel %d)) !",
		    patch, branch, fec, chip, channel));
    return -1;
  }
  return GetPadRow(patch, CodeHWAddress(branch, fec, chip, channel));
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetGlobalPadRow(Int_t patch, Int_t branch, Int_t fec, Int_t chip,
				    Int_t channel) const
{
  // Get Pad Row (for full sector) from the hardware address
  if ( (patch >= fNrcu) || (branch >= fNbranch) || (chip >= fNaltro) || (channel >= fNchannel)
       || (patch < 0) || (branch < 0) || (chip < 0) || (channel < 0) ) {
    AliWarning(Form("Coordinates outside range (patch %d, branch %d, fec %d, chip %d, channel %d)) !",
		    patch, branch, fec, chip, channel));
    return -1;
  }
  if ( patch < 2 ) return GetPadRow(patch, branch, fec, chip, channel);
  else             return GetPadRow(patch, branch, fec, chip, channel) + fNpadrowIROC;
}


//_____________________________________________________________________________
  Int_t AliTPCmapper::GetPad(Int_t patch, Int_t branch, Int_t fec, Int_t chip, Int_t channel) const
{
  // Get pad from hardware coordinates
  if ( (patch >= fNrcu) || (branch >= fNbranch) || (chip >= fNaltro) || (channel >= fNchannel)
       || (patch < 0) || (branch < 0) || (chip < 0) || (channel < 0) ) {
    AliWarning(Form("Coordinates outside range (patch %d, branch %d, fec %d, chip %d, channel %d)) !",
		    patch, branch, fec, chip, channel));
    return -1;
  }
  return GetPad(patch, CodeHWAddress(branch, fec, chip, channel));
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetBranch(Int_t roc, Int_t padrow, Int_t pad) const
{
  // Get the branch to which this pad is connected. The FECs connected to
  // one RCU are divided into two branches: A(=0) and B(=1). This information
  // can be extracted from the hardware address.
  return DecodedHWAddressBranch(GetHWAddress(roc, padrow, pad));
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetBranchSector(Int_t globalpadrow, Int_t pad) const
{
  // Get Branch from pad coordinates, where globalpadrow is counted
  // for a full sector (0 ... 158)
  return DecodedHWAddressBranch(GetHWAddressSector(globalpadrow, pad));
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetFEChw(Int_t roc, Int_t padrow, Int_t pad) const
{
  // Get the FEC number in hardware numbering. The FECs are numbered from 0 (in the
  // center of the partition) to 8 (partition 3, 4, 5), 9 (partition 0, 2), 11
  // (partition 1, branch A) or 12 (partition 1, branch B). This information
  // can be extracted from the hardware address.
  return DecodedHWAddressFECaddr(GetHWAddress(roc, padrow, pad));
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetFEChwSector(Int_t globalpadrow, Int_t pad) const
{
  // Get the FEC number in hardware numbering from pad coordinates, where 
  // globalpadrow is counted for a full sector (0 ... 158)
  return DecodedHWAddressFECaddr(GetHWAddressSector(globalpadrow, pad));
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetFEC(Int_t roc, Int_t padrow, Int_t pad) const
{
  // Get the FEC number in offline-oriented numbering. The FECs are numbered from 0
  // 17 (partition 3, 4, 5), 19 (partition 0, 2) or 24 (partition 1).
  Int_t patch  = GetPatch(roc, padrow, pad);
  Int_t fec    = DecodedHWAddressFECaddr(GetHWAddress(roc, padrow, pad));
  Int_t branch = DecodedHWAddressBranch(GetHWAddress(roc, padrow, pad));
  if ( (fec < 0) || (branch < 0) || (patch < 0) ) return -1;
  return HwToOffline(patch, branch, fec);
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetFECSector(Int_t globalpadrow, Int_t pad) const
{
  // Get the FEC number in offline-oriented numbering. globalpadrow is
  // counted for a full sector (0 ... 158)
  Int_t patch  = GetPatchSector(globalpadrow, pad);
  Int_t fec    = DecodedHWAddressFECaddr(GetHWAddressSector(globalpadrow, pad));
  Int_t branch = DecodedHWAddressBranch(GetHWAddressSector(globalpadrow, pad));
  if ( (fec < 0) || (branch < 0) || (patch < 0) ) return -1;
  return HwToOffline(patch, branch, fec);
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetChip(Int_t roc, Int_t padrow, Int_t pad) const
{
  // Get Chip (ALTRO) index (0 ... 7) from pad coordinates
  return DecodedHWAddressChipaddr(GetHWAddress(roc, padrow, pad));
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetChipSector(Int_t globalpadrow, Int_t pad) const
{
  // Get Chip (ALTRO) index (0 ... 7) from pad coordinates, where 
  // globalpadrow is counted for a full sector (0 ... 158)
  return DecodedHWAddressChipaddr(GetHWAddressSector(globalpadrow, pad));
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetChannel(Int_t roc, Int_t padrow, Int_t pad) const
{
  // Get Channel index (0 ... 15) from pad coordinates
  return DecodedHWAddressChanneladdr(GetHWAddress(roc, padrow, pad));
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetChannelSector(Int_t globalpadrow, Int_t pad) const
{
  // Get Channel index (0 ... 15) from pad coordinates, where 
  // globalpadrow is counted for a full sector (0 ... 158)
  return DecodedHWAddressChanneladdr(GetHWAddressSector(globalpadrow, pad));
}


//_____________________________________________________________________________
Int_t AliTPCmapper::CodeHWAddress(Int_t branch, Int_t fec, Int_t chip, Int_t channel) const
{
  // Get Hardware address from channel, altro, fec and branch coordinates
  return ((branch&1)<<11) + ((fec&0xf)<<7) + ((chip&0x7)<<4) + (channel&0xf);
}


//_____________________________________________________________________________
Int_t AliTPCmapper::DecodedHWAddressBranch(Int_t hwAddress) const
{
  // Get branch index (0, 1) from hardware address
  if ( hwAddress < 0 ) return -1;
  return ((hwAddress>>11)&1);
}


//_____________________________________________________________________________
Int_t AliTPCmapper::DecodedHWAddressFECaddr(Int_t hwAddress) const
{
  // Get FEC index (0 ... 12) from hardware address
  if ( hwAddress < 0 ) return -1;
  return ((hwAddress>>7)&0xf);
}


//_____________________________________________________________________________
Int_t AliTPCmapper::DecodedHWAddressChipaddr(Int_t hwAddress) const
{
  // Get ALTRO index (0 ... 7) from hardware address
  if ( hwAddress < 0 ) return -1;
  return ((hwAddress>>4)&0x7);
}


//_____________________________________________________________________________
Int_t AliTPCmapper::DecodedHWAddressChanneladdr(Int_t hwAddress) const
{
  // Get channel index (0 ... 15) from hardware address
  if ( hwAddress < 0 ) return -1;
  return ((hwAddress&0xf));
}


//______________________________________________________________
Int_t AliTPCmapper::GetNpads(Int_t roc, Int_t padrow) const{
  // Get number of pads in padrow for this ROC.
  AliTPCROC *fROC = AliTPCROC::Instance();
  Int_t retval = fROC->GetNPads((UInt_t)roc, (UInt_t)padrow);
  return retval;
}


//______________________________________________________________
Int_t AliTPCmapper::GetNpads(Int_t globalpadrow) const{
  // Get number of pads in padrow, where globalpadrow is counted for a full sector (0 ... 158)

  if ( globalpadrow >= fNpadrow ) {
    AliWarning(Form("Padrow outside range (globalpadrow %d) !", globalpadrow));
    return -1;
  }
  if ( globalpadrow < fNpadrowIROC ) return GetNpads(0,  globalpadrow);  // IROC
  else  return GetNpads(36, globalpadrow - fNpadrowIROC);                // OROC

  return -1;
}


//______________________________________________________________
Int_t AliTPCmapper::GetNpadrows(Int_t roc) const
{
  // Get number of padrows
  if      (roc < 36) return fNpadrowIROC;
  else if (roc < 72) return fNpadrowOROC;
  return -1;
}


//______________________________________________________________
/*
Double_t AliTPCmapper::GetPadXlocal(Int_t globalpadrow) const
{
  // Get local x coordinate of pad, where globalpadrow is counted for a full sector (0 ... 158)

  if ( globalpadrow >= fNpadrow ) {
    AliWarning(Form("Padrow outside range (globalpadrow %d) !", globalpadrow));
    return -1.0;
  }

  //IROC
  if ( globalpadrow < fNpadrowIROC )
    return (852.25 + 7.5*(Double_t)globalpadrow)/10.; //divide by 10 to get cm

  globalpadrow -= fNpadrowIROC;

  if ( globalpadrow < 64 ) //OROC inner part
    return (10.* globalpadrow + 1351.)/10.;         //divide by 10 to get cm

  //OROC outer part
  return (15.*(globalpadrow - 64) + 1993.5)/10.;    //divide by 10 to get cm
}
*/

//______________________________________________________________
/*
Double_t AliTPCmapper::GetPadYlocal(Int_t globalpadrow, Int_t pad) const
{
  // Get local y coordinate of pad, where globalpadrow is counted for a full sector (0 ... 158)

  if ( globalpadrow >= fNpadrow ) {
    AliWarning(Form("Padrow outside range (globalpadrow %d) !", globalpadrow));
    return -1.0;
  }

  Int_t padsInRow = GetNpads(globalpadrow);
  if ( (padsInRow < 0) || (pad >= padsInRow) ) {
    AliWarning(Form("Pad index outside range (pad %d) !", pad));
    return -1.0;
  }

  //IROC
  if ( globalpadrow < fNpadrowIROC )
    return (2.* padsInRow - 4.*pad - 2.)*1.e-1;  //divide by 10 to get cm

  //OROC
  return (3.* padsInRow -6.*pad - 3.)*1.e-1;  //divide by 10 to get cm
}
*/

//______________________________________________________________
/*
Double_t AliTPCmapper::GetPadXglobal(Int_t globalpadrow, Int_t pad, Int_t sector) const
{
  // Get global x coordinate of pad, where globalpadrow is counted for a full sector (0 ... 158)

  if ( globalpadrow >= fNpadrow ) {
    AliWarning(Form("Padrow outside range (globalpadrow %d) !", globalpadrow));
    return -1.0;
  }

  Int_t padsInRow = GetNpads(globalpadrow);
  if ( (padsInRow < 0) || (pad >= padsInRow) ) {
    AliWarning(Form("Pad index outside range (pad %d) !", pad));
    return -1.0;
  }

  Double_t angle = (Double_t)(( sector * 20. ) + 10. ) * TMath::DegToRad();
  return GetPadXlocal(globalpadrow) * TMath::Cos(angle) -
    GetPadYlocal(globalpadrow, pad) * TMath::Sin(angle);
}
*/

//______________________________________________________________
/*
Double_t AliTPCmapper::GetPadYglobal(Int_t globalpadrow, Int_t pad,Int_t sector) const
{
  // Get global y coordinate of pad, where globalpadrow is counted for a full sector (0 ... 158)

  if ( globalpadrow >= fNpadrow ) {
    AliWarning(Form("Padrow outside range (globalpadrow %d) !", globalpadrow));
    return -1.0;
  }

  Int_t padsInRow = GetNpads(globalpadrow);
  if ( (padsInRow < 0) || (pad >= padsInRow) ) {
    AliWarning(Form("Pad index outside range (pad %d) !", pad));
    return -1.0;
  }

  Double_t angle = (Double_t)(( sector * 20. ) + 10. ) * TMath::DegToRad();
  return GetPadXlocal(globalpadrow) * TMath::Sin(angle) +
    GetPadYlocal(globalpadrow, pad) * TMath::Cos(angle);
}
*/

//______________________________________________________________
/*
Double_t AliTPCmapper::GetPadWidth(Int_t globalpadrow) const
{
  //  Get pad width, where globalpadrow is counted for a full sector (0 ... 158)

  if ( globalpadrow >= fNpadrow ) {
    AliWarning(Form("Padrow outside range (globalpadrow %d) !", globalpadrow));
    return -1.0;
  }

  if (globalpadrow < fNpadrowIROC ) // IROC
    return 0.4;
  return 0.6;
}
*/

//______________________________________________________________
/*
Double_t AliTPCmapper::GetPadLength(Int_t globalpadrow) const
{
  // Get pad length, where globalpadrow is counted for a full sector (0 ... 158)

  if ( globalpadrow >= fNpadrow ) {
    AliWarning(Form("Padrow outside range (globalpadrow %d) !", globalpadrow));
    return -1.0;
  }

  if ( globalpadrow < fNpadrowIROC ) return  0.75;
  if ( globalpadrow < 127 )          return 1.0;
  return 1.5;
}
*/

//_____________________________________________________________________________
Int_t AliTPCmapper::GetNfec(Int_t patch) const
{
  // Get size of readout partition (number of FECs) for this rcu (patch) index (0 ... 5)
  Int_t retval = 0;
  switch(patch){
  case(0):
    retval = 18;
    break;
  case(1):
    retval = 25;
    break;
  case(2):
    retval = 18;
    break;
  case(3):
    retval = 20;
    break;
  case(4):
    retval = 20;
    break;
  case(5):
    retval = 20;
    break;
  };
  return retval;
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetNfec(Int_t patch, Int_t branch) const
{
  // Get size of readout partition (number of FECs) for this branch
  Int_t retval = 0;
  switch(patch){
  case(0):
    retval = 9;
    break;
  case(1):
    retval = 13;
    break;
  case(2):
    retval = 9;
    break;
  case(3):
    retval = 10;
    break;
  case(4):
    retval = 10;
    break;
  case(5):
    retval = 10;
    break;
  };
  
  if( (branch == 1) && (patch == 1) ){
    retval = 12;
  }

  return retval;
}


//_____________________________________________________________________________
Int_t AliTPCmapper::OfflineToHwFec(Int_t patch, Int_t fec) const
{
  // Convert FEC position in offline-like numbering to hardware numbering (fec).

  if ( (patch < 0) || (fec < 0) ) {
    AliWarning(Form("Patch (%d) or Fec number (%d) outside range !", patch, fec));
    return -1;
  }
  if ( (patch > 5) || (fec >= GetNfec(patch)) ) {
    AliWarning(Form("Patch (%d) or Fec number (%d) outside range !", patch, fec));
    return -1;
  }

  Int_t fecsInBranchA = GetNfec(patch, 0);
  if ( fec < fecsInBranchA ) // branch A
    return (fecsInBranchA - 1 - fec);
  else                       // branch B
    return (fec - fecsInBranchA);

  return -1;
}


//_____________________________________________________________________________
Int_t AliTPCmapper::OfflineToHwBranch(Int_t patch, Int_t fec) const
{
  // Convert fec position in offline-like numbering to hardware numbering (branch).

  if ( (patch < 0) || (fec < 0) ) {
    AliWarning(Form("Patch (%d) or Fec number (%d) outside range !", patch, fec));
    return -1;
  }
  if ( (patch > 5) || (fec >= GetNfec(patch)) ) {
    AliWarning(Form("Patch (%d) or Fec number (%d) outside range !", patch, fec));
    return -1;
  }
  if ( fec < GetNfec(patch, 0) ) return 0; // branch A
  else                                         return 1; // branch B

  return -1;
}


//_____________________________________________________________________________
Int_t AliTPCmapper::HwToOffline(Int_t patch, Int_t branch, Int_t fec) const
{
  // Convert hardware FEC position (branch, fec) to the offline-oriented numbering

  if ( (patch < 0) || (fec < 0) || (branch < 0) ) {
    AliWarning(Form("Patch (%d), branch (%d) or Fec number (%d) outside range !", patch, branch, fec));
    return -1;
  }
  if ( (patch > 5) || (branch > 1) || (fec >= GetNfec(patch, branch)) ) {
    AliWarning(Form("Patch (%d), branch (%d) or Fec number (%d) outside range !", patch, branch, fec));
    return -1;
  }
  Int_t fecsInBranchA = GetNfec(patch, 0);
  if ( branch == 0 )  // branch A
    return (fecsInBranchA - 1 - fec);
  else                // branch B
    return (fec + fecsInBranchA);

  return -1;
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetEquipmentID(Int_t roc, Int_t padrow, Int_t pad) const
{
  // Get EqID from pad coordinate. The Roc index is
  // needed as well to determine if it is IROC or OROC. 

  Int_t side = GetSideFromRoc(roc);
  if ( side < 0 ) return -1;
  Int_t sector = GetSectorFromRoc(roc);
  if ( sector < 0 ) return -1;
  Int_t patch = GetPatch(roc, padrow, pad);
  if ( patch < 0 ) return -1;
  return GetEquipmentIDfromPatch(side, sector, patch);
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetEquipmentIDsector(Int_t side, Int_t sector, Int_t globalpadrow, Int_t pad) const
{
  // Get EqID from pad coordinate, where padrow is counted for a full sector (0 ... 158)
  Int_t patch = GetPatchSector(globalpadrow, pad);
  if ( patch < 0 ) return -1;
  Int_t roc = GetRocFromPatch(side, sector, patch);
  if ( roc < 0 ) return -1;

  if ( globalpadrow < fNpadrowIROC )
    return GetEquipmentID(roc, globalpadrow, pad);
  else
    return GetEquipmentID(roc, globalpadrow-fNpadrowIROC, pad);
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetEquipmentIDfromPatch(Int_t side, Int_t sector, Int_t patch) const
{
  // Get EqID from patch (rcu).

  Int_t roc = GetRocFromPatch(side, sector, patch);
  Int_t ddl = 0;
  if (patch < 2)  // IROC
    ddl = roc*2 + patch;
  else            // OROC
    ddl = (roc-36)*4 + 36*2 + (patch-2);
  // Add offset. TPC has detectorID = 3
  return ddl+fTpcDdlOffset;
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetPatchFromEquipmentID(Int_t equipmentID) const
{
  // Get rcu (patch) index (0 ... 5) from equipment ID
  Int_t retval = 0;

  if ( (equipmentID < fTpcDdlOffset) || (equipmentID > 983) ) {
    AliWarning(Form("Equipment ID (%d) outside range !", equipmentID));
    return -1;
  }
  if ( ( (int)equipmentID - 840 ) < 0) retval = (equipmentID-768)%2;
  else                                 retval = (equipmentID-840)%4 + 2;
  return retval;
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetSideFromEquipmentID(Int_t equipmentID) const
{
  // Get side from Eq ID
  if ( equipmentID < fTpcDdlOffset ) {
    AliWarning(Form("Equipment ID (%d) outside range !", equipmentID));
    return -1;
  }
  if      ( equipmentID < 804 ) return 0;
  else if ( equipmentID < 840 ) return 1;
  else if ( equipmentID < 912 ) return 0;
  else if ( equipmentID < 984 ) return 1;
  else {
    AliWarning(Form("Equipment ID (%d) outside range !", equipmentID));
    return -1;
  }
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetSectorFromEquipmentID(Int_t equipmentID) const
{
  // Get sector index (0 ... 17) from equipment ID
  Int_t retval = 0;
  if ( (equipmentID < fTpcDdlOffset) || (equipmentID >= fTpcDdlOffset+216) ) {
    AliWarning(Form("Equipment ID (%d) outside range !", equipmentID));
    return -1;
  }
  Int_t side   = GetSideFromEquipmentID(equipmentID);
  if ( side < 0 ) return -1;

  if ( (equipmentID - 840) < 0 ) { // IROC
    if ( side == 0 ) retval = (equipmentID-fTpcDdlOffset)/2;
    else             retval = (equipmentID-fTpcDdlOffset-18*2)/2;
  } else {                         // OROC
    if ( side == 0 ) retval = (equipmentID-840)/4;
    else             retval = (equipmentID-840-18*4)/4;
  }
  return retval;
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetRocFromEquipmentID(Int_t equipmentID) const
{
  // Get ROC index (0 ... 71) from equipment ID
  Int_t side   = GetSideFromEquipmentID(equipmentID);
  if ( side < 0 ) return -1;
  Int_t sector = GetSectorFromEquipmentID(equipmentID);
  if ( sector < 0 ) return -1;
  Int_t patch  = GetPatchFromEquipmentID(equipmentID);
  if ( patch < 0 ) return -1;

  return GetRocFromPatch(side, sector, patch);
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetSectorFromRoc(Int_t roc) const
{
  // get the sector number (0 ... 17) from the roc number (0 ... 71)

  if ( roc < 0 ) {
    AliWarning(Form("Roc outside range (roc %d) !", roc));
    return -1;
  } else if ( roc < 18 ) {   // inner sector, A side
    return roc;
  } else if ( roc < 36 ) {   // inner sector, C side
    return (roc-18);
  } else if ( roc < 54 ) {   // outer sector, A side
    return (roc-36);
  } else if ( roc < 72 ) {   // outer sector, C side
    return (roc-54);
  } else {
    AliWarning(Form("Roc outside range (roc %d) !", roc));
    return -1;
  }
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetSideFromRoc(Int_t roc) const
{
  // get the side (0, 1) from the roc number (0 ... 71)

  if ( roc < 0 ) {
    AliWarning(Form("Roc outside range (roc %d) !", roc));
    return -1;
  } else if ( roc < 18 ) {   // inner sector, A side
    return 0;
  } else if ( roc < 36 ) {   // inner sector, C side
    return 1;
  } else if ( roc < 54 ) {   // outer sector, A side
    return 0;
  } else if ( roc < 72 ) {   // outer sector, C side
    return 1;
  } else { 
    AliWarning(Form("Roc outside range (roc %d) !", roc));
    return -1;
  } 
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetRocFromPatch(Int_t side, Int_t sector, Int_t patch) const
{
  // Get Roc (0 ... 71) from side (0, 1), sector (0 ... 17) and patch (0 ... 5)

  if ( (side < 0) || (side >= fNside) ) {
    AliWarning(Form("Side outside range (side %d) !", side));
    return -1;
  }
  if ( (sector < 0) || (sector >= fNsector) ) {
    AliWarning(Form("Sector outside range (sector %d) !", sector));
    return -1;
  }
  if ( (patch < 0) || (patch >= fNrcu) ) {
    AliWarning(Form("Patch (rcu) outside range (patch %d) !", patch));
    return -1;
  }

  if ( side == 0 ) { // A side
    if ( patch < 2 ) return sector;         // IROC
    else             return 36+sector;      // OROC
  } else {           // C side
    if ( patch < 2 ) return 18+sector;      // IROC
    else             return 54+sector;      // OROC
  }
}


//_____________________________________________________________________________
Int_t AliTPCmapper::GetRoc(Int_t side, Int_t sector, Int_t globalpadrow, Int_t pad) const
{
  // Get Roc (0 ... 71) from side (0, 1), sector (0 ... 17) and pad coordinates

  Int_t patch = GetPatchSector(globalpadrow, pad);
  if ( patch < 0 ) return -1;
  return GetRocFromPatch(side, sector, patch);
}


//_____________________________________________________________________________
  Bool_t AliTPCmapper::IsIROC(Int_t roc) const
{
  // Is this ROC an IROC?
  if ( roc < 0 ) {
    AliWarning(Form("Roc outside range (roc %d) !", roc));
    return -1;
  } else if ( roc < 36 ) {   // inner sector
    return true;
  } else if ( roc < 72 ) {   // outer sector, C side
    return false;
  } else {
    AliWarning(Form("Roc outside range (roc %d) !", roc));
    return -1;
  }
}


//_____________________________________________________________________________
  Bool_t AliTPCmapper::IsOROC(Int_t roc) const
{
  // Is this ROC an OROC?
  if ( roc < 0 ) {
    AliWarning(Form("Roc outside range (roc %d) !", roc));
    return -1;
  } else if ( roc < 36 ) {   // inner sector
    return false;
  } else if ( roc < 72 ) {   // outer sector, C side
    return true;
  } else {
    AliWarning(Form("Roc outside range (roc %d) !", roc));
    return -1;
  }
}

// EOF
