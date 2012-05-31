#ifndef AliTPCmapper_H
#define AliTPCmapper_H
/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//
// AliTPCmapper class
// Class for all mapping functions (hardware coordinates <-> pad coordinates)
// Author: Christian Lippmann
//       
/////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliTPCAltroMapping;

class AliTPCmapper : public TObject{

public:

  AliTPCmapper();
  AliTPCmapper(const char * dirname);
  virtual ~AliTPCmapper();

  AliTPCmapper& operator = (const AliTPCmapper& mapper);
  AliTPCmapper(const AliTPCmapper& mapper);

  void Init(const char * dirname);
  //
  AliTPCAltroMapping **GetAltroMapping() { return fMapping; };

  // ALTRO mapping functions
  Int_t GetPad(Int_t patch, Int_t hwAddress) const;
  Int_t GetPad(Int_t patch, Int_t branch, Int_t fec, Int_t chip, Int_t channel) const;
  Int_t GetPadRow(Int_t patch, Int_t hwAddress) const;
  Int_t GetPadRow(Int_t patch, Int_t branch, Int_t fec, Int_t chip, Int_t channel) const;

  // ALTRO mapping functions on roc level (padrow = 0 ... kNpadrowIROC, kNpadrowOROC)
  Int_t GetHWAddress(Int_t roc, Int_t padrow, Int_t pad) const;
  Int_t GetRcu(Int_t roc, Int_t padrow, Int_t pad) const;
  Int_t GetPatch(Int_t roc, Int_t padrow, Int_t pad) const;
  Int_t GetBranch(Int_t roc, Int_t padrow, Int_t pad) const;
  Int_t GetFEChw(Int_t roc, Int_t padrow, Int_t pad) const;
  Int_t GetFEC(Int_t roc, Int_t padrow, Int_t pad) const;
  Int_t GetChip(Int_t roc, Int_t padrow, Int_t pad) const;
  Int_t GetChannel(Int_t roc, Int_t padrow, Int_t pad) const;

  // ALTRO mapping functions on sector level (globalpadrow = 0 ... kNpadrow)
  Int_t GetGlobalPadRow(Int_t patch, Int_t hwAddress) const;
  Int_t GetGlobalPadRow(Int_t patch, Int_t branch, Int_t fec, Int_t chip, Int_t channel) const;
  Int_t GetHWAddressSector(Int_t globalpadrow, Int_t pad) const;
  Int_t GetRcuSector(Int_t globalpadrow, Int_t pad) const;
  Int_t GetPatchSector(Int_t globalpadrow, Int_t pad) const;
  Int_t GetBranchSector(Int_t globalpadrow, Int_t pad) const;
  Int_t GetFEChwSector(Int_t globalpadrow, Int_t pad) const;
  Int_t GetFECSector(Int_t globalpadrow, Int_t pad) const;
  Int_t GetChipSector(Int_t globalpadrow, Int_t pad) const;
  Int_t GetChannelSector(Int_t globalpadrow, Int_t pad) const;

  // Coding and decoding of hardware addresses
  Int_t CodeHWAddress(Int_t branch,  Int_t fec,  Int_t chip, Int_t channel) const;
  Int_t DecodedHWAddressBranch(Int_t hwAddress) const;
  Int_t DecodedHWAddressFECaddr(Int_t hwAddress) const;
  Int_t DecodedHWAddressChipaddr(Int_t hwAddress) const;
  Int_t DecodedHWAddressChanneladdr(Int_t hwAddress) const;

  // Pad Geometry on sector level (padrow = 0 ... kNpadrow)
  Int_t    GetNpads(Int_t roc, Int_t padrow) const;
  Int_t    GetNpads(Int_t globalpadrow) const;
  Int_t    GetNpadrows(Int_t roc) const;
  /*
  Double_t GetPadXlocal(Int_t globalpadrow) const;
  Double_t GetPadYlocal(Int_t globalpadrow, Int_t pad) const;
  Double_t GetPadXglobal(Int_t globalpadrow, Int_t pad, Int_t sector) const;
  Double_t GetPadYglobal(Int_t globalpadrow, Int_t pad, Int_t sector) const;
  Double_t GetPadWidth(Int_t globalpadrow) const;
  Double_t GetPadLength(Int_t globalpadrow) const;
  */

  // Conversion between hardware FEC numbering and official numbering
  Int_t HwToOffline(Int_t patch, Int_t branch, Int_t fec) const;
  Int_t OfflineToHwBranch(Int_t patch, Int_t fec) const;
  Int_t OfflineToHwFec(Int_t patch, Int_t fec) const;

  // More mapping functions
  Int_t GetEquipmentID(Int_t roc, Int_t padrow, Int_t pad) const;
  Int_t GetEquipmentIDsector(Int_t side, Int_t sector, Int_t globalpadrow, Int_t pad) const;
  Int_t GetEquipmentIDfromPatch(Int_t side, Int_t sector, Int_t patch) const;
  Int_t GetSectorFromRoc(Int_t roc) const;
  Int_t GetSideFromRoc(Int_t roc) const;
  Int_t GetRocFromPatch(Int_t side, Int_t sector, Int_t patch) const;
  Int_t GetRoc(Int_t side, Int_t sector, Int_t globalpadrow, Int_t pad) const;
  Int_t GetSideFromEquipmentID(Int_t equipmentID) const;
  Int_t GetSectorFromEquipmentID(Int_t equipmentID) const;
  Int_t GetRocFromEquipmentID(Int_t equipmentID) const;
  Int_t GetPatchFromEquipmentID(Int_t equipmentID) const;

  // Even more
  Int_t  GetNfec(Int_t patch, Int_t branch) const;
  Int_t  GetNfec(Int_t patch) const;
  Bool_t IsIROC(Int_t roc) const;
  Bool_t IsOROC(Int_t roc) const;
  
  Int_t  GetTpcDdlOffset() const {return fTpcDdlOffset;}
  Int_t  GetNumDdl() const {return fNside*fNsector*fNrcu; }

 private:

  Int_t fNside;        // TPC has 2 sides
  Int_t fNsector;      // TPC side has 18 sectors
  Int_t fNrcu;         // Sector has 6 RCUs (patches)
  Int_t fNbranch;      // RCU has 2 branches
  Int_t fNaltro;       // FEC has 8 ALTROs
  Int_t fNchannel;     // ALTRO has 16 channels
  Int_t fNpadrow;      // Sector has 159 padrows
  Int_t fNpadrowIROC;  // IROC has 63 padrows
  Int_t fNpadrowOROC;  // OROC has 96 padrows

  Int_t fTpcDdlOffset; // DDL offset for TPC

  AliTPCAltroMapping *fMapping[6];    // The ALTRO mapping for each patch (rcu)

  ClassDef(AliTPCmapper,2)

};

#endif
