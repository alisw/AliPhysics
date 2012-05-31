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
***************************************************************************/

/* 
 * author: Roberto Preghenella (R+), Roberto.Preghenella@bo.infn.it
 */

//////////////////////////////////////////////////////////////////////
//                                                                  //
//                                                                  //
//   This class provides the TOF FEE reader.                        //
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <TSystem.h>
#include "AliTOFFEEReader.h"
#include "AliTOFFEEConfig.h"
#include "AliTOFFEElightConfig.h"
#include "AliTOFRawStream.h"
#include "AliTOFGeometry.h"
#include "AliTOFcalibHisto.h"
#include "AliLog.h"
#include <fstream>

ClassImp(AliTOFFEEReader)

//_______________________________________________________________

AliTOFFEEReader::AliTOFFEEReader() :
  TObject(),
  fFEEConfig(new AliTOFFEEConfig()),
  fFEElightConfig(new AliTOFFEElightConfig()),
  fChannelEnabled(),
  fMatchingWindow(),
  fLatencyWindow()
{
  /* 
   * 
   * default constructor 
   *
   */

  Reset();
}

//_______________________________________________________________

AliTOFFEEReader::AliTOFFEEReader(const AliTOFFEEReader &source) :
  TObject(source),
  fFEEConfig(new AliTOFFEEConfig()),
  fFEElightConfig(new AliTOFFEElightConfig())
{
  /* 
   * 
   * copy constructor 
   *
   */

  Reset();
  memcpy(fFEEConfig, source.fFEEConfig, sizeof(AliTOFFEEConfig));
  memcpy(fFEElightConfig, source.fFEElightConfig, sizeof(AliTOFFEElightConfig));
}

//_______________________________________________________________

AliTOFFEEReader &
AliTOFFEEReader::operator=(const AliTOFFEEReader &source)
{
  /* 
   * 
   * operator = 
   * 
   */

  TObject::operator=(source);
  memcpy(fFEEConfig, source.fFEEConfig, sizeof(AliTOFFEEConfig));
  memcpy(fFEElightConfig, source.fFEElightConfig, sizeof(AliTOFFEElightConfig));
  return *this;
}

//_______________________________________________________________

AliTOFFEEReader::~AliTOFFEEReader()
{
  /* 
   *
   * default destructor 
   *
   */

  delete fFEEConfig;
  delete fFEElightConfig;
}

//_______________________________________________________________

void
AliTOFFEEReader::ResetChannelEnabledArray()
{
  /*
   *
   * reset channel enabled array
   *
   */

  for (Int_t iIndex = 0; iIndex < GetNumberOfIndexes(); iIndex++)
    fChannelEnabled[iIndex] = kFALSE;
}

//_______________________________________________________________

void
AliTOFFEEReader::ResetTriggerMaskArray()
{
  /*
   *
   * reset trigger mask array
   *
   */

  for (Int_t iddl = 0; iddl < GetNumberOfDDLs(); iddl++)
    fTriggerMask[iddl] = 0x0;
}

//_______________________________________________________________

void
AliTOFFEEReader::Reset()
{
  /*
   *
   * reset 
   *
   */

  for (Int_t iIndex = 0; iIndex < GetNumberOfIndexes(); iIndex++) {
    fChannelEnabled[iIndex] = kFALSE;
    fMatchingWindow[iIndex] = 0;
    fLatencyWindow[iIndex] = 0;
  }

  for (Int_t iddl = 0; iddl < GetNumberOfDDLs(); iddl++)
    fTriggerMask[iddl] = 0x0;
}

//_______________________________________________________________

void
AliTOFFEEReader::LoadFEEConfig(const Char_t *FileName) const
{
  /*
   *
   * load FEE config
   *
   */

  Char_t *expandedFileName = gSystem->ExpandPathName(FileName);
  std::ifstream is;
  is.open(expandedFileName, std::ios::binary);
  is.read((Char_t *)fFEEConfig, sizeof(AliTOFFEEConfig));
  is.close();
}

//_______________________________________________________________

void
AliTOFFEEReader::LoadFEElightConfig(const Char_t *FileName) const
{
  /*
   *
   * load FEElight config
   *
   */

  Char_t *expandedFileName = gSystem->ExpandPathName(FileName);
  std::ifstream is;
  is.open(expandedFileName, std::ios::binary);
  is.read((Char_t *)fFEElightConfig, sizeof(AliTOFFEElightConfig));
  is.close();
}

//_______________________________________________________________

Int_t
AliTOFFEEReader::ParseFEEConfig()
{
  /* 
   *
   * parse FEE config
   *
   * loops over all FEE channels, checks whether they are enabled
   * and sets channel enabled 
   *
   */

  AliInfo("parsing TOF FEE config");

  AliTOFRawStream rawStream;
  Int_t nEnabled = 0;
  Int_t volume[5], index;
  Int_t temp;

  Reset();

  /* loop over all FEE channels */
  for (Int_t iDDL = 0; iDDL < GetNumberOfDDLs(); iDDL++)
    for (Int_t iTRM = 0; iTRM < GetNumberOfTRMs(); iTRM++)
      for (Int_t iChain = 0; iChain < GetNumberOfChains(); iChain++)
	for (Int_t iTDC = 0; iTDC < GetNumberOfTDCs(); iTDC++)
	  for (Int_t iChannel = 0; iChannel < GetNumberOfChannels(); iChannel++)
	    /* check whether FEE channel is enabled */
	    if (IsChannelEnabled(iDDL, iTRM + 3, iChain, iTDC, iChannel)) {
	      /* convert FEE channel indexes into detector indexes */
	      rawStream.EquipmentId2VolumeId(iDDL, iTRM + 3, iChain, iTDC, iChannel, volume);
	      /* swap padx and padz to fit AliTOFGeometry::GetIndex behaviour */
	      temp = volume[4]; volume[4] = volume[3]; volume[3] = temp; 
	      /* check if index is ok */
	      if (volume[0] < 0 || volume[0] > 17 ||
		  volume[1] < 0 || volume[1] > 4 ||
		  volume[2] < 0 || volume[2] > 18 ||
		  volume[3] < 0 || volume[3] > 1 ||
		  volume[4] < 0 || volume[4] > 47)
		continue;
	      /* convert detector indexes into calibration index */
	      index = AliTOFGeometry::GetIndex(volume);
	      /* check calibration index */
	      if (index != -1 && index < GetNumberOfIndexes()) {
		/* set calibration channel enabled */
		fChannelEnabled[index] = kTRUE;
		fMatchingWindow[index] = GetMatchingWindow(iDDL, iTRM + 3, iChain, iTDC, iChannel);
		nEnabled++;
	      }
	    }
  return nEnabled;
}

//_______________________________________________________________

Int_t
AliTOFFEEReader::ParseFEElightConfig()
{
  /* 
   *
   * parse FEElight config
   *
   * loops over all FEE channels, checks whether they are enabled
   * and sets channel enabled 
   *
   */

  AliInfo("parsing TOF FEElight config");

  Reset();

  AliTOFcalibHisto calibHisto;
  calibHisto.LoadCalibHisto();

  Int_t nEnabled = 0, index;
  AliTOFFEEchannelConfig *channelConfig = NULL;
  for (Int_t i = 0; i < GetNumberOfIndexesEO(); i++) {
    channelConfig = fFEElightConfig->GetChannelConfig(i);
    if (!channelConfig->IsEnabled()) continue;
    /* get index DO from index EO */
    index = (Int_t)calibHisto.GetCalibMap(AliTOFcalibHisto::kIndex, i);
    if (index == -1) continue;
    nEnabled++;
    fChannelEnabled[index] = channelConfig->IsEnabled();
    fMatchingWindow[index] = channelConfig->GetMatchingWindow();
    fLatencyWindow[index] = channelConfig->GetLatencyWindow();
  }

  AliTOFFEEtriggerConfig *triggerConfig = NULL;
  for (Int_t iddl = 0; iddl < GetNumberOfDDLs(); iddl++) {
    triggerConfig = fFEElightConfig->GetTriggerConfig(iddl);
    fTriggerMask[iddl] = triggerConfig->GetStatusMap();
  }
 
  return nEnabled;
}

//_______________________________________________________________

Bool_t 
AliTOFFEEReader::IsChannelEnabled(Int_t iDDL, Int_t iTRM, Int_t iChain, Int_t iTDC, Int_t iChannel) const
{
  /*
   *
   * is channel enabled
   *
   * checks whether a FEE channel is enabled using the
   * TOF FEE config object.
   *
   */

  AliTOFFEEConfig *feeConfig;
  AliTOFCrateConfig *crateConfig;
  AliTOFTRMConfig *trmConfig;
  Int_t maskPB, maskTDC, activeChip;
  
  /* get and check fee config */
  if (!(feeConfig = GetFEEConfig()))
    return kFALSE;
  
  /* get and check crate config */
  if (!(crateConfig = feeConfig->GetCrateConfig(iDDL)))
    return kFALSE;
  
  /* get and check TRM config */
  if (!(trmConfig = crateConfig->GetTRMConfig(iTRM - 3)))
    return kFALSE;

  /* check DRM enabled */
  if (!crateConfig->IsDRMEnabled())
    return kFALSE;

  /* check TRM enabled */
  if (!crateConfig->IsTRMEnabled(iTRM - 3))
    return kFALSE;

  /* switch chain */
  switch (iChain) {
    /* chain A */
  case 0:
    /* check chain enabled */
    if (trmConfig->GetChainAFlag() != 1)
      return kFALSE;
    /* get active chip mask */
    activeChip = trmConfig->GetActiveChipA();
    /* switch TDC */
    switch (iTDC) {
    case 0: case 1: case 2:
      maskPB = trmConfig->GetMaskPB0();
      break;
    case 3: case 4: case 5:
      maskPB = trmConfig->GetMaskPB1();
      break;
    case 6: case 7: case 8:
      maskPB = trmConfig->GetMaskPB2();
      break;
    case 9: case 10: case 11:
      maskPB = trmConfig->GetMaskPB3();
      break;
    case 12: case 13: case 14:
      maskPB = trmConfig->GetMaskPB4();
      break;
    default:
      return kFALSE;
      break;  
    } /* switch TDC */
    break; /* chain A */
    /* chain B */
  case 1:
    /* check chain enabled */
    if (trmConfig->GetChainBFlag() != 1)
      return kFALSE;
    /* get active chip mask */
    activeChip = trmConfig->GetActiveChipB();
    /* switch TDC */
    switch (iTDC) {
    case 0: case 1: case 2:
      maskPB = trmConfig->GetMaskPB5();
      break;
    case 3: case 4: case 5:
      maskPB = trmConfig->GetMaskPB6();
      break;
    case 6: case 7: case 8:
      maskPB = trmConfig->GetMaskPB7();
      break;
    case 9: case 10: case 11:
      maskPB = trmConfig->GetMaskPB8();
      break;
    case 12: case 13: case 14:
      maskPB = trmConfig->GetMaskPB9();
      break;
    default:
      return kFALSE;
      break;  
    } /* switch TDC */
    break; /* chain B */
  default:
    return kFALSE;
    break;
  } /* switch chain */

  /* check chip enabled */
  if (!(activeChip & (0x1 << iTDC)))
    return kFALSE;

  /* check channel enabled */
  maskTDC = (maskPB & (0xFF << ((iTDC % 3) * 8))) >> ((iTDC % 3) * 8);
  if (maskTDC & (0x1 << iChannel))
    return kTRUE;
  else
    return kFALSE;
  
}

//_______________________________________________________________

Int_t 
AliTOFFEEReader::GetMatchingWindow(Int_t iDDL, Int_t iTRM, Int_t, Int_t, Int_t) const
{
  /*
   *
   * get matching window
   *
   * checks whether a FEE channel is enabled using the
   * TOF FEE config object and return the associated
   * matching window
   *
   */
  
  AliTOFFEEConfig *feeConfig;
  AliTOFCrateConfig *crateConfig;
  AliTOFTRMConfig *trmConfig;

  /* get and check fee config */
  if (!(feeConfig = GetFEEConfig()))
    return 0;
  
  /* get and check crate config */
  if (!(crateConfig = feeConfig->GetCrateConfig(iDDL)))
    return 0;
  
  /* get and check TRM config */
  if (!(trmConfig = crateConfig->GetTRMConfig(iTRM - 3)))
    return 0;

  /* check DRM enabled */
  if (!crateConfig->IsDRMEnabled())
    return 0;

  /* check TRM enabled */
  if (!crateConfig->IsTRMEnabled(iTRM - 3))
    return 0;

  return trmConfig->GetMatchingWindow();
}


void
AliTOFFEEReader::DumpFEEConfig()
{
  /*
   * 
   * dump FEE config
   *
   */

  AliTOFFEEConfig *feeConfig = GetFEEConfig();
  AliTOFCrateConfig *crateConfig;
  AliTOFDRMConfig *drmConfig;
  AliTOFLTMConfig *ltmConfig;
  AliTOFTRMConfig *trmConfig;

  AliInfo("-------------------------------------");
  AliInfo("dumping TOF FEE config");
  AliInfo("-------------------------------------");
  AliInfo(Form("version: %d", feeConfig->GetVersion()));
  AliInfo(Form("dump time: %d", (Int_t)feeConfig->GetDumpTime()));
  AliInfo(Form("run number: %d", feeConfig->GetRunNumber()));
  AliInfo(Form("run type: %d", feeConfig->GetRunType()));
  AliInfo("-------------------------------------");
  
  /* loop over crates */
  for (Int_t iCrate = 0; iCrate < AliTOFFEEConfig::GetNumberOfCrates(); iCrate++) {
    crateConfig = feeConfig->GetCrateConfig(iCrate);
    
    /* check crate config */
    if (!crateConfig)
      continue;
    
    /* check DRM enabled */
    if (!crateConfig->IsDRMEnabled())
    continue;

    AliInfo(Form("crate id: %02d", iCrate));

    /* dump DRM config */
    drmConfig = crateConfig->GetDRMConfig();
    AliInfo(Form("DRM is enabled: drmId=%d, slotMask=%03x", drmConfig->GetDRMId(), drmConfig->GetSlotMask()));

    /* dump LTM config if enabled */
    if (crateConfig->IsLTMEnabled()) {
      ltmConfig = crateConfig->GetLTMConfig();
      AliInfo(Form("LTM is enabled: threshold=%d", ltmConfig->GetThreshold()));
    }
    
    /* dump CPDM config if enabled */
    if (crateConfig->IsCPDMEnabled()) {
      AliInfo(Form("CPDM is enabled"));
    }
    
    /* loop over TRMs */
    for (Int_t iTRM = 0; iTRM < AliTOFCrateConfig::GetNumberOfTRMs(); iTRM++) {

      trmConfig = crateConfig->GetTRMConfig(iTRM);

      /* check TRM config */
      if (!trmConfig)
	continue;
      
      /* check TRM enabled */
      if (!crateConfig->IsTRMEnabled(iTRM))
	continue;

      /* dump TRM config */
      AliInfo(Form("TRM%02d is enabled: matchWin=%d, latWin=%d, packFlag=%d", iTRM + 3, trmConfig->GetMatchingWindow(), trmConfig->GetLatencyWindow(), trmConfig->GetPackingFlag()));
      
      /* check TRM chain A flag */
      if (trmConfig->GetChainAFlag() == 1) {
	AliInfo(Form("TRM%02d chainA is enabled: activeChip=%04X, PB0=%06X, PB1=%06X, PB2=%06X, PB3=%06X, PB4=%06X", iTRM + 3, trmConfig->GetActiveChipA(), trmConfig->GetMaskPB0(), trmConfig->GetMaskPB1(), trmConfig->GetMaskPB2(), trmConfig->GetMaskPB3(), trmConfig->GetMaskPB4()));
      }

      /* check TRM chain B flag */
      if (trmConfig->GetChainBFlag() == 1) {
	AliInfo(Form("TRM%02d chainB is enabled: activeChip=%04X, PB5=%06X, PB6=%06X, PB7=%06X, PB8=%06X, PB9=%06X", iTRM + 3, trmConfig->GetActiveChipB(), trmConfig->GetMaskPB5(), trmConfig->GetMaskPB6(), trmConfig->GetMaskPB7(), trmConfig->GetMaskPB8(), trmConfig->GetMaskPB9()));
      }
      

      
    } /* loop over TRMs */
    AliInfo("-------------------------------------");
  } /* loop over crates */
}

//_______________________________________________________________

void
AliTOFFEEReader::CreateFEElightConfig(const Char_t *filename)
{
  /*
   *
   * create FEElight config 
   *
   */

  AliTOFFEElightConfig lightConfig;

  for (Int_t i = 0; i < GetNumberOfIndexes(); i++) {
    if (fChannelEnabled[i]) {
      lightConfig.GetChannelConfig(i)->SetStatus(AliTOFFEEchannelConfig::kStatusEnabled);
      lightConfig.GetChannelConfig(i)->SetMatchingWindow(fMatchingWindow[i]);
      lightConfig.GetChannelConfig(i)->SetLatencyWindow(fLatencyWindow[i]);
    }
    else {
      lightConfig.GetChannelConfig(i)->SetStatus(0x0);
      lightConfig.GetChannelConfig(i)->SetMatchingWindow(0);
      lightConfig.GetChannelConfig(i)->SetLatencyWindow(0);
    }
  }

  Char_t *expandedFileName = gSystem->ExpandPathName(filename);
  std::ofstream os;
  os.open(expandedFileName, std::ios::binary);
  os.write((Char_t *)&lightConfig, sizeof(AliTOFFEElightConfig));
  os.close();
  
}
