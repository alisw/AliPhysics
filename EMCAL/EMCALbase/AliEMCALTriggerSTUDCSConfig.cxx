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

#include "AliEMCALTriggerSTUDCSConfig.h"

// ROOT system
#include "TClonesArray.h"
#include "TVector2.h"

// Standard libraries
#include <bitset>
#include <iomanip>
#include <iostream>
#include <sstream>

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerSTUDCSConfig) ;
/// \endcond

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerSTUDCSConfig::AliEMCALTriggerSTUTRUErrorCount) ;
/// \endcond

AliEMCALTriggerSTUDCSConfig::AliEMCALTriggerSTUDCSConfig() : TObject(),
fGetRawData(1),
fRegion(0xFFFFFFFF),
fFw(0x2A012),
fPatchSize(0),
fMedian(0)
{
  for (int i = 0; i < 3; i++) 
  {
    for (int j = 0; j < 2; j++) 
    {
      fG[i][j] = 0;
      fJ[i][j] = 0;
    }
  }
  
  memset(fPHOSScale, 0, sizeof(Int_t) * 4);
  memset(fTRUErrorCounts, 0, sizeof(TClonesArray *) * 68);
}

AliEMCALTriggerSTUDCSConfig::AliEMCALTriggerSTUDCSConfig(const AliEMCALTriggerSTUDCSConfig &obj) : TObject(),
fGetRawData(1),
fRegion(0xFFFFFFFF),
fFw(0x2A012),
fPatchSize(0),
fMedian(0)
{
  for (int i = 0; i < 3; i++) 
  {
    for (int j = 0; j < 2; j++) 
    {
      fG[i][j] = obj.GetG(i,j);
      fJ[i][j] = obj.GetJ(i,j);
    }
  }
  
  memset(fPHOSScale, 0, sizeof(Int_t) * 4);
  memset(fTRUErrorCounts, 0, sizeof(TClonesArray *) * 68);
  
  SetRawData(obj.GetRawData());
  SetRegion(obj.GetRegion());
  SetFw(obj.GetFw());
  for (int i = 0; i < 4; i++) {
    SetPHOSScale(i,obj.GetPHOSScale(i));
  }
  for (int i = 0; i < 68 ; i++) {
    TClonesArray * gTRUErrorCounts = obj.GetErrorCountsForTRU(i);
    if (!gTRUErrorCounts) continue;
    for (int j = 0; j < gTRUErrorCounts->GetEntries(); j++) {
      AliEMCALTriggerSTUTRUErrorCount * fErrorCount = (AliEMCALTriggerSTUTRUErrorCount *) gTRUErrorCounts->At(j);
      if (fErrorCount) {
        SetTRUErrorCounts(i,fErrorCount->GetTime(),fErrorCount->GetErrorCount());
      }
    }
  }

  SetPatchSize(obj.GetPatchSize());
  SetMedianMode(obj.GetMedianMode());
}

AliEMCALTriggerSTUDCSConfig::~AliEMCALTriggerSTUDCSConfig()
{
  for(int itru = 0; itru < 68; itru++)
  {
    if(fTRUErrorCounts[itru]) delete fTRUErrorCounts[itru];
  }
}

void AliEMCALTriggerSTUDCSConfig::GetSegmentation(TVector2& v1, TVector2& v2, TVector2& v3, TVector2& v4) const
{
  v1.Set(1., 1.);
  v2.Set(2., 2.);
  v3.Set(4., 4.);
  
  //Double_t js = 2 + (fFw >> 16); // Old method for getting patch size, valid in Run 1
  Double_t js = 2 + fPatchSize; // Patch Size = 0 or 2, from OCDB
  v4.Set(js, js);
}

void  AliEMCALTriggerSTUDCSConfig::SetTRUErrorCounts(Int_t itru, Int_t itime, ULong64_t errorcounts)
{
  if(itru > 67) return;
  
  if(!fTRUErrorCounts[itru])
    fTRUErrorCounts[itru] = new TClonesArray("AliEMCALTriggerSTUDCSConfig::AliEMCALTriggerSTUTRUErrorCount");
  
  AliEMCALTriggerSTUTRUErrorCount test(itime, errorcounts), *found(NULL);
  
  if((found = dynamic_cast<AliEMCALTriggerSTUTRUErrorCount *>(fTRUErrorCounts[itru]->FindObject(&test))))
  {
    found->SetValue(itime, errorcounts);
  } 
  else
  {
    Int_t nErrorCountsTRU = fTRUErrorCounts[itru]->GetEntries();
    new((*(fTRUErrorCounts[itru]))[nErrorCountsTRU]) AliEMCALTriggerSTUTRUErrorCount(itime, errorcounts);
  }
}

TClonesArray *AliEMCALTriggerSTUDCSConfig::GetErrorCountsForTRU(Int_t itru) const
{
  if(itru > 67) return NULL;
  
  return fTRUErrorCounts[itru];
}

///
/// Checks for equalness according to the time stamp.
//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTUDCSConfig::AliEMCALTriggerSTUTRUErrorCount::IsEqual(const TObject *o) const
{
  const AliEMCALTriggerSTUTRUErrorCount *test = dynamic_cast<const AliEMCALTriggerSTUTRUErrorCount *>(o);
  
  if(!test) return false;
  
  return test->fTime == fTime;
}

//
// Compare time-dependent error counts based on the time information.
//_____________________________________________________________________________
Int_t AliEMCALTriggerSTUDCSConfig::AliEMCALTriggerSTUTRUErrorCount::Compare(const TObject *o) const
{
  const AliEMCALTriggerSTUTRUErrorCount *test = dynamic_cast<const AliEMCALTriggerSTUTRUErrorCount *>(o);
  
  if(!test) return 1;
  
  if(fTime > test->fTime) return  1;
  if(fTime < test->fTime) return -1;
  
  return 0;
}

bool AliEMCALTriggerSTUDCSConfig::operator==(const AliEMCALTriggerSTUDCSConfig &other) const {
  return (fGetRawData == other.fGetRawData) && (fRegion == other.fRegion) &&
         (fFw == other.fFw) && (fPatchSize == other.fPatchSize) && (fMedian == other.fMedian) &&
         !memcmp(fPHOSScale, other.fPHOSScale, sizeof(Int_t) * 4) &&
         !memcmp(fG, other.fG, sizeof(Int_t) * 6) &&
         !memcmp(fJ, other.fJ, sizeof(Int_t) * 6);
}

std::ostream &operator<<(std::ostream &stream, const AliEMCALTriggerSTUDCSConfig &config){
  stream << "Gamma High: (" << config.fG[0][0] << ", " << config.fG[1][0] << ", " << config.fG[2][0] << ")" << std::endl;
  stream << "Gamma Low:  (" << config.fG[0][1] << ", " << config.fG[1][1] << ", " << config.fG[2][1] << ")" << std::endl;
  stream << "Jet High:   (" << config.fJ[0][0] << ", " << config.fJ[1][0] << ", " << config.fJ[2][0] << ")" << std::endl;
  stream << "Jet Low:    (" << config.fJ[0][1] << ", " << config.fJ[1][1] << ", " << config.fJ[2][1] << ")" << std::endl;
  stream << "GetRawData: " << config.fGetRawData 
         << ", Region: " << std::hex << config.fRegion << std::dec << "(" << std::bitset<sizeof(config.fRegion) * 8>(config.fRegion) << ")"
         << ", Median: " << config.fMedian 
         << ", Firmware: " << std::hex << config.fFw << std::dec 
         << ", PHOS Scale: (" << config.fPHOSScale[0] << ", " << config.fPHOSScale[1] << ", " << config.fPHOSScale[2] << ", " << config.fPHOSScale[3]
         << ")" << std::endl;
  return stream;
}

std::string AliEMCALTriggerSTUDCSConfig::ToJSON() const {
  std::stringstream jsonstring;
  jsonstring << "{" 
             << "\"fG\":[[" << fG[0][0] << "," << fG[1][0] << "," << fG[2][0] << "],[" << fG[0][1] << "," << fG[1][1] << "," << fG[2][1] <<"]],"
             << "\"fJ\":[[" << fJ[0][0] << "," << fJ[1][0] << "," << fJ[2][0] << "],[" << fJ[0][1] << "," << fJ[1][1] << "," << fJ[2][1] <<"]],"
             << "\"fRawData\":" << fGetRawData << ","
             << "\"fRegion\":" << fRegion << ","
             << "\"fFirmware\":" << fFw << ","
             << "\"fMedian\":" << fMedian << ","
             << "\"fPHOSScale\":[" << fPHOSScale[0] << "," << fPHOSScale[1] << "," << fPHOSScale[2] << "," << fPHOSScale[3] << "]"
             << "}";

  return jsonstring.str();
}
