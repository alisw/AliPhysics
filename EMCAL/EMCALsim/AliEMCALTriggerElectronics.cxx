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

#include <vector>
#include "AliEMCALTriggerElectronics.h"
#include "AliEMCALTriggerTRU.h"
#include "AliEMCALTriggerSTU.h"
#include "AliEMCALGeometry.h"
#include "AliRunLoader.h"
#include "AliCDBManager.h"
#include "AliEMCAL.h" 
#include "AliRun.h" 
#include "AliEMCALTriggerDCSConfig.h"
#include "AliEMCALTriggerData.h"
#include "AliEMCALDigit.h"
#include "AliCaloRawStreamV3.h"
#include "AliEMCALTriggerSTURawStream.h"
#include "AliEMCALDigit.h"
#include "AliEMCALTriggerRawDigit.h"
#include "AliEMCALTriggerPatch.h"
#include "AliEMCALTriggerTRUDCSConfig.h"
#include "AliEMCALTriggerSTUDCSConfig.h"

#include <TVector2.h>

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerElectronics) ;
/// \endcond

///
/// Constructor
//__________________
AliEMCALTriggerElectronics::AliEMCALTriggerElectronics(const AliEMCALTriggerDCSConfig *dcsConf) : TObject(),
fTRU(new TClonesArray("AliEMCALTriggerTRU",52)),
fSTU(0x0),
fGeometry(0)
{
  TVector2 rSize;
	
  rSize.Set(4.,24.);	

  AliRunLoader *rl = AliRunLoader::Instance();
  if (rl->GetAliRun() && rl->GetAliRun()->GetDetector("EMCAL")) {
    AliEMCAL *emcal = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"));
    if (emcal) fGeometry = emcal->GetGeometry();
  }
  
  if (!fGeometry) {	
    fGeometry =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
    AliError("Cannot access geometry, create a new default one!");
  }
    
  fNTRU = fGeometry->GetNTotalTRU() ;

  AliDebug(999,TString::Format("fGeometry <<%s>> has %d TRUs; dcsConf has %d TRUs.\n",fGeometry->GetName(),fNTRU,dcsConf->GetTRUArr()->GetSize()));
  TString fGeometryName = fGeometry->GetName();

  // Update here for future trigger mapping versions
  Int_t iTriggerMapping = 1 + (Int_t) fGeometryName.Contains("DCAL");
  switch (iTriggerMapping ) {
    case 1:
    AliDebug(999,"Trigger Mapping V1.  Will not include DCAL.");
    break;
    case 2:
    AliDebug(999,"Trigger Mapping V2.  Will include DCAL.");
  }


  // STU for EMCAL
//	rSize.Set( 48., 64. );
  rSize.Set( 64., 48. );
  AliEMCALTriggerSTUDCSConfig* stuConf = dcsConf->GetSTUDCSConfig();
  fSTU = new AliEMCALTriggerSTU(stuConf, rSize);
	
  AliCDBManager * cdb = AliCDBManager::Instance();
  Int_t runNumber = cdb->GetRun();

  // Checking Firmware from DCS config to choose algorithm
  fMedianMode = stuConf->GetMedianMode();
  AliEMCALTriggerSTUDCSConfig* stuConfDCal = 0;

  if (iTriggerMapping == 1) {
    // In Run 1, the patch size was stored in the FW version number, not OCDB
    Int_t fFW = stuConf->GetFw();
    if ((fFW >> 16) == 2) {
      stuConf->SetPatchSize(2);
      AliDebug(999,TString::Format("Patch size has been set to 16x16 (value 2) based on the EMCALSTU firmware version. This should only occur for Run 1"));
    }
  }

  if (iTriggerMapping >= 2) {
    rSize.Set( 40., 48. );  // This should be accurate
    stuConfDCal = dcsConf->GetSTUDCSConfig(true);
    if (stuConfDCal) {
      fSTUDCAL = new AliEMCALTriggerSTU(stuConfDCal, rSize);
      AliDebug(999,TString::Format("Found DCAL STU firmware %x.",stuConfDCal->GetFw()));
    } else {
      AliError("No DCS Config found for DCAL, but Run 2 Geometry requested!  Cloning EMCAL STU config.");
      stuConfDCal = new AliEMCALTriggerSTUDCSConfig(*stuConf);
      AliDebug(999,"Manually setting DCAL STU fW object to 0xd007.");
      // Manually setting DCAL STU DCS fW version to 0xd007 (Don't count blank PHOS patches)
      stuConfDCal->SetFw(0xd007);
      // Patchsize for DCAL should be 0 (8x8 jet patch)
      stuConfDCal->SetPatchSize(0);

      fSTUDCAL = new AliEMCALTriggerSTU(stuConfDCal, rSize);
    }
  } else fSTUDCAL = 0;

  // active status for STUs
  if(!stuConf->GetRegion()) fSTU->SetActive(false);
  if(stuConfDCal && fSTUDCAL && !stuConfDCal->GetRegion()) fSTUDCAL->SetActive(false);

  // Checking if the L1 thresholds are missing
  for (int ithr = 0; ithr < 2; ithr++) {
    if (stuConf->GetG(2,ithr) == 0) {
      if(stuConf->GetRegion())
        AliError(Form("EMCAL DCS STU has 0 threshold for EGA %d.  This trigger will not be properly simulated.  Check OCDB!!!",ithr));
      else
        AliWarning(Form("EMCAL DCS STU has 0 threshold for EGA %d. EMCAL STU inactive - this trigger will not be simulated.", ithr));
    }
    if (stuConf->GetJ(2,ithr) == 0) {
      if(stuConf->GetRegion())
        AliError(Form("EMCAL DCS STU has 0 threshold for EJE %d.  This trigger will not be properly simulated.  Check OCDB!!!",ithr));
      else
        AliWarning(Form("EMCAL DCS STU has 0 threshold for EJE %d. EMCAL STU inactive - this trigger will not be simulated.", ithr));
    }

    if (stuConfDCal) {
      if (stuConfDCal->GetG(2,ithr) == 0) {
        if(stuConfDCal->GetRegion())
          AliError(Form("DCAL DCS STU has 0 threshold for EGA %d.  This trigger will not be properly simulated.  Check OCDB!!!",ithr));
        else
          AliWarning(Form("DCAL DCS STU has 0 threshold for EGA %d. DCAL STU inactive - this trigger will not be simulated.", ithr));
      }
      if (stuConfDCal->GetG(2,ithr) == 0) {
        if(stuConfDCal->GetRegion())
          AliError(Form("DCAL DCS STU has 0 threshold for EJE %d.  This trigger will not be properly simulated.  Check OCDB!!!",ithr));
        else
          AliWarning(Form("DCAL DCS STU has 0 threshold for EJE %d. DCAL STU inactive - this trigger will not be simulated.", ithr));
      }
    }
  }


  // Tests
  Int_t iFastORIndex = -1;
  Int_t iSM = -1;
  Int_t iEta = -1;
  Int_t iPhi = -1;

  TString str = "map";

  // 32 TRUs
  Int_t iTRU = 0; // iTRU is an index skipping fake TRUs
  for (Int_t i = 0; i < fNTRU; i++) { // 'i' is the geometric index (0-51)

    // Skip virtual TRUs in PHOS hole
    if (i == 34 || i == 35 ||
        i == 40 || i == 41 ||
        i == 46 || i == 47) {
      continue;
    }
    AliEMCALTriggerTRUDCSConfig *truConf = dcsConf->GetTRUDCSConfig(fGeometry->GetOnlineIndexFromTRUIndex(iTRU));


    Int_t iADC = 0 ; // any iADC will give the right SM
    fGeometry->GetAbsFastORIndexFromTRU(i,iADC,iFastORIndex);
    fGeometry->GetPositionInSMFromAbsFastORIndex(iFastORIndex, iSM, iEta, iPhi);
    Int_t tmpTRU = 0;
    Int_t tmpADC = 0;
    fGeometry->GetTRUFromAbsFastORIndex(iFastORIndex, tmpTRU , tmpADC);

    Int_t iSMType = fGeometry->GetSMType(iSM);
    
    if (truConf) { 		
      switch (iTriggerMapping) {
        case 1:
          rSize.Set(4.,24.);	
          break;
        case 2:
          switch (iSMType) {
            case AliEMCALGeometry::kEMCAL_Standard:
              rSize.Set(12.,8.);
              break;
            case AliEMCALGeometry::kDCAL_Standard:
              rSize.Set(12.,8.);
              break;
            case AliEMCALGeometry::kEMCAL_3rd:
              rSize.Set(4.,24.);
              break;
            case AliEMCALGeometry::kDCAL_Ext:
              rSize.Set(4.,24.);
              break;
            default:
            AliError("TRU Initialization for EMCAL half modules not implemented!");
          }
          break;
        default:
          AliError("TRU Initialization for this Trigger Mapping is not implemented!");
      }
      // AliEMCALTriggerTRU initialization needs Map Type 0 for ASide and 1 for CSide ADC mappings
      Int_t iMapType = 0;

      if (iTriggerMapping == 1) iMapType = iTRU % 2; // Run 1 Geometry (TRUs alternate sides)
      else if (iTriggerMapping == 2) {
        iMapType = iTRU % 2;
        if (i < 30) iMapType = (i % 6) / 3;
        else if (i < 32) iMapType = i % 2;
        else if (i < 50) iMapType = ((i - 2) % 6) / 3;
        else if (i < 52) iMapType = i % 2;
      } else {
        AliFatal("Unknown Trigger Mapping requested. May need implementation.");
      }

      AliEMCALTriggerTRU *mytru = new ((*fTRU)[i]) AliEMCALTriggerTRU(truConf, rSize, iMapType);
      //AliEMCALTriggerTRU *mytru = new ((*fTRU)[i]) AliEMCALTriggerTRU(truConf, rSize, iTRU % 2); // Run 1

      AliDebug(999,Form("Building TRU %d with dimensions %d x %d and MapType %d\n",iTRU,int(rSize.X()),int(rSize.Y()),iMapType));
      if(dcsConf->IsTRUEnabled(fGeometry->GetOnlineIndexFromTRUIndex(iTRU))) {
        if (truConf->GetGTHRL0() <= 1) { // Checking for the null L0 threshold
          AliError(Form("Active TRU %d DCS config is missing L0 threshold.  L0 Trigger will not be simulated for this TRU.",iTRU));
          mytru->SetActive(false);
        } 
      } else {
        AliDebug(999,Form("Inactive TRU %d - L0 trigger will not be simulated for this TRU", iTRU));
        mytru->SetActive(false);
      }

      AliEMCALTriggerTRU *oTRU = static_cast<AliEMCALTriggerTRU*>(fTRU->At(i));
      switch (iSMType) {
        case AliEMCALGeometry::kEMCAL_Standard:
        case AliEMCALGeometry::kEMCAL_Half:
        case AliEMCALGeometry::kEMCAL_3rd:
          // Build in the EMCAL STU
          fSTU->Build(str,i,oTRU->Map(),oTRU->RegionSize(),iTriggerMapping); 
          break;
        case AliEMCALGeometry::kDCAL_Standard:
        case AliEMCALGeometry::kDCAL_Ext:
          // Build in the DCAL STU
          fSTUDCAL->Build(str,i,oTRU->Map(),oTRU->RegionSize(),iTriggerMapping);
          break;
        default:
        AliError("Unrecognized geometry for TRU at STU building.  Update EMCAL/EMCALsim/AliEMCALTriggerElectronics.cxx.");
      }
    }
    else AliError(Form("TRU config not found in OCDB for iTRU = %d\n",iTRU)); 

    iTRU++; // real TRU index
  }
}

///
/// Destructor
//________________
AliEMCALTriggerElectronics::~AliEMCALTriggerElectronics()
{	
  fTRU->Delete();
 	delete fSTU;
	if (fSTUDCAL) delete fSTUDCAL;
}

///
/// Digits to trigger
//__________________
void AliEMCALTriggerElectronics::Digits2Trigger(TClonesArray* digits, const Int_t V0M[], AliEMCALTriggerData* data)
{
  // Digits to trigger
  Int_t pos, px, py, id; 

  Int_t region[104][48], posMap[104][48]; // [phi] [eta]

  for (Int_t i = 0; i < 104; i++) for (Int_t j = 0; j < 48; j++) 
  {
    region[i][j] =  0;
    posMap[i][j] = -1;
  }
  
  // Loading raw digits
  for (Int_t i = 0; i < digits->GetEntriesFast(); i++)
  {
    AliEMCALTriggerRawDigit* digit = (AliEMCALTriggerRawDigit*)digits->At(i);
    
    id = digit->GetId();
    
    Int_t iTRU, iADC;
    
    Bool_t isOK1 = fGeometry->GetTRUFromAbsFastORIndex(id, iTRU, iADC);
    
    if (!isOK1 || iTRU > fNTRU-1) continue;

    for (Int_t j = 0; j < digit->GetNSamples(); j++)
    {
      Int_t time, amp;
      Bool_t isOK2 = digit->GetTimeSample(j, time, amp);
      
      if (isOK1 && isOK2 && amp) {
        AliDebug(999, Form("=== TRU# %2d ADC# %2d time# %2d signal %d ===", iTRU, iADC, time, amp));

        AliEMCALTriggerTRU * etr = (static_cast<AliEMCALTriggerTRU*>(fTRU->At(iTRU)));
        if (etr) {
          if (data->GetMode())
            etr->SetADC(iADC, time, 4 * amp);
          else
            etr->SetADC(iADC, time,     amp);
        }
      }
    }
    
  //	if (fGeometry->GetPositionInEMCALFromAbsFastORIndex(id, px, py)) posMap[px][py] = i;
    if (fGeometry->GetPositionInEMCALFromAbsFastORIndex(id, py, px)) posMap[px][py] = i;
  }

  Int_t iL0 = 0;

  std::vector<int> timeL0(fNTRU);
  Int_t  timeL0min = 999;
  
  for (Int_t truid = 0; truid < fNTRU; truid++) 
  {
    AliEMCALTriggerTRU *currenttru = static_cast<AliEMCALTriggerTRU*>(fTRU->At(truid));
    if (!currenttru) continue;
    
    AliDebug(999, Form("===========< TRU %2d >============", truid));
    
    if (currenttru->L0()) // L0 recomputation: *ALWAYS* done from FALTRO
    {
        iL0++;
      
        timeL0[truid] = currenttru->GetL0Time();
      
        if (!timeL0[truid]) AliWarning(Form("TRU# %d has 0 trigger time",truid));
      
        if (timeL0[truid] < timeL0min) timeL0min = timeL0[truid];
      
        data->SetL0Trigger(AliEMCALTriggerData::L0DecisionOrigin_t::kOriginRecalc, truid, AliEMCALTriggerData::L0TriggerStatus_t::kFired); // TRU# i has issued a L0
    }
    else
      data->SetL0Trigger(AliEMCALTriggerData::L0DecisionOrigin_t::kOriginRecalc, truid, AliEMCALTriggerData::L0TriggerStatus_t::kNotFired);
  }


  AliDebug(999, Form("=== %2d TRU (out of %2d) has issued a L0 / Min L0 time: %d", iL0, fTRU->GetEntriesFast(), timeL0min));
  AliEMCALTriggerRawDigit* dig = 0x0;
  

  if (iL0 && (data->IsSimulationMode() || !fSTU->GetDCSConfig()->GetRawData())) 
  {
    // Update digits after L0 calculation
    for (Int_t i = 0; i < fNTRU; i++)
    {
      AliEMCALTriggerTRU *iTRU = static_cast<AliEMCALTriggerTRU*>(fTRU->At(i));
      if (!iTRU) continue;
      
//			Int_t reg[24][4];
//			for (int j = 0; j < 24; j++) for (int k = 0; k < 4; k++) reg[j][k] = 0;

      TVector2 * rSize = iTRU->RegionSize();
      Int_t nPhi = (Int_t) rSize->X();
      Int_t nEta = (Int_t) rSize->Y();


      Int_t  **reg = (Int_t **) malloc(nPhi * sizeof(Int_t *));
      if (!reg) {
        AliFatal("Failed to allocate memory for region.");
      }
      for (int j = 0; j < nPhi; j++) { 
        reg[j] = (Int_t *) malloc (nEta * sizeof(Int_t));
        if (!reg[j]) AliFatal("Failed to allocate memory for subregion.");
      }
    
      for (int j = 0; j < nPhi; j++) for (int k = 0; k < nEta; k++) reg[j][k] = 0;
//			for (int j = 0; j < nEta; j++) for (int k = 0; k < nPhi; k++) reg[j][k] = 0;
          
      // Loading information from TRU into matrix reg
      iTRU->GetL0Region(timeL0min, reg);

      for (int j = 0; j < iTRU->RegionSize()->X(); j++) // phi
      {  
        for (int k = 0; k < iTRU->RegionSize()->Y(); k++) // eta
        {
          if (reg[j][k]  // (phi,eta)
//					if (true 
            && 
            fGeometry->GetAbsFastORIndexFromPositionInTRU(i, k, j, id)  // (iTRU, Eta, Phi, FastOR)
//						fGeometry->GetAbsFastORIndexFromPositionInTRU(i, j, k, id)
            &&
            fGeometry->GetPositionInEMCALFromAbsFastORIndex(id, py, px)) // (FastOR, Eta, Phi)
  //					fGeometry->GetPositionInEMCALFromAbsFastORIndex(id, px, py))
          {
            pos = posMap[px][py];
                  
            if (pos == -1)
            {
              // Add a new digit
              posMap[px][py] = digits->GetEntriesFast();

              new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
                    
              dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);							
            }
            else
            {
              dig = (AliEMCALTriggerRawDigit*)digits->At(pos);
            }
            
            // 14b to 12b STU time sums
            reg[j][k] >>= 2; 
            
            dig->SetL1TimeSum(reg[j][k]);
          }
        }
      }
    }
  }

  if (iL0 && data->IsSimulationMode()) //Simulation
  {
    for (Int_t i = 0; i < fNTRU; i++)
    {
      AliEMCALTriggerTRU *iTRU = static_cast<AliEMCALTriggerTRU*>(fTRU->At(i));
      if (!iTRU) continue;
      
      AliDebug(999, Form("=== TRU# %2d found %d patches", i, (iTRU->Patches()).GetEntriesFast()));
      
      TIter next(&iTRU->Patches());
      while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)next())
      {
//        p->Position(px, py); // x is phi, y is eta
        p->Position(py, px); // x is phi, y is eta
      
//				if (fGeometry->GetAbsFastORIndexFromPositionInTRU(i, px, py, id) 
//					&& 
//					fGeometry->GetPositionInEMCALFromAbsFastORIndex(id, px, py)) p->SetPosition(px, py);

        // Local 2 Global
//        if (fGeometry->GetAbsFastORIndexFromPositionInTRU(i, py, px, id) 
        if (fGeometry->GetAbsFastORIndexFromPositionInTRU(i, px, py, id) 
          && 
          fGeometry->GetPositionInEMCALFromAbsFastORIndex(id, py, px)) p->SetPosition(px, py);

//				if (fGeometry->GetAbsFastORIndexFromPositionInTRU(i, px, py, id)  //These patches used (x,y) = eta,phi
//					&& 
//					fGeometry->GetPositionInEMCALFromAbsFastORIndex(id, py, px)) p->SetPosition(px, py); // Now using (x,y) = (phi,eta)


        if (AliDebugLevel()) p->Print("");
        
        Int_t peaks = p->Peaks();
                
//        Int_t sizeX = (Int_t) ((iTRU->PatchSize())->X() * (iTRU->SubRegionSize())->X());
//        Int_t sizeY = (Int_t) ((iTRU->PatchSize())->Y() * (iTRU->SubRegionSize())->Y());
        Int_t sizeX = (Int_t) ((iTRU->PatchSize())->Y() * (iTRU->SubRegionSize())->Y());
        Int_t sizeY = (Int_t) ((iTRU->PatchSize())->X() * (iTRU->SubRegionSize())->X());
 
        for (Int_t j = 0; j < sizeX * sizeY; j++)
        {
          if (peaks & (1 << j))
          {
//            pos = posMap[px + j % sizeX][py + j / sizeX];
            pos = posMap[px + j / sizeX][py + j % sizeX];
              
            if (pos == -1) {
              // Add a new digit
//              posMap[px + j % sizeX][py + j / sizeX] = digits->GetEntriesFast();
              posMap[px + j / sizeX][py + j % sizeX] = digits->GetEntriesFast();

              new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
                
              dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);							
            } else {
              dig = (AliEMCALTriggerRawDigit*)digits->At(pos);
            }

            if (dig) dig->SetL0Time(p->Time());
          }
        }
            
        pos = posMap[px][py];
          
        if (pos == -1) {
          // Add a new digit
          posMap[px][py] = digits->GetEntriesFast();

          new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
            
          dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);
          
        }
        else
        {
          dig = (AliEMCALTriggerRawDigit*)digits->At(pos);
        }
          
        if (dig) dig->SetTriggerBit(kL0, 0);
      }
    }
  }

  // Prepare STU for L1 calculation

	Int_t nSTURegionSizeX = fSTU->RegionSize()->X();
	if (fSTUDCAL) nSTURegionSizeX += fSTUDCAL->RegionSize()->X();
  Int_t nSTURegionSizeY = fSTU->RegionSize()->Y();

  for (int i = 0; i < nSTURegionSizeX; i++) {
    for (int j = 0; j < nSTURegionSizeY; j++) {
      
      pos = posMap[i][j];

      if (pos >= 0) {
        
        AliEMCALTriggerRawDigit *digit = (AliEMCALTriggerRawDigit*)digits->At(pos);
    
        if (digit && digit->GetL1TimeSum() > -1) region[i][j] = digit->GetL1TimeSum();
      }
    }
  }

  AliDebug(999,"==================== STU  ====================");
  if (AliDebugLevel() >= 999) {
    fSTU->Scan();
    if (fSTUDCAL) fSTUDCAL->Scan();
  }
  AliDebug(999,"==============================================");
	
	fSTU->SetRegion(region);
	if (fSTUDCAL) fSTUDCAL->SetRegion(&region[64]); //DCAL's region is subset of region
//	if (fSTUDCAL) fSTUDCAL->SetRegion(region); 
  
  if (data->IsRawDataMode()) // Reconstruction
  {
    for (int ithr = 0; ithr < 2; ithr++) {
      AliDebug(999, Form(" THR %d / EGA %d / EJE %d", ithr, data->GetL1GammaThreshold(ithr), data->GetL1JetThreshold(ithr)));
                 
      fSTU->SetThreshold(kL1GammaHigh + ithr, data->GetL1GammaThreshold(ithr));
      fSTU->SetThreshold(kL1JetHigh + ithr, data->GetL1JetThreshold(  ithr));
      if (fSTUDCAL) {
        fSTUDCAL->SetThreshold(kL1GammaHigh + ithr, data->GetL1GammaThreshold(ithr));
        fSTUDCAL->SetThreshold(kL1JetHigh + ithr, data->GetL1JetThreshold(  ithr));
      }
    }
  }
  else  // Simulation
  {
    for (int ithr = 0; ithr < 2; ithr++) {
      //
      fSTU->ComputeThFromV0(kL1GammaHigh + ithr, V0M); 
    // Why is this commented out?
//      data->SetL1GammaThreshold(ithr, fSTU->GetThreshold(kL1GammaHigh + ithr));
      
      fSTU->ComputeThFromV0(kL1JetHigh + ithr,   V0M);
  //    data->SetL1JetThreshold(ithr, fSTU->GetThreshold(kL1JetHigh + ithr)  );
      

      AliCDBManager * cdb = AliCDBManager::Instance();
      Int_t runNumber = cdb->GetRun();

      TString fGeometryName = fGeometry->GetName();
      Int_t iTriggerMapping = 1 + (Int_t) fGeometryName.Contains("DCAL");


      AliDebug(999, Form("STU EMCAL THR %d EGA %d EJE %d", ithr, fSTU->GetThreshold(kL1GammaHigh + ithr), fSTU->GetThreshold(kL1JetHigh + ithr)));
      if (fSTUDCAL) {
        // Setting threshold for DCAL.
        fSTUDCAL->ComputeThFromV0(kL1GammaHigh + ithr, V0M); 
        fSTUDCAL->ComputeThFromV0(kL1JetHigh + ithr,   V0M);

       AliDebug(999, Form("STU DCAL THR %d EGA %d EJE %d", ithr, fSTUDCAL->GetThreshold(kL1GammaHigh + ithr), fSTUDCAL->GetThreshold(kL1JetHigh + ithr)));
      }
  
    }
  }


  if (fSTUDCAL && fMedianMode) { // Execute run 2 algorithm
    AliDebug(999,"Executing event-by-event median background subtraction for trigger");	

    Int_t fBkgRhoEMCAL = fSTU->GetMedianEnergy();
    Int_t fBkgRhoDCAL = fSTUDCAL->GetMedianEnergy();


    AliDebug(999,Form("Found Rho from EMCAL: %d,    Rho from DCAL: %d",fBkgRhoEMCAL,fBkgRhoDCAL));

    fSTU->SetBkgRho(fBkgRhoDCAL);
    fSTUDCAL->SetBkgRho(fBkgRhoEMCAL);

  } else {
    AliDebug(999,"Not doing event-by-event background subtraction for trigger");
  } 
  

  for (int ithr = 0; ithr < 2; ithr++) {
    //
    fSTU->Reset();
  
    fSTU->L1(kL1GammaHigh + ithr);
    
    TIterator* nP = 0x0;
    
    nP = (fSTU->Patches()).MakeIterator();
    
    AliDebug(999, Form("=== EMCAL STU found %d gamma patches", (fSTU->Patches()).GetEntriesFast()));
    
    while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)nP->Next()) 
    {			
      p->Position(px, py);
      
      if (AliDebugLevel()) p->Print("");
      
      if (fGeometry->GetAbsFastORIndexFromPositionInEMCAL(py, px, id))
      {
        if (posMap[px][py] == -1)
        {
          posMap[px][py] = digits->GetEntriesFast();
          
          // Add a new digit
          new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
          
          dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);
        } 
        else
        {

          dig = (AliEMCALTriggerRawDigit*)digits->At(posMap[px][py]);								
        }
        
        if (AliDebugLevel()) dig->Print("");
        
        if (dig) dig->SetTriggerBit(kL1GammaHigh + ithr, 0);
      }
    }

    fSTU->Reset();
    
    fSTU->L1(kL1JetHigh + ithr);
    
    nP = (fSTU->Patches()).MakeIterator();
    
    AliDebug(999, Form("=== EMCAL STU found %d jet patches", (fSTU->Patches()).GetEntriesFast()));
    
    while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)nP->Next()) 
    {			
      p->Position(px, py);
      
      if (AliDebugLevel()) p->Print("");
      
      if (fGeometry->GetAbsFastORIndexFromPositionInEMCAL(py, px, id))
      {
        if (posMap[px][py] == -1)
        {
          posMap[px][py] = digits->GetEntriesFast();
          
          // Add a new digit
          new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
          
          dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);
        } 
        else
        {
          dig = (AliEMCALTriggerRawDigit*)digits->At(posMap[px][py]);
        }
        
        if (AliDebugLevel()) dig->Print("");
        
        if (dig) dig->SetTriggerBit(kL1JetHigh + ithr, 0);
      }
    }

    // DCAL STU
    if (fSTUDCAL) {
      Int_t iDCALOffset = fSTU->RegionSize()->X();

      // Gamma L1

      fSTUDCAL->Reset();
      
      fSTUDCAL->L1(kL1GammaHigh + ithr);
      
      nP = (fSTUDCAL->Patches()).MakeIterator();
      
      AliDebug(999, Form("=== DCAL STU found %d gamma patches", (fSTUDCAL->Patches()).GetEntriesFast()));

      while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)nP->Next()) 
      {			
        p->Position(px, py);
        px += iDCALOffset;
        if (AliDebugLevel()) p->Print("");
        
        if (fGeometry->GetAbsFastORIndexFromPositionInEMCAL(py, px, id))
        {
          if (posMap[px][py] == -1)
          {
            posMap[px][py] = digits->GetEntriesFast();
            
            // Add a new digit
            new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
            
            dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);
          } 
          else
          {
            dig = (AliEMCALTriggerRawDigit*)digits->At(posMap[px][py]);								
          }
          
          if (AliDebugLevel()) dig->Print("");
          
          if (dig) dig->SetTriggerBit(kL1GammaHigh + ithr, 0);
        }
      }
      // Jet L1 

      fSTUDCAL->Reset();
      
      fSTUDCAL->L1(kL1JetHigh + ithr);
      
      nP = (fSTUDCAL->Patches()).MakeIterator();
      
      AliDebug(999, Form("=== DCAL STU found %d jet patches", (fSTUDCAL->Patches()).GetEntriesFast()));

      while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)nP->Next()) 
      {			
        p->Position(px, py);
        px += iDCALOffset;
        if (AliDebugLevel()) p->Print("");
        
        if (fGeometry->GetAbsFastORIndexFromPositionInEMCAL(py, px, id))
        {
          if (posMap[px][py] == -1)
          {
            posMap[px][py] = digits->GetEntriesFast();
            
            // Add a new digit
            new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
            
            dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);
          } 
          else
          {
            dig = (AliEMCALTriggerRawDigit*)digits->At(posMap[px][py]);
          }
          
          if (AliDebugLevel()) dig->Print("");
          

          if (dig) dig->SetTriggerBit(kL1JetHigh + ithr, 0);
        }
      }

    }

  }
  // Now reset the electronics for a fresh start with next event
  Reset();
}

///
/// Reset
//__________________
void AliEMCALTriggerElectronics::Reset()
{	
  TIter nextTRU(fTRU);
  while ( AliEMCALTriggerTRU *TRU = (AliEMCALTriggerTRU*)nextTRU() ) TRU->Reset();
  
  fSTU->Reset();
  if (fSTUDCAL) fSTUDCAL->Reset();
}
