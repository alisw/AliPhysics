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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliGRPManager class                                                    //
// The class can be used in order to access and read the Global Run       //
// Parameters entry from OCDB.                                            //
// It has a methods to set the magnetic field instanton and return        //
// the run and event info objects.                                        //
//                                                                        //
// cvetan.cheshkov@cern.ch 15/06/2009                                     //
//                                                                        //
// Usage:                                                                 //
// AliGRPManager grpMan;                                                  //
// Bool_t status = kTRUE;                                                 //
// status = grpMan.ReadGRPEntry(); // Read the corresponding OCDB entry   //
// status = grpMan.SetMagField();  // Set global field instanton          //
// AliRunInfo *runInfo = grpMan.GetRunInfo();// Get instance of run info  //
//                                                                        //
// Note: CDB manager should be initialized beforehand                     //
////////////////////////////////////////////////////////////////////////////

#include <TGeoGlobalMagField.h>
#include <TPRegexp.h>

#include "AliGRPManager.h"
#include "AliLog.h"
#include "AliRunInfo.h"
#include "AliGRPObject.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliMagF.h"

ClassImp(AliGRPManager)

//_____________________________________________________________________________
AliGRPManager::AliGRPManager() :
  TObject(),
  fGRPData(NULL)
{
  // Default constructor
}

//_____________________________________________________________________________
AliGRPManager::~AliGRPManager()
{
  // Destructor
  if (fGRPData) delete fGRPData;
}

//_____________________________________________________________________________
Bool_t AliGRPManager::ReadGRPEntry()
{
  //------------------------------------
  // Initialization of the GRP entry. 
  // Returns kTRUE in case of success.
  //------------------------------------

  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");

  if (entry) {

    TMap* m = dynamic_cast<TMap*>(entry->GetObject());  // old GRP entry

    if (m) {
       AliInfo("Found a TMap in GRP/GRP/Data, converting it into an AliGRPObject");
       m->Print();
       fGRPData = new AliGRPObject();
       fGRPData->ReadValuesFromMap(m);
    }

    else {
       AliInfo("Found an AliGRPObject in GRP/GRP/Data, reading it");
       fGRPData = dynamic_cast<AliGRPObject*>(entry->GetObject());  // new GRP entry
       entry->SetOwner(0);
    }

    AliCDBManager::Instance()->UnloadFromCache("GRP/GRP/Data");
  }

  if (!fGRPData) {
     AliError("No GRP entry found in OCDB!");
     return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliGRPManager::SetMagField()
{
  // Dealing with the magnetic field map
  // Construct the mag field map from the data in GRP
  // Set the global mag field instance

  if ( TGeoGlobalMagField::Instance()->IsLocked() ) {
    AliInfo("Running with the externally locked B field !");
  }
  else {
    // Construct the field map out of the information retrieved from GRP.
    Bool_t ok = kTRUE;
    // L3
    Float_t l3Current = fGRPData->GetL3Current((AliGRPObject::Stats)0);
    if (l3Current == AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 current !");
      ok = kFALSE;
    }
    
    Char_t l3Polarity = fGRPData->GetL3Polarity();
    if (l3Polarity == AliGRPObject::GetInvalidChar()) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 polarity !");
      ok = kFALSE;
    }

    // Dipole
    Float_t diCurrent = fGRPData->GetDipoleCurrent((AliGRPObject::Stats)0);
    if (diCurrent == AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole current !");
      ok = kFALSE;
    }

    Char_t diPolarity = fGRPData->GetDipolePolarity();
    if (diPolarity == AliGRPObject::GetInvalidChar()) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole polarity !");
      ok = kFALSE;
    }

    TString beamType = fGRPData->GetBeamType();
    if (beamType==AliGRPObject::GetInvalidString()) {
      AliError("GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
      beamType = "UNKNOWN";
      ok = kFALSE;
    }

    Float_t beamEnergy = fGRPData->GetBeamEnergy();
    if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
      beamEnergy = 0;
      ok = kFALSE;
    }
    beamEnergy /= 120E3;       // energy is provided in MeV*120

    // read special bits for the polarity convention and map type
    Int_t  polConvention = fGRPData->IsPolarityConventionLHC() ? AliMagF::kConvLHC : AliMagF::kConvDCS2008;
    Bool_t uniformB = fGRPData->IsUniformBMap();

    if (ok) { 
      if ( !SetFieldMap(l3Current, diCurrent, l3Polarity ? -1:1, diPolarity ? -1:1, 
			polConvention,uniformB,beamEnergy, beamType.Data())) {
	AliError("Failed to create a B field map !");
	ok = kFALSE;
      }
      AliInfo("Running with the B field constructed out of GRP !");
    }
    else {
      AliError("B field is neither set nor constructed from GRP ! Exitig...");
    }
    return ok;
  }

  return kTRUE;
}

//_____________________________________________________________________________
AliRunInfo* AliGRPManager::GetRunInfo()
{
  // Constructs and returns an object
  // containing the run information
  // The user code is the owner of the object

  TString lhcState = fGRPData->GetLHCState();
  if (lhcState==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the LHC state ! Using UNKNOWN");
    lhcState = "UNKNOWN";
  }

  TString beamType = fGRPData->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
    beamType = "UNKNOWN";
  }

  Float_t beamEnergy = fGRPData->GetBeamEnergy();
  if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
    beamEnergy = 0;
  }
  // energy is provided in MeV*120
  beamEnergy /= 120E3;

  TString runType = fGRPData->GetRunType();
  if (runType==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the run type ! Using UNKNOWN");
    runType = "UNKNOWN";
  }

  Int_t activeDetectors = fGRPData->GetDetectorMask();
  if (activeDetectors==AliGRPObject::GetInvalidUInt()) {
    AliError("GRP/GRP/Data entry:  missing value for the detector mask ! Using 1074790399");
    activeDetectors = 1074790399;
  }

  return new AliRunInfo(lhcState, beamType, beamEnergy, runType, activeDetectors);
}

//_____________________________________________________________________________
Bool_t AliGRPManager::SetFieldMap(Float_t l3Cur, Float_t diCur, Float_t l3Pol, 
				  Float_t diPol, Int_t convention, Bool_t uniform,
				  Float_t beamenergy, const Char_t *beamtype, const Char_t *path) 
{
  //------------------------------------------------
  // The magnetic field map, defined externally...
  // L3 current 30000 A  -> 0.5 T
  // L3 current 12000 A  -> 0.2 T
  // dipole current 6000 A
  // The polarities must match the convention (LHC or DCS2008) 
  // unless the special uniform map was used for MC
  //------------------------------------------------
  const Float_t l3NominalCurrent1=30000.; // (A)
  const Float_t l3NominalCurrent2=12000.; // (A)
  const Float_t diNominalCurrent =6000. ; // (A)

  const Float_t tolerance=0.03; // relative current tolerance
  const Float_t zero=77.;       // "zero" current (A)
  //
  AliMagF::BMap_t map;
  double sclL3,sclDip;
  //
  l3Cur = TMath::Abs(l3Cur);
  diCur = TMath::Abs(diCur);
  //
  if (TMath::Abs((sclDip=diCur/diNominalCurrent)-1.) > tolerance && !uniform) {
    if (diCur <= zero) sclDip = 0.; // some small current.. -> Dipole OFF
    else {
      AliError(Form("Wrong dipole current (%f A)!",diCur));
      return kFALSE;
    }
  }
  //
  if (uniform) { 
    // special treatment of special MC with uniform mag field (normalized to 0.5 T)
    // no check for scaling/polarities are done
    map   = AliMagF::k5kGUniform;
    sclL3 = l3Cur/l3NominalCurrent1; 
  }
  else {
    if      (TMath::Abs((sclL3=l3Cur/l3NominalCurrent1)-1.) < tolerance) map  = AliMagF::k5kG;
    else if (TMath::Abs((sclL3=l3Cur/l3NominalCurrent2)-1.) < tolerance) map  = AliMagF::k2kG;
    else if (l3Cur <= zero)                                { sclL3 = 0;  map  = AliMagF::k5kGUniform;}
    else {
      AliError(Form("Wrong L3 current (%f A)!",l3Cur));
      return kFALSE;
    }
  }
  //
  if (sclDip!=0 && (map==AliMagF::k5kG || map==AliMagF::k2kG) &&
      ((convention==AliMagF::kConvLHC     && l3Pol!=diPol) ||
       (convention==AliMagF::kConvDCS2008 && l3Pol==diPol)) ) { 
    AliError(Form("Wrong combination for L3/Dipole polarities (%c/%c) for convention %d",
		  l3Pol>0?'+':'-',diPol>0?'+':'-',AliMagF::GetPolarityConvention()));
    return kFALSE;
  }
  //
  if (l3Pol<0) sclL3  = -sclL3;
  if (diPol<0) sclDip = -sclDip;
  //
  AliMagF::BeamType_t btype = AliMagF::kNoBeamField;
  TString btypestr = beamtype;
  btypestr.ToLower();
  TPRegexp protonBeam("(proton|p)\\s*-?\\s*\\1");
  TPRegexp ionBeam("(lead|pb|ion|a)\\s*-?\\s*\\1");
  if (btypestr.Contains(ionBeam)) btype = AliMagF::kBeamTypeAA;
  else if (btypestr.Contains(protonBeam)) btype = AliMagF::kBeamTypepp;
  else {
    AliInfo(Form("Cannot determine the beam type from %s, assume no LHC magnet field",beamtype));
  }
  char ttl[50];
  sprintf(ttl,"L3: %+5d Dip: %+4d kA; %s",(int)TMath::Sign(l3Cur,float(sclL3)),
	  (int)TMath::Sign(diCur,float(sclDip)),uniform ? " Constant":"");
  AliMagF* fld = new AliMagF("MagneticFieldMap", ttl, 2, sclL3, sclDip, 10., map, path, 
			     btype,beamenergy);
  TGeoGlobalMagField::Instance()->SetField( fld );
  TGeoGlobalMagField::Instance()->Lock();
  //
  return kTRUE;
}
