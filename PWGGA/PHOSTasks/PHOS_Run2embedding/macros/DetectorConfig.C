/*
 * AliDPG - ALICE Experiment Data Preparation Group
 * Detector configuration script
 *
 */

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

enum EDetector_t {
  kDetectorDefault,
  kDetectorMuon,
  kDetectorCustom,
  kNDetectors
};

const Char_t *DetectorName[kNDetectors] = {
  "Default",
  "Muon",
  "Custom"
};

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

enum EMagnet_t {
  kMagnetDefault,
  kMagnetOff,
  kMagnet5kG,
  kNMagnets
};

const Char_t *MagnetName[kNMagnets] = {
  "Default",
  "Off",
  "5kG"
};

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

Int_t iABSO   = 1;
Int_t iACORDE = 1;
Int_t iAD     = 1;
Int_t iDIPO   = 1;
Int_t iEMCAL  = 1;
Int_t iFMD    = 1;
Int_t iFRAME  = 1;
Int_t iHALL   = 1;
Int_t iITS    = 1;
Int_t iMAG    = 1;
Int_t iMUON   = 1;
Int_t iPHOS   = 1;
Int_t iPIPE   = 1;
Int_t iPMD    = 1;
Int_t iHMPID  = 1;
Int_t iSHIL   = 1;
Int_t iT0     = 1;
Int_t iTOF    = 1;
Int_t iTPC    = 1;
Int_t iTRD    = 1;
Int_t iVZERO  = 1;
Int_t iZDC    = 1;
  
/*****************************************************************/

void
DetectorConfig(Int_t tag, Int_t run)
{

  switch (tag) {

    // kDetectorDefault
  case kDetectorDefault:
    DetectorDefault(run);
    break;
    
    // kDetectorMuon
  case kDetectorMuon:
    DetectorMuon(run);
    break;
    
    // kDetectorCustom
  case kDetectorCustom:
    if ((gROOT->LoadMacro("DetectorCustom.C")) != 0) {
      printf("ERROR: cannot find DetectorCustom.C\n");
      abort();
      return;
    }
    DetectorCustom(run);
    break;
    
  }
  
}

/*****************************************************************/

void
DetectorDefault(Int_t run)
{
  /*
   * DetectorDefault
   * configures the detectors to the default 
   * configuration automatically according to run number
   *
   */
  
  gROOT->LoadMacro("Sim/Utils.C");
  Int_t year = RunToYear(run);
  
  iABSO   = 0;
  iACORDE = 0;
  iAD     = year < 2015 ? 0 : 0;
  iDIPO   = 0;
  iEMCAL  = 0;
  iFMD    = 0;
  iFRAME  = 1;
  iHALL   = 0;
  iITS    = 1;
  iMAG    = 1;
  iMUON   = 0;
  iPHOS   = 1;
  iPIPE   = 1;
  iPMD    = 0;
  iHMPID  = 0;
  iSHIL   = 0;
  iT0     = 0;
  iTOF    = 1;
  iTPC    = 1;
  iTRD    = 1;
  iVZERO  = 0;
  iZDC    = 0;

  DetectorInit(run);
}
  
/*****************************************************************/

void
DetectorMuon(Int_t run)
{
  /*
   * DetectorMuon
   * configures the detectors to the Muon 
   * configuration automatically according to run number
   *
   */
  
  gROOT->LoadMacro("Sim/Utils.C");
  Int_t year = RunToYear(run);
  
  iABSO   = 1;
  iACORDE = 0;
  iAD     = year < 2015 ? 0 : 1;
  iDIPO   = 1;
  iEMCAL  = 0;
  iFMD    = 1;
  iFRAME  = 1;
  iHALL   = 1;
  iITS    = 1;
  iMAG    = 1;
  iMUON   = 1;
  iPHOS   = 0;
  iPIPE   = 1;
  iPMD    = 0;
  iHMPID  = 0;
  iSHIL   = 1;
  iT0     = 1;
  iTOF    = 0;
  iTPC    = 0;
  iTRD    = 0;
  iVZERO  = 1;
  iZDC    = 0;

  DetectorInit(run);
}
  
/*****************************************************************/

void
DetectorInit(Int_t run)
{
  /*
   * DetectorInit
   * initialise the detectors to the default 
   * configuration automatically according to run number
   *
   */
  
  gROOT->LoadMacro("Sim/Utils.C");
  Int_t year = RunToYear(run);


  //=================== Alice BODY parameters =============================
  AliBODY *BODY = new AliBODY("BODY", "Alice envelop");
  
  
  if (iMAG)
    {
      //=================== MAG parameters ============================
      // --- Start with Magnet since detector layouts may be depending ---
      // --- on the selected Magnet dimensions ---
      AliMAG *MAG = new AliMAG("MAG", "Magnet");
    }
  

  if (iABSO)
    {
      //=================== ABSO parameters ============================
      AliABSO *ABSO = new AliABSOv3("ABSO", "Muon Absorber");
    }

  if (iDIPO)
    {
      //=================== DIPO parameters ============================

      AliDIPO *DIPO = new AliDIPOv3("DIPO", "Dipole version 3");
    }

  if (iHALL)
    {
      //=================== HALL parameters ============================

      AliHALL *HALL = new AliHALLv3("HALL", "Alice Hall");
    }


  if (iFRAME)
    {
      //=================== FRAME parameters ============================

      if (year < 2015) {
 	AliFRAMEv2 *FRAME = new AliFRAMEv2("FRAME", "Space Frame");
	FRAME->SetHoles(1);
      }
      else {
 	AliFRAMEv3 *FRAME = new AliFRAMEv3("FRAME", "Space Frame");
	FRAME->SetHoles(1);
      }
    }

  if (iSHIL)
    {
      //=================== SHIL parameters ============================

      AliSHIL *SHIL = new AliSHILv3("SHIL", "Shielding Version 3");
    }


  if (iPIPE)
    {
      //=================== PIPE parameters ============================

      AliPIPE *PIPE = new AliPIPEv3("PIPE", "Beam Pipe");
    }
 
  if (iITS)
    {
      //=================== ITS parameters ============================
      
      AliITS *ITS  = new AliITSv11("ITS","ITS v11");
    }

  if (iTPC)
    {
      //============================ TPC parameters =====================

      AliTPC *TPC = new AliTPCv2("TPC", "Default");
    }


  if (iTOF) {
    //=================== TOF parameters ============================

    AliTOF *TOF = new AliTOFv6T0("TOF", "normal TOF");
  }


  if (iHMPID)
    {
      //=================== HMPID parameters ===========================

      AliHMPID *HMPID = new AliHMPIDv3("HMPID", "normal HMPID");
    }


  if (iZDC)
    {
      //=================== ZDC parameters ============================

      if (year == 2010) {
	AliZDC *ZDC = new AliZDCv3("ZDC", "normal ZDC");
	ZDC->SetSpectatorsTrack();
        ZDC->SetLumiLength(0.);
      }
      else if (year < 2015) {
	AliZDC *ZDC = new AliZDCv3("ZDC", "normal ZDC");
	//Collimators aperture
	ZDC->SetVCollSideCAperture(0.85);
	ZDC->SetVCollSideCCentre(0.);
	ZDC->SetVCollSideAAperture(0.75);
	ZDC->SetVCollSideACentre(0.);
	//Detector position
	ZDC->SetYZNC(1.6);
	ZDC->SetYZNA(1.6);
	ZDC->SetYZPC(1.6);
	ZDC->SetYZPA(1.6);
      }
      else {
	AliZDC *ZDC = new AliZDCv4("ZDC", "normal ZDC");
	ZDC->SetLumiLength(0.);
	ZDC->SetVCollSideCAperture(2.8);
	ZDC->SetVCollSideCApertureNeg(2.8);
      }
    }

  if (iTRD)
    {
      //=================== TRD parameters ============================

      AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");
      AliTRDgeometry *geoTRD = TRD->GetGeometry();
      // Partial geometry: modules at 0,1,7,8,9,16,17
      // starting at 3h in positive direction

      if (year == 2010) {
	geoTRD->SetSMstatus(2,0);
	geoTRD->SetSMstatus(3,0);
	geoTRD->SetSMstatus(4,0);
	geoTRD->SetSMstatus(5,0);
	geoTRD->SetSMstatus(6,0);
	geoTRD->SetSMstatus(11,0);
	geoTRD->SetSMstatus(12,0);
	geoTRD->SetSMstatus(13,0);
	geoTRD->SetSMstatus(14,0);
	geoTRD->SetSMstatus(15,0);
	geoTRD->SetSMstatus(16,0);	
      }
      else if (year == 2011) {
	geoTRD->SetSMstatus(2,0);
	geoTRD->SetSMstatus(3,0);
	geoTRD->SetSMstatus(4,0);
	geoTRD->SetSMstatus(5,0);
	geoTRD->SetSMstatus(6,0);
	geoTRD->SetSMstatus(12,0);
	geoTRD->SetSMstatus(13,0);
	geoTRD->SetSMstatus(14,0);
      }
    }
  
  if (iFMD)
    {
      //=================== FMD parameters ============================
      
      AliFMD *FMD = new AliFMDv1("FMD", "normal FMD");
    }

  if (iMUON)
    {
      //=================== MUON parameters ===========================
      // New MUONv1 version (geometry defined via builders)
      AliMUON *MUON = new AliMUONv1("MUON", "default");
      // activate trigger efficiency by cells
      MUON->SetTriggerEffCells(1); // backward compatibility
      MUON->SetTriggerResponseV1(2); // backward compatibility
    }

  if (iPHOS)
    {
      //=================== PHOS parameters ===========================

      if (year < 2015) {
	AliPHOS *PHOS = new AliPHOSv1("PHOS", "noCPV_Modules123");
      }
      else {
	AliPHOS *PHOS = new AliPHOSv1("PHOS", "Run2");
      }
	
    }


  if (iPMD)
    {
      //=================== PMD parameters ============================

      AliPMD *PMD = new AliPMDv1("PMD", "normal PMD");
    }

  if (iT0)
    {
      //=================== T0 parameters ============================
      AliT0 *T0 = new AliT0v1("T0", "T0 Detector");
    }

  if (iEMCAL)
    {
      //=================== EMCAL parameters ============================

      if (year < 2015) {
	AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_FIRSTYEARV1");
      }
      else {
	AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETE12SMV1_DCAL_8SM", kFALSE);
      }
	
    }

  if (iACORDE)
    {
      //=================== ACORDE parameters ============================

      AliACORDE *ACORDE = new AliACORDEv1("ACORDE", "normal ACORDE");
    }

  if (iVZERO)
    {
      //=================== ACORDE parameters ============================
      
      AliVZERO *VZERO = new AliVZEROv7("VZERO", "normal VZERO");
    }  

  if (iAD){
    //=================== AD parameters ============================

    AliAD *AD = new AliADv1("AD", "normal AD");
  }         
  
}

