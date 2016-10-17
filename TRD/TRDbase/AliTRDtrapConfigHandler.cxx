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
//  MCM configuraton handler                                              //
//                                                                        //
//  Author: U. Westerhoff                                                 //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTRDtrapConfigHandler.h"

#include <iostream>
#include <sstream>
#include <iomanip>

#include "AliLog.h"

#include "AliTRDfeeParam.h"
#include "AliTRDtrapConfig.h"
#include "AliTRDmcmSim.h"
#include "AliTRDgeometry.h"
#include "AliTRDcalibDB.h"

#include "TMath.h"
#include "TGeoMatrix.h"
#include "TGraph.h"

#include "AliTRDCalOnlineGainTable.h"
#include "AliTRDCalOnlineGainTableROC.h"
#include "AliTRDCalOnlineGainTableMCM.h"

using namespace std;

AliTRDtrapConfigHandler::AliTRDtrapConfigHandler(AliTRDtrapConfig *cfg) :
     ltuParam()
     , fRestrictiveMask((0x3ffff << 11) | (0x1f << 6) | 0x3f)
     , fTrapConfig(cfg)
     , fGtbl()
{

}


AliTRDtrapConfigHandler::~AliTRDtrapConfigHandler()
{

}

void AliTRDtrapConfigHandler::Init()
{
  if (!fTrapConfig) {
    AliError("No TRAPconfig given");
    return;
  }

  // setup of register allocation
  // I/O configuration which we don't care about
  fTrapConfig->SetTrapRegAlloc(AliTRDtrapConfig::kSEBDOU, AliTRDtrapConfig::kAllocNone);
  // position look-up table by layer
  for (Int_t iBin = 0; iBin < 128; iBin++)
    fTrapConfig->SetTrapRegAlloc((AliTRDtrapConfig::TrapReg_t) (AliTRDtrapConfig::kTPL00 + iBin), AliTRDtrapConfig::kAllocByLayer);
  // ... individual
  fTrapConfig->SetTrapRegAlloc(AliTRDtrapConfig::kC14CPUA, AliTRDtrapConfig::kAllocByMCM);
  fTrapConfig->SetTrapRegAlloc(AliTRDtrapConfig::kC15CPUA, AliTRDtrapConfig::kAllocByMCM);

  // setup of DMEM allocation
  for(Int_t iAddr = AliTRDtrapConfig::fgkDmemStartAddress;
            iAddr < (AliTRDtrapConfig::fgkDmemWords + AliTRDtrapConfig::fgkDmemStartAddress); iAddr++) {

    if(iAddr == AliTRDmcmSim::fgkDmemAddrDeflCorr)
      fTrapConfig->SetDmemAlloc(iAddr, AliTRDtrapConfig::kAllocByMCMinSM);

    else if(iAddr == AliTRDmcmSim::fgkDmemAddrNdrift)
      fTrapConfig->SetDmemAlloc(iAddr, AliTRDtrapConfig::kAllocByDetector);

    else if((iAddr >= AliTRDmcmSim::fgkDmemAddrDeflCutStart) && (iAddr <= AliTRDmcmSim::fgkDmemAddrDeflCutEnd))
      fTrapConfig->SetDmemAlloc(iAddr, AliTRDtrapConfig::kAllocByMCMinSM);

    else if((iAddr >= AliTRDmcmSim::fgkDmemAddrTrackletStart) && (iAddr <= AliTRDmcmSim::fgkDmemAddrTrackletEnd))
      fTrapConfig->SetDmemAlloc(iAddr, AliTRDtrapConfig::kAllocByMCM);

    else if((iAddr >= AliTRDmcmSim::fgkDmemAddrLUTStart) && (iAddr <= AliTRDmcmSim::fgkDmemAddrLUTEnd))
      fTrapConfig->SetDmemAlloc(iAddr, AliTRDtrapConfig::kAllocGlobal);

    else if(iAddr == AliTRDmcmSim::fgkDmemAddrLUTcor0)
      fTrapConfig->SetDmemAlloc(iAddr, AliTRDtrapConfig::kAllocByMCMinSM);

    else if(iAddr == AliTRDmcmSim::fgkDmemAddrLUTcor1)
      fTrapConfig->SetDmemAlloc(iAddr, AliTRDtrapConfig::kAllocByMCMinSM);

    else if(iAddr == AliTRDmcmSim::fgkDmemAddrLUTnbins)
      fTrapConfig->SetDmemAlloc(iAddr, AliTRDtrapConfig::kAllocGlobal);

    else if(iAddr == AliTRDmcmSim::fgkDmemAddrLUTLength)
      fTrapConfig->SetDmemAlloc(iAddr, AliTRDtrapConfig::kAllocGlobal);

    else
      fTrapConfig->SetDmemAlloc(iAddr, AliTRDtrapConfig::kAllocGlobal);
  }
}

void AliTRDtrapConfigHandler::ResetMCMs()
{
   //
   // Reset all MCM registers and DMEM
   //

  if (!fTrapConfig) {
    AliError("No TRAPconfig given");
    return;
  }

   fTrapConfig->ResetRegs();
   fTrapConfig->ResetDmem();
}


Int_t AliTRDtrapConfigHandler::LoadConfig()
{
  // load a default configuration which is suitable for simulation
  // for a detailed description of the registers see the TRAP manual
  // if you want to resimulate tracklets on real data use the appropriate config instead

  if (!fTrapConfig) {
    AliError("No TRAPconfig given");
    return -1;
  }

  // prepare ltuParam
  // ndrift (+ 5 binary digits)
  ltuParam.SetNtimebins(20 << 5);
  // deflection + tilt correction
  ltuParam.SetRawOmegaTau(0.16133);
  // deflection range table
  ltuParam.SetRawPtMin(0.1);
  // magnetic field
  ltuParam.SetRawMagField(0.0);
  // scaling factors for q0, q1
  ltuParam.SetRawScaleQ0(0);
  ltuParam.SetRawScaleQ1(0);
  // disable length correction and tilting correction
  ltuParam.SetRawLengthCorrectionEnable(kFALSE);
  ltuParam.SetRawTiltCorrectionEnable(kFALSE);

  for (Int_t iDet = 0; iDet < 540; iDet++) {
    // HC header configuration bits
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kC15CPUA, 0x2102, iDet); // zs, deh

    // no. of timebins
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kC13CPUA, 24, iDet);

    // pedestal filter
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kFPNP, 4*10, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kFPTC, 0, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kFPBY, 0, iDet); // bypassed!

    // gain filter
    for (Int_t adc = 0; adc < 20; adc++) {
      fTrapConfig->SetTrapReg(AliTRDtrapConfig::TrapReg_t(AliTRDtrapConfig::kFGA0+adc), 40, iDet);
      fTrapConfig->SetTrapReg(AliTRDtrapConfig::TrapReg_t(AliTRDtrapConfig::kFGF0+adc), 15, iDet);
    }
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kFGTA, 20, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kFGTB, 2060, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kFGBY, 0, iDet);  // bypassed!

    // tail cancellation
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kFTAL, 200, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kFTLL, 0, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kFTLS, 200, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kFTBY, 0, iDet);

    // tracklet calculation
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kTPQS0, 5, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kTPQE0, 10, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kTPQS1, 11, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kTPQE1, 20, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kTPFS, 5, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kTPFE, 20, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kTPVBY, 0, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kTPVT, 10, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kTPHT, 150, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kTPFP, 40, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kTPCL, 1, iDet);
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kTPCT, 10, iDet);

    // apply ltuParams
    ConfigureDyCorr(iDet);
    ConfigureDRange(iDet); // deflection range
    ConfigureNTimebins(iDet);  // timebins in the drift region
    ConfigurePIDcorr(iDet);  // scaling parameters for the PID

    // event buffer
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kEBSF, 1, iDet);  // 0: store filtered; 1: store unfiltered

    // zs applied to data stored in event buffer (sel. by EBSF)
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kEBIS, 15, iDet); // single indicator threshold (plus two digits)
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kEBIT, 30, iDet); // sum indicator threshold (plus two digits)
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kEBIL, 0xf0, iDet);   // lookup table
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kEBIN, 0, iDet);      // neighbour sensitivity

    // raw data
    fTrapConfig->SetTrapReg(AliTRDtrapConfig::kNES, (0x0000 << 16) | 0x1000, iDet);
  }

  // ****** hit position LUT

  // now calculate it from PRF
  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();

  Double_t padResponse[3]; // pad response left, central, right
  Double_t padResponseR[3]; // pad response left, central, right
  Double_t padResponseL[3]; // pad response left, central, right

  for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
    TGraph gr(128);
    for (Int_t iBin = 0; iBin < 256*0.5; iBin++) {
      cal->PadResponse(1., iBin*1./256.,    iLayer, padResponse);
      cal->PadResponse(1., iBin*1./256.-1., iLayer, padResponseR);
      cal->PadResponse(1., iBin*1./256.+1., iLayer, padResponseL);
      gr.SetPoint(iBin, (0.5 * (padResponseR[1] - padResponseL[1])/padResponse[1] * 256), iBin);
    }
    for (Int_t iBin = 0; iBin < 128; iBin++) {
      Int_t corr = (Int_t) (gr.Eval(iBin)) - iBin;
      if (corr < 0)
        corr = 0;
      else if (corr > 31)
        corr = 31;
      for (Int_t iStack = 0; iStack < 540/6; iStack++) {
        fTrapConfig->SetTrapReg((AliTRDtrapConfig::TrapReg_t) (AliTRDtrapConfig::kTPL00 + iBin), corr, 6*iStack + iLayer);
      }
    }
  }
  // ****** hit position LUT configuration end

  return 0;
}


Int_t AliTRDtrapConfigHandler::LoadConfig(TString filename)
{
   //
  // load a TRAP configuration from a file
   // The file format is the format created by the standalone
   // command coder scc / show_cfdat but without the running number
   // scc /show_cfdat adds as the first column
   // which are two tools to inspect/export configurations from wingDB
   //

  if (!fTrapConfig) {
    AliError("No TRAPconfig given");
    return -1;
  }

   Int_t ignoredLines=0;
   Int_t ignoredCmds=0;
   Int_t readLines=0;


   AliDebug(5, Form("Processing file %s", filename.Data()));
   std::ifstream infile;
   infile.open(filename.Data(), std::ifstream::in);
   if (!infile.is_open()) {
    AliError(Form("Can not open MCM configuration file %s", filename.Data()));
    return kFALSE;
   }

   UInt_t cmd;
   Int_t extali, addr, data;

   // reset restrictive mask
   fRestrictiveMask = (0x3ffff << 11) | (0x1f << 6) | 0x3f;
   char linebuffer[512];
   istringstream line;

   while(infile.getline(linebuffer, 512) && infile.good()) {
      line.clear();
      line.str(linebuffer);
      cmd=999;
      extali=-1;
      addr=-1;
      data=-1;
      line >> std::skipws >> cmd >> addr >> data >> extali;  // the lines read from config file can contain additional columns.
      // Therefore the detour via istringstream

      if(cmd!=999 && addr != -1 && extali!=-1) {

	 if(cmd==fgkScsnCmdWrite) {
	    for(Int_t det=0; det<AliTRDgeometry::Ndet(); det++) {
	       UInt_t rocpos = (1 << (AliTRDgeometry::GetSector(det)+11)) | (1 << (AliTRDgeometry::GetStack(det)+6)) | (1 << AliTRDgeometry::GetLayer(det));
	       AliDebug(1, Form("checking restriction: mask=0x%08x, rocpos=0x%08x", fRestrictiveMask, rocpos));
	       if ((fRestrictiveMask & rocpos) == rocpos) {
		 AliDebug(1, Form("match: %i %i %i %i", cmd, extali, addr, data));
		  AddValues(det, cmd, extali, addr, data);
	       }
	    }
	 }

	 else if(cmd == fgkScsnLTUparam) {
	    ProcessLTUparam(extali, addr, data);
	 }

	 else if(cmd == fgkScsnCmdRestr) {
	    fRestrictiveMask = data;
	   AliDebug(1, Form("updated restrictive mask to 0x%08x", fRestrictiveMask));
	 }

	 else if((cmd == fgkScsnCmdReset) ||
		 (cmd == fgkScsnCmdRobReset)) {
	   fTrapConfig->ResetRegs();
	 }

	 else if (cmd == fgkScsnCmdSetHC) {
	   Int_t fullVersion = ((data & 0x7F00) >> 1) | (data & 0x7f);

	   for (Int_t iDet = 0; iDet < AliTRDgeometry::Ndet(); iDet++) {
	     Int_t smls = (AliTRDgeometry::GetSector(iDet) << 6) | (AliTRDgeometry::GetLayer(iDet) << 3) | AliTRDgeometry::GetStack(iDet);

	     for (Int_t iRob = 0; iRob < 8; iRob++) {
	       // HC mergers
	       fTrapConfig->SetTrapReg(AliTRDtrapConfig::kC14CPUA, 0xc << 16, iDet, iRob, 17);
	       fTrapConfig->SetTrapReg(AliTRDtrapConfig::kC15CPUA, ((1<<29) | (fullVersion<<15) | (1<<12) | (smls<<1) | (iRob%2)), iDet, iRob, 17);

	       // board mergers
	       fTrapConfig->SetTrapReg(AliTRDtrapConfig::kC14CPUA, 0, iDet, iRob, 16);
	       fTrapConfig->SetTrapReg(AliTRDtrapConfig::kC15CPUA, ((1<<29) | (fullVersion<<15) | (1<<12) | (smls<<1) | (iRob%2)), iDet, iRob, 16);

	       // and now for the others
	       for (Int_t iMcm = 0; iMcm < 16; iMcm++) {
		 fTrapConfig->SetTrapReg(AliTRDtrapConfig::kC14CPUA, iMcm | (iRob << 4) | (3 << 16), iDet, iRob, iMcm);
		 fTrapConfig->SetTrapReg(AliTRDtrapConfig::kC15CPUA, ((1<<29) | (fullVersion<<15) | (1<<12) | (smls<<1) | (iRob%2)), iDet, iRob, iMcm);
	       }
	     }
	   }
	 }

	 else if((cmd == fgkScsnCmdRead) ||
		 (cmd == fgkScsnCmdPause) ||
		 (cmd == fgkScsnCmdPtrg) ||
		 (cmd == fgkScsnCmdHwPtrg) ||
		 (cmd == fgkScsnCmdRobPower) ||
		 (cmd == fgkScsnCmdTtcRx) ||
		 (cmd == fgkScsnCmdMcmTemp) ||
		 (cmd == fgkScsnCmdOri) ||
		 (cmd == fgkScsnCmdPM) ) {
	   AliDebug(2, Form("ignored SCSN command: %i %i %i %i", cmd, addr, data, extali));
	 }

	 else {
            AliWarning(Form("unknown SCSN command: %i %i %i %i", cmd, addr, data, extali));
	    ignoredCmds++;
	 }

	 readLines++;
      }

      else if(!infile.eof() && !infile.good()) {
	 infile.clear();
	 infile.ignore(256, '\n');
	 ignoredLines++;
      }

      if(!infile.eof())
	 infile.clear();
   }

   infile.close();

   AliDebug(5, Form("Ignored lines: %i, ignored cmds: %i", ignoredLines, ignoredCmds));


   if(ignoredLines>readLines)
      AliError(Form("More than 50 %% of the input file could not be processed. Perhaps you should check the input file %s", filename.Data()));


   return kTRUE;
}



Int_t AliTRDtrapConfigHandler::SetGaintable(AliTRDCalOnlineGainTable const &gtbl)
{
   fGtbl=gtbl;
   return 0;
}


void AliTRDtrapConfigHandler::ProcessLTUparam(Int_t dest, Int_t addr, UInt_t data)
{
   //
   // Process the LTU parameters and stores them in internal class variables
   // or transfer the stored values to AliTRDtrapConfig, depending on the dest parameter
   //

   switch (dest) {

   case 0: // set the parameters in AliTRDtrapConfig
      for(Int_t det=0; det<AliTRDgeometry::Ndet(); det++) {
	 ConfigureDyCorr(det);
	 ConfigureDRange(det); // deflection range
	 ConfigureNTimebins(det);  // timebins in the drift region
	 ConfigurePIDcorr(det);  // scaling parameters for the PID
      }
      break;

   case 1: // set variables
      switch (addr) {

      case 0: ltuParam.SetPtMin(data); break; // pt_min in GeV/c (*1000)
      case 1: ltuParam.SetMagField(data); break; // B in T (*1000)
      case 2: ltuParam.SetOmegaTau(data); break; // omega*tau
      case 3: ltuParam.SetNtimebins(data); break;
	// ntimbins: drift time (for 3 cm) in timebins (5 add. bin. digits)
      case 4: ltuParam.SetScaleQ0(data); break;
      case 5: ltuParam.SetScaleQ1(data); break;
      case 6: ltuParam.SetLengthCorrectionEnable(data); break;
      case 7: ltuParam.SetTiltCorrectionEnable(data); break;
      }
      break;

   default:
      AliError(Form("dest %i not implemented", dest));
   }

}


void AliTRDtrapConfigHandler::ConfigureNTimebins(Int_t det)
{
   //
   // Set timebins in the drift region
   //

  if (!fTrapConfig) {
    AliError("No TRAPconfig given");
    return;
  }

  AddValues(det, fgkScsnCmdWrite, 127, AliTRDmcmSim::fgkDmemAddrNdrift, ltuParam.GetNtimebins());
}



void AliTRDtrapConfigHandler::ConfigureDyCorr(Int_t det)
{
   //
   //  Deflection length correction
   //  due to Lorentz angle and tilted pad correction
   //  This correction is in units of padwidth / (256*32)
   //

  if (!fTrapConfig) {
    AliError("No TRAPconfig given");
    return;
  }

   Int_t nRobs = AliTRDgeometry::GetStack(det) == 2 ? 6 : 8;

  for (Int_t r = 0; r < nRobs; r++) {
    for (Int_t m = 0; m < 16; m++) {
      Int_t dest =  1<<10 | r<<7 | m;
      Int_t dyCorrInt = ltuParam.GetDyCorrection(det, r, m);
      AddValues(det, fgkScsnCmdWrite, dest, AliTRDmcmSim::fgkDmemAddrDeflCorr, dyCorrInt);
    }
  }
}





void AliTRDtrapConfigHandler::ConfigureDRange(Int_t det)
{
   //
   // deflection range LUT
   // range calculated according to B-field (in T) and pt_min (in GeV/c)
   // if pt_min < 0.1 GeV/c the maximal allowed range for the tracklet
   // deflection (-64..63) is used
   //

  if (!fTrapConfig) {
    AliError("No TRAPconfig given");
    return;
  }

  Int_t nRobs = AliTRDgeometry::GetStack(det) == 2 ? 6 : 8;

  Int_t dyMinInt;
  Int_t dyMaxInt;

   for (Int_t r = 0; r < nRobs; r++) {
      for (Int_t m = 0; m < 16; m++) {
	 for (Int_t c = 0; c < 18; c++) {

	   // cout << "maxdefl: " << maxDeflAngle << ", localPhi " << localPhi << endl;
	   // cout << "r " << r << ", m" << m << ", c " << c << ", min angle: " << localPhi-maxDeflAngle << ", max: " << localPhi+maxDeflAngle
	   // 	<< ", min int: " << dyMinInt << ", max int: " << dyMaxInt << endl;
	   Int_t dest =  1<<10 | r<<7 | m;
	   Int_t lutAddr = AliTRDmcmSim::fgkDmemAddrDeflCutStart + 2*c;
	   ltuParam.GetDyRange(det, r, m, c, dyMinInt, dyMaxInt);
	   AddValues(det, fgkScsnCmdWrite, dest, lutAddr+0, dyMinInt);
	   AddValues(det, fgkScsnCmdWrite, dest, lutAddr+1, dyMaxInt);
	 }
      }
   }
}

void AliTRDtrapConfigHandler::PrintGeoTest()
{
   //
   // Prints some information about the geometry. Only for debugging
   //

   int sm=0;
   //   for(int sm=0; sm<6; sm++) {
   for(int stack=0; stack<5; stack++) {
      for(int layer=0; layer<6; layer++) {

	 Int_t det = sm*30+stack*6+layer;
	 for (Int_t r = 0; r < 6; r++) {
	    for (Int_t m = 0; m < 16; m++) {
	       for (Int_t c = 7; c < 8; c++) {
		 cout << stack << ";" << layer << ";" << r << ";" << m
		      << ";" << ltuParam.GetX(det, r, m)
		      << ";" << ltuParam.GetLocalY(det, r, m, c)
		      << ";" << ltuParam.GetLocalZ(det, r, m) << endl;
	       }
	    }
	 }
      }
   }
   // }
}


void AliTRDtrapConfigHandler::ConfigurePIDcorr(Int_t det)
{
   //
   // Calculate the MCM individual correction factors for the PID
   // and transfer them to AliTRDtrapConfig
   //

  if (!fTrapConfig) {
    AliError("No TRAPconfig given");
    return;
  }

   static const Int_t addrLUTcor0 = AliTRDmcmSim::fgkDmemAddrLUTcor0;
   static const Int_t addrLUTcor1 = AliTRDmcmSim::fgkDmemAddrLUTcor1;

   UInt_t cor0;
   UInt_t cor1;

   Int_t nRobs = AliTRDgeometry::GetStack(det) == 2 ? 6 : 8;

   for (Int_t r=0; r<nRobs; r++) {
      for(Int_t m=0; m<16; m++) {
	 Int_t dest =  1<<10 | r<<7 | m;
	 if(fGtbl.GetGainTableROC(det) && fGtbl.GetGainTableROC(det)->GetGainTableMCM(r, m))
	    ltuParam.GetCorrectionFactors(det, r, m, 9, cor0, cor1, fGtbl.GetGainTableROC(det)->GetGainTableMCM(r, m)->GetMCMGain());
	 else
	    ltuParam.GetCorrectionFactors(det, r, m, 9, cor0, cor1);
	 AddValues(det, fgkScsnCmdWrite, dest, addrLUTcor0, cor0);
	 AddValues(det, fgkScsnCmdWrite, dest, addrLUTcor1, cor1);
    }
  }
}


Bool_t AliTRDtrapConfigHandler::AddValues(UInt_t det, UInt_t cmd, UInt_t extali, Int_t addr, UInt_t data)
{
   // transfer the informations provided by LoadConfig to the internal class variables

  if (!fTrapConfig) {
    AliError("No TRAPconfig given");
    return kFALSE;
  }

  if(cmd != fgkScsnCmdWrite) {
    AliError(Form("Invalid command received: %i", cmd));
    return kFALSE;
  }

  AliTRDtrapConfig::TrapReg_t mcmReg = fTrapConfig->GetRegByAddress(addr);
  Int_t rocType = AliTRDgeometry::GetStack(det) == 2 ? 0 : 1;

  static const int mcmListSize=40;  // 40 is more or less arbitrary
  Int_t mcmList[mcmListSize];

  // configuration registers
  if(mcmReg >= 0 && mcmReg < AliTRDtrapConfig::kLastReg) {

    for(Int_t linkPair=0; linkPair<fgkMaxLinkPairs; linkPair++) {
      if(AliTRDfeeParam::ExtAliToAli(extali, linkPair, rocType, mcmList, mcmListSize)!=0) {
	Int_t i=0;
        while((i < mcmListSize) && (mcmList[i] != -1)) {
          if(mcmList[i]==127) {
	    AliDebug(1, Form("broadcast write to %s: 0x%08x",
			     fTrapConfig->GetRegName((AliTRDtrapConfig::TrapReg_t) mcmReg), data));
            fTrapConfig->SetTrapReg( (AliTRDtrapConfig::TrapReg_t) mcmReg, data, det);
	  }
          else {
	    AliDebug(1, Form("individual write to %s (%i, %i): 0x%08x",
			     fTrapConfig->GetRegName((AliTRDtrapConfig::TrapReg_t) mcmReg), (mcmList[i]>>7), (mcmList[i]&0x7F), data));
            fTrapConfig->SetTrapReg( (AliTRDtrapConfig::TrapReg_t) mcmReg, data, det, (mcmList[i]>>7)&0x7, (mcmList[i]&0x7F));
	  }
          i++;
        }
      }
    }
    return kTRUE;
  }
  // DMEM
  else if ( (addr >= AliTRDtrapConfig::fgkDmemStartAddress) &&
	    (addr < (AliTRDtrapConfig::fgkDmemStartAddress + AliTRDtrapConfig::fgkDmemWords))) {
    for(Int_t linkPair=0; linkPair<fgkMaxLinkPairs; linkPair++) {
      if(AliTRDfeeParam::ExtAliToAli(extali, linkPair, rocType, mcmList, mcmListSize)!=0) {
        Int_t i=0;
        while(mcmList[i] != -1 && i < mcmListSize) {
          if(mcmList[i] == 127)
	     fTrapConfig->SetDmem(addr, data, det, 0, 127);
          else
	     fTrapConfig->SetDmem(addr, data, det, mcmList[i] >> 7, mcmList[i] & 0x7f);
          i++;
        }
      }
    }
    return kTRUE;
  }
  else if ( (addr >= AliTRDtrapConfig::fgkImemStartAddress) &&
	    (addr < (AliTRDtrapConfig::fgkImemStartAddress + AliTRDtrapConfig::fgkImemWords))) {
    // IMEM is ignored for now
    return kTRUE;
  }
  else if ( (addr >= AliTRDtrapConfig::fgkDbankStartAddress) &&
	    (addr < (AliTRDtrapConfig::fgkDbankStartAddress + AliTRDtrapConfig::fgkImemWords))) {
    // DBANK is ignored for now
    return kTRUE;
  }
  else {
    AliError(Form("Writing to unhandled address 0x%04x", addr));
    return kFALSE;
  }
}
