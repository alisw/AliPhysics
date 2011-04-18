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
#include <iomanip>

#include "AliLog.h"

#include "AliTRDtrapConfig.h"
#include "AliTRDgeometry.h"
#include "AliTRDcalibDB.h"

#include "TMath.h"
#include "TGeoMatrix.h"
#include "TGraph.h"

using namespace std;

ClassImp(AliTRDtrapConfigHandler)

AliTRDtrapConfigHandler::AliTRDtrapConfigHandler() :
     ltuParam()
     , fRestrictiveMask((0x3ffff << 11) | (0x1f << 6) | 0x3f)
{

}


AliTRDtrapConfigHandler::~AliTRDtrapConfigHandler()
{

}

void AliTRDtrapConfigHandler::ResetMCMs()
{
   //
   // Reset all MCM registers and DMEM
   //

   AliTRDtrapConfig *cfg = AliTRDtrapConfig::Instance();
   cfg->ResetRegs();
   cfg->ResetDmem();
}


Int_t AliTRDtrapConfigHandler::LoadConfig()
{
  // load a default configuration which is suitable for simulation
  // for a detailed description of the registers see the TRAP manual
  // if you want to resimulate tracklets on real data use the appropriate config instead

  AliTRDtrapConfig *cfg = AliTRDtrapConfig::Instance();

  // HC header configuration bits
  cfg->SetTrapReg(AliTRDtrapConfig::kC15CPUA, 0x2102); // zs, deh

  // no. of timebins
  cfg->SetTrapReg(AliTRDtrapConfig::kC13CPUA, 24);

  // pedestal filter
  cfg->SetTrapReg(AliTRDtrapConfig::kFPNP, 4*10);
  cfg->SetTrapReg(AliTRDtrapConfig::kFPTC, 0);
  cfg->SetTrapReg(AliTRDtrapConfig::kFPBY, 0); // bypassed!

  // gain filter
  for (Int_t adc = 0; adc < 20; adc++) {
    cfg->SetTrapReg(AliTRDtrapConfig::TrapReg_t(AliTRDtrapConfig::kFGA0+adc), 40);
    cfg->SetTrapReg(AliTRDtrapConfig::TrapReg_t(AliTRDtrapConfig::kFGF0+adc), 15);
  }
  cfg->SetTrapReg(AliTRDtrapConfig::kFGTA, 20);
  cfg->SetTrapReg(AliTRDtrapConfig::kFGTB, 2060);
  cfg->SetTrapReg(AliTRDtrapConfig::kFGBY, 0);  // bypassed!

  // tail cancellation
  cfg->SetTrapReg(AliTRDtrapConfig::kFTAL, 200);
  cfg->SetTrapReg(AliTRDtrapConfig::kFTLL, 0);
  cfg->SetTrapReg(AliTRDtrapConfig::kFTLS, 200);
  cfg->SetTrapReg(AliTRDtrapConfig::kFTBY, 0);

  // tracklet calculation
  cfg->SetTrapReg(AliTRDtrapConfig::kTPQS0, 5);
  cfg->SetTrapReg(AliTRDtrapConfig::kTPQE0, 10);
  cfg->SetTrapReg(AliTRDtrapConfig::kTPQS1, 11);
  cfg->SetTrapReg(AliTRDtrapConfig::kTPQE1, 20);
  cfg->SetTrapReg(AliTRDtrapConfig::kTPFS, 5);
  cfg->SetTrapReg(AliTRDtrapConfig::kTPFE, 20);
  cfg->SetTrapReg(AliTRDtrapConfig::kTPVBY, 0);
  cfg->SetTrapReg(AliTRDtrapConfig::kTPVT, 10);
  cfg->SetTrapReg(AliTRDtrapConfig::kTPHT, 150);
  cfg->SetTrapReg(AliTRDtrapConfig::kTPFP, 40);
  cfg->SetTrapReg(AliTRDtrapConfig::kTPCL, 1);
  cfg->SetTrapReg(AliTRDtrapConfig::kTPCT, 10);

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

  // apply ltuParams to all detectors
  for(Int_t det=0; det<AliTRDgeometry::Ndet(); det++) {
     ConfigureDyCorr(det);
     ConfigureDRange(det); // deflection range
     ConfigureNTimebins(det);  // timebins in the drift region
     ConfigurePIDcorr(det);  // scaling parameters for the PID
  }

  // ****** hit position LUT

  // now calculate it from PRF
  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();

  Double_t padResponse[3]; // pad response left, central, right
  Double_t padResponseR[3]; // pad response left, central, right
  Double_t padResponseL[3]; // pad response left, central, right

  for (Int_t iBin = 0; iBin < 128; iBin++)
    cfg->SetTrapReg((AliTRDtrapConfig::TrapReg_t) (AliTRDtrapConfig::kTPL00 + iBin), 0, 0, 0, 0);

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
        cfg->SetTrapReg((AliTRDtrapConfig::TrapReg_t) (AliTRDtrapConfig::kTPL00 + iBin), corr, 6*iStack + iLayer);
      }
    }
  }
  // ****** hit position LUT configuration end

  // event buffer
  cfg->SetTrapReg(AliTRDtrapConfig::kEBSF, 1);  // 0: store filtered; 1: store unfiltered

  // zs applied to data stored in event buffer (sel. by EBSF)
  cfg->SetTrapReg(AliTRDtrapConfig::kEBIS, 15 << 2); // single indicator threshold (plus two digits)
  cfg->SetTrapReg(AliTRDtrapConfig::kEBIT, 30 << 2); // sum indicator threshold (plus two digits)
  cfg->SetTrapReg(AliTRDtrapConfig::kEBIL, 0xf0);   // lookup table
  cfg->SetTrapReg(AliTRDtrapConfig::kEBIN, 0);      // neighbour sensitivity

  // raw data
  cfg->SetTrapReg(AliTRDtrapConfig::kNES, (0x0000 << 16) | 0x1000);

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

   Int_t ignoredLines=0;
   Int_t ignoredCmds=0;
   Int_t readLines=0;


   AliTRDtrapConfig *cfg = AliTRDtrapConfig::Instance();

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

   while(infile.good()) {
      cmd=999;
      extali=-1;
      addr=-1;
      data=-1;
      infile >> std::skipws >> cmd >> addr >> data >> extali;

      if(cmd!=999 && extali!=-1 && addr != -1 && data!= -1 && extali!=-1) {

	 if(cmd==fgkScsnCmdWrite) {
	    for(Int_t det=0; det<AliTRDgeometry::Ndet(); det++) {
	       UInt_t rocpos = (1 << (AliTRDgeometry::GetSector(det)+11)) | (1 << (AliTRDgeometry::GetStack(det)+6)) | (1 << AliTRDgeometry::GetLayer(det));
	       if ((fRestrictiveMask & rocpos) == rocpos)
		  cfg->AddValues(det, cmd, extali, addr, data);
	    }
	 }

	 else if(cmd == fgkScsnLTUparam) {
	    ProcessLTUparam(extali, addr, data);
	 }

	 else if(cmd == fgkScsnCmdRestr) {
	    fRestrictiveMask = data;
	 }

	 else
	    ignoredCmds++;

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
  AliTRDtrapConfig::Instance()->AddValues(det, fgkScsnCmdWrite, 127, AliTRDtrapConfig::fgkDmemAddrNdrift, ltuParam.GetNtimebins());
}



void AliTRDtrapConfigHandler::ConfigureDyCorr(Int_t det)
{
   //
   //  Deflection length correction
   //  due to Lorentz angle and tilted pad correction
   //  This correction is in units of padwidth / (256*32)
   //

   Int_t nRobs = AliTRDgeometry::GetStack(det) == 2 ? 6 : 8;

  for (Int_t r = 0; r < nRobs; r++) {
    for (Int_t m = 0; m < 16; m++) {
      Int_t dest =  1<<10 | r<<7 | m;
      Int_t dyCorrInt = ltuParam.GetDyCorrection(det, r, m);
      AliTRDtrapConfig::Instance()->AddValues(det, fgkScsnCmdWrite, dest, AliTRDtrapConfig::fgkDmemAddrDeflCorr, dyCorrInt);
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
	   Int_t lutAddr = AliTRDtrapConfig::fgkDmemAddrDeflCutStart + 2*c;
	   ltuParam.GetDyRange(det, r, m, c, dyMinInt, dyMaxInt);
	   AliTRDtrapConfig::Instance()->AddValues(det, fgkScsnCmdWrite, dest, lutAddr+0, dyMinInt);
	   AliTRDtrapConfig::Instance()->AddValues(det, fgkScsnCmdWrite, dest, lutAddr+1, dyMaxInt);
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

   static const Int_t addrLUTcor0 = AliTRDtrapConfig::fgkDmemAddrLUTcor0;
   static const Int_t addrLUTcor1 = AliTRDtrapConfig::fgkDmemAddrLUTcor1;

   UInt_t cor0;
   UInt_t cor1;

   Int_t nRobs = AliTRDgeometry::GetStack(det) == 2 ? 6 : 8;

   for (Int_t r=0; r<nRobs; r++) {
      for(Int_t m=0; m<16; m++) {
	 Int_t dest =  1<<10 | r<<7 | m;
	 ltuParam.GetCorrectionFactors(det, r, m, 9, cor0, cor1);
	 AliTRDtrapConfig::Instance()->AddValues(det, fgkScsnCmdWrite, dest, addrLUTcor0, cor0);
	 AliTRDtrapConfig::Instance()->AddValues(det, fgkScsnCmdWrite, dest, addrLUTcor1, cor1);
    }
  }
}
