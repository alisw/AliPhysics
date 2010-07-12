
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

#include "AliTRDpadPlane.h"
#include "AliTRDtrapConfig.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"
#include "AliLog.h"

#include "AliTRDarrayDictionary.h"

#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliMagWrapCheb.h"
#include "AliTRDgeometry.h"

#include "TMath.h"
#include "TGeoMatrix.h"

using namespace std;

ClassImp(AliTRDtrapConfigHandler)


AliTRDtrapConfigHandler::AliTRDtrapConfigHandler() :
     fGeo(NULL)
     , fDet(0)
     , fBField(0)
     , fOmegaTau(1)
     , fPtMin(0)
     , fNTimebins(0)
     , fScaleQ0(1)
     , fScaleQ1(1)
     , fPidTracklengthCorr(kFALSE)
     , fTiltCorr(kTRUE)
{
   fGeo = new AliTRDgeometry;
}


AliTRDtrapConfigHandler::~AliTRDtrapConfigHandler()
{
   delete fGeo;
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


Int_t AliTRDtrapConfigHandler::LoadConfig(TString filename, Int_t det)
{
   //
  // load a TRAP configuration from a file
   // The file format is the format created by the standalone
   // command coder scc / show_cfdat but without the running number
   // scc /show_cfdat adds as the first column
   // which are two tools to inspect/export configurations from wingDB
   //


   fDet = det;
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

   while(infile.good()) {
      cmd=999;
      extali=-1;
      addr=-1;
      data=-1;
      infile >> std::skipws >> cmd >> addr >> data >> extali;
      //      std::cout << "no: " << no << ", cmd " << cmd << ", extali " << extali << ", addr " << addr << ", data " << data <<  endl;

      if(cmd!=999 && extali!=-1 && addr != -1 && data!= -1 && extali!=-1) {
	 if(cmd==fgkScsnCmdWrite)
	    cfg->AddValues(det, cmd, extali, addr, data);
	 else if(cmd == fgkScsnLTUparam)
	    ProcessLTUparam(extali, addr, data);
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
      ConfigureDyCorr();
      ConfigureDRange(); // deflection range
      ConfigureNTimebins();  // timebins in the drift region
      ConfigurePIDcorr();  // scaling parameters for the PID
      break;

   case 1: // set variables
      switch (addr) {

      case 0: fPtMin =  float(data) / 1000.; break; // pt_min in GeV/c (*1000)
      case 1: fBField = float(data) / 1000. ; break; // B in T (*1000)
      case 2: fOmegaTau = float(data) / 1.0E6  ; break; // omega*tau
      case 3: fNTimebins = data; break;
	// ntimbins: drift time (for 3 cm) in timebins (5 add. bin. digits)
      case 4: fScaleQ0 = data; break;
      case 5: fScaleQ1 = data; break;
      case 6: fPidTracklengthCorr = (bool) data; break;
      case 7: fTiltCorr = (bool) data; break;
      }
      break;

   default:
      AliError(Form("dest %i not implemented", dest));
   }

}


void AliTRDtrapConfigHandler::ConfigureNTimebins()
{
   //
   // Set timebins in the drift region
   //
   AliTRDtrapConfig::Instance()->AddValues(fDet, fgkScsnCmdWrite, 127, AliTRDtrapConfig::fgkDmemAddrNdrift, fNTimebins);
}



void AliTRDtrapConfigHandler::ConfigureDyCorr()
{
   //
   //  Deflection length correction
   //  due to Lorentz angle and tilted pad correction
   //  This correction is in units of padwidth / (256*32)
   //

   Int_t nRobs;
   if(fGeo->GetStack(fDet) == 2)
      nRobs=6;
   else
      nRobs=8;

   Double_t dyLorentz = - fOmegaTau * fGeo->CdrHght();    // CdrHght: Height of the drift region
   Double_t globalPos[3];
   Double_t tiltingAngle = fGeo->GetPadPlane(fDet)->GetTiltingAngle();
   Double_t scalePad = 256. * 32.;
   Int_t dyCorrInt=0;

    for (Int_t r = 0; r < nRobs; r++) {
	for (Int_t m = 0; m < 16; m++) {
	   if(GetPadPosNonRot(r, m, 9, globalPos)==0) {
	      double dyTilt = ( fGeo->CdrHght() *  TMath::Tan(tiltingAngle * TMath::DegToRad()) * globalPos[2]/globalPos[0]);

	      double dyCorr;
	      if(fTiltCorr==kTRUE)
		 dyCorr = dyLorentz + dyTilt;
	      else
		 dyCorr = dyTilt;


	      dyCorrInt =  TMath::Nint(dyCorr * scalePad /  fGeo->GetPadPlane(fDet)->GetWidthIPad());  // The correction is in units of 1/256 of the
	                                                                                               // pad width, including 5 binary digits
	   }
	   else {
	      AliError("No transformation matrix available");
	      dyCorrInt=0;
	   }
	   Int_t dest =  1<<10 | r<<7 | m;
	   AliTRDtrapConfig::Instance()->AddValues(fDet, fgkScsnCmdWrite, dest, AliTRDtrapConfig::fgkDmemAddrDeflCorr, dyCorrInt);
	}
    }
}





void AliTRDtrapConfigHandler::ConfigureDRange()
{
   //
   // deflection range LUT
   // range calculated according to B-field (in T) and pt_min (in GeV/c)
   // if pt_min < 0.1 GeV/c the maximal allowed range for the tracklet
   // deflection (-64..63) is used
   //

   static const int x=0;
   static const int y=1;
   static const Double_t dyBin = 140e-6;

   Int_t dyMinInt=fgkDyMinCut;
   Int_t dyMaxInt=fgkDyMaxCut;
   Double_t mcmPos[3];

   Int_t nRobs=-1;
   if(fGeo->GetStack(fDet) == 2)
      nRobs=6;
   else
      nRobs=8;

   for (Int_t r = 0; r < nRobs; r++) {
      for (Int_t m = 0; m < 16; m++) {
	 for (Int_t c = 0; c < 18; c++) {

	    if(fPtMin<0.1) {
	       dyMinInt=fgkDyMinCut;
	       dyMaxInt=fgkDyMaxCut;
	    }
	    else {
	       if(GetPadPosNonRot(r, m, c, mcmPos)==0) {
		  Double_t radius = fPtMin/(0.3*fBField);

		  double vertexPos[2] = {0,0};

		  Double_t distanceX = (vertexPos[x]-mcmPos[x]) / 100.; // cm -> m
		  Double_t distanceY = (vertexPos[y]-mcmPos[y]) / 100.; // cm -> m

		  Double_t maxDeflTemp = (TMath::Sqrt( Square(distanceX) + Square(distanceY)) / 2) / radius;
		  Double_t localPhi = TMath::ATan2(distanceY, distanceX);

		  Double_t maxDeflAngle=0;
		  if(maxDeflTemp < 1. ) {
		     maxDeflAngle = TMath::ASin(maxDeflTemp);
		     Double_t dyMin = fGeo->CdrHght()/100. * TMath::Tan(localPhi - maxDeflAngle);  // CdrHght: Height of the drift region in cm
		     Double_t dyMax = fGeo->CdrHght()/100. * TMath::Tan(localPhi + maxDeflAngle);

		     dyMinInt = Int_t (dyMin / dyBin);
		     dyMaxInt = Int_t (dyMax / dyBin);

		     if(dyMinInt < fgkDyMinCut)
			dyMinInt = fgkDyMinCut;
		     if(dyMaxInt > fgkDyMaxCut)
			dyMaxInt = fgkDyMaxCut;
		  }
		  else {
		     dyMinInt = fgkDyMinCut;
		     dyMaxInt = fgkDyMaxCut;
		  }

// 		  cout << "maxdefl: " << maxDeflAngle << ", localPhi " << localPhi << endl;
// 		  cout << "r " << r << ", m" << m << ", c " << c << ", min angle: " << localPhi-maxDeflAngle << ", max: " << localPhi+maxDeflAngle
// 		       << ", min int: " << dyMinInt << ", max int: " << dyMaxInt << endl;
	       }
	       else {
		  AliError("No geometry model loaded\n");
	       }
	    }

	    Int_t dest =  1<<10 | r<<7 | m;
 	    Int_t lutAddr = AliTRDtrapConfig::fgkDmemAddrDeflCutStart + 2*c;
	    AliTRDtrapConfig::Instance()->AddValues(fDet, fgkScsnCmdWrite, dest, lutAddr+0, dyMinInt);
	    AliTRDtrapConfig::Instance()->AddValues(fDet, fgkScsnCmdWrite, dest, lutAddr+1, dyMaxInt);
	 }
      }
   }
}

void AliTRDtrapConfigHandler::PrintGeoTest()
{
   //
   // Prints some information about the geometry. Only for debugging
   //

   Double_t mcmPos[3];
   int sm=0;
   //   for(int sm=0; sm<6; sm++) {
   for(int stack=0; stack<5; stack++) {
      for(int layer=0; layer<6; layer++) {

	 fDet = sm*30+stack*6+layer;
	 for (Int_t r = 0; r < 6; r++) {
	    for (Int_t m = 0; m < 16; m++) {
	       for (Int_t c = 7; c < 8; c++) {
		  GetPadPosNonRot(r, m, c, mcmPos);
		  cout << stack << ";" << layer << ";" << r << ";" << m
		       << ";" << mcmPos[0] << ";" << mcmPos[1] << ";" << mcmPos[2] << endl;
	       }
	    }
	 }
      }
   }
   // }
}



void AliTRDtrapConfigHandler::ConfigurePIDcorr()
{
   //
   // Calculate the MCM individual correction factors for the PID
   // and transfer them to AliTRDtrapConfig
   //

   static const Int_t addrLUTcor0 = AliTRDtrapConfig::fgkDmemAddrLUTcor0;
   static const Int_t addrLUTcor1 = AliTRDtrapConfig::fgkDmemAddrLUTcor1;

   UInt_t cor0;
   UInt_t cor1;

   Double_t globalPos[3];

   Int_t nRobs=-1;
   if(fGeo->GetStack(fDet) == 2)
      nRobs=6;
   else
      nRobs=8;

   for (Int_t r=0; r<nRobs; r++) {
      for(Int_t m=0; m<16; m++) {

	 if(GetPadPosNonRot(r, m, 9, globalPos)==0) {
	    Double_t elongation = TMath::Abs(TMath::Sqrt(globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1] + globalPos[2]*globalPos[2]) / globalPos[0]);

	    if(fPidTracklengthCorr==kFALSE) {
	       cor0 = fScaleQ0;
	       cor1 = fScaleQ1;
	    }
	    else {
	       cor0 = Int_t ((1.0*fScaleQ0* (1/elongation) ));
	       cor1 = Int_t ((1.0*fScaleQ1* (1/elongation) ));
	    }

	    Int_t dest =  1<<10 | r<<7 | m;
	    AliTRDtrapConfig::Instance()->AddValues(fDet, fgkScsnCmdWrite, dest, addrLUTcor0, cor0);
	    AliTRDtrapConfig::Instance()->AddValues(fDet, fgkScsnCmdWrite, dest, addrLUTcor1, cor1);
	 }
	 else {
	    AliError("No transformation matrix available");
	 }
      }
   }
}



Int_t AliTRDtrapConfigHandler::GetPadPosNonRot(Int_t rob, Int_t mcm, Int_t channel, Double_t trackCoor[3])
{
   //
   // Calcutate the gobal coordinates for an mcm channel in the supermodule at position -0.5
   //

   Int_t stack = fGeo->GetStack(fDet);
   Int_t layer = fGeo->GetLayer(fDet);

   AliTRDpadPlane *plane = fGeo->GetPadPlane(layer, stack);
   if(plane==NULL) {
      AliError(Form("stack %i, layer %i, det %i", stack, layer, fDet));
      return 1;
   }

   Double_t locYZ[2];
   Double_t loc[3];

   Double_t radialPos = fGeo->AnodePos()-0.83; //
   //   Double_t radialPos = 300.65 + 12.60 * layer; // cm // taken from libTRD/geometry/geometry.cc, probably not the final value

   GetLocalPadPos(plane, rob, mcm, channel, locYZ);

   loc[0] = radialPos;
   loc[1] = locYZ[0];
   loc[2] = locYZ[1];

   //transform from loc[3] - coordinates of the pad
   // Go to tracking coordinates

   TGeoHMatrix *fMatrix  = fGeo->GetClusterMatrix(fDet);
   if(fMatrix==NULL) {
      AliError(Form("stack %i, layer %i, det %i", stack, layer, fDet));
      return 2;
   }
   fMatrix->LocalToMaster(loc, trackCoor);
   return 0;
}


void AliTRDtrapConfigHandler::GetLocalPadPos(AliTRDpadPlane *plane, Int_t rob, Int_t mcm, Int_t channel, Double_t result[2])
{
   //
   // calculate the local coordinates for an mcm channel
   //

   Double_t localY, localZ;

    Int_t padCol;
    if(rob%2 == 0) //side a
       padCol  = (mcm % fgkMCMperROBCol) * fgkPadsPerMCM + channel;
    else
       padCol  = (mcm % fgkMCMperROBCol) * fgkPadsPerMCM + (plane->GetNcols()/2) + channel;

    Int_t padRow = ((Int_t) floor(rob/2.0)) * fgkMCMperROBRow + ((Int_t) floor(mcm/4));

    if(padCol<0 || padCol>= plane->GetNcols())
       AliError(Form("Invalid pad col: %i\n", padCol));

    if(padRow<0 || padRow>= plane->GetNrows())
       AliError(Form("Invalid pad row: %i\n", padRow));

    if(padCol+1 == plane->GetNcols()) // last pad
       localY = plane->GetColPos(padCol) + (plane->GetColEnd()-plane->GetColPos(padCol))/2;
    else
       localY = plane->GetColPos(padCol) + (plane->GetColPos(padCol+1)-plane->GetColPos(padCol))/2;

    if(padRow+1 == plane->GetNrows())
       localZ = plane->GetRowPosROC(padRow) + (plane->GetRowEndROC() - plane->GetRowPosROC(padRow))/2;
    else
       localZ = plane->GetRowPosROC(padRow) + (plane->GetRowPosROC(padRow+1) - plane->GetRowPosROC(padRow))/2;

    //    std::cout << "pad col " << padCol << ", pad row " << padRow << std::endl;
    //    std::cout << "pos y (col) " << localY << ", pos z (row) " << localZ << std::endl;

    result[0]=localY;
    result[1]=localZ;
}


Double_t AliTRDtrapConfigHandler::Square(Double_t val)
{
   //
   // calculate the square of the argument
   //

   return val*val;
}
