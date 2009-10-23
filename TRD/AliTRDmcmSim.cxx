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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD MCM (Multi Chip Module) simulator                                    //
//  which simulated the TRAP processing after the AD-conversion              //
//  The relevant parameters (i.e. configuration registers of the TRAP        //
//  configuration are taken from AliTRDtrapConfig.                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <fstream>  // needed for raw data dump

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TLine.h>
#include <TMath.h>
#include <TRandom.h>
#include <TClonesArray.h>

#include "AliLog.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliTRDdigit.h"

#include "AliTRDfeeParam.h"
#include "AliTRDtrapConfig.h"
#include "AliTRDSimParam.h"
#include "AliTRDgeometry.h"
#include "AliTRDcalibDB.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDarrayADC.h"
#include "AliTRDarrayDictionary.h"
#include "AliTRDpadPlane.h"
#include "AliTRDtrackletMCM.h"
#include "AliTRDmcmSim.h"

#include "AliMagF.h"
#include "TGeoGlobalMagField.h"

ClassImp(AliTRDmcmSim)

Bool_t AliTRDmcmSim::fgApplyCut = kTRUE;

Float_t AliTRDmcmSim::fgChargeNorm = 65000.;
Int_t AliTRDmcmSim::fgAddBaseline = 0;

Int_t AliTRDmcmSim::fgPidNBinsQ0 = 40;
Int_t AliTRDmcmSim::fgPidNBinsQ1 = 50;
Bool_t AliTRDmcmSim::fgPidLutDelete = kFALSE;
Int_t AliTRDmcmSim::fgPidLutDefault[40][50] = {
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  }, 
  { 0, 9, 6, 12, 29, 53, 76, 94, 107, 116, 122, 126, 128, 129, 129, 129, 128, 127, 126, 124, 122, 120, 117, 115, 112, 109, 107, 104, 101, 99, 96, 94, 91, 89, 87, 85, 83, 81, 79, 78, 77, 75, 74, 73, 72, 72, 71, 71, 70, 70  }, 
  { 0, 14, 8, 17, 37, 66, 94, 116, 131, 140, 146, 150, 152, 153, 153, 152, 150, 148, 145, 143, 139, 136, 132, 129, 125, 121, 118, 114, 110, 107, 104, 101, 98, 95, 93, 91, 89, 87, 85, 83, 82, 81, 80, 79, 78, 77, 77, 76, 76, 75  }, 
  { 0, 33, 19, 34, 69, 112, 145, 167, 181, 189, 194, 196, 197, 197, 196, 194, 191, 188, 184, 180, 175, 170, 165, 159, 154, 148, 143, 137, 132, 127, 123, 118, 114, 111, 107, 104, 101, 99, 96, 94, 92, 91, 89, 88, 87, 86, 85, 85, 84, 84  }, 
  { 0, 82, 52, 83, 136, 180, 205, 218, 226, 230, 232, 233, 233, 233, 232, 230, 228, 226, 223, 219, 215, 210, 205, 199, 193, 187, 180, 173, 167, 160, 154, 148, 142, 136, 131, 127, 122, 119, 115, 112, 109, 106, 104, 102, 100, 99, 97, 96, 95, 94  }, 
  { 0, 132, 96, 136, 185, 216, 231, 238, 242, 244, 245, 245, 245, 245, 245, 244, 243, 242, 240, 238, 236, 233, 230, 226, 222, 217, 212, 206, 200, 193, 187, 180, 173, 167, 161, 155, 149, 144, 139, 134, 130, 126, 123, 120, 117, 114, 112, 110, 108, 107  }, 
  { 0, 153, 120, 160, 203, 227, 238, 243, 246, 247, 248, 249, 249, 249, 248, 248, 247, 246, 245, 244, 243, 241, 239, 237, 234, 231, 228, 224, 219, 215, 209, 204, 198, 192, 186, 180, 174, 168, 163, 157, 152, 147, 143, 139, 135, 131, 128, 125, 123, 120  }, 
  { 0, 156, 128, 166, 207, 229, 239, 244, 247, 248, 249, 249, 249, 249, 249, 249, 248, 247, 247, 246, 244, 243, 242, 240, 238, 236, 233, 230, 227, 224, 220, 216, 212, 207, 202, 197, 192, 187, 181, 176, 171, 166, 161, 156, 152, 148, 144, 140, 137, 134  }, 
  { 0, 152, 128, 166, 206, 228, 239, 244, 246, 248, 249, 249, 249, 249, 249, 248, 248, 247, 246, 245, 244, 243, 241, 240, 238, 236, 234, 232, 229, 226, 224, 220, 217, 214, 210, 206, 202, 197, 193, 188, 184, 179, 174, 170, 166, 161, 157, 153, 150, 146  }, 
  { 0, 146, 126, 164, 203, 226, 237, 243, 246, 247, 248, 248, 248, 248, 248, 247, 247, 246, 245, 244, 242, 241, 239, 238, 236, 234, 232, 230, 227, 225, 223, 220, 217, 215, 212, 209, 205, 202, 199, 195, 191, 187, 183, 179, 175, 171, 168, 164, 160, 156  }, 
  { 0, 140, 123, 160, 200, 224, 235, 241, 244, 246, 247, 247, 247, 247, 247, 246, 245, 244, 243, 242, 240, 238, 237, 235, 233, 230, 228, 226, 224, 221, 219, 217, 215, 212, 210, 207, 205, 202, 200, 197, 194, 191, 188, 184, 181, 178, 174, 171, 168, 164  }, 
  { 0, 133, 119, 156, 196, 220, 233, 239, 243, 245, 245, 246, 246, 246, 245, 244, 243, 242, 241, 239, 237, 235, 233, 231, 229, 226, 224, 221, 219, 216, 214, 212, 210, 208, 206, 204, 202, 199, 197, 195, 193, 191, 188, 186, 183, 181, 178, 175, 172, 169  }, 
  { 0, 127, 115, 152, 192, 217, 230, 237, 241, 243, 244, 244, 244, 244, 243, 242, 241, 240, 238, 236, 234, 232, 229, 227, 224, 221, 218, 216, 213, 210, 208, 206, 203, 201, 200, 198, 196, 194, 193, 191, 190, 188, 186, 185, 183, 181, 179, 177, 174, 172  }, 
  { 0, 121, 111, 147, 187, 213, 227, 235, 239, 241, 242, 243, 243, 242, 241, 240, 239, 237, 236, 233, 231, 228, 225, 222, 219, 216, 213, 210, 207, 204, 201, 199, 196, 194, 192, 191, 189, 188, 187, 185, 184, 183, 182, 181, 180, 178, 177, 176, 174, 172  }, 
  { 0, 116, 107, 142, 181, 209, 224, 232, 237, 239, 240, 241, 241, 240, 239, 238, 237, 235, 233, 230, 227, 224, 221, 218, 214, 211, 207, 204, 200, 197, 194, 191, 189, 187, 185, 183, 182, 180, 179, 178, 178, 177, 176, 175, 175, 174, 173, 172, 172, 170  }, 
  { 0, 112, 103, 136, 176, 204, 220, 229, 234, 237, 238, 239, 239, 238, 237, 236, 234, 232, 230, 227, 224, 221, 217, 213, 209, 205, 201, 198, 194, 190, 187, 184, 181, 179, 177, 175, 174, 172, 171, 171, 170, 169, 169, 169, 168, 168, 168, 168, 167, 167  }, 
  { 0, 107, 99, 131, 170, 199, 216, 226, 231, 234, 236, 237, 237, 236, 235, 234, 232, 230, 227, 224, 221, 217, 213, 209, 205, 200, 196, 192, 188, 184, 180, 177, 174, 172, 169, 167, 166, 164, 163, 162, 162, 161, 161, 161, 161, 161, 161, 162, 162, 162  }, 
  { 0, 104, 94, 125, 164, 193, 212, 222, 228, 232, 233, 234, 234, 234, 233, 231, 229, 227, 224, 221, 218, 214, 210, 205, 201, 196, 191, 187, 182, 178, 174, 171, 168, 165, 162, 160, 158, 157, 155, 154, 154, 153, 153, 153, 153, 154, 154, 154, 155, 155  }, 
  { 0, 100, 90, 119, 157, 188, 207, 219, 225, 229, 231, 232, 232, 231, 230, 229, 227, 224, 222, 218, 215, 211, 206, 202, 197, 192, 187, 182, 178, 173, 169, 165, 162, 158, 156, 153, 151, 149, 148, 147, 146, 146, 145, 145, 145, 146, 146, 147, 148, 148  }, 
  { 0, 97, 86, 113, 150, 182, 202, 215, 222, 226, 228, 229, 229, 229, 228, 226, 224, 222, 219, 216, 212, 208, 203, 199, 194, 188, 183, 178, 173, 169, 164, 160, 156, 153, 150, 147, 145, 143, 141, 140, 139, 138, 138, 138, 138, 138, 139, 139, 140, 141  }, 
  { 0, 94, 82, 107, 144, 176, 197, 210, 218, 223, 225, 227, 227, 227, 226, 224, 222, 220, 217, 213, 209, 205, 201, 196, 191, 186, 180, 175, 170, 165, 160, 156, 152, 148, 145, 142, 139, 137, 135, 134, 132, 131, 131, 131, 131, 131, 131, 132, 132, 133  }, 
  { 0, 92, 78, 102, 137, 169, 192, 206, 215, 220, 223, 224, 224, 224, 223, 222, 220, 217, 215, 211, 207, 203, 199, 194, 188, 183, 178, 172, 167, 162, 157, 152, 148, 144, 140, 137, 134, 132, 130, 128, 127, 125, 125, 124, 124, 124, 124, 125, 125, 126  }, 
  { 0, 90, 75, 96, 131, 163, 187, 202, 211, 216, 220, 221, 222, 222, 221, 220, 218, 215, 212, 209, 205, 201, 197, 192, 187, 181, 176, 170, 165, 159, 154, 149, 145, 141, 137, 133, 130, 128, 125, 123, 122, 120, 119, 118, 118, 118, 118, 118, 119, 119  }, 
  { 0, 88, 71, 91, 124, 157, 181, 197, 207, 213, 217, 219, 219, 219, 219, 217, 216, 213, 211, 207, 204, 200, 195, 190, 185, 180, 174, 169, 163, 158, 152, 147, 142, 138, 134, 130, 127, 124, 121, 119, 117, 116, 114, 114, 113, 112, 112, 112, 112, 113  }, 
  { 0, 87, 68, 86, 118, 151, 176, 192, 203, 210, 214, 216, 217, 217, 217, 215, 214, 212, 209, 206, 202, 198, 194, 189, 184, 179, 173, 167, 162, 156, 151, 146, 141, 136, 132, 128, 124, 121, 118, 116, 114, 112, 110, 109, 108, 108, 107, 107, 107, 107  }, 
  { 0, 85, 65, 81, 112, 144, 170, 188, 199, 206, 211, 213, 214, 215, 214, 213, 212, 210, 207, 204, 201, 197, 193, 188, 183, 178, 172, 167, 161, 155, 150, 145, 140, 135, 130, 126, 122, 119, 116, 113, 111, 109, 107, 106, 105, 104, 103, 103, 102, 102  }, 
  { 0, 84, 62, 77, 106, 138, 165, 183, 195, 203, 208, 210, 212, 212, 212, 211, 210, 208, 206, 203, 200, 196, 192, 187, 183, 177, 172, 166, 161, 155, 150, 144, 139, 134, 129, 125, 121, 117, 114, 111, 109, 106, 104, 103, 101, 100, 99, 99, 98, 98  }, 
  { 0, 84, 60, 73, 101, 133, 159, 178, 191, 199, 204, 208, 209, 210, 210, 209, 208, 206, 204, 202, 199, 195, 191, 187, 182, 177, 172, 166, 161, 155, 150, 144, 139, 134, 129, 124, 120, 116, 113, 110, 107, 104, 102, 100, 99, 98, 96, 96, 95, 95  }, 
  { 0, 83, 58, 69, 96, 127, 154, 174, 187, 196, 201, 205, 207, 208, 208, 207, 206, 205, 203, 200, 197, 194, 190, 186, 182, 177, 172, 167, 161, 156, 150, 145, 139, 134, 129, 124, 120, 116, 112, 109, 106, 103, 101, 99, 97, 95, 94, 93, 92, 92  }, 
  { 0, 82, 56, 66, 91, 121, 149, 169, 183, 192, 198, 202, 204, 206, 206, 206, 205, 203, 201, 199, 196, 193, 190, 186, 182, 177, 172, 167, 162, 156, 151, 145, 140, 135, 129, 125, 120, 116, 112, 108, 105, 102, 100, 97, 95, 94, 92, 91, 90, 89  }, 
  { 0, 82, 54, 62, 86, 116, 144, 165, 179, 189, 195, 199, 202, 203, 204, 204, 203, 202, 200, 198, 196, 193, 189, 186, 182, 177, 173, 168, 163, 157, 152, 146, 141, 136, 130, 125, 121, 116, 112, 108, 105, 102, 99, 96, 94, 92, 91, 89, 88, 87  }, 
  { 0, 82, 52, 59, 82, 111, 139, 160, 175, 185, 192, 197, 200, 201, 202, 202, 201, 200, 199, 197, 195, 192, 189, 186, 182, 178, 173, 168, 163, 158, 153, 148, 142, 137, 132, 127, 122, 117, 113, 109, 105, 102, 99, 96, 94, 92, 90, 88, 87, 85  }, 
  { 0, 82, 50, 56, 78, 106, 134, 156, 171, 182, 189, 194, 197, 199, 200, 200, 200, 199, 198, 196, 194, 191, 188, 185, 182, 178, 174, 169, 164, 159, 154, 149, 144, 138, 133, 128, 123, 118, 114, 110, 106, 102, 99, 96, 93, 91, 89, 87, 86, 84  }, 
  { 0, 82, 49, 54, 74, 102, 129, 151, 167, 179, 186, 191, 195, 197, 198, 198, 198, 197, 196, 195, 193, 191, 188, 185, 182, 178, 174, 170, 165, 161, 156, 151, 145, 140, 135, 130, 125, 120, 115, 111, 107, 103, 100, 97, 94, 91, 89, 87, 85, 83  }, 
  { 0, 82, 47, 51, 70, 97, 124, 147, 164, 175, 183, 189, 192, 195, 196, 197, 197, 196, 195, 194, 192, 190, 188, 185, 182, 178, 175, 171, 166, 162, 157, 152, 147, 142, 137, 132, 127, 122, 117, 112, 108, 104, 101, 97, 94, 91, 89, 87, 85, 83  }, 
  { 0, 83, 46, 49, 67, 93, 120, 143, 160, 172, 180, 186, 190, 192, 194, 195, 195, 195, 194, 193, 191, 189, 187, 185, 182, 179, 175, 172, 167, 163, 159, 154, 149, 144, 139, 134, 129, 124, 119, 114, 110, 106, 102, 98, 95, 92, 89, 87, 85, 83  }, 
  { 0, 83, 45, 47, 64, 89, 116, 139, 156, 169, 177, 184, 188, 190, 192, 193, 193, 193, 193, 192, 190, 189, 187, 184, 182, 179, 176, 172, 168, 164, 160, 156, 151, 146, 141, 136, 131, 126, 121, 116, 112, 108, 104, 100, 96, 93, 90, 88, 85, 83  }, 
  { 0, 84, 44, 45, 61, 85, 111, 134, 152, 165, 175, 181, 185, 188, 190, 191, 192, 192, 191, 191, 189, 188, 186, 184, 182, 179, 176, 173, 169, 166, 162, 157, 153, 148, 143, 138, 133, 128, 124, 119, 114, 110, 106, 102, 98, 95, 91, 89, 86, 84  }, 
  { 0, 85, 43, 43, 58, 81, 107, 131, 149, 162, 172, 178, 183, 186, 188, 190, 190, 190, 190, 189, 188, 187, 186, 184, 182, 179, 176, 173, 170, 167, 163, 159, 155, 150, 145, 141, 136, 131, 126, 121, 117, 112, 108, 104, 100, 96, 93, 90, 87, 85  }, 
  { 0, 85, 42, 41, 55, 78, 103, 127, 145, 159, 169, 176, 181, 184, 186, 188, 189, 189, 189, 188, 188, 186, 185, 183, 181, 179, 177, 174, 171, 168, 164, 160, 156, 152, 148, 143, 138, 134, 129, 124, 119, 115, 110, 106, 102, 98, 95, 91, 88, 86  } 
};

Int_t (*AliTRDmcmSim::fgPidLut) = *fgPidLutDefault;

//_____________________________________________________________________________
AliTRDmcmSim::AliTRDmcmSim() : TObject()
  ,fInitialized(kFALSE)
  ,fMaxTracklets(-1) 
  ,fDetector(-1)
  ,fRobPos(-1)
  ,fMcmPos(-1)
  ,fRow (-1)
  ,fNADC(-1)
  ,fNTimeBin(-1)
  ,fADCR(NULL)
  ,fADCF(NULL)
  ,fMCMT(NULL)
  ,fTrackletArray(NULL)      
  ,fZSM(NULL)
  ,fZSM1Dim(NULL)
  ,fFeeParam(NULL)
  ,fTrapConfig(NULL)
  ,fSimParam(NULL)
  ,fCommonParam(NULL)
  ,fCal(NULL)
  ,fGeo(NULL)
  ,fDigitsManager(NULL)
  ,fPedAcc(NULL)
  ,fGainCounterA(NULL)
  ,fGainCounterB(NULL)
  ,fTailAmplLong(NULL)
  ,fTailAmplShort(NULL)
  ,fNHits(0)
  ,fFitReg(NULL)
{
  //
  // AliTRDmcmSim default constructor
  // By default, nothing is initialized.
  // It is necessary to issue Init before use.
}

AliTRDmcmSim::~AliTRDmcmSim() 
{
  //
  // AliTRDmcmSim destructor
  //

  if(fInitialized) {
    for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
      delete [] fADCR[iadc];
      delete [] fADCF[iadc];
      delete [] fZSM [iadc];
    }
    delete [] fADCR;
    delete [] fADCF;
    delete [] fZSM;
    delete [] fZSM1Dim;
    delete [] fMCMT;
 
    delete [] fPedAcc;
    delete [] fGainCounterA;
    delete [] fGainCounterB;
    delete [] fTailAmplLong;
    delete [] fTailAmplShort;
    delete [] fFitReg;
    
    fTrackletArray->Delete();
    delete fTrackletArray;
    delete fGeo;
  }
}

void AliTRDmcmSim::Init( Int_t det, Int_t robPos, Int_t mcmPos, Bool_t /* newEvent */ ) 
{
  //
  // Initialize the class with new geometry information
  // fADC array will be reused with filled by zero
  //
   
  if (!fInitialized) {
    fFeeParam      = AliTRDfeeParam::Instance();
    fTrapConfig    = AliTRDtrapConfig::Instance();
    fSimParam      = AliTRDSimParam::Instance();
    fCommonParam   = AliTRDCommonParam::Instance();
    fCal           = AliTRDcalibDB::Instance();
    fGeo           = new AliTRDgeometry();
  }

  fDetector      = det;
  fRobPos        = robPos;
  fMcmPos        = mcmPos;
  fNADC          = fFeeParam->GetNadcMcm();
  fNTimeBin      = fCal->GetNumberOfTimeBins();
  fRow           = fFeeParam->GetPadRowFromMCM( fRobPos, fMcmPos );
  fMaxTracklets  = fFeeParam->GetMaxNrOfTracklets();  
  
  if (!fInitialized) {
    fADCR    = new Int_t *[fNADC];
    fADCF    = new Int_t *[fNADC];
    fZSM     = new Int_t *[fNADC];
    fZSM1Dim = new Int_t  [fNADC];
    fGainCounterA = new UInt_t[fNADC];
    fGainCounterB = new UInt_t[fNADC];
    for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
      fADCR[iadc] = new Int_t[fNTimeBin];
      fADCF[iadc] = new Int_t[fNTimeBin];
      fZSM [iadc] = new Int_t[fNTimeBin];
    }
    
    // filter registers
    fPedAcc = new UInt_t[fNADC]; // accumulator for pedestal filter
    fTailAmplLong = new UShort_t[fNADC];
    fTailAmplShort = new UShort_t[fNADC];
    
    // tracklet calculation
    fFitReg = new FitReg_t[fNADC]; 
    fTrackletArray = new TClonesArray("AliTRDtrackletMCM", fMaxTracklets);
    
    fMCMT = new UInt_t[fMaxTracklets];
  }

  fInitialized = kTRUE;

  Reset();
}

void AliTRDmcmSim::Reset()
{
  // Resets the data values and internal filter registers
  // by re-initialising them

  for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
    for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {
      fADCR[iadc][it] = 0;
      fADCF[iadc][it] = 0;
      fZSM [iadc][it] = 1;   // Default unread = 1
    }
    fZSM1Dim[iadc] = 1;      // Default unread = 1
    fGainCounterA[iadc] = 0;
    fGainCounterB[iadc] = 0;
  }
  
  for(Int_t i = 0; i < fMaxTracklets; i++) {
    fMCMT[i] = 0;
  }
  
  FilterPedestalInit();
  FilterGainInit();
  FilterTailInit(fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFPNP)); //??? not really correct if gain filter is active
}

void AliTRDmcmSim::SetNTimebins(Int_t ntimebins) 
{
  fNTimeBin = ntimebins;
  for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
    delete fADCR[iadc];
    delete fADCF[iadc];
    delete fZSM[iadc];
    fADCR[iadc] = new Int_t[fNTimeBin];
    fADCF[iadc] = new Int_t[fNTimeBin];
    fZSM [iadc] = new Int_t[fNTimeBin];
  }
}

Bool_t AliTRDmcmSim::LoadMCM(AliRunLoader* const runloader, Int_t det, Int_t rob, Int_t mcm) 
{
  // loads the ADC data as obtained from the digitsManager for the specified MCM

  Init(det, rob, mcm);

  if (!runloader) {
    AliError("No Runloader given");
    return kFALSE;
  }

  AliLoader *trdLoader = runloader->GetLoader("TRDLoader");
  if (!trdLoader) {
    AliError("Could not get TRDLoader");
    return kFALSE;
  }

  Bool_t retval = kTRUE;
  trdLoader->LoadDigits();
  fDigitsManager = 0x0;
  AliTRDdigitsManager *digMgr = new AliTRDdigitsManager();
  digMgr->SetSDigits(0);
  digMgr->CreateArrays();
  digMgr->ReadDigits(trdLoader->TreeD());
  AliTRDarrayADC *digits = (AliTRDarrayADC*) digMgr->GetDigits(det);
  if (digits->HasData()) {
    digits->Expand();

    if (fNTimeBin != digits->GetNtime()) 
      SetNTimebins(digits->GetNtime());

    Int_t padrow = fFeeParam->GetPadRowFromMCM(rob, mcm);
    Int_t padcol = 0;
    for (Int_t ch = 0; ch < fNADC; ch++) {
      padcol = GetCol(ch);
      fZSM1Dim[ch] = 1;
      if (padcol < 0) {
        fZSM1Dim[ch] = 0;
        for (Int_t tb = 0; tb < fNTimeBin; tb++) {
          fADCR[ch][tb] = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFP) + (fgAddBaseline << fgkAddDigits);
          fADCF[ch][tb] = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFP) + (fgAddBaseline << fgkAddDigits);
        }
      }
      else {
        for (Int_t tb = 0; tb < fNTimeBin; tb++) {
          if (digits->GetData(padrow,padcol, tb) < 0) {
            fZSM1Dim[ch] = 0; 
            fADCR[ch][tb] = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFP) + (fgAddBaseline << fgkAddDigits);
            fADCF[ch][tb] = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFP) + (fgAddBaseline << fgkAddDigits);
          }
          else {
            fADCR[ch][tb] = (digits->GetData(padrow, padcol, tb) + fgAddBaseline) << fgkAddDigits;
            fADCF[ch][tb] = (digits->GetData(padrow, padcol, tb) + fgAddBaseline) << fgkAddDigits;
          }
        }
      }
    }
  }
  else 
    retval = kFALSE;
  
  delete digMgr;
  
  return retval;
}

void AliTRDmcmSim::NoiseTest(Int_t nsamples, Int_t mean, Int_t sigma, Int_t inputGain, Int_t inputTail)
{
  // This function can be used to test the filters. 
  // It feeds nsamples of ADC values with a gaussian distribution specified by mean and sigma.
  // The filter chain implemented here consists of:
  // Pedestal -> Gain -> Tail
  // With inputGain and inputTail the input to the gain and tail filter, respectively, 
  // can be chosen where 
  // 0: noise input
  // 1: pedestal output
  // 2: gain output
  // The input has to be chosen from a stage before. 
  // The filter behaviour is controlled by the TRAP parameters from AliTRDtrapConfig in the 
  // same way as in normal simulation.
  // The functions produces four histograms with the values at the different stages.

  TH1F *h   = new TH1F("noise", "Gaussian Noise;sample;ADC count",
                       nsamples, 0, nsamples);
  TH1F *hfp = new TH1F("pedf", "Noise #rightarrow Pedestal filter;sample;ADC count", nsamples, 0, nsamples);
  TH1F *hfg = new TH1F("pedg", "Pedestal #rightarrow Gain;sample;ADC count", nsamples, 0, nsamples);
  TH1F *hft = new TH1F("pedt", "Gain #rightarrow Tail;sample;ADC count", nsamples, 0, nsamples);
  h->SetStats(kFALSE);
  hfp->SetStats(kFALSE);
  hfg->SetStats(kFALSE);
  hft->SetStats(kFALSE);
  
  Int_t value;  // ADC count with noise (10 bit)
  Int_t valuep; // pedestal filter output (12 bit)
  Int_t valueg; // gain filter output (12 bit)
  Int_t valuet; // tail filter value (12 bit)
  
  for (Int_t i = 0; i < nsamples; i++) {
    value = (Int_t) gRandom->Gaus(mean, sigma);  // generate noise with gaussian distribution 
    h->SetBinContent(i, value);

    valuep = FilterPedestalNextSample(1, 0, ((Int_t) value) << 2);
    
    if (inputGain == 0)
      valueg = FilterGainNextSample(1, ((Int_t) value) << 2);
    else 
      valueg = FilterGainNextSample(1, valuep); 
    
    if (inputTail == 0)
      valuet = FilterTailNextSample(1, ((Int_t) value) << 2);
    else if (inputTail == 1)
      valuet = FilterTailNextSample(1, valuep); 
    else
      valuet = FilterTailNextSample(1, valueg); 

    hfp->SetBinContent(i, valuep >> 2);
    hfg->SetBinContent(i, valueg >> 2);
    hft->SetBinContent(i, valuet >> 2);
  }

  TCanvas *c = new TCanvas; 
  c->Divide(2,2);
  c->cd(1);
  h->Draw();
  c->cd(2);
  hfp->Draw();
  c->cd(3);
  hfg->Draw();
  c->cd(4);
  hft->Draw();
}

Bool_t AliTRDmcmSim::CheckInitialized()
{
  //
  // Check whether object is initialized
  //

  if( ! fInitialized ) {
    AliDebug(2, Form ("AliTRDmcmSim is not initialized but function other than Init() is called."));
  }
  return fInitialized;
}

void AliTRDmcmSim::Print(Option_t* const option) const
{
  // Prints the data stored and/or calculated for this MCM.
  // The output is controlled by option which can be a sequence of any of 
  // the following characters:
  // R - prints raw ADC data
  // F - prints filtered data 
  // H - prints detected hits
  // T - prints found tracklets
  // The later stages are only useful when the corresponding calculations 
  // have been performed.

  printf("MCM %i on ROB %i in detector %i\n", fMcmPos, fRobPos, fDetector);

  TString opt = option;
  if (opt.Contains("R")) {
    printf("Raw ADC data (10 bit):\n");
    for (Int_t iTimeBin = 0; iTimeBin < fNTimeBin; iTimeBin++) {
      for (Int_t iChannel = 0; iChannel < fNADC; iChannel++) {
        printf("%5i", fADCR[iChannel][iTimeBin] >> fgkAddDigits);
      }
      printf("\n");
    }
  }

  if (opt.Contains("F")) {
    printf("Filtered data (12 bit):\n");
    for (Int_t iTimeBin = 0; iTimeBin < fNTimeBin; iTimeBin++) {
      for (Int_t iChannel = 0; iChannel < fNADC; iChannel++) {
        printf("%5i", fADCF[iChannel][iTimeBin]);
      }
      printf("\n");
    }
  }

  if (opt.Contains("H")) {
    printf("Found %i hits:\n", fNHits);
    for (Int_t iHit = 0; iHit < fNHits; iHit++) {
      printf("Hit %3i in timebin %2i, ADC %2i has charge %3i and position %3i\n",
             iHit,  fHits[iHit].fTimebin, fHits[iHit].fChannel, fHits[iHit].fQtot, fHits[iHit].fYpos);
    }
  }

  if (opt.Contains("T")) {
    printf("Tracklets:\n");
    for (Int_t iTrkl = 0; iTrkl < fTrackletArray->GetEntriesFast(); iTrkl++) {
      printf("tracklet %i: 0x%08x\n", iTrkl, ((AliTRDtrackletMCM*) (*fTrackletArray)[iTrkl])->GetTrackletWord());
    }
  }
}

void AliTRDmcmSim::Draw(Option_t* const option) 
{
  // Plots the data stored in a 2-dim. timebin vs. ADC channel plot.
  // The option selects what data is plotted and can be a sequence of 
  // the following characters:
  // R - plot raw data (default)
  // F - plot filtered data (meaningless if R is specified)
  // In addition to the ADC values:
  // H - plot hits 
  // T - plot tracklets

  TString opt = option;

  TH2F *hist = new TH2F("mcmdata", Form("Data of MCM %i on ROB %i in detector %i", \
                                        fMcmPos, fRobPos, fDetector), \
                        fNADC, -0.5, fNADC-.5, fNTimeBin, -.5, fNTimeBin-.5);
  hist->GetXaxis()->SetTitle("ADC Channel");
  hist->GetYaxis()->SetTitle("Timebin");
  hist->SetStats(kFALSE);

  if (opt.Contains("R")) {
    for (Int_t iTimeBin = 0; iTimeBin < fNTimeBin; iTimeBin++) {
      for (Int_t iAdc = 0; iAdc < fNADC; iAdc++) {
        hist->SetBinContent(iAdc+1, iTimeBin+1, fADCR[iAdc][iTimeBin] >> fgkAddDigits);
      }
    }
  }
  else {
    for (Int_t iTimeBin = 0; iTimeBin < fNTimeBin; iTimeBin++) {
      for (Int_t iAdc = 0; iAdc < fNADC; iAdc++) {
        hist->SetBinContent(iAdc+1, iTimeBin+1, fADCF[iAdc][iTimeBin] >> fgkAddDigits);
      }
    }
  }
  hist->Draw("colz");

  if (opt.Contains("H")) {
    TGraph *grHits = new TGraph();
    for (Int_t iHit = 0; iHit < fNHits; iHit++) {
      grHits->SetPoint(iHit, 
                       fHits[iHit].fChannel + 1 + fHits[iHit].fYpos/256., 
                       fHits[iHit].fTimebin);
    }
    grHits->Draw("*");
  }

  if (opt.Contains("T")) {
    TLine *trklLines = new TLine[4];
    for (Int_t iTrkl = 0; iTrkl < fTrackletArray->GetEntries(); iTrkl++) {
      AliTRDpadPlane *pp = fGeo->GetPadPlane(fDetector);
      AliTRDtrackletMCM *trkl = (AliTRDtrackletMCM*) (*fTrackletArray)[iTrkl];
      Float_t offset = pp->GetColPos(fFeeParam->GetPadColFromADC(fRobPos, fMcmPos, 19)) + 19 * pp->GetWidthIPad();
      trklLines[iTrkl].SetX1((offset -  trkl->GetY()) / pp->GetWidthIPad());
      trklLines[iTrkl].SetY1(0);
      trklLines[iTrkl].SetX2((offset - (trkl->GetY() + ((Float_t) trkl->GetdY())*140e-4)) / pp->GetWidthIPad());
      trklLines[iTrkl].SetY2(fNTimeBin - 1);
      trklLines[iTrkl].SetLineColor(2);
      trklLines[iTrkl].SetLineWidth(2);
      printf("Tracklet %i: y = %f, dy = %f, offset = %f\n", iTrkl, trkl->GetY(), (trkl->GetdY() * 140e-4), offset);
      trklLines[iTrkl].Draw();
    }
  }
}

void AliTRDmcmSim::SetData( Int_t iadc, Int_t* const adc )
{
  //
  // Store ADC data into array of raw data
  //

  if( !CheckInitialized() ) return;

  if( iadc < 0 || iadc >= fNADC ) {
			//Log (Form ("Error: iadc is out of range (should be 0 to %d).", fNADC-1));
    return;
  }

  for( Int_t it = 0 ;  it < fNTimeBin ; it++ ) {
    fADCR[iadc][it] = (Int_t) (adc[it]) << fgkAddDigits;
    fADCF[iadc][it] = (Int_t) (adc[it]) << fgkAddDigits;
  }
}

void AliTRDmcmSim::SetData( Int_t iadc, Int_t it, Int_t adc )
{
  //
  // Store ADC data into array of raw data
  //

  if( !CheckInitialized() ) return;

  if( iadc < 0 || iadc >= fNADC ) {
    //Log (Form ("Error: iadc is out of range (should be 0 to %d).", fNADC-1));
    return;
  }

  fADCR[iadc][it] = adc << fgkAddDigits;
  fADCF[iadc][it] = adc << fgkAddDigits;
}

void AliTRDmcmSim::SetData(AliTRDarrayADC* const adcArray, AliTRDdigitsManager *digitsManager)
{
  // Set the ADC data from an AliTRDarrayADC

  if (!fInitialized) {
    AliError("Called uninitialized! Nothing done!");
    return;
  }

  fDigitsManager = digitsManager;

  if (fNTimeBin != adcArray->GetNtime())
    SetNTimebins(adcArray->GetNtime());

  Int_t offset = (fMcmPos % 4) * 21 + (fRobPos % 2) * 84;

//  Int_t firstAdc = 0;
//  Int_t lastAdc = fNADC-1;
//
//  while (GetCol(firstAdc) < 0) {
//    for (Int_t iTimeBin = 0; iTimeBin < fNTimeBin; iTimeBin++) {
//      fADCR[firstAdc][iTimeBin] = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFP) + (fgAddBaseline << fgkAddDigits);
//      fADCF[firstAdc][iTimeBin] = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFP) + (fgAddBaseline << fgkAddDigits);
//    }
//    firstAdc++;
//  }
//
//  while (GetCol(lastAdc) < 0) {
//    for (Int_t iTimeBin = 0; iTimeBin < fNTimeBin; iTimeBin++) {
//      fADCR[lastAdc][iTimeBin] = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFP) + (fgAddBaseline << fgkAddDigits);
//      fADCF[lastAdc][iTimeBin] = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFP) + (fgAddBaseline << fgkAddDigits);
//    }
//    lastAdc--;
//  }

  for (Int_t iTimeBin = 0; iTimeBin < fNTimeBin; iTimeBin++) {
    for (Int_t iAdc = 0; iAdc < fNADC; iAdc++) {
      Int_t value = adcArray->GetDataByAdcCol(GetRow(), 20-iAdc + offset, iTimeBin);
      if (value < 0 || (20-iAdc + offset < 1) || (20-iAdc + offset > 165)) {
        fADCR[iAdc][iTimeBin] = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFP) + (fgAddBaseline << fgkAddDigits);
        fADCF[iAdc][iTimeBin] = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFP) + (fgAddBaseline << fgkAddDigits);
      }
      else {
        fADCR[iAdc][iTimeBin] = (adcArray->GetData(GetRow(), GetCol(iAdc), iTimeBin)  + fgAddBaseline) << fgkAddDigits;
        fADCF[iAdc][iTimeBin] = (adcArray->GetData(GetRow(), GetCol(iAdc), iTimeBin)  + fgAddBaseline) << fgkAddDigits;
      }
    }
  }
}

void AliTRDmcmSim::SetDataPedestal( Int_t iadc )
{
  //
  // Store ADC data into array of raw data
  //

  if( !CheckInitialized() ) return;

  if( iadc < 0 || iadc >= fNADC ) {
    return;
  }

  for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {
    fADCR[iadc][it] = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFP) + (fgAddBaseline << fgkAddDigits);
    fADCF[iadc][it] = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFP) + (fgAddBaseline << fgkAddDigits);
  }
}

Int_t AliTRDmcmSim::GetCol( Int_t iadc )
{
  //
  // Return column id of the pad for the given ADC channel
  //

  if( !CheckInitialized() ) 
    return -1;

  Int_t col = fFeeParam->GetPadColFromADC(fRobPos, fMcmPos, iadc);
  if (col < 0 || col >= fFeeParam->GetNcol()) 
    return -1;
  else 
    return col;
}

Int_t AliTRDmcmSim::ProduceRawStream( UInt_t *buf, Int_t maxSize, UInt_t iEv)
{
  //
  // Produce raw data stream from this MCM and put in buf
  // Returns number of words filled, or negative value 
  // with -1 * number of overflowed words
  //

  UInt_t  x;
  Int_t   nw  = 0;  // Number of written words
  Int_t   of  = 0;  // Number of overflowed words
  Int_t   rawVer   = fFeeParam->GetRAWversion();
  Int_t **adc;
  Int_t   nActiveADC = 0;	// number of activated ADC bits in a word

  if( !CheckInitialized() ) return 0;

  if( fFeeParam->GetRAWstoreRaw() ) {
    adc = fADCR;
  } else {
    adc = fADCF;
  }

  // Produce MCM header
  x = (1<<31) | (fRobPos << 28) | (fMcmPos << 24) | ((iEv % 0x100000) << 4) | 0xC;

  if (nw < maxSize) {
    buf[nw++] = x;
    //printf("\nMCM header: %X ",x);
  }
  else {
    of++;
  }

  // Produce ADC mask : nncc cccm mmmm mmmm mmmm mmmm mmmm 1100
  // 				n : unused , c : ADC count, m : selected ADCs
  if( rawVer >= 3 ) {
    x = 0;
    for( Int_t iAdc = 0 ; iAdc < fNADC ; iAdc++ ) {
      if( fZSM1Dim[iAdc] == 0 ) { //  0 means not suppressed
		x = x | (1 << (iAdc+4) );	// last 4 digit reserved for 1100=0xc
		nActiveADC++;		// number of 1 in mmm....m
      }
    }
	x = x | (1 << 30) | ( ( 0x3FFFFFFC ) & (~(nActiveADC) << 25) ) | 0xC;	// nn = 01, ccccc are inverted, 0xc=1100
	//printf("nActiveADC=%d=%08X, inverted=%X ",nActiveADC,nActiveADC,x );

    if (nw < maxSize) {
      buf[nw++] = x;
      //printf("ADC mask: %X nMask=%d ADC data: ",x,nActiveADC);
    }
    else {
      of++;
    }
  }

  // Produce ADC data. 3 timebins are packed into one 32 bits word
  // In this version, different ADC channel will NOT share the same word

  UInt_t aa=0, a1=0, a2=0, a3=0;

  for (Int_t iAdc = 0; iAdc < 21; iAdc++ ) {
    if( rawVer>= 3 && fZSM1Dim[iAdc] != 0 ) continue; // Zero Suppression, 0 means not suppressed
    aa = !(iAdc & 1) + 2;
    for (Int_t iT = 0; iT < fNTimeBin; iT+=3 ) {
      a1 = ((iT    ) < fNTimeBin ) ? adc[iAdc][iT  ] >> fgkAddDigits : 0;
      a2 = ((iT + 1) < fNTimeBin ) ? adc[iAdc][iT+1] >> fgkAddDigits : 0;
      a3 = ((iT + 2) < fNTimeBin ) ? adc[iAdc][iT+2] >> fgkAddDigits : 0;
      x = (a3 << 22) | (a2 << 12) | (a1 << 2) | aa;
      if (nw < maxSize) {
        buf[nw++] = x;
        //printf("%08X ",x);
      }
      else {
        of++;
      }
    }
  }

  if( of != 0 ) return -of; else return nw;
}

Int_t AliTRDmcmSim::ProduceTrackletStream( UInt_t *buf, Int_t maxSize )
{
  //
  // Produce tracklet data stream from this MCM and put in buf
  // Returns number of words filled, or negative value 
  // with -1 * number of overflowed words
  //

  Int_t   nw  = 0;  // Number of written words
  Int_t   of  = 0;  // Number of overflowed words
    
  if( !CheckInitialized() ) return 0;

  // Produce tracklet data. A maximum of four 32 Bit words will be written per MCM 
  // fMCMT is filled continuously until no more tracklet words available

  for (Int_t iTracklet = 0; iTracklet < fTrackletArray->GetEntriesFast(); iTracklet++) {
    if (nw < maxSize) 
      buf[nw++] = ((AliTRDtrackletMCM*) (*fTrackletArray)[iTracklet])->GetTrackletWord();
    else 
      of++;
  }
  
  if( of != 0 ) return -of; else return nw;
}

void AliTRDmcmSim::Filter()
{
  //
  // Filter the raw ADC values. The active filter stages and their
  // parameters are taken from AliTRDtrapConfig.
  // The raw data is stored separate from the filtered data. Thus, 
  // it is possible to run the filters on a set of raw values 
  // sequentially for parameter tuning.
  //

  if( !CheckInitialized() ) {
    AliError("got called before initialization! Nothing done!");
    return;
  }

  // Apply filters sequentially. Bypass is handled by filters
  // since counters and internal registers may be updated even 
  // if the filter is bypassed.
  // The first filter takes the data from fADCR and 
  // outputs to fADCF. 
  
  // Non-linearity filter not implemented.
  FilterPedestal();
  FilterGain();
  FilterTail();
  // Crosstalk filter not implemented.
}

void AliTRDmcmSim::FilterPedestalInit() 
{
  // Initializes the pedestal filter assuming that the input has 
  // been constant for a long time (compared to the time constant).

//  UShort_t    fpnp = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFPNP); // 0..511 -> 0..127.75, pedestal at the output
  UShort_t    fptc = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFPTC); // 0..3, 0 - fastest, 3 - slowest
  UShort_t    shifts[4] = {11, 14, 17, 21}; //??? where to take shifts from?

  for (Int_t iAdc = 0; iAdc < fNADC; iAdc++)
    fPedAcc[iAdc] = (fSimParam->GetADCbaseline() << 2) * (1<<shifts[fptc]);
}

UShort_t AliTRDmcmSim::FilterPedestalNextSample(Int_t adc, Int_t timebin, UShort_t value)
{
  // Returns the output of the pedestal filter given the input value.
  // The output depends on the internal registers and, thus, the 
  // history of the filter.

  UShort_t    fpnp = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFPNP); // 0..511 -> 0..127.75, pedestal at the output
  UShort_t    fptc = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFPTC); // 0..3, 0 - fastest, 3 - slowest
  UShort_t    fpby = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFPBY); // 0..1 the bypass, active low
  UShort_t    shifts[4] = {11, 14, 17, 21}; //??? where to come from

  UShort_t accumulatorShifted;
  Int_t correction;
  UShort_t inpAdd;
  
  inpAdd = value + fpnp;

  if (fpby == 0) //??? before or after update of accumulator
    return value;

  accumulatorShifted = (fPedAcc[adc] >> shifts[fptc]) & 0x3FF;   // 10 bits
  if (timebin == 0) // the accumulator is disabled in the drift time
  {
    correction = (value & 0x3FF) - accumulatorShifted;
    fPedAcc[adc] = (fPedAcc[adc] + correction) & 0x7FFFFFFF;             // 31 bits
  }
  
  if (inpAdd <= accumulatorShifted)
    return 0;
  else
  {
    inpAdd = inpAdd - accumulatorShifted;
    if (inpAdd > 0xFFF) 
      return 0xFFF;
    else 
      return inpAdd;
  }
}

void AliTRDmcmSim::FilterPedestal()
{
  //
  // Apply pedestal filter
  //
  // As the first filter in the chain it reads data from fADCR 
  // and outputs to fADCF. 
  // It has only an effect if previous samples have been fed to 
  // find the pedestal. Currently, the simulation assumes that 
  // the input has been stable for a sufficiently long time.

  for (Int_t iTimeBin = 0; iTimeBin < fNTimeBin; iTimeBin++) {
    for (Int_t iAdc = 0; iAdc < fNADC; iAdc++) {
      fADCF[iAdc][iTimeBin] = FilterPedestalNextSample(iAdc, iTimeBin, fADCR[iAdc][iTimeBin]);
    }
  }
}

void AliTRDmcmSim::FilterGainInit()
{
  // Initializes the gain filter. In this case, only threshold 
  // counters are reset.

  for (Int_t iAdc = 0; iAdc < fNADC; iAdc++) {
    // these are counters which in hardware continue 
    // until maximum or reset
    fGainCounterA[iAdc] = 0;
    fGainCounterB[iAdc] = 0;
  }
}

UShort_t AliTRDmcmSim::FilterGainNextSample(Int_t adc, UShort_t value)
{
  // Apply the gain filter to the given value.
  // BEGIN_LATEX O_{i}(t) = #gamma_{i} * I_{i}(t) + a_{i} END_LATEX
  // The output depends on the internal registers and, thus, the 
  // history of the filter.

  UShort_t    fgby = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFGBY); // bypass, active low
  UShort_t    fgf  = fTrapConfig->GetTrapReg(AliTRDtrapConfig::TrapReg_t(AliTRDtrapConfig::kFGF0 + adc)); // 0x700 + (0 & 0x1ff);
  UShort_t    fga  = fTrapConfig->GetTrapReg(AliTRDtrapConfig::TrapReg_t(AliTRDtrapConfig::kFGA0 + adc)); // 40;
  UShort_t    fgta = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFGTA); // 20;
  UShort_t    fgtb = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFGTB); // 2060;

  UInt_t tmp;

  value &= 0xFFF;
  tmp = (value * fgf) >> 11;
  if (tmp > 0xFFF) tmp = 0xFFF;

  if (fgby == 1)
    value = AddUintClipping(tmp, fga, 12);

  // Update threshold counters 
  // not really useful as they are cleared with every new event
  if ((fGainCounterA[adc] == 0x3FFFFFF) || (fGainCounterB[adc] == 0x3FFFFFF))
  {
    if (value >= fgtb) 
      fGainCounterB[adc]++;
    else if (value >= fgta) 
      fGainCounterA[adc]++;
  }

  return value;
}

void AliTRDmcmSim::FilterGain()
{
  // Read data from fADCF and apply gain filter.

  for (Int_t iAdc = 0; iAdc < fNADC; iAdc++) {
    for (Int_t iTimeBin = 0; iTimeBin < fNTimeBin; iTimeBin++) {
        fADCF[iAdc][iTimeBin] = FilterGainNextSample(iAdc, fADCF[iAdc][iTimeBin]);
    }
  }
}

void AliTRDmcmSim::FilterTailInit(Int_t baseline)
{
  // Initializes the tail filter assuming that the input has 
  // been at the baseline value (configured by FTFP) for a 
  // sufficiently long time.

  // exponents and weight calculated from configuration
  UShort_t    alphaLong = 0x3ff & fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFTAL); // the weight of the long component
  UShort_t    lambdaLong = (1 << 10) | (1 << 9) | (fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFTLL) & 0x1FF); // the multiplier
  UShort_t    lambdaShort = (0 << 10) | (1 << 9) | (fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFTLS) & 0x1FF); // the multiplier

  Float_t lambdaL = lambdaLong  * 1.0 / (1 << 11);
  Float_t lambdaS = lambdaShort * 1.0 / (1 << 11);
  Float_t alphaL  = alphaLong   * 1.0 / (1 << 11);
  Float_t qup, qdn;
  qup = (1 - lambdaL) * (1 - lambdaS);
  qdn = 1 - lambdaS * alphaL - lambdaL * (1 - alphaL);
  Float_t kdc = qup/qdn;

  Float_t kt, ql, qs;
  UShort_t aout;
  
  kt = kdc * baseline;
  aout = baseline - (UShort_t) kt;
  ql = lambdaL * (1 - lambdaS) *      alphaL;
  qs = lambdaS * (1 - lambdaL) * (1 - alphaL);

  for (Int_t iAdc = 0; iAdc < fNADC; iAdc++) {
    fTailAmplLong[iAdc]  = (UShort_t) (aout * ql / (ql + qs));
    fTailAmplShort[iAdc] = (UShort_t) (aout * qs / (ql + qs));
  }
}

UShort_t AliTRDmcmSim::FilterTailNextSample(Int_t adc, UShort_t value)
{
  // Returns the output of the tail filter for the given input value. 
  // The output depends on the internal registers and, thus, the 
  // history of the filter.

  // exponents and weight calculated from configuration
  UShort_t    alphaLong = 0x3ff & fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFTAL); // the weight of the long component
  UShort_t    lambdaLong = (1 << 10) | (1 << 9) | (fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFTLL) & 0x1FF); // the multiplier
  UShort_t    lambdaShort = (0 << 10) | (1 << 9) | (fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFTLS) & 0x1FF); // the multiplier

  Float_t lambdaL = lambdaLong  * 1.0 / (1 << 11);
  Float_t lambdaS = lambdaShort * 1.0 / (1 << 11);
  Float_t alphaL  = alphaLong   * 1.0 / (1 << 11);
  Float_t qup, qdn;
  qup = (1 - lambdaL) * (1 - lambdaS);
  qdn = 1 - lambdaS * alphaL - lambdaL * (1 - alphaL);
//  Float_t kdc = qup/qdn;

  UInt_t aDiff;
  UInt_t alInpv;
  UShort_t aQ;
  UInt_t tmp;
  
  UShort_t inpVolt = value & 0xFFF;    // 12 bits
      
  if (fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFTBY) == 0) // bypass mode, active low
    return value;
  else
  {   
    // add the present generator outputs
    aQ = AddUintClipping(fTailAmplLong[adc], fTailAmplShort[adc], 12);

    // calculate the difference between the input the generated signal
    if (inpVolt > aQ) 
      aDiff = inpVolt - aQ;
    else                
      aDiff = 0;

    // the inputs to the two generators, weighted
    alInpv = (aDiff * alphaLong) >> 11;

    // the new values of the registers, used next time
    // long component
    tmp = AddUintClipping(fTailAmplLong[adc], alInpv, 12);
    tmp =  (tmp * lambdaLong) >> 11;
    fTailAmplLong[adc] = tmp & 0xFFF;
    // short component
    tmp = AddUintClipping(fTailAmplShort[adc], aDiff - alInpv, 12);
    tmp =  (tmp * lambdaShort) >> 11;
    fTailAmplShort[adc] = tmp & 0xFFF;

    // the output of the filter
    return aDiff;
  }
}

void AliTRDmcmSim::FilterTail()
{
  // Apply tail filter

  for (Int_t iTimeBin = 0; iTimeBin < fNTimeBin; iTimeBin++) {
    for (Int_t iAdc = 0; iAdc < fNADC; iAdc++) {
      fADCF[iAdc][iTimeBin] = FilterTailNextSample(iAdc, fADCF[iAdc][iTimeBin]);
    }
  }
}

void AliTRDmcmSim::ZSMapping()
{
  //
  // Zero Suppression Mapping implemented in TRAP chip
  //
  // See detail TRAP manual "Data Indication" section:
  // http://www.kip.uni-heidelberg.de/ti/TRD/doc/trap/TRAP-UserManual.pdf
  //

  //??? values should come from TRAPconfig
  Int_t eBIS = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kEBIS); // TRAP default = 0x4  (Tis=4)
  Int_t eBIT = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kEBIT); // TRAP default = 0x28 (Tit=40)
  Int_t eBIL = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kEBIL); // TRAP default = 0xf0
                                                                 // (lookup table accept (I2,I1,I0)=(111)
                                                                 // or (110) or (101) or (100))
  Int_t eBIN = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kEBIN); // TRAP default = 1 (no neighbor sensitivity)
  Int_t ep   = 0; // fTrapConfig->GetTrapReg(AliTRDtrapConfig::kFPNP); //??? really subtracted here

  Int_t **adc = fADCF;

  if( !CheckInitialized() ) {
    AliError("got called uninitialized! Nothing done!");    
    return;
  }

  for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {
    for( Int_t iadc = 1 ; iadc < fNADC-1; iadc++ ) {

      // Get ADC data currently in filter buffer
      Int_t ap = adc[iadc-1][it] - ep; // previous
      Int_t ac = adc[iadc  ][it] - ep; // current
      Int_t an = adc[iadc+1][it] - ep; // next

      // evaluate three conditions
      Int_t i0 = ( ac >=  ap && ac >=  an ) ? 0 : 1; // peak center detection
      Int_t i1 = ( ap + ac + an > eBIT )    ? 0 : 1; // cluster
      Int_t i2 = ( ac > eBIS )              ? 0 : 1; // absolute large peak

      Int_t i = i2 * 4 + i1 * 2 + i0;    // Bit position in lookup table
      Int_t d = (eBIL >> i) & 1;         // Looking up  (here d=0 means true
                                         // and d=1 means false according to TRAP manual)

      fZSM[iadc][it] &= d;
      if( eBIN == 0 ) {  // turn on neighboring ADCs
	fZSM[iadc-1][it] &= d;
	fZSM[iadc+1][it] &= d;
      }
    }
  }

  // do 1 dim projection
  for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
    for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {
      fZSM1Dim[iadc] &= fZSM[iadc][it];
    }
  }
}

void AliTRDmcmSim::DumpData( const char * const f, const char * const target )
{
  //
  // Dump data stored (for debugging).
  // target should contain one or multiple of the following characters
  //   R   for raw data
  //   F   for filtered data
  //   Z   for zero suppression map
  //   S   Raw dat astream
  // other characters are simply ignored
  //

  UInt_t tempbuf[1024];

  if( !CheckInitialized() ) return;

  std::ofstream of( f, std::ios::out | std::ios::app );
  of << Form("AliTRDmcmSim::DumpData det=%03d sm=%02d stack=%d layer=%d rob=%d mcm=%02d\n",
	     fDetector, fGeo->GetSector(fDetector), fGeo->GetStack(fDetector), 
             fGeo->GetSector(fDetector), fRobPos, fMcmPos );

  for( Int_t t=0 ; target[t] != 0 ; t++ ) {
    switch( target[t] ) {
    case 'R' :
    case 'r' :
      of << Form("fADCR (raw ADC data)\n");
      for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
	of << Form("  ADC %02d: ", iadc);
	for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {
	  of << Form("% 4d",  fADCR[iadc][it]);
	}
	of << Form("\n");
      }
      break;
    case 'F' :
    case 'f' :
      of << Form("fADCF (filtered ADC data)\n");
      for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
	of << Form("  ADC %02d: ", iadc);
	for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {
	  of << Form("% 4d",  fADCF[iadc][it]);
	}
	of << Form("\n");
      }
      break;
    case 'Z' :
    case 'z' :
      of << Form("fZSM and fZSM1Dim (Zero Suppression Map)\n");
      for( Int_t iadc = 0 ; iadc < fNADC; iadc++ ) {
	of << Form("  ADC %02d: ", iadc);
	if( fZSM1Dim[iadc] == 0 ) { of << " R   " ; } else { of << " .   "; } // R:read .:suppressed
	for( Int_t it = 0 ; it < fNTimeBin ; it++ ) {
	  if( fZSM[iadc][it] == 0 ) { of << " R"; } else { of << " ."; } // R:read .:suppressed
	}
	of << Form("\n");
      }
      break;
    case 'S' :
    case 's' :
      Int_t s = ProduceRawStream( tempbuf, 1024 ); 
      of << Form("Stream for Raw Simulation size=%d rawver=%d\n", s, fFeeParam->GetRAWversion());
      of << Form("  address  data\n");
      for( Int_t i = 0 ; i < s ; i++ ) {
	of << Form("  %04x     %08x\n", i, tempbuf[i]);
      }
    }
  }
}

void AliTRDmcmSim::AddHitToFitreg(Int_t adc, UShort_t timebin, UShort_t qtot, Short_t ypos, Int_t label) 
{
  // Add the given hit to the fit register which is lateron used for 
  // the tracklet calculation. 
  // In addition to the fit sums in the fit register MC information 
  // is stored.

  if ((timebin >= fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPQS0)) && 
      (timebin <  fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPQE0)))
    fFitReg[adc].fQ0 += qtot;
  
  if ((timebin >= fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPQS1)) && 
      (timebin <  fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPQE1)))
    fFitReg[adc].fQ1 += qtot;
  
  if ((timebin >= fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFS) ) && 
      (timebin <  fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFE)))
  {
    fFitReg[adc].fSumX  += timebin;
    fFitReg[adc].fSumX2 += timebin*timebin;
    fFitReg[adc].fNhits++;
    fFitReg[adc].fSumY  += ypos;
    fFitReg[adc].fSumY2 += ypos*ypos;
    fFitReg[adc].fSumXY += timebin*ypos;
  }

  // register hits (MC info)
  fHits[fNHits].fChannel = adc;
  fHits[fNHits].fQtot = qtot;
  fHits[fNHits].fYpos = ypos;
  fHits[fNHits].fTimebin = timebin;
  fHits[fNHits].fLabel = label;
  fNHits++;
}

void AliTRDmcmSim::CalcFitreg() 
{
  // Preprocessing.
  // Detect the hits and fill the fit registers.
  // Requires 12-bit data from fADCF which means Filter() 
  // has to be called before even if all filters are bypassed.

  //???
  // TRAP parameters:
  const UShort_t lutPos[128] = {   // move later to some other file
    0,  1,  1,  2,  2,  3,  3,  4,  4,  5,  5,  6,  6,  7,  7,  8,  8,  9,  9, 10, 10, 11, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15,
    16, 16, 16, 17, 17, 18, 18, 19, 19, 19, 20, 20, 20, 21, 21, 22, 22, 22, 23, 23, 23, 24, 24, 24, 24, 25, 25, 25, 26, 26, 26, 26,
    27, 27, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 27, 27, 27, 27, 26,
    26, 26, 26, 25, 25, 25, 24, 24, 23, 23, 22, 22, 21, 21, 20, 20, 19, 18, 18, 17, 17, 16, 15, 14, 13, 12, 11, 10,  9,  8,  7,  7};
  
  //??? to be clarified:
  UInt_t adcMask = 0xffffffff;
  
  UShort_t timebin, adcch, adcLeft, adcCentral, adcRight, hitQual, timebin1, timebin2, qtotTemp;
  Short_t ypos, fromLeft, fromRight, found;
  UShort_t qTotal[19]; // the last is dummy
  UShort_t marked[6], qMarked[6], worse1, worse2;
  
  timebin1 = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFS); 
  if (fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPQS0) 
      < timebin1)
    timebin1 = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPQS0);
  timebin2 = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFE); 
  if (fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPQE1) 
      > timebin2)
    timebin2 = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPQE1);

  // reset the fit registers
  fNHits = 0; 
  for (adcch = 0; adcch < fNADC-2; adcch++) // due to border channels
  {
    fFitReg[adcch].fNhits = 0;
    fFitReg[adcch].fQ0    = 0;
    fFitReg[adcch].fQ1    = 0;
    fFitReg[adcch].fSumX  = 0;
    fFitReg[adcch].fSumY  = 0;
    fFitReg[adcch].fSumX2 = 0;
    fFitReg[adcch].fSumY2 = 0;
    fFitReg[adcch].fSumXY = 0;
  }
  
  for (timebin = timebin1; timebin < timebin2; timebin++)
  {
    // first find the hit candidates and store the total cluster charge in qTotal array
    // in case of not hit store 0 there.
    for (adcch = 0; adcch < fNADC-2; adcch++) {
      if ( ( (adcMask >> adcch) & 7) == 7) //??? all 3 channels are present in case of ZS
      {
        adcLeft  = fADCF[adcch  ][timebin];
        adcCentral  = fADCF[adcch+1][timebin];
        adcRight = fADCF[adcch+2][timebin];
        if (fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPVBY) == 1) 
          hitQual = ( (adcLeft * adcRight) < 
                       (fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPVT) * adcCentral) );
        else            
          hitQual = 1;
        // The accumulated charge is with the pedestal!!!
        qtotTemp = adcLeft + adcCentral + adcRight;
        if ( (hitQual) &&
             (qtotTemp >= fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPHT)) &&
             (adcLeft <= adcCentral) &&
             (adcCentral > adcRight) )
          qTotal[adcch] = qtotTemp;
        else
          qTotal[adcch] = 0;
      }
      else
        qTotal[adcch] = 0; //jkl
      AliDebug(10,Form("ch %2d   qTotal %5d",adcch, qTotal[adcch]));
    }

    fromLeft = -1;
    adcch = 0;
    found = 0;
    marked[4] = 19; // invalid channel
    marked[5] = 19; // invalid channel
    qTotal[19] = 0;
    while ((adcch < 16) && (found < 3))
    {
      if (qTotal[adcch] > 0)
      {
        fromLeft = adcch;
        marked[2*found+1]=adcch;
        found++;
      }
      adcch++;
    }
    
    fromRight = -1;
    adcch = 18;
    found = 0;
    while ((adcch > 2) && (found < 3))
    {
      if (qTotal[adcch] > 0)
      {
        marked[2*found]=adcch;
        found++;
        fromRight = adcch;
      }
      adcch--;
    }

    AliDebug(10,Form("Fromleft=%d, Fromright=%d",fromLeft, fromRight));
    // here mask the hit candidates in the middle, if any
    if ((fromLeft >= 0) && (fromRight >= 0) && (fromLeft < fromRight))
      for (adcch = fromLeft+1; adcch < fromRight; adcch++)
        qTotal[adcch] = 0;
    
    found = 0;
    for (adcch = 0; adcch < 19; adcch++)
      if (qTotal[adcch] > 0) found++;
    // NOT READY

    if (found > 4) // sorting like in the TRAP in case of 5 or 6 candidates!
    {
      if (marked[4] == marked[5]) marked[5] = 19;
      for (found=0; found<6; found++)
      {
        qMarked[found] = qTotal[marked[found]] >> 4;
        AliDebug(10,Form("ch_%d qTotal %d qTotals %d",marked[found],qTotal[marked[found]],qMarked[found]));
      }
      
      Sort6To2Worst(marked[0], marked[3], marked[4], marked[1], marked[2], marked[5],
                    qMarked[0],
                    qMarked[3],
                    qMarked[4],
                    qMarked[1],
                    qMarked[2],
                    qMarked[5],
                    &worse1, &worse2);
      // Now mask the two channels with the smallest charge
      if (worse1 < 19)
      {
        qTotal[worse1] = 0;
        AliDebug(10,Form("Kill ch %d\n",worse1));
      }
      if (worse2 < 19)
      {
        qTotal[worse2] = 0;
        AliDebug(10,Form("Kill ch %d\n",worse2));
      }
    }
    
    for (adcch = 0; adcch < 19; adcch++) {
      if (qTotal[adcch] > 0) // the channel is marked for processing
      {
        adcLeft  = fADCF[adcch  ][timebin];
        adcCentral  = fADCF[adcch+1][timebin];
        adcRight = fADCF[adcch+2][timebin];
        // hit detected, in TRAP we have 4 units and a hit-selection, here we proceed all channels!
        // subtract the pedestal TPFP, clipping instead of wrapping
        
        Int_t regTPFP = fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPFP);
        AliDebug(10, Form("Hit found, time=%d, adcch=%d/%d/%d, adc values=%d/%d/%d, regTPFP=%d, TPHT=%d\n",
               timebin, adcch, adcch+1, adcch+2, adcLeft, adcCentral, adcRight, regTPFP, 
               fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPHT)));

        if (adcLeft  < regTPFP) adcLeft  = 0; else adcLeft  -= regTPFP;
        if (adcCentral  < regTPFP) adcCentral  = 0; else adcCentral  -= regTPFP;
        if (adcRight < regTPFP) adcRight = 0; else adcRight -= regTPFP;

        // Calculate the center of gravity
        // checking for adcCentral != 0 (in case of "bad" configuration)
        if (adcCentral == 0)
          continue;
        ypos = 128*(adcLeft - adcRight) / adcCentral;
        if (ypos < 0) ypos = -ypos;
        // make the correction using the LUT
        ypos = ypos + lutPos[ypos & 0x7F];
        if (adcLeft > adcRight) ypos = -ypos;

        // label calculation
        Int_t mcLabel = -1;
        if (fDigitsManager) {
          Int_t label[9] = { 0 }; // up to 9 different labels possible
          Int_t count[9] = { 0 };
          Int_t maxIdx = -1;
          Int_t maxCount = 0;
          Int_t nLabels = 0;
          Int_t padcol[3]; 
          padcol[0] = fFeeParam->GetPadColFromADC(fRobPos, fMcmPos, adcch);
          padcol[1] = fFeeParam->GetPadColFromADC(fRobPos, fMcmPos, adcch+1);
          padcol[2] = fFeeParam->GetPadColFromADC(fRobPos, fMcmPos, adcch+2);
          Int_t padrow = fFeeParam->GetPadRowFromMCM(fRobPos, fMcmPos);
          for (Int_t iDict = 0; iDict < 3; iDict++) {
            if (!fDigitsManager->UsesDictionaries() || fDigitsManager->GetDictionary(fDetector, iDict) == 0) {
              AliError("Cannot get dictionary");
              continue;
            }
            AliTRDarrayDictionary *dict = (AliTRDarrayDictionary*) fDigitsManager->GetDictionary(fDetector, iDict);
            if (dict->GetDim() == 0) {
              AliError(Form("Dictionary %i of det. %i has dim. 0", fDetector, iDict));
              continue;
            }
            dict->Expand();
            for (Int_t iPad = 0; iPad < 3; iPad++) {
              if (padcol[iPad] < 0) 
                continue;
              Int_t currLabel = dict->GetData(padrow, padcol[iPad], timebin); //fDigitsManager->GetTrack(iDict, padrow, padcol, timebin, fDetector);
	      AliDebug(10, Form("Read label: %4i for det: %3i, row: %i, col: %i, tb: %i\n", currLabel, fDetector, padrow, padcol[iPad], timebin));
              for (Int_t iLabel = 0; iLabel < nLabels; iLabel++) {
                if (currLabel == label[iLabel]) {
                  count[iLabel]++;
                  if (count[iLabel] > maxCount) {
                    maxCount = count[iLabel];
                    maxIdx = iLabel;
                  }
                  currLabel = 0;
                  break;
                }
              } 
              if (currLabel > 0) {
                label[nLabels++] = currLabel;
              }
            }
          }
          if (maxIdx >= 0)
            mcLabel = label[maxIdx];
        }

        // add the hit to the fitregister
        AddHitToFitreg(adcch, timebin, qTotal[adcch], ypos, mcLabel);
      }
    }
  }
}

void AliTRDmcmSim::TrackletSelection() 
{
  // Select up to 4 tracklet candidates from the fit registers  
  // and assign them to the CPUs.

  UShort_t adcIdx, i, j, ntracks, tmp;
  UShort_t trackletCand[18][2]; // store the adcch[0] and number of hits[1] for all tracklet candidates

  ntracks = 0;
  for (adcIdx = 0; adcIdx < 18; adcIdx++) // ADCs
    if ( (fFitReg[adcIdx].fNhits 
          >= fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPCL)) &&
         (fFitReg[adcIdx].fNhits+fFitReg[adcIdx+1].fNhits
          >= fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPCT)))
    {
      trackletCand[ntracks][0] = adcIdx;
      trackletCand[ntracks][1] = fFitReg[adcIdx].fNhits+fFitReg[adcIdx+1].fNhits;
      AliDebug(10,Form("%d  %2d %4d\n", ntracks, trackletCand[ntracks][0], trackletCand[ntracks][1]));
      ntracks++;
    };

  for (i=0; i<ntracks;i++) 
    AliDebug(10,Form("%d %d %d\n",i,trackletCand[i][0], trackletCand[i][1]));

  if (ntracks > 4)
  {
    // primitive sorting according to the number of hits
    for (j = 0; j < (ntracks-1); j++)
    {
      for (i = j+1; i < ntracks; i++)
      {
        if ( (trackletCand[j][1]  < trackletCand[i][1]) ||
             ( (trackletCand[j][1] == trackletCand[i][1]) && (trackletCand[j][0] < trackletCand[i][0]) ) )
        {
          // swap j & i
          tmp = trackletCand[j][1];
          trackletCand[j][1] = trackletCand[i][1];
          trackletCand[i][1] = tmp;
          tmp = trackletCand[j][0];
          trackletCand[j][0] = trackletCand[i][0];
          trackletCand[i][0] = tmp;
        }
      }
    }
    ntracks = 4; // cut the rest, 4 is the max
  }
  // else is not necessary to sort
  
  // now sort, so that the first tracklet going to CPU0 corresponds to the highest adc channel - as in the TRAP
  for (j = 0; j < (ntracks-1); j++)
  {
    for (i = j+1; i < ntracks; i++)
    {
      if (trackletCand[j][0] < trackletCand[i][0])
      {
        // swap j & i
        tmp = trackletCand[j][1];
        trackletCand[j][1] = trackletCand[i][1];
        trackletCand[i][1] = tmp;
        tmp = trackletCand[j][0];
        trackletCand[j][0] = trackletCand[i][0];
        trackletCand[i][0] = tmp;
      }
    }
  }
  for (i = 0; i < ntracks; i++)  // CPUs with tracklets.
    fFitPtr[i] = trackletCand[i][0]; // pointer to the left channel with tracklet for CPU[i]
  for (i = ntracks; i < 4; i++)  // CPUs without tracklets
    fFitPtr[i] = 31;            // pointer to the left channel with tracklet for CPU[i] = 31 (invalid)
  AliDebug(10,Form("found %i tracklet candidates\n", ntracks));
  for (i = 0; i < 4; i++)
    AliDebug(10,Form("fitPtr[%i]: %i\n", i, fFitPtr[i]));
}

void AliTRDmcmSim::FitTracklet()
{
  // Perform the actual tracklet fit based on the fit sums 
  // which have been filled in the fit registers. 

  // parameters in fitred.asm (fit program)
  Int_t decPlaces = 5;
  Int_t rndAdd = 0;
  if (decPlaces >  1) 
    rndAdd = (1 << (decPlaces-1)) + 1;
  else if (decPlaces == 1)
    rndAdd = 1;
  Int_t ndriftDp = 5;  // decimal places for drift time
  Long64_t shift = ((Long64_t) 1 << 32);


  // calculated in fitred.asm
  Int_t padrow = ((fRobPos >> 1) << 2) | (fMcmPos >> 2);
  Int_t yoffs = (((((fRobPos & 0x1) << 2) + (fMcmPos & 0x3)) * 18) << 8) - 
    ((18*4*2 - 18*2 - 1) << 7);
  yoffs = yoffs << decPlaces; // holds position of ADC channel 1
  Int_t layer = fDetector % 6;
  UInt_t scaleY = (UInt_t) ((0.635 + 0.03 * layer)/(256.0 * 160.0e-4) * shift);
  UInt_t scaleD = (UInt_t) ((0.635 + 0.03 * layer)/(256.0 * 140.0e-4) * shift);
  // previously taken from geometry:
  // UInt_t scaleYold = (UInt_t) (shift * (pp->GetWidthIPad() / (256 * 160e-4)));
  // UInt_t scaleDold = (UInt_t) (shift * (pp->GetWidthIPad() / (256 * 140e-4)));


  // should come from trapConfig (DMEM) 
  AliTRDpadPlane *pp = fGeo->GetPadPlane(fDetector);
  Float_t scaleSlope = (256 / pp->GetWidthIPad()) * (1 << decPlaces); // only used for calculation of corrections and cut
  Int_t ndrift   = 20 << ndriftDp; //??? value in simulation?
  Int_t deflCorr = (Int_t) (TMath::Tan(fCommonParam->GetOmegaTau(fCal->GetVdriftAverage(fDetector))) * fGeo->CdrHght() * scaleSlope); // -370;
  Int_t tiltCorr = (Int_t) (pp->GetRowPos(padrow) / fGeo->GetTime0(fDetector % 6) * fGeo->CdrHght() * scaleSlope * 
                            TMath::Tan(pp->GetTiltingAngle() / 180. * TMath::Pi()));
//  printf("vdrift av.: %f\n", fCal->GetVdriftAverage(fDetector));
//  printf("chamber height: %f\n", fGeo->CdrHght());
//  printf("omega tau: %f\n", fCommonParam->GetOmegaTau(fCal->GetVdriftAverage(fDetector)));
//  printf("deflection correction: %i\n", deflCorr);
  Float_t ptcut = 2.3;
  AliMagF* fld = (AliMagF *) TGeoGlobalMagField::Instance()->GetField();
  Double_t bz = 0;
  if (fld) {
    bz       = 0.1 * fld->SolenoidField();   // kGauss -> Tesla
  }
//  printf("Bz: %f\n", bz);
  Float_t x0 = fGeo->GetTime0(fDetector % 6);
  Float_t y0 = pp->GetColPos(fFeeParam->GetPadColFromADC(fRobPos, fMcmPos, 10));
  Float_t alphaMax = TMath::ASin( (TMath::Sqrt(TMath::Power(x0/100., 2) + TMath::Power(y0/100., 2)) * 
                                   0.3 * TMath::Abs(bz) ) / (2 * ptcut));
//  printf("alpha max: %f\n", alphaMax * 180/TMath::Pi());
  Int_t minslope = -1 * (Int_t) (fGeo->CdrHght() * TMath::Tan(TMath::ATan(y0/x0) + alphaMax) / 140.e-4);
  Int_t maxslope = -1 * (Int_t) (fGeo->CdrHght() * TMath::Tan(TMath::ATan(y0/x0) - alphaMax) / 140.e-4);


  // local variables for calculation
  Long64_t mult, temp, denom; //???
  UInt_t q0, q1, qTotal;          // charges in the two windows and total charge
  UShort_t nHits;                 // number of hits
  Int_t slope, offset;            // slope and offset of the tracklet
  Int_t sumX, sumY, sumXY, sumX2; // fit sums from fit registers
  //int32_t SumY2;                // not used in the current TRAP program
  FitReg_t *fit0, *fit1;          // pointers to relevant fit registers
  
//  const uint32_t OneDivN[32] = {  // 2**31/N : exactly like in the TRAP, the simple division here gives the same result!
//      0x00000000, 0x80000000, 0x40000000, 0x2AAAAAA0, 0x20000000, 0x19999990, 0x15555550, 0x12492490,
//      0x10000000, 0x0E38E380, 0x0CCCCCC0, 0x0BA2E8B0, 0x0AAAAAA0, 0x09D89D80, 0x09249240, 0x08888880,
//      0x08000000, 0x07878780, 0x071C71C0, 0x06BCA1A0, 0x06666660, 0x06186180, 0x05D17450, 0x0590B210,
//      0x05555550, 0x051EB850, 0x04EC4EC0, 0x04BDA120, 0x04924920, 0x0469EE50, 0x04444440, 0x04210840};

  for (Int_t cpu = 0; cpu < 4; cpu++) {
    if (fFitPtr[cpu] == 31)
    {
      fMCMT[cpu] = 0x10001000; //??? AliTRDfeeParam::GetTrackletEndmarker(); 
    }
    else
    {
      fit0 = &fFitReg[fFitPtr[cpu]  ];
      fit1 = &fFitReg[fFitPtr[cpu]+1]; // next channel

      mult = 1;
      mult = mult << (32 + decPlaces);
      mult = -mult;

      // Merging
      nHits   = fit0->fNhits + fit1->fNhits; // number of hits
      sumX    = fit0->fSumX  + fit1->fSumX;
      sumX2   = fit0->fSumX2 + fit1->fSumX2;
      denom   = nHits*sumX2 - sumX*sumX;

      mult    = mult / denom; // exactly like in the TRAP program
      q0      = fit0->fQ0    + fit1->fQ0;
      q1      = fit0->fQ1    + fit1->fQ1;
      sumY    = fit0->fSumY  + fit1->fSumY  + 256*fit1->fNhits;
      sumXY   = fit0->fSumXY + fit1->fSumXY + 256*fit1->fSumX;

      slope   = nHits*sumXY - sumX * sumY;
      AliDebug(5, Form("slope from fitreg: %i", slope));
      offset  = sumX2*sumY  - sumX * sumXY;
      temp    = mult * slope;
      slope   = temp >> 32; // take the upper 32 bits
      slope   = -slope;
      temp    = mult * offset;
      offset  = temp >> 32; // take the upper 32 bits

      offset = offset + yoffs;
      AliDebug(5, Form("slope: %i, slope * ndrift: %i, deflCorr: %i, tiltCorr: %i", slope, slope * ndrift, deflCorr, tiltCorr));
      slope  = ((slope * ndrift) >> ndriftDp) + deflCorr + tiltCorr;
      offset = offset - (fFitPtr[cpu] << (8 + decPlaces));
      
      AliDebug(5, Form("Det: %3i, ROB: %i, MCM: %2i: deflection: %i, min: %i, max: %i", fDetector, fRobPos, fMcmPos, slope, minslope, maxslope));
      temp    = slope;
      temp    = temp * scaleD;
      slope   = (temp >> 32);
      AliDebug(5, Form("slope after scaling: %i", slope));

      temp    = offset;
      temp    = temp * scaleY;
      offset  = (temp >> 32);
        
      // rounding, like in the TRAP
      slope   = (slope  + rndAdd) >> decPlaces;
      AliDebug(5, Form("slope after shifting: %i", slope));
      offset  = (offset + rndAdd) >> decPlaces;

      Bool_t rejected = kFALSE;
      if ((slope < minslope) || (slope > maxslope))
        rejected = kTRUE;

      if (rejected && GetApplyCut())
      {
        fMCMT[cpu] = 0x10001000; //??? AliTRDfeeParam::GetTrackletEndmarker();
      }
      else
      {
        if (slope > 63 || slope < -64) { // wrapping in TRAP!
          AliError(Form("Overflow in slope: %i, tracklet discarded!", slope));
          fMCMT[cpu] = 0x10001000;
          continue;
        }

        slope   = slope  &   0x7F; // 7 bit
        
        if (offset > 0xfff || offset < -0xfff) 
          AliWarning("Overflow in offset");
        offset  = offset & 0x1FFF; // 13 bit

	Float_t length = TMath::Sqrt(1 + (pp->GetRowPos(padrow) * pp->GetRowPos(padrow) +
					  (fFeeParam->GetPadColFromADC(fRobPos, fMcmPos, 10) * pp->GetWidthIPad() *
					   fFeeParam->GetPadColFromADC(fRobPos, fMcmPos, 10) * pp->GetWidthIPad())) /
				     (fGeo->GetTime0(fDetector % 6)*fGeo->GetTime0(fDetector % 6)));

	//        qTotal  = (q1 / nHits) >> 1;
	qTotal = GetPID(q0/length/fgChargeNorm, q1/length/fgChargeNorm);
        if (qTotal > 0xff)
          AliWarning("Overflow in charge");
        qTotal  = qTotal & 0xFF; // 8 bit, exactly like in the TRAP program
        
        // assemble and store the tracklet word
        fMCMT[cpu] = (qTotal << 24) | (padrow << 20) | (slope << 13) | offset;

        // calculate MC label
        Int_t mcLabel = -1;
	Int_t nHits0 = 0;
	Int_t nHits1 = 0;
        if (fDigitsManager) {
          Int_t label[30] = {0}; // up to 30 different labels possible
          Int_t count[30] = {0};
          Int_t maxIdx = -1;
          Int_t maxCount = 0;
          Int_t nLabels = 0;
          for (Int_t iHit = 0; iHit < fNHits; iHit++) {
            if ((fHits[iHit].fChannel - fFitPtr[cpu] < 0) ||
                (fHits[iHit].fChannel - fFitPtr[cpu] > 1))
              continue;

	    // counting contributing hits
	    if (fHits[iHit].fTimebin >= fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPQS0) &&
		fHits[iHit].fTimebin <  fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPQE0))
	      nHits0++;
	    if (fHits[iHit].fTimebin >= fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPQS1) &&
		fHits[iHit].fTimebin <  fTrapConfig->GetTrapReg(AliTRDtrapConfig::kTPQE1))
	      nHits1++;

            Int_t currLabel = fHits[iHit].fLabel;
            for (Int_t iLabel = 0; iLabel < nLabels; iLabel++) {
              if (currLabel == label[iLabel]) {
                count[iLabel]++;
                if (count[iLabel] > maxCount) {
                  maxCount = count[iLabel];
                  maxIdx = iLabel;
                }
                currLabel = 0;
                break;
              }
            }
            if (currLabel > 0) {
              label[nLabels++] = currLabel;
            }
          }
          if (maxIdx >= 0)
            mcLabel = label[maxIdx];
        }
        new ((*fTrackletArray)[fTrackletArray->GetEntriesFast()]) AliTRDtrackletMCM((UInt_t) fMCMT[cpu], fDetector*2 + fRobPos%2, fRobPos, fMcmPos);
        ((AliTRDtrackletMCM*) (*fTrackletArray)[fTrackletArray->GetEntriesFast()-1])->SetLabel(mcLabel);

       
        ((AliTRDtrackletMCM*) (*fTrackletArray)[fTrackletArray->GetEntriesFast()-1])->SetNHits(fit0->fNhits + fit1->fNhits);
	((AliTRDtrackletMCM*) (*fTrackletArray)[fTrackletArray->GetEntriesFast()-1])->SetNHits0(nHits0);
        ((AliTRDtrackletMCM*) (*fTrackletArray)[fTrackletArray->GetEntriesFast()-1])->SetNHits1(nHits1);
        ((AliTRDtrackletMCM*) (*fTrackletArray)[fTrackletArray->GetEntriesFast()-1])->SetQ0(q0);
        ((AliTRDtrackletMCM*) (*fTrackletArray)[fTrackletArray->GetEntriesFast()-1])->SetQ1(q1);
      }
    }
  }
}

Int_t AliTRDmcmSim::GetPID(Float_t q0, Float_t q1) 
{
  // get PID from accumulated charges q0 and q1

  Int_t binQ0 = (Int_t) (q0 * fgPidNBinsQ0) + 1;
  Int_t binQ1 = (Int_t) (q1 * fgPidNBinsQ1) + 1;
  binQ0 = binQ0 >= fgPidNBinsQ0 ? fgPidNBinsQ0-1 : binQ0;
  binQ1 = binQ1 >= fgPidNBinsQ0 ? fgPidNBinsQ0-1 : binQ1;

  return fgPidLut[binQ0*fgPidNBinsQ1+binQ1];
}

void AliTRDmcmSim::SetPIDlut(Int_t *lut, Int_t nbinsq0, Int_t nbinsq1)
{
  // set a user-defined PID LUT

  if (fgPidLutDelete)
    delete [] fgPidLut;

  fgPidLutDelete = kFALSE;
  fgPidLut = lut;
  fgPidNBinsQ0 = nbinsq0;
  fgPidNBinsQ1 = nbinsq1;
}

void AliTRDmcmSim::SetPIDlut(TH2F *lut)
{
  // set a user-defined PID LUT from a 2D histogram

  if (fgPidLutDelete)
    delete [] fgPidLut;

  fgPidNBinsQ0 = lut->GetNbinsX();
  fgPidNBinsQ1 = lut->GetNbinsY();

  fgPidLut = new Int_t[fgPidNBinsQ0*fgPidNBinsQ1];

  for (Int_t ix = 0; ix < fgPidNBinsQ0; ix++) {
    for (Int_t iy = 0; iy < fgPidNBinsQ1; iy++) {
      fgPidLut[ix*fgPidNBinsQ1 + iy] = (Int_t) (256. * lut->GetBinContent(ix, iy));
    }
  }

  fgPidLutDelete = kTRUE;
}

void AliTRDmcmSim::SetPIDlutDefault()
{
  // use the default PID LUT

  if (fgPidLutDelete )
    delete [] fgPidLut;

  fgPidLutDelete = kFALSE;
  fgPidLut = *fgPidLutDefault;
  fgPidNBinsQ0 = 40;
  fgPidNBinsQ1 = 50;
}

void AliTRDmcmSim::Tracklet()
{
  // Run the tracklet calculation by calling sequentially:
  // CalcFitreg(); TrackletSelection(); FitTracklet()
  // and store the tracklets 

  if (!fInitialized) {
    AliError("Called uninitialized! Nothing done!");
    return;
  }

  fTrackletArray->Delete();

  CalcFitreg();
  if (fNHits == 0)
    return;
  TrackletSelection();
  FitTracklet();
}

Bool_t AliTRDmcmSim::StoreTracklets() 
{
  // store the found tracklets via the loader

  if (fTrackletArray->GetEntriesFast() == 0) 
    return kTRUE;

  AliRunLoader *rl = AliRunLoader::Instance();
  AliDataLoader *dl = 0x0;
  if (rl)
    dl = rl->GetLoader("TRDLoader")->GetDataLoader("tracklets");
  if (!dl) {
    AliError("Could not get the tracklets data loader!");
    return kFALSE;
  }

  TTree *trackletTree = dl->Tree();
  if (!trackletTree) {
    dl->MakeTree();
    trackletTree = dl->Tree();
  }
  
  AliTRDtrackletMCM *trkl = 0x0;
  TBranch *trkbranch = trackletTree->GetBranch("mcmtrklbranch");
  if (!trkbranch)
    trkbranch = trackletTree->Branch("mcmtrklbranch", "AliTRDtrackletMCM", &trkl, 32000);
  
  for (Int_t iTracklet = 0; iTracklet < fTrackletArray->GetEntriesFast(); iTracklet++) {
    trkl = ((AliTRDtrackletMCM*) (*fTrackletArray)[iTracklet]);
    trkbranch->SetAddress(&trkl);
//      printf("filling tracklet 0x%08x\n", trkl->GetTrackletWord());
    trkbranch->Fill();
  }
  dl->WriteData("OVERWRITE");

  return kTRUE;
}

void AliTRDmcmSim::WriteData(AliTRDarrayADC *digits)
{
  // write back the processed data configured by EBSF
  // EBSF = 1: unfiltered data; EBSF = 0: filtered data
  // zero-suppressed valued are written as -1 to digits

  if (!fInitialized) {
    AliError("Called uninitialized! Nothing done!");
    return;
  }

//  Int_t firstAdc = 0;
//  Int_t lastAdc = fNADC - 1;
//
//  while (GetCol(firstAdc) < 0)
//    firstAdc++;
//
//  while (GetCol(lastAdc) < 0) 
//    lastAdc--;

  Int_t offset = (fMcmPos % 4) * 21 + (fRobPos % 2) * 84;

  if (fTrapConfig->GetTrapReg(AliTRDtrapConfig::kEBSF) != 0) // store unfiltered data
  {
    for (Int_t iAdc = 0; iAdc < fNADC; iAdc++) {
      if (fZSM1Dim[iAdc] == 1) {
        for (Int_t iTimeBin = 0; iTimeBin < fNTimeBin; iTimeBin++) {
          digits->SetDataByAdcCol(GetRow(), 20-iAdc + offset, iTimeBin, -1);
//          printf("suppressed: %i, %i, %i, %i, now: %i\n", fDetector, GetRow(), GetCol(iAdc), iTimeBin, 
//                 digits->GetData(GetRow(), GetCol(iAdc), iTimeBin));
        }
      }
    }
  }
  else {
    for (Int_t iAdc = 0; iAdc < fNADC; iAdc++) {
      if (fZSM1Dim[iAdc] == 0) {
        for (Int_t iTimeBin = 0; iTimeBin < fNTimeBin; iTimeBin++) {
          digits->SetDataByAdcCol(GetRow(), 20-iAdc + offset, iTimeBin, (fADCF[iAdc][iTimeBin] >> fgkAddDigits) - fgAddBaseline);
        }
      }
      else {
        for (Int_t iTimeBin = 0; iTimeBin < fNTimeBin; iTimeBin++) {
          digits->SetDataByAdcCol(GetRow(), 20-iAdc + offset, iTimeBin, -1);
//          printf("suppressed: %i, %i, %i, %i\n", fDetector, GetRow(), GetCol(iAdc), iTimeBin);
        }
      }
    }
  }
}

// help functions, to be cleaned up

UInt_t AliTRDmcmSim::AddUintClipping(UInt_t a, UInt_t b, UInt_t nbits) const
{
  // 
  // This function adds a and b (unsigned) and clips to 
  // the specified number of bits. 
  //  

  UInt_t sum = a + b;
  if (nbits < 32)
  {
    UInt_t maxv = (1 << nbits) - 1;;
    if (sum > maxv) 
      sum = maxv;
  }
  else
  {
    if ((sum < a) || (sum < b)) 
      sum = 0xFFFFFFFF;
  }
  return sum;
}

void AliTRDmcmSim::Sort2(UShort_t  idx1i, UShort_t  idx2i, \
                            UShort_t  val1i, UShort_t  val2i, \
                            UShort_t *idx1o, UShort_t *idx2o, \
                            UShort_t *val1o, UShort_t *val2o) const
{
  // sorting for tracklet selection

    if (val1i > val2i)
    {
        *idx1o = idx1i;
        *idx2o = idx2i;
        *val1o = val1i;
        *val2o = val2i;
    }
    else
    {
        *idx1o = idx2i;
        *idx2o = idx1i;
        *val1o = val2i;
        *val2o = val1i;
    }
}

void AliTRDmcmSim::Sort3(UShort_t  idx1i, UShort_t  idx2i, UShort_t  idx3i, \
                            UShort_t  val1i, UShort_t  val2i, UShort_t  val3i, \
                            UShort_t *idx1o, UShort_t *idx2o, UShort_t *idx3o, \
                            UShort_t *val1o, UShort_t *val2o, UShort_t *val3o)
{
  // sorting for tracklet selection

    Int_t sel;


    if (val1i > val2i) sel=4; else sel=0;
    if (val2i > val3i) sel=sel + 2;
    if (val3i > val1i) sel=sel + 1;
    //printf("input channels %d %d %d, charges %d %d %d sel=%d\n",idx1i, idx2i, idx3i, val1i, val2i, val3i, sel);
    switch(sel)
    {
        case 6 : // 1 >  2  >  3            => 1 2 3
        case 0 : // 1 =  2  =  3            => 1 2 3 : in this case doesn't matter, but so is in hardware!
            *idx1o = idx1i;
            *idx2o = idx2i;
            *idx3o = idx3i;
            *val1o = val1i;
            *val2o = val2i;
            *val3o = val3i;
            break;

        case 4 : // 1 >  2, 2 <= 3, 3 <= 1  => 1 3 2
            *idx1o = idx1i;
            *idx2o = idx3i;
            *idx3o = idx2i;
            *val1o = val1i;
            *val2o = val3i;
            *val3o = val2i;
            break;

        case 2 : // 1 <= 2, 2 > 3, 3 <= 1   => 2 1 3
            *idx1o = idx2i;
            *idx2o = idx1i;
            *idx3o = idx3i;
            *val1o = val2i;
            *val2o = val1i;
            *val3o = val3i;
            break;

        case 3 : // 1 <= 2, 2 > 3, 3  > 1   => 2 3 1
            *idx1o = idx2i;
            *idx2o = idx3i;
            *idx3o = idx1i;
            *val1o = val2i;
            *val2o = val3i;
            *val3o = val1i;
            break;

        case 1 : // 1 <= 2, 2 <= 3, 3 > 1   => 3 2 1
            *idx1o = idx3i;
            *idx2o = idx2i;
            *idx3o = idx1i;
            *val1o = val3i;
            *val2o = val2i;
            *val3o = val1i;
        break;

        case 5 : // 1 > 2, 2 <= 3, 3 >  1   => 3 1 2
            *idx1o = idx3i;
            *idx2o = idx1i;
            *idx3o = idx2i;
            *val1o = val3i;
            *val2o = val1i;
            *val3o = val2i;
        break;

        default: // the rest should NEVER happen!
            AliError("ERROR in Sort3!!!\n");
        break;
    }
//    printf("output channels %d %d %d, charges %d %d %d \n",*idx1o, *idx2o, *idx3o, *val1o, *val2o, *val3o);
}

void AliTRDmcmSim::Sort6To4(UShort_t  idx1i, UShort_t  idx2i, UShort_t  idx3i, UShort_t  idx4i, UShort_t  idx5i, UShort_t  idx6i, \
                               UShort_t  val1i, UShort_t  val2i, UShort_t  val3i, UShort_t  val4i, UShort_t  val5i, UShort_t  val6i, \
                               UShort_t *idx1o, UShort_t *idx2o, UShort_t *idx3o, UShort_t *idx4o, \
                               UShort_t *val1o, UShort_t *val2o, UShort_t *val3o, UShort_t *val4o)
{
  // sorting for tracklet selection

    UShort_t idx21s, idx22s, idx23s, dummy;
    UShort_t val21s, val22s, val23s;
    UShort_t idx23as, idx23bs;
    UShort_t val23as, val23bs;

    Sort3(idx1i, idx2i, idx3i, val1i, val2i, val3i,
                 idx1o, &idx21s, &idx23as,
                 val1o, &val21s, &val23as);

    Sort3(idx4i, idx5i, idx6i, val4i, val5i, val6i,
                 idx2o, &idx22s, &idx23bs,
                 val2o, &val22s, &val23bs);

    Sort2(idx23as, idx23bs, val23as, val23bs, &idx23s, &dummy, &val23s, &dummy);

    Sort3(idx21s, idx22s, idx23s, val21s, val22s, val23s,
                 idx3o, idx4o, &dummy,
                 val3o, val4o, &dummy);

}

void AliTRDmcmSim::Sort6To2Worst(UShort_t  idx1i, UShort_t  idx2i, UShort_t  idx3i, UShort_t  idx4i, UShort_t  idx5i, UShort_t  idx6i, \
                                    UShort_t  val1i, UShort_t  val2i, UShort_t  val3i, UShort_t  val4i, UShort_t  val5i, UShort_t  val6i, \
                                    UShort_t *idx5o, UShort_t *idx6o)
{
  // sorting for tracklet selection

    UShort_t idx21s, idx22s, idx23s, dummy1, dummy2, dummy3, dummy4, dummy5;
    UShort_t val21s, val22s, val23s;
    UShort_t idx23as, idx23bs;
    UShort_t val23as, val23bs;

    Sort3(idx1i, idx2i,   idx3i, val1i, val2i, val3i,
                 &dummy1, &idx21s, &idx23as,
                 &dummy2, &val21s, &val23as);

    Sort3(idx4i, idx5i, idx6i, val4i, val5i, val6i,
                 &dummy1, &idx22s, &idx23bs,
                 &dummy2, &val22s, &val23bs);

    Sort2(idx23as, idx23bs, val23as, val23bs, &idx23s, idx5o, &val23s, &dummy1);

    Sort3(idx21s, idx22s, idx23s, val21s, val22s, val23s,
                 &dummy1, &dummy2, idx6o,
                 &dummy3, &dummy4, &dummy5);
//    printf("idx21s=%d, idx23as=%d, idx22s=%d, idx23bs=%d, idx5o=%d, idx6o=%d\n",
//            idx21s,    idx23as,    idx22s,    idx23bs,    *idx5o,    *idx6o);
}


