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

////////////////////////////////////////////////////////////////////////////
//                                                                        
//  Parameters for Pre-Trigger Simulation                                 
//                                                                        
//  Author: F. Reidt (Felix.Reidt@cern.ch)                                
//               
//  This class controls the parameters used by the pretrigger simulation.
//  A configuration file ca be loaded by calling LoadConfigurationFromFile()
//  The generation of look up tables is also done in this class and has to
//  be done only once after a configuration file was loaded.
//  If no configuration file was loaded, the standard
//  configuration in LoadStandardConfiguration() would be used.
//                                                        
////////////////////////////////////////////////////////////////////////////
#include "TArrayI.h"
#include "TObjString.h"
#include "TString.h"

#include <fstream>

#include "AliLog.h"

#include "AliTRDptrgParam.h"

ClassImp(AliTRDptrgParam)

AliTRDptrgParam *AliTRDptrgParam::fgInstance = 0;


//______________________________________________________________________________
AliTRDptrgParam::AliTRDptrgParam()
  : TObject(),
    fTLMUInputStretch(0),
    fTLMUcmatrices(0x0),
    fTLMUmultiplicity(0x0), 
    fTLMUoutput(0x0),
    fFEBT0Thresholds(0x0),
    fFEBT0Multiplicities(0x0),
    fFEBT0LUTs(0x0),
    fFEBV0Thresholds(0x0),
    fFEBV0Multiplicities(0x0),
    fFEBV0LUTs(0x0),
    fCBLUTs(0x0),
    fPTmasks(AliTRDptrgPTmasks())
{
  // ctor
  
  // initialize coincidence matrices
  this->fTLMUcmatrices = new UInt_t*[3];
  for (UInt_t iMatrix = 0; iMatrix < 3; iMatrix++) {
    this->fTLMUcmatrices[iMatrix] = new UInt_t[18];
    for (UInt_t iSlice = 0; iSlice < 18; iSlice++) {
      this->fTLMUcmatrices[iMatrix][iSlice] = 0;
    }
  }
  
  // initialize multiplicity slices
  this->fTLMUmultiplicity = new UInt_t*[9];
  for (UInt_t iSlice = 0; iSlice < 9; iSlice++) {
    this->fTLMUmultiplicity[iSlice] = new UInt_t[2];
    this->fTLMUmultiplicity[iSlice][0] = 577; // disabled
    this->fTLMUmultiplicity[iSlice][1] = 0;
  }
  
  // initialize output muxer
  this->fTLMUoutput = new Int_t*[8];
  for (UInt_t iBit = 0; iBit < 8; iBit++) {
    this->fTLMUoutput[iBit] = new Int_t[2];
    this->fTLMUoutput[iBit][0] = -1; // cmatrix disabled
    this->fTLMUoutput[iBit][1] = -1; // multslice disabled
  }
  
  // initialize T0 FEB thresholds
  this->fFEBT0Thresholds = new UInt_t*[2];
  this->fFEBT0Thresholds[0] = new UInt_t[12];
  this->fFEBT0Thresholds[1] = new UInt_t[12];
  for (Int_t iChan = 0; iChan < 12; iChan++) {
    this->fFEBT0Thresholds[0][iChan] = 4294967295U; 
    this->fFEBT0Thresholds[1][iChan] = 4294967295U;
    // writing 2^32-1 disables the input because all used adcs have 
    // less than 32 bits
  }
  
  // initialize T0 Multiplicity
  this->fFEBT0Multiplicities = new UInt_t**[2];
  this->fFEBT0Multiplicities[0] = new UInt_t*[2];
  this->fFEBT0Multiplicities[1] = new UInt_t*[2];
  this->fFEBT0Multiplicities[0][0] = new UInt_t[2];
  this->fFEBT0Multiplicities[0][1] = new UInt_t[2];
  this->fFEBT0Multiplicities[1][0] = new UInt_t[2];
  this->fFEBT0Multiplicities[1][1] = new UInt_t[2];
  this->fFEBT0Multiplicities[0][0][0] = 4294967295U;
  this->fFEBT0Multiplicities[0][0][1] = 4294967295U;
  this->fFEBT0Multiplicities[0][1][0] = 4294967295U;
  this->fFEBT0Multiplicities[0][1][1] = 4294967295U;
  this->fFEBT0Multiplicities[1][0][0] = 4294967295U;
  this->fFEBT0Multiplicities[1][0][1] = 4294967295U;
  this->fFEBT0Multiplicities[1][1][0] = 4294967295U;
  this->fFEBT0Multiplicities[1][1][1] = 4294967295U;
  // writing 2^32-1 disables the input because all used adcs have 
  // less than 32 bits

  // initialize T0 LUTs
  // this->fFEBT0LUTs = 0x0; (done in member initialization list)
  // further initialization is done in AliTRDptrgParam::GenerateLUTs()
  

  // initialize V0 FEB Thresholds
  this->fFEBV0Thresholds = new UInt_t**[2];
  for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
    this->fFEBV0Thresholds[iPosition] = new UInt_t*[4];
    for (UInt_t iCard = 0; iCard < 4; iCard++) {
      this->fFEBV0Thresholds[iPosition][iCard] = new UInt_t[8];
      for (UInt_t iChannel = 0; iChannel < 8; iChannel++) {
        this->fFEBV0Thresholds[iPosition][iCard][iChannel] = 4294967295U;
      }
    }
  }
 
  // initialize V0 Multiplicities
  this->fFEBV0Multiplicities = new UInt_t***[2];
  for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
    this->fFEBV0Multiplicities[iPosition] = new UInt_t**[4];
    for (UInt_t iCard = 0; iCard < 4; iCard++) {
      this->fFEBV0Multiplicities[iPosition][iCard] = new UInt_t*[2];
      for (UInt_t iLUT = 0; iLUT < 2; iLUT++) {
        this->fFEBV0Multiplicities[iPosition][iCard][iLUT] = new UInt_t[2];
        this->fFEBV0Multiplicities[iPosition][iCard][iLUT][0] = 4294967295U;
        this->fFEBV0Multiplicities[iPosition][iCard][iLUT][1] = 0x0;      
      }
    }
  }
  
  // initialize V0 LUTs
  //  this->fFEBV0LUTs = 0x0; (done in member initialization list)
  // further initialization is done in AliTRDptrgParam::GenerateLUTs()

  // initialize CB LUTs
  // this->fCBLUTs = 0x0; (done in member initialization list)
  // further initialization is done in AliTRDptrgParam::GenerateLUTs()

  this->LoadStandardConfiguration(); // load standard configuration
}

//______________________________________________________________________________
AliTRDptrgParam::~AliTRDptrgParam() 
{
  // dtor
  
  // delete coincidence matrices
  if (this->fTLMUcmatrices != 0x0) {
    for (UInt_t iMatrix = 0; iMatrix < 3; iMatrix++) {
      if (this->fTLMUcmatrices[iMatrix] != 0x0) {
        delete[] this->fTLMUcmatrices[iMatrix];
        this->fTLMUcmatrices[iMatrix] = 0x0;
      }
    }
    delete[] this->fTLMUcmatrices;
    this->fTLMUcmatrices = 0x0;
  }
  
  // delete multiplicity slices
  if (this->fTLMUmultiplicity != 0x0) {
    for (UInt_t iSlice = 0; iSlice < 9; iSlice++) {
      if (this->fTLMUmultiplicity[iSlice] != 0x0) {
        delete[] this->fTLMUmultiplicity[iSlice];
        this->fTLMUmultiplicity[iSlice] = 0x0;
      }
    }
    delete[] this->fTLMUmultiplicity;
    this->fTLMUmultiplicity = 0x0;
  }

  // delete output mux
  if (this->fTLMUoutput != 0x0) {
    for (UInt_t iBit = 0; iBit < 8; iBit++) {
      if (this->fTLMUoutput[iBit] != 0x0) {
        delete[] this->fTLMUoutput[iBit];
        this->fTLMUoutput[iBit] = 0x0;
      }
    }
    delete[] this->fTLMUoutput;
    this->fTLMUoutput = 0x0;
  }

  // delete T0 FEB thresholds
  if (this->fFEBT0Thresholds != 0x0) {
    if (this->fFEBT0Thresholds[0] != 0x0) {
      delete[] this->fFEBT0Thresholds[0];
      this->fFEBT0Thresholds[0] = 0x0;
    }
    if (this->fFEBT0Thresholds[1] != 0x0) {
      delete[] this->fFEBT0Thresholds[1];
      this->fFEBT0Thresholds[1] = 0x0;
    }
    delete[] this->fFEBT0Thresholds;
    this->fFEBT0Thresholds = 0x0;
  }
 
  // delete T0 multiplicities
  if (this->fFEBT0Multiplicities != 0x0) {
    for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
      if (this->fFEBT0Multiplicities[iPosition] != 0x0) {
        for (UInt_t iLUT = 0; iLUT < 2; iLUT++) {
          if (this->fFEBT0Multiplicities[iPosition][iLUT] != 0x0) {
            delete[] this->fFEBT0Multiplicities[iPosition][iLUT];
            this->fFEBT0Multiplicities[iPosition][iLUT] = 0x0;
          }
	}
        delete[] this->fFEBT0Multiplicities[iPosition];
	this->fFEBT0Multiplicities[iPosition] = 0x0;
      }
    }
    delete[] this->fFEBT0Multiplicities;
    this->fFEBT0Multiplicities = 0x0;
  }  

  // delete T0 LUTs
  if (this->fFEBT0LUTs != 0x0) {
    for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
      if (this->fFEBT0LUTs[iPosition] != 0x0) {
        for (UInt_t iLUT = 0; iLUT < 2; iLUT++) {
          if (this->fFEBT0LUTs[iPosition][iLUT] != 0x0) {
            delete[] this->fFEBT0LUTs[iPosition][iLUT];
            this->fFEBT0LUTs[iPosition][iLUT] = 0x0;
          }
	}
        delete[] this->fFEBT0LUTs[iPosition];
	this->fFEBT0LUTs[iPosition] = 0x0;
      }
    }
    delete[] this->fFEBT0LUTs;
    this->fFEBT0LUTs = 0x0;
  }  

  // delete V0 FEB thresholds
  if (this->fFEBV0Thresholds != 0x0) {
    for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
      if (this->fFEBV0Thresholds[iPosition] != 0x0) {
        for (UInt_t iCard = 0; iCard < 4; iCard++) {
          if (this->fFEBV0Thresholds[iPosition][iCard] != 0x0) {
            delete[] this->fFEBV0Thresholds[iPosition][iCard];
            this->fFEBV0Thresholds[iPosition][iCard] = 0x0;
	  }
	}
        delete[] this->fFEBV0Thresholds[iPosition]; 
        this->fFEBV0Thresholds[iPosition] = 0x0;
      }
    }
    delete[] this->fFEBV0Thresholds;
    this->fFEBV0Thresholds = 0x0;
  }

  // delete V0 multiplicities
  if (this->fFEBV0Multiplicities != 0x0) {
    for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
      if (this->fFEBV0Multiplicities[iPosition] != 0x0) {
        for (UInt_t iCard = 0; iCard < 4; iCard++) {
          if (this->fFEBV0Multiplicities[iPosition][iCard] != 0x0) {
            for (UInt_t iLUT = 0; iLUT < 2; iLUT++) {
              if (this->fFEBV0Multiplicities[iPosition][iCard][iLUT] != 0x0) {
                delete[] this->fFEBV0Multiplicities[iPosition][iCard][iLUT];
                this->fFEBV0Multiplicities[iPosition][iCard][iLUT] = 0x0;
	      }
	    }
            delete[] this->fFEBV0Multiplicities[iPosition][iCard];
            this->fFEBV0Multiplicities[iPosition][iCard] = 0x0;
	  }
	}
        delete[] this->fFEBV0Multiplicities[iPosition]; 
        this->fFEBV0Multiplicities[iPosition] = 0x0;
      }
    }
    delete[] this->fFEBV0Multiplicities;
    this->fFEBV0Multiplicities = 0x0;
  } 

  // delete V0 LUTs
  if (this->fFEBV0LUTs != 0x0) {
    for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
      if (this->fFEBV0LUTs[iPosition] != 0x0) {
        for (UInt_t iCard = 0; iCard < 4; iCard++) {
          if (this->fFEBV0LUTs[iPosition][iCard] != 0x0) {
            for (UInt_t iLUT = 0; iLUT < 2; iLUT++) {
              if (this->fFEBV0LUTs[iPosition][iCard][iLUT] != 0x0) {
                delete[] this->fFEBV0LUTs[iPosition][iCard][iLUT];
                this->fFEBV0LUTs[iPosition][iCard][iLUT] = 0x0;
	      }
	    }
            delete[] this->fFEBV0LUTs[iPosition][iCard];
            this->fFEBV0LUTs[iPosition][iCard] = 0x0;
	  }
	}
        delete[] this->fFEBV0LUTs[iPosition]; 
        this->fFEBV0LUTs[iPosition] = 0x0;
      }
    }
    delete[] this->fFEBV0LUTs;
    this->fFEBV0LUTs = 0x0;
  } 

  // delete CB LUTs
  if (this->fCBLUTs != 0x0) {
    for (UInt_t iCB = 0; iCB < 3; iCB++) {
      if (this->fCBLUTs[iCB] != 0x0) {
        for (UInt_t iLUT = 0; iLUT < 2; iLUT++) {
          if (this->fCBLUTs[iCB][iLUT] != 0x0) {
            delete[] this->fCBLUTs[iCB][iLUT];
            this->fCBLUTs[iCB][iLUT] = 0x0;
          }
	}
        if (iCB == kB) {
          // CB-B has 3 LUTs!
          if (this->fCBLUTs[iCB][2] != 0x0) {
            delete[] this->fCBLUTs[iCB][2];
            this->fCBLUTs[iCB][2] = 0x0;
          }
        }
        delete[] this->fCBLUTs[iCB];
	this->fCBLUTs[iCB] = 0x0;
      }
    }
    delete[] this->fCBLUTs;
    this->fCBLUTs = 0x0;
  }  
}

//______________________________________________________________________________
Int_t AliTRDptrgParam::CheckVariables() const
{
  // checks whether variables are deleted early enough
 
  // check coincidence matrices
  if (this->fTLMUcmatrices != 0x0) {
    for (UInt_t iMatrix = 0; iMatrix < 3; iMatrix++) {
      if (this->fTLMUcmatrices[iMatrix] == 0x0) {
        return -1;
      }
    }
  }
  else {
    return -2;
  }
  
  // check multiplicity slices
  if (this->fTLMUmultiplicity != 0x0) {
    for (UInt_t iSlice = 0; iSlice < 9; iSlice++) {
      if (this->fTLMUmultiplicity[iSlice] == 0x0) {
        return -3;
      }
    }
  }
  else {
    return -4;
  }

  // check output mux
  if (this->fTLMUoutput != 0x0) {
    for (UInt_t iBit = 0; iBit < 8; iBit++) {
      if (this->fTLMUoutput[iBit] == 0x0) {
        return -5;
      }
    }
  }
  else {
    return -6;
  }

  // check T0 FEB thresholds
  if (this->fFEBT0Thresholds != 0x0) {
    if (this->fFEBT0Thresholds[0] == 0x0) {
      return -7;
    }
    if (this->fFEBT0Thresholds[1] == 0x0) {
      return -8;
    }
  }
  else {
    return -9;
  }

  // check T0 multiplicities
  if (this->fFEBT0Multiplicities != 0x0) {
    for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
      if (this->fFEBT0Multiplicities[iPosition] != 0x0) {
        for (UInt_t iLUT = 0; iLUT < 2; iLUT++) {
          if (this->fFEBT0Multiplicities[iPosition][iLUT] == 0x0) {
            return -10;
          }
	}
      }
      else {
        return -11;
      }
    }
  }
  else {
    return -12;
  }
  

  // check T0 LUTs
  if (this->fFEBT0LUTs != 0x0) {
    for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
      if (this->fFEBT0LUTs[iPosition] != 0x0) {
        for (UInt_t iLUT = 0; iLUT < 2; iLUT++) {
          if (this->fFEBT0LUTs[iPosition][iLUT] == 0x0) {
            return -13;
          }
	}
      }
      else {
        return -14;
      }
    }
  }  
  else {
    return -15;
  }

  // check V0 FEB thresholds
  if (this->fFEBV0Thresholds != 0x0) {
    for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
      if (this->fFEBV0Thresholds[iPosition] != 0x0) {
        for (UInt_t iCard = 0; iCard < 4; iCard++) {
          if (this->fFEBV0Thresholds[iPosition][iCard] == 0x0) {
            return -16;
	  }
	}
      }
      else {
        return -17;
      }
    }
  }
  else {
    return -18;
  }

  // check V0 multiplicities
  if (this->fFEBV0Multiplicities != 0x0) {
    for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
      if (this->fFEBV0Multiplicities[iPosition] != 0x0) {
        for (UInt_t iCard = 0; iCard < 4; iCard++) {
          if (this->fFEBV0Multiplicities[iPosition][iCard] != 0x0) {
            for (UInt_t iLUT = 0; iLUT < 2; iLUT++) {
              if (this->fFEBV0Multiplicities[iPosition][iCard][iLUT] == 0x0) {
                return -19;
	      }
	    }
	  }
          else {
            return -20;
	  }
	}
      }
      else {
        return -21;
      }
    }
  } 
  else {
    return -22;
  }

  // check V0 LUTs
  if (this->fFEBV0LUTs != 0x0) {
    for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
      if (this->fFEBV0LUTs[iPosition] != 0x0) {
        for (UInt_t iCard = 0; iCard < 4; iCard++) {
          if (this->fFEBV0LUTs[iPosition][iCard] != 0x0) {
            for (UInt_t iLUT = 0; iLUT < 2; iLUT++) {
              if (this->fFEBV0LUTs[iPosition][iCard][iLUT] == 0x0) {
                return -23;
	      }
	    }
	  }
          else {
            return -24;
	  }
	}
      }
      else {
        return -25;
      }
    }
  } 
  else {
    return -26;
  }

  // check CB LUTs
  if (this->fCBLUTs != 0x0) {
    for (UInt_t iCB = 0; iCB < 3; iCB++) {
      if (this->fCBLUTs[iCB] != 0x0) {
        for (UInt_t iLUT = 0; iLUT < 2; iLUT++) {
          if (this->fCBLUTs[iCB][iLUT] == 0x0) {
            return -27;
          }
	}
        if (iCB == kB) {
          if (this->fCBLUTs[iCB][2] == 0x0) {
            return -28;
          }
        }
      }
      else {
        return -29;
      }
    }
  }  
  else {
    return -30;
  }
  return 0;
}

//______________________________________________________________________________
AliTRDptrgParam* AliTRDptrgParam::Instance() 
{
  // get (or create) the single instance

  if (fgInstance == 0) 
    fgInstance = new AliTRDptrgParam();

  return fgInstance;
}

//______________________________________________________________________________
void AliTRDptrgParam::Terminate() 
{
  // destruct the instance

  if (fgInstance != 0) {
    delete fgInstance;
    fgInstance = 0x0;
  }
}

//______________________________________________________________________________
void AliTRDptrgParam::LoadStandardConfiguration() {
  // loads a standard configuration parameters for testing 

  // TLMU Input Masks
  for (UInt_t iSM = 0; iSM < 18; iSM++) {
    this->fTLMUInputMask[iSM] = 0xFFFFFFFF; // enable all input bits
  }
 
  // TLMU Input Stretch
  this->fTLMUInputStretch = 0; // not used in simulation

  // TLMU Coincidence Matrices
  //
  // Matrix 0: Back-To-Back
  // Matrix 1: Back-To-Back +/-1
  // Matrix 2: Back-To-Back +/-2
  for (UInt_t iMatrix = 0; iMatrix < 3; iMatrix++) {
    for (UInt_t iSlice = 0; iSlice < 18; iSlice++) {
      if (iMatrix == 0) {
        if (iSlice < 9) {
          this->fTLMUcmatrices[iMatrix][iSlice] = 0x201 << iSlice; 
          // Back-To-Back 
          AliDebug(5, Form("fTLMUcmatrices[%d][%d]=0x%x",iMatrix,iSlice,
                           this->fTLMUcmatrices[iMatrix][iSlice]));
	}
        // because of symmetrie the other slices are not necessary
      } 
      else if (iMatrix == 1)  {
        // Back-To-Back +/- 1
        if (iSlice < 8) {
          this->fTLMUcmatrices[iMatrix][iSlice] = 0x381 << iSlice;
        }
        else if (iSlice == 8) {
          this->fTLMUcmatrices[iMatrix][iSlice] = 0x30101;
        }
        else if (iSlice == 9) {
          this->fTLMUcmatrices[iMatrix][iSlice] = 0x20203;
        }
        else {
          this->fTLMUcmatrices[iMatrix][iSlice] = 0x407 << (iSlice - 10);
        } 
        AliDebug(5, Form("fTLMUcmatrices[%d][%d]=0x%x",iMatrix,iSlice,
                         this->fTLMUcmatrices[iMatrix][iSlice])); 
      }
      else if (iMatrix == 2) {
        // Back-To-Back +/-2
        if (iSlice < 7 ) {
          this->fTLMUcmatrices[iMatrix][iSlice] = 0xF81 << iSlice;
        }
        else if (iSlice == 7) {
          this->fTLMUcmatrices[iMatrix][iSlice] = 0x3C081;
        }
        else if (iSlice == 8) {
          this->fTLMUcmatrices[iMatrix][iSlice] = 0x38103;
        }
        else if (iSlice == 9) {
          this->fTLMUcmatrices[iMatrix][iSlice] = 0x30207;
        }
        else if (iSlice == 10) {
          this->fTLMUcmatrices[iMatrix][iSlice] = 0x2040F;
        }
        else {
          this->fTLMUcmatrices[iMatrix][iSlice] = 0x81F << (iSlice - 11);
        } 
        AliDebug(5, Form("fTLMUcmatrices[%d][%d]=0x%x",iMatrix,iSlice,
                         this->fTLMUcmatrices[iMatrix][iSlice]));     
      }
    } 
  }

  // TLMU Mulitplicity
  this->fTLMUmultiplicity[0][0] = 0;
  this->fTLMUmultiplicity[0][1] = 10;
  this->fTLMUmultiplicity[1][0] = 10;
  this->fTLMUmultiplicity[1][1] = 25;
  this->fTLMUmultiplicity[2][0] = 25;
  this->fTLMUmultiplicity[2][1] = 50;
  this->fTLMUmultiplicity[3][0] = 50;
  this->fTLMUmultiplicity[3][1] = 100;
  this->fTLMUmultiplicity[4][0] = 100;
  this->fTLMUmultiplicity[4][1] = 200;
  this->fTLMUmultiplicity[5][0] = 200;
  this->fTLMUmultiplicity[5][1] = 350;
  this->fTLMUmultiplicity[6][0] = 350;
  this->fTLMUmultiplicity[6][1] = 400;
  this->fTLMUmultiplicity[7][0] = 400;
  this->fTLMUmultiplicity[7][1] = 576;
  this->fTLMUmultiplicity[8][0] = 1;
  this->fTLMUmultiplicity[8][1] = 576;
 
  // TLMU output
  this->fTLMUoutput[0][0] = 0;
  this->fTLMUoutput[1][0] = 1;
  this->fTLMUoutput[2][0] = 2;
  this->fTLMUoutput[3][1] = 0;
  this->fTLMUoutput[4][1] = 1;
  this->fTLMUoutput[5][1] = 2;
  this->fTLMUoutput[6][1] = 3;
  this->fTLMUoutput[7][1] = 8;

  // T0 FEB Thresholds
  for (UInt_t iChannel = 0; iChannel < 12; iChannel++) {
    this->fFEBT0Thresholds[0][iChannel] = 10;
    this->fFEBT0Thresholds[1][iChannel] = 10;
  }

  // T0 Multiplicities
  for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
    for (UInt_t iLUT = 0; iLUT < 2; iLUT++) {
      if (iLUT == 0) {
        this->fFEBT0Multiplicities[iPosition][iLUT][0] = 0;
      }
      else {
	this->fFEBT0Multiplicities[iPosition][iLUT][0] = 5;
      }
      this->fFEBT0Multiplicities[iPosition][iLUT][1] = 0xFFF;
    }
  }

  // V0 FEB Thresholds
  for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
    for (UInt_t iCard = 0; iCard < 4; iCard++) {
      for (UInt_t iChannel = 0; iChannel < 8; iChannel++) {
        this->fFEBV0Thresholds[iPosition][iCard][iChannel] = 10;
      }
    }
  }

  // V0 Multiplicities
  for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
    for (UInt_t iCard = 0; iCard < 4; iCard++) {
      for (UInt_t iLUT = 0; iLUT < 2; iLUT++) {
        if (iLUT == 0) {
          this->fFEBV0Multiplicities[iPosition][iCard][iLUT][0] = 0;
        }
        else {
          this->fFEBV0Multiplicities[iPosition][iCard][iLUT][0] = 3;
        }
        this->fFEBV0Multiplicities[iPosition][iCard][iLUT][1] = 0xFF;
      }
    }
  }

  // CB-A LUT equations
  this->fCBALUTequ[0] = "T0_0 || (V0-0_0 || V0-1_0 || V0-2_0 || V0-3_0)";
  this->fCBALUTequ[1] = "!T0_1 && !V0-0_1 && !V0-1_1 && !V0-2_1 && !V0-3_1";

  // CB-C LUT equations
  this->fCBCLUTequ[0] = "T0_0 || ( V0-0_0 || V0-1_0 || V0-2_0 || V0-3_0 )";
  this->fCBCLUTequ[1] = "!T0_1 && !V0-0_1 && !V0-1_1 && !V0-2_1 && !V0-3_1";

  // CB-B LUT equations
  this->fCBBLUTequ[0] = "CB-A_1 && !CB-C_1 && TLMU_7";
  this->fCBBLUTequ[1] = "!CB-A_1 && CB-C_1 && TLMU_7";
  this->fCBBLUTequ[2] = "CB-A_1 && CB-C_1 && TLMU_7";

  // PT output mask
  this->fPTmasks.fLUTs[0] = kTRUE;
  this->fPTmasks.fLUTs[1] = kTRUE;
  this->fPTmasks.fLUTs[2] = kTRUE;
  this->fPTmasks.fCBA[0] = kTRUE;
  this->fPTmasks.fCBC[0] = kTRUE;
  for (Int_t i = 1; i < 7; i++) {
    this->fPTmasks.fTLMU[i] = kTRUE;
  }

  return;  
}

//______________________________________________________________________________
Bool_t AliTRDptrgParam::LoadConfigurationFromFile(TString filename) {
  // Reads pretrigger configuration file and forwards identifiers and values
  // to the corresponding parser functions
  // This method is only checking for certain keywords at the beginning of a
  // line in the config file

  ifstream inputFile;
  inputFile.open(filename.Data());
  TString line;
  TString identifier;
  TString value;
  std::string str;


  if (inputFile.is_open())
  {
    AliDebug(5, "---- Reading configuration file ----");
    while (getline(inputFile, str)) {
      line = str;
   
      AliDebug(5, Form("line: %s\n", line.Data()));
      if (line.Index("TLMU") == 0) {
        this->PrepareLine(line, identifier, value);
        if (!this->ParseTLMU(identifier, value)) {
	  return kFALSE;
        }
      }
      else if (line.Index("FEB") == 0) {
        this->PrepareLine(line, identifier, value);
        if (!this->ParseFEB(identifier, value)) {
	  return kFALSE;
        }
      }
      else if (line.Index("CBB") == 0) {
        this->PrepareLine(line, identifier, value);
        if (!this->ParseCBB(identifier, value)) {
          return kFALSE;
	}
      }
      else if ((line.Index("CBA") == 0) ||
               (line.Index("CBC") == 0)) {
        this->PrepareLine(line, identifier, value);
        if (!this->ParseCBAC(identifier, value)) {
          return kFALSE;
	}
      }
    }
    AliDebug(5, "---- Finished reading configuration file ----");
    inputFile.close();
    return kTRUE;
  }
  else
  {
    AliDebug(5, "Error opening configuration file");
    return kFALSE;
  }
  return kTRUE;
}


//______________________________________________________________________________
Int_t AliTRDptrgParam::GenerateLUTs() {
  // generates all LUTs defined inside this object, this schould called only
  // once, after configuration is loaded in order to save cpu time

  // generation method:
  // walk through address space
  // mask address with input mask =>  get multiplicity of masked value
  // if (multiplicity of masked value) > multiplicity condition
  // write 1 in LUT
  
  // T0
  this->fFEBT0LUTs = new Int_t**[2]; // 2 A + C side
  for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
    // iPosition = 0 -> A, iPosition = 1 -> C
    this->fFEBT0LUTs[iPosition] = new Int_t*[2]; // 2 LUTs per side
    for (UInt_t iLUT = 0; iLUT < 2; iLUT++) {
      // LUT with 12 input channels 2^12=4096
      this->fFEBT0LUTs[iPosition][iLUT] = new Int_t[4096];
      AliDebug(5, Form("Generate FEBT0LUT[%d][%d]: (0x%x)>%d", iPosition, iLUT,
                       this->fFEBT0Multiplicities[iPosition][iLUT][1],
                       this->fFEBT0Multiplicities[iPosition][iLUT][0]));
      for (UInt_t iEntry = 0; iEntry < 4096; iEntry++) {
        // Check whether that entry belongs to a multiplicity exceeding the 
        // threshold
        if (this->GetMultiplicity(iEntry & 
	      this->fFEBT0Multiplicities[iPosition][iLUT][1]) > 
            this->fFEBT0Multiplicities[iPosition][iLUT][0]) {
          this->fFEBT0LUTs[iPosition][iLUT][iEntry] = 1;
	}
        else {
          // initialize LUT (not done before !)
          this->fFEBT0LUTs[iPosition][iLUT][iEntry] = 0;
        }
        AliDebug(10, Form("FEBT0LUTs[%d][%d][0x%x]=%d", iPosition, iLUT, iEntry,
                         this->fFEBT0LUTs[iPosition][iLUT][iEntry])); 
      }
      AliDebug(5, Form("Generated FEBT0LUTs[%d][%d]", iPosition, iLUT));
    }
  }

  // V0
  this->fFEBV0LUTs = new Int_t***[2]; // 2 A + C side
  for (UInt_t iPosition = 0; iPosition < 2; iPosition++) {
    // iPosition = 0 -> A, iPosition = 1 -> C
    this->fFEBV0LUTs[iPosition] = new Int_t**[4]; // 4 FEBs per side
    for (UInt_t iFEB = 0; iFEB < 4; iFEB++) {
      this->fFEBV0LUTs[iPosition][iFEB] = new Int_t*[2]; // 2 LUTs per FEB
      for (UInt_t iLUT = 0; iLUT < 2; iLUT++) {
      // LUT with 10 input channels 2^10=1024
        this->fFEBV0LUTs[iPosition][iFEB][iLUT] = new Int_t[1024];
        AliDebug(5, Form("Generate FEBV0LUT[%d][%d][%d]: (0x%x)>%d", iPosition, 
                         iFEB, iLUT,
                         this->fFEBV0Multiplicities[iPosition][iFEB][iLUT][1],
                         this->fFEBV0Multiplicities[iPosition][iFEB][iLUT][0]));
        for (UInt_t iEntry = 0; iEntry < 1024; iEntry++) {
          // Check whether that entry belongs to a multiplicity exceeding the 
          // threshold
          if (this->GetMultiplicity(iEntry & 
	        this->fFEBV0Multiplicities[iPosition][iFEB][iLUT][1]) > 
              this->fFEBV0Multiplicities[iPosition][iFEB][iLUT][0]) {
            this->fFEBV0LUTs[iPosition][iFEB][iLUT][iEntry] = 1;
	  }
          else {
            // initialize LUT (not done before !)
            this->fFEBV0LUTs[iPosition][iFEB][iLUT][iEntry] = 0;
          }
          AliDebug(10, Form("FEBV0LUTs[%d][%d][%d][0x%x]=%d", iPosition, iFEB,
                            iLUT, iEntry, 
                            this->fFEBV0LUTs[iPosition][iFEB][iLUT][iEntry]));
        }
        AliDebug(5, Form("Generated FEBV0LUTs[%d][%d][%d]", iPosition, iFEB, 
                         iLUT));
      }
    }
  }

  // ControlBoxes (CB-x)
  // initialize LUTs
  this->fCBLUTs = new Int_t**[3];
  for (Int_t iCB = 0; iCB < 3; iCB++) {
    if (iCB == kB) { // B
      fCBLUTs[iCB] = new Int_t*[3];
      this->fCBLUTs[iCB][0] = 0x0;
      this->fCBLUTs[iCB][1] = 0x0;
      this->fCBLUTs[iCB][2] = 0x0;
    }
    else { // A + C
      fCBLUTs[iCB] = new Int_t*[2];
      this->fCBLUTs[iCB][0] = 0x0;
      this->fCBLUTs[iCB][1] = 0x0;
    }
  }
  
  // CB-A (CB = 1 / kA)
  for (Int_t iLUT = 0; iLUT < 2; iLUT++) {
    this->fCBLUTs[1][iLUT] = this->GenerateLUTbasedOnEq(this->fCBALUTequ[iLUT],
                                                        10,1);
    for (Int_t iEntry = 0; iEntry < 1024; iEntry++) {
      AliDebug(10, Form("fCBLUTs[@A][%d][0x%x]=%d", iLUT, iEntry,
                        this->fCBLUTs[1][iLUT][iEntry]));
    }
  }

  // CB-C (CB = 2 / kC)
  for (Int_t iLUT = 0; iLUT < 2; iLUT++) {
    this->fCBLUTs[2][iLUT] = this->GenerateLUTbasedOnEq(this->fCBCLUTequ[iLUT],
                                                        10,1);
    for (Int_t iEntry = 0; iEntry < 1024; iEntry++) {
      AliDebug(6, Form("fCBLUTs[@C][%d][0x%x]=%d", iLUT, iEntry,
                       this->fCBLUTs[2][iLUT][iEntry]));
    }
  }  
 
  // CB-B (CB = 0 / kB)
  for (Int_t iLUT = 0; iLUT < 3; iLUT++) {
    this->fCBLUTs[0][iLUT] = this->GenerateLUTbasedOnEq(this->fCBBLUTequ[iLUT], 
                                                        12,1);
    for (Int_t iEntry = 0; iEntry < 4096; iEntry++) {
      AliDebug(10, Form("fCBLUTs[@B][%d][0x%x]=%d", iLUT, iEntry,
                        this->fCBLUTs[0][iLUT][iEntry]));
    }
  }

  AliDebug(5, "LUTs were generated!");
  return 0;
}

//______________________________________________________________________________
UInt_t AliTRDptrgParam::GetMultiplicity(UInt_t BitVector) const {
  // returns the multiplicity of a given bit vector
  
  UInt_t result = 0;
  UInt_t temp = 0x01;
  for (UInt_t iBit = 0; iBit < 32; iBit++) {
    if (BitVector & temp) {
      result++;
    }
    temp <<= 1;
  }

  return result;
}

//______________________________________________________________________________
UInt_t AliTRDptrgParam::GetMultiplicity(Int_t BitVector) const {
  // returns the multiplicity of a given bit vector
  
  UInt_t result = 0;
  UInt_t temp = 0x01;
  for (UInt_t iBit = 0; iBit < 32; iBit++) {
    if (BitVector & temp) {
      result++;
    }
    temp <<= 1;
  }

  return result;
}

//______________________________________________________________________________
//
//      Configuration file parsing (helper functions)
//______________________________________________________________________________

//______________________________________________________________________________
void AliTRDptrgParam::PrepareLine(TString line, TString& identifier, 
                                  TString& value) {
  // Prepares a line for parsing
  // divide identifier and value 

  // clear identifier and value
  identifier.Clear();
  value.Clear();  

  Int_t iLetter = 0;
  while ((line[iLetter] != ' ') && (line[iLetter] != '\t') && 
         (line[iLetter] != '#') && (iLetter < line.Length())) {
    // read identifier
    identifier += line[iLetter];
    iLetter++;
  }
  while (((line[iLetter] == ' ') || (line[iLetter] == '\t')) &&
	 (iLetter < line.Length()) && (line[iLetter] != '#')) {
    // omit whitespaces and tabs in between
    iLetter++;
  }
  while(iLetter < line.Length()) {
    // read value or equation and remove white spaces and tabs
    //if ((line[iLetter] != ' ') && (line[iLetter] != '\t') &&
    //    (line[iLetter] != '#')) {
    if (line[iLetter] != '#') {
      value += line[iLetter];
    }
    iLetter++;
  }
}

//______________________________________________________________________________
TString AliTRDptrgParam::CleanTString(TString string) {
  // Removes white spaces and tabs

  TString result;
  result.Clear();
  for (Int_t iLetter = 0; iLetter < string.Length(); iLetter++) {
    if ((string[iLetter] != ' ') && (string[iLetter] != '\t')) {
      result += string[iLetter];
    }
  } 
  AliDebug(5, Form("Cleaned string: %s", result.Data()));  
  return result;
}

//______________________________________________________________________________
void AliTRDptrgParam::SplitUpValues(TString value, TObjArray& arr) {
  // splits up multiple values into a TObjArray

  TString temp;
  temp.Clear();
  temp.Resize(0);
  for (Int_t iLetter = 0; iLetter < value.Length(); iLetter++) {
    if ((value[iLetter] != ' ') && (value[iLetter] != '\t')) {
      // add another letter
      temp += value[iLetter];
    }
    else {
      // seperator found: white space or tabs
      if (temp.Length()) {
        TObjString* t = new TObjString(temp.Data());
        arr.Add(t);      
        temp.Clear();
        temp.Resize(0);
      }
    }
  }
  if (temp.Length() != 0) {
    TObjString* t = new TObjString(temp.Data());
    arr.Add(t);         
  }
}

//______________________________________________________________________________
UInt_t AliTRDptrgParam::BinaryTStringToInt(TString number) const {
  // converts a binary TString to an integer
  UInt_t temp = 0x01;
  UInt_t result = 0x0;
  for (Int_t i = number.Length() - 1; i >= 0; i--) {
    if (number[i] == '1') {
      result |= temp;
      temp <<= 1;
    }
    else if (number[i] == '0')  {
      temp <<= 1;    
    }
  }
  return result;
}

//______________________________________________________________________________
Bool_t AliTRDptrgParam::ParseMultiplicityCondition(TString condition, 
                                                   UInt_t* threshold,
                                                   UInt_t* mask) {
  // converts a string formed like "M(1111_1111)>4" to a input mask and a
  // multiplicity threshold

  // check whether condition is starting with "M(
  if ((condition[0] != 'M') || ( condition[1] != '(')) {
    return kFALSE;
  }
  
  TString maskStr = "";
  Int_t iLetter = 0;
  
  // extract input mask
  while (condition[iLetter] != ')') {
    maskStr += condition[iLetter++];
  }
  (*mask) = BinaryTStringToInt(maskStr);
  if ((*mask) == 0) { 
    AliDebug(5, Form("Invalid input mask: %s,[%s]", maskStr.Data(), 
                      condition.Data()));
    return kFALSE;
  }
  
  // ensure that ')' is followed by a '>'
  if (condition[++iLetter] != '>') {
    AliDebug(5, Form("multiplicity condition is incorrectly formed: %s", 
                     condition.Data()));
    return kFALSE;
  }
  iLetter++; // move on to the first digit

  TString thresholdStr = "";
  // gain threshold string
  while (((condition[iLetter] != ' ') || (condition[iLetter] != '\t')) &&
         (iLetter < condition.Length())) {
    thresholdStr += condition[iLetter++];
  }
  (*threshold) = thresholdStr.Atoi(); // convert string to integer
  AliDebug(5, Form("mask: 0x%x, multiplicity threshold: %d", (*mask), 
                   (*threshold)));
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliTRDptrgParam::ParseFEB(TString identifier, TString value) {
  // Parse FEB configuration

  //------------------------------------------
  if (identifier.Index("FEB/T0/A/THR") == 0) {
    // FEB T0 thresholds at A side

    TObjArray arr;
    arr.Clear();
    SplitUpValues(value, arr);
    
    if (arr.GetEntries() != 12) {
      AliError("Wrong FEB T0 Threshold count, it must be 12 channels!");
      return kFALSE; 
    }
    for (Int_t iValue = 0; iValue < arr.GetEntries(); iValue++) {
      this->fFEBT0Thresholds[0][iValue] =  
        (dynamic_cast<TObjString*>(arr[iValue])->GetString()).Atoi();
      AliDebug(5, Form("FEB/T0/A/THR[%d]=%d", iValue, 
                       this->fFEBT0Thresholds[0][iValue])); 
    }
    return kTRUE;
  }

  //-----------------------------------------------
  else if (identifier.Index("FEB/T0/C/THR") == 0) {
    // FEB T0 thresholds at c side

    TObjArray arr;
    arr.Clear();
    SplitUpValues(value, arr);
    
    if (arr.GetEntries() != 12) {
      AliError("Wrong FEB T0 Threshold count, it must be 12 channels!");
      return kFALSE; 
    }
    for (Int_t iValue = 0; iValue < arr.GetEntries(); iValue++) {
      this->fFEBT0Thresholds[1][iValue] =  
        (dynamic_cast<TObjString*>(arr[iValue])->GetString()).Atoi();
      AliDebug(5, Form("FEB/T0/C/THR[%d]=%d", iValue, 
                       this->fFEBT0Thresholds[1][iValue])); 
    }
    return kTRUE;
  }
  
  //--------------------------------------------------
  else if ((identifier.Index("FEB/V0/A0/THR") == 0) ||
      (identifier.Index("FEB/V0/A1/THR") == 0) ||
      (identifier.Index("FEB/V0/A2/THR") == 0) ||
      (identifier.Index("FEB/V0/A3/THR") == 0)) {
    // FEB V0 thresholds at a side (cards 0,1,2,3)
   
    TString cardIDstr = identifier(8, 1);
    Int_t cardID = cardIDstr.Atoi();
 
    TObjArray arr;
    arr.Clear();
    SplitUpValues(value, arr);
    
    if (arr.GetEntries() != 8) {
      AliError("Wrong FEB V0 Threshold count, it must be 8 channels!");
      return kFALSE; 
    }
    for (Int_t iValue = 0; iValue < arr.GetEntries(); iValue++) {
      this->fFEBV0Thresholds[0][cardID][iValue] =  
        (dynamic_cast<TObjString*>(arr[iValue])->GetString()).Atoi();
      AliDebug(5, Form("FEB/V0/A%d/THR[%d]=%d", cardID, iValue, 
                       this->fFEBV0Thresholds[0][cardID][iValue])); 
    }
  }

  //--------------------------------------------------
  else if ((identifier.Index("FEB/V0/C0/THR") == 0) ||
      (identifier.Index("FEB/V0/C1/THR") == 0) ||
      (identifier.Index("FEB/V0/C2/THR") == 0) ||
      (identifier.Index("FEB/V0/C3/THR") == 0)) {
    // FEB V0 thresholds at c side (cards 0,1,2,3)
   
    TString cardIDstr = identifier(8, 1);
    Int_t cardID = cardIDstr.Atoi();
 
    TObjArray arr;
    arr.Clear();
    SplitUpValues(value, arr);
    
    if (arr.GetEntries() != 8) {
      AliError("Wrong FEB V0 Threshold count, it must be 8 channels!");
      return kFALSE; 
    }
    for (Int_t iValue = 0; iValue < arr.GetEntries(); iValue++) {
      this->fFEBV0Thresholds[1][cardID][iValue] =  
        (dynamic_cast<TObjString*>(arr[iValue])->GetString()).Atoi();
      AliDebug(5, Form("FEB/V0/C%d/THR[%d]=%d", cardID, iValue, 
                       this->fFEBV0Thresholds[1][cardID][iValue])); 
    }
  }
  
  //-----------------------------------------------
  else if (identifier.Index("FEB/T0/A/LUT") == 0) {
    // FEB T0 look up tables at A side

    TString lutIDstr = identifier(13, 1);
    Int_t lutID = lutIDstr.Atoi();
    
    UInt_t val = 0;
    UInt_t mask = 0;
    ParseMultiplicityCondition(value, &val, &mask);
    this->fFEBT0Multiplicities[0][lutID][0] = val;
    this->fFEBT0Multiplicities[0][lutID][1] = mask;
    AliDebug(5, Form("FEBT0Multiplicities[0/A][%d][val] = %d", lutID, val));
    AliDebug(5, Form("FEBT0Multiplicities[0/A][%d][mask] = %d", lutID, mask));
    
    return kTRUE;
  }

  //-----------------------------------------------
  else if (identifier.Index("FEB/T0/C/LUT") == 0) {
    // FEB T0 look up tables at C side

    TString lutIDstr = identifier(13, 1);
    Int_t lutID = lutIDstr.Atoi();
    
    UInt_t val = 0;
    UInt_t mask = 0;
    ParseMultiplicityCondition(value, &val, &mask);
    this->fFEBT0Multiplicities[1][lutID][0] = val;
    this->fFEBT0Multiplicities[1][lutID][1] = mask;
    AliDebug(5, Form("FEBT0Multiplicities[1/C][%d][val] = %d", lutID, val));
    AliDebug(5, Form("FEBT0Multiplicities[1/C][%d][mask] = %d", lutID, mask));

    return kTRUE;
  }

  //--------------------------------------------------
  else if ((identifier.Index("FEB/V0/A0/LUT") == 0) ||
           (identifier.Index("FEB/V0/A1/LUT") == 0) ||
           (identifier.Index("FEB/V0/A2/LUT") == 0) ||
           (identifier.Index("FEB/V0/A3/LUT") == 0)) {
    // FEB V0 look up tables at A side

    TString cardIDstr = identifier(8, 1);
    Int_t cardID = cardIDstr.Atoi();
 
    TString lutIDstr = identifier(14, 1);
    Int_t lutID = lutIDstr.Atoi();
     
    UInt_t val = 0;
    UInt_t mask = 0;
    ParseMultiplicityCondition(value, &val, &mask);
    this->fFEBV0Multiplicities[0][cardID][lutID][0] = val;
    this->fFEBV0Multiplicities[0][cardID][lutID][1] = mask;
    AliDebug(5, Form("FEBV0Multiplicities[0/A][%d][%d][val] = %d", cardID, 
                     lutID, val));
    AliDebug(5, Form("FEBV0Multiplicities[0/A][%d][%d][mask] = %d", cardID, 
                     lutID, mask));

    return kTRUE;
  }
   
  //--------------------------------------------------
  else if ((identifier.Index("FEB/V0/C0/LUT") == 0) ||
           (identifier.Index("FEB/V0/C1/LUT") == 0) ||
           (identifier.Index("FEB/V0/C2/LUT") == 0) ||
           (identifier.Index("FEB/V0/C3/LUT") == 0)) {
    // FEB V0 look up tables at C side

    TString cardIDstr = identifier(8, 1);
    Int_t cardID = cardIDstr.Atoi();
 
    TString lutIDstr = identifier(14, 1);
    Int_t lutID = lutIDstr.Atoi();
     
    UInt_t val = 0;
    UInt_t mask = 0;
    ParseMultiplicityCondition(value, &val, &mask);
    this->fFEBV0Multiplicities[1][cardID][lutID][0] = val;
    this->fFEBV0Multiplicities[1][cardID][lutID][1] = mask;
    AliDebug(5, Form("FEBV0Multiplicities[1/C][%d][%d][val] = %d", cardID, 
                     lutID, val));
    AliDebug(5, Form("FEBV0Multiplicities[1/C][%d][%d][mask] = %d", cardID, 
                     lutID, mask));

    return kTRUE;
  }
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliTRDptrgParam::ParseCBAC(TString identifier, TString value) {
  // Parse CB-A and CB-C configuration
 
  if (identifier.Index("CBA/LUT/") == 0) {
    // parse CB-A's logical equations
    
    TString eqIDstr = identifier(8, 1);
    Int_t eqID = eqIDstr.Atoi();
   
    if ((eqID == 0) || (eqID == 1)) {
      this->fCBALUTequ[eqID] = this->CleanTString(value);
      AliDebug(5, Form("fCBALUTequ[%d]=%s", eqID, this->fCBALUTequ[eqID].Data())); 
    }
    return kTRUE;    
  }
    
  else if (identifier.Index("CBC/LUT/") == 0) {
    // parse CB-C's logical equations
    
    TString eqIDstr = identifier(8, 1);
    Int_t eqID = eqIDstr.Atoi();
   
    if ((eqID == 0) || (eqID == 1)) {
      this->fCBCLUTequ[eqID] = this->CleanTString(value);
      AliDebug(5, Form("fCBCLUTequ[%d]=%s", eqID, this->fCBCLUTequ[eqID].Data()));
    }
    return kTRUE;
  }

  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliTRDptrgParam::ParseCBB(TString identifier, TString value) {
  // Parse CBB configuration
  
  if (identifier.Index("CBB/LUT/") == 0) {
    // parse CB-B's logical equations
    
    TString eqIDstr = identifier(8, 1);
    Int_t eqID = eqIDstr.Atoi();
   
    if ((eqID == 0) || (eqID == 1) || (eqID == 2)) {
      this->fCBBLUTequ[eqID] = this->CleanTString(value);
      AliDebug(5, Form("fCBBLUTequ[%d]=%s", eqID, this->fCBBLUTequ[eqID].Data()));
    }
    return kTRUE;
  }
  
  // PT masks 
  else if (identifier.Index("CBB/PT/MASK/CB-A_0") == 0) { // CB-A_0
    if (value.Index("YES") == 0) {
      this->fPTmasks.fCBA[0] = kTRUE;     
    }
    else {
      this->fPTmasks.fCBA[0] = kFALSE;
    }
    AliDebug(5, Form("CBB/PT/MASK/CB-A_0=%d", this->fPTmasks.fCBA[0]));
    return kTRUE;
  }
  else if (identifier.Index("CBB/PT/MASK/CB-A_1") == 0) { // CB-A_1
    if (value.Index("YES") == 0) {
      this->fPTmasks.fCBA[1] = kTRUE;     
    }
    else {
      this->fPTmasks.fCBA[1] = kFALSE;
    }
    AliDebug(5, Form("CBB/PT/MASK/CB-A_1=%d", this->fPTmasks.fCBA[1]));
    return kTRUE;
  } 
  else if (identifier.Index("CBB/PT/MASK/CB-C_0") == 0) { // CB-C_0
    if (value.Index("YES") == 0) {
      this->fPTmasks.fCBC[0] = kTRUE;     
    }
    else {
      this->fPTmasks.fCBC[0] = kFALSE;
    }
    AliDebug(5, Form("CBB/PT/MASK/CB-C_0=%d",this->fPTmasks.fCBC[0]));
    return kTRUE;
  }
  else if (identifier.Index("CBB/PT/MASK/CB-C_1") == 0) { // CB-C_1
    if (value.Index("YES") == 0) {
      this->fPTmasks.fCBC[1] = kTRUE;     
    }
    else {
      this->fPTmasks.fCBC[1] = kFALSE;
    }
    AliDebug(5, Form("CBB/PT/MASK/CB-C_1=%d", this->fPTmasks.fCBC[1]));
    return kTRUE;
  } 
  else if (identifier.Index("CBB/PT/MASK/CB-B_0") == 0) { // CB-B_0
    if (value.Index("YES") == 0) {
      this->fPTmasks.fLUTs[0] = kTRUE;     
    }
    else {
      this->fPTmasks.fLUTs[0] = kFALSE;
    }
    AliDebug(5, Form("CBB/PT/MASK/CB-B_0=%d",this->fPTmasks.fLUTs[0]));
    return kTRUE;
  }
  else if (identifier.Index("CBB/PT/MASK/CB-B_1") == 0) { // CB-B_1
    if (value.Index("YES") == 0) {
      this->fPTmasks.fLUTs[1] = kTRUE;     
    }
    else {
      this->fPTmasks.fLUTs[1] = kFALSE;
    }
    AliDebug(5, Form("CBB/PT/MASK/CB-B_1/=%d", this->fPTmasks.fLUTs[1]));
    return kTRUE;
  }
  else if (identifier.Index("CBB/PT/MASK/CB-B_2") == 0) { // CB-B_2
    if (value.Index("YES") == 0) {
      this->fPTmasks.fLUTs[2] = kTRUE;     
    }
    else {
      this->fPTmasks.fLUTs[2] = kFALSE;
    }
    AliDebug(5, Form("CBB/PT/MASK/CB-B_2/=%d", this->fPTmasks.fLUTs[2]));
    return kTRUE;
  }
  else if (identifier.Index("BB/PT/MASK/TLMU_") == 0) {
    TString indexStr = identifier(16, 1);
    Int_t index = indexStr.Atoi();
    if (value.Index("YES") == 0) {
      this->fPTmasks.fTLMU[index] = kTRUE;
    }
    else {
      this->fPTmasks.fTLMU[index] = kFALSE;
    }
    AliDebug(5, Form("CBB/PT/MASK/TLMU_%d=%d", index, 
                     this->fPTmasks.fTLMU[index]));
    return kTRUE;
  }
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliTRDptrgParam::ParseTLMU(TString identifier, TString value) {
  // Parse TLMU configuration

  if (identifier.Index("TLMU/IMASK/SEC") == 0) {
    // TLMU input masks
    TString indexStr = identifier(14,2);
    Int_t index = indexStr.Atoi();
    if ((index < 0) || (index > 17)) {
      AliDebug(5, "Wrong section index in TLMU input mask");
      return kFALSE;
    }
    this->fTLMUInputMask[index] = BinaryTStringToInt(value);
    AliDebug(5, Form("%d %x\n", index, this->fTLMUInputMask[index]));
    return kTRUE;
  }

  //-----------------------------------------------
  else if (identifier.Index("TLMU/CMATRIX") == 0) {
    // TLMU coincidence matrices

    // matrix index
    TString matrixIndexStr = identifier(12,1);
    Int_t matrixIndex = matrixIndexStr.Atoi();
    // entry index
    TString indexStr = identifier(17,2);
    Int_t index = indexStr.Atoi();
    this->fTLMUcmatrices[matrixIndex][index] = BinaryTStringToInt(value);
    AliDebug(5, Form("%d 0x%x\n", matrixIndex, 
                     this->fTLMUcmatrices[matrixIndex][index]));
    return kTRUE;
  }

  //---------------------------------------------
  else if (identifier.Index("TLMU/MCNTR") == 0) {
    // TLMU multiplicity counter setup
    
    TString indexStr = identifier(10,1);
    Int_t index = indexStr.Atoi();
    TObjArray arr;

    SplitUpValues(value, arr);
    
    TString t0 = (dynamic_cast<TObjString*>(arr[0]))->GetString();
    TString t1 = (dynamic_cast<TObjString*>(arr[1]))->GetString();
  
    this->fTLMUmultiplicity[index][0] = t0.Atoi();
    this->fTLMUmultiplicity[index][1] = t1.Atoi();
 
    AliDebug(5, Form("%d: %d  %d", index, this->fTLMUmultiplicity[index][0], 
                     this->fTLMUmultiplicity[index][1]));      

    return kTRUE;
  }
  
  //----------------------------------------------
  else if (identifier.Index("TLMU/OUTMUX") == 0) {
    // TLMU output signal assignment
    TObjArray arr;
    SplitUpValues(value, arr);
  
    if (arr.GetEntries() > 8) {
      AliError("Too many TLMU output signals assigned");
      return kFALSE;
    } 
  
    for (Int_t iEntry = 0; iEntry < arr.GetEntries(); iEntry++) {
      TString t = (dynamic_cast<TObjString*>(arr[iEntry]))->GetString(); 
      
      TString indexStr = t(2,1);
      if (t.Index("CM") == 0) { // coincidence matrix
        this->fTLMUoutput[iEntry][0] = indexStr.Atoi();
      }
      else if (t.Index("MC") == 0) { // multiplicity
        this->fTLMUoutput[iEntry][1] = indexStr.Atoi();
      }
      AliDebug(5, Form("TLMU output: cm = %d, mc = %d", 
                       this->fTLMUoutput[iEntry][0], 
                       this->fTLMUoutput[iEntry][1]));
    }
    return kTRUE;
  }
  return kTRUE;
}

//______________________________________________________________________________
//
//       Logical Equation to LUT processing (helper functions)
//______________________________________________________________________________

//______________________________________________________________________________
Int_t AliTRDptrgParam::LookUp(TString* const identifier) const {
  // Transforms identifier into look up table address bit
  //
  // this function has to be extended/changed when different identifiers for
  // other equations and destination LUTs should be used

  if (identifier->CompareTo("T0_0", TString::kIgnoreCase) == 0)  
    return 0x001;
  else if (identifier->CompareTo("T0_1", TString::kIgnoreCase) == 0) 
    return 0x002; 
  else if (identifier->CompareTo("V0-0_0", TString::kIgnoreCase) == 0) 
    return 0x004; 
  else if (identifier->CompareTo("V0-0_1", TString::kIgnoreCase) == 0) 
    return 0x008; 
  else if (identifier->CompareTo("V0-1_0", TString::kIgnoreCase) == 0) 
    return 0x010; 
  else if (identifier->CompareTo("V0-1_1", TString::kIgnoreCase) == 0) 
    return 0x020; 
  else if (identifier->CompareTo("V0-2_0", TString::kIgnoreCase) == 0) 
    return 0x040; 
  else if (identifier->CompareTo("V0-2_1", TString::kIgnoreCase) == 0) 
    return 0x080; 
  else if (identifier->CompareTo("V0-3_0", TString::kIgnoreCase) == 0) 
    return 0x100; 
  else if (identifier->CompareTo("V0-3_1", TString::kIgnoreCase) == 0) 
    return 0x200; 
  else if (identifier->CompareTo("CB-A_0", TString::kIgnoreCase) == 0) 
    return 0x001; 
  else if (identifier->CompareTo("CB-A_1", TString::kIgnoreCase) == 0) 
    return 0x002; 
  else if (identifier->CompareTo("CB-C_0", TString::kIgnoreCase) == 0) 
    return 0x004; 
  else if (identifier->CompareTo("CB-C_1", TString::kIgnoreCase) == 0) 
    return 0x008; 
  else if (identifier->CompareTo("TLMU_0", TString::kIgnoreCase) == 0)  
    return 0x010; 
  else if (identifier->CompareTo("TLMU_1", TString::kIgnoreCase) == 0) 
    return 0x020; 
  else if (identifier->CompareTo("TLMU_2", TString::kIgnoreCase) == 0) 
    return 0x040; 
  else if (identifier->CompareTo("TLMU_3", TString::kIgnoreCase) == 0) 
    return 0x080; 
  else if (identifier->CompareTo("TLMU_4", TString::kIgnoreCase) == 0) 
    return 0x100; 
  else if (identifier->CompareTo("TLMU_5", TString::kIgnoreCase) == 0) 
    return 0x200; 
  else if (identifier->CompareTo("TLMU_6", TString::kIgnoreCase) == 0) 
    return 0x400; 
  else if (identifier->CompareTo("TLMU_7", TString::kIgnoreCase) == 0) 
    return 0x800; 
  else return 0x0; // Error
}

//______________________________________________________________________________
void AliTRDptrgParam::MergeResults(TArrayI*& partResult1, TArrayI*& partResult2,
                                   TArrayI*& results,
                                   TArrayI*& signalsInvolved1, 
                                   TArrayI*& signalsInvolved2, 
                                   TArrayI*& signalsInvolved, 
                                   Bool_t useOR) {
  // merges result and signal involved arrays
  // uses logical OR (or=kTRUE) and AND (or==kFALSE) as merging function
  
  // check whether input data is valid
  if ((partResult1 == 0x0) || (partResult2 == 0x0) || 
      (signalsInvolved1 == 0x0) || (signalsInvolved2 == 0x0)) {
    AliError("fatal logical equation processing error!");
  }
 
  // allocate results and signalsInvolved 
  results = new TArrayI(0);    
  signalsInvolved = new TArrayI(0);

  // merge arrays (necessary for OR and AND)
  for (Int_t i = 0; i < partResult1->GetSize(); i++) {
    for (Int_t j = 0; j < partResult2->GetSize(); j++) {
      results->Set(results->GetSize() + 1); // increment size
      (*results)[results->GetSize() - 1] =  // add combination
        (*partResult1)[i] | (*partResult2)[j];
 
      signalsInvolved->Set(signalsInvolved->GetSize() + 1);
      (*signalsInvolved)[signalsInvolved->GetSize() - 1] = 
         (*signalsInvolved1)[i] | (*signalsInvolved2)[j];
    }
  }
  
  if (useOR) { // only necessary for OR
    // add partResult1
    for (Int_t i = 0; i < partResult1->GetSize(); i++) {
      results->Set(results->GetSize() + 1);
      (*results)[results->GetSize() - 1] = (*partResult1)[i];
      
      signalsInvolved->Set(signalsInvolved->GetSize() + 1);
      (*signalsInvolved)[signalsInvolved->GetSize()-1] = (*signalsInvolved1)[i];
    }
    // add partResult2
    for (Int_t i = 0; i < partResult2->GetSize(); i++) {
      results->Set(results->GetSize() + 1);
      (*results)[results->GetSize() - 1] = (*partResult2)[i];
      
      signalsInvolved->Set(signalsInvolved->GetSize() + 1);
      (*signalsInvolved)[signalsInvolved->GetSize()-1] = (*signalsInvolved2)[i];
    }
  }
  
  // debug output
  AliDebug(5, "merging results: ");
  for (Int_t i = 0; i < results->GetSize(); i++) {
    AliDebug(5, Form("0x%x 0x%x", (*results)[i], (*signalsInvolved)[i]));
  }

  // free memory
  delete partResult1;
  partResult1 = 0x0;
  delete partResult2;
  partResult2 = 0x0; 
}

//______________________________________________________________________________
void AliTRDptrgParam::ConvertLogicalEqToBitVectors(TString eq, 
                                                   TArrayI*& results,
                                                   TArrayI*& signalsInvolved) {
  // converts a logical equation to a LUT
  //
  // input string must not contain white spaces or tabs
  // only identifiers, ||, &&, (, ) and ! are allowed
  //
  // neglected signals are assumed to be zero in this function
  // this problem is solved by "void CheckSignalsInvolved(...)"

  AliDebug(5, Form("eq: %s", eq.Data()));

  // temporary variables used before/while merging
  TArrayI* partResult1 = 0x0;
  TArrayI* partResult2 = 0x0;
  TArrayI* partResult3 = 0x0;
  TArrayI* partResult4 = 0x0;
  TArrayI* signalsInvolved1 = 0x0;
  TArrayI* signalsInvolved2 = 0x0;
  TArrayI* signalsInvolved3 = 0x0;
  TArrayI* signalsInvolved4 = 0x0;
 
  Int_t iChar = 0; // counter variable
  
  // variables needed for correct operator order (&& before ||!)
  Int_t foundORbefore = -1; // found an || in that string (-1 = not found)
  Int_t foundAND = -1; // found an &&
  Int_t foundORafter = -1; // found a second OR after &&

  // variables needed for correct bracket processing
  Int_t enteredBrackets = 0; // indicates in which bracket layer the parser is
  Int_t bracketLevelAtZero = -1; // when enteredBrackets = 0 was reached first
  // after it ascended  

  while ((iChar < eq.Length())) { //--------------------------------------------
    // walk through string

    // operators ---------------------------------------------------------------
    if ((enteredBrackets == 0 ) && (eq[iChar] != '(') && (eq[iChar] != ')'))  {
      // '|'
      if (eq[iChar] == '|') {
        if (eq[iChar + 1] == '|') { // ||
          iChar++; // jump to the next charakter
          if (foundAND == -1) {
            foundORbefore = iChar;
          }
          else if ((foundORafter == -1) && (foundAND != -1)) {
            foundORafter = iChar;
          }
        }
        else { // bit-wise and not supported
          AliError(Form("LogicalEquation incorrect: %s", eq.Data()));
          AliError("bit-wise AND (&) not supported for now");
          return;
        }
      }
      // '&' 
      else if (eq[iChar] == '&') {
        if (eq[iChar] == '&') { // ||
          iChar++; // jump to the next charakter
          if (foundAND == -1) {
            foundAND = iChar;
          }
        }
        else { // bit-wise or not supported
          AliError(Form("LogicalEquation incorrect: %s", eq.Data()));
          AliError("bit-wise OR (|) not supported for now");
	  return;
        }
      }
    }
    // brackets ----------------------------------------------------------------
    // '(' 
    if (eq[iChar] == '(') {
      enteredBrackets++;      
    }
    // ')' 
    else if (eq[iChar] == ')') {
      enteredBrackets--;
      if (enteredBrackets < 0) {        
        AliError(Form("LogicalEquation incorrect: %s", eq.Data()));     
        AliError("Too many )s");
      }
      if ((enteredBrackets == 0) && (bracketLevelAtZero == -1) &&
          (foundAND == -1) && (foundORbefore == -1)) {
        // needed to detected equations encapsulated in brackets: (...)
        bracketLevelAtZero = iChar;
      }
    }      
    iChar++; // go on to the next letter/char
  } //--------------------------------------------------------------------------

  if (bracketLevelAtZero == (eq.Length() - 1)) { // strip ( ) and process again
    ConvertLogicalEqToBitVectors(eq(1, eq.Length() -2), results, 
                                 signalsInvolved);
    return;     
  }
  else if (foundAND == -1) { // no AND
    if (foundORbefore != -1) { // only OR / || found and no AND
      ConvertLogicalEqToBitVectors(eq(0, foundORbefore-1), partResult1, 
                                   signalsInvolved1); 
      ConvertLogicalEqToBitVectors(eq(foundORbefore+1, 
                                      eq.Length()-foundORbefore-1),
                                   partResult2, signalsInvolved2);
      
      MergeResults(partResult1, partResult2, results, signalsInvolved1,
                   signalsInvolved2, signalsInvolved, kTRUE);
      return;
    } 
    else { // only identifier remained!
      results = new TArrayI(1);
      signalsInvolved = new TArrayI(1);
      if (eq[0] != '!') { // identifier without negation
        (*results)[0] = LookUp(&eq); 
        (*signalsInvolved)[0] = (*results)[0];
      }
      else { // identifier with negation
        (*results)[0] = 0;
        TString eqNegated = eq(1, eq.Length()-1);
        (*signalsInvolved)[0] = LookUp(&eqNegated);
      } 
      return;
    }
  }
  // found single or multiple AND / && 
  else if ((foundORafter != -1) && (foundORbefore != -1)) { 
    // found: ...||...&&...||...   
    ConvertLogicalEqToBitVectors(eq(0, foundORbefore-1), partResult1, 
                                 signalsInvolved1);
    ConvertLogicalEqToBitVectors(eq(foundORbefore+1, 
                                    foundORafter-foundORbefore-2),
                                 partResult2, signalsInvolved2);
    ConvertLogicalEqToBitVectors(eq(foundORafter+1, eq.Length()-foundORafter-1),
                                 partResult3, signalsInvolved3);

    // merge: 4 = 1 || 2 
    MergeResults(partResult1, partResult2, partResult4, signalsInvolved1,
                 signalsInvolved2, signalsInvolved4, kTRUE);
    // merge results = 3 || 4
    MergeResults(partResult3, partResult4, results, signalsInvolved3,
                 signalsInvolved4, signalsInvolved, kTRUE);
    return;
  } 
  else if (foundORbefore != -1) { 
    // found ...||...&&...
    ConvertLogicalEqToBitVectors(eq(0, foundORbefore - 1), partResult1, 
                                 signalsInvolved1); 
    ConvertLogicalEqToBitVectors(eq(foundORbefore+1, 
                                    eq.Length()-foundORbefore-1),
                                 partResult2, signalsInvolved2);
   
    MergeResults(partResult1, partResult2, results, signalsInvolved1,
                 signalsInvolved2, signalsInvolved, kTRUE);
    return;
  }
  else if (foundORafter != -1) {
    // found  ...&&...||...
    ConvertLogicalEqToBitVectors(eq(0, foundORafter - 1), partResult1, 
                                 signalsInvolved1); 
    ConvertLogicalEqToBitVectors(eq(foundORafter+1, 
                                    eq.Length()-foundORafter-1),
                                 partResult2, signalsInvolved2);
   
    MergeResults(partResult1, partResult2, results, signalsInvolved1,
                 signalsInvolved2, signalsInvolved, kTRUE);
    return;
  }
  else /* if (foundAND != -1)*/ { // found ...&&...
    ConvertLogicalEqToBitVectors(eq(0, foundAND-1), partResult1, 
                                 signalsInvolved1); 
    ConvertLogicalEqToBitVectors(eq(foundAND+1, eq.Length()-foundAND-1), 
                                 partResult2, signalsInvolved2); 

    MergeResults(partResult1, partResult2, results, signalsInvolved1,
                 signalsInvolved2, signalsInvolved, kFALSE);    
    return;
  }
  
  AliError("Logical equation parser error!");
  return;
}

//______________________________________________________________________________
void AliTRDptrgParam::CheckSignalsInvolved(TArrayI*& results, 
                                           TArrayI*& signalsInvolved,
                                           Int_t inputWidth) {
  // checks whether all input signals are taken into account
  //
  // this function is needed to be able to write equations which contain not all
  // possible signals and which are not mentioned in the equation do not effect
  // the result
  // X=B&&C=(A||!A)&&B&&C

  // this routine is quite inefficient but working O((2^inputWidth)^3)

  // generate mask:
  Int_t temp = 0x1;
  Int_t mask = 0x0;
  for (Int_t iSignal = 0; iSignal < inputWidth; iSignal++) {
    mask |= temp;
    temp <<= 1; // move temp to the next bit 
  }
  
  for (Int_t iResult = 0; iResult < results->GetSize(); iResult++) {
    // tricky: size of results increases while loop is iterating
    // that is needed to generate all valid input signal combinations
    if (mask != (*signalsInvolved)[iResult]) {
      // not all input signals are taken into account
      Int_t inputSignal = 0x1;
      for (Int_t iSignal = 0; iSignal < inputWidth; iSignal++) {
        if (!(inputSignal & (*signalsInvolved)[iResult])) {
          Int_t newInvolvedSignalCombination = 
            (*signalsInvolved)[iResult] | inputSignal;
          Int_t newResult = inputSignal | (*results)[iResult];
          Bool_t signalCombinationAlreadyEnlisted = kFALSE;
          for (Int_t iEntry = 0; iEntry < signalsInvolved->GetSize(); iEntry++){
            // this loop is needed to reduce the amount of equal entries in 
            // signalsInvolved
            // maybe a table with all possible input values could reduce the
            // computional effort, but this would consume a lot of ram
            if ((signalsInvolved->At(iEntry) == newInvolvedSignalCombination) &&
	        (results->At(iEntry) == newResult)) {
              signalCombinationAlreadyEnlisted = kTRUE;
              break;
	    }
          }
	  if (!signalCombinationAlreadyEnlisted) {
            results->Set(results->GetSize() + 1);
            (*results)[results->GetSize() - 1] = inputSignal | 
                                                 (*results)[iResult];
            // add variant with active bit, variant with inactive signal
            // is already containt in the results array

            // update signalsInvolved:
            signalsInvolved->Set(signalsInvolved->GetSize() + 1);
            (*signalsInvolved)[signalsInvolved->GetSize() - 1] =
              (*signalsInvolved)[iResult] | inputSignal;
	  }
	}
        inputSignal <<= 1; // move temp to the next input signal      
      }      
    }
  }
  return;
}

//______________________________________________________________________________
Int_t* AliTRDptrgParam::GenerateLUTbasedOnEq(TString eq, Int_t inputWidth, 
                                             Int_t initValue) {
  // Processes the conversion of a logical equation to a look up table
  
  TArrayI* results = 0x0;
  TArrayI* signalsInvolved = 0x0;
 
  ConvertLogicalEqToBitVectors(eq, results, signalsInvolved);
  // generate bit vectors

  CheckSignalsInvolved(results, signalsInvolved, inputWidth);
  // add bit vectors for signals which are not taken into account

  Int_t lutSize =  0x1 << inputWidth; // 2^inputwidth elements
  Int_t* resultingLUT = new Int_t[lutSize]; // create LUT
  for (Int_t iLUTentry = 0; iLUTentry < lutSize; iLUTentry++) { // init LUT
    resultingLUT[iLUTentry] = 0;
  }
  for (Int_t iEntry = 0; iEntry < results->GetSize(); iEntry++) {
    resultingLUT[(*results)[iEntry]] = initValue;
  }
  
  if (results != 0x0) {
    delete results;
    results = 0x0;
  }
  if (signalsInvolved != 0x0) {
    delete signalsInvolved;
    signalsInvolved = 0x0;
  }
  
  return resultingLUT;
}

//______________________________________________________________________________
//___ GETTER FUNCTIONS__________________________________________________________
//______________________________________________________________________________
UInt_t* AliTRDptrgParam::GetFEBT0Thresholds(AliTRDptrgFEBPosition_t FEBposition)
  const
{
  // get T0 FEB Thresholds
  return this->fFEBT0Thresholds[FEBposition - 1];
  // 0 kB, 1= kA, 2, kC => -1 because T0FEBs are only in position kA and kC
}

//_____________________________________________________________________________
Int_t* AliTRDptrgParam::GetFEBT0LUT(AliTRDptrgFEBPosition_t FEBposition, 
                                    Int_t iLUT) {
  // get T0 FEB LUTs
  if (this->fFEBT0LUTs == 0x0) {
    this->GenerateLUTs();
  }
  return this->fFEBT0LUTs[FEBposition - 1][iLUT]; 
  // 0 kB, 1= kA, 2, kC => -1 because T0FEBs are only in position kA and kC
} 

//______________________________________________________________________________
UInt_t* AliTRDptrgParam::GetFEBV0Thresholds(AliTRDptrgFEBPosition_t FEBposition,
                                           Int_t iCard) const {
  // get V0 FEB Thresholds
  return this->fFEBV0Thresholds[FEBposition - 1][iCard];
  // 0 kB, 1= kA, 2, kC => -1 because T0FEBs are only in position kA and kC
}

//______________________________________________________________________________
Int_t* AliTRDptrgParam::GetFEBV0LUT(AliTRDptrgFEBPosition_t FEBposition, 
                                    Int_t iCard, Int_t iLUT) {
  // get V0 FEB LUTs
  if (this->fFEBV0LUTs == 0x0) {
    this->GenerateLUTs();
  }
  return this->fFEBV0LUTs[FEBposition - 1][iCard][iLUT];
  // 0 kB, 1= kA, 2, kC => -1 because T0FEBs are only in position kA and kC
}

//______________________________________________________________________________
Int_t* AliTRDptrgParam::GetCBLUT(UInt_t iCB, Int_t LUTid) {
  // return control box LUT
  // iCB: 0 = B, 1 = A, 2 = C
  if (this->fCBLUTs == 0x0) {
    this->GenerateLUTs();
  }
  return this->fCBLUTs[iCB][LUTid];
}

