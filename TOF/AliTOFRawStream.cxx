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

/*
Revision 0.02  2005/07/28 A. De Caro:
        Update format TOF raw data
               (temporary solution) 
        Correction of few wrong corrispondences
               between 'software' and 'hardware' numberings

Revision 0.01  2005/07/22 A. De Caro
        Implement methods Next()
	                  GetSector(),
	                  GetPlate(),
	                  GetStrip(),
	                  GetPadZ(),
	                  GetPadX()
*/

///////////////////////////////////////////////////////////////////////////////
//
// This class provides access to TOF raw data in DDL files.
//
// It loops over all TOF raw data given by the AliRawReader.
//
///////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliRawReader.h"

#include "AliTOFGeometry.h"
#include "AliTOFRawStream.h"

ClassImp(AliTOFRawStream)


//_____________________________________________________________________________
AliTOFRawStream::AliTOFRawStream(AliRawReader* rawReader)
{
  //
  // create an object to read TOF raw digits
  //

  fRawReader = rawReader;
  fDDL = -1;
  fTRM = -1;
  fTDC = -1;
  fTDCChannel = -1;
  fTof = -1;
  fADC = -1;
  fErrorFlag = -1;
  //fCounter = -1; // v0.01

  fTOFGeometry = new AliTOFGeometry();

  fRawReader->Select(5);

}

//_____________________________________________________________________________
AliTOFRawStream::AliTOFRawStream(const AliTOFRawStream& stream) :
  TObject(stream)
{
  //
  // copy constructor
  //

  fRawReader = NULL;
  fDDL = -1;
  fTRM = -1;
  fTDC = -1;
  fTDCChannel = -1;
  fTof = -1;
  fADC = -1;
  fErrorFlag = -1;
  //fCounter = -1; // v0.01

  fTOFGeometry = new AliTOFGeometry();

}

//_____________________________________________________________________________
AliTOFRawStream& AliTOFRawStream::operator = (const AliTOFRawStream& stream)
{
  //
  // assignment operator
  //

  fRawReader = stream.fRawReader;
  fDDL = stream.fDDL;
  fTRM = stream.fTRM;
  fTDC = stream.fTDC;
  fTDCChannel = stream.fTDCChannel;
  fTof = stream.fTof;
  fADC = stream.fADC;
  fErrorFlag = stream.fErrorFlag;
  //fCounter = stream.fCounter; // v0.01

  fTOFGeometry = stream.fTOFGeometry;

  return *this;

}

//_____________________________________________________________________________
AliTOFRawStream::~AliTOFRawStream()
{
// destructor

  fTOFGeometry = 0;

}


//_____________________________________________________________________________
Bool_t AliTOFRawStream::Next()
{
// read the next raw digit
// returns kFALSE if there is no digit left

  //fCounter++; // v0.01

  UInt_t data;

  /*
  // v0.01
  if (fCounter==0) {
    if (!fRawReader->ReadNextInt(data)) return kFALSE;
    Int_t dummy = data & 0x007F; // first  7 bits
    AliInfo(Form("This is the number of the current DDL file %2i",dummy));
  }
  */

  if (!fRawReader->ReadNextInt(data)) return kFALSE;

  fTRM         =  data        & 0x000F; // first  4 bits
  fTDC         = (data >>  4) & 0x001F; // next   5 bits
  fTDCChannel  = (data >>  9) & 0x0007; // next   3 bits
  fADC         = (data >> 12) & 0x000FFFFF; // last  20 bits // v0.01
  //fADC         = (data >> 12) & 0x00FF; // next   8 bits // v0.02

  if (!fRawReader->ReadNextInt(data)) return kFALSE;

  fErrorFlag   =  data        & 0x00FF; // first 8 bits
  fTof         = (data >>  8) & 0x00FFFFFF; // last 24 bits // v0.01
  //fTof         = (data >>  8) & 0x0FFF; // next 12 bits // v0.02

  fDDL  = fRawReader->GetDDLID();

  //AliInfo(Form("fDDL  %2i,fTRM %2i, fTDC %2i, fTDCChannel %1i",fDDL,fTRM,fTDC,fTDCChannel));

  return kTRUE;
}
//_____________________________________________________________________________

Int_t AliTOFRawStream::GetSector() const
{
  //
  // Transform the Hardware coordinate to Software coordinate
  // and also write it in the digit form
  //

  Int_t iSector = -1;
  iSector = (Int_t)((Float_t)fDDL/AliTOFGeometry::NDDL());

  return iSector;

}
//_____________________________________________________________________________

Int_t AliTOFRawStream::GetPlate() const
{
  //
  // Transform the Hardware coordinate to Software coordinate
  // and also write it in the digit form
  //

  Int_t iPlate = -1;

  Int_t localDDL = fDDL%AliTOFGeometry::NDDL();
  
  Int_t iPadPerSector = fTDCChannel + AliTOFGeometry::NCh()*(fTDC + AliTOFGeometry::NTdc()*(fTRM + AliTOFGeometry::NTRM()*localDDL));

  Int_t iStripPerSector = (Int_t)(iPadPerSector/(Float_t)AliTOFGeometry::NpadX()/(Float_t)AliTOFGeometry::NpadZ());

  if (iStripPerSector     < fTOFGeometry->NStripC())
    iPlate = 0;
  else if (iStripPerSector>=fTOFGeometry->NStripC() &&
	   iStripPerSector< fTOFGeometry->NStripC()+AliTOFGeometry::NStripB())
    iPlate = 1;
  else if (iStripPerSector>=fTOFGeometry->NStripC()+AliTOFGeometry::NStripB() &&
	   iStripPerSector< fTOFGeometry->NStripC()+AliTOFGeometry::NStripB()+AliTOFGeometry::NStripA())
    iPlate = 2;
  else if (iStripPerSector>=fTOFGeometry->NStripC()+AliTOFGeometry::NStripB()+AliTOFGeometry::NStripA() &&
	   iStripPerSector< fTOFGeometry->NStripC()+AliTOFGeometry::NStripB()+AliTOFGeometry::NStripA()+AliTOFGeometry::NStripB())
    iPlate = 3;
  else if (iStripPerSector>=fTOFGeometry->NStripC()+AliTOFGeometry::NStripB()+AliTOFGeometry::NStripA()+AliTOFGeometry::NStripB() &&
	   iStripPerSector< fTOFGeometry->NStripC()+AliTOFGeometry::NStripB()+AliTOFGeometry::NStripA()+AliTOFGeometry::NStripB()+fTOFGeometry->NStripC())
    iPlate = 4;
  else
    iPlate = -1;

  return iPlate;

}
//_____________________________________________________________________________

Int_t AliTOFRawStream::GetStrip() const
{
  //
  // Transform the Hardware coordinate to Software coordinate
  // and also write it in the digit form
  //

  Int_t iStrip = -1;

  Int_t localDDL = fDDL%AliTOFGeometry::NDDL();

  Int_t iPadPerSector = fTDCChannel + AliTOFGeometry::NCh()*(fTDC + AliTOFGeometry::NTdc()*(fTRM + AliTOFGeometry::NTRM()*localDDL));

  Int_t iStripPerSector = (Int_t)(iPadPerSector/(Float_t)AliTOFGeometry::NpadX()/(Float_t)AliTOFGeometry::NpadZ());

  if (iStripPerSector      < fTOFGeometry->NStripC())
    iStrip = iStripPerSector;
  else if (iStripPerSector >=fTOFGeometry->NStripC() &&
	   iStripPerSector < fTOFGeometry->NStripC()+AliTOFGeometry::NStripB())
    iStrip = iStripPerSector-fTOFGeometry->NStripC();
  else if (iStripPerSector >=fTOFGeometry->NStripC()+AliTOFGeometry::NStripB() &&
	   iStripPerSector < fTOFGeometry->NStripC()+AliTOFGeometry::NStripB()+AliTOFGeometry::NStripA()) 
    iStrip = iStripPerSector-fTOFGeometry->NStripC()-AliTOFGeometry::NStripB();
  else if (iStripPerSector >=fTOFGeometry->NStripC()+AliTOFGeometry::NStripB()+AliTOFGeometry::NStripA() &&
	   iStripPerSector < fTOFGeometry->NStripC()+AliTOFGeometry::NStripB()+AliTOFGeometry::NStripA()+AliTOFGeometry::NStripB())
    iStrip = iStripPerSector-fTOFGeometry->NStripC()-AliTOFGeometry::NStripB()-AliTOFGeometry::NStripA();
  else if (iStripPerSector >=fTOFGeometry->NStripC()+AliTOFGeometry::NStripB()+AliTOFGeometry::NStripA()+AliTOFGeometry::NStripB() &&
	   iStripPerSector < fTOFGeometry->NStripC()+AliTOFGeometry::NStripB()+AliTOFGeometry::NStripA()+AliTOFGeometry::NStripB()+fTOFGeometry->NStripC())
    iStrip = iStripPerSector-fTOFGeometry->NStripC()-AliTOFGeometry::NStripB()-AliTOFGeometry::NStripA()-AliTOFGeometry::NStripB();
  else
    iStrip = -1;

  return iStrip;

}
//_____________________________________________________________________________

Int_t AliTOFRawStream::GetPadZ() const
{
  //
  // Transform the Hardware coordinate to Software coordinate
  // and also write it in the digit form
  //

  Int_t iPadZ = -1;

  Int_t localDDL = fDDL%AliTOFGeometry::NDDL();

  Int_t iPadPerSector = fTDCChannel + AliTOFGeometry::NCh()*(fTDC + AliTOFGeometry::NTdc()*(fTRM + AliTOFGeometry::NTRM()*localDDL));
  
  Int_t iStripPerSector = (Int_t)(iPadPerSector/(Float_t)AliTOFGeometry::NpadX()/(Float_t)AliTOFGeometry::NpadZ());

  Int_t iPadPerStrip = iPadPerSector-AliTOFGeometry::NpadX()*AliTOFGeometry::NpadZ()*iStripPerSector;

  iPadZ = iPadPerStrip%AliTOFGeometry::NpadZ();

  return iPadZ;

}
//_____________________________________________________________________________

Int_t AliTOFRawStream::GetPadX() const
{
  //
  // Transform the Hardware coordinate to Software coordinate
  // and also write it in the digit form
  //

  Int_t iPadX = -1;

  Int_t localDDL = fDDL%AliTOFGeometry::NDDL();

  Int_t iPadPerSector = fTDCChannel + AliTOFGeometry::NCh()*(fTDC + AliTOFGeometry::NTdc()*(fTRM + AliTOFGeometry::NTRM()*localDDL));

  Int_t iStripPerSector = (Int_t)(iPadPerSector/(Float_t)AliTOFGeometry::NpadX()/(Float_t)AliTOFGeometry::NpadZ());

  Int_t iPadPerStrip = iPadPerSector-AliTOFGeometry::NpadX()*AliTOFGeometry::NpadZ()*iStripPerSector;

  iPadX = (Int_t)(iPadPerStrip/(Float_t)AliTOFGeometry::NpadZ());

  return iPadX;

}
//_____________________________________________________________________________
