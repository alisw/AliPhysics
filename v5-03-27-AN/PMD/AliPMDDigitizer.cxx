/***************************************************************************
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
//-----------------------------------------------------//
//                                                     //
//  Source File : PMDDigitizer.cxx, Version 00         //
//                                                     //
//  Date   : September 20 2002                         //
//                                                     //
//-----------------------------------------------------//

#include <Riostream.h>
#include <TBRIK.h>
#include <TNode.h>
#include <TTree.h>
#include <TGeometry.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TParticle.h>
#include <TRandom.h>

#include "AliLog.h"
#include "AliRun.h"
#include "AliHit.h"
#include "AliDetector.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliConfig.h"
#include "AliMagF.h"
#include "AliDigitizationInput.h"
#include "AliDigitizer.h"
#include "AliHeader.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliMC.h"

#include "AliPMD.h"
#include "AliPMDhit.h"
#include "AliPMDcell.h"
#include "AliPMDsdigit.h"
#include "AliPMDdigit.h"
#include "AliPMDCalibData.h"
#include "AliPMDPedestal.h"
#include "AliPMDDigitizer.h"


ClassImp(AliPMDDigitizer)

AliPMDDigitizer::AliPMDDigitizer() :
  fRunLoader(0),
  fPMDHit(0),
  fPMD(0),
  fPMDLoader(0),
  fCalibGain(GetCalibGain()),
  fCalibPed(GetCalibPed()),
  fSDigits(0),
  fDigits(0),
  fCPVCell(0),
  fCell(0),
  fNsdigit(0),
  fNdigit(0),
  fDetNo(0),
  fZPos(361.5)   // in units of cm, default position of PMD
{
  // Default Constructor
  //
  for (Int_t i = 0; i < fgkTotUM; i++)
    {
      for (Int_t j = 0; j < fgkRow; j++)
	{
	  for (Int_t k = 0; k < fgkCol; k++)
	    {
	      fCPV[i][j][k]         = 0.;
	      fPRE[i][j][k]         = 0.;
	      fCPVCounter[i][j][k]  =  0; 
	      fPRECounter[i][j][k]  =  0;
	      fCPVTrackNo[i][j][k]  = -1;
	      fPRETrackNo[i][j][k]  = -1;
	      fCPVTrackPid[i][j][k] = -1;
	      fPRETrackPid[i][j][k] = -1;
	    }
	}
    }

}
//____________________________________________________________________________
AliPMDDigitizer::AliPMDDigitizer(const AliPMDDigitizer& digitizer):
  AliDigitizer(digitizer),
  fRunLoader(0),
  fPMDHit(0),
  fPMD(0),
  fPMDLoader(0),
  fCalibGain(GetCalibGain()),
  fCalibPed(GetCalibPed()),
  fSDigits(0),
  fDigits(0),
  fCPVCell(0),
  fCell(0),
  fNsdigit(0),
  fNdigit(0),
  fDetNo(0),
  fZPos(361.5)   // in units of cm, default position of PMD
{
  // copy constructor
  AliError("Copy constructor not allowed ");
  
}
//____________________________________________________________________________
AliPMDDigitizer & AliPMDDigitizer::operator=(const AliPMDDigitizer& /*digitizer*/)
{
  // Assignment operator
  AliError("Assignement operator not allowed ");

  return *this;
}
//____________________________________________________________________________
AliPMDDigitizer::AliPMDDigitizer(AliDigitizationInput* digInput):
  AliDigitizer(digInput),
  fRunLoader(0),
  fPMDHit(0),
  fPMD(0),
  fPMDLoader(0),
  fCalibGain(GetCalibGain()),
  fCalibPed(GetCalibPed()),
  fSDigits(new TClonesArray("AliPMDsdigit", 1000)),
  fDigits(new TClonesArray("AliPMDdigit", 1000)),
  fCPVCell(0),
  fCell(0),
  fNsdigit(0),
  fNdigit(0),
  fDetNo(0),
  fZPos(361.5)// in units of cm, This is the default position of PMD
{
  // ctor which should be used


  for (Int_t i = 0; i < fgkTotUM; i++)
    {
      for (Int_t j = 0; j < fgkRow; j++)
	{
	  for (Int_t k = 0; k < fgkCol; k++)
	    {
	      fCPV[i][j][k]         = 0.;
	      fPRE[i][j][k]         = 0.;
	      fCPVCounter[i][j][k]  =  0; 
	      fPRECounter[i][j][k]  =  0;
	      fCPVTrackNo[i][j][k]  = -1;
	      fPRETrackNo[i][j][k]  = -1;
	      fCPVTrackPid[i][j][k] = -1;
	      fPRETrackPid[i][j][k] = -1;
	    }
	}
    }
}

//____________________________________________________________________________
AliPMDDigitizer::~AliPMDDigitizer()
{
  // Default Destructor
  //
  if (fSDigits) {
    fSDigits->Delete();
    delete fSDigits;
    fSDigits=0;
  }
  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
    fDigits=0;
  }
  fCPVCell.Delete();
  fCell.Delete();
}
//
// Member functions
//
//____________________________________________________________________________
void AliPMDDigitizer::OpengAliceFile(const char *file, Option_t *option)
{
  // Loads galice.root file and corresponding header, kinematics
  // hits and sdigits or digits depending on the option
  //

  TString evfoldname = AliConfig::GetDefaultEventFolderName();
  fRunLoader = AliRunLoader::GetRunLoader(evfoldname);
  if (!fRunLoader)
      fRunLoader = AliRunLoader::Open(file,AliConfig::GetDefaultEventFolderName(), "UPDATE");
  
  if (!fRunLoader)
   {
     AliError(Form("Can not open session for file %s.",file));
   }

  const char *cHS = strstr(option,"HS");
  const char *cHD = strstr(option,"HD");
  const char *cSD = strstr(option,"SD");
  
  if(cHS || cHD)
    {
      if (!fRunLoader->GetAliRun()) fRunLoader->LoadgAlice();
      if (!fRunLoader->TreeE()) fRunLoader->LoadHeader();
      if (!fRunLoader->TreeK()) fRunLoader->LoadKinematics();
  
      gAlice = fRunLoader->GetAliRun();
  
      if (gAlice)
	{
	  AliDebug(1,"Alirun object found");
	}
      else
	{
	  AliError("Could not found Alirun object");
	}
  
      fPMD  = (AliPMD*)gAlice->GetDetector("PMD");
    }

  fPMDLoader = fRunLoader->GetLoader("PMDLoader");
  if (fPMDLoader == 0x0)
    {
      AliError("Can not find PMDLoader");
    }


  if (cHS)
    {
      fPMDLoader->LoadHits("READ");
      fPMDLoader->LoadSDigits("recreate");
    }
  else if (cHD)
    {
      fPMDLoader->LoadHits("READ");
      fPMDLoader->LoadDigits("recreate");
    }
  else if (cSD)
    {
      fPMDLoader->LoadSDigits("READ");
      fPMDLoader->LoadDigits("recreate");
    }
}
//____________________________________________________________________________
void AliPMDDigitizer::Hits2SDigits(Int_t ievt)
{
  // This reads the PMD Hits tree and assigns the right track number
  // to a cell and stores in the summable digits tree
  //

  const Int_t kPi0 = 111;
  const Int_t kGamma = 22;
  Int_t npmd = 0;
  Int_t trackno = 0;
  Int_t smnumber = 0;
  Int_t trackpid = 0;
  Int_t mtrackno = 0;
  Int_t mtrackpid = 0;

  Float_t xPos = 0., yPos = 0., zPos = 0.;
  Int_t xpad = -1, ypad = -1;
  Float_t edep = 0.;
  Float_t vx = -999.0, vy = -999.0, vz = -999.0;


  if (!fSDigits) fSDigits = new TClonesArray("AliPMDsdigit", 1000);
  ResetSDigit();

  AliDebug(1,Form("Event Number = %d",ievt));
  Int_t nparticles = fRunLoader->GetHeader()->GetNtrack();
  AliDebug(1,Form("Number of Particles = %d",nparticles));

  //

  fRunLoader->GetEvent(ievt);
  // ------------------------------------------------------- //
  // Pointer to specific detector hits.
  // Get pointers to Alice detectors and Hits containers

  TTree* treeH = fPMDLoader->TreeH();
  
  Int_t ntracks    = (Int_t) treeH->GetEntries();
  AliDebug(1,Form("Number of Tracks in the TreeH = %d", ntracks));
  TTree* treeS = fPMDLoader->TreeS();
  if (treeS == 0x0)
    {
      fPMDLoader->MakeTree("S");
      treeS = fPMDLoader->TreeS();
    }
  Int_t bufsize = 16000;
  treeS->Branch("PMDSDigit", &fSDigits, bufsize);
  
  TClonesArray* hits = 0;
  if (fPMD) hits = fPMD->Hits();

  // Start loop on tracks in the hits containers

  for (Int_t track=0; track<ntracks;track++)
    {
      gAlice->GetMCApp()->ResetHits();
      treeH->GetEvent(track);
      if (fPMD)
	{
	  npmd = hits->GetEntriesFast();
	  for (Int_t ipmd = 0; ipmd < npmd; ipmd++)
	    {
	      fPMDHit = (AliPMDhit*) hits->UncheckedAt(ipmd);
	      trackno = fPMDHit->GetTrack();
	      //  get kinematics of the particles

	      TParticle* mparticle = gAlice->GetMCApp()->Particle(trackno);
	      trackpid  = mparticle->GetPdgCode();
	      Int_t  ks = mparticle->GetStatusCode();
	      Int_t imo;
	      Int_t tracknoOld=0, trackpidOld=0, statusOld = 0;
	      
	      if (mparticle->GetFirstMother() == -1)
		{
		  tracknoOld  = trackno;
		  trackpidOld = trackpid;
		  statusOld   = -1;
		}
	      Int_t igstatus = 0;

	      Int_t trnotemp = trackno;    // Modified on 25th Nov 2009
	      if(ks==1||(imo = mparticle->GetFirstMother())<0 ){
		vx = mparticle->Vx();
		vy = mparticle->Vy();
		vz = mparticle->Vz();
		
		if(trackpid==kGamma||trackpid==11||trackpid==-11||
		   trackpid==kPi0)igstatus=1;
	      }
	      
	      
	      while(((imo = mparticle->GetFirstMother()) >= 0)&& 
		    (ks = mparticle->GetStatusCode() <1) )
		{
		  mparticle =  gAlice->GetMCApp()->Particle(imo);
		  trackpid = mparticle->GetPdgCode();
		  ks = mparticle->GetStatusCode();
		  vx = mparticle->Vx();
		  vy = mparticle->Vy();
		  vz = mparticle->Vz();
		  
		  // Modified on 25th Nov 2009

		  trnotemp = trackno;
                  if(trackpid == 111)
		    {
		      trackno = trnotemp;
		    }
                  if(trackpid != 111)
		    {
		      trackno=imo;
		    }
		  // end of modification on 25th Nov 2009
		}
	 
	      if(trackpid==kGamma||trackpid==11||trackpid==-11||
		 trackpid==kPi0)igstatus=1;
	      mtrackpid=trackpid;
	      mtrackno=trackno;
	      trackpid=trackpidOld;
	      trackno=tracknoOld;
	      
	      //-----------------end of modification----------------
	      Float_t ptime = fPMDHit->GetTime()*1e6;   // time in microsec
	      if (ptime < 0. || ptime > 1.2) continue;  

	      xPos = fPMDHit->X();
	      yPos = fPMDHit->Y();
	      zPos = fPMDHit->Z();

	      edep       = fPMDHit->GetEnergy();
	      Int_t vol1 = fPMDHit->GetVolume(1); // Column
	      Int_t vol2 = fPMDHit->GetVolume(2); // Row
	      Int_t vol7 = fPMDHit->GetVolume(4); // Serial Module No


	      // -----------------------------------------//
	      // In new geometry after adding electronics //
	      // For Super Module 1 & 2                   //
	      //  nrow = 48, ncol = 96                    //
	      // For Super Module 3 & 4                   //
	      //  nrow = 96, ncol = 48                    //
	      // -----------------------------------------//

	      if (vol7 < 24)
		{
		  smnumber = vol7;
		}
	      else
		{
		  smnumber = vol7 - 24;
		}
	      Int_t vol8 = smnumber/6 + 1;   // fake supermodule

	      if (vol8 == 1 || vol8 == 2)
		{
		  xpad = vol2;
		  ypad = vol1;
		}
	      else if (vol8 == 3 || vol8 == 4)
		{
		  xpad = vol1;
		  ypad = vol2;
		}

	      AliDebug(2,Form("Zposition = %f Edeposition = %f",zPos,edep));

	      if (vol7 < 24)
		{
		  // PRE
		  fDetNo = 0;
		}
	      else
		{
		  // CPV
		  fDetNo = 1;
		}

	      Int_t smn = smnumber;
	      Int_t ixx = xpad     - 1;
	      Int_t iyy = ypad     - 1;
	      if (fDetNo == 0)
		{
		  fPRE[smn][ixx][iyy] += edep;
		  fPRECounter[smn][ixx][iyy]++;

		  AliPMDcell* cell = new AliPMDcell(mtrackno,smn,ixx,iyy,edep);
		  fCell.Add(cell);
		}
	      else if(fDetNo == 1)
		{
		  fCPV[smn][ixx][iyy] += edep;
		  fCPVCounter[smn][ixx][iyy]++;
		  AliPMDcell* cpvcell = new AliPMDcell(mtrackno,smn,ixx,iyy,edep); 
		  fCPVCell.Add(cpvcell);
		}
	    }
	}
    } // Track Loop ended
  TrackAssignment2CPVCell();
  TrackAssignment2Cell();
  ResetCell();

  Float_t deltaE      = 0.;
  Int_t   detno       = 0;
  Int_t   trno        = -1;
  Int_t   trpid       = -99;

  for (Int_t idet = 0; idet < 2; idet++)
    {
      for (Int_t ism = 0; ism < fgkTotUM; ism++)
	{
	  for (Int_t jrow = 0; jrow < fgkRow; jrow++)
	    {
	      for (Int_t kcol = 0; kcol < fgkCol; kcol++)
		{
		  if (idet == 0)
		    {
		      deltaE = fPRE[ism][jrow][kcol];
		      trno   = fPRETrackNo[ism][jrow][kcol];
		      detno = 0;
		    }
		  else if (idet == 1)
		    {
		      deltaE = fCPV[ism][jrow][kcol];
		      trno   = fCPVTrackNo[ism][jrow][kcol];
		      detno = 1;
		    }
		  if (deltaE > 0.)
		    {
		      // Natasha
		      TParticle *mparticle = gAlice->GetMCApp()->Particle(trno);
		      trpid = mparticle->GetPdgCode();
		      AddSDigit(trno,trpid,detno,ism,jrow,kcol,deltaE);
		    }
		}
	    }
	  treeS->Fill();
	  ResetSDigit();
	}
    }
  fPMDLoader->WriteSDigits("OVERWRITE");
  ResetCellADC();
}
//____________________________________________________________________________

void AliPMDDigitizer::Hits2Digits(Int_t ievt)
{
  // This reads the PMD Hits tree and assigns the right track number
  // to a cell and stores in the digits tree
  //
  const Int_t kPi0 = 111;
  const Int_t kGamma = 22;
  Int_t npmd      = 0;
  Int_t trackno   = 0;
  Int_t smnumber  = 0;
  Int_t trackpid  = 0;
  Int_t mtrackno  = 0;
  Int_t mtrackpid = 0;

  Float_t xPos = 0., yPos = 0., zPos = 0.;
  Int_t xpad = -1, ypad = -1;
  Float_t edep = 0.;
  Float_t vx = -999.0, vy = -999.0, vz = -999.0;

  if (!fDigits) fDigits = new TClonesArray("AliPMDdigit", 1000);
  ResetDigit();

  AliDebug(1,Form("Event Number =  %d",ievt));
  Int_t nparticles = fRunLoader->GetHeader()->GetNtrack();
  AliDebug(1,Form("Number of Particles = %d", nparticles));

  fRunLoader->GetEvent(ievt);
  // ------------------------------------------------------- //
  // Pointer to specific detector hits.
  // Get pointers to Alice detectors and Hits containers

  fPMD  = (AliPMD*)gAlice->GetDetector("PMD");
  fPMDLoader = fRunLoader->GetLoader("PMDLoader");

  if (fPMDLoader == 0x0)
    {
      AliError("Can not find PMD or PMDLoader");
    }
  TTree* treeH = fPMDLoader->TreeH();
  Int_t ntracks    = (Int_t) treeH->GetEntries();
  AliDebug(1,Form("Number of Tracks in the TreeH = %d", ntracks));
  fPMDLoader->LoadDigits("recreate");
  TTree* treeD = fPMDLoader->TreeD();
  if (treeD == 0x0)
    {
      fPMDLoader->MakeTree("D");
      treeD = fPMDLoader->TreeD();
    }
  Int_t bufsize = 16000;
  treeD->Branch("PMDDigit", &fDigits, bufsize);
  
  TClonesArray* hits = 0;
  if (fPMD) hits = fPMD->Hits();

  // Start loop on tracks in the hits containers

  for (Int_t track=0; track<ntracks;track++)
    {
      gAlice->GetMCApp()->ResetHits();
      treeH->GetEvent(track);
      
      if (fPMD)
	{
	  npmd = hits->GetEntriesFast();
	  for (Int_t ipmd = 0; ipmd < npmd; ipmd++)
	    {
	      fPMDHit = (AliPMDhit*) hits->UncheckedAt(ipmd);
	      trackno = fPMDHit->GetTrack();
	      
	      //  get kinematics of the particles
	      
	      TParticle* mparticle = gAlice->GetMCApp()->Particle(trackno);
	      trackpid  = mparticle->GetPdgCode();
	      Int_t  ks = mparticle->GetStatusCode();
	      Int_t imo;
	      Int_t tracknoOld=0, trackpidOld=0, statusOld = 0;
	      if (mparticle->GetFirstMother() == -1)
		{
		  tracknoOld  = trackno;
		  trackpidOld = trackpid;
		  statusOld   = -1;
		}

	      Int_t igstatus = 0;

	      Int_t trnotemp = trackno;   // modified on 25th Nov 2009
	      if(ks==1||(imo = mparticle->GetFirstMother())<0 ){
		vx = mparticle->Vx();
		vy = mparticle->Vy();
		vz = mparticle->Vz();
		
		if(trackpid==kGamma||trackpid==11||trackpid==-11||trackpid==kPi0)
		  igstatus=1;
	      }
	      
	      
	      while(((imo = mparticle->GetFirstMother()) >= 0)&& 
		    (ks = mparticle->GetStatusCode() <1) )
		{
		  mparticle =  gAlice->GetMCApp()->Particle(imo);
		  trackpid = mparticle->GetPdgCode();
		  ks = mparticle->GetStatusCode();
		  vx = mparticle->Vx();
		  vy = mparticle->Vy();
		  vz = mparticle->Vz();
		  
		  // Modified on 25th Nov 2009
		  
		  trnotemp = trackno;
                  if(trackpid == 111)
		    {
		      trackno = trnotemp;
		    }
                  if(trackpid != 111)
		    {
		      trackno=imo;
		    }
		}
	 
	      if(trackpid==kGamma||trackpid==11||trackpid==-11||trackpid==kPi0)
		igstatus=1;
	      mtrackpid=trackpid;
	      mtrackno=trackno;
	      trackpid=trackpidOld;
	      trackno=tracknoOld;

	      Float_t ptime = fPMDHit->GetTime()*1e6;
	      if (ptime < 0. || ptime > 1.2) continue;
	      
	      xPos = fPMDHit->X();
	      yPos = fPMDHit->Y();
	      zPos = fPMDHit->Z();
	      edep       = fPMDHit->GetEnergy();
	      Int_t vol1 = fPMDHit->GetVolume(1); // Column
	      Int_t vol2 = fPMDHit->GetVolume(2); // Row
	      Int_t vol7 = fPMDHit->GetVolume(4); // Serial Module No

	      // -----------------------------------------//
	      // In new geometry after adding electronics //
	      // For Super Module 1 & 2                   //
	      //  nrow = 48, ncol = 96                    //
	      // For Super Module 3 & 4                   //
	      //  nrow = 96, ncol = 48                    //
	      // -----------------------------------------//
	      
	      if (vol7 < 24)
		{
		  smnumber = vol7;
		}
	      else
		{
		  smnumber = vol7 - 24;
		}
	      Int_t vol8 = smnumber/6 + 1;    // fake supermodule

	      if (vol8 == 1 || vol8 == 2)
		{
		  xpad = vol2;
		  ypad = vol1;
		}
	      else if (vol8 == 3 || vol8 == 4)
		{
		  xpad = vol1;
		  ypad = vol2;
		}

	      AliDebug(2,Form("ZPosition = %f Edeposition = %f",zPos,edep));

	      if (vol7 < 24)
		{
		  // PRE
		  fDetNo = 0;
		}
	      else
		{
		  fDetNo = 1;
		}

	      Int_t smn = smnumber;
	      Int_t ixx = xpad     - 1;
	      Int_t iyy = ypad     - 1;
	      if (fDetNo == 0)
		{
		  fPRE[smn][ixx][iyy] += edep;
		  fPRECounter[smn][ixx][iyy]++;

		  AliPMDcell* cell = new AliPMDcell(mtrackno,smn,ixx,iyy,edep);

		  fCell.Add(cell);
		}
	      else if(fDetNo == 1)
		{
		  fCPV[smn][ixx][iyy] += edep;
		  fCPVCounter[smn][ixx][iyy]++;
		  AliPMDcell* cpvcell = new AliPMDcell(mtrackno,smn,ixx,iyy,edep); 
		  fCPVCell.Add(cpvcell);
		}
	    }
	}
    } // Track Loop ended
  TrackAssignment2CPVCell();
  TrackAssignment2Cell();
  ResetCell();

  Float_t gain1  = 1.;
  Float_t adc    = 0. ;
  Float_t deltaE = 0.;
  Int_t detno    = 0;
  Int_t trno     = 1;
  Int_t trpid    = -99;

  for (Int_t idet = 0; idet < 2; idet++)
  {
      for (Int_t ism = 0; ism < fgkTotUM; ism++)
      {
	  for (Int_t jrow = 0; jrow < fgkRow; jrow++)
	  {
	      for (Int_t kcol = 0; kcol < fgkCol; kcol++)
	      {
		  if (idet == 0)
		  {
		      deltaE = fPRE[ism][jrow][kcol];
		      trno   = fPRETrackNo[ism][jrow][kcol];
		      detno  = 0;
		  }
		  else if (idet == 1)
		  {
		      deltaE = fCPV[ism][jrow][kcol];
		      trno   = fCPVTrackNo[ism][jrow][kcol];
		      detno  = 1;
		  }
		  if (deltaE > 0.)
		  {
		      MeV2ADC(deltaE,adc);

		      // To decalibrate the adc values
		      //
		      gain1 = Gain(idet,ism,jrow,kcol);
		      if (gain1 != 0.)
		      {
			  Int_t adcDecalib = (Int_t)(adc/gain1);
			  adc = (Float_t) adcDecalib;
		      }
		      else if(gain1 == 0.)
		      {
			  adc = 0.;
		      }

		      // Pedestal Decalibration
		      Int_t   pedmeanrms = 
			  fCalibPed->GetPedMeanRms(idet,ism,jrow,kcol);
		      Int_t   pedrms1    = (Int_t) pedmeanrms%100;
		      Float_t pedrms     = (Float_t)pedrms1/10.;
		      Float_t pedmean    = 
			  (Float_t) (pedmeanrms - pedrms1)/1000.0;
		      if (adc > 0.)
		      {
			  adc += (pedmean + 3.0*pedrms);
			  TParticle *mparticle
			    = gAlice->GetMCApp()->Particle(trno);
			  trpid = mparticle->GetPdgCode();
			  
			  AddDigit(trno,trpid,detno,ism,jrow,kcol,adc);
		      }
		  }
	      } // column loop
	  } // row    loop
	  treeD->Fill();
	  ResetDigit();
      } // supermodule loop
  } // detector loop
  
  fPMDLoader->WriteDigits("OVERWRITE");
  ResetCellADC();

}
//____________________________________________________________________________


void AliPMDDigitizer::SDigits2Digits(Int_t ievt)
{
  // This reads the PMD sdigits tree and converts energy deposition
  // in a cell to ADC and stores in the digits tree
  //

  fRunLoader->GetEvent(ievt);

  TTree* treeS = fPMDLoader->TreeS();
  AliPMDsdigit  *pmdsdigit;
  TBranch *branch = treeS->GetBranch("PMDSDigit");
  if(!branch)
    {
      AliError("PMD Sdigit branch does not exist");
      return;
    }
  if (!fSDigits) fSDigits = new TClonesArray("AliPMDsdigit", 1000);
  branch->SetAddress(&fSDigits);

  TTree* treeD = fPMDLoader->TreeD();
  if (treeD == 0x0)
    {
      fPMDLoader->MakeTree("D");
      treeD = fPMDLoader->TreeD();
    }
  Int_t bufsize = 16000;
  if (!fDigits) fDigits = new TClonesArray("AliPMDdigit", 1000);
  treeD->Branch("PMDDigit", &fDigits, bufsize);

  Int_t   trno = 1, trpid = 0, det = 0, smn = 0;
  Int_t   irow = 0, icol = 0;
  Float_t edep = 0., adc = 0.;

  Int_t nmodules = (Int_t) treeS->GetEntries();
  AliDebug(1,Form("Number of modules = %d",nmodules));

  for (Int_t imodule = 0; imodule < nmodules; imodule++)
    {
      treeS->GetEntry(imodule);
      Int_t nentries = fSDigits->GetLast();
      AliDebug(2,Form("Number of entries per module = %d",nentries+1));
      for (Int_t ient = 0; ient < nentries+1; ient++)
	{
	  pmdsdigit = (AliPMDsdigit*)fSDigits->UncheckedAt(ient);
	  trno   = pmdsdigit->GetTrackNumber();
	  trpid  = pmdsdigit->GetTrackPid();
	  det    = pmdsdigit->GetDetector();
	  smn    = pmdsdigit->GetSMNumber();
	  irow   = pmdsdigit->GetRow();
	  icol   = pmdsdigit->GetColumn();
	  edep   = pmdsdigit->GetCellEdep();

	  MeV2ADC(edep,adc);

	  // To decalibrte the adc values
	  //
	  Float_t gain1 = Gain(det,smn,irow,icol);
	  if (gain1 != 0.)
	  {
	    Int_t adcDecalib = (Int_t)(adc/gain1);
	    adc = (Float_t) adcDecalib;
	  }
	  else if(gain1 == 0.)
	  {
	      adc = 0.;
	  }
	  // Pedestal Decalibration
	  Int_t   pedmeanrms = fCalibPed->GetPedMeanRms(det,smn,irow,icol);
	  Int_t   pedrms1    = (Int_t) pedmeanrms%100;
	  Float_t pedrms     = (Float_t)pedrms1/10.;
	  Float_t pedmean    = (Float_t) (pedmeanrms - pedrms1)/1000.0;
	  if(adc > 0.)
	  {
	      adc += (pedmean + 3.0*pedrms);
	      AddDigit(trno,trpid,det,smn,irow,icol,adc);
	  }

	}
      treeD->Fill();
      ResetDigit();
    }
  fPMDLoader->WriteDigits("OVERWRITE");

}
//____________________________________________________________________________
void AliPMDDigitizer::Digitize(Option_t *option)
{
  // Does the event merging and digitization
  const char *cdeb = strstr(option,"deb");
  if(cdeb)
    {
      AliDebug(100," *** PMD Exec is called ***");
    }

  Int_t ninputs = fDigInput->GetNinputs();
  AliDebug(1,Form("Number of files to be processed = %d",ninputs));
  ResetCellADC();

  for (Int_t i = 0; i < ninputs; i++)
    {
      Int_t troffset = fDigInput->GetMask(i);
      MergeSDigits(i, troffset);
    }

  fRunLoader = AliRunLoader::GetRunLoader(fDigInput->GetOutputFolderName());
  fPMD  = (AliPMD*)gAlice->GetDetector("PMD");
  fPMDLoader = fRunLoader->GetLoader("PMDLoader");
  if (fPMDLoader == 0x0)
    {
      AliError("Can not find PMD or PMDLoader");
    }
  fPMDLoader->LoadDigits("update");
  TTree* treeD = fPMDLoader->TreeD();
  if (treeD == 0x0)
    {
      fPMDLoader->MakeTree("D");
      treeD = fPMDLoader->TreeD();
    }
  Int_t bufsize = 16000;
  if (!fDigits) fDigits = new TClonesArray("AliPMDdigit", 1000);
  treeD->Branch("PMDDigit", &fDigits, bufsize);

  Float_t adc    = 0.;
  Float_t deltaE = 0.;
  Int_t detno    = 0;
  Int_t trno     = 1;
  Int_t trpid    = -99;

  for (Int_t idet = 0; idet < 2; idet++)
    {
      for (Int_t ism = 0; ism < fgkTotUM; ism++)
	{
	  for (Int_t jrow = 0; jrow < fgkRow; jrow++)
	    {
	      for (Int_t kcol = 0; kcol < fgkCol; kcol++)
		{
		  if (idet == 0)
		    {
		      deltaE = fPRE[ism][jrow][kcol];
		      trno   = fPRETrackNo[ism][jrow][kcol];
		      trpid  = fPRETrackPid[ism][jrow][kcol];
		      detno  = 0;
		    }
		  else if (idet == 1)
		    {
		      deltaE = fCPV[ism][jrow][kcol];
		      trno   = fCPVTrackNo[ism][jrow][kcol];
		      trpid  = fCPVTrackPid[ism][jrow][kcol];
		      detno  = 1;
		    }
		  if (deltaE > 0.)
		    {
		      MeV2ADC(deltaE,adc);

                      //
		      // Gain decalibration
		      //
		      Float_t gain1 = Gain(idet,ism,jrow,kcol);

		      if (gain1 != 0.)
		      {
			  Int_t adcDecalib = (Int_t)(adc/gain1);
			  adc = (Float_t) adcDecalib;
		      }
		      else if(gain1 == 0.)
		      {
			  adc = 0.;
		      }
		      // Pedestal Decalibration
		      Int_t   pedmeanrms = 
			  fCalibPed->GetPedMeanRms(idet,ism,jrow,kcol);
		      Int_t   pedrms1    = (Int_t) pedmeanrms%100;
		      Float_t pedrms     = (Float_t)pedrms1/10.;
		      Float_t pedmean    = 
			  (Float_t) (pedmeanrms - pedrms1)/1000.0;
		      if (adc > 0.)
		      {
			  adc += (pedmean + 3.0*pedrms);
			  AddDigit(trno,trpid,detno,ism,jrow,kcol,adc);
		      }

		    }
		} // column loop
	    } // row    loop
	  treeD->Fill();
	  ResetDigit();
	} // supermodule loop
    } // detector loop
  fPMDLoader->WriteDigits("OVERWRITE");
  fPMDLoader->UnloadDigits();
  ResetCellADC();
}
//____________________________________________________________________________
void AliPMDDigitizer::TrackAssignment2CPVCell()
{
  // This block assigns the cell id when there are
  // multiple tracks in a cell according to the
  // energy deposition
  // This method added by Ajay
  Bool_t jsort = false;

  Int_t i = 0, j = 0, k = 0;

  Int_t   *status1;
  Int_t   *status2;
  Int_t   *trnarray;  
  Float_t *fracEdp;
  Float_t *trEdp;
  
  Int_t   ****cpvTrack;
  Float_t ****cpvEdep;

  cpvTrack = new Int_t ***[fgkTotUM];
  cpvEdep  = new Float_t ***[fgkTotUM];
  for (i=0; i<fgkTotUM; i++)
    {
      cpvTrack[i] = new Int_t **[fgkRow];
      cpvEdep[i]  = new Float_t **[fgkRow];
    }

  for (i = 0; i < fgkTotUM; i++)
    {
      for (j = 0; j < fgkRow; j++)
	{
	  cpvTrack[i][j] = new Int_t *[fgkCol];
	  cpvEdep[i][j]  = new Float_t *[fgkCol];
	}
    }
  for (i = 0; i < fgkTotUM; i++)
    {
      for (j = 0; j < fgkRow; j++)
	{
	  for (k = 0; k < fgkCol; k++)
	    {
	      Int_t nn = fCPVCounter[i][j][k];
	      if(nn > 0)
		{
		  cpvTrack[i][j][k] = new Int_t[nn];
		  cpvEdep[i][j][k] = new Float_t[nn];
		}
	      else
		{
		  nn = 1;
		  cpvTrack[i][j][k] = new Int_t[nn];
		  cpvEdep[i][j][k] = new Float_t[nn];
		}		      
	      fCPVCounter[i][j][k] = 0;
	    }
	}
    }


  Int_t nentries = fCPVCell.GetEntries();
 
  Int_t   mtrackno = 0, ism = 0, ixp = 0, iyp = 0;
  Float_t edep = 0.;
  for (i = 0; i < nentries; i++)
    {
      AliPMDcell* cpvcell = (AliPMDcell*)fCPVCell.UncheckedAt(i);
      
      mtrackno = cpvcell->GetTrackNumber();
      ism      = cpvcell->GetSMNumber();
      ixp      = cpvcell->GetX();
      iyp      = cpvcell->GetY();
      edep     = cpvcell->GetEdep();
      Int_t nn = fCPVCounter[ism][ixp][iyp];
      cpvTrack[ism][ixp][iyp][nn] = (Int_t) mtrackno;
      cpvEdep[ism][ixp][iyp][nn] = edep;
      fCPVCounter[ism][ixp][iyp]++;
    }
  
  Int_t iz = 0, il = 0;
  Int_t im = 0, ix = 0, iy = 0;
  Int_t nn = 0;
  for (im=0; im<fgkTotUM; im++)
    {
      for (ix=0; ix<fgkRow; ix++)
	{
	  for (iy=0; iy<fgkCol; iy++)
	    {
	      nn = fCPVCounter[im][ix][iy];
	      if (nn > 1)
		{
		  // This block handles if a cell is fired
		  // many times by many tracks
		  status1  = new Int_t[nn];
		  status2  = new Int_t[2*nn];
		  trnarray = new Int_t[nn];
		  for (iz = 0; iz < nn; iz++)
		    {
		      status1[iz] = cpvTrack[im][ix][iy][iz];
		    }
		  TMath::Sort(nn,status1,status2,jsort);
		  Int_t trackOld = -99999;
		  Int_t track, trCount = 0;
		  for (iz = 0; iz < nn; iz++)
		    {
		      track = status1[status2[iz]];
		      if (trackOld != track)
			{
			  trnarray[trCount] = track;
			  trCount++;
			}			      
		      trackOld = track;
		    }
		  delete [] status1;
		  delete [] status2;
		  Float_t totEdp = 0.;
		  trEdp = new Float_t[trCount];
		  fracEdp = new Float_t[trCount];
		  for (il = 0; il < trCount; il++)
		    {
		      trEdp[il] = 0.;
		      track = trnarray[il];
		      for (iz = 0; iz < nn; iz++)
			{
			  if (track == cpvTrack[im][ix][iy][iz])
			    {
			      trEdp[il] += cpvEdep[im][ix][iy][iz];
			    }
			}
		      totEdp += trEdp[il];
		    }
		  Int_t ilOld = 0;
		  Float_t fracOld = 0.;
		  
		  for (il = 0; il < trCount; il++)
		    {
		      fracEdp[il] = trEdp[il]/totEdp;
		      if (fracOld < fracEdp[il])
			{
			  fracOld = fracEdp[il];
			  ilOld = il;
			}
		    }
		  fCPVTrackNo[im][ix][iy] = trnarray[ilOld];
		  delete [] fracEdp;
		  delete [] trEdp;
		  delete [] trnarray;
		}
	      else if (nn == 1)
		{
		  // This only handles if a cell is fired
		  // by only one track
		  
		  fCPVTrackNo[im][ix][iy] = cpvTrack[im][ix][iy][0];
		  
		}
	      else if (nn ==0)
		{
		  // This is if no cell is fired
		  fCPVTrackNo[im][ix][iy] = -999;
		}
	    } // end of iy
	} // end of ix
    } // end of im
  
  // Delete all the pointers
  
 for (i = 0; i < fgkTotUM; i++)
    {
      for (j = 0; j < fgkRow; j++)
	{
	  for (k = 0; k < fgkCol; k++)
	    {
	      delete []cpvTrack[i][j][k];
	      delete []cpvEdep[i][j][k];
	    }
	}
    }
 
  for (i = 0; i < fgkTotUM; i++)
    {
      for (j = 0; j < fgkRow; j++)
	{
	  delete [] cpvTrack[i][j];
	  delete [] cpvEdep[i][j];
	}
    }
  
  for (i = 0; i < fgkTotUM; i++)
    {
      delete [] cpvTrack[i];
      delete [] cpvEdep[i];
    }
  delete [] cpvTrack;
  delete [] cpvEdep;
  
  // 
  // End of the cell id assignment
  //
}
//____________________________________________________________________________

void AliPMDDigitizer::MergeSDigits(Int_t filenumber, Int_t troffset)
{
  // merging sdigits
  fRunLoader = AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(filenumber));
  fPMDLoader = fRunLoader->GetLoader("PMDLoader");
  fPMDLoader->LoadSDigits("read");
  TTree* treeS = fPMDLoader->TreeS();
  AliPMDsdigit  *pmdsdigit;
  TBranch *branch = treeS->GetBranch("PMDSDigit");
  if (!fSDigits) fSDigits = new TClonesArray("AliPMDsdigit", 1000);
  branch->SetAddress(&fSDigits);

  Int_t   itrackno = 1, itrackpid = 0, idet = 0, ism = 0;
  Int_t   ixp = 0, iyp = 0;
  Float_t edep = 0.;
  Int_t nmodules = (Int_t) treeS->GetEntries();
  AliDebug(1,Form("Number of Modules in the treeS = %d",nmodules));
  AliDebug(1,Form("Track Offset = %d",troffset));
  for (Int_t imodule = 0; imodule < nmodules; imodule++)
    {
      treeS->GetEntry(imodule);
      Int_t nentries = fSDigits->GetLast();
      AliDebug(2,Form("Number of Entries per Module = %d",nentries));
      for (Int_t ient = 0; ient < nentries+1; ient++)
	{
	  pmdsdigit = (AliPMDsdigit*)fSDigits->UncheckedAt(ient);
	  itrackno  = pmdsdigit->GetTrackNumber();
	  itrackpid = pmdsdigit->GetTrackPid();
	  idet      = pmdsdigit->GetDetector();
	  ism       = pmdsdigit->GetSMNumber();
	  ixp       = pmdsdigit->GetRow();
	  iyp       = pmdsdigit->GetColumn();
	  edep      = pmdsdigit->GetCellEdep();
	  if (idet == 0)
	    {
	      if (fPRE[ism][ixp][iyp] < edep)
		{
		  fPRETrackNo[ism][ixp][iyp] = troffset + itrackno;
		  fPRETrackPid[ism][ixp][iyp] = itrackpid;
		}
	      fPRE[ism][ixp][iyp] += edep;
	    }
	  else if (idet == 1)
	    {
	      if (fCPV[ism][ixp][iyp] < edep)
		{
		  fCPVTrackNo[ism][ixp][iyp] = troffset + itrackno;
		  fCPVTrackPid[ism][ixp][iyp] = itrackpid;
		}
	      fCPV[ism][ixp][iyp] += edep;
	    }
	}
    }

}
// ----------------------------------------------------------------------
void AliPMDDigitizer::TrackAssignment2Cell()
{
  // 
  // This block assigns the cell id when there are
  // multiple tracks in a cell according to the
  // energy deposition
  //
  Bool_t jsort = false;

  Int_t i = 0, j = 0, k = 0;

  Int_t   *status1;
  Int_t   *status2;
  Int_t   *trnarray;
  Float_t *fracEdp;
  Float_t *trEdp;
  
  Int_t   ****pmdTrack;
  Float_t ****pmdEdep;

  pmdTrack = new Int_t ***[fgkTotUM];
  pmdEdep  = new Float_t ***[fgkTotUM];
  for (i=0; i<fgkTotUM; i++)
    {
      pmdTrack[i] = new Int_t **[fgkRow];
      pmdEdep[i]  = new Float_t **[fgkRow];
    }

  for (i = 0; i < fgkTotUM; i++)
    {
      for (j = 0; j < fgkRow; j++)
	{
	  pmdTrack[i][j] = new Int_t *[fgkCol];
	  pmdEdep[i][j]  = new Float_t *[fgkCol];
	}
    }
  
  for (i = 0; i < fgkTotUM; i++)
    {
      for (j = 0; j < fgkRow; j++)
	{
	  for (k = 0; k < fgkCol; k++)
	    {
	      Int_t nn = fPRECounter[i][j][k];
	      if(nn > 0)
		{
		  pmdTrack[i][j][k] = new Int_t[nn];
		  pmdEdep[i][j][k] = new Float_t[nn];
		}
	      else
		{
		  nn = 1;
		  pmdTrack[i][j][k] = new Int_t[nn];
		  pmdEdep[i][j][k] = new Float_t[nn];
		}
	      fPRECounter[i][j][k] = 0;
	    }
	}
    }


  Int_t nentries = fCell.GetEntries();

  Int_t   mtrackno, ism, ixp, iyp;
  Float_t edep;

  for (i = 0; i < nentries; i++)
    {
      AliPMDcell* cell = (AliPMDcell*)fCell.UncheckedAt(i);
      
      mtrackno = cell->GetTrackNumber();
      ism      = cell->GetSMNumber();
      ixp      = cell->GetX();
      iyp      = cell->GetY();
      edep     = cell->GetEdep();
      Int_t nn = fPRECounter[ism][ixp][iyp];
      pmdTrack[ism][ixp][iyp][nn] = (Int_t) mtrackno;
      pmdEdep[ism][ixp][iyp][nn] = edep;
      fPRECounter[ism][ixp][iyp]++;
    }
  
  Int_t iz = 0, il = 0;
  Int_t im = 0, ix = 0, iy = 0;
  Int_t nn = 0;
  
  for (im=0; im<fgkTotUM; im++)
    {
      for (ix=0; ix<fgkRow; ix++)
	{
	  for (iy=0; iy<fgkCol; iy++)
	    {
	      nn = fPRECounter[im][ix][iy];
	      if (nn > 1)
		{
		  // This block handles if a cell is fired
		  // many times by many tracks
		  status1  = new Int_t[nn];
		  status2  = new Int_t[2*nn];
		  trnarray = new Int_t[nn];
		  for (iz = 0; iz < nn; iz++)
		    {
		      status1[iz] = pmdTrack[im][ix][iy][iz];
		    }
		  TMath::Sort(nn,status1,status2,jsort);
		  Int_t trackOld = -99999;
		  Int_t track, trCount = 0;
		  for (iz = 0; iz < nn; iz++)
		    {
		      track = status1[status2[iz]];
		      if (trackOld != track)
			{
			  trnarray[trCount] = track;
			  trCount++;
			}
		      trackOld = track;
		    }
		  delete [] status1;
		  delete [] status2;
		  Float_t totEdp = 0.;
		  trEdp = new Float_t[trCount];
		  fracEdp = new Float_t[trCount];
		  for (il = 0; il < trCount; il++)
		    {
		      trEdp[il] = 0.;
		      track = trnarray[il];
		      for (iz = 0; iz < nn; iz++)
			{
			  if (track == pmdTrack[im][ix][iy][iz])
			    {
			      trEdp[il] += pmdEdep[im][ix][iy][iz];
			    }
			}
		      totEdp += trEdp[il];
		    }
		  Int_t ilOld = 0;
		  Float_t fracOld = 0.;
		  
		  for (il = 0; il < trCount; il++)
		    {
		      fracEdp[il] = trEdp[il]/totEdp;
		      if (fracOld < fracEdp[il])
			{
			  fracOld = fracEdp[il];
			  ilOld = il;
			}
		    }
		  fPRETrackNo[im][ix][iy] = trnarray[ilOld];
		  delete [] fracEdp;
		  delete [] trEdp;
		  delete [] trnarray;
		}
	      else if (nn == 1)
		{
		  // This only handles if a cell is fired
		  // by only one track
		  
		  fPRETrackNo[im][ix][iy] = pmdTrack[im][ix][iy][0];
		  
		}
	      else if (nn ==0)
		{
		  // This is if no cell is fired
		  fPRETrackNo[im][ix][iy] = -999;
		}
	    } // end of iy
	} // end of ix
    } // end of im
  
  // Delete all the pointers
  
  for (i = 0; i < fgkTotUM; i++)
    {
      for (j = 0; j < fgkRow; j++)
	{
	  for (k = 0; k < fgkCol; k++)
	    {
	      delete [] pmdTrack[i][j][k];
	      delete [] pmdEdep[i][j][k];
	    }
	}
    }
  
  for (i = 0; i < fgkTotUM; i++)
    {
      for (j = 0; j < fgkRow; j++)
	{
	  delete [] pmdTrack[i][j];
	  delete [] pmdEdep[i][j];
	}
    }
  
  for (i = 0; i < fgkTotUM; i++)
    {
      delete [] pmdTrack[i];
      delete [] pmdEdep[i];
    }
  delete [] pmdTrack;
  delete [] pmdEdep;
  // 
  // End of the cell id assignment
  //
}
//____________________________________________________________________________
void AliPMDDigitizer::MeV2ADC(Float_t mev, Float_t & adc) const
{
  // This converts the simulated edep to ADC according to the
  // Test Beam Data
  // PS Test in June 2010, Voltage @ 1300 V
  // KeV - ADC conversion for 12bit ADC
  // MPV data used for the fit and taken here

  // constants are from Test Beam 2010
  
  const Float_t kConstant   = 0.612796;
  const Float_t kSlope      = 130.158;
  
  Float_t adc12bit = kSlope*mev*0.001 + kConstant;
  if (adc12bit < 0.) adc12bit = 0.;

  //Introducing Readout Resolution for ALICE-PMD

  Float_t sigrr     = 0.605016 - 0.000273*adc12bit + 6.54e-8*adc12bit*adc12bit;
  Float_t adcwithrr = gRandom->Gaus(adc12bit,sigrr);

  if(adcwithrr < 0.)
    {
      adc = 0.;
    }
  else if(adcwithrr >= 0. && adcwithrr < 1600.0)
    {
      adc = adcwithrr;
    }
  else if (adcwithrr >= 1600.0)
    {
      adc = 1600.0;
    }

}
//____________________________________________________________________________
void AliPMDDigitizer::AddSDigit(Int_t trnumber, Int_t trpid, Int_t det,
				Int_t smnumber, Int_t irow, Int_t icol,
				Float_t adc)
{
  // Add SDigit
  //
  if (!fSDigits) fSDigits = new TClonesArray("AliPMDsdigit", 1000);
  TClonesArray &lsdigits = *fSDigits;
  new(lsdigits[fNsdigit++])  AliPMDsdigit(trnumber,trpid,det,smnumber,irow,icol,adc);
}
//____________________________________________________________________________

void AliPMDDigitizer::AddDigit(Int_t trnumber, Int_t trpid, Int_t det,
			       Int_t smnumber, Int_t irow, Int_t icol,
			       Float_t adc)
{
  // Add Digit
  //
  if (!fDigits) fDigits = new TClonesArray("AliPMDdigit", 1000);
  TClonesArray &ldigits = *fDigits;
  new(ldigits[fNdigit++]) AliPMDdigit(trnumber,trpid, det,smnumber,irow,icol,adc);
}
//____________________________________________________________________________

void AliPMDDigitizer::SetZPosition(Float_t zpos)
{
  fZPos = zpos;
}
//____________________________________________________________________________
Float_t AliPMDDigitizer::GetZPosition() const
{
  return fZPos;
}
//____________________________________________________________________________

void AliPMDDigitizer::ResetCell()
{
  // clears the cell array and also the counter
  //  for each cell
  //
  fCPVCell.Delete();
  fCell.Delete();
  for (Int_t i = 0; i < fgkTotUM; i++)
    {
      for (Int_t j = 0; j < fgkRow; j++)
	{
	  for (Int_t k = 0; k < fgkCol; k++)
	    {
	      fCPVCounter[i][j][k] = 0; 
	      fPRECounter[i][j][k] = 0;
	    }
	}
    }
}
//____________________________________________________________________________
void AliPMDDigitizer::ResetSDigit()
{
  // Clears SDigits
  fNsdigit = 0;
  if (fSDigits) fSDigits->Delete();
}
//____________________________________________________________________________
void AliPMDDigitizer::ResetDigit()
{
  // Clears Digits
  fNdigit = 0;
  if (fDigits) fDigits->Delete();
}
//____________________________________________________________________________

void AliPMDDigitizer::ResetCellADC()
{
  // Clears individual cells edep and track number
  for (Int_t i = 0; i < fgkTotUM; i++)
    {
      for (Int_t j = 0; j < fgkRow; j++)
	{
	  for (Int_t k = 0; k < fgkCol; k++)
	    {
	      fCPV[i][j][k]         = 0.;
	      fPRE[i][j][k]         = 0.;
	      fCPVTrackNo[i][j][k]  = 0;
	      fPRETrackNo[i][j][k]  = 0;
	      fCPVTrackPid[i][j][k] = -1;
	      fPRETrackPid[i][j][k] = -1;
	    }
	}
    }
}
//____________________________________________________________________________

void AliPMDDigitizer::UnLoad(Option_t *option)
{
  // Unloads all the root files
  //
  const char *cS = strstr(option,"S");
  const char *cD = strstr(option,"D");

  fRunLoader->UnloadgAlice();
  fRunLoader->UnloadHeader();
  fRunLoader->UnloadKinematics();

  if (cS)
    {
      fPMDLoader->UnloadHits();
    }
  if (cD)
    {
      fPMDLoader->UnloadHits();
      fPMDLoader->UnloadSDigits();
    }
}

//----------------------------------------------------------------------
Float_t AliPMDDigitizer::Gain(Int_t det, Int_t smn, Int_t row, Int_t col) const
{
  // returns of the gain of the cell
  // Added this method by ZA

  //cout<<" I am here in gain "<<fCalibData<< "smn,row, col "<<smn
  //<<" "<<row<<" "<<col<<endl;

  if(!fCalibGain) {
    AliError("No calibration data loaded from CDB!!!");
    return 1;
  }

  Float_t GainFact;
  GainFact = fCalibGain->GetGainFact(det,smn,row,col);
  return GainFact;
}
//----------------------------------------------------------------------
AliPMDCalibData* AliPMDDigitizer::GetCalibGain() const
{
  // The run number will be centralized in AliCDBManager,
  // you don't need to set it here!
  // Added this method by ZA
  // Cleaned up by Alberto
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("PMD/Calib/Gain");
  
  if(!entry) AliFatal("Calibration object retrieval failed!");
  
  AliPMDCalibData *calibdata=0;
  if (entry) calibdata = (AliPMDCalibData*) entry->GetObject();
  
  if (!calibdata)  AliFatal("No calibration data from calibration database !");
  
  return calibdata;
}
//----------------------------------------------------------------------
AliPMDPedestal* AliPMDDigitizer::GetCalibPed() const
{
  // The run number will be centralized in AliCDBManager,
  // you don't need to set it here!

  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("PMD/Calib/Ped");
  
  if(!entry) AliFatal("Pedestal object retrieval failed!");
  
  AliPMDPedestal *pedestal=0;
  if (entry) pedestal = (AliPMDPedestal*) entry->GetObject();
  
  if (!pedestal)  AliFatal("No pedestal data from calibration database !");
  
  return pedestal;
}
