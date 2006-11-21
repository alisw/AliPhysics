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
#include "AliRunDigitizer.h"
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
#include "AliPMDDigitizer.h"


ClassImp(AliPMDDigitizer)

AliPMDDigitizer::AliPMDDigitizer() :
  fRunLoader(0),
  fPMDHit(0),
  fPMD(0),
  fPMDLoader(0),
  fCalibData(GetCalibData()),
  fSDigits(0),
  fDigits(0),
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
	      fCPV[i][j][k] = 0.; 
	      fPRE[i][j][k] = 0.; 
	      fPRECounter[i][j][k] =  0; 
	      fPRETrackNo[i][j][k] = -1; 
	      fCPVTrackNo[i][j][k] = -1; 
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
  fCalibData(GetCalibData()),
  fSDigits(0),
  fDigits(0),
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
AliPMDDigitizer::AliPMDDigitizer(AliRunDigitizer* manager):
  AliDigitizer(manager),
  fRunLoader(0),
  fPMDHit(0),
  fPMD(0),
  fPMDLoader(0),
  fCalibData(GetCalibData()),
  fSDigits(new TClonesArray("AliPMDsdigit", 1000)),
  fDigits(new TClonesArray("AliPMDdigit", 1000)),
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
	      fCPV[i][j][k] = 0.; 
	      fPRE[i][j][k] = 0.; 
	      fPRECounter[i][j][k] =  0; 
	      fPRETrackNo[i][j][k] = -1; 
	      fCPVTrackNo[i][j][k] = -1; 
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
    fRunLoader = AliRunLoader::Open(file,AliConfig::GetDefaultEventFolderName(),
				    "UPDATE");
  
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
  Int_t npmd;
  Int_t trackno;
  Int_t smnumber;
  Int_t trackpid;
  Int_t mtrackno;
  Int_t mtrackpid;

  Float_t xPos, yPos, zPos;
  Int_t xpad = -1, ypad = -1;
  Float_t edep;
  Float_t vx = -999.0, vy = -999.0, vz = -999.0;


  if (!fSDigits) fSDigits = new TClonesArray("AliPMDsdigit", 1000);
  ResetSDigit();

  AliDebug(1,Form("Event Number = %d",ievt));
  Int_t nparticles = fRunLoader->GetHeader()->GetNtrack();
  AliDebug(1,Form("Number of Particles = %d",nparticles));
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
      gAlice->ResetHits();
      treeH->GetEvent(track);
      if (fPMD) 
	{
	  npmd = hits->GetEntriesFast();
	  for (int ipmd = 0; ipmd < npmd; ipmd++) 
	    {
	      fPMDHit = (AliPMDhit*) hits->UncheckedAt(ipmd);
	      trackno = fPMDHit->GetTrack();
	      //  get kinematics of the particles

	      TParticle* mparticle = gAlice->GetMCApp()->Particle(trackno);
	      trackpid  = mparticle->GetPdgCode();

	      Int_t igatr = -999;
	      Int_t ichtr = -999;
	      Int_t igapid = -999;
	      Int_t imo;
	      Int_t igen = 0;
	      Int_t idmo = -999;

	      Int_t tracknoOld=0, trackpidOld=0, statusOld = 0;
	      if (mparticle->GetFirstMother() == -1)
		{
		  tracknoOld  = trackno;
		  trackpidOld = trackpid;
		  statusOld   = -1;
		}
	      Int_t igstatus = 0;
	      while((imo = mparticle->GetFirstMother()) >= 0)
		{
		  igen++;

		  mparticle =  gAlice->GetMCApp()->Particle(imo);
		  idmo = mparticle->GetPdgCode();
		  
		  vx = mparticle->Vx();
		  vy = mparticle->Vy();
		  vz = mparticle->Vz();
		
		  //printf("==> Mother ID %5d %5d %5d Vertex: %13.3f %13.3f %13.3f\n", igen, imo, idmo, vx, vy, vz);
		  //fprintf(fpw1,"==> Mother ID %5d %5d %5d Vertex: %13.3f %13.3f %13.3f\n", igen, imo, idmo, vx, vy, vz);
		  if ((idmo == kGamma || idmo == -11 || idmo == 11) && vx == 0. && vy == 0. && vz == 0.)
		    {
		      igatr = imo;
		      igapid = idmo;
		      igstatus = 1;
		    }
		  if(igstatus == 0)
		    {
		      if (idmo == kPi0 && vx == 0. && vy == 0. && vz == 0.)
			{
			  igatr = imo;
			  igapid = idmo;
			}
		    }
		  ichtr = imo;
		}

	      if (idmo == kPi0 && vx == 0. && vy == 0. && vz == 0.)
		{
		  mtrackno = igatr;
		  mtrackpid = igapid;
		}
	      else
		{
		  mtrackno  = ichtr;
		  mtrackpid = idmo;
		}
	      if (statusOld == -1)
		{
		  mtrackno  = tracknoOld;
		  mtrackpid = trackpidOld;
		}
	      xPos = fPMDHit->X();
	      yPos = fPMDHit->Y();
	      zPos = fPMDHit->Z();
	      
	      edep       = fPMDHit->GetEnergy();
	      Int_t vol1 = fPMDHit->GetVolume(1); // Column
	      Int_t vol2 = fPMDHit->GetVolume(2); // Row
	      Int_t vol7 = fPMDHit->GetVolume(7); // UnitModule
	      Int_t vol8 = fPMDHit->GetVolume(8); // SuperModule


	      // -----------------------------------------//
	      // In new geometry after adding electronics //
	      // For Super Module 1 & 2                   //
	      //  nrow = 48, ncol = 96                    //
	      // For Super Module 3 & 4                   //
	      //  nrow = 96, ncol = 48                    //
	      // -----------------------------------------//


	      
	      smnumber = (vol8-1)*6 + vol7;

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
	      //Float_t zposition = TMath::Abs(zPos);
	      if (zPos < fZPos)
		{
		  // CPV
		  fDetNo = 1;
		}
	      else if (zPos > fZPos)
		{
		  // PMD
		  fDetNo = 0;
		}
	      Int_t smn = smnumber - 1;
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
		  fCPVTrackNo[smn][ixx][iyy] = mtrackno;
		}
	    }
	}
    } // Track Loop ended

  TrackAssignment2Cell();
  ResetCell();

  Float_t deltaE      = 0.;
  Int_t   detno       = 0;
  Int_t   trno        = -1;

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
		      AddSDigit(trno,detno,ism,jrow,kcol,deltaE);
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
  Int_t npmd;
  Int_t trackno;
  Int_t smnumber;
  Int_t trackpid;
  Int_t mtrackno;
  Int_t mtrackpid;

  Float_t xPos, yPos, zPos;
  Int_t xpad = -1, ypad = -1;
  Float_t edep;
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
      gAlice->ResetHits();
      treeH->GetEvent(track);
      
      if (fPMD) 
	{
	  npmd = hits->GetEntriesFast();
	  for (int ipmd = 0; ipmd < npmd; ipmd++) 
	    {
	      fPMDHit = (AliPMDhit*) hits->UncheckedAt(ipmd);
	      trackno = fPMDHit->GetTrack();
	      
	      //  get kinematics of the particles
	      
	      TParticle* mparticle = gAlice->GetMCApp()->Particle(trackno);
	      trackpid  = mparticle->GetPdgCode();

	      Int_t igatr = -999;
	      Int_t ichtr = -999;
	      Int_t igapid = -999;
	      Int_t imo;
	      Int_t igen = 0;
	      Int_t idmo = -999;

	      Int_t tracknoOld=0, trackpidOld=0, statusOld = 0;
	      if (mparticle->GetFirstMother() == -1)
		{
		  tracknoOld  = trackno;
		  trackpidOld = trackpid;
		  statusOld   = -1;
		}

	      Int_t igstatus = 0;
	      while((imo = mparticle->GetFirstMother()) >= 0)
		{
		  igen++;

		  mparticle =  gAlice->GetMCApp()->Particle(imo);
		  idmo = mparticle->GetPdgCode();
		  
		  vx = mparticle->Vx();
		  vy = mparticle->Vy();
		  vz = mparticle->Vz();
		
		  //printf("==> Mother ID %5d %5d %5d Vertex: %13.3f %13.3f %13.3f\n", igen, imo, idmo, vx, vy, vz);
		  //fprintf(fpw1,"==> Mother ID %5d %5d %5d Vertex: %13.3f %13.3f %13.3f\n", igen, imo, idmo, vx, vy, vz);
		  if ((idmo == kGamma || idmo == -11 || idmo == 11) && vx == 0. && vy == 0. && vz == 0.)
		    {
		      igatr = imo;
		      igapid = idmo;
		      igstatus = 1;
		    }
		  if(igstatus == 0)
		    {
		      if (idmo == kPi0 && vx == 0. && vy == 0. && vz == 0.)
			{
			  igatr = imo;
			  igapid = idmo;
			}
		    }
		  ichtr = imo;
		}

	      if (idmo == kPi0 && vx == 0. && vy == 0. && vz == 0.)
		{
		  mtrackno = igatr;
		  mtrackpid = igapid;
		}
	      else
		{
		  mtrackno  = ichtr;
		  mtrackpid = idmo;
		}
	      if (statusOld == -1)
		{
		  mtrackno  = tracknoOld;
		  mtrackpid = trackpidOld;
		}
	      
	      xPos = fPMDHit->X();
	      yPos = fPMDHit->Y();
	      zPos = fPMDHit->Z();
	      edep       = fPMDHit->GetEnergy();
	      Int_t vol1 = fPMDHit->GetVolume(1); // Column
	      Int_t vol2 = fPMDHit->GetVolume(2); // Row
	      Int_t vol7 = fPMDHit->GetVolume(7); // UnitModule
	      Int_t vol8 = fPMDHit->GetVolume(8); // SuperModule


	      // -----------------------------------------//
	      // In new geometry after adding electronics //
	      // For Super Module 1 & 2                   //
	      //  nrow = 48, ncol = 96                    //
	      // For Super Module 3 & 4                   //
	      //  nrow = 96, ncol = 48                    //
	      // -----------------------------------------//
	      
	      smnumber = (vol8-1)*6 + vol7;

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

	      AliDebug(2,Form("ZPosition = %f Edeposition = %d",zPos,edep));
	      //Float_t zposition = TMath::Abs(zPos);

	      if (zPos < fZPos)
		{
		  // CPV
		  fDetNo = 1;
		}
	      else if (zPos > fZPos)
		{
		  // PMD
		  fDetNo = 0;
		}

	      Int_t smn = smnumber - 1;
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
		  fCPVTrackNo[smn][ixx][iyy] = mtrackno;
		}
	    }
	}
    } // Track Loop ended

  TrackAssignment2Cell();
  ResetCell();

  Float_t gain1;
  Float_t adc;
  Float_t deltaE = 0.;
  Int_t detno = 0;
  Int_t trno = 1;
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
		      gain1 = Gain(idet,ism,jrow,kcol);

		      deltaE = fPRE[ism][jrow][kcol]*gain1;
		      trno   = fPRETrackNo[ism][jrow][kcol];
		      detno = 0;
		    }
		  else if (idet == 1)
		    {
		      gain1 = Gain(idet,ism,jrow,kcol);
		      deltaE = fCPV[ism][jrow][kcol]*gain1;
		      trno   = fCPVTrackNo[ism][jrow][kcol];
		      detno = 1;
		    }
		  if (deltaE > 0.)
		    {
		      MeV2ADC(deltaE,adc);
		      AddDigit(trno,detno,ism,jrow,kcol,adc);
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

  Int_t   trno, det, smn;
  Int_t   irow, icol;
  Float_t edep, adc;

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
	  det    = pmdsdigit->GetDetector();
	  smn    = pmdsdigit->GetSMNumber();
	  irow   = pmdsdigit->GetRow();
	  icol   = pmdsdigit->GetColumn();
	  edep   = pmdsdigit->GetCellEdep();

	  MeV2ADC(edep,adc);
	  AddDigit(trno,det,smn,irow,icol,adc);      
	}
      treeD->Fill();
      ResetDigit();
    }
  fPMDLoader->WriteDigits("OVERWRITE");

}
//____________________________________________________________________________
void AliPMDDigitizer::Exec(Option_t *option)
{
  // Does the event merging and digitization
  const char *cdeb = strstr(option,"deb");
  if(cdeb)
    {
      AliDebug(100," *** PMD Exec is called ***");
    }

  Int_t ninputs = fManager->GetNinputs();
  AliDebug(1,Form("Number of files to be processed = %d",ninputs));
  ResetCellADC();

  for (Int_t i = 0; i < ninputs; i++)
    {
      Int_t troffset = fManager->GetMask(i);
      MergeSDigits(i, troffset);
    }

  fRunLoader = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
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

  Float_t adc;
  Float_t deltaE = 0.;
  Int_t detno = 0;
  Int_t trno = 1;

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
		      MeV2ADC(deltaE,adc);
		      AddDigit(trno,detno,ism,jrow,kcol,adc);
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

void AliPMDDigitizer::MergeSDigits(Int_t filenumber, Int_t troffset)
{
  // merging sdigits
  fRunLoader = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(filenumber));
  fPMDLoader = fRunLoader->GetLoader("PMDLoader");
  fPMDLoader->LoadSDigits("read");
  TTree* treeS = fPMDLoader->TreeS();
  AliPMDsdigit  *pmdsdigit;
  TBranch *branch = treeS->GetBranch("PMDSDigit");
  if (!fSDigits) fSDigits = new TClonesArray("AliPMDsdigit", 1000);
  branch->SetAddress(&fSDigits);

  Int_t   itrackno, idet, ism;
  Int_t   ixp, iyp;
  Float_t edep;
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
		}
	      fPRE[ism][ixp][iyp] += edep;
	    }
	  else if (idet == 1)
	    {
	      if (fCPV[ism][ixp][iyp] < edep)
		{
		  fCPVTrackNo[ism][ixp][iyp] = troffset + itrackno;
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

  Int_t i, j, k;

  Float_t *fracEdp;
  Float_t *trEdp;
  Int_t *status1;
  Int_t *status2;
  Int_t *trnarray;
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
  
  Int_t iz, il;
  Int_t im, ix, iy;
  Int_t nn;
  
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
		  status2  = new Int_t[nn];
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
  // To be done
  //

  // PS Test in September 2003
  // MeV - ADC conversion for 10bit ADC

  const Float_t kConstant   = 7.181;
  const Float_t kErConstant = 0.6899;
  const Float_t kSlope      = 35.93;
  const Float_t kErSlope    = 0.306;

  //gRandom->SetSeed();

  Float_t cons   = gRandom->Gaus(kConstant,kErConstant);
  Float_t slop   = gRandom->Gaus(kSlope,kErSlope);

  Float_t adc10bit = slop*mev*0.001 + cons;

  // 12 bit adc

  Int_t adc12bit = (Int_t) (4.0*adc10bit);

  if(adc12bit < 3000)
    {
      adc = (Float_t) adc12bit;
    }
  else if (adc12bit >= 3000)
    {
      adc = 3000.0;
    }

}
//____________________________________________________________________________
void AliPMDDigitizer::AddSDigit(Int_t trnumber, Int_t det, Int_t smnumber, 
  Int_t irow, Int_t icol, Float_t adc)
{
  // Add SDigit
  //
  if (!fSDigits) fSDigits = new TClonesArray("AliPMDsdigit", 1000);
  TClonesArray &lsdigits = *fSDigits;
  new(lsdigits[fNsdigit++])  AliPMDsdigit(trnumber,det,smnumber,irow,icol,adc);
}
//____________________________________________________________________________

void AliPMDDigitizer::AddDigit(Int_t trnumber, Int_t det, Int_t smnumber, 
  Int_t irow, Int_t icol, Float_t adc)
{
  // Add Digit
  //
  if (!fDigits) fDigits = new TClonesArray("AliPMDdigit", 1000);
  TClonesArray &ldigits = *fDigits;
  new(ldigits[fNdigit++]) AliPMDdigit(trnumber,det,smnumber,irow,icol,adc);
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
  fCell.Delete();
  for (Int_t i = 0; i < fgkTotUM; i++)
    {
      for (Int_t j = 0; j < fgkRow; j++)
	{
	  for (Int_t k = 0; k < fgkCol; k++)
	    {
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
	      fCPV[i][j][k] = 0.; 
	      fPRE[i][j][k] = 0.; 
	      fPRETrackNo[i][j][k] = 0;
	      fCPVTrackNo[i][j][k] = 0;
	    }
	}
    }
}
//------------------------------------------------------
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

  if(!fCalibData) {
    AliError("No calibration data loaded from CDB!!!");
    return 1;
  }

  Float_t GainFact;
  GainFact = fCalibData->GetGainFact(det,smn,row,col);
  printf("\t gain=%10.3f\n",GainFact);
  return GainFact;
}
//----------------------------------------------------------------------
AliPMDCalibData* AliPMDDigitizer::GetCalibData() const
{
  // The run number will be centralized in AliCDBManager,
  // you don't need to set it here!
  // Added this method by ZA
  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("PMD/Calib/Data");
  
  if(!entry){
    AliWarning("Calibration object retrieval failed! Dummy calibration will be used.");
    
    // this just remembers the actual default storage. No problem if it is null.
    AliCDBStorage *origStorage = AliCDBManager::Instance()->GetDefaultStorage();
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
    
    entry = AliCDBManager::Instance()->Get("PMD/Calib/Data");
    
    // now reset the original default storage to AliCDBManager...
    AliCDBManager::Instance()->SetDefaultStorage(origStorage);  
  }
  
  AliPMDCalibData *calibdata=0;
  if (entry) calibdata = (AliPMDCalibData*) entry->GetObject();
  
  if (!calibdata)  AliError("No calibration data from calibration database !");
  
  return calibdata;
}
