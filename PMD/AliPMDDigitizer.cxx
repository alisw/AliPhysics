//-----------------------------------------------------//
//                                                     //
//  Source File : PMDDigitizer.cxx, Version 00         //
//                                                     //
//  Date   : September 20 2002                         //
//                                                     //
//-----------------------------------------------------//

#include <TBRIK.h>
#include <TNode.h>
#include <TTree.h>
#include <TGeometry.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TParticle.h>

#include "AliRun.h"
#include "AliPMD.h"
#include "AliHit.h"
#include "AliDetector.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliConfig.h"
#include "AliMagF.h"
#include "AliRunDigitizer.h"
#include "AliHeader.h"

#include "AliPMDcell.h"
#include "AliPMDsdigit.h"
#include "AliPMDdigit.h"
#include "AliPMDDigitizer.h"
#include "AliPMDClustering.h"
#include "AliMC.h"

ClassImp(AliPMDDigitizer)
//
// Constructor
//
AliPMDDigitizer::AliPMDDigitizer()
{
  if (!fSDigits) fSDigits = new TClonesArray("AliPMDsdigit", 1000);  
  fNsdigit = 0;
  if (!fDigits) fDigits = new TClonesArray("AliPMDdigit", 1000);  
  fNdigit = 0;

  for (Int_t i = 0; i < fTotUM; i++)
    {
      for (Int_t j = 0; j < fRow; j++)
	{
	  for (Int_t k = 0; k < fCol; k++)
	    {
	      fCPV[i][j][k] = 0.; 
	      fPMD[i][j][k] = 0.; 
	    }
	}
    }

  if (!fCell) fCell = new TObjArray();

  fZPos = 361.5; // in units of cm, This is the default position of PMD
}
AliPMDDigitizer::~AliPMDDigitizer()
{
  delete fSDigits;
  delete fDigits;
  delete fCell;
}
//
// Member functions
//
void AliPMDDigitizer::OpengAliceFile(Char_t *file, Option_t *option)
{

  fRunLoader = AliRunLoader::Open(file,AliConfig::fgkDefaultEventFolderName,
				  "UPDATE");
  
  if (!fRunLoader)
   {
     Error("Open","Can not open session for file %s.",file);
   }
  
  fRunLoader->LoadgAlice();
  fRunLoader->LoadHeader();
  fRunLoader->LoadKinematics();

  gAlice = fRunLoader->GetAliRun();
  
  if (gAlice)
    {
      printf("<AliPMDdigitizer::Open> ");
      printf("AliRun object found on file.\n");
    }
  else
    {
      printf("<AliPMDdigitizer::Open> ");
      printf("Could not find AliRun object.\n");
    }

  PMD  = (AliPMD*)gAlice->GetDetector("PMD");
  pmdloader = fRunLoader->GetLoader("PMDLoader");
  if (pmdloader == 0x0)
    {
      cerr<<"Hits2Digits : Can not find PMD or PMDLoader\n";
    }

  const char *cHS = strstr(option,"HS");
  const char *cHD = strstr(option,"HD");
  const char *cSD = strstr(option,"SD");

  if (cHS)
    {
      pmdloader->LoadHits("READ");
      pmdloader->LoadSDigits("recreate");
    }
  else if (cHD)
    {
      pmdloader->LoadHits("READ");
      pmdloader->LoadDigits("recreate");
    }
  else if (cSD)
    {
      pmdloader->LoadSDigits("READ");
      pmdloader->LoadDigits("recreate");
    }

}
void AliPMDDigitizer::Hits2SDigits(Int_t ievt)
{
  cout << " -------- Beginning of Hits2SDigits ----------- " << endl;

  Int_t kPi0 = 111;
  Int_t kGamma = 22;
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

  
  ResetSDigit();

  printf("Event Number =  %d \n",ievt); 
  Int_t nparticles = fRunLoader->GetHeader()->GetNtrack();
  printf("Number of Particles = %d \n", nparticles);
  fRunLoader->GetEvent(ievt);
  Particles = gAlice->GetMCApp()->Particles();
  // ------------------------------------------------------- //
  // Pointer to specific detector hits.
  // Get pointers to Alice detectors and Hits containers

  treeH = pmdloader->TreeH();
  
  Int_t ntracks    = (Int_t) treeH->GetEntries();
  printf("Number of Tracks in the TreeH = %d \n", ntracks);

  treeS = pmdloader->TreeS();
  if (treeS == 0x0)
    {
      pmdloader->MakeTree("S");
      treeS = pmdloader->TreeS();
    }
  Int_t bufsize = 16000;
  treeS->Branch("PMDSDigit", &fSDigits, bufsize); 
  
  if (PMD) PMDhits   = PMD->Hits();

  // Start loop on tracks in the hits containers

  
  for (Int_t track=0; track<ntracks;track++) 
    {
      gAlice->ResetHits();
      treeH->GetEvent(track);
      
      if (PMD) 
	{
	  npmd = PMDhits->GetEntriesFast();
	  for (int ipmd = 0; ipmd < npmd; ipmd++) 
	    {
	      pmdHit = (AliPMDhit*) PMDhits->UncheckedAt(ipmd);
	      trackno = pmdHit->GetTrack();

	      //  get kinematics of the particles
	      
	      particle = gAlice->GetMCApp()->Particle(trackno);
	      trackpid  = particle->GetPdgCode();

	      Int_t igatr = -999;
	      Int_t ichtr = -999;
	      Int_t igapid = -999;
	      Int_t imo;
	      Int_t igen = 0;
	      Int_t id_mo = -999;

	      TParticle*  mparticle = particle;
	      Int_t trackno_old=0, trackpid_old=0, status_old = 0;
	      if (mparticle->GetFirstMother() == -1)
		{
		  trackno_old  = trackno;
		  trackpid_old = trackpid;
		  status_old   = -1;
		}

	      Int_t igstatus = 0;
	      while((imo = mparticle->GetFirstMother()) >= 0)
		{
		  igen++;
		  mparticle =  gAlice->GetMCApp()->Particle(imo);
		  id_mo = mparticle->GetPdgCode();
		  
		  vx = mparticle->Vx();
		  vy = mparticle->Vy();
		  vz = mparticle->Vz();
		
		  //printf("==> Mother ID %5d %5d %5d Vertex: %13.3f %13.3f %13.3f\n", igen, imo, id_mo, vx, vy, vz);
		  //fprintf(fpw1,"==> Mother ID %5d %5d %5d Vertex: %13.3f %13.3f %13.3f\n", igen, imo, id_mo, vx, vy, vz);
		  if ((id_mo == kGamma || id_mo == -11 || id_mo == 11) && vx == 0. && vy == 0. && vz == 0.)
		    {
		      igatr = imo;
		      igapid = id_mo;
		      igstatus = 1;
		    }
		  if(igstatus == 0)
		    {
		      if (id_mo == kPi0 && vx == 0. && vy == 0. && vz == 0.)
			{
			  igatr = imo;
			  igapid = id_mo;
			}
		    }
		  ichtr = imo;
		}

	      if (id_mo == kPi0 && vx == 0. && vy == 0. && vz == 0.)
		{
		  mtrackno = igatr;
		  mtrackpid = igapid;
		}
	      else
		{
		  mtrackno  = ichtr;
		  mtrackpid = id_mo;
		}
	      if (status_old == -1)
		{
		  mtrackno  = trackno_old;
		  mtrackpid = trackpid_old;
		}
	      
	      xPos = pmdHit->X();
	      yPos = pmdHit->Y();
	      zPos = pmdHit->Z();
	      edep       = pmdHit->fEnergy;
	      Int_t Vol1 = pmdHit->fVolume[1]; // Column
	      Int_t Vol2 = pmdHit->fVolume[2]; // Row
	      Int_t Vol3 = pmdHit->fVolume[3]; // UnitModule
	      Int_t Vol6 = pmdHit->fVolume[6]; // SuperModule
	      // -----------------------------------------//
	      // For Super Module 1 & 2                   //
	      //  nrow = 96, ncol = 48                    //
	      // For Super Module 3 & 4                   //
	      //  nrow = 48, ncol = 96                    //
	      // -----------------------------------------//
	      
	      smnumber = (Vol6-1)*6 + Vol3;

	      if (Vol6 == 1 || Vol6 == 2)
		{
		  xpad = Vol1;
		  ypad = Vol2;
		}
	      else if (Vol6 == 3 || Vol6 == 4)
		{
		  xpad = Vol2;
		  ypad = Vol1;
		}

	      //cout << "zpos = " << zPos << " edep = " << edep << endl;

	      Float_t zposition = TMath::Abs(zPos);
	      if (zposition < fZPos)
		{
		  // CPV
		  fDetNo = 1;
		}
	      else if (zposition > fZPos)
		{
		  // PMD
		  fDetNo = 0;
		}
	      Int_t smn = smnumber - 1;
	      Int_t ixx = xpad     - 1;
	      Int_t iyy = ypad     - 1;
	      if (fDetNo == 0)
		{
		  fPMD[smn][ixx][iyy] += edep;
		  fPMDCounter[smn][ixx][iyy]++;

		  pmdcell = new AliPMDcell(mtrackno,smn,ixx,iyy,edep);

		  fCell->Add(pmdcell);
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
  Int_t   cellno      = 0;

  for (Int_t idet = 0; idet < 2; idet++)
    {
      for (Int_t ism = 0; ism < fTotUM; ism++)
	{
	  for (Int_t jrow = 0; jrow < fRow; jrow++)
	    {
	      for (Int_t kcol = 0; kcol < fCol; kcol++)
		{
		  cellno = jrow*fCol + kcol;
		  if (idet == 0)
		    {
		      deltaE = fPMD[ism][jrow][kcol];
		      trno   = fPMDTrackNo[ism][jrow][kcol];
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
		      AddSDigit(trno,detno,ism,cellno,deltaE);
		    }
		}
	    }
	  treeS->Fill();
	  ResetSDigit();
	}
    }
  pmdloader->WriteSDigits("OVERWRITE");
  ResetCellADC();

  //  cout << " -------- End of Hits2SDigit ----------- " << endl;
}

void AliPMDDigitizer::Hits2Digits(Int_t ievt)
{
  Int_t kPi0 = 111;
  Int_t kGamma = 22;
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

  
  ResetDigit();

  printf("Event Number =  %d \n",ievt); 

  Int_t nparticles = fRunLoader->GetHeader()->GetNtrack();
  printf("Number of Particles = %d \n", nparticles);
  fRunLoader->GetEvent(ievt);
  Particles = gAlice->GetMCApp()->Particles();
  // ------------------------------------------------------- //
  // Pointer to specific detector hits.
  // Get pointers to Alice detectors and Hits containers

  PMD  = (AliPMD*)gAlice->GetDetector("PMD");
  pmdloader = fRunLoader->GetLoader("PMDLoader");

  if (pmdloader == 0x0)
    {
      cerr<<"Hits2Digits method : Can not find PMD or PMDLoader\n";
    }
  treeH = pmdloader->TreeH();
  Int_t ntracks    = (Int_t) treeH->GetEntries();
  printf("Number of Tracks in the TreeH = %d \n", ntracks);
  pmdloader->LoadDigits("recreate");
  treeD = pmdloader->TreeD();
  if (treeD == 0x0)
    {
      pmdloader->MakeTree("D");
      treeD = pmdloader->TreeD();
    }
  Int_t bufsize = 16000;
  treeD->Branch("PMDDigit", &fDigits, bufsize); 
  
  if (PMD) PMDhits   = PMD->Hits();

  // Start loop on tracks in the hits containers

  for (Int_t track=0; track<ntracks;track++) 
    {
      gAlice->ResetHits();
      treeH->GetEvent(track);
      
      if (PMD) 
	{
	  npmd = PMDhits->GetEntriesFast();
	  for (int ipmd = 0; ipmd < npmd; ipmd++) 
	    {
	      pmdHit = (AliPMDhit*) PMDhits->UncheckedAt(ipmd);
	      trackno = pmdHit->GetTrack();
	      
	      //  get kinematics of the particles
	      
	      particle = gAlice->GetMCApp()->Particle(trackno);
	      trackpid  = particle->GetPdgCode();

	      Int_t igatr = -999;
	      Int_t ichtr = -999;
	      Int_t igapid = -999;
	      Int_t imo;
	      Int_t igen = 0;
	      Int_t id_mo = -999;

	      TParticle*  mparticle = particle;
	      Int_t trackno_old=0, trackpid_old=0, status_old = 0;
	      if (mparticle->GetFirstMother() == -1)
		{
		  trackno_old  = trackno;
		  trackpid_old = trackpid;
		  status_old   = -1;
		}

	      Int_t igstatus = 0;
	      while((imo = mparticle->GetFirstMother()) >= 0)
		{
		  igen++;
		  mparticle =  gAlice->GetMCApp()->Particle(imo);
		  id_mo = mparticle->GetPdgCode();
		  
		  vx = mparticle->Vx();
		  vy = mparticle->Vy();
		  vz = mparticle->Vz();
		
		  //printf("==> Mother ID %5d %5d %5d Vertex: %13.3f %13.3f %13.3f\n", igen, imo, id_mo, vx, vy, vz);
		  //fprintf(fpw1,"==> Mother ID %5d %5d %5d Vertex: %13.3f %13.3f %13.3f\n", igen, imo, id_mo, vx, vy, vz);
		  if ((id_mo == kGamma || id_mo == -11 || id_mo == 11) && vx == 0. && vy == 0. && vz == 0.)
		    {
		      igatr = imo;
		      igapid = id_mo;
		      igstatus = 1;
		    }
		  if(igstatus == 0)
		    {
		      if (id_mo == kPi0 && vx == 0. && vy == 0. && vz == 0.)
			{
			  igatr = imo;
			  igapid = id_mo;
			}
		    }
		  ichtr = imo;
		}

	      if (id_mo == kPi0 && vx == 0. && vy == 0. && vz == 0.)
		{
		  mtrackno = igatr;
		  mtrackpid = igapid;
		}
	      else
		{
		  mtrackno  = ichtr;
		  mtrackpid = id_mo;
		}
	      if (status_old == -1)
		{
		  mtrackno  = trackno_old;
		  mtrackpid = trackpid_old;
		}
	      
	      xPos = pmdHit->X();
	      yPos = pmdHit->Y();
	      zPos = pmdHit->Z();
	      edep       = pmdHit->fEnergy;
	      Int_t Vol1 = pmdHit->fVolume[1]; // Column
	      Int_t Vol2 = pmdHit->fVolume[2]; // Row
	      Int_t Vol3 = pmdHit->fVolume[3]; // UnitModule
	      Int_t Vol6 = pmdHit->fVolume[6]; // SuperModule

	      // -----------------------------------------//
	      // For Super Module 1 & 2                   //
	      //  nrow = 96, ncol = 48                    //
	      // For Super Module 3 & 4                   //
	      //  nrow = 48, ncol = 96                    //
	      // -----------------------------------------//
	      
	      smnumber = (Vol6-1)*6 + Vol3;

	      if (Vol6 == 1 || Vol6 == 2)
		{
		  xpad = Vol1;
		  ypad = Vol2;
		}
	      else if (Vol6 == 3 || Vol6 == 4)
		{
		  xpad = Vol2;
		  ypad = Vol1;
		}

	      //cout << "-zpos = " << -zPos << endl;

	      Float_t zposition = TMath::Abs(zPos);

	      if (zposition < fZPos)
		{
		  // CPV
		  fDetNo = 1;
		}
	      else if (zposition > fZPos)
		{
		  // PMD
		  fDetNo = 0;
		}

	      Int_t smn = smnumber - 1;
	      Int_t ixx = xpad     - 1;
	      Int_t iyy = ypad     - 1;
	      if (fDetNo == 0)
		{
		  fPMD[smn][ixx][iyy] += edep;
		  fPMDCounter[smn][ixx][iyy]++;

		  pmdcell = new AliPMDcell(mtrackno,smn,ixx,iyy,edep);

		  fCell->Add(pmdcell);
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

  Float_t deltaE = 0.;
  Int_t detno = 0;
  Int_t trno = 1;
  Int_t cellno;

  for (Int_t idet = 0; idet < 2; idet++)
    {
      for (Int_t ism = 0; ism < fTotUM; ism++)
	{
	  for (Int_t jrow = 0; jrow < fRow; jrow++)
	    {
	      for (Int_t kcol = 0; kcol < fCol; kcol++)
		{
		  cellno = jrow*fCol + kcol;
		  if (idet == 0)
		    {
		      deltaE = fPMD[ism][jrow][kcol];
		      trno   = fPMDTrackNo[ism][jrow][kcol];
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
		      AddDigit(trno,detno,ism,cellno,deltaE);
		    }
		} // column loop
	    } // row    loop
	} // supermodule loop
      treeD->Fill();
      ResetDigit();
    } // detector loop

  pmdloader->WriteDigits("OVERWRITE");
  ResetCellADC();
  
  //  cout << " -------- End of Hits2Digit ----------- " << endl;
}


void AliPMDDigitizer::SDigits2Digits(Int_t ievt)
{
  //  cout << " -------- Beginning of SDigits2Digit ----------- " << endl;
  fRunLoader->GetEvent(ievt);

  treeS = pmdloader->TreeS();
  AliPMDsdigit  *pmdsdigit;
  TBranch *branch = treeS->GetBranch("PMDSDigit");
  branch->SetAddress(&fSDigits);

  treeD = pmdloader->TreeD();
  if (treeD == 0x0)
    {
      pmdloader->MakeTree("D");
      treeD = pmdloader->TreeD();
    }
  Int_t bufsize = 16000;
  treeD->Branch("PMDDigit", &fDigits, bufsize); 

  Int_t   trno, det, smn;
  Int_t   cellno;
  Float_t edep, adc;

  Int_t nmodules = (Int_t) treeS->GetEntries();

  for (Int_t imodule = 0; imodule < nmodules; imodule++)
    {
      treeS->GetEntry(imodule); 
      Int_t nentries = fSDigits->GetLast();
      //cout << " nentries = " << nentries << endl;
      for (Int_t ient = 0; ient < nentries+1; ient++)
	{
	  pmdsdigit = (AliPMDsdigit*)fSDigits->UncheckedAt(ient);
	  trno   = pmdsdigit->GetTrackNumber();
	  det    = pmdsdigit->GetDetector();
	  smn    = pmdsdigit->GetSMNumber();
	  cellno = pmdsdigit->GetCellNumber();
	  edep   = pmdsdigit->GetCellEdep();

	  MeV2ADC(edep,adc);
	  AddDigit(trno,det,smn,cellno,adc);      
	}
      treeD->Fill();
      ResetDigit();
    }
  pmdloader->WriteDigits("OVERWRITE");
  //  cout << " -------- End of SDigits2Digit ----------- " << endl;
}

void AliPMDDigitizer::TrackAssignment2Cell()
{
  // 
  // This block assigns the cell id when there are
  // multiple tracks in a cell according to the
  // energy deposition
  //
  Bool_t jsort = false;

  Int_t i, j, k;

  Float_t *frac_edp;
  Float_t *tr_edp;
  Int_t *status1;
  Int_t *status2;
  Int_t *trnarray;
  Int_t   ****PMDTrack;
  Float_t ****PMDEdep;

  PMDTrack = new Int_t ***[fTotUM];
  PMDEdep  = new Float_t ***[fTotUM];
  for (i=0; i<fTotUM; i++)
    {
      PMDTrack[i] = new Int_t **[fRow];
      PMDEdep[i]  = new Float_t **[fRow];
    }

  for (i = 0; i < fTotUM; i++)
    {
      for (j = 0; j < fRow; j++)
	{
	  PMDTrack[i][j] = new Int_t *[fCol];
	  PMDEdep[i][j]  = new Float_t *[fCol];
	}
    }
  
  for (i = 0; i < fTotUM; i++)
    {
      for (j = 0; j < fRow; j++)
	{
	  for (k = 0; k < fCol; k++)
	    {
	      Int_t nn = fPMDCounter[i][j][k];
	      if(nn > 0)
		{
		  PMDTrack[i][j][k] = new Int_t[nn];
		  PMDEdep[i][j][k] = new Float_t[nn];
		}
	      else
		{
		  nn = 1;
		  PMDTrack[i][j][k] = new Int_t[nn];
		  PMDEdep[i][j][k] = new Float_t[nn];
		}		      
	      fPMDCounter[i][j][k] = 0;
	    }
	}
    }


  Int_t nentries = fCell->GetEntries();

  Int_t   mtrackno, ism, ixp, iyp;
  Float_t edep;

  for (i = 0; i < nentries; i++)
    {
      pmdcell = (AliPMDcell*)fCell->UncheckedAt(i);
      
      mtrackno = pmdcell->GetTrackNumber();
      ism = pmdcell->GetSMNumber();
      ixp = pmdcell->GetX();
      iyp = pmdcell->GetY();
      edep = pmdcell->GetEdep();
      Int_t nn = fPMDCounter[ism][ixp][iyp];
      //      cout << " nn = " << nn << endl;
      PMDTrack[ism][ixp][iyp][nn] = (Int_t) mtrackno;
      PMDEdep[ism][ixp][iyp][nn] = edep;
      fPMDCounter[ism][ixp][iyp]++;
    }
  
  Int_t iz, il;
  Int_t im, ix, iy;
  Int_t nn;
  
  for (im=0; im<fTotUM; im++)
    {
      for (ix=0; ix<fRow; ix++)
	{
	  for (iy=0; iy<fCol; iy++)
	    {
	      nn = fPMDCounter[im][ix][iy];
	      if (nn > 1)
		{
		  // This block handles if a cell is fired
		  // many times by many tracks
		  status1  = new Int_t[nn];
		  status2  = new Int_t[nn];
		  trnarray = new Int_t[nn];
		  for (iz = 0; iz < nn; iz++)
		    {
		      status1[iz] = PMDTrack[im][ix][iy][iz];
		    }
		  TMath::Sort(nn,status1,status2,jsort);
		  Int_t track_old = -99999;
		  Int_t track, tr_count = 0;
		  for (iz = 0; iz < nn; iz++)
		    {
		      track = status1[status2[iz]];
		      if (track_old != track)
			{
			  trnarray[tr_count] = track;
			  tr_count++;
			}			      
		      track_old = track;
		    }
		  delete status1;
		  delete status2;
		  Float_t tot_edp = 0.;
		  tr_edp = new Float_t[tr_count];
		  frac_edp = new Float_t[tr_count];
		  for (il = 0; il < tr_count; il++)
		    {
		      tr_edp[il] = 0.;
		      track = trnarray[il];
		      for (iz = 0; iz < nn; iz++)
			{
			  if (track == PMDTrack[im][ix][iy][iz])
			    {
			      tr_edp[il] += PMDEdep[im][ix][iy][iz];
			    }
			}
		      tot_edp += tr_edp[il];
		    }
		  Int_t il_old = 0;
		  Float_t frac_old = 0.;
		  
		  for (il = 0; il < tr_count; il++)
		    {
		      frac_edp[il] = tr_edp[il]/tot_edp;
		      if (frac_old < frac_edp[il])
			{
			  frac_old = frac_edp[il];
			  il_old = il;
			}
		    }
		  fPMDTrackNo[im][ix][iy] = trnarray[il_old];
		  delete frac_edp;
		  delete tr_edp;
		  delete trnarray;
		}
	      else if (nn == 1)
		{
		  // This only handles if a cell is fired
		  // by only one track
		  
		  fPMDTrackNo[im][ix][iy] = PMDTrack[im][ix][iy][0];
		  
		}
	      else if (nn ==0)
		{
		  // This is if no cell is fired
		  fPMDTrackNo[im][ix][iy] = -999;
		}
	    } // end of iy
	} // end of ix
    } // end of im
  
  // Delete all the pointers
  
  for (i = 0; i < fTotUM; i++)
    {
      for (j = 0; j < fRow; j++)
	{
	  for (k = 0; k < fCol; k++)
	    {
	      delete [] PMDTrack[i][j][k];
	      delete [] PMDEdep[i][j][k];
	    }
	}
    }
  
  for (i = 0; i < fTotUM; i++)
    {
      for (j = 0; j < fRow; j++)
	{
	  delete [] PMDTrack[i][j];
	  delete [] PMDEdep[i][j];
	}
    }
  
  for (i = 0; i < fTotUM; i++)
    {
      delete [] PMDTrack[i];
      delete [] PMDEdep[i];
    }
  delete PMDTrack;
  delete PMDEdep;
  // 
  // End of the cell id assignment
  //
}


void AliPMDDigitizer::MeV2ADC(Float_t mev, Float_t & adc)
{
  // To be done

  adc = mev*1.;
}
void AliPMDDigitizer::AddSDigit(Int_t trnumber, Int_t det, Int_t smnumber, 
  Int_t cellnumber, Float_t adc)
{
  TClonesArray &lsdigits = *fSDigits;
  AliPMDsdigit *newcell;
  newcell = new AliPMDsdigit(trnumber,det,smnumber,cellnumber,adc);
  new(lsdigits[fNsdigit++]) AliPMDsdigit(newcell);
  delete newcell;
}

void AliPMDDigitizer::AddDigit(Int_t trnumber, Int_t det, Int_t smnumber, 
  Int_t cellnumber, Float_t adc)
{
  TClonesArray &ldigits = *fDigits;
  AliPMDdigit *newcell;
  newcell = new AliPMDdigit(trnumber,det,smnumber,cellnumber,adc);
  new(ldigits[fNdigit++]) AliPMDdigit(newcell);
  delete newcell;
}

Int_t AliPMDDigitizer::Convert2RealSMNumber(Int_t smnumber1)
{
  Int_t smnumber = -999;

  if (smnumber1==1)  smnumber =  1;
  if (smnumber1==2)  smnumber = 10;
  if (smnumber1==3)  smnumber = 19;
  if (smnumber1==4)  smnumber =  1;
  if (smnumber1==5)  smnumber = 10;
  if (smnumber1==6)  smnumber = 19;
  if (smnumber1==7)  smnumber =  2;
  if (smnumber1==8)  smnumber =  3;
  if (smnumber1==9)  smnumber =  4;
  if (smnumber1==10) smnumber =  5;
  if (smnumber1==11) smnumber =  6;
  if (smnumber1==12) smnumber =  7;
  if (smnumber1==13) smnumber =  8;
  if (smnumber1==14) smnumber =  9;
  if (smnumber1==15) smnumber = 11;
  if (smnumber1==16) smnumber = 12;
  if (smnumber1==17) smnumber = 13;
  if (smnumber1==18) smnumber = 14;
  if (smnumber1==19) smnumber = 15;
  if (smnumber1==20) smnumber = 16;
  if (smnumber1==21) smnumber = 17;
  if (smnumber1==22) smnumber = 18;
  if (smnumber1==23) smnumber = 20;
  if (smnumber1==24) smnumber = 21;
  if (smnumber1==25) smnumber = 22;
  if (smnumber1==26) smnumber = 23;
  if (smnumber1==27) smnumber = 24;
  if (smnumber1==28) smnumber = 25;
  if (smnumber1==29) smnumber = 26;
  if (smnumber1==30) smnumber = 27;

  return smnumber;
}
void AliPMDDigitizer::SetZPosition(Float_t zpos)
{
  fZPos = zpos;
}
Float_t AliPMDDigitizer::GetZPosition() const
{
  return fZPos;
}

void AliPMDDigitizer::ResetCell()
{
  fCell->Clear();
  for (Int_t i = 0; i < fTotUM; i++)
    {
      for (Int_t j = 0; j < fRow; j++)
	{
	  for (Int_t k = 0; k < fCol; k++)
	    {
	      fPMDCounter[i][j][k] = 0; 
	    }
	}
    }
}
void AliPMDDigitizer::ResetSDigit()
{
  fNsdigit = 0;
  if (fSDigits) fSDigits->Clear();
}
void AliPMDDigitizer::ResetDigit()
{
  fNdigit = 0;
  if (fDigits) fDigits->Clear();
}

void AliPMDDigitizer::ResetCellADC()
{
  for (Int_t i = 0; i < fTotUM; i++)
    {
      for (Int_t j = 0; j < fRow; j++)
	{
	  for (Int_t k = 0; k < fCol; k++)
	    {
	      fCPV[i][j][k] = 0.; 
	      fPMD[i][j][k] = 0.; 
	    }
	}
    }
}

void AliPMDDigitizer::UnLoad(Option_t *option)
{
  const char *cS = strstr(option,"S");
  const char *cD = strstr(option,"D");

  fRunLoader->UnloadgAlice();
  fRunLoader->UnloadHeader();
  fRunLoader->UnloadKinematics();

  if (cS)
    {
      pmdloader->UnloadHits();
    }
  if (cD)
    {
      pmdloader->UnloadHits();
      pmdloader->UnloadSDigits();
    }
}
