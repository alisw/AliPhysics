#include "AliGenHalo.h"
#include "AliGenMUONlib.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliConst.h"
#include <TDirectory.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TTree.h>
#include <stdlib.h>
 ClassImp(AliGenHalo)
     AliGenHalo::AliGenHalo()
	 :AliGenerator(-1)
{
    fName="Halo";
    fTitle="Halo from LHC Tunnel";
    // Set the default file 
    fFileName="~/marsip/marsip5.mu";
//
//  Read all particles
    fNpart=-1;
    fp=0;
}

AliGenHalo::AliGenHalo(Int_t npart)
    :AliGenerator(npart)
{
    fName="Halo";
    fTitle="Halo from LHC Tunnel";
    // Set the default file 
    fFileName="~/marsip/marsip5.mu";
//
//  Read all particles
    fNpart=-1;
    fp=0;
}

//____________________________________________________________
AliGenHalo::~AliGenHalo()
{
}

//____________________________________________________________
void AliGenHalo::Init() 
{}

//____________________________________________________________
void AliGenHalo::Generate()
{
    FILE *fp = fopen(fFileName,"r");
    if (fp) {
	printf("\n File %s opened for reading ! \n ", fFileName);
    } else {
	printf("\n Opening of file %s failed ! \n ", fFileName);
    }
//
// MARS particle codes
    // const Int_t imars[12]={0,14, 13, 8, 9, 11, 12, 5, 6, 1, 3, 2};
  const Int_t imars[12]={0,kProton,kNeutron,kPiPlus,kPiMinus,kKPlus,kKMinus,
			 kMuonPlus,kMuonMinus,kGamma,kElectron,kPositron};
 
  Float_t polar[3]= {0,0,0};
  Float_t origin[3];
  Float_t p[3], p0;
  Float_t ekin, wgt, tx, ty, tz, txy;
  Float_t amass;
  //
  Int_t ipart, ncols, nt;
  
  Int_t nread=0;
  origin[2]=2650;
  
  while(1) {
      ncols = fscanf(fp,"%i %f %f %f %f %f %f",
		     &ipart, &ekin, &wgt, 
		     &origin[0], &origin[1],
		     &tx, &ty);
      if (ncols < 0) break;
      nread++;
      if (fNpart !=-1 && nread > fNpart) break;
      ipart = imars[ipart];
      amass = TDatabasePDG::Instance()->GetParticle(ipart)->Mass();
      p0=sqrt(ekin*ekin + 2.*amass);
      
      txy=TMath::Sqrt(tx*tx+ty*ty);
      if (txy == 1.) {
	  tz=0;
      } else {
	  tz=-TMath::Sqrt(1.-txy);
      }
      p[0]=p0*tx;
      p[1]=p0*ty;
      p[2]=p0*tz;
      fParentWeight=wgt;
      gAlice->SetTrack(1,-1,ipart,p,origin,polar,0,"Halo+",nt,fParentWeight);
      origin[2]=-origin[2];
      p[2]=-p[2];
      gAlice->SetTrack(1,-1,ipart,p,origin,polar,0,"Halo-",nt,fParentWeight);
      origin[2]=-origin[2];
      p[2]=-p[2];
  }
}
 








