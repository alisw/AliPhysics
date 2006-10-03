//*KEEP,CopyRight,T=C.
/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
//*KEND.

//*    +-------------------------------------------------------------+
//*    |                                                             |
//*    |                                                             |
//*    |                        DPMJET 3.0                           |
//*    |                                                             |
//*    |							     |
//*    |	 S. Roesler+), R. Engel#), J. Ranft*)		     |
//*    |							     |
//*    |	 +) CERN, TIS-RP				     |
//*    |	    CH-1211 Geneva 23, Switzerland		     |
//*    |	    Email: Stefan.Roesler@cern.ch		     |
//*    |							     |
//*    |	 #) University of Delaware, BRI 		     |
//*    |	    Newark, DE 19716, USA			     |
//*    |							     |
//*    |	 *) University of Siegen, Dept. of Physics	     |
//*    |	    D-57068 Siegen, Germany			     |
//*    |							     |
//*    |							     |
//*    |       http://home.cern.ch/sroesler/dpmjet3.html	     |
//*    |							     |
//*    |							     |
//*    |       Monte Carlo models used for event generation:	     |
//*    |	  PHOJET 1.12, JETSET 7.4 and LEPTO 6.5.1	     |
//*    |							     |
//*    +-------------------------------------------------------------+

//*KEEP,TDPMjet.
#include "TDPMjet.h"
//*KEEP,DPMCOMMON.
#include "DPMcommon.h"
//*KEEP,TParticle,T=C++.
#include "TParticle.h"
//*KEND.

//*KEEP,TROOT.
#include "TROOT.h"
//*KEND.

#ifndef WIN32
# define dt_dtuini dt_dtuini_
# define dt_getemu de_getemu_
# define dt_kkinc  dt_kkinc_
# define pho_phist pho_phist_
# define dt_dtuout dt_dtuout_
# define dt_rndm   dt_rndm_
# define dt_rndmst dt_rndmst_
# define dt_rndmin dt_rndmin_
# define dt_rndmou dt_rndmou_
# define type_of_call
#else
# define dt_dtuini DT_DTUINI
# define dt_getemu DT_GETEMU
# define dt_kkinc  DT_KKINC
# define pho_phist PHO_PHIST
# define dt_dtuout DT_DTUOUT
# define dt_rndm   DT_RNDM
# define dt_rndmst DT_RNDMST
# define dt_rndmin DT_RNDMIN
# define dt_rndmou DT_RNDMOU
# define type_of_call _stdcall
#endif

#ifndef WIN32
extern "C" void   type_of_call dt_dtuini(Int_t & , Double_t &, Int_t & , Int_t &, 
			Int_t &, Int_t &, Int_t &, Int_t &);
extern "C" double type_of_call dt_getemu(Int_t &, Int_t &, Int_t &, Int_t &);
extern "C" void   type_of_call dt_kkinc(Int_t &, Int_t &, Int_t &, Int_t &,
				    Int_t &, Double_t &, Int_t &, Int_t &);
extern "C" void   type_of_call pho_phist(Int_t &, Double_t &);
extern "C" void   type_of_call dt_dtuout();
extern "C" void   type_of_call dt_rndm(Int_t &);
extern "C" void   type_of_call dt_rndmst(Int_t &, Int_t &, Int_t &, Int_t &);
extern "C" void   type_of_call dt_rndmin(Int_t &, Int_t &, Int_t &, Int_t &, Int_t &, Int_t &);
extern "C" void   type_of_call dt_rndmou(Int_t &, Int_t &, Int_t &, Int_t &, Int_t &, Int_t &);

#else

#endif

ClassImp(TDPMjet)


//______________________________________________________________________________
    TDPMjet::TDPMjet() : 
	TGenerator("dpmjet","dpmjet"),
	fNEvent(0),
	fIp(0),
	fIpz(0),
	fIt(0),
	fItz(0),
	fEpn(0.),
	fPpn(0.),
	fCMEn(0.),
	fIdp(0),
	fBmin(0.),
	fBmax(0.),
	fFCentr(0),
	fPi0Decay(0),
	fProcess(kDpmMb)
{
// Default Constructor
}

//______________________________________________________________________________
TDPMjet::TDPMjet(DpmProcess_t  iproc, Int_t Ip=208, Int_t Ipz=82, Int_t It=208, Int_t Itz=82, 
		 Double_t Epn=2700., Double_t CMEn=5400.) 
    : TGenerator("dpmjet","dpmjet"),
      fNEvent(0),
      fIp(Ip),
      fIpz(Ipz),
      fIt(It),
      fItz(Itz),
      fEpn(Epn),
      fCMEn(CMEn),
      fIdp(0),
      fBmin(0.),
      fBmax(0.),
      fFCentr(0),
      fPi0Decay(0),
      fProcess(iproc)
{  
    printf("TDPMJet Constructor %d %d %d %d \n", Ip, Ipz, It, Itz);
}


//______________________________________________________________________________
Int_t TDPMjet::ImportParticles(TClonesArray *particles, Option_t *option)
{
//
//  Default primary creation method. It reads the /HEPEVT/ common block which
//  has been filled by the GenerateEvent method. If the event generator does
//  not use the HEPEVT common block, This routine has to be overloaded by
//  the subclasses.
//  The function loops on the generated particles and store them in
//  the TClonesArray pointed by the argument particles.
//  The default action is to store only the stable particles 
//  This can be demanded explicitly by setting the option = "Final"
//  If the option = "All", all the particles are stored.
//
  if(particles==0) return 0;
  TClonesArray &Particles = *particles;
  Particles.Clear();
  Int_t numpart = 0; 	 // Total number of produced particles
  Int_t numStabpart = 0; // Total number of produced stable particles
  Double_t entot = 0; 	 // Total energy in final state (from stable particles)
  
  numpart = DTEVT1.nhkk;
  for(Int_t i=0; i<numpart; i++){
     if(DTEVT1.isthkk[i]==1 || DTEVT1.isthkk[i]==-1 || DTEVT1.isthkk[i]==1001){
        numStabpart++;
	entot += DTEVT1.phkk[i][3]; // PHKK[i][3] <-> PHKK(4,i)
     } 
  }
  //printf("\n TDPMjet: DPMJET stack contains %d particles", numpart);
  // printf("\n TDPMjet: Final not decayed particles: %d",    numStabpart);
  //printf("\n TDPMjet: Total energy: %f GeV          \n",   entot);
  Int_t nump = 0;
  
  if(!strcmp(option,"") || !strcmp(option,"Final")){
      for (Int_t i=0; i < numpart; i++) {
	  
	  if (DTEVT1.isthkk[i] == 1) {
	   //
	   //  Use the common block values for the TParticle constructor
	   //
	      new(Particles[nump]) TParticle(
		  DTEVT1.idhkk[i],
		  DTEVT1.isthkk[i],
		  -1,
		  -1,
		  -1,
		  -1,
		  DTEVT1.phkk[i][0],
		  DTEVT1.phkk[i][1],
		  DTEVT1.phkk[i][2],
		  DTEVT1.phkk[i][3],
		  
		  DTEVT1.vhkk[i][0],
		  DTEVT1.vhkk[i][1],
		  DTEVT1.vhkk[i][2],
		  DTEVT1.vhkk[i][3]);
	      nump++;
	  }
      }
  }
  else if(!strcmp(option,"All")){
      nump = numpart; 
      for (Int_t i=0; i <= numpart; i++){
        
	  // DTEVT1.JMOHKK[i][0] pointer to the entry of the 1st mother of entry i
	  Int_t iParent = DTEVT1.jmohkk[i][0] - 1;
	     
	  if(iParent >= 0){
	      TParticle *mother = (TParticle*) (Particles.UncheckedAt(iParent));
	      mother->SetLastDaughter(i);
	      if(mother->GetFirstDaughter() == -1) mother->SetFirstDaughter(i);
	  } 
	  // --- PDGcode for residual nuclei (idhkk=80000) 
	  // --- 	10000*Z + 10*A
	  // --- DPMJET -> idres = mass #, idxres = charge
	  if(DTEVT1.idhkk[i] == 80000) 
	      DTEVT1.idhkk[i] = 10000*DTEVT2.idxres[i]+10*DTEVT2.idres[i];
/*
	  if(DTEVT2.idxres[i] != 0) 
	      printf("\n	pc#%d -> A = %d, Z = %d -> PDGcode = %d\n",
		     i,DTEVT2.idres[i],DTEVT2.idxres[i],DTEVT1.idhkk[i]);
*/	  
	  new(Particles[i]) TParticle(
	      DTEVT1.idhkk[i],
	      DTEVT1.isthkk[i],
	      iParent,
	      -1,
	      -1,
	      -1,
	      
	      DTEVT1.phkk[i][0],
	      DTEVT1.phkk[i][1],
	      DTEVT1.phkk[i][2],
	      DTEVT1.phkk[i][3],
	      
	      DTEVT1.vhkk[i][0],
	      DTEVT1.vhkk[i][1],
	      DTEVT1.vhkk[i][2],
	      DTEVT1.vhkk[i][3]);
      } // Particle loop  
  }
  return nump;
}


//====================== access to dpmjet subroutines =========================
//______________________________________________________________________________
void TDPMjet::Initialize()
{
//
//  Write standard DPMJET input cards 
//
    FILE* out = fopen("dpmjet.inp","w");
//  Projectile and Target definition 
    fprintf(out, "PROJPAR   %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n", (Float_t) fIp, (Float_t) fIpz,  0., 0., 0., 0.);
    fprintf(out, "TARPAR    %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n", (Float_t) fIt, (Float_t) fItz,  0., 0., 0., 0.);
//  Beam energy and crossing-angle
    fprintf(out, "BEAM      %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n",fEpn, fEpn, 0., 0., 0., 0.);
//  Centrality
    fprintf(out, "CENTRAL   %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n",-1., fBmin, fBmax, 0., 0., 0.);
//  Particle decays
    if (fPi0Decay) 
    fprintf(out, "PARDECAY  %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n", 2., 0., 0., 0., 0., 0.);    
//
//  PHOJET specific
    fprintf(out, "PHOINPUT\n");
    fprintf(out, "DEBUG      0 0 0 \n");
    
    if (fProcess == kDpmMb) {
	fprintf(out, "PROCESS           1 0 1 1 1 1 1 1\n");
    } else if (fProcess == kDpmMbNonDiffr) {
	fprintf(out, "PROCESS           1 0 1 1 0 0 0 1\n");
    }
    
    fprintf(out, "ENDINPUT\n");
//
//  START card
    fprintf(out, "START            1.0       0.0\n");
    fprintf(out, "STOP\n");
    fclose(out);
//
//  Call DPMJET initialisation
    Int_t iemu = 0; // No emulsion (default)
    Dt_Dtuini(1, fEpn, fIp, fIpz, fIt, fItz, fIdp, iemu);

}


//______________________________________________________________________________
void TDPMjet::GenerateEvent()
{
   // Generates one event;
   fNEvent++;
   DTEVNO.nevent=fNEvent;
   Int_t kkmat=-1;
   Float_t Elab = fEpn;
   Int_t irej=0;
   Dt_Kkinc(fIp, fIpz, fIt, fItz, fIdp, Elab, kkmat, irej);
   if(irej!=0) return;
}
//______________________________________________________________________________
void TDPMjet::Dt_Dtuini(int nevts, double epn, int npmass, int npchar, 
   			int ntmass, int ntchar, int idp, int iemu)
{
  // Call dmpjet routine DT_DTUINI passing the parameters 
  // in a way accepted by Fortran routines				   
     

   printf("\n-------------------------------------------\n");
   printf("\n		Dt_Dtuini called with:\n\n");
   printf(" Projectile	-> A = %d, Z = %d \n",npmass, npchar);
   printf(" Target    	-> A = %d, Z = %d \n",ntmass, ntchar);
   printf(" Proj. LAB E	-> E = %f GeV \n",epn);
   printf(" nevts = %d, idp = %d, iemu = %d \n",nevts,idp,iemu);
   printf("\n-------------------------------------------\n");

   dt_dtuini(nevts, epn, npmass, npchar, ntmass, ntchar, idp, iemu);
    
}

//______________________________________________________________________________
void TDPMjet::Dt_Kkinc(int npmass, int npchar, int ntmass, int ntchar, 
   		       int idp, double elab, int kkmat, int irej)
{
  // Call dmpjet routine DT_KKINC passing the parameters 
  // in a way accepted by Fortran routines				   
  dt_kkinc(npmass, npchar, ntmass, ntchar, idp, elab, kkmat, irej);
}

//______________________________________________________________________________
void TDPMjet::Pho_Phist(int imode, double weight)
{
  // Call dmpjet routine PHO_PHIST passing the parameters 
  // in a way accepted by Fortran routines
  
  pho_phist(imode,weight);				   

}

//______________________________________________________________________________
void TDPMjet::Dt_Dtuout()
{
  // Call dmpjet routine DT_DTUOT passing the parameters 
  // in a way accepted by Fortran routines				   
  
  dt_dtuout();

}

//______________________________________________________________________________
Int_t TDPMjet::GetEvNum() const
{
	return DTEVT1.nevhkk;
}
//______________________________________________________________________________
Int_t TDPMjet::GetEntriesNum() const
{
	return DTEVT1.nhkk;
}
//______________________________________________________________________________
Int_t TDPMjet::GetNumStablePc() const
{
	Int_t NumStablePc = 0;
	for(Int_t i=0; i<DTEVT1.nhkk; i++){
	   if(DTEVT1.isthkk[i] == 1) NumStablePc++;
	}
	return NumStablePc;
}

//______________________________________________________________________________
Float_t TDPMjet::GetTotEnergy() const
{
  	Float_t TotEnergy = 0.;
  	for(Int_t i=0; i<DTEVT1.nhkk; i++){
     	  if(DTEVT1.isthkk[i] == 1)
	    TotEnergy += DTEVT1.phkk[i][3]; // PHKK[i][3] <-> PHKK(4,i)
  	}
	return TotEnergy;
}

//______________________________________________________________________________
Int_t TDPMjet::GetStatusCode(Int_t evnum) const 
{
	return DTEVT1.isthkk[evnum];	
}
//______________________________________________________________________________
Int_t TDPMjet::GetPDGCode(Int_t evnum) const   
{
	return DTEVT1.idhkk[evnum];
}
//______________________________________________________________________________
Double_t TDPMjet::Getpx(Int_t evnum) const       
{
	return DTEVT1.phkk[evnum][0];
}
//______________________________________________________________________________
Double_t TDPMjet::Getpy(Int_t evnum) const      
{
	return DTEVT1.phkk[evnum][1];
}
//______________________________________________________________________________
Double_t TDPMjet::Getpz(Int_t evnum) const       
{
	return DTEVT1.phkk[evnum][2];
}
//______________________________________________________________________________
Double_t TDPMjet::GetEnergy(Int_t evnum) const	      
{
	return DTEVT1.phkk[evnum][3];
}
//______________________________________________________________________________
Double_t TDPMjet::GetMass(Int_t evnum) const 	      
{
	return DTEVT1.phkk[evnum][4];
}
//______________________________________________________________________________
Int_t    TDPMjet::GetFragmentA(Int_t evnum) const	
{
	return DTEVT2.idres[evnum];
}
//______________________________________________________________________________
Int_t    TDPMjet::GetFragmentZ(Int_t evnum) const	
{
	return DTEVT2.idxres[evnum];
}
//______________________________________________________________________________
Double_t TDPMjet::GetXSFrac() const 	      
{
	return DTIMPA.xsfrac;
}
//______________________________________________________________________________
Double_t TDPMjet::GetBImpac() const 	      
{
	return DTGLCP.bimpac;
}
//______________________________________________________________________________
Double_t TDPMjet::GetProjRadius() const 	      
{
	return DTGLCP.rproj;
}
//______________________________________________________________________________
Double_t TDPMjet::GetTargRadius() const 	      
{
	return DTGLCP.rtarg;
}
//______________________________________________________________________________
Int_t TDPMjet::GetProjWounded() const
{
	return DTGLCP.nwasam;
}
//______________________________________________________________________________
Int_t TDPMjet::GetTargWounded() const
{
	return DTGLCP.nwbsam;
}
//______________________________________________________________________________
Int_t TDPMjet::GetProjSpectators() const
{
	return DTGLCP.nwtaac;
}
//______________________________________________________________________________
Int_t TDPMjet::GetTargSpectators() const
{
	return DTGLCP.nwtbac;
}

//______________________________________________________________________________
void TDPMjet::Dt_Rndm(int idummy)
{
	dt_rndm(idummy);
}

//______________________________________________________________________________
void TDPMjet::Dt_Rndmst(int na1, int na2, int na3, int nb1)
{
	dt_rndmst(na1, na2, na3, nb1);
}

//______________________________________________________________________________
void TDPMjet::Dt_Rndmin(int u, int c, int cd, int cm, int i, int j)
{
	dt_rndmin(u, c, cd, cm, i, j);
}

//______________________________________________________________________________
void TDPMjet::Dt_Rndmou(int u, int c, int cd, int cm, int i, int j)
{
	dt_rndmou(u, c, cd, cm, i, j);
}

