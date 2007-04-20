////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// THydjet                                                                    //
//                                                                            //
// THydjet is an interface class to fortran version of Hydjet event generator //
//                                                                            //
//      -------------------------------------------------------------         //
//      HYDJET, fast MC code to simulate flow effects, jet production         // 
//      and jet quenching in heavy ion AA collisions at the LHC               //
//      -------------------------------------------------------------         //
//      This code is merging HYDRO (flow effects), PYTHIA6.4 (hard jet        // 
//      production) and PYQUEN (jet quenching)                                //
//      --------------------------------------------------------------        //
//                                                                            //
//      Igor Lokhtin, SINP MSU, Moscow, RU                                    //
//        e-mail: Igor.Lokhtin@cern.ch                                        // 
//                                                                            //
//      Reference for HYDJET:                                                 // 
//      I.P. Lokhtin, A.M. Snigirev,                                          // 
//      Eur. Phys. J. C 46 (2006) 211.                                        //
//                                                                            //
//      References for HYDRO:                                                 // 
//      N.A.Kruglov, I.P.Lokhtin, L.I.Sarycheva, A.M.Snigirev,                // 
//      Z. Phys. C 76 (1997) 99;                                              //  
//      I.P.Lokhtin, L.I.Sarycheva, A.M.Snigirev,                             // 
//      Phys. Lett. B 537 (2002) 261;                                         //   
//    I.P.Lokhtin, A.M.Snigirev, Preprint SINP MSU 2004-14/753,hep-ph/0312204.//
//                                                                            //
//      References for PYQUEN:                                                // 
//      I.P.Lokhtin, A.M.Snigirev, Eur.Phys.J. C16 (2000) 527;                //   
//    I.P.Lokhtin, A.M.Snigirev, Preprint SINP MSU 2004-13/752, hep-ph/0406038.//
//                                                                             //
//      References for PYTHIA:                                                 //
//      T.Sjostrand et al., Comput.Phys.Commun. 135 (2001) 238;                // 
//      T.Sjostrand, S. Mrena and P. Skands, hep-ph/0603175.                   //
//                                                                             //
//      Reference for JETSET event format:                                     //
//      T.Sjostrand, Comput.Phys.Commun. 82 (1994) 74.                         //
//                                                                             // 
//      --------------------------------------------------------------         //
//      Web-page:                                                              //
//      http://cern.ch/lokhtin/hydro                                           //
//      --------------------------------------------------------------         //
//                                                                             //
//**************************************************************************** //

#include "THydjet.h"
#include "TObjArray.h"
#include "HydCommon.h"
#include "TParticle.h"
#include "TROOT.h"

#ifndef WIN32
# define pyinit pyinit_
# define hydro hydro_
# define type_of_call
#else
# define pyinit PYINIT
# define hydro HYDRO
# define type_of_call _stdcall
#endif

extern "C" void type_of_call hydro(float* A, int* ifb, float* bmin,
                                    float* bmax, float* bfix, int* nh);
//extern "C" void type_of_call luedit(Int_t& medit);
#ifndef WIN32
extern "C" void type_of_call pyinit( const char *frame, const char *beam, const char *target,
                                     double *win, Long_t l_frame, Long_t l_beam,
                                     Long_t l_target);
#else
extern "C" void type_of_call pyinit( const char *frame,  Long_t l_frame,
                                    const char *beam,   Long_t l_beam,
                                    const char *target, Long_t l_target,
                                    double *win
                                    );
#endif

#include <TClonesArray.h>

ClassImp(THydjet)

THydjet::THydjet() :
  TGenerator("Hydjet","Hydjet"),
  fEfrm(5500),
  fFrame("CMS"),
  fAw(207),
  fIfb(0),
  fBmin(0.),
  fBmax(1.),
  fBfix(0.),
  fNh(20000)
{
// Default constructor
}

//______________________________________________________________________________
THydjet::THydjet(Float_t efrm, const char *frame="CMS",
		 Float_t aw=207., Int_t ifb=0, Float_t bmin=0, Float_t bmax=1, Float_t bfix=0,
		 Int_t nh=20000) :
  TGenerator("Hydjet","Hydjet"),
  fEfrm(efrm),
  fFrame(frame),
  fAw(aw),
  fIfb(ifb),
  fBmin(bmin),
  fBmax(bmax),
  fBfix(bfix),
  fNh(nh)
{
// THydjet constructor:
}

//______________________________________________________________________________
THydjet::~THydjet()
{
// Destructor
}


TObjArray* THydjet::ImportParticles(Option_t *option)
{
//
//  Default primary creation method. It reads the /LUJETS common block which
//  has been filled by the GenerateEvent method.
//  The function loops on the generated particles and store them in
//  the TClonesArray pointed by the argument particles.
//  The default action is to store only the stable particles (LUJETS.k[0][i] == 1)
//  This can be demanded explicitly by setting the option = "Final"
//  If the option = "All", all the particles are stored.
//
    fParticles->Clear();
    Int_t numpart = LUJETS.n;
    printf("\n THydjet: Hydjet stack contains %d particles.", numpart);
    Int_t nump = 0;
    if(!strcmp(option,"") || !strcmp(option,"Final")) {
        for(Int_t i = 0; i < numpart; i++) {
            if(LUJETS.k[0][i] == 1) {
              //Use the common block values for the TParticle constructor
              nump++;
              TParticle* p = new TParticle(
              LUJETS.k[1][i], LUJETS.k[0][i] ,
              LUJETS.k[2][i], -1, LUJETS.k[3][i], LUJETS.k[4][i],
              LUJETS.p[0][i], LUJETS.p[1][i], LUJETS.p[2][i], LUJETS.p[3][i] ,
              LUJETS.v[0][i], LUJETS.v[1][i], LUJETS.v[2][i], LUJETS.v[3][i]
              );
              fParticles->Add(p);
	         }
	     }
    }
    else if(!strcmp(option,"All")) {
	       nump = numpart;
	       for(Int_t i = 0; i < numpart; i++){
              TParticle* p = new TParticle(
              LUJETS.k[1][i], LUJETS.k[0][i] ,
              LUJETS.k[2][i], -1, LUJETS.k[3][i], LUJETS.k[4][i],
              LUJETS.p[0][i], LUJETS.p[1][i], LUJETS.p[2][i], LUJETS.p[3][i] ,
              LUJETS.v[0][i], LUJETS.v[1][i], LUJETS.v[2][i], LUJETS.v[3][i]
              );
              fParticles->Add(p);
          }
    }
    return fParticles;
}

Int_t THydjet::ImportParticles(TClonesArray *particles, Option_t *option)
{
//
//  Default primary creation method. It reads the /LUJETS common block which
//  has been filled by the GenerateEvent method.
//  The function loops on the generated particles and store them in
//  the TClonesArray pointed by the argument particles.
//  The default action is to store only the stable particles (LUJETS.k[0][i] == 1)
//  This can be demanded explicitly by setting the option = "Final"
//  If the option = "All", all the particles are stored.
//
  if (particles == 0) return 0;
  TClonesArray &particlesR = *particles;
  particlesR.Clear();
  Int_t numpart = LUJETS.n;
  printf("\n THydjet: Hydjet stack contains %d particles.", numpart);
  Int_t nump = 0;
  if(!strcmp(option,"") || !strcmp(option,"Final")) {
        for(Int_t i = 0; i < numpart; i++) {
            if(LUJETS.k[0][i] == 1) {
              //Use the common block values for the TParticle constructor
              nump++;
              new(particlesR[i]) TParticle(
              LUJETS.k[1][i], LUJETS.k[0][i] ,
              LUJETS.k[2][i], -1, LUJETS.k[3][i], LUJETS.k[4][i],
              LUJETS.p[0][i], LUJETS.p[1][i], LUJETS.p[2][i], LUJETS.p[3][i] ,
              LUJETS.v[0][i], LUJETS.v[1][i], LUJETS.v[2][i], LUJETS.v[3][i]
              );
	         }
	     }
    }
    else if(!strcmp(option,"All")){
	       nump = numpart;
	       for(Int_t i = 0; i < numpart; i++){
              new(particlesR[i]) TParticle(
              LUJETS.k[1][i], LUJETS.k[0][i] ,
              LUJETS.k[2][i], -1, LUJETS.k[3][i], LUJETS.k[4][i],
              LUJETS.p[0][i], LUJETS.p[1][i], LUJETS.p[2][i], LUJETS.p[3][i] ,
              LUJETS.v[0][i], LUJETS.v[1][i], LUJETS.v[2][i], LUJETS.v[3][i]
              );
          }
    }
    return nump;
}

//______________________________________________________________________________
void THydjet::SetEfrm(Float_t efrm)
{
// Set the centre of mass (CMS) or lab-energy (LAB)
   fEfrm=efrm;
}
//______________________________________________________________________________
void THydjet::SetFrame(const char* frame)
{
// Set the frame type ("CMS" or "LAB")
   fFrame=frame;
}
//______________________________________________________________________________
/*void THydjet::SetProj(const char* proj)
{
// Set the projectile type
   fProj=proj;
}
//______________________________________________________________________________
void THydjet::SetTarg(const char* targ)
{
// Set the target type
   fTarg=targ;
}
*/
//______________________________________________________________________________
void THydjet::SetAw(Float_t aw)
{
// Set the projectile-targed atomic number
   fAw=aw;
}
//______________________________________________________________________________
void THydjet::SetIfb(Int_t ifb)
{
// flag of type of centrality generation
   fIfb=ifb;
} 
//______________________________________________________________________________
void THydjet::SetBmin(Float_t bmin)
{
// set minimum impact parameter in units of nucleus radius RA
   fBmin=bmin;
}
//______________________________________________________________________________
void THydjet::SetBmax(Float_t bmax)
{
// set maximum impact parameter in units of nucleus radius RA
   fBmax=bmax;
} 
//______________________________________________________________________________
void THydjet::SetBfix(Float_t bfix)
{
// Set fixed impact parameter in units of nucleus radius RA
   fBfix=bfix;
} 
//______________________________________________________________________________
void THydjet::SetNh(Int_t nh)
{
// Set mean soft hadron multiplicity in central Pb-Pb collisions
   fNh=nh;
} 
//______________________________________________________________________________
Float_t THydjet::GetEfrm() const
{
// Get the centre of mass (CMS) or lab-energy (LAB)
   return fEfrm;
}
//______________________________________________________________________________
const char* THydjet::GetFrame() const
{
// Get the frame type ("CMS" or "LAB")
   return fFrame.Data();
}
//______________________________________________________________________________
/*const char* THydjet::GetProj() const
{
// Get the projectile type
   return fProj;
}
//______________________________________________________________________________
const char* THydjet::GetTarg() const
{
// Set the target type
   return fTarg;
}
*/
//______________________________________________________________________________
Float_t THydjet::GetAw() const
{
// Get the projectile atomic number
   return fAw;
}
//______________________________________________________________________________
Int_t THydjet::GetIfb() const
{
// Get flag of type of centrality generation
   return fIfb;
}
//______________________________________________________________________________
Float_t THydjet::GetBmin() const
{
// Get minimum impact parameter in units of nucleus radius RA
   return fBmin;
} 
//______________________________________________________________________________
Float_t THydjet::GetBmax() const
{
// Get maximum impact parameter in units of nucleus radius RA
   return fBmax;
}
//______________________________________________________________________________
Float_t THydjet::GetBfix() const
{
// Get fixed impact parameter in units of nucleus radius RA
   return fBfix;
}
//______________________________________________________________________________
Int_t THydjet::GetNh() const
{
// Get mean soft hadron multiplicity in central Pb-Pb collisions
   return fNh;
} 

//====================== access to common HYFLOW ===============================

//______________________________________________________________________________
const void THydjet::SetYTFL(Float_t value) const
{
   HYFLOW.ytfl=value;
}

//______________________________________________________________________________
Float_t THydjet::GetYTFL() const
{
   return HYFLOW.ytfl;
}

//______________________________________________________________________________
const void THydjet::SetYLFL(Float_t value) const
{
   HYFLOW.ylfl=value;
}

//______________________________________________________________________________
Float_t THydjet::GetYLFL() const
{
   return HYFLOW.ylfl;
}

//______________________________________________________________________________
const void THydjet::SetFPART(Float_t value) const
{
   HYFLOW.fpart=value;
}


//______________________________________________________________________________
Float_t THydjet::GetFPART() const
{
   return HYFLOW.fpart;
}


//====================== access to common HYJPAR ===============================

//______________________________________________________________________________
const void THydjet::SetNHSEL(Int_t value) const
{
   HYJPAR.nhsel=value;
}

//______________________________________________________________________________
Int_t THydjet::GetNHSEL() const
{
   return HYJPAR.nhsel;
}

//______________________________________________________________________________
const void THydjet::SetPTMIN(Float_t value) const
{
   HYJPAR.ptmin=value;
}

//______________________________________________________________________________
Float_t THydjet::GetPTMIN() const
{
   return HYJPAR.ptmin;
}

//______________________________________________________________________________
const void THydjet::SetNJET(Int_t value) const
{
  HYJPAR.njet=value;
}

//______________________________________________________________________________
Int_t  THydjet::GetNJET() const
{
   return HYJPAR.njet;
}

//====================== access to common HYFPAR ===============================

//______________________________________________________________________________
Float_t THydjet::GetBGEN() const
{
   return HYFPAR.bgen;
}

//______________________________________________________________________________
Int_t THydjet::GetNBCOL() const
{
   return HYFPAR.nbcol;
}

//______________________________________________________________________________
Int_t THydjet::GetNPART() const
{
   return HYFPAR.npart;
}

//______________________________________________________________________________
Int_t THydjet::GetNPYT() const
{
   return HYFPAR.npyt;
}

//______________________________________________________________________________
Int_t THydjet::GetNHYD() const
{
   return HYFPAR.nhyd;
}


//====================== access to common LUJETS ===============================

//______________________________________________________________________________
Int_t THydjet::GetN() const
{
   return LUJETS.n;
}

//______________________________________________________________________________
Int_t THydjet::GetK(Int_t key1, Int_t key2) const
{
  // Get Particle codes information
   if ( key1<1 || key1>150000 ) {
      printf("ERROR in THydjet::GetK(key1,key2):\n");
      printf("      key1=%i is out of range [1..150000]\n",key1);
      return 0;
   }

   if ( key2<1 || key2>5 ) {
      printf("ERROR in THydjet::GetK(key1,key2):\n");
      printf("      key2=%i is out of range [1..5]\n",key2);
      return 0;
   }

   return   LUJETS.k[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THydjet::GetP(Int_t key1, Int_t key2) const
{
   // Get Particle four momentum and mass
   if ( key1<1 || key1>150000 ) {
      printf("ERROR in THydjet::GetP(key1,key2):\n");
      printf("      key1=%i is out of range [1..150000]\n",key1);
      return 0;
   }

   if ( key2<1 || key2>5 ) {
      printf("ERROR in THydjet::GetP(key1,key2):\n");
      printf("      key2=%i is out of range [1..5]\n",key2);
      return 0;
   }

   return   LUJETS.p[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THydjet::GetV(Int_t key1, Int_t key2) const
{
   // Get particle vertex, production time and lifetime
   if ( key1<1 || key1>150000 ) {
      printf("ERROR in THydjet::GetV(key1,key2):\n");
      printf("      key1=%i is out of range [1..150000]\n",key1);
      return 0;
   }

   if ( key2<1 || key2>5 ) {
      printf("ERROR in THydjet::GetV(key1,key2):\n");
      printf("      key2=%i is out of range [1..5]\n",key2);
      return 0;
   }

   return   LUJETS.v[key2-1][key1-1];
}

//====================== access to common HYJETS ===============================

//______________________________________________________________________________
Int_t THydjet::GetNL() const
{
   return HYJETS.nl;
}

//______________________________________________________________________________
Int_t THydjet::GetKL(Int_t key1, Int_t key2) const
{
   // Get Particle codes information
   if ( key1<1 || key1>150000 ) {
      printf("ERROR in THydjet::GetKL(key1,key2):\n");
      printf("      key1=%i is out of range [1..150000]\n",key1);
      return 0;
   }

   if ( key2<1 || key2>5 ) {
      printf("ERROR in THydjet::GetKL(key1,key2):\n");
      printf("      key2=%i is out of range [1..5]\n",key2);
      return 0;
   }

   return   HYJETS.kl[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THydjet::GetPL(Int_t key1, Int_t key2) const
{
   // Get Particle four momentum and mass
   if ( key1<1 || key1>150000 ) {
      printf("ERROR in THydjet::GetPL(key1,key2):\n");
      printf("      key1=%i is out of range [1..150000]\n",key1);
      return 0;
   }

   if ( key2<1 || key2>5 ) {
      printf("ERROR in THydjet::GetPL(key1,key2):\n");
      printf("      key2=%i is out of range [1..5]\n",key2);
      return 0;
   }

   return   HYJETS.pl[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THydjet::GetVL(Int_t key1, Int_t key2) const
{
   // Get particle vertex, production time and lifetime
   if ( key1<1 || key1>150000 ) {
      printf("ERROR in THydjet::GetVL(key1,key2):\n");
      printf("      key1=%i is out of range [1..150000]\n",key1);
      return 0;
   }

   if ( key2<1 || key2>5 ) {
      printf("ERROR in THydjet::GetVL(key1,key2):\n");
      printf("      key2=%i is out of range [1..5]\n",key2);
      return 0;
   }

   return   HYJETS.vl[key2-1][key1-1];
}


//====================== access to common PYDAT1 ===============================

//______________________________________________________________________________
void THydjet::SetMSTU(Int_t key, Int_t value)
{
   //Set MSTU in Pythia 
   if ( key<1 || key>200 ) {
      printf("ERROR in THydjet::SetMSTU(key,value):\n");
      printf("      key=%i is out of range [1..200]\n",key);
   }
   PYDAT1.mstu[key-1] = value;
}

//______________________________________________________________________________
void THydjet::SetPARU(Int_t key, Double_t value)
{
   //Set PARU in Pythia 
   if ( key<1 || key>200 ) {
      printf("ERROR in THydjet::SetPARU(key,value):\n");
      printf("      key=%i is out of range [1..200]\n",key);
   }
   PYDAT1.paru[key-1] = value;
}

//______________________________________________________________________________
void THydjet::SetMSTJ(Int_t key, Int_t value)
{
   //Set MSTJ in Pythia 
   if ( key<1 || key>200 ) {
      printf("ERROR in THydjet::SetMSTJ(key,value):\n");
      printf("      key=%i is out of range [1..200]\n",key);
   }
   PYDAT1.mstj[key-1] = value;
}

//______________________________________________________________________________
void THydjet::SetPARJ(Int_t key, Double_t value)
{
   //Set PARJ in Pythia 
   if ( key<1 || key>200 ) {
      printf("ERROR in THydjet::SetPARJ(key,value):\n");
      printf("      key=%i is out of range [1..200]\n",key);
   }
   PYDAT1.parj[key-1] = value;
}


//====================== access to common PYSUBS ===============================

//______________________________________________________________________________
const void THydjet::SetMSEL(Int_t value) const
{
  PYSUBS.msel=value;
}

//______________________________________________________________________________
void THydjet::SetCKIN(Int_t key, Double_t value)
{
   //Set CKIN in Pythia 
   if ( key<1 || key>200 ) {
      printf("ERROR in THydjet::SetCKIN(key,value):\n");
      printf("      key=%i is out of range [1..200]\n",key);
   }
   PYSUBS.ckin[key-1] = value;
}

//====================== access to common PYPARS ===============================

//______________________________________________________________________________
void THydjet::SetMSTP(Int_t key, Int_t value)
{
   //Set MSTP in Pythia 
   if ( key<1 || key>200 ) {
      printf("ERROR in THydjet::SetMSTP(key,value):\n");
      printf("      key=%i is out of range [1..200]\n",key);
   }
   PYPARS.mstp[key-1] = value;
}


//====================== access to Hijing subroutines =========================


//______________________________________________________________________________
void THydjet::Initialize()
{

   // Initialize PYTHIA for hard parton-parton scattering
   if ( (!strcmp(fFrame.Data(), "CMS     "  )) &&
        (!strcmp(fFrame.Data(), "LAB     "  ))){
      printf("WARNING! In THydjet:Initialize():\n");
      printf(" specified frame=%s is neither CMS or LAB\n",fFrame.Data());
      printf(" resetting to default \"CMS\" .");
      fFrame="CMS";
   }
   Int_t nhselflag = GetNHSEL();
   if(nhselflag != 0) {
      Double_t lwin = fEfrm;
      Long_t  s1    = strlen(fFrame);
      Long_t  s2    = strlen("p");
      Long_t  s3    = strlen("p");
#ifndef WIN32
     pyinit(fFrame,"p","p",&lwin,s1,s2,s3);
#else
     pyinit(fFrame, s1, "p" , s2, "p", s3, &lwin);
#endif
   }
}


//______________________________________________________________________________
void THydjet::GenerateEvent()
{
// Generates one event;
  float xbmin = fBmin;
  float xbmax = fBmax;
  float xbfix = fBfix;
  float xAw   = fAw;
  hydro(&xAw,&fIfb,&xbmin,&xbmax,&xbfix,&fNh);

}
//______________________________________________________________________________
void THydjet::Hydro()
{
  // Generates one event;
  float xbmin = fBmin;
  float xbmax = fBmax;
  float xbfix = fBfix;
  float xAw   = fAw;
  hydro(&xAw,&fIfb,&xbmin,&xbmax,&xbfix,&fNh);
}
