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

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// THijing                                                                    //
//                                                                            //
// THijing is an interface class to F77 version of Hijing 1.36                //
// event generator, written by X.N. Wang and  M. Gyulassy.                    //
// For details see http://nta2.lbl.gov/~xnwang                                //
//                                                                            //
//          **************************************************                //
//	    *	  |		_______      /  ------/      *		      //
//	    *	----- ------	 |_____|     /_/     /       *		      //
//	    *	 ||    /	|_____|      /    /	     *		      //
//	    *	 /|   /_/	/_______    /_  /    _       *		      //
//	    *	/ |	/ /	/  /  / |	 -------     *		      //
//	    *	  |    / /	 /  /  |     /     |	     *		      //
//	    *	  |   / /	/  / _|    /   -------       *                //
//	    *						     *		      //
//	    **************************************************		      //
//				  HIJING				      //
//		   Heavy Ion Jet INteraction Generator			      //
//				    by					      //
//		       X. N. Wang  and  M. Gyulassy			      //
//		       Lawrence Berkeley Laboratory                           //
//****************************************************************************//


#include <TClonesArray.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TROOT.h>

#include "Hcommon.h"
#include "THijing.h"

#ifndef WIN32
# define hijset hijset_
# define hijing hijing_
# define profile profile_
# define rluget_hijing rluget_hijing_
# define rluset_hijing rluset_hijing_
# define type_of_call
#else
# define hijset HIJSET
# define hijing HIJING
# define profile PROFILE
# define rluget_hijing RLUGET_HIJING
# define rluset_hijing RLUSET_HIJING
# define type_of_call _stdcall
#endif

#ifndef WIN32
extern "C" void type_of_call hijset(Float_t & , const char *, 
				   const char *, const char *,
				   Int_t & , Int_t &, Int_t &,
				   Int_t &, const int, 
				   const int, const int);
extern "C" float type_of_call profile(Float_t &);
extern "C" void type_of_call hijing(const char *, Float_t  &,
                                   Float_t &, const int);
extern "C" void type_of_call rluget_hijing(Int_t & lfn, Int_t & move);

extern "C" void type_of_call rluset_hijing(Int_t & lfn, Int_t & move);

#else
#endif



ClassImp(THijing)


THijing::THijing(): 
    TGenerator("Hijing","Hijing"),
    fEfrm(5500.),
    fFrame("CMS"),
    fProj("A"),
    fTarg("A"),
    fIap(208),
    fIzp(82),
    fIat(208),
    fIzt(82),
    fBmin(0.),
    fBmax(5.)
{
// Default constructor 
}

//______________________________________________________________________________
THijing::THijing(Float_t efrm, const char *frame="CMS", 
		 const char *proj="A", const char *targ="A", 
		 Int_t iap=207, Int_t izp=82, Int_t iat=207, Int_t izt=82,
		 Float_t bmin=0, Float_t bmax=20): 
    TGenerator("Hijing","Hijing"),
    fEfrm(efrm),
    fFrame(frame),
    fProj(proj),
    fTarg(targ),
    fIap(iap),
    fIzp(izp),
    fIat(iat),
    fIzt(izt),
    fBmin(bmin),
    fBmax(bmax)
{
// THijing constructor: 
// Note that there may be only one functional THijing object
// at a time, so it's not use to create more than one instance of it.
}

//______________________________________________________________________________
THijing::~THijing()
{
// Destructor
}


TObjArray* THijing::ImportParticles(Option_t *option)
{
//
//  Default primary creation method. It reads the /HEPEVT/ common block which
//  has been filled by the GenerateEvent method. If the event generator does
//  not use the HEPEVT common block, This routine has to be overloaded by
//  the subclasses.
//  The function loops on the generated particles and store them in
//  the TClonesArray pointed by the argument particles.
//  The default action is to store only the stable particles (ISTHEP = 1)
//  This can be demanded explicitly by setting the option = "Final"
//  If the option = "All", all the particles are stored.
//
    fParticles->Clear();
    Int_t numpart = HIMAIN1.natt;
    printf("\n THijing: HIJING stack contains %d particles.", numpart);
    printf("\n THijing: Total energy:         %f           ", HIMAIN1.eatt);
    printf("\n THijing: Number of hard scatterings: %d     ", HIMAIN1.jatt);
    Int_t nump = 0;
    if (!strcmp(option,"") || !strcmp(option,"Final")) {
	for (Int_t i = 0; i < numpart; i++) {
	  
	    if (HIMAIN2.katt[3][i] == 1) {
//
//  Use the common block values for the TParticle constructor
//
		nump++;
		TParticle* p = new TParticle(
		    HIMAIN2.katt[0][i], HIMAIN2.katt[3][i] ,
		    -1, -1, -1, -1,
		    HIMAIN2.patt[0][i], HIMAIN2.patt[1][i], HIMAIN2.patt[2][i], HIMAIN2.patt[3][i] ,
		    HIMAIN2.vatt[0][i], HIMAIN2.vatt[1][i], HIMAIN2.vatt[2][i], HIMAIN2.vatt[3][i] 
		  );
		p->SetUniqueID(HIMAIN2.katt[1][i]);
		fParticles->Add(p);
	    }
	}
    }
    else if (!strcmp(option,"All")) {
	nump = numpart; 
	for (Int_t i = 0; i < numpart; i++) {
	    
	    Int_t iParent = HIMAIN2.katt[2][i]-1;
	    
	    if (iParent >= 0) {
		TParticle *mother = (TParticle*) (fParticles->UncheckedAt(iParent));	   
		mother->SetLastDaughter(i);
		if (mother->GetFirstDaughter()==-1)
		    mother->SetFirstDaughter(i);
	    }
	    
	    TParticle* p = new TParticle(
		HIMAIN2.katt[0][i], HIMAIN2.katt[3][i], iParent,
		-1, -1, -1,
		HIMAIN2.patt[0][i], HIMAIN2.patt[1][i], HIMAIN2.patt[2][i], HIMAIN2.patt[3][i] ,
		HIMAIN2.vatt[0][i], HIMAIN2.vatt[1][i], HIMAIN2.vatt[2][i], HIMAIN2.vatt[3][i]
		);
	    p->SetUniqueID(HIMAIN2.katt[1][i]);
	    fParticles->Add(p);
	}
    }
    return fParticles;
}

Int_t THijing::ImportParticles(TClonesArray *particles, Option_t *option)
{
//
//  Default primary creation method. It reads the /HEPEVT/ common block which
//  has been filled by the GenerateEvent method. If the event generator does
//  not use the HEPEVT common block, This routine has to be overloaded by
//  the subclasses.
//  The function loops on the generated particles and store them in
//  the TClonesArray pointed by the argument particles.
//  The default action is to store only the stable particles (ISTHEP = 1)
//  This can be demanded explicitly by setting the option = "Final"
//  If the option = "All", all the particles are stored.
//
  if (particles == 0) return 0;
  TClonesArray &particlesR = *particles;
  particlesR.Clear();
  Int_t numpart = HIMAIN1.natt;
  printf("\n THijing: HIJING stack contains %d particles.", numpart);
  printf("\n THijing: Total energy:         %f           ", HIMAIN1.eatt);
  printf("\n THijing: Number of hard scatterings: %d     ", HIMAIN1.jatt);
  Int_t nump = 0;
  if (!strcmp(option,"") || !strcmp(option,"Final")) {
      for (Int_t i = 0; i < numpart; i++) {
	  
	  if (HIMAIN2.katt[3][i] == 1) {
//
//  Use the common block values for the TParticle constructor
//
	    nump++;
	    new(particlesR[i]) TParticle(
		  HIMAIN2.katt[0][i] ,
		  HIMAIN2.katt[3][i] ,
		  -1 ,
		  -1,
		  -1,
		  -1,
		  
		  HIMAIN2.patt[0][i] ,
		  HIMAIN2.patt[1][i] ,
		  HIMAIN2.patt[2][i] ,
		  HIMAIN2.patt[3][i] ,

		  HIMAIN2.vatt[0][i] ,
		  HIMAIN2.vatt[1][i] ,
		  HIMAIN2.vatt[2][i] ,
		  HIMAIN2.vatt[3][i] 
		  );
	    particlesR[i]->SetUniqueID(HIMAIN2.katt[1][i]);
	  }
      }
  }
  else if (!strcmp(option,"All")) {
      nump = numpart; 
      for (Int_t i = 0; i < numpart; i++) {

	  Int_t iParent = HIMAIN2.katt[2][i]-1;
	  
	  if (iParent >= 0) {
	      TParticle *mother = (TParticle*) (particlesR.UncheckedAt(iParent));	   
	      mother->SetLastDaughter(i);
	      if (mother->GetFirstDaughter()==-1)
		  mother->SetFirstDaughter(i);
	  }

	  new(particlesR[i]) TParticle(
	      HIMAIN2.katt[0][i] ,
	      HIMAIN2.katt[3][i] ,
	      iParent,
	      -1,
	      -1,
	      -1,
	      
	      HIMAIN2.patt[0][i] ,
	      HIMAIN2.patt[1][i] ,
	      HIMAIN2.patt[2][i] ,
	      HIMAIN2.patt[3][i] ,
	      
	      HIMAIN2.vatt[0][i] ,
	      HIMAIN2.vatt[1][i] ,
	      HIMAIN2.vatt[2][i] ,
	      HIMAIN2.vatt[3][i]
	      );
	  particlesR[i]->SetUniqueID(HIMAIN2.katt[1][i]);
      }
  }
  return nump;
}

//______________________________________________________________________________
void THijing::SetEFRM(Float_t efrm)
{
// Set the centre of mass (CMS) or lab-energy (LAB)
   fEfrm=efrm;
} 
//______________________________________________________________________________
void THijing::SetFRAME(const char* frame)
{
// Set the frame type ("CMS" or "LAB")
   fFrame=frame;
} 
//______________________________________________________________________________
void THijing::SetPROJ(const char* proj)
{
// Set the projectile type
   fProj=proj;
} 
//______________________________________________________________________________
void THijing::SetTARG(const char* targ)
{
// Set the target type
   fTarg=targ;
} 
//______________________________________________________________________________
void THijing::SetIAP(Int_t iap)
{
// Set the projectile atomic number
   fIap=iap;
} 
//______________________________________________________________________________
void THijing::SetIZP(Int_t izp)
{
// Set the projectile charge number
   fIzp=izp;
} 
//______________________________________________________________________________
void THijing::SetIAT(Int_t iat)
{
// Set the target atomic number
   fIat=iat;
} 
//______________________________________________________________________________
void THijing::SetIZT(Int_t izt)
{
// Set the target charge number
   fIzt=izt;
} 
//______________________________________________________________________________
void THijing::SetBMIN(Float_t bmin)
{
// Set the minimum impact parameter
   fBmin=bmin;
} 
//______________________________________________________________________________
void THijing::SetBMAX(Float_t bmax)
{
// Set the maximum impact parameter
   fBmax=bmax;
} 
//______________________________________________________________________________
Float_t THijing::GetEFRM() const
{
// Get the centre of mass (CMS) or lab-energy (LAB)
   return fEfrm;
} 
//______________________________________________________________________________
const char* THijing::GetFRAME() const
{
// Get the frame type ("CMS" or "LAB")
   return fFrame.Data();
} 
//______________________________________________________________________________
const char* THijing::GetPROJ() const
{
// Get the projectile type
   return fProj;
} 
//______________________________________________________________________________
const char* THijing::GetTARG() const
{
// Set the target type
   return fTarg;
} 
//______________________________________________________________________________
Int_t THijing::GetIAP() const
{
// Get the projectile atomic number
   return fIap;
} 
//______________________________________________________________________________
Int_t THijing::GetIZP() const
{
// Get the projectile charge number
   return fIzp;
} 
//______________________________________________________________________________
Int_t THijing::GetIAT() const
{
// Get the target atomic number
   return fIat;
} 
//______________________________________________________________________________
Int_t THijing::GetIZT() const
{
// Get the target charge number
   return fIzt;
} 
//______________________________________________________________________________
Float_t THijing::GetBMIN() const
{
// Get the minimum impact parameter
   return fBmin;
} 
//______________________________________________________________________________
Float_t THijing::GetBMAX() const
{
// Get the maximum impact parameter
   return fBmax;
} 

//====================== access to common HIPARNT ===============================

//______________________________________________________________________________
void THijing::SetHIPR1(Int_t key,Float_t value)
{
// Set the values of array HIPR1 in common HIPARNT
   if ( key<1 || key>100 ) {
      printf ("ERROR in THijing:SetHIPR1(key,value): \n ");
      printf ("      key=%i is out of range [1..100]!\n",key);
      return;
   }

   HIPARNT.hipr1[key-1]=value;

}

//______________________________________________________________________________
Float_t THijing::GetHIPR1(Int_t key) const
{
// Get the values of array HIPR1 in common HIPARNT
   if ( key<1 || key>100 ) {
      printf ("ERROR in THijing:GetHIPR1(key): \n ");
      printf ("      key=%i is out of range [1..100]!\n",key);
      return 0;
   }

   return HIPARNT.hipr1[key-1];

}

//______________________________________________________________________________
void THijing::SetIHPR2(Int_t key,Int_t value)
{
// Set the values of array HIPR2 in common HIPARNT
   if ( key<1 || key>50 ) {
      printf ("ERROR in THijing:SetIHPR2(key,value): \n ");
      printf ("      key=%i is out of range [1..50]!\n",key);
      return;
   }

   HIPARNT.ihpr2[key-1]=value;

}

//______________________________________________________________________________
Int_t THijing::GetIHPR2(Int_t key) const
{
// Get the values of array HIPR2 in common HIPARNT
   if ( key<1 || key>50 ) {
      printf ("ERROR in THijing:GetIHPR2(key): \n ");
      printf ("      key=%i is out of range [1..50]!\n",key);
      return 0;
   }

   return HIPARNT.ihpr2[key-1];

}


//______________________________________________________________________________
Float_t THijing::GetHINT1(Int_t key) const
{
// Get the values of array HINT1 in common HIPARNT
   if ( key<1 || key>100 ) {
      printf ("ERROR in THijing:GetHINT1(key): \n ");
      printf ("      key=%i is out of range [1..100]!\n",key);
      return 0;
   }

   return HIPARNT.hint1[key-1];

}


//______________________________________________________________________________
Int_t THijing::GetIHNT2(Int_t key) const
{
// Get the values of array HINT2 in common HIPARNT
   if ( key<1 || key>50 ) {
      printf ("ERROR in THijing:GetIHNT2(key): \n ");
      printf ("      key=%i is out of range [1..50]!\n",key);
      return 0;
   }

   return HIPARNT.ihnt2[key-1];

}


//====================== access to common HIMAIN1 ===============================

//______________________________________________________________________________
Int_t  THijing::GetNATT() const
{
// Get the number of particles produces
   return HIMAIN1.natt;

}

//______________________________________________________________________________
Float_t  THijing::GetEATT() const
{
// Get total energy of particles

   return HIMAIN1.eatt;

}

//______________________________________________________________________________
Int_t  THijing::GetJATT() const
{
// Get number of hard scatterings

   return HIMAIN1.jatt;

}

//______________________________________________________________________________
Int_t  THijing::GetNT() const
{
// Get number of target participants

   return HIMAIN1.nt;

}

//______________________________________________________________________________
Int_t  THijing::GetNP() const
{
// Get number of projectile participants
   return HIMAIN1.np;

}


//______________________________________________________________________________
Int_t  THijing::GetN0() const
{
// Get number of N-N collisions
   return HIMAIN1.n0;

}
//______________________________________________________________________________
Int_t  THijing::GetN01() const
{
// Get number of N-wounded collisions

   return HIMAIN1.n01;

}

//______________________________________________________________________________
Int_t  THijing::GetN10() const
{
// Get number of wounded-N collisions

   return HIMAIN1.n10;

}

//______________________________________________________________________________
Int_t  THijing::GetN11() const
{
// Get number of wounded-wounded collisions

   return HIMAIN1.n11;

}

//______________________________________________________________________________
Float_t  THijing::GetBB() const
{
// Get impact parameter

   return HIMAIN1.bb;

}

//====================== access to common HIMAIN2 ===============================

//______________________________________________________________________________
Int_t THijing::GetKATT(Int_t key1, Int_t key2) const
{
// Get values of array KATT in common HIMAIN2
   if ( key1<1 || key1>200000 ) {
      printf("ERROR in THijing::GetKATT(key1,key2):\n");
      printf("      key1=%i is out of range [1..200000]\n",key1);
      return 0;
   }

   if ( key2<1 || key2>4 ) {
      printf("ERROR in THijing::GetKATT(key1,key2):\n");
      printf("      key2=%i is out of range [1..4]\n",key2);
      return 0;
   }
   
   return   HIMAIN2.katt[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPATT(Int_t key1, Int_t key2) const
{
// Get values of array PATT in common HIMAIN2
   if ( key1<1 || key1>200000 ) {
      printf("ERROR in THijing::GetPATT(key1,key2):\n");
      printf("      key1=%i is out of range [1..130000]\n",key1);
      return 0;
   }

   if ( key2<1 || key2>4 ) {
      printf("ERROR in THijing::GetPATT(key1,key2):\n");
      printf("      key2=%i is out of range [1..4]\n",key2);
      return 0;
   }
   
   return   HIMAIN2.patt[key2-1][key1-1];
}

Float_t THijing::GetVATT(Int_t key1, Int_t key2) const
{
// Get values of array VATT in common HIMAIN2
   if ( key1<1 || key1>200000 ) {
      printf("ERROR in THijing::GetVATT(key1,key2):\n");
      printf("      key1=%i is out of range [1..130000]\n",key1);
      return 0;
   }

   if ( key2<1 || key2>4 ) {
      printf("ERROR in THijing::GetVATT(key1,key2):\n");
      printf("      key2=%i is out of range [1..4]\n",key2);
      return 0;
   }
   
   return   HIMAIN2.vatt[key2-1][key1-1];
}

//====================== access to common HIJJET1 ===============================

//______________________________________________________________________________
Int_t THijing::GetNPJ(Int_t key) const
{
// Get values of array NPJ of common HIJJET1
   if ( key<1 || key>300 ) {
      printf("ERROR in THijing::GetNPJ(key):\n");
      printf("      key=%i is out of range [1..300]\n",key);
      return 0;
   }
   return HIJJET1.npj[key-1];
}

//______________________________________________________________________________
Int_t THijing::GetKFPJ(Int_t key1, Int_t key2) const
{
// Get values of array KFPJ in common HIJJET1
   if ( key1<1 || key1>300 ) {
      printf("ERROR in THijing::GetKFPJ(key1):\n");
      printf("      key1=%i is out of range [1..300]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>500 ) {
      printf("ERROR in THijing::GetKFPJ(key1,key2):\n");
      printf("      key2=%i is out of range [1..500]\n",key2);
      return 0;
   }
   
   return HIJJET1.kfpj[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPJPX(Int_t key1, Int_t key2) const
{
// Get values of array PJPX in common HIJJET1
   if ( key1<1 || key1>300 ) {
      printf("ERROR in THijing::GetPJPX(key1):\n");
      printf("      key1=%i is out of range [1..300]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>500 ) {
      printf("ERROR in THijing::GetPJPX(key1,key2):\n");
      printf("      key2=%i is out of range [1..500]\n",key2);
      return 0;
   }
   
   return HIJJET1.pjpx[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPJPY(Int_t key1, Int_t key2) const
{
// Get values of array PJPY in common HIJJET1
   if ( key1<1 || key1>300 ) {
      printf("ERROR in THijing::GetPJPY(key1):\n");
      printf("      key1=%i is out of range [1..300]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>500 ) {
      printf("ERROR in THijing::GetPJPY(key1,key2):\n");
      printf("      key2=%i is out of range [1..500]\n",key2);
      return 0;
   }
   
   return HIJJET1.pjpy[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPJPZ(Int_t key1, Int_t key2) const
{
// Get values of array PJPZ in common HIJJET1
   if ( key1<1 || key1>300 ) {
      printf("ERROR in THijing::GetPJPZ(key1):\n");
      printf("      key1=%i is out of range [1..300]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>500 ) {
      printf("ERROR in THijing::GetPJPZ(key1,key2):\n");
      printf("      key2=%i is out of range [1..500]\n",key2);
      return 0;
   }
   
   return HIJJET1.pjpz[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPJPE(Int_t key1, Int_t key2) const
{
// Get values of array PJPE in common HIJJET1
   if ( key1<1 || key1>300 ) {
      printf("ERROR in THijing::GetPJPE(key1):\n");
      printf("      key1=%i is out of range [1..300]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>500 ) {
      printf("ERROR in THijing::GetPJPE(key1,key2):\n");
      printf("      key2=%i is out of range [1..500]\n",key2);
      return 0;
   }
   
   return HIJJET1.pjpe[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPJPM(Int_t key1, Int_t key2) const
{
// Get values of array PJPM in common HIJJET1
   if ( key1<1 || key1>300 ) {
      printf("ERROR in THijing::GetPJPM(key1):\n");
      printf("      key1=%i is out of range [1..300]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>500 ) {
      printf("ERROR in THijing::GetPJPM(key1,key2):\n");
      printf("      key2=%i is out of range [1..500]\n",key2);
      return 0;
   }
   
   return HIJJET1.pjpm[key2-1][key1-1];
}

//______________________________________________________________________________
Int_t THijing::GetNTJ(Int_t key) const
{
// Get values of array NTJ in common HIJJET1
   if ( key<1 || key>300 ) {
      printf("ERROR in THijing::GetNTJ(key):\n");
      printf("      key=%i is out of range [1..300]\n",key);
      return 0;
   }
   return HIJJET1.ntj[key-1];
}

//______________________________________________________________________________
Int_t THijing::GetKFTJ(Int_t key1, Int_t key2) const
{
// Get values of array KFTJ in common HIJJET1
   if ( key1<1 || key1>300 ) {
      printf("ERROR in THijing::GetKFTJ(key1):\n");
      printf("      key1=%i is out of range [1..300]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>500 ) {
      printf("ERROR in THijing::GetKFTJ(key1,key2):\n");
      printf("      key2=%i is out of range [1..500]\n",key2);
      return 0;
   }
   
   return HIJJET1.kftj[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPJTX(Int_t key1, Int_t key2) const
{
// Get values of array PJTX in common HIJJET1
   if ( key1<1 || key1>300 ) {
      printf("ERROR in THijing::GetPJTX(key1):\n");
      printf("      key1=%i is out of range [1..300]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>500 ) {
      printf("ERROR in THijing::GetPJTX(key1,key2):\n");
      printf("      key2=%i is out of range [1..500]\n",key2);
      return 0;
   }
   
   return HIJJET1.pjtx[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPJTY(Int_t key1, Int_t key2) const
{
// Get values of array PJTY in common HIJJET1
   if ( key1<1 || key1>300 ) {
      printf("ERROR in THijing::GetPJTY(key1):\n");
      printf("      key1=%i is out of range [1..300]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>500 ) {
      printf("ERROR in THijing::GetPJTY(key1,key2):\n");
      printf("      key2=%i is out of range [1..500]\n",key2);
      return 0;
   }
   
   return HIJJET1.pjty[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPJTZ(Int_t key1, Int_t key2) const
{
// Get values of array PJTZ in common HIJJET1
   if ( key1<1 || key1>300 ) {
      printf("ERROR in THijing::GetPJTZ(key1):\n");
      printf("      key1=%i is out of range [1..300]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>500 ) {
      printf("ERROR in THijing::GetPJTZ(key1,key2):\n");
      printf("      key2=%i is out of range [1..500]\n",key2);
      return 0;
   }
   
   return HIJJET1.pjtz[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPJTE(Int_t key1, Int_t key2) const
{
// Get values of array PJTE in common HIJJET1
   if ( key1<1 || key1>300 ) {
      printf("ERROR in THijing::GetPJTE(key1):\n");
      printf("      key1=%i is out of range [1..300]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>500 ) {
      printf("ERROR in THijing::GetPJTE(key1,key2):\n");
      printf("      key2=%i is out of range [1..500]\n",key2);
      return 0;
   }
   
   return HIJJET1.pjte[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPJTM(Int_t key1, Int_t key2) const
{
// Get values of array PJTM in common HIJJET1
   if ( key1<1 || key1>300 ) {
      printf("ERROR in THijing::GetPJTM(key1):\n");
      printf("      key1=%i is out of range [1..300]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>500 ) {
      printf("ERROR in THijing::GetPJTM(key1,key2):\n");
      printf("      key2=%i is out of range [1..500]\n",key2);
      return 0;
   }
   
   return HIJJET1.pjtm[key2-1][key1-1];
}

//====================== access to common HIJJET1 ===============================

//______________________________________________________________________________
Int_t THijing::GetNSG() const
{
// Get value of NSG in common HIJJET2
   return HIJJET2.nsg;
}

//______________________________________________________________________________
Int_t THijing::GetNJSG(Int_t key) const
{
// Get values of array NJSG in common HIJJET2
   if ( key<1 || key>900 ) {
      printf ("ERROR in THijing:GetNJSG(key): \n ");
      printf ("      key=%i is out of range [1..900]!\n",key);
      return 0;
   }

   return HIJJET2.njsg[key-1];

}

//______________________________________________________________________________
Int_t THijing::GetIASG(Int_t key1, Int_t key2) const
{
// Get values of IASG in common HIJJET2
   if ( key1<1 || key1>900 ) {
      printf("ERROR in THijing::GetIASG(key1):\n");
      printf("      key1=%i is out of range [1..900]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>3 ) {
      printf("ERROR in THijing::GetIASG(key1,key2):\n");
      printf("      key2=%i is out of range [1..3]\n",key2);
      return 0;
   }
   
   return HIJJET2.iasg[key2-1][key1-1];
}

//______________________________________________________________________________
Int_t THijing::GetK1SG(Int_t key1, Int_t key2) const
{
// Get values of K1SG in common HIJJET2
   if ( key1<1 || key1>900 ) {
      printf("ERROR in THijing::GetK1SG(key1):\n");
      printf("      key1=%i is out of range [1..900]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>100 ) {
      printf("ERROR in THijing::GetK1SG(key1,key2):\n");
      printf("      key2=%i is out of range [1..100]\n",key2);
      return 0;
   }
   
   return HIJJET2.k1sg[key2-1][key1-1];
}

//______________________________________________________________________________
Int_t THijing::GetK2SG(Int_t key1, Int_t key2) const
{
// Get values of K2SG in common HIJJET2
   if ( key1<1 || key1>900 ) {
      printf("ERROR in THijing::GetK2SG(key1):\n");
      printf("      key1=%i is out of range [1..900]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>100 ) {
      printf("ERROR in THijing::GetK2SG(key1,key2):\n");
      printf("      key2=%i is out of range [1..100]\n",key2);
      return 0;
   }
   
   return HIJJET2.k2sg[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPXSG(Int_t key1, Int_t key2) const
{
// Get values of PXSG in common HIJJET2
   if ( key1<1 || key1>900 ) {
      printf("ERROR in THijing::GetPXSG(key1):\n");
      printf("      key1=%i is out of range [1..900]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>100 ) {
      printf("ERROR in THijing::GetPXSG(key1,key2):\n");
      printf("      key2=%i is out of range [1..100]\n",key2);
      return 0;
   }
   
   return HIJJET2.pxsg[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPYSG(Int_t key1, Int_t key2) const
{
// Get values of PYSG in common HIJJET2
   if ( key1<1 || key1>900 ) {
      printf("ERROR in THijing::GetPYSG(key1):\n");
      printf("      key1=%i is out of range [1..900]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>100 ) {
      printf("ERROR in THijing::GetPYSG(key1,key2):\n");
      printf("      key2=%i is out of range [1..100]\n",key2);
      return 0;
   }
   
   return HIJJET2.pysg[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPZSG(Int_t key1, Int_t key2) const
{
// Get values of PZSG in common HIJJET2
   if ( key1<1 || key1>900 ) {
      printf("ERROR in THijing::GetPZSG(key1):\n");
      printf("      key1=%i is out of range [1..900]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>100 ) {
      printf("ERROR in THijing::GetPZSG(key1,key2):\n");
      printf("      key2=%i is out of range [1..100]\n",key2);
      return 0;
   }
   
   return HIJJET2.pzsg[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPESG(Int_t key1, Int_t key2) const
{
// Get values of PESG in common HIJJET2
   if ( key1<1 || key1>900 ) {
      printf("ERROR in THijing::GetPESG(key1):\n");
      printf("      key1=%i is out of range [1..900]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>100 ) {
      printf("ERROR in THijing::GetPESG(key1,key2):\n");
      printf("      key2=%i is out of range [1..100]\n",key2);
      return 0;
   }
   
   return HIJJET2.pesg[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPMSG(Int_t key1, Int_t key2) const
{
// Get values of PMSG in common HIJJET2
   if ( key1<1 || key1>900 ) {
      printf("ERROR in THijing::GetPMSG(key1):\n");
      printf("      key1=%i is out of range [1..900]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>100 ) {
      printf("ERROR in THijing::GetPMSG(key1,key2):\n");
      printf("      key2=%i is out of range [1..100]\n",key2);
      return 0;
   }
   
   return HIJJET2.pmsg[key2-1][key1-1];
}

//====================== access to common HISTRNG ===============================

//______________________________________________________________________________
Int_t THijing::GetNFP(Int_t key1, Int_t key2) const
{
// Get values of array NFP in common HISTRNG
   if ( key1<1 || key1>300 ) {
      printf("ERROR in THijing::GetNFP(key1):\n");
      printf("      key1=%i is out of range [1..300]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>15 ) {
      printf("ERROR in THijing::GetNFP(key1,key2):\n");
      printf("      key2=%i is out of range [1..15]\n",key2);
      return 0;
   }
   
   return HISTRNG.nfp[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPP(Int_t key1, Int_t key2) const
{
// Get values of array PP in common HISTRNG
   if ( key1<1 || key1>300 ) {
      printf("ERROR in THijing::GetPP(key1):\n");
      printf("      key1=%i is out of range [1..300]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>15 ) {
      printf("ERROR in THijing::GetPP(key1,key2):\n");
      printf("      key2=%i is out of range [1..15]\n",key2);
      return 0;
   }
   
   return HISTRNG.pp[key2-1][key1-1];
}

//______________________________________________________________________________
Int_t THijing::GetNFT(Int_t key1, Int_t key2) const
{
// Get values of array NFT in common HISTRNG
   if ( key1<1 || key1>300 ) {
      printf("ERROR in THijing::GetNFT(key1):\n");
      printf("      key1=%i is out of range [1..300]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>15 ) {
      printf("ERROR in THijing::GetNFT(key1,key2):\n");
      printf("      key2=%i is out of range [1..15]\n",key2);
      return 0;
   }
   
   return HISTRNG.nft[key2-1][key1-1];
}

//______________________________________________________________________________
Float_t THijing::GetPT(Int_t key1, Int_t key2) const
{
// Get values of array PT in common HISTRNG
   if ( key1<1 || key1>300 ) {
      printf("ERROR in THijing::GetPT(key1):\n");
      printf("      key1=%i is out of range [1..300]\n",key1);
      return 0;
   }
   if ( key2<1 || key2>15 ) {
      printf("ERROR in THijing::GetPT(key1,key2):\n");
      printf("      key2=%i is out of range [1..15]\n",key2);
      return 0;
   }
   
   return HISTRNG.pt[key2-1][key1-1];
}

void THijing::SetPARJ(Int_t key, Float_t parm) 
{
// Set values of array PARJ in common HISTRNG
    if ( key < 1 || key > 200) {
	printf("ERROR in THijing::SetPARJ(key,parm):\n");
	printf("      key=%i is out of range [1..200]\n",key);
    }
    
    LUDAT1_HIJING.parj[key-1] = parm;
}


void THijing::SetMSTJ(Int_t key, Int_t parm) 
{
// Set values of array MSTJ in common HISTRNG
    if ( key < 1 || key > 200) {
	printf("ERROR in THijing::SetMSTJ(key,parm):\n");
	printf("      key=%i is out of range [1..200]\n",key);
    }
    
    LUDAT1_HIJING.mstj[key-1] = parm;
}


//====================== access to Hijing subroutines =========================


//______________________________________________________________________________
void THijing::Initialize()
{
//////////////////////////////////////////////////////////////////////////////////
// Calls Hijset with the either default parameters or the ones set by the user  //
// via SetEFRM, SetFRAME, SetPROJ, SetTARG, SetIAP, SetIZP, SetIAT, SetIZT      //
//////////////////////////////////////////////////////////////////////////////////

   if ( (!strcmp(fFrame.Data(), "CMS     "  )) &&
        (!strcmp(fFrame.Data(), "LAB     "  ))){
      printf("WARNING! In THijing:Initialize():\n");
      printf(" specified frame=%s is neither CMS or LAB\n",fFrame.Data());
      printf(" resetting to default \"CMS\" .");
      fFrame="CMS";
   }

   if ( (!strcmp(fProj.Data(), "A       "     )) &&
        (!strcmp(fProj.Data(), "P       "     )) &&
        (!strcmp(fProj.Data(), "PBAR    "  ))){
      printf("WARNING! In THijing:Initialize():\n");
      printf(" specified projectile=%s is neither A, P or PBAR\n",fProj.Data());
      printf(" resetting to default \"A\" .");
      fProj="A";
   }

   if ( (!strcmp(fTarg.Data(), "A       "     )) &&
        (!strcmp(fTarg.Data(), "P       "     )) &&
        (!strcmp(fTarg.Data(), "PBAR    "  ))){
      printf("WARNING! In THijing:Initialize():\n");
      printf(" specified target=%s is neither A, P or PBAR\n",fTarg.Data());
      printf(" resetting to default \"A\" .");
      fTarg="A";
   }

   printf(" %s-%s at %f GeV \n",fProj.Data(),fTarg.Data(),fEfrm);

   Hijset(fEfrm,fFrame.Data(),fProj.Data(),fTarg.Data(),fIap,fIzp,fIat,fIzt);

   printf(" %s-%s at %f GeV \n",fProj.Data(),fTarg.Data(),fEfrm);
}


//______________________________________________________________________________
void THijing::GenerateEvent()
{
// Generates one event;

   Hijing(fFrame.Data(),fBmin,fBmax);

}
//______________________________________________________________________________
void THijing::Hijset(float efrm, const char *frame, const char *proj,  
		     const char *targ, int iap, int izp, int iat, int izt)
{
// Call HIJING routine HIJSET passing the parameters in a way accepted by
// Fortran routines				   

  int s1 = strlen(frame);
  int s2 = strlen(proj);
  int s3 = strlen(targ);
  printf("s1 = %d s2 = %d s3 = %d\n",s1,s2,s3);
#ifndef WIN32 
  hijset(efrm, frame, proj, targ, iap, izp, iat, izt, s1, s2, s3);
#else
  hijset(efrm, frame, s1, proj, s2, targ, s3, iap, izp, iat, izt);
#endif
}
//______________________________________________________________________________
void THijing::Hijing(const char *frame, float bmin, float bmax)
{
// Call HIJING routine HIJSET passing the parameters in a way accepted by
// Fortran routines				   

  int s1 = strlen(frame);
  
#ifndef WIN32 
  hijing(frame, bmin, bmax, s1);
#else
  hijing(frame, s1, bmin, bmax);
#endif
}


Float_t  THijing::Profile(float b)
{
// Call HIJING routine PROFILE 
  return profile(b);
}


void  THijing::Rluget(Int_t lfn, Int_t move)
{
// write seed to file
  rluget_hijing(lfn, move);
}


void  THijing::Rluset(Int_t lfn, Int_t move)
{
// read seed from file 
  rluset_hijing(lfn, move);
}

