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

//*KEEP,THIJING.
#include "THijing.h"
//*KEEP,HCOMMON.
#include "Hcommon.h"
//*KEEP,TMCParticle.
//#include "TMCParticle.h"
//*KEEP,TParticle,T=C++.
#include "TParticle.h"
//*KEND.

//*KEEP,TCanvas.
//#include "TCanvas.h"
//*KEEP,TView.
//#include "TView.h"
//*KEEP,TROOT.
#include "TROOT.h"
//*KEEP,TPaveText.
//#include "TPaveText.h"
//*KEND.

#ifndef WIN32
# define hijset hijset_
# define hijing hijing_
# define type_of_call
#else
# define hijset HIJSET
# define hijing HIJING
# define type_of_call _stdcall
#endif

#ifndef WIN32
//extern "C" void type_of_call hijset(float &efrm, const char *frame, 
//				   const char *proj, const char *targ,
//				   int &iap, int &izp, int &iat,
//				   int &izt, Long_t l_frame, 
//				   Long_t l_proj, Long_t l_targ);
extern "C" void type_of_call hijset(Float_t & , const char *, 
				   const char *, const char *,
				   Int_t & , Int_t &, Int_t &,
				   Int_t &, const int, 
				   const int, const int);
//extern "C" void type_of_call hijing(const char *frame, float &bmin,
//                                   float &bmax, Long_t l_frame);
extern "C" void type_of_call hijing(const char *, Float_t &,
                                   Float_t &, const int);
#else
//extern "C" void type_of_call hijset(float &efrm, const char *frame, 
//				   Long_t l_frame, const char *proj, 
//				   Long_t l_proj,  const char *targ, 
//				   Long_t l_targ,
//				   int &iap, int &izp, int &iat,
//				   int &izt);
//extern "C" void type_of_call hijing(const char *frame, Long_t l_frame, 
//				   float &bmin, float &bmax);
#endif

ClassImp(THijing)

// void THijing::Streamer(TBuffer &R__b){}
//______________________________________________________________________________
THijing::THijing() : TGenerator("Hijing","Hijing")
{
// THijing constructor: creates a TClonesArray in which it will store all
// particles. Note that there may be only one functional THijing object
// at a time, so it's not use to create more than one instance of it.

//    delete fParticles; // was allocated as TObjArray in TGenerator

//    fParticles = new TClonesArray("TMCParticle",50);

}

//______________________________________________________________________________
THijing::THijing(Float_t efrm, const char *frame="CMS", 
	const char *proj="A", const char *targ="A", Int_t iap=207, 
	Int_t izp=82, Int_t iat=207, Int_t izt=82, Float_t bmin=0, 
	Float_t bmax=20) : TGenerator("Hijing","Hijing")
{
// THijing constructor: creates a TClonesArray in which it will store all
// particles. Note that there may be only one functional THijing object
// at a time, so it's not use to create more than one instance of it.

//    delete fParticles; // was allocated as TObjArray in TGenerator

//    fParticles = new TClonesArray("TMCParticle",50);

      fEfrm=efrm;
      fFrame=frame;
      fProj=proj;
      fTarg=targ;
      fIap=iap;
      fIzp=izp;
      fIat=iat;
      fIzt=izt;
      fBmin=bmin;
      fBmax=bmax;
}

//______________________________________________________________________________
THijing::~THijing()
{
// Destroys the object, deletes and disposes all TMCParticles currently on list.

//    if (fParticles) {
//     fParticles->Delete();
//      delete fParticles;
//      fParticles = 0;
//   }
}

//______________________________________________________________________________
//void THijing::Draw(Option_t *option)
//{
// Event display - not supported for THijing yet.

//   if (!gPad) {
//      if (!gROOT->GetMakeDefCanvas()) return;
//      (gROOT->GetMakeDefCanvas())();
//      gPad->GetCanvas()->SetFillColor(13);
//   }

//   static Float_t rbox = 1000;
//   Float_t rmin[3],rmax[3];
//   TView *view = gPad->GetView();
//   if (!strstr(option,"same")) {
//      if (view) { view->GetRange(rmin,rmax); rbox = rmax[2];}
//      gPad->Clear();
//   }

//   AppendPad(option);

//   view = gPad->GetView();
//   //    compute 3D view
//   if (view) {
//      view->GetRange(rmin,rmax);
//      rbox = rmax[2];
//   } else {
//      view = new TView(1);
//      view->SetRange(-rbox,-rbox,-rbox, rbox,rbox,rbox );
//   }

//   TPaveText *pt = new TPaveText(-0.94,0.85,-0.25,0.98,"br");
//   pt->AddText((char*)GetName());
//   pt->AddText((char*)GetTitle());
//   pt->SetFillColor(42);
//   pt->Draw();
//}

//______________________________________________________________________________
//TObjArray *THijing::ImportParticles(Option_t *)
//{
// Fills TClonesArray fParticles list with particles from common LUJETS.
// Old contents of a list are cleared. This function should be called after
// any change in common LUJETS, however GetParticles() method  calls it
// automatically - user don't need to care about it. In case you make a call
// to LuExec() you must call this method yourself to transfer new data from
// common LUJETS to the fParticles list.
//
//   fParticles->Clear();
//
//   Int_t numpart   = LUJETS.n;
//   TClonesArray &a = *((TClonesArray*)fParticles);
//
//   for (Int_t i = 0; i < numpart; i++) {
//	new(a[i]) TMCParticle(LUJETS.k[0][i] ,
//			      LUJETS.k[1][i] ,
//			      LUJETS.k[2][i] ,
//			      LUJETS.k[3][i] ,
//			      LUJETS.k[4][i] ,
//
//			      LUJETS.p[0][i] ,
//			      LUJETS.p[1][i] ,
//			      LUJETS.p[2][i] ,
//			      LUJETS.p[3][i] ,
//			      LUJETS.p[4][i] ,
//
//			      LUJETS.v[0][i] ,
//			      LUJETS.v[1][i] ,
//			      LUJETS.v[2][i] ,
//			      LUJETS.v[3][i] ,
//			      LUJETS.v[4][i]);

//   }
//   return fParticles;
//}

//______________________________________________________________________________
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
  TClonesArray &Particles = *particles;
  Particles.Clear();
  Int_t numpart = HIMAIN1.natt;
  if (!strcmp(option,"") || !strcmp(option,"Final")) {
    for (Int_t i = 0; i<=numpart; i++) {
	if (HIMAIN2.katt[3][i] == 1) {
//
//  Use the common block values for the TParticle constructor
//
	  new(Particles[i]) TParticle(
			      HIMAIN2.katt[0][i] ,
			      HIMAIN2.katt[1][i] ,
			      HIMAIN2.katt[2][i] ,
			      -1,
			      -1,
			      -1,

			      HIMAIN2.patt[0][i] ,
			      HIMAIN2.patt[1][i] ,
			      HIMAIN2.patt[2][i] ,
			      HIMAIN2.patt[3][i] ,

			      0,
			      0,
			      0,
			      0);
	}
    }
  }
  else if (!strcmp(option,"All")) {
    for (Int_t i = 0; i<=numpart; i++) {
	  new(Particles[i]) TParticle(
			      HIMAIN2.katt[0][i] ,
			      HIMAIN2.katt[1][i] ,
			      HIMAIN2.katt[2][i] ,
			      -1,
			      -1,
			      -1,

			      HIMAIN2.patt[0][i] ,
			      HIMAIN2.patt[1][i] ,
			      HIMAIN2.patt[2][i] ,
			      HIMAIN2.patt[3][i] ,

			      0,
			      0,
			      0,
			      0);
    }
  }
  return numpart;
}

//______________________________________________________________________________
void THijing::SetEFRM(Float_t efrm)
{
   fEfrm=efrm;
} 
//______________________________________________________________________________
void THijing::SetFRAME(const char* frame)
{
   fFrame=frame;
} 
//______________________________________________________________________________
void THijing::SetPROJ(const char* proj)
{
   fProj=proj;
} 
//______________________________________________________________________________
void THijing::SetTARG(const char* targ)
{
   fTarg=targ;
} 
//______________________________________________________________________________
void THijing::SetIAP(Int_t iap)
{
   fIap=iap;
} 
//______________________________________________________________________________
void THijing::SetIZP(Int_t izp)
{
   fIzp=izp;
} 
//______________________________________________________________________________
void THijing::SetIAT(Int_t iat)
{
   fIat=iat;
} 
//______________________________________________________________________________
void THijing::SetIZT(Int_t izt)
{
   fIzt=izt;
} 
//______________________________________________________________________________
void THijing::SetBMIN(Float_t bmin)
{
   fBmin=bmin;
} 
//______________________________________________________________________________
void THijing::SetBMAX(Float_t bmax)
{
   fBmax=bmax;
} 
//______________________________________________________________________________
Float_t THijing::GetEFRM() const
{
   return fEfrm;
} 
//______________________________________________________________________________
const char* THijing::GetFRAME() const
{
   return fFrame.Data();
} 
//______________________________________________________________________________
const char* THijing::GetPROJ() const
{
   return fProj;
} 
//______________________________________________________________________________
const char* THijing::GetTARG() const
{
   return fTarg;
} 
//______________________________________________________________________________
Int_t THijing::GetIAP() const
{
   return fIap;
} 
//______________________________________________________________________________
Int_t THijing::GetIZP() const
{
   return fIzp;
} 
//______________________________________________________________________________
Int_t THijing::GetIAT() const
{
   return fIat;
} 
//______________________________________________________________________________
Int_t THijing::GetIZT() const
{
   return fIzt;
} 
//______________________________________________________________________________
Float_t THijing::GetBMIN() const
{
   return fBmin;
} 
//______________________________________________________________________________
Float_t THijing::GetBMAX() const
{
   return fBmax;
} 

//====================== access to common HIPARNT ===============================

//______________________________________________________________________________
void THijing::SetHIPR1(Int_t key,Float_t value)
{
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

   return HIMAIN1.natt;

}

//______________________________________________________________________________
Float_t  THijing::GetEATT() const
{

   return HIMAIN1.eatt;

}

//______________________________________________________________________________
Int_t  THijing::GetJATT() const
{

   return HIMAIN1.jatt;

}

//______________________________________________________________________________
Int_t  THijing::GetNT() const
{

   return HIMAIN1.nt;

}

//______________________________________________________________________________
Int_t  THijing::GetNP() const
{

   return HIMAIN1.np;

}


//______________________________________________________________________________
Int_t  THijing::GetN0() const
{

   return HIMAIN1.n0;

}
//______________________________________________________________________________
Int_t  THijing::GetN01() const
{

   return HIMAIN1.n01;

}

//______________________________________________________________________________
Int_t  THijing::GetN10() const
{

   return HIMAIN1.n10;

}

//______________________________________________________________________________
Int_t  THijing::GetN11() const
{

   return HIMAIN1.n11;

}

//______________________________________________________________________________
Float_t  THijing::GetBB() const
{

   return HIMAIN1.bb;

}

//====================== access to common HIMAIN2 ===============================

//______________________________________________________________________________
Int_t THijing::GetKATT(Int_t key1, Int_t key2) const
{
   if ( key1<1 || key1>130000 ) {
      printf("ERROR in THijing::GetKATT(key1,key2):\n");
      printf("      key1=%i is out of range [1..130000]\n",key1);
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
   if ( key1<1 || key1>130000 ) {
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

//====================== access to common HIJJET1 ===============================

//______________________________________________________________________________
Int_t THijing::GetNPJ(Int_t key) const
{
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
   return HIJJET2.nsg;
}

//______________________________________________________________________________
Int_t THijing::GetNJSG(Int_t key) const
{
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



//====================== access to Hijing subroutines =========================


//______________________________________________________________________________
void THijing::Initialize()
{
//////////////////////////////////////////////////////////////////////////////////
// Calls Hijset with the either default parameters or the ones set by the user  //
// via SetEFRM, SetFRAME, SetPROJ, SetTARG, SetIAP, SetIZP, SetIAT, SetIZT      //
//////////////////////////////////////////////////////////////////////////////////

   if ( (!strcmp(fFrame.Data(), "CMS"  )) &&
        (!strcmp(fFrame.Data(), "LAB"  ))){
      printf("WARNING! In THijing:Initialize():\n");
      printf(" specified frame=%s is neither CMS or LAB\n",fFrame.Data());
      printf(" resetting to default \"CMS\" .");
      fFrame="CMS";
   }

   if ( (!strcmp(fProj.Data(), "A"     )) &&
        (!strcmp(fProj.Data(), "P"     )) &&
        (!strcmp(fProj.Data(), "PBAR"  ))){
      printf("WARNING! In THijing:Initialize():\n");
      printf(" specified projectile=%s is neither A, P or PBAR\n",fProj.Data());
      printf(" resetting to default \"A\" .");
      fProj="A";
   }

   if ( (!strcmp(fTarg.Data(), "A"     )) &&
        (!strcmp(fTarg.Data(), "P"     )) &&
        (!strcmp(fTarg.Data(), "PBAR"  ))){
      printf("WARNING! In THijing:Initialize():\n");
      printf(" specified target=%s is neither A, P or PBAR\n",fTarg.Data());
      printf(" resetting to default \"A\" .");
      fTarg="A";
   }



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
