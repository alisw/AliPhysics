///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Interface to the HIGZ package for the GEANT drawing package              //
//                                                                           //
//Begin_Html
/*
<img src="gif/THIGZClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "THIGZ.h"
#include "TGeant3.h"
#include "TGraph.h"
#include "TColor.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TPolyMarker.h"
#include "TPaveLabel.h"
#include "TText.h"
#include "Rtypes.h" 
#include "AliCallf77.h" 
#include "TPaveTree.h" 
#include <string.h> 
#include <ctype.h> 
   
#ifndef WIN32
# define hplint hplint_
# define hplend hplend_
# define hplfra hplfra_
# define igmeta igmeta_
# define iselnt iselnt_
# define igiwty igiwty_
# define igqwk  igqwk_

# define iacwk  iacwk_
# define iclrwk iclrwk_
# define iclwk  iclwk_
# define idawk  idawk_
# define ifa    ifa_
# define igbox  igbox_
# define ightor ightor_
# define igpave igpave_
# define igpid  igpid_
# define igq    igq_
# define igrng  igrng_
# define igsa   igsa_
# define igset  igset_
# define igterm igterm_
# define iopwk  iopwk_
# define ipl    ipl_
# define ipm    ipm_
# define irqlc  irqlc_
# define iscr   iscr_
# define isfaci isfaci_
# define isfais isfais_
# define isln   isln_
# define ismk   ismk_
# define islwsc islwsc_
# define isplci isplci_
# define ispmci ispmci_
# define istxci istxci_
# define isvp   isvp_
# define iswn   iswn_
# define itx    itx_
 
#else
# define iacwk  IACWK
# define iclrwk ICLRWK
# define iclwk  ICLWK
# define idawk  IDAWK
# define ifa    IFA
# define igbox  IGBOX
# define ightor IGHTOR
# define igpave IGPAVE
# define igpid  IGPID
# define igq    IGQ
# define igrng  IGRNG
# define igsa   IGSA
# define igset  IGSET
# define igterm IGTERM
# define iopwk  IOPWK
# define ipl    IPL
# define ipm    IPM
# define irqlc  IRQLC
# define iscr   ISCR
# define isfaci ISFACI
# define isfais ISFAIS
# define isln   ISLN
# define ismk   ISMK
# define islwsc ISLWSC
# define isplci ISPLCI
# define ispmci ISPMCI
# define istxci ISTXCI
# define isvp   ISVP
# define iswn   ISWN
# define itx    ITX
#endif

static Int_t npid = 0;
static char cpar[1200];

THIGZ *higz = 0;

ClassImp(THIGZ)
   
//_______________________________________________________________
THIGZ::THIGZ()
      :TCanvas()
{
}
   
//_______________________________________________________________
THIGZ::THIGZ(Int_t size)
      :TCanvas("higz","higz",size,size)
{
   higz = this;
   Reset();
   SetFillColor(10);
}
//_______________________________________________________________
THIGZ::~THIGZ()
{
}
//_______________________________________________________________
Float_t THIGZ::Get(const char *name)
{
   if (!strcmp(name,"FAIS")) return fFAIS;
   if (!strcmp(name,"FASI")) return fFASI;
   if (!strcmp(name,"LTYP")) return fLTYP;
   if (!strcmp(name,"BASL")) return fBASL;
   if (!strcmp(name,"LWID")) return fLWID;
   if (!strcmp(name,"MTYP")) return fMTYP;
   if (!strcmp(name,"MSCF")) return fMSCF;
   if (!strcmp(name,"PLCI")) return fPLCI;
   if (!strcmp(name,"PMCI")) return fPMCI;
   if (!strcmp(name,"FACI")) return fFACI;
   if (!strcmp(name,"TXCI")) return fTXCI;
   if (!strcmp(name,"TXAL")) return fTXAL;
   if (!strcmp(name,"CHHE")) return fCHHE*fY2;
   if (!strcmp(name,"TANG")) return fTANG;
   if (!strcmp(name,"TXFP")) return fTXFP;
   if (!strcmp(name,"BORD")) return fBORD;
   if (!strcmp(name,"NCOL")) return fNCOL;
   if (!strcmp(name,"DRMD")) return fDRMD;
   if (!strcmp(name,"SYNC")) return fSYNC;
   if (!strcmp(name,"CLIP")) return fCLIP;
   if (!strcmp(name,"2BUF")) return f2BUF;
   return -1;
}
//_______________________________________________________________
void THIGZ::Set(const char *name, Float_t val)
{
   if (!strcmp(name,"FAIS")) {fFAIS = Int_t(val); return;}
   if (!strcmp(name,"FASI")) {fFASI = Int_t(val); return;}
   if (!strcmp(name,"LTYP")) {fLTYP = Int_t(val); return;}
   if (!strcmp(name,"BASL")) {fBASL = val; return;}
   if (!strcmp(name,"LWID")) {fLWID = val; if (val <= 0) fLWID=1; return;}
   if (!strcmp(name,"MTYP")) {fMTYP = Int_t(val); return;}
   if (!strcmp(name,"MSCF")) {fMSCF = val; return;}
   if (!strcmp(name,"PLCI")) {fPLCI = Int_t(val); return;}
   if (!strcmp(name,"PMCI")) {fPMCI = Int_t(val); return;}
   if (!strcmp(name,"FACI")) {fFACI = Int_t(val)%1000; return;}
   if (!strcmp(name,"TXCI")) {fTXCI = Int_t(val); return;}
   if (!strcmp(name,"TXAL")) {
      Int_t align = Int_t(val); 
      Int_t aligh = align/10;
      Int_t aligv = align%10;
      Int_t ah = aligh;
      Int_t av = aligv;
      if (aligh == 0) ah = 1;
      if (aligv == 1 || aligh == 2) av = 3;
      if (aligv == 3)               av = 2;
      if (aligv == 0)               av = 1;
      fTXAL = 10*ah + av;
      return;
   }
   if (!strcmp(name,"CHHE")) {fCHHE = val/fY2; if (fCHHE <= 0) fCHHE = 0.03; return;}
   if (!strcmp(name,"TANG")) {fTANG = val; return;}
   if (!strcmp(name,"TXFP")) {fTXFP = Int_t(val); return;}
   if (!strcmp(name,"BORD")) {fBORD = Int_t(val); return;}
   if (!strcmp(name,"NCOL")) {fNCOL = Int_t(val); return;}
   if (!strcmp(name,"DRMD")) {fDRMD = Int_t(val); return;}
   if (!strcmp(name,"SYNC")) {fSYNC = Int_t(val); return;}
   if (!strcmp(name,"CLIP")) {fCLIP = Int_t(val); return;}
   if (!strcmp(name,"2BUF")) {f2BUF = Int_t(val); return;}
}
//_______________________________________________________________
void THIGZ::Reset(Option_t *)
{
   fFAIS = 0;
   fFASI = 1;
   fLTYP = 1;
   fBASL = 1;
   fLWID = 1;
   fMTYP = 1;
   fMSCF = 1;
   fPLCI = 1;
   fPMCI = 1;
   fFACI = 1;
   fTXCI = 1;
   fTXAL = 12;
   fCHHE = 0.03;
   fTANG = 0;
   fTXFP = 62;
   fBORD = 0;
   fNCOL = 8;
   fDRMD = 1;
   fSYNC = 1;
   fCLIP = 0;
   f2BUF = 0;
   fPID  = 0;
   fPname= "";
}
//___________________________________________ 
void THIGZ::Gsatt(const char *name, const char *att, Int_t val)
{ 
//  NAME   Volume name
//  IOPT   Name of the attribute to be set
//  IVAL   Value to which the attribute is to be set
//
//  name= "*" stands for all the volumes.
//  iopt can be chosen among the following :
//     see TGeant3::Gsatt for more details
//

   AliMC* pMC = AliMC::GetMC(); 
   pMC->Gsatt(name,att,val);
}

//___________________________________________ 
void THIGZ::Gdopt(const char *name,const char *value)
{ 
//  NAME   Option name
//  VALUE  Option value
//
//  To set/modify the drawing options.
//     see TGeant3::Gdopt for more details
//
   AliMC* pMC = AliMC::GetMC(); 
   pMC->Gdopt(name,value);
}

 
//___________________________________________ 
void THIGZ::Gdraw(const char *name,Float_t theta, Float_t phi, Float_t psi,Float_t u0,Float_t v0,Float_t ul,Float_t vl)
{ 
//  NAME   Volume name
//  +
//  THETA  Viewing angle theta (for 3D projection)
//  PHI    Viewing angle phi (for 3D projection)
//  PSI    Viewing angle psi (for 2D rotation)
//  U0     U-coord. (horizontal) of volume origin
//  V0     V-coord. (vertical) of volume origin
//  SU     Scale factor for U-coord.
//  SV     Scale factor for V-coord.
//     see TGeant3::Gdraw for more details
//

   AliMC* pMC = AliMC::GetMC(); 
   pMC->Gdraw(name,theta,phi,psi,u0,v0,ul,vl);
}

 
//___________________________________________ 
void THIGZ::Gdrawc(const char *name,Int_t axis, Float_t cut,Float_t u0,Float_t v0,Float_t ul,Float_t vl)
{ 
//  NAME   Volume name
//  CAXIS  Axis value
//  CUTVAL Cut plane distance from the origin along the axis
//  +
//  U0     U-coord. (horizontal) of volume origin
//  V0     V-coord. (vertical) of volume origin
//  SU     Scale factor for U-coord.
//  SV     Scale factor for V-coord.
//     see TGeant3::Gdrawc for more details
//
 
   TGeant3 *geant3=(TGeant3*)AliMC::GetMC();
   geant3->Gdrawc(name,axis,cut,u0,v0,ul,vl);
}

 
//___________________________________________ 
void THIGZ::Gdspec(const char *name)
{ 
//  NAME   Volume name
//
//  Shows 3 views of the volume (two cut-views and a 3D view), together with
//  its geometrical specifications. The 3D drawing will
//  be performed according the current values of the options HIDE and
//  SHAD and according the current CVOL clipping parameters for that
//  volume.
//     see TGeant3::Gdspec for more details
//

   TGeant3 *geant3=(TGeant3*)AliMC::GetMC();
   geant3->Gdspec(name);
}

 
//___________________________________________ 
void THIGZ::Gdtree(const char *name,Int_t levmax, Int_t isel)
{ 
//  NAME   Volume name
//  LEVMAX Depth level
//  ISELT  Options
//
//  This function draws the logical tree,
//  Each volume in the tree is represented by a TPaveTree object.
//  Double-clicking on a TPaveTree draws the specs of the corresponding volume.
//  Use TPaveTree pop-up menu to select:
//    - drawing specs
//    - drawing tree
//    - drawing tree of parent
//

   TGeant3 *geant3=(TGeant3*)AliMC::GetMC();
   geant3->Gdtree(name,levmax,isel);
}

 
//______________________________________________________________________________
void THIGZ::SetBOMB(Float_t boom)
{
//  BOOM  : Exploding factor for volumes position 
// 
//  To 'explode' the detector. If BOOM is positive (values smaller
//  than 1. are suggested, but any value is possible)
//  all the volumes are shifted by a distance
//  proportional to BOOM along the direction between their centre
//  and the origin of the MARS; the volumes which are symmetric
//  with respect to this origin are simply not shown.
//  BOOM equal to 0 resets the normal mode.
//  A negative (greater than -1.) value of
//  BOOM will cause an 'implosion'; for even lower values of BOOM
//  the volumes' positions will be reflected respect to the origin.
//  This command can be useful to improve the 3D effect for very
//  complex detectors. The following commands will make explode the
//  detector:

   TGeant3 *geant3=(TGeant3*)AliMC::GetMC();
   geant3->SetBOMB(boom);
}

//_______________________________________________________________
extern "C" void type_of_call iacwk(Int_t &)
{
//   printf("iacwk called,wid=%d\n",wid);
}

//_______________________________________________________________
extern "C" void type_of_call iclrwk(Int_t &,Int_t &)
{
//   printf("iclrwk called\n");
   higz->Clear();
}

//_______________________________________________________________
extern "C" void type_of_call iclwk(Int_t &,Int_t &)
{
//   printf("iclwk called\n");

}

//_______________________________________________________________
extern "C" void type_of_call idawk(Int_t &)
{
//   printf("idawk called\n");
}

//_______________________________________________________________
extern "C" void type_of_call ifa(Int_t &n,Float_t *x, Float_t *y)
{
//   printf("ifa called, n=%d, pname=%s, pid=%d, x1=%f, y1=%f,x2=%f, y2=%f, x3=%f, y3=%f\n",n,higz->fPname.Data(),higz->fPID,x[0],y[0],x[1],y[1],x[2],y[2]);
   TGraph *gr = new TGraph(n,x,y);
   gr->SetFillColor(higz->fFACI);
   gr->Draw("f");
}

//_______________________________________________________________
extern "C" void type_of_call igbox(Float_t &x1,Float_t &x2,Float_t &y1,Float_t &y2)
{
   printf("igbox called, x1=%f, y1=%f, x2=%f, y2=%f\n",x1,y1,x2,y2);
   TBox *box = new TBox(x1,y1,x2,y2);
   box->SetLineColor(higz->fPLCI);
   box->SetFillColor(higz->fFACI);
   box->SetFillStyle(higz->fFASI);
   box->Draw();
}

//_______________________________________________________________
extern "C" void type_of_call ightor(Float_t &h,Float_t &l,Float_t &s,Float_t &r,Float_t &g,Float_t &b)
{
//   printf("ightor called\n");
   TColor *col1 = gROOT->GetColor(1);
   if (col1) col1->HLStoRGB(h,l,s,r,g,b);
}

//_______________________________________________________________
extern "C" void type_of_call igpave(Float_t &x1,Float_t &x2,Float_t &yy1,Float_t &yy2,Float_t &,Int_t &isbox,Int_t &isfram,const char *, const Int_t)
{
   char text[5];
   strncpy(text,(char*)&higz->fPID,4);
   text[4] = 0;
   Float_t y1 = yy1;
   Float_t y2 = yy2;
   if (y1 > y2) { y1 = yy2; y2 = yy1;}
   Float_t y = 0.5*(y1+y2);
   Float_t dymax = 0.06*(higz->GetY2()-higz->GetY1());
   Float_t dymin = -higz->PixeltoY(12);
   if (y2-y1 > dymax) {
      y1 = y - 0.5*dymax;
      y2 = y + 0.5*dymax;
   }
   if (y2-y1 < dymin) {
      y1 = y - 0.5*dymin;
      y2 = y + 0.5*dymin;
   }
   TPaveTree *pt = new TPaveTree(x1,y1,x2,y2,text);
   pt->SetTextColor(5);
   pt->SetFillColor(isbox%1000);
   pt->SetLineColor(isfram%1000);
   pt->Draw();
//   printf("igpave called, text=%s, Pname=%s, x1=%f, y1=%f, x2=%f, y2=%f, isbox=%d, isfram=%d\n",text,higz->fPname.Data(),x1,y1,x2,y2,isbox,isfram);
}

//_______________________________________________________________
#ifndef WIN32
extern "C" void type_of_call igpid(Int_t &,const char *name,Int_t &pid,const char *, const Int_t l1, const Int_t)
#else
extern "C" void type_of_call igpid(Int_t &,const char *name,const Int_t l1, Int_t &pid,const char *, const Int_t)
#endif
{
   npid++;
//   if(npid&100 == 0) printf("igpid called, npid=%d\n",npid);
   strncpy(cpar,name,l1); cpar[l1] = 0;
   higz->fPname = cpar;
   higz->fPID = pid;
      
//   char text[5];
//   strncpy(text,(char*)&higz->fPID,4);
//   text[4] = 0;
//   printf("igpid called, level=%d, name=%s, pid=%d, cpid=%s\n",level,cpar,pid,text);
}

//_______________________________________________________________
#ifndef WIN32
extern "C" void type_of_call igq(const char *name,Float_t &rval, const Int_t l1)
#else
extern "C" void type_of_call igq(const char *name,const Int_t l1, Float_t &rval)
#endif
{
   strncpy(cpar,name,l1); cpar[l1] = 0;
//   printf("igq called, name=%s\n",cpar);
   rval = higz->Get(cpar);
}

//_______________________________________________________________
extern "C" void type_of_call igrng(Float_t &xsize,Float_t &ysize)
{
//   printf("igrng called, xsize=%f, ysize=%f\n",xsize,ysize);
   higz->Range(0,0,xsize,ysize);
}

//_______________________________________________________________
extern "C" void type_of_call igsa(Int_t &)
{
//   printf("igsa called, iwk=%d\n",iwk);
}

//_______________________________________________________________
#ifndef WIN32
extern "C" void type_of_call igset(const char *name,Float_t &rval, const Int_t l1)
#else
extern "C" void type_of_call igset(const char *name, const Int_t l1,Float_t &rval)
#endif
{
   strncpy(cpar,name,l1); cpar[l1] = 0;
//   printf("igset called, name=%s, rval=%f\n",cpar,rval);
   higz->Set(cpar,rval);
}

//_______________________________________________________________
extern "C" void type_of_call igterm()
{
   printf("igterm called\n");
}

//_______________________________________________________________
extern "C" void type_of_call iopwk(Int_t &iwkid,Int_t &iconid,Int_t &iwtypi)
{
   printf("iopwk called, iwkid=%d, iconid=%d, iwtypi=%d\n",iwkid,iconid,iwtypi);
}

//_______________________________________________________________
extern "C" void type_of_call ipl(Int_t &n,Float_t *x,Float_t *y)
{
//   printf("ipl called, n=%d, x[0]=%f,y[0]=%f, x[1]=%f, y[1]=%f\n",n,x[0],y[0],x[1],y[1]);
   if (n <= 2) {
      TLine *l = new TLine(x[0],y[0],x[1],y[1]);
      l->SetLineColor(higz->fPLCI);
      l->SetLineStyle(higz->fLTYP);
      l->SetLineWidth(Short_t(higz->fLWID));
      l->Draw();
   } else {
      TPolyLine *pl = new TPolyLine(n,x,y);
      pl->SetLineColor(higz->fPLCI);
      pl->SetLineStyle(higz->fLTYP);
      pl->SetLineWidth(Short_t(higz->fLWID));
      pl->Draw();
   }
}

//_______________________________________________________________
extern "C" void type_of_call ipm(Int_t &n,Float_t *x,Float_t *y)
{
   printf("ipm called, n=%d\n",n);
   TPolyMarker *pm = new TPolyMarker(n,x,y);
   pm->SetMarkerColor(higz->fPMCI);
   pm->SetMarkerStyle(higz->fMTYP);
   pm->SetMarkerSize(higz->fMSCF);
   pm->Draw();
}

//_______________________________________________________________
extern "C" void type_of_call irqlc(Int_t&, Int_t&, Int_t&, Int_t&, Float_t&, Float_t&)

{
   printf("irqlc called\n");
}

//_______________________________________________________________
extern "C" void type_of_call iscr(Int_t &,Int_t &ici,Float_t &r,Float_t &g,Float_t &b)
{
//   printf("iscr called, ici=%d, r=%f, g=%f, b=%f\n",ici,r,g,b);
   Int_t color = ici;
//   if (color > 10) color += 35;
   TColor *col = gROOT->GetColor(color);
   if (!col) {
      col = new TColor(color,r,g,b);
      return;
   } else {
      col->SetRGB(r,g,b);
   }
}

//_______________________________________________________________
extern "C" void type_of_call isfaci(Int_t &col)
{
//   printf("isfaci called, col=%d\n",col);
   Int_t color = col%1000;
//   if (color > 10) color += 35;
   higz->fFACI = color;
}

//_______________________________________________________________
extern "C" void type_of_call isfais(Int_t &is)
{
//   printf("isfais called, is=%d\n",is);
   higz->fFAIS = is;
}

//_______________________________________________________________
extern "C" void type_of_call isln(Int_t &ln)
{
//   printf("isln called, ln=%d\n",ln);
   higz->fLTYP = ln;
}

//_______________________________________________________________
extern "C" void type_of_call ismk(Int_t &mk)
{
//   printf("ismk called, mk=%d\n",mk);
   higz->fMTYP = mk;
}

//_______________________________________________________________
extern "C" void type_of_call islwsc(Float_t &wl)
{
//   printf("islwsc called, wl=%f\n",wl);
   higz->fLWID = wl;
}

//_______________________________________________________________
extern "C" void type_of_call isplci(Int_t &col)
{
//   printf("isplci called, col=%d\n",col);
   Int_t color = col%1000;
//   if (color > 10) color += 35;
   higz->fPLCI = color;
}

//_______________________________________________________________
extern "C" void type_of_call ispmci(Int_t &col)
{
//   printf("ispmci called, col=%d\n",col);
   Int_t color = col%1000;
//   if (color > 10) color += 35;
   higz->fPMCI = color;
}

//_______________________________________________________________
extern "C" void type_of_call istxci(Int_t &col)
{
//   printf("istxci called, col=%d\n",col);
   Int_t color = col%1000;
//   if (color > 10) color += 35;
   higz->fTXCI = color;
}

//_______________________________________________________________
extern "C" void type_of_call isvp(Int_t &,Float_t &,Float_t &,Float_t &,Float_t &)
{
//   printf("isvp called, nt=%d, x1=%f, y1=%f, x2=%f, y2=%f\n",nt,x1,y1,x2,y2);
}

//_______________________________________________________________
extern "C" void type_of_call iswn(Int_t &,Float_t &x1,Float_t &x2,Float_t &y1,Float_t &y2)

{
//   printf("iswn called, nt=%d, x1=%f, y1=%f, x2=%f, y2=%f\n",nt,x1,y1,x2,y2);
   higz->Range(x1,y1,x2,y2);
}

//_______________________________________________________________
extern "C" void type_of_call itx(Float_t &x,Float_t &y,const char *ptext, const Int_t l1p)
{
   if (higz->fPname == "Tree") return;
   Int_t l1=l1p;
   strncpy(cpar,ptext,l1); cpar[l1] = 0;
//printf("itx called, x=%f, y=%f, text=%s, l1=%d\n",x,y,cpar,l1);
   while(l1) {
      if (cpar[l1-1] != ' ' && cpar[l1-1] != '$') break;
      l1--;
      cpar[l1] = 0;
   }
   char *small = strstr(cpar,"<");
   while(small && *small) {
      small[0] = tolower(small[1]);
      small++;
   }
   TText *text = new TText(x,y,cpar);
   text->SetTextColor(higz->fTXCI);
   text->SetTextSize(higz->fCHHE);
   text->SetTextFont(higz->fTXFP);
   text->SetTextAlign(higz->fTXAL);
   text->SetTextAngle(higz->fTANG);
   text->Draw();   
}

//_______________________________________________________________
extern "C" void type_of_call hplint(Int_t &)
{
   printf("hplint called\n");
   new THIGZ(600);
}

//_______________________________________________________________
extern "C" void type_of_call hplend()
{
   printf("hplend called\n");
}

//_______________________________________________________________
extern "C" void type_of_call hplfra(Float_t &x1,Float_t &x2,Float_t &y1, Float_t &y2,const char *, const Int_t)
{
//   printf("hplfra called, x1=%f, y1=%f, x2=%f, y2=%f\n",x1,y1,x2,y2);
   higz->Clear();
   higz->Range(x1,y1,x2,y2);
}

//_______________________________________________________________
extern "C" void type_of_call igmeta(Int_t &, Int_t &)
{
   printf("igmeta called\n");
}

//_______________________________________________________________
extern "C" void type_of_call iselnt(Int_t &)
{
//   printf("iselnt called, nt=%d\n",nt);
}

//_______________________________________________________________
extern "C" Int_t type_of_call igiwty(Int_t &)
{
//   printf("igiwty called, wid=%d\n",wid);
   return 1;
}

//_______________________________________________________________
extern "C" void type_of_call igqwk(Int_t &, const char *name, Float_t &rval, const Int_t l1)
{
   strncpy(cpar,name,l1); cpar[l1] = 0;
//   printf("igqwk called, wid=%d, pname=%s\n",wid,cpar);
   rval = higz->Get(cpar);
}
 
