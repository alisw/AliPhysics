///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Interface Class to the Geant3.21 MonteCarlo                              //
//                                                                           //
//Begin_Html
/*
<img src="picts/TGeant3Class.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TGeant3.h" 
#include "TROOT.h" 
#include "THIGZ.h" 
#include "ctype.h" 
#include <TDatabasePDG.h>
#include "AliCallf77.h" 
 
#ifndef WIN32 
# define gzebra  gzebra_ 
# define grfile  grfile_ 
# define gpcxyz  gpcxyz_ 
# define ggclos  ggclos_ 
# define glast   glast_ 
# define ginit   ginit_ 
# define gcinit  gcinit_ 
# define grun    grun_ 
# define gtrig   gtrig_ 
# define gtrigc  gtrigc_ 
# define gtrigi  gtrigi_ 
# define gwork   gwork_ 
# define gzinit  gzinit_ 
# define gfmate  gfmate_ 
# define gfpart  gfpart_ 
# define gftmed  gftmed_ 
# define gmate   gmate_ 
# define gpart   gpart_ 
# define gsdk    gsdk_ 
# define gsmate  gsmate_ 
# define gsmixt  gsmixt_ 
# define gspart  gspart_ 
# define gstmed  gstmed_ 
# define gsckov  gsckov_
# define gstpar  gstpar_ 
# define gfkine  gfkine_ 
# define gfvert  gfvert_ 
# define gskine  gskine_ 
# define gsvert  gsvert_ 
# define gphysi  gphysi_ 
# define gdebug  gdebug_ 
# define gekbin  gekbin_ 
# define gfinds  gfinds_ 
# define gsking  gsking_ 
# define gskpho  gskpho_ 
# define gsstak  gsstak_ 
# define gsxyz   gsxyz_ 
# define gtrack  gtrack_ 
# define gtreve  gtreve_ 
# define gtreve_root  gtreve_root_ 
# define grndm   grndm_ 
# define grndmq  grndmq_ 
# define gdtom   gdtom_ 
# define glmoth  glmoth_ 
# define gmedia  gmedia_ 
# define gmtod   gmtod_ 
# define gsdvn   gsdvn_ 
# define gsdvn2  gsdvn2_ 
# define gsdvs   gsdvs_ 
# define gsdvs2  gsdvs2_ 
# define gsdvt   gsdvt_ 
# define gsdvt2  gsdvt2_
# define gsord   gsord_ 
# define gspos   gspos_ 
# define gsposp  gsposp_ 
# define gsrotm  gsrotm_ 
# define gprotm  gprotm_ 
# define gsvolu  gsvolu_ 
# define gprint  gprint_ 
# define gdinit  gdinit_ 
# define gdopt   gdopt_ 
# define gdraw   gdraw_ 
# define gdrayt  gdrayt_
# define gdrawc  gdrawc_ 
# define gdrawx  gdrawx_ 
# define gdhead  gdhead_ 
# define gdwmn1  gdwmn1_ 
# define gdwmn2  gdwmn2_ 
# define gdwmn3  gdwmn3_ 
# define gdxyz   gdxyz_ 
# define gdcxyz  gdcxyz_ 
# define gdman   gdman_ 
# define gdspec  gdspec_ 
# define gdtree  gdtree_ 
# define gdelet  gdelet_ 
# define gdclos  gdclos_ 
# define gdshow  gdshow_ 
# define gdopen  gdopen_ 
# define dzshow  dzshow_ 
# define gsatt   gsatt_ 
# define gfpara  gfpara_
# define gckpar  gckpar_
# define gckmat  gckmat_
# define geditv  geditv_
# define mzdrop  mzdrop_

# define ertrak  ertrak_
# define ertrgo  ertrgo_
 
# define setbomb setbomb_
# define setclip setclip_
# define gcomad gcomad_

#else 
# define gzebra  GZEBRA 
# define grfile  GRFILE 
# define gpcxyz  GPCXYZ 
# define ggclos  GGCLOS 
# define glast   GLAST 
# define ginit   GINIT 
# define gcinit  GCINIT 
# define grun    GRUN 
# define gtrig   GTRIG 
# define gtrigc  GTRIGC 
# define gtrigi  GTRIGI 
# define gwork   GWORK 
# define gzinit  GZINIT 
# define gfmate  GFMATE 
# define gfpart  GFPART 
# define gftmed  GFTMED 
# define gmate   GMATE 
# define gpart   GPART 
# define gsdk    GSDK 
# define gsmate  GSMATE 
# define gsmixt  GSMIXT 
# define gspart  GSPART 
# define gstmed  GSTMED 
# define gsckov  GSCKOV
# define gstpar  GSTPAR 
# define gfkine  GFKINE 
# define gfvert  GFVERT 
# define gskine  GSKINE 
# define gsvert  GSVERT 
# define gphysi  GPHYSI 
# define gdebug  GDEBUG 
# define gekbin  GEKBIN 
# define gfinds  GFINDS 
# define gsking  GSKING 
# define gskpho  GSKPHO 
# define gsstak  GSSTAK 
# define gsxyz   GSXYZ 
# define gtrack  GTRACK 
# define gtreve  GTREVE 
# define gtreve_root  GTREVE_ROOT
# define grndm   GRNDM
# define grndmq  GRNDMQ
# define gdtom   GDTOM 
# define glmoth  GLMOTH 
# define gmedia  GMEDIA 
# define gmtod   GMTOD 
# define gsdvn   GSDVN 
# define gsdvn2  GSDVN2 
# define gsdvs   GSDVS 
# define gsdvs2  GSDVS2 
# define gsdvt   GSDVT 
# define gsdvt2  GSDVT2
# define gsord   GSORD 
# define gspos   GSPOS 
# define gsposp  GSPOSP 
# define gsrotm  GSROTM 
# define gprotm  GPROTM 
# define gsvolu  GSVOLU 
# define gprint  GPRINT 
# define gdinit  GDINIT
# define gdopt   GDOPT 
# define gdraw   GDRAW
# define gdrayt  GDRAYT
# define gdrawc  GDRAWC
# define gdrawx  GDRAWX 
# define gdhead  GDHEAD
# define gdwmn1  GDWMN1
# define gdwmn2  GDWMN2
# define gdwmn3  GDWMN3
# define gdxyz   GDXYZ
# define gdcxyz  GDCXYZ
# define gdman   GDMAN
# define gdfspc  GDFSPC
# define gdspec  GDSPEC
# define gdtree  GDTREE
# define gdelet  GDELET
# define gdclos  GDCLOS
# define gdshow  GDSHOW
# define gdopen  GDOPEN
# define dzshow  DZSHOW 
# define gsatt   GSATT 
# define gfpara  GFPARA
# define gckpar  GCKPAR
# define gckmat  GCKMAT
# define geditv  GEDITV
# define mzdrop  MZDROP 

# define ertrak  ERTRAK
# define ertrgo  ERTRGO
 
# define setbomb SETBOMB
# define setclip SETCLIP
# define gcomad  GCOMAD
 
#endif 

//____________________________________________________________________________ 
extern "C" 
{
  //
  // Prototypes for GEANT functions
  //
  void type_of_call gzebra(const int&); 

  void type_of_call gpcxyz(); 

  void type_of_call ggclos(); 

  void type_of_call glast(); 

  void type_of_call ginit(); 

  void type_of_call gcinit(); 

  void type_of_call grun(); 

  void type_of_call gtrig(); 

  void type_of_call gtrigc(); 

  void type_of_call gtrigi(); 

  void type_of_call gwork(const int&); 

  void type_of_call gzinit(); 

  void type_of_call gmate(); 

  void type_of_call gpart(); 

  void type_of_call gsdk(Int_t &, Float_t *, Int_t *); 

  void type_of_call gfkine(Int_t &, Float_t *, Float_t *, Int_t &,
			   Int_t &, Float_t *, Int_t &); 

  void type_of_call gfvert(Int_t &, Float_t *, Int_t &, Int_t &, 
			   Float_t &, Float_t *, Int_t &); 

  void type_of_call gskine(Float_t *,Int_t &, Int_t &, Float_t *,
			   Int_t &, Int_t &); 

  void type_of_call gsvert(Float_t *,Int_t &, Int_t &, Float_t *,
			   Int_t &, Int_t &); 

  void type_of_call gphysi(); 

  void type_of_call gdebug(); 

  void type_of_call gekbin(); 

  void type_of_call gfinds(); 

  void type_of_call gsking(Int_t &); 

  void type_of_call gskpho(Int_t &); 

  void type_of_call gsstak(Int_t &); 

  void type_of_call gsxyz(); 

  void type_of_call gtrack(); 

  void type_of_call gtreve(); 

  void type_of_call gtreve_root(); 

  void type_of_call grndm(Float_t *, const Int_t &); 

  void type_of_call grndmq(Int_t &, Int_t &, const Int_t &,
			   DEFCHARD DEFCHARL); 

  void type_of_call gdtom(Float_t *, Float_t *, Int_t &); 

  void type_of_call glmoth(DEFCHARD, Int_t &, Int_t &, Int_t *,
			   Int_t *, Int_t * DEFCHARL); 

  void type_of_call gmedia(Float_t *, Int_t &); 

  void type_of_call gmtod(Float_t *, Float_t *, Int_t &); 

  void type_of_call gsrotm(const Int_t &, const Float_t &, const Float_t &,
			   const Float_t &, const Float_t &, const Float_t &,
			   const Float_t &); 

  void type_of_call gprotm(const Int_t &); 

  void type_of_call grfile(const Int_t&, DEFCHARD, 
			   DEFCHARD DEFCHARL DEFCHARL); 

  void type_of_call gfmate(const Int_t&, DEFCHARD, Float_t &, Float_t &,
			   Float_t &, Float_t &, Float_t &, Float_t *,
			   Int_t& DEFCHARL); 

  void type_of_call gfpart(const Int_t&, DEFCHARD, Int_t &, Float_t &,
			   Float_t &, Float_t &, Float_t *, Int_t & DEFCHARL); 

  void type_of_call gftmed(const Int_t&, DEFCHARD, Int_t &, Int_t &, Int_t &,
			   Float_t &, Float_t &, Float_t &, Float_t &,
			   Float_t &, Float_t &, Float_t *, Int_t * DEFCHARL); 

  void type_of_call gsmate(const Int_t&, DEFCHARD, Float_t &, Float_t &,
			   Float_t &, Float_t &, Float_t &, Float_t *,
			   Int_t & DEFCHARL); 

  void type_of_call gsmixt(const Int_t&, DEFCHARD, Float_t *, Float_t *,
			   Float_t &, Int_t &, Float_t * DEFCHARL); 

  void type_of_call gspart(const Int_t&, DEFCHARD, Int_t &, Float_t &,
			   Float_t &, Float_t &, Float_t *, Int_t & DEFCHARL); 


  void type_of_call gstmed(const Int_t&, DEFCHARD, Int_t &, Int_t &, Int_t &,
			   Float_t &, Float_t &, Float_t &, Float_t &,
			   Float_t &, Float_t &, Float_t *, Int_t & DEFCHARL); 

  void type_of_call gsckov(Int_t &itmed, Int_t &npckov, Float_t *ppckov,
			   Float_t *absco, Float_t *effic, Float_t *rindex);
  void type_of_call gstpar(const Int_t&, DEFCHARD, Float_t & DEFCHARL); 

  void type_of_call gsdvn(DEFCHARD,DEFCHARD, Int_t &, Int_t &
			  DEFCHARL DEFCHARL); 

  void type_of_call gsdvn2(DEFCHARD,DEFCHARD, Int_t &, Int_t &, Float_t &,
			   Int_t & DEFCHARL DEFCHARL); 

  void type_of_call gsdvs(DEFCHARD,DEFCHARD, Float_t &, Int_t &, Int_t &
			  DEFCHARL DEFCHARL); 

  void type_of_call gsdvs2(DEFCHARD,DEFCHARD, Float_t &, Int_t &, Float_t &,
			   Int_t & DEFCHARL DEFCHARL); 

  void type_of_call gsdvt(DEFCHARD,DEFCHARD, Float_t &, Int_t &, Int_t &,
			  Int_t & DEFCHARL DEFCHARL); 

  void type_of_call gsdvt2(DEFCHARD,DEFCHARD, Float_t &, Int_t &, Float_t&,
			   Int_t &, Int_t & DEFCHARL DEFCHARL); 

  void type_of_call gsord(DEFCHARD, Int_t & DEFCHARL); 

  void type_of_call gspos(DEFCHARD, Int_t &, DEFCHARD, Float_t &, Float_t &,
			  Float_t &, Int_t &, DEFCHARD DEFCHARL DEFCHARL
			  DEFCHARL); 

  void type_of_call gsposp(DEFCHARD, Int_t &, DEFCHARD, Float_t &, Float_t &,
			   Float_t &, Int_t &, DEFCHARD,  
			   Float_t *, Int_t & DEFCHARL DEFCHARL DEFCHARL); 

  void type_of_call gsvolu(DEFCHARD, DEFCHARD, Int_t &, Float_t *, Int_t &,
			   Int_t & DEFCHARL DEFCHARL); 

  void type_of_call gsatt(DEFCHARD, DEFCHARD, Int_t & DEFCHARL DEFCHARL); 

  void type_of_call gfpara(DEFCHARD , Int_t&, Int_t&, Int_t&, Int_t&, Float_t*,
			   Float_t* DEFCHARL);

  void type_of_call gckpar(Int_t&, Int_t&, Float_t*);

  void type_of_call gckmat(Int_t&, DEFCHARD DEFCHARL);

  void type_of_call gprint(DEFCHARD,const int& DEFCHARL); 

  void type_of_call gdinit(); 

  void type_of_call gdopt(DEFCHARD,DEFCHARD DEFCHARL DEFCHARL); 
  
  void type_of_call gdraw(DEFCHARD,Float_t &,Float_t &, Float_t &,Float_t &,
			  Float_t &, Float_t &, Float_t & DEFCHARL); 
  void type_of_call gdrayt(DEFCHARD,Float_t &,Float_t &, Float_t &,Float_t &,
			   Float_t &, Float_t &, Float_t & DEFCHARL); 
  void type_of_call gdrawc(DEFCHARD,Int_t &, Float_t &, Float_t &, Float_t &,
                          Float_t &, Float_t & DEFCHARL); 
  void type_of_call gdrawx(DEFCHARD,Float_t &, Float_t &, Float_t &, Float_t &,
			   Float_t &, Float_t &, Float_t &, Float_t &,
			   Float_t & DEFCHARL); 
  void type_of_call gdhead(Int_t &,DEFCHARD, Float_t & DEFCHARL);
  void type_of_call gdxyz(Int_t &);
  void type_of_call gdcxyz();
  void type_of_call gdman(Float_t &, Float_t &);
  void type_of_call gdwmn1(Float_t &, Float_t &);
  void type_of_call gdwmn2(Float_t &, Float_t &);
  void type_of_call gdwmn3(Float_t &, Float_t &);
  void type_of_call gdspec(DEFCHARD DEFCHARL);
  void type_of_call gdfspc(DEFCHARD, Int_t &, Int_t & DEFCHARL) {;}
  void type_of_call gdtree(DEFCHARD, Int_t &, Int_t & DEFCHARL);

  void type_of_call gdopen(Int_t &);
  void type_of_call gdclos();
  void type_of_call gdelet(Int_t &);
  void type_of_call gdshow(Int_t &);
  void type_of_call geditv(Int_t &) {;}


  void type_of_call dzshow(DEFCHARD,const int&,const int&,DEFCHARD,const int&,
			   const int&, const int&, const int& DEFCHARL
			   DEFCHARL); 

  void type_of_call mzdrop(Int_t&, Int_t&, DEFCHARD DEFCHARL);

  void type_of_call setbomb(Float_t &);
  void type_of_call setclip(DEFCHARD, Float_t &,Float_t &,Float_t &,Float_t &,
			    Float_t &, Float_t & DEFCHARL); 
  void type_of_call gcomad(DEFCHARD, Int_t*& DEFCHARL); 

  void type_of_call ertrak(const Float_t *const x1, const Float_t *const p1,
			   const Float_t *x2, const Float_t *p2,
			   const Int_t &ipa, DEFCHARD DEFCHARL);

  void type_of_call ertrgo();
}

//
// Geant3 global pointer
//
static Int_t defSize = 600;

ClassImp(TGeant3) 
 
//____________________________________________________________________________ 
TGeant3::TGeant3()
{ 
  //
  // Default constructor
  //
} 
 
//____________________________________________________________________________ 
TGeant3::TGeant3(const char *title, Int_t nwgeant) 
       :AliMC("TGeant3",title) 
{
  //
  // Standard constructor for TGeant3 with ZEBRA initialisation
  // 
   
  if(nwgeant) {
    gzebra(nwgeant); 
    ginit(); 
    gzinit();
  } else {
    gcinit();
  }
  //
  // Load Address of Geant3 commons    
  LoadAddress(); 
  //
  // Zero number of particles
  fNPDGCodes=0;
} 

//____________________________________________________________________________ 
Int_t TGeant3::CurrentMaterial(Float_t &a, Float_t &z, Float_t &dens,
			       Float_t &radl, Float_t &absl) const
{
  //
  // Return the parameters of the current material during transport
  //
  z     = fGcmate->z;
  a     = fGcmate->a;
  dens  = fGcmate->dens;
  radl  = fGcmate->radl;
  absl  = fGcmate->absl;
  return 1;  //this could be the number of elements in mixture
}
   
//____________________________________________________________________________ 
void TGeant3::DefaultRange()
{ 
  //
  // Set range of current drawing pad to 20x20 cm
  //
  if (!higz) {
    new THIGZ(defSize); 
    gdinit();
  }
  higz->Range(0,0,20,20);
}

//____________________________________________________________________________ 
void TGeant3::InitHIGZ() 
{ 
  //
  // Initialise HIGZ
  //
  if (!higz) {
    new THIGZ(defSize); 
    gdinit();
  }
}
 
//____________________________________________________________________________ 
void TGeant3::LoadAddress() 
{ 
  //
  // Assigns the address of the GEANT common blocks to the structures
  // that allow their access from C++
  //
  Int_t *addr;
  gcomad(PASSCHARD("QUEST"), (int*&) fQuest PASSCHARL("QUEST"));
  gcomad(PASSCHARD("GCBANK"),(int*&) fGcbank  PASSCHARL("GCBANK"));
  gcomad(PASSCHARD("GCLINK"),(int*&) fGclink  PASSCHARL("GCLINK"));
  gcomad(PASSCHARD("GCCUTS"),(int*&) fGccuts  PASSCHARL("GCCUTS"));
  gcomad(PASSCHARD("GCFLAG"),(int*&) fGcflag  PASSCHARL("GCFLAG"));
  gcomad(PASSCHARD("GCKINE"),(int*&) fGckine  PASSCHARL("GCKINE"));
  gcomad(PASSCHARD("GCKING"),(int*&) fGcking  PASSCHARL("GCKING"));
  gcomad(PASSCHARD("GCKIN2"),(int*&) fGckin2  PASSCHARL("GCKIN2"));
  gcomad(PASSCHARD("GCKIN3"),(int*&) fGckin3  PASSCHARL("GCKIN3"));
  gcomad(PASSCHARD("GCMATE"),(int*&) fGcmate  PASSCHARL("GCMATE"));
  gcomad(PASSCHARD("GCTMED"),(int*&) fGctmed  PASSCHARL("GCTMED"));
  gcomad(PASSCHARD("GCTRAK"),(int*&) fGctrak  PASSCHARL("GCTRAK"));
  gcomad(PASSCHARD("GCTPOL"),(int*&) fGctpol  PASSCHARL("GCTPOL"));
  gcomad(PASSCHARD("GCVOLU"),(int*&) fGcvolu  PASSCHARL("GCVOLU"));
  gcomad(PASSCHARD("GCNUM"), (int*&) fGcnum   PASSCHARL("GCNUM"));
  gcomad(PASSCHARD("GCSETS"),(int*&) fGcsets  PASSCHARL("GCSETS"));
  gcomad(PASSCHARD("GCPHYS"),(int*&) fGcphys  PASSCHARL("GCPHYS"));
  gcomad(PASSCHARD("GCOPTI"),(int*&) fGcopti  PASSCHARL("GCOPTI"));
  gcomad(PASSCHARD("GCTLIT"),(int*&) fGctlit  PASSCHARL("GCTLIT"));
  gcomad(PASSCHARD("GCVDMA"),(int*&) fGcvdma  PASSCHARL("GCVDMA"));

  gcomad(PASSCHARD("IQ"), addr  PASSCHARL("IQ"));
  fZiq = addr;
  gcomad(PASSCHARD("LQ"), addr  PASSCHARL("LQ"));
  fZlq = addr;
  fZq       = (float*)fZiq; 
} 

//_____________________________________________________________________________
void TGeant3::GeomIter()
{
  //
  // Geometry iterator for moving upward in the geometry tree
  // Initialise the iterator
  //
  fNextVol=fGcvolu->nlevel;
}

//____________________________________________________________________________ 
Int_t TGeant3::NextVolUp(Text_t *name, Int_t &copy)
{
  //
  // Geometry iterator for moving upward in the geometry tree
  // Return next volume up
  //
  Int_t i, gname;
  fNextVol--;
  if(fNextVol>=0) {
    gname=fGcvolu->names[fNextVol];
    strncpy(name,(char *) &gname, 4);
    name[4]='\0';
    copy=fGcvolu->number[fNextVol];
    i=fGcvolu->lvolum[fNextVol];
    if(gname == fZiq[fGclink->jvolum+i]) return i;
    else printf("GeomTree: Volume %s not found in bank\n",name);
  }
  return 0;
}

//_____________________________________________________________________________
Int_t TGeant3::CurrentVolID(Int_t &copy) const
{
  //
  // Returns the current volume ID and copy number
  //
  Int_t i, gname;
  if( (i=fGcvolu->nlevel-1) < 0 ) {
    Warning("CurrentVolID","Stack depth only %d\n",fGcvolu->nlevel);
  } else {
    gname=fGcvolu->names[i];
    copy=fGcvolu->number[i];
    i=fGcvolu->lvolum[i];   
    if(gname == fZiq[fGclink->jvolum+i]) return i;
    else Warning("CurrentVolID","Volume %4s not found\n",(char*)&gname);
  }
  return 0;
}

//_____________________________________________________________________________
Int_t TGeant3::CurrentVolOffID(Int_t off, Int_t &copy) const
{
  //
  // Return the current volume "off" upward in the geometrical tree 
  // ID and copy number
  //
  Int_t i, gname;
  if( (i=fGcvolu->nlevel-off-1) < 0 ) {
    Warning("CurrentVolOffID","Offset requested %d but stack depth %d\n",
	    off,fGcvolu->nlevel);
  } else {
    gname=fGcvolu->names[i];
    copy=fGcvolu->number[i];          
    i=fGcvolu->lvolum[i];    
    if(gname == fZiq[fGclink->jvolum+i]) return i;
    else Warning("CurrentVolOffID","Volume %4s not found\n",(char*)&gname);
  }
  return 0;
}

//_____________________________________________________________________________
const char* TGeant3::CurrentVolName() const
{
  //
  // Returns the current volume name
  //
  Int_t i, gname;
  char *name;
  if( (i=fGcvolu->nlevel-1) < 0 ) {
    Warning("CurrentVolName","Stack depth %d\n",fGcvolu->nlevel);
  } else {
    gname=fGcvolu->names[i];
    name = new char[5];
    strncpy(name,(char *) &gname, 4);
    name[4]='\0';
    i=fGcvolu->lvolum[i];   
    if(gname == fZiq[fGclink->jvolum+i]) return name;
    else Warning("CurrentVolName","Volume %4s not found\n",name);
  }
  return 0;
}

//_____________________________________________________________________________
const char* TGeant3::CurrentVolOffName(Int_t off) const
{
  //
  // Return the current volume "off" upward in the geometrical tree 
  // ID, name and copy number
  // if name=0 no name is returned
  //
  Int_t i, gname;
  char *name;
  if( (i=fGcvolu->nlevel-off-1) < 0 ) {
    Warning("CurrentVolOffName",
	    "Offset requested %d but stack depth %d\n",off,fGcvolu->nlevel);
  } else {
    gname=fGcvolu->names[i];
    name = new char[5];
    strncpy(name,(char *) &gname, 4);
    name[4]='\0';
    i=fGcvolu->lvolum[i];    
    if(gname == fZiq[fGclink->jvolum+i]) return name;
    else Warning("CurrentVolOffName","Volume %4s not found\n",name);
  }
  return 0;
}

//_____________________________________________________________________________
Int_t TGeant3::IdFromPDG(Int_t pdg) const 
{
  //
  // Return Geant3 code from PDG and pseudo ENDF code

  for(Int_t i=0;i<fNPDGCodes;++i)
    if(pdg==fPDGCode[i]) return i;
  return -1;
}

//_____________________________________________________________________________
Int_t TGeant3::PDGFromId(Int_t id) const 
{
  if(id>0 && id<fNPDGCodes) return fPDGCode[id];
  else return -1;
}

//_____________________________________________________________________________
void TGeant3::DefineParticles() 
{
  //
  // Define standard Geant 3 particles
  Gpart();
  //
  // Load standard numbers for GEANT particles and PDG conversion
  fPDGCode[fNPDGCodes++]=-99;   //  0 = unused location
  fPDGCode[fNPDGCodes++]=22;    //  1 = photon
  fPDGCode[fNPDGCodes++]=-11;   //  2 = positron
  fPDGCode[fNPDGCodes++]=11;    //  3 = electron
  fPDGCode[fNPDGCodes++]=12;    //  4 = neutrino e
  fPDGCode[fNPDGCodes++]=-13;   //  5 = muon +
  fPDGCode[fNPDGCodes++]=13;    //  6 = muon -
  fPDGCode[fNPDGCodes++]=111;   //  7 = pi0
  fPDGCode[fNPDGCodes++]=211;   //  8 = pi+
  fPDGCode[fNPDGCodes++]=-211;  //  9 = pi-
  fPDGCode[fNPDGCodes++]=130;   // 10 = Kaon Long
  fPDGCode[fNPDGCodes++]=321;   // 11 = Kaon +
  fPDGCode[fNPDGCodes++]=-321;  // 12 = Kaon -
  fPDGCode[fNPDGCodes++]=2112;  // 13 = Neutron
  fPDGCode[fNPDGCodes++]=2212;  // 14 = Proton
  fPDGCode[fNPDGCodes++]=-2212; // 15 = Anti Proton
  fPDGCode[fNPDGCodes++]=310;   // 16 = Kaon Short
  fPDGCode[fNPDGCodes++]=221;   // 17 = Eta
  fPDGCode[fNPDGCodes++]=3122;  // 18 = Lambda
  fPDGCode[fNPDGCodes++]=3222;  // 19 = Sigma +
  fPDGCode[fNPDGCodes++]=3212;  // 20 = Sigma 0
  fPDGCode[fNPDGCodes++]=3112;  // 21 = Sigma -
  fPDGCode[fNPDGCodes++]=3322;  // 22 = Xi0
  fPDGCode[fNPDGCodes++]=3312;  // 23 = Xi-
  fPDGCode[fNPDGCodes++]=3334;  // 24 = Omega-
  fPDGCode[fNPDGCodes++]=-2112; // 25 = Anti Proton
  fPDGCode[fNPDGCodes++]=-3122; // 26 = Anti Proton
  fPDGCode[fNPDGCodes++]=-3222; // 27 = Anti Sigma -
  fPDGCode[fNPDGCodes++]=-3212; // 28 = Anti Sigma 0
  fPDGCode[fNPDGCodes++]=-3112; // 29 = Anti Sigma 0
  fPDGCode[fNPDGCodes++]=-3322; // 30 = Anti Xi 0
  fPDGCode[fNPDGCodes++]=-3312; // 31 = Anti Xi +
  fPDGCode[fNPDGCodes++]=-3334; // 32 = Anti Omega +


  Int_t mode[6];
  Int_t kz, ipa;
  Float_t bratio[6];

  /* --- Define additional particles */
  Gspart(33, "OMEGA(782)", 3, 0.782, 0., 7.836e-23);
  fPDGCode[fNPDGCodes++]=223;   // 33 = Omega(782)
  
  Gspart(34, "PHI(1020)", 3, 1.019, 0., 1.486e-22);
  fPDGCode[fNPDGCodes++]=333;   // 34 = PHI (1020)

  Gspart(35, "D +", 4, 1.87, 1., 1.066e-12);
  fPDGCode[fNPDGCodes++]=411;   // 35 = D+

  Gspart(36, "D -", 4, 1.87, -1., 1.066e-12);
  fPDGCode[fNPDGCodes++]=-411;  // 36 = D-

  Gspart(37, "D 0", 3, 1.865, 0., 4.2e-13);
  fPDGCode[fNPDGCodes++]=421;   // 37 = D0

  Gspart(38, "ANTI D 0", 3, 1.865, 0., 4.2e-13);
  fPDGCode[fNPDGCodes++]=-421;  // 38 = D0 bar

  fPDGCode[fNPDGCodes++]=-99;  // 39 = unassigned

  fPDGCode[fNPDGCodes++]=-99;  // 40 = unassigned

  fPDGCode[fNPDGCodes++]=-99;  // 41 = unassigned

  Gspart(42, "RHO +", 4, 0.768, 1., 4.353e-24);
  fPDGCode[fNPDGCodes++]=213;   // 42 = RHO+

  Gspart(43, "RHO -", 4, 0.768, -1., 4.353e-24);
  fPDGCode[fNPDGCodes++]=-213;   // 40 = RHO-

  Gspart(44, "RHO 0", 3, 0.768, 0., 4.353e-24);
  fPDGCode[fNPDGCodes++]=113;   // 37 = D0

  //
  // Use ENDF-6 mapping for ions = 10000*z+10*a+iso
  // and add 1 000 000
  // and numbers above 5 000 000 for special applications
  //

  const Int_t kion=10000000;

  const Int_t kspe=50000000;

  TDatabasePDG *pdgDB = TDatabasePDG::Instance();

  const Double_t autogev=0.9314943228;
  const Double_t hslash = 1.0545726663e-27;
  const Double_t erggev = 1/1.6021773349e-3;
  const Double_t hshgev = hslash*erggev;
  const Double_t yearstosec = 3600*24*365.25;


  pdgDB->AddParticle("Deuteron","Deuteron",2*autogev+8.071e-3,kTRUE,
		     0,1,"Ion",kion+10020);
  fPDGCode[fNPDGCodes++]=kion+10020;   // 45 = Deuteron

  pdgDB->AddParticle("Triton","Triton",3*autogev+14.931e-3,kFALSE,
		     hshgev/(12.33*yearstosec),1,"Ion",kion+10030);
  fPDGCode[fNPDGCodes++]=kion+10030;   // 46 = Triton

  pdgDB->AddParticle("Alpha","Alpha",4*autogev+2.424e-3,kTRUE,
		     hshgev/(12.33*yearstosec),2,"Ion",kion+20040);
  fPDGCode[fNPDGCodes++]=kion+20040;   // 47 = Alpha

  fPDGCode[fNPDGCodes++]=0;   // 48 = geantino mapped to rootino

  pdgDB->AddParticle("HE3","HE3",3*autogev+14.931e-3,kFALSE,
		     0,2,"Ion",kion+20030);
  fPDGCode[fNPDGCodes++]=kion+20030;   // 49 = HE3

  pdgDB->AddParticle("Cherenkov","Cherenkov",0,kFALSE,
		     0,0,"Special",kspe+50);
  fPDGCode[fNPDGCodes++]=kspe+50;   // 50 = Cherenkov

/* --- Define additional decay modes --- */
/* --- omega(783) --- */
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 33;
    bratio[0] = 89.;
    bratio[1] = 8.5;
    bratio[2] = 2.5;
    mode[0] = 70809;
    mode[1] = 107;
    mode[2] = 908;
    Gsdk(ipa, bratio, mode);
/* --- phi(1020) --- */
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 34;
    bratio[0] = 49.;
    bratio[1] = 34.4;
    bratio[2] = 12.9;
    bratio[3] = 2.4;
    bratio[4] = 1.3;
    mode[0] = 1112;
    mode[1] = 1610;
    mode[2] = 4407;
    mode[3] = 90807;
    mode[4] = 1701;
    Gsdk(ipa, bratio, mode);
/* --- D+ --- */
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 35;
    bratio[0] = 25.;
    bratio[1] = 25.;
    bratio[2] = 25.;
    bratio[3] = 25.;
    mode[0] = 80809;
    mode[1] = 120808;
    mode[2] = 111208;
    mode[3] = 110809;
    Gsdk(ipa, bratio, mode);
/* --- D- --- */
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 36;
    bratio[0] = 25.;
    bratio[1] = 25.;
    bratio[2] = 25.;
    bratio[3] = 25.;
    mode[0] = 90908;
    mode[1] = 110909;
    mode[2] = 121109;
    mode[3] = 120908;
    Gsdk(ipa, bratio, mode);
/* --- D0 --- */
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 37;
    bratio[0] = 33.;
    bratio[1] = 33.;
    bratio[2] = 33.;
    mode[0] = 809;
    mode[1] = 1208;
    mode[2] = 1112;
    Gsdk(ipa, bratio, mode);
/* --- Anti D0 --- */
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 38;
    bratio[0] = 33.;
    bratio[1] = 33.;
    bratio[2] = 33.;
    mode[0] = 809;
    mode[1] = 1109;
    mode[2] = 1112;
    Gsdk(ipa, bratio, mode);
/* --- rho+ --- */
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 42;
    bratio[0] = 100.;
    mode[0] = 807;
    Gsdk(ipa, bratio, mode);
/* --- rho- --- */
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 43;
    bratio[0] = 100.;
    mode[0] = 907;
    Gsdk(ipa, bratio, mode);
/* --- rho0 --- */
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 44;
    bratio[0] = 100.;
    mode[0] = 707;
    Gsdk(ipa, bratio, mode);
    /*
// --- jpsi ---
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 113;
    bratio[0] = 50.;
    bratio[1] = 50.;
    mode[0] = 506;
    mode[1] = 605;
    Gsdk(ipa, bratio, mode);
// --- upsilon --- 
    ipa = 114;
    Gsdk(ipa, bratio, mode);
// --- phi --- 
    ipa = 115;
    Gsdk(ipa, bratio, mode);
    */

}

//_____________________________________________________________________________
Int_t TGeant3::VolId(Text_t *name) const
{
  //
  // Return the unique numeric identifier for volume name
  //
  Int_t gname, i;
  strncpy((char *) &gname, name, 4);
  for(i=1; i<=fGcnum->nvolum; i++)
    if(gname == fZiq[fGclink->jvolum+i]) return i;
  printf("VolId: Volume %s not found\n",name);
  return 0;
}

//_____________________________________________________________________________
Int_t TGeant3::NofVolumes() const 
{
  //
  // Return total number of volumes in the geometry
  //
  return fGcnum->nvolum;
}

//_____________________________________________________________________________
const char* TGeant3::VolName(Int_t id) const
{
  //
  // Return the volume name given the volume identifier
  //
  static char name[5];
  if(id<1 || id > fGcnum->nvolum || fGclink->jvolum<=0) 
    strcpy(name,"NULL");
  else
    strncpy(name,(char *)&fZiq[fGclink->jvolum+id],4);
  name[4]='\0';
  return name;
}

//_____________________________________________________________________________
void TGeant3::TrackPosition(TLorentzVector &xyz) const
{
  //
  // Return the current position in the master reference frame of the
  // track being transported
  //
  xyz[0]=fGctrak->vect[0];
  xyz[1]=fGctrak->vect[1];
  xyz[2]=fGctrak->vect[2];
  xyz[3]=fGctrak->tofg;
}

//_____________________________________________________________________________
Float_t TGeant3::TrackTime() const
{
  //
  // Return the current time of flight of the track being transported
  //
  return fGctrak->tofg;
}

//_____________________________________________________________________________
void TGeant3::TrackMomentum(TLorentzVector &xyz) const
{
  //
  // Return the direction and the momentum (GeV/c) of the track
  // currently being transported
  //
  Double_t ptot=fGctrak->vect[6];
  xyz[0]=fGctrak->vect[3]*ptot;
  xyz[1]=fGctrak->vect[4]*ptot;
  xyz[2]=fGctrak->vect[5]*ptot;
  xyz[3]=fGctrak->getot;
}

//_____________________________________________________________________________
Float_t TGeant3::TrackCharge() const
{
  //
  // Return charge of the track currently transported
  //
  return fGckine->charge;
}

//_____________________________________________________________________________
Float_t TGeant3::TrackMass() const
{
  //
  // Return the mass of the track currently transported
  //
  return fGckine->amass;
}

//_____________________________________________________________________________
Int_t TGeant3::TrackPid() const
{
  //
  // Return the id of the particle transported
  //
  return PDGFromId(fGckine->ipart);
}

//_____________________________________________________________________________
Float_t TGeant3::TrackStep() const
{
  //
  // Return the length in centimeters of the current step
  //
  return fGctrak->step;
}

//_____________________________________________________________________________
Float_t TGeant3::TrackLength() const
{
  //
  // Return the length of the current track from its origin
  //
  return fGctrak->sleng;
}

//_____________________________________________________________________________
Bool_t TGeant3::IsTrackInside() const
{
  //
  // True if the track is not at the boundary of the current volume
  //
  return (fGctrak->inwvol==0);
}

//_____________________________________________________________________________
Bool_t TGeant3::IsTrackEntering() const
{
  //
  // True if this is the first step of the track in the current volume
  //
  return (fGctrak->inwvol==1);
}

//_____________________________________________________________________________
Bool_t TGeant3::IsTrackExiting() const
{
  //
  // True if this is the last step of the track in the current volume
  //
  return (fGctrak->inwvol==2);
}

//_____________________________________________________________________________
Bool_t TGeant3::IsTrackOut() const
{
  //
  // True if the track is out of the setup
  //
  return (fGctrak->inwvol==3);
}

//_____________________________________________________________________________
Bool_t TGeant3::IsTrackStop() const
{
  //
  // True if the track energy has fallen below the threshold 
  //
  return (fGctrak->istop==2);
}

//_____________________________________________________________________________
Int_t   TGeant3::NSecondaries() const
{
  //
  // Number of secondary particles generated in the current step
  //
  return fGcking->ngkine;
}

//_____________________________________________________________________________
Int_t   TGeant3::CurrentEvent() const
{
  //
  // Number of the current event
  //
  return fGcflag->idevt;
}

//_____________________________________________________________________________
void    TGeant3::ProdProcess(char* proc) const
{
  //
  // Name of the process that has produced the secondary particles
  // in the current step
  //
  const Int_t ipmec[13] = { 5,6,7,8,9,10,11,12,21,23,25,105,108 };
  Int_t mec, km, im;
  //
  if(fGcking->ngkine>0) {
    for (km = 0; km < fGctrak->nmec; ++km) {
      for (im = 0; im < 13; ++im) {
	if (fGctrak->lmec[km] == ipmec[im]) {
	  mec = fGctrak->lmec[km];
	  if (0 < mec && mec < 31) {
	    strncpy(proc,(char *)&fGctrak->namec[mec - 1],4);
	  } else if (mec - 100 <= 30 && mec - 100 > 0) {
	    strncpy(proc,(char *)&fGctpol->namec1[mec - 101],4);
	  }
	  proc[4]='\0';
	  return;
	}
      }
    }
    strcpy(proc,"UNKN");
  } else strcpy(proc,"NONE");
}

//_____________________________________________________________________________
void    TGeant3::GetSecondary(Int_t isec, Int_t& ipart, Float_t* x, Float_t* p)
{
  //
  // Get the parameters of the secondary track number isec produced
  // in the current step
  //
  Int_t i;
  if(-1<isec && isec<fGcking->ngkine) {
    ipart=Int_t (fGcking->gkin[isec][4] +0.5);
    for(i=0;i<3;i++) {
      x[i]=fGckin3->gpos[isec][i];
      p[i]=fGcking->gkin[isec][i];
    }
    x[3]=fGcking->tofd[isec];
    p[3]=fGcking->gkin[isec][3];
  } else {
    printf(" * TGeant3::GetSecondary * Secondary %d does not exist\n",isec);
    x[0]=x[1]=x[2]=x[3]=p[0]=p[1]=p[2]=p[3]=0;
    ipart=0;
  }
}

//_____________________________________________________________________________
void TGeant3::InitLego()
{
  SetSWIT(4,0);
  SetDEBU(0,0,0);  //do not print a message 
}

//_____________________________________________________________________________
Bool_t TGeant3::IsTrackDisappeared() const
{
  //
  // True if the current particle has disappered
  // either because it decayed or because it underwent
  // an inelastic collision
  //
  return (fGctrak->istop==1);
}

//_____________________________________________________________________________
Bool_t TGeant3::IsTrackAlive() const
{
  //
  // True if the current particle is alive and will continue to be
  // transported
  //
  return (fGctrak->istop==0);
}

//_____________________________________________________________________________
void TGeant3::StopTrack()
{
  //
  // Stop the transport of the current particle and skip to the next
  //
  fGctrak->istop=1;
}

//_____________________________________________________________________________
void TGeant3::StopEvent()
{
  //
  // Stop simulation of the current event and skip to the next
  //
  fGcflag->ieotri=1;
}

//_____________________________________________________________________________
Float_t TGeant3::MaxStep() const
{
  //
  // Return the maximum step length in the current medium
  //
  return fGctmed->stemax;
}

//_____________________________________________________________________________
void TGeant3::SetColors()
{
  //
  // Set the colors for all the volumes
  // this is done sequentially for all volumes
  // based on the number of their medium
  //
  Int_t kv, icol;
  Int_t jvolum=fGclink->jvolum;
  //Int_t jtmed=fGclink->jtmed;
  //Int_t jmate=fGclink->jmate;
  Int_t nvolum=fGcnum->nvolum;
  char name[5];
  //
  //    Now for all the volumes
  for(kv=1;kv<=nvolum;kv++) {
    //     Get the tracking medium
    Int_t itm=Int_t (fZq[fZlq[jvolum-kv]+4]);
    //     Get the material
    //Int_t ima=Int_t (fZq[fZlq[jtmed-itm]+6]);
    //     Get z
    //Float_t z=fZq[fZlq[jmate-ima]+7];
    //     Find color number
    //icol = Int_t(z)%6+2;
    //icol = 17+Int_t(z*150./92.);
    //icol = kv%6+2;
    icol = itm%6+2;
    strncpy(name,(char*)&fZiq[jvolum+kv],4);
    name[4]='\0';
    Gsatt(name,"COLO",icol);
  }
}

//_____________________________________________________________________________
void TGeant3::SetMaxStep(Float_t maxstep)
{
  //
  // Set the maximum step allowed till the particle is in the current medium
  //
  fGctmed->stemax=maxstep;
}

//_____________________________________________________________________________
void TGeant3::SetMaxNStep(Int_t maxnstp)
{
  //
  // Set the maximum number of steps till the particle is in the current medium
  //
  fGctrak->maxnst=maxnstp;
}

//_____________________________________________________________________________
Int_t TGeant3::GetMaxNStep() const
{
  //
  // Maximum number of steps allowed in current medium
  //
  return fGctrak->maxnst;
}

//_____________________________________________________________________________
void TGeant3::Material(Int_t& kmat, const char* name, Float_t a, Float_t z,
		       Float_t dens, Float_t radl, Float_t absl, Float_t* buf,
		       Int_t nwbuf)
{
  //
  // Defines a Material
  // 
  //  kmat               number assigned to the material
  //  name               material name
  //  a                  atomic mass in au
  //  z                  atomic number
  //  dens               density in g/cm3
  //  absl               absorbtion length in cm
  //                     if >=0 it is ignored and the program 
  //                     calculates it, if <0. -absl is taken
  //  radl               radiation length in cm
  //                     if >=0 it is ignored and the program 
  //                     calculates it, if <0. -radl is taken
  //  buf                pointer to an array of user words
  //  nbuf               number of user words
  //
  Int_t jmate=fGclink->jmate;
  kmat=1;
  Int_t ns, i;
  if(jmate>0) {
    ns=fZiq[jmate-2];
    kmat=ns+1;
    for(i=1; i<=ns; i++) {
      if(fZlq[jmate-i]==0) {
	kmat=i;
	break;
      }
    }
  }
  gsmate(kmat,PASSCHARD(name), a, z, dens, radl, absl, buf,
	 nwbuf PASSCHARL(name)); 
}

//_____________________________________________________________________________
void TGeant3::Mixture(Int_t& kmat, const char* name, Float_t* a, Float_t* z, 
		      Float_t dens, Int_t nlmat, Float_t* wmat)
{
  //
  // Defines mixture OR COMPOUND IMAT as composed by 
  // THE BASIC NLMAT materials defined by arrays A,Z and WMAT
  // 
  // If NLMAT > 0 then wmat contains the proportion by
  // weights of each basic material in the mixture. 
  // 
  // If nlmat < 0 then WMAT contains the number of atoms 
  // of a given kind into the molecule of the COMPOUND
  // In this case, WMAT in output is changed to relative
  // weigths.
  //
  Int_t jmate=fGclink->jmate;
  kmat=1;
  Int_t ns, i;
  if(jmate>0) {
    ns=fZiq[jmate-2];
    kmat=ns+1;
    for(i=1; i<=ns; i++) {
      if(fZlq[jmate-i]==0) {
	kmat=i;
	break;
      }
    }
  }
  gsmixt(kmat,PASSCHARD(name), a, z,dens, nlmat,wmat PASSCHARL(name)); 
}

//_____________________________________________________________________________
void TGeant3::Medium(Int_t& kmed, const char* name, Int_t nmat, Int_t isvol,
		     Int_t ifield, Float_t fieldm, Float_t tmaxfd,
		     Float_t stemax, Float_t deemax, Float_t epsil,
		     Float_t stmin, Float_t* ubuf, Int_t nbuf)
{
  //
  //  kmed      tracking medium number assigned
  //  name      tracking medium name
  //  nmat      material number
  //  isvol     sensitive volume flag
  //  ifield    magnetic field
  //  fieldm    max. field value (kilogauss)
  //  tmaxfd    max. angle due to field (deg/step)
  //  stemax    max. step allowed
  //  deemax    max. fraction of energy lost in a step
  //  epsil     tracking precision (cm)
  //  stmin     min. step due to continuos processes (cm)
  //
  //  ifield = 0 if no magnetic field; ifield = -1 if user decision in guswim;
  //  ifield = 1 if tracking performed with grkuta; ifield = 2 if tracking
  //  performed with ghelix; ifield = 3 if tracking performed with ghelx3.
  //  
  Int_t jtmed=fGclink->jtmed;
  kmed=1;
  Int_t ns, i;
  if(jtmed>0) {
    ns=fZiq[jtmed-2];
    kmed=ns+1;
    for(i=1; i<=ns; i++) {
      if(fZlq[jtmed-i]==0) {
	kmed=i;
	break;
      }
    }
  }
  gstmed(kmed, PASSCHARD(name), nmat, isvol, ifield, fieldm, tmaxfd, stemax,
	 deemax, epsil, stmin, ubuf, nbuf PASSCHARL(name)); 
}

//_____________________________________________________________________________
void TGeant3::Matrix(Int_t& krot, Float_t thex, Float_t phix, Float_t they,
		     Float_t phiy, Float_t thez, Float_t phiz)
{
  //
  //  krot     rotation matrix number assigned
  //  theta1   polar angle for axis i
  //  phi1     azimuthal angle for axis i
  //  theta2   polar angle for axis ii
  //  phi2     azimuthal angle for axis ii
  //  theta3   polar angle for axis iii
  //  phi3     azimuthal angle for axis iii
  //
  //  it defines the rotation matrix number irot.
  //  
  Int_t jrotm=fGclink->jrotm;
  krot=1;
  Int_t ns, i;
  if(jrotm>0) {
    ns=fZiq[jrotm-2];
    krot=ns+1;
    for(i=1; i<=ns; i++) {
      if(fZlq[jrotm-i]==0) {
	krot=i;
	break;
      }
    }
  }
  gsrotm(krot, thex, phix, they, phiy, thez, phiz);
}

//_____________________________________________________________________________
Int_t TGeant3::GetMedium() const
{
  //
  // Return the number of the current medium
  //
  return fGctmed->numed;
}

//_____________________________________________________________________________
Float_t TGeant3::Edep() const
{
  //
  // Return the energy lost in the current step
  //
  return fGctrak->destep;
}

//_____________________________________________________________________________
Float_t TGeant3::Etot() const
{
  //
  // Return the total energy of the current track
  //
  return fGctrak->getot;
}

//_____________________________________________________________________________
void TGeant3::Rndm(Float_t* r, const Int_t n) const
{
  //
  // Return an array of n random numbers uniformly distributed 
  // between 0 and 1 not included
  //
  Grndm(r,n);
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GBASE
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//____________________________________________________________________________ 
void  TGeant3::Gfile(const char *filename, const char *option) 
{ 
  //
  //    Routine to open a GEANT/RZ data base. 
  //
  //    LUN logical unit number associated to the file 
  //
  //    CHFILE RZ file name   
  //  
  //    CHOPT is a character string which may be  
  //        N  To create a new file 
  //        U  to open an existing file for update 
  //       " " to open an existing file for read only
  //        Q  The initial allocation (default 1000 records) 
  //           is given in IQUEST(10)
  //        X  Open the file in exchange format
  //        I  Read all data structures from file to memory 
  //        O  Write all data structures from memory to file 
  // 
  // Note:
  //      If options "I"  or "O" all data structures are read or
  //         written from/to file and the file is closed. 
  //      See routine GRMDIR to create subdirectories  
  //      See routines GROUT,GRIN to write,read objects 
  //  
  grfile(21, PASSCHARD(filename), PASSCHARD(option) PASSCHARL(filename)
	 PASSCHARL(option)); 
} 
 
//____________________________________________________________________________ 
void  TGeant3::Gpcxyz() 
{ 
  //
  //    Print track and volume parameters at current point
  //
  gpcxyz(); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Ggclos() 
{ 
  //
  //   Closes off the geometry setting.
  //   Initializes the search list for the contents of each
  //   volume following the order they have been positioned, and
  //   inserting the content '0' when a call to GSNEXT (-1) has
  //   been required by the user.
  //   Performs the development of the JVOLUM structure for all 
  //   volumes with variable parameters, by calling GGDVLP. 
  //   Interprets the user calls to GSORD, through GGORD.
  //   Computes and stores in a bank (next to JVOLUM mother bank)
  //   the number of levels in the geometrical tree and the
  //   maximum number of contents per level, by calling GGNLEV.
  //   Sets status bit for CONCAVE volumes, through GGCAVE.
  //   Completes the JSET structure with the list of volume names 
  //   which identify uniquely a given physical detector, the
  //   list of bit numbers to pack the corresponding volume copy 
  //   numbers, and the generic path(s) in the JVOLUM tree, 
  //   through the routine GHCLOS. 
  //
  ggclos(); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Glast() 
{ 
  //
  // Finish a Geant run
  //
  glast(); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gprint(const char *name) 
{ 
  //
  // Routine to print data structures
  // CHNAME   name of a data structure
  // 
  char vname[5];
  Vname(name,vname);
  gprint(PASSCHARD(vname),0 PASSCHARL(vname)); 
} 

//_____________________________________________________________________________
void  TGeant3::Grun() 
{ 
  //
  // Steering function to process one run
  //
  grun(); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gtrig() 
{ 
  //
  // Steering function to process one event
  //
  gtrig(); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gtrigc() 
{ 
  //
  // Clear event partition
  //
  gtrigc(); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gtrigi() 
{ 
  //
  // Initialises event partition
  //
  gtrigi(); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gwork(Int_t nwork) 
{ 
  //
  // Allocates workspace in ZEBRA memory
  //
  gwork(nwork); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gzinit() 
{ 
  //
  // To initialise GEANT/ZEBRA data structures
  //
  gzinit(); 
} 
 
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GCONS
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 
//_____________________________________________________________________________
void  TGeant3::Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z,  
		      Float_t &dens, Float_t &radl, Float_t &absl,
		      Float_t* ubuf, Int_t& nbuf) 
{ 
  //
  // Return parameters for material IMAT 
  //
  gfmate(imat, PASSCHARD(name), a, z, dens, radl, absl, ubuf, nbuf
	 PASSCHARL(name)); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gfpart(Int_t ipart, char *name, Int_t &itrtyp,  
		   Float_t &amass, Float_t &charge, Float_t &tlife) 
{ 
  //
  // Return parameters for particle of type IPART
  //
  Float_t *ubuf=0; 
  Int_t   nbuf; 
  Int_t igpart = IdFromPDG(ipart);
  gfpart(igpart, PASSCHARD(name), itrtyp, amass, charge, tlife, ubuf, nbuf
	 PASSCHARL(name)); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gftmed(Int_t numed, char *name, Int_t &nmat, Int_t &isvol,  
		   Int_t &ifield, Float_t &fieldm, Float_t &tmaxfd, 
		    Float_t &stemax, Float_t &deemax, Float_t &epsil, 
		    Float_t &stmin, Float_t *ubuf, Int_t *nbuf) 
{ 
  //
  // Return parameters for tracking medium NUMED
  //
  gftmed(numed, PASSCHARD(name), nmat, isvol, ifield, fieldm, tmaxfd, stemax,  
         deemax, epsil, stmin, ubuf, nbuf PASSCHARL(name)); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gmate() 
{ 
  //
  // Define standard GEANT materials
  //
  gmate(); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gpart() 
{ 
  //
  //  Define standard GEANT particles plus selected decay modes
  //  and branching ratios.
  //
  gpart(); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gsdk(Int_t ipart, Float_t *bratio, Int_t *mode) 
{ 
//  Defines branching ratios and decay modes for standard
//  GEANT particles.
   gsdk(ipart,bratio,mode); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gsmate(Int_t imat, const char *name, Float_t a, Float_t z,  
		   Float_t dens, Float_t radl, Float_t absl) 
{ 
  //
  // Defines a Material
  // 
  //  kmat               number assigned to the material
  //  name               material name
  //  a                  atomic mass in au
  //  z                  atomic number
  //  dens               density in g/cm3
  //  absl               absorbtion length in cm
  //                     if >=0 it is ignored and the program 
  //                     calculates it, if <0. -absl is taken
  //  radl               radiation length in cm
  //                     if >=0 it is ignored and the program 
  //                     calculates it, if <0. -radl is taken
  //  buf                pointer to an array of user words
  //  nbuf               number of user words
  //
  Float_t *ubuf=0; 
  Int_t   nbuf=0; 
  gsmate(imat,PASSCHARD(name), a, z, dens, radl, absl, ubuf, nbuf
	 PASSCHARL(name)); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gsmixt(Int_t imat, const char *name, Float_t *a, Float_t *z,  
		   Float_t dens, Int_t nlmat, Float_t *wmat) 
{ 
  //
  //       Defines mixture OR COMPOUND IMAT as composed by 
  //       THE BASIC NLMAT materials defined by arrays A,Z and WMAT
  // 
  //       If NLMAT.GT.0 then WMAT contains the PROPORTION BY
  //       WEIGTHS OF EACH BASIC MATERIAL IN THE MIXTURE. 
  // 
  //       If NLMAT.LT.0 then WMAT contains the number of atoms 
  //       of a given kind into the molecule of the COMPOUND
  //       In this case, WMAT in output is changed to relative
  //       weigths.
  //
  gsmixt(imat,PASSCHARD(name), a, z,dens, nlmat,wmat PASSCHARL(name)); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gspart(Int_t ipart, const char *name, Int_t itrtyp,  
		   Float_t amass, Float_t charge, Float_t tlife) 
{ 
  //
  // Store particle parameters
  //
  // ipart           particle code
  // name            particle name
  // itrtyp          transport method (see GEANT manual)
  // amass           mass in GeV/c2
  // charge          charge in electron units
  // tlife           lifetime in seconds
  //
  Float_t *ubuf=0; 
  Int_t   nbuf=0; 
  gspart(ipart,PASSCHARD(name), itrtyp, amass, charge, tlife, ubuf, nbuf
	 PASSCHARL(name)); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gstmed(Int_t numed, const char *name, Int_t nmat, Int_t isvol,  
		      Int_t ifield, Float_t fieldm, Float_t tmaxfd,
		      Float_t stemax, Float_t deemax, Float_t epsil,
		      Float_t stmin) 
{ 
  //
  //  NTMED  Tracking medium number
  //  NAME   Tracking medium name
  //  NMAT   Material number
  //  ISVOL  Sensitive volume flag
  //  IFIELD Magnetic field
  //  FIELDM Max. field value (Kilogauss)
  //  TMAXFD Max. angle due to field (deg/step)
  //  STEMAX Max. step allowed
  //  DEEMAX Max. fraction of energy lost in a step
  //  EPSIL  Tracking precision (cm)
  //  STMIN  Min. step due to continuos processes (cm)
  //
  //  IFIELD = 0 if no magnetic field; IFIELD = -1 if user decision in GUSWIM;
  //  IFIELD = 1 if tracking performed with GRKUTA; IFIELD = 2 if tracking
  //  performed with GHELIX; IFIELD = 3 if tracking performed with GHELX3.
  //  
  Float_t *ubuf=0; 
  Int_t   nbuf=0; 
  gstmed(numed,PASSCHARD(name), nmat, isvol, ifield, fieldm, tmaxfd, stemax,
	 deemax, epsil, stmin, ubuf, nbuf PASSCHARL(name)); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gsckov(Int_t itmed, Int_t npckov, Float_t *ppckov,
		      Float_t *absco, Float_t *effic, Float_t *rindex)
{ 
  //
  //    Stores the tables for UV photon tracking in medium ITMED 
  //    Please note that it is the user's responsability to 
  //    provide all the coefficients:
  //
  //
  //       ITMED       Tracking medium number
  //       NPCKOV      Number of bins of each table
  //       PPCKOV      Value of photon momentum (in GeV)
  //       ABSCO       Absorbtion coefficients 
  //                   dielectric: absorbtion length in cm
  //                   metals    : absorbtion fraction (0<=x<=1)
  //       EFFIC       Detection efficiency for UV photons 
  //       RINDEX      Refraction index (if=0 metal)
  //
  gsckov(itmed,npckov,ppckov,absco,effic,rindex);
}

//_____________________________________________________________________________
void  TGeant3::Gstpar(Int_t itmed, const char *param, Float_t parval) 
{ 
  //
  //  To change the value of cut  or mechanism "CHPAR"
  //      to a new value PARVAL  for tracking medium ITMED
  //    The  data   structure  JTMED   contains  the   standard  tracking
  //  parameters (CUTS and flags to control the physics processes)  which
  //  are used  by default  for all  tracking media.   It is  possible to
  //  redefine individually  with GSTPAR  any of  these parameters  for a
  //  given tracking medium. 
  //  ITMED     tracking medium number 
  //  CHPAR     is a character string (variable name) 
  //  PARVAL    must be given as a floating point.
  //
  gstpar(itmed,PASSCHARD(param), parval PASSCHARL(param)); 
} 
 
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GCONS
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 
//_____________________________________________________________________________
void  TGeant3::Gfkine(Int_t itra, Float_t *vert, Float_t *pvert, Int_t &ipart,
		      Int_t &nvert) 
{ 
  //           Storing/Retrieving Vertex and Track parameters
  //           ---------------------------------------------- 
  //
  //  Stores vertex parameters. 
  //  VERT      array of (x,y,z) position of the vertex 
  //  NTBEAM    beam track number origin of the vertex 
  //            =0 if none exists  
  //  NTTARG    target track number origin of the vertex
  //  UBUF      user array of NUBUF floating point numbers
  //  NUBUF       
  //  NVTX      new vertex number (=0 in case of error). 
  //  Prints vertex parameters.
  //  IVTX      for vertex IVTX.
  //            (For all vertices if IVTX=0) 
  //  Stores long life track parameters.
  //  PLAB      components of momentum 
  //  IPART     type of particle (see GSPART)
  //  NV        vertex number origin of track
  //  UBUF      array of NUBUF floating point user parameters 
  //  NUBUF
  //  NT        track number (if=0 error).
  //  Retrieves long life track parameters.
  //  ITRA      track number for which parameters are requested
  //  VERT      vector origin of the track  
  //  PVERT     4 momentum components at the track origin 
  //  IPART     particle type (=0 if track ITRA does not exist)
  //  NVERT     vertex number origin of the track 
  //  UBUF      user words stored in GSKINE. 
  //  Prints initial track parameters. 
  //  ITRA      for track ITRA 
  //            (For all tracks if ITRA=0) 
  //
  Float_t *ubuf=0; 
  Int_t   nbuf; 
  gfkine(itra,vert,pvert,ipart,nvert,ubuf,nbuf); 
} 

//_____________________________________________________________________________
void  TGeant3::Gfvert(Int_t nvtx, Float_t *v, Int_t &ntbeam, Int_t &nttarg,
		      Float_t &tofg) 
{ 
  //
  //       Retrieves the parameter of a vertex bank
  //       Vertex is generated from tracks NTBEAM NTTARG
  //       NVTX is the new vertex number 
  //
  Float_t *ubuf=0; 
  Int_t   nbuf; 
  gfvert(nvtx,v,ntbeam,nttarg,tofg,ubuf,nbuf); 
} 
 
//_____________________________________________________________________________
Int_t TGeant3::Gskine(Float_t *plab, Int_t ipart, Int_t nv, Float_t *buf,
		      Int_t nwbuf) 
{ 
  //
  //       Store kinematics of track NT into data structure
  //       Track is coming from vertex NV
  //
  Int_t nt = 0; 
  gskine(plab, ipart, nv, buf, nwbuf, nt); 
  return nt; 
} 
 
//_____________________________________________________________________________
Int_t TGeant3::Gsvert(Float_t *v, Int_t ntbeam, Int_t nttarg, Float_t *ubuf,
		      Int_t nwbuf) 
{ 
  //
  //       Creates a new vertex bank 
  //       Vertex is generated from tracks NTBEAM NTTARG 
  //       NVTX is the new vertex number
  //
  Int_t nwtx = 0; 
  gsvert(v, ntbeam, nttarg, ubuf, nwbuf, nwtx); 
  return nwtx; 
} 
 
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GPHYS
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//_____________________________________________________________________________
void  TGeant3::Gphysi() 
{ 
  //
  //       Initialise material constants for all the physics
  //       mechanisms used by GEANT
  //
  gphysi(); 
} 
 
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GTRAK
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 
//_____________________________________________________________________________
void  TGeant3::Gdebug() 
{ 
  //
  // Debug the current step
  //
  gdebug(); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gekbin() 
{ 
  //
  //       To find bin number in kinetic energy table
  //       stored in ELOW(NEKBIN)
  //
  gekbin(); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gfinds() 
{ 
  //
  //       Returns the set/volume parameters corresponding to 
  //       the current space point in /GCTRAK/
  //       and fill common /GCSETS/
  // 
  //       IHSET  user set identifier 
  //       IHDET  user detector identifier 
  //       ISET set number in JSET  
  //       IDET   detector number in JS=LQ(JSET-ISET) 
  //       IDTYPE detector type (1,2)  
  //       NUMBV  detector volume numbers (array of length NVNAME)
  //       NVNAME number of volume levels
  //
  gfinds(); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gsking(Int_t igk) 
{ 
  //
  //   Stores in stack JSTAK either the IGKth track of /GCKING/,
  //    or the NGKINE tracks when IGK is 0.
  //
  gsking(igk); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gskpho(Int_t igk) 
{ 
  //
  //  Stores in stack JSTAK either the IGKth Cherenkov photon of  
  //  /GCKIN2/, or the NPHOT tracks when IGK is 0.                
  //
  gskpho(igk); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gsstak(Int_t iflag) 
{ 
  //
  //   Stores in auxiliary stack JSTAK the particle currently 
  //    described in common /GCKINE/. 
  // 
  //   On request, creates also an entry in structure JKINE :
  //    IFLAG =
  //     0 : No entry in JKINE structure required (user) 
  //     1 : New entry in JVERTX / JKINE structures required (user)
  //    <0 : New entry in JKINE structure at vertex -IFLAG (user)
  //     2 : Entry in JKINE structure exists already (from GTREVE)
  //
  gsstak(iflag); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gsxyz() 
{ 
  //
  //   Store space point VECT in banks JXYZ 
  //
  gsxyz(); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gtrack() 
{ 
  //
  //   Controls tracking of current particle 
  //
  gtrack(); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gtreve() 
{ 
  //
  //   Controls tracking of all particles belonging to the current event
  //
  gtreve(); 
} 

//_____________________________________________________________________________
void  TGeant3::Gtreve_root() 
{ 
  //
  //   Controls tracking of all particles belonging to the current event
  //
  gtreve_root(); 
} 

//_____________________________________________________________________________
void  TGeant3::Grndm(Float_t *rvec, const Int_t len) const
{
  //
  //   To generate a vector RVECV of LEN random numbers 
  //   Copy of the CERN Library routine RANECU 
  grndm(rvec,len);
}

//_____________________________________________________________________________
void  TGeant3::Grndmq(Int_t &is1, Int_t &is2, const Int_t iseq,
		      const Text_t *chopt)
{
  //
  //  To set/retrieve the seed of the random number generator
  //
  grndmq(is1,is2,iseq,PASSCHARD(chopt) PASSCHARL(chopt));
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GDRAW
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//_____________________________________________________________________________
void  TGeant3::Gdxyz(Int_t it)
{
  //
  // Draw the points stored with Gsxyz relative to track it
  //
  gdxyz(it);
}

//_____________________________________________________________________________
void  TGeant3::Gdcxyz()
{
  //
  // Draw the position of the current track
  //
  gdcxyz();
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GGEOM
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//_____________________________________________________________________________
void  TGeant3::Gdtom(Float_t *xd, Float_t *xm, Int_t iflag) 
{ 
  //
  //  Computes coordinates XM (Master Reference System
  //  knowing the coordinates XD (Detector Ref System)
  //  The local reference system can be initialized by
  //    - the tracking routines and GDTOM used in GUSTEP
  //    - a call to GSCMED(NLEVEL,NAMES,NUMBER)
  //        (inverse routine is GMTOD)
  // 
  //   If IFLAG=1  convert coordinates
  //      IFLAG=2  convert direction cosinus
  //
  gdtom(xd, xm, iflag); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Glmoth(const char* iudet, Int_t iunum, Int_t &nlev, Int_t *lvols,
		      Int_t *lindx) 
{ 
  //
  //   Loads the top part of the Volume tree in LVOLS (IVO's),
  //   LINDX (IN indices) for a given volume defined through
  //   its name IUDET and number IUNUM.
  // 
  //   The routine stores only upto the last level where JVOLUM
  //   data structure is developed. If there is no development
  //   above the current level, it returns NLEV zero.
  Int_t *idum=0; 
  glmoth(PASSCHARD(iudet), iunum, nlev, lvols, lindx, idum PASSCHARL(iudet)); 
} 

//_____________________________________________________________________________
void  TGeant3::Gmedia(Float_t *x, Int_t &numed) 
{ 
  //
  //   Finds in which volume/medium the point X is, and updates the
  //    common /GCVOLU/ and the structure JGPAR accordingly. 
  // 
  //   NUMED returns the tracking medium number, or 0 if point is
  //         outside the experimental setup.
  //
  gmedia(x,numed); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gmtod(Float_t *xm, Float_t *xd, Int_t iflag) 
{ 
  //
  //       Computes coordinates XD (in DRS) 
  //       from known coordinates XM in MRS 
  //       The local reference system can be initialized by
  //         - the tracking routines and GMTOD used in GUSTEP
  //         - a call to GMEDIA(XM,NUMED)
  //         - a call to GLVOLU(NLEVEL,NAMES,NUMBER,IER) 
  //             (inverse routine is GDTOM) 
  //
  //        If IFLAG=1  convert coordinates 
  //           IFLAG=2  convert direction cosinus
  //
  gmtod(xm, xd, iflag); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gsdvn(const char *name, const char *mother, Int_t ndiv,
		     Int_t iaxis) 
{ 
  //
  // Create a new volume by dividing an existing one
  // 
  //  NAME   Volume name
  //  MOTHER Mother volume name
  //  NDIV   Number of divisions
  //  IAXIS  Axis value
  //
  //  X,Y,Z of CAXIS will be translated to 1,2,3 for IAXIS.
  //  It divides a previously defined volume.
  //  
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);
  gsdvn(PASSCHARD(vname), PASSCHARD(vmother), ndiv, iaxis PASSCHARL(vname)
	PASSCHARL(vmother)); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gsdvn2(const char *name, const char *mother, Int_t ndiv,
		      Int_t iaxis, Float_t c0i, Int_t numed) 
{ 
  //
  // Create a new volume by dividing an existing one
  // 
  // Divides mother into ndiv divisions called name
  // along axis iaxis starting at coordinate value c0.
  // the new volume created will be medium number numed.
  //
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);
  gsdvn2(PASSCHARD(vname), PASSCHARD(vmother), ndiv, iaxis, c0i, numed
	 PASSCHARL(vname) PASSCHARL(vmother)); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gsdvs(const char *name, const char *mother, Float_t step,
		     Int_t iaxis, Int_t numed) 
{ 
  //
  // Create a new volume by dividing an existing one
  // 
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);
  gsdvs(PASSCHARD(vname), PASSCHARD(vmother), step, iaxis, numed
	PASSCHARL(vname) PASSCHARL(vmother)); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gsdvs2(const char *name, const char *mother, Float_t step,
		      Int_t iaxis, Float_t c0, Int_t numed) 
{ 
  //
  // Create a new volume by dividing an existing one
  // 
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);
  gsdvs2(PASSCHARD(vname), PASSCHARD(vmother), step, iaxis, c0, numed
	 PASSCHARL(vname) PASSCHARL(vmother)); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gsdvt(const char *name, const char *mother, Float_t step,
		     Int_t iaxis, Int_t numed, Int_t ndvmx) 
{ 
  //
  // Create a new volume by dividing an existing one
  // 
  //       Divides MOTHER into divisions called NAME along
  //       axis IAXIS in steps of STEP. If not exactly divisible 
  //       will make as many as possible and will centre them 
  //       with respect to the mother. Divisions will have medium 
  //       number NUMED. If NUMED is 0, NUMED of MOTHER is taken.
  //       NDVMX is the expected maximum number of divisions
  //          (If 0, no protection tests are performed) 
  //
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);
  gsdvt(PASSCHARD(vname), PASSCHARD(vmother), step, iaxis, numed, ndvmx
	PASSCHARL(vname) PASSCHARL(vmother)); 
} 

//_____________________________________________________________________________
void  TGeant3::Gsdvt2(const char *name, const char *mother, Float_t step,
		      Int_t iaxis, Float_t c0, Int_t numed, Int_t ndvmx) 
{ 
  //
  // Create a new volume by dividing an existing one
  //                                                                    
  //           Divides MOTHER into divisions called NAME along          
  //            axis IAXIS starting at coordinate value C0 with step    
  //            size STEP.                                              
  //           The new volume created will have medium number NUMED.    
  //           If NUMED is 0, NUMED of mother is taken.                 
  //           NDVMX is the expected maximum number of divisions        
  //             (If 0, no protection tests are performed)              
  //
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);
  gsdvt2(PASSCHARD(vname), PASSCHARD(vmother), step, iaxis, c0,
	 numed, ndvmx PASSCHARL(vname) PASSCHARL(vmother)); 
} 

//_____________________________________________________________________________
void  TGeant3::Gsord(const char *name, Int_t iax) 
{ 
  //
  //    Flags volume CHNAME whose contents will have to be ordered 
  //    along axis IAX, by setting the search flag to -IAX
  //           IAX = 1    X axis 
  //           IAX = 2    Y axis 
  //           IAX = 3    Z axis 
  //           IAX = 4    Rxy (static ordering only  -> GTMEDI)
  //           IAX = 14   Rxy (also dynamic ordering -> GTNEXT)
  //           IAX = 5    Rxyz (static ordering only -> GTMEDI)
  //           IAX = 15   Rxyz (also dynamic ordering -> GTNEXT)
  //           IAX = 6    PHI   (PHI=0 => X axis)
  //           IAX = 7    THETA (THETA=0 => Z axis)
  //
  char vname[5];
  Vname(name,vname);
  gsord(PASSCHARD(vname), iax PASSCHARL(vname)); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gspos(const char *name, Int_t nr, const char *mother, Float_t x,
		     Float_t y, Float_t z, Int_t irot, const char *konly) 
{ 
  //
  // Position a volume into an existing one
  //
  //  NAME   Volume name
  //  NUMBER Copy number of the volume
  //  MOTHER Mother volume name
  //  X      X coord. of the volume in mother ref. sys.
  //  Y      Y coord. of the volume in mother ref. sys.
  //  Z      Z coord. of the volume in mother ref. sys.
  //  IROT   Rotation matrix number w.r.t. mother ref. sys.
  //  ONLY   ONLY/MANY flag
  //
  //  It positions a previously defined volume in the mother.
  //  
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);
  gspos(PASSCHARD(vname), nr, PASSCHARD(vmother), x, y, z, irot,
	PASSCHARD(konly) PASSCHARL(vname) PASSCHARL(vmother)
	PASSCHARL(konly)); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gsposp(const char *name, Int_t nr, const char *mother,  
		   Float_t x, Float_t y, Float_t z, Int_t irot,
		      const char *konly, Float_t *upar, Int_t np ) 
{ 
  //
  //      Place a copy of generic volume NAME with user number
  //      NR inside MOTHER, with its parameters UPAR(1..NP)
  //
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);
  gsposp(PASSCHARD(vname), nr, PASSCHARD(vmother), x, y, z, irot,
	 PASSCHARD(konly), upar, np PASSCHARL(vname) PASSCHARL(vmother)
	 PASSCHARL(konly)); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gsrotm(Int_t nmat, Float_t theta1, Float_t phi1, Float_t theta2,
		      Float_t phi2, Float_t theta3, Float_t phi3) 
{ 
  //
  //  nmat   Rotation matrix number
  //  THETA1 Polar angle for axis I
  //  PHI1   Azimuthal angle for axis I
  //  THETA2 Polar angle for axis II
  //  PHI2   Azimuthal angle for axis II
  //  THETA3 Polar angle for axis III
  //  PHI3   Azimuthal angle for axis III
  //
  //  It defines the rotation matrix number IROT.
  //  
  gsrotm(nmat, theta1, phi1, theta2, phi2, theta3, phi3); 
} 
 
//_____________________________________________________________________________
void  TGeant3::Gprotm(Int_t nmat) 
{ 
  //
  //    To print rotation matrices structure JROTM
  //     nmat     Rotation matrix number
  //
  gprotm(nmat); 
} 
 
//_____________________________________________________________________________
Int_t TGeant3::Gsvolu(const char *name, const char *shape, Int_t nmed,  
		      Float_t *upar, Int_t npar) 
{ 
  //
  //  NAME   Volume name
  //  SHAPE  Volume type
  //  NUMED  Tracking medium number
  //  NPAR   Number of shape parameters
  //  UPAR   Vector containing shape parameters
  //
  //  It creates a new volume in the JVOLUM data structure.
  //  
  Int_t ivolu = 0; 
  char vname[5];
  Vname(name,vname);
  char vshape[5];
  Vname(shape,vshape);
  gsvolu(PASSCHARD(vname), PASSCHARD(vshape), nmed, upar, npar, ivolu
	 PASSCHARL(vname) PASSCHARL(vshape)); 
  return ivolu; 
} 
 
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//           T H E    D R A W I N G   P A C K A G E
//           ======================================
//  Drawing functions. These functions allow the visualization in several ways
//  of the volumes defined in the geometrical data structure. It is possible
//  to draw the logical tree of volumes belonging to the detector (DTREE),
//  to show their geometrical specification (DSPEC,DFSPC), to draw them
//  and their cut views (DRAW, DCUT). Moreover, it is possible to execute
//  these commands when the hidden line removal option is activated; in
//  this case, the volumes can be also either translated in the space
//  (SHIFT), or clipped by boolean operation (CVOL). In addition, it is
//  possible to fill the surfaces of the volumes
//  with solid colours when the shading option (SHAD) is activated.
//  Several tools (ZOOM, LENS) have been developed to zoom detailed parts
//  of the detectors or to scan physical events as well.
//  Finally, the command MOVE will allow the rotation, translation and zooming
//  on real time parts of the detectors or tracks and hits of a simulated event.
//  Ray-tracing commands. In case the command (DOPT RAYT ON) is executed,
//  the drawing is performed by the Geant ray-tracing;
//  automatically, the color is assigned according to the tracking medium of each
//  volume and the volumes with a density lower/equal than the air are considered
//  transparent; if the option (USER) is set (ON) (again via the command (DOPT)),
//  the user can set color and visibility for the desired volumes via the command
//  (SATT), as usual, relatively to the attributes (COLO) and (SEEN).
//  The resolution can be set via the command (SATT * FILL VALUE), where (VALUE)
//  is the ratio between the number of pixels drawn and 20 (user coordinates).
//  Parallel view and perspective view are possible (DOPT PROJ PARA/PERS); in the
//  first case, we assume that the first mother volume of the tree is a box with
//  dimensions 10000 X 10000 X 10000 cm and the view point (infinetely far) is
//  5000 cm far from the origin along the Z axis of the user coordinates; in the
//  second case, the distance between the observer and the origin of the world
//  reference system is set in cm by the command (PERSP NAME VALUE); grand-angle
//  or telescopic effects can be achieved changing the scale factors in the command
//  (DRAW). When the final picture does not occupy the full window,
//  mapping the space before tracing can speed up the drawing, but can also
//  produce less precise results; values from 1 to 4 are allowed in the command
//  (DOPT MAPP VALUE), the mapping being more precise for increasing (VALUE); for
//  (VALUE = 0) no mapping is performed (therefore max precision and lowest speed).
//  The command (VALCUT) allows the cutting of the detector by three planes
//  ortogonal to the x,y,z axis. The attribute (LSTY) can be set by the command
//  SATT for any desired volume and can assume values from 0 to 7; it determines
//  the different light processing to be performed for different materials:
//  0 = dark-matt, 1 = bright-matt, 2 = plastic, 3 = ceramic, 4 = rough-metals,
//  5 = shiny-metals, 6 = glass, 7 = mirror. The detector is assumed to be in the
//  dark, the ambient light luminosity is 0.2 for each basic hue (the saturation
//  is 0.9) and the observer is assumed to have a light source (therefore he will
//  produce parallel light in the case of parallel view and point-like-source
//  light in the case of perspective view).
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//_____________________________________________________________________________
void TGeant3::Gsatt(const char *name, const char *att, Int_t val)
{ 
  //
  //  NAME   Volume name
  //  IOPT   Name of the attribute to be set
  //  IVAL   Value to which the attribute is to be set
  //
  //  name= "*" stands for all the volumes.
  //  iopt can be chosen among the following :
  //  
  //     WORK   0=volume name is inactive for the tracking
  //            1=volume name is active for the tracking (default)
  //  
  //     SEEN   0=volume name is invisible
  //            1=volume name is visible (default)
  //           -1=volume invisible with all its descendants in the tree
  //           -2=volume visible but not its descendants in the tree
  //  
  //     LSTY   line style 1,2,3,... (default=1)
  //            LSTY=7 will produce a very precise approximation for
  //            revolution bodies.
  //  
  //     LWID   line width -7,...,1,2,3,..7 (default=1)
  //            LWID<0 will act as abs(LWID) was set for the volume
  //            and for all the levels below it. When SHAD is 'ON', LWID
  //            represent the linewidth of the scan lines filling the surfaces
  //            (whereas the FILL value represent their number). Therefore
  //            tuning this parameter will help to obtain the desired
  //            quality/performance ratio.
  //  
  //     COLO   colour code -166,...,1,2,..166 (default=1)
  //            n=1=black
  //            n=2=red;    n=17+m, m=0,25, increasing luminosity according to 'm';
  //            n=3=green;  n=67+m, m=0,25, increasing luminosity according to 'm';
  //            n=4=blue;   n=117+m, m=0,25, increasing luminosity according to 'm';
  //            n=5=yellow; n=42+m, m=0,25, increasing luminosity according to 'm';
  //            n=6=violet; n=142+m, m=0,25, increasing luminosity according to 'm';
  //            n=7=lightblue; n=92+m, m=0,25, increasing luminosity according to 'm';
  //            colour=n*10+m, m=1,2,...9, will produce the same colour
  //            as 'n', but with increasing luminosity according to 'm';
  //            COLO<0 will act as if abs(COLO) was set for the volume
  //            and for all the levels below it.
  //            When for a volume the attribute FILL is > 1 (and the
  //            option SHAD is on), the ABS of its colour code must be < 8
  //            because an automatic shading of its faces will be
  //            performed.
  //  
  //     FILL  (1992) fill area  -7,...,0,1,...7 (default=0)
  //            when option SHAD is "on" the FILL attribute of any
  //            volume can be set different from 0 (normal drawing);
  //            if it is set to 1, the faces of such volume will be filled
  //            with solid colours; if ABS(FILL) is > 1, then a light
  //            source is placed along the observer line, and the faces of
  //            such volumes will be painted by colours whose luminosity
  //            will depend on the amount of light reflected;
  //            if ABS(FILL) = 1, then it is possible to use all the 166
  //            colours of the colour table, becouse the automatic shading
  //            is not performed;
  //            for increasing values of FILL the drawing will be performed
  //            with higher and higher resolution improving the quality (the
  //            number of scan lines used to fill the faces increases with FILL);
  //            it is possible to set different values of FILL
  //            for different volumes, in order to optimize at the same time
  //            the performance and the quality of the picture;
  //            FILL<0 will act as if abs(FILL) was set for the volume
  //            and for all the levels below it.
  //            This kind of drawing can be saved in 'picture files'
  //            or in view banks.
  //            0=drawing without fill area
  //            1=faces filled with solid colours and resolution = 6
  //            2=lowest resolution (very fast)
  //            3=default resolution
  //            4=.................
  //            5=.................
  //            6=.................
  //            7=max resolution
  //            Finally, if a coloured background is desired, the FILL
  //            attribute for the first volume of the tree must be set
  //            equal to -abs(colo), colo being >0 and <166.
  //  
  //     SET   set number associated to volume name
  //     DET   detector number associated to volume name
  //     DTYP  detector type (1,2)
  //  
  InitHIGZ();
  char vname[5];
  Vname(name,vname);
  char vatt[5];
  Vname(att,vatt);
  gsatt(PASSCHARD(vname), PASSCHARD(vatt), val PASSCHARL(vname)
	PASSCHARL(vatt)); 
} 

//_____________________________________________________________________________
void TGeant3::Gfpara(const char *name, Int_t number, Int_t intext, Int_t& npar,
			 Int_t& natt, Float_t* par, Float_t* att)
{
  //
  // Find the parameters of a volume
  //
  gfpara(PASSCHARD(name), number, intext, npar, natt, par, att
	 PASSCHARL(name));
}

//_____________________________________________________________________________
void TGeant3::Gckpar(Int_t ish, Int_t npar, Float_t* par)
{
  //
  // Check the parameters of a shape
  //
  gckpar(ish,npar,par);
}

//_____________________________________________________________________________
void TGeant3::Gckmat(Int_t itmed, char* natmed)
{
  //
  // Check the parameters of a tracking medium
  //
  gckmat(itmed, PASSCHARD(natmed) PASSCHARL(natmed));
}

//_____________________________________________________________________________
void TGeant3::Gdelete(Int_t iview)
{ 
  //
  //  IVIEW  View number
  //
  //  It deletes a view bank from memory.
  //
  gdelet(iview);
}
 
//_____________________________________________________________________________
void TGeant3::Gdopen(Int_t iview)
{ 
  //
  //  IVIEW  View number
  //
  //  When a drawing is very complex and requires a long time to be
  //  executed, it can be useful to store it in a view bank: after a
  //  call to DOPEN and the execution of the drawing (nothing will
  //  appear on the screen), and after a necessary call to DCLOSE,
  //  the contents of the bank can be displayed in a very fast way
  //  through a call to DSHOW; therefore, the detector can be easily
  //  zoomed many times in different ways. Please note that the pictures
  //  with solid colours can now be stored in a view bank or in 'PICTURE FILES'
  //
  InitHIGZ();
  higz->Clear();
  gdopen(iview);
}
 
//_____________________________________________________________________________
void TGeant3::Gdclose()
{ 
  //
  //  It closes the currently open view bank; it must be called after the
  //  end of the drawing to be stored.
  //
  gdclos();
}
 
//_____________________________________________________________________________
void TGeant3::Gdshow(Int_t iview)
{ 
  //
  //  IVIEW  View number
  //
  //  It shows on the screen the contents of a view bank. It
  //  can be called after a view bank has been closed.
  //
  gdshow(iview);
} 

//_____________________________________________________________________________
void TGeant3::Gdopt(const char *name,const char *value)
{ 
  //
  //  NAME   Option name
  //  VALUE  Option value
  //
  //  To set/modify the drawing options.
  //     IOPT   IVAL      Action
  //  
  //     THRZ    ON       Draw tracks in R vs Z
  //             OFF (D)  Draw tracks in X,Y,Z
  //             180
  //             360
  //     PROJ    PARA (D) Parallel projection
  //             PERS     Perspective
  //     TRAK    LINE (D) Trajectory drawn with lines
  //             POIN       " " with markers
  //     HIDE    ON       Hidden line removal using the CG package
  //             OFF (D)  No hidden line removal
  //     SHAD    ON       Fill area and shading of surfaces.
  //             OFF (D)  Normal hidden line removal.
  //     RAYT    ON       Ray-tracing on.
  //             OFF (D)  Ray-tracing off.
  //     EDGE    OFF      Does not draw contours when shad is on.
  //             ON  (D)  Normal shading.
  //     MAPP    1,2,3,4  Mapping before ray-tracing.
  //             0   (D)  No mapping.
  //     USER    ON       User graphics options in the raytracing.
  //             OFF (D)  Automatic graphics options.
  //  
  InitHIGZ();
  char vname[5];
  Vname(name,vname);
  char vvalue[5];
  Vname(value,vvalue);
  gdopt(PASSCHARD(vname), PASSCHARD(vvalue) PASSCHARL(vname)
	PASSCHARL(vvalue)); 
} 
 
//_____________________________________________________________________________
void TGeant3::Gdraw(const char *name,Float_t theta, Float_t phi, Float_t psi,
		    Float_t u0,Float_t v0,Float_t ul,Float_t vl)
{ 
  //
  //  NAME   Volume name
  //  +
  //  THETA  Viewing angle theta (for 3D projection)
  //  PHI    Viewing angle phi (for 3D projection)
  //  PSI    Viewing angle psi (for 2D rotation)
  //  U0     U-coord. (horizontal) of volume origin
  //  V0     V-coord. (vertical) of volume origin
  //  SU     Scale factor for U-coord.
  //  SV     Scale factor for V-coord.
  //
  //  This function will draw the volumes,
  //  selected with their graphical attributes, set by the Gsatt
  //  facility. The drawing may be performed with hidden line removal
  //  and with shading effects according to the value of the options HIDE
  //  and SHAD; if the option SHAD is ON, the contour's edges can be
  //  drawn or not. If the option HIDE is ON, the detector can be
  //  exploded (BOMB), clipped with different shapes (CVOL), and some
  //  of its parts can be shifted from their original
  //  position (SHIFT). When HIDE is ON, if
  //  the drawing requires more than the available memory, the program
  //  will evaluate and display the number of missing words
  //  (so that the user can increase the
  //  size of its ZEBRA store). Finally, at the end of each drawing (with HIDE on),
  //  the program will print messages about the memory used and
  //  statistics on the volumes' visibility.
  //  The following commands will produce the drawing of a green
  //  volume, specified by NAME, without using the hidden line removal
  //  technique, using the hidden line removal technique,
  //  with different linewidth and colour (red), with
  //  solid colour, with shading of surfaces, and without edges.
  //  Finally, some examples are given for the ray-tracing. (A possible
  //  string for the NAME of the volume can be found using the command DTREE).
  //
  InitHIGZ();
  higz->Clear();
  char vname[5];
  Vname(name,vname);
  if (fGcvdma->raytra != 1) {
    gdraw(PASSCHARD(vname), theta,phi,psi,u0,v0,ul,vl PASSCHARL(vname)); 
  } else {
    gdrayt(PASSCHARD(vname), theta,phi,psi,u0,v0,ul,vl PASSCHARL(vname)); 
  }
} 
 
//_____________________________________________________________________________
void TGeant3::Gdrawc(const char *name,Int_t axis, Float_t cut,Float_t u0,
		     Float_t v0,Float_t ul,Float_t vl)
{ 
  //
  //  NAME   Volume name
  //  CAXIS  Axis value
  //  CUTVAL Cut plane distance from the origin along the axis
  //  +
  //  U0     U-coord. (horizontal) of volume origin
  //  V0     V-coord. (vertical) of volume origin
  //  SU     Scale factor for U-coord.
  //  SV     Scale factor for V-coord.
  //
  //  The cut plane is normal to caxis (X,Y,Z), corresponding to iaxis (1,2,3),
  //  and placed at the distance cutval from the origin.
  //  The resulting picture is seen from the the same axis.
  //  When HIDE Mode is ON, it is possible to get the same effect with
  //  the CVOL/BOX function.
  //  
  InitHIGZ();
  higz->Clear();
  char vname[5];
  Vname(name,vname);
  gdrawc(PASSCHARD(vname), axis,cut,u0,v0,ul,vl PASSCHARL(vname)); 
} 
 
//_____________________________________________________________________________
void TGeant3::Gdrawx(const char *name,Float_t cutthe, Float_t cutphi,
		     Float_t cutval, Float_t theta, Float_t phi, Float_t u0,
		     Float_t v0,Float_t ul,Float_t vl)
{ 
  //
  //  NAME   Volume name
  //  CUTTHE Theta angle of the line normal to cut plane
  //  CUTPHI Phi angle of the line normal to cut plane
  //  CUTVAL Cut plane distance from the origin along the axis
  //  +
  //  THETA  Viewing angle theta (for 3D projection)
  //  PHI    Viewing angle phi (for 3D projection)
  //  U0     U-coord. (horizontal) of volume origin
  //  V0     V-coord. (vertical) of volume origin
  //  SU     Scale factor for U-coord.
  //  SV     Scale factor for V-coord.
  //
  //  The cut plane is normal to the line given by the cut angles
  //  cutthe and cutphi and placed at the distance cutval from the origin.
  //  The resulting picture is seen from the viewing angles theta,phi.
  //
  InitHIGZ();
  higz->Clear();
  char vname[5];
  Vname(name,vname);
  gdrawx(PASSCHARD(vname), cutthe,cutphi,cutval,theta,phi,u0,v0,ul,vl
	 PASSCHARL(vname)); 
}
 
//_____________________________________________________________________________
void TGeant3::Gdhead(Int_t isel, const char *name, Float_t chrsiz)
{ 
  //
  //  Parameters
  //  +
  //  ISEL   Option flag  D=111110
  //  NAME   Title
  //  CHRSIZ Character size (cm) of title NAME D=0.6
  //
  //  ISEL =
  //   0      to have only the header lines
  //   xxxxx1 to add the text name centered on top of header
  //   xxxx1x to add global detector name (first volume) on left
  //   xxx1xx to add date on right
  //   xx1xxx to select thick characters for text on top of header
  //   x1xxxx to add the text 'EVENT NR x' on top of header
  //   1xxxxx to add the text 'RUN NR x' on top of header
  //  NOTE that ISEL=x1xxx1 or ISEL=1xxxx1 are illegal choices,
  //  i.e. they generate overwritten text.
  //
  gdhead(isel,PASSCHARD(name),chrsiz PASSCHARL(name));
}

//_____________________________________________________________________________
void TGeant3::Gdman(Float_t u, Float_t v, const char *type)
{ 
  //
  //  Draw a 2D-man at position (U0,V0)
  //  Parameters
  //  U      U-coord. (horizontal) of the centre of man' R
  //  V      V-coord. (vertical) of the centre of man' R
  //  TYPE   D='MAN' possible values: 'MAN,WM1,WM2,WM3'
  // 
  //   CALL GDMAN(u,v),CALL GDWMN1(u,v),CALL GDWMN2(u,v),CALL GDWMN2(u,v)
  //  It superimposes the picure of a man or of a woman, chosen among
  //  three different ones, with the same scale factors as the detector
  //  in the current drawing.
  //
  TString opt = type;
   if (opt.Contains("WM1")) {
     gdwmn1(u,v);
   } else if (opt.Contains("WM3")) {
     gdwmn3(u,v);
   } else if (opt.Contains("WM2")) {
     gdwmn2(u,v);
   } else {
     gdman(u,v);
   }
}
 
//_____________________________________________________________________________
void TGeant3::Gdspec(const char *name)
{ 
  //
  //  NAME   Volume name
  //
  //  Shows 3 views of the volume (two cut-views and a 3D view), together with
  //  its geometrical specifications. The 3D drawing will
  //  be performed according the current values of the options HIDE and
  //  SHAD and according the current SetClipBox clipping parameters for that
  //  volume.
  //  
  InitHIGZ();
  higz->Clear();
  char vname[5];
  Vname(name,vname);
  gdspec(PASSCHARD(vname) PASSCHARL(vname)); 
} 
 
//_____________________________________________________________________________
void TGeant3::DrawOneSpec(const char *name)
{ 
  //
  //  Function called when one double-clicks on a volume name
  //  in a TPavelabel drawn by Gdtree.
  //
  THIGZ *higzSave = higz;
  higzSave->SetName("higzSave");
  THIGZ *higzSpec = (THIGZ*)gROOT->FindObject("higzSpec");
  //printf("DrawOneSpec, higz=%x, higzSpec=%x\n",higz,higzSpec);
  if (higzSpec) higz     = higzSpec;
  else          higzSpec = new THIGZ(defSize);
  higzSpec->SetName("higzSpec");
  higzSpec->cd();
  higzSpec->Clear();
  char vname[5];
  Vname(name,vname);
  gdspec(PASSCHARD(vname) PASSCHARL(vname)); 
  higzSpec->Update();
  higzSave->cd();
  higzSave->SetName("higz");
  higz = higzSave;
} 

//_____________________________________________________________________________
void TGeant3::Gdtree(const char *name,Int_t levmax, Int_t isel)
{ 
  //
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
  InitHIGZ();
  higz->Clear();
  char vname[5];
  Vname(name,vname);
  gdtree(PASSCHARD(vname), levmax, isel PASSCHARL(vname)); 
  higz->fPname = "";
} 

//_____________________________________________________________________________
void TGeant3::GdtreeParent(const char *name,Int_t levmax, Int_t isel)
{ 
  //
  //  NAME   Volume name
  //  LEVMAX Depth level
  //  ISELT  Options
  //
  //  This function draws the logical tree of the parent of name.
  //  
  InitHIGZ();
  higz->Clear();
  // Scan list of volumes in JVOLUM
  char vname[5];
  Int_t gname, i, jvo, in, nin, jin, num;
  strncpy((char *) &gname, name, 4);
  for(i=1; i<=fGcnum->nvolum; i++) {
    jvo = fZlq[fGclink->jvolum-i];
    nin = Int_t(fZq[jvo+3]);
    if (nin == -1) nin = 1;
    for (in=1;in<=nin;in++) {
      jin = fZlq[jvo-in];
      num = Int_t(fZq[jin+2]);
      if(gname == fZiq[fGclink->jvolum+num]) {
	strncpy(vname,(char*)&fZiq[fGclink->jvolum+i],4);
	vname[4] = 0;           
	gdtree(PASSCHARD(vname), levmax, isel PASSCHARL(vname)); 
	higz->fPname = "";
	return;
      }
    }
  }
} 
 
//_____________________________________________________________________________
void TGeant3::SetABAN(Int_t par)
{
  //
  // par = 1 particles will be stopped according to their residual
  //         range if they are not in a sensitive material and are
  //         far enough from the boundary
  //       0 particles are transported normally
  //
  fGcphys->dphys1 = par;
}
 
 
//_____________________________________________________________________________
void TGeant3::SetANNI(Int_t par)
{
  //
  //   To control positron annihilation.
  //    par =0 no annihilation
  //        =1 annihilation. Decays processed.
  //        =2 annihilation. No decay products stored.
  //
  fGcphys->ianni = par;
}
 
 
//_____________________________________________________________________________
void TGeant3::SetAUTO(Int_t par)
{
  //
  //  To control automatic calculation of tracking medium parameters:
  //   par =0 no automatic calculation;
  //       =1 automati calculation.
  //  
  fGctrak->igauto = par;
}
 
 
//_____________________________________________________________________________
void TGeant3::SetBOMB(Float_t boom)
{
  //
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
  //
  InitHIGZ();
  setbomb(boom);
}
 
//_____________________________________________________________________________
void TGeant3::SetBREM(Int_t par)
{
  //
  //  To control bremstrahlung.
  //   par =0 no bremstrahlung
  //       =1 bremstrahlung. Photon processed.
  //       =2 bremstrahlung. No photon stored.
  //  
  fGcphys->ibrem = par;
}
 
 
//_____________________________________________________________________________
void TGeant3::SetCKOV(Int_t par)
{
  //
  //  To control Cerenkov production
  //   par =0 no Cerenkov;
  //       =1 Cerenkov;
  //       =2 Cerenkov with primary stopped at each step.
  //  
  fGctlit->itckov = par;
}
 
 
//_____________________________________________________________________________
void  TGeant3::SetClipBox(const char *name,Float_t xmin,Float_t xmax,
			  Float_t ymin,Float_t ymax,Float_t zmin,Float_t zmax)
{
  //
  //  The hidden line removal technique is necessary to visualize properly
  //  very complex detectors. At the same time, it can be useful to visualize
  //  the inner elements of a detector in detail. This function allows
  //  subtractions (via boolean operation) of BOX shape from any part of
  //  the detector, therefore showing its inner contents.
  //  If "*" is given as the name of the
  //  volume to be clipped, all volumes are clipped by the given box.
  //  A volume can be clipped at most twice.
  //  if a volume is explicitely clipped twice,
  //  the "*" will not act on it anymore. Giving "." as the name
  //  of the volume to be clipped will reset the clipping.
  //  Parameters
  //  NAME   Name of volume to be clipped 
  //  +
  //  XMIN   Lower limit of the Shape X coordinate
  //  XMAX   Upper limit of the Shape X coordinate
  //  YMIN   Lower limit of the Shape Y coordinate
  //  YMAX   Upper limit of the Shape Y coordinate
  //  ZMIN   Lower limit of the Shape Z coordinate
  //  ZMAX   Upper limit of the Shape Z coordinate
  //
  //  This function performs a boolean subtraction between the volume
  //  NAME and a box placed in the MARS according the values of the given
  //  coordinates.
  
  InitHIGZ();
  char vname[5];
  Vname(name,vname);
  setclip(PASSCHARD(vname),xmin,xmax,ymin,ymax,zmin,zmax PASSCHARL(vname));   
} 

//_____________________________________________________________________________
void TGeant3::SetCOMP(Int_t par)
{
  //
  //  To control Compton scattering
  //   par =0 no Compton
  //       =1 Compton. Electron processed.
  //       =2 Compton. No electron stored.
  //  
  //
  fGcphys->icomp = par;
}
  
//_____________________________________________________________________________
void TGeant3::SetCUTS(Float_t cutgam,Float_t cutele,Float_t cutneu,
		      Float_t cuthad,Float_t cutmuo ,Float_t bcute ,
		      Float_t bcutm ,Float_t dcute ,Float_t dcutm ,
		      Float_t ppcutm, Float_t tofmax)
{
  //
  //  CUTGAM   Cut for gammas              D=0.001
  //  CUTELE   Cut for electrons           D=0.001
  //  CUTHAD   Cut for charged hadrons     D=0.01
  //  CUTNEU   Cut for neutral hadrons     D=0.01
  //  CUTMUO   Cut for muons               D=0.01
  //  BCUTE    Cut for electron brems.     D=-1.
  //  BCUTM    Cut for muon brems.         D=-1.
  //  DCUTE    Cut for electron delta-rays D=-1.
  //  DCUTM    Cut for muon delta-rays     D=-1.
  //  PPCUTM   Cut for e+e- pairs by muons D=0.01
  //  TOFMAX   Time of flight cut          D=1.E+10
  //
  //   If the default values (-1.) for       BCUTE ,BCUTM ,DCUTE ,DCUTM
  //   are not modified, they will be set to CUTGAM,CUTGAM,CUTELE,CUTELE
  //   respectively.
  //  If one of the parameters from CUTGAM to PPCUTM included
  //  is modified, cross-sections and energy loss tables must be
  //  recomputed via the function Gphysi.
  //
  fGccuts->cutgam = cutgam;
  fGccuts->cutele = cutele;
  fGccuts->cutneu = cutneu;
  fGccuts->cuthad = cuthad;
  fGccuts->cutmuo = cutmuo;
  fGccuts->bcute  = bcute;
  fGccuts->bcutm  = bcutm;
  fGccuts->dcute  = dcute;
  fGccuts->dcutm  = dcutm;
  fGccuts->ppcutm = ppcutm;
  fGccuts->tofmax = tofmax;   
}

//_____________________________________________________________________________
void TGeant3::SetDCAY(Int_t par)
{
  //
  //  To control Decay mechanism.
  //   par =0 no decays.
  //       =1 Decays. secondaries processed.
  //       =2 Decays. No secondaries stored.
  //  
  fGcphys->idcay = par;
}
 
 
//_____________________________________________________________________________
void TGeant3::SetDEBU(Int_t emin, Int_t emax, Int_t emod)
{
  //
  // Set the debug flag and frequency
  // Selected debug output will be printed from
  // event emin to even emax each emod event
  //
  fGcflag->idemin = emin;
  fGcflag->idemax = emax;
  fGcflag->itest  = emod;
}
 
 
//_____________________________________________________________________________
void TGeant3::SetDRAY(Int_t par)
{
  //
  //  To control delta rays mechanism.
  //   par =0 no delta rays.
  //       =1 Delta rays. secondaries processed.
  //       =2 Delta rays. No secondaries stored.
  //  
  fGcphys->idray = par;
}
 
//_____________________________________________________________________________
void TGeant3::SetHADR(Int_t par)
{
  //
  //  To control hadronic interactions.
  //   par =0 no hadronic interactions.
  //       =1 Hadronic interactions. secondaries processed.
  //       =2 Hadronic interactions. No secondaries stored.
  //  
  fGcphys->ihadr = par;
}
 
//_____________________________________________________________________________
void TGeant3::SetKINE(Int_t kine, Float_t xk1, Float_t xk2, Float_t xk3,
		      Float_t xk4, Float_t xk5, Float_t xk6, Float_t xk7,
		      Float_t xk8, Float_t xk9, Float_t xk10)
{
  //
  // Set the variables in /GCFLAG/ IKINE, PKINE(10)
  // Their meaning is user defined
  //
  fGckine->ikine    = kine;
  fGckine->pkine[0] = xk1;
  fGckine->pkine[1] = xk2;
  fGckine->pkine[2] = xk3;
  fGckine->pkine[3] = xk4;
  fGckine->pkine[4] = xk5;
  fGckine->pkine[5] = xk6;
  fGckine->pkine[6] = xk7;
  fGckine->pkine[7] = xk8;
  fGckine->pkine[8] = xk9;
  fGckine->pkine[9] = xk10;
}
 
//_____________________________________________________________________________
void TGeant3::SetLOSS(Int_t par)
{
  //
  //  To control energy loss.
  //   par =0 no energy loss;
  //       =1 restricted energy loss fluctuations;
  //       =2 complete energy loss fluctuations;
  //       =3 same as 1;
  //       =4 no energy loss fluctuations.
  //  If the value ILOSS is changed, then cross-sections and energy loss
  //  tables must be recomputed via the command 'PHYSI'.
  //  
  fGcphys->iloss = par;
}
 
 
//_____________________________________________________________________________
void TGeant3::SetMULS(Int_t par)
{
  //
  //  To control multiple scattering.
  //   par =0 no multiple scattering.
  //       =1 Moliere or Coulomb scattering.
  //       =2 Moliere or Coulomb scattering.
  //       =3 Gaussian scattering.
  //  
  fGcphys->imuls = par;
}
 
 
//_____________________________________________________________________________
void TGeant3::SetMUNU(Int_t par)
{
  //
  //  To control muon nuclear interactions.
  //   par =0 no muon-nuclear interactions.
  //       =1 Nuclear interactions. Secondaries processed.
  //       =2 Nuclear interactions. Secondaries not processed.
  //  
  fGcphys->imunu = par;
}
 
//_____________________________________________________________________________
void TGeant3::SetOPTI(Int_t par)
{
  //
  //  This flag controls the tracking optimisation performed via the
  //  GSORD routine:
  //      1 no optimisation at all; GSORD calls disabled;
  //      0 no optimisation; only user calls to GSORD kept;
  //      1 all non-GSORDered volumes are ordered along the best axis;
  //      2 all volumes are ordered along the best axis.
  //  
  fGcopti->ioptim = par;
}
 
//_____________________________________________________________________________
void TGeant3::SetPAIR(Int_t par)
{
  //
  //  To control pair production mechanism.
  //   par =0 no pair production.
  //       =1 Pair production. secondaries processed.
  //       =2 Pair production. No secondaries stored.
  //  
  fGcphys->ipair = par;
}
 
 
//_____________________________________________________________________________
void TGeant3::SetPFIS(Int_t par)
{
  //
  //  To control photo fission mechanism.
  //   par =0 no photo fission.
  //       =1 Photo fission. secondaries processed.
  //       =2 Photo fission. No secondaries stored.
  //  
  fGcphys->ipfis = par;
}
  
//_____________________________________________________________________________
void TGeant3::SetPHOT(Int_t par)
{
  //
  //  To control Photo effect.
  //   par =0 no photo electric effect.
  //       =1 Photo effect. Electron processed.
  //       =2 Photo effect. No electron stored.
  //  
  fGcphys->iphot = par;
}
 
//_____________________________________________________________________________
void TGeant3::SetRAYL(Int_t par)
{
  //
  //  To control Rayleigh scattering.
  //   par =0 no Rayleigh scattering.
  //       =1 Rayleigh.
  //  
  fGcphys->irayl = par;
}
 
//_____________________________________________________________________________
void TGeant3::SetSWIT(Int_t sw, Int_t val)
{
  //
  //  sw    Switch number
  //  val   New switch value
  //
  //  Change one element of array ISWIT(10) in /GCFLAG/
  //  
  if (sw <= 0 || sw > 10) return;
  fGcflag->iswit[sw-1] = val;
}
 
 
//_____________________________________________________________________________
void TGeant3::SetTRIG(Int_t nevents)
{
  //
  // Set number of events to be run
  //
  fGcflag->nevent = nevents;
}
 
//_____________________________________________________________________________
void TGeant3::SetUserDecay(Int_t pdg)
{
  //
  // Force the decays of particles to be done with Pythia
  // and not with the Geant routines. 
  // just kill pointers doing mzdrop
  //
  Int_t ipart = IdFromPDG(pdg);
  if(ipart<0) {
    printf("Particle %d not in geant\n",pdg);
    return;
  }
  Int_t jpart=fGclink->jpart;
  Int_t jpa=fZlq[jpart-ipart];
  //
  if(jpart && jpa) {
    Int_t jpa1=fZlq[jpa-1];
    if(jpa1)
      mzdrop(fGcbank->ixcons,jpa1,PASSCHARD(" ") PASSCHARL(" "));
    Int_t jpa2=fZlq[jpa-2];
    if(jpa2)
      mzdrop(fGcbank->ixcons,jpa2,PASSCHARD(" ") PASSCHARL(" "));
  }
}

//______________________________________________________________________________
void TGeant3::Vname(const char *name, char *vname)
{
  //
  //  convert name to upper case. Make vname at least 4 chars
  //
  Int_t l = strlen(name);
  Int_t i;
  l = l < 4 ? l : 4;
  for (i=0;i<l;i++) vname[i] = toupper(name[i]);
  for (i=l;i<4;i++) vname[i] = ' ';
  vname[4] = 0;      
}
 
//______________________________________________________________________________
void TGeant3::Ertrgo()
{
  ertrgo();
}

//______________________________________________________________________________
void TGeant3::Ertrak(const Float_t *const x1, const Float_t *const p1, 
			const Float_t *x2, const Float_t *p2,
			Int_t ipa,  Option_t *chopt)
{
  ertrak(x1,p1,x2,p2,ipa,PASSCHARD(chopt) PASSCHARL(chopt));
}
        
//_____________________________________________________________________________
void TGeant3::WriteEuclid(const char* filnam, const char* topvol,
			  Int_t number, Int_t nlevel)
{
  //
  //
  //     ******************************************************************
  //     *                                                                *
  //     *  Write out the geometry of the detector in EUCLID file format  *
  //     *                                                                *
  //     *       filnam : will be with the extension .euc                 *
  //     *       topvol : volume name of the starting node                *
  //     *       number : copy number of topvol (relevant for gsposp)     *
  //     *       nlevel : number of  levels in the tree structure         *
  //     *                to be written out, starting from topvol         *
  //     *                                                                *
  //     *       Author : M. Maire                                        *
  //     *                                                                *
  //     ******************************************************************
  //
  //     File filnam.tme is written out with the definitions of tracking
  //     medias and materials.
  //     As to restore original numbers for materials and medias, program
  //     searches in the file euc_medi.dat and comparing main parameters of
  //     the mat. defined inside geant and the one in file recognizes them
  //     and is able to take number from file. If for any material or medium,
  //     this procedure fails, ordering starts from 1.
  //     Arrays IOTMED and IOMATE are used for this procedure
  //
  const char shape[][5]={"BOX ","TRD1","TRD2","TRAP","TUBE","TUBS","CONE",
			 "CONS","SPHE","PARA","PGON","PCON","ELTU","HYPE",
			 "GTRA","CTUB"};
  Int_t i, end, itm, irm, jrm, k, nmed;
  Int_t imxtmed=0;
  Int_t imxmate=0;
  FILE *lun;
  char *filext, *filetme;
  char natmed[21], namate[21];
  char natmedc[21], namatec[21];
  char key[5], name[5], mother[5], konly[5];
  char card[133];
  Int_t iadvol, iadtmd, iadrot, nwtot, iret;
  Int_t mlevel, numbr, natt, numed, nin, ndata;
  Int_t iname, ivo, ish, jvo, nvstak, ivstak;
  Int_t jdiv, ivin, in, jin, jvin, irot;
  Int_t jtm, imat, jma, flag=0, imatc;
  Float_t az, dens, radl, absl, a, step, x, y, z;
  Int_t npar, ndvmx, left;
  Float_t zc, densc, radlc, abslc, c0, tmaxfd;
  Int_t nparc, numb;
  Int_t iomate[100], iotmed[100];
  Float_t par[50], att[20], ubuf[50];
  Float_t *qws;
  Int_t   *iws;
  Int_t level, ndiv, iaxe;
  Int_t itmedc, nmatc, isvolc, ifieldc, nwbufc, isvol, nmat, ifield, nwbuf;
  Float_t fieldmc, tmaxfdc, stemaxc, deemaxc, epsilc, stminc, fieldm;
  Float_t tmaxf, stemax, deemax, epsil, stmin;
  const char *f10000="!\n%s\n!\n";
  //Open the input file
  end=strlen(filnam);
  for(i=0;i<end;i++) if(filnam[i]=='.') {
    end=i;
    break;
  }
  filext=new char[end+4];
  filetme=new char[end+4];
  strncpy(filext,filnam,end);
  strncpy(filetme,filnam,end);
  //
  // *** The output filnam name will be with extension '.euc'
  strcpy(&filext[end],".euc");
  strcpy(&filetme[end],".tme");
  lun=fopen(filext,"w");
  //
  // *** Initialisation of the working space
  iadvol=fGcnum->nvolum;
  iadtmd=iadvol+fGcnum->nvolum;
  iadrot=iadtmd+fGcnum->ntmed;
  if(fGclink->jrotm) {
    fGcnum->nrotm=fZiq[fGclink->jrotm-2];
  } else {
    fGcnum->nrotm=0;
  }
  nwtot=iadrot+fGcnum->nrotm;
  qws = new float[nwtot+1];
  for (i=0;i<nwtot+1;i++) qws[i]=0;
  iws = (Int_t*) qws;
  mlevel=nlevel;
  if(nlevel==0) mlevel=20;
  //
  // *** find the top volume and put it in the stak
  numbr = number>0 ? number : 1;
  Gfpara(topvol,numbr,1,npar,natt,par,att);
  if(npar <= 0) {
    printf(" *** GWEUCL *** top volume : %s number : %3d can not be a valid root\n",
	   topvol, numbr);
    return;
  }
  //
  // ***  authorized shape ?
  strncpy((char *)&iname, topvol, 4);
  ivo=0;
  for(i=1; i<=fGcnum->nvolum; i++) if(fZiq[fGclink->jvolum+i]==iname) {
    ivo=i;
    break;
  }
  jvo = fZlq[fGclink->jvolum-ivo];
  ish = Int_t (fZq[jvo+2]);
  if(ish > 12) {
    printf(" *** GWEUCL *** top volume : %s number : %3d can not be a valid root\n",
	   topvol, numbr);
  }
  //
  level = 1;
  nvstak = 1;
  iws[nvstak]     = ivo;
  iws[iadvol+ivo] = level;
  ivstak = 0;
  //
  //*** flag all volumes and fill the stak
  //
 L10:
  //
  //    pick the next volume in stak
  ivstak += 1;
  ivo   = TMath::Abs(iws[ivstak]);
  jvo   = fZlq[fGclink->jvolum - ivo];
  //
  //     flag the tracking medium
  numed =  Int_t (fZq[jvo + 4]);
  iws[iadtmd + numed] = 1;
  //
  //    get the daughters ...
  level = iws[iadvol+ivo];
  if (level < mlevel) {
    level +=  1;
    nin = Int_t (fZq[jvo + 3]);
    //
    //       from division ...
    if (nin < 0) {
      jdiv = fZlq[jvo  - 1];
      ivin =  Int_t (fZq[jdiv + 2]);
      nvstak += 1;
      iws[nvstak]      = -ivin;
      iws[iadvol+ivin] =  level;
      //
      //       from position ...
    } else if (nin > 0) {
      for(in=1; in<=nin; in++) {
	jin  = fZlq[jvo - in];
	ivin =  Int_t (fZq[jin + 2 ]);
	jvin = fZlq[fGclink->jvolum - ivin];
	ish  =  Int_t (fZq[jvin + 2]);
	//              authorized shape ?
	if (ish <= 12) {
	  //                 not yet flagged ?
	  if (iws[iadvol+ivin]==0) {
	    nvstak += 1;
	    iws[nvstak]      = ivin;
	    iws[iadvol+ivin] = level;
	  }
	  //                 flag the rotation matrix
	  irot = Int_t ( fZq[jin + 4 ]);
	  if (irot > 0) iws[iadrot+irot] = 1;
	}
      }
    }
  }
  //
  //     next volume in stak ?
  if (ivstak < nvstak) goto L10;
  //
  // *** restore original material and media numbers
  // file euc_medi.dat is needed to compare materials and medias
  //
  FILE* luncor=fopen("euc_medi.dat","r");
  //
  if(luncor) {
    for(itm=1; itm<=fGcnum->ntmed; itm++) {
      if (iws[iadtmd+itm] > 0) {
	jtm = fZlq[fGclink->jtmed-itm];
	strncpy(natmed,(char *)&fZiq[jtm+1],20);
	imat =  Int_t (fZq[jtm+6]);
	jma  = fZlq[fGclink->jmate-imat];
	if (jma <= 0) {
	  printf(" *** GWEUCL *** material not defined for tracking medium %5i %s\n",itm,natmed);
	  flag=1;
	} else {
	  strncpy(namate,(char *)&fZiq[jma+1],20);
	}
	//*
	//** find the material original number
	rewind(luncor);
      L23:
	iret=fscanf(luncor,"%4s,%130s",key,card);
	if(iret<=0) goto L26;
	flag=0;
	if(!strcmp(key,"MATE")) {
	  sscanf(card,"%d %s %f %f %f %f %f %d",&imatc,namatec,&az,&zc,&densc,&radlc,&abslc,&nparc);
	  Gfmate(imat,namate,a,z,dens,radl,absl,par,npar);
	  if(!strcmp(namatec,namate)) {
	    if(az==a && zc==z && densc==dens && radlc==radl 
	       && abslc==absl && nparc==nparc) {
	      iomate[imat]=imatc;
	      flag=1;
	      printf("*** GWEUCL *** material : %3d '%s' restored as %3d\n",imat,namate,imatc);
	    } else {
	      printf("*** GWEUCL *** different definitions for material: %s\n",namate);
	    }
	  }
	}
	if(strcmp(key,"END") && !flag) goto L23;
	if (!flag) {
	  printf("*** GWEUCL *** cannot restore original number for material: %s\n",namate);
	}
	//*
	//*
	//***  restore original tracking medium number
	rewind(luncor);
      L24:
	iret=fscanf(luncor,"%4s,%130s",key,card);
	if(iret<=0) goto L26;
	flag=0;
	if (!strcmp(key,"TMED")) {
	  sscanf(card,"%d %s %d %d %d %f %f %f %f %f %f %d\n",
		 &itmedc,natmedc,&nmatc,&isvolc,&ifieldc,&fieldmc,
		 &tmaxfdc,&stemaxc,&deemaxc,&epsilc,&stminc,&nwbufc);
	  Gftmed(itm,natmed,nmat,isvol,ifield,fieldm,tmaxf,stemax,deemax,
			epsil,stmin,ubuf,&nwbuf);
	  if(!strcmp(natmedc,natmed)) {
	    if (iomate[nmat]==nmatc && nwbuf==nwbufc) {
	      iotmed[itm]=itmedc;
	      flag=1;
	      printf("*** GWEUCL *** medium   : %3d '%20s' restored as %3d\n",
		     itm,natmed,itmedc);
	    } else {
	      printf("*** GWEUCL *** different definitions for tracking medium: %s\n",natmed);
	    }
	  }
	}
	if(strcmp(key,"END") && !flag) goto L24;
	if(!flag) {
	  printf("cannot restore original number for medium : %s\n",natmed);
	  goto L27;
	}
      }
    }
    goto L29;
    //*
  }
 L26:   printf("*** GWEUCL *** cannot read the data file\n");
 L27:   flag=2;
 L29:   if(luncor) fclose (luncor);
  //
  //
  // *** write down the tracking medium definition
  //
  strcpy(card,"!       Tracking medium");
  fprintf(lun,f10000,card);
  //
  for(itm=1;itm<=fGcnum->ntmed;itm++) {
    if (iws[iadtmd+itm]>0) {
      jtm  = fZlq[fGclink->jtmed-itm];
      strncpy(natmed,(char *)&fZiq[jtm+1],20);
      natmed[20]='\0';
      imat =  Int_t (fZq[jtm+6]);
      jma  = fZlq[fGclink->jmate-imat];
      //*  order media from one, if comparing with database failed
      if (flag==2) {
	iotmed[itm]=++imxtmed;
	iomate[imat]=++imxmate;
      }
      //*
      if(jma<=0) {
	strcpy(namate,"                  ");
	printf(" *** GWEUCL *** material not defined for tracking medium %5d %s\n",
	       itm,natmed);
      } else {
	strncpy(namate,(char *)&fZiq[jma+1],20);
	namate[20]='\0';
      }
      fprintf(lun,"TMED %3d '%20s' %3d '%20s'\n",iotmed[itm],natmed,iomate[imat],namate);
    }
  }
  //*
      //* *** write down the rotation matrix
  //*
  strcpy(card,"!       Reperes");
  fprintf(lun,f10000,card);
  //
  for(irm=1;irm<=fGcnum->nrotm;irm++) {
    if (iws[iadrot+irm]>0) {
      jrm  = fZlq[fGclink->jrotm-irm];
      fprintf(lun,"ROTM %3d",irm);
      for(k=11;k<=16;k++) fprintf(lun," %8.3f",fZq[jrm+k]);
      fprintf(lun,"\n");
    }
  }
  //*
  //* *** write down the volume definition
  //*
  strcpy(card,"!       Volumes");
  fprintf(lun,f10000,card);
  //*
  for(ivstak=1;ivstak<=nvstak;ivstak++) {
    ivo = iws[ivstak];
    if (ivo>0) {
      strncpy(name,(char *)&fZiq[fGclink->jvolum+ivo],4);
      name[4]='\0';
      jvo  = fZlq[fGclink->jvolum-ivo];
      ish   = Int_t (fZq[jvo+2]);
      nmed  = Int_t (fZq[jvo+4]);
      npar  = Int_t (fZq[jvo+5]);
      if (npar>0) {
	if (ivstak>1) for(i=0;i<npar;i++) par[i]=fZq[jvo+7+i];
	Gckpar (ish,npar,par);
	fprintf(lun,"VOLU '%4s' '%4s' %3d %3d\n",name,shape[ish-1],iotmed[nmed],npar);
	for(i=0;i<(npar-1)/6+1;i++) {
	  fprintf(lun,"     ");
	  left=npar-i*6;
	  for(k=0;k<(left<6?left:6);k++) fprintf(lun," %11.5f",par[i*6+k]);
	  fprintf(lun,"\n");
	}
      } else {
	fprintf(lun,"VOLU '%4s' '%4s' %3d %3d\n",name,shape[ish-1],iotmed[nmed],npar);
      }
    }
  }
  //*
  //* *** write down the division of volumes
  //*
  fprintf(lun,f10000,"!       Divisions");
  for(ivstak=1;ivstak<=nvstak;ivstak++) {
    ivo = TMath::Abs(iws[ivstak]);
    jvo  = fZlq[fGclink->jvolum-ivo];
    ish  = Int_t (fZq[jvo+2]);
    nin  = Int_t (fZq[jvo+3]);
    //*        this volume is divided ...
    if (nin<0) {
      jdiv = fZlq[jvo-1];
      iaxe = Int_t ( fZq[jdiv+1]);
      ivin = Int_t ( fZq[jdiv+2]);
      ndiv = Int_t ( fZq[jdiv+3]);
      c0   =  fZq[jdiv+4];
      step =  fZq[jdiv+5];
      jvin = fZlq[fGclink->jvolum-ivin];
      nmed = Int_t ( fZq[jvin+4]);
      strncpy(mother,(char *)&fZiq[fGclink->jvolum+ivo ],4);
      mother[4]='\0';
      strncpy(name,(char *)&fZiq[fGclink->jvolum+ivin],4);
      name[4]='\0';
      if ((step<=0.)||(ish>=11)) {
	//*              volume with negative parameter or gsposp or pgon ...
	fprintf(lun,"DIVN '%4s' '%4s' %3d %3d\n",name,mother,ndiv,iaxe);
      } else if ((ndiv<=0)||(ish==10)) {
	//*              volume with negative parameter or gsposp or para ...
	ndvmx = TMath::Abs(ndiv);
	fprintf(lun,"DIVT '%4s' '%4s' %11.5f %3d %3d %3d\n",
		name,mother,step,iaxe,iotmed[nmed],ndvmx);
      } else {
	//*              normal volume : all kind of division are equivalent
	fprintf(lun,"DVT2 '%4s' '%4s' %11.5f %3d %11.5f %3d %3d\n",
		name,mother,step,iaxe,c0,iotmed[nmed],ndiv);
      }
    }
  }
  //*
  //* *** write down the the positionnement of volumes
  //*
  fprintf(lun,f10000,"!       Positionnements\n");
  //
  for(ivstak = 1;ivstak<=nvstak;ivstak++) {
    ivo = TMath::Abs(iws[ivstak]);
    strncpy(mother,(char*)&fZiq[fGclink->jvolum+ivo ],4);
    mother[4]='\0';
    jvo  = fZlq[fGclink->jvolum-ivo];
    nin  = Int_t( fZq[jvo+3]);
    //*        this volume has daughters ...
    if (nin>0) {
      for (in=1;in<=nin;in++) {
	jin  = fZlq[jvo-in];
	ivin =  Int_t (fZq[jin +2]);
	numb =  Int_t (fZq[jin +3]);
	irot =  Int_t (fZq[jin +4]);
	x    =  fZq[jin +5];
	y    =  fZq[jin +6];
	z    =  fZq[jin +7];
	strcpy(konly,"ONLY");
	if (fZq[jin+8]!=1.) strcpy(konly,"MANY");
	strncpy(name,(char*)&fZiq[fGclink->jvolum+ivin],4);
	name[4]='\0';
	jvin = fZlq[fGclink->jvolum-ivin];
	ish  = Int_t (fZq[jvin+2]);
	//*              gspos or gsposp ?
	ndata = fZiq[jin-1];
	if (ndata==8) {
	  fprintf(lun,"POSI '%4s' %4d '%4s' %11.5f %11.5f %11.5f %3d '%4s'\n",
		  name,numb,mother,x,y,z,irot,konly);
	} else {
	  npar =  Int_t (fZq[jin+9]);
	  for(i=0;i<npar;i++) par[i]=fZq[jin+10+i];
	  Gckpar (ish,npar,par);
	  fprintf(lun,"POSP '%4s' %4d '%4s' %11.5f %11.5f %11.5f %3d '%4s' %3d\n",
		  name,numb,mother,x,y,z,irot,konly,npar);
	  fprintf(lun,"     ");
	  for(i=0;i<npar;i++) fprintf(lun," %11.5f",par[i]);
	  fprintf(lun,"\n");
	}
      }
    }
  }
  //*
  fprintf(lun,"END\n");
  fclose(lun);
  //*
  //****** write down the materials and medias *****
  //*
  lun=fopen(filetme,"w");
  //*
  for(itm=1;itm<=fGcnum->ntmed;itm++) {
    if (iws[iadtmd+itm]>0) {
      jtm  = fZlq[fGclink->jtmed-itm];
      strncpy(natmed,(char*)&fZiq[jtm+1],4);
      imat =  Int_t (fZq[jtm+6]);
      jma  =  Int_t (fZlq[fGclink->jmate-imat]);
      //*  material
      Gfmate (imat,namate,a,z,dens,radl,absl,par,npar);
      fprintf(lun,"MATE %4d '%20s'%11.5E %11.5E %11.5E %11.5E %11.5E %3d\n",
	     iomate[imat],namate,a,z,dens,radl,absl,npar);
      //*
      if (npar>0) {
	  fprintf(lun,"     ");
	  for(i=0;i<npar;i++) fprintf(lun," %11.5f",par[i]);
	  fprintf(lun,"\n");
      }
      //*  medium
      Gftmed(itm,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin,par,&npar);
      fprintf(lun,"TMED %4d '%20s' %3d %1d %3d %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %3d\n",
	      iotmed[itm],natmed,iomate[nmat],isvol,ifield,
	      fieldm,tmaxfd,stemax,deemax,epsil,stmin,npar);
      //*
      if (npar>0) {
	  fprintf(lun,"     ");
	  for(i=0;i<npar;i++) fprintf(lun," %11.5f",par[i]);
	  fprintf(lun,"\n");
      }
      
    }
  }
  fprintf(lun,"END\n");
  printf(" *** GWEUCL *** file: %s is now written out\n",filext);
  printf(" *** GWEUCL *** file: %s is now written out\n",filetme);
  // Clean up
  delete [] filext;
  delete [] filetme;
  delete [] qws;
  iws=0;
  return;
}

//_____________________________________________________________________________
void TGeant3::Streamer(TBuffer &R__b)
{
  //
  // Stream an object of class TGeant3.
  //
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    AliMC::Streamer(R__b);
    R__b >> fNextVol;
    R__b >> fNPDGCodes;
    R__b.ReadStaticArray(fPDGCode);
  } else {
    R__b.WriteVersion(TGeant3::IsA());
    AliMC::Streamer(R__b);
    R__b << fNextVol;
    R__b << fNPDGCodes;
    R__b.WriteArray(fPDGCode, fNPDGCodes);
  }
}


