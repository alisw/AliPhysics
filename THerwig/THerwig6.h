#ifndef THERWIG6_H
#define THERWIG6_H

// declaration of c++ Class THerwig6 to be used in ROOT
// this is a c++ interface to the F77 Herwig6 program
// author: j. g. contreras jgcn@moni.mda.cinvestav.mx
// date: december 22, 2000

/*

 Class THerwig6 is an interface to the Herwig program

C-----------------------------------------------------------------------
C                           H E R W I G
C
C            a Monte Carlo event generator for simulating
C        +---------------------------------------------------+
C        | Hadron Emission Reactions With Interfering Gluons |
C        +---------------------------------------------------+
C I.G. Knowles(*), G. Marchesini(+), M.H. Seymour($) and B.R. Webber(#)
C-----------------------------------------------------------------------
C with Minimal Supersymmetric Standard Model Matrix Elements by
C                  S. Moretti($) and K. Odagiri($)
C-----------------------------------------------------------------------
C R parity violating Supersymmetric Decays and Matrix Elements by
C                          P. Richardson(&)
C-----------------------------------------------------------------------
C matrix element corrections to top decay and Drell-Yan type processes
C                         by G. Corcella(+)
C-----------------------------------------------------------------------
C Deep Inelastic Scattering and Heavy Flavour Electroproduction by
C                  G. Abbiendi(@) and L. Stanco(%)
C-----------------------------------------------------------------------
C and Jet Photoproduction in Lepton-Hadron Collisions by J. Chyla(~)
C-----------------------------------------------------------------------
C(*)  Department of Physics & Astronomy, University of Edinburgh
C(+)  Dipartimento di Fisica, Universita di Milano
C($)  Rutherford Appleton Laboratory
C(#)  Cavendish Laboratory, Cambridge
C(&)  Department of Physics, University of Oxford
C(@)  Dipartimento di Fisica, Universita di Bologna
C(%)  Dipartimento di Fisica, Universita di Padova
C(~)  Institute of Physics, Prague
C-----------------------------------------------------------------------
C                  Version 6.100 - 16th December 1999
C-----------------------------------------------------------------------
C Main reference:
C    G.Marchesini,  B.R.Webber,  G.Abbiendi,  I.G.Knowles,  M.H.Seymour,
C    and L.Stanco, Computer Physics Communications 67 (1992) 465.
C-----------------------------------------------------------------------
C Please send e-mail about  this program  to one of the  authors at the
C following Internet addresses:
C    I.Knowles@ed.ac.uk        Giuseppe.Marchesini@mi.infn.it
C    M.Seymour@rl.ac.uk        webber@hep.phy.cam.ac.uk
C-----------------------------------------------------------------------
*/

/* declarations from ROOT */
#include "TGenerator.h"

typedef enum
{
   kHwCharm         =  1704,
   kHwBeauty        =  1705,
   kHwCharmMCATNLO  = -1704,
   kHwBeautyMCATNLO = -1705,
   kHwJetsMCATNLO   = -1396
} Process_t;

class TObjArray;




/* THerwig6 class declaration */
class THerwig6 : public TGenerator {
//----------------------------------------------------------------------------
//  functions:
//----------------------------------------------------------------------------
public:
				// ****** constructors and destructor
  THerwig6();
  THerwig6(const THerwig6 & source);
  THerwig6 & operator=(const THerwig6 & /*source*/) {
    Fatal("THerwig6","Assignment operator not implemented yet");
    return *this;
  }
  virtual ~THerwig6();
  static THerwig6 *Instance();

  // acces to hep common block
  int         GetNEVHEP        () const;
  int         GetNhep          () const;
  int         GetISTHEP    (int i)const;
  int         GetIDHEP     (int i)const;
  int         GetJMOHEP (int i, int j) const;
  int         GetJDAHEP (int i, int j) const;
  double      GetPHEP   (int i, int j) const;
  double      GetVHEP   (int i, int j) const;
  int         GetIPART1        () const;
  int         GetIPART2        () const;
  char*       GetPART1         () const;
  char*       GetPART2         () const;
  double      GetEBEAM1        () const;
  double      GetEBEAM2        () const;
  double      GetPBEAM1        () const;
  double      GetPBEAM2        () const;
  int         GetIPROC         () const;
  int         GetMAXEV         () const;
  double      GetQCDLAM        () const;    
  void        SetQCDLAM   (double q) const;       
  double      GetVQCUT         () const;    
  void        SetVQCUT    (double v) const;       
  double      GetVGCUT         () const;    
  void        SetVGCUT    (double v) const;       
  double      GetVPCUT         () const;    
  void        SetVPCUT    (double v) const;
  double      GetCLMAX         () const;    
  void        SetCLMAX    (double c) const;       
  double      GetCLPOW         () const;    
  void        SetCLPOW    (double c) const;       
  double      GetPSPLT    (int i) const;    
  void        SetPSPLT    (int i, double p) const;
  double      GetQDIQK         () const;
  void        SetQDIQK    (double q) const;
  double      GetPDIQK         () const;
  void        SetPDIQK    (double p) const;   
  double      GetQSPAC         () const;
  void        SetQSPAC    (double q) const;   
  double      GetPTRMS         () const;
  void        SetPTRMS    (double p) const;   
  double      GetENSOF         () const;
  void        SetENSOF    (double e) const;   
  int         GetIPRINT        () const;
  void        SetIPRINT   (int i) const;      
  int         GetMODPDF   (int i) const;
  void        SetMODPDF   (int i, int j) const;
  int         GetNSTRU         () const; 
  void        SetNSTRU    (int i) const;       
  char*       GetAUTPDF     (int i) const;         
  void        SetAUTPDF(int i,const char* s) const;
  char*       GetBDECAY        () const;           
  double      GetAVWGT         () const;
  int         GetMAXPR         () const;
  void        SetMAXPR    (int i) const;      
  int         GetMAXER         () const;
  void        SetMAXER    (int i) const;      
  int         GetNRN      (int i) const;
  void        SetNRN    (int i, int j) const; 
  double      GetEVWGT         () const;

  int         GetIDHW     (int i) const;

  int         GetIERROR        () const;

  // /HWHARD/
  double      GetPTMIN         () const;
  void        SetPTMIN    (double d) const;
  double      GetPTMAX         () const;
  void        SetPTMAX    (double d) const;
  double      GetPTPOW         () const;
  void        SetPTPOW    (double d) const;
  double      GetYJMIN         () const;
  void        SetYJMIN    (double d) const;
  double      GetYJMAX         () const;
  void        SetYJMAX    (double d) const;
  double      GetQ2MIN         () const;
  void        SetQ2MIN    (double d) const;
  double      GetQ2MAX         () const;
  void        SetQ2MAX    (double d) const;
  double      GetYBMIN         () const;
  void        SetYBMIN    (double d) const;
  double      GetYBMAX         () const;
  void        SetYBMAX    (double d) const;
  double      GetZJMAX        ()  const;
  void        SetZJMAX    (double d) const;
  int         GetIHPRO        () const;
  
  // /HWPROP/
  double      GetRMASS      (int i) const;
  void        SetRMASS    (int i, double r) const;


  void        GetRNAME (int i, char a[9]) const;

  // Herwig6 routines
  // the user would call
  //   Initialize
  //   change by himself the parameters s/he wants
  //   Hwusta to make stable the particles s/he wants
  //   PrepareRun
  //   GenerateEvent as many times as wished
  // An example is given in SetupTest

  void             GenerateEvent();
  void             Initialize(const char *beam, const char *target, double pbeam1, double pbeam2, int iproc);
  void             InitializeJimmy(const char *beam, const char *target, double pbeam1, double pbeam2, int iproc);
  void             PrepareRun();
  void             PrepareRunJimmy();
  void             OpenFortranFile(int lun, char* name);
  void             CloseFortranFile(int lun);
  Int_t            ImportParticles(TClonesArray *particles, Option_t *option="");
  TObjArray       *ImportParticles(Option_t *option="");
  TObjArray       *Particles() const { return fParticles; }
  void             Hwigin();
  void             Hwuinc();
  void             Hwusta(const char * name);
  void             Hweini();
  void             Hwuine();
  void             Hwepro();
  void             Hwbgen();
  void             Hwdhob();
  void             Hwcfor();
  void             Hwcdec();
  void             Hwdhad();
  void             Hwdhvy();
  void             Hwmevt();
  void             Hwufne();
  void             Hwefin();
  void             Hwiodk(int iopt);
  void             SetupTest();
  void             PrintEvt();
 // Jimmy subroutines:
  void             Jminit();
  void             Jimmin();
  void             Jmefin();
protected:
  static   THerwig6 *fgInstance; // Singlet instance
  ClassDef(THerwig6,0)  //Interface to Herwig6.1 Event Generator
};

#endif
