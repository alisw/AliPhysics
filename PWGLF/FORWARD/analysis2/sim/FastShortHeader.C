#ifndef FASTSHORTHEADER_C
#define FASTSHORTHEADER_C
#ifndef __CINT__
# include <TString.h>
#endif

//====================================================================
/** 
 * Header structure 
 */
struct FastShortHeader {
  UInt_t   fRunNo;
  UInt_t   fEventId;
  UInt_t   fNtgt;
  UInt_t   fNproj;
  UInt_t   fNbin;
  UInt_t   fType;
  Double_t fIpX;
  Double_t fIpY;
  Double_t fIpZ;
  Double_t fB;
  Double_t fC;
  Double_t fPhiR;
  UInt_t   fNSpecNproj;  // # of spectator neutrons in projectile
  UInt_t   fNSpecNtgt;   // # of spectator neutrons in target 
  UInt_t   fNSpecPproj;  // # of spectator protons in projectile
  UInt_t   fNSpecPtgt;   // # of spectator protons in target
  
  void Print()
  {
    Printf(" Run #/Event:          %9d/%9d", fRunNo, fEventId);
    Printf(" Participants/binary:  %4d/%4d/%3d", fNtgt, fNproj, fNbin);
    Printf(" Event type:           %7s%12s",(fType==1?"Non":
					     fType==2?"Single":
					     "Double"), "-diffractive");
    Printf(" IP:                   (%-5.1f,%-5.1f,%-5.1f)",fIpX,fIpY,fIpZ);
    Printf(" Impact par./cent.:    (%13f/%-3d)", fB, Int_t(fC));
    Printf(" Reaction plane:       %19f", fPhiR);
    Printf(" Specs (Nt,Np,Pt,Pp):  %4d/%4d/%4d/%4d",
	   fNSpecNtgt, fNSpecNproj, fNSpecPtgt, fNSpecPproj);
  }
  void Clear(Option_t* option="")
  {
    Reset(0,0);
  }
  void Reset(UInt_t runNo, UInt_t eventNo)
  {
    fRunNo      = runNo;
    fEventId    = eventNo;
    fIpX        = 1024;
    fIpY        = 1024;
    fIpZ        = 1024;
    fNtgt       = -1;
    fNproj      = -1;
    fNbin       = -1;
    fPhiR       = -1;
    fB          = -1;
    fC          = -1;
    fNSpecNtgt  = -1;
    fNSpecNproj = -1;
    fNSpecPtgt  = -1;
    fNSpecPproj = -1;
    fEG         = kUnknown;
  }

  enum EG_t {
    kUnknown = 0,
    kPythia,
    kHijing,
    kDPMJet,
    kEPOS,
  } fEG;
};
#endif
//
// EOF
//
