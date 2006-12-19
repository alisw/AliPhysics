#ifdef COMPILING
#include "ITS/AliITSvPPRasymmFMD.h"	// ITS
#include "STRUCT/AliDIPOv2.h"		// DIPO
#include "STRUCT/AliFRAMEv2.h"		// FRAME
#include "STRUCT/AliSHILv2.h"		// SHIL
#include "TPC/AliTPCv2.h"		// TPC
#include "TOF/AliTOFv4T0.h"		// TOF
#include "HMPID/AliHMPIDv1.h"		// HMPID
#include "ZDC/AliZDCv2.h"		// ZDC
#include "TRD/AliTRDv1.h"		// TRD
#include "MUON/AliMUONv1.h"		// MUON
#include "PHOS/AliPHOSv1.h"		// PHOS
#include "PMD/AliPMDv1.h"		// PMD
#include "T0/AliT0v1.h"			// T0
#include "EMCAL/AliEMCALv1.h"		// EMCAL
#include "VZERO/AliVZEROv3.h"		// VZERO

template <typename T>
struct Dummy : public T 
{
  Dummy() : T() {}
  Dummy(const char* n) : T(n, Form("%s dummy", n)) {}
  void StepManager() {}
  ClassDef(Dummy, 1);
};

typedef Dummy<AliDIPOv2>          DummyDIPO;
typedef Dummy<AliFRAMEv2>         DummyFRAME;
typedef Dummy<AliSHILv2>          DummySHIL;
typedef Dummy<AliITSvPPRasymmFMD> DummyITS;
typedef Dummy<AliTPCv2>           DummyTPC;
typedef Dummy<AliTOFv4T0>	  DummyTOF;
typedef Dummy<AliHMPIDv1>         DummyHMPID;
typedef Dummy<AliZDCv2>           DummyZDC;
typedef Dummy<AliTRDv1>           DummyTRD;
typedef Dummy<AliMUONv1>          DummyMUON;
typedef Dummy<AliPHOSv1>          DummyPHOS;
typedef Dummy<AliPMDv1>           DummyPMD;
typedef Dummy<AliT0v1>            DummyT0;
typedef Dummy<AliEMCALv1>         DummyEMCAL;
typedef Dummy<AliVZEROv3>         DummyVZERO;
#endif


  
//
// EOF
//
