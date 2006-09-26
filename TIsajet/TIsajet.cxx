#include "TParticle.h"
#include "TSystem.h"
#include "TIsajet.h"
#include "Icommon.h"
#include "Riostream.h"
#include "math.h"
#include "TROOT.h"
#include "TMath.h"

ClassImp(TIsajet)

static TRandom * sRandom;

/**************************************************************************/

TIsajet::TIsajet() : TGenerator("Isajet", "Isajet")
{
//  Default constructor        
//  Set random number
    if (!sRandom) sRandom=fRandom;

// Initialising equivalence structures in FRGPAR :
// EQUIVALENCE (PMIX1(1,1),PMIXX1(1))
// EQUIVALENCE (PMIX2(1,1),PMIXX2(1))
// EQUIVALENCE (FRPAR(1),PUD)

    FRGPAR.pmixx1[0] = &FRGPAR.pmix1[0][0];
    FRGPAR.pmixx2[0] = &FRGPAR.pmix2[0][0];
    FRGPAR.frpar[0] = &FRGPAR.pud;

    for (Int_t i = 1; i < 6; i++) {
	FRGPAR.pmixx1[i] = FRGPAR.pmixx1[i-1] + 1;
	FRGPAR.pmixx1[i] = FRGPAR.pmixx1[i-1] + 1;
    }
    
    for (Int_t i = 1; i < 32; i++) {
	FRGPAR.frpar[i] = FRGPAR.frpar[i-1] +1;
    }

// Internal logical flags ensuring that theta and y limits are not simultaneously set :

    DYLIM.ywset = DYLIM.thwset = false;

// Defaults :

    title = ("Default title.");
    jobtype = ("TWOJET");

    RestoreDefaults();
}

/**************************************************************************/

TIsajet::~TIsajet()
{
    // No housekeeping required at the moment.
}

/**************************************************************************/

void TIsajet::Initialise() 
{
    
// Writes parameter file and stores common block variables in
// data members, according to booleans.
// If TIsajet is being used in online-control mode, the parameter file is
// unnecessary, hence the output is echoed to the screen instead. 

    const char *fname =  gSystem->ExpandPathName("$ALICE_ROOT/ISAJET/data/myjob.par");
    ofstream Write(fname, ios::out);

    ostream *Writer = &Write;

    if (online) Writer = &cout;

    *Writer << title << '\n';
    *Writer << PRIMAR.ecm << ',' << PRIMAR.nevent << ",1,1/\n";
    *Writer << jobtype << '\n';

    center_energy = PRIMAR.ecm;

    if (setBeams) {
	beam1_type = PRIMAR.idin[0];
	beam2_type = PRIMAR.idin[1];

	*Writer << "BEAMS\n";
	if (PRIMAR.idin[0] == -1120) *Writer << "AP,";
	else if (PRIMAR.idin[0] == 1220) *Writer << "N,";
	else if (PRIMAR.idin[0] == -1220) *Writer << "AN,";
	else *Writer << "P,";

	if (PRIMAR.idin[1] == -1120) *Writer << "AP/\n";
	else if (PRIMAR.idin[1] == 1220) *Writer << "N/\n";
	else if (PRIMAR.idin[1] == -1220) *Writer << "AN/\n";
	else *Writer << "P/\n";
    }

    if (setCutjet) {
	cutoff_mass = QCDPAR.cutjet;
	*Writer << "CUTJET\n" << QCDPAR.cutjet << "/\n";
    }
    
    if (setFragment) {
	for (Int_t i = 0; i < 32; i++) frag_params[i] = *FRGPAR.frpar[i];

	*Writer << "FRAGMENT\n";
	for (Int_t i = 0; i < 31; i++) *Writer << FRGPAR.frpar[i] << ',';
	*Writer << FRGPAR.frpar[31] << "/\n";
    }
    
    if (setJettype1) {
	for (Int_t i = 0; i < TYPES.njttyp[0]; i++) jet1_type[i] = XTYPES.jetyp[0][i];
	num_jet_type[0] = TYPES.njttyp[0];

	*Writer << "JETTYPE1\n";
	for (Int_t i = 0; i < TYPES.njttyp[0]-1; i++) *Writer << XTYPES.jetyp[0][i] << ',';
	*Writer << XTYPES.jetyp[0][TYPES.njttyp[0]-1] << "/\n";
    }
    if (setJettype2) {
	for (Int_t i = 0; i < TYPES.njttyp[1]; i++) jet2_type[i] = XTYPES.jetyp[1][i];
	num_jet_type[0] = TYPES.njttyp[0];

	*Writer << "JETTYPE2\n";
	for (Int_t i = 0; i < TYPES.njttyp[1]-1; i++) *Writer << XTYPES.jetyp[1][i] << ',';
	*Writer << XTYPES.jetyp[1][TYPES.njttyp[1]-1] << "/\n";
    }
    if (setJettype3) {
	for (Int_t i = 0; i < TYPES.njttyp[2]; i++) jet3_type[i] = XTYPES.jetyp[2][i];
	num_jet_type[0] = TYPES.njttyp[0];

	*Writer << "JETTYPE3\n";
	for (Int_t i = 0; i < TYPES.njttyp[2]-1; i++) *Writer << XTYPES.jetyp[2][i] << ',';
	*Writer << XTYPES.jetyp[2][TYPES.njttyp[2]-1] << "/\n";
    }
    
	
    if (setLambda) {
	qcd_lambda = QCDPAR.alam;
	*Writer << "LAMBDA\n" << QCDPAR.alam << "/\n";
    }
    

    if (setNodcay) {
	forbid_decay = NODCAY.nodcay;

	*Writer << "NODCAY\n";
	if (NODCAY.nodcay) *Writer << "TRUE/\n";
	else *Writer << "FALSE/\n";
    }

    if (setNoeta) {
	forbid_eta = NODCAY.noeta;

	*Writer << "NOETA\n";
	if (NODCAY.noeta) *Writer << "TRUE/\n";
	else *Writer << "FALSE/\n";
    }

    if (setNoevolve) {
	forbid_evolve = NODCAY.noevol;
	
	*Writer << "NOEVOLVE\n";
	if (NODCAY.noevol) *Writer << "TRUE/\n";
	else *Writer << "FALSE/\n";
    }

    if (setNohadron) {
	forbid_hadron = NODCAY.nohadr;
	
	*Writer << "NOHADRON\n";
	if (NODCAY.nohadr) *Writer << "TRUE/\n";
	else *Writer << "FALSE/\n";
    }

    if (setNopi0) {
	forbid_pi0 = NODCAY.nopi0;
	
	*Writer << "NOPI0\n";
	if (NODCAY.nopi0) *Writer << "TRUE/\n";
	else *Writer << "FALSE/\n";
    }
	
    if (setNsigma) {
	generate_sigma = PRIMAR.nsigma;
	*Writer << "NSIGMA\n" << PRIMAR.nsigma << "/\n";
    }    
    
    if (setP) {
	for (Int_t i = 0; i < 3; i++) {
	    p_limits[2 * i] = JETLIM.pmin[i];
	    p_limits[2 * i + 1] = JETLIM.pmax[i];
	}	

	*Writer << "P\n";
	*Writer << JETLIM.pmin[0] << ',' << JETLIM.pmax[0] << ',';
	*Writer << JETLIM.pmin[1] << ',' << JETLIM.pmax[1] << ',';
	*Writer << JETLIM.pmin[2] << ',' << JETLIM.pmax[2] << "/\n";
    }

    if (setPhi) {
	for (Int_t i = 0; i < 3; i++) {
	    phi_limits[2 * i] = JETLIM.phimin[i];
	    phi_limits[2 * i + 1] = JETLIM.phimax[i];
	}

	*Writer << "PHI\n";
	*Writer << JETLIM.phimin[0] << ',' << JETLIM.phimax[0] << ',';
	*Writer << JETLIM.phimin[1] << ',' << JETLIM.phimax[1] << ',';
	*Writer << JETLIM.phimin[2] << ',' << JETLIM.phimax[2] << "/\n";
    }
    
    if (setPt) {
	for (Int_t i = 0; i < 3; i++) {
	    pt_limits[2 * i] = JETLIM.ptmin[i];
	    pt_limits[2 * i + 1] = JETLIM.ptmax[i];
	}

	*Writer << "PT\n";
	*Writer << JETLIM.ptmin[0] << ',' << JETLIM.ptmax[0] << ',';
	*Writer << JETLIM.ptmin[1] << ',' << JETLIM.ptmax[1] << ',';
	*Writer << JETLIM.ptmin[2] << ',' << JETLIM.ptmax[2] << "/\n";
    }

    if (setTheta) {
	for (Int_t i = 0; i < 3; i++) {
	    theta_limits[2 * i] = JETLIM.thmin[i];
	    theta_limits[2 * i + 1] = JETLIM.thmax[i];
	}

	*Writer << "THETA\n";
	*Writer << JETLIM.thmin[0] << ',' << JETLIM.thmax[0] << ',';
	*Writer << JETLIM.thmin[1] << ',' << JETLIM.thmax[1] << ',';
	*Writer << JETLIM.thmin[2] << ',' << JETLIM.thmax[2] << "/\n";
    }

    if (setX) {
	for (Int_t i = 0; i < 3; i++) {
	    x_limits[2 * i] = JETLIM.xjmin[i];
	    x_limits[2 * i + 1] = JETLIM.xjmax[i];
	}

	*Writer << "X\n";
	*Writer << JETLIM.xjmin[0] << ',' << JETLIM.xjmax[0] << ',';
	*Writer << JETLIM.xjmin[1] << ',' << JETLIM.xjmax[1] << ',';
	*Writer << JETLIM.xjmin[2] << ',' << JETLIM.xjmax[2] << "/\n";
    }

    if (setY) {
	for (Int_t i = 0; i < 3; i++) {
	    y_limits[2 * i] = JETLIM.yjmin[i];
	    y_limits[2 * i + 1] = JETLIM.yjmax[i];
	}

	*Writer << "Y\n";
	*Writer << JETLIM.yjmin[0] << ',' << JETLIM.yjmax[0] << ',';
	*Writer << JETLIM.yjmin[1] << ',' << JETLIM.yjmax[1] << ',';
	*Writer << JETLIM.yjmin[2] << ',' << JETLIM.yjmax[2] << "/\n";
    }

    if (setXgen) {
	for (Int_t i = 0; i < 8; i++) peter_jet_frag[i] = FRGPAR.xgen[i];

	*Writer << "XGEN\n";
	for (Int_t i = 0; i < 7; i++) *Writer << FRGPAR.xgen[i] << ',';
	*Writer << FRGPAR.xgen[7] << "/\n";
    }

    if (setPdf) {
	*Writer << "PDFLIB\n";
	for (Int_t i = 0; i < num_Pdf; i++) *Writer << "\'" << pdfpar[i] << "\'" << ',' << pdfval[i] << "/\n";
    }
    

    *Writer << "END\n";
    *Writer << "STOP\n";
    Write.close();

//  Stuff for online-control mode :

    if (online) {
	KEYS.reac = jobtype;
	KEYS.keyon = false;
	for (Int_t i = 0; i < KEYS.mxkeys; i++) KEYS.keys[i] = false;
	
	if (!strcmp(KEYS.reac, "TWOJET")) {
	    KEYS.keys[0] = true;
	    KEYS.ikey = 1;
	    PRIMAR.njet = 2;
	}
	else if (!strcmp(KEYS.reac, "MINBIAS")) {
	    KEYS.keys[3] = true;
	    KEYS.ikey = 4;
	    PRIMAR.njet = 0;
	}
	else {
	    printf("Error in TIsajet::Initialise :\n");
	    printf("Invalid job type %s.\n", KEYS.reac);
	    printf("Only TWOJET and MINBIAS are currently supported for online mode.\n");
	    return;
	}

	if (setPdf) {
//	    PDFinit();
	}
    }
}

/**************************************************************************/

void TIsajet::Reload() 
{
//
// Sets the common block variables to the data member values.
//

    SetECM(center_energy);

    if (setBeams) {
	SetIDIN(0, beam1_type);
	SetIDIN(1, beam2_type);
    }

    if (setCutjet) SetCUTJET(cutoff_mass);
    
    if (setFragment) SetAllFRPAR(frag_params, 32);
    
    if (setJettype1) for (Int_t i = 0; i < num_jet_type[0]; i++) SetJETYP(0, jet1_type[i]);
    
    if (setJettype2) for (Int_t i = 0; i < num_jet_type[1]; i++) SetJETYP(1, jet2_type[i]);
    
    if (setJettype3) for (Int_t i = 0; i < num_jet_type[2]; i++) SetJETYP(2, jet3_type[i]);
	
    if (setLambda) SetALAM(qcd_lambda);

    if (setNodcay) SetNODCAY(forbid_decay);

    if (setNoeta) SetNOETA(forbid_eta);

    if (setNoevolve) SetNOEVOL(forbid_evolve);

    if (setNohadron) SetNOHADR(forbid_hadron);

    if (setNopi0) SetNOPI0(forbid_pi0);
	
    if (setNsigma) SetNSIGMA(generate_sigma);
    
    if (setP) {
	for (Int_t i = 0; i < 3; i++) {
	    SetPMIN(p_limits[2 * i], i);
	    SetPMAX(p_limits[2 * i + 1], i);
	}	
    }

    if (setPhi) {
	for (Int_t i = 0; i < 3; i++) {
	    SetPHIMIN(phi_limits[2 * i], i);
	    SetPHIMAX(phi_limits[2 * i + 1], i);
	}	
    }
    
    if (setPt) {
	for (Int_t i = 0; i < 3; i++) {
	    SetPTMIN(pt_limits[2 * i], i);
	    SetPTMAX(pt_limits[2 * i + 1], i);
	}	
    }

    if (setTheta) {
	for (Int_t i = 0; i < 3; i++) {
	    SetTHMIN(theta_limits[2 * i], i);
	    SetTHMAX(theta_limits[2 * i + 1], i);
	}	
    }

    if (setX) {
	for (Int_t i = 0; i < 3; i++) {
	    SetXJMIN(x_limits[2 * i], i);
	    SetXJMAX(x_limits[2 * i + 1], i);
	}	
    }

    if (setY) {
	for (Int_t i = 0; i < 3; i++) {
	    SetYJMIN(y_limits[2 * i], i);
	    SetYJMAX(y_limits[2 * i + 1], i);
	}	
    }

    if (setXgen) SetAllXGEN(peter_jet_frag, 8);
}

/**************************************************************************/

void TIsajet::RestoreDefaults() 
{
// Booleans indicating which keywords should be written into the parameter file.

    setBeams = setCutjet = setFragment = setJettype1 = false;
    setJettype2 = setJettype3 = setLambda = setNodcay = false;
    setNoeta = setNoevolve = setNohadron = setNopi0 = false;
    setNsigma = setP = setPhi = setPt  = setTheta = false;
    setX = setXgen = setY = setPdf = false;
    num_Pdf = 0;

// Calling on FORTRAN for initialisation of variables

    Openfiles();
    Int_t a, b, c, d, e;
    
    a = -54;
    b = 0;
    c = 51;
    d = 53;

    Isaini(a, b, c, d);    
    e = 0;
    Isabeg(e);
}

/**************************************************************************/

Int_t TIsajet::ImportParticles(TClonesArray *particles, Option_t *option)
{
//
//  Default primary creation method. It reads the /HEPEVT/ common block which
//  has been filled by the GenerateEvent method. If the event generator does
//  not use the HEPEVT common block, this routine has to be overloaded by
//  the subclasses.
//  The function loops on the generated particles and stores them in
//  the TClonesArray pointed by the argument particles.
//  The default action is to store only the stable particles (ISTHEP = 1)
//  This can be demanded explicitly by setting the option = "Final"
//  If the option = "All", all the particles are stored.
//
  if (particles == 0) return 0;
  TClonesArray &Particles = *particles;
  Particles.Clear();
  TDatabasePDG* converter = TDatabasePDG::Instance();
  Int_t numpart = PARTCL.nptcl;
  printf("\n TIsajet: ISAJET stack contains %d particles.", numpart);
  printf("\n TIsajet: Total energy:         %f           ", PRIMAR.ecm);
  Int_t nump = 0;
  if ((!strcmp(option,"")) || (!strcmp(option,"Final"))) {
      for (Int_t i = 0; i < numpart; i++) {
	  
	  if (PARTCL.idcay[i] == 0) {  // Check whether particle is stable.
//  
//  Use the common block values for the TParticle constructor
//
	    nump++;
	    new(Particles[i]) TParticle(
		  converter->ConvertIsajetToPdg(PARTCL.ident[i]) , // PDG code
		  0 , // Status - currently a default
		  
		  -1, // Mothers and daughters - not used for stable particles
		  -1,
		  -1,
		  -1,
		  
		  PARTCL.pptcl[i][0] ,  // x, y, z and 0 momenta
		  PARTCL.pptcl[i][1] ,
		  PARTCL.pptcl[i][2] ,
		  PARTCL.pptcl[i][3] ,
		  
		  0, // Velocities - currently not used.
		  0,
		  0,
		  0);
	  }
      }
  }
  else if (!strcmp(option,"All")) {
      nump=numpart; 
      for (Int_t i = 0; i < numpart; i++) {

	  // Determine mother particle. Set to -1 if the particle originates from
	  // a parton or is a beam particle.

	  Int_t origin = PARTCL.iorig[i];
	  Int_t jet = origin / PARTCL.ipack;
	  origin = origin - (jet * PARTCL.ipack);
	  
	  if (origin < 0) origin = 0;

	  // Determine first and last decay products. Both are -1 if the particle is stable.
	  // Note this means they are set to 0, because one is subtracted after decoding;
	  // this avoid off-by-one errors relative to the FORTRAN.

	  Int_t first_Daughter = 0;
	  Int_t last_Daughter = 0;
	  
	  if (PARTCL.idcay[i] != 0) {
	      first_Daughter = PARTCL.idcay[i] / PARTCL.ipack;
	      last_Daughter = PARTCL.idcay[i] - (first_Daughter * PARTCL.ipack);
	  }	  
	  new(Particles[i]) TParticle(
	      converter->ConvertIsajetToPdg(PARTCL.ident[i]) ,
	      0,

	      origin - 1, 
	      -1,
	      first_Daughter - 1,
	      last_Daughter - 1,
	      
	      PARTCL.pptcl[i][0] ,
	      PARTCL.pptcl[i][1] ,
	      PARTCL.pptcl[i][2] ,
	      PARTCL.pptcl[i][3] ,
	      
	      0,
	      0,
	      0,
	      0);
      }
  }
  return nump;
}

/**************************************************************************/

void TIsajet::GenerateEvent() 
{
    Int_t e, ok, done;
 
//    e = 0;

//    if (online) Isabg2(e);
//    else Isabeg(e);

    e = 1;
    Isaevt(e, ok, done);
}

/**************************************************************************/

void TIsajet::SetJobtype(Char_t *val) 
{
    if ((!strcmp(val, "TWOJET")) || (!strcmp(val, "E+E-")) ||
    (!strcmp(val, "DRELLYAN")) || (!strcmp(val, "MINBIAS")) ||
    (!strcmp(val, "SUSY")) || (!strcmp(val, "WPAIR")) ||
    (!strcmp(val, "HIGGS")) || (!strcmp(val, "PHOTON")) ||
    (!strcmp(val, "TCOLOR")) || (!strcmp(val, "WHIGGS")) ||
    (!strcmp(val, "EXTRADIM")) || (!strcmp(val, "ZJJ"))) {
	    jobtype = val;
    }
    else {
	printf("Error in TIsajet::SetJobtype :\n");
	printf("Invalid reaction keyword %s.\n", val);
	printf("Valid keywords are : TWOJET, E+E-, DRELLYAN,\n");
	printf("MINBIAS, SUSY, WPAIR, HIGGS, PHOTON, TCOLOR,\n");
	printf("WHIGGS, EXTRADIM and ZJJ.\n");
    }    
}

/**************************************************************************/

void TIsajet::GetJobtype() const 
{
    printf ("Current job type is %s.\n", jobtype);
}

/**************************************************************************/

void TIsajet::SetOnline(Bool_t val) 
{
    online = val;
}

/**************************************************************************/

Bool_t TIsajet::GetOnline() const
{
    return online;
}

/**************************************************************************/

void TIsajet::SetPDF(Char_t *name, Float_t val)
{
    if (num_Pdf < 19) {
	pdfpar[num_Pdf] = name;
	pdfval[num_Pdf] = val;
	num_Pdf++;
	setPdf = true;
    }
    else {
	printf ("Error in TIsajet::SetPDF :\n");
	printf ("Maximum of twenty PDF parameters may be set.\n");
    }
}

/**************************************************************************/

// Access routines for common blocks.
// Begins DYLIM access.

/**************************************************************************/

void TIsajet::SetQMIN(Float_t val)
{
    if (val > DYLIM.qmax) {
	printf("Error in TIsajet::SetQMIN : \n");
	printf("You may not set QMIN to a value larger than QMAX = %f.\n", DYLIM.qmax);
	return;
    }

    DYLIM.qmin = val;
    if (!DYLIM.ywset) SetYWLIMS();
}

/**************************************************************************/

Float_t TIsajet::GetQMIN() const
{
    return DYLIM.qmin;
}

/**************************************************************************/

void TIsajet::SetQMAX(Float_t val)
{
    if (val < DYLIM.qmin) {
	printf("Error in TIsajet::SetQMAX : \n");
	printf("You may not set QMAX to a value less than QMIN = %f.\n", DYLIM.qmin);
	return;
    }

    DYLIM.qmax = val;
}

/**************************************************************************/

Float_t TIsajet::GetQMAX() const
{
    return DYLIM.qmax;
}

/**************************************************************************/

void TIsajet::SetQTMIN(Float_t val)
{
    if (val > DYLIM.qtmax) {
	printf("Error in TIsajet::SetQTMIN : \n");
	printf("You may not set QTMIN to a value larger than QTMAX = %f.\n", DYLIM.qtmax);
	return;
    }
    DYLIM.qtmin = val;
    if (!DYLIM.ywset) SetYWLIMS();
}

/**************************************************************************/

Float_t TIsajet::GetQTMIN() const
{
    return DYLIM.qtmin;
}

/**************************************************************************/

void TIsajet::SetQTMAX(Float_t val)
{
    if (val < DYLIM.qtmin) {
	printf("Error in TIsajet::SetQTMAX : \n");
	printf("You may not set QTMAX to a value less than QTMIN = %f.\n", DYLIM.qtmin);
	return;
    }

    DYLIM.qtmax = val;
    if (!DYLIM.ywset) SetYWLIMS();
}

/**************************************************************************/

Float_t TIsajet::GetQTMAX() const
{
    return DYLIM.qtmax;
}

/**************************************************************************/

void TIsajet::SetYWMIN(Float_t val)
{
    if (val > DYLIM.ywmax) {
	printf("Error in TIsajet::SetYWMIN : \n");
	printf("You may not set YWMIN to a value larger than YWMAX = %f.\n", DYLIM.ywmax);
	return;
    }

    if (DYLIM.thwset) {
	printf("Error in TIsajet::SetYWMIN :\n");
	printf("May not set both theta and y limits. Use SetTHWLIMS, then set YWMIN.\n");
    }
    else {
	DYLIM.ywset = true;
	DYLIM.ywmin = val;
    }
}

/**************************************************************************/

Float_t TIsajet::GetYWMIN() const
{
    return DYLIM.ywmin;
}

/**************************************************************************/

void TIsajet::SetYWMAX(Float_t val)
{
    if (val < DYLIM.ywmin) {
	printf("Error in TIsajet::SetYWMAX : \n");
	printf("You may not set YWMAX to a value less than YWMIN = %f.\n", DYLIM.ywmin);
	return;
    }

    if (DYLIM.thwset) {
	printf("Error in TIsajet::SetYWMAX :\n");
	printf("May not set both theta and y limits. Use SetTHWLIMS, then set YWMAX.\n");
    }
    else {
	DYLIM.ywset = true;
	DYLIM.ywmax = val;
    }
}

/**************************************************************************/

Float_t TIsajet::GetYWMAX() const
{
    return DYLIM.ywmax;
}

/**************************************************************************/

void TIsajet::SetYWLIMS() 
{
    Float_t rot = sqrt(DYLIM.qmin * DYLIM.qmin + DYLIM.qtmin * DYLIM.qtmin);
    DYLIM.ywmax = acosh(PRIMAR.halfe / rot);
    DYLIM.ywmin = -DYLIM.ywmax;
    DYLIM.ywset = false;
}

/**************************************************************************/

void TIsajet::SetXWMIN(Float_t val)
{
    if (val > DYLIM.xwmax) {
	printf("Error in TIsajet::SetXWMIN : \n");
	printf("You may not set XWMIN to a value larger than XWMAX = %f.\n", DYLIM.xwmax);
	return;
    }
    DYLIM.xwmin = val;
}

/**************************************************************************/

Float_t TIsajet::GetXWMIN() const
{
    return DYLIM.xwmin;
}

/**************************************************************************/

void TIsajet::SetXWMAX(Float_t val)
{
    if (val < DYLIM.xwmin) {
	printf("Error in TIsajet::SetXWMAX : \n");
	printf("You may not set XWMAX to a value less than XWMIN = %f.\n", DYLIM.xwmin);
	return;
    }

    DYLIM.xwmax = val;
}

/**************************************************************************/

Float_t TIsajet::GetXWMAX() const
{
    return DYLIM.xwmax;
}

/**************************************************************************/

void TIsajet::SetTHWMIN(Float_t val)
{
    if (val > DYLIM.thwmax) {
	printf("Error in TIsajet::SetTHWMIN : \n");
	printf("You may not set THWMIN to a value larger than THWMAX = %f.\n", DYLIM.thwmax);
	return;
    }

    if (DYLIM.ywset) {
	printf("Error in TIsajet::SetTHWMIN :\n");
	printf("May not set both theta and y limits. Use SetYWLIMS, then set THWMIN.\n");
    }
    else {
	DYLIM.thwset = true;
        DYLIM.thwmin = val;
    }
}

/**************************************************************************/

Float_t TIsajet::GetTHWMIN() const
{
    return DYLIM.thwmin;
}

/**************************************************************************/

void TIsajet::SetTHWMAX(Float_t val)
{
    if (val < DYLIM.thwmin) {
	printf("Error in TIsajet::SetTHWMAX : \n");
	printf("You may not set THWMAX to a value less than THWMIN = %f.\n", DYLIM.thwmin);
	return;
    }

    if (DYLIM.ywset) {
	printf("Error in TIsajet::SetTHWMAX :\n");
	printf("May not set both theta and y limits. Use SetYWLIMS, then set THWMAX.\n");
    }
    else {
	DYLIM.thwset = true;
	DYLIM.thwmax = val;
    }
}

/**************************************************************************/

Float_t TIsajet::GetTHWMAX() const
{
    return DYLIM.thwmax;
}

/**************************************************************************/

void TIsajet::SetTHWLIMS() 
{
    DYLIM.thwmin = 0;
    DYLIM.thwmax = TMath::Pi();
    DYLIM.thwset = false;
}

/**************************************************************************/

void TIsajet::SetPHWMIN(Float_t val)
{
    if (val > DYLIM.phwmax) {
	printf("Error in TIsajet::SetPHWMIN : \n");
	printf("You may not set PHWMIN to a value larger than PHWMAX = %f.\n", DYLIM.phwmax);
	return;
    }
    DYLIM.phwmin = val;
}

/**************************************************************************/

Float_t TIsajet::GetPHWMIN() const
{
    return DYLIM.phwmin;
}

/**************************************************************************/

void TIsajet::SetPHWMAX(Float_t val)
{
    if (val < DYLIM.phwmin) {
	printf("Error in TIsajet::SetPHWMAX : \n");
	printf("You may not set PHWMAX to a value less than PHWMIN = %f.\n", DYLIM.phwmin);
	return;
    }

    DYLIM.phwmax = val;
}

/**************************************************************************/

Float_t TIsajet::GetPHWMAX() const
{
    return DYLIM.phwmax;
}

/**************************************************************************/

Bool_t TIsajet::GetSETLMQ(Int_t index) const
{
    Int_t length = (sizeof DYLIM.setlmq / sizeof DYLIM.setlmq[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetSETLMQ : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return DYLIM.setlmq[index];
}

/**************************************************************************/

// End of DYLIM access. 
// Begins EEPAR access.

/**************************************************************************/

void TIsajet::SetPLEP(Float_t val)
{
    EEPAR.plep = val;
}

/**************************************************************************/

Float_t TIsajet::GetPLEP() const
{
    return EEPAR.plep;
}

/**************************************************************************/

void TIsajet::SetPLEM(Float_t val)
{
    EEPAR.plem = val;
}

/**************************************************************************/

Float_t TIsajet::GetPLEM() const
{
    return EEPAR.plem;
}

/**************************************************************************/

void TIsajet::SetRSHMIN(Float_t val)
{
    if (val > EEPAR.rshmax) {
	printf("Error in TIsajet::SetRSHMIN : \n");
	printf("You may not set RSHMIN to a value larger than RSHMAX = %f.\n", EEPAR.rshmax);
	return;
    }
    EEPAR.rshmin = val;
}

/**************************************************************************/

Float_t TIsajet::GetRSHMIN() const
{
    return EEPAR.rshmin;
}

/**************************************************************************/

void TIsajet::SetRSHMAX(Float_t val)
{
    if (val < EEPAR.rshmin) {
	printf("Error in TIsajet::SetRSHMAX : \n");
	printf("You may not set RSHMAX to a value less than RSHMIN = %f.\n", EEPAR.rshmin);
	return;
    }

    EEPAR.rshmax = val;
}

/**************************************************************************/

Float_t TIsajet::GetRSHMAX() const
{
    return EEPAR.rshmax;
}

/**************************************************************************/

void TIsajet::SetUPSLON(Float_t val)
{
    EEPAR.upslon = val;
}

/**************************************************************************/

Float_t TIsajet::GetUPSLON() const
{
    return EEPAR.upslon;
}

/**************************************************************************/

void TIsajet::SetSIGZ(Float_t val)
{
    EEPAR.sigz = val;
}

/**************************************************************************/

Float_t TIsajet::GetSIGZ() const
{
    return EEPAR.sigz;
}

/**************************************************************************/

Bool_t TIsajet::GetIBREM() const
{
    return EEPAR.ibrem;
}

/**************************************************************************/

Bool_t TIsajet::GetIBEAM() const
{
    return EEPAR.ibeam;
}

/**************************************************************************/

Float_t TIsajet::GetSGMXEE() const
{
    return EEPAR.sgmxee;
}

/**************************************************************************/

// End of EEPAR access.
// Begins FORCE access.

/**************************************************************************/

Int_t TIsajet::GetNFORCE() const
{
    return FORCE.nforce;
}

/**************************************************************************/

void TIsajet::SetIFORCE(const Int_t val[], Int_t arraySize, Bool_t anti)
{
    if (GetNFORCE() >= FORCE.mxforc - anti) {
	printf ("ERROR in TIsajet::SetIFORCE :\n");
	printf ("Cannot have more than %d forced decays.\n", FORCE.mxforc );
	return;
    }

    if ((arraySize > 6) || (arraySize < 2)) {
	printf ("Error in TIsajet::SetIFORCE : \n");
	printf ("Invalid array size %d; must be 2-6.\n", arraySize);
	return;
    }
    
    for (Int_t i = 0; i < FORCE.nforce; i++) {
	if (FORCE.iforce[i] == val[0]) {
	    printf ("Error in TIsajet::SetIFORCE : \n");
	    printf ("Particle %d has already been forced, index %d.\n", val[0], i);
	    return;
	}
    }
    

    FORCE.iforce[FORCE.nforce] = val[0];
    for (Int_t i = 1; i < arraySize; i++) {
	FORCE.mforce[FORCE.nforce][i-1] = val[i];
    }

    FORCE.nforce++;

    printf ("Decay channel %d -> ", val[0]);
    for (Int_t i = 1; i < arraySize; i++) {
	printf ("%d, ", val[i]);
    }
    printf ("set. \n");

    if (anti) {
	Int_t* antivals = new Int_t[arraySize];
	for (Int_t i = 0; i < arraySize; i++){
	    antivals[i] = (0 - val[i]);
	}
	SetIFORCE(antivals, arraySize, false);
    }
}

/**************************************************************************/

void TIsajet::UnForce(Int_t index, Bool_t anti)
{
    if (FORCE.nforce == 0) {
	printf ("Error in TIsajet::UnForce : \n");
	printf ("No decays have been forced.\n");
	return;
    }

    if ((index < 0) || (index >= FORCE.nforce)) {
	printf ("Error in TIsajet::UnForce : \n");
	printf ("Invalid decay index %d; range is 0-%d.\n", index, FORCE.nforce-1);
	return;
    }
    
    Int_t particle_ID = FORCE.iforce[index];

    for (Int_t i = index; i < FORCE.mxforc - 1; i++) {
	FORCE.iforce[i] = FORCE.iforce[i+1];
	for (Int_t j = 0; j < 5; j++) {
	    FORCE.mforce[i][j] = FORCE.mforce[i+1][j];
	}
    }
    FORCE.iforce[FORCE.mxforc - 1] = 0;
    for (Int_t j = 0; j < 5; j++) {
	FORCE.mforce[FORCE.mxforc - 1][j] = 0;
    }
    
    FORCE.nforce--;

    printf ("Decay of %d unforced.\n", particle_ID);

    if (anti) UnForceID(-particle_ID, false);
}

/**************************************************************************/

void TIsajet::UnForceID(Int_t particle_ID, Bool_t anti) 
{
    if (FORCE.nforce == 0) {
	printf ("Error in TIsajet::UnForceID : \n");
	printf ("No decays have been forced.\n");
	return;
    }

    for (Int_t i = 0; i < FORCE.nforce; i++) {
	if (FORCE.iforce[i] == particle_ID) {
	    UnForce(i, anti);
	    return;
	}
    }

    printf ("Error in TIsajet::UnForceID : \n");
    printf ("Cannot find particle %d.\n", particle_ID);
}

/**************************************************************************/

Int_t* TIsajet::GetIFORCE(Int_t index) const
{
    if (FORCE.nforce == 0) {
	printf ("Error in TIsajet::GetIFORCE : \n");
	printf ("No decays have been forced.\n");
	return 0;
    }
    
    if ((index < 0) || (index >= FORCE.nforce)) {
	printf ("Error in TIsajet::GetIFORCE : \n");
	printf ("Invalid decay index %d; range is 0-%d.\n", index, FORCE.nforce-1);
	return 0;
    }

    Int_t*  decay = new Int_t[6];
    decay[0] = FORCE.iforce[index];

    for (Int_t i = 1; i < 6; i++) {
	decay[i] = FORCE.mforce[index][i-1];
    }
    
    return decay;
}

/**************************************************************************/

Int_t TIsajet::GetMEFORC(Int_t index) const
{
    if (FORCE.nforce == 0) {
	printf ("Error in TIsajet::GetMEFORCE : \n");
	printf ("No decays have been forced.\n");
	return 0;
    }

    if ((index < 0) || (index >= FORCE.nforce)) {
	printf ("Error in TIsajet::GetMEFORC : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, FORCE.nforce-1);
	return 0;
    }

    return FORCE.meforc[index];
}

/**************************************************************************/

// End of FORCE access.
// Begins FRGPAR access.

/**************************************************************************/

void TIsajet::SetFRPAR(Float_t val, Int_t index)
{
    Int_t length = (sizeof FRGPAR.frpar / sizeof FRGPAR.frpar[0]);
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetFRPAR : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }
    
    *FRGPAR.frpar[index] = val;
    setFragment = true;
}

/**************************************************************************/

void TIsajet::SetAllFRPAR(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof FRGPAR.frpar / sizeof FRGPAR.frpar[0]);
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllFRPAR : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetFRPAR(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetFRPAR(Int_t index) const 
{
    Int_t length = (sizeof FRGPAR.frpar / sizeof FRGPAR.frpar[0]);
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetFRPAR : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return *FRGPAR.frpar[index];
}

/**************************************************************************/

void TIsajet::SetPUD(Float_t val) 
{
    SetFRPAR(val, 0);
}

/**************************************************************************/

Float_t TIsajet::GetPUD() const 
{
    return GetFRPAR(0);
}

/**************************************************************************/

void TIsajet::SetPBARY(Float_t val) 
{
    SetFRPAR(val, 1);    
}

/**************************************************************************/

Float_t TIsajet::GetPBARY() const 
{
    return GetFRPAR(1);
}

/**************************************************************************/

void TIsajet::SetSIGQT(Float_t val) 
{
    SetFRPAR(val, 2);    
}

/**************************************************************************/

Float_t TIsajet::GetSIGQT() const 
{
    return GetFRPAR(2);
}

/**************************************************************************/

void TIsajet::SetPEND(Float_t val) 
{
    SetFRPAR(val, 3);    
}

/**************************************************************************/

Float_t TIsajet::GetPEND() const 
{
    return GetFRPAR(3);
}

/**************************************************************************/

void TIsajet::SetXGEN(Float_t val, Int_t index)
{
    Int_t length = (sizeof FRGPAR.xgen / sizeof FRGPAR.xgen[0]);
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetXGEN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }
    SetFRPAR(val, index + 4);
    setXgen = true;
}

/**************************************************************************/

void TIsajet::SetAllXGEN(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof FRGPAR.xgen / sizeof FRGPAR.xgen[0]);
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllXGEN : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetXGEN(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetXGEN(Int_t index) const 
{
    Int_t length = (sizeof FRGPAR.xgen / sizeof FRGPAR.xgen[0]);
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetXGEN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return GetFRPAR(index + 4);
}

/**************************************************************************/

void TIsajet::SetPSPIN1(Float_t val, Int_t index)
{
    Int_t length = (sizeof FRGPAR.pspin1 / sizeof FRGPAR.pspin1[0]);
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetPSPIN1 : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }
    
    SetFRPAR(val, index + 12);
}

/**************************************************************************/

void TIsajet::SetAllPSPIN1(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof FRGPAR.pspin1 / sizeof FRGPAR.pspin1[0]);
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllPSPIN1 : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetPSPIN1(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetPSPIN1(Int_t index) const 
{
    Int_t length = (sizeof FRGPAR.xgen / sizeof FRGPAR.xgen[0]);
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetPSPIN1 : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return GetFRPAR(index + 12);
}

/**************************************************************************/

void TIsajet::SetPMIX1(Float_t val, Int_t index1, Int_t index2)
{
    Int_t col_num = (sizeof FRGPAR.pmix1[0] / sizeof FRGPAR.pmix1[0][0]);
    Int_t row_num = (sizeof FRGPAR.pmix1 / (sizeof FRGPAR.pmix1[0][0] * col_num));
    
    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::SetPMIX1 : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::SetPMIX1 : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return;
    }
    
    FRGPAR.pmix1[index1][index2] = val;
    setFragment = true;
}

/**************************************************************************/

void TIsajet::SetAllPMIX1(const Float_t val[2][3]) 
{
    for (Int_t i = 0; i < 2; i++) {
	for (Int_t j = 0; j < 3; j++) {
	    SetPMIX1(val[i][j], i, j);
	}
    }
}

/**************************************************************************/

void TIsajet::SetColumnPMIX1(const Float_t val[], Int_t col) 
{
    Int_t col_num = (sizeof FRGPAR.pmix1[0] / sizeof FRGPAR.pmix1[0][0]);
    Int_t row_num = (sizeof FRGPAR.pmix1 / (sizeof FRGPAR.pmix1[0][0] * col_num));
    
    if ((col < 0) || (col >= col_num)) {
	printf ("Error in TIsajet::SetColumnPMIX1 : \n");
	printf ("Invalid column index %d, range is 0-%d.\n", col, col_num-1);
	return;
    }
    
    for (Int_t i = 0; i < row_num; i++) {
	SetPMIX1(val[i], i, col);
    }
}

/**************************************************************************/

Float_t TIsajet::GetPMIX1(Int_t index1, Int_t index2) const 
{
    Int_t col_num = (sizeof FRGPAR.pmix1[0] / sizeof FRGPAR.pmix1[0][0]);
    Int_t row_num = (sizeof FRGPAR.pmix1 / (sizeof FRGPAR.pmix1[0][0] * col_num));

    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::GetPMIX1 : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return 0;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::GetPMIX1 : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return 0;
    }

    return FRGPAR.pmix1[index1][index2];
}

/**************************************************************************/

void TIsajet::SetPMIX2(Float_t val, Int_t index1, Int_t index2)
{
    Int_t col_num = (sizeof FRGPAR.pmix2[0] / sizeof FRGPAR.pmix2[0][0]);
    Int_t row_num = (sizeof FRGPAR.pmix2 / (sizeof FRGPAR.pmix2[0][0] * col_num));
    
    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::SetPMIX2 : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::SetPMIX2 : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return;
    }
    
    FRGPAR.pmix2[index1][index2] = val;
    setFragment = true;
}

/**************************************************************************/

void TIsajet::SetAllPMIX2(const Float_t val[2][3]) 
{
    for (Int_t i = 0; i < 2; i++) {
	for (Int_t j = 0; j < 3; j++) {
	    SetPMIX2(val[i][j], i, j);
	}
    }
}

/**************************************************************************/

void TIsajet::SetColumnPMIX2(const Float_t val[], Int_t col) 
{
    Int_t col_num = (sizeof FRGPAR.pmix2[0] / sizeof FRGPAR.pmix2[0][0]);
    Int_t row_num = (sizeof FRGPAR.pmix2 / (sizeof FRGPAR.pmix2[0][0] * col_num));
    
    if ((col < 0) || (col >= col_num)) {
	printf ("Error in TIsajet::SetColumnPMIX2 : \n");
	printf ("Invalid column index %d, range is 0-%d.\n", col, col_num-1);
	return;
    }
    
    for (Int_t i = 0; i < row_num; i++) {
	SetPMIX2(val[i], i, col);
    }
}

/**************************************************************************/

Float_t TIsajet::GetPMIX2(Int_t index1, Int_t index2) const 
{
    Int_t col_num = (sizeof FRGPAR.pmix2[0] / sizeof FRGPAR.pmix2[0][0]);
    Int_t row_num = (sizeof FRGPAR.pmix2 / (sizeof FRGPAR.pmix2[0][0] * col_num));

    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::GetPMIX2 : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return 0;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::GetPMIX2 : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return 0;
    }

    return FRGPAR.pmix2[index1][index2];
}

/**************************************************************************/

void TIsajet::SetPMIXX1(Float_t val, Int_t index)
{
    Int_t length = (sizeof FRGPAR.pmixx1 / sizeof FRGPAR.pmixx1[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetPMIXX1 : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }
    
    *FRGPAR.pmixx1[index] = val;
    setFragment = true;
}

/**************************************************************************/

void TIsajet::SetAllPMIXX1(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof FRGPAR.pmixx1 / sizeof FRGPAR.pmixx1[0]);    
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllPMIXX1 : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetPMIXX1(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetPMIXX1(Int_t index) const 
{
    Int_t length = (sizeof FRGPAR.pmixx1 / sizeof FRGPAR.pmixx1[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetPMIXX1 : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return *FRGPAR.pmixx1[index];
}

/**************************************************************************/

void TIsajet::SetPMIXX2(Float_t val, Int_t index)
{
    Int_t length = (sizeof FRGPAR.pmixx2 / sizeof FRGPAR.pmixx2[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetPMIXX2 : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }
    
    *FRGPAR.pmixx2[index] = val;
    setFragment = true;
}

/**************************************************************************/

void TIsajet::SetAllPMIXX2(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof FRGPAR.pmixx2 / sizeof FRGPAR.pmixx2[0]);    
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllPMIXX2 : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetPMIXX2(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetPMIXX2(Int_t index) const 
{
    Int_t length = (sizeof FRGPAR.pmixx2 / sizeof FRGPAR.pmixx2[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetPMIXX2 : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return *FRGPAR.pmixx2[index];
}

/**************************************************************************/

void TIsajet::SetXGENSS(Float_t val, Int_t index)
{
    Int_t length = (sizeof FRGPAR.xgenss / sizeof FRGPAR.xgenss[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetXGENSS : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }
    
     FRGPAR.xgenss[index] = val;
}

/**************************************************************************/

void TIsajet::SetAllXGENSS(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof FRGPAR.xgenss / sizeof FRGPAR.xgenss[0]);    
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllXGENSS : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetXGENSS(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetXGENSS(Int_t index) const 
{
    Int_t length = (sizeof FRGPAR.xgenss / sizeof FRGPAR.xgenss[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetXGENSS : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return FRGPAR.xgenss[index];
}

/**************************************************************************/

// End of FRGPAR access.
// Begins HCON access.

/**************************************************************************/

Float_t TIsajet::GetANWWWW(Int_t index1, Int_t index2, Int_t index3) const
{
    Int_t elem_Size = sizeof HCON.anwwww[0][0][0];
    Int_t thd_Dim_Length = (sizeof HCON.anwwww[0][0] / elem_Size);
    Int_t sec_Dim_Length = (sizeof HCON.anwwww[0] / (elem_Size * thd_Dim_Length));
    Int_t fst_Dim_Length = (sizeof HCON.anwwww / (elem_Size * thd_Dim_Length * sec_Dim_Length));
    
    if ((index1 < 0) || (index1 >= fst_Dim_Length)) {
	printf ("Error in TIsajet::GetANWWWW : \n");
	printf ("Invalid first index %d; range is 0-%d.\n", index1, fst_Dim_Length-1);
	return 0;
    }
    
    if ((index2 < 0) || (index2 >= sec_Dim_Length)) {
	printf ("Error in TIsajet::GetANWWWW : \n");
	printf ("Invalid second index %d; range is 0-%d.\n", index2, sec_Dim_Length-1);
	return 0;
    }

    if ((index3 < 0) || (index3 >= thd_Dim_Length)) {
	printf ("Error in TIsajet::GetANWWWW : \n");
	printf ("Invalid third index %d; range is 0-%d.\n", index3, thd_Dim_Length-1);
	return 0;
    }

    return HCON.anwwww[index1][index2][index3];
}

/**************************************************************************/

Float_t TIsajet::GetADWWWW(Int_t index1, Int_t index2) const
{
    Int_t elem_Size = sizeof HCON.adwwww[0][0];
    Int_t sec_Dim_Length = (sizeof HCON.adwwww[0] / elem_Size);
    Int_t fst_Dim_Length = (sizeof HCON.adwwww / (elem_Size * sec_Dim_Length));
    
    if ((index1 < 0) || (index1 >= fst_Dim_Length)) {
	printf ("Error in TIsajet::GetADWWWW : \n");
	printf ("Invalid first index %d; range is 0-%d.\n", index1, fst_Dim_Length-1);
	return 0;
    }
    
    if ((index2 < 0) || (index2 >= sec_Dim_Length)) {
	printf ("Error in TIsajet::GetADWWWW : \n");
	printf ("Invalid second index %d; range is 0-%d.\n", index2, sec_Dim_Length-1);
	return 0;
    }

    return HCON.adwwww[index1][index2];
}

/**************************************************************************/

Float_t TIsajet::GetAIWWWW(Int_t index) const 
{
    Int_t length = (sizeof HCON.aiwwww / sizeof HCON.aiwwww[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetAIWWWW : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return HCON.aiwwww[index];
}

/**************************************************************************/

Float_t TIsajet::GetHMASS() const 
{
    return HCON.hmass;
}

/**************************************************************************/

Float_t TIsajet::GetHGAM() const 
{
    return HCON.hgam;
}

/**************************************************************************/

Float_t TIsajet::GetHGAMS(Int_t index) const 
{
    Int_t length = (sizeof HCON.hgams / sizeof HCON.hgams[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetHGAMS : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return HCON.hgams[index];
}

/**************************************************************************/

Float_t TIsajet::GetETAHGG() const 
{
    return HCON.etahgg;
}

/**************************************************************************/

Int_t TIsajet::GetMATCHH(Int_t index) const 
{
    Int_t length = (sizeof HCON.matchh / sizeof HCON.matchh[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetMATCHH : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return HCON.matchh[index];
}

/**************************************************************************/

Float_t TIsajet::GetZSTARS(Int_t index1, Int_t index2) const 
{
    Int_t elem_Size = sizeof HCON.zstars[0][0];
    Int_t sec_Dim_Length = (sizeof HCON.zstars[0] / elem_Size);
    Int_t fst_Dim_Length = (sizeof HCON.zstars / (elem_Size * sec_Dim_Length));
    
    if ((index1 < 0) || (index1 >= fst_Dim_Length)) {
	printf ("Error in TIsajet::GetZSTARS : \n");
	printf ("Invalid first index %d; range is 0-%d.\n", index1, fst_Dim_Length-1);
	return 0;
    }
    
    if ((index2 < 0) || (index2 >= sec_Dim_Length)) {
	printf ("Error in TIsajet::GetZSTARS : \n");
	printf ("Invalid second index %d; range is 0-%d.\n", index2, sec_Dim_Length-1);
	return 0;
    }

    return HCON.zstars[index1][index2];
}

/**************************************************************************/

void TIsajet::SetIHTYPE(Int_t val) 
{
    if ((val < 82) || (val > 84)) {
	printf ("Error in TIsajet::SetIHTYPE : \n");
	printf ("Invalid input value %d. Possible values are 82, 83, 84.\n", val);
	return;
    }
    
    HCON.ihtype = val;
}

/**************************************************************************/

void TIsajet::SetIHTYPE(Char_t val[])
{
    if (!strcmp("HL0", val)) {
	HCON.ihtype = 82;
    }
    else if (!strcmp("HH0", val)) {
	HCON.ihtype = 83;
    }
    else if (!strcmp("HA0", val)){
	HCON.ihtype = 84;
    }
    else {
	printf ("Error in TIsajet::SetIHTYPE : \n");
	printf ("Invalid input string %s. Possible strings are HL0, HH0, HA0.\n", val);
    }
}

/**************************************************************************/

Int_t TIsajet::GetIHTYPE() const 
{
    return HCON.ihtype;
}

/**************************************************************************/

Float_t TIsajet::GetHGAMSS(Int_t index1, Int_t index2) const 
{
    Int_t elem_Size = sizeof HCON.hgamss[0][0];
    Int_t sec_Dim_Length = (sizeof HCON.hgamss[0] / elem_Size);
    Int_t fst_Dim_Length = (sizeof HCON.hgamss / (elem_Size * sec_Dim_Length));
    
    if ((index1 < 0) || (index1 >= fst_Dim_Length)) {
	printf ("Error in TIsajet::GetHGAMSS : \n");
	printf ("Invalid first index %d; range is 0-%d.\n", index1, fst_Dim_Length-1);
	return 0;
    }
    
    if ((index2 < 0) || (index2 >= sec_Dim_Length)) {
	printf ("Error in TIsajet::GetHGAMSS : \n");
	printf ("Invalid second index %d; range is 0-%d.\n", index2, sec_Dim_Length-1);
	return 0;
    }

    return HCON.hgamss[index1][index2];
}

/**************************************************************************/

// End of HCON access
// Begins JETLIM access

/**************************************************************************/

void TIsajet::SetPMIN(Float_t val, Int_t index)
{
    Int_t length = (sizeof JETLIM.pmin / sizeof JETLIM.pmin[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetPMIN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }

    if (val > JETLIM.pmax[index]) {
	printf("Error in TIsajet::SetPMIN : \n");
	printf("You may not set PMIN to a value larger than PMAX = %f.\n", JETLIM.pmax[index]);
	return;
    }

    JETLIM.pmin[index] = val;
    setP = true;
}

/**************************************************************************/

void TIsajet::SetAllPMIN(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof JETLIM.pmin / sizeof JETLIM.pmin[0]);    
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllPMIN : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetPMIN(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetPMIN(Int_t index) const 
{
    Int_t length = (sizeof JETLIM.pmin / sizeof JETLIM.pmin[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetPMIN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETLIM.pmin[index];
}

/**************************************************************************/

void TIsajet::SetPMAX(Float_t val, Int_t index)
{
    Int_t length = (sizeof JETLIM.pmax / sizeof JETLIM.pmax[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetPMAX : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }
    
    if (val < JETLIM.pmin[index]) {
	printf("Error in TIsajet::SetPMAX : \n");
	printf("You may not set PMAX to a value larger than PMIN = %f.\n", JETLIM.pmin[index]);
	return;
    }

    JETLIM.pmax[index] = val;
    setP = true;
}

/**************************************************************************/

void TIsajet::SetAllPMAX(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof JETLIM.pmax / sizeof JETLIM.pmax[0]);    
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllPMAX : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetPMAX(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetPMAX(Int_t index) const 
{
    Int_t length = (sizeof JETLIM.pmax / sizeof JETLIM.pmax[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetPMAX : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETLIM.pmax[index];
}

/**************************************************************************/

void TIsajet::SetPTMIN(Float_t val, Int_t index)
{
    Int_t length = (sizeof JETLIM.ptmin / sizeof JETLIM.ptmin[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetPTMIN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }
/* andreas 7/8/2001
    if (val > JETLIM.ptmax[index]) {
	printf("Error in TIsajet::SetPTMIN : \n");
	printf("You may not set PTMIN to a value larger than PTMAX = %f.\n", JETLIM.ptmax[index]);
	return;
    }
*/    
     JETLIM.ptmin[index] = val;
//     if (!setY) SetYJLIMS();
//     if (!setTheta) SetTHLIMS();
     setPt = true;
}

/**************************************************************************/

void TIsajet::SetAllPTMIN(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof JETLIM.ptmin / sizeof JETLIM.ptmin[0]);    
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllPTMIN : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetPTMIN(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetPTMIN(Int_t index) const 
{
    Int_t length = (sizeof JETLIM.ptmin / sizeof JETLIM.ptmin[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetPTMIN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETLIM.ptmin[index];
}

/**************************************************************************/

void TIsajet::SetPTMAX(Float_t val, Int_t index)
{
    Int_t length = (sizeof JETLIM.ptmax / sizeof JETLIM.ptmax[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetPTMAX : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }

    if (val < JETLIM.ptmin[index]) {
	printf("Error in TIsajet::SetPTMAX : \n");
	printf("You may not set PTMAX to a value larger than PTMIN = %f.\n", JETLIM.ptmin[index]);
	return;
    }

    JETLIM.ptmax[index] = val;
    setPt = true;
}

/**************************************************************************/

void TIsajet::SetAllPTMAX(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof JETLIM.ptmax / sizeof JETLIM.ptmax[0]);    
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllPTMAX : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetPTMAX(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetPTMAX(Int_t index) const 
{
    Int_t length = (sizeof JETLIM.ptmax / sizeof JETLIM.ptmax[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetPTMAX : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETLIM.ptmax[index];
}

/**************************************************************************/

void TIsajet::SetYJMIN(Float_t val, Int_t index)
{
    Int_t length = (sizeof JETLIM.yjmin / sizeof JETLIM.yjmin[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetYJMIN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }

    if (val > JETLIM.yjmax[index]) {
	printf("Error in TIsajet::SetYJMIN : \n");
	printf("You may not set YJMIN to a value larger than YJMAX = %f.\n", JETLIM.yjmax[index]);
	return;
    }

    if (setTheta) {
	printf("Error in TIsajet::SetYJMIN :\n");
	printf("May not set both theta and y limits. Use SetTHLIMS, then set YJMIN.\n");
	return;
    }

    setY =  true;
    JETLIM.yjmin[index] = val;
}

/**************************************************************************/

void TIsajet::SetAllYJMIN(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof JETLIM.yjmin / sizeof JETLIM.yjmin[0]);    
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllYJMIN : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetYJMIN(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetYJMIN(Int_t index) const 
{
    Int_t length = (sizeof JETLIM.yjmin / sizeof JETLIM.yjmin[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetYJMIN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETLIM.yjmin[index];
}

/**************************************************************************/

void TIsajet::SetYJMAX(Float_t val, Int_t index)
{
    Int_t length = (sizeof JETLIM.yjmax / sizeof JETLIM.yjmax[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetYJMAX : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }

    if (val < JETLIM.yjmin[index]) {
	printf("Error in TIsajet::SetYJMAX : \n");
	printf("You may not set YJMAX to a value larger than YJMIN = %f.\n", JETLIM.yjmin[index]);
	return;
    }

    if (setTheta) {
	printf("Error in TIsajet::SetYJMAX :\n");
	printf("May not set both theta and y limits. Use SetTHLIMS, then set YJMAX.\n");
	return;
    }

    setY = true;
    JETLIM.yjmax[index] = val;
}

/**************************************************************************/

void TIsajet::SetAllYJMAX(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof JETLIM.yjmax / sizeof JETLIM.yjmax[0]);    
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllYJMAX : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetYJMAX(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetYJMAX(Int_t index) const 
{
    Int_t length = (sizeof JETLIM.yjmax / sizeof JETLIM.yjmax[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetYJMAX : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETLIM.yjmax[index];
}

/**************************************************************************/

void TIsajet::SetYJLIMS() 
{
    for (Int_t i = 0; i < JETLIM.mxlim; i++) {
	JETLIM.yjmax[i] = acosh(PRIMAR.halfe / JETLIM.ptmin[i]);
	JETLIM.yjmax[i] = -JETLIM.yjmin[i];
    }
    setY = false;
}

/**************************************************************************/

void TIsajet::SetPHIMIN(Float_t val, Int_t index)
{
    Int_t length = (sizeof JETLIM.phimin / sizeof JETLIM.phimin[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetPHIMIN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }

    if (val > JETLIM.phimax[index]) {
	printf("Error in TIsajet::SetPHIMIN : \n");
	printf("You may not set PHIMIN to a value larger than PHIMAX = %f.\n", JETLIM.phimax[index]);
	return;
    }

    JETLIM.phimin[index] = val;
    setPhi = true;
}

/**************************************************************************/

void TIsajet::SetAllPHIMIN(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof JETLIM.phimin / sizeof JETLIM.phimin[0]);    
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllPHIMIN : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetPHIMIN(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetPHIMIN(Int_t index) const 
{
    Int_t length = (sizeof JETLIM.phimin / sizeof JETLIM.phimin[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetPHIMIN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETLIM.phimin[index];
}

/**************************************************************************/

void TIsajet::SetPHIMAX(Float_t val, Int_t index)
{
    Int_t length = (sizeof JETLIM.phimax / sizeof JETLIM.phimax[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetPHIMAX : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }

    if (val < JETLIM.phimin[index]) {
	printf("Error in TIsajet::SetPHIMAX : \n");
	printf("You may not set PHIMAX to a value larger than PHIMIN = %f.\n", JETLIM.phimin[index]);
	return;
    }

    JETLIM.phimax[index] = val;
    setPhi = true;
}

/**************************************************************************/

void TIsajet::SetAllPHIMAX(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof JETLIM.phimax / sizeof JETLIM.phimax[0]);    
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllPHIMAX : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetPHIMAX(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetPHIMAX(Int_t index) const 
{
    Int_t length = (sizeof JETLIM.phimax / sizeof JETLIM.phimax[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetPHIMAX : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETLIM.phimax[index];
}

/**************************************************************************/

void TIsajet::SetXJMIN(Float_t val, Int_t index)
{
    Int_t length = (sizeof JETLIM.xjmin / sizeof JETLIM.xjmin[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetXJMIN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }
    if (val > JETLIM.xjmax[index]) {
	printf("Error in TIsajet::SetXJMIN : \n");
	printf("You may not set XJMIN to a value larger than XJMAX = %f.\n", JETLIM.xjmax[index]);
	return;
    }

    JETLIM.xjmin[index] = val;
    setX = true;
}

/**************************************************************************/

void TIsajet::SetAllXJMIN(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof JETLIM.xjmin / sizeof JETLIM.xjmin[0]);    
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllXJMIN : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetXJMIN(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetXJMIN(Int_t index) const 
{
    Int_t length = (sizeof JETLIM.xjmin / sizeof JETLIM.xjmin[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetXJMIN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETLIM.xjmin[index];
}

/**************************************************************************/

void TIsajet::SetXJMAX(Float_t val, Int_t index)
{
    Int_t length = (sizeof JETLIM.xjmax / sizeof JETLIM.xjmax[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetXJMAX : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }

    if (val < JETLIM.xjmin[index]) {
	printf("Error in TIsajet::SetXJMAX : \n");
	printf("You may not set XJMAX to a value larger than XJMIN = %f.\n", JETLIM.xjmin[index]);
	return;
    }

    JETLIM.xjmax[index] = val;
    setX = true;
}

/**************************************************************************/

void TIsajet::SetAllXJMAX(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof JETLIM.xjmax / sizeof JETLIM.xjmax[0]);    
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllXJMAX : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetXJMAX(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetXJMAX(Int_t index) const 
{
    Int_t length = (sizeof JETLIM.xjmax / sizeof JETLIM.xjmax[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetXJMAX : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETLIM.xjmax[index];
}

/**************************************************************************/

void TIsajet::SetTHMIN(Float_t val, Int_t index)
{
    Int_t length = (sizeof JETLIM.thmin / sizeof JETLIM.thmin[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetTHMIN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }

    if (val > JETLIM.thmax[index]) {
	printf("Error in TIsajet::SetTHMIN : \n");
	printf("You may not set THMIN to a value larger than THMAX = %f.\n", JETLIM.thmax[index]);
	return;
    }

    if (setY) {
	printf("Error in TIsajet::SetTHMIN :\n");
	printf("May not set both theta and y limits. Use SetYJLIMS, then set THMIN.\n");
	return;
    }

    setTheta = true;
    JETLIM.thmin[index] = val;
    
}

/**************************************************************************/

void TIsajet::SetAllTHMIN(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof JETLIM.thmin / sizeof JETLIM.thmin[0]);    
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllTHMIN : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetTHMIN(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetTHMIN(Int_t index) const 
{
    Int_t length = (sizeof JETLIM.thmin / sizeof JETLIM.thmin[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetTHMIN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETLIM.thmin[index];
}

/**************************************************************************/

void TIsajet::SetTHMAX(Float_t val, Int_t index)
{
    Int_t length = (sizeof JETLIM.thmax / sizeof JETLIM.thmax[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetTHMAX : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }

    if (val < JETLIM.thmin[index]) {
	printf("Error in TIsajet::SetTHMAX : \n");
	printf("You may not set THMAX to a value larger than THMIN = %f.\n", JETLIM.thmin[index]);
	return;
    }

    if (setY) {
	printf("Error in TIsajet::SetTHMAX :\n");
	printf("May not set both theta and y limits. Use SetYJLIMS, then set THMAX.\n");
	return;
    }

    setTheta = true;
    JETLIM.thmax[index] = val;
}

/**************************************************************************/

void TIsajet::SetAllTHMAX(const Float_t val[], Int_t arraySize) 
{
    Int_t length = (sizeof JETLIM.thmax / sizeof JETLIM.thmax[0]);    
    if (arraySize != length) {
	printf ("Error in TIsajet::SetAllTHMAX : \n");
	printf ("Array must have %d elements.\n", length);
	return;
    }
    
    for (Int_t i = 0; i < arraySize; i++) {
	SetTHMAX(val[i], i);
    }
}

/**************************************************************************/

Float_t TIsajet::GetTHMAX(Int_t index) const 
{
    Int_t length = (sizeof JETLIM.thmax / sizeof JETLIM.thmax[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetTHMAX : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETLIM.thmax[index];
}

/**************************************************************************/

void TIsajet::SetTHLIMS() 
{
    Float_t tmin;
    for (Int_t i = 0; i < JETLIM.mxlim; i++) {
	tmin = acosh(PRIMAR.halfe / JETLIM.ptmin[i]);
	JETLIM.thmin[i] = 2*atan(exp(tmin));
	JETLIM.thmax[i] = 2*atan(exp(-tmin));
    }
    setTheta = false;
}

/**************************************************************************/

Bool_t TIsajet::GetSETLMJ(Int_t index) const
{
    Int_t length = (sizeof JETLIM.setlmj / sizeof JETLIM.setlmj[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetSETLMJ : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETLIM.setlmj[index];
}

/**************************************************************************/

// Ends JETLIM access.
// Begins JETPAR access.

/**************************************************************************/

Float_t TIsajet::GetP(Int_t index) const
{
    Int_t length = (sizeof JETPAR.p / sizeof JETPAR.p[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetP : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETPAR.p[index];
}

/**************************************************************************/

Float_t TIsajet::GetPT(Int_t index) const
{
    Int_t length = (sizeof JETPAR.pt / sizeof JETPAR.pt[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetPT : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETPAR.pt[index];
}

/**************************************************************************/

Float_t TIsajet::GetYJ(Int_t index) const
{
    Int_t length = (sizeof JETPAR.yj / sizeof JETPAR.yj[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetYJ : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETPAR.yj[index];
}

/**************************************************************************/

Float_t TIsajet::GetPHI(Int_t index) const
{
    Int_t length = (sizeof JETPAR.phi / sizeof JETPAR.phi[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetPHI : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETPAR.phi[index];
}

/**************************************************************************/

Float_t TIsajet::GetXJ(Int_t index) const
{
    Int_t length = (sizeof JETPAR.xj / sizeof JETPAR.xj[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetXJ : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETPAR.xj[index];
}

/**************************************************************************/

Float_t TIsajet::GetTH(Int_t index) const
{
    Int_t length = (sizeof JETPAR.th / sizeof JETPAR.th[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetTH : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETPAR.th[index];
}

/**************************************************************************/

Float_t TIsajet::GetCTH(Int_t index) const
{
    Int_t length = (sizeof JETPAR.cth / sizeof JETPAR.cth[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetCTH : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETPAR.cth[index];
}

/**************************************************************************/

Float_t TIsajet::GetSTH(Int_t index) const
{
    Int_t length = (sizeof JETPAR.sth / sizeof JETPAR.sth[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetSTH : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETPAR.sth[index];
}

/**************************************************************************/

Int_t TIsajet::GetJETTYP(Int_t index) const
{
    Int_t length = (sizeof JETPAR.jettyp / sizeof JETPAR.jettyp[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetJETTYP : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETPAR.jettyp[index];
}

/**************************************************************************/

Float_t TIsajet::GetSHAT() const 
{
    return JETPAR.shat;
}

/**************************************************************************/

Float_t TIsajet::GetTHAT() const 
{
    return JETPAR.that;
}

/**************************************************************************/

Float_t TIsajet::GetUHAT() const 
{
    return JETPAR.uhat;
}

/**************************************************************************/

Float_t TIsajet::GetQSQ() const 
{
    return JETPAR.qsq;
}

/**************************************************************************/

Float_t TIsajet::GetX1() const 
{
    return JETPAR.x1;
}

/**************************************************************************/

Float_t TIsajet::GetX2() const 
{
    return JETPAR.x2;
}

/**************************************************************************/

Float_t TIsajet::GetPBEAM(Int_t index) const
{
    Int_t length = (sizeof JETPAR.pbeam / sizeof JETPAR.pbeam[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetPBEAM : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETPAR.pbeam[index];
}

/**************************************************************************/

Float_t TIsajet::GetQMW() const 
{
    return JETPAR.qmw;
}

/**************************************************************************/

Float_t TIsajet::GetQW() const 
{
    return JETPAR.qw;
}

/**************************************************************************/

Float_t TIsajet::GetQTW() const 
{
    return JETPAR.qtw;
}

/**************************************************************************/

Float_t TIsajet::GetYW() const 
{
    return JETPAR.yw;
}

/**************************************************************************/

Float_t TIsajet::GetXW() const 
{
    return JETPAR.xw;
}

/**************************************************************************/

Float_t TIsajet::GetTHW() const 
{
    return JETPAR.thw;
}

/**************************************************************************/

Float_t TIsajet::GetQTMW() const 
{
    return JETPAR.qtmw;
}

/**************************************************************************/

Float_t TIsajet::GetPHIW() const 
{
    return JETPAR.phiw;
}

/**************************************************************************/

Float_t TIsajet::GetSHAT1() const 
{
    return JETPAR.shat1;
}

/**************************************************************************/

Float_t TIsajet::GetTHAT1() const 
{
    return JETPAR.that1;
}

/**************************************************************************/

Float_t TIsajet::GetUHAT1() const 
{
    return JETPAR.uhat1;
}

/**************************************************************************/

void TIsajet::SetJWTYP(Int_t val) 
{
    if ((val < 1) || (val > 4) || (val == 2))
    {
	printf ("Error in TIsajet::SetJWTYP : \n");
	printf ("Invalid value  %d; range is 1, 3, and 4.\n", val);
	return;
    }
    
    JETPAR.jwtyp = val;
}

/**************************************************************************/

void TIsajet::SetJWTYP(Char_t val[]) 
{
    Int_t value;
    
    if (!strcmp(val, "GM")) value = 1;
    else if (!strcmp(val, "W+")) value = 3;    
    else if (!strcmp(val, "W-")) value = 3;    
    else if (!strcmp(val, "Z0")) value = 4;
    else 
    {
	printf ("Error in TIsajet::SetJWTYP : \n");
	printf ("Invalid value  %s; possible are GM, Z0, W+ and W-.\n", val);
	return;      
    }

    
    JETPAR.jwtyp = value;
}

/**************************************************************************/

Int_t TIsajet::GetJWTYP() const 
{
    return JETPAR.jwtyp;
}

/**************************************************************************/

Float_t TIsajet::GetALFQSQ() const 
{
    return JETPAR.alfqsq;
}

/**************************************************************************/

Float_t TIsajet::GetCTHW() const 
{
    return JETPAR.cthw;
}

/**************************************************************************/

Float_t TIsajet::GetSTHW() const 
{
    return JETPAR.sthw;
}

/**************************************************************************/

Float_t TIsajet::GetQ0W() const 
{
    return JETPAR.q0w;
}

/**************************************************************************/

Int_t TIsajet::GetINITYP(Int_t index) const
{
    Int_t length = (sizeof JETPAR.inityp / sizeof JETPAR.inityp[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetINITYP : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETPAR.inityp[index];
}

/**************************************************************************/

Int_t TIsajet::GetISIGS() const 
{
    return JETPAR.isigs;
}

/**************************************************************************/

Float_t TIsajet::GetPBEAMS(Int_t index) const
{
    Int_t length = (sizeof JETPAR.pbeams / sizeof JETPAR.pbeams[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetPBEAMS : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return JETPAR.pbeams[index];
}

/**************************************************************************/

// Ends JETPAR access.
// Begins KKGRAV access.

/**************************************************************************/

void TIsajet::SetNEXTRAD(Int_t val) 
{
    KKGRAV.nextrad = val;
}

/**************************************************************************/

Int_t TIsajet::GetNEXTRAD() const 
{
    return KKGRAV.nextrad;
}

/**************************************************************************/

void TIsajet::SetMASSD(Float_t val) 
{
    KKGRAV.massd = val;
}

/**************************************************************************/

Float_t TIsajet::GetMASSD() const 
{
    return KKGRAV.massd;
}

/**************************************************************************/

Float_t TIsajet::GetKKGSD() const 
{
    return KKGRAV.kkgsd;
}

/**************************************************************************/

Float_t TIsajet::GetSURFD() const 
{
    return KKGRAV.surfd;
}

/**************************************************************************/

void TIsajet::SetUVCUT(Bool_t val) 
{
    KKGRAV.uvcut = val;
}

/**************************************************************************/

Bool_t TIsajet::GetUVCUT() const 
{
    return KKGRAV.uvcut;
}

/**************************************************************************/

// Ends KKGRAV access.
// Begins MBGEN access.

/**************************************************************************/

Float_t TIsajet::GetPOMWT(Int_t index) const
{
    Int_t length = (sizeof MBGEN.pomwt / sizeof MBGEN.pomwt[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetPOMWT : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return MBGEN.pomwt[index];
}

/**************************************************************************/

Float_t TIsajet::GetPOMGEN(Int_t index) const
{
    Int_t length = (sizeof MBGEN.pomgen / sizeof MBGEN.pomgen[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetPOMGEN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return MBGEN.pomgen[index];
}

/**************************************************************************/

void TIsajet::SetMNPOM(Int_t val) 
{
    if (val > MBGEN.mxpom) {
	printf("Error in TIsajet::SetMNPOM : \n");
	printf("You may not set MNPOM to a value larger than MXPOM = %d.\n", MBGEN.mxpom);
	return;
    }

    MBGEN.mnpom = val;
}

/**************************************************************************/

Int_t TIsajet::GetMNPOM() const 
{
    return MBGEN.mnpom;
}

/**************************************************************************/

void TIsajet::SetMXPOM(Int_t val) 
{
    if (val < MBGEN.mnpom) {
	printf("Error in TIsajet::SetMXPOM : \n");
	printf("You may not set MXPOM to a value less than MNPOM = %d.\n", MBGEN.mnpom);
	return;
    }

    MBGEN.mxpom = val;
}

/**************************************************************************/

Int_t TIsajet::GetMXPOM() const 
{
    return MBGEN.mxpom;
}

/**************************************************************************/

Float_t TIsajet::GetPDIFFR() const 
{
    return MBGEN.pdiffr;
}

/**************************************************************************/

Int_t TIsajet::GetNPOM() const 
{
    return MBGEN.npom;
}

/**************************************************************************/

Float_t TIsajet::GetXBARY(Int_t index) const
{
    Int_t length = (sizeof MBGEN.xbary / sizeof MBGEN.xbary[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetXBARY : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return MBGEN.xbary[index];
}

/**************************************************************************/

Float_t TIsajet::GetDXBARY(Int_t index) const
{
    Int_t length = (sizeof MBGEN.dxbary / sizeof MBGEN.dxbary[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetDXBARY : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return MBGEN.dxbary[index];
}

/**************************************************************************/

Float_t TIsajet::GetXPOM(Int_t index1, Int_t index2) const 
{
    Int_t elem_Size = sizeof MBGEN.xpom[0][0];
    Int_t sec_Dim_Length = (sizeof MBGEN.xpom[0] / elem_Size);
    Int_t fst_Dim_Length = (sizeof MBGEN.xpom / (elem_Size * sec_Dim_Length));
    
    if ((index1 < 0) || (index1 >= fst_Dim_Length)) {
	printf ("Error in TIsajet::GetXPOM : \n");
	printf ("Invalid first index %d; range is 0-%d.\n", index1, fst_Dim_Length-1);
	return 0;
    }
    
    if ((index2 < 0) || (index2 >= sec_Dim_Length)) {
	printf ("Error in TIsajet::GetXPOM : \n");
	printf ("Invalid second index %d; range is 0-%d.\n", index2, sec_Dim_Length-1);
	return 0;
    }

    return MBGEN.xpom[index1][index2];
}

/**************************************************************************/

// Ends MBGEN access.
// Begins MGLIMS access.

/**************************************************************************/

void TIsajet::SetEHMGMN(Float_t val) 
{
    if (val > MGLIMS.ehmgmx) {
	printf("Error in TIsajet::SetEHMGMN : \n");
	printf("You may not set EHMGMN to a value larger than EHMGMX = %f.\n", MGLIMS.ehmgmx);
	return;
    }

    MGLIMS.ehmgmn = val;
}

/**************************************************************************/

Float_t TIsajet::GetEHMGMN() const 
{
    return MGLIMS.ehmgmn;
}

/**************************************************************************/

void TIsajet::SetEHMGMX(Float_t val) 
{
    if (val < MGLIMS.ehmgmn) {
	printf("Error in TIsajet::SetEHMGMX : \n");
	printf("You may not set EHMGMX to a value less than EHMGMN = %f.\n", MGLIMS.ehmgmn);
	return;
    }

    MGLIMS.ehmgmx = val;
}

/**************************************************************************/

Float_t TIsajet::GetEHMGMX() const 
{
    return MGLIMS.ehmgmx;
}

/**************************************************************************/

Float_t TIsajet::GetYHMGMN() const 
{
    return MGLIMS.yhmgmn;
}

/**************************************************************************/

Float_t TIsajet::GetYHMGMX() const 
{
    return MGLIMS.yhmgmx;
}

/**************************************************************************/

void TIsajet::SetAMIJMN(Float_t val, Int_t index1, Int_t index2)
{
    Int_t col_num = (sizeof MGLIMS.amijmn[0] / sizeof MGLIMS.amijmn[0][0]);
    Int_t row_num = (sizeof MGLIMS.amijmn / (sizeof MGLIMS.amijmn[0][0] * col_num));
    
    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::SetAMIJMN : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::SetAMIJMN : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return;
    }

    if (val > MGLIMS.amijmx[index1][index2]) {
	printf("Error in TIsajet::SetAMIJMN : \n");
	printf("You may not set AMIJMN to a value larger than AMIJMX = %f.\n", MGLIMS.amijmx[index1][index2]);
	return;
    }
    
    MGLIMS.amijmn[index1][index2] = val;
}

/**************************************************************************/

void TIsajet::SetAllAMIJMN(const Float_t val[MGLIMS.mxlim][MGLIMS.mxlim]) 
{
    for (Int_t i = 0; i < 2; i++) {
	for (Int_t j = 0; j < 3; j++) {
	    SetAMIJMN(val[i][j], i, j);
	}
    }
}

/**************************************************************************/

void TIsajet::SetColumnAMIJMN(const Float_t val[], Int_t col) 
{
    Int_t col_num = (sizeof MGLIMS.amijmn[0] / sizeof MGLIMS.amijmn[0][0]);
    Int_t row_num = (sizeof MGLIMS.amijmn / (sizeof MGLIMS.amijmn[0][0] * col_num));
    
    if ((col < 0) || (col >= col_num)) {
	printf ("Error in TIsajet::SetColumnAMIJMN : \n");
	printf ("Invalid column index %d, range is 0-%d.\n", col, col_num-1);
	return;
    }
    
    for (Int_t i = 0; i < row_num; i++) {
	SetAMIJMN(val[i], i, col);
    }
}

/**************************************************************************/

Float_t TIsajet::GetAMIJMN(Int_t index1, Int_t index2) const 
{
    Int_t col_num = (sizeof MGLIMS.amijmn[0] / sizeof MGLIMS.amijmn[0][0]);
    Int_t row_num = (sizeof MGLIMS.amijmn / (sizeof MGLIMS.amijmn[0][0] * col_num));

    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::GetAMIJMN : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return 0;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::GetAMIJMN : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return 0;
    }

    return MGLIMS.amijmn[index1][index2];
}

/**************************************************************************/

void TIsajet::SetAMIJMX(Float_t val, Int_t index1, Int_t index2)
{
    Int_t col_num = (sizeof MGLIMS.amijmx[0] / sizeof MGLIMS.amijmx[0][0]);
    Int_t row_num = (sizeof MGLIMS.amijmx / (sizeof MGLIMS.amijmx[0][0] * col_num));
    
    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::SetAMIJMX : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::SetAMIJMX : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return;
    }

    if (val < MGLIMS.amijmn[index1][index2]) {
	printf("Error in TIsajet::SetAMIJMX : \n");
	printf("You may not set AMIJMX to a value less than AMIJMN = %f.\n", MGLIMS.amijmn[index1][index2]);
	return;
    }
    
    MGLIMS.amijmx[index1][index2] = val;
}

/**************************************************************************/

void TIsajet::SetAllAMIJMX(const Float_t val[MGLIMS.mxlim][MGLIMS.mxlim]) 
{
    for (Int_t i = 0; i < 2; i++) {
	for (Int_t j = 0; j < 3; j++) {
	    SetAMIJMX(val[i][j], i, j);
	}
    }
}

/**************************************************************************/

void TIsajet::SetColumnAMIJMX(const Float_t val[], Int_t col) 
{
    Int_t col_num = (sizeof MGLIMS.amijmx[0] / sizeof MGLIMS.amijmx[0][0]);
    Int_t row_num = (sizeof MGLIMS.amijmx / (sizeof MGLIMS.amijmx[0][0] * col_num));
    
    if ((col < 0) || (col >= col_num)) {
	printf ("Error in TIsajet::SetColumnAMIJMX : \n");
	printf ("Invalid column index %d, range is 0-%d.\n", col, col_num-1);
	return;
    }
    
    for (Int_t i = 0; i < row_num; i++) {
	SetAMIJMX(val[i], i, col);
    }
}

/**************************************************************************/

Float_t TIsajet::GetAMIJMX(Int_t index1, Int_t index2) const 
{
    Int_t col_num = (sizeof MGLIMS.amijmx[0] / sizeof MGLIMS.amijmx[0][0]);
    Int_t row_num = (sizeof MGLIMS.amijmx / (sizeof MGLIMS.amijmx[0][0] * col_num));

    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::GetAMIJMX : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return 0;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::GetAMIJMX : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return 0;
    }

    return MGLIMS.amijmx[index1][index2];
}

/**************************************************************************/

Bool_t TIsajet::GetFIXMIJ(Int_t index1, Int_t index2) const 
{
    Int_t elem_Size = sizeof MGLIMS.fixmij[0][0];
    Int_t sec_Dim_Length = (sizeof MGLIMS.fixmij[0] / elem_Size);
    Int_t fst_Dim_Length = (sizeof MGLIMS.fixmij / (elem_Size * sec_Dim_Length));
    
    if ((index1 < 0) || (index1 >= fst_Dim_Length)) {
	printf ("Error in TIsajet::GetFIXMIJ : \n");
	printf ("Invalid first index %d; range is 0-%d.\n", index1, fst_Dim_Length-1);
	return 0;
    }
    
    if ((index2 < 0) || (index2 >= sec_Dim_Length)) {
	printf ("Error in TIsajet::GetFIXMIJ : \n");
	printf ("Invalid second index %d; range is 0-%d.\n", index2, sec_Dim_Length-1);
	return 0;
    }

    return MGLIMS.fixmij[index1][index2];
}

/**************************************************************************/

// Ends MGLIMS access.
// Begins NODCAY access.

/**************************************************************************/

void TIsajet::SetNODCAY(Bool_t val) 
{
    NODCAY.nodcay = val;
    setNodcay = true;
}

/**************************************************************************/

Bool_t TIsajet::GetNODCAY() const 
{
    return NODCAY.nodcay;
}

/**************************************************************************/

void TIsajet::SetNOETA(Bool_t val) 
{
    NODCAY.noeta = val;
    setNoeta = true;
}

/**************************************************************************/

Bool_t TIsajet::GetNOETA() const 
{
    return NODCAY.noeta;
}

/**************************************************************************/

void TIsajet::SetNOPI0(Bool_t val) 
{
    NODCAY.nopi0 = val;
    setNopi0 = true;
}

/**************************************************************************/

Bool_t TIsajet::GetNOPI0() const 
{
    return NODCAY.nopi0;
}

/**************************************************************************/

void TIsajet::SetNONUNU(Bool_t val) 
{
    NODCAY.nonunu = val;
}

/**************************************************************************/

Bool_t TIsajet::GetNONUNU() const 
{
    return NODCAY.nonunu;
}

/**************************************************************************/

void TIsajet::SetNOEVOL(Bool_t val) 
{
    NODCAY.noevol = val;
    setNoevolve = true;
}

/**************************************************************************/

Bool_t TIsajet::GetNOEVOL() const 
{
    return NODCAY.noevol;
}

/**************************************************************************/

void TIsajet::SetNOHADR(Bool_t val) 
{
    NODCAY.nohadr = val;
    setNohadron = true;
}

/**************************************************************************/

Bool_t TIsajet::GetNOHADR() const 
{
    return NODCAY.nohadr;
}

/**************************************************************************/

void TIsajet::SetNOGRAV(Bool_t val) 
{
    NODCAY.nograv = val;
}

/**************************************************************************/

Bool_t TIsajet::GetNOGRAV() const 
{
    return NODCAY.nograv;
}

/**************************************************************************/

// Ends NODCAY access.
// Begins PARTCL access.

/**************************************************************************/

Int_t TIsajet::GetNPTCL() const 
{
    return PARTCL.nptcl;
}

/**************************************************************************/

Float_t TIsajet::GetPX(Int_t index) const 
{
    if ((index < 0) || (index >= PARTCL.nptcl)) {
	printf ("Error in TIsajet::GetPX : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, PARTCL.nptcl-1);
	return 0;
    }

    return PARTCL.pptcl[index][0];
}

/**************************************************************************/

Float_t TIsajet::GetPY(Int_t index) const 
{
    if ((index < 0) || (index >= PARTCL.nptcl)) {
	printf ("Error in TIsajet::GetPY : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, PARTCL.nptcl-1);
	return 0;
    }

    return PARTCL.pptcl[index][1];
}

/**************************************************************************/

Float_t TIsajet::GetPZ(Int_t index) const 
{
    if ((index < 0) || (index >= PARTCL.nptcl)) {
	printf ("Error in TIsajet::GetPZ : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, PARTCL.nptcl-1);
	return 0;
    }

    return PARTCL.pptcl[index][2];
}

/**************************************************************************/

Float_t TIsajet::GetP0(Int_t index) const 
{
    if ((index < 0) || (index >= PARTCL.nptcl)) {
	printf ("Error in TIsajet::GetP0 : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, PARTCL.nptcl-1);
	return 0;
    }

    return PARTCL.pptcl[index][3];
}

/**************************************************************************/

Float_t TIsajet::GetMASS(Int_t index) const 
{
    if ((index < 0) || (index >= PARTCL.nptcl)) {
	printf ("Error in TIsajet::GetMASS : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, PARTCL.nptcl-1);
	return 0;
    }

    return PARTCL.pptcl[index][4];
}

/**************************************************************************/

Float_t TIsajet::GetORIG(Int_t index) const 
{
    if ((index < 0) || (index >= PARTCL.nptcl)) {
	printf ("Error in TIsajet::GetORIG : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, PARTCL.nptcl-1);
	return 0;
    }

    return PARTCL.iorig[index];
}

/**************************************************************************/

Float_t TIsajet::GetIDENT(Int_t index) const 
{
    if ((index < 0) || (index >= PARTCL.nptcl)) {
	printf ("Error in TIsajet::GetIDENT : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, PARTCL.nptcl-1);
	return 0;
    }

    return PARTCL.ident[index];
}

/**************************************************************************/

Float_t TIsajet::GetIDCAY(Int_t index) const 
{
    if ((index < 0) || (index >= PARTCL.nptcl)) {
	printf ("Error in TIsajet::GetIDCAY : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, PARTCL.nptcl-1);
	return 0;
    }

    return PARTCL.idcay[index];
}

/**************************************************************************/

// Ends PARTCL access.
// Begins PRIMAR access.

/**************************************************************************/

Int_t TIsajet::GetNJET() const 
{
    return PRIMAR.njet;
}

/**************************************************************************/

Float_t TIsajet::GetSCM() const 
{
    return PRIMAR.scm;
}

/**************************************************************************/

Float_t TIsajet::GetHALFE() const 
{
    return PRIMAR.halfe;
}

/**************************************************************************/

void TIsajet::SetECM(Float_t val) 
{
    if (val < 0) {
	printf ("Error in TIsajet::SetECM :\n");
	printf ("Cannot set energy to a negative value.\n");
	return;
    }
    
    PRIMAR.ecm = val;
    PRIMAR.scm = val*val;
    PRIMAR.halfe = val / 2;
}

/**************************************************************************/

Float_t TIsajet::GetECM() const 
{
    return PRIMAR.ecm;
}

/**************************************************************************/

void TIsajet::SetIDIN(Int_t val, Int_t index)
{
    Int_t length = (sizeof PRIMAR.idin / sizeof PRIMAR.idin[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetIDIN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }

    if ((val = 1120) || (val = 1220) || (val = -1120) || (val = -1220)) {
	PRIMAR.idin[index] = val;
    }
    else {
	printf ("Error in TIsajet::SetIDIN : \n");
	printf ("Invalid input value %d. Possible values are 1120, 1220, -1120, -1220.\n", val);
	return;
    }

    setBeams = true;
}

/**************************************************************************/

void TIsajet::SetIDIN(const Char_t val[], Int_t index)
{
    Int_t length = (sizeof PRIMAR.idin / sizeof PRIMAR.idin[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetIDIN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }

    if (!strcmp("P", val)) {
	PRIMAR.idin[index] = 1120;
    }
    else if (!strcmp("AP", val)) {
	PRIMAR.idin[index] = -1120;
    }
    else if (!strcmp("N", val)) {
	PRIMAR.idin[index] = 1220;
    }
    else if (!strcmp("AN", val)) {
	PRIMAR.idin[index] = -1220;
    }
    else {
	printf ("Error in TIsajet::SetIDIN : \n");
	printf ("Invalid input string %s. Possible strings are P, AP, N, and AN.\n", val);
	return;
    }
}

/**************************************************************************/

Int_t TIsajet::GetIDIN(Int_t index) const 
{
    Int_t length = (sizeof PRIMAR.idin / sizeof PRIMAR.idin[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetIDIN : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return PRIMAR.idin[index];
}

/**************************************************************************/

Int_t TIsajet::GetNEVENT() const 
{
    return PRIMAR.nevent;
}

/**************************************************************************/

void TIsajet::SetNTRIES(Int_t val) 
{
    PRIMAR.ntries = val;
}

/**************************************************************************/

Int_t TIsajet::GetNTRIES() const 
{
    return PRIMAR.ntries;
}

/**************************************************************************/

void TIsajet::SetNSIGMA(Int_t val) 
{
    PRIMAR.nsigma = val;
    setNsigma = true;
}

/**************************************************************************/

Int_t TIsajet::GetNSIGMA() const 
{
    return PRIMAR.nsigma;
}

/**************************************************************************/

// Ends PRIMAR access.
// Begins QCDPAR access.

/**************************************************************************/

void TIsajet::SetALAM(Float_t val) 
{
    QCDPAR.alam = val;
    QCDPAR.alam2 = val*val;
    setLambda = true;
}

/**************************************************************************/

Float_t TIsajet::GetALAM() const 
{
    return QCDPAR.alam;
}

/**************************************************************************/

Float_t TIsajet::GetALAM2() const 
{
    return QCDPAR.alam2;
}

/**************************************************************************/

void TIsajet::SetCUTJET(Float_t val) 
{
    QCDPAR.cutjet = val;
    setCutjet = true;
}

/**************************************************************************/

Float_t TIsajet::GetCUTJET() const 
{
    return QCDPAR.cutjet;
}

/**************************************************************************/

void TIsajet::SetISTRUC(Int_t val)
{
    if ((val < 1) || (val > 6)) {
	printf ("Error in TIsajet::SetISTRUC : \n");
	printf ("Invalid input value %d. Possible values are 1 through 6.\n", val);
	return;
    }
    QCDPAR.istruc = val;
}

/**************************************************************************/

void TIsajet::SetISTRUC(const Char_t val[])
{
    if (!strcmp("OWENS", val)) {
	QCDPAR.istruc = 1;
    }
    else if (!strcmp("BAIER", val)) {
	QCDPAR.istruc = 2;
    }
    else if ((!strcmp("EICHTEN", val)) || (!strcmp("EHLQ", val))) {
	QCDPAR.istruc = 3;
    }
    else if ((!strcmp("DUKE", val)) || (!strcmp("DO", val))) {    
	QCDPAR.istruc = 4;
    }
    else if (!strcmp("CTEQ2L", val)) {
	QCDPAR.istruc = 5;
    }
    else if ((!strcmp("CTEQ", val)) || (!strcmp("CTEQ3L", val))) {    
	QCDPAR.istruc = 6;
    }
    else {
	printf ("Error in TIsajet::SetISTRUC : \n");
	printf ("Invalid input string %s. Possible strings are OWENS, BAIER, EICHTEN, \n", val);
        printf ("EHLQ, DUKE, DO, CTEQ2L, CTEQ, and CTEQ3L.\n");
	return;
    }
}

/**************************************************************************/

Int_t TIsajet::GetISTRUC() const 
{
    return QCDPAR.istruc;
}

/**************************************************************************/

// Ends QCDPAR access.
// Begins QLMASS access.

/**************************************************************************/

void TIsajet::SetAMLEP(Float_t val, Int_t index)
{
    Int_t length = (sizeof QLMASS.amlep / sizeof QLMASS.amlep[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetAMLEP : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }
    
    if (((index < 5) && (index > 7)) &&
	((index < 21) && (index > 26)) &&
	((index < 29) && (index > 36)) &&
	((index < 39) && (index > 40)) &&
	((index < 63) && (index > 71))) 
	{
	    printf ("Error in TIsajet::SetAMLEP : \n");
	    printf ("Index %d may not be set by the user. Valid indices are : \n", index);
	    printf ("5-7, 21-26, 29-36, 39-40 and 63-71.\n");
	    return;
	}
	    
    QLMASS.amlep[index] = val;
}

/**************************************************************************/

Float_t TIsajet::GetAMLEP(Int_t index) const 
{
    Int_t length = (sizeof QLMASS.amlep / sizeof QLMASS.amlep[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetAMLEP : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }
    if   ((index < 5) || (index > 71) ||
	  ((index > 7)  && (index < 21)) || 
	  ((index > 26)  && (index < 29)) || 
	  ((index > 36)  && (index < 39)) || 
	  ((index > 40)  && (index < 63)))   
	{
	    printf ("Possible error in TIsajet::GetAMLEP : \n");
	    printf ("Index %d may not be set by the user. Valid indices are : \n", index);
	    printf ("5-7, 21-26, 29-36, 39-40 and 63-71. \n");
	    printf ("To return the value of this index, use GetAnyAMLEP(Int_t index).\n");
	    return 0;
	}


    return QLMASS.amlep[index];
}

/**************************************************************************/

Float_t TIsajet::GetAnyAMLEP(Int_t index) const 
{
    Int_t length = (sizeof QLMASS.amlep / sizeof QLMASS.amlep[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetAnyAMLEP : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return QLMASS.amlep[index];
}

/**************************************************************************/

void TIsajet::SetTquarkMass(Float_t val) 
{
    QLMASS.amlep[5] = val;
}

/**************************************************************************/

Float_t TIsajet::GetTquarkMass() const 
{
    return QLMASS.amlep[5];
}

/**************************************************************************/

void TIsajet::SetXquarkMass(Float_t val) 
{
    QLMASS.amlep[6] = val;
}

/**************************************************************************/

Float_t TIsajet::GetXquarkMass() const 
{
    return QLMASS.amlep[6];
}

/**************************************************************************/

void TIsajet::SetYquarkMass(Float_t val) 
{
    QLMASS.amlep[7] = val;
}

/**************************************************************************/

Float_t TIsajet::GetYquarkMass() const 
{
    return QLMASS.amlep[7];
}

/**************************************************************************/

void TIsajet::SetUtildeMass(Float_t val) 
{
    QLMASS.amlep[21] = val;
}

/**************************************************************************/

Float_t TIsajet::GetUtildeMass() const 
{
    return QLMASS.amlep[21];
}

/**************************************************************************/

void TIsajet::SetDtildeMass(Float_t val) 
{
    QLMASS.amlep[22] = val;
}

/**************************************************************************/

Float_t TIsajet::GetDtildeMass() const 
{
    return QLMASS.amlep[22];
}

/**************************************************************************/

void TIsajet::SetStildeMass(Float_t val) 
{
    QLMASS.amlep[23] = val;
}

/**************************************************************************/

Float_t TIsajet::GetStildeMass() const 
{
    return QLMASS.amlep[23];
}

/**************************************************************************/

void TIsajet::SetCtildeMass(Float_t val) 
{
    QLMASS.amlep[24] = val;
}

/**************************************************************************/

Float_t TIsajet::GetCtildeMass() const 
{
    return QLMASS.amlep[24];
}

/**************************************************************************/

void TIsajet::SetBtildeMass(Float_t val) 
{
    QLMASS.amlep[25] = val;
}

/**************************************************************************/

Float_t TIsajet::GetBtildeMass() const 
{
    return QLMASS.amlep[25];
}

/**************************************************************************/

void TIsajet::SetTtildeMass(Float_t val) 
{
    QLMASS.amlep[26] = val;
}

/**************************************************************************/

Float_t TIsajet::GetTtildeMass() const 
{
    return QLMASS.amlep[26];
}

/**************************************************************************/

void TIsajet::SetGtildeMass(Float_t val) 
{
    QLMASS.amlep[29] = val;
}

/**************************************************************************/

Float_t TIsajet::GetGtildeMass() const 
{
    return QLMASS.amlep[29];
}

/**************************************************************************/

void TIsajet::SetGammatildeMass(Float_t val) 
{
    QLMASS.amlep[30] = val;
}

/**************************************************************************/

Float_t TIsajet::GetGammatildeMass() const 
{
    return QLMASS.amlep[30];
}

/**************************************************************************/

void TIsajet::SetNuEtildeMass(Float_t val) 
{
    QLMASS.amlep[31] = val;
}

/**************************************************************************/

Float_t TIsajet::GetNuEtildeMass() const 
{
    return QLMASS.amlep[31];
}

/**************************************************************************/

void TIsajet::SetEtildeMass(Float_t val) 
{
    QLMASS.amlep[32] = val;
}

/**************************************************************************/

Float_t TIsajet::GetEtildeMass() const 
{
    return QLMASS.amlep[32];
}

/**************************************************************************/

void TIsajet::SetNuMutildeMass(Float_t val) 
{
    QLMASS.amlep[33] = val;
}

/**************************************************************************/

Float_t TIsajet::GetNuMutildeMass() const 
{
    return QLMASS.amlep[33];
}

/**************************************************************************/

void TIsajet::SetMutildeMass(Float_t val) 
{
    QLMASS.amlep[34] = val;
}

/**************************************************************************/

Float_t TIsajet::GetMutildeMass() const 
{
    return QLMASS.amlep[34];
}

/**************************************************************************/

void TIsajet::SetNuTautildeMass(Float_t val) 
{
    QLMASS.amlep[35] = val;
}

/**************************************************************************/

Float_t TIsajet::GetNuTautildeMass() const 
{
    return QLMASS.amlep[35];
}

/**************************************************************************/

void TIsajet::SetTautildeMass(Float_t val) 
{
    QLMASS.amlep[36] = val;
}

/**************************************************************************/

Float_t TIsajet::GetTautildeMass() const 
{
    return QLMASS.amlep[36];
}

/**************************************************************************/

void TIsajet::SetWplustildeMass(Float_t val) 
{
    QLMASS.amlep[39] = val;
}

/**************************************************************************/

Float_t TIsajet::GetWplustildeMass() const 
{
    return QLMASS.amlep[39];
}

/**************************************************************************/

void TIsajet::SetZ0tildeMass(Float_t val) 
{
    QLMASS.amlep[40] = val;
}

/**************************************************************************/

Float_t TIsajet::GetZ0tildeMass() const 
{
    return QLMASS.amlep[40];
}

/**************************************************************************/

void TIsajet::SetHiggsMesonMass(Float_t val, Int_t index)
{
    if ((index < 1) || (index > 9)) {
	printf ("Error in TIsajet::SetHiggsMesonMass : \n");
	printf ("Invalid Higgs meson number index %d; range is 1-9.\n", index);
	return;
    }
    
    QLMASS.amlep[62 + index] = val;
}

/**************************************************************************/

Float_t TIsajet::GetHiggsMesonMass(Int_t index) const
{
    if ((index < 1) || (index > 9)) {
	printf ("Error in TIsajet::GetHiggsMesonMass : \n");
	printf ("Invalid Higgs meson number index %d; range is 1-9.\n", index);
	return 0;
    }
    
    return QLMASS.amlep[62 + index];
}

/**************************************************************************/

Int_t TIsajet::GetNQLEP() const
{
    return QLMASS.nqlep;
}

/**************************************************************************/

Int_t TIsajet::GetNMES() const
{
    return QLMASS.nmes;
}

/**************************************************************************/

Int_t TIsajet::GetNBARY() const
{
    return QLMASS.nbary;
}

/**************************************************************************/

// Ends QLMASS access.
// Begins SEED access.

/**************************************************************************/

void TIsajet::SetSEED(const Char_t val[24]) 
{
    Int_t length = (sizeof SEED.xseed / sizeof SEED.xseed[0]);    
    
    for (Int_t i = 0; i < length; i++) {
	SEED.xseed[i] = val[i];
    }
}

/**************************************************************************/

Char_t* TIsajet::GetSEED() const 
{
    return SEED.xseed;
}

/**************************************************************************/

// Ends SEED access - short and sweet, wasn't it?
// Begins SUGNU access, an entirely different business.

/**************************************************************************/

void TIsajet::SetXNUSUG(Float_t val, Int_t index)
{
    Int_t length = (sizeof SUGNU.xnusug / sizeof SUGNU.xnusug[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::SetXNUSUG : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return;
    }
    
    SUGNU.xnusug[index] = val;
}

/**************************************************************************/

Float_t TIsajet::GetXNUSUG(Int_t index) const 
{
    Int_t length = (sizeof SUGNU.xnusug / sizeof SUGNU.xnusug[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetXNUSUG : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return SUGNU.xnusug[index];
}

/**************************************************************************/

void TIsajet::SetGauginoMass(Float_t val, Int_t index)
{
    if ((index < 1) || (index > 9)) {
	printf ("Error in TIsajet::SetGauginoMass : \n");
	printf ("Invalid gaugino number %d; range is 1-3.\n", index);
	return;
    }
    
    SUGNU.xnusug[index-1] = val;
}

/**************************************************************************/
 
Float_t TIsajet::GetGauginoMass(Int_t index) const
{
    if ((index < 1) || (index > 3)) {
	printf ("Error in TIsajet::GetGauginoMass : \n");
	printf ("Invalid gaugino number %d; range is 1-3.\n", index);
	return 0;
    }
    
    return SUGNU.xnusug[index-1];
}

/**************************************************************************/

void TIsajet::SetAtau(Float_t val) 
{
    SUGNU.xnusug[3] = val;
}

/**************************************************************************/

Float_t TIsajet::GetAtau() const 
{
    return SUGNU.xnusug[3];
}

/**************************************************************************/

void TIsajet::SetAb(Float_t val) 
{
    SUGNU.xnusug[4] = val;
}

/**************************************************************************/

Float_t TIsajet::GetAb() const 
{
    return SUGNU.xnusug[4];
}

/**************************************************************************/

void TIsajet::SetAt(Float_t val) 
{
    SUGNU.xnusug[5] = val;
}

/**************************************************************************/

Float_t TIsajet::GetAt() const 
{
    return SUGNU.xnusug[5];
}

/**************************************************************************/

void TIsajet::SetHiggsDmass(Float_t val) 
{
    SUGNU.xnusug[6] = val;
}

/**************************************************************************/

Float_t TIsajet::GetHiggsDmass() const 
{
    return SUGNU.xnusug[6];
}

/**************************************************************************/

void TIsajet::SetHiggsUmass(Float_t val) 
{
    SUGNU.xnusug[7] = val;
}

/**************************************************************************/

Float_t TIsajet::GetHiggsUmass() const 
{
    return SUGNU.xnusug[7];
}

/**************************************************************************/

void TIsajet::SetERmass(Float_t val) 
{
    SUGNU.xnusug[8] = val;
}

/**************************************************************************/

Float_t TIsajet::GetERmass() const 
{
    return SUGNU.xnusug[8];
}

/**************************************************************************/

void TIsajet::SetELmass(Float_t val) 
{
    SUGNU.xnusug[9] = val;
}

/**************************************************************************/

Float_t TIsajet::GetELmass() const 
{
    return SUGNU.xnusug[9];
}

/**************************************************************************/

void TIsajet::SetDRmass(Float_t val) 
{
    SUGNU.xnusug[10] = val;
}

/**************************************************************************/

Float_t TIsajet::GetDRmass() const 
{
    return SUGNU.xnusug[10];
}

/**************************************************************************/

void TIsajet::SetURmass(Float_t val) 
{
    SUGNU.xnusug[11] = val;
}

/**************************************************************************/

Float_t TIsajet::GetURmass() const 
{
    return SUGNU.xnusug[11];
}

/**************************************************************************/

void TIsajet::SetULmass(Float_t val) 
{
    SUGNU.xnusug[12] = val;
}

/**************************************************************************/

Float_t TIsajet::GetULmass() const 
{
    return SUGNU.xnusug[12];
}

/**************************************************************************/

void TIsajet::SetTauRmass(Float_t val) 
{
    SUGNU.xnusug[13] = val;
}

/**************************************************************************/

Float_t TIsajet::GetTauRmass() const 
{
    return SUGNU.xnusug[13];
}

/**************************************************************************/

void TIsajet::SetTauLmass(Float_t val) 
{
    SUGNU.xnusug[14] = val;
}

/**************************************************************************/

Float_t TIsajet::GetTauLmass() const 
{
    return SUGNU.xnusug[14];
}

/**************************************************************************/

void TIsajet::SetBRmass(Float_t val) 
{
    SUGNU.xnusug[15] = val;
}

/**************************************************************************/

Float_t TIsajet::GetBRmass() const 
{
    return SUGNU.xnusug[15];
}

/**************************************************************************/

void TIsajet::SetTRmass(Float_t val) 
{
    SUGNU.xnusug[16] = val;
}

/**************************************************************************/

Float_t TIsajet::GetTRmass() const 
{
    return SUGNU.xnusug[16];
}

/**************************************************************************/

void TIsajet::SetTLmass(Float_t val) 
{
    SUGNU.xnusug[17] = val;
}

/**************************************************************************/

Float_t TIsajet::GetTLmass() const 
{
    return SUGNU.xnusug[17];
}

/**************************************************************************/

// Ends XNUSUG access.
// Begins TCPAR access.

/**************************************************************************/

void TIsajet::SetTCMRHO(Float_t val) 
{
    TCPAR.tcmrho = val;
}

/**************************************************************************/

Float_t TIsajet::GetTCMRHO() const 
{
    return TCPAR.tcmrho;
}

/**************************************************************************/

void TIsajet::SetTCGRHO(Float_t val) 
{
    TCPAR.tcgrho = val;
}

/**************************************************************************/

Float_t TIsajet::GetTCGRHO() const 
{
    return TCPAR.tcgrho;
}

/**************************************************************************/

// Ends TCPAR access.
// Begins TYPES access.

/**************************************************************************/

Int_t TIsajet::GetLOC(Int_t index) const 
{
    Int_t length = (sizeof TYPES.loc / sizeof TYPES.loc[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetLOC : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return TYPES.loc[index];
}

/**************************************************************************/

Int_t TIsajet::GetNTYP() const 
{
    return TYPES.ntyp;
}

/**************************************************************************/

Int_t TIsajet::GetNJTTYP(Int_t index) const 
{
    Int_t length = (sizeof TYPES.njttyp / sizeof TYPES.njttyp[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetNJTYP : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return TYPES.njttyp[index];
}

/**************************************************************************/

Int_t TIsajet::GetNWWTYP(Int_t index) const 
{
    Int_t length = (sizeof TYPES.nwwtyp / sizeof TYPES.nwwtyp[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetNWWTYP : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return TYPES.nwwtyp[index];
}

/**************************************************************************/

Int_t TIsajet::GetNWMODE(Int_t index) const 
{
    Int_t length = (sizeof TYPES.nwmode / sizeof TYPES.nwmode[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetNWMODE : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return TYPES.nwmode[index];
}

/**************************************************************************/

// Ends TYPES access.
// Begins XMSSM access.

/**************************************************************************/

Bool_t TIsajet::GetGOMSSM() const
{
    return XMSSM.gomssm;
}

/**************************************************************************/

Bool_t TIsajet::GetGOSUG() const
{
    return XMSSM.gosug;
}

/**************************************************************************/

Bool_t TIsajet::GetGOGMSB() const
{
    return XMSSM.gogmsb;
}

/**************************************************************************/

Bool_t TIsajet::GetGOAMSB() const
{
    return XMSSM.goamsb;
}

/**************************************************************************/

Bool_t TIsajet::GetAL3UNI() const
{
    return XMSSM.al3uni;
}

/**************************************************************************/

void TIsajet::SetXGLSS(Float_t val) 
{
    XMSSM.xglss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXGLSS() const 
{
    return XMSSM.xglss;
}

/**************************************************************************/

void TIsajet::SetXMUSS(Float_t val) 
{
    XMSSM.xmuss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXMUSS() const 
{
    return XMSSM.xmuss;
}

/**************************************************************************/

void TIsajet::SetXHASS(Float_t val) 
{
    XMSSM.xhass = val;
}

/**************************************************************************/

Float_t TIsajet::GetXHASS() const 
{
    return XMSSM.xhass;
}

/**************************************************************************/

void TIsajet::SetXTBSS(Float_t val) 
{
    XMSSM.xtbss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXTBSS() const 
{
    return XMSSM.xtbss;
}

/**************************************************************************/

void TIsajet::SetXQ1SS(Float_t val) 
{
    XMSSM.xq1ss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXQ1SS() const 
{
    return XMSSM.xq1ss;
}

/**************************************************************************/

void TIsajet::SetXDRSS(Float_t val) 
{
    XMSSM.xdrss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXDRSS() const 
{
    return XMSSM.xdrss;
}

/**************************************************************************/

void TIsajet::SetXURSS(Float_t val) 
{
    XMSSM.xurss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXURSS() const 
{
    return XMSSM.xurss;
}

/**************************************************************************/

void TIsajet::SetXL1SS(Float_t val) 
{
    XMSSM.xl1ss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXL1SS() const 
{
    return XMSSM.xl1ss;
}

/**************************************************************************/

void TIsajet::SetXERSS(Float_t val) 
{
    XMSSM.xerss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXERSS() const 
{
    return XMSSM.xerss;
}

/**************************************************************************/

void TIsajet::SetXQ2SS(Float_t val) 
{
    XMSSM.xq2ss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXQ2SS() const 
{
    return XMSSM.xq2ss;
}

/**************************************************************************/

void TIsajet::SetXSRSS(Float_t val) 
{
    XMSSM.xsrss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXSRSS() const 
{
    return XMSSM.xsrss;
}

/**************************************************************************/

void TIsajet::SetXCRSS(Float_t val) 
{
    XMSSM.xcrss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXCRSS() const 
{
    return XMSSM.xcrss;
}

/**************************************************************************/

void TIsajet::SetXL2SS(Float_t val) 
{
    XMSSM.xl2ss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXL2SS() const 
{
    return XMSSM.xl2ss;
}

/**************************************************************************/

void TIsajet::SetXMRSS(Float_t val) 
{
    XMSSM.xmrss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXMRSS() const 
{
    return XMSSM.xmrss;
}

/**************************************************************************/

void TIsajet::SetXQ3SS(Float_t val) 
{
    XMSSM.xq3ss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXQ3SS() const 
{
    return XMSSM.xq3ss;
}

/**************************************************************************/

void TIsajet::SetXBRSS(Float_t val) 
{
    XMSSM.xbrss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXBRSS() const 
{
    return XMSSM.xbrss;
}

/**************************************************************************/

void TIsajet::SetXTRSS(Float_t val) 
{
    XMSSM.xtrss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXTRSS() const 
{
    return XMSSM.xtrss;
}

/**************************************************************************/

void TIsajet::SetXL3SS(Float_t val) 
{
    XMSSM.xl3ss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXL3SS() const 
{
    return XMSSM.xl3ss;
}

/**************************************************************************/

void TIsajet::SetXTARSS(Float_t val) 
{
    XMSSM.xtarss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXTARSS() const 
{
    return XMSSM.xtarss;
}

/**************************************************************************/

void TIsajet::SetXATSS(Float_t val) 
{
    XMSSM.xatss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXATSS() const 
{
    return XMSSM.xatss;
}

/**************************************************************************/

void TIsajet::SetXABSS(Float_t val) 
{
    XMSSM.xabss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXABSS() const 
{
    return XMSSM.xabss;
}

/**************************************************************************/

void TIsajet::SetXATASS(Float_t val) 
{
    XMSSM.xatass = val;
}

/**************************************************************************/

Float_t TIsajet::GetXATASS() const 
{
    return XMSSM.xatass;
}

/**************************************************************************/

void TIsajet::SetXM1SS(Float_t val) 
{
    XMSSM.xm1ss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXM1SS() const 
{
    return XMSSM.xm1ss;
}

/**************************************************************************/

void TIsajet::SetXM2SS(Float_t val) 
{
    XMSSM.xm2ss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXM2SS() const 
{
    return XMSSM.xm2ss;
}

/**************************************************************************/

void TIsajet::SetXM0SU(Float_t val) 
{
    XMSSM.xm0su = val;
}

/**************************************************************************/

Float_t TIsajet::GetXM0SU() const 
{
    return XMSSM.xm0su;
}

/**************************************************************************/

void TIsajet::SetXMHSU(Float_t val) 
{
    XMSSM.xmhsu = val;
}

/**************************************************************************/

Float_t TIsajet::GetXMHSU() const 
{
    return XMSSM.xmhsu;
}

/**************************************************************************/

void TIsajet::SetXA0SU(Float_t val) 
{
    XMSSM.xa0su = val;
}

/**************************************************************************/

Float_t TIsajet::GetXA0SU() const 
{
    return XMSSM.xa0su;
}

/**************************************************************************/

void TIsajet::SetXTGBSU(Float_t val) 
{
    XMSSM.xtgbsu = val;
}

/**************************************************************************/

Float_t TIsajet::GetXTGBSU() const 
{
    return XMSSM.xtgbsu;
}

/**************************************************************************/

void TIsajet::SetXSMUSU(Float_t val) 
{
    XMSSM.xsmusu = val;
}

/**************************************************************************/

Float_t TIsajet::GetXSMUSU() const 
{
    return XMSSM.xsmusu;
}

/**************************************************************************/

void TIsajet::SetXLAMGM(Float_t val) 
{
    XMSSM.xlamgm = val;
}

/**************************************************************************/

Float_t TIsajet::GetXLAMGM() const 
{
    return XMSSM.xlamgm;
}

/**************************************************************************/

void TIsajet::SetXMESGM(Float_t val) 
{
    XMSSM.xmesgm = val;
}

/**************************************************************************/

Float_t TIsajet::GetXMESGM() const 
{
    return XMSSM.xmesgm;
}

/**************************************************************************/

void TIsajet::SetXN5GM(Float_t val) 
{
    XMSSM.xn5gm = val;
}

/**************************************************************************/

Float_t TIsajet::GetXN5GM() const 
{
    return XMSSM.xn5gm;
}

/**************************************************************************/

void TIsajet::SetXCMGV(Float_t val) 
{
    XMSSM.xcmgv = val;
}

/**************************************************************************/

Float_t TIsajet::GetXCMGV() const 
{
    return XMSSM.xcmgv;
}

/**************************************************************************/

void TIsajet::SetMGVTO(Float_t val) 
{
    XMSSM.mgvto = val;
}

/**************************************************************************/

Float_t TIsajet::GetMGVTO() const 
{
    return XMSSM.mgvto;
}

/**************************************************************************/

void TIsajet::SetXRSLGM(Float_t val) 
{
    XMSSM.xrslgm = val;
}

/**************************************************************************/

Float_t TIsajet::GetXRSLGM() const 
{
    return XMSSM.xrslgm;
}

/**************************************************************************/

void TIsajet::SetXDHDGM(Float_t val) 
{
    XMSSM.xdhdgm = val;
}

/**************************************************************************/

Float_t TIsajet::GetXDHDGM() const 
{
    return XMSSM.xdhdgm;
}

/**************************************************************************/

void TIsajet::SetXDHUGM(Float_t val) 
{
    XMSSM.xdhugm = val;
}

/**************************************************************************/

Float_t TIsajet::GetXDHUGM() const 
{
    return XMSSM.xdhugm;
}

/**************************************************************************/

void TIsajet::SetXDYGM(Float_t val) 
{
    XMSSM.xdygm = val;
}

/**************************************************************************/

Float_t TIsajet::GetXDYGM() const 
{
    return XMSSM.xdygm;
}

/**************************************************************************/

void TIsajet::SetXN51GM(Float_t val) 
{
    XMSSM.xn51gm = val;
}

/**************************************************************************/

Float_t TIsajet::GetXN51GM() const 
{
    return XMSSM.xn51gm;
}

/**************************************************************************/

void TIsajet::SetXN52GM(Float_t val) 
{
    XMSSM.xn52gm = val;
}

/**************************************************************************/

Float_t TIsajet::GetXN52GM() const 
{
    return XMSSM.xn52gm;
}

/**************************************************************************/

void TIsajet::SetXN53GM(Float_t val) 
{
    XMSSM.xn53gm = val;
}

/**************************************************************************/

Float_t TIsajet::GetXN53GM() const 
{
    return XMSSM.xn53gm;
}

/**************************************************************************/

void TIsajet::SetXMN3NR(Float_t val) 
{
    XMSSM.xmn3nr = val;
}

/**************************************************************************/

Float_t TIsajet::GetXMN3NR() const 
{
    return XMSSM.xmn3nr;
}

/**************************************************************************/

void TIsajet::SetXMAJNR(Float_t val) 
{
    XMSSM.xmajnr = val;
}

/**************************************************************************/

Float_t TIsajet::GetXMAJNR() const 
{
    return XMSSM.xmajnr;
}

/**************************************************************************/

void TIsajet::SetXANSS(Float_t val) 
{
    XMSSM.xanss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXANSS() const 
{
    return XMSSM.xanss;
}

/**************************************************************************/

void TIsajet::SetXNRSS(Float_t val) 
{
    XMSSM.xnrss = val;
}

/**************************************************************************/

Float_t TIsajet::GetXNRSS() const 
{
    return XMSSM.xnrss;
}

/**************************************************************************/

void TIsajet::SetXSBCS(Float_t val) 
{
    XMSSM.xsbcs = val;
}

/**************************************************************************/

Float_t TIsajet::GetXSBCS() const 
{
    return XMSSM.xsbcs;
}

/**************************************************************************/

// Ends XMSSM access.
// Begins XTYPES access.

/**************************************************************************/

Char_t* TIsajet::GetPARTYP(Int_t index) const 
{
    Int_t length = (sizeof XTYPES.partyp / sizeof XTYPES.partyp[0]);    
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetPARTYP : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return XTYPES.partyp[index];
}

/**************************************************************************/

void TIsajet::SetTITLE(Char_t *val) 
{
    title = XTYPES.title = val;
}

/**************************************************************************/
Char_t* TIsajet::GetTITLE() const 
{
    return XTYPES.title;
}

/**************************************************************************/

void TIsajet::SetJETYP(Int_t index, Char_t val[]) 
{
    Int_t col_num = (sizeof XTYPES.jetyp[0] / sizeof XTYPES.jetyp[0][0]);
    Int_t row_num = (sizeof XTYPES.jetyp / (sizeof XTYPES.jetyp[0][0] * col_num));

    if ((index < 0) || (index >= row_num)) {
	printf ("Error in TIsajet::SetJETYP : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index, row_num-1);
	return;
    }

    if (TYPES.njttyp[index] >= col_num) {
	printf ("Error in TIsajet::SetJETYP : \n");
	printf ("Cannot set more than %d jet types.\n", col_num-1);
	return;
    }
    
    if ((!strcmp(val, "ALL")) || (!strcmp(val, "GL")) ||
	(!strcmp(val, "UP")) || (!strcmp(val, "UB")) || 
	(!strcmp(val, "DN")) || (!strcmp(val, "DB")) || 
	(!strcmp(val, "ST")) || (!strcmp(val, "SB")) || 
	(!strcmp(val, "CH")) || (!strcmp(val, "CB")) || 
	(!strcmp(val, "BT")) || (!strcmp(val, "BB")) || 
	(!strcmp(val, "TP")) || (!strcmp(val, "TB")) || 
	(!strcmp(val, "X")) || (!strcmp(val, "XB")) || 
	(!strcmp(val, "Y")) || (!strcmp(val, "YB")) || 
	(!strcmp(val, "E-")) || (!strcmp(val, "E+")) || 
	(!strcmp(val, "MU-")) || (!strcmp(val, "MU+")) || 
	(!strcmp(val, "TAU-")) || (!strcmp(val, "TAU+")) || 
	(!strcmp(val, "NUS")) || (!strcmp(val, "GM")) || 
	(!strcmp(val, "W+")) || (!strcmp(val, "W-")) || 
	(!strcmp(val, "Z0"))) {
	
	XTYPES.jetyp[index][TYPES.njttyp[index]++] = val;
    }
    else {
	printf ("Error in TIsajet::SetJETYP : \n");
	printf ("Invalid jet type %s; valid types are\n", val);
	printf ("ALL, GL, UP, UB, DN, DB, ST, SB,\n");
	printf ("CH, CB, BT, BB, TP, TB, X, XB, Y, YB,\n");
	printf ("E-, E+, MU-, MU+, TAU-, TAU+, NUS, GM,\n");
	printf ("W+, W- and Z0.\n");
	return;
    }

    if (index == 0) setJettype1 = true;
    else if (index == 1) setJettype2 = true;
    else if (index == 2) setJettype3 = true;
}

/**************************************************************************/

Char_t* TIsajet::GetJETYP(Int_t index1, Int_t index2) const 
{
    Int_t col_num = (sizeof XTYPES.jetyp[0] / sizeof XTYPES.jetyp[0][0]);
    Int_t row_num = (sizeof XTYPES.jetyp / (sizeof XTYPES.jetyp[0][0] * col_num));

    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::GetJETYP : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return 0;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::GetJETYP : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return 0;
    }

    return XTYPES.jetyp[index1][index2];
}

/**************************************************************************/

void TIsajet::SetWWTYP(Char_t val[], Int_t index1, Int_t index2)
{
    Int_t col_num = (sizeof XTYPES.wwtyp[0] / sizeof XTYPES.wwtyp[0][0]);
    Int_t row_num = (sizeof XTYPES.wwtyp / (sizeof XTYPES.wwtyp[0][0] * col_num));
    
    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::SetWWTYP : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::SetWWTYP : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return;
    }
    
    XTYPES.wwtyp[index1][index2] = val;
}

/**************************************************************************/

void TIsajet::SetAllWWTYP(Char_t* val[2][30]) 
{
    for (Int_t i = 0; i < 2; i++) {
	for (Int_t j = 0; j < 30; j++) {
	    SetWWTYP(val[i][j], i, j);
	}
    }
}

/**************************************************************************/

void TIsajet::SetColumnWWTYP(Char_t* val[], Int_t col) 
{
    Int_t col_num = (sizeof XTYPES.wwtyp[0] / sizeof XTYPES.wwtyp[0][0]);
    Int_t row_num = (sizeof XTYPES.wwtyp / (sizeof XTYPES.wwtyp[0][0] * col_num));
    
    if ((col < 0) || (col >= col_num)) {
	printf ("Error in TIsajet::SetColumnWWTYP : \n");
	printf ("Invalid column index %d, range is 0-%d.\n", col, col_num-1);
	return;
    }
    
    for (Int_t i = 0; i < row_num; i++) {
	SetWWTYP(val[i], i, col);
    }
}

/**************************************************************************/

Char_t* TIsajet::GetWWTYP(Int_t index1, Int_t index2) const 
{
    Int_t col_num = (sizeof XTYPES.wwtyp[0] / sizeof XTYPES.wwtyp[0][0]);
    Int_t row_num = (sizeof XTYPES.wwtyp / (sizeof XTYPES.wwtyp[0][0] * col_num));

    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::GetWWTYP : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return 0;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::GetWWTYP : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return 0;
    }

    return XTYPES.wwtyp[index1][index2];
}

/**************************************************************************/

void TIsajet::SetWMODES(Char_t val[], Int_t index1, Int_t index2)
{
    Int_t col_num = (sizeof XTYPES.wmodes[0] / sizeof XTYPES.wmodes[0][0]);
    Int_t row_num = (sizeof XTYPES.wmodes / (sizeof XTYPES.wmodes[0][0] * col_num));
    
    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::SetWMODES : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::SetWMODES : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return;
    }
    
    XTYPES.wmodes[index1][index2] = val;
}

/**************************************************************************/

void TIsajet::SetAllWMODES(Char_t* val[2][30]) 
{
    for (Int_t i = 0; i < 2; i++) {
	for (Int_t j = 0; j < 30; j++) {
	    SetWMODES(val[i][j], i, j);
	}
    }
}

/**************************************************************************/

void TIsajet::SetColumnWMODES(Char_t* val[], Int_t col) 
{
    Int_t col_num = (sizeof XTYPES.wmodes[0] / sizeof XTYPES.wmodes[0][0]);
    Int_t row_num = (sizeof XTYPES.wmodes / (sizeof XTYPES.wmodes[0][0] * col_num));
    
    if ((col < 0) || (col >= col_num)) {
	printf ("Error in TIsajet::SetColumnWMODES : \n");
	printf ("Invalid column index %d, range is 0-%d.\n", col, col_num-1);
	return;
    }
    
    for (Int_t i = 0; i < row_num; i++) {
	SetWMODES(val[i], i, col);
    }
}

/**************************************************************************/

Char_t* TIsajet::GetWMODES(Int_t index1, Int_t index2) const 
{
    Int_t col_num = (sizeof XTYPES.wmodes[0] / sizeof XTYPES.wmodes[0][0]);
    Int_t row_num = (sizeof XTYPES.wmodes / (sizeof XTYPES.wmodes[0][0] * col_num));

    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::GetWMODES : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return 0;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::GetWMODES : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return 0;
    }

    return XTYPES.wmodes[index1][index2];
}

/**************************************************************************/

// Ends XTYPES access.
// Begins WCON access.

/**************************************************************************/

void TIsajet::SetSIN2W(Float_t val)
{
    WCON.sin2w = val;
}

/**************************************************************************/

Float_t TIsajet::GetSIN2W() const
{
    return WCON.sin2w;
}

/**************************************************************************/

void TIsajet::SetWMASS(Float_t w, Float_t z)
{

// This is how the FORTRAN does it. Don't ask me why.

    WCON.wmass[0] = 0;
    WCON.wmass[1] = WCON.wmass[2] = w;
    WCON.wmass[3] = z;
}

/**************************************************************************/

Float_t TIsajet::GetWMASS(Int_t index) const 
{
    Int_t length = (sizeof WCON.wmass / sizeof WCON.wmass[0]);
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetWMASS : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return WCON.wmass[index];
}

/**************************************************************************/

void TIsajet::SetWMass(Float_t val)
{
    WCON.wmass[1] = WCON.wmass[2] = val;
}

/**************************************************************************/

void TIsajet::SetZMass(Float_t val)
{
    WCON.wmass[3]  = val;
}

/**************************************************************************/

Float_t TIsajet::GetWGAM(Int_t index) const 
{
    Int_t length = (sizeof WCON.wgam / sizeof WCON.wgam[0]);
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetWGAM : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return WCON.wgam[index];
}

/**************************************************************************/

Float_t TIsajet::GetAQ(Int_t index1, Int_t index2) const 
{
    Int_t col_num = (sizeof WCON.aq[0] / sizeof WCON.aq[0][0]);
    Int_t row_num = (sizeof WCON.aq / (sizeof WCON.aq[0][0] * col_num));

    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::GetAQ : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return 0;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::GetAQ : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return 0;
    }

    return WCON.aq[index1][index2];
}

/**************************************************************************/

Float_t TIsajet::GetBQ(Int_t index1, Int_t index2) const 
{
    Int_t col_num = (sizeof WCON.bq[0] / sizeof WCON.bq[0][0]);
    Int_t row_num = (sizeof WCON.bq / (sizeof WCON.bq[0][0] * col_num));

    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::GetBQ : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return 0;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::GetBQ : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return 0;
    }

    return WCON.bq[index1][index2];
}

/**************************************************************************/

Float_t TIsajet::GetCOUT(Int_t index) const 
{
    Int_t length = (sizeof WCON.cout / sizeof WCON.cout[0]);
    if ((index < 0) || (index >= length)) {
	printf ("Error in TIsajet::GetCOUT : \n");
	printf ("Invalid array index %d; range is 0-%d.\n", index, length-1);
	return 0;
    }

    return WCON.cout[index];
}

/**************************************************************************/

Int_t TIsajet::GetMATCH() const
{
    return WCON.match;
}

/**************************************************************************/

Float_t TIsajet::GetWCBR(Int_t index1, Int_t index2) const 
{
    Int_t col_num = (sizeof WCON.wcbr[0] / sizeof WCON.wcbr[0][0]);
    Int_t row_num = (sizeof WCON.wcbr / (sizeof WCON.wcbr[0][0] * col_num));

    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::GetWCBR : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return 0;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::GetWCBR : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return 0;
    }

    return WCON.wcbr[index1][index2];
}

/**************************************************************************/

void TIsajet::SetCUTOFF(Float_t val)
{
    WCON.cutoff = val;
}

/**************************************************************************/

Float_t TIsajet::GetCUTOFF() const
{
    return WCON.cutoff;
}

/**************************************************************************/

void TIsajet::SetCUTPOW(Float_t val)
{
    WCON.cutpow = val;
}

/**************************************************************************/

Float_t TIsajet::GetCUTPOW() const
{
    return WCON.cutpow;
}

/**************************************************************************/

Float_t TIsajet::GetTBRWW(Int_t index1, Int_t index2) const 
{
    Int_t col_num = (sizeof WCON.tbrww[0] / sizeof WCON.tbrww[0][0]);
    Int_t row_num = (sizeof WCON.tbrww / (sizeof WCON.tbrww[0][0] * col_num));

    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::GetTBRWW : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return 0;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::GetTBRWW : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return 0;
    }

    return WCON.tbrww[index1][index2];
}

/**************************************************************************/

Float_t TIsajet::GetRBRWW(Int_t index1, Int_t index2, Int_t index3) const 
{
    Int_t elem_Size = sizeof WCON.rbrww[0][0][0];
    Int_t thd_Dim_Length = (sizeof WCON.rbrww[0][0] / elem_Size);
    Int_t sec_Dim_Length = (sizeof WCON.rbrww[0] / (elem_Size * thd_Dim_Length));
    Int_t fst_Dim_Length = (sizeof WCON.rbrww / (elem_Size * thd_Dim_Length * sec_Dim_Length));
    
    if ((index1 < 0) || (index1 >= fst_Dim_Length)) {
	printf ("Error in TIsajet::GetRBRWW : \n");
	printf ("Invalid first index %d; range is 0-%d.\n", index1, fst_Dim_Length-1);
	return 0;
    }
    
    if ((index2 < 0) || (index2 >= sec_Dim_Length)) {
	printf ("Error in TIsajet::GetRBRWW : \n");
	printf ("Invalid second index %d; range is 0-%d.\n", index2, sec_Dim_Length-1);
	return 0;
    }

    if ((index3 < 0) || (index3 >= thd_Dim_Length)) {
	printf ("Error in TIsajet::GetRBRWW : \n");
	printf ("Invalid third index %d; range is 0-%d.\n", index3, thd_Dim_Length-1);
	return 0;
    }

    return WCON.rbrww[index1][index2][index3];

}

/**************************************************************************/

Float_t TIsajet::GetEZ() const
{
    return WCON.ez;
}

/**************************************************************************/

Float_t TIsajet::GetAQDP(Int_t index1, Int_t index2) const 
{
    Int_t col_num = (sizeof WCON.aqdp[0] / sizeof WCON.aqdp[0][0]);
    Int_t row_num = (sizeof WCON.aqdp / (sizeof WCON.aqdp[0][0] * col_num));

    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::GetAQDP : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return 0;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::GetAQDP : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return 0;
    }

    return WCON.aqdp[index1][index2];
}

/**************************************************************************/

Float_t TIsajet::GetBQDP(Int_t index1, Int_t index2) const 
{
    Int_t col_num = (sizeof WCON.bqdp[0] / sizeof WCON.bqdp[0][0]);
    Int_t row_num = (sizeof WCON.bqdp / (sizeof WCON.bqdp[0][0] * col_num));

    if ((index1 < 0) || (index1 >= row_num)) {
	printf ("Error in TIsajet::GetBQDP : \n");
	printf ("Invalid array row index %d; range is 0-%d.\n", index1, row_num-1);
	return 0;
    }

    if ((index2 < 0) || (index2 >= col_num)) {
	printf ("Error in TIsajet::GetBQDP : \n");
	printf ("Invalid array column index %d; range is 0-%d.\n", index2, col_num-1);
	return 0;
    }

    return WCON.bqdp[index1][index2];
}

/**************************************************************************/

Float_t TIsajet::GetEZDP() const
{
    return WCON.ezdp;
}

/**************************************************************************/

void TIsajet::SetWFUDGE(Float_t val)
{
    WCON.wfudge = val;
}

/**************************************************************************/

Float_t TIsajet::GetWFUDGE() const
{
    return WCON.wfudge;
}

/**************************************************************************/

// Ends WCON access.

#ifndef WIN32
# define isaini  isaini_
# define isaevt  isaevt_
# define isabeg  isabeg_
# define isabg2  isabg2_
# define openfiles openfiles_
# define pdfinit pdfinit_
# define ranf    ranf_
# define type_of_call
#else
# define isaini  ISAINI
# define isaevt  ISAEVT
# define isabeg  ISABEG
# define isabg2  ISABG2
# define openfiles OPENFILES
# define pdfinit PDFINIT
# define ranf    RANF
# define type_of_call _stdcall
#endif

extern "C" void type_of_call isaini(Int_t& j, Int_t& k, Int_t& m, Int_t& n);
extern "C" void type_of_call isaevt(Int_t& j, Int_t& k, Int_t& m);
extern "C" void type_of_call openfiles();
extern "C" void type_of_call pdfinit(Int_t &, Int_t &);
extern "C" void type_of_call isabeg(Int_t& ifl);
extern "C" void type_of_call isabg2(Int_t& ifl);

void TIsajet::Isaini(Int_t& j, Int_t& k, Int_t& m, Int_t& n) 
{
    isaini(j, k, m, n);
}

void TIsajet::Isaevt(Int_t& j, Int_t& k, Int_t& m) 
{
    isaevt(j, k, m);
}

void TIsajet::Openfiles() 
{
    openfiles();
}

void TIsajet::PDFinit(Int_t& pdfpar, Int_t& pdfval ) 
{
    pdfinit(pdfpar, pdfval);
}

void TIsajet::Isabeg(Int_t& ifl) 
{
    isabeg(ifl);
}

void TIsajet::Isabg2(Int_t& ifl) 
{
    isabg2(ifl);
}

extern "C" {
    Double_t type_of_call ranf(Int_t & /*idum*/) 
    {
	Float_t r;
	do r=sRandom->Rndm(); while(0 >= r || r >= 1);
	return r;
    }
}







    



