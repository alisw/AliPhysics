//
//###################################################################
//#        EPOS 1.67     K. WERNER, T. PIEROG, S. PORTEBOEUF.       #
//#                      Contact: werner@subatech.in2p3.fr          #
//###################################################################
//
// TEpos.cxx
// 
// Wraper class for interfacing EPOS model, derived from ROOT's TGenerator.
// It generates temporary input file for the model, providing user with
// ability to add his/hers own lines to the input.
// Output is read directly from common blocks.
//
//      Author: Piotr Ostrowski, postrow@if.pw.edu.pl
//


#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "eposproc.h"
#include "EPOScommon.h"
#include "AliGenEposIsajetToPdgConverter.h"
#include "TEpos.h"

using namespace std;

ClassImp(TEpos)

TEpos::TEpos()  : TGenerator("EPOS", "Epos event generator"),
		fLaproj(0),
		fMaproj(0),
		fLatarg(0),
		fMatarg(0),
		fBminim(0.0),
		fBmaxim(10000.0),
		fPhimin(0.0),
		fPhimax(TMath::TwoPi()),
		fEcms(-1),
		fSplitting(kFALSE),
		fNoDecays(),
		fExtraInputLines(),
		fIdConverter()
{
	fIdConverter=new AliGenEposIsajetToPdgConverter();
}

TEpos::TEpos(const TEpos&)  : TGenerator("EPOS", "Epos event generator"),
		fLaproj(0),
		fMaproj(0),
		fLatarg(0),
		fMatarg(0),
		fBminim(0.0),
		fBmaxim(10000.0),
		fPhimin(0.0),
		fPhimax(TMath::TwoPi()),
		fEcms(-1),
		fSplitting(kFALSE),
		fNoDecays(),
		fExtraInputLines(),
		fIdConverter()
{
	fIdConverter = new AliGenEposIsajetToPdgConverter();
}

TEpos::~TEpos() {
	delete fIdConverter;
}

TEpos& TEpos::operator=(const TEpos& epos) {
  //operator=
	if (!fIdConverter && (this != &epos)) {
		fIdConverter = new AliGenEposIsajetToPdgConverter();
	}
	return *this;
}

void TEpos::Initialize() {
  // Generates input file and prepares EPOS to read from it.
	Int_t nopeno = 0;
	GenerateInputFile();
	aaset(0);
	atitle();
	xiniall();

	const char *inputFileName = GetInputFileName();
	Int_t nameLength = strlen(inputFileName);
	setinp(inputFileName, nameLength, nameLength);
	aread();
	while(copen.nopen == -1) {
		copen.nopen=nopeno;
	    prnt1.iecho=1;
		xiniall();
	    aread();
	}
	Int_t ishini;
	utpri("aamain",prnt1.ish,ishini,4,6);
	if(appli.model != 1)
		IniModel(appli.model);
	ebin.nrebin = 1;
	ainit();
	aseed(2);
}

void TEpos::GenerateEvent() {
//	cseed.seedj = gRandom->Rndm() * 1e10;
	Int_t n = 1;
	evgen(n);
}

Int_t TEpos::ImportParticles(TClonesArray *particles, Option_t *) {
  //Fills provided ClonesArray with generated particles
	particles->Clear();
	if (!cevt.nevt) return 0;
	Int_t numpart = cptl.nptl;
	printf("%d particles generated\n", numpart);
	for (Int_t iPart=0; iPart<numpart; iPart++) {
		Int_t tFather = cptl.iorptl[iPart] - 1;
		tFather = tFather < -1 ? -1 : tFather;
		if (tFather> -1) {
			TParticle *mother = (TParticle*) (particles->UncheckedAt(tFather));
      			mother->SetLastDaughter(iPart);
      			if (mother->GetFirstDaughter()==-1)
 				mother->SetFirstDaughter(iPart);
    		}
		Int_t status = cptl.istptl[iPart] + 1;
		Int_t pdg = fIdConverter->ConvertIsajetToPdg(cptl.idptl[iPart]);
    		new((*particles)[iPart]) TParticle(pdg, status,
				 tFather, -1, -1, -1,
				 cptl.pptl[iPart][0], cptl.pptl[iPart][1],cptl.pptl[iPart][2],cptl.pptl[iPart][3],
				 cptl.xorptl[iPart][0]*1.e-12, cptl.xorptl[iPart][1]*1.e-12, cptl.xorptl[iPart][2]*1.e-12, cptl.xorptl[iPart][3]*1e-12);
    		(*particles)[iPart]->SetUniqueID(iPart);
  	}

  	return numpart;
}

TObjArray*  TEpos::ImportParticles(Option_t *) {
  //Creates new particle array
	fParticles->Clear();
	if (!cevt.nevt) return NULL;
	Int_t numpart = cptl.nptl;
	printf("%d particles generated\n", numpart);
	for (Int_t iPart=0; iPart<numpart; iPart++) {
		Int_t tFather = cptl.iorptl[iPart] - 1;
		tFather = tFather < -1 ? -1 : tFather;
		if (tFather> -1) {
			TParticle *mother = (TParticle*) (fParticles->UncheckedAt(tFather));
      			mother->SetLastDaughter(iPart);
      			if (mother->GetFirstDaughter()==-1)
 				mother->SetFirstDaughter(iPart);
    		}
		Int_t status = cptl.istptl[iPart] + 1;
		Int_t pdg = fIdConverter->ConvertIsajetToPdg(cptl.idptl[iPart]);
    		TParticle* p = new TParticle(pdg, status,
				 tFather, -1, -1, -1,
				 cptl.pptl[iPart][0], cptl.pptl[iPart][1],cptl.pptl[iPart][2],cptl.pptl[iPart][3],
				 cptl.xorptl[iPart][0]*1.e-12, cptl.xorptl[iPart][1]*1.e-12, cptl.xorptl[iPart][2]*1.e-12, cptl.xorptl[iPart][3]*1e-12);
    		p->SetUniqueID(iPart);
    		fParticles->Add(p);
  	}

  	return fParticles;
}

void TEpos::AddNoDecay(Int_t nodecay) {
	fNoDecays.push_back(nodecay);
}

void TEpos::AddExtraInputLine(const char *line) {
	fExtraInputLines.push_back(line);
}

void TEpos::GenerateInputFile() {
  // Generate input file in EPOS format
	ofstream file(GetInputFileName(), ios_base::out | ios_base::trunc);
	char epo[256];
	char *epoEnv = getenv("EPO");
	if (epoEnv) {
		strncpy(epo, epoEnv, 255);
	} else {
		strncpy(epo, getenv("ALICE_ROOT"), 255);
	}
	strncat(epo, "/EPOS/epos167", 255);

	file << "fname pathnx " << epo << "/" << endl;
	file << "fname histo none" << endl;
	file << "fname copy none" << endl;
	file << "fname log none" << endl;
	file << "fname check none" << endl;
//	file << "fname data /tmp/epos.out" << endl;
	file << "fname initl " << epo << "/epos.initl" << endl;
	file << "fname inidi " << epo << "/epos.inidi" << endl;
	file << "fname inidr " << epo << "/epos.inidr" << endl;
	file << "fname iniev " << epo << "/epos.iniev" << endl;
	file << "fname inirj " << epo << "/epos.inirj" << endl;
	file << "fname inics " << epo << "/epos.inics" << endl;
	file << "fname inigrv " << epo << "/epos.inigrv" << endl;
	file << "fqgsjet dat " << epo << "/qgsjet/qgsjet.dat" << endl;
	file << "fqgsjet ncs " << epo << "/qgsjet/qgsjet.ncs" << endl;
	file << "fqgsjetII dat " << epo << "/qgsjetII/qgsdat-II-03" << endl;
	file << "fqgsjetII ncs " << epo << "/qgsjetII/sectnu-II-03" << endl;
	file << "nodecay 120" << endl;
	file << "nodecay -120" << endl;
	file << "nodecay 130" << endl;
	file << "nodecay -130" << endl;
	file << "nodecay -20" << endl;
	file << "nodecay 20" << endl;
	file << "nodecay 14" << endl;
	file << "nodecay -14" << endl;
	file << "set ndecay 1111110" << endl;
	file << "echo on" << endl;

//	.optns file
	int precision = file.precision();
	file.precision(15);
	file << "set seedj " << (gRandom->Rndm() * 1e14) << endl;
	file.precision(precision);
	file << "application hadron" << endl;
	file << "set laproj " << fLaproj << endl;
	file << "set maproj " << fMaproj << endl;
	file << "set latarg " << fLatarg << endl;
	file << "set matarg " << fMatarg << endl;
	file << "set bminim " << fBminim << endl;
	file << "set bmaxim " << fBmaxim << endl;
	file << "set phimin " << fPhimin << endl;
	file << "set phimax " << fPhimax << endl;
	file << "set ecms " << fEcms << endl;

	for(unsigned int i = 0; i < fNoDecays.size(); ++i) {
		file << "nodecay " << fNoDecays[i] << endl;
	}

	file << "switch splitting " << (fSplitting ? "on" : "off") << endl;

	file << "frame nucleon-nucleon" << endl;

	for(unsigned int i = 0; i < fExtraInputLines.size(); ++i) {
		file << fExtraInputLines[i] << endl;
	}

//    file << "output epos" << endl;
//    file << "record event nevt nptl b endrecord" << endl;
//    file << "record particle i id fa mo c1 c2 st endrecord" << endl;

	file << "input " << epo << "/epos.param" << endl;
	file << "runprogram" << endl;
	file.close();
}

Float_t TEpos::GetBimevt() const { return cevt.bimevt; }
Float_t TEpos::GetPhievt() const { return cevt.phievt; }
Int_t TEpos::GetKolevt() const { return cevt.kolevt; }
Int_t TEpos::GetKoievt() const { return cevt.koievt; }
Float_t TEpos::GetPmxevt() const { return cevt.pmxevt; }
Float_t TEpos::GetEgyevt() const { return cevt.egyevt; }
Int_t TEpos::GetNpjevt() const { return cevt.npjevt; }
Int_t TEpos::GetNtgevt() const { return cevt.ntgevt; }
Int_t TEpos::GetNpnevt() const { return cevt.npnevt; }
Int_t TEpos::GetNppevt() const { return cevt.nppevt; }
Int_t TEpos::GetNtnevt() const { return cevt.ntnevt; }
Int_t TEpos::GetNtpevt() const { return cevt.ntpevt; }
Int_t TEpos::GetJpnevt() const { return cevt.jpnevt; }
Int_t TEpos::GetJppevt() const { return cevt.jppevt; }
Int_t TEpos::GetJtnevt() const { return cevt.jtnevt; }
Int_t TEpos::GetJtpevt() const { return cevt.jtpevt; }
Float_t TEpos::GetXbjevt() const { return cevt.xbjevt; }
Float_t TEpos::GetQsqevt() const { return cevt.qsqevt; }
Int_t TEpos::GetNglevt() const { return cevt.nglevt; }
Float_t TEpos::GetZppevt() const { return cevt.zppevt; }
Float_t TEpos::GetZptevt() const { return cevt.zptevt; }

