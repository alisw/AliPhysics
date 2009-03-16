/*
 *###################################################################
 *#        EPOS 1.67     K. WERNER, T. PIEROG, S. PORTEBOEUF.       #
 *#                      Contact: werner@subatech.in2p3.fr          #
 *###################################################################
 *
 * TEpos.cxx
 * 
 * Wraper class for interfacing EPOS model, derived from ROOT's TGenerator.
 * It generates temporary input file for the model, providing user with
 * ability to add his/hers own lines to the input.
 * Output is read directly from common blocks.
 *
 *      Author: Piotr Ostrowski, postrow@if.pw.edu.pl
 */


#include <TClonesArray.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TRandom.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "eposproc.h"
#include "EPOScommon.h"
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
		fPhimax(2*3.1415927),
		fEcms(-1),
		fSplitting(kFALSE),
		fNoDecays(),
		fExtraInputLines()
{
}

TEpos::~TEpos() {}

void TEpos::Initialize() {
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
}

void TEpos::GenerateEvent() {
	cseed.seedj = gRandom->Rndm() * 1e10;
	aseed(2);
	Int_t n = 1;
	evgen(n);
}

Int_t TEpos::ImportParticles(TClonesArray *particles, Option_t *) {
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
		TDatabasePDG *pdgDb = TDatabasePDG::Instance();
		Int_t status = cptl.istptl[iPart] + 1;
		Int_t pdg = pdgDb->ConvertIsajetToPdg(cptl.idptl[iPart]);
		if (pdg == 0) {
			printf("TEpos: Warning, unknown particle, index: %d, ISAJET: %d\n",iPart,cptl.idptl[iPart]);
		}
    		new((*particles)[iPart]) TParticle(pdg, status,
				 tFather, -1, -1, -1,
				 cptl.pptl[iPart][0], cptl.pptl[iPart][1],cptl.pptl[iPart][2],cptl.pptl[iPart][3],
				 cptl.xorptl[iPart][0]*1.e-12, cptl.xorptl[iPart][1]*1.e-12, cptl.xorptl[iPart][2]*1.e-12, cptl.xorptl[iPart][3]*1e-12);
    		(*particles)[iPart]->SetUniqueID(iPart);
  	}

  	return numpart;
}

TObjArray*  TEpos::ImportParticles(Option_t *) {
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
		TDatabasePDG *pdgDb = TDatabasePDG::Instance();
		Int_t status = cptl.istptl[iPart] + 1;
		Int_t pdg = pdgDb->ConvertIsajetToPdg(cptl.idptl[iPart]);
		if (pdg == 0) {
			printf("TEpos: Warning, unknown particle, index: %d, ISAJET: %d\n",iPart,cptl.idptl[iPart]);
		}
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
	ofstream file(GetInputFileName(), ios_base::out | ios_base::trunc);
	char epo[256];
	strcpy(epo, getenv("ALICE_ROOT"));
	strcat(epo, "/EPOS/epos167");

	file << "fname pathnx " << epo << "/" << endl;
	file << "fname histo none" << endl;
	file << "fname copy none" << endl;
	file << "fname log none" << endl;
	file << "fname check none" << endl;
	file << "fname data /tmp/epos.out" << endl;
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

    file << "output epos" << endl;
    file << "record event nevt nptl b endrecord" << endl;
    file << "record particle i id fa mo c1 c2 st endrecord" << endl;

	file << "input " << epo << "/epos.param" << endl;
	file << "runprogram" << endl;
	file.close();
}

Float_t TEpos::GetBimevt() { return cevt.bimevt; }
Float_t TEpos::GetPhievt() { return cevt.phievt; }
Int_t TEpos::GetKolevt() { return cevt.kolevt; }
Int_t TEpos::GetKoievt() { return cevt.koievt; }
Float_t TEpos::GetPmxevt() { return cevt.pmxevt; }
Float_t TEpos::GetEgyevt() { return cevt.egyevt; }
Int_t TEpos::GetNpjevt() { return cevt.npjevt; }
Int_t TEpos::GetNtgevt() { return cevt.ntgevt; }
Int_t TEpos::GetNpnevt() { return cevt.npnevt; }
Int_t TEpos::GetNppevt() { return cevt.nppevt; }
Int_t TEpos::GetNtnevt() { return cevt.ntnevt; }
Int_t TEpos::GetNtpevt() { return cevt.ntpevt; }
Int_t TEpos::GetJpnevt() { return cevt.jpnevt; }
Int_t TEpos::GetJppevt() { return cevt.jppevt; }
Int_t TEpos::GetJtnevt() { return cevt.jtnevt; }
Int_t TEpos::GetJtpevt() { return cevt.jtpevt; }
Float_t TEpos::GetXbjevt() { return cevt.xbjevt; }
Float_t TEpos::GetQsqevt() { return cevt.qsqevt; }
Int_t TEpos::GetNglevt() { return cevt.nglevt; }
Float_t TEpos::GetZppevt() { return cevt.zppevt; }
Float_t TEpos::GetZptevt() { return cevt.zptevt; }

