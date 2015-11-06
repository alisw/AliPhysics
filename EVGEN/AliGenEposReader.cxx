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
//
// Realisation of the AliGenReader interface to be used with AliGenExFile.
// NextEvent() loops over events
// and NextParticle() loops over particles.
// This implementation reads EPOS v3.111 output format (from ROOT trees)
// Author: Igor Lakomov <Igor.Lakomov@cern.ch>

#include "AliGenEposReader.h"
#include "AliGenEpos3EventHeader.h"

using namespace std;

// Arrays for PDG codes
// EPOS codes-->PDG
// a=269
 const Int_t a=269;
 Int_t eposid[a]={1, 2, 3, 4, 5, 6, 10, 9, 12, -12, 11, -11, 14, -14, 13, -13, 16, 15, 110, 120, -120, 220, 130, -130, 230, -230, 20, -20, 330, 111, 121, -121, 221, 131, -131, 231, -231, 331, -140, 240, 1120, 1220, 2130, 1130, 1230, 2230, 1330, 2330, 1111, 1121, 1221, 2221, 1131, 2231, 1331, 2331, 3331, 2140, 17, 18, 19, 0, 99, 1112, 1113, 1114, 2222, 2223, 2224, 2224, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1233, 1234, 1235, 1236, 1237, 1238, 1239, 1132, 1133, 1134, 2232, 2233, 2234, 90, 80, 81, 85, 86, 82, 83, 84, 1200, 2300, 1300, 2400, 1400, 3400, 2500, 1500, 3500, 4500, 2200, 1200, 2300, 1300, 3300, 2400, 1400, 3400, 4400, 2500, 1500, 3500, 4500, 5500, 800000091, 800000092, 800000093, 800000094, -340, 340, -241, 241, -141, 141, -341, 341, 250, 150, 350, 450, 251, 151, 351, 451, 440, 441, 550, 551, 2240, 1240, 1140, 2241, 1241, 3240, 2340, 3140, 1340, 3340, 2341, 1341, 3341, 2440, 2441, 1440, 1441, 3440, 3441, 4441, 2250, 2150, 3250, 4250, 1250, 1150, 3150, 4150, 2350, 1350, 3350, 4350, 2450, 1450, 3450, 4450, 2550, 1550, 3550, 2251, 1151, 2351, 1351, 3351, 2451, 1451, 3451, 4451, 2551, 1551, 3551, 4551, 5551, 123, 122, 233, 232, 133, 132, 143, 132, 243, 242, 343, 342, 223, 222, 113, 112, 333, 332, 443, 442, 444, 253, 252, 153, 152, 353, 352, 453, 452, 553, 552, 124, 125, 234, 235, 134, 135, 144, 135, 244, 245, 344, 345, 224, 225, 114, 115, 334, 335, 444, 445, 254, 255, 154, 155, 354, 355, 454, 455, 554, 555, 11099, 12099, 22099, 33099, 44099, 112099, 122099, 800000110, 800000990};
 Int_t pdgid[a]={2, 1, 3, 4, 5, 6, 22, 21, 11, -11, 12, -12, 13, -13, 14, -14, 15, 16, 111, 211, -211, 221, 321, -321, 311, -311, 310, -310, 331, 113, 213, -213, 223, 323, -323, 313, -313, 333, 421, -411, 2212, 2112, 3122, 3222, 3212, 3112, 3322, 3312, 2224, 2214, 2114, 1114, 3224, 3114, 3324, 3314, 3334, 4122, 99, 99, 99, 99, 99, 32224, 12224, 12222, 31114, 11114, 11112, 21114, 12212, 2124, 32214, 2216, 12214, 22124, 11212, 13122, 3124, 23122, 13212, 23212, 53122, 13216, 13222, 23222, 13226, 13112, 23112, 13116, 23, 24, 25, 32, 33, 35, 36, 37, 2101, 3101, 3201, 4101, 4201, 4301, 5101, 5201, 5301, 5401, 1103, 1103, 3103, 3203, 3303, 4103, 4203, 4303, 4403, 5103, 5203, 5303, 5403, 5503, 91, 92, 93, 94, 431, -431, 413, -413, 423, -423, 433, -433, 511, 521, 531, 541, 513, 523, 533, 543, 441, 443, 551, 553, 4112, 4212, 4222, 4114, 4224, 4132, 4312, 4232, 4322, 4332, 4314, 4324, 4334, 4412, 4414, 4422, 4424, 4432, 4434, 4444, 5112, 5122, 5132, 5142, 5212, 5222, 5232, 5242, 5312, 5322, 5332, 5342, 5412, 5422, 5432, 5442, 5512, 5522, 5532, 5114, 5224, 5314, 5324, 5334, 5414, 5424, 5434, 5444, 5524, 5524, 5534, 5544, 5554, 10213, 10211, 10313, 10311, 10323, 10321, 10423, 10421, 10413, 10411, 10433, 10431, 10113, 10111, 10223, 10221, 10333, 10331, 10443, 10441, 10443, 10513, 10511, 10523, 10521, 10533, 10531, 10543, 10541, 10553, 10551, 20213, 215, 20313, 315, 20323, 325, 20423, 425, 20413, 415, 20433, 435, 20113, 115, 20223, 225, 20333, 335, 20443, 445, 20513, 515, 20523, 525, 20533, 535, 20543, 545, 20553, 555, 9900110, 9900210, 9900220, 9900330, 9900440, 9902210, 9902110, 110, 990};
	
ClassImp(AliGenEposReader)

AliGenEposReader::AliGenEposReader():
  AliGenReader(),
  fNcurrent(0),
  fNparticle(0),
  fTreeNtuple(0),
  fTreeHeader(0),
  fCurrentEvent(0),
  fCurrentParticle(0),
  fGenEventHeader(0),
  fIversn(0),
  fLaproj(0),
  fMaproj(0),
  fLatarg(0),
  fMatarg(0),
  fEngy(0.),
  fNfull(0),
  fNfreeze(0),
  bim(0)
   
{
  //default constructor
}
    
AliGenEposReader::AliGenEposReader(const AliGenEposReader &reader):
  AliGenReader(reader),
  fNcurrent(0),
  fNparticle(0),
  fTreeNtuple(0),
  fTreeHeader(0),
  fCurrentEvent(0),
  fCurrentParticle(0),
  fGenEventHeader(0),
  fIversn(0),
  fLaproj(0),
  fMaproj(0),
  fLatarg(0),
  fMatarg(0),
  fEngy(0.),
  fNfull(0),
  fNfreeze(0),
  bim(0)
  
{
  reader.Copy(*this);
}
    
AliGenEposReader::~AliGenEposReader(){ delete fTreeNtuple; }

void AliGenEposReader::Init()
{
  //reset file and open a new root file
  TFile *file=0;
  if (!file) {
    file = new TFile(fFileName);
    AliInfo(Form("File %s opened", fFileName));
    file->cd();
  }
  else {
    AliError(Form("Couldn't open input file: %s", fFileName));
  }
 
  // Initialization
  fTreeHeader = (TTree*)gDirectory->Get("teposhead");
  fTreeHeader->SetMakeClass(1);
  fTreeHeader->SetBranchAddress("iversn",&fIversn);
  fTreeHeader->SetBranchAddress("laproj",&fLaproj);
  fTreeHeader->SetBranchAddress("maproj",&fMaproj);
  fTreeHeader->SetBranchAddress("latarg",&fLatarg);
  fTreeHeader->SetBranchAddress("matarg",&fMatarg);
  fTreeHeader->SetBranchAddress("engy",&fEngy);
  fTreeHeader->SetBranchAddress("nfull",&fNfull);
  fTreeHeader->SetBranchAddress("nfreeze",&fNfull);
  fTreeHeader->GetEvent(0);

  fTreeNtuple = (TTree*)gDirectory->Get("teposevent");
  fTreeNtuple->SetMakeClass(1);
  const Int_t maxnp=fTreeNtuple->GetMaximum("np");
  zus.resize(maxnp);
  px.resize(maxnp);
  py.resize(maxnp);
  pz.resize(maxnp);
  e.resize(maxnp);
  x.resize(maxnp);
  y.resize(maxnp);
  z.resize(maxnp);
  t.resize(maxnp);
  id.resize(maxnp);
  ist.resize(maxnp);
  ity.resize(maxnp);
  ior.resize(maxnp);
  jor.resize(maxnp);
  fTreeNtuple->SetBranchAddress("np",&np);
  fTreeNtuple->SetBranchAddress("bim",&bim);
  fTreeNtuple->SetBranchAddress("zus",&zus[0]);
  fTreeNtuple->SetBranchAddress("px",&px[0]);
  fTreeNtuple->SetBranchAddress("py",&py[0]);
  fTreeNtuple->SetBranchAddress("pz",&pz[0]);
  fTreeNtuple->SetBranchAddress("e",&e[0]);
  fTreeNtuple->SetBranchAddress("x",&x[0]);
  fTreeNtuple->SetBranchAddress("y",&y[0]);
  fTreeNtuple->SetBranchAddress("z",&z[0]);
  fTreeNtuple->SetBranchAddress("t",&t[0]);
  fTreeNtuple->SetBranchAddress("id",&id[0]);
  fTreeNtuple->SetBranchAddress("ist",&ist[0]);
  fTreeNtuple->SetBranchAddress("ity",&ity[0]);
  fTreeNtuple->SetBranchAddress("ior",&ior[0]);
  fTreeNtuple->SetBranchAddress("jor",&jor[0]);
}

Int_t  AliGenEposReader::NextEvent()
{	
  fCurrentParticle = 0;
  fMothersMap.clear();
  Int_t nentries = (Int_t) fTreeNtuple->GetEntries();
  if(fCurrentEvent < nentries){
	fTreeNtuple->GetEvent(fCurrentEvent);
        if (np>(Int_t)zus.size()) Init();
        fGenEventHeader = new AliGenEpos3EventHeader();
        ((AliGenEpos3EventHeader*)fGenEventHeader)->SetIversn(fIversn);
        ((AliGenEpos3EventHeader*)fGenEventHeader)->SetLaproj(fLaproj);
        ((AliGenEpos3EventHeader*)fGenEventHeader)->SetMaproj(fMaproj);
        ((AliGenEpos3EventHeader*)fGenEventHeader)->SetLatarg(fLatarg);
        ((AliGenEpos3EventHeader*)fGenEventHeader)->SetMatarg(fMatarg);
        ((AliGenEpos3EventHeader*)fGenEventHeader)->SetEngy(fEngy);
        ((AliGenEpos3EventHeader*)fGenEventHeader)->SetNfull(fNfull);
        ((AliGenEpos3EventHeader*)fGenEventHeader)->SetNfreeze(fNfull);
        ((AliGenEpos3EventHeader*)fGenEventHeader)->SetBim(bim);
        for (Int_t i=0; i<np; i++) fMothersMap.insert( std::pair<Int_t,Int_t>((zus[i]>0&&ior[i]!=i&&ior[i]>=0)?ior[i]:-1,i));
	fCurrentEvent++;
	return np;
  }	
  AliError("No more events in the file.");
  return 0;
}

TParticle*  AliGenEposReader::NextParticle()
{
  std::pair <std::multimap<Int_t,Int_t>::iterator, std::multimap<Int_t,Int_t>::iterator> ret = fMothersMap.equal_range(fCurrentParticle);
  TParticle* particle =
     new TParticle(EposToPdg(id[fCurrentParticle]),
		   (ist[fCurrentParticle]==0) ? 1 : 0,
		   (zus[fCurrentParticle]>0&&ior[fCurrentParticle]!=fCurrentParticle) ? ior[fCurrentParticle] : -1,
		   -1,
		   ret.first->second,
		   ret.first->second+fMothersMap.count(fCurrentParticle)-1,
		   px[fCurrentParticle],
		   py[fCurrentParticle],
		   pz[fCurrentParticle],
		   e[fCurrentParticle],
		   x[fCurrentParticle],
		   y[fCurrentParticle],
		   z[fCurrentParticle],
		   t[fCurrentParticle]);
    
  if (particle && particle->GetStatusCode()==1) {
      particle->SetBit(kTransportBit);
  }
  fCurrentParticle++;
  return particle;
}

Int_t AliGenEposReader::EposToPdg(Int_t code)
{
  for(Int_t i=0;i<a;i++) {
    if (eposid[i]==code) return pdgid[i];
  }
//  cout<<"Particle id cannot be find in PDG! return 0; \n";
  return 0;	
}

void  AliGenEposReader::RewindEvent()
{
  fCurrentParticle = 0; 
}

AliGenEposReader&  AliGenEposReader::operator=(const  AliGenEposReader& rhs)
{
  rhs.Copy(*this);
  return(*this);
}

void AliGenEposReader::Copy(TObject&) const
{
  Fatal("Copy","Not implemented!\n");
}
