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

#include <stdlib.h>

#include <TROOT.h> 
#include <TArrayI.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TParticleClassPDG.h>
#include <TPDGCode.h> 

#include "AliIonPDGCodes.h" 

ClassImp(AliIonPDGCodes)
 
//______________________________________________________________________________
    AliIonPDGCodes::AliIonPDGCodes():
	TObject(),
	fNIon(200)
{ 
    for(Int_t i=0; i<fNIon; i++) fPDGCode[i]=0;
}

AliIonPDGCodes::AliIonPDGCodes(const AliIonPDGCodes &/*PDGCodes*/)
    :TObject(),
     fNIon(200)
{
  for(Int_t i=0; i<fNIon; i++) fPDGCode[i]=0;
}

AliIonPDGCodes&  AliIonPDGCodes::operator=(const AliIonPDGCodes &pdg)
{
  for(Int_t i=0; i<fNIon; i++) fPDGCode[i]=0;
  if (this != &pdg) {
      TObject::operator=(pdg);
      fNIon = pdg.fNIon;
        for(Int_t i=0; i<fNIon; i++) fPDGCode[i]=pdg.fPDGCode[i];
  }
  return *this;
}

//______________________________________________________________________________
void AliIonPDGCodes::AddParticlesToPdgDataBase()
{

  // --- PDG codes for ions are added to the PDGDatabase
  // --- PDGCode = Z*10000 + A + Offset(to avoid multiple definitions)
  
  TDatabasePDG *pdgDB = TDatabasePDG::Instance();
  
  const Int_t kOffset=10000000;
  
  pdgDB->AddParticle("H","H",1.00794,kTRUE,0,1,"Ion",1*10000+1*10+kOffset,-1,8);
  //pdgDB->AddParticle("Helium","Helium",4.0026,kTRUE,0,2,"Ion",2*10000+4*10+kOffset,-1,8);
  pdgDB->AddParticle("Li7","Li7",6.941,kTRUE,0,3,"Ion",3*10000+7*10+kOffset,-1,8);
  pdgDB->AddParticle("Be","Be",9.01218,kTRUE,0,4,"Ion",4*10000+9*10+kOffset,-1,8);
  pdgDB->AddParticle("B","B",10.811,kTRUE,0,5,"Ion",5*10000+11*10+kOffset,-1,8);
  pdgDB->AddParticle("C12","C12",12.0107,kTRUE,0,6,"Ion",6*10000+12*10+kOffset,-1,8);
  pdgDB->AddParticle("N","N",14.00674,kTRUE,0,7,"Ion",7*10000+14*10+kOffset,-1,8);
  pdgDB->AddParticle("O","O",15.9994,kTRUE,0,8,"Ion",8*10000+16*10+kOffset,-1,8);
  pdgDB->AddParticle("F","F",18.9984,kTRUE,0,9,"Ion",9*10000+19*10+kOffset,-1,8);
  pdgDB->AddParticle("Ne","Ne",20.1797,kTRUE,0,10,"Ion",10*10000+20*10+kOffset,-1,8);
  pdgDB->AddParticle("Na","Na",22.98977,kTRUE,0,11,"Ion",11*10000+23*10+kOffset,-1,8);
  pdgDB->AddParticle("Mg","Mg",24.3050,kTRUE,0,12,"Ion",12*10000+24*10+kOffset,-1,8);
  pdgDB->AddParticle("Al","Al",26.9815,kTRUE,0,13,"Ion",13*10000+27*10+kOffset,-1,8);
  pdgDB->AddParticle("Si","Si",28.0855,kTRUE,0,14,"Ion",14*10000+28*10+kOffset,-1,8);
  pdgDB->AddParticle("P","P",30.97376,kTRUE,0,15,"Ion",15*10000+31*10+kOffset,-1,8);
  pdgDB->AddParticle("S","S",32.066,kTRUE,0,16,"Ion",16*10000+32*10+kOffset,-1,8);
  pdgDB->AddParticle("Cl","Cl",35.4527,kTRUE,0,17,"Ion",17*10000+35*10+kOffset,-1,8);
  pdgDB->AddParticle("Ar","Ar",39.948,kTRUE,0,18,"Ion",18*10000+40*10+kOffset,-1,8);
  pdgDB->AddParticle("K","K",39.0983,kTRUE,0,19,"Ion",19*10000+39*10+kOffset,-1,8);
  pdgDB->AddParticle("Ca","Ca",40.078,kTRUE,0,20,"Ion",20*10000+40*10+kOffset,-1,8);
  pdgDB->AddParticle("Sc","Sc",44.95591,kTRUE,0,21,"Ion",21*10000+45*10+kOffset,-1,8);
  pdgDB->AddParticle("Ti","Ti",47.867,kTRUE,0,22,"Ion",22*10000+48*10+kOffset,-1,8);
  pdgDB->AddParticle("V","V",50.9415,kTRUE,0,23,"Ion",23*10000+51*10+kOffset,-1,8);
  pdgDB->AddParticle("Cr","Cr",51.9961,kTRUE,0,24,"Ion",24*10000+52*10+kOffset,-1,8);
  pdgDB->AddParticle("Mn","Mn",54.938,kTRUE,0,25,"Ion",25*10000+55*10+kOffset,-1,8);
  pdgDB->AddParticle("Fe","Fe",55.845,kTRUE,0,26,"Ion",26*10000+56*10+kOffset,-1,8);
  pdgDB->AddParticle("Co","Co",58.933200,kTRUE,0,27,"Ion",27*10000+59*10+kOffset,-1,8);
  pdgDB->AddParticle("Ni","Ni",58.6934,kTRUE,0,28,"Ion",28*10000+59*10+kOffset,-1,8);
  pdgDB->AddParticle("Cu","Cu",63.546,kTRUE,0,29,"Ion",29*10000+63*10+kOffset,-1,8);
  pdgDB->AddParticle("Zn","Zn",65.39,kTRUE,0,30,"Ion",30*10000+65*10+kOffset,-1,8);
  pdgDB->AddParticle("Ga","Ga",69.723,kTRUE,0,31,"Ion",31*10000+70*10+kOffset,-1,8);
  pdgDB->AddParticle("Ge","Ge",72.61,kTRUE,0,32,"Ion",32*10000+73*10+kOffset,-1,8);
  pdgDB->AddParticle("As","As",74.9216,kTRUE,0,33,"Ion",33*10000+75*10+kOffset,-1,8);
  pdgDB->AddParticle("Se","Se",78.96,kTRUE,0,34,"Ion",34*10000+79*10+kOffset,-1,8);
  pdgDB->AddParticle("Br","Br",79.904,kTRUE,0,35,"Ion",35*10000+80*10+kOffset,-1,8);
  pdgDB->AddParticle("Kr","Kr",83.80,kTRUE,0,36,"Ion",36*10000+84*10+kOffset,-1,8);
  pdgDB->AddParticle("Rb","Rb",85.4678,kTRUE,0,37,"Ion",37*10000+85*10+kOffset,-1,8);
  pdgDB->AddParticle("Sr","Sr",87.62,kTRUE,0,38,"Ion",38*10000+88*10+kOffset,-1,8);
  pdgDB->AddParticle("Y","Y",88.90585,kTRUE,0,39,"Ion",39*10000+89*10+kOffset,-1,8);
  pdgDB->AddParticle("Zr","Zr",91.224,kTRUE,0,40,"Ion",40*10000+91*10+kOffset,-1,8);
  pdgDB->AddParticle("Nb","Nb",92.90638,kTRUE,0,41,"Ion",41*10000+93*10+kOffset,-1,8);
  pdgDB->AddParticle("Mo","Mo",95.94,kTRUE,0,42,"Ion",42*10000+96*10+kOffset,-1,8);
  pdgDB->AddParticle("Tc","Tc",97.907,kTRUE,0,43,"Ion",43*10000+98*10+kOffset,-1,8);
  pdgDB->AddParticle("Ru","Rut",101.07,kTRUE,0,44,"Ion",44*10000+101*10+kOffset,-1,8);
  pdgDB->AddParticle("Rh","Rh",102.9055,kTRUE,0,45,"Ion",45*10000+103*10+kOffset,-1,8);
  pdgDB->AddParticle("Pd","Pd",106.42,kTRUE,0,46,"Ion",46*10000+106*10+kOffset,-1,8);
  pdgDB->AddParticle("Ag","Ag",107.8682,kTRUE,0,47,"Ion",47*10000+108*10+kOffset,-1,8);
  pdgDB->AddParticle("Cd","Cd",112.411,kTRUE,0,48,"Ion",48*10000+112*10+kOffset,-1,8);
  pdgDB->AddParticle("In","In",114.818,kTRUE,0,49,"Ion",49*10000+115*10+kOffset,-1,8);
  pdgDB->AddParticle("Sn","Sn",118.710,kTRUE,0,50,"Ion",50*10000+119*10+kOffset,-1,8);
  pdgDB->AddParticle("Sb","Sb",121.760,kTRUE,0,51,"Ion",51*10000+122*10+kOffset,-1,8);
  pdgDB->AddParticle("Te","Te",127.60,kTRUE,0,52,"Ion",52*10000+128*10+kOffset,-1,8);
  pdgDB->AddParticle("I","I",126.90447,kTRUE,0,53,"Ion",53*10000+127*10+kOffset,-1,8);
  pdgDB->AddParticle("Xe","Xe",131.29,kTRUE,0,54,"Ion",54*10000+131*10+kOffset,-1,8);
  pdgDB->AddParticle("Cs","Cs",132.90545,kTRUE,0,55,"Ion",55*10000+133*10+kOffset,-1,8);
  pdgDB->AddParticle("Ba","Ba",137.327,kTRUE,0,56,"Ion",56*10000+137*10+kOffset,-1,8);
  pdgDB->AddParticle("La","La",138.906,kTRUE,0,57,"Ion",57*10000+139*10+kOffset,-1,8);
  pdgDB->AddParticle("Ce","Ce",140.116,kTRUE,0,58,"Ion",58*10000+140*10+kOffset,-1,8);
  pdgDB->AddParticle("Pr","Pr",140.908,kTRUE,0,59,"Ion",59*10000+141*10+kOffset,-1,8);
  pdgDB->AddParticle("Nd","Nd",144.240,kTRUE,0,60,"Ion",60*10000+144*10+kOffset,-1,8);
  pdgDB->AddParticle("Pm","Pm",144.913,kTRUE,0,61,"Ion",61*10000+145*10+kOffset,-1,8);
  pdgDB->AddParticle("Sm","Sm",150.360,kTRUE,0,62,"Ion",62*10000+150*10+kOffset,-1,8);
  pdgDB->AddParticle("Eu","Eu",151.964,kTRUE,0,63,"Ion",63*10000+152*10+kOffset,-1,8);
  pdgDB->AddParticle("Gd","Gd",157.250,kTRUE,0,64,"Ion",64*10000+157*10+kOffset,-1,8);
  pdgDB->AddParticle("Tb","Tb",158.925,kTRUE,0,65,"Ion",65*10000+159*10+kOffset,-1,8);
  pdgDB->AddParticle("Dy","Dy",162.500,kTRUE,0,66,"Ion",66*10000+162*10+kOffset,-1,8);
  pdgDB->AddParticle("Ho","Ho",164.930,kTRUE,0,67,"Ion",67*10000+165*10+kOffset,-1,8);
  pdgDB->AddParticle("Er","Er",167.260,kTRUE,0,68,"Ion",68*10000+167*10+kOffset,-1,8);
  pdgDB->AddParticle("Tm","Tm",168.934,kTRUE,0,69,"Ion",69*10000+169*10+kOffset,-1,8);
  pdgDB->AddParticle("Yb","Yb",173.040,kTRUE,0,70,"Ion",70*10000+173*10+kOffset,-1,8);
  pdgDB->AddParticle("Lu","Lu",174.967,kTRUE,0,71,"Ion",71*10000+175*10+kOffset,-1,8);
  pdgDB->AddParticle("Hf","Hf",178.49,kTRUE,0,72,"Ion",72*10000+178*10+kOffset,-1,8);
  pdgDB->AddParticle("Ta","Ta",180.9479,kTRUE,0,73,"Ion",73*10000+181*10+kOffset,-1,8);
  pdgDB->AddParticle("W","W",183.84,kTRUE,0,74,"Ion",74*10000+184*10+kOffset,-1,8);
  pdgDB->AddParticle("Re","Re",186.207,kTRUE,0,75,"Ion",75*10000+186*10+kOffset,-1,8);
  pdgDB->AddParticle("Os","Os",190.23,kTRUE,0,76,"Ion",76*10000+190*10+kOffset,-1,8);
  pdgDB->AddParticle("Ir","Ir",192.217,kTRUE,0,77,"Ion",77*10000+192*10+kOffset,-1,8);
  pdgDB->AddParticle("Pt","Pt",195.078,kTRUE,0,78,"Ion",78*10000+195*10+kOffset,-1,8);
  pdgDB->AddParticle("Au","Au",196.96655,kTRUE,0,79,"Ion",79*10000+197*10+kOffset,-1,8);
  pdgDB->AddParticle("Hg","Hg",200.59,kTRUE,0,80,"Ion",80*10000+200*10+kOffset,-1,8);
  pdgDB->AddParticle("Tl","Tl",204.38,kTRUE,0,81,"Ion",81*10000+204*10+kOffset,-1,8);
  pdgDB->AddParticle("Pb","Pb",207.2,kTRUE,0,82,"Ion",82*10000+207*10+kOffset,-1,8);
  pdgDB->AddParticle("Bi","Bi",208.98038,kTRUE,0,83,"Ion",83*10000+209*10+kOffset,-1,8);
  // --- Ions added to run DPMJET (Chiara)
  pdgDB->AddParticle("Mg25","Mg25",25.,kTRUE,0,12,"Ion",12*10000+25*10+kOffset,-1,8);
  pdgDB->AddParticle("Co58","Co58",58.,kTRUE,0,27,"Ion",27*10000+58*10+kOffset,-1,8);
  pdgDB->AddParticle("Ti46","Ti46",46.,kTRUE,0,22,"Ion",22*10000+46*10+kOffset,-1,8);
  pdgDB->AddParticle("Ni58","Ni58",58.,kTRUE,0,28,"Ion",28*10000+58*10+kOffset,-1,8);
  pdgDB->AddParticle("Zn30","Zn30",67.,kTRUE,0,30,"Ion",30*10000+67*10+kOffset,-1,8);
  pdgDB->AddParticle("Se74","Se74",74.,kTRUE,0,34,"Ion",34*10000+74*10+kOffset,-1,8);
  pdgDB->AddParticle("Sr84","Sr84",84.,kTRUE,0,38,"Ion",38*10000+84*10+kOffset,-1,8);
  pdgDB->AddParticle("Sb147","Sb147",147.,kTRUE,0,51,"Ion",51*10000+147*10+kOffset,-1,8);
  pdgDB->AddParticle("I131","I131",131.,kTRUE,0,53,"Ion",53*10000+131*10+kOffset,-1,8);
  pdgDB->AddParticle("Pr132","Pr132",132.,kTRUE,0,59,"Ion",59*10000+132*10+kOffset,-1,8);
  pdgDB->AddParticle("Nd134","Nd134",134.,kTRUE,0,60,"Ion",60*10000+134*10+kOffset,-1,8);
  pdgDB->AddParticle("Dy144","Dy144",144.,kTRUE,0,66,"Ion",66*10000+144*10+kOffset,-1,8);
  pdgDB->AddParticle("Tm156","Tm156",156.,kTRUE,0,69,"Ion",69*10000+156*10+kOffset,-1,8);
  pdgDB->AddParticle("Yb177","Yb177",177.,kTRUE,0,70,"Ion",70*10000+177*10+kOffset,-1,8);
  pdgDB->AddParticle("Ir173","Ir173",173.,kTRUE,0,77,"Ion",77*10000+173*10+kOffset,-1,8);
  pdgDB->AddParticle("Os185","Os185",185.,kTRUE,0,76,"Ion",76*10000+185*10+kOffset,-1,8);
  
  printf("\n\n	AliIonPDGCodes -> definition of PDG codes for ions\n\n");

}

//_____________________________________________________________________________
Int_t AliIonPDGCodes::IdFromPDG(Int_t pdg) const 
{
  //
  // Return Geant3 code from PDG and pseudo ENDF code
  //
  const Int_t kOffset=10000000;
  for(Int_t i=0; i<fNIon; ++i)
    if(pdg==fPDGCode[i]) return i+kOffset;
  return -1;
}

//_____________________________________________________________________________
Int_t AliIonPDGCodes::PDGFromId(Int_t id) const 
{
  //
  // Return PDG code and pseudo ENDF code from Geant3 code
  //
  const Int_t kOffset=10000000;
  id -= kOffset;
  if(id >= 0 && id < fNIon) return fPDGCode[id];
  else return -1;
}

