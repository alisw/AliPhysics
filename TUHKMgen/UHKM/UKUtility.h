#ifndef UKUTILITY_INCLUDED
#define UKUTILITY_INCLUDED
 
/*                                                                            
                                                                            
        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
                           November. 2, 2005                                

*/

class TLorentzVector;
class TVector3;
class TH1F;

class Particle;

void IsotropicR3(Double_t r, Double_t *pX, Double_t *pY, Double_t *pZ);
void IsotropicR3(Double_t r, TVector3 &pos);
void MomAntiMom(TLorentzVector &mom, Double_t mass, TLorentzVector &antiMom,
		Double_t antiMass, Double_t initialMass);

extern const Double_t GeV;
extern const Double_t MeV;
extern const Double_t fermi;
extern const Double_t mbarn;
extern const Double_t hbarc; 
extern const Double_t w;
extern const Double_t hbarc_squared; 

#endif
