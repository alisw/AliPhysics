//                                                                            
//                                                                            
//        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
//      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
//                           November. 2, 2005                                
//
//

#ifndef UKUTILITY_H
#define UKUTILITY_H
 
class TLorentzVector;
class TVector3;

void IsotropicR3(Double_t r, Double_t *pX, Double_t *pY, Double_t *pZ);
void IsotropicR3(Double_t r, TVector3 &pos);
void MomAntiMom(TLorentzVector &mom, Double_t mass, TLorentzVector &antiMom,
		Double_t antiMass, Double_t initialMass);

const Double_t kGeV = 1.;
const Double_t kFermi = 1.;
const Double_t kHbarc = 0.197 * kGeV * kFermi; 

#endif
