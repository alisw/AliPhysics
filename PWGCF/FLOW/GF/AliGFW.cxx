/*
Author: Vytautas Vislavicius
Extention of Generic Flow (https://arxiv.org/abs/1312.3572 by A. Bilandzic et al.)
Class steers the initialization and calculation of n-particle correlations. Uses recursive function, all terms are calculated only once.
Latest version includes the calculation of any number of gaps and any combination of harmonics (including eg symmetric cumulants, etc.)
If used, modified, or distributed, please aknowledge the author of this code.
*/

#include "AliGFW.h"
AliGFW::AliGFW():
  fInitialized(kFALSE)
{
};

AliGFW::~AliGFW() {
  for(auto pItr = fCumulants.begin(); pItr != fCumulants.end(); ++pItr)
    pItr->DestroyComplexVectorArray();
};

void AliGFW::AddRegion(string refName, Int_t lNhar, Int_t lNpar, Double_t lEtaMin, Double_t lEtaMax, Int_t lNpT, Int_t BitMask) {
  if(lNpT < 1) {
    printf("Number of pT bins cannot be less than 1! Not adding anything.\n");
    return;
  };
  if(lEtaMin >= lEtaMax) {
    printf("Eta min. cannot be more than eta max! Not adding...\n");
    return;
  };
  if(refName=="") {
    printf("Region must have a name!\n");
    return;
  };
  Region lOneRegion;
  lOneRegion.Nhar = lNhar; //Number of harmonics
  lOneRegion.Npar = lNpar; //Number of powers
  lOneRegion.NparVec = vector<Int_t> {}; //if powers defined, then set this to empty vector
  lOneRegion.EtaMin = lEtaMin; //Min. eta
  lOneRegion.EtaMax = lEtaMax; //Max. eta
  lOneRegion.NpT = lNpT; //Number of pT bins
  lOneRegion.rName = refName; //Name of the region
  lOneRegion.BitMask = BitMask; //Bit mask
  AddRegion(lOneRegion);
};
void AliGFW::AddRegion(string refName, Int_t lNhar, Int_t *lNparVec, Double_t lEtaMin, Double_t lEtaMax, Int_t lNpT, Int_t BitMask) {
  if(lNpT < 1) {
    printf("Number of pT bins cannot be less than 1! Not adding anything.\n");
    return;
  };
  if(lEtaMin >= lEtaMax) {
    printf("Eta min. cannot be more than eta max! Not adding...\n");
    return;
  };
  if(refName=="") {
    printf("Region must have a name!\n");
    return;
  };
  Region lOneRegion;
  lOneRegion.Nhar = lNhar; //Number of harmonics
  lOneRegion.Npar = 0; //If vector with powers defined, set this to zero
  lOneRegion.NparVec = vector<Int_t> {};//lNparVec; //vector with powers for each harmonic
  for(Int_t i=0;i<lNhar;i++) lOneRegion.NparVec.push_back(lNparVec[i]);
  lOneRegion.EtaMin = lEtaMin; //Min. eta
  lOneRegion.EtaMax = lEtaMax; //Max. eta
  lOneRegion.NpT = lNpT; //Number of pT bins
  lOneRegion.rName = refName; //Name of the region
  lOneRegion.BitMask = BitMask; //Bit mask
  AddRegion(lOneRegion);
};
Int_t AliGFW::CreateRegions() {
  if(fRegions.size()<1) {
    printf("No regions set. Skipping...\n");
    return 0;
  };
  Int_t nRegions=0;
  for(auto pItr=fRegions.begin(); pItr!=fRegions.end(); pItr++) {
    AliGFWCumulant *lCumulant = new AliGFWCumulant();
    if(pItr->NparVec.size()) {
      lCumulant->CreateComplexVectorArrayVarPower(pItr->Nhar, pItr->NparVec, pItr->NpT);
    } else {
      lCumulant->CreateComplexVectorArray(pItr->Nhar, pItr->Npar, pItr->NpT);
    };
    fCumulants.push_back(*lCumulant);
    ++nRegions;
  };
  if(nRegions) fInitialized=kTRUE;
  return nRegions;
};
void AliGFW::Fill(Double_t eta, Int_t ptin, Double_t phi, Double_t weight, Int_t mask, Double_t SecondWeight) {
  if(!fInitialized) CreateRegions();
  if(!fInitialized) return;
  for(Int_t i=0;i<(Int_t)fRegions.size();++i) {
    if(fRegions.at(i).EtaMin<eta && fRegions.at(i).EtaMax>eta && (fRegions.at(i).BitMask&mask))
      fCumulants.at(i).FillArray(ptin,phi,weight,SecondWeight);
  };
};
complex<Double_t> AliGFW::TwoRec(Int_t n1, Int_t n2, Int_t p1, Int_t p2, Int_t ptbin, AliGFWCumulant *r1, AliGFWCumulant *r2, AliGFWCumulant *r3) {
  complex<Double_t> part1 = r1->Vec(n1,p1,ptbin);
  complex<Double_t> part2 = r2->Vec(n2,p2,ptbin);
  complex<Double_t> part3 = r3?r3->Vec(n1+n2,p1+p2,ptbin):complex<Double_t>(0.,0.);
  complex<Double_t> formula = part1*part2-part3;
  return formula;
};
complex<Double_t> AliGFW::RecursiveCorr(AliGFWCumulant *qpoi, AliGFWCumulant *qref, AliGFWCumulant *qol, Int_t ptbin, vector<Int_t> &hars) {
  vector<Int_t> pows;
  for(Int_t i=0; i<(Int_t)hars.size(); i++)
    pows.push_back(1);
  return RecursiveCorr(qpoi, qref, qol, ptbin, hars, pows);
};

complex<Double_t> AliGFW::RecursiveCorr(AliGFWCumulant *qpoi, AliGFWCumulant *qref, AliGFWCumulant *qol, Int_t ptbin, vector<Int_t> &hars, vector<Int_t> &pows) {
  if((pows.at(0)!=1) && qol) qpoi=qol; //if the power of POI is not unity, then always use overlap (if defined).
  //Only valid for 1 particle of interest though!
  if(hars.size()<2) return qpoi->Vec(hars.at(0),pows.at(0),ptbin);
  if(hars.size()<3) return TwoRec(hars.at(0), hars.at(1),pows.at(0),pows.at(1), ptbin, qpoi, qref, qol);
  Int_t harlast=hars.at(hars.size()-1);
  Int_t powlast=pows.at(pows.size()-1);
  hars.erase(hars.end()-1);
  pows.erase(pows.end()-1);
  complex<Double_t> formula = RecursiveCorr(qpoi, qref, qol, ptbin, hars, pows)*qref->Vec(harlast,powlast);
  Int_t lDegeneracy=1;
  Int_t harSize = (Int_t)hars.size();
  for(Int_t i=harSize-1;i>=0;i--) {
  //checking if current configuration is a permutation of the next one.
  //Need to have more than 2 harmonics though, otherwise it doesn't make sense.
    if(i>2) { //only makes sense when we have more than two harmonics remaining
      if(hars.at(i) == hars.at(i-1) && pows.at(i) == pows.at(i-1)) {//if it is a permutation, then increase degeneracy and continue;
        lDegeneracy++;
        continue;
      };
    }
    hars.at(i)+=harlast;
    pows.at(i)+=powlast;
    //The issue is here. In principle, if i=0 (dif), then the overlap is only qpoi (0, if no overlap);
    //Otherwise, if we are not working with the 1st entry (dif.), then overlap will always be from qref
    //One should thus (probably) make a check if i=0, then qovl=qpoi, otherwise qovl=qref. But need to think more
    //-- This is not aplicable anymore, since the overlap is explicitly specified
    complex<Double_t> subtractVal = RecursiveCorr(qpoi, qref, qol, ptbin, hars, pows);
    if(lDegeneracy>1) { subtractVal *= lDegeneracy; lDegeneracy=1; };
    formula-=subtractVal;
    hars.at(i)-=harlast;
    pows.at(i)-=powlast;

  };
  hars.push_back(harlast);
  pows.push_back(powlast);
  return formula;
};
void AliGFW::Clear() {
  for(auto ptr = fCumulants.begin(); ptr!=fCumulants.end(); ++ptr) ptr->ResetQs();
};
AliGFW::CorrConfig AliGFW::GetCorrelatorConfig(string config, string head, Bool_t ptdif) {
  //First remove all ; and ,:
  s_replace_all(config,","," ");
  s_replace_all(config,";"," ");
  s_replace_all(config,"| ","|");
  //If pT-bin is provided, then look for & remove space before "(" (so that it's clean afterwards)
  while(s_index(config," (")>-1) s_replace_all(config," (","(");
  //Then make sure we don't have any double-spaces:
  while(s_index(config,"  ")>-1) s_replace_all(config,"  "," ");
  vector<Int_t> regs;
  vector<Int_t> hars;
  Int_t sz1=0;
  Int_t szend=0;
  string ts, ts2;
  CorrConfig ReturnConfig;
  //Fetch region descriptor
  if(!s_tokenize(config,ts,szend,"{")) {
    printf("Could not find any harmonics!\n");
    return ReturnConfig;
  };
  szend=0;
  Int_t counter=0;
  while(s_tokenize(config,ts,szend,"{")) {
    counter++;
    ReturnConfig.Regs.push_back(vector<Int_t> {});
    ReturnConfig.Hars.push_back(vector<Int_t> {});
    // ReturnConfig.ptInd.push_back(vector<Int_t> {});
    ReturnConfig.Overlap.push_back(-1); //initially, assume no overlap
    //Check if there's a particular pT bin I should be using here. If so, store it (otherwise, it's bin 0)
    Int_t ptbin=-1;
    Int_t sz2=0;
    if(s_contains(ts,"(")) {
      if(!s_contains(ts,")")) {printf("Missing \")\" in the configurator. Returning...\n"); return ReturnConfig; };
      sz2 = s_index(ts,"(");
      sz1=sz2+1;
      s_tokenize(ts,ts2,sz1,")");
      ptbin=stoi(ts2);
      ts.erase(sz2,(sz1-sz2+1));
      szend-=(sz1-sz2); //szend also becomes shorter
      //also need to remove this from config now:
      sz2 = s_index(config,"(");
      sz1 = s_index(config,")");
      config.erase(sz2,sz1-sz2+1);
    };
    ReturnConfig.ptInd.push_back(ptbin);
    sz1=0;
    //Fetch regions
    while(s_tokenize(ts,ts2,sz1," ")) {
      if(sz1>=szend) break;
      Bool_t isOverlap=s_contains(ts2,"|");
      if(isOverlap) ts2.erase(0,1); //If overlap, remove the delimiter |
      Int_t ind=FindRegionByName(ts2);
      if(ts2 == " " || ts2 == "") continue;
      if(ind<0) {
        printf("Could not find region named %s!\n",ts2.c_str());
        break;
      };
      if(!isOverlap)
        ReturnConfig.Regs.at(counter-1).push_back(ind);
      else ReturnConfig.Overlap.at((int)ReturnConfig.Overlap.size()-1) = ind;
    };
    string harstr;
    s_tokenize(config,harstr,szend,"}");
    Int_t dummys=0;
    //Fetch harmonics
    while(s_tokenize(harstr,ts,dummys," ")) ReturnConfig.Hars.at(counter-1).push_back(stoi(ts));
  };
  ReturnConfig.Head = head;
  ReturnConfig.pTDif = ptdif;
  // ReturnConfig.pTbin = ptbin;
  return ReturnConfig;
};

complex<Double_t> AliGFW::Calculate(Int_t poi, Int_t ref, vector<Int_t> hars, Int_t ptbin) {
  AliGFWCumulant *qref = &fCumulants.at(ref);
  AliGFWCumulant *qpoi = &fCumulants.at(poi);
  AliGFWCumulant *qovl = qpoi;
  return RecursiveCorr(qpoi, qref, qovl, ptbin, hars);
};
// complex<Double_t> AliGFW::Calculate(CorrConfig corconf, Int_t ptbin, Bool_t SetHarmsToZero, Bool_t DisableOverlap) {
//    vector<Int_t> ptbins;
//    for(Int_t i=0;i<(Int_t)corconf.size();i++) ptbins.push_back(ptbin);
//    return Calculate(corconf,ptbins,SetHarmsToZero,DisableOverlap);
// }
complex<Double_t> AliGFW::Calculate(CorrConfig corconf, Int_t ptbin, Bool_t SetHarmsToZero, Bool_t DisableOverlap) {
  if(corconf.Regs.size()==0) return complex<Double_t>(0,0); //Check if we have any regions at all
  // if(ptbins.size()!=corconf.Regs.size()) {printf("Number of pT-bins is not the same as number of subevents!\n"); return complex<Double_t>(0,0); };
  complex<Double_t> retval(1,0);
  Int_t ptInd;
  for(Int_t i=0;i<(Int_t)corconf.Regs.size();i++) { //looping over all regions
    if(corconf.Regs.at(i).size()==0)  return complex<Double_t>(0,0); //again, if no regions in the current subevent, then quit immediatelly
    ptInd = corconf.ptInd.at(i); //for i=0 (potentially, POI)
    if(ptInd<0) ptInd = ptbin;
    // Int_t ptbin = ptbins.at(i);
    //picking up the indecies of regions...
    Int_t poi = corconf.Regs.at(i).at(0);
    Int_t ref = (corconf.Regs.at(i).size()>1)?corconf.Regs.at(i).at(1):corconf.Regs.at(i).at(0);
    Int_t ovl = corconf.Overlap.at(i);
    //and regions themselves
    AliGFWCumulant *qref = &fCumulants.at(ref);
    AliGFWCumulant *qpoi = &fCumulants.at(poi);
    if(!qref->IsPtBinFilled(ptInd)) return complex<Double_t>(0,0); //if REF is not filled, don't even continue. Could be redundant, but should save little CPU time
    if(!qpoi->IsPtBinFilled(ptInd)) return complex<Double_t>(0,0);//if POI is not filled, don't even continue. Could be redundant, but should save little CPU time
    AliGFWCumulant *qovl=0;
    //Check if in the ref. region we have enough particles (no. of particles in the region >= no of harmonics for subevent)
    Int_t sz1 = corconf.Hars.at(i).size();
    if(poi!=ref) sz1--;
    if(qref->GetN() < sz1) return complex<Double_t>(0,0);
    //Then, figure the overlap
    if(ovl > -1) //if overlap is defined, then (unless it's explicitly disabled)
      qovl = DisableOverlap?0:&fCumulants.at(ovl);
    else if(ref==poi) qovl = qref; //If ref and poi are the same, then the same is for overlap. Only, when OL not explicitly defined
    if(SetHarmsToZero) for(Int_t j=0;j<(Int_t)corconf.Hars.at(i).size();j++) corconf.Hars.at(i).at(j) = 0;
    retval *= RecursiveCorr(qpoi, qref, qovl, ptInd, corconf.Hars.at(i));
  }
  return retval;
};

complex<Double_t> AliGFW::Calculate(Int_t poi, vector<Int_t> hars) {
  AliGFWCumulant *qpoi = &fCumulants.at(poi);
  return RecursiveCorr(qpoi, qpoi, qpoi, 0, hars);
};
Int_t AliGFW::FindRegionByName(string refName) {
  for(Int_t i=0;i<(Int_t)fRegions.size();i++) if(fRegions.at(i).rName == refName) return i;
  return -1;
};
//String processing:
Int_t AliGFW::s_index(string &instr, const string &pattern, const Int_t &spos) {
  return instr.find(pattern,spos);
};
Bool_t AliGFW::s_contains(string &instr, const string &pattern) {
  return (s_index(instr,pattern)>-1);
};
void AliGFW::s_replace(string &instr, const string &pattern1, const string &pattern2, const Int_t &spos) {
  Int_t lpos = s_index(instr,pattern1,spos);
  if(lpos<0) return;
  instr.replace(lpos,pattern1.size(),pattern2);
};
void AliGFW::s_replace_all(string &instr, const string &pattern1, const string &pattern2) {
  Int_t lpos=s_index(instr,pattern1);
  while(lpos>-1) { s_replace(instr,pattern1,pattern2,lpos); lpos=s_index(instr,pattern1,lpos); };
};
Bool_t AliGFW::s_tokenize(string &instr, string &subs, Int_t &spos, const string &delim) {
  if(spos<0 || spos>=(Int_t)instr.size()) {spos=-1; subs=""; return kFALSE;};
  Int_t lpos = s_index(instr,delim,spos);
  if(lpos<0) lpos=instr.size();
  subs = instr.substr(spos,lpos-spos);
  spos=lpos+1;
  return kTRUE;
}
