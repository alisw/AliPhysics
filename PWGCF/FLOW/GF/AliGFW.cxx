/*
Author: Vytautas Vislavicius
Extention of Generic Flow (https://arxiv.org/abs/1312.3572)
*/
#include "AliGFW.h"
/*TODOs:
need to add flags to have control over what is added, e.g. what happens, when I have several overlapping regions of different types: reference, pT-diff unID and pT-diff. ID?
*/
AliGFW::AliGFW():
  fInitialized(kFALSE)
{
};

AliGFW::~AliGFW() {
  for(auto pItr = fCumulants.begin(); pItr != fCumulants.end(); ++pItr)
    pItr->DestroyComplexVectorArray();
};

void AliGFW::AddRegion(TString refName, Int_t lNhar, Int_t lNpar, Double_t lEtaMin, Double_t lEtaMax, Int_t lNpT, Int_t BitMask) {
  if(lNpT < 1) {
    printf("Number of pT bins cannot be less than 1! Not adding anything.\n");
    return;
  };
  if(lEtaMin >= lEtaMax) {
    printf("Eta min. cannot be more than eta max! Not adding...\n");
    return;
  };
  if(refName.EqualTo("")) {
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
void AliGFW::AddRegion(TString refName, Int_t lNhar, Int_t *lNparVec, Double_t lEtaMin, Double_t lEtaMax, Int_t lNpT, Int_t BitMask) {
  if(lNpT < 1) {
    printf("Number of pT bins cannot be less than 1! Not adding anything.\n");
    return;
  };
  if(lEtaMin >= lEtaMax) {
    printf("Eta min. cannot be more than eta max! Not adding...\n");
    return;
  };
  if(refName.EqualTo("")) {
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
void AliGFW::SplitRegions() {
  //Simple case. Will not look for overlaps, etc. Everything is left for end-used
};

Int_t AliGFW::CreateRegions() {
  if(fRegions.size()<1) {
    printf("No regions set. Skipping...\n");
    return 0;
  };
  SplitRegions();
  //for(auto pitr = fRegions.begin(); pitr!=fRegions.end(); pitr++) pitr->PrintStructure();
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
      fCumulants.at(i).FillArray(eta,ptin,phi,weight,SecondWeight);
  };
};
TComplex AliGFW::TwoRec(Int_t n1, Int_t n2, Int_t p1, Int_t p2, Int_t ptbin, AliGFWCumulant *r1, AliGFWCumulant *r2, AliGFWCumulant *r3) {
  TComplex part1 = r1->Vec(n1,p1,ptbin);
  TComplex part2 = r2->Vec(n2,p2,ptbin);
  TComplex part3 = r3?r3->Vec(n1+n2,p1+p2,ptbin):TComplex(0,0);
  TComplex formula = part1*part2-part3;
  return formula;
};
TComplex AliGFW::RecursiveCorr(AliGFWCumulant *qpoi, AliGFWCumulant *qref, AliGFWCumulant *qol, Int_t ptbin, vector<Int_t> hars, vector<Int_t> pows) {
  if(pows.size()==0) //if powers are not initialized, initialize them to 1
    for(Int_t i=0; i<(Int_t)hars.size(); i++)
      pows.push_back(1);
  if(hars.size()<2) return qpoi->Vec(hars.at(0),pows.at(0),ptbin);
  if(hars.size()<3) return TwoRec(hars.at(0), hars.at(1),pows.at(0),pows.at(1), ptbin, qpoi, qref, qol);
  Int_t harlast=hars.at(hars.size()-1);
  Int_t powlast=pows.at(pows.size()-1);
  hars.erase(hars.end()-1);
  pows.erase(pows.end()-1);
  TComplex formula = RecursiveCorr(qpoi, qref, qol, ptbin, hars, pows)*qref->Vec(harlast,powlast);
  for(Int_t i=0;i<(Int_t)hars.size();i++) {
    vector<Int_t> lhars = hars;
    vector<Int_t> lpows = pows;
    lhars.at(i)+=harlast;
    lpows.at(i)+=powlast;
    //The issue is here. In principle, if i=0 (dif), then the overlap is only qpoi (0, if no overlap);
    //Otherwise, if we are not working with the 1st entry (dif.), then overlap will always be from qref
    //One should thus (probably) make a check if i=0, then qovl=qpoi, otherwise qovl=qref. But need to think more
    formula-=RecursiveCorr(qpoi, qref, qol, ptbin, lhars, lpows);
  };
  return formula;
};
void AliGFW::Clear() {
  for(auto ptr = fCumulants.begin(); ptr!=fCumulants.end(); ++ptr) ptr->ResetQs();
  fCalculatedNames.clear();
  fCalculatedQs.clear();
};
TComplex AliGFW::Calculate(TString config, Bool_t SetHarmsToZero) {
  if(config.EqualTo("")) {
    printf("Configuration empty!\n");
    return TComplex(0,0);
  };
  TString tmp;
  Ssiz_t sz1=0;
  TComplex ret(1,0);
  while(config.Tokenize(tmp,sz1,"}")) {
    if(SetHarmsToZero) SetHarmonicsToZero(tmp);
    /*Int_t ind=FindCalculated(tmp);
    if(ind>=0) {
      ret*=fCalculatedQs.at(ind);
      continue;
    };*/
    TComplex val=CalculateSingle(tmp);
    ret*=val;
    fCalculatedQs.push_back(val);
    fCalculatedNames.push_back(tmp);
  };
  return ret;
};
TComplex AliGFW::CalculateSingle(TString config) {
  //First remove all ; and ,:
  config.ReplaceAll(","," ");
  config.ReplaceAll(";"," ");
  //Then make sure we don't have any double-spaces:
  while(config.Index("  ")>-1) config.ReplaceAll("  "," ");
  vector<Int_t> regs;
  vector<Int_t> hars;
  Int_t ptbin=0;
  Ssiz_t sz1=0;
  Ssiz_t szend=0;
  TString ts, ts2;
  //find the pT-bin:
  if(config.Tokenize(ts,sz1,"(")) {
    config.Tokenize(ts,sz1,")");
    ptbin=ts.Atoi();
  };
  //Fetch region descriptor
  if(sz1<0) sz1=0;
  if(!config.Tokenize(ts,szend,"{")) {
    printf("Could not find harmonics!\n");
    return TComplex(0,0);
  };
  //Fetch regions
  while(ts.Tokenize(ts2,sz1," ")) {
    if(sz1>=szend) break;
    Int_t ind=FindRegionByName(ts2);
    if(ts2.EqualTo(" ") || ts2.EqualTo("")) continue;
    if(ind<0) {
      printf("Could not find region named %s!\n",ts2.Data());
      break;
    };
    regs.push_back(ind);
  };
  //Fetch harmonics
  while(config.Tokenize(ts,szend," ")) hars.push_back(ts.Atoi());
  if(regs.size()==1) return Calculate(regs.at(0),hars);
  return Calculate(regs.at(0),regs.at(1),hars,ptbin);
};
AliGFW::CorrConfig AliGFW::GetCorrelatorConfig(TString config, TString head, Bool_t ptdif) {
  //First remove all ; and ,:
  config.ReplaceAll(","," ");
  config.ReplaceAll(";"," ");
  config.ReplaceAll("| ","|");
  //Then make sure we don't have any double-spaces:
  while(config.Index("  ")>-1) config.ReplaceAll("  "," ");
  vector<Int_t> regs;
  vector<Int_t> hars;
  Int_t ptbin=0;
  Ssiz_t sz1=0;
  Ssiz_t szend=0;
  TString ts, ts2;
  CorrConfig ReturnConfig;
  //find the pT-bin:
  if(config.Tokenize(ts,sz1,"(")) {
    config.Tokenize(ts,sz1,")");
    ptbin=ts.Atoi();
  };
  //Fetch region descriptor
  if(sz1<0) sz1=0;

  if(!config.Tokenize(ts,szend,"{")) {
    printf("Could not find harmonics!\n");
    return ReturnConfig;
  };
  szend=0;
  Int_t counter=0;
  while(config.Tokenize(ts,szend,"{")) {
    counter++;
    sz1=0;
    //Fetch regions
    while(ts.Tokenize(ts2,sz1," ")) {
      if(sz1>=szend) break;
      Bool_t isOverlap=ts2.Contains("|");
      if(isOverlap) ts2.Remove(0,1); //If overlap, remove the delimiter |
      Int_t ind=FindRegionByName(ts2);
      if(ts2.EqualTo(" ") || ts2.EqualTo("")) continue;
      if(ind<0) {
        printf("Could not find region named %s!\n",ts2.Data());
        break;
      };
      if(counter==1) {
        if(!isOverlap)
          ReturnConfig.Regs.push_back(ind);
        else ReturnConfig.Overlap1 = ind;
      } else {
        if(isOverlap)
          ReturnConfig.Regs2.push_back(ind);
        else ReturnConfig.Overlap2 = ind;
      }
    };
    TString harstr;
    config.Tokenize(harstr,szend,"}");
    Ssiz_t dummys=0;
    //Fetch harmonics
    if(counter==1)
      while(harstr.Tokenize(ts,dummys," ")) ReturnConfig.Hars.push_back(ts.Atoi());
    else
      while(harstr.Tokenize(ts,dummys," ")) ReturnConfig.Hars2.push_back(ts.Atoi());
  };
  ReturnConfig.Head = head;
  ReturnConfig.pTDif = ptdif;
  return ReturnConfig;
};

TComplex AliGFW::Calculate(Int_t poi, Int_t ref, vector<Int_t> hars, Int_t ptbin) {
  AliGFWCumulant *qref = &fCumulants.at(ref);
  AliGFWCumulant *qpoi = &fCumulants.at(poi);
  AliGFWCumulant *qovl = qpoi;
  return RecursiveCorr(qpoi, qref, qovl, ptbin, hars);
};
TComplex AliGFW::Calculate(CorrConfig corconf, Int_t ptbin, Bool_t SetHarmsToZero, Bool_t DisableOverlap) {
  if(corconf.Regs.size()==0) return TComplex(0,0);
  Int_t poi = corconf.Regs.at(0);
  Int_t ref = (corconf.Regs.size()>1)?corconf.Regs.at(1):corconf.Regs.at(0);
  AliGFWCumulant *qref = &fCumulants.at(ref);
  AliGFWCumulant *qpoi = &fCumulants.at(poi);
  AliGFWCumulant *qovl=0;
  if(corconf.Overlap1 > -1)
    qovl = DisableOverlap?0:(&fCumulants.at(corconf.Overlap1));//;DisableOverlap?0:qpoi;
  else if(ref==poi) qovl = qref; //If ref and poi are the same, then the same is for overlap. Only, when OL not explicitly defined
  if(!qpoi->IsPtBinFilled(ptbin)) return TComplex(0,0);
  //if(!qref->IsPtBinFilled(ptbin)) return TComplex(0,0);
  if(SetHarmsToZero) for(Int_t i=0;i<(Int_t)corconf.Hars.size();i++) corconf.Hars.at(i) = 0;
  TComplex retval = RecursiveCorr(qpoi, qref, qovl, ptbin, corconf.Hars);
  if(corconf.Regs2.size()==0) return retval;
  poi = corconf.Regs2.at(0);
  ref = (corconf.Regs2.size()>1)?corconf.Regs2.at(1):corconf.Regs2.at(0);
  qref = &fCumulants.at(ref);
  qpoi = &fCumulants.at(poi);
  if(corconf.Overlap2 > -1)
    qovl = DisableOverlap?0:(&fCumulants.at(corconf.Overlap2));//;DisableOverlap?0:qpoi;
  else if(ref==poi) qovl = qref; //Only when OL is not explicitly defined, then set it to ref/POI if they are the same
  if(SetHarmsToZero) for(Int_t i=0;i<(Int_t)corconf.Hars2.size();i++) corconf.Hars2.at(i) = 0;
  retval*=RecursiveCorr(qpoi, qref, qovl, 0, corconf.Hars2);
  return retval;
};

TComplex AliGFW::Calculate(Int_t poi, vector<Int_t> hars) {
  AliGFWCumulant *qpoi = &fCumulants.at(poi);
  return RecursiveCorr(qpoi, qpoi, qpoi, 0, hars);
};
Int_t AliGFW::FindRegionByName(TString refName) {
  for(Int_t i=0;i<(Int_t)fRegions.size();i++) if(fRegions.at(i).rName.EqualTo(refName)) return i;
  return -1;
};
Int_t AliGFW::FindCalculated(TString identifier) {
  if(fCalculatedNames.size()==0) return -1;
  for(Int_t i=0;i<(Int_t)fCalculatedNames.size();i++) {
    if(fCalculatedNames.at(i).EqualTo(identifier)) return i;
  };
  return -1;
};
Bool_t AliGFW::SetHarmonicsToZero(TString &instr) {
  TString tmp;
  Ssiz_t sz1=0, sz2;
  if(!instr.Tokenize(tmp,sz1,"{")) {
    printf("AliGFW::SetHarmonicsToZero: could not find \"{\" token in %s\n",instr.Data());
    return kFALSE;
  };
  sz2=sz1;
  Int_t indc=0;
  while(instr.Tokenize(tmp,sz1," ")) indc++;
  if(!indc) {
    printf("AliGFW::SetHarmonicsToZero: did not find any harmonics in %s, nothing to replace\n",instr.Data());
    return kFALSE;
  };
  instr.Remove(sz2);
  for(Int_t i=0;i<indc;i++) instr.Append("0 ");
  return kTRUE;
};
