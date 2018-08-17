#include "TInfo.h"
#include <TDatime.h>
#include <TFile.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH2F.h>
#include <TMap.h>
#include <TNtuple.h>

Float_t TInfo::AbsMinT(Int_t type) const
{
  Float_t min=1e9;
  for (Int_t i=0;i<NSensors();++i) {
    if (!IsValid(i))
      continue;
    Float_t val = T(i,type);
    if (val<min)
      min = val;
  }
  return min;
}

Float_t TInfo::AbsMaxT(Int_t type) const
{
  Float_t max=-1e9;
  for (Int_t i=0;i<NSensors();++i) {
    if (!IsValid(i))
      continue;
    Float_t val = T(i,type);
    if (val>max)
      max = val;
  }
  return max;
}

Float_t TInfo::AvgTempSM(Int_t sm, Int_t type) const
{
  Int_t nActiveSensors = 8;
  Double_t temp = 0;
  for (Int_t i=sm*8;i<(sm+1)*8;++i) {
    if (T(i,type) <= 15 || T(i,type) > 27 ) // take out sensors which are nonesense
      nActiveSensors--;
    else
      temp += T(i,3);
  }
  if (nActiveSensors != 0)
    temp /= nActiveSensors;
  return temp;
}

TH2 *TInfo::GetHist(Int_t type) const
{
  TH2F *h = new TH2F(Form("h%dT%d_%s",fRunNo,type,Type(type)),";col;rows",8,-0.5,7.5,20,-0.5,19.5);
  h->SetTitle(Form("%s for run %d (frac=%.1f%%, length=%.1fh)",Type(type),fRunNo,Fraction()*100,(fLastTime-fFirstTime)/3600.));
  h->SetDirectory(0);
  h->SetStats(0);
  Double_t min=AbsMinT(type);
  Double_t max=AbsMaxT(type);
  h->SetMinimum(min);
  h->SetMaximum(max);
  for (Int_t i=0;i<NSensors();++i) {
    if (!IsValid(i))
      continue;
    Double_t val = T(i,type);
    Double_t bin = GetBin(i);
    h->SetBinContent(bin,val);
  }
  return h;
}

void TInfo::Print(Option_t *option) const
{
  cout << "Runno: " << fRunNo << " with average time " << fAvTime << " and " << Nvalid() << " entries" << endl;
  for (Int_t i=0;i<NSensors();++i) {
    if (IsValid(i))
      cout << "  " << i << " avgT=" << fAvgT.At(i) << ", rmsT=" << fRmsT.At(i) << ", minT=" << fMinT.At(i) << ", maxT=" << fMaxT.At(i) << ", diff=" << Diff(i) << endl;
  }
}

Float_t TInfo::T(Int_t ns, Int_t type) const
{
  Double_t val = 0;
  if (type==1)
    val = MinT(ns);
  else if (type==2)
    val = MaxT(ns);
  else if (type==3)
    val = AvgT(ns);
  else if (type==4)
    val = RmsT(ns);
  else if (type==5)
    val = (MinT(ns)+MaxT(ns))/2;
  else
    val = Diff(ns);
  return val;
}

/*static*/
Int_t TInfo::GetBin(Int_t ns)
{
  static TH2F *h=0;
  if (h==0) {
    h = new TH2F("hdummy",";col;rows",8,-0.5,7.5,20,-0.5,19.5);
    h->SetDirectory(0);
  }

  switch (ns) {
  case  0: return h->GetBin(1,19); break;
  case  1: return h->GetBin(1,20); break;
  case  2: return h->GetBin(2,19); break;
  case  3: return h->GetBin(2,20); break;
  case  4: return h->GetBin(3,19); break;
  case  5: return h->GetBin(3,20); break;
  case  6: return h->GetBin(4,19); break;
  case  7: return h->GetBin(4,20); break;
  case  8: return h->GetBin(8,20); break;
  case  9: return h->GetBin(8,19); break;
  case 10: return h->GetBin(7,20); break;
  case 11: return h->GetBin(7,19); break;
  case 12: return h->GetBin(6,20); break;
  case 13: return h->GetBin(6,19); break;
  case 14: return h->GetBin(5,20); break;
  case 15: return h->GetBin(5,19); break;

  case 16: return h->GetBin(1,17); break;
  case 17: return h->GetBin(1,18); break;
  case 18: return h->GetBin(2,17); break;
  case 19: return h->GetBin(2,18); break;
  case 20: return h->GetBin(3,17); break;
  case 21: return h->GetBin(3,18); break;
  case 22: return h->GetBin(4,17); break;
  case 23: return h->GetBin(4,18); break;
  case 24: return h->GetBin(8,18); break;
  case 25: return h->GetBin(8,17); break;
  case 26: return h->GetBin(7,18); break;
  case 27: return h->GetBin(7,17); break;
  case 28: return h->GetBin(6,18); break;
  case 29: return h->GetBin(6,17); break;
  case 30: return h->GetBin(5,18); break;
  case 31: return h->GetBin(5,17); break;

  case 32: return h->GetBin(1,15); break;
  case 33: return h->GetBin(1,16); break;
  case 34: return h->GetBin(2,15); break;
  case 35: return h->GetBin(2,16); break;
  case 36: return h->GetBin(3,15); break;
  case 37: return h->GetBin(3,16); break;
  case 38: return h->GetBin(4,15); break;
  case 39: return h->GetBin(4,16); break;
  case 40: return h->GetBin(8,16); break;
  case 41: return h->GetBin(8,15); break;
  case 42: return h->GetBin(7,16); break;
  case 43: return h->GetBin(7,15); break;
  case 44: return h->GetBin(6,16); break;
  case 45: return h->GetBin(6,15); break;
  case 46: return h->GetBin(5,16); break;
  case 47: return h->GetBin(5,15); break;

  case 48: return h->GetBin(1,13); break;
  case 49: return h->GetBin(1,14); break;
  case 50: return h->GetBin(2,13); break;
  case 51: return h->GetBin(2,14); break;
  case 52: return h->GetBin(3,13); break;
  case 53: return h->GetBin(3,14); break;
  case 54: return h->GetBin(4,13); break;
  case 55: return h->GetBin(4,14); break;
  case 56: return h->GetBin(8,14); break;
  case 57: return h->GetBin(8,13); break;
  case 58: return h->GetBin(7,14); break;
  case 59: return h->GetBin(7,13); break;
  case 60: return h->GetBin(6,14); break;
  case 61: return h->GetBin(6,13); break;
  case 62: return h->GetBin(5,14); break;
  case 63: return h->GetBin(5,13); break;

  case 64: return h->GetBin(1,11); break;
  case 65: return h->GetBin(1,12); break;
  case 66: return h->GetBin(2,11); break;
  case 67: return h->GetBin(2,12); break;
  case 68: return h->GetBin(3,11); break;
  case 69: return h->GetBin(3,12); break;
  case 70: return h->GetBin(4,11); break;
  case 71: return h->GetBin(4,12); break;
  case 72: return h->GetBin(8,12); break;
  case 73: return h->GetBin(8,11); break;
  case 74: return h->GetBin(7,12); break;
  case 75: return h->GetBin(7,11); break;
  case 76: return h->GetBin(6,12); break;
  case 77: return h->GetBin(6,11); break;
  case 78: return h->GetBin(5,12); break;
  case 79: return h->GetBin(5,11); break;

  case 80: return h->GetBin(1,9);  break;
  case 81: return h->GetBin(1,10); break;
  case 82: return h->GetBin(2,9);  break;
  case 83: return h->GetBin(2,10); break;
  case 84: return h->GetBin(3,9); break;
  case 85: return h->GetBin(3,10); break;
  case 86: return h->GetBin(4,9); break;
  case 87: return h->GetBin(4,10); break;
  case 88: return h->GetBin(8,10); break;
  case 89: return h->GetBin(8,9);  break;
  case 90: return h->GetBin(7,10); break;
  case 91: return h->GetBin(7,9);  break;
  case 92: return h->GetBin(6,10); break;
  case 93: return h->GetBin(6,9);  break;
  case 94: return h->GetBin(5,10); break;
  case 95: return h->GetBin(5,9);  break;

  case  96: return h->GetBin(1,7); break;
  case  97: return h->GetBin(1,8); break;
  case  98: return h->GetBin(2,7); break;
  case  99: return h->GetBin(2,8); break;
  case 100: return h->GetBin(3,7); break;
  case 101: return h->GetBin(3,8); break;
  case 102: return h->GetBin(4,7); break;
  case 103: return h->GetBin(4,8); break;
  case 104: return h->GetBin(8,8); break;
  case 105: return h->GetBin(8,7); break;
  case 106: return h->GetBin(7,8); break;
  case 107: return h->GetBin(7,7); break;
  case 108: return h->GetBin(6,8); break;
  case 109: return h->GetBin(6,7); break;
  case 110: return h->GetBin(5,8); break;
  case 111: return h->GetBin(5,7); break;

  case 112: return h->GetBin(1,5); break;
  case 113: return h->GetBin(1,6); break;
  case 114: return h->GetBin(2,5); break;
  case 115: return h->GetBin(2,6); break;
  case 116: return h->GetBin(3,5); break;
  case 117: return h->GetBin(3,6); break;
  case 118: return h->GetBin(4,5); break;
  case 119: return h->GetBin(4,6); break;
  case 120: return h->GetBin(8,6); break;
  case 121: return h->GetBin(8,5); break;
  case 122: return h->GetBin(7,6); break;
  case 123: return h->GetBin(7,5); break;
  case 124: return h->GetBin(6,6); break;
  case 125: return h->GetBin(6,5); break;
  case 126: return h->GetBin(5,6); break;
  case 127: return h->GetBin(5,5); break;

  case 128: return h->GetBin(1,3); break;
  case 129: return h->GetBin(1,4); break;
  case 130: return h->GetBin(2,3); break;
  case 131: return h->GetBin(2,4); break;
  case 132: return h->GetBin(3,3); break;
  case 133: return h->GetBin(3,4); break;
  case 134: return h->GetBin(4,3); break;
  case 135: return h->GetBin(4,4); break;
  case 136: return h->GetBin(8,4); break;
  case 137: return h->GetBin(8,3); break;
  case 138: return h->GetBin(7,4); break;
  case 139: return h->GetBin(7,3); break;
  case 140: return h->GetBin(6,4); break;
  case 141: return h->GetBin(6,3); break;
  case 142: return h->GetBin(5,4); break;
  case 143: return h->GetBin(5,3); break;

  case 144: return h->GetBin(1,1); break;
  case 145: return h->GetBin(1,2); break;
  case 146: return h->GetBin(2,1); break;
  case 147: return h->GetBin(2,2); break;
  case 148: return h->GetBin(3,1); break;
  case 149: return h->GetBin(3,2); break;
  case 150: return h->GetBin(4,1); break;
  case 151: return h->GetBin(4,2); break;
  case 152: return h->GetBin(8,2); break;
  case 153: return h->GetBin(8,1); break;
  case 154: return h->GetBin(7,2); break;
  case 155: return h->GetBin(7,1); break;
  case 156: return h->GetBin(6,2); break;
  case 157: return h->GetBin(6,1); break;
  case 158: return h->GetBin(5,2); break;
  case 159: return h->GetBin(5,1); break;
  }
  return 0;
}

Int_t TInfo::SensId(Int_t sm, Int_t row, Int_t col)
{
  Int_t ret = 0;

  Int_t nrows=12;
  Int_t ncols=12;
  if (sm>11 && sm<18){
    ncols=16;
  }
  if (sm==10||sm==11||sm==18||sm==19) {
    ncols=24;
    nrows=4;
  }

  Bool_t oddsm = sm%2;
  if (sm==11) {
    ret = 91;
  } else if (sm==19) {
    ret = 155;
  } else if (!oddsm) {
    ret = sm*8;
  } else {
    ret = (sm+1)*8-1;
  }

  if (!oddsm) {
    Int_t x = col/ncols;
    Int_t b = 0;
    if (row<nrows)
      b = 1;
    ret += 2*x + b;
  } else {
    Int_t x = col/ncols;
    Int_t b = 0;
    if (row<nrows)
      b = 1;
    ret -= 2*x + b;
  }

  return ret;
}

const char *TInfo::Type(Int_t type)
{
  TString title("MaxT-MinT");
  if (type==1)
    title="MinT";
  else if (type==2)
    title="MaxT";
  else if (type==3)
    title="AvgT";
  else if (type==4)
    title="RmsT";
  else if (type==5)
    title="AvgMinMaxT";
  return Form("%s",title.Data());
}
