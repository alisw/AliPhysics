
//Author:        Anders Strand Vestbo
//Author:        Uli Frankenfeld
//Last Modified: 05.12.2000

#include "AliL3Logging.h"
#include "AliL3Transform.h"
//#include <TFile.h>
#include <math.h>
//___________________________
// AliL3Transform
//
// Transformation class for ALICE TPC.

ClassImp(AliL3Transform);


AliL3Transform::AliL3Transform(){
  //constructor
  Init();
}


AliL3Transform::~AliL3Transform(){
}

void AliL3Transform::Init(){
  //sector:
  fNTimeBins = 512;  //?uli
  fNRowLow = 55;
  fNRowUp = 119;
  fNSectorLow = 36;
  fNSectorUp = 36;
  fNSector = 72;
  fPadPitchWidthLow = 0.400000;
  fPadPitchWidthUp = 0.600000;
  fZWidth = 0.56599998474121093750;
  fZSigma = 0.22880849748219134199;

  //slices:
  fNSlice = 36;
  fNRow = 174;
  fPi = 3.14159265358979323846;
  for(Int_t i=0;i<36;i++){
    fCos[i] = cos(fPi/9*i);
    fSin[i] = sin(fPi/9*i);
  }

  fX[0] = 88.850006103515625;
  fX[1] = 89.600006103515625;
  fX[2] = 90.350006103515625;
  fX[3] = 91.100006103515625;
  fX[4] = 91.850006103515625;
  fX[5] = 92.600006103515625;
  fX[6] = 93.350006103515625;
  fX[7] = 94.100006103515625;
  fX[8] = 94.850006103515625;
  fX[9] = 95.600006103515625;
  fX[10] = 96.350006103515625;
  fX[11] = 97.100006103515625;
  fX[12] = 97.850006103515625;
  fX[13] = 98.600006103515625;
  fX[14] = 99.350006103515625;
  fX[15] = 100.100006103515625;
  fX[16] = 100.850006103515625;
  fX[17] = 101.600006103515625;
  fX[18] = 102.350006103515625;
  fX[19] = 103.100006103515625;
  fX[20] = 103.850006103515625;
  fX[21] = 104.600006103515625;
  fX[22] = 105.350006103515625;
  fX[23] = 106.100006103515625;
  fX[24] = 106.850006103515625;
  fX[25] = 107.600006103515625;
  fX[26] = 108.350006103515625;
  fX[27] = 109.100006103515625;
  fX[28] = 109.850006103515625;
  fX[29] = 110.600006103515625;
  fX[30] = 111.350006103515625;
  fX[31] = 112.100006103515625;
  fX[32] = 112.850006103515625;
  fX[33] = 113.600006103515625;
  fX[34] = 114.350006103515625;
  fX[35] = 115.100006103515625;
  fX[36] = 115.850006103515625;
  fX[37] = 116.600006103515625;
  fX[38] = 117.350006103515625;
  fX[39] = 118.100006103515625;
  fX[40] = 118.850006103515625;
  fX[41] = 119.600006103515625;
  fX[42] = 120.350006103515625;
  fX[43] = 121.100006103515625;
  fX[44] = 121.850006103515625;
  fX[45] = 122.600006103515625;
  fX[46] = 123.350006103515625;
  fX[47] = 124.100006103515625;
  fX[48] = 124.850006103515625;
  fX[49] = 125.600006103515625;
  fX[50] = 126.350006103515625;
  fX[51] = 127.100006103515625;
  fX[52] = 127.850006103515625;
  fX[53] = 128.600006103515625;
  fX[54] = 129.350006103515625;
  fX[55] = 132.574996948242188;
  fX[56] = 133.574996948242188;
  fX[57] = 134.574996948242188;
  fX[58] = 135.574996948242188;
  fX[59] = 136.574996948242188;
  fX[60] = 137.574996948242188;
  fX[61] = 138.574996948242188;
  fX[62] = 139.574996948242188;
  fX[63] = 140.574996948242188;
  fX[64] = 141.574996948242188;
  fX[65] = 142.574996948242188;
  fX[66] = 143.574996948242188;
  fX[67] = 144.574996948242188;
  fX[68] = 145.574996948242188;
  fX[69] = 146.574996948242188;
  fX[70] = 147.574996948242188;
  fX[71] = 148.574996948242188;
  fX[72] = 149.574996948242188;
  fX[73] = 150.574996948242188;
  fX[74] = 151.574996948242188;
  fX[75] = 152.574996948242188;
  fX[76] = 153.574996948242188;
  fX[77] = 154.574996948242188;
  fX[78] = 155.574996948242188;
  fX[79] = 156.574996948242188;
  fX[80] = 157.574996948242188;
  fX[81] = 158.574996948242188;
  fX[82] = 159.574996948242188;
  fX[83] = 160.574996948242188;
  fX[84] = 161.574996948242188;
  fX[85] = 162.574996948242188;
  fX[86] = 163.574996948242188;
  fX[87] = 164.574996948242188;
  fX[88] = 165.574996948242188;
  fX[89] = 166.574996948242188;
  fX[90] = 167.574996948242188;
  fX[91] = 168.574996948242188;
  fX[92] = 169.574996948242188;
  fX[93] = 170.574996948242188;
  fX[94] = 171.574996948242188;
  fX[95] = 172.574996948242188;
  fX[96] = 173.574996948242188;
  fX[97] = 174.574996948242188;
  fX[98] = 175.574996948242188;
  fX[99] = 176.574996948242188;
  fX[100] = 177.574996948242188;
  fX[101] = 178.574996948242188;
  fX[102] = 179.574996948242188;
  fX[103] = 180.574996948242188;
  fX[104] = 181.574996948242188;
  fX[105] = 182.574996948242188;
  fX[106] = 183.574996948242188;
  fX[107] = 184.574996948242188;
  fX[108] = 185.574996948242188;
  fX[109] = 186.574996948242188;
  fX[110] = 187.574996948242188;
  fX[111] = 188.574996948242188;
  fX[112] = 189.574996948242188;
  fX[113] = 190.574996948242188;
  fX[114] = 191.574996948242188;
  fX[115] = 192.574996948242188;
  fX[116] = 193.574996948242188;
  fX[117] = 194.574996948242188;
  fX[118] = 195.574996948242188;
  fX[119] = 196.574996948242188;
  fX[120] = 197.574996948242188;
  fX[121] = 198.574996948242188;
  fX[122] = 199.574996948242188;
  fX[123] = 200.574996948242188;
  fX[124] = 201.574996948242188;
  fX[125] = 202.574996948242188;
  fX[126] = 203.574996948242188;
  fX[127] = 204.574996948242188;
  fX[128] = 205.574996948242188;
  fX[129] = 206.574996948242188;
  fX[130] = 207.574996948242188;
  fX[131] = 208.574996948242188;
  fX[132] = 209.574996948242188;
  fX[133] = 210.574996948242188;
  fX[134] = 211.574996948242188;
  fX[135] = 212.574996948242188;
  fX[136] = 213.574996948242188;
  fX[137] = 214.574996948242188;
  fX[138] = 215.574996948242188;
  fX[139] = 216.574996948242188;
  fX[140] = 217.574996948242188;
  fX[141] = 218.574996948242188;
  fX[142] = 219.574996948242188;
  fX[143] = 220.574996948242188;
  fX[144] = 221.574996948242188;
  fX[145] = 222.574996948242188;
  fX[146] = 223.574996948242188;
  fX[147] = 224.574996948242188;
  fX[148] = 225.574996948242188;
  fX[149] = 226.574996948242188;
  fX[150] = 227.574996948242188;
  fX[151] = 228.574996948242188;
  fX[152] = 229.574996948242188;
  fX[153] = 230.574996948242188;
  fX[154] = 231.574996948242188;
  fX[155] = 232.574996948242188;
  fX[156] = 233.574996948242188;
  fX[157] = 234.574996948242188;
  fX[158] = 235.574996948242188;
  fX[159] = 236.574996948242188;
  fX[160] = 237.574996948242188;
  fX[161] = 238.574996948242188;
  fX[162] = 239.574996948242188;
  fX[163] = 240.574996948242188;
  fX[164] = 241.574996948242188;
  fX[165] = 242.574996948242188;
  fX[166] = 243.574996948242188;
  fX[167] = 244.574996948242188;
  fX[168] = 245.574996948242188;
  fX[169] = 246.574996948242188;
  fX[170] = 247.574996948242188;
  fX[171] = 248.574996948242188;
  fX[172] = 249.574996948242188;
  fX[173] = 250.574996948242188;
  fNPads[0] = 71;
  fNPads[1] = 71;
  fNPads[2] = 71;
  fNPads[3] = 73;
  fNPads[4] = 73;
  fNPads[5] = 73;
  fNPads[6] = 75;
  fNPads[7] = 75;
  fNPads[8] = 75;
  fNPads[9] = 77;
  fNPads[10] = 77;
  fNPads[11] = 77;
  fNPads[12] = 79;
  fNPads[13] = 79;
  fNPads[14] = 79;
  fNPads[15] = 81;
  fNPads[16] = 81;
  fNPads[17] = 81;
  fNPads[18] = 83;
  fNPads[19] = 83;
  fNPads[20] = 83;
  fNPads[21] = 85;
  fNPads[22] = 85;
  fNPads[23] = 85;
  fNPads[24] = 87;
  fNPads[25] = 87;
  fNPads[26] = 87;
  fNPads[27] = 89;
  fNPads[28] = 89;
  fNPads[29] = 89;
  fNPads[30] = 89;
  fNPads[31] = 91;
  fNPads[32] = 91;
  fNPads[33] = 91;
  fNPads[34] = 93;
  fNPads[35] = 93;
  fNPads[36] = 93;
  fNPads[37] = 95;
  fNPads[38] = 95;
  fNPads[39] = 95;
  fNPads[40] = 97;
  fNPads[41] = 97;
  fNPads[42] = 97;
  fNPads[43] = 99;
  fNPads[44] = 99;
  fNPads[45] = 99;
  fNPads[46] = 101;
  fNPads[47] = 101;
  fNPads[48] = 101;
  fNPads[49] = 103;
  fNPads[50] = 103;
  fNPads[51] = 103;
  fNPads[52] = 105;
  fNPads[53] = 105;
  fNPads[54] = 105;
  fNPads[55] = 73;
  fNPads[56] = 73;
  fNPads[57] = 73;
  fNPads[58] = 75;
  fNPads[59] = 75;
  fNPads[60] = 75;
  fNPads[61] = 75;
  fNPads[62] = 77;
  fNPads[63] = 77;
  fNPads[64] = 77;
  fNPads[65] = 79;
  fNPads[66] = 79;
  fNPads[67] = 79;
  fNPads[68] = 81;
  fNPads[69] = 81;
  fNPads[70] = 81;
  fNPads[71] = 81;
  fNPads[72] = 83;
  fNPads[73] = 83;
  fNPads[74] = 83;
  fNPads[75] = 85;
  fNPads[76] = 85;
  fNPads[77] = 85;
  fNPads[78] = 85;
  fNPads[79] = 87;
  fNPads[80] = 87;
  fNPads[81] = 87;
  fNPads[82] = 89;
  fNPads[83] = 89;
  fNPads[84] = 89;
  fNPads[85] = 91;
  fNPads[86] = 91;
  fNPads[87] = 91;
  fNPads[88] = 91;
  fNPads[89] = 93;
  fNPads[90] = 93;
  fNPads[91] = 93;
  fNPads[92] = 95;
  fNPads[93] = 95;
  fNPads[94] = 95;
  fNPads[95] = 95;
  fNPads[96] = 97;
  fNPads[97] = 97;
  fNPads[98] = 97;
  fNPads[99] = 99;
  fNPads[100] = 99;
  fNPads[101] = 99;
  fNPads[102] = 101;
  fNPads[103] = 101;
  fNPads[104] = 101;
  fNPads[105] = 101;
  fNPads[106] = 103;
  fNPads[107] = 103;
  fNPads[108] = 103;
  fNPads[109] = 105;
  fNPads[110] = 105;
  fNPads[111] = 105;
  fNPads[112] = 105;
  fNPads[113] = 107;
  fNPads[114] = 107;
  fNPads[115] = 107;
  fNPads[116] = 109;
  fNPads[117] = 109;
  fNPads[118] = 109;
  fNPads[119] = 111;
  fNPads[120] = 111;
  fNPads[121] = 111;
  fNPads[122] = 111;
  fNPads[123] = 113;
  fNPads[124] = 113;
  fNPads[125] = 113;
  fNPads[126] = 115;
  fNPads[127] = 115;
  fNPads[128] = 115;
  fNPads[129] = 115;
  fNPads[130] = 117;
  fNPads[131] = 117;
  fNPads[132] = 117;
  fNPads[133] = 119;
  fNPads[134] = 119;
  fNPads[135] = 119;
  fNPads[136] = 121;
  fNPads[137] = 121;
  fNPads[138] = 121;
  fNPads[139] = 121;
  fNPads[140] = 123;
  fNPads[141] = 123;
  fNPads[142] = 123;
  fNPads[143] = 125;
  fNPads[144] = 125;
  fNPads[145] = 125;
  fNPads[146] = 125;
  fNPads[147] = 127;
  fNPads[148] = 127;
  fNPads[149] = 127;
  fNPads[150] = 129;
  fNPads[151] = 129;
  fNPads[152] = 129;
  fNPads[153] = 129;
  fNPads[154] = 131;
  fNPads[155] = 131;
  fNPads[156] = 131;
  fNPads[157] = 133;
  fNPads[158] = 133;
  fNPads[159] = 133;
  fNPads[160] = 135;
  fNPads[161] = 135;
  fNPads[162] = 135;
  fNPads[163] = 135;
  fNPads[164] = 137;
  fNPads[165] = 137;
  fNPads[166] = 137;
  fNPads[167] = 139;
  fNPads[168] = 139;
  fNPads[169] = 139;
  fNPads[170] = 139;
  fNPads[171] = 141;
  fNPads[172] = 141;
  fNPads[173] = 141;
}


Double_t AliL3Transform::GetEta(Float_t *xyz)
{
  Double_t r3 = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
  Double_t eta = 0.5 * log((r3+xyz[2])/(r3-xyz[2]));
  return eta;
}

Double_t AliL3Transform::GetPhi(Float_t *xyz)
{
  
  Double_t phi = atan2(xyz[1],xyz[0]);
  //if(phi<0) phi=phi+2*TMath::Pi();
  return phi;
}


Bool_t AliL3Transform::Slice2Sector(Int_t slice, Int_t slicerow, Int_t & sector, Int_t &row) const{
  if(slicerow<0&&slicerow>=fNRow) return kFALSE;
  if(slice<0||slice>=fNSlice) return kFALSE;

  if(slicerow<fNRowLow){
    sector = slice;
    row    = slicerow;
  }
  else {
    sector = slice+fNSlice;
    row    = slicerow-fNRowLow;
  }
  return kTRUE;
}

Bool_t AliL3Transform::Sector2Slice(Int_t & slice, Int_t  sector) const{
  if(sector<0||sector>=fNSector) return kFALSE;
  if(sector<fNSectorLow) slice = sector;
  else          slice = sector - fNSectorLow;
  return kTRUE;
}

Bool_t AliL3Transform::Sector2Slice(Int_t & slice, Int_t & slicerow,Int_t  sector, Int_t row) const{
  if(sector<0||sector>=fNSector||row<0) return kFALSE;
  if(sector<fNSectorLow){
    if(row>=fNRowLow) return kFALSE;
    slice = sector;
    slicerow = row;
  }
  else{
    if(row>=fNRowUp) return kFALSE;
    slice = sector - fNSectorLow;
    slicerow = row + fNRowLow;
  }
  return kTRUE;
}

Double_t AliL3Transform::Row2X(Int_t slicerow){
  if(slicerow<0||slicerow>=fNRow) return 0;
  return fX[slicerow];
}

void AliL3Transform::Local2Global(Float_t *xyz,Int_t slice)
{
  //Transformation to global coordinate system
  Float_t x0 = xyz[0];
  Float_t y0 = xyz[1];
  Float_t cs,sn;
  cs = fCos[slice];
  sn = fSin[slice];
  xyz[0]=x0*cs-y0*sn;
  xyz[1]=x0*sn+y0*cs;
  xyz[2]=xyz[2];//global z=local z
}

void AliL3Transform::Local2GlobalAngle(Float_t *angle,Int_t slice){
  angle[0] = fmod(angle[0]+slice*fPi/9,2*fPi);
}

void AliL3Transform::Global2LocalAngle(Float_t *angle,Int_t slice){
  angle[0] = angle[0]-slice*fPi/9;
  if(angle[0]<0) angle[0]+=2*fPi;
}

void AliL3Transform::Raw2Local(Float_t *xyz,Int_t sector,Int_t row,Float_t pad,Float_t time)
{
  //Transformation from rawdata to local coordinate system
  
  Int_t slice,slicerow;
  Sector2Slice(slice, slicerow, sector, row);  

  xyz[0]=Row2X(slicerow); 
  Int_t npads= fNPads[slicerow];
  if(sector<fNSectorLow)
    xyz[1]=(pad-0.5*(npads-1))*fPadPitchWidthLow;
  else
    xyz[1]=(pad-0.5*(npads-1))*fPadPitchWidthUp;
  xyz[2]=fZWidth*time-3.*fZSigma;
  Int_t sign=-1;
  Int_t nis=fNSectorLow;
  Int_t nos=fNSectorUp;
  
  if((sector<nis)/2 || ((sector-nis)<nos/2)) sign=1;
  xyz[2]=sign*(250.-xyz[2]);

}

void AliL3Transform::Local2Global(Float_t *xyz,Int_t sector,Int_t row)
{
  //Transformation to global coordinate system
  Int_t slice,slicerow;
  Sector2Slice(slice, slicerow, sector, row);  
  Float_t r=Row2X(slicerow);
  Float_t cs = fCos[slice];
  Float_t sn = fSin[slice];

  xyz[0]=r*cs-xyz[1]*sn;
  xyz[1]=r*sn+xyz[1]*cs;
  xyz[2]=xyz[2];//global z=local z
}

Double_t AliL3Transform::GetMaxY(Int_t slicerow)
{
 
 if(slicerow < fNRowLow)
     return fPadPitchWidthLow*fNPads[slicerow]/2; 
 
 else
     return fPadPitchWidthUp*fNPads[slicerow]/2;
}

void AliL3Transform::Global2Local(Float_t *xyz,Int_t sector)
{
  Int_t slice;
  Sector2Slice(slice, sector);  
  Float_t cs = fCos[slice];
  Float_t sn = fSin[slice];
  Float_t x1 = xyz[0]*cs + xyz[1]*sn;
  Float_t y1 = -xyz[0]*sn + xyz[1]*cs;
  xyz[0] = x1;
  xyz[1] = y1;
}

void AliL3Transform::Raw2Global(Float_t *xyz,Int_t sector,Int_t row,Float_t pad,Float_t time)
{
  //Transformation from raw to global coordinates
  
  Raw2Local(xyz,sector,row,pad,time);
  Local2Global(xyz,sector,row);
}

void AliL3Transform::Local2Raw(Float_t *xyz,Int_t sector,Int_t row)
{
  //Transformation from local coordinates to raw
  
  Int_t slice,slicerow;
  Sector2Slice(slice, slicerow, sector, row);  
   
  if(sector<fNSectorLow)
    xyz[1]=xyz[1]/fPadPitchWidthLow+0.5*(fNPads[slicerow]-1);
  else
    xyz[1]=xyz[1]/fPadPitchWidthUp+0.5*(fNPads[slicerow]-1);
  Int_t sign=-1;
  Int_t nis=fNSectorLow;
  Int_t nos=fNSectorUp;
 
  if ((sector<nis/2) || ((sector-nis)<nos/2)) sign=1; 
  xyz[2]=250-sign*xyz[2];
  xyz[2]=(xyz[2]+3.*fZSigma)/fZWidth;
}

void AliL3Transform::Global2Raw(Float_t *xyz,Int_t sector,Int_t row)
{
  //Transformation from global coordinates to raw. 

  Global2Local(xyz,sector);
  Local2Raw(xyz,sector,row);

}
