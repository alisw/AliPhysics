// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>, Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ASV

#include "AliL3Logging.h"
#include "AliL3Transform.h"
//#include <TFile.h>
#include <math.h>

//_____________________________________________________________
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
  fNTimeBins = 446;
  fNRowLow = 64;
  fNRowUp = 112;
  fNSectorLow = 36;
  fNSectorUp = 36;
  fNSector = 72;
  fPadPitchWidthLow = 0.400000;
  fPadPitchWidthUp = 0.600000;
  fZWidth = 0.56599998474121093750;
  fZSigma = 0.22880849748219134199;

  //slices:
  fNSlice = 36;
  fNRow = 176;
  fPi = 3.141592653589793;
/*
  for(Int_t i=0;i<36;i++){
    fCos[i] = cos(fPi/9*i);
    fSin[i] = sin(fPi/9*i);
  }
*/
  for(Int_t i=0;i<36;i++){
    fCos[i] = cos( (2*fPi/18) * (i+0.5) );
    fSin[i] = sin( (2*fPi/18) * (i+0.5) );
  }

  fX[0] = 84.570007324218750;
  fX[1] = 85.320007324218750;
  fX[2] = 86.070007324218750;
  fX[3] = 86.820007324218750;
  fX[4] = 87.570007324218750;
  fX[5] = 88.320007324218750;
  fX[6] = 89.070007324218750;
  fX[7] = 89.820007324218750;
  fX[8] = 90.570007324218750;
  fX[9] = 91.320007324218750;
  fX[10] = 92.070007324218750;
  fX[11] = 92.820007324218750;
  fX[12] = 93.570007324218750;
  fX[13] = 94.320007324218750;
  fX[14] = 95.070007324218750;
  fX[15] = 95.820007324218750;
  fX[16] = 96.570007324218750;
  fX[17] = 97.320007324218750;
  fX[18] = 98.070007324218750;
  fX[19] = 98.820007324218750;
  fX[20] = 99.570007324218750;
  fX[21] = 100.320007324218750;
  fX[22] = 101.070007324218750;
  fX[23] = 101.820007324218750;
  fX[24] = 102.570007324218750;
  fX[25] = 103.320007324218750;
  fX[26] = 104.070007324218750;
  fX[27] = 104.820007324218750;
  fX[28] = 105.570007324218750;
  fX[29] = 106.320007324218750;
  fX[30] = 107.070007324218750;
  fX[31] = 107.820007324218750;
  fX[32] = 108.570007324218750;
  fX[33] = 109.320007324218750;
  fX[34] = 110.070007324218750;
  fX[35] = 110.820007324218750;
  fX[36] = 111.570007324218750;
  fX[37] = 112.320007324218750;
  fX[38] = 113.070007324218750;
  fX[39] = 113.820007324218750;
  fX[40] = 114.570007324218750;
  fX[41] = 115.320007324218750;
  fX[42] = 116.070007324218750;
  fX[43] = 116.820007324218750;
  fX[44] = 117.570007324218750;
  fX[45] = 118.320007324218750;
  fX[46] = 119.070007324218750;
  fX[47] = 119.820007324218750;
  fX[48] = 120.570007324218750;
  fX[49] = 121.320007324218750;
  fX[50] = 122.070007324218750;
  fX[51] = 122.820007324218750;
  fX[52] = 123.570007324218750;
  fX[53] = 124.320007324218750;
  fX[54] = 125.070007324218750;
  fX[55] = 125.820007324218750;
  fX[56] = 126.570007324218750;
  fX[57] = 127.320007324218750;
  fX[58] = 128.070007324218750;
  fX[59] = 128.820007324218750;
  fX[60] = 129.570007324218750;
  fX[61] = 130.320007324218750;
  fX[62] = 131.070007324218750;
  fX[63] = 131.820007324218750;
  fX[64] = 135.054992675781250;
  fX[65] = 136.054992675781250;
  fX[66] = 137.054992675781250;
  fX[67] = 138.054992675781250;
  fX[68] = 139.054992675781250;
  fX[69] = 140.054992675781250;
  fX[70] = 141.054992675781250;
  fX[71] = 142.054992675781250;
  fX[72] = 143.054992675781250;
  fX[73] = 144.054992675781250;
  fX[74] = 145.054992675781250;
  fX[75] = 146.054992675781250;
  fX[76] = 147.054992675781250;
  fX[77] = 148.054992675781250;
  fX[78] = 149.054992675781250;
  fX[79] = 150.054992675781250;
  fX[80] = 151.054992675781250;
  fX[81] = 152.054992675781250;
  fX[82] = 153.054992675781250;
  fX[83] = 154.054992675781250;
  fX[84] = 155.054992675781250;
  fX[85] = 156.054992675781250;
  fX[86] = 157.054992675781250;
  fX[87] = 158.054992675781250;
  fX[88] = 159.054992675781250;
  fX[89] = 160.054992675781250;
  fX[90] = 161.054992675781250;
  fX[91] = 162.054992675781250;
  fX[92] = 163.054992675781250;
  fX[93] = 164.054992675781250;
  fX[94] = 165.054992675781250;
  fX[95] = 166.054992675781250;
  fX[96] = 167.054992675781250;
  fX[97] = 168.054992675781250;
  fX[98] = 169.054992675781250;
  fX[99] = 170.054992675781250;
  fX[100] = 171.054992675781250;
  fX[101] = 172.054992675781250;
  fX[102] = 173.054992675781250;
  fX[103] = 174.054992675781250;
  fX[104] = 175.054992675781250;
  fX[105] = 176.054992675781250;
  fX[106] = 177.054992675781250;
  fX[107] = 178.054992675781250;
  fX[108] = 179.054992675781250;
  fX[109] = 180.054992675781250;
  fX[110] = 181.054992675781250;
  fX[111] = 182.054992675781250;
  fX[112] = 183.054992675781250;
  fX[113] = 184.054992675781250;
  fX[114] = 185.054992675781250;
  fX[115] = 186.054992675781250;
  fX[116] = 187.054992675781250;
  fX[117] = 188.054992675781250;
  fX[118] = 189.054992675781250;
  fX[119] = 190.054992675781250;
  fX[120] = 191.054992675781250;
  fX[121] = 192.054992675781250;
  fX[122] = 193.054992675781250;
  fX[123] = 194.054992675781250;
  fX[124] = 195.054992675781250;
  fX[125] = 196.054992675781250;
  fX[126] = 197.054992675781250;
  fX[127] = 198.054992675781250;
  fX[128] = 199.054992675781250;
  fX[129] = 200.054992675781250;
  fX[130] = 201.054992675781250;
  fX[131] = 202.054992675781250;
  fX[132] = 203.054992675781250;
  fX[133] = 204.054992675781250;
  fX[134] = 205.054992675781250;
  fX[135] = 206.054992675781250;
  fX[136] = 207.054992675781250;
  fX[137] = 208.054992675781250;
  fX[138] = 209.054992675781250;
  fX[139] = 210.054992675781250;
  fX[140] = 211.054992675781250;
  fX[141] = 212.054992675781250;
  fX[142] = 213.054992675781250;
  fX[143] = 214.054992675781250;
  fX[144] = 215.054992675781250;
  fX[145] = 216.054992675781250;
  fX[146] = 217.054992675781250;
  fX[147] = 218.054992675781250;
  fX[148] = 219.054992675781250;
  fX[149] = 220.054992675781250;
  fX[150] = 221.054992675781250;
  fX[151] = 222.054992675781250;
  fX[152] = 223.054992675781250;
  fX[153] = 224.054992675781250;
  fX[154] = 225.054992675781250;
  fX[155] = 226.054992675781250;
  fX[156] = 227.054992675781250;
  fX[157] = 228.054992675781250;
  fX[158] = 229.054992675781250;
  fX[159] = 230.054992675781250;
  fX[160] = 231.054992675781250;
  fX[161] = 232.054992675781250;
  fX[162] = 233.054992675781250;
  fX[163] = 234.054992675781250;
  fX[164] = 235.054992675781250;
  fX[165] = 236.054992675781250;
  fX[166] = 237.054992675781250;
  fX[167] = 238.054992675781250;
  fX[168] = 239.054992675781250;
  fX[169] = 240.054992675781250;
  fX[170] = 241.054992675781250;
  fX[171] = 242.054992675781250;
  fX[172] = 243.054992675781250;
  fX[173] = 244.054992675781250;
  fX[174] = 245.054992675781250;
  fX[175] = 246.054992675781250;
  fNPads[0] = 67;
  fNPads[1] = 67;
  fNPads[2] = 67;
  fNPads[3] = 69;
  fNPads[4] = 69;
  fNPads[5] = 69;
  fNPads[6] = 71;
  fNPads[7] = 71;
  fNPads[8] = 71;
  fNPads[9] = 73;
  fNPads[10] = 73;
  fNPads[11] = 73;
  fNPads[12] = 75;
  fNPads[13] = 75;
  fNPads[14] = 75;
  fNPads[15] = 77;
  fNPads[16] = 77;
  fNPads[17] = 77;
  fNPads[18] = 79;
  fNPads[19] = 79;
  fNPads[20] = 79;
  fNPads[21] = 81;
  fNPads[22] = 81;
  fNPads[23] = 81;
  fNPads[24] = 83;
  fNPads[25] = 83;
  fNPads[26] = 83;
  fNPads[27] = 85;
  fNPads[28] = 85;
  fNPads[29] = 85;
  fNPads[30] = 87;
  fNPads[31] = 87;
  fNPads[32] = 87;
  fNPads[33] = 89;
  fNPads[34] = 89;
  fNPads[35] = 89;
  fNPads[36] = 91;
  fNPads[37] = 91;
  fNPads[38] = 91;
  fNPads[39] = 93;
  fNPads[40] = 93;
  fNPads[41] = 93;
  fNPads[42] = 95;
  fNPads[43] = 95;
  fNPads[44] = 95;
  fNPads[45] = 97;
  fNPads[46] = 97;
  fNPads[47] = 97;
  fNPads[48] = 99;
  fNPads[49] = 99;
  fNPads[50] = 99;
  fNPads[51] = 101;
  fNPads[52] = 101;
  fNPads[53] = 101;
  fNPads[54] = 103;
  fNPads[55] = 103;
  fNPads[56] = 103;
  fNPads[57] = 105;
  fNPads[58] = 105;
  fNPads[59] = 105;
  fNPads[60] = 107;
  fNPads[61] = 107;
  fNPads[62] = 107;
  fNPads[63] = 109;
  fNPads[64] = 73;
  fNPads[65] = 75;
  fNPads[66] = 75;
  fNPads[67] = 75;
  fNPads[68] = 77;
  fNPads[69] = 77;
  fNPads[70] = 77;
  fNPads[71] = 77;
  fNPads[72] = 79;
  fNPads[73] = 79;
  fNPads[74] = 79;
  fNPads[75] = 81;
  fNPads[76] = 81;
  fNPads[77] = 81;
  fNPads[78] = 83;
  fNPads[79] = 83;
  fNPads[80] = 83;
  fNPads[81] = 83;
  fNPads[82] = 85;
  fNPads[83] = 85;
  fNPads[84] = 85;
  fNPads[85] = 87;
  fNPads[86] = 87;
  fNPads[87] = 87;
  fNPads[88] = 87;
  fNPads[89] = 89;
  fNPads[90] = 89;
  fNPads[91] = 89;
  fNPads[92] = 91;
  fNPads[93] = 91;
  fNPads[94] = 91;
  fNPads[95] = 93;
  fNPads[96] = 93;
  fNPads[97] = 93;
  fNPads[98] = 93;
  fNPads[99] = 95;
  fNPads[100] = 95;
  fNPads[101] = 95;
  fNPads[102] = 97;
  fNPads[103] = 97;
  fNPads[104] = 97;
  fNPads[105] = 97;
  fNPads[106] = 99;
  fNPads[107] = 99;
  fNPads[108] = 99;
  fNPads[109] = 101;
  fNPads[110] = 101;
  fNPads[111] = 101;
  fNPads[112] = 103;
  fNPads[113] = 103;
  fNPads[114] = 103;
  fNPads[115] = 103;
  fNPads[116] = 105;
  fNPads[117] = 105;
  fNPads[118] = 105;
  fNPads[119] = 107;
  fNPads[120] = 107;
  fNPads[121] = 107;
  fNPads[122] = 107;
  fNPads[123] = 109;
  fNPads[124] = 109;
  fNPads[125] = 109;
  fNPads[126] = 111;
  fNPads[127] = 111;
  fNPads[128] = 111;
  fNPads[129] = 113;
  fNPads[130] = 113;
  fNPads[131] = 113;
  fNPads[132] = 113;
  fNPads[133] = 115;
  fNPads[134] = 115;
  fNPads[135] = 115;
  fNPads[136] = 117;
  fNPads[137] = 117;
  fNPads[138] = 117;
  fNPads[139] = 117;
  fNPads[140] = 119;
  fNPads[141] = 119;
  fNPads[142] = 119;
  fNPads[143] = 121;
  fNPads[144] = 121;
  fNPads[145] = 121;
  fNPads[146] = 123;
  fNPads[147] = 123;
  fNPads[148] = 123;
  fNPads[149] = 123;
  fNPads[150] = 125;
  fNPads[151] = 125;
  fNPads[152] = 125;
  fNPads[153] = 127;
  fNPads[154] = 127;
  fNPads[155] = 127;
  fNPads[156] = 127;
  fNPads[157] = 129;
  fNPads[158] = 129;
  fNPads[159] = 129;
  fNPads[160] = 131;
  fNPads[161] = 131;
  fNPads[162] = 131;
  fNPads[163] = 133;
  fNPads[164] = 133;
  fNPads[165] = 133;
  fNPads[166] = 133;
  fNPads[167] = 135;
  fNPads[168] = 135;
  fNPads[169] = 135;
  fNPads[170] = 137;
  fNPads[171] = 137;
  fNPads[172] = 137;
  fNPads[173] = 137;
  fNPads[174] = 139;
  fNPads[175] = 139;
}

Double_t AliL3Transform::GetEta(Float_t *xyz)
{
  Double_t r3 = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
  Double_t eta = 0.5 * log((r3+xyz[2])/(r3-xyz[2]));
  return eta;
}

Double_t AliL3Transform::GetEta(Int_t padrow,Int_t pad,Int_t time)
{
  Float_t xyz[3];
  Int_t sector,row;
  Slice2Sector(0,padrow,sector,row);
  Raw2Local(xyz,sector,row,pad,time);
  
  return GetEta(xyz);
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
  angle[0] = fmod(angle[0]+(slice+0.5)*(2*fPi/18),2*fPi);
}

void AliL3Transform::Global2LocalAngle(Float_t *angle,Int_t slice){
  angle[0] = angle[0]-(slice+0.5)*(2*fPi/18);
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

void AliL3Transform::Global2Local(Float_t *xyz,Int_t sector,Bool_t isSlice)
{
  
  Int_t slice;
  if(!isSlice)
    Sector2Slice(slice, sector);  
  else
    slice = sector;
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
