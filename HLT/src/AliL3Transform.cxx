/* $Id$
 Author: Anders Vestbo <mailto:vestbo@fi.uib.no>,          
         Uli Frankenfeld <mailto:franken@fi.uib.no>
 changes done by Constantin Loizides <mailto:loizides@ikf.physik.uni-frankfurt.de>
 -- Copyright&copy ASV
*/

#include "AliL3Logging.h"
#include "AliL3Transform.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/** \class AliL3Transform 
//<pre>
//_____________________________________________________________
// AliL3Transform
//
// Transformation class for ALICE TPC.
//
// Class which contains all detector specific parameters for the TPC,
// and different useful functions for coordinate transforms.
//
// The class is completely static, which means that no object needs
// to be instantiated. Function calls should then be done like:
//
// AliL3Transform::GetEta(xyz);
//
// IMPORTANT: If used as is, default detector parameters will be used,
//            and you really have to make sure that these correspond to
//            the AliROOT version you are currently working on!!
//            You should therefore always initialize the parameters by
//
//            AliL3Transform::Init(path);
// 
//            where path is a char*, giving the path to where file containing
//            the detector parameter is located. This file should be called
//            "l3transform.config", and can be created with macro exa/Make_Init.C.
//</pre>
*/

ClassImp(AliL3Transform)

Double_t AliL3Transform::fBField = 0.2;
Int_t AliL3Transform::fBFieldFactor = 1;
Int_t AliL3Transform::fVersion = 0;
Int_t AliL3Transform::fNTimeBins = 446 ;
Int_t AliL3Transform::fNRowLow = 64 ;
Int_t AliL3Transform::fNRowUp = 112 ;
Int_t AliL3Transform::fNSectorLow = 36 ;
Int_t AliL3Transform::fNSectorUp = 36 ;
Int_t AliL3Transform::fNSector = 72 ;
Double_t AliL3Transform::fPadPitchWidthLow = 0.400000 ;
Double_t AliL3Transform::fPadPitchWidthUp = 0.600000 ;
Double_t AliL3Transform::fZWidth = 0.56599998474121093750 ;
Double_t AliL3Transform::fZSigma = 0.22880849748219134199 ;
Double_t AliL3Transform::fZOffset = 0.68642549244657402596;
Double_t AliL3Transform::fDiffT = 0.0219999998807907104;
Double_t AliL3Transform::fDiffL = 0.0219999998807907104;
Double_t AliL3Transform::fAnodeWireSpacing = 0.25;
Double_t AliL3Transform::fInnerPadLength = 0.75;
Double_t AliL3Transform::fOuterPadLength = 1.;
Double_t AliL3Transform::fInnerPRFSigma = 0.203811079263687134;
Double_t AliL3Transform::fOuterPRFSigma = 0.299324512481689453;
Double_t AliL3Transform::fTimeSigma = 0.228808626532554626;
Double_t AliL3Transform::fZLength = 250.;
Int_t AliL3Transform::fNSlice = 36 ;
Int_t AliL3Transform::fNRow = 176 ;
Double_t AliL3Transform::fNRotShift = 0.5 ;
Double_t AliL3Transform::fPi = 3.141592653589793 ;
Int_t AliL3Transform::fNPatches = 6;
Int_t AliL3Transform::fRows[6][2] = {{0,31},{32,63},{64,91},{92,119},{120,143},{144,175}};
Int_t AliL3Transform::fNRows[6] = {32,32,28,28,24,32};
Double_t AliL3Transform::fX[176] = {84.570007324218750,
                                    85.320007324218750,
                                    86.070007324218750,
                                    86.820007324218750,
                                    87.570007324218750,
                                    88.320007324218750,
                                    89.070007324218750,
                                    89.820007324218750,
                                    90.570007324218750,
                                    91.320007324218750,
                                    92.070007324218750,
                                    92.820007324218750,
                                    93.570007324218750,
                                    94.320007324218750,
                                    95.070007324218750,
                                    95.820007324218750,
                                    96.570007324218750,
                                    97.320007324218750,
                                    98.070007324218750,
                                    98.820007324218750,
                                    99.570007324218750,
                                    100.320007324218750,
                                    101.070007324218750,
                                    101.820007324218750,
                                    102.570007324218750,
                                    103.320007324218750,
                                    104.070007324218750,
                                    104.820007324218750,
                                    105.570007324218750,
                                    106.320007324218750,
                                    107.070007324218750,
                                    107.820007324218750,
                                    108.570007324218750,
                                    109.320007324218750,
                                    110.070007324218750,
                                    110.820007324218750,
                                    111.570007324218750,
                                    112.320007324218750,
                                    113.070007324218750,
                                    113.820007324218750,
                                    114.570007324218750,
                                    115.320007324218750,
                                    116.070007324218750,
                                    116.820007324218750,
                                    117.570007324218750,
                                    118.320007324218750,
                                    119.070007324218750,
                                    119.820007324218750,
                                    120.570007324218750,
                                    121.320007324218750,
                                    122.070007324218750,
                                    122.820007324218750,
                                    123.570007324218750,
                                    124.320007324218750,
                                    125.070007324218750,
                                    125.820007324218750,
                                    126.570007324218750,
                                    127.320007324218750,
                                    128.070007324218750,
                                    128.820007324218750,
                                    129.570007324218750,
                                    130.320007324218750,
                                    131.070007324218750,
                                    131.820007324218750,
                                    135.054992675781250,
                                    136.054992675781250,
                                    137.054992675781250,
                                    138.054992675781250,
                                    139.054992675781250,
                                    140.054992675781250,
                                    141.054992675781250,
                                    142.054992675781250,
                                    143.054992675781250,
                                    144.054992675781250,
                                    145.054992675781250,
                                    146.054992675781250,
                                    147.054992675781250,
                                    148.054992675781250,
                                    149.054992675781250,
                                    150.054992675781250,
                                    151.054992675781250,
                                    152.054992675781250,
                                    153.054992675781250,
                                    154.054992675781250,
                                    155.054992675781250,
                                    156.054992675781250,
                                    157.054992675781250,
                                    158.054992675781250,
                                    159.054992675781250,
                                    160.054992675781250,
                                    161.054992675781250,
                                    162.054992675781250,
                                    163.054992675781250,
                                    164.054992675781250,
                                    165.054992675781250,
                                    166.054992675781250,
                                    167.054992675781250,
                                    168.054992675781250,
                                    169.054992675781250,
                                    170.054992675781250,
                                    171.054992675781250,
                                    172.054992675781250,
                                    173.054992675781250,
                                    174.054992675781250,
                                    175.054992675781250,
                                    176.054992675781250,
                                    177.054992675781250,
                                    178.054992675781250,
                                    179.054992675781250,
                                    180.054992675781250,
                                    181.054992675781250,
                                    182.054992675781250,
                                    183.054992675781250,
                                    184.054992675781250,
                                    185.054992675781250,
                                    186.054992675781250,
                                    187.054992675781250,
                                    188.054992675781250,
                                    189.054992675781250,
                                    190.054992675781250,
                                    191.054992675781250,
                                    192.054992675781250,
                                    193.054992675781250,
                                    194.054992675781250,
                                    195.054992675781250,
                                    196.054992675781250,
                                    197.054992675781250,
                                    198.054992675781250,
                                    199.054992675781250,
                                    200.054992675781250,
                                    201.054992675781250,
                                    202.054992675781250,
                                    203.054992675781250,
                                    204.054992675781250,
                                    205.054992675781250,
                                    206.054992675781250,
                                    207.054992675781250,
                                    208.054992675781250,
                                    209.054992675781250,
                                    210.054992675781250,
                                    211.054992675781250,
                                    212.054992675781250,
                                    213.054992675781250,
                                    214.054992675781250,
                                    215.054992675781250,
                                    216.054992675781250,
                                    217.054992675781250,
                                    218.054992675781250,
                                    219.054992675781250,
                                    220.054992675781250,
                                    221.054992675781250,
                                    222.054992675781250,
                                    223.054992675781250,
                                    224.054992675781250,
                                    225.054992675781250,
                                    226.054992675781250,
                                    227.054992675781250,
                                    228.054992675781250,
                                    229.054992675781250,
                                    230.054992675781250,
                                    231.054992675781250,
                                    232.054992675781250,
                                    233.054992675781250,
                                    234.054992675781250,
                                    235.054992675781250,
                                    236.054992675781250,
                                    237.054992675781250,
                                    238.054992675781250,
                                    239.054992675781250,
                                    240.054992675781250,
                                    241.054992675781250,
                                    242.054992675781250,
                                    243.054992675781250,
                                    244.054992675781250,
                                    245.054992675781250,
                                    246.054992675781250,
};

Int_t AliL3Transform::fNPads[176] = {67,
                                     67,
                                     67,
                                     69,
                                     69,
                                     69,
                                     71,
                                     71,
                                     71,
                                     73,
                                     73,
                                     73,
                                     75,
                                     75,
                                     75,
                                     77,
                                     77,
                                     77,
                                     79,
                                     79,
                                     79,
                                     81,
                                     81,
                                     81,
                                     83,
                                     83,
                                     83,
                                     85,
                                     85,
                                     85,
                                     87,
                                     87,
                                     87,
                                     89,
                                     89,
                                     89,
                                     91,
                                     91,
                                     91,
                                     93,
                                     93,
                                     93,
                                     95,
                                     95,
                                     95,
                                     97,
                                     97,
                                     97,
                                     99,
                                     99,
                                     99,
                                     101,
                                     101,
                                     101,
                                     103,
                                     103,
                                     103,
                                     105,
                                     105,
                                     105,
                                     107,
                                     107,
                                     107,
                                     109,
                                     73,
                                     75,
                                     75,
                                     75,
                                     77,
                                     77,
                                     77,
                                     77,
                                     79,
                                     79,
                                     79,
                                     81,
                                     81,
                                     81,
                                     83,
                                     83,
                                     83,
                                     83,
                                     85,
                                     85,
                                     85,
                                     87,
                                     87,
                                     87,
                                     87,
                                     89,
                                     89,
                                     89,
                                     91,
                                     91,
                                     91,
                                     93,
                                     93,
                                     93,
                                     93,
                                     95,
                                     95,
                                     95,
                                     97,
                                     97,
                                     97,
                                     97,
                                     99,
                                     99,
                                     99,
                                     101,
                                     101,
                                     101,
                                     103,
                                     103,
                                     103,
                                     103,
                                     105,
                                     105,
                                     105,
                                     107,
                                     107,
                                     107,
                                     107,
                                     109,
                                     109,
                                     109,
                                     111,
                                     111,
                                     111,
                                     113,
                                     113,
                                     113,
                                     113,
                                     115,
                                     115,
                                     115,
                                     117,
                                     117,
                                     117,
                                     117,
                                     119,
                                     119,
                                     119,
                                     121,
                                     121,
                                     121,
                                     123,
                                     123,
                                     123,
                                     123,
                                     125,
                                     125,
                                     125,
                                     127,
                                     127,
                                     127,
                                     127,
                                     129,
                                     129,
                                     129,
                                     131,
                                     131,
                                     131,
                                     133,
                                     133,
                                     133,
                                     133,
                                     135,
                                     135,
                                     135,
                                     137,
                                     137,
                                     137,
                                     137,
                                     139,
                                     139,
};



void AliL3Transform::Init(const Char_t* path)
{
  //Overwrite the parameters with values stored in file "l3transform.config" in path.
  //If file does not exist, old default values will be used.
  
  Char_t *pathname=new Char_t[1024];
  strcpy(pathname,path);
  strcat(pathname,"/l3transform.config");
  
  FILE *fptr=fopen(pathname,"r");
  if(!fptr){
    LOG(AliL3Log::kWarning,"AliL3Transform::Init","File Open")
      <<"Pointer to Config File \""<<pathname<<"\" 0x0!"<<ENDLOG;
    return;
  }

  Char_t d1[250], d2[100], d3[100];
  Int_t dummy=0;
  Double_t ddummy=0.0;

  while(!feof(fptr)) {
    fscanf(fptr,"%s",d1);

    if(strcmp(d1,"fBFieldFactor")==0){fscanf(fptr,"%s %d %s",d2,&dummy,d3);fBFieldFactor=(Int_t)dummy;fBField=fBFieldFactor*0.2;}
    else if(strcmp(d1,"fNTimeBins")==0){fscanf(fptr,"%s %d %s",d2,&dummy,d3);fNTimeBins=(Int_t)dummy;}
    else if(strcmp(d1,"fNRowLow")==0){fscanf(fptr,"%s %d %s",d2,&dummy,d3);fNRowLow=(Int_t)dummy;}    
    if(fNRowLow != 64)
      LOG(AliL3Log::kError,"AliL3Transform::Init","Overflow")
	<<"Number of inner PadRows should be 64! Check and fgrep the code for 64 to see the consequences of this major change!"<<ENDLOG;
    else if(strcmp(d1,"fNRowUp")==0){fscanf(fptr,"%s %d %s",d2,&dummy,d3);fNRowUp=(Int_t)dummy;}
    else if(strcmp(d1,"fNSectorLow")==0){fscanf(fptr,"%s %d %s",d2,&dummy,d3);fNSectorLow=(Int_t)dummy;}
    else if(strcmp(d1,"fNSectorUp")==0){fscanf(fptr,"%s %d %s",d2,&dummy,d3);fNSectorUp=(Int_t)dummy;}
    else if(strcmp(d1,"fNSector")==0){fscanf(fptr,"%s %d %s",d2,&dummy,d3);fNSector=(Int_t)dummy;}
    else if(strcmp(d1,"fPadPitchWidthLow")==0){fscanf(fptr,"%s %lf %s",d2,&ddummy,d3);fPadPitchWidthLow=(Double_t)ddummy;}
    else if(strcmp(d1,"fPadPitchWidthUp")==0){fscanf(fptr,"%s %lf %s",d2,&ddummy,d3);fPadPitchWidthUp=(Double_t)ddummy;}
    else if(strcmp(d1,"fZWidth")==0){fscanf(fptr,"%s %lf %s",d2,&ddummy,d3);fZWidth=(Double_t)ddummy;}
    else if(strcmp(d1,"fZSigma")==0){fscanf(fptr,"%s %lf %s",d2,&ddummy,d3);fZSigma=(Double_t)ddummy;}
    else if(strcmp(d1,"fZLength")==0){fscanf(fptr,"%s %lf %s",d2,&ddummy,d3);fZLength=(Double_t)ddummy;}
    else if(strcmp(d1,"fZOffset")==0){fscanf(fptr,"%s %lf %s",d2,&ddummy,d3);fZOffset=(Double_t)ddummy;}
    else if(strcmp(d1,"fNSlice")==0){fscanf(fptr,"%s %d %s",d2,&dummy,d3);fNSlice=(Int_t)dummy;}
    else if(strcmp(d1,"fDiffT")==0){fscanf(fptr,"%s %lf %s",d2,&ddummy,d3);fDiffT=(Double_t)ddummy;}
    else if(strcmp(d1,"fDiffL")==0){fscanf(fptr,"%s %lf %s",d2,&ddummy,d3);fDiffL=(Double_t)ddummy;}
    else if(strcmp(d1,"fInnerPadLength")==0){fscanf(fptr,"%s %lf %s",d2,&ddummy,d3);fInnerPadLength=(Double_t)ddummy;}
    else if(strcmp(d1,"fOuterPadLength")==0){fscanf(fptr,"%s %lf %s",d2,&ddummy,d3);fOuterPadLength=(Double_t)ddummy;}
    else if(strcmp(d1,"fInnerPRFSigma")==0){fscanf(fptr,"%s %lf %s",d2,&ddummy,d3);fInnerPRFSigma=(Double_t)ddummy;}
    else if(strcmp(d1,"fOuterPRFSigma")==0){fscanf(fptr,"%s %lf %s",d2,&ddummy,d3);fOuterPRFSigma=(Double_t)ddummy;}
    else if(strcmp(d1,"fTimeSigma")==0){fscanf(fptr,"%s %lf %s",d2,&ddummy,d3);fTimeSigma=(Double_t)ddummy;}
    else if(strcmp(d1,"fNRow")==0){
      fscanf(fptr,"%s %d %s",d2,&dummy,d3);fNRow=(Int_t)dummy;
      if(fNRow!=176){
	LOG(AliL3Log::kError,"AliL3Transform::Init","Overflow")<<"Number of PadRows should be 176! Check and fgrep the code for 176 to see the consequences of this major change!"<<ENDLOG;
      }
    }
    else if(strcmp(d1,"fNRotShift")==0){fscanf(fptr,"%s %lf %s",d2,&ddummy,d3);fNRotShift=(Double_t)ddummy;}
    else if(strcmp(d1,"fPi")==0){fscanf(fptr,"%s %lf %s",d2,&ddummy,d3);fPi=(Double_t)ddummy;}
    else if(strcmp(d1,"fX[0]")==0){
      fscanf(fptr,"%s %lf %s",d2,&ddummy,d3);fX[0]=(Double_t)ddummy;
      for(Int_t i=1;i<fNRow;i++){fscanf(fptr,"%s %s %lf %s",d1,d2,&ddummy,d3);fX[i]=(Double_t)ddummy;}
    }
    else if(strcmp(d1,"fNPads[0]")==0){
      fscanf(fptr,"%s %d %s",d2,&dummy,d3);fNPads[0]=(Int_t)dummy;
      for(Int_t i=1;i<fNRow;i++){fscanf(fptr,"%s %s %d %s",d1,d2,&dummy,d3);fNPads[i]=(Int_t)dummy;}
    }
  }

  fclose(fptr);
  delete pathname;
  fVersion=1; //new version

}

Double_t AliL3Transform::GetEta(Float_t *xyz)
{
  Double_t r3 = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
  Double_t eta = 0.5 * log((r3+xyz[2])/(r3-xyz[2]));
  return eta;
}

void AliL3Transform::XYZtoRPhiEta(Float_t *rpe, Float_t *xyz)
{
  rpe[0] = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
  rpe[1] = atan2(xyz[1],xyz[0]);
  rpe[2] = 0.5 * log((rpe[0]+xyz[2])/(rpe[0]-xyz[2]));
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

Bool_t AliL3Transform::Slice2Sector(Int_t slice, Int_t slicerow, Int_t & sector, Int_t &row)
{
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

Bool_t AliL3Transform::Sector2Slice(Int_t & slice, Int_t  sector)
{
  if(sector<0||sector>=fNSector) return kFALSE;
  if(sector<fNSectorLow) slice = sector;
  else          slice = sector - fNSectorLow;
  return kTRUE;
}

Bool_t AliL3Transform::Sector2Slice(Int_t & slice, Int_t & slicerow,Int_t  sector, Int_t row)
{
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
  cs = cos( (2*fPi/18) * (slice+fNRotShift) );
  sn = sin( (2*fPi/18) * (slice+fNRotShift) );
  xyz[0]=x0*cs-y0*sn;
  xyz[1]=x0*sn+y0*cs;
  xyz[2]=xyz[2];//global z=local z
}

void AliL3Transform::Local2GlobalAngle(Float_t *angle,Int_t slice){
  angle[0] = fmod(angle[0]+(slice+fNRotShift)*(2*fPi/18),2*fPi);
}

void AliL3Transform::Global2LocalAngle(Float_t *angle,Int_t slice){
  angle[0] = angle[0]-(slice+fNRotShift)*(2*fPi/18);
  if(angle[0]<0) angle[0]+=2*fPi;
}

void AliL3Transform::Raw2Local(Float_t *xyz,Int_t sector,Int_t row,Float_t pad,Float_t time)
{
  //Transformation from rawdata to local coordinate system
  
  Int_t slice,slicerow;
  Sector2Slice(slice, slicerow, sector, row);  

  //X-Value
  xyz[0]=Row2X(slicerow); 

  //Y-Value
  Int_t npads= fNPads[slicerow];
  if(sector<fNSectorLow)
    xyz[1]=(pad-0.5*(npads-1))*fPadPitchWidthLow;
  else
    xyz[1]=(pad-0.5*(npads-1))*fPadPitchWidthUp;

  //Z-Value (remember PULSA Delay)
  //xyz[2]=fZWidth*time-3.*fZSigma;
  xyz[2]=fZWidth*time-fZOffset;
  if(slice < 18)
    xyz[2]=fZLength-xyz[2];
  else
    xyz[2]=xyz[2]-fZLength;
  
}

void AliL3Transform::Local2Global(Float_t *xyz,Int_t sector,Int_t row)
{
  //Transformation to global coordinate system
  Int_t slice,slicerow;
  Sector2Slice(slice, slicerow, sector, row);  
  Float_t r=Row2X(slicerow);
  Float_t cs = cos( (2*fPi/18) * (slice+fNRotShift) );
  Float_t sn = sin( (2*fPi/18) * (slice+fNRotShift) );

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
  Float_t cs = cos( (2*fPi/18) * (slice+fNRotShift) );
  Float_t sn = sin( (2*fPi/18) * (slice+fNRotShift) );
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
  xyz[2]=fZLength-sign*xyz[2];
  xyz[2]=(xyz[2]+fZOffset)/fZWidth;

}

void AliL3Transform::Global2Raw(Float_t *xyz,Int_t sector,Int_t row)
{
  //Transformation from global coordinates to raw. 

  Global2Local(xyz,sector);
  Local2Raw(xyz,sector,row);

}
