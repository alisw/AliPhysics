/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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
#include <iostream.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <AliITS.h>
#include <AliITSgeom.h>
#include <AliITSdigit.h>
#include "AliITSDDLRawData.h"

ClassImp(AliITSDDLRawData)

////////////////////////////////////////////////////////////////////////////////////////
AliITSDDLRawData::AliITSDDLRawData(){
  fIndex=-1;
  fHalfStaveModule=-1;
}

////////////////////////////////////////////////////////////////////////////////////////

AliITSDDLRawData::AliITSDDLRawData(const AliITSDDLRawData &source){
  // Copy Constructor
  this->fIndex=source.fIndex;
  this->fHalfStaveModule=source.fHalfStaveModule;
  return;
}

////////////////////////////////////////////////////////////////////////////////////////

AliITSDDLRawData& AliITSDDLRawData::operator=(const AliITSDDLRawData &source){
  //Assigment operator
  this->fIndex=source.fIndex;
  this->fHalfStaveModule=source.fHalfStaveModule;
  return *this;
}

////////////////////////////////////////////////////////////////////////////////////////
//STRIP 
//

void AliITSDDLRawData::GetDigitsSSD(TClonesArray *ITSdigits,Int_t mod, ULong_t *buf){
  Int_t ix;
  Int_t iz;
  Int_t is;
  ULong_t Word;
  ULong_t BaseWord;
  Int_t ndigits = ITSdigits->GetEntries();
  AliITSdigit *digs;
  if(ndigits){
    for (Int_t digit=0;digit<ndigits;digit++) {
      digs = (AliITSdigit*)ITSdigits->UncheckedAt(digit);
      iz=digs->fCoord1;  // If iz==0, N side and if iz=1 P side
      ix=digs->fCoord2;  // Strip Numbar
      is=digs->fSignal;  // ADC Signal
      // cout<<" Module:"<<mod-500<<" N/P side:"<<iz<<" Strip Number:"<<ix<<" Amplidute:"<<is-1<<endl;
      BaseWord=0;
      Word=is-1;
      PackWord(BaseWord,Word,0,9);//ADC data
      Word=ix;
      PackWord(BaseWord,Word,10,19);//Strip Number
      Word=iz;      
      PackWord(BaseWord,Word,20,20);//ADC Channel ID (N or P side)
      Word=mod;
      PackWord(BaseWord,Word,21,31);//ADC module ID
      fIndex++;
      buf[fIndex]=BaseWord;
    }//end for
  }//end if
  return;
}//end GetDigitsSSD

////////////////////////////////////////////////////////////////////////////////////////
//Silicon Drift Detector
//

void AliITSDDLRawData::GetDigitsSDD(TClonesArray *ITSdigits,Int_t mod, ULong_t *buf){  
  Int_t ix;
  Int_t iz;
  Int_t is;
  ULong_t Word;
  ULong_t BaseWord;
  Int_t ndigits = ITSdigits->GetEntries();
  AliITSdigit *digs;
  if(ndigits){
    //cout<<"Mudule "<<mod<<" number of digits "<<ndigits<<endl;
    for (Int_t digit=0;digit<ndigits;digit++) {
      digs = (AliITSdigit*)ITSdigits->UncheckedAt(digit);
      iz=digs->fCoord1;  // Anode
      ix=digs->fCoord2;  // Time
      is=digs->fSignal;  // ADC Signal
      //      cout<<"Amplitude value:"<<is<<" Time Bucket:"<<ix<<" Anode:"<<iz<<endl;
      if (is>255){cout<<"WARNING (!) bits words is needed)!!!\n";}
      BaseWord=0;
      /*
      //10 bits words for amplitude value
      Word=is;
      PackWord(BaseWord,Word,0,9);//ADC data
      Word=ix;
      PackWord(BaseWord,Word,10,17);//Time bucket
      Word=iz;
      PackWord(BaseWord,Word,18,26);//Anode Number
      Word=mod;
      PackWord(BaseWord,Word,27,31);//Module number
      */
      
      //8bits words for amplitude value
      Word=is;
      PackWord(BaseWord,Word,0,7);//ADC data
      Word=ix;
      PackWord(BaseWord,Word,8,15);//Time bucket
      Word=iz;
      PackWord(BaseWord,Word,16,24);//Anode Number
      Word=mod;
      PackWord(BaseWord,Word,25,31);//Module number
     
      fIndex++;
      buf[fIndex]=BaseWord;
    }//end for
  }//end if
  return;
}//end GetDigitsSDD

////////////////////////////////////////////////////////////////////////////////////////
//PIXEL 
//

void AliITSDDLRawData::GetDigitsSPD(TClonesArray *ITSdigits,Int_t mod, ULong_t *buf){
  Int_t ix;
  Int_t iz;
  Int_t ChipNo=0;
  ULong_t BaseWord=0;
  ULong_t HitRow=0;
  Int_t ChipHitCount=0;  //Number of Hit in the current chip
  Int_t PreviousChip=-1; //Previuos chip respect to the actual aone
  Int_t ndigits = ITSdigits->GetEntries(); //number of digits in the current module
  //cout<<"      Number of digits in the current module:"<<ndigits<<" module:"<<mod<<endl;
  AliITSdigit *digs;
  fHalfStaveModule++;    //It's a private variable used to distinguish between the firs  
                         //and the second module of an Half Stave Module
  if(ndigits){
    //loop over digits
    for (Int_t digit=0;digit<ndigits;digit++){
      digs = (AliITSdigit*)ITSdigits->UncheckedAt(digit);
      /*---------------------------------------------------------------------------
       *     Each module contains 5 read out chips of 256 rows and 32 columns.
       *     So, the cell number in Z direction varies from 1 to 160.  Therefore,
       *     to get the chip address (0 to 4), we need to divide column number by 32.
       *     ---------------------------------------------------------------------*/
      iz=digs->fCoord1;  // Cell number in Z direction 
      ix=digs->fCoord2;  // Cell number in X direction
      ChipNo=iz/32;
      HitRow=iz-ChipNo*32;
      if(fHalfStaveModule){
	ChipNo+=5;
	fHalfStaveModule=-1;
      }//end if
      //cout<<"Chip number of the current digit:"<<ChipNo<<" Row:"<<HitRow<<" Column:"<<ix<<endl;
      if(PreviousChip==-1){
	//loop over chip without digits 
	//Even if there aren't digits for a given chip 
	//the chip header and the chip trailer are stored
	for(Int_t i=0;i<(iz/32);i++){
	  if(ChipNo>4)
	    WriteChipHeader(i+5,(mod/2),BaseWord);
	  else
	    WriteChipHeader(i,(mod/2),BaseWord);
	  WriteChipTrailer(buf,ChipHitCount,BaseWord);
	}//end for
	PreviousChip=ChipNo;
	WriteChipHeader(ChipNo,(mod/2),BaseWord);
	ChipHitCount++;
	WriteHit(buf,ix,HitRow,BaseWord);
      }//end if
      else{
	ChipHitCount++;
	if(PreviousChip!=ChipNo){
	  WriteChipTrailer(buf,ChipHitCount-1,BaseWord);
	  for(Int_t i=PreviousChip+1;i<ChipNo;i++){
	    WriteChipHeader(i,(mod/2),BaseWord);
	    WriteChipTrailer(buf,0,BaseWord);
	  }//end for
	  WriteChipHeader(ChipNo,(mod/2),BaseWord);
	  PreviousChip=ChipNo;
	}//end if
	WriteHit(buf,ix,HitRow,BaseWord);
      }//end else
    }//end for
    //Even if there aren't digits for a given chip 
    //the chip header and the chip trailer are stored
    Int_t End=4;
    if(ChipNo>4)End+=5;
    WriteChipTrailer(buf,ChipHitCount,BaseWord);
    for(Int_t i=ChipNo+1;i<=End;i++){
      WriteChipHeader(i,(mod/2),BaseWord);
      WriteChipTrailer(buf,0,BaseWord);
    }//end for
  }//end if
  else{
    //In this module there aren't digits but
    //the chip header and chip trailer are store anyway
    if(fHalfStaveModule){
      ChipNo=5;
      fHalfStaveModule=-1;
    }//end if
    for(Int_t i=0;i<5;i++){
      WriteChipHeader(ChipNo+i,(mod/2),BaseWord);
      WriteChipTrailer(buf,ChipHitCount,BaseWord);
    }//end for
  }//end else 
  return;
}//end GetDigitsSPD

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AliITSDDLRawData::PackWord(ULong_t &BaseWord, ULong_t Word, Int_t StartBit, Int_t StopBit){
  ULong_t DummyWord,OffSet;
  Int_t   Length;
  ULong_t Sum;
  //The BaseWord is being filled with 1 from StartBit to StopBit
  Length=StopBit-StartBit+1;
  Sum=(ULong_t)pow(2,Length)-1;
  if(Word > Sum){
    cout<<"WARNING::Word to be filled is not within desired length"<<endl;
    cout<<"Word:"<<Word<<" Start bit:"<<StartBit<<" Stop Bit:"<<StopBit<<endl;
    exit(-1);
  }
  OffSet=Sum;
  OffSet<<=StartBit;
  BaseWord=BaseWord|OffSet;
  //The Word to be filled is shifted to the position StartBit
  //and the remaining  Left and Right bits are filled with 1
  Sum=(ULong_t)pow(2,StartBit)-1;
  DummyWord=0xFFFFFFFF<<Length;
  DummyWord +=Word;
  DummyWord<<=StartBit;
  DummyWord+=Sum;
  BaseWord=BaseWord&DummyWord;
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AliITSDDLRawData::UnpackWord(ULong_t PackedWord, Int_t StartBit, Int_t StopBit, ULong_t &Word){ 	
  ULong_t OffSet;
  Int_t Length;
  Length=StopBit-StartBit+1;
  OffSet=(ULong_t)pow(2,Length)-1;
  OffSet<<=StartBit;
  Word=PackedWord&OffSet;
  Word>>=StartBit;
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliITSDDLRawData::RawDataSPD(AliITS *ITS,TTree *TD ,Int_t LDCsNumber){

  //Silicon Pixel Detector
  const Int_t DDLNumber=20;       // Number of DDL in SPD
  const Int_t ModulePerDDL=12;    // Number of modules in each DDL 
  //Row    ==> DDLs
  //Column ==> Modules
  Int_t SPDMap[DDLNumber][ModulePerDDL]={{ 0, 1, 4, 5, 80, 81, 84, 85, 88, 89, 92, 93},
					 { 2, 3, 6, 7, 82, 83, 86, 87, 90, 91, 94, 95},
					 { 8, 9,12,13, 96, 97,100,101,104,105,108,109},
					 {10,11,14,15, 98, 99,102,103,106,107,110,111},
					 {16,17,20,21,112,113,116,117,120,121,124,125},
					 {18,19,22,23,114,115,118,119,122,123,126,127},
					 {24,25,28,29,128,129,132,133,136,137,140,141},
					 {26,27,30,31,130,131,134,135,138,139,142,143},
					 {32,33,36,37,144,145,148,149,152,153,156,157},
					 {34,35,38,39,146,147,150,151,154,155,158,159},
					 {40,41,44,45,160,161,164,165,168,169,172,173},
					 {42,43,46,47,162,163,166,167,170,171,174,175},
					 {48,47,50,51,176,177,180,181,184,185,188,189},
					 {50,51,54,55,178,179,182,183,186,187,190,191},
					 {56,57,60,61,192,193,196,197,200,201,204,205},
					 {58,59,62,63,194,195,198,199,202,203,206,207},
					 {64,65,68,69,208,209,212,213,216,217,220,221},
					 {66,67,70,71,210,211,214,215,218,219,222,223},
					 {72,73,76,77,224,225,228,229,232,233,236,237},
					 {74,75,78,79,226,227,230,231,234,235,238,239}};
  Int_t DDLPerFile=DDLNumber/LDCsNumber;
  if(DDLNumber%LDCsNumber)DDLPerFile++;
  cout<<"Number of DDL per File: "<<DDLPerFile<<endl;
  Int_t subd=0;
  Int_t CountDDL=0;
  Int_t SliceNumber=1;
  const Int_t size=21000; //256*32*5=40960 max number of digits per module
  ULong_t buf[size];      //One buffer cell can contain 2 digits 
  fIndex=-1;
  TClonesArray *ITSdigits  = ITS->DigitsAddress(subd);

  Int_t nbytes = 0; 
  char fileName[15];

  ULong_t MiniHeaderPosition=0; //variable used to store the position of the Mini Header inside a file
  ofstream outfile;         // logical name of the output file 
  Int_t Flag=0;             // 0==> Uncompressed data 1==>Compressed Data
  Int_t Detector=1;         // 1==>ITS (Pixel) 2==>ITS (Drift) 3==>ITS (Strip) 0==>TPC ......
  ULong_t Size=0;           // Size of the data block that follows the mini header
  Int_t MagicWord=0x123456;  // Magic word used to distinguish between data and garbage
  sprintf(fileName,"SPDslice%d",SliceNumber); //The name of  the output file. There are as many slides as the number of LDC
  outfile.open(fileName,ios::binary);
  ULong_t MiniHeader[3];
  Int_t MiniHeaderSize=sizeof(ULong_t)*3;
  Int_t Version=1;          //Version of the mini header 
  //loop over DDLs
  for(Int_t i=0;i<DDLNumber;i++){
    CountDDL++;
    //write Dummy MINI HEADER
    MiniHeader[0]=Size;
    PackWord(MiniHeader[1],MagicWord,8,31);
    PackWord(MiniHeader[1],Detector,0,7);
    PackWord(MiniHeader[2],i,16,31);
    PackWord(MiniHeader[2],Flag,8,15);
    PackWord(MiniHeader[2],Version,0,7);
    MiniHeaderPosition=outfile.tellp();
    outfile.write((char*)(MiniHeader),MiniHeaderSize);
    //Loops over Modules of a particular DDL
    for (Int_t mod=0; mod<ModulePerDDL; mod++){
      ITS->ResetDigits();
      nbytes += TD->GetEvent(SPDMap[i][mod]);
      //For each Module, buf contains the array of data words in Binary format	  
      //fIndex gives the number of 32 bits words in the buffer for each module
      GetDigitsSPD(ITSdigits,SPDMap[i][mod],buf);
      outfile.write((char *)buf,((fIndex+1)*sizeof(ULong_t)));
      for(Int_t i=0;i<(fIndex+1);i++){
	buf[i]=0;
      }//end for
      fIndex=-1;
    }//end for
    
    //Write REAL MINI HEADER
    ULong_t CurrentFilePosition=outfile.tellp();
    outfile.seekp(MiniHeaderPosition);
    Size=CurrentFilePosition-MiniHeaderPosition-MiniHeaderSize;
    
    outfile.write((char*)(&Size),sizeof(ULong_t));
    outfile.seekp(CurrentFilePosition);
    if(CountDDL==DDLPerFile){
      outfile.close();
      SliceNumber++;
      sprintf(fileName,"SPDslice%d",SliceNumber); 
      if(i!=(DDLNumber-1))
	outfile.open(fileName,ios::binary);
      CountDDL=0;
    }//end if
  }//end for
  outfile.close();
  return 0;  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliITSDDLRawData::RawDataSSD(AliITS *ITS,TTree *TD ,Int_t LDCsNumber){
  //Strip detector
  const Int_t DDLNumber=16;        // Number of DDL in SSD
  const Int_t ModulePerDDL=109;    // Number of modules in each DDL 
  //DDL from 32 to 47 (16 DDL)
  //Row    ==> DDLs
  //Column ==> Modules
  Int_t SSDMap[DDLNumber][ModulePerDDL]={
    //104
    //DDL[32][]=
    { 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510,
      522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532,
      1204,1205,1206,1207,1208,1209,1210,1211,1212,1213,1214,
      1226,1227,1228,1229,1230,1231,1232,1233,1234,1235,1236,
      1248,1249,1250,1251,1252,1253,1254,1255,1256,1257,1258,1259,
      2098,2099,2100,2101,2102,2103,2104,2105,2106,2107,2108,2109,
      2123,2124,2125,2126,2127,2128,2129,2130,2131,2132,2133,2134,
      2148,2149,2150,2151,2152,2153,2154,2155,2156,2157,2158,2159,
      2173,2174,2175,2176,2177,2178,2179,2180,2181,2182,2183,2184,-1,-1,-1,-1,-1},    
    //93
    //DDL[33][]=
    { 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576,
      588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598,
      610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620,
      1273,1274,1275,1276,1277,1278,1279,1280,1281,1282,1283,1284,
      1298,1299,1300,1301,1302,1303,1304,1305,1306,1307,1308,1309,
      1323,1324,1325,1326,1327,1328,1329,1330,1331,1332,1333,1334,
      1348,1349,1350,1351,1352,1353,1354,1355,1356,1357,1358,1359,
      1373,1374,1375,1376,1377,1378,1379,1380,1381,1382,1383,1384,
      -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,-1,-1,-1,-1},
    //103
    //DDL[34][]=
    { 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642,
      654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664,
      676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686,
      698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708,
      720, 721, 722, 723, 724, 725, 726, 727, 728, 729, 730,
      1398,1399,1400,1401,1402,1403,1404,1405,1406,1407,1408,1409,
      1423,1424,1425,1426,1427,1428,1429,1430,1431,1432,1433,1434,
      1448,1449,1450,1451,1452,1453,1454,1455,1456,1457,1458,1459,
      1473,1474,1475,1476,1477,1478,1479,1480,1481,1482,1483,1484,-1,-1,-1,-1,-1,-1},
    //104
    //DDL[35][]=
    { 742, 743, 744, 745, 746, 747, 748, 749, 750, 751, 752,
      764, 765, 766, 767, 768, 769, 770, 771, 772, 773, 774,
      786, 787, 788, 789, 790, 791, 792, 793, 794, 795, 796,
      808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818,
      1498,1499,1500,1501,1502,1503,1504,1505,1506,1507,1508,1509,
      1523,1524,1525,1526,1527,1528,1529,1530,1531,1532,1533,1534,
      1548,1549,1550,1551,1552,1553,1554,1555,1556,1557,1558,1559,
      1573,1574,1575,1576,1577,1578,1579,1580,1581,1582,1583,1584,
      1598,1599,1600,1601,1602,1603,1604,1605,1606,1607,1608,1609,-1,-1,-1,-1,-1},
    //104
    //DDL[36[]=
    { 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840,
      852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862,
      874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884,
      896, 897, 898, 899, 900, 901, 902, 903, 904, 905, 906,
      1623,1624,1625,1626,1627,1628,1629,1630,1631,1632,1633,1634,
      1648,1649,1650,1651,1652,1653,1654,1655,1656,1657,1658,1659,
      1673,1674,1675,1676,1677,1678,1679,1680,1681,1682,1682,1684,
      1698,1699,1700,1701,1702,1703,1704,1705,1706,1707,1708,1709,
      1723,1724,1725,1726,1727,1728,1729,1730,1731,1732,1733,1734,-1,-1,-1,-1,-1},
    //104
    //DDL[37][]=
    { 918, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928,
      940, 941, 942, 943, 944, 945, 946, 947, 948, 949, 950,
      962, 963, 964, 965, 966, 967, 968, 969, 970, 971, 972,
      984, 985, 986, 987, 988, 989, 990, 991, 992, 993, 994,
      1748,1749,1750,1751,1752,1753,1754,1755,1756,1757,1758,1759,
      1773,1774,1775,1776,1777,1778,1779,1780,1781,1782,1783,1784,
      1798,1799,1800,1801,1802,1803,1804,1805,1806,1807,1808,1809,
      1823,1824,1825,1826,1827,1828,1829,1830,1831,1832,1833,1834,
      1848,1849,1850,1851,1852,1853,1854,1855,1856,1857,1858,1859,-1,-1,-1,-1,-1},
    //103
    //DDL[38][]=
    {1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,
     1028,1029,1030,1031,1032,1033,1034,1035,1036,1037,1038,
     1050,1051,1052,1053,1054,1055,1056,1057,1058,1059,1060,
     1072,1073,1074,1075,1076,1077,1078,1079,1080,1081,1082,
     1094,1095,1096,1097,1098,1099,1100,1101,1102,1103,1104,
     1873,1874,1875,1876,1877,1878,1879,1880,1881,1882,1883,1884,
     1898,1899,1900,1901,1902,1903,1904,1905,1906,1907,1908,1909,
     1923,1924,1925,1926,1927,1928,1929,1930,1931,1932,1933,1934,
     1948,1949,1950,1951,1952,1953,1954,1955,1956,1957,1958,1959,-1,-1,-1,-1,-1,-1},
    //104
    //DDL[39][]=
    {1116,1117,1118,1119,1120,1121,1122,1123,1124,1125,1126,
     1138,1139,1140,1141,1142,1143,1144,1145,1146,1147,1148,
     1160,1161,1162,1163,1164,1165,1166,1167,1168,1169,1170,
     1182,1183,1184,1185,1186,1187,1188,1189,1190,1191,1192,
     1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,
     1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,
     2023,2024,2025,2026,2027,2028,2029,2030,2031,2032,2033,2034,
     2048,2049,2050,2051,2052,2053,2054,2055,2056,2057,2058,2059,
     2073,2074,2075,2076,2077,2078,2079,2080,2081,2082,2083,2084,-1,-1,-1,-1,-1},
    //109
    //DDL[40][]=
    { 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521,
      533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543,
      1215,1216,1217,1218,1219,1220,1221,1222,1223,1224,1225,
      1237,1238,1239,1240,1241,1242,1243,1244,1245,1246,1247,
      1260,1261,1262,1263,1264,1265,1266,1267,1268,1269,1270,1271,1272,
      2110,2111,2112,2113,2114,2115,2116,2117,2118,2119,2120,2121,2122,
      2135,2136,2137,2138,2139,2140,2141,2142,2143,2144,2145,2146,2147,
      2160,2161,2162,2163,2164,2165,2166,2167,2168,2169,2170,2171,2172,
      2185,2186,2187,2188,2189,2190,2191,2192,2193,2194,2195,2196,2197},
    //109
    //DDL[41][]=
    { 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565,
      577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587,
      599, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609,
      621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631,
      1285,1286,1287,1288,1289,1290,1291,1292,1293,1294,1295,1296,1297,
      1310,1311,1312,1313,1314,1315,1316,1317,1318,1319,1320,1321,1322,
      1335,1336,1337,1338,1339,1340,1341,1342,1443,1344,1345,1346,1347,
      1360,1361,1362,1363,1364,1365,1366,1367,1368,1369,1370,1371,1372,
      1385,1386,1387,1388,1389,1390,1391,1392,1393,1394,1395,1396,1397},
    //107
    //DDL[42][]=
    { 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653,
      665, 666, 667, 668, 669, 670, 671, 672, 673, 674, 675,
      687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697,
      709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719,
      731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741,
      1410,1411,1412,1413,1414,1415,1416,1417,1418,1419,1420,1421,1422,
      1435,1436,1437,1438,1439,1440,1441,1442,1443,1444,1445,1446,1447,
      1460,1461,1462,1463,1464,1465,1466,1467,1468,1469,1470,1471,1472,
      1485,1486,1487,1488,1489,1490,1491,1492,1493,1494,1495,1496,1497,-1,-1},
    //109
    //DDL[43][]=
    { 753, 754, 755, 756, 757, 758, 759, 760, 761, 762, 763,
      775, 776, 777, 778, 779, 780, 781, 782, 783, 784, 785,
      797, 798, 799, 800, 801, 802, 803, 804, 805, 806, 807,
      819, 820, 821, 822, 823, 824, 825, 826, 827, 828, 829,
      1510,1511,1512,1513,1514,1515,1516,1517,1518,1519,1520,1521,1522,
      1535,1536,1537,1538,1539,1540,1541,1542,1543,1544,1545,1546,1547,
      1560,1561,1562,1563,1564,1565,1566,1567,1568,1569,1570,1571,1572,
      1585,1586,1587,1588,1589,1590,1591,1592,1593,1584,1595,1596,1597,
      1610,1611,1612,1613,1614,1615,1616,1617,1618,1619,1620,1621,1622},
    //109
    //DDL[44][]=
    { 841, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851,
      863, 864, 865, 866, 867, 868, 869, 870, 871, 872, 873,
      885, 886, 887, 888, 889, 890, 891, 892, 893, 894, 895,
      907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917,
      1635,1636,1637,1638,1639,1640,1641,1642,1643,1644,1645,1646,1647,
      1660,1661,1662,1663,1664,1665,1666,1667,1668,1669,1670,1671,1672,
      1685,1686,1687,1688,1689,1690,1691,1692,1693,1694,1695,1696,1697,
      1710,1711,1712,1713,1714,1715,1716,1717,1718,1719,1720,1721,1722,
      1735,1736,1737,1738,1739,1740,1741,1742,1743,1744,1745,1746,1747},
    //109
    //DDL[45][]=
    {929, 930, 931, 932, 933, 934, 935, 936, 937, 938, 939,
     951, 952, 953, 954, 955, 956, 957, 958, 959, 960, 961,
     973, 974, 975, 976, 977, 978, 979, 980, 981, 982, 983,
     995, 996, 997, 998, 999,1000,1001,1002,1003,1004,1005,
     1760,1761,1762,1763,1764,1765,1766,1767,1768,1769,1770,1771,1772,
     1785,1786,1787,1788,1789,1790,1791,1792,1793,1794,1795,1796,1797,
     1810,1811,1812,1813,1814,1815,1816,1817,1818,1819,1820,1821,1822,
     1835,1836,1837,1838,1839,1840,1841,1842,1843,1844,1845,1846,1847,
     1860,1861,1862,1863,1864,1865,1866,1867,1868,1869,1870,1871,1872},
    //109
    //DDL[46][]=
    {1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,
     1039,1040,1041,1042,1043,1044,1045,1046,1047,1048,1049,
     1061,1062,1063,1064,1065,1066,1067,1068,1069,1070,1071,
     1083,1084,1085,1086,1087,1088,1089,1090,1091,1092,1093,
     1105,1106,1107,1108,1109,1110,1111,1112,1113,1114,1115,
     1885,1886,1887,1888,1889,1890,1891,1892,1893,1894,1895,1896,1897,
     1910,1911,1912,1913,1914,1915,1916,1917,1918,1919,1920,1921,1922,
     1935,1936,1937,1938,1939,1940,1941,1942,1943,1944,1945,1946,1947,
     1960,1961,1962,1963,1964,1965,1966,1967,1968,1969,1970,1971,1972},
    //109
    //DDL[47][]=
    {1127,1128,1129,1130,1131,1132,1133,1134,1135,1136,1137,
     1149,1150,1151,1152,1153,1154,1155,1156,1157,1158,1159,
     1171,1172,1173,1174,1175,1176,1177,1178,1179,1180,1181,
     1193,1194,1195,1196,1197,1198,1199,1200,1201,1202,1203,
     1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,
     2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,
     2035,2036,2037,2038,2039,2040,2041,2042,2043,2044,2045,2046,2047,
     2060,2061,2062,2063,2064,2065,2066,2067,2068,2069,2070,2071,2072,
     2085,2086,2087,2088,2089,2090,2091,2092,2093,2094,2095,2096,2097} };


  Int_t DDLPerFile=DDLNumber/LDCsNumber;
  if(20%LDCsNumber)DDLPerFile++;
  cout<<"Number of DDL per File: "<<DDLPerFile<<endl;
  Int_t subd=2;          //SSD
  Int_t CountDDL=0;
  Int_t SliceNumber=1;
  const Int_t size=1536;//768*2 Number of stripe * number of sides(N and P)
  ULong_t buf[size];      
  fIndex=-1;
  Int_t nbytes = 0; 
  TClonesArray *ITSdigits  = ITS->DigitsAddress(subd);
  char fileName[15];

  ULong_t MiniHeaderPosition=0; //variable used to store the position of the Mini Header inside the file
  ofstream outfile;             // logical name of the output file 
  Int_t Flag=0;                 // 0==> Uncompressed data 1==>Compressed Data
  Int_t Detector=3;             // 1==>ITS (Pixel) 2==>ITS (Drift) 3==>ITS (Strip) 0==>TPC ......
  ULong_t Size=0;               // Size of the data block that follows the mini header
  Int_t MagicWord=0x123456;     // Magic word used to distinguish between data and garbage
  sprintf(fileName,"SSDslice%d",SliceNumber); //The name of  the output file. There are as many slides as the number of LDC
  outfile.open(fileName,ios::binary);
  ULong_t MiniHeader[3];
  Int_t MiniHeaderSize=sizeof(ULong_t)*3;
  Int_t Version=1;              //Version of the mini header 
  //loop over DDLs
  
  for(Int_t i=0;i<DDLNumber;i++){
    CountDDL++;
    //write Dummy MINI HEADER
    MiniHeader[0]=Size;
    PackWord(MiniHeader[1],MagicWord,8,31);
    PackWord(MiniHeader[1],Detector,0,7);
    PackWord(MiniHeader[2],i,16,31);
    PackWord(MiniHeader[2],Flag,8,15);
    PackWord(MiniHeader[2],Version,0,7);
    MiniHeaderPosition=outfile.tellp();
    outfile.write((char*)(MiniHeader),MiniHeaderSize);
    
    //Loops over Modules of a particular DDL
    for (Int_t mod=0; mod<ModulePerDDL; mod++){
      if(SSDMap[i][mod]!=-1){
	ITS->ResetDigits();
	nbytes += TD->GetEvent(SSDMap[i][mod]);
	//For each Module, buf contains the array of data words in Binary format	  
	//fIndex gives the number of 32 bits words in the buffer for each module
	GetDigitsSSD(ITSdigits,mod,buf);
	outfile.write((char *)buf,((fIndex+1)*sizeof(ULong_t)));
	for(Int_t i=0;i<(fIndex+1);i++){
	  buf[i]=0;
	}//end for
	fIndex=-1;
      }//end if
    }//end for
    //Write REAL MINI HEADER
    ULong_t CurrentFilePosition=outfile.tellp();
    outfile.seekp(MiniHeaderPosition);
    Size=CurrentFilePosition-MiniHeaderPosition-MiniHeaderSize;
    
    outfile.write((char*)(&Size),sizeof(ULong_t));
    outfile.seekp(CurrentFilePosition);
    if(CountDDL==DDLPerFile){
      outfile.close();
      SliceNumber++;
      sprintf(fileName,"SSDslice%d",SliceNumber); 
      if(i!=(DDLNumber-1))
	outfile.open(fileName,ios::binary);
      CountDDL=0;
    }//end if
  }//end for
  outfile.close();
  return 0;  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliITSDDLRawData::RawDataSDD(AliITS *ITS,TTree *TD ,Int_t LDCsNumber){
  //Silicon Drift Detector
  const Int_t DDLNumber=12;       // Number of DDL in SPD
  const Int_t ModulePerDDL=22;    // Number of modules in each DDL 
  //Row    ==> DDLs
  //Column ==> Modules  
  Int_t SDDMap[DDLNumber][ModulePerDDL]= {{240,241,242,246,247,248,252,253,254,258,259,260,264,265,266,270,271,272,276,277,278,-1},
					  {243,244,245,249,250,251,255,256,257,261,262,263,267,268,269,273,274,275,279,280,281,-1},
					  {282,283,284,288,289,290,294,295,296,300,301,302,306,307,308,312,313,314,318,319,320,-1},
					  {285,286,287,291,292,293,297,298,299,303,304,305,309,310,311,315,316,317,321,322,323,-1},
					  {324,325,326,327,332,333,334,335,340,341,342,343,348,349,350,351,356,357,358,359,364,365},
					  {328,329,330,331,336,337,338,339,344,345,346,347,352,353,354,355,360,361,362,363,368,369},
					  {366,367,372,373,374,375,380,381,382,383,388,389,390,391,396,397,398,399,404,405,406,407},
					  {370,371,376,377,378,379,384,385,386,387,392,393,394,395,400,401,402,403,408,409,410,411},
					  {412,413,414,415,420,421,422,423,428,429,430,431,436,436,438,439,444,445,446,447,452,453},
					  {416,417,418,419,424,425,426,427,432,433,434,435,440,441,442,443,448,449,450,451,456,457},
					  {454,455,460,461,462,463,468,469,470,471,476,477,478,479,484,485,486,487,492,493,494,495},
					  {458,459,464,465,466,467,472,473,474,475,480,481,482,483,488,489,490,491,496,497,498,499}};
  
  Int_t DDLPerFile=DDLNumber/LDCsNumber;
  if(DDLNumber%LDCsNumber)DDLPerFile++;
  cout<<"Number of DDL per File: "<<DDLPerFile<<endl;
  Int_t subd=1;
  Int_t CountDDL=0;
  Int_t SliceNumber=1;
  const Int_t size=131072; //256*512
  ULong_t buf[size];      
  fIndex=-1;
  Int_t nbytes = 0; 
  TClonesArray *ITSdigits  = ITS->DigitsAddress(subd);
  char fileName[15];

  ULong_t MiniHeaderPosition=0; // variable used to store the position of the Mini Header inside a file
  ofstream outfile;             // logical name of the output file 
  Int_t Flag=0;                 // 0==> Uncompressed data 1==>Compressed Data
  Int_t Detector=2;             // 1==>ITS (Pixel) 2==>ITS (Drift) 3==>ITS (Strip) 0==>TPC ......
  ULong_t Size=0;               // Size of the data block that follows the mini header
  Int_t MagicWord=0x123456;     // Magic word used to distinguish between data and garbage
  sprintf(fileName,"SDDslice%d",SliceNumber); //The name of  the output file. There are as many slides as the number of LDC
  outfile.open(fileName,ios::binary);
  ULong_t MiniHeader[3];
  Int_t MiniHeaderSize=sizeof(ULong_t)*3;
  Int_t Version=1;             //Version of the mini header 
  //loop over DDLs
  for(Int_t i=0;i<DDLNumber;i++){
    CountDDL++;
    //write Dummy MINI HEADER
    MiniHeader[0]=Size;
    PackWord(MiniHeader[1],MagicWord,8,31);
    PackWord(MiniHeader[1],Detector,0,7);
    PackWord(MiniHeader[2],i,16,31);
    PackWord(MiniHeader[2],Flag,8,15);
    PackWord(MiniHeader[2],Version,0,7);
    MiniHeaderPosition=outfile.tellp();
    outfile.write((char*)(MiniHeader),MiniHeaderSize);
    //Loops over Modules of a particular DDL
    for (Int_t mod=0; mod<ModulePerDDL; mod++){
      if(SDDMap[i][mod]!=-1){
	ITS->ResetDigits();
	nbytes += TD->GetEvent(SDDMap[i][mod]);
	//For each Module, buf contains the array of data words in Binary format	  
	//fIndex gives the number of 32 bits words in the buffer for each module
	//	cout<<"MODULE NUMBER:"<<SDDMap[i][mod]<<endl;
	GetDigitsSDD(ITSdigits,mod,buf);
	outfile.write((char *)buf,((fIndex+1)*sizeof(ULong_t)));
	for(Int_t i=0;i<(fIndex+1);i++){
	  buf[i]=0;
	}//end for
	fIndex=-1;
      }//end if
    }//end for
    
    //Write REAL MINI HEADER
    ULong_t CurrentFilePosition=outfile.tellp();
    outfile.seekp(MiniHeaderPosition);
    Size=CurrentFilePosition-MiniHeaderPosition-MiniHeaderSize;
    
    outfile.write((char*)(&Size),sizeof(ULong_t));
    outfile.seekp(CurrentFilePosition);
    if(CountDDL==DDLPerFile){
      outfile.close();
      SliceNumber++;
      sprintf(fileName,"SDDslice%d",SliceNumber); 
      if(i!=(DDLNumber-1))
	outfile.open(fileName,ios::binary);
      CountDDL=0;
    }//end if
  }//end for
  outfile.close();
  return 0;  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AliITSDDLRawData::WriteChipHeader(Int_t ChipAddr,Int_t EventCnt,ULong_t &BaseWord){
  //cout<<"Chip: "<<ChipAddr<<" Half Stave module:"<<EventCnt<<endl;
  BaseWord=0;
  PackWord(BaseWord,ChipAddr,0,3);
  PackWord(BaseWord,EventCnt,4,10);
  PackWord(BaseWord,0x7,11,13);
  PackWord(BaseWord,0x1,14,15);
  return;
}//end WriteChipHeader

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AliITSDDLRawData::ReadChipHeader(Int_t &ChipAddr,Int_t &EventCnt,ULong_t BaseWord){
  ULong_t temp=0;
  UnpackWord(BaseWord,0,3,temp);
  ChipAddr=(Int_t)temp;
  UnpackWord(BaseWord,4,10,temp);
  EventCnt=(Int_t)temp;
  cout<<"Chip: "<<ChipAddr<<" Half Stave module:"<<EventCnt<<endl;
  return;
}//end ReadChipHeader

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void  AliITSDDLRawData::WriteChipTrailer(ULong_t *buf,Int_t ChipHitCount,ULong_t &BaseWord){
  //pixel fill word
  if((ChipHitCount%2)!=0){
    PackWord(BaseWord,0xFECD,0,15);
  }
  PackWord(BaseWord,ChipHitCount,16,28);
  PackWord(BaseWord,0x0,30,31);
  fIndex++;
  buf[fIndex]=BaseWord;
  BaseWord=0;
  return;
}//end WriteChipTrailer

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void  AliITSDDLRawData::ReadChipTrailer(Int_t &ChipHitCount,ULong_t BaseWord){
  ULong_t temp=0;
  UnpackWord(BaseWord,16,28,temp);
  ChipHitCount=(Int_t)temp;
  return;
}//end ReadChipTrailer

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void  AliITSDDLRawData::WriteHit(ULong_t *buf,Int_t RowAddr,Int_t HitAddr,ULong_t &BaseWord){
  if(!BaseWord){
    PackWord(BaseWord,HitAddr,0,4);
    PackWord(BaseWord,RowAddr,5,12);
    PackWord(BaseWord,2,14,15);
  }//end if
  else{
    PackWord(BaseWord,HitAddr,16,20);
    PackWord(BaseWord,RowAddr,21,28);
    PackWord(BaseWord,2,30,31);
    fIndex++;
    buf[fIndex]=BaseWord;
    BaseWord=0;
  }//end else
  return;
}//end WriteHit

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AliITSDDLRawData::TestFormat(){
  ifstream f;
  Int_t  LDCsNumber=2;
  ofstream ftxt("File2.txt");
  ULong_t Size=0;
  char filename[15];
  ULong_t DDLNumber=0;
  ULong_t MiniHeader[3];
  Int_t MiniHeaderSize=sizeof(ULong_t)*3;
  for(Int_t i=1;i<=LDCsNumber;i++){
    sprintf(filename,"SPDslice%d",i);  
    f.open(filename,ios::binary|ios::in);
    if(!f){exit(1);}
    //loop over the DDL block 
    //here the Mini Header is read
    while( f.read((char*)(MiniHeader),MiniHeaderSize)){
      //cout<<"Block Size: "<<Size<<endl;
      Size=MiniHeader[0];
      UnpackWord(MiniHeader[2],16,31,DDLNumber);
      ftxt<<"DDL NUMBER:"<<DDLNumber<<endl;
      ULong_t Word=0;
      ULong_t Code=0;
      ULong_t Decoded1,Decoded2=0;
      for(ULong_t j=0;j<(Size/4);j++){
	f.read((char*)(&Word),sizeof(Word)); //32 bits word
	Code=0;
	UnpackWord(Word,14,15,Code);
	DecodeWord(Code,Word,0,Decoded1,Decoded2);
	switch (Code){
	case 0://trailer
	  ftxt<<"Number of Hit:"<<Decoded1<<endl;
	  break;
	case 1://header
	  ftxt<<"Half Stave Number:"<<Decoded1<<" Chip Number:"<<Decoded2<<endl;
	  break;
	case 2://hit
	  ftxt<<"Row:"<<Decoded1<<" Column:"<<Decoded2<<endl;
	  break;
	case 3://fill word
	  break;
	}//end switch
	Code=0;
	UnpackWord(Word,30,31,Code);
	DecodeWord(Code,Word,1,Decoded1,Decoded2);
	switch (Code){
	case 0://trailer
	  ftxt<<"Number of Hit:"<<Decoded1<<endl;
	  break;
	case 1://header
	  ftxt<<"Half Stave Number:"<<Decoded1<<" Chip Number:"<<Decoded2<<endl;
	  break;
	case 2://hit
	  ftxt<<"Row:"<<Decoded1<<" Column:"<<Decoded2<<endl;
	  break;
	case 3://fill word
	  break;
	}//end switch
      }//end for
    }//end while
    f.close();
  }//end for
  ftxt.close();
  return;
}

void AliITSDDLRawData::DecodeWord(ULong_t Code,ULong_t BaseWord,Int_t FirstHalf,ULong_t &Decoded1,ULong_t &Decoded2){
  //FirstHalf=0 ==>bits from 0 to 15
  //FirstHalf=1 ==>bits from 16 to 31
  if(!FirstHalf){
    switch (Code){
    case 0://trailer
      UnpackWord(BaseWord,0,12,Decoded1);
      break;
    case 1://header
      UnpackWord(BaseWord,4,10,Decoded1);
      UnpackWord(BaseWord,0,3,Decoded2);
      break;
    case 2://hit
      UnpackWord(BaseWord,5,12,Decoded1);
      UnpackWord(BaseWord,0,4,Decoded2);
      break;//fill word
    case 3:
      UnpackWord(BaseWord,0,13,Decoded1);
      break;
    }//end switch
  }
  else{
    switch (Code){
    case 0://trailer
      UnpackWord(BaseWord,16,28,Decoded1);
      break;
    case 1://header
      UnpackWord(BaseWord,20,26,Decoded1);
      UnpackWord(BaseWord,16,19,Decoded2);
      break;
    case 2://hit
      UnpackWord(BaseWord,21,28,Decoded1);
      UnpackWord(BaseWord,16,20,Decoded2);
      break;
    case 3://fill word
      UnpackWord(BaseWord,16,29,Decoded1);
      break;
    }//end switch
  }
  return;
}
