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

#include <TObjArray.h>
#include "Riostream.h"
#include <math.h>
#include "AliTPCCompression.h"
#include "AliTPCBuffer160.h"
#include "AliTPCHuffman.h"

ClassImp(AliTPCCompression)
//////////////////////////////////////////////////////////////////////////////////////////////////
AliTPCCompression::AliTPCCompression(){
  fDimBuffer=sizeof(ULong_t)*8;
  fFreeBitsBuffer=fDimBuffer;
  fReadBits=0;
  fPos=0;
  fBuffer=0;
  fVerbose=0;
  fFillWords=0;
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////
AliTPCCompression::AliTPCCompression(const AliTPCCompression &source){
  this->fDimBuffer=source.fDimBuffer;
  this->fFreeBitsBuffer=source.fFreeBitsBuffer;
  this->fReadBits=source.fReadBits;
  this->fPos=source.fPos;
  this->fBuffer=source.fBuffer;
  this->fVerbose=source.fVerbose;
  this->fFillWords=source.fFillWords;
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////
AliTPCCompression&  AliTPCCompression::operator=(const AliTPCCompression &source){
  this->fDimBuffer=source.fDimBuffer;
  this->fFreeBitsBuffer=source.fFreeBitsBuffer;
  this->fReadBits=source.fReadBits;
  this->fPos=source.fPos;
  this->fBuffer=source.fBuffer;
  this->fVerbose=source.fVerbose;
  this->fFillWords=source.fFillWords;
  return *this;
} 
//////////////////////////////////////////////////////////////////////////////////////////////////
void AliTPCCompression::NextTable(Int_t Val,Int_t &NextTableType,Int_t &BunchLen,Int_t &Count){
  /*
    Table index:
    0==> Bunch length value     
    1==> Time Bin value 
    2==> 1-samples bunch
    3==> Central samples
    4==> Border samples
  */  
  switch (NextTableType){
  case 0:{
    BunchLen=Val-2;
    NextTableType=1;
    break;
  }//end case 0
  case 1:{
    if (BunchLen==1)NextTableType=2;
    else{
      NextTableType=4;
      Count=1;
    }
    break;
  }//end case 1
  case 2:{
    NextTableType=0;
    break;
  }//end case 2
  case 3:{
    Count++;
    if (Count==(BunchLen-1)){
      NextTableType=4;
    }
    break;
  }//end case 3
  case 4:{
    if (Count==1){
      if (BunchLen>2)
	NextTableType=3;
      else
	Count++;
    }
    else
      NextTableType=0;
    break;
  }//end case 4
  }//end switch
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliTPCCompression::FillTables(const char* fSource,AliTPCHTable* table[],const Int_t NumTables){
  //This method is used to compute the frequencies of the symbols in the source file
  AliTPCBuffer160 Buff(fSource,0);
  ULong_t CountWords=0;
  ULong_t CountTrailer=0;
  Int_t NumWords,PadNum,RowNum,SecNum=0;
  Int_t Value=0;
  ULong_t Stat[5]={0,0,0,0,0};
  Int_t EndFill=0;
  Int_t End=1;
  while(Buff.ReadTrailerBackward(NumWords,PadNum,RowNum,SecNum) !=-1 ){
    if(End){
      EndFill=Buff.GetFillWordsNum();
      End=0;
    }//endif
    CountTrailer++;
    if (NumWords%4){
      fFillWords+=4-NumWords%4;
      for(Int_t j=0;j<(4-NumWords%4);j++){
	Value=Buff.GetNextBackWord();
      }//end for
    }//end if
    
    Int_t Packet[1024];
    Int_t TimePos[345];
    Int_t Tp=0;
    for(Int_t i=0;i<345;i++)TimePos[i]=0;
    for(Int_t i=0;i<1024;i++)Packet[i]=0;
    
    Int_t NextTableType=0;
    Int_t BunchLen=0;
    Int_t Count=0;
    for(Int_t i=0;i<NumWords;i++){
      Value=Buff.GetNextBackWord();
      Packet[i]=Value;
      if(NextTableType==1){
	TimePos[Tp]=i;
	Tp++;
      }
      NextTable(Value,NextTableType,BunchLen,Count);
    }//end for
    //computing the Time gap between two bunches
    Int_t temp=0;
    Tp--;
    Int_t PreviousTime=Packet[TimePos[Tp]];
    for(Int_t i=Tp-1;i>=0;i--){
      Int_t TimPos=TimePos[i];
      Int_t BunchLen=Packet[TimPos-1]-2;
      temp=Packet[TimPos];
      Packet[TimPos]=Packet[TimPos]-PreviousTime-BunchLen;
      PreviousTime=temp;
    }
    NextTableType=0;
    Count=0;
    BunchLen=0;
    for(Int_t i=0;i<NumWords;i++){
      Value=Packet[i];
      table[NextTableType]->SetFrequency(Value);
      Stat[NextTableType]++;
      NextTable(Value,NextTableType,BunchLen,Count);
      CountWords++;
    }//end for
  }//end while
  cout<<"Number of words:       "<<CountWords<<endl;
  cout<<"Number of trailers:    "<<CountTrailer<<endl;
  cout<<"Number of fill words   "<<fFillWords+EndFill<<endl;
  cout<<"Total number of words: "<<CountWords+CountTrailer*4+fFillWords<<endl;
  //STATISTICS  
  stat.open("Statistics");
  stat<<"Number of words:..........................................."<<CountWords<<endl;
  stat<<"Number of trailers (4 10 bits words in each one)..........."<<CountTrailer<<endl;
  stat<<"Number of fill words:......................................"<<fFillWords+EndFill<<endl;
  stat<<"Total number of words:....................................."<<CountWords+CountTrailer*4+fFillWords+EndFill<<endl;
  stat<<"-----------------------------------------"<<endl;
  stat<<"Number of Bunches............."<<Stat[0]<<endl;
  stat<<"Number of Time bin............"<<Stat[1]<<endl;
  stat<<"Number of One Samples Bunch..."<<Stat[2]<<endl;
  stat<<"Number of Central Samples....."<<Stat[3]<<endl;
  stat<<"Number of Border Samples......"<<Stat[4]<<endl;
  stat<<"-----------------------------------------"<<endl;
  ULong_t FileDimension=(ULong_t)ceil(double((CountTrailer*4+CountWords+fFillWords+EndFill)*10/8));
  stat<<"Total file Size in bytes.."<<FileDimension<<endl;
  Double_t Percentage=ceil((fFillWords+EndFill)*125)/FileDimension;
  stat<<"Fill Words................"<<(ULong_t)ceil((fFillWords+EndFill)*10/8)<<" bytes   "<<Percentage<<"%"<<endl;  
  Percentage=(Double_t)CountTrailer*500/FileDimension;
  stat<<"Trailer..................."<<CountTrailer*5<<" bytes   "<<Percentage<<"%"<<endl;

  Percentage=(Double_t)((Stat[0]+Stat[1]+Stat[2]+Stat[3]+Stat[4])) *125/FileDimension;
  stat<<"Data......................"<<(ULong_t)ceil((Stat[0]+Stat[1]+Stat[2]+Stat[3]+Stat[4])*10/8)<<" bytes   "<<Percentage<<"%"<<endl;

  Percentage=(Double_t)(Stat[0]*125)/FileDimension;
  stat<<"Bunch....................."<<(ULong_t)ceil(Stat[0]*10/8)<<" bytes  "<<Percentage<<"%"<<endl;  //  
  Percentage=(Double_t)(Stat[1]*125)/FileDimension;
  stat<<"Time......................"<<(ULong_t)ceil(Stat[1]*10/8)<<" bytes  "<<Percentage<<"%"<<endl;  //  


  Percentage=(Double_t)((Stat[2]+Stat[3]+Stat[4])) *125/FileDimension;
  stat<<"Amplitude values.........."<<(ULong_t)ceil((Stat[2]+Stat[3]+Stat[4])*10/8)<<" bytes  "<<Percentage<<"%"<<endl;
  Percentage=(Double_t)(Stat[2]*125)/FileDimension;
  stat<<"     One Samples..............."<<(ULong_t)ceil(Stat[2]*10/8)<<" bytes  "<<Percentage<<"%"<<endl;  //  
  Percentage=(Double_t)(Stat[3]*125)/FileDimension;
  stat<<"     Central Samples..........."<<(ULong_t)ceil(Stat[3]*10/8)<<" bytes  "<<Percentage<<"%"<<endl;  //  
  Percentage=(Double_t)(Stat[4]*125)/FileDimension;
  stat<<"     Border Samples............"<<(ULong_t)ceil(Stat[4]*10/8)<<" bytes  "<<Percentage<<"%"<<endl;  //  
  stat.close();
  return 0;
}
////////////////////////////////////////////////////////////////////////////////////////
Int_t AliTPCCompression::StoreTables(AliTPCHTable* table[],const Int_t NumTable){
  char filename[15];
  ofstream fTable;
  for(Int_t k=0;k<NumTable;k++){
    sprintf(filename,"Table%d.dat",k); 
    fTable.open(filename,ios::binary);
    Int_t dim=table[k]->Size();
    //Table dimension is written into a file
    fTable.write((char*)(&dim),sizeof(Int_t));
    //One table is written into a file
    for(Int_t i=0;i<dim;i++){
      UChar_t CodeLen=table[k]->CodeLen()[i];
      //      ULong_t Code=(ULong_t)table[k]->Code()[i];
      Double_t Code=table[k]->Code()[i];
      fTable.write((char*)(&CodeLen),sizeof(UChar_t));
      //fTable.write((char*)(&Code),sizeof(ULong_t));
      fTable.write((char*)(&Code),sizeof(Double_t));
    } //end for
    fTable.close();
  }//end for
  return 0;
}
////////////////////////////////////////////////////////////////////////////////////////
Int_t AliTPCCompression::CreateTables(const char* fSource,const Int_t NumTables){
    Int_t n=10;// 10 bits per symbol 
  /*
    Table index:
    0==> Bunch length values     
    1==> Time Bin values 
    2==> 1-samples bunch
    3==> Central samples
    4==> Border samples
  */
  AliTPCHTable ** table = new AliTPCHTable*[NumTables];
  //The table is inizialized with the rigth number of rows 
  for(Int_t i=0;i<NumTables;i++){table[i]=new  AliTPCHTable((Int_t)(pow(2,n)));}
  //The frequencies are calculated and the tables are filled
  if (fVerbose)
    cout<<"Filling tables...\n";
  //The get the frequencies 
  FillTables(fSource,table,NumTables);

  //This part will be used in the table optimization phase
  /*
  for(Int_t i=0;i<NumTables;i++){
    table[i]->CompleteTable(i);
  }
  */
  if(fVerbose){
    cout<<"Entropy of Bunch length table........."<<table[0]->GetEntropy()<<endl;
    cout<<"Entropy of Time bin table............."<<table[1]->GetEntropy()<<endl;
    cout<<"Entropy of one Sample bunch table....."<<table[2]->GetEntropy()<<endl;
    cout<<"Entropy of Central Sample table......."<<table[3]->GetEntropy()<<endl;
    cout<<"Entropy Border Samples table.........."<<table[4]->GetEntropy()<<endl;
  }
  stat.open("Statistics",ios::app);
  stat<<endl;
  stat<<"----------------- ENTROPY for castomized tables --------------------------"<<endl;
  stat<<"Entropy of Bunch length table......."<<table[0]->GetEntropy()<<endl;
  stat<<"Entropy of Time bin table..........."<<table[1]->GetEntropy()<<endl;
  stat<<"Entropy of one Sample bunch table..."<<table[2]->GetEntropy()<<endl;
  stat<<"Entropy of Central Sample table....."<<table[3]->GetEntropy()<<endl;
  stat<<"Entropy Border Samples table........"<<table[4]->GetEntropy()<<endl;
  stat.close();
 
  if (fVerbose)
    cout<<"Tables filled \n";
  //Tables are saved in a sequence of text file and using the macro Histo.C is it possible to get
  //a series of histograms rappresenting the frequency distribution
  table[0]->StoreFrequencies("BunchLenFreq.txt");
  table[1]->StoreFrequencies("TimeFreq.txt");
  table[2]->StoreFrequencies("Sample1Freq.txt");
  table[3]->StoreFrequencies("SCentralFreq.txt");
  table[4]->StoreFrequencies("SBorderFreq.txt");
  if (fVerbose)
    cout<<"Creating Tables..\n";
  //One Huffman tree is created for each table starting from the frequencies of the symbols
  for(Int_t i=0;i<NumTables;i++){
    table[i]->BuildHTable();
    if (fVerbose==2){
      cout<<"Number of elements inside the table:"<<table[i]->GetWordsNumber();
      switch(i){
      case 0:{
	cout<<" (Bunch Length)"<<endl;
	break;
      }
      case 1:{
	cout<<" (Time Bin)"<<endl;
	break;
      }
      case 2:{
	cout<<" (1 Samples Bunch)"<<endl;
	break;
      }
      case 3:{
	cout<<" (Central Samples)"<<endl;
	break;
      }
      case 4:{
	cout<<" (Border Samples)"<<endl;
	break;
      }
      }//end switch
      table[i]->PrintTable();
    }
  }
  //The tables are saved ad binary files
  StoreTables(table,NumTables);
  //The tables stored in memory are deleted; 
  for(Int_t i=0;i<NumTables;i++)delete table[i];
  delete [] table;
  return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Int_t AliTPCCompression::RetrieveTables(AliTPCHTable* table[],Int_t NumTable){
  if (fVerbose)
    cout<<"Retrieving tables from files \n";
  //  ULong_t Code;
  Double_t Code;
  UChar_t CodeLen;
  ifstream fTable;  
  char filename[15];
  //The following for loop is used to generate the Huffman trees acording to the tables
  for(Int_t k=0;k<NumTable;k++){
    Int_t Dim;//this variable contains the table dimension
    sprintf(filename,"Table%d.dat",k); 
    fTable.open(filename,ios::binary);
    fTable.read((char*)(&Dim),sizeof(Int_t));
    if (fVerbose)
      cout<<"Table dimension: "<<Dim<<endl;
    table[k]=new AliTPCHTable(Dim);
    for(Int_t i=0;i<Dim;i++){
      fTable.read((char*)(&CodeLen),sizeof(UChar_t));
      table[k]->SetCodeLen(CodeLen,i);
      //      fTable.read((char*)(&Code),sizeof(ULong_t));
      fTable.read((char*)(&Code),sizeof(Double_t));
      table[k]->SetCode(Mirror((ULong_t)Code,CodeLen),i);
    }//end for 
    fTable.close();
  }//end for 
  if (fVerbose)
    cout<<"Trees generated \n";
  //At this point the trees are been built
  return 0;
}
////////////////////////////////////////////////////////////////////////////////////////
/*                               COMPRESSION                                          */
////////////////////////////////////////////////////////////////////////////////////////

void AliTPCCompression::StoreValue(ULong_t val,UChar_t len){
  if (len<=fFreeBitsBuffer){           // val is not splitted in two buffer
    fFreeBitsBuffer-=len;
    fBuffer=fBuffer<<len;
    fBuffer=fBuffer|val;    
    if(!fFreeBitsBuffer){              // if the buffer is full it is written into a file 
      f.write((char*)(&fBuffer),sizeof(ULong_t));	
      fFreeBitsBuffer=fDimBuffer;
      fBuffer=0;
    }
  }//end if
  else{                               //val has to be splitted in two buffers
    fBuffer=fBuffer<<fFreeBitsBuffer;
    ULong_t temp;
    temp=val;
    temp=temp>>(len-fFreeBitsBuffer);
    fBuffer=fBuffer|temp;
    f.write((char*)(&fBuffer),sizeof(ULong_t));
    fFreeBitsBuffer=fDimBuffer-(len-fFreeBitsBuffer);
    val=val<<fFreeBitsBuffer;
    val=val>>fFreeBitsBuffer;
    fBuffer=val;
  }//end else
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////
void AliTPCCompression::Flush(){
  //The last buffen cannot be completely full
  if(fFreeBitsBuffer<fDimBuffer){
    fBuffer=fBuffer<<fFreeBitsBuffer;
    f.write((char*)(&fBuffer),sizeof(ULong_t));	 
  }//end if
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////
ULong_t AliTPCCompression::Mirror(ULong_t val,UChar_t len){
  ULong_t specular=0;
  ULong_t Mask=0x1;
  ULong_t bit;
  for(Int_t i=0;i<len;i++){
    bit=val&Mask;
    bit=bit>>i;
    specular=specular<<1;
    specular=specular|bit;
    Mask=Mask<<1;
  }
  return specular;
}
//////////////////////////////////////////////////////////////////////////////////////////////////
Int_t AliTPCCompression::CompressData(AliTPCHTable* table[],Int_t NumTable,const char* fSource,const char* fDest){
  cout<<" COMPRESSION "<<endl;
  cout<<"compression of the file "<<fSource<<" Output File: "<<fDest<<endl;
  //the output file is open
  f.open(fDest,ios::binary|ios::out);
  //Tables are written into the output file
  for(Int_t k=0;k<NumTable;k++){
    Int_t dim=table[k]->Size();
    //Table dimension is written into a file
    f.write((char*)(&dim),sizeof(Int_t));
    //One table is written into a file
    for(Int_t i=0;i<dim;i++){
      UChar_t CodeLen=table[k]->CodeLen()[i];
      ULong_t Code=(ULong_t)table[k]->Code()[i];
      f.write((char*)(&CodeLen),sizeof(UChar_t));
      f.write((char*)(&Code),sizeof(ULong_t));
    } //end for
  }//end for

  // Source file is open
  AliTPCBuffer160 Buff(fSource,0);
  //coded words are written into the output file
  Int_t NumWords,PadNum,RowNum,SecNum=0;
  ULong_t StoredWords=0;
  Int_t Value=0;
  ULong_t NumPacket=0;
  while(Buff.ReadTrailerBackward(NumWords,PadNum,RowNum,SecNum) !=-1 ){
    NumPacket++;
    if (NumWords%4){
      for(Int_t j=0;j<(4-NumWords%4);j++){
	Value=Buff.GetNextBackWord();
      }//end for
    }//end if

    Int_t Packet[1024];
    Int_t TimePos[345];
    Int_t Tp=0;
    for(Int_t i=0;i<345;i++)TimePos[i]=0;
    for(Int_t i=0;i<1024;i++)Packet[i]=0;

    Int_t NextTableType=0;
    Int_t BunchLen=0;
    Int_t Count=0;
    for(Int_t i=0;i<NumWords;i++){
      Value=Buff.GetNextBackWord();
      Packet[i]=Value;
      if(NextTableType==1){
	TimePos[Tp]=i;
	Tp++;
      }
      NextTable(Value,NextTableType,BunchLen,Count);
    }//end for
    //computing the Time gap between two bunches
    Int_t temp=0;
    Tp--;
    Int_t PreviousTime=Packet[TimePos[Tp]];
    for(Int_t i=Tp-1;i>=0;i--){
      Int_t TimPos=TimePos[i];
      Int_t BunchLen=Packet[TimPos-1]-2;
      temp=Packet[TimPos];
      Packet[TimPos]=Packet[TimPos]-PreviousTime-BunchLen;
      PreviousTime=temp;
    }//end for
    NextTableType=0;
    Count=0;
    BunchLen=0;
    Int_t TimeBin=0;
    //All the words for one pad are compressed and stored in the compress file
    for(Int_t i=0;i<NumWords;i++){
      Value=Packet[i];
      if(NextTableType==1)TimeBin=Value;
      if(NextTableType>1){
	//	ULong_t val=(ULong_t)table[NextTableType]->Code()[Value];     // val is the code
	Double_t val=table[NextTableType]->Code()[Value];     // val is the code
	UChar_t len=table[NextTableType]->CodeLen()[Value];  // len is the length (number of bits)of val
	StoreValue(Mirror((ULong_t)val,len),len);
	StoredWords++;
      }//end if
      NextTable(Value,NextTableType,BunchLen,Count);
      if(NextTableType==0){
	//	ULong_t val=(ULong_t)table[1]->Code()[TimeBin];     // val is the code
	Double_t val=table[1]->Code()[TimeBin];     // val is the code
	UChar_t len=table[1]->CodeLen()[TimeBin];  // len is the length (number of bits)of val
	StoreValue(Mirror((ULong_t)val,len),len);
	//val=(ULong_t)table[NextTableType]->Code()[(BunchLen+2)];     // val is the code
	val=table[NextTableType]->Code()[(BunchLen+2)];     // val is the code
	len=table[NextTableType]->CodeLen()[(BunchLen+2)];  // len is the length (number of bits)of val
	StoreValue(Mirror((ULong_t)val,len),len);
	StoredWords+=2;
      }
    }//end for
    //Trailer
    StoreValue(NumWords,10);
    StoreValue(PadNum,10);
    StoreValue(RowNum,10);
    StoreValue(SecNum,9);
    StoreValue(1,1);
    StoredWords+=4;
  }//end  while
  StoreValue(NumPacket,32);
  cout<<"Number of strored packet: "<<NumPacket<<endl;
  StoreValue(1,1);
  //The last buffen cannot be completely full
  Flush();
  cout<<"Number of stored words: "<<StoredWords<<endl;
  f.close();
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Int_t AliTPCCompression::CompressDataOptTables(Int_t NumTable,const char* fSource,const char* fDest){
  if (fVerbose){
    cout<<" BackWord COMPRESSION "<<endl;
    cout<<"compression of the file "<<fSource<<" Output File: "<<fDest<<endl;
  }
  //Tables are read from the files (Each codeword has been "Mirrored")
  AliTPCHTable ** table = new AliTPCHTable*[NumTable];
  RetrieveTables(table,NumTable);
  //the output file is open
  f.open(fDest,ios::binary|ios::out);
  // Source file is open
  AliTPCBuffer160 Buff(fSource,0);
  //coded words are written into a file
  Int_t NumWords,PadNum,RowNum,SecNum=0;
  ULong_t  StoredWords=0;
  Int_t    Value=0;
  ULong_t  NumPacket=0;
  Double_t Stat[5]={0.,0.,0.,0.,0.};
  ULong_t  TrailerNumber=0;
  Double_t  NumElem[5]={0,0,0,0,0};
  Double_t FillWords=0.;
  stat.open("Statistics",ios::app);
  stat<<endl;
  stat<<"-------------------COMPRESSION STATISTICS----------"<<endl;
  Int_t End=1;
  while(Buff.ReadTrailerBackward(NumWords,PadNum,RowNum,SecNum) !=-1 ){
    if(End){
      FillWords=Buff.GetFillWordsNum();
      End=0;
    }//endif

    NumPacket++;
    if (NumWords%4){
      FillWords+=4-NumWords%4;
      for(Int_t j=0;j<(4-NumWords%4);j++){
	Value=Buff.GetNextBackWord();
      }//end for
    }//end if

    Int_t Packet[1024];
    Int_t TimePos[345];
    Int_t Tp=0;
    for(Int_t i=0;i<345;i++)TimePos[i]=0;
    for(Int_t i=0;i<1024;i++)Packet[i]=0;

    Int_t NextTableType=0;
    Int_t BunchLen=0;
    Int_t Count=0;
    for(Int_t i=0;i<NumWords;i++){
      Value=Buff.GetNextBackWord();
      Packet[i]=Value;
      if(NextTableType==1){
	TimePos[Tp]=i;
	Tp++;
      }
      NextTable(Value,NextTableType,BunchLen,Count);
    }//end for
    //computing the Time gap between two bunches
    Int_t temp=0;
    Tp--;
    Int_t PreviousTime=Packet[TimePos[Tp]];
    for(Int_t i=Tp-1;i>=0;i--){
      Int_t TimPos=TimePos[i];
      Int_t BunchLen=Packet[TimPos-1]-2;
      temp=Packet[TimPos];
      Packet[TimPos]=Packet[TimPos]-PreviousTime-BunchLen;
      PreviousTime=temp;
    }//end for

    NextTableType=0;
    Count=0;
    BunchLen=0;
    Int_t TimeBin=0;
    for(Int_t i=0;i<NumWords;i++){
      Value=Packet[i];
      if(NextTableType==1)TimeBin=Value;
      if(NextTableType>1){
	//ULong_t val=(ULong_t)table[NextTableType]->Code()[Value];     // val is the code
	Double_t val=table[NextTableType]->Code()[Value];     // val is the code
	UChar_t len=table[NextTableType]->CodeLen()[Value];  // len is the length (number of bits)of val
	Stat[NextTableType]+=len;
	NumElem[NextTableType]++;
	StoreValue((ULong_t)val,len);
	StoredWords++;
      }//end if
      NextTable(Value,NextTableType,BunchLen,Count);
      if(NextTableType==0){
	//	ULong_t val=(ULong_t)table[1]->Code()[TimeBin];     // val is the code
	Double_t val=table[1]->Code()[TimeBin];     // val is the code
	UChar_t len=table[1]->CodeLen()[TimeBin];  // len is the length (number of bits)of val
	Stat[1]+=len;
	NumElem[1]++;
	StoreValue((ULong_t)val,len);
	//	val=(ULong_t)table[NextTableType]->Code()[(BunchLen+2)];     // val is the code
	val=table[NextTableType]->Code()[(BunchLen+2)];     // val is the code
	len=table[NextTableType]->CodeLen()[(BunchLen+2)];  // len is the length (number of bits)of val
	StoreValue((ULong_t)val,len);
	Stat[NextTableType]+=len;
	NumElem[NextTableType]++;
	StoredWords+=2;
      }
    }//end for
    //Trailer
    StoreValue(NumWords,10);
    StoreValue(PadNum,10);
    StoreValue(RowNum,10);
    StoreValue(SecNum,9);
    StoreValue(1,1);
    StoredWords+=4;
    TrailerNumber++;
  }//end  while
  StoreValue(NumPacket,32);
  if(fVerbose)
    cout<<"Number of strored packet: "<<NumPacket<<endl;
  StoreValue(1,1);
  //The last buffen cannot be completely full
  Flush();
  if(fVerbose)
    cout<<"Number of stored words: "<<StoredWords<<endl;
  f.close();
  //Tables are deleted
  for(Int_t i=0;i<NumTable;i++){
    delete table[i];
  }//end for
  delete [] table;
  Double_t dimension=(ULong_t)ceil((Stat[0]+Stat[1]+Stat[2]+Stat[3]+Stat[4])/8)+TrailerNumber*5;
  stat<<"Trailer Dimension in bytes......"<<TrailerNumber*5<<endl;
  stat<<"Data Dimension in bytes........."<<(ULong_t)ceil((Stat[0]+Stat[1]+Stat[2]+Stat[3]+Stat[4])/8)<<endl;
  stat<<"Compressed file dimension......."<<(ULong_t)dimension<<endl;
  /*
  stat<<(ULong_t)TrailerNumber<<endl;
  stat<<(ULong_t)FillWords<<endl;
  stat<<(ULong_t)NumElem[0]<<endl;
  stat<<(ULong_t)NumElem[1]<<endl;
  stat<<(ULong_t)NumElem[2]<<endl;
  stat<<(ULong_t)NumElem[3]<<endl;
  stat<<(ULong_t)NumElem[4]<<endl;
  */
  FillWords=(FillWords+NumElem[0]+NumElem[1]+NumElem[2]+NumElem[3]+NumElem[4]+TrailerNumber*4)*10/8;
  stat<<"Original file dimension........."<<(ULong_t)FillWords<<endl;

  Double_t ratio=(dimension/FillWords)*100;
  stat<<"Compression ratio (Compressed/Uncompressed)..."<<ratio<<"%"<<endl;
  stat<<endl;
  stat<<"Bunch length size in bytes......"<<(ULong_t)ceil(Stat[0]/8)<<" Comppression.."<<(Stat[0]/NumElem[0])*10<<"%"<<endl;
  
  stat<<"Time gap size in bytes.........."<<(ULong_t)ceil(Stat[1]/8)<<" Comppression.."<<(Stat[1]/NumElem[1])*10<<"%"<<endl;
  stat<<"Amplitude values in bytes......."<<(ULong_t)ceil((Stat[2]+Stat[3]+Stat[4])/8)<<" Comppression.."<<
    ((Stat[2]+Stat[3]+Stat[4])/(NumElem[2]+NumElem[3]+NumElem[4]))*10<<"%"<<endl;
  stat<<"     One Samples in bytes............"<<(ULong_t)ceil(Stat[2]/8)<<" Comppression.."<<(Stat[2]/NumElem[2])*10<<"%"<<endl;
  stat<<"     Central Samples size in bytes..."<<(ULong_t)ceil(Stat[3]/8)<<" Comppression.."<<(Stat[3]/NumElem[3])*10<<"%"<<endl;
  stat<<"     Border Samples size in bytes...."<<(ULong_t)ceil(Stat[4]/8)<<" Comppression.."<<(Stat[4]/NumElem[4])*10<<"%"<<endl;
  stat<<endl;
  stat<<"Average number of bits per word"<<endl;
  stat<<"Bunch length ......"<<Stat[0]/NumElem[0]<<endl;
  stat<<"Time gap .........."<<Stat[1]/NumElem[1]<<endl;
  stat<<"One Samples........"<<Stat[2]/NumElem[2]<<endl;
  stat<<"Central Samples ..."<<Stat[3]/NumElem[3]<<endl;
  stat<<"Border Samples....."<<Stat[4]/NumElem[4]<<endl;
  stat.close();
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
/*                               DECOMPRESSION                                        */
////////////////////////////////////////////////////////////////////////////////////////
void AliTPCCompression::CreateTrees(AliTPCHNode *RootNode[],const Int_t NumTables){
  //The first part of the compressed file cotains the tables
  //The following for loop is used to generate the Huffman trees acording to the tables
  if(fVerbose)
    cout<<"Creating the Huffman trees \n";
  AliTPCHNode *node=0;
  //  ULong_t Code;
  Double_t Code;
  UChar_t CodeLen;
  //loop over the numbero of tables
  for(Int_t k=0;k<NumTables;k++){
    RootNode[k]=new AliTPCHNode(); //RootNode is the root of the tree
    Int_t Dim;//this variable contains the table dimension
    f.read((char*)(&Dim),sizeof(Int_t));
    if (fVerbose)
      cout<<"Table dimension: "<<Dim<<endl;
    //loop over the words of a table
    for(Int_t i=0;i<Dim;i++){
      f.read((char*)(&CodeLen),sizeof(UChar_t));
      //f.read((char*)(&Code),sizeof(ULong_t));
      f.read((char*)(&Code),sizeof(Double_t));
      node=RootNode[k];
      for(Int_t j=1;j<=CodeLen;j++){
	ULong_t bit,val=0;
	val=(ULong_t)pow(2,CodeLen-j);
	bit=(ULong_t)Code&val; 
	AliTPCHNode *temp=node;
	if(bit){
	  node=node->GetRight();
	  if(!node){
	    node=new AliTPCHNode();
	    temp->SetRight(node);
	  }//end if
	}//end if
	else{
	  node=node->GetLeft();
	  if(!node){
	    node=new AliTPCHNode();
	    temp->SetLeft(node);
	  }//end if
	}//end else
      }//end for
      if(CodeLen){
	node->SetSymbol(i);
	node->SetFrequency(CodeLen);
      }//end if
      //cout<<node->GetSymbol()<<"  "<<(Int_t)node->GetFrequency()<<endl;
    }//end for 
  }//end for 
  if (fVerbose)
    cout<<"Trees generated \n";
  //At this point the trees are been built
}
//////////////////////////////////////////////////////////////////////////////////////////////////
void AliTPCCompression::CreateTreesFromFile(AliTPCHNode *RootNode[],const Int_t NumTables){
  if(fVerbose)
    cout<<"Creating the Huffman trees \n";
  AliTPCHNode *node=0;
  // ULong_t Code;
  Double_t Code;
  UChar_t CodeLen;
  ifstream fTable;  
  char filename[15];
  //The following for loop is used to generate the Huffman trees acording to the tables
  //loop over the tables
  for(Int_t k=0;k<NumTables;k++){
    RootNode[k]=new AliTPCHNode(); //RootNode is the root of the tree
    Int_t Dim=0;//this variable contains the table dimension
    sprintf(filename,"Table%d.dat",k); 
    fTable.open(filename,ios::binary);
    fTable.read((char*)(&Dim),sizeof(Int_t));
    if (fVerbose)
      cout<<"Table dimension: "<<Dim<<endl;
    //loop over the words of one table
    for(Int_t i=0;i<Dim;i++){
      fTable.read((char*)(&CodeLen),sizeof(UChar_t));
      //fTable.read((char*)(&Code),sizeof(ULong_t));
      fTable.read((char*)(&Code),sizeof(Double_t));
      node=RootNode[k];
      for(Int_t j=1;j<=CodeLen;j++){
	ULong_t bit,val=0;
	val=(ULong_t)pow(2,CodeLen-j);
	bit=(ULong_t)Code&val; 
	AliTPCHNode *temp=node;
	if(bit){
	  node=node->GetRight();
	  if(!node){
	    node=new AliTPCHNode();
	    temp->SetRight(node);
	  }//end if 
	}//end if
	else{
	  node=node->GetLeft();
	  if(!node){
	    node=new AliTPCHNode();
	    temp->SetLeft(node);
	  }//end if
	}//end else
      }//end for
      if(CodeLen){
	node->SetSymbol(i);
	node->SetFrequency(CodeLen);
      }//end if
    }//end for 
    fTable.close();
  }//end for 
  if (fVerbose)
    cout<<"Trees generated \n";
  //At this point the trees are been built
}
//////////////////////////////////////////////////////////////////////////////////////////////////
void AliTPCCompression::DeleteHuffmanTree(AliTPCHNode* node){
  //This function deletes all the nodes of an Huffman tree
  //In an Huffman tree any internal node has always two children
  if (node){
    DeleteHuffmanTree(node->GetLeft());
    DeleteHuffmanTree(node->GetRight());
    //    cout<<node->GetSymbol()<<"  "<<(Int_t)node->GetFrequency()<<endl;
    delete node;
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////////
void AliTPCCompression::VisitHuffmanTree(AliTPCHNode* node){
  //This function realizes an in order visit of a binary tree
  if (node){
    cout<<node->GetSymbol()<<" "<<node->GetFrequency()<<endl;
    VisitHuffmanTree(node->GetLeft());
    VisitHuffmanTree(node->GetRight());
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////////
ULong_t AliTPCCompression::ReadWord(Int_t NumberOfBit){
  ULong_t Result=0;
  ULong_t bit=0;
  for (Int_t i=0;i<NumberOfBit;i++){
    if (fReadBits==32){
      fPos-=sizeof(ULong_t);
      f.seekg(fPos);
      f.read((char*)(&fBuffer),sizeof(ULong_t));
      fReadBits=0;
    }//end if
    ULong_t mask=0;
    mask=(ULong_t)pow(2,fReadBits);
    bit=fBuffer&mask;
    bit=bit>>fReadBits;
    fReadBits++;
    bit=bit<<i;
    Result=Result|bit;
  }//end for
  return Result;
}
//////////////////////////////////////////////////////////////////////////////////////////////////
void AliTPCCompression::ReadTrailer(Int_t &WordsNumber,Int_t &PadNumber,Int_t &RowNumber,Int_t &SecNumber){
  ReadWord(1);
  SecNumber=ReadWord(9);
  RowNumber=ReadWord(10);
  PadNumber=ReadWord(10);
  WordsNumber=ReadWord(10);
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////
ULong_t AliTPCCompression::GetDecodedWord(AliTPCHNode* root){
  AliTPCHNode *node=root;
  ULong_t symbol=0;
  Bool_t decoded=0;
  while(!decoded){
    ULong_t bit=ReadWord(1);
    if(bit)
      node=node->GetRight();
    else
      node=node->GetLeft();
    if (!(node->GetLeft())){
      symbol=node->GetSymbol();
      decoded=1;
    }
  }//end while
  return symbol;
}
//////////////////////////////////////////////////////////////////////////////////////////////////
Int_t AliTPCCompression::DecompressData(Int_t NumTables,const char* fname,char* fDest){
  cout<<"   DECOMPRESSION:"<<endl;
  cout<<"Source File "<<fname<<" Destination File "<<fDest<<endl; 
  f.open(fname,ios::binary|ios::in);
  if(!f){cout<<"File doesn't exist\n";return -1;}
  AliTPCHNode ** RootNode = new AliTPCHNode*[NumTables];
  //Creation of the Huffman trees
  CreateTrees(RootNode,NumTables);
  //to go to the end of the file
  f.seekg(0,ios::end);
  //to get the file dimension in byte
  fPos=f.tellg();
  fPos-=sizeof(ULong_t);
  f.seekg(fPos);
  fReadBits=0;
  fBuffer=0;
  f.read((char*)(&fBuffer),sizeof(ULong_t));
  Int_t bit=0;
  ULong_t Mask=0x1;
  while(!bit){
    bit=fBuffer&Mask;
    Mask=Mask<<1;
    fReadBits++;
  }
  ULong_t PacketNumber=ReadWord(sizeof(ULong_t)*8);
  cout<<"Number of Packect: "<<PacketNumber<<endl;
  AliTPCBuffer160 BufferFile(fDest,1);
  ULong_t k=0;
  ULong_t WordsRead=0; //number of read coded words 
  while(k<PacketNumber){
    Int_t NumWords,PadNumber,RowNumber,SecNumber=0;
    ReadTrailer(NumWords,PadNumber,RowNumber,SecNumber);
    k++;
    WordsRead+=4;
    Int_t PreviousTime=-1;
    Int_t Time=0;
    Int_t NextTableType=0;
    Int_t BunchLen=0;
    Int_t Count=0;
    for(Int_t i=0;i<NumWords;i++){
      ULong_t symbol=GetDecodedWord(RootNode[NextTableType]);
      WordsRead++;
      //Time reconstruction
      if (NextTableType==1){
	if (PreviousTime!=-1){
	  PreviousTime=symbol+PreviousTime+BunchLen;
	}
	else PreviousTime=symbol;
	Time=PreviousTime;
      }
      if(NextTableType>1)
	BufferFile.FillBuffer(symbol);
      NextTable(symbol,NextTableType,BunchLen,Count); 
      if(NextTableType==0){
	BufferFile.FillBuffer(Time);
	BufferFile.FillBuffer(BunchLen+2);
	BunchLen=0;
      }
    }//end for
    BufferFile.WriteTrailer(NumWords,PadNumber,RowNumber,SecNumber);
  }//end while
  cout<<"Number of decoded words:"<<WordsRead<<endl;
  f.close();
  //The trees are deleted 
  for(Int_t j=0;j<NumTables;j++){
      DeleteHuffmanTree(RootNode[j]);
  }//end for
  delete [] RootNode; 
  return 0; 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t AliTPCCompression::DecompressDataOptTables(Int_t NumTables,const char* fname,char* fDest){
  if(fVerbose){
    cout<<"   DECOMPRESSION:"<<endl;
    cout<<"Source File "<<fname<<" Destination File "<<fDest<<endl; 
  }
  AliTPCHNode ** RootNode = new AliTPCHNode*[NumTables];
  //Creation of the Huffman trees
  CreateTreesFromFile(RootNode,NumTables);
  f.open(fname,ios::binary|ios::in);
  if(!f){cout<<"File doesn't exist\n";return -1;}
  //to go to the end of the file
  f.seekg(0,ios::end);
  //to get the file dimension in byte
  fPos=f.tellg();
  fPos-=sizeof(ULong_t);
  f.seekg(fPos);
  fReadBits=0;
  fBuffer=0;
  f.read((char*)(&fBuffer),sizeof(ULong_t));
  Int_t bit=0;
  ULong_t Mask=0x1;
  while(!bit){
    bit=fBuffer&Mask;
    Mask=Mask<<1;
    fReadBits++;
  }
  ULong_t PacketNumber=ReadWord(sizeof(ULong_t)*8);
  if(fVerbose){
    cout<<"Number of Packect: "<<PacketNumber<<endl;
  }
  AliTPCBuffer160 BufferFile(fDest,1);
  ULong_t k=0;
  ULong_t WordsRead=0; //number of read coded words 
  while(k<PacketNumber){
    Int_t NumWords,PadNumber,RowNumber,SecNumber=0;
    ReadTrailer(NumWords,PadNumber,RowNumber,SecNumber);
    k++;
    WordsRead+=4;
    Int_t PreviousTime=-1;
    Int_t Time=0;
    Int_t NextTableType=0;
    Int_t BunchLen=0;
    Int_t Count=0;
    for(Int_t i=0;i<NumWords;i++){
      ULong_t symbol=GetDecodedWord(RootNode[NextTableType]);
      WordsRead++;
      //Time reconstruction
      if (NextTableType==1){
	if (PreviousTime!=-1){
	  PreviousTime=symbol+PreviousTime+BunchLen;
	}
	else PreviousTime=symbol;
	Time=PreviousTime;
      }
      if(NextTableType>1)
	BufferFile.FillBuffer(symbol);
      NextTable(symbol,NextTableType,BunchLen,Count); 
      if(NextTableType==0){
	BufferFile.FillBuffer(Time);
	BufferFile.FillBuffer(BunchLen+2);
	BunchLen=0;
      }
    }//end for
    BufferFile.WriteTrailer(NumWords,PadNumber,RowNumber,SecNumber);
  }//end while
  if(fVerbose){
    cout<<"Number of decoded words:"<<WordsRead<<endl;
  }
  f.close();
  //The trees are deleted 
  for(Int_t j=0;j<NumTables;j++){
      DeleteHuffmanTree(RootNode[j]);
  }//end for
  delete [] RootNode;
  return 0; 
}

///////////////////////////////////////////////////////////////////////////////////////////

void AliTPCCompression::ReadAltroFormat(char* fileOut,char* fileIn){
  ofstream ftxt(fileOut);
  AliTPCBuffer160 Buff(fileIn,0);
  Int_t NumWords,PadNum,RowNum,SecNum=0;
  Int_t Value=0;
  while(Buff.ReadTrailerBackward(NumWords,PadNum,RowNum,SecNum) !=-1 ){
    ftxt<<"W:"<<NumWords<<" P:"<<PadNum<<" R:"<<RowNum<<" S:"<<SecNum<<endl;
    if (NumWords%4){
      for(Int_t j=0;j<(4-NumWords%4);j++){
	Value=Buff.GetNextBackWord();
      }//end for
    }//end if
    for(Int_t i=0;i<NumWords;i++){
      Value=Buff.GetNextBackWord();
      ftxt<<Value<<endl;
    }//end for
  }//end while
  ftxt.close();
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////
