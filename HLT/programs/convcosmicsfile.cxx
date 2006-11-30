// $Id$

// Author: Constantin Loizides <loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

/**
 Program to convert big <-> little endian cosmics data of 02/2003.
*/

#include "AliHLTStandardIncludes.h"
#include "AliHLTRootTypes.h"

#if __GNUC__ == 3
using namespace std;
#endif

//#define DEBUG
//#define NOCONV

//this flag is for the cosmics/pubsub test
//#define ADDFILESIZE

#ifdef ADDFILESIZE
#include <sys/stat.h>
struct stat stat_results;
#endif

Int_t Convert4(Int_t i)
{ //BigEndian i0i1i2i3 -> LittleEndian i3i2i1i0
#ifdef NOCONV
  return i;
#else
  Char_t *p=(Char_t*)&i;

  Char_t temp[4];
  temp[0]=p[3];
  temp[1]=p[2];
  temp[2]=p[1];
  temp[3]=p[0];

  return (*(Int_t*)temp);
#endif
}

Short_t Convert2(Short_t s)
{ //BigEndian i0i1 -> LittleEndian i1i0
#ifdef NOCONV
  return s;
#else
  Char_t *p=(Char_t*)&s;

  Char_t temp[2];
  temp[0]=p[1];
  temp[1]=p[0];

  return (*(Short_t*)temp);
#endif
}

int main(Int_t argc, Char_t **argv)
{
  //p1 -> filename
  //p2 -> file is in little endian
  Int_t islittle=0;

  if((sizeof(Int_t) != 4) || (sizeof(Short_t) != 2)) {
    cerr << "Check architecture to run the conversion on! Int_t should be 32 and Short_t should be 16 bit." << endl;
    exit(1);
  }

  Char_t fname[1024];
  if(argc==1)
    sprintf(fname,"%s","/home/loizides/tmp/cosmics/run0013/evt0000");
  else if(argc==2)
    strcpy(fname,argv[1]);
  else if(argc>2){
    strcpy(fname,argv[1]);
    islittle=1;
  }

  ifstream *in = new ifstream();
  in->open(fname,fstream::binary);
  if(!in->is_open()){
    cerr << "Error opening input file " << fname << endl;
    exit(1);
  }

  Char_t sname[1024];
#ifdef ADDFILESIZE
  sprintf(sname,"%s.added",fname);
#else
  sprintf(sname,"%s.conv",fname);
#endif
  ofstream *out = new ofstream();
  out->open(sname,fstream::in | fstream::out | fstream::binary | fstream::trunc);
  if(!in->is_open()){
    cerr << "Error opening output file " << sname << endl;
    exit(1);
  }


  Int_t dummy4;
  Short_t dummy2;
#ifdef ADDFILESIZE
  if (stat(fname, &stat_results) == 0){
    dummy4=stat_results.st_size/4;
    cout << "Add file size: " << dummy4 << endl;
  } else {
    cerr << "Error stating input file " << fname << endl;
    exit(1);
  }
  out->write((Char_t*)&dummy4,sizeof(dummy4));
#endif

  in->read((Char_t*)&dummy4,sizeof(dummy4));
  const Int_t knumofChannels = Convert4(dummy4);
#ifdef DEBUG
  cout << knumofChannels << endl;
#endif
  out->write((Char_t*)&knumofChannels,sizeof(knumofChannels));

  Int_t channelloop=knumofChannels;
  if(islittle) channelloop=dummy4;
  for(Int_t i = 0 ; i < channelloop ; i++){
    in->read((Char_t*)&dummy2,sizeof(dummy2));
    Short_t channel = Convert2(dummy2);
    out->write((Char_t*)&channel,sizeof(channel));
#ifdef DEBUG
    cout << channel << " ";
#endif
  }
#ifdef DEBUG
  cout << endl;
#endif

  in->read((Char_t*)&dummy4,sizeof(dummy4));
  const Int_t numofChannelsTest = Convert4(dummy4);
#ifdef DEBUG
  cout << numofChannelsTest << endl;
#endif
  out->write((Char_t*)&numofChannelsTest,sizeof(numofChannelsTest));

  if (knumofChannels != numofChannelsTest){
    cerr << "Error in File format: \"channels " << knumofChannels << " must be channels-test " << numofChannelsTest << "\" "<< endl;
  }

  in->read((Char_t*)&dummy4,sizeof(dummy4));
  const Int_t knumofTimebins = Convert4(dummy4);
#ifdef DEBUG
  cout << knumofTimebins << endl;
#endif
  out->write((Char_t*)&knumofTimebins,sizeof(knumofTimebins));

  Int_t timeloop=knumofTimebins;
  if(islittle) timeloop=dummy4;

  cout << "Input:  " << fname << endl;
  cout << "Output: " << sname << endl;
  cout << "Channels: " << knumofChannels << endl;
  cout << "Timebins: " << knumofTimebins << endl;

  for(int channel = 0; channel < channelloop; channel++){
    for(int timebin = 0 ; timebin < timeloop ; timebin++){
      in->read((Char_t*)&dummy2,sizeof(dummy2));
      Short_t charge = Convert2(dummy2);
#ifdef DEBUG
      cout << charge << " ";
#endif
      out->write((Char_t*)&charge,sizeof(charge));
    }
#ifdef DEBUG
    cout << endl;
#endif
  }

  in->close();
  out->close();
}



/*
  Data format as of 28/01/2003 R. Bramm explained,
  sorry for the German text, too lazy to translate now

  Die Files beinhalten 512 Kanaele (Altrokanäle also 32 Altros)
  mit je 999 Timebins. -> depends on detailed setting of the test

  Format ist
  Anzahl der Kaenaele als int32
  Liste der Kanaele je int16
  Anzahl der Kaenaele als int32
  Anzahl der Timebins

  und dann vollkommenformatlos
  Anzahl der Kaenaele * Anzahl der Timebins (512*1000 values)

  Bezogen auf das Mapping:
----------------------------
  Die Kanalnummer, die die nehmen, ist Spalte 0 oder 9 also der 
  Altrochannel. Die GSI nimmt aber einfach die "FEC Channels" 
  was einfach die Anschlusse des Kabels durchgezaehlt sind ... 
  dummerweise ist FECChannel != Altrochannel

  FEC ist Front End Card
  FEC Channel ist die Kanalnummer ... 
  Altro channel ist durchnummeriert vom 0 bis 511 also fuer 32 Altros. 
  (insofern ist der Name falsch ... weils mehrere altros sind.)

  es geht aber mit 
  8,9,10,11,12,13,14,15,7,6,5,4,3,2,1,0 (alles + 16)
  weiter weil der Altro auf der Unterseite liegt

  pad und row kann man umrechnen (steht nicht die formel im .doc ?) in x,y 
  timebin ist "komplizierter" denn da brauchtm an driftgeschw. + samplingfreq
  ja, aber in meiner tabelle mit den Altrochannels steht alles mit drin.

  es existiert noch
  http://137.138.45.89:8080/svn/repos/TPCMapping/trunk/MappingData/mappingRowPad.data
  dadrin steht die anzahl der pads pro row, rownummer, liste der pads
  die pads sind in "altrochannel" codiert ... war einfacher das so einzuwurschteln 
  ... brauche ich aber nur fuer die Row view und ich hatte keinen Bock, das 
  jedes mal aus dem Ulifile zu erzeugen ...

  wenn du in http://mactpc01.cern.ch:8080
  /svn/repos/TPCMapping/trunk/MappingData/mapping32313029.data
  schaust, dann ist das erste die "altrochannelnummer" und in spalte 7 
  (c array zaehlweise) 
  steht der FEC channel. in den Events steht also immer 999 mal n adc value 
  fuer jeden "altrochannel" in genau der reihenfolge wie in dem obigen file.

  Im OM.C loope ich von 
  Zeile 124 bis 175 uebers file 
  interessant ist nur das einlesen
  von 124 bis 127
  dann berechne ich die baseline
  128 - 132
  dann kommt maxadc
  134 - 143
  dann rudimentaeres Clusterfindeing abe das ist unwichtig
  143 - 163
  dann kommt das fuellen des TPC Topviews, wo auch das "mapping " drinsteckt
  171-174

  du solltest lesen 
  512 als int32
  0 als int16
  1 als int16
  2 als int16
  3 als int16
  ...
  511 als int16
  512 als int32
  999 als int32
  dann int16 ...


  So zur frage wer ordert:
  Luciano und ich haben gebraeinstormed und kamen zum ergebnis, dass es 
  auch kein Problem waere die RCU pad by pad auslesen zu lassen ... 
  also "komplett korrekt" gemapped.
  man braucht nur nen lookuptable. Da ist keine Latenz, die es problematisch 
  macht, das man 
  fee 1 channel 1 = pad 1 
  fee 1 channel 2 = pad 2
  fee 1 channel 4 = pad 3
  fee 1 channel 3 = pad 4
  fee 2 channel 1 = pad 5
  ...
*/
