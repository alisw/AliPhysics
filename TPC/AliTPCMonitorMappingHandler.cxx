/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/*
  $Log$
  Revision 1.2  2007/10/12 13:36:27  cvetan
  Coding convention fixes from Stefan

  Revision 1.1  2007/09/17 10:23:31  cvetan
  New TPC monitoring package from Stefan Kniege. The monitoring package can be started by running TPCMonitor.C macro located in macros folder.

*/ 

////////////////////////////////////////////////////////////////////////
////
//// AliTPCMonitorMappingHandler class
////
//// Class for handling mapping information TPC  
////  
//// The mapping information for the TPC front end electornics (pads, front end cards) 
//// are handled by this class.
//// The information from the design mapping and from the hardware address can be 
//// cross checked in the TPCMonitor.C. 
//// Malfunctioning front end cards can be identified by looking at single channels 
//// displayed with  the TPCMonitor.  
////   
//// 
//// Authors: Roland Bramm, 
////          Stefan Kniege, IKF, Frankfurt
////       
/////////////////////////////////////////////////////////////////////////


#include <cstdlib>
#include "AliTPCMonitorMappingHandler.h"
#include "TH1.h"
#include "TLegend.h"
#include "AliLog.h"
#include <Riostream.h>

ClassImp(AliTPCMonitorMappingHandler)

//_____________________________________________________________________________________________
AliTPCMonitorMappingHandler::AliTPCMonitorMappingHandler(const Char_t* name, const Char_t* title): 
  TNamed(name,title),
  fnumofChannels(0),
  fmaxHWAdress(0),
  fsizeofArray(0),
  fmapping(new Short_t*[24000]),
  fmappingChannelinRow(new Int_t*[160]),
  fu2ftestmapping(new Short_t*[7000]),
  fMapHwFECglobal(new Int_t*[24000]),
  fecGainMap(new Float_t*[7000])
{
  // Constructor : Initialize mapping arrays
 
  for(Int_t in = 0; in<160; in++)
    {
      Int_t* hold = new Int_t[150]; 
      for(Int_t jn = 0; jn<150;jn++)  hold[jn]=0;
      fmappingChannelinRow[in]= hold;
    }
  
  for(Int_t i = 0; i<7000; i++)
    {
      Short_t* hold = new Short_t[8]; 
      for(Int_t j = 0; j<8;j++)  hold[j]=0;
      fu2ftestmapping[i]= hold;
    }
  
  for(Int_t i = 0; i<24000; i++)
    {
      Int_t* hold = new Int_t[2]; 
      for(Int_t j = 0; j<2;j++)  hold[j]=0;
      fMapHwFECglobal[i]= hold;
    }
  
  
  for(Int_t i = 0; i<7000; i++)
    {
      Float_t* hold = new Float_t[128]; 
      for(Int_t j = 0; j<128;j++)  hold[j]=0;
      fecGainMap[i]= hold;
    }
}

//____________________________________________________________________________
AliTPCMonitorMappingHandler::AliTPCMonitorMappingHandler(const  AliTPCMonitorMappingHandler &maphand) :
  TNamed(maphand.GetName(),maphand.GetTitle()),
  fnumofChannels(maphand.fnumofChannels),
  fmaxHWAdress(maphand.fmaxHWAdress),
  fsizeofArray(maphand.fsizeofArray),
  fmapping(new Short_t*[24000]),
  fmappingChannelinRow(new Int_t*[160]),
  fu2ftestmapping(new Short_t*[7000]),
  fMapHwFECglobal(new Int_t*[24000]),
  fecGainMap(new Float_t*[7000])
{
  // copy constructor
 
 
  for(Int_t in = 0; in<160; in++)
    {
      Int_t* hold = new Int_t[150]; 
      for(Int_t jn = 0; jn<150;jn++)  hold[jn]=maphand.fmappingChannelinRow[in][jn];
      fmappingChannelinRow[in]= hold;
    }
  
  for(Int_t i = 0; i<7000; i++)
    {
      Short_t* hold = new Short_t[8]; 
      for(Int_t j = 0; j<8;j++)  hold[j]=maphand.fu2ftestmapping[i][j];
      fu2ftestmapping[i]= hold;
    }

  for(Int_t i = 0; i<24000; i++)
    {
      Int_t* hold = new Int_t[2]; 
      for(Int_t j = 0; j<2;j++)  hold[j]=maphand.fMapHwFECglobal[i][j];
      fMapHwFECglobal[i]= hold;
    }
  
  for(Int_t i = 0; i<7000; i++)
    {
      Float_t* hold = new Float_t[128]; 
      for(Int_t j = 0; j<128;j++)  hold[j]=maphand.fecGainMap[i][j];
      fecGainMap[i]= hold;
    }

}


//____________________________________________________________________________
AliTPCMonitorMappingHandler &AliTPCMonitorMappingHandler:: operator= (const AliTPCMonitorMappingHandler& maphand)
{
  // assignment operator
  if(this!=&maphand)
    {
      fnumofChannels=maphand.fnumofChannels;
      fmaxHWAdress=maphand.fmaxHWAdress;
      fsizeofArray=maphand.fsizeofArray;

      fmapping = new Short_t*[24000]; // empty
  
      
      fmappingChannelinRow = new Int_t*[160];
      for(Int_t in = 0; in<160; in++)
	{
	  Int_t* hold = new Int_t[150]; 
	  for(Int_t jn = 0; jn<150;jn++)  hold[jn]=maphand.fmappingChannelinRow[in][jn];
	  fmappingChannelinRow[in]= hold;
	}
      
      fu2ftestmapping      = new Short_t*[7000];
      for(Int_t i = 0; i<7000; i++)
	{
	  Short_t* hold = new Short_t[8]; 
	  for(Int_t j = 0; j<8;j++)  hold[j]=maphand.fu2ftestmapping[i][j];
	  fu2ftestmapping[i]= hold;
	}
      
      fMapHwFECglobal = new Int_t*[24000];
      for(Int_t i = 0; i<24000; i++)
	{
	  Int_t* hold = new Int_t[2]; 
	  for(Int_t j = 0; j<2;j++)  hold[j]=maphand.fMapHwFECglobal[i][j];
	  fMapHwFECglobal[i]= hold;
	}
      
      fecGainMap = new Float_t*[7000];
      for(Int_t i = 0; i<7000; i++)
	{
	  Float_t* hold = new Float_t[128]; 
	  for(Int_t j = 0; j<128;j++)  hold[j]=maphand.fecGainMap[i][j];
	  fecGainMap[i]= hold;
	}
      
      
    }
  return *this;
}


//_____________________________________________________________________________________________
AliTPCMonitorMappingHandler::~AliTPCMonitorMappingHandler() 
{
  // Destructor
  Short_t* temp;
  for(Int_t i = 0; i < 24000 ; i++) {
    temp = (Short_t *) fmapping[i];
    if(temp[1] != -1)
      delete temp;
  }
  temp = (Short_t *) fmapping;
  delete temp;
  
  
  Int_t* temp1;
  for(Int_t i = 0; i < 24000 ; i++) {
    temp1 = (Int_t *) fMapHwFECglobal[i];
    delete temp1;
  }
  temp1 = (Int_t *) fMapHwFECglobal;
  delete temp1;
 
  for(Int_t i = 0; i < 7000 ; i++) {
    temp = (Short_t *) fu2ftestmapping[i];
    delete temp;
  }
 
  temp = (Short_t *) fu2ftestmapping;
  delete temp;
 	 
  for(Int_t i = 0; i < 150 ; i++) {
    temp1= (Int_t *) fmappingChannelinRow[i];
    delete temp1;
  }
  temp1 = (Int_t *) fmappingChannelinRow;
  delete temp1;

  Float_t* temp2;
  for(Int_t i = 0; i<7000; i++)
    {
      temp2= (Float_t *) fecGainMap[i];
      delete temp2;
    }
  temp2 = (Float_t *) fecGainMap;
  delete temp2;

	
}

//_____________________________________________________________________________________________
void AliTPCMonitorMappingHandler::ReadMapping(char* mapfile)
{
  // Read global Mapping file
  // Format of data in mapping file:
  // column 0:  hadrware address
  // column 1:  readout index in IROC/OROC
  // column 2:  global pad row number (0-158)
  // column 3   pad number in row
  // column 4:  connector
  // column 5:  pin
  // column 6:  fec number in IROC/OROC
  // column 7:  fec channel
  // column 8:  fec connector
  // column 9:  altro channel
  // column 10: altro chip
 
  // Mapping data for a given hardware address are placed at the 
  // index corresponding to the value of the hardware address in 
  // the fmapping array..
  // The mapping information for the hardware address 'hwaddr'
  // can hence be found in fmapping[hwaddr]


  Short_t*  mappingRow; 
  char readcarry[255];
  Int_t version = -1;
  Int_t actPos  = 0;
  Int_t oldPos  = 0;

  ifstream *in = new ifstream();
  in->open(mapfile);
  *in >> readcarry;
  version = atoi(readcarry);
  
  *in >> readcarry;
  fnumofChannels = atoi(readcarry);
  *in >> readcarry;
  fmaxHWAdress = atoi(readcarry);
  fsizeofArray = fmaxHWAdress;
	
  Short_t *fmappingEmptyRow = new Short_t[11];
  for(Int_t i = 0; i < 11 ; i++) {
    fmappingEmptyRow[i] = 0;
  }
  fmappingEmptyRow[1] = -1;
  for(Int_t i = 0; i < fnumofChannels ; i++) {
    mappingRow = new Short_t[11];
    for(Int_t j = 0 ; j < 11 ; j++) {
      *in >> readcarry;
      mappingRow[j] = atoi(readcarry);
    }
    actPos = mappingRow[0];
    fmapping[actPos] = mappingRow;
    if( (actPos - oldPos) > 1) {
      for(Int_t j = (oldPos+1); j < actPos; j++) {
        fmapping[j] = fmappingEmptyRow;
      }
    }
    oldPos = actPos;
    fmapping[actPos][0] = i;
  }
  
  in->close();
  delete in;

  delete [] fmappingEmptyRow;
}

//_____________________________________________________________________________________________
Int_t  AliTPCMonitorMappingHandler::ReadFECMapping(char* u2ftestfile)
{
  // Read in Mapping of global FEC numbers to branches, rcu patches etc.
  // Format in Mapping file
  // column 0: global card number (set to 0 if card number does not exist)
  // column 1: side of the TPC
  // column 2: sector on the side   (0-17)
  // column 3: rcu in the sector    (0-5)
  // column 4: fec number in rcu    (0-24)
  // column 5: fec number in branch (0-12)
  // column 6: branch number        (0,1)
  // column 7: local hardware address of first channel (not used)

  // Order of data is kept in fu2ftestmapping
  
  ifstream datin(u2ftestfile);
  
  Int_t carry     = 0;
  Int_t ncards    = 0;
 
  for(Int_t ind = 0; ind<7000; ind++)
    {
      for(Int_t entr = 0; entr<8; entr++)
	{
	  datin >> carry ;
	  fu2ftestmapping[ind][entr] = carry  ;
	  
	  if(entr==0 && carry!=0) ncards++;
	}
    }
  return ncards ;
}

 

//_____________________________________________________________________________________________
void AliTPCMonitorMappingHandler::ReadfecHwMap(Int_t sector)
{
  // Create mapping of global FEC numbers to hardware addresses for a given sector

  Int_t  fside           =0;
  Int_t  fsector        = 0;
  Int_t  fec            = 0;
  Int_t  branch         = 0;
  Int_t  rcupatch      = 0;
  Int_t  altrchann      = 0;
  Int_t  altrchip       = 0;
  Int_t  nextHwAddress  = 0;
  Int_t  nfecs          = 0;

  if(sector/18==0) fside =65;
  else             fside =67;
  
  if(sector>18)    fsector= sector-18;
  else             fsector= sector   ;

  for(Int_t ind = 0; ind<7000; ind++)
    { 
      if((Int_t)U2fGetSide(ind)==fside && U2fGetSector(ind)==fsector) 
	{
	  nfecs++;
	  fec            = U2fGetFECinBranch(ind);
	  branch         = U2fGetBranch(ind);
	  rcupatch      = U2fGetRCU(ind);
	  	  
	  for(Int_t ch = 0; ch<128; ch++)
	    {
	      altrchann      = ch%16;
	      altrchip       = ch/16;
	      
	      nextHwAddress  = (   ((branch&1)<<11) + (fec<<7) + (altrchip<<4) + (altrchann)  + ((rcupatch-1)<<12) ); 
	      
	      fMapHwFECglobal[nextHwAddress][0] = ind;
	      fMapHwFECglobal[nextHwAddress][1] = ch ;
	    }
	}
    }
}
//_____________________________________________________________________________________________
void AliTPCMonitorMappingHandler::ReadfecGainMap(char* fecgainmap)
{
  // Read global gain calibration pap
  // Format in file :
  // colummn 0      : FEC number 
  // colums  1-128  : gain calibration factors 
  ifstream datin(fecgainmap);
  
  Int_t   fecnr  = 0;
  Float_t val    = 0.0 ;
  
  while(!datin.eof())
  {
    datin >> fecnr ;
    for(Int_t in = 0; in<128; in++)
    {
      datin >> val ;
      fecGainMap[fecnr][in] = val;
    }
  }
}

//_____________________________________________________________________________________________
void  AliTPCMonitorMappingHandler::ReadRowMappingGlob(char* fpathtoMappingRowfile) 
{
  // Read mapping of hardware addresses in rows
  // Format of file:
  // column 0:        global row number (0-158)
  // column 1:        number of pads in this row (npads)
  // column 2-npads:  hardware addresses for these pads
 
  char readcarry[256];
  char readcarry2[256];
  ifstream *in = new ifstream();
  in->open(fpathtoMappingRowfile);
  
  for(Int_t i = 0; i < 159 ; i++) {
    *in >> readcarry;   // row number 
    *in >> readcarry2;  // numof pads
    fmappingChannelinRow[i][0] = atoi(readcarry2);
    fmappingChannelinRow[i][1] = atoi(readcarry);
    
    for(Int_t j = 2 ; j < fmappingChannelinRow[i][0]+2 ; j++) {
      *in >> readcarry;
      fmappingChannelinRow[i][j] = atoi(readcarry);
    }
  }
  in->close();
}



//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::GetNumOfChannels() const
{
  // Return number of channels
  return fnumofChannels;
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::GetSizeofArray() const
{
  // Return sise of global mapping fmapping array.
  // Value orresponds to max value of hardware addresses
  return fsizeofArray;
}


//_____________________________________________________________________________________________
Short_t* AliTPCMonitorMappingHandler::GetLine(Int_t hwaddr)const
{
  // Return pointer to mapping array for the hardware address hwaddr 
  Short_t* retval;
  if(hwaddr <= fsizeofArray)
    retval = fmapping[hwaddr];
  else
    retval = 0;
  return retval;
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::GetIndex(Int_t hwaddr) const
{
  // Return readout index for the hardware address hwaddr
  Int_t retval;
  if(hwaddr <= fsizeofArray)
    retval = fmapping[hwaddr][1];
  else
    retval = 0;
  return retval;
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::GetPadRow(Int_t hwaddr) const
{
  // Return global pad row  for the hardware address hwaddr 
  Int_t retval;
  if(hwaddr <= fsizeofArray)
    retval = fmapping[hwaddr][2];
  else
    retval = 0;
  return retval;
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::GetPad(Int_t hwaddr) const
{
  // Return pad number in row for the hardware address hwaddr 
  Int_t retval;
  if(hwaddr < fsizeofArray)
    retval = fmapping[hwaddr][3];
  else
    retval = 0;
  return retval;
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::GetConnector(Int_t hwaddr) const
{
  // Return connector for the hardware address hwaddr 
  Int_t retval;
  if(hwaddr <= fsizeofArray)
    retval = fmapping[hwaddr][4];
  else
    retval = 0;
  return retval;
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::GetPin(Int_t hwaddr) const
{
  // Return pin for the hardware address hwaddr 
  Int_t retval;
  if(hwaddr <= fsizeofArray)
    retval = fmapping[hwaddr][5];
  else
    retval = 0;
  return retval;
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::GetFEC(Int_t hwaddr) const
{
  // Return fec number in IROC/OROC  for the hardware address hwaddr 
  Int_t retval;
  if(hwaddr <= fsizeofArray)
    retval = fmapping[hwaddr][6];
  else
    retval = 0;
  return retval;
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::GetFECchannel(Int_t hwaddr) const
{
  // Return FEC channel for the hardware address hwaddr 
  Int_t retval;
  if(hwaddr < fsizeofArray)
    retval = fmapping[hwaddr][7];
  else
    retval = 0;
  return retval;
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::GetFECconnector(Int_t hwaddr) const
{
  // Return FEC connector for the hardware address hwaddr 
  Int_t retval;
  if(hwaddr <= fsizeofArray)
    retval = fmapping[hwaddr][8];
  else
    retval = 0;
  return retval;
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::GetAltroChannel(Int_t hwaddr)  const
{
  // Return Altro channel for the hardware address hwaddr 
  Int_t retval;
  if(hwaddr <= fsizeofArray)
    retval = fmapping[hwaddr][9];
  else
    retval = 0;
  return retval;
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::GetAltro(Int_t hwaddr)  const
{
  // Return Altro chip number in FEC for the hardware address hwaddr 
  Int_t retval;
  if(hwaddr <= fsizeofArray)
    retval = fmapping[hwaddr][10];
  else
    retval = 0;
  return retval;
}


//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::GetNumofPads(Int_t row) 
{
  // Return number of pads in row
  if(row<159)
    return fmappingChannelinRow[row][0];
  else  {
    AliError("Row number to high");
    return 0;
  }
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::GetPadAddInRow(Int_t row,Int_t pad )
{
  // Return hardware address for given pad in row
  if(row<159)
    return fmappingChannelinRow[row][pad+2];
  else {
    AliError("Row number to high");
    return 0;
  }
}



//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::U2fGetFECnr(Int_t index) const
{
  // Return FEC number for index  (FEC number should be equal to index)
  return  fu2ftestmapping[index][0];
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::U2fGetSide(Int_t fecnr) const
{
  // Return side on which FEC is installed
  return  fu2ftestmapping[fecnr][1];
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::U2fGetSector(Int_t fecnr) const
{
  // Return sector in which FEC is installed
  return  fu2ftestmapping[fecnr][2];
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::U2fGetRCU(Int_t fecnr) const
{
  // Rerurn rcu patch in which FEC is installed
  return  fu2ftestmapping[fecnr][3];
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::U2fGetFECinRCU(Int_t fecnr) const
{
  // Return index of FEC in RCU (0-25)
  return  fu2ftestmapping[fecnr][4];
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::U2fGetFECinBranch(Int_t fecnr) const
{
  // Return index of FEC in branch (0-12)
  return  fu2ftestmapping[fecnr][5];
}

//_____________________________________________________________________________________________
Int_t AliTPCMonitorMappingHandler::U2fGetBranch(Int_t fecnr) const
{
  // Return branch in which FEC is installed (0,1)
  return  fu2ftestmapping[fecnr][6];
    
}
