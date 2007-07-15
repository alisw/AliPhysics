#ifndef ALIHLTPHOSMAPPER_H
#define ALIHLTPHOSMAPPER_H

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2006                                       *
 *                                                                        * 
 * Author: Per Thomas Hille perthi@fys.uio.no for the ALICE DCS Project.  *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                * 
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//#include "PhosFeeClient.h"

//#include "stdio.h"
//#include <iostream>
#include <cstdlib>
#include <assert.h>
//#include "AliHLTPHOSCommonDefs.h"
//#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSBase.h"
    
//            PhosHLTConst
using namespace std;
using namespace PhosHLTConst;

class AliHLTPHOSMapper : public AliHLTPHOSBase
{
 public:
  AliHLTPHOSMapper();
  virtual ~AliHLTPHOSMapper();
  void InitAltroMapping(); 
  

  struct altromap{ 
    //    int mod;
    int row;
    int col;
    int gain;
    //    int rcu;
    //    int branch;
    //    int card;
    //    int chip;
    //    int chan;
    //    int csp;
    //    int num;
    //    int hid;
  };


altromap *hw2geomapPtr;


  /*
  struct FEE_CSP_MAP{ 
    int row;
    int col;
    int gain;
    int csp;
    int num;
  };
  */

  /*
  struct ALTRO_GEO_MAP{ 
    int mod;
    int gain;
    int row;
    int col;
    int rcu;
    int branch;
    int card;
    int chip;
    int chan;
    int csp;
    int num;
    int hid;
  };
  */

  //  void AddCsp(int csp, int chip, int chHi, int chLo, int numHi, int numLo);

 
  //Function to generate Active Channel List (ACL)for each of the four readout partitions
  //Of the Phos Module. The ACL register is 256x16 bit big.
  //  void GenerateACL(int startZ, int endZ, int startX, int endX, int mID, int acl[RCUS_PER_MODULE][256], unsigned long int afl[RCUS_PER_MODULE]);

  //  void InitAltroCspMapping();

  /*
  inline int Geo2hid(int mod, int gain, int row, int col);
  inline int Hid2mod(int hid);
  inline int Hid2gain(int hid);
  inline int Hid2row(int hid);
  inline int Hid2col(int hid);
  inline int ExtractHid(char *objName);  

  inline void InitAltroMapping(int saveMapping);
  inline void PrintHistMapInfo(char *objName);
  */

  /*
  int Geo2hid(int mod, int gain, int row, int col);
  int Hid2mod(int hid);
  int Hid2gain(int hid);
  int Hid2row(int hid);
  int Hid2col(int hid);
  int ExtractHid(char *objName);  

  void InitAltroMapping(int saveMapping);
  void PrintHistMapInfo(char *objName);

  FEE_CSP_MAP CSP_MAP[N_ALTROS][N_ALTROCHANNELS];
  ALTRO_GEO_MAP ALTRO_MAP[N_MODULES*N_RCUS*N_BRANCHES*N_FEECS*N_ALTROS*N_ALTROCHANNELS];  
  int hdw2geo[N_MODULES][N_RCUS][N_BRANCHES][N_FEECS][N_ALTROS][N_ALTROCHANNELS];
  int geo2hdw[N_MODULES][N_GAINS][ N_XCOLUMNS_MOD][N_ZROWS_MOD]; 
  */

 private:



};

#endif
