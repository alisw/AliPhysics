#ifndef ALIHLTPHOSMIPCLUSTERMANAGER_H
#define ALIHLTPHOSMIPCLUSTERMANAGER_H

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
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
#include "AliHLTPHOSConstants.h"
#include "Rtypes.h"

using namespace PhosHLTConst;

#define maxMipsCandidates  5000
#define Z_CENTER_RANGE 2   //center +/-2 (=5)
#define X_CENTER_RANGE 2   //center +/-2 (=5)


class AliHLTPHOSMipCluster;
class AliHLTPHOSDigit;

class  AliHLTPHOSMipClusterManager
{
public:
  AliHLTPHOSMipClusterManager();
  virtual ~AliHLTPHOSMipClusterManager();

  void AddDigit(AliHLTPHOSDigit *digit);
  void PrintDigitInfo(AliHLTPHOSDigit *digit);
  void PrintClusters(int minEntries = 0);
  void ResetMipManager(); 
  AliHLTPHOSMipCluster* GetCluster(int z, int x);

  
private:
  void ResetClusterFarm();
  void AddToMipTable(AliHLTPHOSDigit *digit);
  bool IsMipCandidate(AliHLTPHOSDigit *digit);
  AliHLTPHOSMipCluster* MoveCluster(int zFrom, int xFrom, int zTo, int xTo );
  void NewMipTableEntry(AliHLTPHOSDigit *digit, int z, int x); 
  AliHLTPHOSMipCluster* GetNeighbour(int zRange, int xRange, int zIn, int xIn);
  //  AliHLTPHOSMipCluster  *fMipClusterFarm[maxMipsCandidates];
  AliHLTPHOSMipCluster *fMipClusterTable[N_ZROWS_MOD][N_XCOLUMNS_MOD];
  AliHLTPHOSMipCluster *fMipClusterFarm[N_ZROWS_MOD][N_XCOLUMNS_MOD];
  int fMipCnt;

  bool fIsMoreClusters;
  
  

};

#endif

