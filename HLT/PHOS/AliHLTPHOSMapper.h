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

//#include "AliHLTPHOSBase.h"

//using namespace PhosHLTConst;

//class AliHLTPHOSMapper : public AliHLTPHOSBase
class AliHLTPHOSMapper 
{
 public:
  AliHLTPHOSMapper();
  virtual ~AliHLTPHOSMapper();
  void InitAltroMapping(); 
  bool GetIsInitializedMapping();
  char* GetFilePath();

  struct fAltromap{ 
    int fZRow; // Coordinate in Z direction (beam direction) relatve too one RCU
    int fXCol; // Coordinate in X direction (perpendicular too beam direction an parallell to ground) relatve too one RCU
    int fGain; // Gain (high gain = 1, low gain = 0)
  };

  fAltromap *fHw2geomapPtr; //pointer to structure holding information about geometrical address 
  
  char fFilepath[1024];

 private:
  bool fIsInitializedMapping;
  AliHLTPHOSMapper(const AliHLTPHOSMapper & );
  AliHLTPHOSMapper & operator = (const AliHLTPHOSMapper &);

};

#endif
