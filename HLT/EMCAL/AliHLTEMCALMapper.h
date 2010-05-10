#ifndef ALIHLTEMCALMAPPER_H
#define ALIHLTEMCALMAPPER_H

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


#include  "AliHLTCaloMapper.h"


class  AliHLTEMCALMapper : public AliHLTCaloMapper
{
 public:
  AliHLTEMCALMapper(  const unsigned long  specifiaction );
  virtual ~AliHLTEMCALMapper();
  virtual Bool_t InitAltroMapping( const unsigned long specification ); 
  virtual void InitDDLSpecificationMapping();
  virtual void GetLocalCoord(const int channelId, Float_t* localCoord) const; 
  
 private:
  AliHLTEMCALMapper();
  AliHLTEMCALMapper(const AliHLTEMCALMapper & );
  AliHLTEMCALMapper & operator = (const AliHLTEMCALMapper &);
  const char* DDL2RcuMapFileName(const int ddlid) const;
  
};

#endif
