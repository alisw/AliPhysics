#ifndef ALIHLTEMCALRAWANALYZERCOMPONENT_H
#define ALIHLTEMCALRAWANALYZERCOMPONENT_H

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

#include "AliHLTCaloRawAnalyzerComponentv3.h"
 
//class AliHLTCaloMapper;

class  AliHLTEMCALRawAnalyzerComponent : public AliHLTCaloRawAnalyzerComponentv3
{
 public:
  AliHLTEMCALRawAnalyzerComponent();
  virtual ~AliHLTEMCALRawAnalyzerComponent();
  //  virtual int DoInit(int argc =0, const char** argv  = 0);
  // virtual int Deinit() {} ;
  //  virtual const char* GetComponentID() = 0;
  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>& list);

  virtual const char* GetComponentID() = 0;
  /** interface function, see @ref AliHLTComponent for description */
  
  //  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>& list);
  //  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>& list);
 
 /** interface function, see @ref AliHLTComponent for description */
  //  virtual AliHLTComponentDataType GetOutputDataType();

  /** interface function, see @ref AliHLTComponent for description */
 
  //  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

  /** interface function, see @ref AliHLTComponent for description */
  virtual AliHLTComponent* Spawn() = 0; 

 protected:
  // virtual bool CheckInputDataType(const AliHLTComponentDataType &datatype) = 0;
  virtual bool CheckInputDataType(const AliHLTComponentDataType &datatype);
  

 private:
  AliHLTEMCALRawAnalyzerComponent(const AliHLTEMCALRawAnalyzerComponent & );

  /** Keep the assignement operator private since it should not be used */
  AliHLTEMCALRawAnalyzerComponent & operator = (const AliHLTEMCALRawAnalyzerComponent  &);

  virtual void InitMapping();

  //  AliHLTCaloMapper *fMapperPtr;

};

#endif
