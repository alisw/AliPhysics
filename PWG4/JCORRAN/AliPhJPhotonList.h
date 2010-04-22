// $Id: AliPhJPhotonList.h,v 1.4 2008/05/08 13:44:45 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file AliPhJPhotonList.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.4 $
  \date $Date: 2008/05/08 13:44:45 $
*/
////////////////////////////////////////////////////

#ifndef ALIPHJPHOTONLIST_H
#define ALIPHJPHOTONLIST_H


#include "TClonesArray.h"
#include <iostream>
#include <stdlib.h>

#include "JConst.h"
#include "AliPhJPhoton.h"
#include "AliJPhoton.h"

//class TClonesArray;
//class AliJPhoton;
//class AliPhJPhoton;

class AliPhJPhotonList : public TObject {

public:
  AliPhJPhotonList();
  AliPhJPhotonList(expName exp);
  AliPhJPhotonList(const AliPhJPhotonList& a);
  virtual ~AliPhJPhotonList();

  void Reset();

  //getters
  unsigned short GetNPhotons() const { return fPhotons; }
  AliPhJPhoton*  GetPhoton(const unsigned int iph); 
  AliJPhoton*    GetAliJPhoton(const unsigned int iph);    // ALICE getter
  //setters
  void SetNPhotons(const unsigned short nph) { fPhotons = nph; }
  int  SetTClonesArraySize(const unsigned int nph);
  // add Photon
  void AddAliJPhoton(const unsigned int iph);     // ALICE add

  AliPhJPhotonList& operator=(const AliPhJPhotonList& list);

protected:
  TClonesArray *GetList() const { return fPhotonList; }
  TClonesArray *fPhotonList; //photon  list
  unsigned short fPhotons;   //number of photons in the list

private:
  ClassDef(AliPhJPhotonList,1);

};

#endif
