// $Id: AliPhJHeaderList.h,v 1.4 2008/05/08 13:44:45 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file AliPhJHeaderList.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.4 $
  \date $Date: 2008/05/08 13:44:45 $
*/
////////////////////////////////////////////////////

#ifndef ALIPHJHEADERLIST_H
#define ALIPHJHEADERLIST_H

#include "TClonesArray.h"
#include <iostream>

#include "AliPhJBaseHeader.h"
#include "AliJHeader.h"
#include "JConst.h"

//class AliJHeader;
//class AliPhJBaseHeader;
//class TClonesArray;


class AliPhJHeaderList : public TObject {

public:
  AliPhJHeaderList(); // default constructor
  AliPhJHeaderList(expName exp);
  AliPhJHeaderList(const AliPhJHeaderList& a);
  virtual ~AliPhJHeaderList();

  void Reset();

  //getters
  unsigned short GetNHeaders() const { return fHeaders; }
  AliPhJBaseHeader*       GetHeader(const unsigned int ihdr); 
  AliJHeader*    GetAliJHeader(const unsigned int ihdr);    // ALICE getter
  //setters
  void SetNHeaders(const unsigned short nhdr) { fHeaders = nhdr; }
  int  SetTClonesArraySize(const unsigned int nhdr);
  // add header
  void AddAliJHeader(const unsigned int ihdr);     // ALICE add

  AliPhJHeaderList& operator=(const AliPhJHeaderList& list);  

protected:
  TClonesArray *GetList() const { return fHeaderList; }
  TClonesArray *fHeaderList;   //header list
  unsigned short fHeaders;//number of headers
 
  

private:
  ClassDef(AliPhJHeaderList,1)

};

#endif
