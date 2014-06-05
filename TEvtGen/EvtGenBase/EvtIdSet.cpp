//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtId.cc
//
// Description: Class for particle Id used in EvtGen.
//
// Modification history:
//
//    DJL     Jan 4, 2000        Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include <iostream>
#include <string>

EvtIdSet::EvtIdSet(const EvtId name1) {
  _numInList=1;
  _list=new EvtId[_numInList];

  _list[0]=name1;
}

EvtIdSet::EvtIdSet(const std::string name1){
  _numInList=1;
  _list=new EvtId[_numInList];

  _list[0]=EvtPDL::getId(name1);
}


EvtIdSet::EvtIdSet(const EvtId name1,
		   const EvtId name2){
  _numInList=2;
  _list=new EvtId[_numInList];

  _list[0]=name1;
  _list[1]=name2;
}


EvtIdSet::EvtIdSet(const std::string name1,
		   const std::string name2){
  _numInList=2;
  _list=new EvtId[_numInList];

  _list[0]=EvtPDL::getId(name1);
  _list[1]=EvtPDL::getId(name2);
}


EvtIdSet::EvtIdSet(const EvtId name1,
		   const EvtId name2,
		   const EvtId name3){
  _numInList=3;
  _list=new EvtId[_numInList];

  _list[0]=name1;
  _list[1]=name2;
  _list[2]=name3;

}


EvtIdSet::EvtIdSet(const std::string name1,
		   const std::string name2,
		   const std::string name3){
  _numInList=3;
  _list=new EvtId[_numInList];

  _list[0]=EvtPDL::getId(name1);
  _list[1]=EvtPDL::getId(name2);
  _list[2]=EvtPDL::getId(name3);
}


EvtIdSet::EvtIdSet(const EvtId name1,
		   const EvtId name2,
		   const EvtId name3,
		   const EvtId name4){
  _numInList=4;
  _list=new EvtId[_numInList];

  _list[0]=name1;
  _list[1]=name2;
  _list[2]=name3;
  _list[3]=name4;
}


EvtIdSet::EvtIdSet(const std::string name1,
		   const std::string name2,
		   const std::string name3,
		   const std::string name4){
  _numInList=4;
  _list=new EvtId[_numInList];

  _list[0]=EvtPDL::getId(name1);
  _list[1]=EvtPDL::getId(name2);
  _list[2]=EvtPDL::getId(name3);
  _list[3]=EvtPDL::getId(name4);
}


EvtIdSet::EvtIdSet(const EvtId name1,
		   const EvtId name2,
		   const EvtId name3,
		   const EvtId name4,
		   const EvtId name5){
  _numInList=5;
  _list=new EvtId[_numInList];

  _list[0]=name1;
  _list[1]=name2;
  _list[2]=name3;
  _list[3]=name4;
  _list[4]=name5;
}


EvtIdSet::EvtIdSet(const std::string name1,
		   const std::string name2,
		   const std::string name3,
		   const std::string name4,
		   const std::string name5){
  _numInList=5;
  _list=new EvtId[_numInList];

  _list[0]=EvtPDL::getId(name1);
  _list[1]=EvtPDL::getId(name2);
  _list[2]=EvtPDL::getId(name3);
  _list[3]=EvtPDL::getId(name4);
  _list[4]=EvtPDL::getId(name5);
}


EvtIdSet::EvtIdSet(const EvtId name1,
		   const EvtId name2,
		   const EvtId name3,
		   const EvtId name4,
		   const EvtId name5,
		   const EvtId name6){
  _numInList=6;
  _list=new EvtId[_numInList];

  _list[0]=name1;
  _list[1]=name2;
  _list[2]=name3;
  _list[3]=name4;
  _list[4]=name5;
  _list[5]=name6;
}


EvtIdSet::EvtIdSet(const std::string name1,
		   const std::string name2,
		   const std::string name3,
		   const std::string name4,
		   const std::string name5,
		   const std::string name6){
  _numInList=6;
  _list=new EvtId[_numInList];

  _list[0]=EvtPDL::getId(name1);
  _list[1]=EvtPDL::getId(name2);
  _list[2]=EvtPDL::getId(name3);
  _list[3]=EvtPDL::getId(name4);
  _list[4]=EvtPDL::getId(name5);
  _list[5]=EvtPDL::getId(name6);
}


EvtIdSet::EvtIdSet(const EvtId name1,
		   const EvtId name2,
		   const EvtId name3,
		   const EvtId name4,
		   const EvtId name5,
		   const EvtId name6,
		   const EvtId name7){
  _numInList=7;
  _list=new EvtId[_numInList];

  _list[0]=name1;
  _list[1]=name2;
  _list[2]=name3;
  _list[3]=name4;
  _list[4]=name5;
  _list[5]=name6;
  _list[6]=name7;
}


EvtIdSet::EvtIdSet(const std::string name1,
		   const std::string name2,
		   const std::string name3,
		   const std::string name4,
		   const std::string name5,
		   const std::string name6,
		   const std::string name7){
  _numInList=7;
  _list=new EvtId[_numInList];

  _list[0]=EvtPDL::getId(name1);
  _list[1]=EvtPDL::getId(name2);
  _list[2]=EvtPDL::getId(name3);
  _list[3]=EvtPDL::getId(name4);
  _list[4]=EvtPDL::getId(name5);
  _list[5]=EvtPDL::getId(name6);
  _list[6]=EvtPDL::getId(name7);
}


EvtIdSet::EvtIdSet(const EvtId name1,
		   const EvtId name2,
		   const EvtId name3,
		   const EvtId name4,
		   const EvtId name5,
		   const EvtId name6,
		   const EvtId name7,
		   const EvtId name8){
  _numInList=8;
  _list=new EvtId[_numInList];

  _list[0]=name1;
  _list[1]=name2;
  _list[2]=name3;
  _list[3]=name4;
  _list[4]=name5;
  _list[5]=name6;
  _list[6]=name7;
  _list[7]=name8;
}


EvtIdSet::EvtIdSet(const std::string name1,
		   const std::string name2,
		   const std::string name3,
		   const std::string name4,
		   const std::string name5,
		   const std::string name6,
		   const std::string name7,
		   const std::string name8){
  _numInList=8;
  _list=new EvtId[_numInList];

  _list[0]=EvtPDL::getId(name1);
  _list[1]=EvtPDL::getId(name2);
  _list[2]=EvtPDL::getId(name3);
  _list[3]=EvtPDL::getId(name4);
  _list[4]=EvtPDL::getId(name5);
  _list[5]=EvtPDL::getId(name6);
  _list[6]=EvtPDL::getId(name7);
  _list[7]=EvtPDL::getId(name8);
}


EvtIdSet::EvtIdSet(const EvtId name1,
		   const EvtId name2,
		   const EvtId name3,
		   const EvtId name4,
		   const EvtId name5,
		   const EvtId name6,
		   const EvtId name7,
		   const EvtId name8,
		   const EvtId name9){
  _numInList=9;
  _list=new EvtId[_numInList];

  _list[0]=name1;
  _list[1]=name2;
  _list[2]=name3;
  _list[3]=name4;
  _list[4]=name5;
  _list[5]=name6;
  _list[6]=name7;
  _list[7]=name8;
  _list[8]=name9;
}


EvtIdSet::EvtIdSet(const std::string name1,
		   const std::string name2,
		   const std::string name3,
		   const std::string name4,
		   const std::string name5,
		   const std::string name6,
		   const std::string name7,
		   const std::string name8,
		   const std::string name9){
  _numInList=9;
  _list=new EvtId[_numInList];

  _list[0]=EvtPDL::getId(name1);
  _list[1]=EvtPDL::getId(name2);
  _list[2]=EvtPDL::getId(name3);
  _list[3]=EvtPDL::getId(name4);
  _list[4]=EvtPDL::getId(name5);
  _list[5]=EvtPDL::getId(name6);
  _list[6]=EvtPDL::getId(name7);
  _list[7]=EvtPDL::getId(name8);
  _list[8]=EvtPDL::getId(name9);
}


EvtIdSet::EvtIdSet(const EvtId name1,
		   const EvtId name2,
		   const EvtId name3,
		   const EvtId name4,
		   const EvtId name5,
		   const EvtId name6,
		   const EvtId name7,
		   const EvtId name8,
		   const EvtId name9,
		   const EvtId name10){
  _numInList=10;
  _list=new EvtId[_numInList];

  _list[0]=name1;
  _list[1]=name2;
  _list[2]=name3;
  _list[3]=name4;
  _list[4]=name5;
  _list[5]=name6;
  _list[6]=name7;
  _list[7]=name8;
  _list[8]=name9;
  _list[9]=name10;
}


EvtIdSet::EvtIdSet(const std::string name1,
		   const std::string name2,
		   const std::string name3,
		   const std::string name4,
		   const std::string name5,
		   const std::string name6,
		   const std::string name7,
		   const std::string name8,
		   const std::string name9,
		   const std::string name10){
  _numInList=10;
  _list=new EvtId[_numInList];

  _list[0]=EvtPDL::getId(name1);
  _list[1]=EvtPDL::getId(name2);
  _list[2]=EvtPDL::getId(name3);
  _list[3]=EvtPDL::getId(name4);
  _list[4]=EvtPDL::getId(name5);
  _list[5]=EvtPDL::getId(name6);
  _list[6]=EvtPDL::getId(name7);
  _list[7]=EvtPDL::getId(name8);
  _list[8]=EvtPDL::getId(name9);
  _list[9]=EvtPDL::getId(name10);
}


EvtIdSet::EvtIdSet(const EvtId name1,
		   const EvtId name2,
		   const EvtId name3,
		   const EvtId name4,
		   const EvtId name5,
		   const EvtId name6,
		   const EvtId name7,
		   const EvtId name8,
		   const EvtId name9,
		   const EvtId name10,
		   const EvtId name11){
  _numInList=11;
  _list=new EvtId[_numInList];

  _list[0]=name1;
  _list[1]=name2;
  _list[2]=name3;
  _list[3]=name4;
  _list[4]=name5;
  _list[5]=name6;
  _list[6]=name7;
  _list[7]=name8;
  _list[8]=name9;
  _list[9]=name10;
  _list[10]=name11;
}


EvtIdSet::EvtIdSet(const std::string name1,
		   const std::string name2,
		   const std::string name3,
		   const std::string name4,
		   const std::string name5,
		   const std::string name6,
		   const std::string name7,
		   const std::string name8,
		   const std::string name9,
		   const std::string name10,
		   const std::string name11){
  _numInList=11;
  _list=new EvtId[_numInList];

  _list[0]=EvtPDL::getId(name1);
  _list[1]=EvtPDL::getId(name2);
  _list[2]=EvtPDL::getId(name3);
  _list[3]=EvtPDL::getId(name4);
  _list[4]=EvtPDL::getId(name5);
  _list[5]=EvtPDL::getId(name6);
  _list[6]=EvtPDL::getId(name7);
  _list[7]=EvtPDL::getId(name8);
  _list[8]=EvtPDL::getId(name9);
  _list[9]=EvtPDL::getId(name10);
  _list[10]=EvtPDL::getId(name11);
}


EvtIdSet::EvtIdSet(const EvtId name1,
		   const EvtId name2,
		   const EvtId name3,
		   const EvtId name4,
		   const EvtId name5,
		   const EvtId name6,
		   const EvtId name7,
		   const EvtId name8,
		   const EvtId name9,
		   const EvtId name10,
		   const EvtId name11,
		   const EvtId name12){
  _numInList=12;
  _list=new EvtId[_numInList];

  _list[0]=name1;
  _list[1]=name2;
  _list[2]=name3;
  _list[3]=name4;
  _list[4]=name5;
  _list[5]=name6;
  _list[6]=name7;
  _list[7]=name8;
  _list[8]=name9;
  _list[9]=name10;
  _list[10]=name11;
  _list[11]=name12;
}


EvtIdSet::EvtIdSet(const std::string name1,
		   const std::string name2,
		   const std::string name3,
		   const std::string name4,
		   const std::string name5,
		   const std::string name6,
		   const std::string name7,
		   const std::string name8,
		   const std::string name9,
		   const std::string name10,
		   const std::string name11,
		   const std::string name12){
  _numInList=12;
  _list=new EvtId[_numInList];

  _list[0]=EvtPDL::getId(name1);
  _list[1]=EvtPDL::getId(name2);
  _list[2]=EvtPDL::getId(name3);
  _list[3]=EvtPDL::getId(name4);
  _list[4]=EvtPDL::getId(name5);
  _list[5]=EvtPDL::getId(name6);
  _list[6]=EvtPDL::getId(name7);
  _list[7]=EvtPDL::getId(name8);
  _list[8]=EvtPDL::getId(name9);
  _list[9]=EvtPDL::getId(name10);
  _list[10]=EvtPDL::getId(name11);
  _list[11]=EvtPDL::getId(name12);
}


EvtIdSet::EvtIdSet(const EvtIdSet& set1){

  _numInList=set1.sizeOfSet();
  _list=new EvtId[_numInList];
  int i;
  for (i=0;i<_numInList;i++){
    _list[i]=set1.getElem(i);
  }

}
EvtIdSet::EvtIdSet(const EvtIdSet& set1, const EvtIdSet& set2){

  _numInList=set1.sizeOfSet();
  _list=new EvtId[_numInList];
  int i;
  for (i=0;i<_numInList;i++){
    _list[i]=set1.getElem(i);
  }
  //then just append the second list.
  this->append(set2);

}

int EvtIdSet::contains(const EvtId id){

  int i;
  for (i=0;i<_numInList;i++){
    if ( _list[i] == id ) return 1;
  }
  
  return 0;
}

int EvtIdSet::contains(const std::string nm){

  int i;
  for (i=0;i<_numInList;i++){
    if ( _list[i] == EvtPDL::getId(nm) ) return 1;
  }
  
  return 0;
}


void EvtIdSet::append(const EvtIdSet set1){

  int combLen=_numInList+set1.sizeOfSet();
  int uniqueLen=0;
  EvtId *combSet;

  combSet=new EvtId[combLen];

  int i;
  for (i=0;i<combLen;i++){
    if ( i>=_numInList ) {

      //check that there are no overlaps between lists
      int j;
      int isUnique=1;
      for (j=0;j<_numInList;j++){
	if ( _list[j]==set1.getElem(i-_numInList) ) {
	  isUnique=0;
	}
      }
      if ( isUnique==1 ) {
	combSet[uniqueLen]=set1.getElem(i-_numInList);
	uniqueLen+=1;
      }
    }
    else{
      combSet[uniqueLen]=_list[i];
      uniqueLen+=1;
    }

    delete _list;
    _list=new EvtId[uniqueLen];

    _numInList=uniqueLen;
    for (i=0;i<_numInList;i++){
      _list[i]=combSet[i];
    }

    delete combSet;

  }
}

int EvtIdSet::sizeOfSet() const { return _numInList;}

EvtId EvtIdSet::getElem(const int i) const { return _list[i];}








