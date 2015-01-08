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
// Module: EvtGen/EvtSpinDensity.hh
//
// Description: Templated class to implement a hash table on string
//
// Modification history:
//
//    RYD     Aug 2, 2006         Module created
//
//------------------------------------------------------------------------


#ifndef EVTSTRINGHASH_HH
#define EVTSTRINGHASH_HH

#include <string>


template<class T> 
class EvtStringHash{

public:
  inline EvtStringHash(int size);
  inline void add(const std::string& str,T* data);
  inline T* get(const std::string& str);
  inline ~EvtStringHash();

private:
  EvtStringHash();
  int _size;
  inline int hash(const std::string& str);
  std::string*** _strings;
  T*** _data;
  int* _entries;

};


template<class T> 
EvtStringHash<T>::EvtStringHash(int size){

  _size=size;

  typedef std::string** EvtStringPtrPtr;
  typedef T** TPtrPtr;

  _strings=new EvtStringPtrPtr[_size];
  _data=new TPtrPtr[_size];
  _entries=new int[_size];

  int i;

  for(i=0;i<_size;i++){
    _entries[i]=0;
  }

}

template<class T> 
EvtStringHash<T>::~EvtStringHash(){
  
  int i;
  for(i=0;i<_size;i++){
    int j;
    for(j=0;j<_entries[i];j++){
      delete _strings[i][j];
    }    
    if (_entries[i]>0){
      delete [] _strings[i];
      delete [] _data[i];
    }
  }

  delete [] _strings;
  delete [] _data;
  delete [] _entries;

}

template<class T> 
void EvtStringHash<T>::add(const std::string& str,T* data){

  int ihash=hash(str);

  typedef std::string* EvtStringPtr;
  typedef T* TPtr;

  std::string** newstrings=new EvtStringPtr[_entries[ihash]+1];
  T** newdata=new TPtr[_entries[ihash]+1];
				
  int i;

  for(i=0;i<_entries[ihash];i++){
    newstrings[i]=_strings[ihash][i];
    newdata[i]=_data[ihash][i];
  }

  newstrings[_entries[ihash]]=new std::string;
  *(newstrings[_entries[ihash]])=str;
  newdata[_entries[ihash]]=data;
  

  if(_entries[ihash]!=0){
    delete [] _strings[ihash];
    delete [] _data[ihash];
  }
  
  _entries[ihash]++;

  _strings[ihash]=newstrings;
  _data[ihash]=newdata;

}

template<class T> 
T* EvtStringHash<T>::get(const std::string& str){

  int ihash=hash(str);
  
  int i;

  for (i=0;i<_entries[ihash];i++){
    if (*(_strings[ihash][i])==str) return _data[ihash][i];
  }
  
  return 0;


}


template<class T> 
int EvtStringHash<T>::hash(const std::string& str){

  const char* cstr=str.c_str();

  int i=0;

  int value=0;

  while(cstr[i]!=0){
    value+=(int)cstr[i];
    i++;
  }

  return value%_size;

}


#endif





