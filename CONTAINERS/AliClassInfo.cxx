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
Revision 1.1  2000/11/01 16:01:22  kowal2
Classes for handling the new hits structures

*/
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  AliClassInfo                                                             //
//                                                                           //
//Defined to make unified interface to primitive types (in Root described by //
// TDataType) and TClass.                                                    //
// Additional virtual function (comparing to ROOT  TClass) neccesary         //
//                        for AliContainers                                  //
//   virtual void CTORBuffer(void * pointer, UInt_t size=1)                  //
//       should construct buffer of size =size  objects at position pointer  //
//   virtual void DTORBuffer(void * pointer, UInt_t size=1)                  //
//         should destruct buffer of size =size  objects at position pointer //
//   virtual void StreamBuffer(TBuffer& b, const void *object, UInt_t size)  //
//         stream buffer of objects   `
//Begin_Html
//<img src="../gif/AliClassInfo.gif">
//End_Html
//                                                                           //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////


#include "AliClassInfo.h"
#include "TMath.h"
#include "TClass.h"

#include "TROOT.h"
#include "iostream.h"
#include "AliDataType.h"

ClassImp(AliClassInfo)


TList AliClassInfo::fgList; // list of loaded class

AliClassInfo * AliClassInfo::FindClassInfo(const char * name)
{
  //
  TIter      next(&fgList);
  TString sname(name);
  AliClassInfo * info;
  while ((info = (AliClassInfo *) next())) 
    if (info->GetName()==sname) return info;
  return 0;
}



AliClassInfo * AliClassInfo::GenerClassInfo(const char * classname)
{
   //
  // Set class information fClassInfo  according class name
  //  
  char name[100];
  AliClassInfo *info=0;
  sprintf(name,"AliClassType%s",classname);
  info = AliClassInfo::FindClassInfo(classname);
  if (info) return info;
  //
  if ( (!info) &&  (gROOT->GetType(classname),kTRUE)){ 
    //if data type information exist
    //    char line[100];
    //    sprintf(line,"(*((AliClassInfo**)%p))= new AliDataType(\"%s\");",
    //	    &info,classname);
    // sprintf(line,"new AliDataType(\"%s\");",
    //    classname);

    //cout<<line<<"\n";
    //   gROOT->ProcessLine(line);
    new AliDataType(classname);
    info = AliClassInfo::FindClassInfo(classname);
  }   
  if (info) return info;

  TClass * cl = gROOT->GetClass(classname);
  // if exist root class information 
  if ( (!info) && (gROOT->GetClass(classname))){  //is it class?           
    char chinter[1000];
    sprintf(chinter,"%s.C",classname);
    GenerClassInfoCode(classname,kTRUE,cl->GetDeclFileName(),chinter);
    info = AliClassInfo::FindClassInfo(classname);
  }
  if (!info){ 
    TClass * cl = new TClass(classname,0);
    if (cl->Size()>0){ //if root information doesn't exist 
      //but exist cint information
      char chinclude[1000];
      sprintf(chinclude,"%s.h",classname);
      char chinter[1000];
      sprintf(chinter,"%s.C",classname);
      GenerClassInfoCode(classname,kTRUE,chinclude,chinter);
    }
  }
  
  return info;    
}



void AliClassInfo::GenerClassInfoCode(const char * clname, Bool_t load,
				      const char *incpath, const char *outfile)
{

  // gener temporary file - name

  FILE *fout = fopen(outfile,"w");
  if (!clname){
    cerr<<"Class not  specified\n";
    return ;
  }
  char buff[1000];
  const char *pchar =incpath;   
  //replace  with /0
  char * pchar2 =buff; 
  fprintf(fout,"#include \"AliClassInfo.h\"\n");          
  
  if (incpath==0) {}
  else{
    // proces headers - header separated by     
    pchar =incpath;       
    pchar2 =buff; 
    // 
    while (*pchar==' ') pchar++;
    while (*pchar) {
      if (*pchar!=' ') *pchar2++ = *pchar;
      else
	if (*(pchar2-1)!=0) *pchar2++=0;	 
      pchar++;
    }
    *pchar2=0;
    Int_t index = pchar2-buff;   
    for (Int_t i=0;i<index;i++) 
      if ( (i==0) ||(buff[i-1]==0))  
	fprintf(fout,"#include \"%s\"\n",&buff[i]);          
  }
  //process classes
  pchar =clname;
  pchar2 =buff; 
  while (*pchar==' ') pchar++;
  while (*pchar) {
    if (*pchar!=' ') *pchar2++ = *pchar;
    else
      if (*(pchar2-1)!=0) *pchar2++=0;	 
    pchar++;
  }
  *pchar2=0;
  Int_t index = pchar2-buff;   
  for (Int_t i=0;i<index;i++) 
    if ( (i==0) ||(buff[i-1]==0))  
      GenerClassInterface(&buff[i],fout);

  fclose(fout);
  //
  
  if (load) {
    char line[100];
    // gSystem->Rename("/tmp/root_tmpinterinter"
    sprintf(line,".L %s+",outfile);
    cout<<line<<"\nGenerating class Interface \n";
    cout<<line<<"\n*****************************\n";   
    gROOT->ProcessLine(line);
    cout<<line<<"\n*****************************\n";
    cout<<line<<"\nClass Interface generated \n";
  }
}



Bool_t  AliClassInfo::GenerClassInterface(const char * clname, FILE * fout)
{
  //  TClass * cl = gROOT->GetClass("AliLHit",kTRUE);
  fprintf(fout,"\n/************************************************/\n");
  fprintf(fout,"/* Automaticaly generated interface for class     \n");
  fprintf(fout,"                 %s                                \n",clname);
  fprintf(fout,"**************************************************/\n");
  fprintf(fout,"\n\n");
  //constructor
  fprintf(fout,"class AliClass%s : public AliClassInfo {\n",clname);
  fprintf(fout,"public:\n");
  fprintf(fout,"\tAliClass%s(){\n",clname);
  fprintf(fout,"\t  SetName(\"%s\");\n",clname);
  fprintf(fout,"\t  SetTitle(\"Interface for %s class \");\n",clname);
  fprintf(fout,"\t  fgList.Add(this);\n");
  fprintf(fout,"\t  fSize = sizeof(%s);\n\t}\n",clname);
  //
  fprintf(fout,"\tconst char * GetClassName(){ return \"%s\";}\n",clname); 
  //
  fprintf(fout,"\tvirtual TClass* GetClass(){return %s::Class();}\n",clname);
 //placement constructor interface
  fprintf(fout,"\tvoid CTORBuffer(void * pointer, UInt_t size=1)\n\t{\n");
  fprintf(fout,"\t  %s * last = &((%s*)pointer)[size];\n",clname,clname);
  fprintf(fout,"\t  %s * p = (%s*)pointer;\n",clname,clname);
  fprintf(fout,"\t  while (p!=last) new (p++)%s;\n\t}\n",clname);
  //placement destructor interface
  fprintf(fout,"\tvoid DTORBuffer(void * pointer, UInt_t size=1)\n\t{\n");
  fprintf(fout,"\t  %s * last = &((%s*)pointer)[size];\n",clname,clname);
  fprintf(fout,"\t  %s * p = (%s*)pointer;\n",clname,clname);
  fprintf(fout,"\t  while (p!=last) (p++)->~%s();\n\t}\n",clname);
  //streamer interface
  fprintf(fout,"\tvoid StreamBuffer(TBuffer &b,const void * pointer, UInt_t size=1)\n\t{\n");
  fprintf(fout,"\t  for (UInt_t i=0;i<size;i++) ((%s*)pointer)[i].Streamer(b);\n\t}\n",clname);
  //
  fprintf(fout,"\t  void ObjectDump(void *p) {((%s*)p)->Dump();}\n",clname);
  fprintf(fout,"};\n");
  //make instance of the class
  fprintf(fout,"AliClass%s galiclass____%s; \n",clname,clname);  
  return kTRUE;
}








 
