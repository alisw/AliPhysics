#ifndef ALIHBTWeightSPID_H
#define ALIHBTWeightSPID_H
/////////////////////////////////////////////////////////////
//
//This class introduces the weights calculated according 
//with functions of efficiency of identification (TPC+TOF) 
//(calculated by B.V. Batyunia).
//
//Author: Ludmila Malinina, JINR (malinina@sunhe.jinr.ru)
//
/////////////////////////////////////////////////////////////

#include <TObject.h>
class TF1;
class TH1;
class AliHBTPair;
class AliHBTWeightsPID: public TObject
 {
 public:
   AliHBTWeightsPID();
   AliHBTWeightsPID(const AliHBTWeightsPID &source):TObject(source) {
     //Copy ctor needed by the coding conventions but not used
     Fatal("AliHBTWeightsPID","copy ctor not implemented");
   }
   AliHBTWeightsPID & operator=(const AliHBTWeightsPID &/*source*/) {
     //Assignment operator needed by the coding conventions but not used
     Fatal("AliHBTWeightsPID","assignment operator not implemented");
     return * this;
   }
   virtual ~AliHBTWeightsPID(){;}
   static AliHBTWeightsPID* Instance();
   
   Double_t GetWeightPID(const AliHBTPair* trackpair); //get weight calculated Batyunia's  algorithm
   
 protected:
   Float_t fEfficTPC1; // ...?
   Float_t fEfficTPC2; // ...?
   Float_t fEfficTOF1; // ...?
   Float_t fEfficTOF2; // ...?
   
   static AliHBTWeightsPID *fgWeightsPID;// pointer to wrapper of Fortran Lednicky code   
   TH1 *fPtK;   //comment?
   TH1 *fPtKefftpc;//comment?
   TH1 *fPtKefftpcboth;//comment?
   TF1 *fEffic1pol;//comment?
   TF1 *fEffic2pol;//comment?
   TF1 *fEffic3pol;//comment?
   TF1 *fEffic4pol;//comment?
   
   TF1 *fEffic1polTOF;//comment?
   TF1 *fEffic2polTOF;//comment?
   TF1 *fEffic3polTOF;//comment?
   TF1 *fEffic4polTOF;//comment?
   
   ClassDef(AliHBTWeightsPID,2)
 };

#endif  
