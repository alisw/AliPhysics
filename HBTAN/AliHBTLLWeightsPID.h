#ifndef ALIHBTLLWEIGHTSPID_H
#define ALIHBTLLWEIGHTSPID_H
/////////////////////////////////////////////////////////////
//
//This class introduces the weights calculated according 
//with functions of efficiency of identification (TPC+TOF) 
//(calculated by B.V. Batyunia).
//
//Author: Ludmila Malinina, JINR (malinina@sunhe.jinr.ru)
//
/////////////////////////////////////////////////////////////

#include "TH1.h"                                                                                  
#include "TH2.h"                                                                                  
#include <TF1.h>                                                                                  
#include "TCanvas.h"                                                                              
#include "TPad.h"                                                                                 

#include <TObject.h>

class AliHBTPair;
class AliHBTLLWeightsPID: public TObject
 {
   public:
     AliHBTLLWeightsPID();
     virtual ~AliHBTLLWeightsPID(){;}
     static AliHBTLLWeightsPID* Instance();

     Double_t GetWeightPID(const AliHBTPair* trackpair); //get weight calculated Batyunia's  algorithm
     Float_t efficTPC1,efficTPC2,efficTOF1,efficTOF2;     
    
   protected:
                                                                                             
   static AliHBTLLWeightsPID *fWeightsPID;// pointer to wrapper of Fortran Lednicky code   
   TH1 *ptK;   //comment?
   TH1 *ptKefftpc;//comment?
   TH1 *ptKefftpcboth;//comment?
   TF1 *effic1pol;//comment?
   TF1 *effic2pol;//comment?
   TF1 *effic3pol;//comment?
   TF1 *effic4pol;//comment?
   
   TF1 *effic1polTOF;//comment?
   TF1 *effic2polTOF;//comment?
   TF1 *effic3polTOF;//comment?
   TF1 *effic4polTOF;//comment?
	   
   private:

   public:
     ClassDef(AliHBTLLWeightsPID,1)
 };

#endif  
