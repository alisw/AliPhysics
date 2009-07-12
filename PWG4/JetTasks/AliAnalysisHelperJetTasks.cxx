
#include "TROOT.h"
#include "TList.h"
#include "AliMCEvent.h"
#include "AliAODJet.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include <fstream>
#include <iostream>
#include "AliAnalysisHelperJetTasks.h"


ClassImp(AliAnalysisHelperJetTasks)



 
AliGenPythiaEventHeader*  AliAnalysisHelperJetTasks::GetPythiaEventHeader(AliMCEvent *mcEvent){
  
  AliGenEventHeader* genHeader = mcEvent->GenEventHeader();
  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
  if(!pythiaGenHeader){
    // cocktail ??
    AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(genHeader);
    
    if (!genCocktailHeader) {
      Printf("%s %d: Unknown header type (not Pythia or Cocktail)",(char*)__FILE__,__LINE__);
      return 0;
    }
    TList* headerList = genCocktailHeader->GetHeaders();
    for (Int_t i=0; i<headerList->GetEntries(); i++) {
      pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(headerList->At(i));
      if (pythiaGenHeader)
        break;
    }
    if(!pythiaGenHeader){
      Printf("%s %d: PythiaHeader not found!",(char*)__FILE__,__LINE__);
      return 0;
    }
  }
  return pythiaGenHeader;

}


void AliAnalysisHelperJetTasks::PrintStack(AliMCEvent *mcEvent,Int_t iFirst,Int_t iLast,Int_t iMaxPrint){

  AliStack *stack = mcEvent->Stack();
  if(!stack){
    Printf("%s%d No Stack available",(char*)__FILE__,__LINE__);
    return;
  }

  static Int_t iCount = 0;
  if(iCount>iMaxPrint)return;
  Int_t nStack = stack->GetNtrack();
  if(iLast == 0)iLast = nStack;
  else if(iLast > nStack)iLast = nStack;


  Printf("####################################################################");
  for(Int_t np = iFirst;np<iLast;++np){
    TParticle *p = stack->Particle(np);
    Printf("Nr.%d --- Status %d ---- Mother1 %d Mother2 %d Daughter1 %d Daughter2 %d ",
	   np,p->GetStatusCode(),p->GetMother(0),p->GetMother(1),p->GetDaughter(0),p->GetDaughter(1));
    Printf("Eta %3.3f Phi %3.3f ",p->Eta(),p->Phi()); 
    p->Print();    
    Printf("---------------------------------------");
  }
  iCount++;
}




void AliAnalysisHelperJetTasks::GetClosestJets(AliAODJet *genJets,const Int_t &nGenJets,
					       AliAODJet *recJets,const Int_t &nRecJets,
					       Int_t *iGenIndex,Int_t *iRecIndex,
					       Int_t iDebug,Float_t maxDist){

  //
  // Relate the two input jet Arrays
  //

  //
  // The association has to be unique
  // So check in two directions
  // find the closest rec to a gen
  // and check if there is no other rec which is closer
  // Caveat: Close low energy/split jets may disturb this correlation

  // Idea: search in two directions generated e.g (a--e) and rec (1--3)
  // Fill a matrix with Flags (1 for closest rec jet, 2 for closest rec jet
  // in the end we have something like this
  //    1   2   3
  // ------------
  // a| 3   2   0
  // b| 0   1   0
  // c| 0   0   3
  // d| 0   0   1
  // e| 0   0   1
  // Topology
  //   1     2
  //     a         b        
  //
  //  d      c
  //        3     e
  // Only entries with "3" match from both sides
  
  const int kMode = 3;

  

  if(nRecJets==0||nGenJets==0)return;

  for(int i = 0;i < nGenJets;++i)iGenIndex[i] = -1;
  for(int j = 0;j < nRecJets;++j)iRecIndex[j] = -1;

  UShort_t *iFlag = new UShort_t[nGenJets*nRecJets];
  for(int i = 0;i < nGenJets;++i){
    for(int j = 0;j < nRecJets;++j){
      iFlag[i*nGenJets+j] = 0;
    }
  }



  // find the closest distance to the generated
  for(int ig = 0;ig<nGenJets;++ig){
    Float_t dist = maxDist;
    if(iDebug>1)Printf("Gen (%d) p_T %3.3f eta %3.3f ph %3.3f ",ig,genJets[ig].Pt(),genJets[ig].Eta(),genJets[ig].Phi());
    for(int ir = 0;ir<nRecJets;++ir){
      Double_t dR = genJets[ig].DeltaR(&recJets[ir]);
      if(iDebug>1)Printf("Rec (%d) p_T %3.3f eta %3.3f ph %3.3f ",ir,recJets[ir].Pt(),recJets[ir].Eta(),recJets[ir].Phi());
      if(iDebug>1)Printf("Distance (%d)--(%d) %3.3f ",ig,ir,dR);
      if(dR<dist){
	iRecIndex[ig] = ir;
	dist = dR;
      }	
    }
    if(iRecIndex[ig]>=0)iFlag[ig*nGenJets+iRecIndex[ig]]+=1;
    // reset...
    iRecIndex[ig] = -1;
  }
  // other way around
  for(int ir = 0;ir<nRecJets;++ir){
    Float_t dist = maxDist;
    for(int ig = 0;ig<nGenJets;++ig){
      Double_t dR = genJets[ig].DeltaR(&recJets[ir]);
      if(dR<dist){
	iGenIndex[ir] = ig;
	dist = dR;
      }	
    }
    if(iGenIndex[ir]>=0)iFlag[iGenIndex[ir]*nGenJets+ir]+=2;
    // reset...
    iGenIndex[ir] = -1;
  }

  // check for "true" correlations

  if(iDebug>1)Printf(">>>>>> Matrix");

  for(int ig = 0;ig<nGenJets;++ig){
    for(int ir = 0;ir<nRecJets;++ir){
      // Print
      if(iDebug>1)printf("Flag[%d][%d] %d ",ig,ir,iFlag[ig*nGenJets+ir]);

      if(kMode==3){
	// we have a uniqie correlation
	if(iFlag[ig*nGenJets+ir]==3){
	  iGenIndex[ir] = ig;
	  iRecIndex[ig] = ir;
	}
      }
      else{
	// we just take the correlation from on side
	if((iFlag[ig*nGenJets+ir]&2)==2){
	  iGenIndex[ir] = ig;
	}
	if((iFlag[ig*nGenJets+ir]&1)==1){
	  iRecIndex[ig] = ir;
	}
      }
    }
    if(iDebug>1)printf("\n");
  }

  delete [] iFlag;

}





