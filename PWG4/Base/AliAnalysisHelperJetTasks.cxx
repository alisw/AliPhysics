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

//
// Helper Class that contains a lot of 
// usefull static functions jet matchin pythia access etc.
//
// Author: C. Klein-Boesing IKP Muenster 

#include "TROOT.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TList.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TSeqCollection.h"
#include "TMethodCall.h"
#include "TFile.h"
#include "TString.h"
#include "TArrayI.h"
#include "TArrayS.h"
#include "TArrayF.h"

#include "AliMCEvent.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliAODJet.h"
#include "AliAODEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include <AliGenDPMjetEventHeader.h>
#include "AliGenPythiaEventHeader.h"
#include <fstream>
#include <iostream>
#include "AliAnalysisHelperJetTasks.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector.h"

ClassImp(AliAnalysisHelperJetTasks)

Int_t AliAnalysisHelperJetTasks::fgLastProcessType = -1;

 
AliGenPythiaEventHeader*  AliAnalysisHelperJetTasks::GetPythiaEventHeader(const AliMCEvent *mcEvent){
  //
  // Fetch  the pythia event header
  // 
  if(!mcEvent)return 0;
  AliGenEventHeader* genHeader = mcEvent->GenEventHeader();
  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
  if(!pythiaGenHeader){
    // cocktail ??
    AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(genHeader);
    
    if (!genCocktailHeader) {
      AliWarningGeneral(Form(" %s:%d",(char*)__FILE__,__LINE__),"Unknown header type (not Pythia or Cocktail)");
      //      AliWarning(Form("%s %d: Unknown header type (not Pythia or Cocktail)",(char*)__FILE__,__LINE__));
      return 0;
    }
    TList* headerList = genCocktailHeader->GetHeaders();
    for (Int_t i=0; i<headerList->GetEntries(); i++) {
      pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(headerList->At(i));
      if (pythiaGenHeader)
        break;
    }
    if(!pythiaGenHeader){
      AliWarningGeneral(Form(" %s:%d",(char*)__FILE__,__LINE__),"Pythia event header not found");
      return 0;
    }
  }
  return pythiaGenHeader;

}


void AliAnalysisHelperJetTasks::PrintStack(AliMCEvent *mcEvent,Int_t iFirst,Int_t iLast,Int_t iMaxPrint){
  // 
  // Print the stack informatuin up to the iMaxPrint event
  //

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




void AliAnalysisHelperJetTasks::GetClosestJets(AliAODJet *genJets,const Int_t &kGenJets,
					       const AliAODJet *recJets,const Int_t &kRecJets,
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

  // In case we have more jets than kmaxjets only the 
  // first kmaxjets are searched
  // all other are -1
  // use kMaxJets for a test not to fragemnt the memory...

  for(int i = 0;i < kRecJets;++i)iGenIndex[i] = -1;
  for(int j = 0;j < kGenJets;++j)iRecIndex[j] = -1;


  
  const int kMode = 3;

  const Int_t nGenJets = TMath::Min(kMaxJets,kGenJets);
  const Int_t nRecJets = TMath::Min(kMaxJets,kRecJets);

  if(nRecJets==0||nGenJets==0)return;

  // UShort_t *iFlag = new UShort_t[nGenJets*nRecJets];
  UShort_t iFlag[kMaxJets*kMaxJets] = {0}; // all values set to zero
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
      if(iDebug>1){
	printf("Rec (%d) ",ir);
	Printf("p_T %3.3f eta %3.3f ph %3.3f ",recJets[ir].Pt(),recJets[ir].Eta(),recJets[ir].Phi());
      }    
      Double_t dR = genJets[ig].DeltaR(&recJets[ir]);
      if(iDebug>1)Printf("Distance (%d)--(%d) %3.3f ",ig,ir,dR);
      if(dR<dist){
	iRecIndex[ig] = ir;
	dist = dR;
      }	
    }
    if(iRecIndex[ig]>=0)iFlag[ig*nRecJets+iRecIndex[ig]]+=1;
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
    if(iGenIndex[ir]>=0)iFlag[iGenIndex[ir]*nRecJets+ir]+=2;
    // reset...
    iGenIndex[ir] = -1;
  }

  // check for "true" correlations

  if(iDebug>1)Printf(">>>>>> Matrix");

  for(int ig = 0;ig<nGenJets;++ig){
    for(int ir = 0;ir<nRecJets;++ir){
      // Print
      if(iDebug>1)printf("Flag[%d][%d] %d ",ig,ir,iFlag[ig*nRecJets+ir]);

      if(kMode==3){
	// we have a uniqie correlation
	if(iFlag[ig*nRecJets+ir]==3){
	  iGenIndex[ir] = ig;
	  iRecIndex[ig] = ir;
	}
      }
      else{
	// we just take the correlation from on side
	if((iFlag[ig*nRecJets+ir]&2)==2){
	  iGenIndex[ir] = ig;
	}
	if((iFlag[ig*nRecJets+ir]&1)==1){
	  iRecIndex[ig] = ir;
	}
      }
    }
    if(iDebug>1)printf("\n");
  }
}



void AliAnalysisHelperJetTasks::GetClosestJets(const TList *genJetsList,const Int_t &kGenJets,
					       const TList *recJetsList,const Int_t &kRecJets,
					       TArrayI &iGenIndex,TArrayI &iRecIndex,
					       Int_t iDebug,Float_t maxDist){

  // Size indepnedendentt Implemenation of jet matching
  // Thepassed TArrayI should be static in the user function an only increased if needed

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

  iGenIndex.Reset(-1);
  iRecIndex.Reset(-1);
  
  const int kMode = 3;
  const Int_t nGenJets = TMath::Min(genJetsList->GetEntries(),kGenJets);
  const Int_t nRecJets = TMath::Min(recJetsList->GetEntries(),kRecJets);
  if(nRecJets==0||nGenJets==0)return;

  static TArrayS iFlag(nGenJets*nRecJets);
  if(iFlag.GetSize()<(nGenJets*nRecJets)){
    iFlag.Set(nGenJets*nRecJets);
  }
  iFlag.Reset(0);

  // find the closest distance to the generated
  for(int ig = 0;ig<nGenJets;++ig){
    AliAODJet *genJet = (AliAODJet*)genJetsList->At(ig); 
    if(!genJet)continue;

    Float_t dist = maxDist;
    if(iDebug>1)Printf("Gen (%d) p_T %3.3f eta %3.3f ph %3.3f ",ig,genJet->Pt(),genJet->Eta(),genJet->Phi());
    for(int ir = 0;ir<nRecJets;++ir){
      AliAODJet *recJet = (AliAODJet*)recJetsList->At(ir); 
      if(!recJet)continue;
      if(iDebug>1){
	printf("Rec (%d) ",ir);
	Printf("p_T %3.3f eta %3.3f ph %3.3f ",recJet->Pt(),recJet->Eta(),recJet->Phi());
      }    
      Double_t dR = genJet->DeltaR(recJet);
      if(iDebug>1)Printf("Distance (%d)--(%d) %3.3f ",ig,ir,dR);
      if(dR<dist){
	iRecIndex[ig] = ir;
	dist = dR;
      }	
    }
    if(iRecIndex[ig]>=0)iFlag[ig*nRecJets+iRecIndex[ig]]+=1;
    // reset...
    iRecIndex[ig] = -1;
  }
  // other way around

  for(int ir = 0;ir<nRecJets;++ir){
      AliAODJet *recJet = (AliAODJet*)recJetsList->At(ir); 
      if(!recJet)continue;
      Float_t dist = maxDist;
      for(int ig = 0;ig<nGenJets;++ig){
	AliAODJet *genJet = (AliAODJet*)genJetsList->At(ig); 
	if(!genJet)continue;
	Double_t dR = genJet->DeltaR(recJet);
	if(dR<dist){
	iGenIndex[ir] = ig;
	dist = dR;
      }	
    }
    if(iGenIndex[ir]>=0)iFlag[iGenIndex[ir]*nRecJets+ir]+=2;
    // reset...
    iGenIndex[ir] = -1;
  }

  // check for "true" correlations

  if(iDebug>1)Printf(">>>>>> Matrix Size %d",iFlag.GetSize());

  for(int ig = 0;ig<nGenJets;++ig){
    for(int ir = 0;ir<nRecJets;++ir){
      // Print
      if(iDebug>1)printf("Flag2[%d][%d] %d ",ig,ir,iFlag[ig*nRecJets+ir]);

      if(kMode==3){
	// we have a uniqie correlation
	if(iFlag[ig*nRecJets+ir]==3){
	  iGenIndex[ir] = ig;
	  iRecIndex[ig] = ir;
	}
      }
      else{
	// we just take the correlation from on side
	if((iFlag[ig*nRecJets+ir]&2)==2){
	  iGenIndex[ir] = ig;
	}
	if((iFlag[ig*nRecJets+ir]&1)==1){
	  iRecIndex[ig] = ir;
	}
      }
    }
    if(iDebug>1)printf("\n");
  }
}

void AliAnalysisHelperJetTasks::GetJetMatching(const TList *genJetsList, const Int_t &kGenJets,
                                               const TList *recJetsList, const Int_t &kRecJets,
                                               TArrayI &iMatchIndex, TArrayF &fPtFraction,
                                               Int_t iDebug, Float_t maxDist, Int_t mode){

                                            
    // Matching jets from two lists
    // Start with highest energetic jet from first list (generated/embedded)
    // Calculate distance (\Delta R) to all jets from second list (reconstructed)
    // Select N closest jets = kClosestJetsN
    // Check energy fraction from jets from first list in jets from second list
    // Matched jets = jet with largest energy fraction
    // Store index of matched jet in TArrayI iMatchIndex
                                            
    // reset index
    iMatchIndex.Reset(-1);
    fPtFraction.Reset(-1.);
    
    // N closest jets: store list with index and \Delta R
    const Int_t kClosestJetsN = 4; 
    Double_t closestJets[kClosestJetsN][2]; //[][0] = index, [][1] = \Delta R
        
    const Int_t nGenJets = TMath::Min(genJetsList->GetEntries(),kGenJets);
    const Int_t nRecJets = TMath::Min(recJetsList->GetEntries(),kRecJets);
    if(nRecJets==0||nGenJets==0) return;
    
    AliAODJet *genJet = 0x0;
    AliAODJet *recJet = 0x0;
    
    // loop over generated/embedded jets
    for(Int_t ig=0; ig<nGenJets; ++ig){

        for(Int_t i=0; i<kClosestJetsN; ++i){
            closestJets[i][0] = -1;   // index
            closestJets[i][1] = 1e6;  // delta R
        }

        genJet = (AliAODJet*)genJetsList->At(ig);
        //if(!genJet || !JetSelected(genJet)) continue;
        if(!genJet) continue;
        
        // find N closest reconstructed jets
        Double_t deltaR = 0.;
        for(Int_t ir=0; ir<nRecJets; ++ir){
            recJet = (AliAODJet*)recJetsList->At(ir);
            //if(!recJet || !JetSelected(recJet)) continue;
            if(!recJet) continue;
            
            deltaR = genJet->DeltaR(recJet);
            
            Int_t i=kClosestJetsN-1;
            if(deltaR<closestJets[i][1] && deltaR<maxDist){
                closestJets[i][0] = (Double_t) ir; // index
                closestJets[i][1] = deltaR;
                
                // sort array (closest at the beginning)
                while(i>=1 && closestJets[i][1]<closestJets[i-1][1]){
                    Double_t tmpArr[2];
                    for(Int_t j=0; j<2; j++){
                       tmpArr[j] = closestJets[i-1][j];
                       closestJets[i-1][j]   = closestJets[i][j];
                       closestJets[i][j] = tmpArr[j];
                    }
                    i--;
                }
            } 
        } // end: loop over reconstructed jets
        
        // calculate fraction for the N closest jets
        Double_t maxFraction = -1.; // maximum found fraction in one jets
        Double_t cumFraction = 0.; // cummulated fraction of closest jets (for break condition)
        Double_t fraction = 0.;
        Int_t ir = -1;  // index of close reconstruced jet
        
        for(Int_t irc=0; irc<kClosestJetsN; irc++){
            ir = (Int_t)(closestJets[irc][0]);
			if(ir<0 || ir>nRecJets-1) continue;
            recJet = (AliAODJet*)recJetsList->At(ir);
            if(!(recJet)) continue;
            
            fraction = GetFractionOfJet(recJet,genJet,mode);
            
            cumFraction += fraction;
            
            // check if jet fullfills current matching condition
            if(fraction>maxFraction){
                // avoid multiple links
                for(Int_t ij=0; ij<ig; ++ij){
                    if(iMatchIndex[ij]==ir) continue;
                }
                // set index
                maxFraction = fraction;
                fPtFraction[ig] = fraction;                
                iMatchIndex[ig] = ir;
            }
            // break condition: less energy left as already found in one jet or
            // as required for positiv matching
            if(1-cumFraction<maxFraction) break;
        } // end: loop over closest jets
        
        if(iMatchIndex[ig]<0){
            if(iDebug) Printf("Matching failed for (gen) jet #%d", ig);
        }
    }
}

Double_t AliAnalysisHelperJetTasks::GetFractionOfJet(const AliAODJet *recJet, const AliAODJet *genJet, Int_t mode){
  //
  // get the fraction of hte signal jet in the full jt
  //
    Double_t ptGen = genJet->Pt();
    if(ptGen==0.) return 999.;
    
    Double_t ptAssocTracks = 0.; // sum of pT of tracks found in both jets
    
    // look at tracks inside jet
    TRefArray *genTrackList = genJet->GetRefTracks();
    TRefArray *recTrackList = recJet->GetRefTracks();
    Int_t nTracksGenJet = genTrackList->GetEntriesFast();
    Int_t nTracksRecJet = recTrackList->GetEntriesFast();
    
    AliAODTrack* recTrack;
    AliAODTrack* genTrack;
    for(Int_t ir=0; ir<nTracksRecJet; ++ir){
        recTrack = (AliAODTrack*)(recTrackList->At(ir));
        if(!recTrack) continue;
        
        for(Int_t ig=0; ig<nTracksGenJet; ++ig){
            genTrack = (AliAODTrack*)(genTrackList->At(ig));
            if(!genTrack) continue;
            
            // look if it points to the same track
            if( (mode&1)!=0 && genTrack==recTrack){
                ptAssocTracks += genTrack->Pt();
                break;
            }
 
            if( (mode&2)!=0 
                    && genTrack->GetLabel()>-1
                    && recTrack->GetLabel()>-1
                    && genTrack->GetLabel()==recTrack->GetLabel()){

                ptAssocTracks += genTrack->Pt();
                break; 
            }
        }
    }
    
    // calculate fraction
    Double_t fraction = ptAssocTracks/ptGen;
    
    return fraction;
}


void  AliAnalysisHelperJetTasks::MergeOutputDirs(const char* cFiles,const char* cPattern,const char *cOutFile,Bool_t bUpdate){
  // Routine to merge only directories containing the pattern
  //
  TString outFile(cOutFile);
  if(outFile.Length()==0)outFile = Form("%s.root",cPattern);
  ifstream in1;
  in1.open(cFiles);
  // open all files
  TList fileList;
  char cFile[240];
  while(in1>>cFile){// only open the first file
    Printf("Adding %s",cFile);
    TFile *f1 = TFile::Open(cFile);
    fileList.Add(f1);
  }

  TFile *fOut = 0;
  if(fileList.GetEntries()){// open the first file
    TFile* fIn = dynamic_cast<TFile*>(fileList.At(0));
    if(!fIn){
      Printf("Input File not Found");
      return;
    }
    // fetch the keys for the directories
    TList *ldKeys = fIn->GetListOfKeys();
    for(int id = 0;id < ldKeys->GetEntries();id++){
      // check if it is a directory
      TKey *dKey = (TKey*)ldKeys->At(id);
      TDirectory *dir = dynamic_cast<TDirectory*>(dKey->ReadObj());
      if(dir){
	TString dName(dir->GetName());
	if(dName.Contains(cPattern)){
	  // open new file if first match
	  if(!fOut){
	    if(bUpdate)fOut = new TFile(outFile.Data(),"UPDATE");
	    else fOut = new TFile(outFile.Data(),"RECREATE");
	  }
	  // merge all the stuff that is in this directory
	  TList *llKeys = dir->GetListOfKeys();
	  TList *tmplist;
	  TMethodCall callEnv;

	  fOut->cd();
	  TDirectory *dOut = fOut->mkdir(dir->GetName());

	  for(int il = 0;il < llKeys->GetEntries();il++){
	    TKey *lKey = (TKey*)llKeys->At(il);
	    TObject *object = dynamic_cast<TObject*>(lKey->ReadObj());
	    //  from TSeqCollections::Merge
	    if(!object)continue;
	    // If current object is not mergeable just skip it
	    if (!object->IsA()) {
	      continue;
	    }
	    callEnv.InitWithPrototype(object->IsA(), "Merge", "TCollection*");
	    if (!callEnv.IsValid()) {
	      continue;
	    }
	    // Current object mergeable - get corresponding objects in input lists
	    tmplist = new TList();
	    for(int i = 1; i < fileList.GetEntries();i++){
	      TDirectory *fInTmp = dynamic_cast<TDirectory*>(fileList.At(i)); 
	      if(!fInTmp){
		Printf("File %d not found",i);
		continue;
	      }
	      TDirectory *dTmp = (TDirectory*)fInTmp->Get(dName.Data());
	      if(!dTmp){
		Printf("Directory %s not found",dName.Data());
		continue;
	      }
	      TObject* oTmp = (TObject*)dTmp->Get(object->GetName());
	      if(!oTmp){
		Printf("Object %s not found",object->GetName());
		continue;
	      }
	      tmplist->Add(oTmp);
	    }
	    callEnv.SetParam((Long_t) tmplist);
	    callEnv.Execute(object);
	    delete tmplist;tmplist = 0;
	    dOut->cd();
	    object->Write(object->GetName(),TObject::kSingleKey);
	  }
	}
      }
    }
    if(fOut){
      fOut->Close();
    }
  }
}

void  AliAnalysisHelperJetTasks::MergeOutput(char* cFiles, char* cDir, char *cList,char *cOutFile,Bool_t bUpdate){

  // This is used to merge the analysis-output from different 
  // data samples/pt_hard bins
  // in case the eventweigth was set to xsection/ntrials already, this
  // is not needed. Both methods only work in case we do not mix different 
  // pt_hard bins, and do not have overlapping bins

  const Int_t nMaxBins = 12;
  // LHC08q jetjet100: Mean = 1.42483e-03, RMS = 6.642e-05
  // LHC08r jetjet50: Mean = 2.44068e-02, RMS = 1.144e-03
  // LHC08v jetjet15-50: Mean = 2.168291 , RMS = 7.119e-02
  // const Float_t xsection[nBins] = {2.168291,2.44068e-02};

  Float_t xsection[nMaxBins];
  Float_t nTrials[nMaxBins];
  Float_t sf[nMaxBins];
  TList *lIn[nMaxBins];
  TDirectory *dIn[nMaxBins];
  TFile *fIn[nMaxBins];

  ifstream in1;
  in1.open(cFiles);

  char cFile[120];
  Int_t ibTotal = 0;
  while(in1>>cFile){
    fIn[ibTotal] = TFile::Open(cFile);
    if(strlen(cDir)==0){
      dIn[ibTotal] = gDirectory;
    }
    else{
      dIn[ibTotal] = (TDirectory*)fIn[ibTotal]->Get(cDir);
    }
    if(!dIn[ibTotal]){
      Printf("%s:%d No directory %s found, exiting...",__FILE__,__LINE__,cDir);
      fIn[ibTotal]->ls();
      return;
    }

    lIn[ibTotal] = (TList*)dIn[ibTotal]->Get(cList);
    Printf("Merging file %s %s",cFile, cDir);
    if(!lIn[ibTotal]){
      Printf("%s:%d No list %s found, exiting...",__FILE__,__LINE__,cList);
      fIn[ibTotal]->ls();
      return;
    }
    TH1* hTrials = (TH1F*)lIn[ibTotal]->FindObject("fh1Trials");
    if(!hTrials){
      Printf("%s:%d fh1PtHard_Trials not found in list, exiting...",__FILE__,__LINE__);
      return;
    }
    TProfile* hXsec = (TProfile*)lIn[ibTotal]->FindObject("fh1Xsec");
    if(!hXsec){
      Printf("%s:%d fh1Xsec  not found in list, exiting...",__FILE__,__LINE__);
      return;
    }
    xsection[ibTotal] = hXsec->GetBinContent(1);
    nTrials[ibTotal] = hTrials->Integral();
    sf[ibTotal] = xsection[ibTotal]/ nTrials[ibTotal];
    ibTotal++;
  }

  if(ibTotal==0){
    Printf("%s:%d No files found for mergin, exiting",__FILE__,__LINE__);
    return;
  }

  TFile *fOut = 0;
  if(bUpdate)fOut = new TFile(cOutFile,"UPDATE");
  else fOut = new TFile(cOutFile,"RECREATE");
  TDirectory *dOut = fOut->mkdir(dIn[0]->GetName());
  dOut->cd();
  TList *lOut = new TList();
  lOut->SetName(lIn[0]->GetName());

  // for the start scale all...
  for(int ie = 0; ie < lIn[0]->GetEntries();++ie){
    TH1 *h1Add = 0;
    THnSparse *hnAdd = 0;
    for(int ib = 0;ib < ibTotal;++ib){
      // dynamic cast does not work with cint
      TObject *h = lIn[ib]->At(ie);
      if(h->InheritsFrom("TH1")){
	TH1 *h1 = (TH1*)h;
	if(ib==0){
	  h1Add = (TH1*)h1->Clone(h1->GetName());
	  h1Add->Scale(sf[ib]);
	}
	else{
	  h1Add->Add(h1,sf[ib]);
	}
      }
      else if(h->InheritsFrom("THnSparse")){
	THnSparse *hn = (THnSparse*)h;
	if(ib==0){
	  hnAdd = (THnSparse*)hn->Clone(hn->GetName());
	  hnAdd->Scale(sf[ib]);
	}
	else{
	  hnAdd->Add(hn,sf[ib]);
	}
      }
      

    }// ib
    if(h1Add)lOut->Add(h1Add);
    else if(hnAdd)lOut->Add(hnAdd);
  }
  dOut->cd();
  lOut->Write(lOut->GetName(),TObject::kSingleKey);
  fOut->Close();
}

Bool_t AliAnalysisHelperJetTasks::PythiaInfoFromFile(const char* currFile,Float_t &fXsec,Float_t &fTrials){
  //
  // get the cross section and the trails either from pyxsec.root or from pysec_hists.root
  // This is to called in Notify and should provide the path to the AOD/ESD file

  TString file(currFile);  
  fXsec = 0;
  fTrials = 1;

  if(file.Contains("root_archive.zip#")){
    Ssiz_t pos1 = file.Index("root_archive",12,TString::kExact);
    Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
    file.Replace(pos+1,20,"");
  }
  else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  Printf("%s",file.Data());
  
 
   

  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); // problem that we cannot really test the existance of a file in a archive so we have to lvie with open error message from root
  if(!fxsec){
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if(!fxsec){
	// not a severe condition but inciate that we have no information
      return kFALSE;
    }
    else{
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0); 
      if(!key){
	fxsec->Close();
	return kFALSE;
      }
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if(!list){
	fxsec->Close();
	return kFALSE;
      }
      fXsec = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
      fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } // no tree pyxsec.root
  else {
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
    if(!xtree){
      fxsec->Close();
      return kFALSE;
    }
    UInt_t   ntrials  = 0;
    Double_t  xsection  = 0;
    xtree->SetBranchAddress("xsection",&xsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
    fTrials = ntrials;
    fXsec = xsection;
    fxsec->Close();
  }
  return kTRUE;
}

Bool_t AliAnalysisHelperJetTasks::PrintDirectorySize(const char* currFile,Int_t iDetail){

  //
  // Print the size on disk and in memory occuppied by a directory 
  //

  TFile *fDummy = 0;
  if(iDetail>=0)fDummy = TFile::Open("/tmp/dummy.root","RECREATE");

  TFile *fIn = TFile::Open(currFile);
  if(!fIn){
    // not a severe condition but inciate that we have no information
    return kFALSE;
  }
  // find the tlists we want to be independtent of the name so use the Tkey
  TList* keyList = fIn->GetListOfKeys();
  Float_t memorySize = 0;
  Float_t diskSize = 0;

  for(int i = 0;i < keyList->GetEntries();i++){
    TKey* ikey = (TKey*)keyList->At(i); 
    
    //    TList *list = dynamic_cast<TList*>(key->ReadObj());
    //    TNamed *name = dynamic_cast<TNamed*>(ikey->ReadObj());
    TDirectory *dir =  dynamic_cast<TDirectory*>(ikey->ReadObj());
    



    if(dir){
      Printf("%03d    : %60s %8d %8d ",i,dir->GetName(),ikey->GetObjlen(),ikey->GetNbytes());
      TList * dirKeyList = dir->GetListOfKeys();
      for(int j = 0;j<dirKeyList->GetEntries();j++){
	TKey* jkey = (TKey*)dirKeyList->At(j); 
	TList *list =  dynamic_cast<TList*>(jkey->ReadObj());

	memorySize += (Float_t)jkey->GetObjlen()/1024./1024.;
	diskSize +=  (Float_t)jkey->GetNbytes()/1024./1024.;
	if(list){
	  Printf("%03d/%03d: %60s %5.2f MB %5.2f MB",i,j,list->GetName(),(Float_t)jkey->GetObjlen()/1024./1024.,(Float_t)jkey->GetNbytes()/1024./1024.);
	  if(iDetail==i){
	    for(int il = 0;il<list->GetEntries();il++){
	      TObject *ob = list->At(il);
	      if(fDummy){
		fDummy->cd();
		Int_t nBytesWrite = ob->Write();
		Printf("%03d/%03d/%03d: %60s  %5.2f kB",i,j,il,ob->GetName(),(Float_t)nBytesWrite/1024.);
	      }
	    }
	  }
	}
	else{
	  Printf("%03d/%03d: %60s %5.2f MB %5.2f MB",i,j,jkey->GetName(),(Float_t)jkey->GetObjlen()/1024./1024.,(Float_t)jkey->GetNbytes()/1024./1024.);
	}
      }
    }
  }
  Printf("Total %5.2f MB %5.2f MB",memorySize,diskSize);
  return kTRUE;
}


Bool_t  AliAnalysisHelperJetTasks::Selected(Bool_t bSet,Bool_t bNew){
  // 
  // Static helper task, (re)set event by event
  //


  static Bool_t bSelected = kTRUE; // if service task is not run we acccpet all
  if(bSet){
    bSelected = bNew;
  }
  return bSelected;
}

Bool_t  AliAnalysisHelperJetTasks::IsCosmic(){
  return ((SelectInfo()&kIsCosmic)==kIsCosmic);
}

Bool_t  AliAnalysisHelperJetTasks::IsPileUp(){
  return ((SelectInfo()&kIsPileUp)==kIsPileUp);
}


Bool_t  AliAnalysisHelperJetTasks::TestSelectInfo(UInt_t iMask){
  return ((SelectInfo()&iMask)==iMask);
}


Bool_t  AliAnalysisHelperJetTasks::TestEventClass(Int_t iMask){
  if(iMask==0)return kTRUE;
  return (EventClass()==iMask);
}


UInt_t  AliAnalysisHelperJetTasks::SelectInfo(Bool_t bSet,UInt_t iNew){
  // 
  // set event by event the slection info
  // 

  static UInt_t iSelectInfo = 0; //
  if(bSet){
    iSelectInfo = iNew;
  }
  return iSelectInfo;
}


Int_t  AliAnalysisHelperJetTasks::EventClass(Bool_t bSet,Int_t iNew){
  //
  // gloab storage/access to event class
  //

  static Int_t iEventClass = 0; //
  if(bSet){
    iEventClass = iNew;
  }
  return iEventClass;
}




//___________________________________________________________________________________________________________

Bool_t AliAnalysisHelperJetTasks::GetEventShapes(TVector3 &n01,const  TVector3 * pTrack, Int_t nTracks, Double_t * eventShapes)
{       
  // ***
  // Event shape calculation
  // sona.pochybova@cern.ch

  const Int_t kTracks = 1000;
  if(nTracks>kTracks)return kFALSE;

  //variables for thrust calculation
  TVector3 pTrackPerp[kTracks];
  Double_t psum2 = 0;

  TVector3 psum;
  TVector3 psum02;
  TVector3 psum03;

  Double_t psum1 = 0;
  Double_t psum102 = 0;
  Double_t psum103 = 0;

  Double_t thrust[kTracks] = {0.0};
  Double_t th = -3;
  Double_t thrust02[kTracks] = {0.0};
  Double_t th02 = -4;
  Double_t thrust03[kTracks] = {0.0};
  Double_t th03 = -5;

  //Sphericity calculation variables
  TMatrixDSym m(3);
  Double_t s00 = 0;
  Double_t s01 = 0;
  Double_t s02 = 0;
  
  Double_t s10 = 0;
  Double_t s11 = 0;
  Double_t s12 = 0;
  
  Double_t s20 = 0;
  Double_t s21 = 0;
  Double_t s22 = 0;
  
  Double_t ptot = 0;
  
  Double_t c = -10;

//
//loop for thrust calculation  
//
    
  for(Int_t i = 0; i < nTracks; i++)
    {
      pTrackPerp[i].SetXYZ(pTrack[i].X(), pTrack[i].Y(), 0);
      psum2 += pTrackPerp[i].Mag();
    }

  //additional starting axis    
  TVector3 n02;
  n02 = pTrack[1].Unit();
  n02.SetZ(0.);   
  TVector3 n03;
  n03 = pTrack[2].Unit();
  n03.SetZ(0.);

  //switches for calculating thrust for different starting points
  Int_t switch1 = 1;
  Int_t switch2 = 1;
  Int_t switch3 = 1;

  //indexes for iteration of different starting points
  Int_t l1 = 0;
  Int_t l2 = 0;
  Int_t l3 = 0;

  //maximal number of iterations
  //  Int_t nMaxIter = 100;
 
  for(Int_t k = 0; k < nTracks; k++)
    {  
      
      if(switch1 == 1){
	psum.SetXYZ(0., 0., 0.);
	psum1 = 0;
	for(Int_t i = 0; i < nTracks; i++)
	  {
	    psum1 += (TMath::Abs(n01.Dot(pTrackPerp[i])));
	    if (n01.Dot(pTrackPerp[i]) > 0) psum += pTrackPerp[i];
	    if (n01.Dot(pTrackPerp[i]) < 0) psum -= pTrackPerp[i];
	  }
	thrust[l1] = psum1/psum2;
      }

      if(switch2 == 1){
	psum02.SetXYZ(0., 0., 0.);
	psum102 = 0;
	for(Int_t i = 0; i < nTracks; i++)
	  {
	    psum102 += (TMath::Abs(n02.Dot(pTrackPerp[i])));
	    if (n02.Dot(pTrackPerp[i]) > 0) psum02 += pTrackPerp[i];
	    if (n02.Dot(pTrackPerp[i]) < 0) psum02 -= pTrackPerp[i];
	  }
	thrust02[l2] = psum102/psum2;
      }
      
      if(switch3 == 1){
	psum03.SetXYZ(0., 0., 0.);
	psum103 = 0;
	for(Int_t i = 0; i < nTracks; i++)
	  {
	    psum103 += (TMath::Abs(n03.Dot(pTrackPerp[i])));
	    if (n03.Dot(pTrackPerp[i]) > 0) psum03 += pTrackPerp[i];
	    if (n03.Dot(pTrackPerp[i]) < 0) psum03 -= pTrackPerp[i];
	  }
	thrust03[l3] = psum103/psum2;
      }

      //check whether thrust value converged    
      if(TMath::Abs(th-thrust[l1]) < 10e-7){
      	switch1 = 0;
      }
      
      if(TMath::Abs(th02-thrust02[l2]) < 10e-7){
	switch2 = 0;
      }

      if(TMath::Abs(th03-thrust03[l3]) < 10e-7){
	switch3 = 0;
      }

      //if it didn't, continue with the calculation
      if(switch1 == 1){
	th = thrust[l1]; 
	n01 = psum.Unit();
	l1++;
      }

      if(switch2 == 1){
	th02 = thrust02[l2];
	n02 = psum02.Unit();
	l2++;
      }

      if(switch3 == 1){
	th03 = thrust03[l3];
	n03 = psum03.Unit();
	l3++;
      }

      //if thrust values for all starting direction converged check if to the same value
      if(switch2 == 0 && switch1 == 0 && switch3 == 0){
	if(TMath::Abs(th-th02) < 10e-7 && TMath::Abs(th-th03) < 10e-7 && TMath::Abs(th02-th03) < 10e-7){
	  eventShapes[0] = th;
	  AliInfoGeneral(Form(" %s:%d",(char*)__FILE__,__LINE__),Form("===== THRUST VALUE FOUND AT %d :: %f\n", k, th));
	  break;
	}
	//if they did not, reset switches
	else{
	  switch1 = 1;
	  //	  th = -1.;
	  switch2 = 1;
	  //	  th02 = -2.;
	  switch3 = 1;
	  //	  th03 = -4.;
	}
      }

      //      Printf("========== %d +++ th :: %f=============\n", l1, th);
      //      Printf("========== %d +++ th2 :: %f=============\n", l2, th02);
      //      Printf("========== %d +++ th3 :: %f=============\n", l3, th03);
      
    }

  //if no common limitng value was found, take the maximum and take the corresponding thrust axis
  if(switch1 == 1 && switch2 == 1 && switch3 == 1){
    eventShapes[0] = TMath::Max(thrust[l1-1], thrust02[l2-1]);
    eventShapes[0] = TMath::Max(eventShapes[0], thrust03[l3-1]);
    if(TMath::Abs(eventShapes[0]-thrust[l1-1]) < 10e-7)
      n01 = n01;
    if(TMath::Abs(eventShapes[0]-thrust02[l2-1]) < 10e-7)
      n01 = n02;
    if(TMath::Abs(eventShapes[0]-thrust03[l3-1]) < 10e-7)
      n01 = n03;
    Printf("NO LIMITING VALUE FOUND :: MAXIMUM = %f\n", eventShapes[0]);
  }

//
//other event shapes variables
//
  for(Int_t j = 0; j < nTracks; j++)
    {  
      s00 = s00 + (pTrack[j].Px()*pTrack[j].Px())/pTrack[j].Mag();
      s01 = s01 + (pTrack[j].Px()*pTrack[j].Py())/pTrack[j].Mag();
      s02 = s02 + (pTrack[j].Px()*pTrack[j].Pz())/pTrack[j].Mag();
      
      s10 = s10 + (pTrack[j].Py()*pTrack[j].Px())/pTrack[j].Mag();
      s11 = s11 + (pTrack[j].Py()*pTrack[j].Py())/pTrack[j].Mag();
      s12 = s12 + (pTrack[j].Py()*pTrack[j].Pz())/pTrack[j].Mag();
      
      s20 = s20 + (pTrack[j].Pz()*pTrack[j].Px())/pTrack[j].Mag();
      s21 = s21 + (pTrack[j].Pz()*pTrack[j].Py())/pTrack[j].Mag();
      s22 = s22 + (pTrack[j].Pz()*pTrack[j].Pz())/pTrack[j].Mag();
      
      ptot += pTrack[j].Mag();
    }

  if(ptot > 0.)
    {
      m(0,0) = s00/ptot;
      m(0,1) = s01/ptot;
      m(0,2) = s02/ptot;

      m(1,0) = s10/ptot;
      m(1,1) = s11/ptot;
      m(1,2) = s12/ptot;

      m(2,0) = s20/ptot;
      m(2,1) = s21/ptot;
      m(2,2) = s22/ptot;

      TMatrixDSymEigen eigen(m);
      TVectorD eigenVal = eigen.GetEigenValues();

      Double_t sphericity = (3/2)*(eigenVal(2)+eigenVal(1));
      eventShapes[1] = sphericity;

      Double_t aplanarity = (3/2)*(eigenVal(2));
      eventShapes[2] = aplanarity;

      c = 3*(eigenVal(0)*eigenVal(1)+eigenVal(0)*eigenVal(2)+eigenVal(1)*eigenVal(2));
      eventShapes[3] = c;
    }
  return kTRUE;
}
  


 //__________________________________________________________________________________________________________________________
// Trigger Decisions copid from PWG0/AliTriggerAnalysis


Bool_t AliAnalysisHelperJetTasks::IsTriggerFired(const AliVEvent* aEv, Trigger trigger)
{
  // checks if an event has been triggered
  // no usage of ofline trigger here yet
  
  // here we do a dirty hack to take also into account the
  // missing trigger bits and Bunch crossing paatern for real data 


  if(aEv->InheritsFrom("AliESDEvent")){
    const AliESDEvent *esd = (AliESDEvent*)aEv;
    switch (trigger)
      {
      case kAcceptAll:
	{
	  return kTRUE;
	  break;
	}
      case kMB1:
	{
	  if(esd->GetFiredTriggerClasses().Contains("CINT1B-"))return kTRUE;
	  // does the same but without or'ed V0s
	  if(esd->GetFiredTriggerClasses().Contains("CSMBB"))return kTRUE;  
	  if(esd->GetFiredTriggerClasses().Contains("CINT6B-"))return kTRUE; 
	  // this is for simulated data
	  if(esd->GetFiredTriggerClasses().Contains("MB1"))return kTRUE;   
	  break;
	}
      case kMB2:
	{
	  if(esd->GetFiredTriggerClasses().Contains("MB2"))return kTRUE;   
	  break;
	}
      case kMB3:
	{
	  if(esd->GetFiredTriggerClasses().Contains("MB3"))return kTRUE;   
	  break;
	}
      case kSPDGFO:
	{
	  if(esd->GetFiredTriggerClasses().Contains("CSMBB"))return kTRUE;
	  if(esd->GetFiredTriggerClasses().Contains("CINT6B-"))return kTRUE; 
	  if(esd->GetFiredTriggerClasses().Contains("GFO"))return kTRUE;
	  break;
	}
      default:
	{
	  Printf("IsEventTriggered: ERROR: Trigger type %d not implemented in this method", (Int_t) trigger);
	  break;
	}
      }
  }
  else if(aEv->InheritsFrom("AliAODEvent")){
    const AliAODEvent *aod = (AliAODEvent*)aEv;
    switch (trigger)
      {
      case kAcceptAll:
	{
	  return kTRUE;
	  break;
	}
      case kMB1:
	{
	  if(aod->GetFiredTriggerClasses().Contains("CINT1B"))return kTRUE;
	  // does the same but without or'ed V0s
	  if(aod->GetFiredTriggerClasses().Contains("CSMBB"))return kTRUE;   
	  // for sim data
	  if(aod->GetFiredTriggerClasses().Contains("MB1"))return kTRUE;   
	  break;
	}
      case kMB2:
	{
	  if(aod->GetFiredTriggerClasses().Contains("MB2"))return kTRUE;   
	  break;
	}
      case kMB3:
	{
	  if(aod->GetFiredTriggerClasses().Contains("MB3"))return kTRUE;   
	  break;
	}
      case kSPDGFO:
	{
	  if(aod->GetFiredTriggerClasses().Contains("CSMBB"))return kTRUE;	  
	  if(aod->GetFiredTriggerClasses().Contains("GFO"))return kTRUE;   
	  break;
	}
      default:
	{
	  Printf("IsEventTriggered: ERROR: Trigger type %d not implemented in this method", (Int_t) trigger);
	  break;
	}
      }
  }
    return kFALSE;
}


 AliAnalysisHelperJetTasks::MCProcessType  AliAnalysisHelperJetTasks::GetPythiaEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug) {
   //
   // fetch the process type from the mc header
   //


  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(aHeader);

  if (!pythiaGenHeader) {
    //    printf(" AliAnalysisHelperJetTasks::GetProcessType : Unknown gen Header type). \n");
    return kInvalidProcess;
  }


  Int_t pythiaType = pythiaGenHeader->ProcessType();
  fgLastProcessType = pythiaType;
  MCProcessType globalType = kInvalidProcess;  


  if (adebug) {
    printf(" AliAnalysisHelperJetTasks::GetProcessType : Pythia process type found: %d \n",pythiaType);
  }


  if(pythiaType==92||pythiaType==93){
    globalType = kSD;
  }
  else if(pythiaType==94){
    globalType = kDD;
  }
  //else if(pythiaType != 91){ // also exclude elastic to be sure... CKB??
  else {
    globalType = kND;
  }
  return globalType;
}


 AliAnalysisHelperJetTasks::MCProcessType  AliAnalysisHelperJetTasks::GetDPMjetEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug) {
  //
  // get the process type of the event.
  //

  // can only read pythia headers, either directly or from cocktalil header
  AliGenDPMjetEventHeader* dpmJetGenHeader = dynamic_cast<AliGenDPMjetEventHeader*>(aHeader);

  if (!dpmJetGenHeader) {
    printf(" AliAnalysisHelperJetTasks::GetDPMjetProcessType : Unknown header type (not DPMjet or). \n");
    return kInvalidProcess;
  }

  Int_t dpmJetType = dpmJetGenHeader->ProcessType();
  fgLastProcessType = dpmJetType;
  MCProcessType globalType = kInvalidProcess;  


  if (adebug) {
    printf(" AliAnalysisHelperJetTasks::GetDPMJetProcessType : DPMJet process type found: %d \n",dpmJetType);
  }


  if (dpmJetType == 1 || dpmJetType == 4) { // explicitly inelastic plus central diffraction
    globalType = kND;
  }  
  else if (dpmJetType==5 || dpmJetType==6) {
    globalType = kSD;
  }
  else if (dpmJetType==7) {
    globalType = kDD;
  }
  return globalType;
}


Int_t AliAnalysisHelperJetTasks::GetPhiBin(Double_t phi,Int_t fNRPBins)
{
    Int_t phibin=-1;
    if(!(TMath::Abs(phi)<=2*TMath::Pi()))return -1;
    Double_t phiwrtrp=TMath::ACos(TMath::Abs(TMath::Cos(phi)));
    phibin=Int_t(fNRPBins*phiwrtrp/(0.5*TMath::Pi()));
    return phibin;
}

Double_t  AliAnalysisHelperJetTasks::ReactionPlane(Bool_t bSet,Double_t fNew){
  // 
  // Static helper task, (re)set event by event
  //


  static Double_t fRP = 0; // if service task is not run we acccpet all
  if(bSet){
    fRP = fNew;
  }
  return fRP;
}
