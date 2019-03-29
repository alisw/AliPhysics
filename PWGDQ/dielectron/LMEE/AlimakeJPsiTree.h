#ifndef AlimakeJPsiTree_H
#define AlimakeJPsiTree_H

//########################################################################
//#                                                                      #
//#   Class to produce charm cocktail (analysis task)                    #
//#                                                                      #
//#  Authors:                                                            #
//#   Raphaelle Bailhache, Uni Frankfurt / Raphaelle.Bailhache@cern.ch   #
//#   Carsten Klein, Uni Frankfurt / Carsten.Klein@cern.ch               #
//#   Jerome Jung, Uni Frankfurt / s8286523@uni-frankfurt.de             #
//#   Sebastian Scheid, Uni Frankfurt / s.scheid@cern.ch                 #
//#                                                                      #
//########################################################################

#include "AliAnalysisTaskSE.h"
// local files
class TList;
class TTree;
class AliAODMCParticle;
class AliMCEvent;
class AliInputEventHandler;

class AlimakeJPsiTree : public AliAnalysisTaskSE {
  
public:
  AlimakeJPsiTree(); ///< default constructor probably needed for AnalysisManager or such...
  AlimakeJPsiTree(const Char_t* name); ///< named constructor which also creates input and output objects.
  virtual ~AlimakeJPsiTree();
  
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit(){Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  
  
private:

  
protected:
  AliMCEvent*              fMcEvent;    //! MC event
  AliInputEventHandler*    fMcHandler;  //! MC EventHandler
  // Tree
  TTree *ftree;                         //
  AliAODMCParticle *fmother;                   //
  AliAODMCParticle *fdaughter1;                //
  AliAODMCParticle *fdaughter2;                //

  TList     *fOutputList; 

  AlimakeJPsiTree(const AlimakeJPsiTree &c); // not implemented
  AlimakeJPsiTree& operator= (const AlimakeJPsiTree &c); // not implemented
  
  ClassDef(AlimakeJPsiTree,1)
};

#endif

