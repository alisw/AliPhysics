#include "AliITSAlignMille2ConstrArray.h"
#include "AliITSAlignMille2Module.h"

//========================================================================================================
ClassImp(AliITSAlignMille2ConstrArray)

//________________________________________________________________________________________________________
AliITSAlignMille2ConstrArray::AliITSAlignMille2ConstrArray() :
AliITSAlignMille2Constraint(),
fModuleIDs(0),
fModulePatt(0),
fCoeffs(0),
fError(0)
{}

//________________________________________________________________________________________________________
AliITSAlignMille2ConstrArray::AliITSAlignMille2ConstrArray(const Char_t* name,Double_t *parcf,Int_t npar,Double_t val,Double_t err) :
  AliITSAlignMille2Constraint(name,kTypeGaussian,-100,val,0),
  fModuleIDs(0),
  fModulePatt(0),
  fCoeffs(npar),
  fError(err)
{
  // create module
  for (int i=0;i<npar;i++) {fCoeffs[i] = parcf[i]; if (parcf[i]!=0) SetBit(0x1<<i);}
}

//________________________________________________________________________________________________________
AliITSAlignMille2ConstrArray::AliITSAlignMille2ConstrArray(const AliITSAlignMille2ConstrArray& src) :
  AliITSAlignMille2Constraint(src),
  fModuleIDs(0),
  fModulePatt(0),
  fCoeffs(0),
  fError(0)
{/* DUMMY */} 

//________________________________________________________________________________________________________
void AliITSAlignMille2ConstrArray::AddModule(AliITSAlignMille2Module* mod, Bool_t needGeom)
{
  // add module to constraint
  int nmd = GetNModules();
  // check if its already not there
  for (int im=nmd;im--;) {if (mod->GetUniqueID() == (UInt_t)fModuleIDs[im]) return; }
  fModuleIDs.Set(nmd+1);
  fModulePatt.Set(nmd+1);
  fModuleIDs[nmd] = mod->GetUniqueID();
  if (needGeom) { // this is geometrical constraint
    double jacobian[AliITSAlignMille2Module::kMaxParGeom][AliITSAlignMille2Module::kMaxParGeom];
    if (mod->GeomParamsGlobal()) mod->CalcDerivLocGlo(&jacobian[0][0]);
    //
    Short_t patt =  GetPattern();
    if (mod->GeomParamsGlobal()) {
      // the constraint is defined in the module's local frame. If the alignment of geom params is
      // done in the global frame, we need to set the real parameter involved
      for (int i=AliITSAlignMille2Module::kMaxParGeom;i--;) patt &= ~BIT(i);  // reset the geometry parameters
      for (int i=0;i<AliITSAlignMille2Module::kMaxParGeom;i++) {
	if (!IncludesParam(i)) continue;
	for (int j=0;j<AliITSAlignMille2Module::kMaxParGeom;j++) if (jacobian[i][j]!=0) patt |= BIT(j);
      }
    }
    fModulePatt[nmd] = patt;
  }
}

//________________________________________________________________________________________________________
Bool_t AliITSAlignMille2ConstrArray::IncludesModule(Int_t id) const
{
  // is this module mentioned in the list?
  int nmd = GetNModules();
  for (int i=nmd;i--;) if (fModuleIDs[i]==id) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________________________________________
Bool_t AliITSAlignMille2ConstrArray::IncludesModPar(Int_t id,Int_t par) const
{
  // is this module/parameter mentioned in the list?
  int nmd = GetNModules();
  for (int i=nmd;i--;) {
    if (fModuleIDs[i]!=id) continue;
    if (fModulePatt[i] & (0x1<<par)) return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________________________________________
Bool_t AliITSAlignMille2ConstrArray::IncludesModPar(const AliITSAlignMille2Module* mod, Int_t par) const
{
  // is this module/parameter mentioned in the list?
  return IncludesModPar(mod->GetUniqueID(), par);
}


//________________________________________________________________________________________________________
void AliITSAlignMille2ConstrArray::Print(Option_t* ) const
{
  // print data
  printf("#%3d Constraint %s of type %d | Value=%+e Error=%+e\n",GetConstraintID(),GetName(),GetType(),GetValue(),GetError());
  printf("Weights on params: "); for (int i=0;i<GetNCoeffs();i++) printf("%+.3e ",GetCoeff(i)); 
  printf("\nModules involved: \n");
  int nmd = GetNModules();
  for (int i=0;i<nmd;i++) {printf("%4d |",GetModuleID(i)); if ( ((i+1)%14)==0 ) printf("\n");}
  if ( nmd%14 ) printf("\n");
  //
}

