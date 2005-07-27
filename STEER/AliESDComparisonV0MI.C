/*

.L $ALICE_ROOT/STEER/AliGenInfo.C+
.L $ALICE_ROOT/STEER/AliESDComparisonMI.C+
.L AliESDComparisonV0MI.C
.x TDR_style.C

AliKalmanTrack::SetConvConst(1);
AliESDComparisonDraw compv0;  
AliESDComparisonDraw comp;  


TFile fv0("cmpESDTracks.root");
TTree * treev0 = (TTree*)fv0.Get("ESDcmpV0");
TTree * tree = (TTree*)fv0.Get("ESDcmpTracks");
compv0.fTree = treev0;
comp.fTree = tree;
MakeV0Aliases(&compv0);
MakeAliases(comp);

//
//
TChain * chain1 = new TChain("ESDcmpV0","ESDcmpV0");
chain1->AddFile("k0merge/cmpESDTracks.root");
chain1->AddFile("k0merge1000/cmpESDTracks.root");
chain1->AddFile("strangemerge/newv1/cmpESDTracks.root");
compv0.fTree = chain1;
MakeV0Aliases(&compv0);
//
//
TChain * chain3 = new TChain("ESDcmpV0","ESDcmpV0");
//
chain3->AddFile("750801/cmpESDTracks.root");
chain3->AddFile("751760/cmpESDTracks.root");
chain3->AddFile("751764/cmpESDTracks.root");
chain3->AddFile("751765/cmpESDTracks.root");
chain3->AddFile("751766/cmpESDTracks.root");
chain3->AddFile("751767/cmpESDTracks.root");
chain3->AddFile("751768/cmpESDTracks.root");
chain3->AddFile("751769/cmpESDTracks.root");
chain3->AddFile("751772/cmpESDTracks.root");
//
chain3->AddFile("751774/cmpESDTracks.root");
chain3->AddFile("751775/cmpESDTracks.root");
chain3->AddFile("751776/cmpESDTracks.root");
chain3->AddFile("751778/cmpESDTracks.root");
chain3->AddFile("751779/cmpESDTracks.root");
chain3->AddFile("751781/cmpESDTracks.root");
chain3->AddFile("751782/cmpESDTracks.root");
chain3->AddFile("751784/cmpESDTracks.root");
//
chain3->AddFile("751788/cmpESDTracks.root");
chain3->AddFile("751789/cmpESDTracks.root");
chain3->AddFile("751790/cmpESDTracks.root");
chain3->AddFile("751791/cmpESDTracks.root");
chain3->AddFile("751792/cmpESDTracks.root");
chain3->AddFile("751793/cmpESDTracks.root");
chain3->AddFile("751794/cmpESDTracks.root");
chain3->AddFile("751796/cmpESDTracks.root");

compv0.fTree = chain3;
MakeV0Aliases(&compv0);

TChain * chain3 = new TChain("ESDcmpTracks","ESDcmpTracks");
chain3->AddFile("751790/cmpESDTracks.root");
chain3->AddFile("751791/cmpESDTracks.root");
chain3->AddFile("751792/cmpESDTracks.root");
chain3->AddFile("751793/cmpESDTracks.root");
chain3->AddFile("751794/cmpESDTracks.root");
chain3->AddFile("751796/cmpESDTracks.root");
comp.fTree = chain3;
MakeAliases(comp);
//
TChain * chain2 = new TChain("ESDcmpV0","ESDcmpV0");
chain2->AddFile("run50_100/newv6/cmpESDTracks.root");
chain2->AddFile("run20_100/newv6/cmpESDTracks.root");
compv0.fTree = chain2;
MakeV0Aliases(&compv0);

//
//
TChain * chain3 = new TChain("ESDcmpV0","ESDcmpV0");
chain3->AddFile("~/data/HEAD0205/cent1_0/newv7cmpESDTracks.root");

*/

TH1F hrel("hrel","hrel",100,-1,1);
TCut cbase("cbase","RC.fV0Status>-1");
TCut cr("cr","RC.fV0Status==3");
TCut cp("cp","MC.fMotherP.fVx<0.2&&MC.fMotherP.fVy<0.2");
TCut clambda("clambda","MC.fMotherP.fPdgCode==3122");
TCut calambda("calambda","MC.fMotherP.fPdgCode==-3122");
TCut clam("clam","abs(MC.fMotherP.fPdgCode)==3122");
TCut cgamma("cgamma","abs(MC.fMotherP.fPdgCode)==22");
TCut ck0s("ck0s","MC.fMotherP.fPdgCode==313");
TCut ck0b("ck0b","MC.fMotherP.fPdgCode==-313");
TCut cks0("cks0","MC.fMotherP.fPdgCode==310");
TCut crho0("crho0","MC.fMotherP.fPdgCode==113");
TCut cpi("cpi","abs(MC.fMotherP.fPdgCode)==211");
TCut cpa("cpa","MC.fPointAngle>0.999");
TCut cint("cint","abs(MC.fMotherP.fPdgCode)==3122||abs(MC.fMotherP.fPdgCode)==22||MC.fMotherP.fPdgCode==310"+cpa);
TH1F hpull("hpull","hpull",100,-5,5); 
TF1 fim("fim","[2]*exp(-(x-[0])**2/(2*[1]**2))+abs([3])+[4]*(x-0.495)",0,1);


void MakeV0Aliases(AliESDComparisonDraw * compv0)
{
  //
  // Make base, resol and cut aliases
  //
  MakeAliasBase(compv0);
  MakeResolAlias(compv0);
  MakeCutsAliases(compv0);
  MakeCutsAliasesPID(compv0);
}

void MakeAliasBase(AliESDComparisonDraw * compv0){
  //
  //  aliases based on MC data
  //  and on Reconstructed data
  //
  //  comp->fTree->SetAlias("locphi","(100*(180+180*atan2(MC.fParticle.fPy,MC.fParticle.fPx)/pi)%1000.)*0.01");
  compv0->fTree->SetAlias("minlab","min(abs(RC.fV0rec.fLab[0]),abs(RC.fV0rec.fLab[1]))");    
  compv0->fTree->SetAlias("maxlab","max(abs(RC.fV0rec.fLab[0]),abs(RC.fV0rec.fLab[1]))");    
  compv0->fTree->SetAlias("MCE","sqrt((MC.fMCm.fParticle.fPx+MC.fMCd.fParticle.fPx)**2+(MC.fMCm.fParticle.fPy+MC.fMCd.fParticle.fPy)**2+(MC.fMCm.fParticle.fPz+MC.fMCd.fParticle.fPz)**2)");  
  compv0->fTree->SetAlias("Ptpc","sqrt((RC.fV0tpc.fPM[0]+RC.fV0tpc.fPP[0])**2+(RC.fV0tpc.fPM[1]+RC.fV0tpc.fPP[1])**2+(RC.fV0tpc.fPM[2]+RC.fV0tpc.fPP[2])**2)");  
  compv0->fTree->SetAlias("Pits","sqrt((RC.fV0its.fPM[0]+RC.fV0its.fPP[0])**2+(RC.fV0its.fPM[1]+RC.fV0its.fPP[1])**2+(RC.fV0its.fPM[2]+RC.fV0its.fPP[2])**2)");
  compv0->fTree->SetAlias("Prec","sqrt((RC.fV0rec.fPM[0]+RC.fV0rec.fPP[0])**2+(RC.fV0rec.fPM[1]+RC.fV0rec.fPP[1])**2+(RC.fV0rec.fPM[2]+RC.fV0rec.fPP[2])**2)");

  compv0->fTree->SetAlias("PTtpc","sqrt((RC.fV0tpc.fPM[0]+RC.fV0tpc.fPP[0])**2+(RC.fV0tpc.fPM[1]+RC.fV0tpc.fPP[1])**2)");  
  compv0->fTree->SetAlias("PTits","sqrt((RC.fV0its.fPM[0]+RC.fV0its.fPP[0])**2+(RC.fV0its.fPM[1]+RC.fV0its.fPP[1])**2)");
  compv0->fTree->SetAlias("PTrec","sqrt((RC.fV0rec.fPM[0]+RC.fV0rec.fPP[0])**2+(RC.fV0rec.fPM[1]+RC.fV0rec.fPP[1])**2)");
  compv0->fTree->SetAlias("PhiRec","atan2(RC.fV0rec.fPP[1]+RC.fV0rec.fPM[1],RC.fV0rec.fPP[0]+RC.fV0rec.fPM[0])");
  compv0->fTree->SetAlias("PhiTpc","atan2(RC.fV0tpc.fPP[1]+RC.fV0tpc.fPM[1],RC.fV0tpc.fPP[0]+RC.fV0tpc.fPM[0])");
  compv0->fTree->SetAlias("PhiIts","atan2(RC.fV0its.fPP[1]+RC.fV0its.fPM[1],RC.fV0its.fPP[0]+RC.fV0its.fPM[0])");
  compv0->fTree->SetAlias("PhiMC","atan2(MC.fMotherP.fPy,MC.fMotherP.fPx)");

  compv0->fTree->SetAlias("minpt","min(abs(MC.fMCm.fParticle.Pt()),abs(MC.fMCd.fParticle.Pt()))");
  compv0->fTree->SetAlias("minpt2","min(abs(MC.fMCm.fTrackRef.Pt()),abs(MC.fMCd.fTrackRef.Pt()))");
  compv0->fTree->SetAlias("minrows","min(MC.fMCm.fRowsWithDigits,MC.fMCd.fRowsWithDigits)");
  compv0->fTree->SetAlias("maxtan","tan(max(abs(MC.fMCm.fParticle.Theta()-1.578),abs(MC.fMCd.fParticle.Theta()-1.578)))");
  compv0->fTree->SetAlias("maxtanr","max(abs(RC.fV0rec.fParamP.fParam[3]),abs(RC.fV0rec.fParamM.fParam[3]))");
  
  compv0->fTree->SetAlias("findable","RC.fRecStatus>=0&&MC.fMCR>0.8&&MC.fMCR<120&&minrows>100&&minpt>0.2&&maxtan<1.0&&sqrt(MC.fMotherP.fVx**2+MC.fMotherP.fVy**2)<0.00001&&MC.fPointAngle>0.9999");
  
  compv0->fTree->SetAlias("findable2","RC.fRecStatus>=0&&MC.fMCR>0.6&&MC.fMCR<120&&MC.fPointAngle>0.9999&&min(MC.fMCm.fRowsWithDigits,MC.fMCd.fRowsWithDigits)>120&&max(abs(MC.fMCm.fTrackRef.fZ),abs(MC.fMCd.fTrackRef.fZ))<90&&min(abs(MC.fMCm.fParticle.Pt()),abs(MC.fMCd.fParticle.Pt()))>0.2");
  
  compv0->fTree->SetAlias("EMC2","(MC.fMCd.fParticle.fE+MC.fMCm.fParticle.fE)*(MC.fMCd.fParticle.fE+MC.fMCm.fParticle.fE)");
  compv0->fTree->SetAlias("PMC2","(MC.fMCd.fParticle.fPx+MC.fMCm.fParticle.fPx)**2+(MC.fMCd.fParticle.fPy+MC.fMCm.fParticle.fPy)**2+(MC.fMCd.fParticle.fPz+MC.fMCm.fParticle.fPz)**2");
}  

////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////

void MakeResolAlias(AliESDComparisonDraw * compv0){
  //
  // Resolution estimates 
  // SigmaAP0 and SigmaD0 - uses Covariance matrix
  // SigmaAP2 and SigmaD2 - "Restricted Sigma"  in ( 0.5  to 1.5) window from effective sigma parameterization
  //                       
  compv0->fTree->SetAlias("CausalityB","(1-RC.fV0rec.fCausality[0])*(1-RC.fV0rec.fCausality[1])");
  compv0->fTree->SetAlias("CausalityA","sqrt((min(RC.fV0rec.fCausality[2],0.7))*(min(RC.fV0rec.fCausality[3],0.7)))");
  compv0->fTree->SetAlias("p12","sqrt(RC.fV0rec.fParamP.fParam[4]**2+RC.fV0rec.fParamM.fParam[4]**2)");
  compv0->fTree->SetAlias("ptminrec","1./max(abs(RC.fV0rec.fParamP.fParam[4]),abs(RC.fV0rec.fParamM.fParam[4]))");
  compv0->fTree->SetAlias("P0","(exp(-RC.fV0rec.fDistNorm)+0.1)*(RC.fV0rec.fDist2<0.5)*(RC.fV0rec.fDistNorm<5)");
  //
  // resolution estimates in Y and Z
  //
  compv0->fTree->SetAlias("SigmaY","sqrt(RC.fV0rec.fParamP.fCovar[0]+RC.fV0rec.fParamM.fCovar[0]+(RC.fV0rec.fParamP.fCovar[5])*(RC.fV0rec.fParamP.fX-RC.fV0rec.fRr)**2+(RC.fV0rec.fParamM.fCovar[5])*(RC.fV0rec.fParamM.fX-RC.fV0rec.fRr)**2)");  
  compv0->fTree->SetAlias("SigmaZ","sqrt(RC.fV0rec.fParamP.fCovar[2]+RC.fV0rec.fParamM.fCovar[2]+(RC.fV0rec.fParamP.fCovar[9])*(RC.fV0rec.fParamP.fX-RC.fV0rec.fRr)**2+(RC.fV0rec.fParamM.fCovar[9])*(RC.fV0rec.fParamM.fX-RC.fV0rec.fRr)**2)");
  //
  // Normalize distance sigma
  compv0->fTree->SetAlias("SigmaD00","RC.fV0rec.fParamP.fCovar[0]+RC.fV0rec.fParamM.fCovar[0]+RC.fV0rec.fParamP.fCovar[2]+RC.fV0rec.fParamM.fCovar[2]");

  compv0->fTree->SetAlias("SigmaD01","min(RC.fV0rec.fRr-RC.fV0rec.fParamP.fX,0)**2*(RC.fV0rec.fParamP.fCovar[5]+RC.fV0rec.fParamP.fCovar[9])+min(RC.fV0rec.fRr-RC.fV0rec.fParamM.fX,0)**2*(RC.fV0rec.fParamM.fCovar[5]+RC.fV0rec.fParamM.fCovar[9])");
  compv0->fTree->SetAlias("SigmaD0","sqrt(SigmaD00+SigmaD01+0.03**2)");
  //
  // Point Angle sigma
  compv0->fTree->SetAlias("Norm0","sqrt(RC.fV0rec.fPM[0]**2+RC.fV0rec.fPM[1]**2+RC.fV0rec.fPM[2]**2)/Prec");
  compv0->fTree->SetAlias("Norm1","sqrt(RC.fV0rec.fPP[0]**2+RC.fV0rec.fPP[1]**2+RC.fV0rec.fPP[2]**2)/Prec");
  compv0->fTree->SetAlias("SigmaA","(RC.fV0rec.fParamP.fCovar[5]+RC.fV0rec.fParamP.fCovar[9])*(Norm1**2)+(RC.fV0rec.fParamP.fCovar[5]+RC.fV0rec.fParamM.fCovar[9])*(Norm0**2)");
  compv0->fTree->SetAlias("SigmaAP0","sqrt(0.005**2+0.5*((SigmaD0/RC.fV0rec.fRr)**2+SigmaA))");
  //compv0->fTree->SetAlias("SigmaAP0","sqrt(0.0001**2+0.5*((SigmaD0/RC.fV0rec.fRr)**2+SigmaA))");
  //
  // Effective sigma definition
  //  
  compv0->fTree->SetAlias("SigmaAPE0","min((0.005+0.02/(0.1+RC.fV0rec.fRr))*0.4*(0.7+0.3*p12),0.06)");
  compv0->fTree->SetAlias("SigmaAP2","max(min((SigmaAP0+SigmaAPE0)*0.5,1.5*SigmaAPE0),0.5*SigmaAPE0+0.003)");
  compv0->fTree->SetAlias("SigmaDE0","min(sqrt( ((0.02*max(sqrt(RC.fV0rec.fRr)-2.7,0))*p12**2)**2**2+0.06**2),0.5)");
  compv0->fTree->SetAlias("SigmaD2","max(min((SigmaD0+SigmaDE0)*0.5,1.5*SigmaDE0),0.5*SigmaDE0)");

}

void MakeCutsAliases(AliESDComparisonDraw * compv0){
  //
  // Defines likelihood variables
  // LikeAP   - likelihood acoording point angle
  // LikeD    - likelihood according DCA
  // LCausal  - causality likelihood  
  // and Cuts 
  // CLow1 cut - cut on log(LikeAP*LogD)
  // CHigh cut - Causality likehood also used  - cut on log(LikeAP*LikeD*LCausal)
  //
  compv0->fTree->SetAlias("APNorm0","acos(RC.fV0rec.fPointAngle)/SigmaAP0");
  compv0->fTree->SetAlias("DNorm0", "RC.fV0rec.fDist2/SigmaD0");
  compv0->fTree->SetAlias("APNorm2","acos(RC.fV0rec.fPointAngle)/SigmaAP2");
  compv0->fTree->SetAlias("DNorm2", "RC.fV0rec.fDist2/SigmaD2");
  //
  compv0->fTree->SetAlias("LikeAP","(10^-50+exp(-APNorm2**2/2)+0.5*exp(-APNorm2**2/4)+0.25*exp(-APNorm2**2/8))/1.75");
  compv0->fTree->SetAlias("LikeD","10^-50+(exp(-DNorm2*2)+0.5*exp(-DNorm2)+0.25*exp(-DNorm2/2))/1.75");
  //
  compv0->fTree->SetAlias("LikeAP0","10^-50+exp(-APNorm0**2/2)");
  compv0->fTree->SetAlias("LikeD0","10^-50+exp(-DNorm0*2.0)");
  //
  compv0->fTree->SetAlias("LikeAP1","10^-50+(exp(-APNorm0**2/2)+0.5*exp(-APNorm0**2/4))/1.5");
  compv0->fTree->SetAlias("LikeD1","10^-50+(exp(-DNorm0*2.0)+0.3*exp(-DNorm0))/1.3");
  //
  compv0->fTree->SetAlias("LikeAP3","10^-50+exp(-APNorm2**2/2)");
  compv0->fTree->SetAlias("LikeD3","10^-50+exp(-DNorm2*2.0)");
  //
  //
  compv0->fTree->SetAlias("MinCausal","min(RC.fV0rec.fCausality[0],RC.fV0rec.fCausality[1])");
  compv0->fTree->SetAlias("MeanCausal","0.5*(RC.fV0rec.fCausality[0]+RC.fV0rec.fCausality[1])");
  compv0->fTree->SetAlias("MaxCausal","max(RC.fV0rec.fCausality[0],RC.fV0rec.fCausality[1])");
  //  compv0->fTree->SetAlias("LCausal","(1.05-(3*(0.8-exp(-MaxCausal))+(0.8-exp(-MinCausal)))/2)**4");
  compv0->fTree->SetAlias("LCausal","(1.05-(2*(0.8-exp(-max(RC.fV0rec.fCausality[0],RC.fV0rec.fCausality[1])))+2*(0.8-exp(-min(RC.fV0rec.fCausality[0],RC.fV0rec.fCausality[1]))))/2)**4");
  //  compv0->fTree->SetAlias("LCausal","(1.05-(2*(0.8-exp(-max(RC.fV0rec.fCausality[0],RC.fV0rec.fCausality[1])))))**4");
  //  
  //
  compv0->fTree->SetAlias("CLow1","RC.fV0rec.fPointAngle>0.99&&RC.fV0rec.fRr<210&&log(LikeD*LikeAP)>-5.");
  compv0->fTree->SetAlias("CHigh","CLow1&&log(LikeD*LikeAP*LCausal)>-5.0");
  compv0->fTree->SetAlias("CHigh4","CLow1&&log(LikeD*LikeAP*LCausal)>-4.0");
  compv0->fTree->SetAlias("CHigh3","CLow1&&log(LikeD*LikeAP*LCausal)>-3.0");
  compv0->fTree->SetAlias("CHigh2","CLow1&&log(LikeD*LikeAP*LCausal)>-2.0");
  compv0->fTree->SetAlias("CHigh1","CLow1&&log(LikeD*LikeAP*LCausal)>-1.0");
  //  
  compv0->fTree->SetAlias("minCB","min(RC.fV0rec.fCausality[0],RC.fV0rec.fCausality[1])<0.5"); 
}


void MakeCutsAliasesPID(AliESDComparisonDraw * compv0){
  //
  // Cuts optimized for different types of particle
  // use PID of daughter particles as additional criteria 
  //
  compv0->fTree->SetAlias("Spid0","RC.fT1.fESDTrack.fTPCr[0]+RC.fT1.fESDTrack.fTPCr[2]+RC.fT1.fESDTrack.fTPCr[3]+RC.fT1.fESDTrack.fTPCr[4]");
  compv0->fTree->SetAlias("Spid1","RC.fT2.fESDTrack.fTPCr[0]+RC.fT2.fESDTrack.fTPCr[2]+RC.fT2.fESDTrack.fTPCr[3]+RC.fT2.fESDTrack.fTPCr[4]");
  // K0S CUT
  compv0->fTree->SetAlias("Crecks0_0","RC.fT1.fESDTrack.fTPCr[0]/Spid0<0.9&&RC.fT1.fESDTrack.fTPCr[4]/Spid0<0.8&RC.fT1.fESDTrack.fTPCr[3]/Spid0<0.9&&RC.fT2.fESDTrack.fTPCr[0]/Spid1<0.9&&RC.fT2.fESDTrack.fTPCr[4]/Spid1<0.8&&RC.fT2.fESDTrack.fTPCr[3]/Spid1<0.9");
  compv0->fTree->SetAlias("Crecks0_1","RC.fV0rec.fNormDCAPrim[0]>4&&RC.fV0rec.fNormDCAPrim[0]>4");
  compv0->fTree->SetAlias("Crecks0","Crecks0_0&&Crecks0_1&&RC.fV0rec.fRr<25");
  // 
  // LAMBDA CUTS
  compv0->fTree->SetAlias("CrecLambda_0","RC.fT1.fESDTrack.fTPCr[4]/Spid0>0.19&&RC.fT2.fESDTrack.fTPCr[2]/Spid0>0.19");
  compv0->fTree->SetAlias("CrecLambda_1","RC.fT1.fESDTrack.fTPCr[0]/Spid0<0.9&&RC.fT2.fESDTrack.fTPCr[0]/Spid0<0.9");
  compv0->fTree->SetAlias("CrecLambda","CrecLambda_0&&CrecLambda_1");
  //
  // GAMMA CONVERSION CUTS
  compv0->fTree->SetAlias("CrecGamma_0","RC.fT1.fESDTrack.fTPCr[0]/Spid0+RC.fT2.fESDTrack.fTPCr[0]/Spid1>0.4");
  compv0->fTree->SetAlias("CrecGamma_1","RC.fT1.fESDTrack.fTPCr[4]/Spid0<0.8&&RC.fT2.fESDTrack.fTPCr[4]/Spid1<0.8");
  compv0->fTree->SetAlias("CrecGamma_2","RC.fT1.fESDTrack.fTPCr[2]/Spid0<0.9&&RC.fT2.fESDTrack.fTPCr[2]/Spid1<0.9");
  compv0->fTree->SetAlias("CrecGamma","CrecGamma_0&&CrecGamma_1&&CrecGamma_2");
}

//
//----------------------------------------------------------------------------------------------
TH1F * GetInvMass(AliESDComparisonDraw* compv0, TCut cut, Int_t pdg=310, Bool_t all=kFALSE){
  //
  //   return invariant mass histogram
  //   pdg - specify type of particle - (mass window)
  //   all - kTRUE  - all particles shown
  //         kFALSE - only particle with given pdg
  if (pdg==310){
    TH1F * himk0 = new TH1F("himk0","himk0",100,0.46,0.54);
    if (all) {
      compv0->fTree->Draw("RC.fV0rec.GetEffMass(2,2)>>himk0",cut);
    }else{
      compv0->fTree->Draw("RC.fV0rec.GetEffMass(2,2)>>himk0",cut+cks0);
    }
    himk0->SetXTitle("Invariant Mass (GeV)");
    fim.SetParameter(0,0.4975);
    fim.SetParameter(1,0.003);
    fim.SetParameter(2,himk0->GetMaximum());
    fim.SetParameter(3,(himk0->GetBinContent(1)+himk0->GetBinContent(99))*0.5);
    fim.SetParameter(4,0);
    himk0->Fit(&fim);
    himk0->Draw();
    return himk0;
  }
  if (pdg==3122){
    TH1F * hlam0 = new TH1F("hlam0","hlam0",100,1.1,1.13);
    if (all){
      compv0->fTree->Draw("RC.fV0rec.GetEffMass(4,2)>>hlam0",cut);
    }else{
      compv0->fTree->Draw("RC.fV0rec.GetEffMass(4,2)>>hlam0",cut+clambda);
    }
    hlam0->SetXTitle("Invariant Mass (GeV)");
    fim.SetParameter(0,1.11568);
    fim.SetParameter(1,0.0011);
    fim.SetParameter(2,himk0->GetMaximum());
    fim.SetParameter(3,(himk0->GetBinContent(1)+himk0->GetBinContent(99))*0.5);
    fim.SetParameter(4,0);
    hlam0->Fit(&fim);
    hlam0->Draw();
    return hlam0;
  }
  if (pdg==22){
    TH1F * hgamma = new TH1F("hgamma","hgamma",100,0.0,0.2);
    if (all){
      compv0->fTree->Draw("RC.fV0rec.GetEffMass(0,0)>>hgamma",cut);
    }else{
      compv0->fTree->Draw("RC.fV0rec.GetEffMass(0,0)>>hgamma",cut+cgamma);
    }
    hgamma->SetXTitle("Invariant Mass (GeV)");
    
    hgamma->Fit("gaus");
    hgamma->Draw();
    return hgamma;
  }
}

TH1F * GetInvMassPtMin(AliESDComparisonDraw* compv0, TCut cut="RC.fRecStatus==1&&CHigh", Int_t pdg=310)
{
  //
  // effective mass resolution vs. particle momenta
  //
  if (pdg==310){
    compv0->DrawXY("minpt","1000*(RC.fV0rec.GetEffMass(2,2)-0.49767)",cut,"RC.fRecStatus==1",5,0.2,1,-10,10);
  }
  if (pdg==3122){
    compv0->DrawXY("minpt","1000*(RC.fV0rec.GetEffMass(4,2)-1.11568)",cut,"RC.fRecStatus==1",5,0.2,0.5,-10,10);
  }
  TH1F * his = (TH1F*)compv0->fRes->Clone();
  his->SetXTitle("Pt(GeV)");
  his->SetYTitle("Invariant mass sigma (MeV)");
  his->Draw();
  return his;
}

TH1F * GetLogLikelihood(AliESDComparisonDraw* compv0, char * var0="LCausal*LikeAP*LikeD", 
			TCut cut="RC.fRecStatus==1&&CHigh", Int_t min=-10)
{
  //
  // Likelihood
  //
  //char *var0 = "LCausal";
  TH1F * his = new TH1F("loglikelihood","loglikelihood",50,min,0);
  char var[1000];
  sprintf(var,"log(%s)>>loglikelihood",var0);
  compv0->fTree->Draw(var,cut);
  his->SetXTitle("Logarithm of likelihood");
  his->Draw();
  return his;
}

//
//----------------------------------------------------------------------------------------------
TH1F * GetDeltaR(AliESDComparisonDraw* compv0, TCut cut){
  //
  // resolution of vertex position in radial direction 
  //
  TH1F * hdr = new TH1F("hdr","hdr",100,-0.5,0.5);
  compv0->fTree->Draw("RC.fV0rec.fRr-MC.fMCR>>hdr",cut);
  hdr->SetXTitle("Delta r (cm)");
  hdr->Fit("gaus");
  hdr->Draw();
  return hdr;
}

//
//----------------------------------------------------------------------------------------------
TH1F* GetResolP2(AliESDComparisonDraw* compv0, TCut cut, Float_t pmin=0.6, Float_t pmax=2.){
  //
  //
  //
  compv0->DrawXY("MC.fMotherP.P()","100*(MC.fMotherP.P()-Prec)/MC.fMotherP.P()",cut,"RC.fRecStatus==1",5,pmin,pmax,-6,6);
  TH1F * hdp0 = (TH1F*)compv0->fRes->Clone();
  compv0->DrawXY("MC.fMotherP.P()","100*(MC.fMotherP.P()-Pits)/MC.fMotherP.P()",cut,"RC.fRecStatus==1",5,pmin,pmax,-6,6);
  TH1F * hdp1 = (TH1F*)compv0->fRes->Clone();
  compv0->DrawXY("MC.fMotherP.P()","100*(MC.fMotherP.P()-Ptpc)/MC.fMotherP.P()",cut,"RC.fRecStatus==1",5,pmin,pmax,-6,6);
  TH1F * hdp2 = (TH1F*)compv0->fRes->Clone();
  hdp0->SetXTitle("P(GeV)");
  hdp0->SetYTitle("dP/P(%)");
  hdp0->SetMarkerStyle(26);
  hdp1->SetMarkerStyle(27);
  hdp2->SetMarkerStyle(28);  
  hdp0->Draw();
  hdp1->Draw("same");
  hdp2->Draw("same");
  //
  TLegend *legend = new TLegend(0.35,0.15,0.85,0.3, "Momentum resolution");
  legend->SetBorderSize(1);
  legend->AddEntry(hdp0,"Combined V0 finder");
  legend->AddEntry(hdp1,"ESD V0 finder");
  legend->AddEntry(hdp2,"TPC V0 finder");
  //
  legend->Draw();
  return hdp0;
}
//
//----------------------------------------------------------------------------------------------
TH1F* GetResolPhi(AliESDComparisonDraw* compv0, TCut cut, Float_t pmin=0.6, Float_t pmax=2.){
  //
  //
  compv0->DrawXY("MC.fMotherP.P()","1000*(PhiMC-PhiRec)",cut,"RC.fRecStatus==1",5,pmin,pmax,-15,15,20);
  TH1F * hdp0 = (TH1F*)compv0->fRes->Clone();
  compv0->DrawXY("MC.fMotherP.P()","1000*(PhiMC-PhiIts)",cut,"RC.fRecStatus==1",5,pmin,pmax,-15,15,20);
  TH1F * hdp1 = (TH1F*)compv0->fRes->Clone();
  compv0->DrawXY("MC.fMotherP.P()","1000*(PhiMC-PhiTpc)",cut,"RC.fRecStatus==1",5,pmin,pmax,-15,15,20);
  TH1F * hdp2 = (TH1F*)compv0->fRes->Clone();
  hdp0->SetXTitle("P(GeV)");
  hdp0->SetYTitle("d#varphi(mrad)");
  hdp0->SetMarkerStyle(26);
  hdp1->SetMarkerStyle(27);
  hdp2->SetMarkerStyle(28);  
  hdp0->Draw();
  hdp1->Draw("same");
  hdp2->Draw("same");
  //
  TLegend *legend = new TLegend(0.2,0.15,0.55,0.3, "#varphi resolution");
  legend->SetBorderSize(1);
  legend->AddEntry(hdp0,"Combined V0 finder");
  legend->AddEntry(hdp1,"ESD V0 finder");
  legend->AddEntry(hdp2,"TPC V0 finder");
  //
  legend->Draw();
  return hdp0;
}

//
//-----------------------------------------------------------------------------------------------------
TH1F * EffR(AliESDComparisonDraw* compv0, TCut cut, Float_t minr=0.8, Float_t maxr=30, Int_t ndiv=10){  
  //
  // make efficiency plots
  //  
  compv0->Eff("MC.fMCR",cut,"RC.fV0Status==3",ndiv,minr,maxr);
  TH1F * heff0 = (TH1F*)compv0->fRes->Clone();  
  compv0->Eff("MC.fMCR",cut,"RC.fRecStatus==1",ndiv,minr,maxr);
  TH1F * heff1 = (TH1F*)compv0->fRes->Clone();  
  compv0->Eff("MC.fMCR",cut,"RC.fRecStatus==1&&CLow1",ndiv,minr,maxr);
  TH1F * heff2 = (TH1F*)compv0->fRes->Clone();  
  compv0->Eff("MC.fMCR",cut,"RC.fRecStatus==1&&CHigh",ndiv,minr,maxr);
  TH1F * heff3 = (TH1F*)compv0->fRes->Clone();   
  heff0->SetXTitle("Decay radii (cm)");
  heff0->SetMarkerStyle(25);
  heff1->SetMarkerStyle(26);
  heff2->SetMarkerStyle(27);
  heff3->SetMarkerStyle(28);
  heff0->Draw(); heff1->Draw("same");heff2->Draw("same");heff3->Draw("same");
  TLegend *legend = new TLegend(0.2,0.15,0.55,0.3, "Efficiency");
  legend->SetBorderSize(1);
  legend->AddEntry(heff0,"Efficiency to find both decay product");
  legend->AddEntry(heff1,"After rough cuts during tracking");
  legend->AddEntry(heff2,"After resolution cuts");
  legend->AddEntry(heff3,"After resolution and causality cuts");
  //
  legend->Draw();
  return heff3;
}

//
//-----------------------------------------------------------------------------------------------------
TH1F * EffPt(AliESDComparisonDraw* compv0, TCut cut, Float_t min=0.5, Float_t max=1., Int_t ndiv=10){  
  //
  // make efficiency plots
  //  
  compv0->Eff("MC.fMotherP.Pt()",cut,"RC.fV0Status==3",ndiv,min,max);
  TH1F * heff0 = (TH1F*)compv0->fRes->Clone();  
  compv0->Eff("MC.fMotherP.Pt()",cut,"RC.fRecStatus==1",ndiv,min,max);
  TH1F * heff1 = (TH1F*)compv0->fRes->Clone();  
  compv0->Eff("MC.fMotherP.Pt()",cut,"RC.fRecStatus==1&&CLow1",ndiv,min,max);
  TH1F * heff2 = (TH1F*)compv0->fRes->Clone();  
  compv0->Eff("MC.fMotherP.Pt()",cut,"RC.fRecStatus==1&&CHigh",ndiv,min,max);
  TH1F * heff3 = (TH1F*)compv0->fRes->Clone();   
  heff0->SetXTitle("Pt(GeV)");
  heff0->SetMarkerStyle(25);
  heff1->SetMarkerStyle(26);
  heff2->SetMarkerStyle(27);
  heff3->SetMarkerStyle(28);
  heff0->Draw(); heff1->Draw("same");heff2->Draw("same");heff3->Draw("same");
  TLegend *legend = new TLegend(0.2,0.15,0.55,0.3, "Efficiency");
  legend->SetBorderSize(1);
  legend->AddEntry(heff0,"Efficiency to find both decay product");
  legend->AddEntry(heff1,"After rough cuts during tracking");
  legend->AddEntry(heff2,"After resolution cuts");
  legend->AddEntry(heff3,"After resolution and causality cuts");
  //
  legend->Draw();
  return heff3;
}

//
//-----------------------------------------------------------------------------------------------------
TH1F * EffMinPt(AliESDComparisonDraw* compv0, TCut cut, Float_t min=0.5, Float_t max=1., Int_t ndiv=10){  
  //
  // make efficiency plots
  //  
  compv0->Eff("min(MC.fMCm.fParticle.Pt(),MC.fMCd.fParticle.Pt())",cut,"RC.fV0Status==3",ndiv,min,max);
  TH1F * heff0 = (TH1F*)compv0->fRes->Clone();  
  compv0->Eff("min(MC.fMCm.fParticle.Pt(),MC.fMCd.fParticle.Pt())",cut,"RC.fRecStatus==1",ndiv,min,max);
  TH1F * heff1 = (TH1F*)compv0->fRes->Clone();  
  compv0->Eff("min(MC.fMCm.fParticle.Pt(),MC.fMCd.fParticle.Pt())",cut,"RC.fRecStatus==1&&CLow1",ndiv,min,max);
  TH1F * heff2 = (TH1F*)compv0->fRes->Clone();  
  compv0->Eff("min(MC.fMCm.fParticle.Pt(),MC.fMCd.fParticle.Pt())",cut,"RC.fRecStatus==1&&CHigh",ndiv,min,max);
  TH1F * heff3 = (TH1F*)compv0->fRes->Clone();   
  heff0->SetXTitle("Pt(GeV)");
  heff0->SetMarkerStyle(25);
  heff1->SetMarkerStyle(26);
  heff2->SetMarkerStyle(27);
  heff3->SetMarkerStyle(28);
  heff0->Draw(); heff1->Draw("same");heff2->Draw("same");heff3->Draw("same");
  TLegend *legend = new TLegend(0.2,0.15,0.55,0.3, "Efficiency");
  legend->SetBorderSize(1);
  legend->AddEntry(heff0,"Efficiency to find both decay product");
  legend->AddEntry(heff1,"After rough cuts during tracking");
  legend->AddEntry(heff2,"After resolution cuts");
  legend->AddEntry(heff3,"After resolution and causality cuts");
  //
  legend->Draw();
  return heff3;
}


TH1F * GetCumulative(AliESDComparisonDraw* compv0, char * var0,TCut cut = "RC.fRecStatus==1", Float_t min=0, Float_t max=30, Int_t ndiv=20)
{
  //
  // Make cumulative function
  //
  TH1F his("his","his",ndiv,min,max);
  char var[1000];
  sprintf(var,"%s>>his",var0);
  compv0->fTree->Draw(var,cut);
  TH1F *hcumul = new TH1F("hcumul","hcumul",ndiv,min,max);
  Float_t all = (Float_t)his.GetSum();
  {
    for (Int_t i=0;i<=ndiv;i++){
      hcumul->SetBinContent(i,his.Integral(0,i)/all);
    }
  }
  hcumul->SetFillColor(0);
  hcumul->SetXTitle("Threshold (unit)");
  hcumul->SetYTitle("Cumulative function");
  hcumul->Draw();
  return hcumul;
}

TH1F * CompareLikeCumulative(AliESDComparisonDraw* compv0,TCut cut = "RC.fRecStatus==1", Float_t min=0, Float_t max=30, Int_t ndiv=20)
{
  //
  //
  //
  TH1F * h0 = GetCumulative(compv0,"-log(LikeAP*LikeD*LCausal)",cut,min,max,ndiv);
  TH1F * h1 = GetCumulative(compv0,"-log(LikeAP1*LikeD1*LCausal)",cut,min,max,ndiv);
  TH1F * h3 = GetCumulative(compv0,"-log(LikeAP3*LikeD3*LCausal)",cut,min,max,ndiv);
  h0->SetXTitle("Threshold (unit)");
  h0->SetMaximum(1.1);
  h0->SetLineColor(1); h0->SetLineStyle(1);
  h1->SetLineColor(2); h1->SetLineStyle(2);
  h3->SetLineColor(3); h3->SetLineStyle(3);
  h0->Draw();
  h1->Draw("same");
  h3->Draw("same");
  TLegend *legend = new TLegend(0.4,0.12,0.88,0.5, "Cumulative as function of threshold");
  legend->SetBorderSize(1);
  legend->AddEntry(h0,"Likelihhod (P0 parameterization) with three components");
  legend->AddEntry(h1,"Likelihood (P1 parameterization) with one component");
  legend->AddEntry(h3,"Likelihood (P0 parameterization) with one component");
  //
  legend->Draw();
  return h0;
}


TH1F * CompareLikeCumulative2(AliESDComparisonDraw* compv0,TCut cut = "RC.fRecStatus==1", Float_t min=0, Float_t max=30, Int_t ndiv=20)
{
  //
  //
  //
  TH1F * h0 = GetCumulative(compv0,"-log(LikeAP*LikeD*LCausal)",cut,min,max,ndiv);
  TH1F * h1 = GetCumulative(compv0,"-log(LikeAP1*LikeD1*LCausal)",cut,min,max,ndiv);
  TH1F * h3 = GetCumulative(compv0,"-log(LikeAP3*LikeD3*LCausal)",cut,min,max,ndiv);
  TH1F * h4 = GetCumulative(compv0,"-log(RC.fV0rec.GetLikelihoodC(0,0)*RC.fV0rec.GetLikelihoodAP(2,2)*RC.fV0rec.GetLikelihoodD(2,2))",cut,min,max,ndiv);
  h0->SetXTitle("Threshold (unit)");
  h0->SetMaximum(1.1);
  h0->SetLineColor(1); h0->SetLineStyle(1);
  h1->SetLineColor(2); h1->SetLineStyle(2);
  h3->SetLineColor(3); h3->SetLineStyle(3);
  h4->SetLineColor(4); h4->SetLineStyle(4);

  h0->Draw();
  h1->Draw("same");
  h3->Draw("same");
  h4->Draw("same");
  TLegend *legend = new TLegend(0.4,0.12,0.88,0.5, "Cumulative as function of threshold");
  legend->SetBorderSize(1);
  legend->AddEntry(h0,"Likelihhod (P0 parameterization) with three components");
  legend->AddEntry(h1,"Likelihood (P1 parameterization) with one component");
  legend->AddEntry(h3,"Likelihood (P0 parameterization) with one component");
  legend->AddEntry(h4,"Likelihood default");
  //
  legend->Draw();
  return h0;
}



TH1F * CompareCumulativeSec(AliESDComparisonDraw* compv0, Float_t min=0, Float_t max=30, Int_t ndiv=20)
{
  //
  //
  //
  // Float_t min=0; Float_t max=20; Int_t ndiv=1000;
  TH1F * hks0 = GetCumulative(&compv0,"abs(RC.fT1.fITStrack.fP0/sqrt(RC.fT1.fITStrack.fC00))","RC.fV0Status==3&&MC.fMCR>0.5"+cint+cks0,min,max,ndiv);
  TH1F * hks1 = GetCumulative(&compv0,"sqrt((RC.fT1.fITStrack.fP0/sqrt(RC.fT1.fITStrack.fC00))**2+((RC.fT1.fITStrack.fP1-MC.fMotherP.fVz)/sqrt(RC.fT1.fITStrack.fC11))**2)","RC.fV0Status==3&&MC.fMCR>0.5"+cint+cks0,min,max,ndiv);

  TH1F * hlam0 = GetCumulative(&compv0,"abs(RC.fT1.fITStrack.fP0/sqrt(RC.fT1.fITStrack.fC00))","RC.fV0Status==3&&MC.fMCR>0.5"+cint+clam,min,max,ndiv);
  TH1F * hlam1 = GetCumulative(&compv0,"sqrt((RC.fT1.fITStrack.fP0/sqrt(RC.fT1.fITStrack.fC00))**2+((RC.fT1.fITStrack.fP1-MC.fMotherP.fVz)/sqrt(RC.fT1.fITStrack.fC11))**2)","RC.fV0Status==3&&MC.fMCR>0.5"+cint+clam,min,max,ndiv);
  
  TH1F * hgamma0 = GetCumulative(&compv0,"abs(RC.fT1.fITStrack.fP0/sqrt(RC.fT1.fITStrack.fC00))","RC.fV0Status==3&&MC.fMCR>0.5"+cint+cgamma,min,max,ndiv);
  TH1F * hgamma1 = GetCumulative(&compv0,"sqrt((RC.fT1.fITStrack.fP0/sqrt(RC.fT1.fITStrack.fC00))**2+((RC.fT1.fITStrack.fP1-MC.fMotherP.fVz)/sqrt(RC.fT1.fITStrack.fC11))**2)","RC.fV0Status==3&&MC.fMCR>0.5"+cint+cgamma,min,max,ndiv);
  
  hks0->SetXTitle("Threshold (sigma)");
  hks0->SetYTitle("Cumulative probability");
  hks0->SetMaximum(0.8);
  hks0->SetLineColor(1); hks0->SetLineStyle(1); hks0->SetLineWidth(1)
  hks1->SetLineColor(1); hks1->SetLineStyle(2); hks1->SetLineWidth(2);
  hlam0->SetLineColor(1); hlam0->SetLineStyle(3); hlam0->SetLineWidth(1);
  hlam1->SetLineColor(1); hlam1->SetLineStyle(4); hlam1->SetLineWidth(2)
  hgamma0->SetLineColor(1); hgamma0->SetLineStyle(5); hgamma0->SetLineWidth(1);
  hgamma1->SetLineColor(1); hgamma1->SetLineStyle(6);hgamma1->SetLineWidth(2);

  hks0->Draw();
  hks1->Draw("same");
  hlam0->Draw("same");
  hlam1->Draw("same");
  hgamma0->Draw("same");
  hgamma1->Draw("same");

  TLegend *legend = new TLegend(0.4,0.12,0.88,0.5, "Cumulative probability");
  legend->SetBorderSize(1);
  legend->AddEntry(hks0,"K^{0}_{s} decay - normalized DCA in r#varphi");
  legend->AddEntry(hks1,"K^{0}_{s} decay - normalized DCA 3D");
  legend->AddEntry(hlam0,"#Lambda decay - normalized DCA in r#varphi");
  legend->AddEntry(hlam1,"#Lambda decay - normalized DCA 3D");
  legend->AddEntry(hgamma0,"#gamma conversion - normalized DCA in r#varphi");
  legend->AddEntry(hgamma1,"#gamma conversion - normalized DCA 3D");
  //
  legend->Draw();
  return h0;
}




