# **************************************************************************
# * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
# *                                                                        *
# * Author: The ALICE Off-line Project.                                    *
# * Contributors are mentioned in the code where appropriate.              *
# *                                                                        *
# * Permission to use, copy, modify and distribute this software and its   *
# * documentation strictly for non-commercial purposes is hereby granted   *
# * without fee, provided that the above copyright notice appears in all   *
# * copies and that both the copyright notice and this permission notice   *
# * appear in the supporting documentation. The authors make no claims     *
# * about the suitability of this software for any purpose. It is          *
# * provided "as is" without express or implied warranty.                  *
# **************************************************************************/

#Module
set(MODULE EvtGen)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/TEvtGen/${MODULE})

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/TEvtGen
)

# Sources - alphabetical order
set(SRCS
  EvtGen.cpp 
  EvtGenBase/Evt3Rank3C.cpp
  EvtGenBase/EvtAbsLineShape.cpp
  EvtGenBase/EvtAmp.cpp
  EvtGenBase/EvtAmpIndex.cpp
  EvtGenBase/EvtAmpSubIndex.cpp
  EvtGenBase/EvtBlattWeisskopf.cpp
  EvtGenBase/EvtBreitWignerPdf.cpp
  EvtGenBase/EvtCGCoefSingle.cpp
  EvtGenBase/EvtComplex.cpp
  EvtGenBase/EvtConst.cpp
  EvtGenBase/EvtCPUtil.cpp
  EvtGenBase/EvtCyclic3.cpp
  EvtGenBase/EvtDalitzCoord.cpp
  EvtGenBase/EvtDalitzFlatPdf.cpp
  EvtGenBase/EvtDalitzPlot.cpp
  EvtGenBase/EvtDalitzPoint.cpp
  EvtGenBase/EvtDalitzReso.cpp
  EvtGenBase/EvtDalitzResPdf.cpp
  EvtGenBase/EvtDecayAmp.cpp
  EvtGenBase/EvtDecayBase.cpp
  EvtGenBase/EvtDecayIncoherent.cpp
  EvtGenBase/EvtDecayMode.cpp
  EvtGenBase/EvtDecayParm.cpp
  EvtGenBase/EvtDecayProb.cpp
  EvtGenBase/EvtDecayTable.cpp
  EvtGenBase/EvtdFunction.cpp
  EvtGenBase/EvtdFunctionSingle.cpp
  EvtGenBase/EvtDiLog.cpp
  EvtGenBase/EvtDiracParticle.cpp
  EvtGenBase/EvtDiracSpinor.cpp
  EvtGenBase/EvtEvalHelAmp.cpp
  EvtGenBase/EvtExtGeneratorCommandsTable.cpp
  EvtGenBase/EvtFlatLineShape.cpp
  EvtGenBase/EvtFlatte.cpp
  EvtGenBase/EvtGammaMatrix.cpp
  EvtGenBase/EvtGenKine.cpp
  EvtGenBase/EvtHepMCEvent.cpp
  EvtGenBase/EvtHighSpinParticle.cpp
  EvtGenBase/EvtId.cpp
  EvtGenBase/EvtIdSet.cpp
  EvtGenBase/EvtIncoherentMixing.cpp
  EvtGenBase/EvtIntegPdf1D.cpp
  EvtGenBase/EvtIntervalFlatPdf.cpp
  EvtGenBase/EvtKine.cpp
  EvtGenBase/EvtLASSAmp.cpp
  EvtGenBase/EvtManyDeltaFuncLineShape.cpp
  EvtGenBase/EvtMassAmp.cpp
  EvtGenBase/EvtMBreitWigner.cpp
  EvtGenBase/EvtMHelAmp.cpp
  EvtGenBase/EvtMNode.cpp
  EvtGenBase/EvtModelAlias.cpp
  EvtGenBase/EvtModel.cpp
  EvtGenBase/EvtMParticle.cpp
  EvtGenBase/EvtMRes.cpp
  EvtGenBase/EvtMTree.cpp
  EvtGenBase/EvtMTrivialLS.cpp
  EvtGenBase/EvtMultiChannelParser.cpp
  EvtGenBase/EvtNeutrinoParticle.cpp
  EvtGenBase/EvtNonresonantAmp.cpp
  EvtGenBase/EvtOrthogVector.cpp
  EvtGenBase/EvtParser.cpp
  EvtGenBase/EvtParserXml.cpp
  EvtGenBase/EvtParticle.cpp
  EvtGenBase/EvtParticleDecay.cpp
  EvtGenBase/EvtParticleDecayList.cpp
  EvtGenBase/EvtParticleFactory.cpp
  EvtGenBase/EvtPartProp.cpp
  EvtGenBase/EvtPDL.cpp
  EvtGenBase/EvtPhotonParticle.cpp
  EvtGenBase/EvtPoint1D.cpp
  EvtGenBase/EvtPropBreitWigner.cpp
  EvtGenBase/EvtPropBreitWignerRel.cpp
  EvtGenBase/EvtPropFlatte.cpp
  EvtGenBase/EvtPropGounarisSakurai.cpp
  EvtGenBase/EvtPto3PAmp.cpp
  EvtGenBase/EvtPto3PAmpFactory.cpp
  EvtGenBase/EvtPto3PAmpSmpResolution.cpp
  EvtGenBase/EvtRadCorr.cpp
  EvtGenBase/EvtRandom.cpp
  EvtGenBase/EvtRanf.cpp
  EvtGenBase/EvtRanFor.cpp
  EvtGenBase/EvtRaritaSchwinger.cpp
  EvtGenBase/EvtRaritaSchwingerParticle.cpp
  EvtGenBase/EvtRelBreitWignerBarrierFact.cpp
  EvtGenBase/EvtReport.cpp
  EvtGenBase/EvtResonance2.cpp
  EvtGenBase/EvtResonance.cpp
  EvtGenBase/EvtScalarParticle.cpp
  EvtGenBase/EvtSecondary.cpp
  EvtGenBase/EvtSemiLeptonicAmp.cpp
  EvtGenBase/EvtSemiLeptonicBaryonAmp.cpp
  EvtGenBase/EvtSemiLeptonicScalarAmp.cpp
  EvtGenBase/EvtSemiLeptonicTensorAmp.cpp
  EvtGenBase/EvtSemiLeptonicVectorAmp.cpp
  EvtGenBase/EvtSimpleRandomEngine.cpp
  EvtGenBase/EvtSpinAmp.cpp
  EvtGenBase/EvtSpinDensity.cpp
  EvtGenBase/EvtSpinType.cpp
  EvtGenBase/EvtStdHep.cpp
  EvtGenBase/EvtStringParticle.cpp
  EvtGenBase/EvtSymTable.cpp
  EvtGenBase/EvtTensor3C.cpp
  EvtGenBase/EvtTensor4C.cpp
  EvtGenBase/EvtTensorParticle.cpp
  EvtGenBase/EvtTwoBodyKine.cpp
  EvtGenBase/EvtTwoBodyVertex.cpp
  EvtGenBase/EvtValError.cpp
  EvtGenBase/EvtVector3C.cpp
  EvtGenBase/EvtVector3R.cpp
  EvtGenBase/EvtVector4C.cpp
  EvtGenBase/EvtVector4R.cpp
  EvtGenBase/EvtVectorParticle.cpp
  EvtGenModels/EvtBaryonPCR.cpp
  EvtGenModels/EvtBaryonPCRFF.cpp
  EvtGenModels/EvtBBScalar.cpp
  EvtGenModels/EvtBcBsNPi.cpp
  EvtGenModels/EvtBcBsStarNPi.cpp
  EvtGenModels/EvtBcPsiNPi.cpp
  EvtGenModels/EvtBcToNPi.cpp
  EvtGenModels/EvtBCVFF.cpp
  EvtGenModels/EvtBcVMuNu.cpp
  EvtGenModels/EvtBcVNpi.cpp
  EvtGenModels/EvtBHadronic.cpp
  EvtGenModels/EvtBsquark.cpp
  EvtGenModels/EvtBto2piCPiso.cpp
  EvtGenModels/EvtBTo3piCP.cpp
  EvtGenModels/EvtBTo4piCP.cpp
  EvtGenModels/EvtBToDDalitzCPK.cpp
  EvtGenModels/EvtBtoKD3P.cpp
  EvtGenModels/EvtBtoKpiCPiso.cpp
  EvtGenModels/EvtBToKpipiCP.cpp
  EvtGenModels/EvtBToPlnuBK.cpp
  EvtGenModels/EvtBToPlnuBKFF.cpp
  EvtGenModels/EvtbTosllAli.cpp
  EvtGenModels/EvtbTosllAliFF.cpp
  EvtGenModels/EvtbTosllAmp.cpp
  EvtGenModels/EvtbTosllBall.cpp
  EvtGenModels/EvtbTosllBallFF.cpp
  EvtGenModels/EvtbTosllScalarAmp.cpp
  EvtGenModels/EvtbTosllVectorAmp.cpp
  EvtGenModels/EvtBToVlnuBall.cpp
  EvtGenModels/EvtBToVlnuBallFF.cpp
  EvtGenModels/EvtBtoXsEtap.cpp
  EvtGenModels/EvtBtoXsgammaAbsModel.cpp
  EvtGenModels/EvtBtoXsgammaAliGreub.cpp
  EvtGenModels/EvtBtoXsgamma.cpp
  EvtGenModels/EvtBtoXsgammaFermiUtil.cpp
  EvtGenModels/EvtBtoXsgammaFixedMass.cpp
  EvtGenModels/EvtBtoXsgammaFlatEnergy.cpp
  EvtGenModels/EvtBtoXsgammaKagan.cpp
  EvtGenModels/EvtBtoXsgammaRootFinder.cpp
  EvtGenModels/EvtBtoXsll.cpp
  EvtGenModels/EvtBtoXsllUtil.cpp
  EvtGenModels/EvtCBTo3piMPP.cpp
  EvtGenModels/EvtCBTo3piP00.cpp
  EvtGenModels/EvtD0gammaDalitz.cpp
  EvtGenModels/EvtD0mixDalitz.cpp
  EvtGenModels/EvtDalitzTable.cpp
  EvtGenModels/EvtDDalitz.cpp
  EvtGenModels/EvtDMix.cpp
  EvtGenModels/EvtEta2MuMuGamma.cpp
  EvtGenModels/EvtEtaDalitz.cpp
  EvtGenModels/EvtFlatQ2.cpp
  EvtGenModels/EvtGenericDalitz.cpp
  EvtGenModels/EvtGoityRoberts.cpp
  EvtGenModels/EvtHelAmp.cpp
  EvtGenModels/EvtHQET2.cpp
  EvtGenModels/EvtHQET2FF.cpp
  EvtGenModels/EvtHQET.cpp
  EvtGenModels/EvtHQETFF.cpp
  EvtGenModels/EvtHypNonLepton.cpp
  EvtGenModels/EvtISGW2.cpp
  EvtGenModels/EvtISGW2FF.cpp
  EvtGenModels/EvtISGW.cpp
  EvtGenModels/EvtISGWFF.cpp
  EvtGenModels/EvtItgAbsFunction.cpp
  EvtGenModels/EvtItgAbsIntegrator.cpp
  EvtGenModels/EvtItgFourCoeffFcn.cpp
  EvtGenModels/EvtItgFunction.cpp
  EvtGenModels/EvtItgPtrFunction.cpp
  EvtGenModels/EvtItgSimpsonIntegrator.cpp
  EvtGenModels/EvtItgThreeCoeffFcn.cpp
  EvtGenModels/EvtItgTwoCoeffFcn.cpp
  EvtGenModels/EvtKKLambdaC.cpp
  EvtGenModels/EvtKKLambdaCFF.cpp
  EvtGenModels/EvtKstarnunu.cpp
  EvtGenModels/EvtKstarstargamma.cpp
  EvtGenModels/EvtLambdaB2LambdaV.cpp
  EvtGenModels/EvtLambdaP_BarGamma.cpp
  EvtGenModels/EvtLb2Lll.cpp
  EvtGenModels/EvtLNuGamma.cpp
  EvtGenModels/EvtMelikhov.cpp
  EvtGenModels/EvtMelikhovFF.cpp
  EvtGenModels/EvtModelReg.cpp
  EvtGenModels/EvtMultibody.cpp
  EvtGenModels/EvtOmegaDalitz.cpp
  EvtGenModels/EvtPartWave.cpp
  EvtGenModels/EvtPFermi.cpp
  EvtGenModels/EvtPhiDalitz.cpp
  EvtGenModels/EvtPhsp.cpp
  EvtGenModels/EvtPi0Dalitz.cpp
  EvtGenModels/EvtPropSLPole.cpp
  EvtGenModels/EvtPto3P.cpp
  EvtGenModels/EvtPVVCPLH.cpp
  EvtGenModels/EvtSingleParticle.cpp
  EvtGenModels/EvtSLBKPole.cpp
  EvtGenModels/EvtSLBKPoleFF.cpp
  EvtGenModels/EvtSll.cpp
  EvtGenModels/EvtSLN.cpp
  EvtGenModels/EvtSLPole.cpp
  EvtGenModels/EvtSLPoleFF.cpp
  EvtGenModels/EvtSSDCP.cpp
  EvtGenModels/EvtSSD_DirectCP.cpp
  EvtGenModels/EvtSSSCP.cpp
  EvtGenModels/EvtSSSCPpng.cpp
  EvtGenModels/EvtSSSCPT.cpp
  EvtGenModels/EvtSTSCP.cpp
  EvtGenModels/EvtSTS.cpp
  EvtGenModels/EvtSVPCP.cpp
  EvtGenModels/EvtSVP.cpp
  EvtGenModels/EvtSVPHelAmp.cpp
  EvtGenModels/EvtSVSCP.cpp
  EvtGenModels/EvtSVSCPiso.cpp
  EvtGenModels/EvtSVSCPLH.cpp
  EvtGenModels/EvtSVS.cpp
  EvtGenModels/EvtSVSNONCPEIGEN.cpp
  EvtGenModels/EvtSVVCP.cpp
  EvtGenModels/EvtSVVCPLH.cpp
  EvtGenModels/EvtSVVHelAmp.cpp
  EvtGenModels/EvtSVVHelCPMix.cpp
  EvtGenModels/EvtSVVNONCPEIGEN.cpp
  EvtGenModels/EvtTauHadnu.cpp
  EvtGenModels/EvtTaulnunu.cpp
  EvtGenModels/EvtTauScalarnu.cpp
  EvtGenModels/EvtTauVectornu.cpp
  EvtGenModels/EvtTSS.cpp
  EvtGenModels/EvtTVP.cpp
  EvtGenModels/EvtTVSPwave.cpp
  EvtGenModels/EvtVectorIsr.cpp
  EvtGenModels/EvtVll.cpp
  EvtGenModels/EvtVPHOtoV.cpp
  EvtGenModels/EvtVPHOtoVISR.cpp
  EvtGenModels/EvtVPHOtoVISRHi.cpp
  EvtGenModels/EvtVSPPwave.cpp
  EvtGenModels/EvtVSSBMixCPT.cpp
  EvtGenModels/EvtVSS.cpp
  EvtGenModels/EvtVSSMix.cpp
  EvtGenModels/EvtVtoSll.cpp
  EvtGenModels/EvtVubAC.cpp
  EvtGenModels/EvtVubBLNP.cpp
  EvtGenModels/EvtVubBLNPHybrid.cpp
  EvtGenModels/EvtVub.cpp
  EvtGenModels/EvtVubdGamma.cpp
  EvtGenModels/EvtVubHybrid.cpp
  EvtGenModels/EvtVubNLO.cpp
  EvtGenModels/EvtVVP.cpp
  EvtGenModels/EvtVVpipi.cpp
  EvtGenModels/EvtVVPIPI_WEIGHTED.cpp
  EvtGenModels/EvtVVSPwave.cpp
  EvtGenModels/EvtWilsonCoefficients.cpp
  EvtGenModels/EvtWnPi.cpp
  EvtGenModels/EvtXPsiGamma.cpp
  EvtGenModels/EvtY3SToY1SpipiMoxhay.cpp
  EvtGenModels/EvtYmSToYnSpipiCLEO.cpp
)

set ( FSRCS
  EvtGenModels/EvtBTo3pi.F
  EvtGenModels/EvtBToKpipi.F
  EvtGenModels/EvtBTo3piMPP.F
  EvtGenModels/EvtLi2Spence.F
  EvtGenModels/EvtBTo3piP00.F
)

# No dictionary

# Add a library to the project using the specified source files
add_library(${MODULE} SHARED ${SRCS} ${FSRCS})

# Do not generate the ROOT map
# Dependencies
set(LIBDEPS HepMC)

# Additional compilation flags
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "")
target_link_libraries(${MODULE} )

# Linking the library
target_link_libraries(${MODULE} ${LIBDEPS})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE} 
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)
install( FILES DECAY.DEC evt.pdl DESTINATION TEvtGen/EvtGen )
install( DIRECTORY DecayTable DESTINATION TEvtGen/EvtGen )
