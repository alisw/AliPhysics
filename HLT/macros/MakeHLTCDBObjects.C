#include "HLT/global/physics/macros/makeConfigurationObjectMultiplicityCorrelations.C"
#include "HLT/global/physics/macros/makeConfigurationObjectdNdPtAnalysis.C"
#include "HLT/exa/makeComponentConfigurationObject.C"
#include "HLT/trigger/macros/makeTriggerConfigurationObject.C"
#include "HLT/TPCLib/HWCFemulator/macro/makeConfigurationObjectTPCHWClusterFinder.C"
#include "HLT/VZERO/macros/makeConfigurationObjectVZEROReconstruction.C"
#include "HLT/ZDC/macros/makeConfigurationObjectZDCReconstruction.C"
#include "HLT/global/macros/makeGlobalHistoConfigObject.C"
#include "HLT/MUON/macros/rootlogon.C"
#include "HLT/MUON/macros/CreateDefaultCDBEntries.C"

void MakeHLTCDBObjects()
{
  //makeConfigurationObjectMultiplicityCorrelations("-addTrigger CPBI1 -addTrigger CPBI2", "", 0, AliCDBRunRange::Infinity(), "");
  makeConfigurationObjectdNdPtAnalysis();             // HLT/ConfigAnalysis/dNdPtAnalysis                           
  makeComponentConfigurationObject("HLT/ConfigHLT/BarrelHighMultiplicity", "-mintracks 100", "", 0, 999999999); // HLT/ConfigHLT/BarrelHighMultiplicity
  makeComponentConfigurationObject("HLT/ConfigHLT/HighMultiplicityTrigger", "-mintracks 10", "", 0, 999999999); // HLT/ConfigHLT/BarrelMultiplicityTrigger
  makeComponentConfigurationObject("HLT/ConfigHLT/BarrelPt_v01", "-mintracks 1 -minpt 0.5", "", 0, 999999999);  // HLT/ConfigHLT/BarrelPt_v01
  makeComponentConfigurationObject("HLT/ConfigHLT/BarrelPt_v02", "-mintracks 1 -minpt 1.0", "", 0, 999999999);  // HLT/ConfigHLT/BarrelPt_v02
  makeComponentConfigurationObject("HLT/ConfigHLT/BarrelPt_v03", "-mintracks 1 -minpt 5.0", "", 0, 999999999);  // HLT/ConfigHLT/BarrelPt_v03
  makeComponentConfigurationObject("HLT/ConfigHLT/CosmicsTrigger", "", "", 0, 999999999);                       // HLT/ConfigHLT/CosmicsTrigger
  makeComponentConfigurationObject("HLT/ConfigHLT/D0Trigger", "-pt 0.7 -dca 0.03 -invmass 0.033 -costhetastar 0.8 -d0 100000 -d0d0 -0.00005 -cospoint 0.8", "", 0, 999999999); // HLT/ConfigHLT/D0Trigger
  makeComponentConfigurationObject("HLT/ConfigHLT/EmcalClusterEnergyTrigger", "-energy 1.0", "", 0, 999999999); // HLT/ConfigHLT/EmcalClusterEnergyTrigger
  makeComponentConfigurationObject("HLT/ConfigHLT/EmcalElectronTrigger", "-energy 1. -minEoverP .8 -maxEoverP 1.2 -dEta 0.05 -dPhi 0.05", "", 0, 999999999);  // HLT/ConfigHLT/EmcalElectronTrigger
  makeComponentConfigurationObject("HLT/ConfigHLT/EmcalJetTrigger", "-energy 20.", "", 0, 999999999);           // HLT/ConfigHLT/EmcalJetTrigger
  makeComponentConfigurationObject("HLT/ConfigHLT/GammaConversionTrigger", "", "", 0, 999999999);               // HLT/ConfigHLT/GammaConversionTrigger
  makeComponentConfigurationObject("HLT/ConfigHLT/GlobalOfflineVertexer", "", "", 0, 999999999);                // HLT/ConfigHLT/GlobalOfflineVertexer
  makeComponentConfigurationObject("HLT/ConfigHLT/GlobalTrackMatcher", "-method 1", "", 0, 999999999);          // HLT/ConfigHLT/GlobalTrackMatcher
  makeComponentConfigurationObject("HLT/ConfigHLT/ITSMultiplicityTrigger", "-nclusters 4", "", 0, 999999999);   // HLT/ConfigHLT/ITSMultiplicityTrigger
  makeComponentConfigurationObject("HLT/ConfigHLT/PhosClusterEnergyTrigger", "-energy 1.0", "", 0, 999999999);  // HLT/ConfigHLT/PhosClusterEnergyTrigger
  makeComponentConfigurationObject("HLT/ConfigHLT/PhosMipTrigger", "-emin 0.080 -emax 0.250 -ncells 4", "", 0, 999999999);  // HLT/ConfigHLT/PhosMipTrigger
  makeComponentConfigurationObject("HLT/ConfigHLT/PrimaryVertexFinder", "", "", 0, 999999999);                  // HLT/ConfigHLT/PrimaryVertexFinder
  makeComponentConfigurationObject("HLT/ConfigHLT/V0Finder", "", "", 0, 999999999);                             // HLT/ConfigHLT/V0Finder
  makeComponentConfigurationObject("HLT/ConfigHLT/V0Histo", "", "", 0, 999999999);                              // HLT/ConfigHLT/V0Histo
  makeComponentConfigurationObject("HLT/ConfigITS/ITSClusterFinderSDD", "", "", 0, 999999999);                  // HLT/ConfigITS/ITSClusterFinderSDD
  makeComponentConfigurationObject("HLT/ConfigITS/ITSClusterFinderSPD", "", "", 0, 999999999);                  // HLT/ConfigITS/ITSClusterFinderSPD
  makeComponentConfigurationObject("HLT/ConfigITS/ITSClusterFinderSSD", "", "", 0, 999999999);                  // HLT/ConfigITS/ITSClusterFinderSSD
  makeComponentConfigurationObject("HLT/ConfigITS/ITSTracker", "", "", 0, 999999999);                           // HLT/ConfigITS/ITSTracker
  makeComponentConfigurationObject("HLT/ConfigPHOS/PHOSClusterizer", "-digitthreshold 0.005 -recpointthreshold 0.1", "", 0, 999999999);   // HLT/ConfigPHOS/PHOSClusterizer
  makeComponentConfigurationObject("HLT/ConfigPHOS/PHOSDigitMaker", "-sethighgainfactor 0.005 -setlowgainfactor 0.08", "", 0, 999999999); // HLT/ConfigPHOS/PHOSDigitMaker
  makeComponentConfigurationObject("HLT/ConfigPHOS/PHOSRawAnalyzer", "", "", 0, 999999999);                     // HLT/ConfigPHOS/PHOSRawAnalyzer
  makeComponentConfigurationObject("HLT/ConfigSample/SampleCalibration", "", "", 0, 999999999);                 // HLT/ConfigSample/SampleCalibration
  makeComponentConfigurationObject("HLT/ConfigSample/SampleComponent1", "-config1 config-param -config2", "", 0, 999999999);   // HLT/ConfigSample/SampleComponent1
  makeComponentConfigurationObject("HLT/ConfigSample/SampleESDAnalysis", "", "", 0, 999999999);                 // HLT/ConfigSample/SampleESDAnalysis
  makeComponentConfigurationObject("HLT/ConfigSample/SampleRawAnalysis", "", "", 0, 999999999);                 // HLT/ConfigSample/SampleRawAnalysis
  makeComponentConfigurationObject("HLT/ConfigTPC/TPCCAGlobalMerger", "", "", 0, 999999999);                    // HLT/ConfigTPC/TPCCAGlobalMerger
  makeComponentConfigurationObject("HLT/ConfigTPC/TPCCalibTime", "-output-size 50000", "", 0, 999999999);       // HLT/ConfigTPC/TPCCalibTime
  makeComponentConfigurationObject("HLT/ConfigTPC/TPCCATracker", "", "", 0, 999999999);                         // HLT/ConfigTPC/TPCCATracker
  makeComponentConfigurationObject("HLT/ConfigTPC/TPCCFComparison", "", "", 0, 999999999);                      // HLT/ConfigTPC/TPCCFComparison
  makeComponentConfigurationObject("HLT/ConfigTPC/TPCClusterFinderDecoder", "", "", 0, 999999999);              // HLT/ConfigTPC/TPCClusterFinderDecoder
  makeComponentConfigurationObject("HLT/ConfigTPC/TPCClusterFinderPacked", "", "", 0, 999999999);               // HLT/ConfigTPC/TPCClusterFinderPacked
  makeComponentConfigurationObject("HLT/ConfigTPC/TPCDataCompressor", "-deflater-mode 2 -mode 1", "", 0, 999999999);  // HLT/ConfigTPC/TPCDataCompressor
  makeComponentConfigurationObject("HLT/ConfigTPC/TPCdEdxMonitoring", "", "", 0, 999999999);                    // HLT/ConfigTPC/TPCdEdxMonitoring
  makeComponentConfigurationObject("HLT/ConfigTPC/TPCGlobalMerger", "", "", 0, 999999999);                      // HLT/ConfigTPC/TPCGlobalMerger
  makeComponentConfigurationObject("HLT/ConfigTPC/TPCHWClusterTransform", "-charge-threshold 0", "", 0, 999999999);   // HLT/ConfigTPC/TPCHWClusterTransform
  makeComponentConfigurationObject("HLT/ConfigTPC/TPCSliceTracker", "", "", 0, 999999999);                      // HLT/ConfigTPC/TPCSliceTracker
  makeComponentConfigurationObject("HLT/ConfigTPC/TPCTrackHisto", "-event-modulo 20 -buffer-size 5000", "", 0, 999999999);  // HLT/ConfigTPC/TPCTrackHisto
  makeComponentConfigurationObject("HLT/ConfigTPC/TPCVertexFinder", "", "", 0, 999999999);                      // HLT/ConfigTPC/TPCVertexFinder
  makeComponentConfigurationObject("HLT/ConfigTRD/ClusterizerComponent", "output_percentage 100 -lowflux -experiment -tailcancellation -faststreamer -yPosMethod LUT", "", 0, 999999999); // HLT/ConfigTRD/ClusterizerComponent
  makeTriggerConfigurationObject("H-Barrel_pT_Single-V0001.001"); // HLT/ConfigHLT/H_._Barrel_pT_Single_._V0001.001              
  makeTriggerConfigurationObject("H-Barrel_pT_Single-V0002.001"); // HLT/ConfigHLT/H_._Barrel_pT_Single_._V0002.001              
  makeTriggerConfigurationObject("H-Barrel_pT_Single-V0003.001"); // HLT/ConfigHLT/H_._Barrel_pT_Single_._V0003.001              
  makeComponentConfigurationObject("HLT/ConfigTRD/TrackerV1Component", "output_percentage 100 -lowflux -PIDmethod NN", "", 0, 999999999);   // HLT/ConfigTRD/TrackerV1Component
  rootlogon();  // HLT/ConfigMUON/DecisionComponent, FieldIntegrals, HitReconstructor, MansoTrackerFSM, TriggerReconstructor
  CreateDefaultCDBEntries();             
  makeConfigurationObjectTPCHWClusterFinder();               // HLT/ConfigTPC/TPCHWClusterFinder                             
  makeConfigurationObjectVZEROReconstruction();              // HLT/ConfigVZERO/VZEROReconstruction                          
  makeConfigurationObjectZDCReconstruction();                // HLT/ConfigZDC/ZDCESDReco                                     
  //the following for HLT/ConfigHLT/GlobalHisto 
  TString s("-max-track-count 8000 ");
  s+="-histogram TrackPt(100,0,20) -size 1000 -expression Track_pt -title p_{T}_[GeV/c] -cut Track_TPCclus>0 ";
  s+="-histogram TrackPhi(180,0,360) -size 1000 -expression Track_phi -title #phi_(deg) -cut Track_TPCclus>0 ";
  s+="-histogram TrackMultiplicity(250,0,5000) -size 1000 -expression trackcount -title TrackMultiplicity ";
  s+="-histogram TrackEta(100,-2,2) -size 1000 -expression Track_eta -title #eta -cut Track_TPCclus>0&&Track_DCAr<7&&Track_DCAr>-7&&Track_pt>0.3&&Track_eta<0.9&&Track_eta>-0.9 ";
  s+="-histogram TrackTPCclus(200,0,200) -size 1000 -expression Track_TPCclus -title TPC_clusters/track -cut Track_TPCclus>0&&Track_DCAr<7&&Track_DCAr>-7&&Track_pt>0.3&&Track_eta<0.9&&Track_eta>-0.9 ";
  s+="-histogram TrackITSclus(7,0,7) -size 1000 -expression Track_ITSclus -title ITS_clusters/track ";
  s+="-histogram TrackTheta(90,0,180) -size 1000 -expression Track_theta -title #theta_(deg) -cut Track_TPCclus>0&&Track_DCAr<7&&Track_DCAr>-7&&Track_pt>0.3&&Track_eta<0.9&&Track_eta>-0.9 ";
  s+="-histogram TrackDCAr(100,-50,50) -size 1000 -expression Track_DCAr -title DCAr_[cm] -cut Track_TPCclus>0 ";
  s+="-histogram TrackCharge -size 1000 -expression Track_charge -title Polarity -cut Track_TPCclus>0 ";
  s+="-histogram VertexXY -size 1000 -expression vertexY:vertexX -title y:x_[cm] -cut nContributors>3 -opt colz ";
  s+="-histogram VertexX(50,-5,5)  -size 1000 -expression vertexX -title x_[cm] -cut nContributors>3 ";
  s+="-histogram VertexY(50,-5,5)  -size 1000 -expression vertexY -title y_[cm] -cut nContributors>3 ";
  s+="-histogram VertexZ(200,-50,50)  -size 1000 -expression vertexZ -title x_[cm] -cut nContributors>3 ";
  s+="-histogram VertexTrendX -size 1000 -expression vertexX:event -title vertexX_vs_event -cut nContributors>3 ";
  s+="-histogram VertexTrendY -size 1000 -expression vertexY:event -title vertexY_vs_event -cut nContributors>3 ";
  makeComponentConfigurationObject("HLT/ConfigHLT/GlobalHisto", s.Data(), "", 0, 999999999);
}

//HLT/Calib/esdLayout                                         HLT/exa/sampleAliHLTOUTHandlerEsdBranch.C ?
//HLT/Calib/RecoParam
//HLT/Calib/StreamerInfo
//HLT/ConfigHLT/esdLayout
//HLT/ConfigHLT/HLTGlobalTrigger
//HLT/ConfigTPC/TPCClusterFinder32Bit
//HLT/ConfigTPC/TPCClusterFinderUnpacked
//HLT/ConfigTPC/TPCDataCompressorHuffmanTables
//can be used to modify an existing object: HLT/programs/adjustOCDBObject.C
