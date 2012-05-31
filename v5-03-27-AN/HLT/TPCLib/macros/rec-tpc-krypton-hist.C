// $Id$
/**
This is only temporary and will be nicer for the user, I just dont have time just now.

To run this: aliroot -b -q rec-tpc-krypton-hist.C
You will have to change the directory for the two clusterdumps. and it have to exist!

You will have to start it in the directory where the raw folders are, and you will also need the Kinematics.root and galice.root at the moment.

It does:
AliRawReaderPublisher
KryptonClusterFinder1 (input: altro data fra AliRawReaderPublisher)
ClusterDump(input: KryptonClusterFinder1)
AltroSelectionComponent(input: altro data fra AliRawReaderPublisher og liste med hardware addresser fra KryptonClusterFinderen )
KryptonClusterFinder2(input: altrodata fra AltroSelectionComponent)
ClusterDump(input: KryptonClusterFinder2)
ClusterHisto(input: KryptonClusterFinder)
HistogramHandler(input: ClusterHisto) 


Sorry that it is not nicer, it will be soon.

Kenneth

 */
{
  AliReconstruction rec;
  //rec.SetInput("zeroSup/");
  rec.SetInput("raw0/");
  //rec.SetInput("../noiseComp/file.root");
  //rec.SetInput("./");
  //rec.SetInput("/scratch/peter/07000015067011.10.root");
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunLocalReconstruction("HLT");
  rec.SetRunReconstruction("");
  rec.SetRunTracking("");
  rec.SetLoadAlignFromCDB(0);
  rec.SetFillTriggerESD(kFALSE);
  rec.SetFillESD("HLT");
  rec.SetOption("HLT", "libAliHLTTPC.so  loglevel=0x7c config=$ALICE_ROOT/HLT/TPCLib/macros/conf-tpc-krypton-hist.C chains=sink1");
  rec.Run();
}
