void
overlaps() {

	gAlice->Init("./Config.C");
	gGeoManager->CheckOverlaps();
	gGeoManager->PrintOverlaps();
}
