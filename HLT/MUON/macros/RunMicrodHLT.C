////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

/* Sample macro indicating how to run the micro dHLT version.
 */

void RunMicrodHLT(Bool_t rundebug = kFALSE)
{
	if (rundebug)
		gSystem->Load("../lib/Linux-debug/libMUONHLT.so");
	else
		gSystem->Load("../lib/Linux/libMUONHLT.so");
	
	cout << "dHLT Version: " << AliMUONHLT::Version() << endl;

	AliMUONDataInterface data;
	data.SetFile();
	
	// Load the trigger data.
	AliMUONHLT::TriggerSource ts;
	ts.DataToUse(AliMUONHLT::TriggerSource::FromHits);
	ts.FillFrom(&data);

	// Load the cluster data.
	AliMUONHLT::ClusterSource cs;
	cs.DataToUse(AliMUONHLT::ClusterSource::FromHits);
	cs.FillFrom(&data);

	// Create a Track output object.
	AliMUONHLT::TrackSink trackout;

	// Create the dHLT process object and connect it to input and output.
	AliMUONHLT::MicrodHLT dhlt;
	dhlt.SetTriggerSource(&ts);
	dhlt.SetClusterSource(&cs);
	dhlt.SetTrackSink(&trackout);
	
	// run the dHLT
	cout << "running..." << endl;
	dhlt.Run();
	
	// Display the output:
	if ( ! trackout.GetFirstEvent())
		cout << "No tracks found." << endl;

	while (trackout.MoreEvents())
	{
		cout << "================= Event: " << trackout.CurrentEvent() << " =================" << endl;
		Int_t blocknum = 0;

		if ( ! trackout.GetFirstBlock() )
			cout << "No track blocks found." << endl;
		
		while (trackout.MoreBlocks())
		{
			cout << "Block: " << blocknum++ << endl;
			cout << "\tReal\tReal\tHLT\tHLT\tHLT\tL0\tL0" << endl;
			cout << "\tP\tPt\tSign\tP\tPt\tSign\tPt" << endl;

			if ( trackout.GetFirstTrack() == NULL )
				cout << "\tNo tracks found." << endl;

			while (trackout.MoreTracks())
			{
				AliMUONHLT::Track* track = trackout.GetTrack();

				// Find the corresponding trigger record.
				AliMUONHLT::TriggerRecord* trigrec = NULL;
				ts.GetEvent(trackout.CurrentEvent());
				
				if (ts.GetFirstBlock())
				{
					do
					{
						trigrec = ts.GetTrigger( track->TriggerID() );
						if (trigrec != NULL) break;
					} while (ts.GetNextBlock());
				};
				if (trigrec == NULL)
				{
					cerr << "Error: could not find corresponding trigger record!" << endl;
					continue;
				};

				// Find the corresponding particle in Kine.
				// Note: track->fTriggerID happens to be the TreeH track number when
				// filling from Hit objects.
				data.GetEvent(trackout.CurrentEvent());
				AliMUONHit* hit = data.Hit(track->TriggerID(), 0);
				if (hit == NULL)
				{
					cerr << "Warning: Could not find hit 0 for track number = "
					<< track->TriggerID() << endl;
					continue;
				};

				// Need to decode the particle index.
				Int_t num_primaries = data.NumberOfTracks();
				Int_t num_tracks = data.NumberOfParticles();
				Int_t particleindex;
				if (hit->Track() < num_primaries)
					particleindex = hit->Track() + num_tracks - num_primaries;
				else
					particleindex =  hit->Track()- num_primaries;

				TParticle* p = data.Particle(particleindex);

				// All data found, now write output:
				cout << "\t" << p->P()
					<< "\t" << p->Pt()
					<< "\t" << track->ParticleSign()
					<< "\t" << track->P()
					<< "\t" << track->Pt()
					<< "\t" << trigrec->ParticleSign()
					<< "\t" << trigrec->Pt()
					<< endl;

				trackout.GetNextTrack();
			};
			trackout.GetNextBlock();
		};
		trackout.GetNextEvent();
	};
};
