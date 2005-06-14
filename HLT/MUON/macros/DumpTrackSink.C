////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

/* Dumps the contents of the TrackSink 'ts' to the screen.
 */
void DumpTrackSink(AliMUONHLT::TrackSink* ts)
{
	cout << "File  : " << ts->FileName() << endl;
	cout << "Folder: " << ts->FolderName() << endl;
	if ( ! ts->GetFirstEvent() )
	{
		cout << "No events found." << endl;
		return;
	};
	
	while (ts->MoreEvents())
	{
		cout << "================= Event: " << ts->CurrentEvent() << " =================" << endl;
		Int_t blocknum = 0;

		if ( ! ts->GetFirstBlock() )
		{
			cout << "No blocks found." << endl;
			ts->GetNextEvent();
			continue;
		};

		while (ts->MoreBlocks())
		{
			cout << "Block: " << blocknum++ << endl;
			if ( ts->GetFirstTrack() == NULL )
			{
				cout << "\tNo tracks found." << endl;
				ts->GetNextBlock();
				continue;
			};
			
			while (ts->MoreTracks())
			{
				const AliMUONHLT::Track* data = ts->GetTrack();

				cout << "\tTrigger ID = " << data->TriggerID()
					<< ", Sign = " << data->ParticleSign()
					<< ", P = " << data->P()
					<< ", Pt = " << data->Pt()
					<< endl;
				cout << "\t\tX\tY\tLeft\tRight\tBottom\tTop" << endl;
				for (Int_t i = 0; i < 10; i++)
				{
					cout << "\t\t" << data->Hit(i).fX
						<< "\t" << data->Hit(i).fY
						<< "\t" << data->RegionOfInterest(i).Left()
						<< "\t" << data->RegionOfInterest(i).Right()
						<< "\t" << data->RegionOfInterest(i).Bottom()
						<< "\t" << data->RegionOfInterest(i).Top()
						<< endl;
				};

				ts->GetNextTrack();
			};
			ts->GetNextBlock();
		};
		ts->GetNextEvent();
	};
};

