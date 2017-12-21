/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

#ifndef ALIJTRACKCUT_H
#define ALIJTRACKCUT_H

class AliJTrack;

class AliJTrackCut {
    public:
        enum { kJTPCOnly, kJRaa, kJGlobalTightDCA, kJGlobalDCA, kJGlobalSDD , kJHybrid, kJHybridAOD86, kJHybridStep1, kJNTrackCuts };

        AliJTrackCut();
        ~AliJTrackCut(){;}
        bool IsSelected ( AliJTrack *track, int icut ) const;
        bool SetMomentum( AliJTrack *track, int icut ) const;
        int GetNCut() const { return kJNTrackCuts; } //TODO
        static AliJTrackCut& GetSpecialInstance();
        static const AliJTrackCut& GetInstance();
    private:

};

#endif
