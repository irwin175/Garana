R__ADD_LIBRARY_PATH(/home/burke/Research/garana_work/garana_build/lib)
R__LOAD_LIBRARY(libGaranaAccessors.so)
R__LOAD_LIBRARY(libGaranaProducts.so)
R__LOAD_LIBRARY(libGaranaDict.so)
R__LOAD_LIBRARY(libGaranaUtility.so)
R__ADD_INCLUDE_PATH(/home/burke/Research/garana_work/garana)

#include "garana/Accessors/TreeManager.h"
#include "garana/Utility/Backtracker.h"
#include <TH1F.h>
#include <iostream>
#include <math.h>

using namespace garana;

vector<double> linspace(double min, double max, UInt_t size)
{
        vector<double> p;
        double step = (max-min)/size;
        for(UInt_t a=0; a<size; a++)
        {
                p.push_back(min+a*step);
        }
        p.push_back(max);
        return p;
}

double MaxP(G4Tree* g4, TreeManager* treeman)
{
        double maxp = 0;
        double p;
        for(UInt_t b=0; b<treeman->NEntries(); b++) //Loops over all of the events in the .root file
        {
                treeman->GetEntry(b);
                for(UInt_t c=0; c<g4->NSim(); c++)
                {
                        if(g4->IsPrimary(c) && g4->PDG(c) == 2212)
                        {
                                p = g4->SimMomEnter(c,0)->P();
                                if(p>maxp)
                                {
                                        maxp=p;
                                }
                        }
                }
        }
        return maxp;
}

pair<float,float> Direction(G4Tree* g4,int ig4p)
{
        double DriftA = acos(g4->SimMomEnter(ig4p,0)->Px()/g4->SimMomEnter(ig4p,0)->P());
        double DriftAzi = atan2(g4->SimMomEnter(ig4p,0)->Py(),g4->SimMomEnter(ig4p,0)->Pz());
        pair<float,float> direction;
        direction = make_pair(DriftA,DriftAzi);
        return direction;
}

//Developing a way to collect a vector of pairs for track ID and percent Energy donated for FS Protons
pair<int,float> TrackDeposition(G4Tree* g4,int ig4p, RecoTree* reco, Backtracker* bt)
{
        pair<int,float> MaxDeposit;
        double max = 0.0;
        for(UInt_t itrk = 0; itrk<reco->NTrack(); itrk++)  //Loop through all tracks in event
        {
                auto templist = bt->TrackToG4ParticlesDeposits(itrk);
                for(auto iter = templist.begin(); iter != templist.end(); iter++)  //Loop through all responsible
                                                                                   //particles for a track
                {
			if(ig4p==iter->first && iter->second>max)       //associate g4 FS proton with list of tracks and fractional energy
                                                //deposition
                        {
                                MaxDeposit = make_pair(itrk,iter->second/g4->SimMomEnter(ig4p,0)->E());
                        }
                }
        }

        return MaxDeposit;
}

