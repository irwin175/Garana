R__ADD_LIBRARY_PATH(/home/burke/Research/garana_work/garana_build/lib)
R__LOAD_LIBRARY(libGaranaAccessors.so)
R__LOAD_LIBRARY(libGaranaProducts.so)
R__LOAD_LIBRARY(libGaranaDict.so)
R__LOAD_LIBRARY(libGaranaUtility.so)
R__ADD_INCLUDE_PATH(/home/burke/Research/garana_work/garana)


/* Macro: FS Proton Reco Efficiency
 * Purpose: calculate the ratio of true final state protons that are matched with at least one track to
 * final state protons that are identified using Backtracker on tracks in the reconstruction tree*/

#include "garana/Accessors/TreeManager.h"
#include "garana/Utility/Backtracker.h"
#include <TH1F.h>
#include <iostream>
#include <math.h>

using namespace garana;


void FirstGaranaMacro()
{

std::cout << "loaded libraries" << std::endl;

string infile = "structuredtree.root";

TreeManager* treeman = new TreeManager(infile);
Backtracker* bt = new Backtracker(treeman);


GenTree* gen = treeman->GetGenTree();
G4Tree*  g4  = treeman->GetG4Tree();
RecoTree* reco = treeman->GetRecoTree();

std::cout<<"Found trees with "<<treeman->NEntries()<<" entries"<<std::endl;

UInt_t TruthCount = 0;
UInt_t RecoCount = 0;
const vector<UInt_t>* particleID = 0;
double Energy[20];
double Efficiency[20];
for(UInt_t ei=0; ei<20; ei++)
{
	for(UInt_t ientry=0; ientry<treeman->NEntries(); ientry++) //Loops over all of the events in the .root file
	{
		treeman->GetEntry(ientry);
		bt->FillMaps();

		for(UInt_t ig4p=0; ig4p<g4->NSim(); ig4p++)	//Finds the number of True FS Protons that are matched to at least
								//one reconstructed track
		{
			if(g4->IsPrimary(ig4p) && g4->PDG(ig4p) == 2212 && bt->G4ParticleToTracks(ig4p)->size()>0 && g4->SimMomEnter(ig4p)->E()>.938)
			{
				std::cout<<g4->SimMomEnter(ig4p)->E()<<std::endl;
				TruthCount += 1;
			}
		}

		for(UInt_t itrk=0; itrk<reco->NTrack(); itrk++)	//Find the number of tracks associated with Final State protons
		{
			particleID = bt->TrackToG4Particles(itrk);
			for(UInt_t bti=0; bti<particleID->size(); bti++)
			{
				if(g4->IsPrimary(particleID->at(bti)) && g4->PDG(particleID->at(bti)) == 2212)
				{
					RecoCount += 1;
				}
			}
		}
		e->Fill(
	}

	Efficiency[ei] = (double)TruthCount/(double)RecoCount;
}
std::cout<<"Truth Count: "<<TruthCount<<std::endl;
std::cout<<"Reconstructed Count: "<<RecoCount<<std::endl;
std::cout<<"Reconstruction Efficiency: "<<(double)TruthCount/(double)RecoCount<<std::endl;


TCanvas* c1 = new TCanvas("Efficiency","Efficiency",600,600);
TGraph* e = new TGraph(20,Energy,Efficiency);
}

vector<double> range(double min, double max, size_t N) 
{
    vector<double> range;
    double delta = (max-min)/double(N-1);
    for(int i=0; i<N; i++) {
        range.push_back(min + i*delta);
    }
    return range;
}
