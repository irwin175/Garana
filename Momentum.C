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

void Momentum()
{

std::cout << "loaded libraries" << std::endl;

string infile = "structuredtree.root";

TreeManager* treeman = new TreeManager(infile);
Backtracker* bt = new Backtracker(treeman);


GenTree* gen = treeman->GetGenTree();
G4Tree*  g4  = treeman->GetG4Tree();
RecoTree* reco = treeman->GetRecoTree();

std::cout<<"Found trees with "<<treeman->NEntries()<<" entries"<<std::endl;

TH1D* t = new TH1D("t","True Protons",100,0,MaxP(g4,treeman));
TH1D* r = new TH1D("r","Reconstructed Protons",100,0,MaxP(g4,treeman));

const vector<UInt_t>* particleID = 0;

for(UInt_t ientry=0; ientry<treeman->NEntries(); ientry++)
	{
	treeman->GetEntry(ientry);
	bt->FillMaps();

	for(UInt_t ig4p=0; ig4p<g4->NSim(); ig4p++)
		{
		if(g4->IsPrimary(ig4p) && g4->PDG(ig4p) == 2212)
			{
			t->Fill(g4->SimMomEnter(ig4p,0)->P());
			}
		}
	for(UInt_t itrk=0; itrk<reco->NTrack(); itrk++)
		{
		particleID = bt->TrackToG4Particles(itrk);
		for(UInt_t bti=0; bti<particleID->size(); bti++)
			{
			if(g4->IsPrimary(particleID->at(bti)) && g4->PDG(particleID->at(bti)) == 2212)
				{
				r->Fill(g4->SimMomEnter(particleID->at(bti),0)->P());
				}
			}
		}
	}

TCanvas* c1 = new TCanvas("Momentum","Momentum",1000,500);
c1->Divide(2,1);

c1->cd(1);
t->Draw("E");

c1->cd(2);
r->Draw("E");

TCanvas* c2 = new TCanvas("Efficiency","Efficiency",1000,1000);
TH1D* e =(TH1D*) r->Clone("e");

e->SetTitle("Garana Reconstruction Efficiency; Momentum (GeV); Efficiency");

e->Divide(r,t);
e->Draw("E");
}
