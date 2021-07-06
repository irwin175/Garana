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


string infile = "structuredtree.root";
TreeManager* treeman = new TreeManager(infile);
G4Tree*  g4  = treeman->GetG4Tree();

double max(G4Tree* g4, TreeManager* treeman)
{
        double maxp = 0;
	double p;
        for(UInt_t b=0; b<treeman->NEntries(); b++) //Loops over all of the events in the .root file
        {treeman->GetEntry(b);
                for(UInt_t c=0; c<g4->NSim(); c++)
                {if(g4->IsPrimary(c) && g4->PDG(c) == 2212)
                        {p = g4->SimMomEnter(c,0)->P();
                                if(p>maxp)
                                {maxp=p;}}}}
        return maxp;}

void MaxP()
{
	std::cout<<max(g4,treeman)<<std::endl;
}
