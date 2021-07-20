R__ADD_LIBRARY_PATH(/home/burke/Research/garana_work/garana_build/lib)
R__LOAD_LIBRARY(libGaranaAccessors.so)
R__LOAD_LIBRARY(libGaranaProducts.so)
R__LOAD_LIBRARY(libGaranaDict.so)
R__LOAD_LIBRARY(libGaranaUtility.so)
R__ADD_INCLUDE_PATH(/home/burke/Research/garana_work/garana)


/* Macro: FS Proton Reco Efficiency
 * Purpose 1: calculate the ratio of true final state protons that are matched with at least one track to
 * final state protons that are identified using Backtracker on tracks in the reconstruction tree (function of momentum)
 * Purpose 2: same calculation but as a function of angle (angle from drift direction, angle around drift direction)*/

#include "garana/Accessors/TreeManager.h"
#include "garana/Utility/Backtracker.h"
#include <TH1F.h>
#include <iostream>
#include <math.h>
#include <TStyle.h>
#include "functions.h"

using namespace garana;

void Eff()
{
string infile = "structuredtree.root";

TreeManager* treeman = new TreeManager(infile);
Backtracker* bt = new Backtracker(treeman);
GenTree* gen = treeman->GetGenTree();
G4Tree*  g4  = treeman->GetG4Tree();
RecoTree* reco = treeman->GetRecoTree();

vector<int> particle_types= {2212,13,211,321};
double binspace [54] = {0.,.02,.04,.06,.08,.1,.12,.14,.16,.18,.2,.22,.24,.26,.28,.3,.32,.34,.36,.38,.4,.5,.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,4.,5.,6.,7.,8.,9.,10.};


//Declare histograms for all particle types
//Proton
TH1F* ptp = new TH1F("ptp","True Proton Momentum",53,binspace);
TH2F* pta = new TH2F("pta","True Proton Drift Angle",5,0,M_PI,10,-M_PI,M_PI);
TH1F* prp = new TH1F("prp","True Momentum for Reconstructed Protons using new Method",53,binspace);
TH2F* pra = new TH2F("pra","Reconstructed Proton Drift Angles",5,0,M_PI,10,-M_PI,M_PI);

//Muons
TH1F* mtp = new TH1F("mtp","True Muon Momentum",53,binspace);
TH2F* mta = new TH2F("mta","True Muon Drift Angle",5,0,M_PI,10,-M_PI,M_PI);
TH1F* mrp = new TH1F("mrp","True Momentum for Reconstructed Muons",53,binspace);
TH2F* mra = new TH2F("mra","Reconstructed Muon Drift Angles",5,0,M_PI,10,-M_PI,M_PI);

//Charged Pions
TH1F* cptp = new TH1F("cptp","True Charged Pion Momentum",53,binspace);
TH2F* cpta = new TH2F("cpta","True Charged Pion Drift Angle",5,0,M_PI,10,-M_PI,M_PI);
TH1F* cprp = new TH1F("cprp","True Momentum for Reconstructed Charged Pions",53,binspace);
TH2F* cpra = new TH2F("cpra","Reconstructed Charged Pion Drift Angles",5,0,M_PI,10,-M_PI,M_PI);

//Charged Kaons
TH1F* ktp = new TH1F("ktp","True Charged Kaon Momentum",53,binspace);
TH2F* kta = new TH2F("kta","True Charged Kaon Drift Angles",5,0,M_PI,10,-M_PI,M_PI);
TH1F* krp = new TH1F("krp","True Momentum for Reconstructed Charged Kaons",53,binspace);
TH2F* kra = new TH2F("kra","Reconstructed Charged Kaon Drift Angles",5,0,M_PI,10,-M_PI,M_PI);

vector<TH1F*> oneD_histos = {ptp,prp,mtp,mrp,cptp,cprp,ktp,krp};
vector<TH2F*> twoD_histos = {pta,pra,mta,mra,cpta,cpra,kta,kra};

for(UInt_t ipdg = 0; ipdg<particle_types.size(); ipdg++)
{
	for(UInt_t ientry=0; ientry<treeman->NEntries(); ientry++)		//Loop through all events
		{
		treeman->GetEntry(ientry);
		bt->FillMaps();

		for(UInt_t ig4p=0; ig4p<g4->NSim(); ig4p++)			//Loop through true particles
			{
			if(g4->IsPrimary(ig4p) && g4->PDG(ig4p) == particle_types.at(ipdg))
				{
					int association = 0;
					if(-220.<gen->NuVertex(bt->G4ParticleToGTruth(ig4p))->X()<220. && sqrt(pow(gen->NuVertex(bt->G4ParticleToGTruth(ig4p))->Y(),2)+pow(gen->NuVertex(bt->G4ParticleToGTruth(ig4p))->Z(),2))<230.)
					{
						auto truth_angles = Direction(g4,ig4p);
						oneD_histos.at(2*ipdg)->Fill(g4->SimMomEnter(ig4p,0)->P());
						twoD_histos.at(2*ipdg)->Fill(truth_angles.first,truth_angles.second);	
						for(UInt_t itrk=0; itrk<reco->NTrack(); itrk++)	//Loop through reconstructed tracks
						{
							auto PrimeParticle = bt->TrackToG4ParticleDeposit(itrk);
							if(ig4p==PrimeParticle.first)
							{
								association += 1;
							}	
						}
					}
				if(association != 0)
				{
					oneD_histos.at(2*ipdg+1)->Fill(g4->SimMomEnter(ig4p,0)->P());
					auto reco_angles = Direction(g4,ig4p);
                                        twoD_histos.at(2*ipdg+1)->Fill(reco_angles.first,reco_angles.second);
				}
				}
			}	
		}
	TCanvas* c1 = new TCanvas("Momentum","Momentum",1500,750);
	c1->Divide(2,1);

	c1->cd(1);
	oneD_histos.at(2*ipdg)->Sumw2();
	oneD_histos.at(2*ipdg)->GetXaxis()->SetRangeUser(0.,3.);
	oneD_histos.at(2*ipdg)->Draw("E");

	c1->cd(2);
	oneD_histos.at(2*ipdg+1)->Sumw2();
	oneD_histos.at(2*ipdg+1)->GetXaxis()->SetRangeUser(0.,3.);
	oneD_histos.at(2*ipdg+1)->Draw("E");
//	c1->Print("mommu.png");

	TCanvas* c3 = new TCanvas("Drift Angle Efficiency","Drift Angle Efficiency",2000,1000);
	c3->Divide(2,1);

	c3->cd(1);
	twoD_histos.at(2*ipdg)->SetTitle("True Angular Distribution about the Drift Direction");
	twoD_histos.at(2*ipdg)->GetXaxis()->SetTitle("Angle from Drift Direction");
	twoD_histos.at(2*ipdg)->GetYaxis()->SetTitle("Angle about Drift Direction");
	twoD_histos.at(2*ipdg)->SetMarkerStyle(kFullDotMedium);
	twoD_histos.at(2*ipdg)->Draw();

	c3->cd(2);
	twoD_histos.at(2*ipdg+1)->SetTitle("Reco Angular Distribution about the Drift Direction");
	twoD_histos.at(2*ipdg+1)->GetXaxis()->SetTitle("Angle from Drift Direction");
	twoD_histos.at(2*ipdg+1)->GetYaxis()->SetTitle("Angle about Drift Direction");
	twoD_histos.at(2*ipdg+1)->SetMarkerStyle(kFullDotMedium);
	twoD_histos.at(2*ipdg+1)->Draw();
//	c3->Print("angmu.png");
}

//Efficiency Plots
//Protons
TCanvas* cp1 = new TCanvas("Proton Efficiency vs. Momentum","Efficiency vs. Momentum",2000,1000);
cp1->Divide(2,1);

cp1->cd(1);
TH1D* epFull =(TH1D*) oneD_histos.at(1)->Clone("epFull");
epFull->Divide(oneD_histos.at(1),oneD_histos.at(0));
TH1D* ep400  =(TH1D*) epFull->Clone("ep400");
epFull->SetTitle("FS Proton Reco Efficiency vs. Momentum");
epFull->GetXaxis()->SetRangeUser(0.,3.);
epFull->GetYaxis()->SetRangeUser(0.,1.5);
epFull->Draw("E");

cp1->cd(2);
ep400->SetTitle("Low Momentum Efficiency");
ep400->GetXaxis()->SetRangeUser(0.,.4);
ep400->GetYaxis()->SetRangeUser(0.,1.0);
ep400->Draw("E");

cp1->Print("EffwFVcut.png");

TCanvas* cp2 = new TCanvas("Garana Proton Reco Efficiency vs. Drift Angle","Garana Reco Efficiency vs. Drift Angle",1000,1000);
TH2D* epa =(TH2D*) twoD_histos.at(1)->Clone("epa");
epa->Divide(twoD_histos.at(1),twoD_histos.at(0));
epa->SetTitle("Proton Reco Efficiency vs. Drift Angle");
epa->GetXaxis()->SetTitle("Angle from Drift Direction");
epa->GetYaxis()->SetTitle("Angle about Drift Direction");
epa->GetZaxis()->SetRangeUser(0.,1.);
epa->Print("DriftAngleEff.png");
epa->Draw("colz");

//cp2->Print("angeffpro.png");

//Muons
TCanvas* cm1 = new TCanvas("Muon Efficiency vs. Momentum","Efficiency vs. Momentum",2000,1000);
cm1->Divide(2,1);

cm1->cd(1);
TH1D* emFull =(TH1D*) oneD_histos.at(3)->Clone("emFull");
emFull->Divide(oneD_histos.at(3),oneD_histos.at(2));
TH1D* em400  =(TH1D*) emFull->Clone("em400");
emFull->SetTitle("FS Muon Reco Efficiency vs. Momentum");
emFull->GetXaxis()->SetRangeUser(0.,3.);
emFull->GetYaxis()->SetRangeUser(0.,1.5);
emFull->Draw("E");

cm1->cd(2);
em400->SetTitle("Low Momentum Efficiency");
em400->GetXaxis()->SetRangeUser(0.,.4);
em400->GetYaxis()->SetRangeUser(0.,1.0);
em400->Draw("E");
//cm1->Print("momeffmu.png");

TCanvas* cm2 = new TCanvas("Garana Muon Reco Efficiency vs. Drift Angle","Garana Reco Efficiency vs. Drift Angle",1000,1000);
TH2D* ema =(TH2D*) twoD_histos.at(3)->Clone("ema");
ema->Divide(twoD_histos.at(3),twoD_histos.at(2));
ema->SetTitle("Muon Reco Efficiency vs. Drift Angle");
ema->GetXaxis()->SetTitle("Angle from Drift Direction");
ema->GetYaxis()->SetTitle("Angle about Drift Direction");
ema->GetZaxis()->SetRangeUser(0.,1.);
ema->Print("DriftAngleEff.png");
ema->Draw("colz");
//cm2->Print("angeffmu.png");


//Charged Pions
TCanvas* ccp1 = new TCanvas("Charged Pion Efficiency vs. Momentum","Efficiency vs. Momentum",2000,1000);
ccp1->Divide(2,1);

ccp1->cd(1);
TH1D* ecpFull =(TH1D*) oneD_histos.at(5)->Clone("ecpFull");
ecpFull->Divide(oneD_histos.at(5),oneD_histos.at(4));
TH1D* ecp400  =(TH1D*) ecpFull->Clone("ecp400");
ecpFull->SetTitle("FS Charged Pion Reco Efficiency vs. Momentum");
ecpFull->GetXaxis()->SetRangeUser(0.,3.);
ecpFull->GetYaxis()->SetRangeUser(0.,1.5);
ecpFull->Draw("E");

ccp1->cd(2);
ecp400->SetTitle("Low Momentum Efficiency");
ecp400->GetXaxis()->SetRangeUser(0.,.4);
ecp400->GetYaxis()->SetRangeUser(0.,1.0);
ecp400->Draw("E");
//ccp1->Print("momeffcp.png");

TCanvas* ccp2 = new TCanvas("Garana Charged Pion Reco Efficiency vs. Drift Angle","Garana Reco Efficiency vs. Drift Angle",1000,1000);
TH2D* ecpa =(TH2D*) twoD_histos.at(5)->Clone("ecpa");
ecpa->Divide(twoD_histos.at(5),twoD_histos.at(4));
ecpa->SetTitle("Charged Pion Reco Efficiency vs. Drift Angle");
ecpa->GetXaxis()->SetTitle("Angle from Drift Direction");
ecpa->GetYaxis()->SetTitle("Angle about Drift Direction");
ecpa->GetZaxis()->SetRangeUser(0.,1.);
ecpa->Print("DriftAngleEff.png");
ecpa->Draw("colz");
//ccp2->Print("angeffcp.png");


//Charged Kaons
TCanvas* ck1 = new TCanvas("Charged Kaon Efficiency vs. Momentum","Charged Kaon Efficiency vs. Momentum",2000,1000);
ck1->Divide(2,1);

ck1->cd(1);
TH1D* ekFull =(TH1D*) oneD_histos.at(7)->Clone("ekFull");
ekFull->Divide(oneD_histos.at(7),oneD_histos.at(6));
TH1D* ek400  =(TH1D*) ekFull->Clone("ek400");
ekFull->SetTitle("FS Charged Kaon Reco Efficiency vs. Momentum");
ekFull->GetXaxis()->SetRangeUser(0.,3.);
ekFull->GetYaxis()->SetRangeUser(0.,1.5);
ekFull->Draw("E");

ck1->cd(2);
ek400->SetTitle("Low Momentum Efficiency");
ek400->GetXaxis()->SetRangeUser(0.,.4);
ek400->GetYaxis()->SetRangeUser(0.,1.0);
ek400->Draw("E");
//ck1->Print("momeffk.png");

TCanvas* ck2 = new TCanvas("Garana Charged Kaon Reco Efficiency vs. Drift Angle","Garana Charged Kaon Reco Efficiency vs. Drift Angle",1000,1000);
TH2D* eka =(TH2D*) twoD_histos.at(7)->Clone("eka");
eka->Divide(twoD_histos.at(7),twoD_histos.at(6));
eka->SetTitle("Charged Kaon Reco Efficiency vs. Drift Angle");
eka->GetXaxis()->SetTitle("Angle from Drift Direction");
eka->GetYaxis()->SetTitle("Angle about Drift Direction");
eka->GetZaxis()->SetRangeUser(0.,1.);
eka->Print("DriftAngleEff.png");
eka->Draw("colz");
//ck2->Print("angeffk.png");

}
