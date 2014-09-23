#include "Common/PileupWeight.h"

#include <stdexcept>

PileupWeight::PileupWeight(const char* data_era)
{
	if(strcmp(data_era, "fall11") == 0 || strncmp(data_era, "fall11-", 7) == 0)
	{
		pileupDataFile = new TFile(("/scratch/hh/dust/naf/cms/user/aburgmei/tautau-fullyhadronic/WH_Analysis/WH_Analysis/Lumi/pileup/Pileup_2011.root"), "READ");
		pileupMCFile = new TFile(("/scratch/hh/dust/naf/cms/user/aburgmei/tautau-fullyhadronic/WH_Analysis/WH_Analysis/Lumi/pileup/Pileup_Fall11.root"), "READ");
	}
	else if(strcmp(data_era, "summer12") == 0 || strncmp(data_era, "summer12-", 9) == 0)
	{
		pileupDataFile = new TFile(("/scratch/hh/dust/naf/cms/user/aburgmei/tautau-fullyhadronic/WH_Analysis/WH_Analysis/Lumi/pileup/Pileup_2012.root"), "READ");
		pileupMCFile = new TFile(("/scratch/hh/dust/naf/cms/user/aburgmei/tautau-fullyhadronic/WH_Analysis/WH_Analysis/Lumi/pileup/Pileup_Summer12_53x.root"), "READ");
	}
	else if(strcmp(data_era, "2011") == 0 || strcmp(data_era, "2012") == 0)
	{
		pileupDataFile = NULL;
		pileupMCFile = NULL;
	}
	else
		throw std::runtime_error("No such data era: " + std::string(data_era));

	if(pileupDataFile != NULL && pileupMCFile != NULL)
	{
		if(pileupDataFile->IsZombie() || pileupMCFile->IsZombie())
			throw std::runtime_error("At least one Pileup file does not exist");
		pileupData = dynamic_cast<TH1D*>(pileupDataFile->Get("pileup"));
		pileupMC = dynamic_cast<TH1D*>(pileupMCFile->Get("pileup"));
		if(!pileupData || !pileupMC)
			throw std::runtime_error("Pileup files do not contain pileup histograms");

		pileupWeight = dynamic_cast<TH1D*>(pileupData->Clone("pileupWeight"));
		pileupWeight->Scale(1.0/pileupWeight->Integral());

		TH1D* temp = dynamic_cast<TH1D*>(pileupMC->Clone("temp"));
		temp->Scale(1.0/temp->Integral());

		pileupWeight->Divide(temp);
		delete temp;
	}
	else
	{
		pileupData = NULL;
		pileupMC = NULL;
		pileupWeight = NULL;
	}
}

float PileupWeight::getWeight(float numTrueInteractions) const
{
	if(!pileupWeight) return 1.0f;
	return pileupWeight->GetBinContent(pileupWeight->FindBin(numTrueInteractions));
}
