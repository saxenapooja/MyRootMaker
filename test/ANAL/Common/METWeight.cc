#include "Common/METWeight.h"

#include <cassert>
#include <stdexcept>

METWeight::METWeight(const char* data_era, const char* channel)
{
  // TODO: Replace this by proper SS weight files (not used...)!
  if(std::string(data_era) != "fall11" && std::string(data_era) != "summer12" && std::string(data_era) != "fall11-1jet" && std::string(data_era) != "fall11-2jets" && std::string(data_era) != "fall11-3jets" && std::string(data_era) != "fall11-4jets" && std::string(data_era) != "summer12-1jet" && std::string(data_era) != "summer12-2jets" && std::string(data_era) != "summer12-3jets" && std::string(data_era) != "summer12-4jets")
  {
    fileMETWeight = new TFile("/nfs/dust/cms/user/aburgm/tautau-fullyhadronic/WH_Analysis/WH_Analysis/CorrFactors/METWeight.root");
    if(fileMETWeight->IsZombie())
      throw std::runtime_error("Failed to load MET weights!");

    fMETWeight = dynamic_cast<TF1*>(fileMETWeight->Get(("fit" + std::string(channel) + std::string(data_era)).c_str()));
    if(!fMETWeight)
      throw std::runtime_error("MET Weight file does not contain function!");
  }
  else
  {
    fMETWeight = NULL;
  }
}

double METWeight::getWeight(double met) const
{
  assert(met >= 0.);
  const double truncatedMet = std::min(met, 100.);

  if(fMETWeight)
    return fMETWeight->Eval(truncatedMet);
  else
    return 1.0;
}
