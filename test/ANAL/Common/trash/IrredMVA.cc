
#include <fstream>

#include "Common/IrredMVA.h"

IrredMVA::IrredMVA()
{
  // If the MVA couldn't be booked than there is no trained MVA yet. Don't bail though --
  // we'll keep running and create trees to be used to train the MVA. In a second run we
  // can then apply the MVA.
  reader = new TMVA::Reader("Silent");

  // Register the variables with the reader
  reader->AddVariable("muPt", &muPt);
  reader->AddVariable("tau1Pt", &tauOSPt);
  reader->AddVariable("tau2Pt", &tauSSPt);
  reader->AddVariable("met", &met);
 
  try
  {
    const std::string weightsFile = "/nfs/dust/cms/user/aburgm/tautau-fullyhadronic/WH_Analysis/WH_Analysis/MVA/atanu/MLP.weights.xml";
    std::ifstream stream(weightsFile.c_str());
    if(!stream)
      throw std::runtime_error("Weight file does not exist");

    // Note BookMVA crashes the program if it fails...
    reader->BookMVA("MLP", weightsFile.c_str());
    std::cout << "Loaded irreducible MVA training from " << weightsFile << std::endl;
  }
  catch(...)
  {
    delete reader;
    reader = NULL;

    std::cout << "Irreducible MVA not available, could not book weights file" << std::endl;
  }
}

double IrredMVA::mvaValue(float muPt, float tauOSPt, float tauSSPt, float met) const
{
  if(!reader) return -1.;

  this->muPt = muPt;
  this->tauOSPt = tauOSPt;
  this->tauSSPt = tauSSPt;
  this->met = met;

  const double value = reader->EvaluateMVA("MLP");
  return value;
}
