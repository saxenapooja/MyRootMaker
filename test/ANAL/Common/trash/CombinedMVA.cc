
#include <fstream>

#include "Common/CombinedMVA.h"

CombinedMVA::CombinedMVA()
{
  // If the MVA couldn't be booked than there is no trained MVA yet. Don't bail though --
  // we'll keep running and create trees to be used to train the MVA. In a second run we
  // can then apply the MVA.
  reader = new TMVA::Reader("Silent");

  // Register the variables with the reader
  reader->AddVariable("BDT8", &BDT8Output);
  reader->AddVariable("MLP", &MLPOutput);
 
  try
  {
    const std::string weightsFile = "/nfs/dust/cms/user/aburgm/tautau-fullyhadronic/WH_Analysis/WH_Analysis/MVA/claudio/WH_MLPBFGS.weights.xml";
    std::ifstream stream(weightsFile.c_str());
    if(!stream)
      throw std::runtime_error("Weight file does not exist");

    // Note BookMVA crashes the program if it fails...
    reader->BookMVA("MLPBFGS", weightsFile.c_str());
    std::cout << "Loaded combined MVA training from " << weightsFile << std::endl;
  }
  catch(...)
  {
    delete reader;
    reader = NULL;

    std::cout << "Combined MVA not available, could not book weights file" << std::endl;
  }
}

double CombinedMVA::mvaValue(float BDT8Output, float MLPOutput) const
{
  if(!reader) return -1.;

  this->BDT8Output = BDT8Output;
  this->MLPOutput = MLPOutput;

  const double value = reader->EvaluateMVA("MLPBFGS");
  return value;
}
