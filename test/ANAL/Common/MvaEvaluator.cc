#include <fstream>
#include "Common/MvaEvaluator.h"

MvaEvaluator::MvaEvaluator(const char* training, const char* method, const char* variables):
  m_method(method)
{
  // If the MVA couldn't be booked than there is no trained MVA yet. Don't bail though --
  // we'll keep running and create trees to be used to train the MVA. In a second run we
  // can then apply the MVA.
  reader = new TMVA::Reader("Silent");
  
  // Register the variables with the reader
  reader->AddVariable("LeadPt", &mvaVars.leadPt);
  reader->AddVariable("SubPt", &mvaVars.subPt);
  reader->AddVariable("Met", &mvaVars.Met);
  if(strcmp(training, "v3-std") == 0)
    {
      reader->AddVariable("abs(DeltaPhiDiTau)", &absDeltaPhiDiTau);
      reader->AddVariable("TMath::Sqrt(DeltaRDiTau*DeltaRDiTau - DeltaPhiDiTau*DeltaPhiDiTau)", &absDeltaEtaDiTau);
    }
  else
    {
      reader->AddVariable("DeltaRDiTau", &mvaVars.deltaRDiTau);
    }
  reader->AddVariable("PtRatio", &mvaVars.ptRatio);
  if(strcmp(variables, "sane-plus-mt") == 0) reader->AddVariable("TransMassLepMET", &mvaVars.transMassLepMET);
  
  try
    {
      const std::string weightsFile = "/nfs/dust/cms/user/aburgm/tautau-fullyhadronic/WH_Analysis/WH_Analysis/MVA/weights-" + std::string(training) + "-" + std::string(variables) + "-2011-2012-data-tau12fake/WH_" + std::string(method) + ".weights.xml";
      std::ifstream stream(weightsFile.c_str());
      if(!stream)
	throw std::runtime_error("Weight file does not exist");
      
      // Note BookMVA crashes the program if it fails...
      reader->BookMVA(method, weightsFile.c_str());
      std::cout << "Loaded MVA training from " << weightsFile << std::endl;
    }
  catch(...)
    {
      delete reader;
      reader = NULL;
      
      std::cout << "MVA " << method << " not available, could not book weights file" << std::endl;
    }
}

bool MvaEvaluator::isAvailable() const
{
  return reader != NULL;
}

double MvaEvaluator::mvaValue(const MvaTree::Variables::MvaVars& vars) const
{
  mvaVars = vars;
  absDeltaPhiDiTau = fabs(vars.deltaPhiDiTau);
  absDeltaEtaDiTau = TMath::Sqrt(vars.deltaRDiTau*vars.deltaRDiTau - vars.deltaPhiDiTau*vars.deltaPhiDiTau);
  const double value = reader->EvaluateMVA(m_method.c_str());
  return value;
}
