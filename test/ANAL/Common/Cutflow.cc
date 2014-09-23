#include "Common/Cutflow.h"
#include <iostream>
#include <stdexcept>
#include <cassert>

Cutflow::Cut::Cut(unsigned int index): index(index), pass(false)
{
}



Cutflow::Cutflow(const char* name): name_(name)
{
  hCutflow = new TH1D(("h" + name_).c_str(), "Cutflow;;N_{events}", 100, -0.5, 99.5);
  hCutflow->GetXaxis()->SetBinLabel(1, "Pre-Skim"); // First bin, to be set by analysis
}


Cutflow::~Cutflow()
{
  for(unsigned int i = 0; i < cuts.size(); ++i)  delete cuts[i];
}


const Cutflow::Cut* Cutflow::AddCut(const char* name)
{
  cuts.push_back(new Cut(cuts.size()));
  const Cut* cut = cuts.back();
  hCutflow->GetXaxis()->SetBinLabel(cut->index + 2, name);
  //  std::cout<<"@@@@@@ cut->index :"<< name<< ", "<< cut->index << std::endl;
  return cut;
}


void Cutflow::Pass(const Cut* cut)
{
  // std::cout<<" cut->index  :"<< cut->index  <<" size  :"<<cuts.size() << std::endl;
  assert( cut->index < cuts.size());
  assert( cut->index == 0 || cuts[cut->index - 1]->pass ) ;
  
  Cut* nonconst_cut = cuts[cut->index];
  nonconst_cut->pass = true;
}


void Cutflow::FinishEvent(double weight)
{
  for(unsigned int i = 0; i < cuts.size(); ++i)
    {
      if(cuts[i]->pass)
	{
	  hCutflow->SetBinContent(i+2, hCutflow->GetBinContent(i+2) + weight);
	  cuts[i]->pass = false;
	}
    }
}


void Cutflow::Print(std::ostream& out) const
{
  unsigned int i;
  for(i = 0; i < cuts.size(); ++i)
    {
      if(!cuts[i]->pass)
	{
	  out << "Cutflow fails at cut " << i << ": " << hCutflow->GetXaxis()->GetBinLabel(i+2) << std::endl;
	  break;
	}
    }
  
  if(i == cuts.size())  out << "Survives the full cutflow" << std::endl;
}



void Cutflow::Print() const
{
  unsigned int i;
  for(i = 0; i < cuts.size(); ++i)
    {
      if(!cuts[i]->pass)
	{
	  std::cout << "Cutflow fails at cut " << i << ": " << hCutflow->GetXaxis()->GetBinLabel(i+2) << std::endl;
	  break;
	}
    }
  
  if(i == cuts.size())  std::cout << "Survives the full cutflow" << std::endl;
}


CutflowRegions::CutflowRegions()
{
}


CutflowRegions::~CutflowRegions()
{
  for(unsigned int i = 0; i < cutflows.size(); ++i)
    delete cutflows[i];
  for(unsigned int i = 0; i < cuts.size(); ++i)
    delete cuts[i];
}


int CutflowRegions::AddRegion(const char* name)
{
  assert(cuts.empty());
  const int size = cutflows.size();
  cutflows.push_back(new Cutflow(name));
  return 1 << size;
}



int CutflowRegions::AllRegions() const
{
  //  std::cout<<"cutflows.size() :"<< cutflows.size() << std::endl;
  //  std::cout<<"return :"<< (1 << cutflows.size())<< std::endl;
  return (1 << cutflows.size()) - 1;
}


const CutflowRegions::Cut* CutflowRegions::AddCut(const char* name)
{
  Cut* cut = new Cut;
  cut->cuts.resize(cutflows.size());

  for(unsigned int i = 0; i < cutflows.size(); ++i)
    cut->cuts[i] = cutflows[i]->AddCut(name);
  cuts.push_back(cut);
  return cut;
}


void CutflowRegions::Pass(const Cut* cut, int regions)
{
  for(unsigned int i = 0; i < cutflows.size(); ++i)
    {
      if(regions & (1 << i))
	cutflows[i]->Pass(cut->cuts[i]);
    }
}



void CutflowRegions::FinishEvent()
{
  for(unsigned int i = 0; i < cutflows.size(); ++i)
    cutflows[i]->FinishEvent(1.);
}



int CutflowRegions::GetCutflowMask(const char* region) const
{
  for(unsigned int i = 0; i < cutflows.size(); ++i)
    if(cutflows[i]->GetName() == region)
      return 1 << i;
  throw std::runtime_error(std::string("No such region: ") + region);
}



const Cutflow* CutflowRegions::GetCutflow(const char* region) const
{
  for(unsigned int i = 0; i < cutflows.size(); ++i)
    if(cutflows[i]->GetName() == region)
      return cutflows[i];
  throw std::runtime_error(std::string("No such region: ") + region);
}
