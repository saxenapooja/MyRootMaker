#ifndef _WHANALYSIS_CUTFLOW_H_
#define _WHANALYSIS_CUTFLOW_H_

#include <TH1D.h>

//cutflow 
class Cutflow
{
 public:

  class Cut
  {
  public:
    Cut(unsigned int index);
    const unsigned int index;
    bool pass;
  };
  
  Cutflow(const char* name);
  ~Cutflow();
  
  const Cut* AddCut(const char* name);
  void Pass(const Cut* cut);
  void FinishEvent(double weight);
  
  const std::string& GetName() const { return name_;    }
  TH1D* GetCutflowHistogram()  const { return hCutflow; }
  
  void Print() const;
  void Print(std::ostream& out) const;
  
 private:
  const std::string name_;
  std::vector<Cut*> cuts;
  TH1D* hCutflow;
};

// Cutflow by region
class CutflowRegions
{
 public:
  class Cut // cut for all regions
  {
  public:
    std::vector<const Cutflow::Cut*> cuts;
  };
  
  CutflowRegions();
  ~CutflowRegions();

  int AddRegion(const char* name);
  int AllRegions() const;

  const Cut* AddCut(const char* name);
  void Pass(const Cut* cut, int regions);
  void FinishEvent();

  int GetCutflowMask(const char* region) const;
  const Cutflow* GetCutflow(const char* region) const;

  TH1D* GetCutflowHistogram(const char* region) const { 
    return GetCutflow(region)->GetCutflowHistogram(); 
  }
  
 private:
  std::vector<Cutflow*> cutflows;
  std::vector<Cut*> cuts;
};

#endif //_WHANALYSIS_CUTFLOW_H_
