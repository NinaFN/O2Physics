
#ifndef PWGCF_GENERICFRAMEWORK_CORE_CALCULATOR_H_
#define PWGCF_GENERICFRAMEWORK_CORE_CALCULATOR_H_

#include <string>
#include <functional>
#include <iostream>
#include <memory>
#include <map>
#include <string>
#include "TH1.h"
#include "TString.h"
#include "FlowContainer.h"
#include "FlowPtContainer.h"

class Calculator
{
 public:
  using ObsFunc = TH1D* (Calculator::*)(int);
  Calculator();
  Calculator(FlowContainer* cont1, FlowPtContainer* cont2);
  ~Calculator();
  void setDataContainers(FlowContainer* cont1, FlowPtContainer* cont2)
  {
    flowcont = cont1;
    flowptcont = cont2;
  }
  TH1D* calculate(ObsFunc func, int i);

  // Flow functions
  TH1D* getVNN(int nhar, int npar, std::string id, int i);
  TH1D* getV22(int i) { return getVNN(2, 2, "ChGap", i); }
  TH1D* getV32(int i) { return getVNN(3, 2, "ChGap", i); }
  TH1D* getV42(int i) { return getVNN(4, 2, "ChGap", i); }
  TH1D* getNN(int npar, int nhar, std::string id, int i);
  TH1D* get22(int i) { return getNN(2, 2, "ChGap", i); }
  TH1D* get24(int i) { return getNN(2, 4, "ChGap", i); }
  TH1D* get32(int i) { return getNN(3, 2, "ChGap", i); }
  TH1D* get42(int i) { return getNN(4, 2, "ChGap", i); }
  // Flow-pt cumulants
  TH1D* getCVN2pt(int nhar, int i);
  TH1D* getCV22pt(int i) { return getCVN2pt(2, i); }
  TH1D* getCV32pt(int i) { return getCVN2pt(3, i); }
  TH1D* getCV42pt(int i) { return getCVN2pt(4, i); }
  // Could in theory be generalised to CVN2ptX, but will we ever need those? In any case trivial to change
  TH1D* getCV22pt2(int i);
  TH1D* getCV22pt3(int i);
  TH1D* getCV22pt4(int i);
  TH1D* getCV24pt(int i);
  TH1D* getCV24pt2(int i);
  // Flow-pt normalised cumulants

  // Pt-functions
  TH1D* getMpt(int i) { return dynamic_cast<TH1D*>(flowptcont->getCorrHist(i, 1)); }

 private:
  // Data members
  FlowContainer* flowcont; //!
  FlowPtContainer* flowptcont; //!
};
#endif // PWGCF_GENERICFRAMEWORK_CORE_CALCULATOR_H_
