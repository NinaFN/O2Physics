#include "Calculator.h"

Calculator::Calculator()
{
}
Calculator::Calculator(FlowContainer* cont1, FlowPtContainer* cont2) : flowcont(cont1), flowptcont(cont2)
{
}
Calculator::~Calculator() {};
TH1D* Calculator::calculate(ObsFunc func, int i)
{
  return (this->*func)(i);
}
TH1D* Calculator::getCVN2pt(int nhar, int i)
{
  if (flowptcont->usesCentralMoments()) {
    LOGF(info,"----------------------in getCVN2pt with CM");
    LOGF(info,"harmonic = %i",nhar);
    auto vn2_mpt_profile = std::make_unique<BootstrapProfile>(*dynamic_cast<BootstrapProfile*>(flowptcont->GetCovList()->FindObject(Form("ChGap%i2pt1_Mpt0", nhar)))); // ChGap%i2pt1 ChGap%i2_mpt1
    TH1D* reth = dynamic_cast<TH1D*>(vn2_mpt_profile->getHist(i));
    
    auto vn2_hist = getNN(nhar, 2, "ChGap", i);
    auto mpt_hist(getMpt(i));
    vn2_hist->Multiply(mpt_hist);
    reth->Add(vn2_hist, -1);
    //LOGF(info,"add is done"); 
    return reth;
  }

  else{
    LOGF(info,"----------------------in getCVN2pt");
    LOGF(info,"harmonic = %i",nhar);
    auto vn2_mpt_profile = std::make_unique<BootstrapProfile>(*dynamic_cast<BootstrapProfile*>(flowptcont->GetCovList()->FindObject(Form("ChGap%i2pt1", nhar)))); // ChGap%i2pt1 ChGap%i2_mpt1
    TH1D* reth = dynamic_cast<TH1D*>(vn2_mpt_profile->getHist(i));
    //LOGF(info,"all good");
    //TH1D* reth = dynamic_cast<TH1D*>(vn2_mpt_profile->getHist(i));

   //LOGF(info,"getNN:");
    auto vn2_hist = getNN(nhar, 2, "ChGap", i);
    //LOGF(info,"getMpt:");
    auto mpt_hist(getMpt(i));
    //LOGF(info,"getMpt finished");
    vn2_hist->Multiply(mpt_hist);
    //LOGF(info,"multiply is done"); 
    reth->Add(vn2_hist, -1);
    //LOGF(info,"add is done"); 
    //LOGF(info,"reth name: %s",reth->GetName());
    return reth;
  }
  return nullptr;
}

TH1D* Calculator::getCV22pt2(int i)
{
  if (flowptcont->usesCentralMoments()) {
    auto v22pt2_Mpt0_profile = dynamic_cast<BootstrapProfile*>(flowptcont->GetCovList()->FindObject("ChGap22pt2_Mpt0"));
    // First get the <v2^2*deltapt^2> term
    TH1D* reth = dynamic_cast<TH1D*>(v22pt2_Mpt0_profile->getHist(i));
    auto mpt_hist = getMpt(i);
    auto v22pt2_Mpt1_profile = dynamic_cast<BootstrapProfile*>(flowptcont->GetCovList()->FindObject("ChGap22pt2_Mpt1"));
    auto v22pt_pt = dynamic_cast<TH1D*>(v22pt2_Mpt1_profile->getHist(i));
    v22pt_pt->Multiply(mpt_hist);
    reth->Add(v22pt_pt);
    auto v22_hist = getNN(2, 2, "ChGap", i);
    auto mpt_sq_hist = dynamic_cast<TH1D*>(mpt_hist->Clone("mpt_sq_hist"));
    mpt_sq_hist->Multiply(mpt_hist);
    auto v22_mptsq_hist = dynamic_cast<TH1D*>(v22_hist->Clone("v22_mptsq_hist"));
    v22_mptsq_hist->Multiply(mpt_sq_hist);
    reth->Add(v22_mptsq_hist);
    // Then subtract <v2^2><deltapt^2>
    auto v22_cm2_hist = dynamic_cast<TH1D*>(v22_hist->Clone("v22_cm2_hist"));
    v22_cm2_hist->Multiply(flowptcont->getCentralMomentHist(i, 2));
    reth->Add(v22_cm2_hist, -1);
    return reth;
  }

  return nullptr;
}

TH1D* Calculator::getCV22pt3(int i)
{
  if (flowptcont->usesCentralMoments()) {
    // First get the <v2^2*deltapt^3> term
    //START
    auto v22pt3_Mpt0_profile = dynamic_cast<BootstrapProfile*>(flowptcont->GetCovList()->FindObject("ChGap22pt3_Mpt0"));
    auto v22pt3_Mpt1_profile = dynamic_cast<BootstrapProfile*>(flowptcont->GetCovList()->FindObject("ChGap22pt3_Mpt1"));
    auto v22pt3_Mpt2_profile = dynamic_cast<BootstrapProfile*>(flowptcont->GetCovList()->FindObject("ChGap22pt3_Mpt2"));
    //LOGF(info,"finish access to CovList");
    TH1D* reth = dynamic_cast<TH1D*>(v22pt3_Mpt0_profile->getHist(i));
    //LOGF(info,"getMpt");
    auto mpt_hist = getMpt(i);
    auto v22pt2_pt = dynamic_cast<TH1D*>(v22pt3_Mpt1_profile->getHist(i));
    v22pt2_pt->Multiply(mpt_hist);
    mpt_hist->Multiply(mpt_hist);
    auto v22pt_pt2 = dynamic_cast<TH1D*>(v22pt3_Mpt2_profile->getHist(i));
    v22pt_pt2->Multiply(mpt_hist);
    mpt_hist->Multiply(mpt_hist);
    //LOGF(info,"getNN");
    auto v22_pt3 = getNN(2, 2, "ChGap", i);
    v22_pt3->Multiply(mpt_hist);
    reth->Add(v22pt2_pt, -1);
    reth->Add(v22pt_pt2, -1);
    reth->Add(v22_pt3, -1);
    //END


    // Then subtract the 3 x <v2^2*deltapt>*<deltapt^2> term
    //START
    auto v22deltapt_deltapt2 = getCV22pt(i);
    LOGF(info,"Multiplying (i,2)");
    v22deltapt_deltapt2->Multiply(flowptcont->getCentralMomentHist(i, 2));
    //LOGF(info,"Multiplying done");
    reth->Add(v22deltapt_deltapt2, -3);
    //END


    // And then subtract the <v2^2>*<deltapt^3> term
    //START
    auto v22_deltapt3 = getNN(2, 2, "ChGap", i);
    LOGF(info,"Multiplying (i,3)");
    v22_deltapt3->Multiply(flowptcont->getCentralMomentHist(i, 3));
    reth->Add(v22_deltapt3, -1);
    //END

    return reth;
  }
  return nullptr;
}
TH1D* Calculator::getCV22pt4(int i)
{
  if (flowptcont->usesCentralMoments()) {
    // First get the <v2^2*deltapt^4> term
    auto v22pt4_Mpt0_profile = dynamic_cast<BootstrapProfile*>(flowptcont->GetCovList()->FindObject("ChGap22pt4_Mpt0"));
    auto v22pt4_Mpt1_profile = dynamic_cast<BootstrapProfile*>(flowptcont->GetCovList()->FindObject("ChGap22pt4_Mpt1"));
    auto v22pt4_Mpt2_profile = dynamic_cast<BootstrapProfile*>(flowptcont->GetCovList()->FindObject("ChGap22pt4_Mpt2"));
    auto v22pt4_Mpt3_profile = dynamic_cast<BootstrapProfile*>(flowptcont->GetCovList()->FindObject("ChGap22pt4_Mpt3"));
    TH1D* reth = dynamic_cast<TH1D*>(v22pt4_Mpt0_profile->getHist(i));
    auto mpt_hist = getMpt(i);
    auto v22pt3_pt = dynamic_cast<TH1D*>(v22pt4_Mpt1_profile->getHist(i));
    v22pt3_pt->Multiply(mpt_hist);
    mpt_hist->Multiply(mpt_hist);
    auto v22pt2_pt2 = dynamic_cast<TH1D*>(v22pt4_Mpt2_profile->getHist(i));
    v22pt2_pt2->Multiply(mpt_hist);
    auto v22pt_pt3 = dynamic_cast<TH1D*>(v22pt4_Mpt3_profile->getHist(i));
    v22pt_pt3->Multiply(mpt_hist);
    mpt_hist->Multiply(mpt_hist);
    auto v22_pt4 = getNN(2, 2, "ChGap", i);
    v22_pt4->Multiply(mpt_hist);
    reth->Add(v22pt3_pt, -1);
    reth->Add(v22pt2_pt2, -1);
    reth->Add(v22pt_pt3, -1);
    reth->Add(v22_pt4, -1);
    // Then subtract the 6 x <v2^2*deltapt2>*<deltapt^2> term
    auto v22deltapt2_deltapt2 = getCV22pt2(i);
    v22deltapt2_deltapt2->Multiply(flowptcont->getCentralMomentHist(i, 2));
    reth->Add(v22deltapt2_deltapt2, -6);
    // Then subtract the 4 x <v2^2*deltapt>*<deltapt^3> term
    auto v22deltapt_deltapt3 = getCV22pt(i);
    v22deltapt_deltapt3->Multiply(flowptcont->getCentralMomentHist(i, 3));
    reth->Add(v22deltapt_deltapt3, -4);
    // Then add the 6 x <v2^2>*<deltapt^2>*<deltapt^2> term
    auto v22_deltapt2_deltapt2 = getNN(2, 2, "ChGap", i);
    v22_deltapt2_deltapt2->Multiply(flowptcont->getCentralMomentHist(i, 2));
    v22_deltapt2_deltapt2->Multiply(flowptcont->getCentralMomentHist(i, 2));
    reth->Add(v22_deltapt2_deltapt2, -4);
    // Then subtract the <v2^2>*<deltapt^4> term
    auto v22_deltapt4 = getNN(2, 2, "ChGap", i);
    v22_deltapt4->Multiply(flowptcont->getCentralMomentHist(i, 4));
    reth->Add(v22_deltapt4, 1);
    return reth;
  }
  return nullptr;
}
TH1D* Calculator::getCV24pt(int i)
{
  LOGF(info,"----------------------in getCV24pt");
  if (flowptcont->usesCentralMoments()) {
    // First get the <v2^4*deltapt> term
    BootstrapProfile* v24pt_Mpt0_profile = dynamic_cast<BootstrapProfile*>(flowptcont->GetCovList()->FindObject("ChGap24pt1_Mpt0"));
    TH1D* reth = dynamic_cast<TH1D*>(v24pt_Mpt0_profile->getHist(i));
    auto v24_pt = getNN(2, 4, "ChGap", i);
    auto mpt_hist = getMpt(i);
    v24_pt->Multiply(mpt_hist);
    reth->Add(v24_pt, -1);
    // Then subtract the 4 x <v2^2*deltapt>*<v2^2> term
    auto v22pt_v22 = getCV22pt(i);
    LOGF(info,"getCV22pt is done"); 
    v22pt_v22->Multiply(getNN(2, 2, "ChGap", i));
    reth->Add(v22pt_v22, -4);
    return reth;
  }

  else {
    // First get the <v2^4*deltapt> term
    BootstrapProfile* v24pt_Mpt0_profile = dynamic_cast<BootstrapProfile*>(flowptcont->GetCovList()->FindObject("ChGap24pt1"));
    TH1D* reth = dynamic_cast<TH1D*>(v24pt_Mpt0_profile->getHist(i));
    auto v24_pt = getNN(2, 4, "ChGap", i);
    auto mpt_hist = getMpt(i);
    v24_pt->Multiply(mpt_hist);
    reth->Add(v24_pt, -1);
    // Then subtract the 4 x <v2^2*deltapt>*<v2^2> term
    auto v22pt_v22 = getCV22pt(i);
    LOGF(info,"getCV22pt is done"); 
    v22pt_v22->Multiply(getNN(2, 2, "ChGap", i));
    reth->Add(v22pt_v22, -4);
    return reth;
  }

  return nullptr;
}
TH1D* Calculator::getCV24pt2(int i)
{
  if (flowptcont->usesCentralMoments()) {
    // First get the <v2^4*deltapt^2> term
    BootstrapProfile* v24pt2_Mpt0_profile = dynamic_cast<BootstrapProfile*>(flowptcont->GetCovList()->FindObject("ChGap24pt2_Mpt0"));
    BootstrapProfile* v24pt2_Mpt1_profile = dynamic_cast<BootstrapProfile*>(flowptcont->GetCovList()->FindObject("ChGap24pt2_Mpt1"));
    TH1D* reth = dynamic_cast<TH1D*>(v24pt2_Mpt0_profile->getHist(i));
    auto v24pt_pt = dynamic_cast<TH1D*>(v24pt2_Mpt1_profile->getHist(i));

    auto mpt_hist = getMpt(i);
    v24pt_pt->Multiply(mpt_hist);
    auto v24_pt2 = getNN(2, 4, "ChGap", i);
    mpt_hist->Multiply(mpt_hist);
    v24_pt2->Multiply(mpt_hist);
    reth->Add(v24pt_pt);
    reth->Add(v24_pt2);
    // Then subtract the 4 x <v2^2>*<vn^2*deltapt^2> term

    // Then subtract the 4 x <v2^2*deltapt>*<v2^2*deltapt> term

    // Then add the 4 x <v2^2>*<v2^2>*<deltapt^2> term

    // Then subtract the <v2^4>*<deltapt^2> term
    return reth;
  }
  return nullptr;
}
TH1D* Calculator::getVNN(int npar, int nhar, std::string id, int i)
{
  TH1D* reth = nullptr;
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  // Should be set via input config (json?)
  // if(fPropagateFCerrors) clone->SetPropagateErrors(kTRUE);
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
  };
  clone->SetIDName(id.c_str());
  switch (npar) {
    case 2:
      reth = dynamic_cast<TH1D*>(clone->GetVN2VsX(nhar, kFALSE, 0));
      break;
    case 4:
      reth = dynamic_cast<TH1D*>(clone->GetVN4VsX(nhar, kFALSE, 0));
      break;
    case 6:
      reth = dynamic_cast<TH1D*>(clone->GetVN6VsX(nhar, kFALSE, 0));
      break;
    case 8:
      reth = dynamic_cast<TH1D*>(clone->GetVN8VsX(nhar, kFALSE, 0));
      break;
    default:
      printf("Number of particles in correlation should be multiple of 2 and only up to 8 (for now)\n");
      reth = nullptr;
      break;
  }
  return reth;
}
TH1D* Calculator::getNN(int npar, int nhar, std::string id, int i)
{
  auto clone = (i < 0) ? flowcont : dynamic_cast<FlowContainer*>(flowcont->Clone("TempClone"));
  // Should be set via input config (json?)
  // if(fPropagateFCerrors) clone->SetPropagateErrors(kTRUE);
  if (i >= 0) {
    clone->OverrideMainWithSub(i, kFALSE);
    Int_t nxb;
    Double_t* xbs = flowcont->GetMultiRebin(nxb);
    clone->SetMultiRebin(nxb, xbs);
  };
  clone->SetIDName(id.c_str());
  //LOGF(info,"end of getNN pre-return");
  return dynamic_cast<TH1D*>(clone->GetHistCorrXXVsMulti(Form("%i%u", nhar, npar), 0));
}
