// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "OutputHandler.h"

ClassImp(OutputHandler);

OutputHandler::OutputHandler() : TNamed("", ""),
                                 fileName("")
{
}
OutputHandler::OutputHandler(std::string filename) : TNamed("", ""),
                                                     fileName(filename)
{
}
OutputHandler::~OutputHandler() {
    delete stathandler;
    delete flowcont;
    delete flowptcont;
};

void OutputHandler::initialise()
{
  auto inputfile = TFile::Open(fileName.c_str());
  std::string dirname = "generic-framework";
  flowcont = dynamic_cast<FlowContainer*>(inputfile->Get(Form("%s/FlowContainer", dirname.c_str())));
  if (!flowcont) {
    printf("Flowcontainer not found\n");
    inputfile->Close();
    return;
  }
  flowptcont = dynamic_cast<FlowPtContainer*>(inputfile->Get(Form("%s/FlowPtContainer", dirname.c_str())));
  if (!flowptcont) {
    printf("FlowPtcontainer not found\n");
    inputfile->Close();
    return;
  }
  inputfile->Close();
  calculator.setDataContainers(flowcont, flowptcont);
  stathandler = new StatHandler();
  // if(!statHandler) { }
  // sysHandler = new SysHandler();
  printf("Outputhandler initialised succesfully\n");
  return;
}
TH1D* OutputHandler::getHist(Calculator::ObsFunc func, StatHandler::StatMethod method){ //, int counter
    // std::vector<TH1D*> samplehists;
    // for(int i = -1;i<10;++i) { //counter was originally 10
    //   LOGF(info,"calculator %i",i);
    //   samplehists.push_back(calculator.calculate(func,i));
    // }
    // auto reth = stathandler->getNominal(samplehists,method);
    //sysHandler->applySyst(reth);

    auto reth = calculator.calculate(func,-1);

    return reth;
}