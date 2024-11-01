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

#ifndef PWGCF_GENERICFRAMEWORK_CORE_OUTPUTHANDLER_H_
#define PWGCF_GENERICFRAMEWORK_CORE_OUTPUTHANDLER_H_

#include <functional>
#include <string>
#include <vector>
#include "TH1.h"
#include "FlowContainer.h"
#include "FlowPtContainer.h"
#include "Calculator.h"
#include "StatHandler.h"

class OutputHandler : public TNamed
{
 public:
  OutputHandler();
  OutputHandler(std::string filename);
  ~OutputHandler();
  /*   OutputHandler(const OutputHandler& other);
    OutputHandler& operator=(const OutputHandler& other)
    {
      if (this != &other) {
      }
      return *this;
    }
    OutputHandler(const OutputHandler&& other) noexcept
    {
    }
    OutputHandler& operator=(const OutputHandler&& other) noexcept
    {
      if (this != &other) {
      }
      return *this;
    } */

  void initialise();

  TH1D* getHist(Calculator::ObsFunc func, StatHandler::StatMethod method = StatHandler::kBootstrap);

 private:
  FlowContainer* flowcont; //!
  FlowPtContainer* flowptcont; //!
  StatHandler* stathandler; //!
  //  SysHandler* sysHandler;
  std::string fileName;
  Calculator calculator; //!

  ClassDef(OutputHandler, 1);
};

#endif // PWGCF_GENERICFRAMEWORK_CORE_OUTPUTHANDLER_H_
