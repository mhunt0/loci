//#############################################################################
//#
//# Copyright 2025, Loci Consulting Services, LLC
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef DOCREPORT_H
#define DOCREPORT_H

#include <rule.h>
#include <iostream>

namespace Loci {
  void describeInputs(rule_db &rdb, std::ostream &o) ;
  inline void describeInputs(rule_db &rdb) { describeInputs(rdb,std::cout) ; }
  void ruleDocumentation(rule_db &rdb, variableSet vars,
                         std::ostream &o) ;
  inline void ruleDocumentation(rule_db &rdb, variableSet vars) {
    ruleDocumentation(rdb,vars,std::cout) ;
  }
}

#endif
