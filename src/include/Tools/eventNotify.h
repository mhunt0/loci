//#############################################################################
//#
//# Copyright 2008, 2015, Mississippi State University
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
#ifndef EVENTNOTIFY_H
#define EVENTNOTIFY_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>

#include <Tools/debug.h>
#include <Tools/identity.h>
#include <Tools/lmutex.h>

#include <list>

namespace Loci {
    
  class eventNotify : public Identity {
  public:
    eventNotify() { }
    virtual ~eventNotify() ;
    virtual void notification()  = 0 ;
  } ;


  class eventDispatcher {
    typedef std::list<eventNotify *> notify_list ;
    notify_list notify_group ;
    lmutex mutex ;
  public:
    eventDispatcher() ;
    ~eventDispatcher() ;
    void engage(eventNotify *) ;
    void disengage(eventNotify *) ;
    void dispatch_notify() ;
  } ;
}

#endif
