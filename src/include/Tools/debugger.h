#ifndef DEBUGGER_H
#define DEBUGGER_H
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace Loci {
  void setup_debugger(const char *execname,const char *gdb, const char *hostname) ;
  void debugger_() ;
  void chopsigs_() ;
}



#endif
