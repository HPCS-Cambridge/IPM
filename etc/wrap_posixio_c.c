
/** HEADER_BEGIN **/

#define _GNU_SOURCE
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "ipm_time.h"
#include "hashtable.h"
#include "perfdata.h"
#include "mod_posixio.h"
#include "ipm_core.h"
//AT
#include "report.h"

#ifdef HAVE_DYNLOAD
#include <execinfo.h>
#include <dlfcn.h>
#endif 

#ifdef HAVE_CALLPATH
#include "mod_callpath.h"
#endif

#ifdef HAVE_KEYHIST
#include "mod_keyhist.h"
#endif

#include <GEN.calltable_posixio.h>

#define MPI3CONST const

/** HEADER_END **/


__CRET__ __real___CFNAME__(__CPARAMS__);


/* ---- wrapping __CFNAME__ ---- */
/*
 * strings in the form __IDENT__ are replaced by the wrapper script
 *
 * CRET         __CRET__
 * CFNAME       __CFNAME__ 
 * CPARAMS      __CPARAMS__
 * CARGS        __CARGS__
 * CARGFMT      __CARGFMT__
 * CRETFMT      __CRETFMT__
 * GET_SSIZE    __GET_SSIZE__
 * GET_RSIZE    __GET_RSIZE__
 * GET_RANK     __GET_RANK__
 * GET_BYTES    __GET_BYTES__
 * RETURN_VALUE __RETURN_VALUE__
 */

#ifdef HAVE_DYNLOAD
__CRET__ __CFNAME__(__CPARAMS__)
#else
__CRET__ __wrap___CFNAME__(__CPARAMS__)
#endif
{
  static int loaded=0;
  static __CRET__ (*__CFNAME___real)(__CPARAMS__);
  int bytes, irank;
  int csite, idx, idx2, regid; 
  double tstart, tstop, t;
  IPM_KEY_TYPE key;

#if __RETURN_VALUE__
  __CRET__ rv;
#endif
  
#ifdef HAVE_DYNLOAD
  if(!loaded) {
    __CFNAME___real=0;
    __CFNAME___real=(__CRET__ (*)(__CPARAMS__)) dlsym(RTLD_NEXT, 
						      "__CFNAME__");
    
    if(!dlerror()) loaded=1;
    else {
      fprintf(stderr, "Error loading __CFNAME__ \n");
      /* handle error */
    }
  }
#endif /* HAVE_DYNLOAD */

  if( ipm_state==STATE_NOTINIT ) {
#ifndef HAVE_MPI
    ipm_init(0);
#endif
  }

  IPM_TIMESTAMP(tstart);

#if __RETURN_VALUE__
#ifdef HAVE_DYNLOAD
  rv=__CFNAME___real(__CARGS__);
#else
  rv=__real___CFNAME__(__CARGS__); 
#endif
#else
#ifdef HAVE_DYNLOAD
  __CFNAME___real(__CARGS__);
#else
  __real___CFNAME__(__CARGS__); 
#endif
#endif 
  IPM_TIMESTAMP(tstop); t=tstop-tstart;

  if( ipm_state!=STATE_ACTIVE || 
      modules[IPM_MODULE_POSIXIO].state!=STATE_ACTIVE ) { 
#if __RETURN_VALUE__
    return rv;
#else
    return;
#endif
  }
  
  __GET_BYTES__(bytes);

  /* Signal on bandwidth > threshold
	 * TODO - Configurable threshold/timestep
	 * TODO - Signal in a different way than printf to tty
	 */
#define BW_THRESH 130
#define BW_TIMESTEP 5

	static double t_bw = -1.0;
	static double mbytes = 0.0;
	double t_diff, bw;

	/* init static vars */
	if (t_bw < 0) {
		t_bw = tstop;
	}

	mbytes += bytes/1048574.0;

	t_diff = tstop - t_bw;

	if (t_diff > BW_TIMESTEP) {
		bw = mbytes / t_diff;
		if (bw > BW_THRESH) {
			printf("__CFNAME__ BW over %d MB/s (rolling %f)\n", BW_THRESH, bw);
		}
		t_bw = tstop;
		mbytes = 0.0;
	}

  csite=0;

#ifdef HAVE_CALLPATH
  /* csite=get_callsite_id(); */
#endif 

  regid=ipm_rstackptr->id;

#ifdef HAVE_POSIXIO_TRACE
  if( task.tracestate && task.tracefile ) {
    
    ipm_state=STATE_NOTACTIVE;
        
#if __RETURN_VALUE__
    fprintf(task.tracefile, "%10.6f __CRETFMT__ %s(__CARGFMT__) %u bytes %g sec\n", tstart, rv, "__CFNAME__", __CARGS__, bytes, t );
#else
    fprintf(task.tracefile, "%10.6f %s(__CARGFMT__) %u bytes %g sec\n", 
	    tstart, "__CFNAME__", __CARGS__, bytes, t);
#endif
    fflush(task.tracefile);    

    ipm_state=STATE_ACTIVE;

  }

#endif
    
  IPM_POSIXIO_KEY(key, __CFID___GLOBAL, 0, bytes, regid, csite);

  IPM_HASH_HKEY(ipm_htable, key, idx);

#ifdef HAVE_KEYHIST
  IPM_XHASH_HKEY(ipm_xhtable,last_hkey,key,idx2);
  ipm_xhtable[idx2].t_tot+=(tstart-last_tstamp);
  ipm_xhtable[idx2].count++;
  KEY_ASSIGN(last_hkey,key);
  last_tstamp=tstop;
#endif
  
  IPM_HASHTABLE_ADD(idx,t,tstart);

  // AT - Interval report

  if( task.flags & FLAG_REPORT_INTERVAL ) {
    // TODO - faster alt for gettimeofday?
    struct timeval tv;
    gettimeofday(&tv, 0);
    ipm_call_count++;

    if( (task.flags & FLAG_INTERVAL_CALL &&
          ipm_call_count == IPM_CALL_INTERVAL) ||
        (task.flags & FLAG_INTERVAL_TIME &&
         IPM_TIMEVAL(tv) - t_interval >= IPM_TIME_INTERVAL) ) {

      /* Need to reset these values before we dump the output -
       * report_xml_atinterval calls POSIX IO, meaning we'll otherwise have an
       * infinite loop here if we track those calls. This in turn causes a stack
       * overflow / too many open files.
       *
       * TODO: find a way to avoid this issue?
       */
      t_interval = IPM_TIMEVAL(tv);
      ipm_call_count = 0;

      report_xml_atinterval();
      for(int i = 0; i < MAXSIZE_HASH; i++) {
        ipm_htable[i].it_count = 0;
      }

      // Compensate for time spent dumping output - TODO necessary?
      gettimeofday(&tv, 0);
      t_interval = IPM_TIMEVAL(tv);
    }
  }

  // END AT
  
#ifdef HAVE_SNAP
 IPM_SNAP;
#endif

#if __RETURN_VALUE__
    return rv;
#else
    return;
#endif
}

