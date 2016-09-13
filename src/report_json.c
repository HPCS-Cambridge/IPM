/* Modified version of report_xml.c to output in JSON */
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <errno.h>

#include "ipm.h"
#include "ipm_types.h"
#include "hashtable.h"
#include "calltable.h"
#include "ipm_time.h"
#include "report.h"
#include "regstack.h"

#ifdef HAVE_MPI
#include <mpi.h>
#include "mod_mpi.h"
#include "GEN.calltable_mpi.h"
#else
#define IS_P2P_CALL(call_) 0
#endif

#ifdef HAVE_PAPI
#include "mod_papi.h"
#endif

#ifdef HAVE_CLUSTERING
#include "mod_clustering.h"
#endif

#ifdef HAVE_OMPTRACEPOINTS
#include "mod_omptracepoints.h"
#endif

#ifdef HAVE_POSIXIO
#include "mod_posixio.h"
#include "GEN.calltable_posixio.h"
#endif

#ifdef HAVE_MPIIO
#include "mod_mpiio.h"
#include "GEN.calltable_mpiio.h"
#endif

#ifdef HAVE_CUDA
#include "mod_cuda.h"
#endif

#define htable_scan(a_, b_, c_) htable_scan_named(a_, b_, c_, __FUNCTION__)

// AT - necessary prototypes
void json_regions(void *ptr, taskdata_t *t, struct region *reg, ipm_hent_t *htab);
void json_noregion(void *ptr, taskdata_t *t, region_t *reg, ipm_hent_t *htab);

char json_logfname[MAXSIZE_FILENAME] = {'\0'}; // AT

/* AT - manage JSON indent depth*/
int json_depth = 0;

/* map internal region IDs to IDs used in the xml file
   this is neccessary to handle the case where we only
   print the first level of regions */
int internal2xml[MAXNUM_REGIONS];

#define PRINT_TO_FILE      0
#define PRINT_TO_BUFFER    1

static int print_selector=PRINT_TO_FILE;
static int print_offset=0;
static unsigned long print_flags=0;

void json_cmdline(void *ptr, taskdata_t *t)
{
  ipm_printf(ptr, "%*c\"cmdline\": {\n", json_depth++, ' ');
  ipm_printf(ptr, "%*c\"md5sum\": \"%s\",\n", json_depth, ' ', t->exec_md5sum);
  ipm_printf(ptr, "%*c\"realpath\": \"%s\",\n", json_depth, ' ', t->exec_realpath);
  ipm_printf(ptr, "%*c\"cmd\": \"%s\"\n", json_depth, ' ', t->exec_cmdline);
  ipm_printf(ptr, "%*c},\n", --json_depth, ' ');
}

void json_job(void *ptr, taskdata_t *t)
{
  ipm_printf(ptr, "%*c\"job-env\": {\n", json_depth++, ' ');
  ipm_printf(ptr, "%*c\"jobid\": %s,\n", json_depth, ' ', t->jobid);
  ipm_printf(ptr, "%*c\"user\": \"%s\",\n", json_depth, ' ', t->user);
  ipm_printf(ptr, "%*c\"allocation\": \"%s\",\n", json_depth, ' ', t->allocation);
  ipm_printf(ptr, "%*c\"ipm\": \"%s\",\n", json_depth, ' ', IPM_VERSION);
  ipm_printf(ptr, "%*c\"t_start\": %d,\n", json_depth, ' ', t->t_start.tv_sec);
  ipm_printf(ptr, "%*c\"t_end\": %d,\n", json_depth, ' ', t->t_stop.tv_sec);
  ipm_printf(ptr, "%*c\"t_total\": %d,\n", json_depth, ' ', t->t_stop.tv_sec - t->t_start.tv_sec);
  ipm_printf(ptr, "%*c\"hosts\": %d,\n", json_depth, ' ', t->nhosts);
  ipm_printf(ptr, "%*c\"tasks\": %d\n", json_depth, ' ', t->ntasks);
  ipm_printf(ptr, "%*c},\n", --json_depth, ' ');
}

void json_profile_header(void *ptr, taskdata_t *t)
{
  ipm_printf(ptr, "{\n");
  json_depth++;
  json_cmdline(ptr, t);
  json_job(ptr, t);
}

void json_profile_footer(void *ptr)
{
  ipm_printf(ptr, "%*c}\n", --json_depth, ' ');
  ipm_printf(ptr, "}");
}

void json_task_header(void *ptr, taskdata_t *t)
{
  ipm_printf(ptr, "%*c\"%d\": {\n", json_depth++, ' ', t->taskid);
  ipm_printf(ptr, "%*c\"hostname\": \"%s\",\n", json_depth, ' ', t->hostname);
  ipm_printf(ptr, "%*c\"pid\": %d,\n", json_depth, ' ', t->pid);
  ipm_printf(ptr, "%*c\"t_init\": %.6f,\n", json_depth, ' ', IPM_TIMEVAL(t->t_start));
  ipm_printf(ptr, "%*c\"t_fini\": %.6f,\n", json_depth, ' ', IPM_TIMEVAL(t->t_stop));
}

void json_task_footer(void *ptr)
{
  ipm_printf(ptr, "%*c},\n", --json_depth, ' ');
}

void json_perf(void *ptr, taskdata_t *t)
{
  double procmem;
  double gflops;
  region_t *reg;
  
  reg = &ipm_app;

  procmem = task.procmem;
#ifdef HAVE_PAPI
  gflops = ipm_papi_gflops(reg->ctr, reg->wtime);
#else
  gflops = 0.0;
#endif

  ipm_printf(ptr, "%*c\"perf\": {\n", json_depth++, ' ');
  ipm_printf(ptr, "%*c\"wtime\": %f,\n", json_depth, ' ', IPM_TIMEVAL(t->t_stop)-IPM_TIMEVAL(t->t_start));
  ipm_printf(ptr, "%*c\"utime\": %f,\n", json_depth, ' ', t->utime);
  ipm_printf(ptr, "%*c\"stime\": %f,\n", json_depth, ' ', t->stime);
  ipm_printf(ptr, "%*c\"mtime\": %f,\n", json_depth, ' ', t->mtime);
  ipm_printf(ptr, "%*c\"gflop\": %f,\n", json_depth, ' ', gflops);
  ipm_printf(ptr, "%*c\"gbyte\": %f\n", json_depth, ' ', procmem);
  ipm_printf(ptr, "%*c},\n", --json_depth, ' ');
}

void json_print_instances(void *ptr, taskdata_t *t, ipm_hent_t *htab, int actv, int reg);

void json_func(void *ptr, taskdata_t *t, region_t *reg, ipm_hent_t *htab, int actv)
{
  int nkey;
  scanspec_t spec;
  scanstats_t stats;
  region_t *tmp;
  int id, res=0;

  /* if needed, walk up region tree and assign
     same id as closest parent has */
  if( internal2xml[reg->id]<0 ) {
    tmp=reg->parent;
    while(tmp) {
      id=internal2xml[tmp->id];
      if( id>=0 ) {
        internal2xml[reg->id]=id;
        break;
      }
      tmp=tmp->parent;
    }
  }

  scanspec_unrestrict_all(&spec);
  scanspec_restrict_activity(&spec, actv, actv);
  scanspec_restrict_region(&spec, reg->id, reg->id);
  //fprintf(stderr, "JSON_FUNC REGION: %d\n", reg->id);

  nkey = htable_scan( htab, &stats, spec );
  if( nkey>0 ) {
    //fprintf(stderr, "NAME: %s\n", ipm_calltable[actv].name);
    ipm_printf(ptr, "%*c\"%s\": {\n", json_depth++, ' ', ipm_calltable[actv].name);
    ipm_printf(ptr, "%*c\"bytes\": %f,\n", json_depth, ' ', stats.bytesum);
    ipm_printf(ptr, "%*c\"count\": %d,\n", json_depth, ' ', stats.hent.count);
    ipm_printf(ptr, "%*c\"instances\": {\n", json_depth++, ' ');
    json_print_instances(ptr, t, htab, actv, reg->id);
    ipm_printf(ptr, "%*c},\n", --json_depth, ' ');
    ipm_printf(ptr, "%*c\"time\": %f\n", json_depth, ' ', stats.hent.t_tot);
    ipm_printf(ptr, "%*c},\n", --json_depth, ' ');
  }

  //if( !(reg->flags)&FLAG_PRINT_EXCLUSIVE) {
  //  /* also print the func entries for all sub-regions, recursively.
  //     -> this makes the listed <func> entries inclusive */
  //  tmp=reg->child;
  //  while(tmp) {
  //    json_func(ptr, t, tmp, htab, actv);
  //    tmp=tmp->next;
  //  }
  //}

  if( stats.hent.count && stats.hent.timestamps ) {
    free(stats.hent.timestamps);
  }
}

void json_region(void *ptr, taskdata_t *t, region_t *reg, ipm_hent_t *htab)
{
  int i, j;
  int offs, range;

  /* print the <func> entries for all modules */

  for( i=0; i<MAXNUM_MODULES; i++ ) {
    offs  = modules[i].ct_offs;
    range = modules[i].ct_range;

    if( !(modules[i].name) || range==0 )
      continue;

    for( j=0; j<range; j++ ) {
      if( !(ipm_calltable[offs+j].name) )
        continue;

      json_func(ptr, t, reg, htab, offs+j);
    }
  }

  // TODO bad & needs fixing & checking
  fseek(ptr, -2, SEEK_CUR);
  ipm_printf(ptr, "\n");

  /* print child regions, only if we are not
     restricting to 1st level regions only */
  if( ((t->flags)&FLAG_NESTED_REGIONS) &&
      !(reg->flags&FLAG_PRINT_EXCLUSIVE) && reg->child ) {
    json_regions(ptr, t, reg, htab);
  }


}

void json_regions(void *ptr, taskdata_t *t, struct region *reg, ipm_hent_t *htab)
{
  struct region *tmp;
  int nreg, res=0;
  int nextid = 0;

  nreg = 0;
  tmp = reg->child;

  while(tmp) {
    nreg++;
    tmp = tmp->next;
  }

  if (reg == t->rstack->child) {
    nreg++; //noregion
  }

  ipm_printf(ptr, "%*c\"Region: %d\": {\n", json_depth++, ' ', reg->id);
  tmp = reg->child;
  while(tmp) {

    if((t->flags)&FLAG_NESTED_REGIONS) {
      internal2xml[tmp->id] = (tmp->id)-1;
    }
    else {
      internal2xml[tmp->id] = ++nextid;
    }

    json_region(ptr, t, tmp, htab);
    tmp = tmp->next;
  }


  if(reg == t->rstack->child) {
    json_noregion(ptr, t, reg, htab);
  }
  ipm_printf(ptr, "%*c}\n", --json_depth, ' ');

}

void json_noregion(void *ptr, taskdata_t *t, region_t *reg, ipm_hent_t *htab)
{
  double wtime, utime, stime, mtime;
  region_t noregion, *tmp;
  int i, res=0;

  rstack_clear_region(&noregion);
  noregion.id=1;
  noregion.nexecs=reg->nexecs;
  sprintf(noregion.name, "ipm_noregion");
  noregion.flags |= FLAG_PRINT_EXCLUSIVE;
  noregion.child = reg->child;

  wtime = reg->wtime;
  utime = reg->utime;
  stime = reg->stime;
  mtime = reg->mtime;

#ifdef HAVE_PAPI
  for( i=0; i<MAXNUM_PAPI_EVENTS; i++ ) {
    noregion.ctr[i]=reg->ctr[i];
  }
#endif

  tmp = reg->child;
  while(tmp) {
    wtime -= tmp->wtime;
    utime -= tmp->utime;
    stime -= tmp->stime;
    mtime -= tmp->mtime;

#ifdef HAVE_PAPI
    for( i=0; i<MAXNUM_PAPI_EVENTS; i++ ) {
      noregion.ctr[i] -= tmp->ctr[i];
    }
#endif

    tmp = tmp->next;
  }

  noregion.wtime=wtime;
  noregion.utime=utime;
  noregion.stime=stime;
  noregion.mtime=mtime;

  json_region(ptr, t, &noregion, htab);
}

void json_print_instances(void *ptr, taskdata_t *t, ipm_hent_t *htab, int actv, int reg)
{
  int i, j;
  int slct, call, bytes,/* reg,*/ csite;
  int op, dtype;
  int inst_count;
  IPM_COUNT_TYPE count;
  IPM_RANK_TYPE orank;
  double tmi, tma, tto;
  int nkey, tid;
  scanstats_t stats;
  char buf[80];


  inst_count = 0;
  nkey = 0;

//#ifdef HAVE_MPI
//  nkey += htable_scan_activity( htab, &stats,
//      MPI_MINID_GLOBAL, MPI_MAXID_GLOBAL);
//  count = stats.hent.count;
//
//  if( stats.hent.count && stats.hent.timestamps ) {
//    free(stats.hent.timestamps);
//  }
//#endif
//
//#ifdef HAVE_OMPTRACEPOINTS
//  nkey += htable_scan_activity( htab, &stats,
//      OMP_MINID_GLOBAL, OMP_MAXID_GLOBAL);
//  count += stats.hent.count;
//
//  if( stats.hent.count && stats.hent.timestamps ) {
//    free(stats.hent.timestamps);
//  }
//#endif
//#ifdef HAVE_POSIXIO
//  nkey += htable_scan_activity( htab, &stats,
//      POSIXIO_MINID_GLOBAL, POSIXIO_MAXID_GLOBAL);
//  count += stats.hent.count;
//
//  if( stats.hent.count && stats.hent.timestamps ) {
//    free(stats.hent.timestamps);
//  }
//#endif
//
//#ifdef HAVE_MPIIO
//  nkey += htable_scan_activity( htab, &stats,
//      MPIIO_MINID_GLOBAL, MPIIO_MAXID_GLOBAL);
//  count += stats.hent.count;
//
//  if( stats.hent.count && stats.hent.timestamps ) {
//    free(stats.hent.timestamps);
//  }
//#endif /* HAVE_MPIIO */
//
//
  for( i=0; i<MAXSIZE_HASH; i++ )
  {
    if( htab[i].count==0 || KEY_GET_ACTIVITY(htab[i].key) != actv || (t->flags & FLAG_NESTED_REGIONS && (KEY_GET_REGION(htab[i].key) != reg)))
      continue;

    slct  = KEY_GET_SELECT(htab[i].key);

    if( slct==IPM_RESOURCE_BYTES_AND_RANK ) {
      bytes = KEY_GET_BYTES(htab[i].key);
    } else {
      bytes = 0;
    }

    tmi = htab[i].t_min;
    tma = htab[i].t_max;
    tto = htab[i].t_tot;

    ipm_printf(ptr, "%*c\"%d\": {\n", json_depth++, ' ', inst_count++);

    ipm_printf(ptr, "%*c\"region\": %d,\n", json_depth, ' ', reg);
    ipm_printf(ptr, "%*c\"bytes\": %d,\n", json_depth, ' ', bytes);
    ipm_printf(ptr, "%*c\"count\": %d,\n", json_depth, ' ', htab[i].count);

    ipm_printf(ptr, "%*c\"t_min\": %f,\n", json_depth, ' ', tmi);
    ipm_printf(ptr, "%*c\"t_max\": %f,\n", json_depth, ' ', tma);


    if (task.flags&FLAG_REPORT_TIMESTAMPS) {
      ipm_printf(ptr, "%*c\"t_tot\": %f,\n", json_depth, ' ', tto);

      ipm_printf(ptr, "%*c\"timestamps\": [\n", json_depth++, ' ');

      for ( j = 0; j < htab[i].count -1; j++ ) {
        ipm_printf(ptr, "%*c%.6f,\n", json_depth, ' ', htab[i].timestamps[j] - IPM_TIMEVAL(t->t_start));
      }
      //Don't print a comma at the last instance
      for (; j < htab[i].count; j++ ) {
        ipm_printf(ptr, "%*c%.6f\n", json_depth, ' ', htab[i].timestamps[j] - IPM_TIMEVAL(t->t_start));
      }

      ipm_printf(ptr, "%*c]\n", --json_depth, ' ');
    }
    else {
      ipm_printf(ptr, "%*c\"t_tot\": %f\n", json_depth, ' ', tto);
    }

    ipm_printf(ptr, "%*c},\n", --json_depth, ' ');
  }

  ipm_printf(ptr, "%*c\"instance_count\": %d\n", json_depth, ' ', inst_count);

}

void json_task(void *ptr, taskdata_t *td, ipm_hent_t *htab)
{
  region_t *ipm_main;
  int i;

  for( i=0; i<MAXNUM_REGIONS; i++ ) {
    internal2xml[i]=-1;
  }

  /* td->rstack->child is ipm_main */
  ipm_main = td->rstack->child;
  internal2xml[ipm_main->id]=0;

  json_task_header(ptr, td);
  json_perf(ptr, td);
  json_regions(ptr, td, ipm_main, htab);
  json_task_footer(ptr);
}

void json_report_set_filename()
{
  struct stat fstat;

  /* see if logdir is available on task 0 */
  if( task.taskid==0 ) {
    if( task.flags&FLAG_OUTFILE ) {
      strncpy(json_logfname, task.fname, MAXSIZE_FILENAME);

    } else {
      if( stat(task.logdir, &fstat) ) {
        IPMERR("logdir %s unavailable, using '.'\n", task.logdir);
        sprintf(task.logdir, ".");
      }
      sprintf(json_logfname, "%s/%s.ipm.json",
          task.logdir, task.fname);
    }
  }
}


#ifdef HAVE_MPI
int report_json_atroot(unsigned long flags) // AT
{
  FILE *f;
  char buf[80];
  int i, j, nreg;
  MPI_Status stat;
  /* other task, hashtable and regions */
  taskdata_t otask;
  ipm_hent_t ohash[MAXSIZE_HASH];
  region_t oregions[MAXNUM_REGIONS];
  region_t *ostack;
  int err;

#ifdef HAVE_CLUSTERING
  procstats_t ostats;
#endif

  print_selector=PRINT_TO_FILE;
  print_flags=flags;
  err=0;

  if( task.taskid==0 ) {

    f=fopen(json_logfname, "w");
    if(!f) {
      IPMERR("Could not open IPM log file: '%s'\n", json_logfname);
      return IPM_EOTHER;
    }

    json_profile_header(f, &task);
    fflush(f);

    ipm_printf(f, "%*c\"tasks\": {\n", json_depth++, ' ');
    json_task(f, &task, ipm_htable);

    fflush(f);

    for( i=1; i<task.ntasks; i++ )
    {
#ifdef HAVE_CLUSTERING
      if( flags&XML_CLUSTERED ) {

        IPM_RECV( &ostats, sizeof(procstats_t),
            MPI_BYTE, i, 36, MPI_COMM_WORLD, &stat );

        /* if rank i is not a new cluster center, we can just
           continue with the next rank. if it is, we receive all the
           information we need to output the task profile */
        if( ostats.clrank != ostats.myrank ) {
          xml_taskcopy(f, &ostats);
          continue;
        }
      }
#endif

      IPM_RECV( ohash, sizeof(ipm_hent_t)*MAXSIZE_HASH,
          MPI_BYTE, i, 33, MPI_COMM_WORLD, &stat);
      // AT - Get timestamp memory from sister processes. TODO - fix tag.
      for( j = 0; j < MAXSIZE_HASH; j++ ) {
        if( ohash[j].count ) {
          ohash[j].timestamps = malloc(sizeof(double)*ohash[j].count);
          if( !ohash[j].timestamps ) {
            MPI_Abort(MPI_COMM_WORLD, ENOMEM);
          }
          IPM_RECV(ohash[j].timestamps, ohash[j].count,
              MPI_DOUBLE, i, 1337, MPI_COMM_WORLD, &stat);
        }
      }


#ifdef HAVE_CUDA
      IPM_RECV( &cudaptr, sizeof(cudaptr), MPI_BYTE, i, 39,
          MPI_COMM_WORLD, &stat );
#endif

      IPM_RECV( &otask, sizeof(taskdata_t), MPI_BYTE, i, 34,
          MPI_COMM_WORLD, &stat );

      IPM_RECV( oregions, sizeof(region_t)*MAXNUM_REGIONS, MPI_BYTE,
          i, 35, MPI_COMM_WORLD, &stat );

      ostack = rstack_unpack(MAXNUM_REGIONS, oregions);

      assert(ostack);
      assert(ostack->child);

      otask.rstack = ostack;

      json_task(f, &otask, ohash);

      // AT - free timestamp memory
      for( j = 0; j < MAXSIZE_HASH; j++ ) {
        if( ohash[j].count ) {
          free(ohash[j].timestamps);
        }
      }

      rstack_cleanup(ostack);
      if( ostack ) IPM_FREE(ostack);

      fflush(f);
    }

    // TODO bad & needs fixing & checking
    fseek(f, -2, SEEK_CUR);
    ipm_printf(f, "\n");

  // Here?? TODO
  //ipm_printf(f, "%*c},\n", --json_depth, ' ');

    json_profile_footer(f);
    chmod(json_logfname, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
    fclose(f);

  } else {

    /* all other ranks send their information to rank 0 */

#ifdef HAVE_CLUSTERING
    if( flags&XML_CLUSTERED ) {

      IPM_SEND( &mystats, sizeof(procstats_t),
          MPI_BYTE, 0, 36, MPI_COMM_WORLD );

      /* check if this rank is a cluster center, if not
         we're done and don't need to send anything more.
         If it is a cluster center, we need to send all the other
         data structures to rank 0 so that it can output the task
         profile */
      if( mystats.clrank!=mystats.myrank )
        return IPM_OK;
    }
#endif

    IPM_SEND( ipm_htable, sizeof(ipm_hent_t)*MAXSIZE_HASH,
        MPI_BYTE, 0, 33, MPI_COMM_WORLD );
    // AT - Send timestamp data to rank 0. Tag 1337 because why not (TODO -
    // normal tag).
    for( j = 0; j < MAXSIZE_HASH; j++ ) {
      if( ipm_htable[j].count ) {
        IPM_SEND(ipm_htable[j].timestamps, ipm_htable[j].count,
            MPI_DOUBLE, 0, 1337, MPI_COMM_WORLD);
      }
    }

#ifdef HAVE_CUDA
    IPM_SEND( &cudaptr, sizeof(cudaptr), MPI_BYTE, 0, 39,
        MPI_COMM_WORLD );
#endif

    IPM_SEND( &task, sizeof(taskdata_t), MPI_BYTE, 0, 34,
        MPI_COMM_WORLD );

    memset( oregions, 0, sizeof(region_t)*MAXNUM_REGIONS );
    rstack_pack( ipm_rstack, MAXNUM_REGIONS, oregions);

    IPM_SEND( oregions, sizeof(region_t)*MAXNUM_REGIONS, MPI_BYTE,
        0, 35, MPI_COMM_WORLD );
  }


  return IPM_OK;
}
#endif  /* HAVE_MPI */
