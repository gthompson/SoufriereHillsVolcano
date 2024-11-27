/*
 *   THIS FILE IS UNDER RCS - DO NOT MODIFY UNLESS YOU HAVE
 *   CHECKED IT OUT USING THE COMMAND CHECKOUT.
 *
 *    $Id$
 *
 *    Revision history:
 *     $Log$
 *     Revision 1.9  2002/01/08 17:02:29  cjbryan
 *     In order to accommodate 24-bit broadband data, traces of data of longer than long are now clipped to +/-2147483648
 *
 *     Revision 1.8  2001/04/12 03:47:36  lombard
 *     Major reorganization of SUDSPA_next to improve efficiency
 *     Reformatted code; cleaned up lots of comments and logit calls
 *     Removed mapping of NC component codes (only needed at UW)
 *     Now traces of data larger than short in are clipped to +/- 32767
 *     counts; before they were byte-truncated.
 *     Initializes SUDS structures to zero before filling them.
 *
 *     Revision 1.7  2001/03/27 22:22:51  cjbryan
 *     ensured that time portion of suds format is hhmmss and not hhmmss.ss as carried by EventTime variable
 *
 *     Revision 1.6  2001/03/22 20:55:12  cjbryan
 *     deleted now obsolete EventInst variable
 *
 *     Revision 1.5  2001/03/21 23:15:23  cjbryan
 *     integrated CVO subnet into output filename
 *
 *     Revision 1.4  2001/03/21 16:41:07  alex
 *     added EventSubnet variable to carry subnet name through to
 *     filename; deleted now extraneous EventMod variable
 *
 *     Revision 1.3  2001/03/21 02:23:03  alex
 *     *** empty log message ***
 *
 *     Revision 1.2  2000/03/10 23:23:43  davidk
 *     changed the sudsputaway routines to match the new PA_XXX
 *     routines in putaway.c.
 *     Moved the SUDSPA_XXX function prototypes to pa_other_head.h
 *
 *     Revision 1.1  2000/02/14 18:51:48  lucky
 *     Initial revision
 *
 *
 */

/* sudsputaway.c 
   
        Routines for writing trace data in PC-SUDS format 


  Thu Jul 22 14:56:37 MDT 1999 Lucky Vidmar
    Major structural changes to accomodate making all putaway routines
    into library calls so that they can be shared between trig2disk
    and wave2disk. Eliminated global variables set and allocated outside
    of this file. We should only have static globals from now on, with
    everything else being passed in via the putaway.c calls.
*/

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <earthworm.h>
#include <time.h>
#include <trace_buf.h>
#include <swap.h>
#include <ws_clientII.h>
#include <sudshead.h>
#include <pa_subs.h>

#define TAG_FILE        '.tag'        /* file containing the last known file tag */
#define MAXTXT           150


FILE    *SUDSfp;                      /* file pointer for the SUDS file */

static  long   *SudsBuffer;           /* write out SUDS data as long integers */
static  short  *SudsBufferShort;      /* write out SUDS data as short integers */
static  char SudsOutputFormat[MAXTXT];

/* Internal Function Prototypes */
static int StructMakeLocal (void *, int, char, int);
static int SwapDo (void *, int);

/************************************************************************
* Initialization function,                                              *
*       This is the Put Away startup intializer. This is called when    *
*       the system first comes up. Here is a chance to look around      *
*       and see if it's possible to do business, and to complain        *
*       if not ,BEFORE an event has to be processed.                    *
*                                                                       *
*       For PCSUDS, all we want to do is to make sure that the          *
*       directory where files should be written exists.                 *
*************************************************************************/
int SUDSPA_init (int OutBufferLen, char *OutDir, char *OutputFormat, 
                 int debug)
{
  /** Allocate SudsBuffer and SudsBufferShort 
       We waste RAM by allocating both the long and short buffers here
	 at the beginning of the code because some fluke (feature?) of NT 
	 which I don't understand becomes unhappy if we do the allocation 
	later. Win2000, of course, behaves differently, and is quite happy
	with buffer allocation after we have determined the format of the
	incoming data */
  if ((SudsBuffer = (long *) malloc (OutBufferLen * sizeof (char))) == NULL)
  {
    logit ("et", "SUDSPA_init: couldn't malloc SudsBuffer\n");
    return EW_FAILURE;
  }
  if ((SudsBufferShort = (short *) malloc (OutBufferLen * sizeof (char))) == NULL)
  {
    logit ("et", "SUDSPA_init: couldn't malloc SudsBufferShort\n");
    return EW_FAILURE;
  }

  /* Make sure that the top level output directory exists */
  if (CreateDir (OutDir) != EW_SUCCESS)
  {
    logit ("e", "SUDSPA_init: Call to CreateDir failed\n");
    return EW_FAILURE;
  }

  if(strlen(OutputFormat) >= sizeof(SudsOutputFormat))
  {
    logit("","SUDSPA_init: Error: OutputFormat(%s) is too long! Quitting!\n",
          OutputFormat);
    return(EW_FAILURE);
  }
  else
  {
    strcpy(SudsOutputFormat,OutputFormat);
  }
  return EW_SUCCESS; 
}

/************************************************************************
*   This is the Put Away event initializer. It's called when a snippet  *
*   has been received, and is about to be processed.                    *
*   It gets to see the pointer to the TraceRequest array,               *
*   and the number of loaded trace structures.                          *
*                                                                       *
*   For PCSUDS, we need to make sure that the target directory          *
*   exists, create it if it does not, then open the SUDS file           *
*   for writing.                                                        *
*************************************************************************/
int SUDSPA_next_ev (char *EventID, TRACE_REQ *ptrReq, int nReq, 
                    char *OutDir, char *EventDate, char *EventTime,
                    char *EventSubnet, int debug)

{
  char    SUDSFile[4*MAXTXT];
  char    hhmmss[7];

  /* Changed by Eugene Lublinsky, 3/31/Y2K */
  /* There are 2 modes of behavior now: the default one when insmod */
  /* is numeric or the second one when insmod is alphanumeric */
  char  tmpEventID[EVENTID_SIZE + 3];

  /* We treat suds file slightly differently - by request
   * from Tom Murray. The files are placed directly into the
   * top level suds directory, with names as follows:
   *
   *    yyyymmdd_hhmmss_insmod_evid.dmx
   */

  /* Changed by Eugene Lublinsky, 3/31/Y2K */
  /* choose which way to go */

  /* Build the file name */
  /* added by murray (and it shows) to set the eventid at 3 characters */
  if (strlen(EventID) < 3)
    sprintf (tmpEventID,"000%s",EventID) ;
  else
    sprintf (tmpEventID,"%s",EventID) ;

  /* ensure that time string is of the form hhmmss (revised EventTime 
     allows a string of the form hhmmss.ss, so truncate it if necessary) */
  if (strlen(EventTime) > 6)
  {
    strncpy(hhmmss, EventTime, 6);
    hhmmss[6] = '\0';
  }
  else
    strcpy(hhmmss, EventTime);

  /* changed by Carol 3/21/01: if no subnet, use network name
     in filename */
  if (EventSubnet[0] != '\0')
    sprintf (SUDSFile, "%s/%s_%s_%s_%s.dmx", OutDir,
             EventDate, hhmmss, EventSubnet, 
             &tmpEventID[strlen(tmpEventID)-3]);
  else
    sprintf (SUDSFile, "%s/%s_%s_%s_%s.dmx", OutDir,
             EventDate, hhmmss, ptrReq->net, 
             &tmpEventID[strlen(tmpEventID)-3]);
        
  /* end of changes */

  if (debug == 1)
    logit ("t", "Opening SUDS file %s\n", SUDSFile);

  /* open file */
  if ((SUDSfp = fopen (SUDSFile, "wb")) == NULL)
  {
    logit ("e", "SUDSPA_next_ev: unable to open file %s: %s\n", 
           SUDSFile, strerror(errno));
    return EW_FAILURE;
  }

  return (EW_SUCCESS);
}

/************************************************************************
* This is the working entry point into the disposal system. This        *
* routine gets called for each trace snippet which has been recovered.  *
* It gets to see the corresponding SNIPPET structure, and the event id  *
*                                                                       *
* For PCSUDS, this routine writes to the SUDS file, pointed to by the   *
* SUDSfp pointer, all of the received trace data in SUDS format:        *
*                                                                       *
*      1. SUDS tag - indicating what follows                            *
*      2. SUDS_STATIONCOMP struct - describe the station                *
*      3. SUDS tag - indicating what follows                            *
*      4. SUDS_DESCRIPTRACE struct - describe the trace data            *
*      5. trace data                                                    *
*                                                                       *
*  One bit of complexity is that we need to write the files in the      *
*  correct byte-order. Based on the OutputFormat parameter, determine   *
*  whether or not to swap bytes before writing the suds file.           *
*                                                                       *
* WARNING: we clip trace data to -2147483648 - +2147483647 so it will   *
*  fit in a long int. Any incoming data that is longer than 32 bits     *
*  will be CLIPPED. cjb 5/18/2001                                       *
*************************************************************************/
/* Process one channel of data */
int SUDSPA_next (TRACE_REQ *getThis, double GapThresh,
                 long OutBufferLen, int debug)
{
  TRACE_HEADER *wf;
  short  *s_data;
  long   *l_data;
  float  *f_data;
  char   *msg_p;        /* pointer into tracebuf data */
  char    datatype;
  int     data_size;
  int     j;
  int     gap_count = 0;
  long    nsamp, nfill;
  long    nfill_max = 0l;
  long    nsamp_this_scn = 0l;
  long    this_size;
  double  begintime, starttime, endtime;
  double  samprate;
  long    fill = 0l;
  long    min, max;
  int     total;
  SUDS_STRUCTTAG          tag;
  SUDS_DESCRIPTRACE       dt;
  SUDS_STATIONCOMP        sc;
  
  /* Check arguments */
  if (getThis == NULL)
  {
    logit ("e", "SUDSPA_next: invalid argument passed in.\n");
    return EW_FAILURE;
  }

  /* Start with a clean slate */
  memset(&tag, 0, sizeof(tag));
  memset(&dt, 0, sizeof(dt));
  memset(&sc, 0, sizeof(sc));
  
  /* Used for computing trace statistics */
  max = 4096;
  min = 0;
  total = 0;

  if ( (msg_p = getThis->pBuf) == NULL)   /* pointer to first message */
  {
    logit ("e", "SUDSPA_next: Message buffer is NULL.\n");
    return EW_FAILURE;
  }
  wf = (TRACE_HEADER *) msg_p;
  
  /* Look at the first TRACE_HEADER and get set up of action */
  if (WaveMsgMakeLocal(wf) < 0)
  {
    logit("e", "SUDSPA_next: unknown trace data type: %s\n",
          wf->datatype);
    return( EW_FAILURE );
  }

  nsamp = wf->nsamp;
  starttime = wf->starttime;
  endtime = wf->endtime;
  samprate = wf->samprate;
  if (samprate < 0.01)
  {
    logit("et", "unreasonable samplerate (%f) for <%s.%s.%s>\n",
          samprate, wf->sta, wf->chan, wf->net);
    return( EW_FAILURE );
  }
  begintime = starttime;
  datatype = 'n';
  if (wf->datatype[0] == 's' || wf->datatype[0] == 'i')
  {
    if (wf->datatype[1] == '2') datatype = 's';
    else if (wf->datatype[1] == '4') datatype = 'l';
  }
  else if (wf->datatype[0] == 't' || wf->datatype[0] == 'f')
  {
    if (wf->datatype[1] == '4') datatype = 'f';
  }
  if (datatype == 'n')
  {
    logit("et", "SUDSPA_next: unsupported datatype: %s\n", datatype);
    return( EW_FAILURE );
  }

  if (debug == 1)
    logit("et", "SUDSPA_next: working on <%s/%s/%s> datatype: %c \n",
			wf->sta, wf->chan, wf->net, datatype);

  /* loop through all the messages for this s-c-n */
  while (1) 
  {
    /* advance message pointer to the data */
    msg_p += sizeof(TRACE_HEADER);
        
    /* check for sufficient memory in output buffer */
    this_size = (nsamp_this_scn + nsamp ) * sizeof(long);
    if ( OutBufferLen < this_size )
    {
      logit( "e", "out of space for <%s.%s.%s>; saving long trace.\n",
             wf->sta, wf->chan, wf->net);
      break;
    }
  
	/* if data are floats, clip to longs cjb 5/18/2001 */
    switch( datatype )
    {
    case 's':
      s_data = (short *)msg_p;
      for ( j = 0; j < nsamp ; j++, nsamp_this_scn++ )
        SudsBuffer[nsamp_this_scn] = (long) s_data[j];
      msg_p += sizeof(short) * nsamp;
      break;
    case 'l':
      l_data = (long *)msg_p;
      for ( j = 0; j < nsamp; j++, nsamp_this_scn++ )
          SudsBuffer[nsamp_this_scn] = l_data[j];
      msg_p += sizeof(long) * nsamp;
      break;
    case 'f':
      f_data = (float *)msg_p;
      /* CLIP the data to long int */
      for ( j = 0; j < nsamp; j++, nsamp_this_scn++ )
      {
        if (l_data[j] < (float)LONG_MIN)
          SudsBuffer[nsamp_this_scn] = LONG_MIN;
        else if (l_data[j] > (float) LONG_MAX)
          SudsBuffer[nsamp_this_scn] = LONG_MAX;
        else
          SudsBuffer[nsamp_this_scn] = (long) l_data[j];
      }
      msg_p += sizeof(float) * nsamp;
      break;
    }
  
    /* End-check based on length of snippet buffer */
    if ((long) msg_p >= ((long) getThis->actLen + (long) getThis->pBuf)) 
    {
      if (debug == 1)
        logit ("", "Setting done for <%s.%s.%s>\n", wf->sta, wf->chan, 
               wf->net);
      break; /* Break out of the while(1) loop 'cuz we're done */
    }
  
    /* msg_p has been advanced to the next TRACE_BUF; localize bytes *
     * and check for gaps.                                           */
    wf = (TRACE_HEADER *) msg_p;
    if (WaveMsgMakeLocal(wf) < 0)
    {
      logit("e", "SUDSPA_next: unknown trace data type: %s\n",
            wf->datatype);
      return( EW_FAILURE );
    }
    nsamp = wf->nsamp;
    starttime = wf->starttime; 
    /* starttime is set for new packet; endtime is still set for old packet */
    if ( endtime + ( 1.0/samprate ) * GapThresh < starttime ) 
    {
      /* there's a gap, so fill it */
      logit("e", "gap in %s.%s.%s: %lf: %lf\n", wf->sta, wf->chan, wf->net,
            endtime, starttime - endtime);
      nfill = (long) (samprate * (starttime - endtime) - 1);
      if ( (nsamp_this_scn + nfill) * (long)sizeof(long) > OutBufferLen ) 
      {
        logit("e", 
              "bogus gap (%d); skipping\n", nfill);
        return(EW_FAILURE);
      }
      /* do the filling */
      for ( j = 0; j < nfill; j++, nsamp_this_scn ++ ) 
        SudsBuffer[nsamp_this_scn] = fill;
      /* keep track of how many gaps and the largest one */
      gap_count++;
      if (nfill_max < nfill) 
        nfill_max = nfill;
    }
    /* Advance endtime to the new packet;        *
     * process this packet in the next iteration */
    endtime = wf->endtime;
  } /* while(1) */
      
  /* figure out min, max, and "noise" */
  for (j = 0; j < 200 && j < nsamp_this_scn; j++)
  {
    if (SudsBuffer[j] > max)
      max = SudsBuffer[j];
    if (SudsBuffer[j] < min)
      min = SudsBuffer[j];
    total += SudsBuffer[j];
  }
  dt.avenoise = ((float)total)/ (float)j;  /* Mean of the first 200 samples */
  for (; j < nsamp_this_scn; j++)
  {
    if (SudsBuffer[j] > max)
      max = SudsBuffer[j];
    if (SudsBuffer[j] < min)
      min = SudsBuffer[j];
  }

  /* If the incoming data were originally short integers, copy the values
     back to an array of shorts; disk space saving feature requested
	 by Gabriel Reyes cjb 6/11/01 */
  if (datatype == 's') 
  {
    for (j = 0; j < nsamp_this_scn; j++)
		SudsBufferShort[j] = (short)SudsBuffer[j];
  }

  /* Convert to the appropriate output format */
#if defined (_INTEL)
  /* we are on intel, data will be read on sparc */
  if (strcmp (SudsOutputFormat, "sparc") == 0)
    for (j = 0; j < nsamp_this_scn; j++)
    {
		if (datatype == 's')
			SwapShort(&SudsBufferShort[j]);
		else
			SwapLong(&SudsBuffer[j]);
	}
#elif defined (_SPARC)
  /* we are on sparc, data will be read on intel */
  if (strcmp (SudsOutputFormat, "intel") == 0)
    for (j = 0; j < nsamp_this_scn; j++)
    {
		if (datatype == 's')
			SwapShort(&SudsBufferShort[j]);
		else
			SwapLong(&SudsBuffer[j]);
	}
#else
  logit ("e", "SUDSPA_next: Can't determine my platform - please compile with either _INTEL or _SPARC set\n");
      return EW_FAILURE;
#endif
               
  /* Write out to the SUDS file */
  /* Fill and write TAG for the STATIONCOMP struct */

  tag.id_struct = STATIONCOMP; /* what follows is STATIONCOMP */
  tag.len_struct = sizeof (SUDS_STATIONCOMP); 
  tag.len_data = 0;
  tag.sync = 'S';

  if (debug == 1)
    logit ("", "Writing tag for %d (%d)\n", tag.id_struct, tag.len_struct);

  if (strcmp (SudsOutputFormat, "sparc") == 0)
    tag.machine = '1';
  else if (strcmp (SudsOutputFormat, "intel") == 0)
    tag.machine = '6';

  if (StructMakeLocal ((void *) &tag, STRUCTTAG, tag.machine, debug) 
      != EW_SUCCESS)
  {
    logit ("et", "SUDSPA_next: Call to StructMakeLocal failed. \n");
    return EW_FAILURE;
  }

  if (fwrite ((void *) &tag, sizeof (SUDS_STRUCTTAG), 1, SUDSfp) != 1)
  {
    logit ("et", "SUDSPA_next: error writing SUDS tag. \n");
    return EW_FAILURE;
  }

  /* Fill and write STATIONCOMP struct */
  /* in SUDS_STATIDENT structure, char st_name[5],
	char network[4], char component where component = v,n,e
	cjb 5/18/2001 */
  if (strlen(wf->sta) > 5) {
	strncpy(sc.sc_name.st_name, wf->sta, 4);
	sc.sc_name.st_name[4] = '\0';
  }
  else
	strcpy (sc.sc_name.st_name, wf->sta); 

 if (strlen(wf->net) > 4) {
	strncpy(sc.sc_name.network, wf->net, 3);
	sc.sc_name.network[3] = '\0';
  }
  else
	strcpy (sc.sc_name.network, wf->net); 

  if (strlen(wf->chan) >= 3)
	sc.sc_name.component = wf->chan[2];
  else 
	sc.sc_name.component = wf->chan[0];

  if ((sc.sc_name.component == 'Z') || (sc.sc_name.component == 'z'))
	  sc.sc_name.component = 'V';

  if ((sc.sc_name.component != 'N') && (sc.sc_name.component != 'n')
	  && (sc.sc_name.component != 'E') && (sc.sc_name.component != 'e')
	  && (sc.sc_name.component != 'V') && (sc.sc_name.component != 'v')
	  && (sc.sc_name.component != 'T') && (sc.sc_name.component != 't'))
	  logit("et", "SUDSPA_next: unknown station component %c \n",
		sc.sc_name.component);

  if (datatype == 's')
	sc.data_type = 's';
  else
	  sc.data_type = 'l';
  sc.data_units = 'd';

  if (debug == 1)
    logit ("", "Writing STATIONCOMP struct for %s.%s.%c\n", 
           sc.sc_name.st_name,
           sc.sc_name.network,     
           sc.sc_name.component);

  if (StructMakeLocal ((void *) &sc, STATIONCOMP, tag.machine, debug) 
      != EW_SUCCESS)
  {
    logit ("et", "SUDSPA_next: Call to StructMakeLocal failed. \n");
    return EW_FAILURE;
  }

  if (fwrite ((void *) &sc, sizeof (SUDS_STATIONCOMP), 1, SUDSfp) != 1)
  {
    logit ("et", "SUDSPA_next: error writing SUDS_STATIONCOMP struct. \n");
    return EW_FAILURE;
  }

  /* Fill and write TAG for the DESCRIPTRACE struct */
  if (datatype == 's')
	data_size = nsamp_this_scn * sizeof (short);
  else
	data_size = nsamp_this_scn * sizeof (long);
  tag.id_struct = DESCRIPTRACE; /* what follows is DESCRIPTRACE */
  tag.len_struct = sizeof (SUDS_DESCRIPTRACE); 
  tag.len_data = (long) data_size;

  if (debug == 1)
    logit ("", "Writing tag for %d (%d)\n", tag.id_struct, tag.len_struct);

  if (StructMakeLocal ((void *) &tag, STRUCTTAG, tag.machine, debug) != EW_SUCCESS)
  {
    logit ("et", "SUDSPA_next: Call to StructMakeLocal failed. \n");
    return EW_FAILURE;
  }

  if (fwrite ((void *) &tag, sizeof (SUDS_STRUCTTAG), 1, SUDSfp) != 1)
  {
    logit ("et", "SUDSPA_next: error writing SUDS tag. \n");
    return EW_FAILURE;
  }

  /* Fill and write DESCRIPTRACE struct */
  /*strcpy (dt.dt_name.st_name, wf->sta); 
  strcpy (dt.dt_name.network, wf->net); 
  dt.dt_name.component = wf->chan[0];
  */
  /* in SUDS_STATIDENT structure, char st_name[5],
	char network[4], char component where component = v,n,e
	cjb 5/18/2001 */
  if (strlen(wf->sta) > 5) {
	strncpy(dt.dt_name.st_name, wf->sta, 4);
	dt.dt_name.st_name[4] = '\0';
  }
  else
	strcpy (dt.dt_name.st_name, wf->sta); 

 if (strlen(wf->net) > 4) {
	strncpy(dt.dt_name.network, wf->net, 3);
	dt.dt_name.network[3] = '\0';
  }
  else
	strcpy (dt.dt_name.network, wf->net); 

  if (strlen(wf->chan) >= 3)
    dt.dt_name.component = wf->chan[2];
  else 
    dt.dt_name.component = wf->chan[0];

  if ((dt.dt_name.component == 'Z') || (dt.dt_name.component == 'z'))
	  dt.dt_name.component = 'V';

  if ((dt.dt_name.component != 'N') && (dt.dt_name.component != 'n')
	  && (dt.dt_name.component != 'E') && (dt.dt_name.component != 'e')
	  && (dt.dt_name.component != 'V') && (dt.dt_name.component != 'v')
	  && (dt.dt_name.component != 'T') && (dt.dt_name.component != 't'))
	  logit("et", "SUDSPA_next: error with station component %c \n",
		dt.dt_name.component);


  if (datatype == 's')
	  dt.datatype = 's';
  else
	  dt.datatype = 'l';
  dt.begintime = begintime;
  dt.length = nsamp_this_scn;
  dt.rate = (float) wf->samprate;
  dt.mindata = (float) min;
  dt.maxdata = (float) max;

  /* Ignore the rest for now - see how it works */

  if (debug == 1)
    logit ("", "Writing DESCRIPTRACE - %d samples (%d,%d) \n", 
           nsamp_this_scn, min, max);

  if (StructMakeLocal ((void *) &dt, DESCRIPTRACE, tag.machine, debug) != EW_SUCCESS)
  {
    logit ("et", "SUDSPA_next: Call to StructMakeLocal failed. \n");
    return EW_FAILURE;
  }

  if (fwrite ((void *) &dt, sizeof (SUDS_DESCRIPTRACE), 1, SUDSfp) != 1)
  {
    logit ("et", "SUDSPA_next: error writing DESCRIPTRACE struct. \n");
    return EW_FAILURE;
  }

  /* write TRACE data - SudsBuffer holds long data;
	  SudsBufferShort holds short data */
  if (debug == 1)
    logit ("", "Writing %d bytes of DESCRIPTRACE data\n", data_size);
  if (datatype == 's')
  {
	if ((long)fwrite ((void *) SudsBufferShort, sizeof (char), data_size, SUDSfp)
		!= data_size)
	{
		logit ("et", "SUDSPA_next: error writing short TRACE data. \n");
		return EW_FAILURE;
	}
  }
  else
  {
	if ((long)fwrite ((void *) SudsBuffer, sizeof (char), data_size, SUDSfp)
		!= data_size)
	{
		logit ("et", "SUDSPA_next: error writing long TRACE data. \n");
		return EW_FAILURE;
	}
  }
  return EW_SUCCESS;
}


/************************************************************************
* This is the Put Away end event routine. It's called after we've       *
* finished processing one event                                         *
*                                                                       *
* For PC-SUDS - close the SUDS file, pointed to by SUDSfp               *
*               free SudsBuffer memory to be nice to everyone else      *
*************************************************************************/
int SUDSPA_end_ev(int debug)
{
  fclose (SUDSfp);
        
  if (debug == 1)
    logit("t", "Closing SUDS file \n");

  return( EW_SUCCESS );
}


/************************************************************************
*       This is the Put Away close routine. It's called after when      *
*       we're being shut down.                                          *
*************************************************************************/
int SUDSPA_close(int debug)
{

  free ((char *) SudsBufferShort);
  free ((char *) SudsBuffer);
  return( EW_SUCCESS );
}


/*
 *
 *  Byte swapping functions
 */

/* swap bytes if necessary, depending on what machine */
/* we have been compiled, and what machine the data   */
/* is being read from/being written to                */

static int StructMakeLocal(void *ptr, int struct_type, char data_type, 
                           int debug)
{
  if (ptr == NULL)
  {
    logit ("e", "StructMakeLocal: NULL pointer passed in\n");
    return EW_FAILURE;
  }

#if defined (_INTEL)
  if (data_type != '6')
  {
    if (debug == 1)
      logit ("", "Swapping from Intel to Sun \n");

    /* if we are on intel, and target data is not intel */
    if (SwapDo (ptr, struct_type) != EW_SUCCESS)
    {
      logit ("e", "StructMakeLocal: Call to SwapDo failed \n");
      return EW_FAILURE;
    }
  }

#elif defined (_SPARC)
  if (data_type == '6')
  {
    if (debug == 1)
      logit ("", "Swapping from Sun to Intel \n");

    /* if we are on solaris, and target data is not solaris */
    if (SwapDo (ptr, struct_type) != EW_SUCCESS)
    {
      logit ("e", "StructMakeLocal: Call to SwapDo failed \n");
      return EW_FAILURE;
    }
  }

#endif
  return EW_SUCCESS;
}


/* Do the dirty work by swapping the innards of SUDS */
/* structures that are of interest to us             */
static int SwapDo( void *ptr, int struct_type)
{
  SUDS_STRUCTTAG *tag;
  SUDS_STATIONCOMP *sc;
  SUDS_DESCRIPTRACE *dt;
  SUDS_TRIGGERS *tr;
  SUDS_DETECTOR *de;
  SUDS_TIMECORRECTION *tc;

  if (ptr == NULL)
  {
    logit ("e", "SwapDo: NULL pointer passed in.\n");
    return EW_FAILURE;
  }

  switch( struct_type ) 
  {
  case STRUCTTAG:
    tag = (SUDS_STRUCTTAG *) ptr;
    SwapShort (&tag->id_struct);
    SwapLong (&tag->len_struct);
    SwapLong (&tag->len_data);
    break;
  case STATIONCOMP:
    sc = (SUDS_STATIONCOMP *) ptr;
    SwapShort (&sc->sc_name.inst_type);
    SwapShort (&sc->azim);
    SwapShort (&sc->incid);
    SwapDouble (&sc->st_lat);
    SwapDouble (&sc->st_long);
    SwapFloat (&sc->elev);
    SwapShort (&sc->rocktype);
    SwapFloat (&sc->max_gain);
    SwapFloat (&sc->clip_value);
    SwapFloat (&sc->con_mvolts);
    SwapShort (&sc->channel);
    SwapShort (&sc->atod_gain);
    SwapLong (&sc->effective);
    SwapFloat (&sc->clock_correct);
    SwapFloat (&sc->station_delay);
    break;
  case DESCRIPTRACE:
    dt = (SUDS_DESCRIPTRACE *) ptr;
    SwapShort (&dt->dt_name.inst_type);
    SwapDouble (&dt->begintime);
    SwapShort (&dt->localtime);
    SwapShort (&dt->digi_by);
    SwapShort (&dt->processed);
    SwapLong (&dt->length);
    SwapFloat (&dt->rate);
    SwapFloat (&dt->mindata);
    SwapFloat (&dt->maxdata);
    SwapFloat (&dt->avenoise);
    SwapLong (&dt->numclip);
    SwapDouble (&dt->time_correct);
    SwapFloat (&dt->rate_correct);
    break;
  case TRIGGERS:
    tr = (SUDS_TRIGGERS *) ptr;
    SwapShort (&tr->tr_name.inst_type);
    SwapShort (&tr->sta);
    SwapShort (&tr->lta);
    SwapShort (&tr->abs_sta);
    SwapShort (&tr->abs_lta);
    SwapShort (&tr->trig_value);
    SwapShort (&tr->num_triggers);
    SwapDouble (&tr->trig_time);
    break;
  case DETECTOR:
    de = (SUDS_DETECTOR *) ptr;
    SwapFloat (&de->versionnum);
    SwapLong (&de->event_number);
    SwapLong (&de->spareL);
    break;
  case TIMECORRECTION:
    tc = (SUDS_TIMECORRECTION *) ptr;
    SwapShort (&tc->tm_name.inst_type);
    SwapDouble (&tc->time_correct);
    SwapFloat (&tc->rate_correct);
    SwapLong (&tc->effective_time);
    SwapShort (&tc->spareM);
    break;
  default:
    logit ("e", "SwapDo: Don't know about type %d\n", struct_type);
    return EW_FAILURE;
  }

  return EW_SUCCESS;
}

