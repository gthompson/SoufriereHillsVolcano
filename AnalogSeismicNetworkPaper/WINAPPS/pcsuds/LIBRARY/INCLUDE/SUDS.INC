!---  SUDS structure definitions and library constants
!---  09-Dec-1993, RB

!---  PC-SUDS version 1.44

!---  Include this file in any program unit that needs access to SUDS.

!---  Note that all strings (CHARACTER*(>1) variables) are NULL (char(0))
!---  terminated.

!----------------------------------------------------------------------

      ! SUDS_STATIDENT: Station identification
      STRUCTURE /suds_statident/
         CHARACTER*4 network     ! Network name
         CHARACTER*5 st_name     ! Station name
         CHARACTER*1 component   ! Component identifier
         INTEGER*2   inst_type   ! Instrument type identifier
      END STRUCTURE

      ! SUDS_ATODINFO: Information on the A to D converter
      STRUCTURE /suds_atodinfo/
         INTEGER*2   base_address     ! base I/O address of this device
         INTEGER*2   device_id        ! device identifier
         INTEGER*2   device_flags     ! device flags
         INTEGER*2   extended_bufs    ! number of extended buffers used
         INTEGER*2   external_mux     ! AtoD external mux control word
         CHARACTER*1 timing_source    ! AtoD timing source: i=internal, e=external
         CHARACTER*1 trigger_source   ! AtoD trigger source: i=internal, e=external
      END STRUCTURE

      ! SUDS_CALIBRATION: Calibration information for a station component
      STRUCTURE /suds_calibr/
         COMPLEX*8   pole             ! pole
         COMPLEX*8   zero             ! zero
      END STRUCTURE
      STRUCTURE /suds_calibration/
         RECORD /suds_statident/ ca_name ! station component identification
         REAL*4      maxgain          ! maximum gain of calibration curve
         REAL*4      normaliz         ! factor to multiply standard calib by to make
                                      ! peak at given frequency=1
         RECORD /suds_calibr/ cal(30) ! calibration info
         INTEGER*4   begint           ! time this calibration becomes effective
         INTEGER*4   endt             ! time this calibration is no longer effective
      END STRUCTURE

      ! SUDS_COMMENT:  Comment tag to be followed by the bytes of comment
      STRUCTURE /suds_comment/
         INTEGER*2   refer           ! structure identifier comment refers to
         INTEGER*2   item            ! item in structure comment refers to
         INTEGER*2   length          ! number of bytes in comment
         INTEGER*2   unused
      END STRUCTURE

      ! SUDS_CHANSET:  Associate station/components into sets.
      STRUCTURE /suds_chanset/
         INTEGER*2   type            ! Set type; 0=single channel(s), 1=orthogonal vector
         INTEGER*2   entries         ! Number of entries in set (these follow as data)
         CHARACTER*4 network         ! Network name
         CHARACTER*5 name            ! Set name
         INTEGER*4   active          ! Set is defined after this time
         INTEGER*4   inactive        ! Set is not defined after this time
      END STRUCTURE

      ! Entries of this form follow SUDS_CHANSET struct.
      STRUCTURE /chansetentry/
         INTEGER*4  inst_num         ! Instrument serial number
         INTEGER*2  stream_num       ! Stream of instrument
         INTEGER*2  chan_num         ! Channel of stream
         RECORD /suds_statident/ st  ! Station/component identification
      END STRUCTURE

      ! SUDS_DESCRIPTRACE:  Descriptive information about a seismic trace.
      !                     Normally followed by waveform.
      STRUCTURE /suds_descriptrace/
         RECORD /suds_statident/ dt_name  ! station component identification
         REAL*8      begintime        ! time of first data sample
         INTEGER*2   localtime        ! minutes to add to GMT to get local time
         CHARACTER*1 datatype         ! s = 12 bit unsigned stored as short int, 0 to 4096,
                                      ! q = 12 bit signed stored as short int, -2048 to 2048,
                                      ! u = 16 bit unsigned stored as short int, 0 to 65536
                                      ! i = 16 bit signed stored as short int, -32767 to 32767,
                                      ! r = 12 bit data, 4 lsb time stored as short int,
                                      ! l = long (32 bit signed integer),
                                      ! f = float (32 bit IEEE real),
                                      ! d = double (64 bit IEEE real),
                                      ! c = complex,
                                      ! v = vector,
                                      ! t = tensor
         CHARACTER*1 descriptor       ! g=good, t=telemetry noise, c=calibration, etc
         INTEGER*2   digi_by          ! agency code who digitized record 0=original
         INTEGER*2   processed        ! processing done on this waveform
         INTEGER*4   length           ! number of samples in trace
         REAL*4      rate             ! samples per second
         REAL*4      mindata          ! minimum value of data (type s,l,f only)
         REAL*4      maxdata          ! maximum value of data (type s,l,f only)
         REAL*4      avenoise         ! average value of first 200 samples (type s,l,f only)
         INTEGER*4   numclip          ! number of clipped datapoints
         REAL*8      time_correct     ! time correction to be added to begintime
         REAL*4      rate_correct     ! rate correction to be added to rate
      END STRUCTURE

      ! SUDS_DETECTOR:  Information on detector program being used
      STRUCTURE /suds_detector/
         CHARACTER*1 dalgorithm       ! triggering algorithm: x=xdetect, m=mdetect
                                      ! e=eqdetect
         CHARACTER*1 event_type       ! c=calibration, e=earthquake, E=explosion,
                                      ! f=free run, n=noise, etc.
         CHARACTER*10 net_node_id     ! network node identification
         REAL*4      versionnum       ! software version number
         INTEGER*4   event_number     ! unique event number assigned locally.
         INTEGER*4   spareL           ! spare
      END STRUCTURE

      ! SUDS_EQUIPMENT:  Equipment making up a station/component. Primarily used for
      !                  maintenance but may be referenced by researcher. One or more
      !                  structures exist for each piece of equipment making up a
      !                  station/component.
      STRUCTURE /suds_equipment/
         RECORD /suds_statident/ this      ! identifier of this piece of equipment
         RECORD /suds_statident/ previous  ! next piece of equipment toward sensor
         RECORD /suds_statident/ next      ! next piece of equipment toward recorder
         CHARACTER*8 serial           ! serial number
         INTEGER*2   model            ! model such as L4, HS10, etc.
         INTEGER*2   knob1            ! knob setting or series resistor value of Lpad
         INTEGER*2   knob2            ! knob setting or shunt  resistor value of Lpad
         INTEGER*2   reason           ! reason change was made
         REAL*4      frequency        ! sensor corner frequency, vco freq, transmitter
                                      ! frequency, etc.
         INTEGER*4   effective        ! date/time these values became effective
      END STRUCTURE

      ! SUDS_ERROR:  Error matrix
      STRUCTURE /suds_error/
         REAL*4      covarr(10)       ! covariance matrix
      END STRUCTURE

      ! SUDS_EVENT:  General information about an event.
      STRUCTURE /suds_event/
         INTEGER*2   authority        ! organization processing the data
         INTEGER*4   number           ! unique event number assigned by organization
         INTEGER*2   felt             ! number of felt reports
         CHARACTER*1 mintensity       ! maximum Modified Mercali Intensity
         CHARACTER*1 ev_type          ! e=earthquake, E=explosion, n=nuclear,
                                      ! i=icequake, b=b_type, n=net, r=regional,
                                      ! t=teleseism, c=calibration, n=noise
         CHARACTER*1 tectonism        ! observed u=uplift, s=subsidence, S=strikeslip
                                      ! faulting, N=normal faulting, T=thrust
         CHARACTER*1 waterwave        ! seiche, tsunami, etc.
         CHARACTER*1 mechanism        ! t=thrust, s=strike-slip, n=normal, e=explosive
         CHARACTER*1 medium           ! medium containing explosion or event
         REAL*4      size             ! magnitude or pounds TNT for explosions
      END STRUCTURE

      ! SUDS_EVENTSETTING:  Settings for earthquake trigger system
      STRUCTURE /suds_eventsetting/
         CHARACTER*4 netwname         ! network name
         REAL*8      beginttime       ! time these values in effect
         INTEGER*2   const1           ! trigger constant 1
         INTEGER*2   const2           ! trigger constant 2
         INTEGER*2   threshold        ! trigger threshold
         INTEGER*2   const3           ! trigger constant 3
         REAL*4      minduration      ! minimum duration for event
         REAL*4      maxduration      ! maximum duration for event
         CHARACTER*1 algorithm        ! triggering algorithm: x=xdetect, m=mdetect
                                      ! e=eqdetect
         CHARACTER*1 spareK           ! spare
         INTEGER*2   spareI           ! spare
      END STRUCTURE

      ! SUDS_EVDESCR:  Descriptive information about an event typically used for
      !                major, destructive earthquakes. This structure is typically
      !                associated with EVENT structure.
      STRUCTURE /suds_evdescr/
         CHARACTER*20 eqname     ! Popular name used to refer to this earthquake
         CHARACTER*16 country    ! country of earthquake
         CHARACTER*16 state      ! state, province or other political subdivision
         INTEGER*2   localtime   ! hours to add to GMT to get local time
         INTEGER*2   spareB
      END STRUCTURE

      ! SUDS_FEATURE:  Observed phase arrival time, amplitude, and period.
      STRUCTURE /suds_feature/
         RECORD /suds_statident/ fe_name ! station component identification
         INTEGER*2   obs_phase        ! observed phase code
         CHARACTER*1 onset            ! wave onset descriptor, i or e
         CHARACTER*1 direction        ! first motion: U,D,+,-
         INTEGER*2   sig_noise        ! ratio ampl. of first peak or trough to noise
         CHARACTER*1 data_source      ! i=interactive,a=automatic,r=rtp, or user code
         CHARACTER*1 tim_qual         ! timing quality given by analyst: 0-4, etc.
                                      !      n=ignore timing
         CHARACTER*1 amp_qual         ! amplitude quality given by analyst: 0-4, etc.
                                      !      n=ignor amplitude information
         CHARACTER*1 ampunits         ! units amplitude measured in: d=digital counts
                                      !      m=mm on develocorder, etc.
         INTEGER*2   gain_range       ! 1 or gain multiplier if gain range in effect
         REAL*8      time             ! phase time, x value where pick was made
         REAL*4      amplitude        ! peak-to-peak amplitude of phase
         REAL*4      period           ! period of waveform measured
         INTEGER*4   time_of_pick     ! time this pick was made
         INTEGER*2   pick_authority   ! organization processing the data
         INTEGER*2   pick_reader      ! person making this pick
      END STRUCTURE

      ! SUDS_FOCALMECH:  General information about a focal mechanism.
      STRUCTURE /suds_focalmech/
         REAL*4      astrike          ! strike of plane a
         REAL*4      adip             ! dip of plane a
         REAL*4      arake            ! rake of plane a
         REAL*4      bstrike          ! strike of plane b
         REAL*4      bdip             ! dip of plane b
         REAL*4      brake            ! rake of plane b
         CHARACTER*1 prefplane        ! preferred plane a or b or blank
         CHARACTER*3 spare
      END STRUCTURE

      ! SUDS_INSTRUMENT: Instrument hardware settings, mainly PADS related
      !                  added by R. Banfill, Jan 1991
      STRUCTURE /suds_instrument/
         RECORD /suds_statident/ in_name ! Station component identification
         INTEGER*2   in_serial        ! Instrument serial number
         INTEGER*2   comps            ! Number of components recorded by instrument
         INTEGER*2   channel          ! Actual channel number on recorder
         CHARACTER*1 sens_type        ! Sensor type a=accel, v=vel, d=disp...
         CHARACTER*1 datatype         ! see SUDS_DESCRIPTRACE.datatype
         INTEGER*4   void_samp        ! Invalid or void sample value
         REAL*4      dig_con          ! Digitizing constant (counts / volt)
         REAL*4      aa_corner        ! Anti-alias filter corner frequency (Hz)
         REAL*4      aa_poles         ! Anti-alias filter poles
         REAL*4      nat_freq         ! Transducer natural frequency (Hz)
         REAL*4      damping          ! Transducer damping coeff.
         REAL*4      mot_con          ! Transducer motion constant (volts / GMU)
         REAL*4      gain             ! Amplifier gain (dB)
         REAL*4      local_x          ! Local coordinate X (meters)
         REAL*4      local_y          ! Local coordinate Y (meters)
         REAL*4      local_z          ! Local coordinate Z (meters)
         INTEGER*4   effective        ! Time these setting took effect
         REAL*4      pre_event        ! Pre-event length (IST+pre_event=trigger time)
         INTEGER*2   trig_num         ! Trigger number on instrument
         CHARACTER*6 study            ! Study name, used to insure unique station names
         INTEGER*2   sn_serial        ! Sensor serial number
      END STRUCTURE

      ! SUDS_LAYERS:  Velocity layers.
      STRUCTURE /suds_layers/
         REAL*4      thickness        ! thickness in kilometers
         REAL*4      pveltop          ! p velocity at top of layer
         REAL*4      pvelbase         ! p velocity at base of layer
         REAL*4      sveltop          ! s velocity at top of layer
         REAL*4      svelbase         ! s velocity at base of layer
         INTEGER*2   function         ! velocity function in layer: 0=constant,
                                      !    1=linear, 2=exponential, etc.
         INTEGER*2   spare
      END STRUCTURE

      ! SUDS_LOCTRACE:  Location of trace.
      STRUCTURE /suds_loctrace/
         RECORD /suds_statident/ lt_name ! station component identification
         INTEGER*4   fileloc          ! pointer to pathname in file system
         INTEGER*4   tapeloc          ! pointer to name of tape or offline storage
         INTEGER*4   beginloc         ! bytes from begining of file to trace
      END STRUCTURE

      ! SUDS_MOMENT:  Moment tensor information.
      STRUCTURE /suds_moment/
         CHARACTER*1 datatypes        ! sum of: 1=polarities, 2=amplitudes,
                                      ! 4=waveforms, etc.
         CHARACTER*1 constraints      ! solution constrained: d=deviatoric,
                                      ! c=double couple
         CHARACTER*2 spareD
         REAL*4      sc_moment        ! scalar moment
         REAL*4      norm_ten(6)      ! normalized moment tensor
      END STRUCTURE

      ! SUDS_MUXDATA:  Header for multiplexed data
      STRUCTURE /suds_muxdata/
         CHARACTER*4 netname          ! network name
         REAL*8      begintime        ! time of first data sample
         INTEGER*2   loctime          ! minutes to add to GMT to get local time
         INTEGER*2   numchans         ! number of channels: if !=1 then multiplexed
         REAL*4      dig_rate         ! samples per second
         CHARACTER*1 typedata         ! see SUDS_DESCRIPTRACE.datatype
         CHARACTER*1 descript         ! g=good, t=telemetry noise, c=calibration, etc
         INTEGER*2   spareG           ! spare
         INTEGER*4   numsamps         ! number of sample sweeps. Typically not known
                                      ! when header is written, but can be added later
         INTEGER*4   blocksize        ! number of demultiplexed samples per channel if
                                      ! data is partially demultiplexed, otherwise=0
      END STRUCTURE

      ! SUDS_ORIGIN: Information about a specific solution for a given event
      STRUCTURE /suds_origin/
         INTEGER*4   number           ! unique event number assigned by organization
         INTEGER*2   authority        ! organization processing the data
         CHARACTER*1 version          ! version of solution within organization
         CHARACTER*1 or_status        ! processing status: f=final, a=automatic, etc
         CHARACTER*1 preferred        ! p=preferred location
         CHARACTER*1 program          ! name of processing program  h=hypo71,
                                      ! l=hypolayer, i=isc, c=centroid, etc.
         CHARACTER*1 depcontrl        ! depth control: f=fixed, etc.
         CHARACTER*1 convergence      ! hypocentral convergence character
         INTEGER*4   region           ! geographic region code assigned locally
         REAL*8      orgtime          ! origin time
         REAL*8      or_lat           ! latitude, north is plus
         REAL*8      or_long          ! longitude, east is plus
         REAL*4      depth            ! depth in kilometers, + down
         REAL*4      err_horiz        ! horizontal error in km
         REAL*4      err_depth        ! vertical error in km
         REAL*4      res_rms          ! rms of residuals
         CHARACTER*6 crustmodel       ! code for model used in this location
         INTEGER*2   gap              ! azimuthal gap in degrees
         REAL*4      nearstat         ! distance in km to nearest station
         INTEGER*2   num_stats        ! number of stations reporting phases
         INTEGER*2   rep_p            ! number of p phases reported
         INTEGER*2   used_p           ! number of p times used in the solution
         INTEGER*2   rep_s            ! number of s phases reported
         INTEGER*2   used_s           ! number of s times used in the solution
         INTEGER*2   mag_type         ! magnitude type: coda,tau,xmag ml,mb,ms,mw
         INTEGER*2   rep_m            ! number of magnitude readings reported
         INTEGER*2   used_m           ! number of magnitude readings used
         REAL*4      magnitude        ! magnitude value
         REAL*4      weight           ! average magnitude weight
         REAL*4      mag_rms          ! rms of magnitude
         INTEGER*4   effective        ! time this solution was calculated
      END STRUCTURE

      ! SUDS_PROFILE:  Grouping of shotgathers by profile.
      STRUCTURE /suds_profile/
         INTEGER*2   junk1            ! What is your suggestion?
      END STRUCTURE

      ! SUDS_RESIDUAL:  Calculated residuals for arrival times, magnitudes, etc.
      STRUCTURE /suds_residual/
         INTEGER*4   event_num        ! unique event number
         RECORD /suds_statident/ re_name ! station component identification
         INTEGER*2   set_phase        ! phase code set for this solution
         CHARACTER*1 set_tim_qual     ! timing quality assigned for this soln: 0-4
         CHARACTER*1 set_amp_qual     ! amplitude quality assigned for this soln: 0-4
         REAL*4      residual         ! traveltime residual or phase magnitude
         REAL*4      weight_used      ! weight used in this solution
         REAL*4      delay            ! delay time or station correction used
         REAL*4      azimuth          ! azimuth event to station, 0 north
         REAL*4      distance         ! distance in km event to station
         REAL*4      emergence        ! angle of emergence from source, 0=down,180=up
      END STRUCTURE

      ! SUDS_SHOTGATHER:  Grouping of waveforms by source event
      STRUCTURE /suds_shotgather/
         INTEGER*2   junk2            ! What is your suggestion?
      END STRUCTURE

      ! SUDS_STATIONCOMP:  Generic station component information
      STRUCTURE /suds_stationcomp/
         RECORD /suds_statident/ sc_name ! station component identification
         INTEGER*2   azim             ! component azimuth clockwise from north
         INTEGER*2   incid            ! component angle of incidence from vertical
                                      !      0 is vertical, 90 is horizontal
         REAL*8      st_lat           ! latitude, north is plus
         REAL*8      st_long          ! longitude, east is plus
         REAL*4      elev             ! elevation in meters
         CHARACTER*1 enclosure        ! d=dam, n=nuclear power plant, v=underground
                                      !      vault, b=buried, s=on surface, etc.
         CHARACTER*1 annotation       ! annotated comment code
         CHARACTER*1 recorder         ! type device data recorded on
         CHARACTER*1 rockclass        ! i=igneous, m=metamorphic, s=sedimentary
         INTEGER*2   rocktype         ! code for type of rock
         CHARACTER*1 sitecondition    ! p=permafrost, etc.
         CHARACTER*1 sensor_type      ! sensor type: d=displacement, v=velocity,
                                      ! a=acceleration, t=time code
         CHARACTER*1 data_type        ! see SUDS_DESCRIPTRACE.datatype
         CHARACTER*1 data_units       ! data units: d=digital counts, v=millivolts,
                                      ! n=nanometers (/sec or /sec/sec)
         CHARACTER*1 polarity         ! n=normal, r=reversed
         CHARACTER*1 st_status        ! d=dead, g=good
         REAL*4      max_gain         ! maximum gain of the amplifier
         REAL*4      clip_value       ! +-value of data where clipping begins
         REAL*4      con_mvolts       ! conversion factor to millivolts: mv per counts
                                      !      0 means not defined or not appropriate
                                      ! max_ground_motion=digital_sample*con_mvolts*
                                      ! max_gain
         INTEGER*2   channel          ! a2d channel number
         INTEGER*2   atod_gain        ! gain of analog to digital converter
         INTEGER*4   effective        ! date/time these values became effective
         REAL*4      clock_correct    ! clock correction in seconds.
         REAL*4      station_delay    ! seismological station delay.
      END STRUCTURE

      ! SUDS_TERMINATOR:  Structure to end a sequence of related structures when
      !                   loaded in a serial file or on a serial device.
      STRUCTURE /suds_terminator/
         INTEGER*2   structid         ! id for structure at beginning of this sequence
         INTEGER*2   spareA
      END STRUCTURE

      ! SUDS_TIMECORRECTION:  Time correction information.
      STRUCTURE /suds_timecorrection/
         RECORD /suds_statident/ tm_name ! time trace station id used to determine
                                         ! correction.
         REAL*8      time_correct     ! time correction to be added to begintime
         REAL*4      rate_correct     ! rate correction to be added to rate
         CHARACTER*1 sync_code        ! synchronization code as follows:
                                      !   0 = total failure,   1 = 1 second synch,
                                      !   2 = 10 second synch, 3 = minute synch,
                                      !   4, 5 = successful decode.
         CHARACTER*1 program          ! program used to decode time:
                                      !   e = irige, c = irigc
         INTEGER*4   effective_time   ! time this correction was calculated
         INTEGER*2   spareM
      END STRUCTURE

      ! SUDS_TRIGGERS:  Earthquake detector trigger statistics
      STRUCTURE /suds_triggers/
         RECORD /suds_statident/ tr_name ! station component identification
         INTEGER*2   sta              ! short term average
         INTEGER*2   lta              ! long term average pre_lta for xdetect
         INTEGER*2   abs_sta          ! short term absolute average
         INTEGER*2   abs_lta          ! long term absolute average
         INTEGER*2   trig_value       ! value of trigger level (eta)
         INTEGER*2   num_triggers     ! number of times triggered during this event
         REAL*8      trig_time        ! time of first trigger
      END STRUCTURE

      ! SUDS_TRIGSETTING:  Settings for earthquake trigger system
      STRUCTURE /suds_trigsetting/
         CHARACTER*4 netwname         ! network name
         REAL*8      beginttime       ! time these values in effect
         INTEGER*2   const1           ! trigger constant 1
         INTEGER*2   const2           ! trigger constant 2
         INTEGER*2   threshold        ! trigger threshold
         INTEGER*2   const3           ! trigger constant 3
         INTEGER*2   const4           ! trigger constant 4
         INTEGER*2   wav_inc          ! weighted average increment
         REAL*4      sweep            ! trigger sweep time in seconds
         REAL*4      aperture         ! seconds for coincident station triggers
         CHARACTER*1 algorithm        ! triggering algorithm: x=xdetect, m=mdetect
                                      ! e=eqdetect
         CHARACTER*1 spareJ           ! spare
         INTEGER*2   spareI           ! spare
      END STRUCTURE

      ! SUDS_VELMODEL:  Velocity model
      STRUCTURE /suds_velmodel/
         CHARACTER*4 netname          ! network name
         CHARACTER*6 modelname        ! model name
         CHARACTER*1 spareE
         CHARACTER*1 modeltype        ! p=profile A to B, a=area within corners A B
         REAL*8      latA             ! latitude of point A, north is plus
         REAL*8      longA            ! longitude of point A, east is plus
         REAL*8      latB             ! latitude of point B, north is plus
         REAL*8      longB            ! longitude of point B, east is plus
         INTEGER*4   time_effective   ! time this model was created
      END STRUCTURE

!----------------------------------------------------------------------

      ! Constants for use with the library functions
      REAL*4 SUDS_VERSION
      PARAMETER( SUDS_VERSION = 1.44 )

      INTEGER*2 SUDS_EOF, SUDS_READONLY, SUDS_READWRITE, SUDS_APPEND,
     +          SUDS_CREATE, TRUE, FALSE
      PARAMETER( SUDS_EOF = -1,
     +           SUDS_READONLY = 1,
     +           SUDS_READWRITE = 2,
     +           SUDS_APPEND = 3,
     +           SUDS_CREATE = 4,
     +           TRUE = 1,
     +           FALSE = 0 )

      ! Void time value
      INTEGER*4 NOTIME
      PARAMETER( NOTIME = -2147472000 ) ! Dec 14 1901 00:00:00 GMT

      ! Structure type identifier constants
      INTEGER*2 NO_STRUCT, STAT_IDENT, STRUCTTAG, TERMINATOR, EQUIPMENT,
     +          STATIONCOMP, MUXDATA, DESCRIPTRACE, LOCTRACE,
     +          CALIBRATION, FEATURE, RESIDUAL, EVENT, EV_DESCRIPT,
     +          ORIGIN, S_ERROR, FOCALMECH, MOMENT, VELMODEL, LAYERS,
     +          COMMENT, PROFILE, SHOTGATHER, CALIB, COMPLEX, TRIGGERS,
     +          TRIGSETTING, EVENTSETTING, DETECTOR, ATODINFO,
     +          TIMECORRECTION, INSTRUMENT, CHANSET
      PARAMETER( NO_STRUCT = 0,
     +           STAT_IDENT = 1,
     +           STRUCTTAG = 2,
     +           TERMINATOR = 3,
     +           EQUIPMENT = 4,
     +           STATIONCOMP = 5,
     +           MUXDATA = 6,
     +           DESCRIPTRACE = 7,
     +           LOCTRACE = 8,
     +           CALIBRATION = 9,
     +           FEATURE = 10,
     +           RESIDUAL = 11,
     +           EVENT = 12,
     +           EV_DESCRIPT = 13,
     +           ORIGIN = 14,
     +           S_ERROR = 15,
     +           FOCALMECH = 16,
     +           MOMENT = 17,
     +           VELMODEL = 18,
     +           LAYERS = 19,
     +           COMMENT = 20,
     +           PROFILE = 21,
     +           SHOTGATHER = 22,
     +           CALIB = 23,
     +           COMPLEX = 24,
     +           TRIGGERS = 25,
     +           TRIGSETTING = 26,
     +           EVENTSETTING = 27,
     +           DETECTOR = 28,
     +           ATODINFO = 29,
     +           TIMECORRECTION = 30,
     +           INSTRUMENT = 31,
     +           CHANSET = 32 )

!----------------------------------------------------------------------

      ! SUDS structure used in library function calls
      ! Union provides access to each SUDS structure type
      STRUCTURE /suds_structs/
         INTEGER*2 type        ! Structure type
         INTEGER*4 data_len    ! Length of data following structure (bytes)
         UNION
            MAP
               RECORD /suds_atodinfo/ ad
            END MAP
            MAP
               RECORD /suds_calibration/ ca
            END MAP
            MAP
               RECORD /suds_comment/ co
            END MAP
            MAP
               RECORD /suds_chanset/ cs
            END MAP
            MAP
               RECORD /suds_descriptrace/ dt
            END MAP
            MAP
               RECORD /suds_detector/ de
            END MAP
            MAP
               RECORD /suds_equipment/ eq
            END MAP
            MAP
               RECORD /suds_error/ er
            END MAP
            MAP
               RECORD /suds_event/ ev
            END MAP
            MAP
               RECORD /suds_eventsetting/ es
            END MAP
            MAP
               RECORD /suds_evdescr/ ed
            END MAP
            MAP
               RECORD /suds_feature/ fe
            END MAP
            MAP
               RECORD /suds_focalmech/ fo
            END MAP
            MAP
               RECORD /suds_instrument/ in
            END MAP
            MAP
               RECORD /suds_layers/ la
            END MAP
            MAP
               RECORD /suds_loctrace/ lo
            END MAP
            MAP
               RECORD /suds_moment/ mo
            END MAP
            MAP
               RECORD /suds_muxdata/ mu
            END MAP
            MAP
               RECORD /suds_origin/ or
            END MAP
            MAP
               RECORD /suds_profile/ pr
            END MAP
            MAP
               RECORD /suds_residual/ re
            END MAP
            MAP
               RECORD /suds_shotgather/ sh
            END MAP
            MAP
               RECORD /suds_stationcomp/ sc
            END MAP
            MAP
               RECORD /suds_terminator/ te
            END MAP
            MAP
               RECORD /suds_timecorrection/ tc
            END MAP
            MAP
               RECORD /suds_triggers/ tr
            END MAP
            MAP
               RECORD /suds_trigsetting/ ts
            END MAP
            MAP
               RECORD /suds_velmodel/ vm
            END MAP
         END UNION
      END STRUCTURE
