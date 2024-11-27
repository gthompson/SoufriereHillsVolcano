// SUDS form definition and declaration

// 24-Aug-1993, RLB

#ifndef _SUDS_FORM_
#define _SUDS_FORM_

#include <stdio.h>
#include "suds.h"

// special formats
#define STIFM   ""
#define STTFM   "%ld"
#define MSTFM   "%lf"
#define LTTFM   "%lf"
#define CALFM   ""
#define COMFM   ""

SUDS_FORM _suds_form[] = {
   NO_STRUCT,   "junk",                 STR,    1,    0,     NOSTRG, "%s",     0,
   STAT_IDENT,  "network",              STR,    4,    0,     "unk",  "%-4.4s", 0,
   STAT_IDENT,  "station name",         STR,    5,    4,     NOSTRG, "%-5.5s", 0,
   STAT_IDENT,  "component",            CHR,    1,    9,     NOCHAR, "%c",     0,
   STAT_IDENT,  "instrument type",      SHT,    1,    10,    NOLIST, "%d",     0,
   STRUCTTAG,   "sync char",            CHR,    1,    0,     "S",    "%c",     0,
   STRUCTTAG,   "machine type",         CHR,    1,    1,  "MACHINE", "%c",     0,
   STRUCTTAG,   "struct number",        SHT,    1,    2,     NODATA, "%d",     0,
   STRUCTTAG,   "struct length",        LNG,    1,    4,     NODATA, "%ld",    0,
   STRUCTTAG,   "data length",          LNG,    1,    8,     NODATA, "%ld",    0,
//10
   TERMINATOR,  "beginning struct",     SHT,    1,    0,     NODATA, "%d",     0,
   TERMINATOR,  "unused",               SHT,    1,    2,     NODATA, "%d",     0,
   EQUIPMENT,   "this station/comp",    STI,    1,    0,     NOSTRG, STIFM,    STAT_IDENT,
   EQUIPMENT,   "last station/comp",    STI,    1,    8,     NOSTRG, STIFM,    STAT_IDENT,
   EQUIPMENT,   "next station/comp",    STI,    1,    16,    NOSTRG, STIFM,    STAT_IDENT,
   EQUIPMENT,   "serial number",        STR,    8,    24,    NOSTRG, "%-8.8s", 0,
   EQUIPMENT,   "model",                SHT,    1,    32,    NOLIST, "%d",     0,
   EQUIPMENT,   "knob 1 setting",       SHT,    1,    34,    NODATA, "%d",     0,
   EQUIPMENT,   "knob 2 setting",       SHT,    1,    36,    NODATA, "%d",     0,
   EQUIPMENT,   "repair reason",        SHT,    1,    38,    NOLIST, "%d",     0,
//20
   EQUIPMENT,   "frequency",            FLT,    1,    40,    NODATA, "%f",     0,
   EQUIPMENT,   "date effective",       STT,    1,    44,    NOTIME, STTFM,    0,
   STATIONCOMP, "station/component",    STI,    1,    0,     NOSTRG, STIFM,    STAT_IDENT,
   STATIONCOMP, "component azimuth",    SHT,    1,    12,    NODATA, "%d",     0,
   STATIONCOMP, "component incidence",  SHT,    1,    14,    NODATA, "%d",     0,
   STATIONCOMP, "latitude",             LLT,    1,    16,    NODATA, LTTFM,    0,
   STATIONCOMP, "longitude",            LLT,    1,    24,    NODATA, LTTFM,    0,
   STATIONCOMP, "elevation, meters",    FLT,    1,    32,    NODATA, "%f",     0,
   STATIONCOMP, "enclosure",            CHR,    1,    36,    NOCHAR, "%c",     0,
   STATIONCOMP, "annotated comment",    MIN,    1,    37,    NOLIST, "%d",     0,
//30
   STATIONCOMP, "recorder type",        CHR,    1,    38,    NOCHAR, "%c",     0,
   STATIONCOMP, "rock class",           CHR,    1,    39,    NOCHAR, "%c",     0,
   STATIONCOMP, "rock type",            SHT,    1,    40,    NOLIST, "%d",     0,
   STATIONCOMP, "site condition",       CHR,    1,    42,    NOCHAR, "%c",     0,
   STATIONCOMP, "sensor type",          CHR,    1,    43,    NOCHAR, "%c",     0,
   STATIONCOMP, "data type",            CHR,    1,    44,    NOCHAR, "%c",     0,
   STATIONCOMP, "data units",           CHR,    1,    45,    NOCHAR, "%c",     0,
   STATIONCOMP, "polarity",             CHR,    1,    46,    NOCHAR, "%c",     0,
   STATIONCOMP, "status",               MIN,    1,    47,    NOLIST, "%d",     0,
   STATIONCOMP, "maximum gain",         FLT,    1,    48,    NODATA, "%f",     0,
//40
   STATIONCOMP, "clipping value",       FLT,    1,    52,    NODATA, "%f",     0,
   STATIONCOMP, "conversion to mvolts", FLT,    1,    56,    NODATA, "%f",     0,
   STATIONCOMP, "channel",              SHT,    1,    60,    NODATA, "%d",     0,
   STATIONCOMP, "atod gain",            SHT,    1,    62,    NODATA, "%d",     0,
   STATIONCOMP, "effective date",       STT,    1,    64,    NOTIME, STTFM,    0,
   STATIONCOMP, "clock correction",     FLT,    1,    68,    "0.0",  "%+f",    0,
   STATIONCOMP, "station delay",        FLT,    1,    72,    "0.0",  "%+f",     0,
   MUXDATA,     "network",              STR,    4,    0,     "unk",  "%-4.4s", 0,
   MUXDATA,     "beginning time",       MST,    1,    4,     NOTIME, MSTFM,    0,
   MUXDATA,     "minutes west of GMT",  SHT,    1,    12,    NODATA, "%d",     0,
//50
   MUXDATA,     "number of channels",   SHT,    1,    14,    NODATA, "%d",     0,
   MUXDATA,     "samples per second",   FLT,    1,    16,    NODATA, "%f",     0,
   MUXDATA,     "type of data",         CHR,    1,    20,    NOCHAR, "%c",     0,
   MUXDATA,     "data descriptor",      CHR,    1,    21,    NOCHAR, "%c",     0,
   MUXDATA,     "spareG",               SHT,    1,    22,    NODATA, "%d",     0,
   MUXDATA,     "number of samples",    LNG,    1,    24,    NODATA, "%ld",    0,
   MUXDATA,     "demultiplexed samples",LNG,    1,    28,    NODATA, "%ld",    0,
   DESCRIPTRACE,"station/component",    STI,    1,    0,     NOSTRG, STIFM,    STAT_IDENT,
   DESCRIPTRACE,"initial sample time",  MST,    1,    12,    NOTIME, MSTFM,    0,
   DESCRIPTRACE,"local time diff",      SHT,    1,    20,    NODATA, "%d",     0,
//60
   DESCRIPTRACE,"data type",            CHR,    1,    22,    NOCHAR, "%c",     0,
   DESCRIPTRACE,"data descriptor",      CHR,    1,    23,    NOCHAR, "%c",     0,
   DESCRIPTRACE,"digitized by",         SHT,    1,    24,    NOLIST, "%d",     0,
   DESCRIPTRACE,"processed by",         SHT,    1,    26,    NOLIST, "%d",     0,
   DESCRIPTRACE,"number of samples",    LNG,    1,    28,    NODATA, "%ld",    0,
   DESCRIPTRACE,"samples per second",   FLT,    1,    32,    NODATA, "%f",     0,
   DESCRIPTRACE,"minimum data value",   FLT,    1,    36,    NODATA, "%f",     0,
   DESCRIPTRACE,"maximum data value",   FLT,    1,    40,    NODATA, "%f",     0,
   DESCRIPTRACE,"average noise",        FLT,    1,    44,    NODATA, "%f",     0,
   DESCRIPTRACE,"num clipped samples",  LNG,    1,    48,    NODATA, "%ld",    0,
//70
   DESCRIPTRACE,"time correction",      DBL,    1,    52,    "0.0",  "%+lf",   0,
   DESCRIPTRACE,"rate correction",      FLT,    1,    60,    "0.0",  "%+f",    0,
   LOCTRACE,    "station/component",    STI,    1,    0,     NOSTRG, STIFM,    STAT_IDENT,
   LOCTRACE,    "pathname",             LNG,    1,    8,     NODATA, "%ld",    0,
   LOCTRACE,    "tape name",            LNG,    1,    12,    NODATA, "%ld",    0,
   LOCTRACE,    "offset",               LNG,    1,    16,    NODATA, "%ld",    0,
   CALIBRATION, "station/component",    STI,    1,    0,     NOSTRG, STIFM,    STAT_IDENT,
   CALIBRATION, "normalization factor", FLT,    1,    8,     NODATA, "%f",     0,
   CALIBRATION, "calibration",          CAL,NOCALPTS, 12,    NODATA, CALFM,    CALIB,
   CALIBRATION, "beginning time",       STT,    1,    512,   NOTIME, STTFM,    0,
//80
   CALIBRATION, "ending time",          STT,    1,    516,   NOTIME, STTFM,    0,
   FEATURE,     "station/component",    STI,    1,    0,     NOSTRG, STIFM,    STAT_IDENT,
   FEATURE,     "observed phase",       SHT,    1,    12,    NOLIST, "%d",     0,
   FEATURE,     "onset descriptor",     CHR,    1,    14,    NOCHAR, "%c",     0,
   FEATURE,     "first motion",         CHR,    1,    15,    NOCHAR, "%c",     0,
   FEATURE,     "signal/noise",         SHT,    1,    16,    NODATA, "%d",     0,
   FEATURE,     "data source",          CHR,    1,    18,    NOCHAR, "%c",     0,
   FEATURE,     "timing quality",       CHR,    1,    19,    NOCHAR, "%c",     0,
   FEATURE,     "amplitude quality",    CHR,    1,    20,    NOCHAR, "%c",     0,
   FEATURE,     "amplitude units",      CHR,    1,    21,    NOCHAR, "%c",     0,
//90
   FEATURE,     "gain range factor",    SHT,    1,    22,    NODATA, "%d",     0,
   FEATURE,     "phase time",           MST,    1,    24,    NOTIME, MSTFM,    0,
   FEATURE,     "phase amplitude",      FLT,    1,    32,    NODATA, "%f",     0,
   FEATURE,     "phase period",         FLT,    1,    36,    NODATA, "%f",     0,
   FEATURE,     "time of pick",         STT,    1,    40,    NOTIME, STTFM,    0,
   FEATURE,     "pick authority",       SHT,    1,    44,    NOLIST, "%d",     0,
   FEATURE,     "pick reader",          SHT,    1,    46,    NODATA, "%d",     0,
   RESIDUAL,    "event number",         LNG,    1,    0,     NODATA, "%ld",    0,
   RESIDUAL,    "station/component",    STI,    1,    4,     NOSTRG, STIFM,    STAT_IDENT,
   RESIDUAL,    "phase code",           SHT,    1,    12,    NOLIST, "%d",     0,
//100
   RESIDUAL,    "timing quality",       CHR,    1,    14,    NOCHAR, "%c",     0,
   RESIDUAL,    "amplitude quality",    CHR,    1,    15,    NOCHAR, "%c",     0,
   RESIDUAL,    "traveltime residual",  FLT,    1,    16,    NODATA, "%f",     0,
   RESIDUAL,    "solution weight",      FLT,    1,    20,    NODATA, "%f",     0,
   RESIDUAL,    "delay time",           FLT,    1,    24,    NODATA, "%f",     0,
   RESIDUAL,    "azimuth to station",   FLT,    1,    28,    NODATA, "%f",     0,
   RESIDUAL,    "distance to station",  FLT,    1,    32,    NODATA, "%f",     0,
   RESIDUAL,    "angle of emergence",   FLT,    1,    36,    NODATA, "%f",     0,
   EVENT,       "organization",         SHT,    1,    0,     NOLIST, "%d",     0,
   EVENT,       "event number",         LNG,    1,    2,     NODATA, "%ld",    0,
//110
   EVENT,       "number felt reports",  SHT,    1,    6,     NODATA, "%d",     0,
   EVENT,       "MM intensity",         CHR,    1,    8,     NOCHAR, "%c",     0,
   EVENT,       "event type",           CHR,    1,    9,     NOCHAR, "%c",     0,
   EVENT,       "tectonism",            CHR,    1,    10,    NOCHAR, "%c",     0,
   EVENT,       "waterwaves",           CHR,    1,    11,    NOCHAR, "%c",     0,
   EVENT,       "mechanism type",       CHR,    1,    12,    NOCHAR, "%c",     0,
   EVENT,       "explosive medium",     CHR,    1,    13,    NOCHAR, "%c",     0,
   EVENT,       "magnitude or lbs TNT", FLT,    1,    14,    NODATA, "%f",     0,
   EV_DESCRIPT, "earthquake name",      STR,    20,   0,     NOSTRG, "%-20.20s",0,
   EV_DESCRIPT, "county",               STR,    16,   20,    NOSTRG, "%16.16s", 0,
//120
   EV_DESCRIPT, "state",                STR,    16,   36,    NOSTRG, "%-16.16s",0,
   EV_DESCRIPT, "local time diff",      SHT,    1,    52,    NODATA, "%d",     0,
   EV_DESCRIPT, "unused",               SHT,    1,    54,    NODATA, "%d",     0,
   ORIGIN,      "event number",         LNG,    1,    0,     NODATA, "%ld",    0,
   ORIGIN,      "organization",         SHT,    1,    4,     NOLIST, "%d",     0,
   ORIGIN,      "version",              CHR,    1,    6,     NOCHAR, "%c",     0,
   ORIGIN,      "processing status",    CHR,    1,    7,     NOCHAR, "%c",     0,
   ORIGIN,      "preferred location",   CHR,    1,    8,     NOCHAR, "%c",     0,
   ORIGIN,      "processing program",   CHR,    1,    9,     NOCHAR, "%c",     0,
   ORIGIN,      "depth control",        CHR,    1,    10,    NOCHAR, "%c",     0,
//130
   ORIGIN,      "convergence",          CHR,    1,    11,    NOCHAR, "%c",     0,
   ORIGIN,      "region",               LNG,    1,    12,    NODATA, "%ld",    0,
   ORIGIN,      "origin time",          MST,    1,    16,    NOTIME, MSTFM,    0,
   ORIGIN,      "latitude",             LLT,    1,    24,    NODATA, LTTFM,    0,
   ORIGIN,      "longitude",            LLT,    1,    32,    NODATA, LTTFM,    0,
   ORIGIN,      "depth, km",            FLT,    1,    40,    NODATA, "%f",     0,
   ORIGIN,      "horizontalerror",      FLT,    1,    44,    NODATA, "%f",     0,
   ORIGIN,      "error in depth",       FLT,    1,    48,    NODATA, "%f",     0,
   ORIGIN,      "rms of residuals",     FLT,    1,    52,    NODATA, "%f",     0,
   ORIGIN,      "crustal model used",   STR,    6,    56,    NOSTRG, "%-6.6s", 0,
//140
   ORIGIN,      "azimuthal gap",        SHT,    1,    62,    NODATA, "%d",     0,
   ORIGIN,      "closest station, km",  FLT,    1,    64,    NODATA, "%f",     0,
   ORIGIN,      "number of stations",   SHT,    1,    68,    NODATA, "%d",     0,
   ORIGIN,      "number p phases rept.",SHT,    1,    70,    NODATA, "%d",     0,
   ORIGIN,      "number p phases used", SHT,    1,    72,    NODATA, "%d",     0,
   ORIGIN,      "number s phases rept", SHT,    1,    74,    NODATA, "%d",     0,
   ORIGIN,      "number s phases used", SHT,    1,    76,    NODATA, "%d",     0,
   ORIGIN,      "magnitude type",       SHT,    1,    78,    NODATA, "%d",     0,
   ORIGIN,      "number mag readings",  SHT,    1,    80,    NODATA, "%d",     0,
   ORIGIN,      "number mag s used",    SHT,    1,    82,    NODATA, "%d",     0,
//150
   ORIGIN,      "magnitude",            FLT,    1,    84,    NODATA, "%f",     0,
   ORIGIN,      "magnitude weight",     FLT,    1,    88,    NODATA, "%f",     0,
   ORIGIN,      "magnitude rms",        FLT,    1,    92,    NODATA, "%f",     0,
   ORIGIN,      "date solution done",   STT,    1,    96,    NOTIME, STTFM,    0,
   ERROR,       "covariance matrix",    FLT,    10,   0,     NODATA, "%f",     0,
   FOCALMECH,   "strike plane a",       FLT,    1,    0,     NODATA, "%f",     0,
   FOCALMECH,   "dip plane a",          FLT,    1,    4,     NODATA, "%f",     0,
   FOCALMECH,   "rake plane a",         FLT,    1,    8,     NODATA, "%f",     0,
   FOCALMECH,   "strike plane b",       FLT,    1,    12,    NODATA, "%f",     0,
   FOCALMECH,   "dip plane b",          FLT,    1,    16,    NODATA, "%f",     0,
//160
   FOCALMECH,   "rake plane b",         FLT,    1,    20,    NODATA, "%f",     0,
   FOCALMECH,   "preferred plane",      CHR,    1,    24,    NOCHAR, "%c",     0,
   FOCALMECH,   "unused",               CHR,    3,    25,    NOCHAR, "%c",     0,
   MOMENT,      "datatypes",            BTS,    1,    0,     "0",    "%x",     0,
   MOMENT,      "constraints",          CHR,    1,    1,     NOCHAR, "%c",     0,
   MOMENT,      "unused",               CHR,    2,    2,     NOCHAR, "%c",     0,
   MOMENT,      "scalar moment",        FLT,    1,    4,     NODATA, "%f",     0,
   MOMENT,      "moment tensor",        FLT,    6,    8,     NODATA, "%f",     0,
   VELMODEL,    "network name",         STR,    4,    0,     NOSTRG, "%-4.4s", 0,
   VELMODEL,    "model name",           STR,    6,    4,     NOSTRG, "%-6.6s", 0,
//170
   VELMODEL,    "unused",               CHR,    1,    10,    NOCHAR, "%c",     0,
   VELMODEL,    "type",                 CHR,    1,    11,    NOCHAR, "%c",     0,
   VELMODEL,    "latitude pt A",        LLT,    1,    12,    NODATA, LTTFM,    0,
   VELMODEL,    "longitude pt A",       LLT,    1,    20,    NODATA, LTTFM,    0,
   VELMODEL,    "latitude pt B",        LLT,    1,    28,    NODATA, LTTFM,    0,
   VELMODEL,    "longitude pt B",       LLT,    1,    36,    NODATA, LTTFM,    0,
   VELMODEL,    "time effective",       STT,    1,    44,    NOTIME, STTFM,    0,
   LAYERS,      "thickness, km",        FLT,    1,    0,     NODATA, "%f",     0,
   LAYERS,      "p velocity at top",    FLT,    1,    4,     NODATA, "%f",     0,
   LAYERS,      "p velocity at base",   FLT,    1,    8,     NODATA, "%f",     0,
//180
   LAYERS,      "s velocity at top",    FLT,    1,    12,    NODATA, "%f",     0,
   LAYERS,      "s velocity at base",   FLT,    1,    16,    NODATA, "%f",     0,
   LAYERS,      "velocity function",    SHT,    1,    20,    NODATA, "%d",     0,
   LAYERS,      "spare",                SHT,    1,    24,    NODATA, "%d",     0,
   COMMENT,     "structure identifier", SHT,    1,    0,     NODATA, "%d",     0,
   COMMENT,     "field number",         SHT,    1,    2,     NODATA, "%d",     0,
   COMMENT,     "comment length",       SHT,    1,    4,     NODATA, "%d",     0,
   COMMENT,     "unused",               SHT,    1,    6,     NODATA, "%d",     0,
   PROFILE,     "not right",            SHT,    1,    6,     NODATA, "%d",     0,
   SHOTGATHER,  "not right",            SHT,    1,    6,     NODATA, "%d",     0,
//190
   CALIB,       "pole",                 CPX,    1,    0,     NODATA, COMFM,    COMPLEX,
   CALIB,       "zero",                 CPX,    1,    4,     NODATA, COMFM,    COMPLEX,
   COMPLEX,     "real",                 FLT,    1,    0,     NODATA, "%f",     0,
   COMPLEX,     "imaginary",            FLT,    1,    4,     NODATA, "%f",     0,
   TRIGGERS,    "station/component",    STI,    1,    0,     NOSTRG, STIFM,    STAT_IDENT,
   TRIGGERS,    "short-term average",   SHT,    1,    12,    NODATA, "%d",     0,
   TRIGGERS,    "long-term average",    SHT,    1,    14,    NODATA, "%d",     0,
   TRIGGERS,    "sh-term abs average",  SHT,    1,    16,    NODATA, "%d",     0,
   TRIGGERS,    "lg-term abs average",  SHT,    1,    18,    NODATA, "%d",     0,
   TRIGGERS,    "trigger level",        SHT,    1,    20,    NODATA, "%d",     0,
//200
   TRIGGERS,    "number of triggers",   SHT,    1,    22,    NODATA, "%d",     0,
   TRIGGERS,    "time of first trigger",MST,    1,    24,    NOTIME, MSTFM,    0,
   TRIGSETTING, "network name",         STR,    4,    0,     "unk",  "%-4.4s", 0,
   TRIGSETTING, "last mod time",        MST,    1,    4,     NOTIME, MSTFM,    0,
   TRIGSETTING, "const1",               SHT,    1,    12,    NODATA, "%d",     0,
   TRIGSETTING, "const2",               SHT,    1,    14,    NODATA, "%d",     0,
   TRIGSETTING, "threshold",            SHT,    1,    16,    NODATA, "%d",     0,
   TRIGSETTING, "const3",               SHT,    1,    18,    NODATA, "%d",     0,
   TRIGSETTING, "const4",               SHT,    1,    20,    NODATA, "%d",     0,
   TRIGSETTING, "weighted ave increment",SHT,   1,    22,    NODATA, "%d",     0,
//210
   TRIGSETTING, "sweep in seconds",     FLT,    1,    24,    NODATA, "%f",     0,
   TRIGSETTING, "aperture",             FLT,    1,    28,    NODATA, "%f",     0,
   TRIGSETTING, "algorithm",            CHR,    1,    32,    NOCHAR, "%c",     0,
   TRIGSETTING, "spare",                CHR,    1,    33,    NOCHAR, "%c",     0,
   TRIGSETTING, "spare",                SHT,    1,    34,    NODATA, "%d",    0,
   EVENTSETTING, "network name",        STR,    4,    0,     "unk",  "%-4.4s", 0,
   EVENTSETTING, "last mod time",       MST,    1,    4,     NOTIME, MSTFM,    0,
   EVENTSETTING, "const1",              SHT,    1,    12,    NODATA, "%d",     0,
   EVENTSETTING, "const2",              SHT,    1,    14,    NODATA, "%d",     0,
   EVENTSETTING, "threshold",           SHT,    1,    16,    NODATA, "%d",     0,
//220
   EVENTSETTING, "const3",              SHT,    1,    18,    NODATA, "%d",     0,
   EVENTSETTING, "minduration",         FLT,    1,    20,    NODATA, "%f",     0,
   EVENTSETTING, "maxduration",         FLT,    1,    24,    NODATA, "%f",     0,
   EVENTSETTING, "algorithm",           CHR,    1,    28,    NOCHAR, "%c",     0,
   EVENTSETTING, "spare",               CHR,    1,    29,    NOCHAR, "%c",     0,
   EVENTSETTING, "spare",               SHT,    1,    30,    NODATA, "%d",    0,
   DETECTOR,    "algorithm",            CHR,    1,    0,     NOCHAR, "%c",     0,
   DETECTOR,    "event type",           CHR,    1,    1,     NOCHAR, "%c",     0,
   DETECTOR,    "network node id",      STR,    10,   2,     "unk",  "%-10.10s",0,
   DETECTOR,    "version number",       FLT,    1,    12,    NODATA, "%f",     0,
//230
   DETECTOR,    "event number",         LNG,    1,    16,    NODATA, "%ld",    0,
   DETECTOR,    "spare",                LNG,    1,    20,    NODATA, "%ld",    0,
   ATODINFO,    "base address",         SHT,    1,    0,     NODATA, "%d",     0,
   ATODINFO,    "device id",            SHT,    1,    2,     NODATA, "%d",     0,
   ATODINFO,    "device flags",         BTW,    1,    4,     "0",    "%u",     0,
   ATODINFO,    "# of extended bufs",   SHT,    1,    6,     NODATA, "%d",     0,
   ATODINFO,    "external mux control word",SHT,1,    8,     NODATA, "%d",     0,
   ATODINFO,    "AtoD timing source",   CHR,    1,    10,    NOCHAR, "%c",     0,
   ATODINFO,    "AtoD trigger source",  CHR,    1,    11,    NOCHAR, "%c",     0,
   TIMECORRECTION,"station/component",  STI,    1,    0,     NOSTRG, STIFM,    STAT_IDENT,
//240
   TIMECORRECTION,"time correction",    DBL,    1,    12,    "0.0",  "%+lf",   0,
   TIMECORRECTION,"rate correction",    FLT,    1,    20,    "0.0",  "%+f",    0,
   TIMECORRECTION,"synchronization code", CHR,  1,    24,    NOCHAR, "%c",     0,
   TIMECORRECTION,"program code",       CHR,    1,    25,    NOCHAR, "%c",     0,
   TIMECORRECTION,"effective time",     STT,    1,    26,    NOTIME, STTFM,    0,
   INSTRUMENT,  "station/component",    STI,    1,    0,     NOSTRG, STIFM,    STAT_IDENT,
   INSTRUMENT,  "inst. serial number",  SHT,    1,    12,    NODATA, "%d",     0,
   INSTRUMENT,  "number of components", SHT,    1,    14,    NODATA, "%d",     0,
   INSTRUMENT,  "channel number",       SHT,    1,    16,    NODATA, "%d",     0,
   INSTRUMENT,  "sensor type",          CHR,    1,    18,    NOCHAR, "%c",     0,
//250
   INSTRUMENT,  "data type",            CHR,    1,    19,    NOCHAR, "%c",     0,
   INSTRUMENT,  "void sample value",    LNG,    1,    20,    NODATA, "%ld",    0,
   INSTRUMENT,  "digitizing constant",  FLT,    1,    24,    NODATA, "%f",     0,
   INSTRUMENT,  "AAF corner freq. (Hz)",FLT,    1,    28,    NODATA, "%f",     0,
   INSTRUMENT,  "AAF poles",            FLT,    1,    32,    NODATA, "%f",     0,
   INSTRUMENT,  "trans natural freq. (Hz)",FLT, 1,    36,    NODATA, "%f",     0,
   INSTRUMENT,  "trans damping coeff.", FLT,    1,    40,    NODATA, "%f",     0,
   INSTRUMENT,  "trans motion constant",FLT,    1,    44,    NODATA, "%f",     0,
   INSTRUMENT,  "amplifier gain (dB)",  FLT,    1,    48,    NODATA, "%f",     0,
   INSTRUMENT,  "local X coord. (meters)",FLT,  1,    52,    NODATA, "%f",     0,
//260
   INSTRUMENT,  "local Y coord. (meters)",FLT,  1,    56,    NODATA, "%f",     0,
   INSTRUMENT,  "local Z coord. (meters)",FLT,  1,    60,    NODATA, "%f",     0,
   INSTRUMENT,  "effective time",       STT,    1,    64,    NOTIME, STTFM,    0,
   INSTRUMENT,  "pre-event memory (secs)",FLT,  1,    68,    NODATA, "%f",     0,
   INSTRUMENT,  "trigger number",       SHT,    1,    72,    NODATA, "%d",     0,
   INSTRUMENT,  "study ID",             STR,    6,    74,    NOSTRG, "%s",     0,
   INSTRUMENT,  "sensor serial number", SHT,    1,    80,    NODATA, "%d",     0,
   CHANSET,     "set type",             SHT,    1,    0,     "0",    "%d",     0,
   CHANSET,     "entries in set",       SHT,    1,    2,     "0",    "%d",     0,
   CHANSET,     "network name",         STR,    4,    4,     NOSTRG, "%s",     0,
//270
   CHANSET,     "set name",             STR,    5,    8,     NOSTRG, "%s",     0,
   CHANSET,     "active time",          STT,    1,    14,    NOTIME, STTFM,    0,
   CHANSET,     "inactive time",        STT,    1,    18,    NOTIME, STTFM,    0
};

// Indicies in _suds_form for beginning of each structure type
int _suds_beg_struct[] = {
   0,
   1,    // SUDS_STATIDENT
   5,    // SUDS_STRUCTTAG
   10,   // SUDS_TERMINATOR
   12,   // SUDS_EQUIPMENT
   22,   // SUDS_STATIONCOMP
   47,   // SUDS_MUXDATA
   57,   // SUDS_DESCRIPTRACE
   72,   // SUDS_LOCTRACE
   76,   // SUDS_CALIBRATION
   81,   // SUDS_FEATURE
   97,   // SUDS_RESIDUAL
   108,  // SUDS_EVENT
   118,  // SUDS_EVDESCR
   123,  // SUDS_ORIGIN
   154,  // SUDS_ERROR
   155,  // SUDS_FOCALMECH
   163,  // SUDS_MOMENT
   168,  // SUDS_VELMODEL
   177,  // SUDS_LAYERS
   184,  // SUDS_COMMENT
   188,  // SUDS_PROFILE
   189,  // SUDS_SHOTGATHER
   190,  // SUDS_CALIBR
   192,  // COMPLEXX
   194,  // SUDS_TRIGGERS
   202,  // SUDS_TRIGSETTING
   215,  // SUDS_EVENTSETTING
   226,  // SUDS_DETECTOR
   232,  // SUDS_ATODINFO
   239,  // SUDS_TIMECORRECTION
   245,  // SUDS_INSTRUMENT
   267   // SUDS_CHANSET
};

int _suds_size_struct[] = {
   0,
   sizeof( SUDS_STATIDENT ),
   sizeof( SUDS_STRUCTTAG ),
   sizeof( SUDS_TERMINATOR ),
   sizeof( SUDS_EQUIPMENT ),
   sizeof( SUDS_STATIONCOMP ),
   sizeof( SUDS_MUXDATA ),
   sizeof( SUDS_DESCRIPTRACE ),
   sizeof( SUDS_LOCTRACE ),
   sizeof( SUDS_CALIBRATION ),
   sizeof( SUDS_FEATURE ),
   sizeof( SUDS_RESIDUAL ),
   sizeof( SUDS_EVENT ),
   sizeof( SUDS_EVDESCR ),
   sizeof( SUDS_ORIGIN ),
   sizeof( SUDS_ERROR ),
   sizeof( SUDS_FOCALMECH ),
   sizeof( SUDS_MOMENT ),
   sizeof( SUDS_VELMODEL ),
   sizeof( SUDS_LAYERS ),
   sizeof( SUDS_COMMENT ),
   sizeof( SUDS_PROFILE ),
   sizeof( SUDS_SHOTGATHER ),
   sizeof( SUDS_CALIBR ),
   sizeof( COMPLEXX ),
   sizeof( SUDS_TRIGGERS ),
   sizeof( SUDS_TRIGSETTING ),
   sizeof( SUDS_EVENTSETTING ),
   sizeof( SUDS_DETECTOR ),
   sizeof( SUDS_ATODINFO ),
   sizeof( SUDS_TIMECORRECTION ),
   sizeof( SUDS_INSTRUMENT ),
   sizeof( SUDS_CHANSET )
};

SUDS_CODES _suds_struct_names[]={
   0,   "No Struct",
   1,   "StationIdentifier",
   2,   "Structure Tag",
   3,   "Terminator",
   4,   "Equipment",
   5,   "StationComp",
   6,   "MuxData",
   7,   "DescripTrace",
   8,   "LocTrace",
   9,   "Calibration",
   10,  "Feature",
   11,  "Residual",
   12,  "Event",
   13,  "EventDescription",
   14,  "Origin",
   15,  "Error",
   16,  "FocalMech",
   17,  "Moment",
   18,  "VelModel",
   19,  "Layers",
   20,  "Comment",
   21,  "Profile",
   22,  "ShotGather",
   23,  "Calib/Points",
   24,  "Complex Number",
   25,  "Triggers",
   26,  "TrigSetting",
   27,  "EventSetting",
   28,  "Detector",
   29,  "AtoDinfo",
   30,  "TimeCorrection",
   31,  "Instrument",
   32,  "ChannelSet",
   0,   0
};

#endif
