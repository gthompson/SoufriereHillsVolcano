                          CONVSEIS v6.0

    Conversion program package for seismological data on PCs.

                 by Lani Oncescu & Mike Rizescu

                            (May, 2001)


Converters:

ASAC2EQ  : from ASCII-SAC into PCEQ format, including despiking and
           decimation;
ASAC2SGA : from ASCII-SAC into SeisGram-ASCII format;
ASAC2SUD : from ASCII-SAC into PC-SUDS demultiplexed format;
ASC2ASAC : from PITSA ASCII into ASCII-SAC format;
ASC2SUDS : from PITSA ASCII into PC-SUDS demultiplexed format;
ASC2V1   : from a general ASCII format into KINEMETRICS - VOL1 format;
CSS2GSE  : from binary CSS3.0 into GSE1.0 (uncompressed) format;
CSS2SGA  : from CSS3.0 into SeisGram-ASCII format;
CSS2SUD  : from CSS3.0 into PC-SUDS demultiplexed format;
EQ2ASC   : from PCEQ into a general ASCII format (with option for GRAPHER
           or GMT plotting program);
EQ2GSE   : from PCEQ into GSE1.0 (uncompressed) format;
EQ2SEGY  : from PCEQ into PC-SEGY format (used by RecSec program);
EQ2SEIS  : from PCEQ into SEIS89 input format;
EQ2SG    : from PCEQ into SeisGram-binary format;
EQ2SUDS  : from PCEQ into PC-SUDS demultiplexed format;
GSE2SG   : from GSE1.0 (uncompressed) into SeisGram-binary format;
GSE2SUDS : from GSE1.0 (uncompressed) into PC-SUDS demultiplexed format;
K22SUDS  : from K2-ASCII and K2-header formats (output of KINEMETRICS-
           KW2ASC) into PC-SUDS format; batch K2SUDS combines both programs;
SDP2EQ   : from Teledyne Geotech SDP (on PDP minicomputers) into PCEQ
           format;
SEGY2ASC : from PC-SEGY into a general ASCII format (with option for
           GRAPHER plotting program);
SG2EQ    : from SeisGram-binary into PCEQ format;
SG2SEIS  : from SeisGram-binary format into SEIS89 input format;
SG2SUDS  : from SeisGram format into PC-SUDS format;
SGA2ASC  : from SeisGram-ASCII into a general ASCII format (with option for
           GRAPHER plotting program);
SGA2EQ   : from SeisGram-ASCII into PCEQ format;
SUD2ASAC : from demultiplexed PC-SUDS format into ASCII-SAC format;
SUDS2CSS : from demultiplexed PC-SUDS format into CSS3.0 format;
SUDS2EQ  : from demultiplexed PC-SUDS format into PCEQ format;
SUDS2GSE : from demultiplexed PC-SUDS into GSE1.0 (uncompressed) format;
SUDS2SGA : from demultiplexed PC-SUDS into SeisGram-ASCII format;
V12SGA   : from KINEMETRICS -VOL1 into SeisGram-ASCII format;
V22EQ    : from KINEMETRICS -VOL2 into PCEQ format;
WV2ASC   : from IASPEI-MDETECT or XDETECT into a general ASCII
           format (with option for GRAPHER plotting program);
WV2EQ    : from IASPEI-MDETECT or XDETEXT  into PCEQ format;
WV2GSE   : from IASPEI-MDETECT or XDETECT into GSE1.0 (uncompressed)
           format.



Utilities:


EVL2SUDS : writing results of location program HYPOPLUS into the
           SUDS_ORIGIN structure of a PC-SUDS file;
FILLEQ   : filling one missing component with zeroes and reentering header
           information (on PCEQ format);
FILTSUDS : filtering, decimating, integrating and skipping channels (on PC-
           SUDS format).
FIXFNAME : coding/decoding filenames produced by K22SUDS converter;
GSE2GSE  : compressing or decompressing GSE1.0 format, converting from
           GSE1.0 to GSE2.0 and vice-versa;
HYPO2GSE : writing phase data or HYPOPLUS location results in GSE2.0 format;
RENCHGSE : renaming channels on GSE1.0 or GSE2.0 files;
RFEQ     : filtering and resampling (on PCEQ format);
SG2SG    : flipping one component, rotating horizontal components and
           computing the product between vertical and radial components
           (on SeisGram-binary format);
SGA2SGA  : summing up SeisGram-ASCII files;
SGY2SGY  : summing up PC-SEGY files;
SORTEQ   : sorting out given channels (on PCEQ format);
SUDSINST : extracting the instrument related information from a PC-SUDS
           file;
SUDSRW   : extracting data from a PC-SUDS file;
SUDSTOOL : listing and interactive changing of channel name, starting time
           and  sampling rate on a PC-SUDS file; it also lists sensor type,
           motion type, data units and content of the SUDS_COMMENT
           structure;
SUDXTR   : extracting phase information (P, S, coda, amplitude and period)
           from a PC-SUDS file;
SUMEQ    : summing up PCEQ files;
TIM2SEGY : generating travel time curves for an arbitrarily layered
           structure, with the output in PC-SEGY format.
V22V2    : extracting a portion of a KINEMETRICS - VOL2 file.



** The output of FILTSUDS (in PC-SUDS format) is readable by SUDS v2.8;
                it now supports real*4 data **

** GSE2GSE supports now conversion into and from GSE2.0 format **

Started with a question-mark on the command line, the programs will produce
a help message !












