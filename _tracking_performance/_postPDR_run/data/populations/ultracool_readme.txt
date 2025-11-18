ultracool.out

From Quinn querying keck to find companions to ultra cool sheet targets that came from Adam burgasser:
""Here’s the result of running this through the Keck LGS AO tip/tilt finder IDL code.  An example of the output is here:

SDSS J000013    00 00 13.500 +25 54 19.80 2000.0 kmag=14.83   lgs=1
   1159-0000095 00 00 14.451 +25 54 52.63 2000.0 rmag=16.0 sep=35.2 b-v=0.41 b-r=0.95 S=0.24
   1159-0000085 00 00 13.028 +25 54 34.93 2000.0 rmag=16.9 sep=16.4 b-v=-0.12 b-r=-0.21 S=0.23

The ultra cool objects are those that have the “lgs=1” flag at the end of the row.  Underneath them will be up to 3 suitable TT stars that were found in the USNO-B1.0 catalog.  If there are no sources underneath a target with an “lgs=1” it means no stars were found with R<18 with a separation <60”.  The tool is not smart enough to know when the TT star is the same as the science target.  Usually if a TT star has a separation <a few arc seconds, it’s probably the same source as the target itself.  

For the identified TT stars, most of the info after the USNO B1.0 identifier is hopefully self-explanatory. The “S” at the end is the predicted Strehl that will be achieved for the target given the properties of the TT star.  These are older estimates and I don’t know how HAKA will factor in to changing the Strehl.

Hope this is usable!""
