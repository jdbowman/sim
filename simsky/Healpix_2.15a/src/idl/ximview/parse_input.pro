; -----------------------------------------------------------------------------
;
;  Copyright (C) 2007-2008   J. P. Leahy
;
;
;  This file is part of Ximview and of HEALPix
;
;  Ximview and HEALPix are free software; you can redistribute them and/or modify
;  them under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2 of the License, or
;  (at your option) any later version.
;
;  Ximview and HEALPix are distributed in the hope that they will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with HEALPix; if not, write to the Free Software
;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;
; -----------------------------------------------------------------------------
; Module PARSE_INPUT: decodes input parameter and converts to image
;
; J. P. Leahy 2008
;
; Contents:
;:
; parse_input_event  Exists solely to delete the optional plot.
; bad2nan            Converts HP bad values to NaN
; hp_check           Sees if data is HEALPix
; get_image          Wrapper that reads both primary HDU and extensions
; get_fits           Opens, decodes and reads a FITS file
; parse_input        Main routine
;
; Standard IDL documentation block follows declaration of parse_input
;
PRO parse_input_event, event
WIDGET_CONTROL, event.TOP, /DESTROY
END

PRO bad2nan, array
; Sets values of array that have HEALPix bad values to NaN
;
COMPILE_OPT IDL2, HIDDEN
hpx_sbadval = -1.6375E30
hpx_dbadval = -1.6375D30

SWITCH SIZE(array, /TYPE) OF
    4:                          ; Single precision
    6: BEGIN                    ; complex single
        flag = hpx_sbadval
        nan  = !values.F_NAN
        BREAK
    END
    5:                          ; Double precision
    9: BEGIN                    ; complex double
        flag = hpx_dbadval
        nan  = !values.D_NAN
        BREAK
    END
    ELSE: RETURN
ENDSWITCH

bad = WHERE(array EQ flag)
IF bad[0] NE -1 THEN array[bad] = nan

END

PRO hp_check, order, header, dsize, healpix, cut4, nside
; Checks header and data size to see if image could be one or several
; Healpix datasets (including CUT4 format
;
; Inputs:
;    Order:  Parse_Input input parameter: if set implies HEALPix.
;    header: Image header
;    dsize:  SIZE structure for dataset
;
; Outputs:
;    Healpix: True if dataset is healpix (including CUT4)
;    Cut4:    True if cut4 format
;    Nside:   HEALPix size parameter
;
COMPILE_OPT IDL2, HIDDEN

healpix = KEYWORD_SET(order)
cut4 = 0B

pixtype = SXPAR(header, 'PIXTYPE', COUNT = count)
hphead  = STRCMP(pixtype, 'HEALPIX', 7, /FOLD_CASE)
healpix = healpix || hphead
IF (count GT 0) && healpix && ~hphead THEN MESSAGE, /INFORMATIONAL, $
          'Assuming HEALPix per input parameters, despite header'

IF healpix || count EQ 0 THEN BEGIN
    nside = SXPAR(header, 'NSIDE')
    nsh = NPIX2NSIDE(dsize.DIMENSIONS[0])
    IF nside EQ 0 THEN nside = nsh
ENDIF
IF nside NE nsh THEN BEGIN
    IF healpix THEN BEGIN       ; Check for cut4 file
                                ; Check if we have required columns
        test = SXPAR(header, 'TTYPE1')
        cut4 = STRCMP(test, 'PIXEL', 5, /FOLD_CASE)
        test = SXPAR(header, 'TTYPE2')
        cut4 = cut4 && STRCMP(test,'SIGNAL', 5, /FOLD_CASE)
                                ; Check for required values of
                                ; optional keywords (if present)
        test = SXPAR(header, 'GRAIN', COUNT = count)
        IF count GT 0 THEN cut4 = cut4 && test EQ 1
        test = SXPAR(header, 'INDXSCHM', COUNT = count) ;
        IF count GT 0 THEN cut4 = $
          cut4 && STRCMP(test, 'EXPLICIT', 8, /FOLD_CASE)
        test = SXPAR(header, 'OBS_NPIX')
        IF test NE 0 THEN cut4 = cut4 && dsize.DIMENSIONS[0] EQ test
        IF ~cut4 THEN nside = -1
    ENDIF ELSE nside = -1       ; could be gridded HEALPix
ENDIF
IF healpix && nside EQ -1 THEN MESSAGE, $
  'Image was claimed to be HEALPix but number of pixels is invalid'

healpix = healpix || (count EQ 0 && nside GT 0)

END

PRO get_image, file, header, naxis, col, verbose, start, exten, data
;
; Wrapper to read image: covers both primary HDU and extensions.
;
COMPILE_OPT IDL2, HIDDEN

nax = SXPAR(header,'NAXIS*')
ntot = 1
FOR i = 2,naxis-1 DO ntot *= nax[i]

ncol = N_ELEMENTS(col)
IF ncol EQ 0 THEN BEGIN
    ncol = ntot
    col = INDGEN(ntot) + 1
ENDIF
                                ; Check that requested plane(s) exists
IF ntot LT MAX(col) THEN MESSAGE, 'Requested image slice(s) not found'

doall = naxis LE 2

IF doall THEN data = READFITS(file, header, EXTEN_NO = exten) ELSE BEGIN
                                ; Read first slice:
    data = READFITS(file, header, NSLICE = col[0]-1, EXTEN_NO = exten)
    IF ncol GT 1 THEN BEGIN
                                ; Use a pointer array for the slices
        temp = TEMPORARY(data)
        data = PTRARR(ncol, /ALLOCATE_HEAP)
        *data[0] = TEMPORARY(temp)
                                ; This is a poor way to access the
                                ; data but at least avoids unnecessary
                                ; virtual memory usage. Really
                                ; lower-level I/O would be better.
        FOR i=1,ncol-1 DO $
          *data[i] = READFITS(file, header, NSLICE = col[i]-1, EXTEN_NO=exten)
    ENDIF
ENDELSE

text = exten GT 0 ? 'Image extension' : 'file'
IF verbose THEN MESSAGE, /INFORMATIONAL, text + ' read at' + $
  STRING(SYSTIME(1) - start)
IF verbose THEN HELP, /MEMORY

END

PRO get_fits, name, col, exten, file, data, header, howto, xtype, $
              verbose, start
;
; Reads a file named "name", "name.fits" or "name.FITS" and extracts
; either a set of slices from an image or a set of columns from a
; table, according to which it encounters first. It only looks in the
; main HDU and the first extension at present.
;
; Inputs:
;    name:   file name or root
; Input/Output
;    col:    list of required image planes or data columns
;            set to list in file if unset on input
;    exten:  FITS extension. If unset look in 0, then 1
; Outputs
;    file:    completed and trimmed file name
;    data:    binary data, or pointer(s) to it
;    header:  FITS header constructed from primary header + extension
;             header (if appropriate).
;    howto:   data format... see header for parse_input
;    xtype:   string containing either "IMAGE" or "BINTABLE"
;    verbose: print diagnostics
;    start:   time (for diagnostics)
;
COMPILE_OPT IDL2, HIDDEN

ncol = N_ELEMENTS(col)
col_set = ncol GT 0

; see if there is a file of this name:
file = STRTRIM(name,2)
info = FILE_INFO(file)
IF ~info.EXISTS THEN BEGIN
    file = file+'.fits'
    info = FILE_INFO(file)
    IF ~info.EXISTS THEN BEGIN
        file = STRTRIM(name,2)+'.FITS'
        info = FILE_INFO(file)
        IF ~info.EXISTS THEN MESSAGE, "Can't find file named " + name
    ENDIF
ENDIF

;error_status = 0                ; For debugging
CATCH, error_status ; In case it's not a fits file
IF error_status THEN BEGIN
    CATCH, /CANCEL
    HELP, /LAST_MESSAGE
    IF N_ELEMENTS(unit) GT 0 THEN FREE_LUN, unit
    MESSAGE, 'Cannot read '+file+'. May not be FITS'
ENDIF

OPENR,  unit, file, /GET_LUN, ERROR = err
IF err NE 0 THEN MESSAGE, 'Cannot open file ' + file

MESSAGE, 'Reading '+file, /INFORMATIONAL
howto = 1

FXHREAD, unit, primhdr, status
IF status NE 0 THEN MESSAGE, 'Problem reading primary FITS header'

; At this point exten may not be set

exten_set = N_ELEMENTS(exten) && exten GT 0
IF ~exten_set THEN BEGIN
                                ; Is there more than just a header?
    naxis = SXPAR(primhdr,'NAXIS')
    IF  naxis GT 0 THEN BEGIN   ; Data is in primary HDU
        get_image, file, primhdr, naxis, col, verbose, start, 0, data
        header = primhdr
        IF SIZE(data,/TYPE) EQ 10 THEN howto = 3

        RETURN                  ; all done!
    ENDIF ELSE IF exten_set THEN MESSAGE, 'No data in primary HDU'
ENDIF

; If we reach here the primary HDU was empty or unwanted. Read the
; extension

IF ~exten_set THEN exten = 1

IF exten GT 1 THEN BEGIN
    status = FXMOVE(unit, exten-1)
    IF status NE 0 THEN MESSAGE, $
      STRING(exten, FORMAT = "('Can''t find extension number',I3)")
ENDIF
FXHREAD, unit, xtnhdr, status
FREE_LUN, unit
IF status NE 0 THEN MESSAGE, 'Problem reading extension header'

xtype = SXPAR(xtnhdr,'XTENSION')

IF xtype EQ 'BINTABLE' OR xtype EQ 'A3DTABLE' THEN BEGIN
                                ; binary table
    xtype = 'BINTABLE'          ; replace synonym
    FXBOPEN, unit, file, exten, xtnhdr
                                ; Check that requested column exists
    msg = "Requested column(s) not found"
    nfields = SXPAR(xtnhdr, 'TFIELDS')

    ttype = SXPAR(xtnhdr,'TTYPE*', COUNT = count)
                                ; Check for first column = pixels
                                ; (doesn't count)
    IF count GE 1 THEN cut4 = STRCMP(ttype[0], 'PIXEL', 5, /FOLD_CASE)

    IF ~col_set THEN col = INDGEN(nfields - cut4) + 1 + cut4 ELSE $
      IF SIZE(col,/TYPE) EQ 7 THEN BEGIN ; We have a list of names...
        ttype = STRTRIM(STRUPCASE(ttype), 2)
        col   = STRTRIM(STRUPCASE(col), 2)
        match, ttype, col, colnum, COUNT= count
        IF count LT ncol THEN MESSAGE, msg
        col = colnum + 1  ; From now on use numbers not names
    ENDIF ELSE BEGIN
        IF nfields - cut4 LT MAX(col) THEN MESSAGE, msg
        IF cut4 then col += 1
    ENDELSE
    IF cut4 AND col[0] NE 1 THEN col = [1, col]
    ncol = N_ELEMENTS(col)

    IF ncol EQ 1 THEN FXBREADM, unit, col, data  ELSE BEGIN
        data = PTRARR(ncol)
        FXBREADM, unit, col, POINTERS=data, PASS_METHOD='POINTER'
        SXADDPAR, header, 'POINTER', 'T'
                                ; Convert to 1-D arrays since the only
                                ; option here (so far) is HEALPix arrays.
        npix = N_ELEMENTS(*data[0])
        FOR i=0,ncol-1 DO *data[i] = REFORM(*data[i], npix, /OVERWRITE)

        howto = 3
    ENDELSE
    IF verbose THEN MESSAGE, /INFORMATIONAL, $
      STRING('Binary extension read at', SYSTIME(1) - start)
    IF verbose THEN HELP, /MEMORY
    FXBCLOSE, unit

    dsize = SIZE(data)
    IF dsize[0] GT 1 && ncol EQ 1 THEN $
      data = REFORM(data,dsize[dsize[0]+2],/OVERWRITE)
                                ; We only read one column
    nmand = 8
ENDIF ELSE IF xtype EQ 'IMAGE   ' THEN BEGIN ; Image extension
    naxis = SXPAR(xtnhdr,'NAXIS')
    get_image, file, xtnhdr, naxis, col, verbose, start, 1, data
    nmand = 5+naxis
    IF SIZE(data,/TYPE) EQ 10 THEN howto = 3
ENDIF ELSE MESSAGE, 'Cannot process extension type '+xtype

; Remove SIMPLE, BITPIX, NAXIS and END cards from primary header and
; add rest to extension header immediately after the mandatory
; keywords.
lp = N_ELEMENTS(primhdr) - 2
primhdr = primhdr[3:lp]
; Strip keywords from primary that recur in extension:
keys = STRMID(primhdr, 0, 8)
pactive = WHERE (keys NE 'COMMENT ' AND keys NE 'HISTORY ' AND $
                 keys NE 'HIERARCH', nprim)

IF nprim EQ 0 THEN header = TEMPORARY(xtnhdr) ELSE BEGIN
    pkeys = keys[pactive]
    keys = STRMID(xtnhdr[nmand:*], 0, 8)
    eactive = WHERE (keys NE 'COMMENT ' AND keys NE 'HISTORY ' AND $
                     keys NE 'HIERARCH', next)
    ekeys = keys[eactive]
    idx = -1L                   ; dummy starter value
    FOR i=0,nprim-1 DO BEGIN
        null = WHERE(pkeys[i] EQ ekeys, hit)
        IF hit EQ 0 THEN idx = [idx,pactive[i]]
    ENDFOR
    IF N_ELEMENTS(idx) GT 1 THEN $
      header = [xtnhdr[0:nmand-1], primhdr[idx[1:*]], xtnhdr[nmand:*]] $
    ELSE header = TEMPORARY(xtnhdr)
ENDELSE
IF verbose THEN MESSAGE, /INFORMATIONAL, $
  STRING('Headers merged at', SYSTIME(1) - start)

END

PRO parse_input, thing, ORDER = order, PROJ = proj, TEMPORARY = temporary, $
                 data, header, name, howto, colab, scale_par, namestr, $
                 polcode, GET_SCALE = get_scale, COLUMN = col, $
                 EXTENSION = exten, PLOT = do_plot, VERBOSE = verbose
;+
; NAME:
;       PARSE_INPUT
;
; PURPOSE:
;       Decodes input parameter (thing): is it
;
;       (a) A stack of 2-D images
;       (b) A stack of HEALPix (HP) arrays
;       (c) A file name for a FITS file containing images or HP arrays
;       (d) A file name minus the ".fits" or ".FITS" ending
;       (e) A structure containing header + set of 2-D images (one per tag)
;       (f) A structure containing header + a set of HP arrays (one per tag)
;       (g) A pointer to any of the above
;       (h) An array of pointers to 2-D images or HP arrays.
;
;       Headers are assumed to be more or less FITS.
;       Returns extracted data, header, and information from the
;       header. Optionally gets robust statistics of the data for use
;       in scaling. Optionally, plots histograms to check statistics.
;
; CATEGORY:
;       Input/Output, Widget
;
; CALLING SEQUENCE:
;
;       PARSE_INPUT, Thing, Data, Header, Name, Howto, Colab, $
;                    Scale_par, Namestr, Polcode
;
; INPUTS:
;       Thing:     An IDL variable that somehow specifies the wanted data
;
; KEYWORD PARAMETERS:
;       ORDER:     'RING' or 'NESTED' (default: header value or 'RING')
;
;       PROJ:      'GRID', 'NPOLE', 'SPOLE' if this is a healpix array
;                  (default: 'GRID'), or if it is already gridded but
;                  header is missing.
;
;       COLUMN:    Requested table column(s) or image plane(s). (Below
;                  referred to simply as maps). Default: all
;
;       EXTENSION: FITS extension containing the data
;                  Default: look in 0, then 1.
;
;       GET_SCALE: If set, get statistics of the data
;
;       PLOT:      Set to launch a widget containing plots of the
;                  histogram of each map around the noise level (only
;                  if /GET_SCALE is also set).
;
;       TEMPORARY: Input array should be overwritten to save space
;
;       VERBOSE:   Print diagnostics
;
; OUTPUTS:
;       Data:      2+D image or array of pointers to 2+D image
;                  (ie 2D or 3D if more than one map requested).
;
;       Header:    FITS-like header describing data (may not fully
;                  conform to FITS standard)
;
;       Name:      Text string to use as suggested title
;
;       Howto:     = 0: use the original variable (= 2+D image)
;                  = 1: data is the 2+D image
;                  = 2: data is a pointer to a 2D image
;                  = 3: data is an array of pointers to 2D images.
;
; OPTIONAL OUTPUTS:
;      Colab:      Array of strings usable as labels for individual maps
;
;      Scale_par:  [4, N_map] array containing min, max, mode, standard
;                  deviation in each map. (Set if /GET_SCALE is set).
;
;      Namestr:    Structure containing header elements suitable for
;                  constructing the "Name" and "Colab" strings.
;
;      Polcode:    Array of code numbers describing polarization each
;                  map: standard FITS "STOKES" codes plus:
;                     5 = linearly polarized intensity
;                     6 = Fractional linear polarization
;                     7 = Polarization angle
;
; SIDE EFFECTS:
;      IF /GET_SCALE and /PLOT are set, a widget is launched showing
;      the histograms of each map around the mode, overplotted with
;      the parabolic fit used to find the mode.
;
;      Thing may be overwritten if /TEMPORARY is set.
;
;      If Thing is a HEALPix array or an array of pointers to HEALPix
;      arrays, any HEALPix bad values will be changed to NaNs. (even
;      if /TEMPORARY is not set).
;
; EXAMPLE:
;               dummy = FINDGEN(100,100)
;               PARSE_INPUT, dummy, data, header, name, howto
;
;      Returns howto = 0, data is undefined, header is a minimal FITS
;      header, and name = "Online data:"
;
;
; MODIFICATION HISTORY:
;       Written by:     J. P. Leahy, January 2008
;       March 2008:     Added CUT4 support
;       April 2008:     Bug fixes
;-
COMPILE_OPT IDL2
ON_ERROR, 2

start = SYSTIME(1)
ispointer = 0
path = ''  &  file = ''
xtype = 'IMAGE'
do_plot = KEYWORD_SET(do_plot)  ; Flag to plot histogram of data.

temporary = KEYWORD_SET(temporary)
; Set temporary: automatically true if there is nothing to return it to...
intype = SIZE(thing, /TYPE)      ; 10 = pointer
temporary = temporary || (intype NE 10 && ~ARG_PRESENT(thing))

verbose = KEYWORD_SET(verbose)
tsize = SIZE(thing,/STRUCTURE)
ncol = N_ELEMENTS(col)
col_set = ncol GT 0
IF ncol GT 1 THEN BEGIN
    col = col[UNIQ(col, BSORT(col))] ; sort & eliminate duplicates
    ncol = N_ELEMENTS(col)
ENDIF

REDO:

CASE tsize.TYPE OF
    0: MESSAGE, 'Undefined image parameter'

    7: BEGIN ; String... should be a file name
        name = ispointer ? (*thing)[0] : thing[0]

        get_fits, name, col, exten, file, data, header, howto, xtype, $
          verbose, start
                        ; Strip directory info from filename (if any)
        sep = PATH_SEP()
        IF sep EQ '\' THEN sep = '\\'
        firstchar = STRSPLIT(file, '.*'+sep, /REGEX)
        path = STRMID(file, 0, firstchar-1)
        file = STRMID(file, firstchar)
        temporary = 1
        ncol = N_ELEMENTS(col)
        col_set = ncol GT 0
        IF ~col_set THEN MESSAGE, 'Internal error: no columns read from file'
    END

    8: BEGIN                    ; Structure or pointer to structure
        IF tsize.N_DIMENSIONS GT 1 THEN $
          MESSAGE, 'Expected image but found array of structures'

                                ; Check it has header and (potential) data:
        header = ispointer ? *thing.(0) : thing.(0)
        hdtype = SIZE(header, /TYPE)
        badhd = hdtype NE 7     ; Should be a FITS header, ie strings.
        IF ~badhd THEN BEGIN    ; check it is a FITS-ish header
            simple = SXPAR(header,'SIMPLE')
            IF ~simple THEN BEGIN ; Not primary header
                xtype = SXPAR(header,'XTENSION',COUNT=count)
                badhd = count EQ 0
            ENDIF
        ENDIF
        ntags = N_TAGS(ispointer ? *thing : thing)
        IF badhd || ntags LT 2 THEN $
          MESSAGE, 'Unknown structure type supplied as input'
        nfields = ntags - 1     ; Number of tags with data
        tagnames = TAG_NAMES(ispointer ? *thing : thing)
                                ; Check for cut4 file
        IF nfields GE 2 THEN $
          cut4 = STRCMP(tagnames[1],  'PIXEL', 5, /FOLD_CASE) && $
                 STRCMP(tagnames[2], 'SIGNAL', 6, /FOLD_CASE) ELSE cut4 = 0

                                ; Have a look at first item to see
                                ; what type it is:
        data = ispointer ? *thing.(1) : thing.(1)
        dsize = SIZE(data, /STRUCTURE)
        datatype = dsize.TYPE
        IF datatype EQ 10 THEN BEGIN
            SXADDPAR, header, 'POINTER', 'T'
            howto = 2
            dsize = SIZE(*data, /STRUCTURE)
        ENDIF ELSE howto = 1
        wasone =  howto EQ 1

        is_image = 0B
        IF ntags EQ 2 && dsize.N_DIMENSIONS GT 1 THEN BEGIN
            is_image = 1B
                                ; See if this is a healpix dataset:
            hp_check, order, header, dsize, healpix, cut4, nside
            npix = dsize.N_ELEMENTS
            nax  = dsize.DIMENSIONS
            mappix = healpix ? nax[0] : nax[0] * nax[1]
            nfields =  npix / mappix
        ENDIF

                                ; Potential message if things go wrong:
        msg = 'Input structure does not contain all requested columns'
        IF ~col_set THEN BEGIN
            col = INDGEN(nfields) + 1
            col_set = 1B
        ENDIF ELSE IF ~is_image && SIZE(col, /TYPE) EQ 7 THEN BEGIN
                                ; list of names...
            IF count LT ncol - cut4 THEN MESSAGE, msg

            tagnames = STRTRIM(STRUPCASE(tagnames), 2)
            col      = STRTRIM(STRUPCASE(col), 2)
            match, tagnames[1:*], col, colnum, COUNT = count
            IF count LT ncol THEN MESSAGE, msg
            col = colnum + 1    ; From now on use numbers not names
        ENDIF ELSE BEGIN
            IF nfields - cut4 LT MAX(col) THEN MESSAGE, msg
            IF cut4 then col += 1
        ENDELSE
        IF cut4 AND col[0] NE 1 THEN col = [1, col]
        ncol = N_ELEMENTS(col)

        IF nfields LT MAX(col) THEN MESSAGE, msg

        IF ncol GT 1 THEN BEGIN ; Store data via a pointer array
            tmp = TEMPORARY(data)
            IF is_image THEN BEGIN ; one 2+D image
                                ; Delete original now, as we need yet
                                ; another copy...
                IF temporary THEN BEGIN
                    IF ispointer THEN PTR_FREE, thing
                    ispointer = 0B
                ENDIF
                IF howto EQ 2 THEN BEGIN ; swap to non-heap variable
                    tmp2 = temporary ? TEMPORARY(*tmp) : *tmp
                    IF temporary THEN PTR_FREE, tmp
                    tmp = TEMPORARY(tmp2)
                ENDIF

                tmp = REFORM(tmp, mappix, nfields, /OVERWRITE)
                data = PTRARR(ncol, /ALLOCATE_HEAP)
                FOR i=0,ncol-1 DO BEGIN
                    tmp2 = tmp[*,col[i]-1]
                    IF ~healpix THEN tmp2 = REFORM(tmp2, nax[0], nax[1], $
                                                   /OVERWRITE)
                    *data[i] = TEMPORARY(tmp2)
                ENDFOR
                tmp = 0
            ENDIF ELSE BEGIN ; One col per tab
                IF howto EQ 2 THEN BEGIN
                    data = PTRARR(ncol)
                    data[0] = tmp
                    FOR i=1,ncol-1 DO $
                      data[i]  = ispointer ? *thing.(col[i]) : thing.(col[i])
                ENDIF ELSE BEGIN
                    data = PTRARR(ncol, /ALLOCATE_HEAP)
                    *data[0] = TEMPORARY(tmp)
                    FOR i=1,ncol-1 DO $
                      *data[i] = ispointer ? *thing.(col[i]) : thing.(col[i])
                    SXADDPAR, header, 'POINTER', 'T'
                ENDELSE
            ENDELSE
            howto = 3
        ENDIF
        IF temporary THEN BEGIN ; Delete thing
; NB: if structure contained pointers, their heap variables have been
; copied to data, unless is_image in which case any pointer was freed above.
            IF ispointer THEN PTR_FREE, thing
            thing = 0
        ENDIF
        IF wasone THEN temporary = 1B ; we can delete data heap variable(s)
    END

    10: BEGIN ; Pointer
        IF ispointer THEN MESSAGE, 'Cannot parse chains of pointers'
        ndims = tsize.N_DIMENSIONS
        CASE ndims OF
            0: BEGIN ; Pointer to something else: find out what
                tsize = SIZE(*thing,/STRUCTURE)
                ispointer = 1
                GOTO, REDO
            END
            1: BEGIN ; Should be array of pointers to 1- or 2-D datasets
                ncol = tsize.DIMENSIONS[0]
                testcol = INDGEN(ncol) + 1
                IF col_set && ARRAY_EQUAL(col,testcol) EQ 0B THEN MESSAGE, $
                  'Please trim data array to required size'
                col = testcol
                col_set = 1B
                T0 = SIZE(*thing[0])
                SWITCH T0[0] OF
                    0: MESSAGE, $
                      'Pointer arrays should point directly to bulk data!'
                    1: ; Same as 2:
                    2: BEGIN
                        FOR i=1,tsize.N_DIMENSIONS DO BEGIN
                            T = SIZE(*thing[i])
                            IF ARRAY_EQUAL(T0[0:T0[0]], T[0:T[0]]) EQ 0B $
                              THEN MESSAGE, 'Mismatched array sizes in input'
                        ENDFOR
                        tsize = SIZE(*thing[0], /STRUCTURE)
                        MKHDR, header, tsize.TYPE, $
                          [tsize.DIMENSIONS[0:T0[0]-1], ndims]
                        SXADDPAR, header, 'POINTER', 'T'
                        data = thing
                        howto = 3
                        BREAK
                    END
                    ELSE: MESSAGE, 'Each pointer should point to 1 or 2-D data'
                ENDSWITCH
            END
            ELSE: MESSAGE, $
              'Expected image but found multi-dimensional array of pointers'
        ENDCASE
    END
    ELSE: BEGIN                 ; Array (hopefully) of numbers
                ; Make sure we have only been passed the planes needed

        MKHDR, header, tsize.TYPE, tsize.DIMENSIONS[0:tsize.N_DIMENSIONS-1]
                                ; See if this is a healpix dataset:
        hp_check, order, header, tsize, healpix, cut4, nside
        npix  = tsize.N_ELEMENTS
        nax   = tsize.DIMENSIONS
        mappix = healpix ? nax[0] : nax[0] * nax[1]
        IF mappix EQ 0 THEN MESSAGE, $
          'Cannot process input: 1-D and does not seem to be HEALPix'

        nblock =  npix / mappix
        testcol = LINDGEN(nblock) + 1
        IF col_set && ~ARRAY_EQUAL(col, testcol) THEN MESSAGE, $
                'Please trim data array to required size'
        col = testcol           ; Assume that all these data passed are wanted.
        ncol = nblock
        col_set = 1B
        IF ncol GT 1 && ~healpix THEN BEGIN
                                ; swap to non-heap variable
            IF temporary THEN BEGIN
                tmp = TEMPORARY(ispointer ? *thing : thing)
                IF ispointer THEN PTR_FREE, thing
            ENDIF ELSE tmp = ispointer ? *thing : thing
            tmp = REFORM(tmp, mappix, nblock, /OVERWRITE)
            data = PTRARR(ncol, /ALLOCATE_HEAP)
            FOR i=0,ncol-1 DO BEGIN
                tmp2 = tmp[*,col[i]-1]
                IF ~healpix THEN tmp2 = REFORM(tmp2, nax[0], nax[1], $
                                               /OVERWRITE)
                *data[i] = TEMPORARY(tmp2)
            ENDFOR
            tmp = 0
            SXADDPAR, header, 'POINTER', 'T'
            howto = 3
        ENDIF ELSE IF ispointer THEN BEGIN
            SXADDPAR, header, 'POINTER', 'T'
            data = thing
            howto = 2
        ENDIF ELSE IF temporary THEN BEGIN
            data = TEMPORARY(thing)
            howto = 1
        ENDIF ELSE howto = 0
    END
ENDCASE

IF verbose THEN MESSAGE, /INFORMATIONAL, $
  STRING('Got '+xtype+' data at', SYSTIME(1) - start)
IF verbose THEN HELP, /MEMORY

datapt = SXPAR(header,'POINTER')
naxis  = SXPAR(header,'NAXIS')

CASE howto OF
    0: dsize = SIZE(thing, /STRUCTURE)
    1: dsize = SIZE(data, /STRUCTURE)
    2: dsize = SIZE(*data, /STRUCTURE)
    3: dsize = SIZE(*data[0], /STRUCTURE)
ENDCASE
IF dsize.N_DIMENSIONS EQ 0 THEN $
  MESSAGE, 'Found scalar when image expected'

; What kind of data do we have?
type = dsize.TYPE
IF howto EQ 3 THEN BEGIN
    type = INTARR(ncol)
    FOR i=0,ncol-1 DO type[i] = SIZE(*data[i], /TYPE)
ENDIF
utype = UNIQ(type[SORT(type)])
match, utype, [7, 8, 10, 11], suba, subb
badtypes = ['STRING', 'STRUCTURE', 'POINTER', 'OBJREF']
;void = WHERE([7, 8, 10, 11] - dsize.TYPE EQ 0, baddata)
IF subb[0] NE -1 THEN MESSAGE, 'Cannot process data type ' + badtypes[subb]

; Find safe data type for means etc:
match, utype, [4, 5, 6, 9], suba, subb
IF subb NE -1 THEN stype = MAX(utype[suba]) ELSE stype = 5

; See if this is a healpix dataset, if not already done:
IF N_ELEMENTS(healpix) EQ 0 THEN $
  hp_check, order, header, dsize, healpix, cut4, nside

IF KEYWORD_SET(get_scale) THEN BEGIN
                                ; Analyse each channel for scaling purposes:
    scale_par = MAKE_ARRAY(4, ncol, TYPE = stype)
    nblock = ncol - cut4

    IF do_plot THEN BEGIN ; Launch widget for plot
        ttit = 'Data Histogram'
        IF nblock GT 1 THEN ttit += 's'
        tlb = WIDGET_BASE(/COLUMN, TITLE = ttit)
        row  = WIDGET_BASE(tlb, /ROW)
        draw   = LONARR(nblock)
        FOR imap = 0, nblock-1 DO BEGIN
            icol = col[imap + cut4]
            base = WIDGET_BASE(row, /COLUMN)
            void = WIDGET_LABEL(base, VALUE = $
                                'Column '+STRTRIM(STRING(icol),2) )
            draw[imap] = WIDGET_DRAW(base, XSIZE=300, YSIZE=300, RETAIN = 2)
        ENDFOR
        done = WIDGET_BUTTON(tlb, VALUE = 'Done')
        WIDGET_CONTROL, tlb, /REALIZE

        swap_lut, dummy, old_graph
    ENDIF

    FOR imap = 0, nblock-1 DO BEGIN
        IF do_plot THEN BEGIN
            WIDGET_CONTROL, draw[imap], GET_VALUE = window
            WSET, window
        ENDIF
        scaling_params, thing, data, imap+cut4, healpix, howto, col, $
          ar, mode, sdev, PLOT = do_plot
        PRINT, imap, ar, mode, sdev, FORMAT = $
       "('Column',I2,': Min & MAX:',2E11.3,', Mode:',E11.3,', est RMS:',E10.3)"
        scale_par[0,imap] = [ar, mode, sdev]
    ENDFOR

    IF verbose THEN MESSAGE, /INFORMATIONAL, $
      STRING('Scaling parameters derived', SYSTIME(1) - start)
    IF verbose THEN HELP, /MEMORY
    IF do_plot THEN BEGIN
        restore_lut, old_graph
                                ; Allow user to get rid of plot!
        XMANAGER, 'parse_input', tlb, /NO_BLOCK
    ENDIF
ENDIF

IF healpix THEN BEGIN    ; seems to be HEALPix: convert to grid.
                                ; First convert HEALPix bad values to NaN
    IF howto GE 2 THEN FOR icol = 0L, ncol-1 DO BEGIN
        dptr = data[icol]
        bad2nan, *dptr
    ENDFOR ELSE IF howto EQ 1 THEN bad2nan, data ELSE bad2nan, thing

    IF cut4 THEN BEGIN          ; must be howto = 3
        tmp = cut4grid(*data[0], data[1:*], header, order, proj, $
                       VERBOSE = verbose)
        IF temporary THEN PTR_FREE, data
        data = TEMPORARY(tmp)
        col = col[1:*] - 1      ; We have eliminated col 0 = PIXEL
        ncol = ncol - 1
    ENDIF ELSE CASE howto OF
        0: data = hpgrid(thing, header, order, proj, VERBOSE = verbose)
        1: data = hpgrid(TEMPORARY(data), header, order, proj, VERBOSE=verbose)
        2: data = temporary ? hpgrid(TEMPORARY(*data), header, order, proj, $
                                     VERBOSE = verbose) $
                            : hpgrid(*data, header, order, proj, $
                                     VERBOSE = verbose)
        3: BEGIN
            tmp = hpgrid(data, header, order, proj, VERBOSE = verbose)
            IF temporary THEN PTR_FREE, data
            data = TEMPORARY(tmp)
        END
    ENDCASE

    IF dsize.N_DIMENSIONS GT 1 THEN howto = 3 ; we get a pointer array

    IF howto LT 3 THEN BEGIN
        dsize = SIZE(data, /STRUCTURE)
        howto = 1
    ENDIF ELSE dsize = SIZE(*data[0], /STRUCTURE)

    temporary = 1B

    IF verbose THEN MESSAGE, /INFORMATIONAL, $
      STRING('Converted to HEALPix at', SYSTIME(1) - start)

ENDIF ELSE BEGIN                ; definitely not HEALPix array
    IF dsize.N_DIMENSIONS EQ 1 THEN $
      MESSAGE, 'Found 1-D array when image expected'
ENDELSE

IF verbose THEN HELP, /MEMORY

ndims = dsize.N_DIMENSIONS
extras = ndims GT 2

; If we are allowed to, strip data down to 2D (but this should never happen)
IF ~col_set THEN BEGIN
    MESSAGE, /INFORMATIONAL, 'Internal problem: col unset'
    CASE howto OF
        0:                      ; Do nothing: never temporary
        1: IF extras THEN data = data[*,*]
        2: IF temporary && extras THEN *data = *data[*,*]
    ENDCASE
    IF temporary OR howto EQ 1 THEN extras = 0
ENDIF

; Interpret & maybe modify header:

types = SXPAR(header,'CTYPE*',COUNT = count)
IF count LT naxis THEN BEGIN
    extra = STRING(INDGEN(naxis-count)+count, FORMAT = "('AXIS',I2)")
    types = count GT 0 ? [types,extra] : extra
ENDIF ELSE IF count GT naxis THEN naxis = count
types = STRTRIM(types,2)

IF ~healpix && KEYWORD_SET(proj)  THEN BEGIN
                                ; May be already gridded HEALPix
    ngrid = dsize.DIMENSIONS[0]
    nside = STRCMP(proj,'GRID', 4, /FOLD_CASE) ? ngrid/5 : ngrid/4
    npix = NSIDE2NPIX(nside)
    IF npix EQ -1 OR ngrid NE dsize.DIMENSIONS[1] THEN BEGIN
        MESSAGE, /INFORMATIONAL, 'Projection set but image wrong size'
    ENDIF ELSE BEGIN            ; Update header
        header = extras ? grid_header(header, nside, proj, ndims-2, $
                                      dsize.DIMENSIONS[2:ndims-1]) $
                        : grid_header(header, nside, proj, 0, 0)
    ENDELSE
ENDIF

IF verbose THEN MESSAGE, /INFORMATIONAL, $
  STRING('Data ready', SYSTIME(1) - start)
IF verbose THEN HELP,/MEMORY

; Form title string:
scodes = ['YX','XY','YY','XX', 'LR', 'RL', 'LL','RR','','I','Q','U','V', $
          'PPOL', 'FPOL', 'PANG']
polcode = REPLICATE(0, ncol)

telescope = SXPAR(header, 'TELESCOP', COUNT = count, /SILENT)
IF count EQ 0 THEN telescope = ''
instrument = SXPAR(header, 'INSTRUME', COUNT = count, /SILENT)
IF count EQ 0 THEN instrument = ''
IF instrument EQ telescope THEN instrument = ''

object = SXPAR(header, 'OBJECT', COUNT = count, /SILENT)
IF count EQ 0 THEN object = ''
creator = SXPAR(header, 'CREATOR', COUNT = count, /SILENT)
IF count EQ 0 THEN creator = ''

IF xtype EQ 'BINTABLE' THEN BEGIN
    ttype = (SXPAR(header, 'TTYPE*', COUNT = count))[col-1]
                                ; Remove underscores & shorten
    FOR i=0,ncol-1 DO BEGIN
        ttype[i] = STRJOIN(STRSPLIT(ttype[i],'_',/EXTRACT),' ')
        ipol  = STREGEX(ttype[i], 'polari(s|z)ation|stokes', /FOLD_CASE, $
                        LENGTH = len)
        IF ipol GE 0 THEN BEGIN
            first = STRMID(ttype[i],0,ipol)
            last  = STRMID(ttype[i],ipol+len)
            icode = WHERE(scodes EQ STRTRIM(first,2))
            IF icode EQ -1 THEN  icode = WHERE(scodes EQ STRTRIM(last,2))
            IF icode NE -1 THEN polcode[i] = icode - 8
;            IF ipol EQ 0 THEN ttype[i] = first + ' Pol'
;            IF ipol GE 0 THEN ttype[i] = 'Pol ' + last
            ttype[i] = first + 'Pol' + last
        ENDIF ELSE BEGIN
            itemp = STREGEX(ttype[i], 'temperature', /FOLD_CASE, LENGTH = len)
            IF itemp GE 0 THEN ttype[i] = STRMID(ttype[i],0,itemp) + 'Temp' $
              + STRMID(ttype[i],itemp+len)
        ENDELSE
    ENDFOR
    colab = STRTRIM(ttype,2)
ENDIF ELSE BEGIN
    colab = 'Map' + STRING(INDGEN(ncol)+1, FORMAT = "(I2)")
ENDELSE

fval = 0
stokes = ''

IF xtype EQ 'IMAGE' && naxis GT 2 THEN BEGIN  ; Process extra axes in header
    nax = SXPAR(header,'NAXIS*', COUNT = count)
    IF count LT naxis THEN nax = [nax, REPLICATE(1,naxis-count)]
                                ; Find planes requested
    colcoord = INTARR(ncol,naxis-2)
    active = BYTARR(naxis)
    block = REPLICATE(1,naxis-2)
    FOR i = 2,naxis-2 DO block[i-1] *= nax[i]*block[i-2]
    list = col-1
    FOR i=naxis-1,2,-1 DO BEGIN
        colcoord[0,i-2] = list / block[i-2]
        list = list - colcoord[*,i-2]*block[i-2]
        active[i] = MIN(colcoord[*,i-2]) LT MAX(colcoord[*,i-2])
    ENDFOR

    rvals = SXPAR(header,'CRVAL*',COUNT=count)
    IF count LT naxis THEN rvals = count GT 0 $
      ? [rvals,REPLICATE(0.,naxis-count)] : REPLICATE(0.,naxis)
    rpixs = SXPAR(header,'CRPIX*',COUNT=count)
    IF count LT naxis THEN rpixs = count GT 0 $
      ? [rpixs,REPLICATE(0.,naxis-count)] : REPLICATE(0,naxis)
    delts = SXPAR(header,'CDELT*',COUNT=count)
    IF count LT naxis THEN delts = count GT 0 $
      ? [delts,REPLICATE(0.,naxis-count)] : REPLICATE(1.,naxis)
    units = SXPAR(header,'CUNIT*',COUNT=count)
    IF count LT naxis THEN units = count GT 0 $
      ? [units,REPLICATE(0.,naxis-count)] : REPLICATE('',naxis)

    IF ncol GT 1 THEN BEGIN     ; At least one dimension > 2 is active
        active = WHERE(active EQ 1B, nactive)
        colco = TRANSPOSE(colcoord[*, active])
        colab = REPLICATE('',ncol)
        da = REBIN(delts[active],nactive,ncol)
        ra = REBIN(rpixs[active],nactive,ncol)
        va = REBIN(rvals[active],nactive,ncol)
        vals = da*(colco - (ra - 1)) + va
        FOR i = 0,ncol-1 DO BEGIN
            FOR j = 0,nactive-1 DO BEGIN
                IF STRCMP(types[active[j]],'STOKES',6) THEN BEGIN
                    istokes = ROUND(vals[j,i])
                    polcode[i] = istokes
                    colab[i] = colab[i] + scodes[istokes+8] + ' Pol'
                ENDIF ELSE colab[i] = colab[i] + ' ' + types[active[j]] + $
                  ': ' + numunit(vals[j,i], units[j])
            ENDFOR
            colab[i] = STRTRIM( STRCOMPRESS(colab[i]), 2)
        ENDFOR
    ENDIF
                                ; Search for frequency info
    fcodes = ['FREQ', 'WAVE', 'AWAV',  'WAVN', 'ENER']
    defun  = ['Hz', 'm', 'm', 'm-1', 'J']
    FOR i=2,naxis-1 DO BEGIN
        match = WHERE(STRCMP(types[i],fcodes,4),hit)
        IF hit GT 0 && ~active[i] THEN BEGIN
            ft = fcodes[match]
            fval  = delts[i]*(colcoord[0,i-2]-rpixs[i]+1) + rvals[i]
            funit = units[i] NE '' ? units[i] : defun[match]
            BREAK
        ENDIF
    ENDFOR
                                 ; Search for stokes info
    match = WHERE(STRCMP(types,'STOKES',6),hit)
    IF hit GT 0 && ~active[match] THEN BEGIN
        match = match[hit-1]
        istokes = ROUND( delts[match]*(colcoord[0,match-2] - $
                                       rpixs[match]+1) + rvals[match] )
        stokes = scodes[istokes+8]+' Pol'
        polcode[*] = istokes
    ENDIF
ENDIF

IF fval GT 0. THEN freq = numunit(fval,funit) $
ELSE BEGIN                      ; Try non-standard keywords:
    freq = SXPAR(header,'FREQ', COUNT = hit)
    IF hit THEN BEGIN
        ftype = SIZE(freq,/TYPE)
        IF ftype NE 7 THEN BEGIN ; not a string: try to guess units
            IF freq LT 1E6 THEN freq = numunit(freq,' GHz')+' (?)' $
            ELSE freq = numunit(freq,' Hz')+' (?)'
        ENDIF
    ENDIF ELSE freq = ''
ENDELSE

namestr = {path: path, file: file, telescope: telescope, $
           instrument: instrument, $
           creator: creator, object: object, freq: freq, stokes: stokes}

namestr.FILE       = file       ? file : 'Online data'
namestr.TELESCOPE  = telescope  ?  ' ' + STRTRIM(telescope, 2)  : ''
namestr.INSTRUMENT = instrument ?  ' ' + STRTRIM(instrument, 2) : ''
namestr.CREATOR    = creator    ?  ' Created by ' + STRTRIM(creator, 2) : ''
namestr.OBJECT     = object     ?  ' ' + STRTRIM(object, 2)  : ''
namestr.FREQ       = freq       ? '  ' + STRTRIM(freq, 2)   : ''
namestr.STOKES     = stokes     ? '  ' + STRTRIM(stokes, 2) : ''

name = namestr.FILE + ': '
FOR i=2,7 DO name += namestr.(i)

END
