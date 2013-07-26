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
; Module XIMVIEW_FITSLOAD
;
; J. P. Leahy February 2008
;
; Contains the event handlers for loading fits files on the fly in Ximview.
;
PRO ximview_fitsload_event, event
; Processes events from the "FITS load" dialog
;
COMPILE_OPT IDL2, HIDDEN
ON_ERROR, 1

WIDGET_CONTROL, event.TOP, GET_UVALUE = info

tag = TAG_NAMES(event, /STRUCTURE_NAME)
CASE STRMID(tag,0,12) OF
    'FILESEL_EVEN': BEGIN
        CASE event.DONE OF
            0: RETURN
            1: BEGIN
                WIDGET_CONTROL, info.EXTID,  GET_VALUE = exten
                WIDGET_CONTROL, info.COLID,  GET_VALUE = column
                WIDGET_CONTROL, info.AUTOID, GET_VALUE = no_auto
                WIDGET_CONTROL, info.PLOTID, GET_VALUE = no_plot
                WIDGET_CONTROL, info.FTEXT,  GET_VALUE = name
                order  = WIDGET_INFO(info.SORTID, /DROPLIST_SELECT)
                IF SIZE(info.HPX,/type) EQ 10 THEN BEGIN
                    proj   = WIDGET_INFO(info.hpx.PROJID, /DROPLIST_SELECT)
                    WIDGET_CONTROL, info.hpx.ROLLID, GET_VALUE = no_roll
                ENDIF ELSE BEGIN
                    proj = ''  &  no_roll = 0B
                ENDELSE

                result = {cancel: 0B, exten: exten, file: event.VALUE, $
                          column: column, auto: 1 - no_auto, $
                          plot: 1 - no_plot, name: name[0], $
                          proj: proj, order: order, roll: 1- no_roll}
            END
            2: result = {cancel: 1B}
        ENDCASE
        WIDGET_CONTROL, info.RETID, SET_UVALUE = result
        WIDGET_CONTROL, event.TOP, /DESTROY
    END
    'WIDGET_TEXT_':             ;
    'WIDGET_DROPL':             ; no action needed
    '': IF event.ID EQ info.AUTOID THEN $ ; Enable plot if autoscale
      WIDGET_CONTROL, info.PLOTID, SENSITIVE = 1S - event.VALUE
    ELSE: MESSAGE, /INFORMATIONAL, 'Received unexpected event ' + tag
ENDCASE

END

PRO ximview_fitsload, event
; Gets file info from user via dialog box and reads in a file
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global
ON_ERROR, 1

start = SYSTIME(1)
WIDGET_CONTROL, event.TOP, GET_UVALUE = state
WIDGET_CONTROL, state.TABS, GET_UVALUE = mode

tabarr = state.TABARR
ntab = N_ELEMENTS(*tabarr)
verbose = state.VERBOSE

first = ntab EQ 0

query = WIDGET_BASE(GROUP_LEADER = event.TOP, /MODAL, $
                    /COLUMN, TITLE = 'Load FITS data')
IF ~first THEN void  = WIDGET_LABEL(query, VALUE = $
                      '!! New image must align with previous ones !!')
base1  = WIDGET_BASE(query, /ROW)
left   = WIDGET_BASE(base1, /COLUMN)
base   = WIDGET_BASE(left,/ROW)
void   = WIDGET_LABEL(base, VALUE = 'FITS extension:', $
                      /ALIGN_LEFT)
extID  = WIDGET_TEXT(base, /EDITABLE, XSIZE = 5)
void   = WIDGET_LABEL(left, /ALIGN_LEFT, VALUE = $
                      '(Leave blank for first data in file)')

base   = WIDGET_BASE(left,/ROW)
void   = WIDGET_LABEL(base, VALUE = 'Columns/slices (list):', $
                      /ALIGN_LEFT)
colID  = WIDGET_TEXT(base, /EDITABLE, XSIZE = 12, VALUE = '*')

autoID = CW_BGROUP(left, ['Yes', 'No'], /ROW, /EXCLUSIVE, $
                   LABEL_LEFT = 'Auto-scale input?', SET_VALUE = 1)
plotID = CW_BGROUP(left, ['Yes', 'No'], /ROW, /EXCLUSIV, $
                   LABEL_LEFT = 'Plot histogram?  ', SET_VALUE = 1)
WIDGET_CONTROL, plotID, SENSITIVE = 0

IF first THEN BEGIN
    base   = WIDGET_BASE(left, /ROW)
    void   = WIDGET_LABEL(base,  VALUE = 'Title:')
    ftext  = WIDGET_TEXT(base, /EDITABLE, XSIZE = 40)
    void   = WIDGET_LABEL(left, /ALIGN_LEFT, VALUE = $
             '                 (Leave blank to use FITS header info)')
ENDIF ELSE BEGIN
    base  = WIDGET_BASE(left, /ROW)
    void  = WIDGET_LABEL(base, VALUE = 'Label for tab:')
    ftext = WIDGET_TEXT(base, /EDITABLE, XSIZE = 10)
ENDELSE
geom   = WIDGET_INFO(base,/GEOMETRY)
void   = WIDGET_BASE(left, /ROW, XSIZE = geom.XSIZE, FRAME = 1) ; spacer
void   = WIDGET_LABEL(left, VALUE = 'For HEALPix datasets only:', $
                      /ALIGN_LEFT)
sortID = WIDGET_DROPLIST(left, TITLE = 'Pixel order:', $
                         VALUE = ['Take from header', 'RING', 'NESTED'])

IF first THEN BEGIN
    projID = WIDGET_DROPLIST(left, VALUE = ['GRID', 'NPOLE', 'SPOLE'], $
                             TITLE = 'Projection:')

    rollID = CW_BGROUP(left, ['Yes', 'No'], /ROW, /EXCLUSIVE, $
                       LABEL_LEFT = 'Interpret image as HEALPix grid?', $
                       SET_VALUE = 1)

    hpx = {projID: projID, rollID: roll_ID}
ENDIF ELSE hpx = 0B

right = WIDGET_BASE(base1, /COLUMN, /ALIGN_BOTTOM)

;void = WIDGET_LABEL(right, VALUE = 'Choose file:', /ALIGN_LEFT)
chooser = CW_FILESEL (right, FILTER=['.fits'], PATH = state.PATH)

                                ; Structure with return ID, internal
                                ; IDs, and defaults
info = {retid: event.ID, extID: extID, colID: colID, sortID: sortID, $
        autoID: autoID, plotID: plotID, ftext: ftext, hpx: hpx, tabarr: tabarr}
WIDGET_CONTROL, query, SET_UVALUE = info, /NO_COPY
WIDGET_CONTROL, event.ID,  GET_UVALUE = save, /NO_COPY

; Fire off widget and collect result
WIDGET_CONTROL, query, /REALIZE
XMANAGER, 'ximview_fitsload', query

WIDGET_CONTROL, event.ID, GET_UVALUE = result, /NO_COPY
WIDGET_CONTROL, event.ID, SET_UVALUE = save, /NO_COPY


                                ; Interpret result here
IF  result.CANCEL THEN RETURN

file = result.FILE
; Convert text to numbers, checking for blank strings:
exten  = STREGEX(result.EXTEN,  '^ *$', /BOOLEAN) ? -1 : FIX(result.EXTEN)

IF ~STREGEX(result.COLUMN, '^ *\*? *$', /BOOLEAN) THEN $
  column = FIX(STRSPLIT(result.COLUMN, ', ',/EXTRACT))

auto_scale = result.AUTO
do_plot = result.PLOT
name = STRTRIM(result.NAME,2)

str = (*tabarr)[0]

IF first THEN BEGIN
    proj = (['GRID', 'NPOLE', 'SPOLE'])[result.PROJ]
    IF result.ORDER NE 0 THEN order = (['','RING', 'NESTED'])[result.ORDER]
    roll = result.ROLL
ENDIF ELSE BEGIN
    proj = state.PROJ
    roll = state.ROLL
ENDELSE
                                ; Turn off blinking, make tab 0 current:
prep_screen, state, mode, tab_arr, oldgraph
                                ; set graphics focus
IF ~state.FOCUS THEN BEGIN
    state.FOCUS = 1
    WIDGET_CONTROL, event.TOP, SET_UVALUE = state
    WIDGET_CONTROL, str.DRAW, /INPUT_FOCUS
    *state.OLD_GRAPH = oldgraph
ENDIF
WIDGET_CONTROL, /HOURGLASS

CATCH, error
IF error NE 0 THEN BEGIN
    CATCH, /CANCEL
    ndat = N_ELEMENTS(data)
    IF ndat NE 0 && SIZE(data,/TYPE) EQ 10 THEN  $
      FOR i=0,ndat-1 DO IF PTR_VALID(data[i]) THEN PTR_FREE, data[i]

    RETURN
ENDIF
                                ; Read in file
parse_input, file, ORDER = order, PROJ = proj, COLUMN = column, $
  data, header, newname, howto, tablab, scale_pars, namestr, polcodes, $
  GET_SCALE = auto_scale, PLOT = do_plot, EXTENSION = exten, VERBOSE = verbose

CATCH, /CANCEL

; create widget tab(s), and stash data in widget uservalues/heap
make_tabs, dummy, proj, column, roll, name, 1B, state, mode, $
  str, first, start, ntab, data, header, newname, howto, $
  tablab, namestr, polcodes, ncol, line, title, extradims, mismatch

IF mismatch THEN RETURN
IF auto_scale THEN BEGIN
                                ; Scale data to byte arrays:
    scale_tabs, ncol, ntab, column, state, auto_scale, scale_pars, howto, $
      extradims, input, range, start
                                ; Draw initial screens
    fill_screens, ncol, ntab, tabarr, mode, first, state, start
ENDIF ELSE BEGIN
                                ; Request manual scaling (ximview_scale)
    newevent = {ID: 0L, TOP: 0L, HANDLER: 0L, VALUE: LINDGEN(ncol)+ntab}
    WIDGET_CONTROL, state.PAD1, SEND_EVENT = newevent
ENDELSE

; Enable blinking if possible
IF ntab + ncol GE 2 THEN BEGIN
    WIDGET_CONTROL, state.BLINK,  /SENSITIVE
    WIDGET_CONTROL, state.FRAMES, /SENSITIVE
    *state.BLINK_SEQ = INDGEN(ntab+ncol)
ENDIF
                                ; Update readout label
title_string = title.HEAD + form_unit((*tabarr)[ntab].UNIT) + title.TAIL
WIDGET_CONTROL, state.READLAB, SET_VALUE = title_string

IF first THEN BEGIN
    PRINTF, state.LOGLUN, line, FORMAT="(A)"
    PRINT, line, FORMAT="(A)"
ENDIF

WIDGET_CONTROL, state.TABS, /SENSITIVE
END
