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
FUNCTION ximget, tab, NO_COPY = no_copy
;+
; NAME: 
;       XIMGET
;
; PURPOSE: 
;       Retrieves an image displayed on a Ximview tab, making it
;       available on the IDL command line.
;
; CATEGORY:
;       Widget helpers
;
; CALLING SEQUENCE: 
;
;        Result = XIMGET(Tab)
;
; INPUTS:
;       tab:     The tab containing the required image. Specify as
;                either the index number (starting with 0 on the left)
;                or the tab title. 
;
; KEYWORD PARAMETERS:
;       NO_COPY: Set to save memory by transfering the image rather
;                than copying it. 
;
; OUTPUTS:
;       Result is a 2-D array, of whatever size and type is stored in
;       XIMVIEW. 
;
; SIDE EFFECTS:
;       If NO_COPY is set, the specified tab on the XIMVIEW widget is
;       deleted.  
;
; RESTRICTIONS:
;       Specifying /NO_COPY when there is only one tab is ignored,
;       since the last tab cannot be deleted. 
;
; EXAMPLE:
;       Display some WMAP map data with XIMVIEW and copy the total
;       intensity image to the command line environment:
;
;            XIMVIEW, 'wmap_band_iqumap_r9_3yr_K_v2', '*', COL=[1,2,3]
;            stokes_I = XIMGET(0)
;
; MODIFICATION HISTORY:
;       Written by:      J. P. Leahy, Feb 2008
;-
COMPILE_OPT IDL2
ON_ERROR, 2

; Find Ximview

test = WIDGET_INFO(/MANAGED)
IF test[0] EQ 0 THEN MESSAGE, 'No widgets currently being managed'

FOR i = 0,N_ELEMENTS(test)-1 DO BEGIN
    uname = WIDGET_INFO(test[i], /UNAME)
    IF uname EQ 'XIMVIEW' THEN GOTO, GOTIT
ENDFOR

MESSAGE, 'Ximview is not currently being managed'

GOTIT:

top = test[i]
WIDGET_CONTROL, top, GET_UVALUE = state

; Find the appropriate tab
IF SIZE(tab,/TYPE) EQ 7 THEN BEGIN ; we have a tab name not number
    labels = STRTRIM(get_tab_uvals(state.TABS), 2)
    itab = WHERE(STRCMP(labels, STRTRIM(tab, 2), /FOLD_CASE), ntab)
    IF itab[0] EQ -1 THEN MESSAGE, 'Tab ' + tab + ' not found'
    IF ntab GT 1 THEN $
        MESSAGE, /INFORMATIONAL, 'Ambiguous tab name, returning first'
    itab = itab[0]
ENDIF ELSE itab = tab

screens = (*state.TABARR).SCREEN
iscreen = WHERE(itab[0] EQ screens)
str = (*state.TABARR)[iscreen]

no_copy = KEYWORD_SET(no_copy)
IF no_copy  && N_ELEMENTS(screens) LE 1 THEN BEGIN
    no_copy = 0B
    MESSAGE, /INFORMATIONAL, 'Overriding /NO_COPY to preserve last tab'
ENDIF

IF ~str.IM_PTR THEN $
  WIDGET_CONTROL, state.LABEL, GET_UVALUE = data, NO_COPY = no_copy $
ELSE IF no_copy THEN BEGIN
    ptr = str.IM_PTR
    data = TEMPORARY(*ptr)
ENDIF ELSE data = *str.IM_PTR

; Delete tab if NO_COPY set
IF no_copy THEN BEGIN
    newevent = {ID: 0L, TOP: 0L, HANDLER: 0L, TAB: itab}
    WIDGET_CONTROL, state.PAD2, SEND_EVENT = newevent
ENDIF

RETURN, data

END

