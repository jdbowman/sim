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
PRO swap_lut, tabstr, old_graph
;
; Saves current graphics state and swaps it for the one set up in the
; tab control structure tabstr.
;
; Inputs:
;   tabstr:      Structure describing current tab

; Inputs in common gr_global
;   windev:      Device used for widget graphics ('X' or 'WIN')
;   colmap:      TRUE if device supports color maps (it should!)
;   redraw_req:  TRUE if a redraw is required to change colour scheme.
;                (will be false for visual classes "PseudoColor" or
;                "DirectColor"
; Output:
;   old_graph:   Structure containing details of former graphics
;   state.
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global, windev, redraw_req, colmap, badcol, syscol

; Set common if needed
IF N_ELEMENTS(windev) EQ 0 THEN BEGIN
    system = STRUPCASE(!VERSION.os_family)
    windev = STRCMP(system, 'WINDOW', 6) ? 'WIN' : 'X'
ENDIF

old_device = !D.name

; Remember previous colour setting
TVLCT, old_red, old_green, old_blue, /GET
old_colour = !P.color
old_background = !P.background

IF old_device NE windev THEN BEGIN
                                ; Switch to windows
    old_x = !X  &  old_y = !Y  &  old_z = !Z
    SET_PLOT, windev
ENDIF ELSE BEGIN                ; Save some space:
    old_x = 0B  &  old_y = 0B  &  old_z = 0B
ENDELSE

; These parameters only apply for windows so get them *after* SET_PLOT
old_window = !D.window

DEVICE, GET_DECOMPOSED = old_decomposed

old_graph = {red: old_red, green: old_green, blue: old_blue, $
             colour: old_colour, background: old_background, $
             x: old_x, y: old_y, z: old_z, window: old_window, $
             decomposed: old_decomposed, device: old_device}

; Install requested graphics state if there is one:
IF SIZE(tabstr, /TYPE) EQ 8 THEN BEGIN
    IF ~colmap THEN BEGIN
        MESSAGE, /INFORMATIONAL, $
          'Apparently this device does not support colour tables'
        RETURN
    ENDIF
    IF ~PTR_VALID(tabstr.LUT) THEN MESSAGE, $
      'Internal error: missing LUT structure on screen ' + $
      STRTRIM( STRING(tabstr.SCREEN), 2)

    lut = *tabstr.LUT
    IF redraw_req THEN DEVICE, DECOMPOSED = 0B $
                  ELSE DEVICE, DECOMPOSED = tabstr.DECOMPOSED

    TVLCT, lut.R, lut.G, lut.B

    !P.color      = lut.LINE
    !P.background = lut.ABSENT
ENDIF

END
