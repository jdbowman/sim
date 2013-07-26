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
PRO restore_lut, scale
;
; Restore old graphics state
;
; Inputs:
;    scale: structure containing old details created by swap_lut
;
COMPILE_OPT IDL2, HIDDEN
COMMON gr_global

IF scale.WINDOW EQ -1 THEN WSET, -1 ELSE BEGIN
    DEVICE, WINDOW_STATE = ws
    IF ws[scale.WINDOW] EQ 1 THEN WSET, scale.WINDOW ELSE BEGIN
        WSET, -1
        IF scale.WINDOW NE 0 THEN MESSAGE, /INFORMATIONAL, $
          'Cannot find old window #'+STRTRIM(STRING(scale.WINDOW),2)
    ENDELSE
ENDELSE

DEVICE, DECOMPOSED = scale.DECOMPOSED
IF windev NE scale.DEVICE THEN BEGIN
    SET_PLOT, scale.DEVICE
    !X = scale.X  &  !Y = scale.Y  &  !Z = scale.Z
ENDIF
TVLCT, scale.RED, scale.GREEN, scale.BLUE

!P.background = scale.background
!P.color      = scale.COLOUR

END
