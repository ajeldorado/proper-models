; Copyright 2020, by the California Institute of Technology. ALL RIGHTS
; RESERVED. United States Government Sponsorship acknowledged. Any
; commercial use must be negotiated with the Office of Technology Transfer
; at the California Institute of Technology.
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;+
; NAME:
;	TAG_EXISTS
; PURPOSE:		
;	Determine if a specified structure tag (member name) exists in
;	a provided structure
; CALLING SEQUENCE:	
;	result = tag_exists( tagname, structure )
; INPUTS:
;	tagname  :	String containing the name of the tag to test for
;	structure :	The variable that contains the struture to test
; RETURNS:
;	Zero (0) if the tag does not exist in the structure, one (1) if it does
; PROCEDURE:
;	Example:
;	a = {x:0, y:0}
;	print, tag_exists('x', a)
;	         1
;	print, tag_exists('z', a)
;	         0
;-
function tag_exists, tag, struct

names = tag_names( struct )
if ( total(names eq strupcase(tag)) gt 0 ) then return, 1 else return, 0

end  
