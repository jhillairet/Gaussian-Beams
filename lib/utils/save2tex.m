## Copyright (C) 1996, 1997  Andreas Weingessel 2005 Michael Creel <michael.creel@uab.es>
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details. 
## 
## You should have received a copy of the GNU General Public License
## along with this file.  If not, write to the Free Software Foundation,
## 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

## usage:  save2tex (X, file, cname, rname, prec, align, lines)
##
## Saves the data of the matrix X in a latex tabular environment to
## file. cname and rname are matrices whose rows are the column
## respectively row headings, if c/rname == "" (default), no headings
## are printed. prec gives the number of digits printed after the
## comma. If prec<0 (default) the number is printed unformatted. align
## gives the alignment of the data, default is "r" (flush right). If
## lines == 0 or == "none", no lines are printed in the
## tabular; if lines == 2 or == "full" there are lines around each cell,
## otherwise (default) there are lines at the border and between the headings and
## the data.
##
## Example:
## a = rand(3,2);
## cnames = str2mat("c1, "c2");
## rnames = str2mat("r1, "r2", "r3");
## save2tex(a, "test.tex", cnames, rnames); 

## Author:  AW <Andreas.Weingessel@ci.tuwien.ac.at>
## Author:  MC <Michael.Creel@uab.es> (minor format modifications)

## Description:  Save to a file in a LaTeX tabular environment

function save2tex(tabledata, file, cname, rname, prec, align, lines)

  	if ((nargin == 1) || (nargin > 7))
		printf("usage: save2tex (X, file, cname, rname, prec, align, lines)\n");
  	endif

	if (nargin < 7) lines = 1; endif
	if (nargin < 6) align = "r"; endif
	if (nargin < 5) prec = -1; endif
	if (nargin < 4) rname = ""; endif
	if (nargin < 3) cname = ""; endif

	if (isstr (lines))
		if (strcmp (lines, "full"))
			lines = 2;
		elseif (strcmp (lines, "none"))
			lines = 0;
		else
			lines = 1;
		endif
	endif

	nr = rows(tabledata);
	nc = columns(tabledata);

	is_rn = columns(rname);
	is_cn = columns(cname);

	if ((is_rn) && (rows(rname) != nr))
		error ("save2tex: Numbers of rows and row names do not match.");
	endif

	if ((is_cn) && (rows(cname) != nc))
		error ("save2tex: Numbers of columns and column names do not match.");
	endif


	### open output file
	FN = fopen (file, "w");
	if (FN < 0)
		error ("save2tex: Can not open File %s", file);
	endif

	### create format line
	LSTR = "";
	LFSTR = "";
# 	if (lines)
# 		LSTR = "|";
# 		if (lines == 2)
# 			LFSTR = "|";
# 		endif
# 	endif
	STR = ["\\begin{tabular}{", LSTR];
	if (is_rn)
		STR = [STR, "l", LSTR];
	endif

	for i=1:nc
		STR = [STR, align, LFSTR];
	endfor
	
	if (lines == 2)
		STR = [STR, "}\n"];
	else
		STR = [STR, LSTR, "}\n"];
	endif
	fprintf (FN, STR);
	
	if (lines)
		fprintf (FN, "\\hline\n");
	endif

	
	### print column headers
	if (is_cn)
		if (is_rn)
			STR = "  & ";
		else
			STR = "  ";
		endif
		
		for i = 1:nc
			STR = [STR, cname(i,:), " & "];
		endfor
		
		les=length(STR);
		fprintf(FN, [STR(1:les-2), "\\\\\n"]);
		
		if (lines)
			fprintf (FN, "\\hline\\hline\n");
		endif
	endif

	### print data
	for i = 1:nr
		if (is_rn)
			STR = [rname(i,:), "  & "];
		else
			STR = "  ";
		endif
		
		if (prec >= 0)
			form = sprintf("%%.%df & ", prec);
		else
			form = "%f & ";
		endif
		
		for j = 1:nc
			STR = [STR, sprintf(form,tabledata(i,j))];
		endfor
		
		les=length(STR);
		
		if (i < nr)
			if (lines == 2)
				fprintf(FN, [STR(1:les-2), "\\\\\n\\hline\n"]);
			else
				fprintf(FN, [STR(1:les-2), "\\\\\n"]);
			endif
		else
			if (lines)
				fprintf(FN, [STR(1:les-2), "\\\\\n\\hline\n"]);
			else
				fprintf(FN, [STR(1:les-2), "\n"]);
			endif
		endif
	endfor
	
	### finish output
	fprintf(FN, "\\end{tabular}\n");
	fclose(FN);

endfunction
