vTop = {
	'% !TEX program = xelatexmk'
	'\documentclass[final]{report}'
	''
	'\usepackage{graphbox}'
	'\graphicspath{{Figures/}}'
	''
	'\usepackage{epstopdf}'
	'\usepackage{longtable}'
	''
	'\newcommand\marginsize{2cm}'
	'\usepackage[letterpaper, margin=\marginsize]{geometry}'
	''
	'\begin{document}'
	'%%%% pagewidth 17.59cm'
	''
	'\begin{center}'
	'\begin{longtable}{|p{0.5cm}|p{1.2cm}|p{11.89cm}|r|}'
	'No. & Platonic Glyph & {Lines (Glyphs per line)} & {Total} \\'
	'\hline'
};

% 	'\includegraphics[align=t,width=0.8cm]{0001} & 1 (3), 2 (4), 3 (1), & Total \\'
% 	''
% 	'\hline'

vBtm = {
	'\end{longtable}'
	'\end{center}'
	'\end{document}'
};

strTop = sprintf('%s\n', vTop{:});
strBtm = sprintf('%s\n', vBtm{:});