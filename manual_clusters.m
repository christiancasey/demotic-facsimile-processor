
%% Script Initialization
switchToCD;
clc;	clear;	close all;	
figure; figure(gcf);	whitebg('w');	colormap lapham;
plottools('off')

% Supress the warning raised by imshow() because why does this even exist?
warning('off','images:imshow:magnificationMustBeFitForDockedFigure');
set(0,'DefaultFigureWindowStyle','docked');

iFig = 0;

%% Load working data from previous step

load glyph_data02;


vGlyphs = vGlyphPolygons(:,2); % Just use the simple 100-pt polygons


%% Save each glyph as svg with glyph index as filename

figure(1);
clf;

vColor = lapham(nGlyphs);
for i = 1:nGlyphs
	
	figure(1);
	hold on;
	fill( vGlyphPolygons{i,1}(:,1), -vGlyphPolygons{i,1}(:,2), vColor(i,:), 'EdgeColor', 'none');
	axis equal;
	axis off;
	hold off;
	
	figure(i+100);
	clf;
	hold on;
	fill(vGlyphs{i}(:,1), -vGlyphs{i}(:,2), [0 0 0], 'EdgeColor', 'none');
	sprintf('glyphs%s%04i.eps',filesep,i)
	axis equal;
	axis off;
	getframe;
	
	
	print(sprintf('glyphs%s%04i.eps',filesep,i), '-depsc', '-painters');
end



%%

manual_cluster_assignments;
vManualClusterAssignments = vManualClusterAssignments(:,2);

X = mGlyphPolygons;

c = [ 1 2 1 1 1 1 1 2 2 2 2 3 1 3 1 3 3 ]';
v = [ 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 ]';

v = [ 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 ;
	  1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 4 ]';
  
  
  
  
  
  





