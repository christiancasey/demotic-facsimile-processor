function strFilename = imagefilename(strName)

	strFilename = sprintf('wip_images%s%s %s.eps', ...
		filesep, strName, datetime('now','TimeZone','UTC','Format','y-MM-dd HH-mm-ss'));
