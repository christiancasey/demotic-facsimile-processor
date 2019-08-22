function switchToCD(bAddToPath)

	if(~exist('bAddToPath'))
		bAddToPath = false;
	end
	
	sCaller = dbstack(1,'-completenames');
	if( length( sCaller ) > 0 )
		strPath = fileparts(sCaller(1).file);
		cd(strPath);
		if( bAddToPath )
			strPD = pathdef;
			if(isempty( strfind( strPD, [ strPath ';' ]) ) )

				fin = fopen([ matlabroot '\toolbox\local\pathdef.m' ],'r');
				strFile = char(fread(fin,inf,'char'))';
				fclose(fin);
				vPDFile = get_tokens(strFile, '%%% BEGIN ENTRIES %%%', 1);
				
				fout = fopen([ matlabroot '\toolbox\local\pathdef.m' ],'w');
				fprintf(fout, '%s', vPDFile{1});
				fprintf(fout, '%s', '%%% BEGIN ENTRIES %%%');
				fprintf(fout, '\n\t''%s;'', ...', strPath);
				fprintf(fout, '%s', vPDFile{2});
				fclose(fout);
				
				path(path, strPath);
			end
			
		end
	
	end
		