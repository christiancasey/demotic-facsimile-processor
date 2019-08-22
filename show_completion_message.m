function show_completion_message()
% Outputs a message to show user that the current script has finished

	sCaller = dbstack(1,'-completenames');
	if( length( sCaller ) > 0 )
		[~, strCurrentFile] = fileparts(sCaller(1).file);
		disp([ 'Done: ' strCurrentFile ]);
	end