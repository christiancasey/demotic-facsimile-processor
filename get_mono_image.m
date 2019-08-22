function img = get_mono_image(img)
% Returns a logical matrix the size of the RGB image: img
% Converts all grayscale values in img to 0 and colors to 1

	if(size(img,3) ~= 3)
		img = img(:,:,1);
		
		if(~all(img(:)==0) && ~all(img(:)==1))
			img = (img > mean(img(:)));
		end
		
		return;
	end
	
	img = (img(:,:,1) ~= img(:,:,2)) | (img(:,:,2) ~= img(:,:,3));