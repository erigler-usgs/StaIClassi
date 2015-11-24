%%
%% Usage:
%%
%% [img_out, hdr_out] = img_resize_keep_aspect(img_in, pix_dim, hdr_in);
%%
%% This function resizes an image while maintaining the aspect ratio of each pixel. 
%% It may be useful when combining images of a common scene taken from different
%% instruments.
%%
%% INPUTS:
%%
%% The input img_in is the image to be resized, and whose pixel aspect ratio
%% (usually square) defines the aspect ratio of the output image, regardless
%% of the actual dimensions of the output image.
%%
%% The input pix_dim is a two-element vector describing the desired dimension
%% of the output image.  The smaller of the ratios of pix_dim to the dimensions
%% of img_in will define the actual scaling to be performed.  The dimension
%% corresponding to the larger of these ratios will have its edges padded with
%% nans in order to ensure that img_in remains centered.
%%
%% A simple check is performed if the optional input hdr_in is passed, just to
%% see if it is consistent with img_in.  Otherwise, it simply passes through
%% this function after its appropriate keyword values have been changed to 
%% reflect the new size of img_out.
%%
%%
%% OUTPUTS:
%%
%% The output img_out is the resized image.
%%
%% The output hdr_out is the modified header that should be used to hold meta
%% data for img_out.
%%
%%
%% AUTHOR(S):     E. Joshua Rigler (EJR)
%%
%%
%% CHANGELOG:     20 June, 2011 - First public release (EJR)
%%
%%

function [img_out, hdr_out] = img_resize_keep_aspect(img_in, pix_dim, hdr_in)
   
   % Get size of input image
   [nr,nc] = size(img_in);
  
   % If hdr_in is passed, assume it is a typical solar FITS header and
   % perform simple check to see if it actually belongs to img_in
   if (exist('hdr_in'))
      nx = hdr_in{strcmp(hdr_in(1:end,1), 'NAXIS1'), 2};
      ny = hdr_in{strcmp(hdr_in(1:end,1), 'NAXIS2'), 2};
      if (nx ~= nc || ny ~= nr)
         error('img_resize_keep_aspect: img_in dimensions do not match hdr_in dimensions');
      endif
   endif
   
   
   % The following mess essentially does:
   %  1) check to see which dimension is going to be grown faster;
   %  2) use the slower-growing dimension to define the step size for 
   %     remapping pixel coordinates;
   %  3) determines offsets from the stepsize and assumption that pixel
   %     coordinates correspond to pixel centers;
   %  4) centers new image and fills pixels beyond the original image's
   %     boundaries with nans
   if (pix_dim(2)/nc > pix_dim(1)/nr)
      noverp = nr/pix_dim(1);
      povern = pix_dim(1)/nr;
      pre1_postend = noverp*(povern-1)/2;
      
      [xi, yi] = meshgrid([1-pre1_postend:nr/(pix_dim(1)):pix_dim(2)/pix_dim(1)*nc+pre1_postend] - ...
                          (pix_dim(2) - pix_dim(1)) / ((pix_dim(1))/nc*2), ...
                          [1-pre1_postend:nr/(pix_dim(1)):nr+pre1_postend]);
   elseif (pix_dim(2)/nc < pix_dim(1)/nr)
      
      noverp = nc/pix_dim(2);
      povern = pix_dim(2)/nc;
      pre1_postend = noverp*(povern-1)/2;
      [xi, yi] = meshgrid([1-pre1_postend:nc/(pix_dim(2)):nc+pre1_postend],...
                          [1-pre1_postend:nc/(pix_dim(2)):pix_dim(1)/pix_dim(2)*nr+pre1_postend] - ...
                          (pix_dim(1) - pix_dim(2)) / ((pix_dim(2))/nr*2) );
   else
      noverp = nr/pix_dim(1);
      povern = pix_dim(1)/nr;
      pre1_postend = noverp*(povern-1)/2;
      [xi, yi] = meshgrid([1-pre1_postend:nc/(pix_dim(2)):nc+pre1_postend],...
                          [1-pre1_postend:nr/(pix_dim(1)):nr+pre1_postend]);
   endif
   
   [img_out] = imremap(img_in, xi, yi, 'bilinear', nan);
   
   
   % If hdr_in was passed, update standard solar FITS image header values
   if (exist('hdr_in'))
      
      hdr_out = hdr_in;
      
      hdr_out{strcmp(hdr_out(1:end,1), 'NAXIS1'), 2} = pix_dim(2);
      hdr_out{strcmp(hdr_out(1:end,1), 'NAXIS2'), 2} = pix_dim(1);
      
      
      % cannot assume CRPIX* are already pointing to center of image
# FIRST ATTEMPT WAS COMPLETELY WRONG, NOT SURE WHAT I WAS THINKING
#       hdr_out{strcmp(hdr_out(1:end,1), 'CRPIX1'), 2} = povern * ...
#                                                      (hdr_out{strcmp(hdr_out(1:end,1), 'CRPIX1'), 2} + ...
#                                                       max((pix_dim(1) - pix_dim(2)) / ((pix_dim(2))/nr*2), 0) );
#       hdr_out{strcmp(hdr_out(1:end,1), 'CRPIX2'), 2} = povern * ...
#                                                      (hdr_out{strcmp(hdr_out(1:end,1), 'CRPIX2'), 2} + ...
#                                                       max((pix_dim(2) - pix_dim(1)) / ((pix_dim(1))/nc*2), 0) );
      hdr_out{strcmp(hdr_out(1:end,1), 'CRPIX1'), 2} = ((hdr_out{strcmp(hdr_out(1:end,1), 'CRPIX1'), 2} * 2 -...
                                                         1) * povern + 1) / 2 + ...
                                                       povern * max( (pix_dim(1) - pix_dim(2)) / ...
                                                                     ((pix_dim(2))/nr*2), 0);
      hdr_out{strcmp(hdr_out(1:end,1), 'CRPIX2'), 2} = ((hdr_out{strcmp(hdr_out(1:end,1), 'CRPIX2'), 2} * 2 -...
                                                         1) * povern + 1) / 2 + ...
                                                       povern * max( (pix_dim(2) - pix_dim(1)) / ...
                                                                     ((pix_dim(1))/nc*2), 0);
      
      
      if any(strcmp(hdr_in(1:end,1), 'SOLAR_R'))
      
         hdr_out{strcmp(hdr_out(1:end,1), 'SOLAR_R'), 2} = povern * ...
                                                         hdr_out{strcmp(hdr_out(1:end,1), 'SOLAR_R'), 2};
      elseif any(strcmp(hdr_in(1:end,1), 'R_SUN'))

         hdr_out{strcmp(hdr_out(1:end,1), 'R_SUN'), 2} = povern * ...
                                                         hdr_out{strcmp(hdr_out(1:end,1), 'R_SUN'), 2};
      else
         error('img_resize_keep_aspect: pixel-radius of sun never specified in header');
      endif
      
      
      hdr_out{strcmp(hdr_out(1:end,1), 'CDELT1'), 2} = noverp * ...
                                                     hdr_out{strcmp(hdr_out(1:end,1), 'CDELT1'), 2};
      hdr_out{strcmp(hdr_out(1:end,1), 'CDELT2'), 2} = noverp * ...
                                                     hdr_out{strcmp(hdr_out(1:end,1), 'CDELT2'), 2};
      
      
      
   endif
   
endfunction