%%
%% Usage:
%%
%% [xsub, ysub] = get_coord_img(img_in[, xsub, ysub])
%%
%% This is a semi-gui utility function for collecting pixel coordinates from
%% an input image.  It may be used to select training data for classifying
%% image pixels.
%%
%%
%% INPUTS:
%%
%% The input argument img_in is a stack of 2D image arrays that contains desired
%% pixels.  Unselected pixels will be displayed in grey scale, while selected
%% pixels will be displayed in a color scale derived from whatever the colormap
%% was when the function was called (be certain the colormap set when this is
%% called is not grayscale too).  Pixels are selected by clicking on the
%% plotted image:
%%
%% * A left-button mouse click will select a vertex to a polygon;
%% * A middle-button mouse click will select the last vertex of a polygon,
%%   then include all pixels to the interior of the polygon in the pixel mask;
%%   if the middle button is clicked with no vertices selected, it will select
%%   a single pixel; if the middle button is clicked as the 2nd vertex, the
%%   vertices will define the min/max corners of a rectangular region, and NOT
%%   a line (the latter may not be the ideal behavior, since I have often found
%%    that the ability to draw a line would be more convenient that drawing a
%%    box -EJR 3/2012)
%% * The first right-button mouse click will select the first corner of a 
%%   rectangular region to zoom in to, and the second right-button mouse click 
%%   will select the second corner of a rectangular region to zoom into.  
%%
%% There are also some important keyboard 'clicks':
%%  
%%  r     - exactly like a middle-button click, except the selected pixels
%%          are removed from the pixel mask;
%%  z     - zoom ~10% centered on the pixel position;
%%  u     - unzoom ~10% centered on the pixel position;
%%  p     - go to previous zoom level;
%%  n     - go to next zoom level;
%%  h     - return to original, or 'home', zoom level;
%%  j     - shift the current pixel value range down (brighten);
%%  l     - shift the current pixel value range up (darken);
%%  i     - expand the current pixel value range (reduce contrast);
%%  k     - contract the current pixel value range (enhance contrast);
%%  1->0  - switch layer of image stack being displayed (zoom, brightness,
%%          and contrast levels are remembered for each layer)
%%  Enter - end selection process and return pixel coordinates;
%%
%% FIXME_01:  this function relies on gnuplot, the default plotting program
%%            for Octave, and a standard Matlab/Octave function called ginput
%%            to read mouse clicks when the mouse cursor is positioned over a 
%%            plot window. It provides some basic real-time instructions by 
%%            printing text to the terminal, but it would be nice if 
%%            instructions could be printed to the plot window instead.
%%
%% The optional inputs xsub and ysub are the x and y subscripts of pixels
%% that have already been selected.  If only xsub is passed, and is a vector, 
%% get_coord assumes that these are linear indices into the 2D image.  These
%% are interpreted in a column major manner, which is standard for Fortran,
%% Matlab, and Octave, but not for C, C++, or IDL.  If only xsub is passed,
%% but it is a 2-column array, assume the first column holds x coordinates,
%% and the second column holds y coordinates for selected pixels in img_in.
%%
%%
%% OUTPUTS:
%%
%% Output arrays xsub and ysub are the x and y subscripts of pixels were
%% selected.  If only 1 output array is requested, the x and y subscripts
%% are converted to linear indices.
%% FIXME_02:  this should also be consistent with how optional inputs
%%            xsub and ysub were passed.
%%
%%
%% AUTHOR(S):     E. Joshua Rigler (EJR)
%%
%%
%% CHANGELOG:     17 Jun, 2011 - First public release (EJR)
%%                19 Jan, 2012 - Major overhaul to allow selecting/removing
%%                               polygons; also zoom levels are tracked, thus
%%                               allowing ALL previous/next zooms to be visited;
%%                               finally, the name was changed to get_coord_img.m
%%                               to differentiate from get_coord_ts.m for 1-D 
%%                               time series...need to update other routines 
%%                               to exploit new version.
%%                11 Mar. 2012 - Modifications to allow ithresh input; this
%%                               is really just a usability enhancement
%%                15 Mar, 2012 - Allow 3D img_in, and introduce controls to
%%                               select channel. Also introduce controls
%%                               for 'contrast' and 'brightness'. Removed
%%                               ability to specify colormap or color axis,
%%                               since these really just complicated the code
%%                               and were better controlled interactively anyway.
%%

function [xsub, ysub] = get_coord_img(img_in, xsub, ysub)

   % the following is NOT ML compatible, but helps with a number of issues that
   % crop up in this interactive function:
   % 1) pressing 'q' closes a plot window in gnuplot, but the x11 terminal has
   %    a way to prevent this with the -ctrlq option...unfortunately this doesn't
   %    work with other terminals
   % 2) the wxt terminal doesn't align characters very well, which results in
   %    offset '+' marks when specifying polygon vertices
   % 3) some global, and even outside environment variables get modified, so 
   %    the unwind_protect block is used to reset these if there is an error.
   unwind_protect
   gplot_term_orig = getenv('GNUTERM');
   [gplot_bin_orig, gplot_args_orig] = gnuplot_binary;
   if (strfind(tolower(computer()), 'linux'))
    putenv('GNUTERM','x11');
    gnuplot_binary(gplot_bin_orig, gplot_args_orig{:}, "-ctrlq");
   elseif (strfind(tolower(computer()), 'mingw'))
     putenv('GNUTERM','windows');
    gnuplot_binary(gplot_bin_orig, gplot_args_orig{:}, "");
   endif
   
      
   % get image dimensions
   [sy, sx, sz] = size(img_in);
   
   % get current colormap to restore on exit
   cmap_restore = colormap;
   
   
   % check inputs
   if (nargin < 1)
      error('Not enough inputs');
   endif

   
   % initialize xsub and ysub if no initial *sub arrays are passed in
   if (nargin == 1)
   
      xsub = ysub = [];
      
   endif
   
   
   % if only xsub is passed... 
   if (nargin == 2)
   
      if (isvector(xsub))
         % ...and it is a vector, assume it is
         % comprised of linear indices
         [ysub, xsub] = ind2sub([sy,sx], xsub);
         
      elseif (size(xsub,2) == 2)
         % ...or else it is a 2-column matrix
         % holding x-y pairs already (NOT i-j pairs)
         % (this is probably a bad idea)
         ysub = xsub(:,2);
         xsub(:,2) = [];
         
      else
         error('Invalid index or subscript array');
      endif
            
   endif
   
   
   % xsub and ysub must be vectors of similar length
   if (nargin == 3)
      
      if (~isvector(xsub) || ~isvector(ysub) || (size(xsub) ~= size(ysub)) )
         error('Invalid index or subscript array');
      endif
      
   endif

   
   if (nargin > 3)
      error('Too many inputs');
   endif
   

   
   %
   % start processing
   %
   
   
   % interpolate current colormap to 240 levels (leaves 16 black for a clean division
   % between the 256 B/W levels and 256 color levels of the image)
   cmap = colormap;
   cmap_240 = interp1([1:rows(cmap)], cmap, [1:(rows(cmap)-1)/(240-1):rows(cmap)]);
   colormap([zeros(16,3);gray(240);zeros(16,3); cmap_240]);
   
   
   % determine a reasonable initial intensity range for each channel
   % (this can be adjusted with ik controls)
   channel_means = mean(reshape(img_in,sx*sy,sz));
   channel_stds = std(reshape(img_in,sx*sy,sz));
   
   channel_mins = channel_means - 1*channel_stds;
   channel_maxes = channel_means + 1*channel_stds;
   
   
#    # this seemed necessary at one point, but I cannot remember why;
#    # delete it if recent mods don't break anything   
#    if (pix_min > min(img_in(:,:,1)(:)))
#       pix_min = max(img_in(:,:,1)(img_in(:,:,1) <= pix_min)(:));
#    endif
#    
#    if (pix_max < max(img_in(:,:,1)(:)))
#       pix_max = min(img_in(:,:,1)(img_in(:,:,1) >= pix_max)(:));
#    endif
   
   
   % always start by displaying first channel...it's easy enough to change using 
   # interactive keyboard controls
   channel = 1;
   
   
   % trim/truncate intensities to min/max range for this channel
   % (FIXME: the following block of code is repeated several times below;
   %         it would be good to write a subfunction, and call it when
   %         necessary, just to keep the code clean and maintainable)
   img_in_tmp = img_in(:,:,channel);
   img_in_tmp(img_in_tmp > channel_maxes(channel)) = channel_maxes(channel);
   img_in_tmp(img_in_tmp < channel_mins(channel)) = channel_mins(channel);
   
   % scale intensities to [0-255], all integer indices into bottom half of colormap
   img_sc = round([(img_in_tmp - channel_mins(channel)) / ...
                   (channel_maxes(channel) - channel_mins(channel))] * 255) + 1;
   
   % an obscure bug seems to exist in Gnuplot wherein integer indices do not
   % align perfectly with the corresonding row of a color table. As a work-
   % around, subtract 1 here, then convert all zeros to 1. Revisit someday to
   % see if this has been fixed.
   img_sc = img_sc - 1;
   img_sc(img_sc <= 0) = 1;
   
   
   % initialize some important variables outside the loop(s)
   px=0;
   py=0;
   button=0;
   channel=1;
   vertex_x = [];
   vertex_y = [];
   plot_refresh = 0;
   xax = [0 sx];
   yax = [0 sy];
   ax_idx = 1;
   xax_stack = xax;
   yax_stack = yax;
   img_mask = zeros(sy,sx);
   
   
   % prompt/instruct the user
   printf('left-click = next vertex; middle-click = close & select polygon; r = close & remove polygon \n');
   printf('right-click = zoom box; z = zoom; u = unzoom; p = previous zoom; n = next zoom; h = home zoom\n');
   printf('i|k = contrast; j|l=brightness; 1->0=channel displayed; d=debug mode...then dbcont|dbquit');
   printf('Enter/Return = done\n');
   fflush(stdout);

   while (1) 
      
      % change pixel values in img_mask for preview plot
      img_mask(:) = 0;
      img_mask(sub2ind([sy,sx], ysub, xsub)) = 256;
      

      % more responsive if only region of interest is plotted, not the entire image,
      % followed by a change in the axes.
      xpix_coords = [max(1, xax(1)):min(sx, xax(2))];
      ypix_coords = [max(1, yax(1)):min(sy, yax(2))];
      image(xpix_coords, ypix_coords, [img_sc + img_mask](ypix_coords, xpix_coords) );
      axis([xax yax], 'xy', 'equal');
      
      % replotting the entire image every time a vertex is added to a polygon is not
      % necessary, and actually quite distracting due to flickering
      if plot_refresh
        plot_refresh = 0;
        text(vertex_x, vertex_y, '+', 'color', [1 1 0], 'horizontalalignment', 'center', 'verticalalignment', 'middle');
      else
        printf('\rSelect pixel, corner of rectangle, or vertex of polygon...\n')
      end
      
      ## This was part of amisguided attempt to clear lines of text so that output didn't 
      ## just scroll off the top of the terminal. The following bit of code is handy for
      ## deleting a line of text though, so I'm leaving it commented out for possible future
      ## use.
      # printf('                                                                                '); % 80 spaces
      # printf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); % 40 backspaces
      # printf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); % 40 backspaces
      fflush(stdout);
      
      
      # a separate loop for pixel selection prevents refresh every time a new
      # vertex pixel is selected...unless it is required for a rezoom
      while (2)
        
        % read in a single mouse/keyboard event using ginput
        [px, py, button] = ginput(1);
        px = round(px);
        py = round(py);
        
        % determine which button was pressed, and act accordingly
        switch button
          
          
          case 2 % finish aggregating 'clicks' and generate next mask
              
              vertex_x = [vertex_x, px];
              vertex_y = [vertex_y, py];
              
              printf('Vertex {%d, %d} selected; generating mask...\n', px, py);
              fflush(stdout);
              
              switch numel(vertex_x)
                
                case 1 % simply add a single pixel coordinate to mask coordinates
                  
                  xsub = [xsub; px];
                  ysub = [ysub; py];
                  
                  xadd = px;
                  yadd = py;

                case 2 % draw a box where 2 pixels define min/max corners, add box's
                       % interior to mask coordinates...this could be done using the
                       % poly2mask function, but since the version distributed via
                       % Octave-Forge is broken IMHO, let's create our own mask for
                       % rectangular boxes
                  
                  xrange = sort(vertex_x);
                  yrange = sort(vertex_y);
                  [xadd, yadd] = meshgrid([xrange(1):xrange(2)], [yrange(1):yrange(2)]);
                  
                  
                otherwise % draw a polygon using more than 2 vertex coordinates, add 
                          % polygon's interior to mask coordinates...this uses the
                          % poly2mask function from the Octave-Forge Image toolbox.
                          % EJR made a minor "fix" to this function, which he uses,
                          % but the official version really isn't so broken that it
                          % is unusable.
                  
                  mask_tmp = poly2mask(vertex_x, vertex_y, size(img_mask,1), size(img_mask,2));
                  [yadd, xadd] = find(mask_tmp);
                  
                
              end % switch numel(vertex_x)
              
              
              % update xsub and ysub
              xsub = [xsub; xadd(:)];
              ysub = [ysub; yadd(:)];
              
              
              % clear vertex coordinates
              vertex_x = [];
              vertex_y = [];

              
              break;
              
              
          case 1 % select next vertex

              vertex_x = [vertex_x, px];
              vertex_y = [vertex_y, py];
              
              printf('Vertex {%d, %d} selected; choose next vertex...\n', px, py);
              fflush(stdout);
              
              % plot a marker for each vertex...these will be erased after last vertex is selected
              text(px, py, '+', 'color', [1 1 0], 'horizontalalignment', 'center', 'verticalalignment', 'middle');
              
          case toascii('r') % remove pixels from list
              
              vertex_x = [vertex_x, px];
              vertex_y = [vertex_y, py];
              
              printf('Vertex {%d, %d} selected; removing mask...\n', px, py);
              fflush(stdout);
              
              switch numel(vertex_x)
                
                case 1 % simply add a single pixel coordinate to mask coordinates
                  
                  xrem = px;
                  yrem = py;
                  
                case 2 % draw a box where 2 pixels define min/max corners, add box's
                       % interior to mask coordinates...this could be done using the
                       % poly2mask function, but since the version distributed via
                       % Octave-Forge is broken IMHO, let's create our own mask for
                       % rectangular boxes
                  
                  xrange = sort(vertex_x);
                  yrange = sort(vertex_y);
                  [xrem, yrem] = meshgrid([xrange(1):xrange(2)], [yrange(1):yrange(2)]);
                  
                otherwise % draw a polygon using more than 2 vertex coordinates, add 
                          % polygon's interior to mask coordinates...this uses the
                          % poly2mask function from the Octave-Forge Image toolbox.
                          % EJR made a minor "fix" to this function, which he uses,
                          % but the official version really isn't so broken that it
                          % is unusable.
                  
                  mask_tmp = poly2mask(vertex_x, vertex_y, size(img_mask,1), size(img_mask,2));
                  [yrem, xrem] = find(mask_tmp);
                
              end % switch numel(vertex_x)
              
              
              % remove xbox from xsub using setdiff.m 
              xsub = setdiff([xsub, ysub], [xrem(:), yrem(:)], 'rows');
              ysub = xsub(:,2);
              xsub(:,2) = [];
              
              
              % reset vertex coordinates
              vertex_x = [];
              vertex_y = [];
      
              break;
              
              
          case 3 % save current axis, and 'zoom' to new axis

              px_f = px;
              py_f = py;
              
              printf('Zoom box corner {%d, %d} selected; now choose opposite corner...\n', px, py);
              fflush(stdout);
              [px, py, button] = ginput(1);
              px = round(px);
              py = round(py);
              
              if (button == 3)
                
                % remove any axes on stack after current axis
                xax_stack(ax_idx+1:end,:) = [];
                yax_stack(ax_idx+1:end,:) = [];
                
                % generate new axis
                xax = sort([px_f px]);
                yax = sort([py_f py]);
                
                % append new axis
                xax_stack = [xax_stack; xax];
                yax_stack = [yax_stack; yax];
                
                % increment ax_idx
                ax_idx = ax_idx + 1;
                
              
                plot_refresh = 1 % force replot
                break;
                
              else
                printf('\nInvalid key or button, try again...\n');
              endif

          
          case toascii('z') % zoom in ~10%, centered on {px,py}
              
              % remove any axes on stack after current axis
              xax_stack(ax_idx+1:end,:) = [];
              yax_stack(ax_idx+1:end,:) = [];
              
              % generate new axis
              xax = round([px - 1/1.1 * diff(xax)/2, px + 1/1.1 * diff(xax)/2]);
              yax = round([py - 1/1.1 * diff(yax)/2, py + 1/1.1 * diff(yax)/2]);
              
              % append new axis
              xax_stack = [xax_stack; xax];
              yax_stack = [yax_stack; yax];
              
              
              % increment ax_idx
              ax_idx = ax_idx + 1;
              
              
              plot_refresh = 1; % force replot
              break;
                
              
          case toascii('u') % zoom out ~10%, centered on {px,py}
              
              % remove any axes on stack after current axis
              xax_stack(ax_idx+1:end,:) = [];
              yax_stack(ax_idx+1:end,:) = [];
              
              % generate new axis
              xax = round([px - 1.1 * diff(xax)/2, px + 1.1 * diff(xax)/2]);
              yax = round([py - 1.1 * diff(yax)/2, py + 1.1 * diff(yax)/2]);
              
              % append new axis
              xax_stack = [xax_stack; xax];
              yax_stack = [yax_stack; yax];
              
              
              % increment ax_idx
              ax_idx = ax_idx + 1;
              
              
              plot_refresh = 1; % force replot
              break;
                
                            
          case toascii('p') % switch to previous axis on stack
              
              % update ax_idx
              ax_idx = max(1, ax_idx - 1);
              
              
              % move to previous axis on stack
              xax = xax_stack(ax_idx,:);
              yax = yax_stack(ax_idx,:);
              
              plot_refresh = 1; % force replot
              break;
              
          
          case toascii('n') % switch to next axis on stack
              
              % update ax_idx
              ax_idx = min(size(xax_stack,1), ax_idx + 1);
              
              
              % move to next axis on stack
              xax = xax_stack(ax_idx,:);
              yax = yax_stack(ax_idx,:);
              
              plot_refresh = 1; % force replot
              break;
              
          
          case toascii('h') % switch to original/default zoom
              
              ax_idx = 1;
              
              
              % move to first axis on stack
              xax = xax_stack(ax_idx,:);
              yax = yax_stack(ax_idx,:);
              
              plot_refresh = 1; % force replot
              break;
              
          
          case toascii('i') % expand color axis (reduce contrast)
            
            channel_range = channel_maxes(channel) - channel_mins(channel);
            channel_mean = mean([channel_maxes(channel), channel_mins(channel)]);
            channel_maxes(channel) = channel_mean + channel_range*(3/4);
            channel_mins(channel) = channel_mean - channel_range*(3/4);
            
            img_in_tmp = img_in(:,:,channel); % gotta do this or you get saturation craziness
            
            % trim/truncate intensities to min/max range for this channel
            img_in_tmp(img_in_tmp > channel_maxes(channel)) = channel_maxes(channel);
            img_in_tmp(img_in_tmp < channel_mins(channel)) = channel_mins(channel);
            
            % scale intensities to [0-255], all integer indices into bottom half of colormap
            img_sc = round([(img_in_tmp - channel_mins(channel)) / ...
                            (channel_maxes(channel) - channel_mins(channel))] * 255) + 1;
            
            % an obscure bug seems to exist in Gnuplot wherein integer indices do not
            % align perfectly with the corresonding row of a color table. As a work-
            % around, subtract 1 here, then convert all zeros to 1. Revisit someday to
            % see if this has been fixed.
            img_sc = img_sc - 1;
            img_sc(img_sc <= 0) = 1;

            plot_refresh = 1; % force replot
            break;

            
          case toascii('k') % contract color axis (increase contrast)
            
            channel_range = channel_maxes(channel) - channel_mins(channel);
            channel_mean = mean([channel_maxes(channel), channel_mins(channel)]);
            channel_maxes(channel) = channel_mean + channel_range*(1/3);
            channel_mins(channel) = channel_mean - channel_range*(1/3);
            
            img_in_tmp = img_in(:,:,channel); % gotta do this or you get saturation craziness
            
            % trim/truncate intensities to min/max range for this channel
            img_in_tmp(img_in_tmp > channel_maxes(channel)) = channel_maxes(channel);
            img_in_tmp(img_in_tmp < channel_mins(channel)) = channel_mins(channel);
            
            % scale intensities to [0-255], all integer indices into bottom half of colormap
            img_sc = round([(img_in_tmp - channel_mins(channel)) / ...
                            (channel_maxes(channel) - channel_mins(channel))] * 255) + 1;
            
            % an obscure bug seems to exist in Gnuplot wherein integer indices do not
            % align perfectly with the corresonding row of a color table. As a work-
            % around, subtract 1 here, then convert all zeros to 1. Revisit someday to
            % see if this has been fixed.
            img_sc = img_sc - 1;
            img_sc(img_sc <= 0) = 1;

            plot_refresh = 1; % force replot
            break;
            
          case toascii('l') % shift color axis up 1/4 current color range
            
            channel_range = channel_maxes(channel) - channel_mins(channel);
            channel_maxes(channel) = channel_maxes(channel) + channel_range/4;
            channel_mins(channel) = channel_mins(channel) + channel_range/4;
            
            img_in_tmp = img_in(:,:,channel); % gotta do this or you get saturation craziness
            
            % trim/truncate intensities to min/max range for this channel
            img_in_tmp(img_in_tmp > channel_maxes(channel)) = channel_maxes(channel);
            img_in_tmp(img_in_tmp < channel_mins(channel)) = channel_mins(channel);
            
            % scale intensities to [0-255], all integer indices into bottom half of colormap
            img_sc = round([(img_in_tmp - channel_mins(channel)) / ...
                            (channel_maxes(channel) - channel_mins(channel))] * 255) + 1;
            
            % an obscure bug seems to exist in Gnuplot wherein integer indices do not
            % align perfectly with the corresonding row of a color table. As a work-
            % around, subtract 1 here, then convert all zeros to 1. Revisit someday to
            % see if this has been fixed.
            img_sc = img_sc - 1;
            img_sc(img_sc <= 0) = 1;

            plot_refresh = 1; % force replot
            break;

            
          case toascii('j') % shift color axis down 1/4 current color range
            
            channel_range = channel_maxes(channel) - channel_mins(channel);
            channel_maxes(channel) = channel_maxes(channel) - channel_range/4;
            channel_mins(channel) = channel_mins(channel) - channel_range/4;
            
            img_in_tmp = img_in(:,:,channel); % gotta do this or you get saturation craziness
            
            % trim/truncate intensities to min/max range for this channel
            img_in_tmp(img_in_tmp > channel_maxes(channel)) = channel_maxes(channel);
            img_in_tmp(img_in_tmp < channel_mins(channel)) = channel_mins(channel);
            
            % scale intensities to [0-255], all integer indices into bottom half of colormap
            img_sc = round([(img_in_tmp - channel_mins(channel)) / ...
                            (channel_maxes(channel) - channel_mins(channel))] * 255) + 1;
            
            % an obscure bug seems to exist in Gnuplot wherein integer indices do not
            % align perfectly with the corresonding row of a color table. As a work-
            % around, subtract 1 here, then convert all zeros to 1. Revisit someday to
            % see if this has been fixed.
            img_sc = img_sc - 1;
            img_sc(img_sc <= 0) = 1;

            plot_refresh = 1; % force replot
            break;

            
          
          case {toascii('1') toascii('2') toascii('3') toascii('4') toascii('5') ...
                toascii('6') toascii('7') toascii('8') toascii('9') toascii('0') }
            
            if (mod(channel,10) == str2num(char(button)))
              channel = channel+10;
            else
              channel = str2num(char(button));
            endif
            
            if (channel > size(img_in,3))
              printf("\nInvalid channel specified, try again...\n");
              fflush(stdout);
              continue;
            else
              img_in_tmp = img_in(:,:,channel);
            endif
            
            % trim/truncate intensities to min/max range for this channel
            img_in_tmp(img_in_tmp > channel_maxes(channel)) = channel_maxes(channel);
            img_in_tmp(img_in_tmp < channel_mins(channel)) = channel_mins(channel);
            
            % scale intensities to [0-255], all integer indices into bottom half of colormap
            img_sc = round([(img_in_tmp - channel_mins(channel)) / ...
                            (channel_maxes(channel) - channel_mins(channel))] * 255) + 1;
            
            % an obscure bug seems to exist in Gnuplot wherein integer indices do not
            % align perfectly with the corresonding row of a color table. As a work-
            % around, subtract 1 here, then convert all zeros to 1. Revisit someday to
            % see if this has been fixed.
            img_sc = img_sc - 1;
            img_sc(img_sc <= 0) = 1;

            plot_refresh = 1; % force replot
            break;
            
          case 13 % 13 is the retrun/enter key...'q' would be better, but by default
                  % gnuplot closes the plot window when 'q' is pressed
              break;
          
          case toascii('d') % enter debug mode..in Linux, this could be done using
                                    % ctrl-c and setting debug_on_interrupt, but ctrl-c is  
                                    % not reliable in Windows
              keyboard;
          
          otherwise
              
              printf('\nInvalid key or button, try again...\n');
              
        end % switch button
        fflush(stdout);
        
        
      end % while(2)
      
      
      % the following conditionals are a consequence of the fact that there is no
      % way to break out of more than one loop at a time in Octave.
      
      % if the plot needs to be updated, skip to the end of this loop
      if plot_refresh
        continue;
      endif
      
      % if the return/enter key was hit inside the inner loop, a break was issued...
      % this only breaks the inner loop, so we must break the outer loop now.
      if (button == 13) 
        break; 
      endif
      
      
      % remove duplicate coordinates inside the loop...this is not as
      % efficient, but it is necessary for plotting the updated image;
      % also, make sure xsub actually exists
      if (ismatrix(xsub))
        xsub = unique([xsub,ysub], 'rows');
        ysub = xsub(:,2);
        xsub(:,2) = [];
      endif
      
   end % while(1)

            
   printf('\n');
   fflush(stdout);


   % make sure xsub and ysub are row vectors, NOT column vectors
   ysub = reshape(ysub, 1, numel(ysub));
   xsub = reshape(xsub, 1, numel(xsub));


   % convert to linear indices if only one or no output arguments are requested
   if (nargout < 2)
      xsub = sub2ind(size(img_in), ysub, xsub);
   endif

   unwind_protect_cleanup
   
   colormap(cmap_restore);
   
   putenv('GNUTERM', gplot_term_orig);
   gnuplot_binary(gplot_bin_orig, gplot_args_orig{:})
   
   end_unwind_protect
   
end
