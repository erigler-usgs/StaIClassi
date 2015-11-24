%%
%% Usage:
%%
%% [ts_idx] = get_coord_ts(ts_in[, ts_idx, rgb_color])
%%
%% This is a utility function for collecting time indices from a time series.
%%
%%
%% INPUTS:
%%
%% The input argument ts_in is a time series that contains some desired values.
%% Unselected values will be displayed normally as a black (or rgb_colored) line
%% on a white background, while selected values are displayed using the inverse.
%% Time ranges are selected by clicking on the plotted image:
%%
%% * A left-button mouse click will select the first of two time inteval boundaries;
%% * A middle-button mouse click will select the second of two time interval boundaries;
%% * The first right-button mouse click will select the first boundary of a
%%   time interval zoom in on, and the second right-button mouse click
%%   will select the second boundary. The y axis is unaffected by this.
%%
%% There are also a few important keyboard 'clicks':
%%
%%  r     - exactly like a middle-button click, except the selected interval
%%          is removed from the pixel mask;
%%  z     - zoom ~10% centered on the time position;
%%  u     - unzoom ~10% centered on the time position;
%%  p     - go to previous zoom level;
%%  n     - go to next zoom level;
%%  h     - return to original, or 'home', zoom level;
%%  j     - shift the current baseline value range down;
%%  l     - shift the current baseline value range up;
%%  i     - expand the current yaxis;
%%  k     - contract the current yaxis;
%%  1->0  - switch channel of the time series stack being displayed
%%  Enter - end selection process and return time indices;
%%
%% FIXME_01:  this function relies on gnuplot, the default plotting program
%%            for Octave, and a standard Matlab/Octave function called ginput
%%            to read mouse clicks when the mouse cursor is positioned over a
%%            plot window. It provides some basic real-time instructions by
%%            printing text to the terminal, but it would be nice if
%%            instructions could be printed to the plot window instead.
%%
%% The optional input ts_idx contains time indices that have already been
%% selected, thus allowing this routine to be called repeatedly.
%%
%%
%% OUTPUTS:
%%
%% Output array tsub is an array of the time array indices, relative to ts_in
%% only, that were selected.
%%
%%
%% AUTHOR(S):     E. Joshua Rigler (EJR)
%%
%%
%% CHANGELOG:     16 Dec, 2011 - First public release (EJR)
%%                22 Mar, 2012 - Major update to make it compatible with
%%                               get_coord_img.m. This included:
%%                               - zoom stack/memory
%%                               - switch between multiple channels
%%                               - scale and shift yaxis dynamically
%%                               - several additional keyboard 'clicks'
%%

function [ts_idx] = get_coord_ts(ts_in, ts_idx, rgb_color)

   % check inputs
   if (nargin < 1)
      error('Not enough inputs');
   endif


   % make certain ts_in is a vector (or vectors stacked along 3rd dimension)
   [si, sj, sk] = size(ts_in);
   if ~(si==1 || sj==1)
     error('One of the first two dimensions must be singular');
   endif

   % force ts_in to be column vector(s)
   ts_in = reshape(ts_in, max(si,sj), 1, sk);
   st = rows(ts_in); % get length of time series


   if (nargin == 1)
      % initialize ts_idx and rgb_color if they are not passed
      ts_idx = [];
      rgb_color = [0 0 0];
   endif


   if (nargin == 2)
    % initialize rgb_color
    rgb_color = [0 0 0];
   endif


   if (nargin == 3)
     if ~(isvector(rgb_color) && numel(rgb_color) == 3)
       error('get_coord_ts.m: rgb_color must be a 3-element vector');
     endif
   endif


   if (nargin > 3)
      error('Too many inputs');
   endif


   %
   % start processing
   %

   % initialize some important variables outside the loop(s)
   px=0;
   py=0;
   button=0;
   channel=1;
   vertex_x = [];
   vertex_y = [];
   plot_refresh = 0;
   xax = [0 st+1];
   xax_stack = xax;
   for i=1:size(ts_in,3)
     yax_stack(1,:,i) = [min(ts_in(:,:,i)) max(ts_in(:,:,i))];
   endfor
   yax = yax_stack(1,:,channel);
   ax_idx = 1;
   ts_mask = zeros(st,1);
   bound1 = [];


   % prompt/instruct the user
   printf('left-click = 1st time boundary; middle-click = single index, or 2nd time boundary; r = unselect/remove range\n');
   printf('right-click = rect. zoom; z = zoom; u = unzoom; p = previous zoom; n = next zoom; h = home zoom\n');
   printf('i|k = y range; j|l=y baseline; 1->0=channel displayed; d=debug mode...then dbcont|dbquit\n');
   printf('Enter/Return = done\n');
   fflush(stdout);


   % associate pixels with predefined classes
   while (1)


      % change pixel values in img_mask for preview plot
      ts_mask(:) = nan;
      ts_mask(ts_idx) = 1;


      % this plots faster if/when you zoom in on features of interest
      plot_mask = ([1:numel(ts_in(:,1,channel))] >= xax(1) & ...
                   [1:numel(ts_in(:,1,channel))] <= xax(2));

      plot(find(plot_mask), ts_in(:,1,channel)(plot_mask), ...
           'color', rgb_color, ...
           'linewidth', 1.5);
      axis([xax yax]);
      hold on;


      % don't draw a patch for each index, but rather ranges of indices
      first_first = min(min(0,ts_idx)); % this is a total kludge to ensure an empty matrix or 0
      last_last = min(min(0,ts_idx)) + length(ts_idx); % this is a continuation of total kludge
      first_idx = [first_first; find(diff(ts_idx) ~= 1)] + 1;
      last_idx = [find(diff(ts_idx) ~= 1); last_last];


      for i=1:numel(first_idx)
        if ~(ts_idx(last_idx(i)) < xax(1) || ts_idx(first_idx(i)) > xax(2))
          patch([ts_idx(first_idx(i))-.5 ts_idx(last_idx(i))+.5 ts_idx(last_idx(i))+.5 ts_idx(first_idx(i))-.5], ...
                [yax(1)+.01*diff(yax) yax(1)+.01*diff(yax) yax(2)-.01*diff(yax) yax(2)-.01*diff(yax)], ...
                rgb_color, 'linewidth', 0);
        endif
      endfor



      plot(find(plot_mask), [ts_in(:,:,channel) .* ts_mask](plot_mask), 'color', 1-rgb_color, 'linewidth', 1.5);
      hold off;


#       printf('\rSelect pixel or first corner of rectangle...')
#       printf('                                                                                '); % 80 spaces
#       printf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); % 40 backspaces
#       printf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); % 40 backspaces
#       fflush(stdout);


      # a separate loop for pixel selection prevents refresh every time a new
      # vertex pixel is selected...unless it is required for a rezoom
      while (2)

        % read in a single event using ginput
        [px, py, button] = ginput(1);
        px = round(px);
        py = round(py);

        % check which button was pressed, and act accordingly
        if isempty(button)
          % cannot use isempty() in case statement, so just check for the
          % return/empty key prior to the switch-case block
          break;
        else

          switch button


            case 1 % set first time span boundary

                bound1 = px;

                printf('First boundary {%d} set; now select other boundary...\n', px);
                fflush(stdout);


            case 2 % set second time boundary, or return a single time index

                if (numel(bound1) == 0)

                  ts_idx = [ts_idx; px];

                else

                  printf('');
                  xrange = sort([bound1,px]);
                  tbox = [xrange(1):xrange(2)];
                  ts_idx = [ts_idx; tbox(:)];

                  % reset bound1
                  bound1 = [];
                endif


                plot_refresh = 1; % force replot
                break;


            case toascii('r') % remove pixel range from list
                              % (would be nice to remove a single pixel too)

                if (numel(bound1) == 0)

                  ts_idx = setdiff(ts_idx, px, 'rows');

                else

                  printf('');
                  xrange = sort([bound1,px]);
                  tbox = [xrange(1):xrange(2)];
                  ts_idx = setdiff(ts_idx, tbox(:), 'rows');

                  % reset bound1
                  bound1 = [];

                  plot_refresh = 1; % force replot
                  break;

                endif


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
                  yax_stack(ax_idx+1:end,:,:) = [];

                  % generate new x axis; keep all old y axes
                  xax = sort([px_f px]);
                  yax = yax;

                  % append new axes
                  xax_stack = [xax_stack; xax];
                  yax_tmp = yax_stack(end,:,:);
                  yax_tmp(1,:,channel) = yax;
                  yax_stack = cat(1, yax_stack, yax_tmp);

                  % increment ax_idx
                  ax_idx = ax_idx + 1;


                  plot_refresh = 1; % force replot
                  break;

                else
                  printf('\nInvalid key or button, try again...\n');
                endif


            case toascii('z') % zoom in ~10%, centered on {px}

                % remove any axes on stack after current axis
                xax_stack(ax_idx+1:end,:) = [];
                yax_stack(ax_idx+1:end,:,:) = [];

                % generate new x axis; keep all old y axes
                xax = round([px - 1/1.1 * diff(xax)/2, px + 1/1.1 * diff(xax)/2]);
                yax = yax;

                % append new axes
                xax_stack = [xax_stack; xax];
                yax_tmp = yax_stack(end,:,:);
                yax_tmp(1,:,channel) = yax;
                yax_stack = cat(1, yax_stack, yax_tmp);


                % increment ax_idx
                ax_idx = ax_idx + 1;


                plot_refresh = 1; % force replot
                break;


            case toascii('u') % zoom out ~10%, centered on {px}

                % remove any axes on stack after current axis
                xax_stack(ax_idx+1:end,:) = [];
                yax_stack(ax_idx+1:end,:,:) = [];

                % generate new x axis; keep old y axis
                xax = round([px - 1.1 * diff(xax)/2, px + 1.1 * diff(xax)/2]);
                yax = yax;

                % append new axes
                xax_stack = [xax_stack; xax];
                yax_tmp = yax_stack(end,:,:);
                yax_tmp(1,:,channel) = yax;
                yax_stack = cat(1, yax_stack, yax_tmp);


                % increment ax_idx
                ax_idx = ax_idx + 1;


                plot_refresh = 1; % force replot
                break;


            case toascii('p') % switch to previous axis on stack

                % update ax_idx
                ax_idx = max(1, ax_idx - 1);


                % move to previous axis on stack
                xax = xax_stack(ax_idx,:);
                yax = yax_stack(ax_idx,:,channel);

                plot_refresh = 1; % force replot
                break;


            case toascii('n') % switch to next axis on stack

                % update ax_idx
                ax_idx = min(size(xax_stack,1), ax_idx + 1);


                % move to next axis on stack
                xax = xax_stack(ax_idx,:);
                yax = yax_stack(ax_idx,:,channel);

                plot_refresh = 1; % force replot
                break;


            case toascii('h') % switch to original/default zoom

                ax_idx = 1;


                % move to first axis on stack
                xax = xax_stack(ax_idx,:);
                yax = yax_stack(ax_idx,:,channel);

                plot_refresh = 1; % force replot
                break;


            case toascii('i') % expand y axis by 50%


                % remove any axes on stack after current axis
                xax_stack(ax_idx+1:end,:) = [];
                yax_stack(ax_idx+1:end,:,:) = [];

                % generate new y axis; keep old x axis
                xax = xax_stack(end,:);
                #yax = round([mean(yax) - diff(yax)*(3/4), mean(yax) + diff(yax)*(3/4)]);
                yax = ([mean(yax) - diff(yax)*(3/4), mean(yax) + diff(yax)*(3/4)]);

                % append new axes
                xax_stack = [xax_stack; xax];
                yax_tmp = yax_stack(end,:,:);
                yax_tmp(1,:,channel) = yax;
                  yax_stack = cat(1, yax_stack, yax_tmp);


                % increment ax_idx
                ax_idx = ax_idx + 1;


                plot_refresh = 1; % force replot
                break;


            case toascii('k') % contract y axis by 33%


                % remove any axes on stack after current axis
                xax_stack(ax_idx+1:end,:) = [];
                yax_stack(ax_idx+1:end,:,:) = [];

                % generate new y axis; keep old x axis
                xax = xax_stack(end,:);
                #yax = round([mean(yax) - diff(yax)*(1/3), mean(yax) + diff(yax)*(1/3)]);
                yax = ([mean(yax) - diff(yax)*(1/3), mean(yax) + diff(yax)*(1/3)]);

                % append new axes
                xax_stack = [xax_stack; xax];
                yax_tmp = yax_stack(end,:,:);
                yax_tmp(1,:,channel) = yax;
                yax_stack = cat(1, yax_stack, yax_tmp);


                % increment ax_idx
                ax_idx = ax_idx + 1;


                plot_refresh = 1; % force replot
                break;


            case toascii('l') % shift y axis up by 25% of y range


                % remove any axes on stack after current axis
                xax_stack(ax_idx+1:end,:) = [];
                yax_stack(ax_idx+1:end,:,:) = [];

                % generate new y axis; keep old x axis
                xax = xax_stack(end,:);
                #yax = round([yax(1) + diff(yax)*(1/4), yax(2) + diff(yax)*(1/4)]);
                yax = ([yax(1) + diff(yax)*(1/4), yax(2) + diff(yax)*(1/4)]);

                % append new axes
                xax_stack = [xax_stack; xax];
                yax_tmp = yax_stack(end,:,:);
                yax_tmp(1,:,channel) = yax;
                  yax_stack = cat(1, yax_stack, yax_tmp);


                % increment ax_idx
                ax_idx = ax_idx + 1;


                plot_refresh = 1; % force replot
                break;


            case toascii('j') % shift y axis down by 25% of y range


                % remove any axes on stack after current axis
                xax_stack(ax_idx+1:end,:) = [];
                yax_stack(ax_idx+1:end,:,:) = [];

                % generate new y axis; keep old x axis
                xax = xax_stack(end,:);
                #yax = round([yax(1) - diff(yax)*(1/4), yax(2) - diff(yax)*(1/4)]);
                yax = ([yax(1) - diff(yax)*(1/4), yax(2) - diff(yax)*(1/4)]);

                % append new axes
                xax_stack = [xax_stack; xax];
                yax_tmp = yax_stack(end,:,:);
                yax_tmp(1,:,channel) = yax;
                yax_stack = cat(1, yax_stack, yax_tmp);


                % increment ax_idx
                ax_idx = ax_idx + 1;


                plot_refresh = 1; % force replot
                break;


            case {toascii('1') toascii('2') toascii('3') toascii('4') ...
                  toascii('5') toascii('6') toascii('7') toascii('8') ...
                  toascii('9') toascii('0') }

                channel_orig = channel;

                if (mod(channel,10) == str2num(char(button)))
                  channel = channel+10;
                else
                  channel = str2num(char(button));
                endif

                if (channel > size(ts_in,3))
                  channel = channel_orig;
                  printf("\nInvalid channel specified, try again...\n");
                  fflush(stdout);
                  continue;
                endif

                xax = xax;
                yax = yax_stack(ax_idx,:,channel);

                plot_refresh = 1; % force replot
                break;


            case toascii('d') % enter debug mode..in Linux, this could be done using
                                      % ctrl-c and setting debug_on_interrupt, but ctrl-c is
                                      % not reliable in Windows
                keyboard;

              otherwise

                  printf('\nInvalid key or button, try again...\n');

            end % endswitch button
          end % endif isempty(button)
          fflush(stdout);


      end % while (2)


      % if the break command is used in a switch block, it only breaks the
      % switch block, not this loop
      if isempty(button) break; endif


      % remove duplicate coordinates inside the loop...this is not as
      % efficient, but it is necessary for plotting the updated image;
      % also, make sure ts_idx actually exists
      if (ismatrix(ts_idx))
         ts_idx = unique(ts_idx);
      endif


   end % endwhile button ~= toascii('q')


   printf('\n');
   fflush(stdout);


   % make sure ts_idx is row vector, NOT column vectors...why again? -EJR 3/2012
   ts_idx = reshape(ts_idx, 1, numel(ts_idx));


end
