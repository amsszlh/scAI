function colors = generateColors(n_colors,format)
% Example:
%   c = generateColors(26);
%   figure
%   image(reshape(c,[1 size(c)]))
if ~exist('format','var') | isempty(format)
    format = 'RGB';  % (format,'HEX')
end
colorSpace = [228,26,28 % red
    55,126,184 % blue
    77,175,74 % green
    152,78,163 % purple
    242 148 3 % orange
    247,129,191 % pink
    188 157 204; % light purple
    166,86,40 % brown
    84 176 228 % light blue
    34 47 117 % dark blue
    27,158,119 % dark green
    178,223,138 % light green
    227 190 0 % yellow%
    %255,255,51 % yellow
    251,154,153 % light red
    231,41,138 % % bright pink
    145 2 65 % dark red
    0 205 209; % C11
    %     255 0 255 % bright pink 
    166,206,227
    206 18 97
    94,79,162 % dark purple
    140 167 123 % light green
    0,68,27 % dark green
    222 220 0 % light yellow
    179,222,105 % light green
    141,211,199 %
    153,153,153 % grey
    ]/255;
if n_colors <= 26
    colors = colorSpace(1:n_colors,:);
end
if strcmpi(format,'HEX')
    colors = rgb2hex(colors);
end


function [ hex ] = rgb2hex(rgb)
% rgb2hex converts rgb color values to hex color format.
%
% This function assumes rgb values are in [r g b] format on the 0 to 1
% scale.  If, however, any value r, g, or b exceed 1, the function assumes
% [r g b] are scaled between 0 and 255.
%
% * * * * * * * * * * * * * * * * * * * *
% SYNTAX:
% hex = rgb2hex(rgb) returns the hexadecimal color value of the n x 3 rgb
%                    values. rgb can be an array.
%
% * * * * * * * * * * * * * * * * * * * *
% EXAMPLES:
%
% myhexvalue = rgb2hex([0 1 0])
%    = #00FF00
%
% myhexvalue = rgb2hex([0 255 0])
%    = #00FF00
%
% myrgbvalues = [.2 .3 .4;
%                .5 .6 .7;
%                .8 .6 .2;
%                .2 .2 .9];
% myhexvalues = rgb2hex(myrgbvalues)
%    = #334D66
%      #8099B3
%      #CC9933
%      #3333E6
%
% * * * * * * * * * * * * * * * * * * * *
% Chad A. Greene, April 2014
%
% Updated August 2014: Functionality remains exactly the same, but it's a
% little more efficient and more robust. Thanks to Stephen Cobeldick for
% his suggestions.
%
% * * * * * * * * * * * * * * * * * * * *
% See also hex2rgb, dec2hex, hex2num, and ColorSpec.

%% Check inputs:

assert(nargin==1,'This function requires an RGB input.')
assert(isnumeric(rgb)==1,'Function input must be numeric.')

sizergb = size(rgb);
assert(sizergb(2)==3,'rgb value must have three components in the form [r g b].')
assert(max(rgb(:))<=255& min(rgb(:))>=0,'rgb values must be on a scale of 0 to 1 or 0 to 255')

%% If no value in RGB exceeds unity, scale from 0 to 255:
if max(rgb(:))<=1
    rgb = round(rgb*255);
else
    rgb = round(rgb);
end

%% Convert (Thanks to Stephen Cobeldick for this clever, efficient solution):

hex(:,2:7) = reshape(sprintf('%02X',rgb.'),6,[]).';
hex(:,1) = '#';



