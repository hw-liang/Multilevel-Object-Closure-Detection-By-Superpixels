% [max_pb,phase] = pbWrapper(I)
% [max_pb,phase] = pbWrapper(I,pb_timing)
%
% Wrap Martin-Fowlkes pb code to do phase calcuations as well.
function [max_pb,phase] = pbWrapper(I,varargin)
if nargin >=2
  pb_timing = varargin{1};
else
  pb_timing=0;
end

[r,c,k] = size(I);

% any missing parameter is substituted by a default value
par = [8,1,21,3];
par(end+1:4)=0;
j = (par>0);
% make the filter size an odd number so that the responses are not skewed
if mod(par(3),2)==0, par(3) = par(3)+1; end
j = num2cell(par);
[n_filter,n_scale,winsz,enlong] = deal(j{:});
n = ceil(winsz/2);

if c < n+1 | r < n+1
  oe = ones(r,c);
  zcrs = zeros(r,c);
  return;
end

if pb_timing, st=clock; end
if (nargin >= 3)
    load(varargin{2});
    max_pb = image_data.pb;
    theta = image_data.pb_theta;
else
    [max_pb,theta] = pbCGTG(I);
end
if pb_timing, fprintf('pb took %.2f minutes\n',etime(clock,st)/60); end

% filter to get phase info
FBo = make_filterbank_odd2(8,1,21,3);

% filter max_pb with FBo, this gets us 2nd derivative at each point/orientation
f = [fliplr(max_pb(:,2:n+1)), max_pb, fliplr(max_pb(:,c-n:c-1))];
f = [flipud(f(2:n+1,:)); f; flipud(f(r-n:r-1,:))];
Fpbo = fft_filt_2(f,FBo,1); 
Fpbo = Fpbo(n+[1:r],n+[1:c],:);

% select orientation using theta map
P_map = zeros(r,c);
for t_i=1:n_filter
  P_map = P_map + (theta==(t_i-1)*pi/n_filter).*Fpbo(:,:,t_i);
end

% threshold to get phase
phase = (P_map >= 0) - (P_map<0);
