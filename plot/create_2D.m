function [out] = create_2D(data,xdim,ydim,zdim,zlow,zhigh,show,xpose,gridsize,minv,maxv,cb)
%function [out] = create_2D(data,xdim,ydim,zdim,zlow,zhigh,show,xpose,gridsize,minv,maxv)

if (~exist('xdim','var')),
   xdim = size(data,1);
   ydim = size(data,2);
   zdim = size(data,3);
   zlow = 1;
   zhigh = zdim;
end;

data = squeeze(data);
if (~exist('xdim','var')),
   if length(size(data)) == 3,
      xdim = size(data,1);ydim = size(data,3);zdim = size(data,3);
      zlow = 1;zhigh = zdim;
   else,
      disp('Error:must be a 3D vector or specify dimensions');
   end;
end;
if (~exist('zlow','var')),
   zlow = 1;
   zhigh = zdim;
end;
if (~exist('cb','var')),
   cb = 1;%colorbar?
end;


if (~exist('show','var')),
   show = 1;
end;
if (~exist('gridsize','var')),
   gridsize = -999;
end;
if (~exist('xpose','var')),
   xpose = 1;
end;

if (~exist('minv','var')),
   minv = -maxN(abs(data));
   %minv=2;
end;

if (~exist('maxv','var')),
   maxv = maxN(abs(data));
   %maxv=2;
end;

data = reshape(data,xdim,ydim,zdim);
data = data(:,:,1:(end-0));
data = reshape(data,xdim,ydim,1,zdim);
a1b = zeros(size(data,2),size(data,1),size(data,3),size(data,4));
%a1b = zeros(size(data,1),size(data,2),size(data,3),size(data,4));

for j = 1:(zdim),
   a1b(:,:,1,end-j+1) = (data(:,end:-1:1,1,j)');
   %a1b(:,:,1,end-j+1) = data(:,end:-1:1,1,j);
end;

if (xpose),
   out = (get_montage(a1b(:,:,:,zlow:zhigh),gridsize)');
else,
   out = (get_montage(a1b(:,:,:,zlow:zhigh),gridsize));
end;
if (show),
   %figure;
   colortables
   %scl = maxN(abs(out));
   % imagesc(out',[minv maxv]);axis off;axis image;colormap(coldhot);
   imagesc(out',[minv maxv]);axis image; colormap(coldhot);
   if (cb),
      colorbar;
   end;
%   imagesc(out,[0 maxN(out)]);%.*(abs(a1)>3*std(flatrow(a1))));
end;