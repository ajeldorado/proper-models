function [pt,rvec,qvec] = polarTransform(input_image, center_vec, rmax, numRadPts, numAngles,method)
%[pt,rvec,qvec] = polarTransform(input_image, center_vec, rmax, numRadPts, numAngles,method)
%   Computes the polar transform of input_image.
%
%   Inputs:
%       input_image: 2D array
%       center_vec: coordinates of origin in the image
%       rmax: max radius to compute polar transform 
%       numRadPts: number of points in the radial direction
%       numAngles: number of polar angles
%       method: interpolation method. See interp2 docs
%   Outputs:
%       pt: 2D array with polar transform
%       rvec: vector of radial positions 
%       qvec: vector of polar angles 

    
    rvec = linspace(0,rmax,numRadPts);% make array of radial points
    qvec = linspace(0,2*pi-2*pi/numAngles,numAngles);% make array of azimuthal points
    radialComp = repmat(rvec', 1, numAngles);% 2D array of radial points 
    angleComp = repmat(qvec, numRadPts, 1);% 2D array of azimuthal points 
    [xComp,yComp] = pol2cart(angleComp,radialComp);% Transform desired polar coords into cartesian coords
    
    % Make image coordinates 
    [rows,cols] = size(input_image);
    
    xvals = (1:rows) - center_vec(1);
    yvals = (1:cols) - center_vec(2);
    
    % compute polar transform
    pt = interp2(xvals,yvals,double(input_image),xComp,yComp,method);

end

