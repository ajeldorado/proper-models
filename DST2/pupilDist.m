function [z2] = pupilDist(f1,f2,z1)
%This function gives the distance a pupil plane is behind a pair of lenses
%(with focal lengths f1 and f2) that are each separated by their focal 
%lengths given an object pupil a distance z1 in front of lens (denoted f1)


z2 = f2*(f1*(f1+f2)-f2*z1)/f1^2;

end

