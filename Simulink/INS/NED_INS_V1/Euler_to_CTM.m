function C = Euler_to_CTM(eul)
%Euler_to_CTM - Converts a set of Euler angles to the corresponding
%coordinate transformation matrix
%
% Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition.
%
% This function created 1/4/2012 by Paul Groves
%
% Inputs:
%   eul     Euler angles describing rotation from beta to alpha in the 
%           order roll, pitch, yaw(rad)
%
% Outputs:
%   C       coordinate transformation matrix describing transformation from
%           beta to alpha

% Copyright 2012, Paul Groves
% License: BSD; see license.txt for details

% Begins

% Precalculate sines and cosines of the Euler angles
sin_phi = sin(eul(1));
cos_phi = cos(eul(1));
sin_theta = sin(eul(2));
cos_theta = cos(eul(2));
sin_psi = sin(eul(3));
cos_psi = cos(eul(3));

% Calculate coordinate transformation matrix using (2.22)
C11 = cos_theta * cos_psi;
C12 = cos_theta * sin_psi;
C13 = -sin_theta;
C21 = -cos_phi * sin_psi + sin_phi * sin_theta * cos_psi;
C22 = cos_phi * cos_psi + sin_phi * sin_theta * sin_psi;
C23 = sin_phi * cos_theta;
C31 = sin_phi * sin_psi + cos_phi * sin_theta * cos_psi;
C32 = -sin_phi * cos_psi + cos_phi * sin_theta * sin_psi;
C33 = cos_phi * cos_theta;

C = [C11,C12,C13;C21,C22,C23;C31,C32,C33];

% Ends