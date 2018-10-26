function [dcvldt] = VanLeer1D(c,v,dz,dt)

% 
% MATLAB Implementation of 1D FVM Van Leer I
%  dcdt for interior cells (iz=2:nz+1)
%

% Internal slopes
slopes = zeros(size(c));
s_mask = ...
    find([0 ((c(2:end-1) > c(3:end) & c(2:end-1) <= c(1:end-2)) | ...
    (c(2:end-1) < c(3:end) & c(2:end-1) >= c(1:end-2))) 0] > 0);
slopes(s_mask) = 2 .* ...
    (c(s_mask+1) - c(s_mask)) .* (c(s_mask) - c(s_mask-1)) ./ ...
    (c(s_mask+1) - c(s_mask-1));

% Time-integrated fluxes
fz_s = dt * v(1:end-1) / dz;
fz_n = dt * v(2:end) / dz;

% Van Leer face-centered values
cvl_n = ...
    (fz_n > 0) .* (c(2:end-1) + 0.5 .* (1 - fz_n) .* slopes(2:end-1)) + ...
    (fz_n <= 0) .* (c(3:end) - 0.5 .* (1 + fz_n) .* slopes(3:end));
cvl_s = ...
    (fz_s > 0) .* (c(1:end-2) + 0.5 .* (1 - fz_s) .* slopes(1:end-2)) + ...
    (fz_s <= 0) .* (c(2:end-1) - 0.5 .* (1 + fz_s) .* slopes(2:end-1));

% Assemble Van Leer flux terms
dcvldt = (cvl_s .* v(1:end-1) - cvl_n .* v(2:end)) / dz; 

end
