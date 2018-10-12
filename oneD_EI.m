%%% 1D simulation for Edge Illumination X-Ray Phase Contrast Imaging %%%

clear;
tic;
%%%   INPUT   %%%
% GEOMETRY %
z1 = 1.6; % source to sample distance
z2 = 0.4; % sample to detector distance
M = (z1+z2)/z1; % magnification

% SAMPLING%
dx = .05e-6; % sampling step along x (m)

% SOURCE %
ener = 21e3; % photons energy (eV)
source_width = 70e-6; % source FWHM along x(m)
la = 1.24e-6./ener; % wavelength (m)
sigma_s_ini = source_width/(2*sqrt(2*log(2))); % source sigma along x(m)

% DETECTOR MASK %
D_num_ap = 21; % D_num_ap <= pixel_num
D_step_ini = 50e-6; % period (m)
D_aper_ini = 20e-6; % apertures (m)
D_beta = 1; % absorbing septa beta
D_delta = 0; % absorbing septa delta
D_thick = 1; % absorbing septa thickness (m)
D_off_ini = 0; % detector mask displacement (m)

% SAMPLE MASK %
S_num_ap = D_num_ap; % S_num_ap <= pixel_num
S_step = D_step_ini/M; % period (m)
S_aper = 12e-6; % aperture size (m)
S_beta = 1; % absorbing septa beta
S_delta = 0; % absorbing septa delta
S_thick = 1; % absorbing septa thickness (m)
S_off = -D_aper_ini/(2*M); % sample mask displacement (m)

% DETECTOR %
pixel_size_ini = D_step_ini; % pixel size along x (m)
pixel_num = D_num_ap; % pixel number
psf_fwhm = 0e-6;

% SAMPLE %
O_beta = 1e-9; % sample beta
O_delta = 1e-6; % sample delta
O_rad = 150e-6; % sample radius (m)
H_beta = 0; % host beta 
H_delta = 0; % host delta
H_thick = 0; % host thickness (m)
num_dit = 8; % number of dithering positions
dit_step = D_step_ini/(M*num_dit); % dithering step (m)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%   PARALLEL GEOMETRY   %%%

z_eff = z2/M;
sigma_s = sigma_s_ini*(M-1)/M;
sigma_psf = psf_fwhm/(M*2*sqrt(2*log(2)));
pixel_size = pixel_size_ini/M;
D_step = D_step_ini/M;
D_aper = D_aper_ini/M;
D_off = D_off_ini/M;

%%%   PARAMETERS   %%%

Lx = (S_num_ap+1)*S_step + abs(S_off) + pixel_size*(pixel_num+1) + 9*(sigma_s+sigma_psf);
if ( z_eff*la/dx ) < ( Lx )%z_eff*la < Lx*dx
    Lx = (z_eff*la/(2*dx) + Lx/2);
    Lx = max(Lx, S_num_ap*S_step + abs(S_off) + z_eff*la/dx);
    Nx = 2*round(Lx/(2*dx)); Lx = Nx*dx;
    x =(-(Nx/2):(Nx/2-1))*dx;
    dkx = 1/(Nx*dx);
    kx = dkx*((-Nx/2):(Nx/2-1));
    Hf = fftshift(exp(-1i*pi*la*((kx).^2)*z_eff));
    fprintf('   Angular spectrum\n');
else
    Nx = 2*round(Lx/(2*dx)); Lx = Nx*dx;
    x =(-(Nx/2):(Nx/2-1))*dx;
    dkx = 1/(Nx*dx);
    kx = dkx*((-Nx/2):(Nx/2-1));
    Hf = dx*fft(fftshift(exp(1i*pi*x.^2./(z_eff*la))))./sqrt(1i*z_eff*la);
    fprintf('   Other\n');
end


[~,n_zero] = min(abs(x));
fprintf('2^%i sampling points\n',round(log2(Nx)));

if round(log2(Nx))>23
    ris = input('Number of sampling points higher than 2^23. Do you want to continue? (Y/N) ','s');
    if (ris == 'n')||(ris == 'N')
        error('Execution stopped');
    end
end

%%%   SIMULATION   %%%
    
ft_Pr = ifftshift(exp(-2*pi*pi*sigma_s^2*kx.^2));

ft_Psf = ifftshift(exp(-2*pi*pi*sigma_psf^2*kx.^2));

Pixel = heaviside(pixel_size/2-abs(x));
Pixel = ifftshift(Pixel)/sum(Pixel(:));
ft_Pixel = fft(Pixel);

Dtra = exp(-1i*2*pi*D_delta*D_thick/la-2*pi*D_beta*D_thick/la);
Stra = exp(-1i*2*pi*S_delta*S_thick/la-2*pi*S_beta*S_thick/la);
Ds = zeros(num_dit,pixel_num);
Dr = zeros(num_dit,pixel_num);
pixel_cen = floor((-(pixel_num-1)/2):((pixel_num-1)/2))*pixel_size;
pixel_ind = round(pixel_cen/dx) + n_zero;
        
Dmask = mask(D_num_ap, D_step, D_aper, x, D_off, Dtra );
Smask = mask(S_num_ap, S_step, S_aper, x, S_off, Stra );
for idx = 1:num_dit
    dit_cen = (-(num_dit-1)/2+(idx-1))*dit_step + S_off;
    sample = wire( O_delta, O_beta, O_rad, H_delta, H_beta, H_thick, x, dit_cen, la);
    Is = abs(ifft(fft(ifftshift(Smask.*sample)).*Hf)).^2;
    Is = fftshift(ifft(ft_Pr.*fft(Is))).*(abs(Dmask).^2);
    Is = abs(fftshift(ifft(ft_Psf.*ft_Pixel.*fft(ifftshift(Is)))));
    Ds(idx,:) = Is(pixel_ind);
end
        
Is = abs(ifft(fft(ifftshift(Smask)).*Hf)).^2;
Is = fftshift(ifft(ft_Pr.*fft(Is))).*(abs(Dmask).^2);
Is = fftshift(ifft(ft_Psf.*ft_Pixel.*fft(ifftshift(Is))));
for idx = 1:num_dit
    Dr(idx,:) = Is(pixel_ind);
end
IN = Ds./Dr;
toc;

figure;plot((1:(pixel_num*num_dit))*dit_step*1e6,IN(:));
xlabel('Position (\mum)');
ylabel('Normalized intensity');