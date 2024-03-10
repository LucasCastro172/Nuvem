clear all;

ny = 1152;
nx = 2048;
FF = zeros(ny,nx,3);
FF2= zeros(ny/2,nx/2,3);

FF(:,:,1) = load('Red.dat');
FF(:,:,2) = load('Green.dat');
FF(:,:,3) = load('Blue.dat');

FF2(:,:,1) = load('Red2.dat');
FF2(:,:,2) = load('Green2.dat');
FF2(:,:,3) = load('Blue2.dat');

figure('color','k')
image(FF)
axis image

figure('color','k')
image(FF2)
axis image

