clear all;

R = load('Red.dat');
G = load('Green.dat');
B = load('Blue.dat');

[ny,nx] = size(R);

FF = zeros(ny,nx,3);
FF(:,:,1) = R;
FF(:,:,2) = G;
FF(:,:,3) = B;

figure('color','k')
image(FF)
axis image
