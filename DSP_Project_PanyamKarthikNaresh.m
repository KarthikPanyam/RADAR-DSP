clc;
clear all;
close all;

%% Target 1 Distance in m, speed in kmph, angle in degrees
target1d =100;
target1speed = -50*5/18;
target1angle = 45;


%% Target 2 Distance in m, speed in kmph, angle in degrees
target2d = 200;
target2speed = 75*5/18;
target2angle = -45;

%% Parameters 
fc = 10*10^9;                     % Carrier frequency
c = 3*10^8;                       % Speed of light const.
lambdac = c/fc;                   % Carrier wavelength

Vmax = 300*5/18;                   % Max velocity that can be detected
Vres = 1.4 *5/18;                 % Velocity resolution
Dmax = 300;                       % Max distance 
Dres = 0.3;                       % Distance resolution



Mr = 16;                           % Number of antenna
Dant = lambdac/2;                 % Distance between the antenna 

%% Chirp Parametrs
bchirp = c/(2*Dres);                                        % BW of the chirp;
nb_chirp = (2*Vmax)/Vres;                                   % Number of chirps
Tchirp = c/(fc*2*nb_chirp*Vres);                            % time period of each chirp
Toverall = nb_chirp*Tchirp;                                 % Overall tx time
kchirp = bchirp/Tchirp;                                     % modulation index
fbmax = abs(2*Vmax/lambdac) + kchirp*2*Dmax/c;              % Max beat frequency

%% Samples
fs = 2*fbmax;                                   % Sample frequency 
Ts = 1/fs;                                      % Sample Period
samples = Dmax/Dres;

ybeat = zeros(floor(nb_chirp),samples,Mr);



%% generating synthetic signal

for  k = 1:Mr
    for n = 1:samples
        for l = 1:nb_chirp

            ybeat(l,n,k) = ybeat(l,n,k)+exp(2*1*j*pi*(2*bchirp*target1d*(n)*Ts/(Tchirp*c)+2*target1speed*(l-1)*Tchirp/lambdac+Dant*sind(target1angle)*(k-1)/lambdac)) + exp(2*1j*pi*(2*bchirp*target2d*(n)*Ts/(Tchirp*c)+2*target2speed*(l-1)*Tchirp/lambdac+Dant*sind(target2angle)*(k-1)/lambdac));
            
        end
    end
end
%% Calculationg the ffts 

y1 = ybeat(:,:,1);
y1fft = fft2(y1);
y1fft_mag = abs(fftshift(y1fft));

%% Finding the peaks
largest1 = max(y1fft_mag(:));
[row1,col1] = find(y1fft_mag == largest1);

temp = max(y1fft_mag);
temp2 = sort(temp,"descend");

largest2 = temp2(1,2);
[row2,col2] = find(y1fft_mag == largest2);

%% Plotting the Range-Doppler plot

xaxis = linspace(-0.5,0.5,samples);
xscalingfactor = Tchirp*c*fs/(2*bchirp);
xscaled = floor(xaxis*xscalingfactor);          % X-Axis---------> Range

yaxis = linspace(-0.5,0.5,nb_chirp);
yscalingfactor = 3.6*lambdac/(2*Tchirp);
yscaled = yaxis*yscalingfactor;                 % Y-Axis---------> Doppler


% Setting the Range-Doppler plot
figure(1);
imagesc([-0.5*xscalingfactor,0.5*xscalingfactor],[-0.5*yscalingfactor,0.5*yscalingfactor],y1fft_mag);
title(" Range-Doppler Plot")
set(gca,'Ydir','normal');
colorbar;
colormap jet;
hold on
grid on;

%plot(xscaled,y1fft_mag(row1,col1)); % Range doppler
xlabel("Range");
ylabel("Doppler");

hold off;

%% 1D Slices
figure(2)

subplot(2,1,1);
plot(yscaled,y1fft_mag(:,col1)); % Velocity
title("Doppler")
xlabel("Doppler (Kmph)");
ylabel(" ");
hold on;
plot(yscaled,y1fft_mag(:,col2));
hold off;

subplot(2,1,2)
plot(xscaled,y1fft_mag(row1,:));% Range
title("Range")
xlabel("Range (Km)")
ylabel(" ")
hold on;
plot(xscaled,y1fft_mag(row2,:));
hold off;

%% Performing the 3D FFT

thirdfft = fftshift(fftn(ybeat));
new1 = abs(thirdfft);

magnitude = new1(:,:,2);



largest3 = max(magnitude(:));
[row3,col3] = find(magnitude == largest3);

temp3 = max(magnitude);
temp4 = sort(temp3,"descend");

largest4 = temp4(1,2);
[row4,col4] = find(magnitude == largest4);
t1Mr = new1(row3,col3,:);
t1Mrreshape = reshape(t1Mr,[1,Mr]);

t2Mr = new1(row4,col4,:);
t2Mrreshape = reshape(t2Mr,[1,Mr]);

poi = linspace(-01,01,Mr);

azi = asind(poi);
figure(3);
subplot(2,1,1);
plot(azi,t1Mrreshape);
title("azimuth");
hold on;
plot(azi,t2Mrreshape);
hold off;
subplot(2,1,2)
plot(poi,t1Mrreshape);
title("Sin(azimuth) ");
hold on;
plot(poi,t2Mrreshape);
hold off;

%% Range-Azimuth plane (Considered the first chirp):

rngazi = new1(1,:,:);
rngazi = reshape(rngazi,[samples,Mr]);
figure(4);
imagesc([-0.5*xscalingfactor,0.5*xscalingfactor],[-1,1],rngazi);
title(" Range-Azimuth Plot")
set(gca,'Ydir','normal');
colorbar;
colormap jet;
hold on
grid on;


xlabel("Range");
ylabel("Sine(azimuth)");

hold off;

%% 1D slices

figure(5)
subplot(2,1,1)
plot(xscaled,y1fft_mag(row1,:));            % Range
title("Range")
xlabel("Range (Km)")
ylabel(" ")
hold on;
plot(xscaled,y1fft_mag(row2,:));
hold off;
subplot(2,1,2)

plot(poi,t1Mrreshape);                         % Sin(azimuth)
title("Sin(azimuth) ");
hold on;
plot(poi,t2Mrreshape);
hold off;

%% Bird's eye view
figure(6)
xyz = zeros(Mr*Mr/2,3);
xyz(:,3)=xyz(:,3)/max(xyz(:,3));  % here "xyz" has size (FFT_azimuth_length*(FFT_range_length/2),3), with the size-3 dimension used for x_coordinates, y_coordinates, and magnitude information 
colors=repmat([1,1,1],length(xyz(:,3)),1); % define a custom color scheme for scatter plot
for i=1:1:length(xyz(:,3))
    colors(i,:)=colors(i,:)*(1.0-xyz(i,3));
end
s=5;
scatter3(xyz(:,1),xyz(:,2),xyz(:,3),s,colors,'o','filled')
view(0,90) 
