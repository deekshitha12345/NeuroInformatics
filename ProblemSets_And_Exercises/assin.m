% Kernel size (odd for symmetry)
kSize = 9;  

%% 1. Inverted U kernel (parabolic shape)
x = linspace(-1, 1, kSize);             % normalized x-axis
invU = -(x.^2) + 1;                     % parabola opening downward
invU = invU / sum(abs(invU));           % normalize

%% 2. Decay kernel (exponential decay)
lambda = 2;                             % decay rate
x = linspace(0, 2, kSize);              % only positive side
decay = exp(-lambda * x);               % exponential decay
decay = decay / sum(decay);             % normalize

%% Display kernels
disp('Inverted U kernel:');
disp(invU);

disp('Decay kernel:');
disp(decay);

%% Example: Apply convolution to a signal
signal = [zeros(1,10), ones(1,10), zeros(1,10)];
convU = conv(signal, invU, 'same');
convDecay = conv(signal, decay, 'same');

figure;
subplot(3,1,1);
stem(signal, 'filled'); title('Original Signal');

subplot(3,1,2);
plot(convU, 'r'); title('Convolution with Inverted U kernel');

subplot(3,1,3);
plot(convDecay, 'b'); title('Convolution with Decay kernel');
