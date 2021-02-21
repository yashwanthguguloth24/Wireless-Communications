clear all;
clc;

% Date 2/7/2021
% Yashwanth


% Jakes Channel model
% It assumes isotropic scattering
% R = R_I + j*R_Q;
% R_I -- inphase component
% R_Q -- quadrature phase component

% Params
N = 30;
M = 0.5*((N/2)-1);
Fs = 1000; % 1000 samples per second
Ts = 1/Fs; % Sampling time
fD = [1,10,100];
Duration = 5;

% Initiazing the matrices
R_I = zeros(M,(Duration*Fs)+1);
R_Q = zeros(M,(Duration*Fs)+1);
alpha_n = zeros(1,M);
beta_n = zeros(1,M);
f_n = zeros(1,M);
R_I_alpha = zeros(1,(Duration*Fs)+1);

% Taking (Duration*Fs)+1 samples for each i for i = 1 to M
for i = 1:1:3
    for n = 1:1:M
        k = 1;
        for t = 0:Ts:Duration
            alpha(n) = ((2*pi*n)/N);
            f_n(n) = fD(i)*cos(alpha(n));
            beta_n(n) = (pi*n)/M;
            % Inphase component
            R_I(n,k) = 2*cos(beta_n(n))*cos(2*pi*f_n(n)*t);
            % Quadrature phase component 
            R_Q(n,k) = 2*sin(beta_n(n))*cos(2*pi*f_n(n)*t);
            % R_I_alpha = cos(2*pi*fd*t)*sqrt(2)
            % R_Q_alpha = 0
            R_I_alpha(k) = 2*cos(2*pi*fD(i)*t)/sqrt(2);
            k = k + 1;
        end
    end
    R_I_sum = sum(R_I)+R_I_alpha;
    R_Q_sum = sum(R_Q);
    R_sum = sqrt(R_I_sum.^2+R_Q_sum.^2);
    avg = sum(R_sum)/((Duration*Fs)+1);
    
    % Plot for different fD
    subplot(3,1,i);
    time = 0:Ts:Duration;
    plot(time,(10*log10(R_sum)-10*log10(avg)));
    ti = ['fD = ' int2str(fD(i)) 'Hz'];
    title(ti);
    xlabel('Time(s)')
    ylabel('Recieved Envelope(dB)');

        
end      
