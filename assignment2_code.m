% eee431 hw2 - Emir Bayer

clear; clc; close all;
rng(88);


%% part 1

% part 1a
T =1;
Ts =T/20;
t =0:Ts:T-Ts;


p = t .* (T - t);
Ep= sum(p.^2)*Ts;

phi1 = p / sqrt(Ep);

figure;
subplot(2,1,1);
plot(t, p, 'LineWidth', 1.5); hold on;
plot(t, -p, 'r--', 'LineWidth', 1.5);
xlabel('time (s)'); ylabel('amplitude');
legend('p(t)', '-p(t)');
title('part 1(a): signal waveforms'); grid on;

subplot(2,1,2);
plot(t, phi1, 'k', 'LineWidth', 1.5);

xlabel('time (s)'); ylabel('amplitude');
title('part 1(a): basis function'); grid on;



% part 1 b-e
n_sim = 1e6;
gam_db= linspace(-10, 9, 10);

Pe_theo=zeros(size(gam_db));

Pe_sim=zeros(size(gam_db));

fprintf('part 1(d) equal prob\n');

for k=1:length(gam_db)

    snr_lin= 10^(gam_db(k)/10);
    Pe_theo(k) =qfunc(sqrt(2*snr_lin));
    
    N0 = Ep/ snr_lin;
   
    bits = rand(n_sim, 1)> 0.5;
    
    tx_val = zeros(n_sim, 1);

    tx_val(bits==0) =sqrt(Ep);
    tx_val(bits==1) =-sqrt(Ep);
    
    noise = sqrt(N0/2)*randn(n_sim, 1);

    y = tx_val + noise;

    
    decisions= y<0;
    
    errs = sum(bits ~=decisions);
    Pe_sim(k) = errs /n_sim;
    
    fprintf('SNR=%.2f dB, Pe_theo=%.5e, Pe_sim=%.5e\n', gam_db(k),Pe_theo(k), Pe_sim(k));
end

figure;
plot(gam_db, Pe_sim, '-o', 'LineWidth',1.2);

grid on; xlabel('\gamma_s (dB)');ylabel('Pe');
title('part 1(d): estimated SEP');

figure;

semilogy(gam_db,Pe_theo, 'k-', 'LineWidth',1.5); hold on;
semilogy(gam_db, Pe_sim, 'ro','LineWidth', 1.2);

grid on; xlabel('\gamma_s (dB)'); ylabel('Pe');
legend('Theory','Simulated');

title('part 1(e): SEP vs SNR');


% part 1 f-h
p1= 0.25; p0= 0.75;
Pe_opt =zeros(size(gam_db));

Pe_sub = zeros(size(gam_db));



fprintf('\npart 1(g)-(h) unequal priors\n');

for k=1:length(gam_db)


    snr_lin= 10^(gam_db(k)/10);
    N0 = Ep/ snr_lin;
    
    bits =rand(n_sim, 1) < p1;
    
    tx_val =zeros(n_sim, 1);
    tx_val(bits==0)= sqrt(Ep);


    tx_val(bits==1)= -sqrt(Ep);
    
    y = tx_val + sqrt(N0/2)* randn(n_sim, 1);
    
    eta =(N0/ (4*sqrt(Ep)))* log(p1/p0);
    dec_opt= y < eta;
    Pe_opt(k)= sum(bits ~= dec_opt) / n_sim;
    
    dec_sub = y< 0;
    Pe_sub(k) =sum(bits ~= dec_sub)/ n_sim;

    
    fprintf('SNR=%.2f dB, Pe_opt=%.5e, Pe_sub=%.5e\n', gam_db(k),Pe_opt(k),Pe_sub(k));
end

figure;

plot(gam_db, Pe_opt,'-o', 'LineWidth',1.2);
grid on; xlabel('\gamma_s (dB)'); ylabel('Pe');
title('part 1(g): MAP estimated SEP');

figure;

semilogy(gam_db, Pe_opt,'b-o', 'LineWidth', 1.2); hold on;

semilogy(gam_db, Pe_sub,'r--s','LineWidth', 1.2);
grid on; legend('MAP', 'ML');
xlabel('\gamma_s (dB)'); ylabel('Pe');
title('part 1(h): unequal priors comparison');


%% part 2

% part 2 a-c
gam_db2 =linspace(-4, 18, 12); 

n_sim2 =1e6;

A1= sqrt(4/7);

syms1 =[A1 A1; A1 -A1; -A1 -A1; -A1 0]'; 

Pe_sim1= zeros(size(gam_db2));

Pe_bound1= zeros(size(gam_db2));

fprintf('\npart 2 constellation 1\n');

for k=1:length(gam_db2)
    snr_lin = 10^(gam_db2(k)/10);
    N0 = 1/snr_lin;
    
    u_sum =0;
    for i=1:4

        for j=1:4
            if i~= j
                d =norm(syms1(:,i) - syms1( :,j));
                u_sum = u_sum+ qfunc(d /sqrt(2*N0));
            end

        end
    end
    Pe_bound1(k) = u_sum/ 4;
    
    tx_inds = randi([1 4],n_sim2, 1);

    tx_vecs = syms1(:,tx_inds);

    rx_vecs = tx_vecs+ sqrt(N0/2)*randn(2, n_sim2);
    


    errs = 0;
    for m=1:n_sim2
        [~, dec] =min(sum((syms1 -rx_vecs(:,m)).^2,1));
        if dec ~=tx_inds(m),errs=errs+1; end
    end
    Pe_sim1(k) =errs/n_sim2;
    
    fprintf('SNR=%.2f, Bound=%.4e, Sim=%.4e\n', gam_db2(k),Pe_bound1(k),Pe_sim1(k));
end

figure;
plot(gam_db2,Pe_sim1,'-o'); grid on;

xlabel('\gamma_s (dB)'); ylabel('Pe');
title('part 2(b): const 1 sim only');

figure;
semilogy(gam_db2, Pe_bound1, 'k-', 'LineWidth', 1.5); hold on;

semilogy(gam_db2,Pe_sim1, 'ro');
grid on; legend('Bound','Simulated');
xlabel('\gamma_s (dB)'); ylabel('Pe');

title('part 2(c): const 1 bound vs sim');




% part 2 d (constellation 2)
A2 =sqrt(1/2);
syms2=[A2 A2; A2 -A2; -A2 -A2; -A2 A2]';

Pe_sim2 = zeros(size(gam_db2));
Pe_bound2 =zeros(size(gam_db2));

fprintf('\npart 2 constellation 2\n');

for k=1:length(gam_db2)
    snr_lin =10^(gam_db2(k)/10);
    N0 =1/snr_lin;
    
    u_sum =0;
    for i=1:4

        for j=1:4
            if i~=j

                d = norm(syms2(:,i) -syms2(:,j));
                u_sum = u_sum +qfunc(d/sqrt(2*N0));
            end

        end
    end
    Pe_bound2(k) =u_sum / 4;
    
    tx_inds= randi([1 4], n_sim2, 1);
    tx_vecs =syms2(:,tx_inds);
    rx_vecs = tx_vecs + sqrt(N0/2)*randn(2, n_sim2);
    
    errs =0;
    for m=1:n_sim2
        [~, dec] = min(sum((syms2 -rx_vecs(:,m)).^2, 1));
        if dec ~= tx_inds(m),errs=errs+1; end
    end
    Pe_sim2(k) =errs/n_sim2;
    
    fprintf('SNR=%.2f, Bound=%.4e, Sim=%.4e\n', gam_db2(k),Pe_bound2(k), Pe_sim2(k));
end

figure;
semilogy(gam_db2, Pe_sim1,'b-o'); hold on;
semilogy(gam_db2,Pe_sim2, 'r-s');

grid on; legend('Const 1', 'Const 2');
xlabel('\gamma_s (dB)'); ylabel('Pe');
title('part 2(d): const comparison');

figure;
semilogy(gam_db2, Pe_bound2, 'k-', 'LineWidth',1.5); hold on;

semilogy(gam_db2, Pe_sim2, 'ro');
grid on; legend('Bound','Simulated');
xlabel('\gamma_s (dB)'); ylabel('Pe');

title('part 2(d): const 2 bound vs sim');


% part 2 e-f
map_nat =[0 0; 0 1; 1 0; 1 1];

map_gray= [0 0; 1 0; 1 1; 0 1];

gam_b_db =linspace(-4, 9, 10);

Ber_nat=zeros(size(gam_b_db));
Ber_gray= zeros(size(gam_b_db));

fprintf('\npart 2(f) BER\n');

for k=1:length(gam_b_db)
    snr_b_lin = 10^(gam_b_db(k)/10);
    snr_s_lin=2 * snr_b_lin; 

    N0 = 1/snr_s_lin;
    
    tx_inds =randi([1 4], n_sim2, 1);
    tx_vecs =syms2(:, tx_inds);

    rx_vecs =tx_vecs + sqrt(N0/2)*randn(2, n_sim2);
    
    dec_inds= zeros(n_sim2, 1);
    for m=1:n_sim2

        [~, dec_inds(m)]= min(sum((syms2- rx_vecs(:,m)).^2,1));
    end
    
    bit_err_n =sum(sum(map_nat(tx_inds,:)~=map_nat(dec_inds,:)));
    Ber_nat(k) = bit_err_n /(2*n_sim2);
    
    bit_err_g= sum(sum(map_gray(tx_inds,:)~=map_gray(dec_inds,:)));

    Ber_gray(k)= bit_err_g/ (2*n_sim2);
    
    fprintf('Eb/N0=%.2f, Nat=%.4e, Gray=%.4e\n',gam_b_db(k), Ber_nat(k),Ber_gray(k));
end

figure;
semilogy(gam_b_db,Ber_nat, 'b-o', 'LineWidth',1.2);

grid on; xlabel('\gamma_b (dB)'); ylabel('BER');
title('part 2(e): Natural Coding BER');

figure;
semilogy(gam_b_db, Ber_nat,'b-o', 'LineWidth',1.2); hold on;

semilogy(gam_b_db, Ber_gray,'r--s', 'LineWidth',1.2);
grid on; legend('Natural', 'Gray');

xlabel('\gamma_b (dB)'); ylabel('BER');
title('part 2(f): BER comparison');